export module jpeg:data.tables;

// data_tables.cpp:
// Defines tables used in JPEG compression, and containers for those tables.

import std;
import msg;
import :data.matrix;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

#define HUFFMAN_EXTEND false

export namespace jpeg {
// A generic container for various types of tables involved in JPEG
// compression. TableMaps have 4 "slots" into which tables may be installed.
// The state of these slots may be saved and loaded, and any tables
// installed to this map will be maintained even when overwritten.
template<typename TableType>
class TableMap {
private:
	static const inline unsigned int n_slots = 4;
	static const inline std::array<int8_t, n_slots> no_tables { -1, -1, -1, -1 };

	// The tables stored in this map.
	std::vector<TableType> tables;

	// A 4-tuple mapping a jpeg table destination value to an index
	// in `tables`. This value may not be accessed directly. Once all
	// desired tables are installed, `save_dest` must be called at least
	// once. The tables are then accessible through `get(0, dest_id)`.
	std::array<int8_t, n_slots> table_dest { -1, -1, -1, -1 };

	// Saved values of table_dest.
	std::vector<std::array<int8_t, n_slots>> table_dest_history;

	// Bounds check for `dest_id` parameters.
	static void check_dest_id(uint8_t dest_id) {
		if (dest_id >= n_slots) {
			throw std::out_of_range(std::format(
				"dest_id >= {}", n_slots
			));
		}
	}

	// Bounds check for `hist_id` parameters.
	void check_hist_id(size_t hist_id) const {
		if (hist_id >= table_dest_history.size()) {
			throw std::out_of_range("hist_id out of range");
		}
	}

	// Gets a history entry. If `hist_id` is negative, returns an empty
	// history entry.
	auto& hist_get_extended(int hist_id) const {
		if (hist_id < 0) {
			return no_tables;
		} else {
			return table_dest_history[hist_id];
		}
	}
public:
	const int hist_last() const {
		return static_cast<int>(table_dest_history.size()) - 1;
	}

	// Gets the table at the given `dest_id`, at `hist_id` point in
	// history.
	const TableType& get(size_t hist_id, uint8_t dest_id) const {
		check_hist_id(hist_id);
		check_dest_id(dest_id);

		int8_t table_idx = table_dest_history[hist_id][dest_id];

		if (table_idx < 0) {
			throw std::runtime_error(std::format(
				"No table installed for (hist, dest) = ({}, {})",
				hist_id, dest_id
			));
		}

		if (table_idx >= tables.size()) {
			throw std::out_of_range(std::format(
				"Table index {} out of valid range!", table_idx
			));
		}

		return tables[table_idx];
	}

	// Gets a vector of const ptrs to every new table that was introduced
	// between `hist_id - 1` and `hist_id`. If `hist_id` is 0, then this
	// just returns every set table for the first history entry.
	auto get_differential(size_t hist_id) const {
		check_hist_id(hist_id);

		std::vector<std::tuple<const TableType*, unsigned int>> v;

		auto& prev_dest = hist_get_extended(static_cast<int>(hist_id) - 1);
		auto& cur_dest = table_dest_history[hist_id];

		for (unsigned int i = 0; i < n_slots; i++) {
			if (cur_dest[i] != prev_dest[i] && cur_dest[i] >= 0) {
				const auto& tbl = tables[cur_dest[i]];
				v.emplace(v.end(), &tbl, i);
			}
		}

		return v;
	}

	// Checks if this TableMap has any tables defined in table_dest.
	bool any_defined_in_dest() {
		for (const auto& v : table_dest) {
			if (v != -1) return true;
		}

		return false;
	}

	// Saves the current state of table_dest.
	void save_dest() {
		table_dest_history.push_back(table_dest);
	}

	// Installs a new table to this TableMap. Previously installed tables
	// will be maintained.
	template<typename... Args>
	void install_table(uint8_t dest_id, Args&&... args) {
		check_dest_id(dest_id);
		if (table_dest[dest_id] >= 0) {
			msg::debug("Replacing table at dest {}", dest_id);
		}

		// Construct new table in place and add a reference to it.
		tables.emplace(tables.end(), std::forward<Args>(args)...);
		table_dest[dest_id] = tables.size() - 1;
	}

	explicit operator std::string() const {
		std::stringstream out;

		out << "TABLES:\n";

		if (tables.size() == 0) {
			out << "None set.\n";
		} else {
			for (unsigned int i = 0; i < tables.size(); i++) {
				const auto& tbl = tables[i];
				out << std::format("TABLE{}:\n{}\n", i, std::string(tbl));
			}
		}

		out << "HISTMAP:\n";
		if (table_dest_history.size() == 0) {
			out << "Empty.\n";
		} else {
			for (unsigned int i = 0; i < table_dest_history.size(); i++) {
				const auto& hist = table_dest_history[i];
				out << std::format(
					"{}: ({}, {}, {}, {})\n",
					i, hist[0], hist[1], hist[2], hist[3]
				);
			}
		}

		out << '\n';

		return out.str();
	}
};

struct QuantizationTable {
	Matrix<uint16_t, 8, 8> data;

	QuantizationTable() = default;
	QuantizationTable(const std::initializer_list<uint16_t>& ilist) : data{ilist} {}
	QuantizationTable(const Matrix<uint16_t, 8, 8>& mat) : data{mat} {}
	QuantizationTable(Matrix<uint16_t, 8, 8>&& mat) : data{std::move(mat)} {}

	// Determines the precision value for this table. 0 indicates an 8-bit
	// table, and 1 indicates a 16-bit table. This value is only used when
	// reading/writing the table from/to a JPEG file.
	uint8_t precision() const {
		if (data.max() <= 255) {
			return 0;
		} else {
			return 1;
		}
	}

	explicit operator std::string() const {
		return std::string(data);
	}
};

enum class TableClass : uint8_t {
	dc = 0, ac = 1
};

const std::map<TableClass, std::string> str_from_table_class = {
	{TableClass::dc, "DC"},
	{TableClass::ac, "AC"}
};

// Struct representing a single Huffman code.
struct HuffmanCode {
	uint16_t code = 0; // huffman code
	uint8_t bits = 0;  // number of bits for this code
	uint8_t value = 0; // value of the code

	// Format the code as a binary string. This differs from the typical
	// formatting; while `code` is 16 bits, only the first `bits` bits
	// are important, so only those are included in the returned string.
	std::string code_str() const {
		if (bits == 0) return "BITS=0";

		auto out = std::string(bits, '\0');
		int c = code;
		for (int b = bits - 1; b >= 0; --b) {
			out[b] = bool(c & 0x01) ? '1' : '0';
			c = c >> 1;
		}

		return std::move(out);
	}

	explicit operator std::string() const {
		return std::format("{:2X}: {}", value, code_str());
	}
};

// Default value for HuffmanTable::lookup_value
const HuffmanCode huffman_zero { 0, 0, 0 };

// Struct representing a Huffman code table.
struct HuffmanTable {
	// This vector contains the three tables HUFFSIZE, HUFFCODE and HUFFVAL,
	// as specified in section C.
	std::vector<HuffmanCode> codes;

	// Tables for increased code lookup speed.
	std::array<uint16_t, 16> size_ptrs {0};
	std::array<uint8_t, 16> size_amts {0};

	// Value lookup table. A little big, but faster.
	std::array<const HuffmanCode*, 255> value_ptrs { &huffman_zero };

	HuffmanTable() = default;
	HuffmanTable(const std::initializer_list<HuffmanCode> ilist) {
		codes.insert(codes.end(), ilist.begin(), ilist.end());
		populate_codes();
	}

	// Generates the table HUFFCODE, and populates additional lookup
	// tables. Requires HUFFSIZE and HUFFVAL to be populated.
	void populate_codes() {
		if (codes.empty()) {
			msg::warn("HuffmanTable: Attempt to populate_codes() on empty table");
			return;
		}

		uint16_t code = 0;
		int cur_size = codes[0].bits;

		for (int k = 0; k < codes.size(); k++) {
			HuffmanCode& entry = codes[k];

			value_ptrs[entry.value] = &entry;

			entry.code = code;
			code++;

			size_amts[cur_size-1]++;

			if (k == codes.size() - 1) continue; // guard against out of range error

			// lookahead and handle potential increase in code size
			int next_size = codes[k+1].bits;
			if (next_size == cur_size) continue;

			size_ptrs[next_size-1] = k+1;

			// next_size will be larger than cur_size, since codes
			// are organized in order of code size.
			code = code << (next_size - cur_size);
			cur_size = next_size;
		}
	}

	// Looks up the Huffman code for a given value. If this value isn't
	// present in the table, then a reference to a "zero" Huffman code is
	// returned. This may be checked by testing if `bits` is 0:
	//
	// if (table.lookup_value(val).bits == 0) {
	//      // handle missing value, valid huffman codes never have 0 bits
	// }
	const HuffmanCode& lookup_value(uint8_t value) const {
		return *value_ptrs[value];
	}

	// Looks up the value of a Huffman code. If the code could not be
	// found, returns -1. If the code is found, its value will be an 8-bit
	// unsigned integer fitted in the lower half of the returned int16_t.
	int16_t lookup_code(uint16_t lookup_code, uint8_t bits) const {
		// Optimization: `codes` is sorted by code length, smallest
		// first. We take advantage of this by iterating only
		// over the codes that are N bits long, eliminating the
		// need to check the size manually.
		//
		// In addition, if there are no codes of size N, then start=0
		// and end=0; this skips the loop altogether.
		//
		// In practice this results in a significant speedup, often up 
		// to 2x overall execution speed for decoding compared to the
		// previous implementation.
		uint16_t start = size_ptrs[bits-1];
		uint16_t end = start + size_amts[bits-1];
		for (int k = start; k < end; k++) {
			const auto& code = codes[k];
			if (lookup_code == code.code) return code.value;
		}

		return -1;
	}

	explicit operator std::string() const {
		if (HUFFMAN_EXTEND) {
			std::stringstream out;
			out << std::format("VALUES...\n");

			for (const auto& code : codes) {
				out << std::string(code) << '\n';
			}

			return out.str();
		} else {
			return std::format("{} VALUES OMITTED", codes.size());
		}
	}

	// A bit hacky... but most baseline process JPEG files I've seen in the
	// wild all use the tables suggested by the JPEG standard. To save some
	// time, this function will spit out the codes formatted as an
	// initializer list that can be used in C++ code.
	// TODO Implement custom tables?
	std::string as_init_list() const {
		std::stringstream out;

		out << "{\n";

		for (const auto& code : codes) {
			// Note: Codes are set to zero because our
			// algorithm can just regenerate them anyway.
			out << std::format("\t{{ 0, {}, 0x{:X} }},\n", code.bits, code.value);
		}

		out << "}";

		return out.str();
	}
};
}
