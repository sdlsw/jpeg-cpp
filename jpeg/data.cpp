export module jpeg:data;

import std;
import msg;

#define HUFFMAN_EXTEND true

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

// For JPEG, we only encounter 8x8 blocks for DCT/quantization/etc.
// The standard document offers a table, so just hardcode it.
const std::array<uint8_t, 64> zigzag_order_from_row_order = {
	0,   1,  5,  6, 14, 15, 27, 28,
	2,   4,  7, 13, 16, 26, 29, 42,
	3,   8, 12, 17, 25, 30, 41, 43,
	9,  11, 18, 24, 31, 40, 44, 53,
	10, 19, 23, 32, 39, 45, 52, 54,
	20, 22, 33, 38, 46, 51, 55, 60,
	21, 34, 37, 47, 50, 56, 59, 61,
	35, 36, 48, 49, 57, 58, 62, 63
};

// However, we can use the hardcoded table to generate the inverse one:
constexpr auto generate_order_inverse(const std::array<uint8_t, 64>& order) {
	std::array<uint8_t, 64> inverse;

	for (int i = 0; i < order.size(); i++) {
		inverse[order[i]] = i;
	}

	return std::move(inverse);
}

const std::array<uint8_t, 64> row_order_from_zigzag_order = generate_order_inverse(zigzag_order_from_row_order);

export namespace jpeg {
enum class StartOfFrameType : uint8_t {
	baseline_seq = 0,
	extended_seq = 1,
	progressive =  2,
	lossless =     3 // always sequential
};

std::map<StartOfFrameType, std::string> str_from_sof_type {
	{StartOfFrameType::baseline_seq, "BASELINE_SEQUENTIAL"},
	{StartOfFrameType::extended_seq, "EXTENDED_SEQUENTIAL"},
	{StartOfFrameType::progressive,  "PROGRESSIVE"},
	{StartOfFrameType::lossless,     "LOSSLESS"}
};

enum class CodingType : uint8_t {
	huffman =    0,
	arithmetic = 1
};

std::map<CodingType, std::string> str_from_coding_type {
	{CodingType::huffman,    "HUFFMAN"},
	{CodingType::arithmetic, "ARITHMETIC"}
};

struct StartOfFrameInfo {
	StartOfFrameType type = StartOfFrameType::baseline_seq;
	CodingType coding = CodingType::huffman;
	bool is_differential = false;

	StartOfFrameInfo() = default;
	StartOfFrameInfo(StartOfFrameType sof_t, CodingType ct, bool diff) :
		type{sof_t}, coding{ct}, is_differential{diff} {}

	// From marker char. In this case, the enum definitions
	// we're using line up with the values given in the standard,
	// so just use bit magic to get them.
	StartOfFrameInfo(unsigned char m) :
		type{StartOfFrameType(m & 0x03)},
		coding{CodingType((m & 0x08) >> 3)},
		is_differential{bool(m & 0x04)} {}

	bool is_baseline() const {
		return type == StartOfFrameType::baseline_seq;
	}

	uint8_t marker_value() const {
		return (
			0xC0 |
			((uint8_t) type) |
			(((uint8_t) coding) << 3) |
			(((uint8_t) is_differential) << 2)
		);
	}

	explicit operator std::string() const {
		return std::format(
			"SOF {} | {}DIFFERENTIAL | {}-CODED",
			str_from_sof_type[type],
			(is_differential) ? "" : "NON-",
			str_from_coding_type[coding]
		);
	}
};

// See Table B.1 - Marker code assignments
enum class MarkerSpecial : uint8_t {
	tem_temporary_use =                   0x01,

	dht_define_huff_table =               0xC4,
	jpg_reserved =                        0xC8,
	dac_define_arith_coding =             0xCC,

	// Other markers
	soi_start_of_img =                    0xD8,
	eoi_end_of_img =                      0xD9,
	sos_start_of_scan =                   0xDA,
	dqt_define_quant_table =              0xDB,
	dnl_define_num_lines =                0xDC,
	dri_define_rst_interval =             0xDD,
	dhp_define_hierarchical_progression = 0xDE,
	exp_expand_reference_components =     0xDF,
	com_comment =                         0xFE
};

std::map<MarkerSpecial, std::string> symbol_from_special_marker {
	{MarkerSpecial::tem_temporary_use,                   "TEM*"},
	{MarkerSpecial::dht_define_huff_table,               "DHT"},
	{MarkerSpecial::jpg_reserved,                        "JPG"},
	{MarkerSpecial::dac_define_arith_coding,             "DAC"},
	{MarkerSpecial::soi_start_of_img,                    "SOI*"},
	{MarkerSpecial::eoi_end_of_img,                      "EOI*"},
	{MarkerSpecial::sos_start_of_scan,                   "SOS"},
	{MarkerSpecial::dqt_define_quant_table,              "DQT"},
	{MarkerSpecial::dnl_define_num_lines,                "DNL"},
	{MarkerSpecial::dri_define_rst_interval,             "DRI"},
	{MarkerSpecial::dhp_define_hierarchical_progression, "DHP"},
	{MarkerSpecial::exp_expand_reference_components,     "EXP"},
	{MarkerSpecial::com_comment,                         "COM"}
};

// Class for encapsulating all the magic numbers/bitmasks involved in
// decoding the marker values.
struct Marker {
	// Note: Technically markers are two bytes, but the first
	// byte is always 0xFF and doesn't contain any interesting
	// information, so leave it out of our representiation.
	unsigned char marker;

	bool is_sof() const {
		// SOF Markers are 0b1100xxxx, except for some special
		// markers stuffed in the middle of the SOF range.
		return (
			(marker & 0xF0) == 0xC0 &&
			marker != (unsigned char)(MarkerSpecial::dht_define_huff_table) &&
			marker != (unsigned char)(MarkerSpecial::jpg_reserved) &&
			marker != (unsigned char)(MarkerSpecial::dac_define_arith_coding)
		);
	}

	StartOfFrameInfo parse_sof() const {
		// SOF markers are basically bitfields, so
		// delegate parsing to another struct...
		return StartOfFrameInfo(marker);
	}

	static Marker sof(const StartOfFrameInfo& info) {
		return Marker(info.marker_value());
	}

	bool is_rstn() const {
		// RSTn Markers are 0b11010xxx
		return (marker & 0xF8) == 0xD0;
	}

	// Get the 'n' of RSTn.
	unsigned char parse_rstn() const {
		return marker & 0x07;
	}

	static Marker rstn(uint8_t n) {
		return Marker(0xD0 | (0x07 & n));
	}

	bool is_appn() const {
		// APPn Markers are 0b1110xxxx
		return (marker & 0xF0) == 0xE0;
	}

	// Get the 'n' of APPn.
	unsigned char parse_appn() const {
		return marker & 0x0F;
	}

	static Marker appn(uint8_t n) {
		return Marker(0xE0 | (0x0F & n));
	}

	// Markers reserved specifically for JPG extensions.
	bool is_jpgn() const {
		return marker >= 0xF0 && marker <= 0xFD;
	}

	// Get the 'n' of JPGn.
	unsigned char parse_jpgn() const {
		return marker & 0x0F;
	}

	// RES markers should be ignored.
	bool is_res() const {
		return marker >= 0x02 && marker <= 0xBF;
	}

	bool is_special() const {
		return !(
			is_sof() ||
			is_rstn() ||
			is_appn() ||
			is_jpgn() ||
			is_res()
		);
	}

	MarkerSpecial parse_special() const {
		return MarkerSpecial { marker };
	}

	bool is_soi() const {
		return MarkerSpecial(marker) == MarkerSpecial::soi_start_of_img;
	}

	bool is_eoi() const {
		return MarkerSpecial(marker) == MarkerSpecial::eoi_end_of_img;
	}

	// Test if a marker is a valid value. The standard specifically marks
	// two values as "invalid". They are used to encode 0xFF in the entropy
	// coded segments.
	bool is_valid() const {
		return marker != 0x00 && marker != 0xFF;
	}

	// If this is true, then this marker does not precede a
	// marker segment.
	bool is_standalone() const {
		if (is_rstn()) return true;

		switch (marker) {
			case (unsigned char)(MarkerSpecial::soi_start_of_img):
			case (unsigned char)(MarkerSpecial::eoi_end_of_img):
			case (unsigned char)(MarkerSpecial::tem_temporary_use):
				return true;
			default:
				return false;
		}
	}

	operator int() const { return marker; };

	explicit operator std::string() const {
		if (!is_valid()) {
			return std::format("INVALID_MARKER 0x{:X}", int(marker));
		}

		if (is_sof()) {
			return std::string(parse_sof());
		} else if (is_rstn()) {
			return std::format("RST{}*", int(parse_rstn()));
		} else if (is_appn()) {
			return std::format("APP{}", int(parse_appn()));
		} else if (is_jpgn()) {
			return std::format("JPG{}", int(parse_jpgn()));
		} else if (is_res()) {
			return std::format("RESERVED 0x{:X}", int(marker));
		} else {
			return symbol_from_special_marker[parse_special()];
		}
	}
};

const std::string hexchars = "0123456789ABCDEF";
std::string string_from_byte_vec(const std::vector<uint8_t>& v) {
	std::string s;
	s.reserve(
		2 * v.size() + // data chars
		v.size() / 4   // spaces
	);

	int spacercount = -1;
	for (int i = 0; i < v.size(); i++) {
		uint8_t x = v[i];
		uint8_t high = (0xF0 & x) >> 4;
		uint8_t low = 0xF & x;

		if (spacercount == 3) {
			s.push_back(' ');
			spacercount = 0;
		} else {
			spacercount++;
		}

		s.push_back(hexchars[high]);
		s.push_back(hexchars[low]);
	}

	return std::move(s);
}

template<typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

template<arithmetic T, size_t W=8, size_t H=8>
class Matrix {
	// Every instance of Matrix must be friends with every other instance.
	template <arithmetic OTHER_T, size_t OTHER_W, size_t OTHER_H>
	friend class Matrix;

private:
	// FIXME big problem here - std::array isn't movable so all of the
	// arithmetic operations are copying entire matrices
	std::array<T, W*H> data {0};

	constexpr void bounds_check(size_t x, size_t y) const {
		if (x >= W) throw std::out_of_range("x coordinate out of range");
		if (y >= H) throw std::out_of_range("y coordinate out of range");
	}

	template<arithmetic RESULT_T, arithmetic OTHER_T, size_t OTHER_W>
	constexpr RESULT_T calc_matmul_val(
		const Matrix<OTHER_T, OTHER_W, W>& other,
		size_t x,
		size_t y
	) const {
		// for self, start at leftmost col of row y and move
		// horizontally
		size_t self_cursor = y*W;

		// for other, start at top row of col x and move vertically
		size_t other_cursor = x;

		RESULT_T sum = 0;
		for (int i = 0; i < W; i++) {
			sum += other.data[other_cursor] * data[self_cursor];
			self_cursor++;
			other_cursor += W;
		}

		return sum;
	}

public:
	// for external access of dimensions
	static const size_t width = W;
	static const size_t height = H;
	static const unsigned int str_float_precision = 5;

	Matrix() = default;
	Matrix(const std::initializer_list<T>& ilist) {
		if (ilist.size() != data.size()) {
			auto s = std::format(
				"Matrix: bad initializer list size {} (expected {})",
				ilist.size(),
				data.size()
			);
			throw std::invalid_argument(s);
		}

		auto iarray = ilist.begin();
		for (size_t i = 0; i < ilist.size(); i++) {
			data[i] = iarray[i];
		}
	}

	decltype(auto) operator[](this auto& self, size_t x, size_t y) {
		self.bounds_check(x, y);
		return self.data[y*W + x];
	}

	// Shortcut access 'operator' for zig-zag ordering. Used only for 8x8
	// matrices.
	decltype(auto) zz(this auto& self, size_t z) {
		if (W != 8 || H != 8) throw std::runtime_error("zz only valid on 8x8 matrix");
		if (z >= 64) throw std::out_of_range("zz: out of range");

		return self.data[row_order_from_zigzag_order[z]];
	}

	// Typical matrix multiplication. By default the number type of the
	// result will be the same as the left operand (A in A*B).
	// In order to manually specify the type of the resulting matrix, you
	// can specify with `A.operator*<desired_type>(B)`.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T, size_t OTHER_W>
	Matrix<RESULT_T, OTHER_W, H> operator*(const Matrix<OTHER_T, OTHER_W, W>& other) const {
		Matrix<RESULT_T, OTHER_W, H> result;

		for (size_t x = 0; x < OTHER_W; x++) {
			for (size_t y = 0; y < H; y++) {
				result[x, y] = calc_matmul_val<RESULT_T>(other, x, y);
			}
		}

		return std::move(result);
	}

	// Creates a transposed copy of this matrix.
	Matrix<T, H, W> transpose() const {
		Matrix<T, H, W> result;

		for (size_t x = 0; x < W; x++) {
			for (size_t y = 0; y < H; y++) {
				result[y, x] = (*this)[x, y];
			}
		}

		return std::move(result);
	}

	// Multiplies each point in this matrix by the corresponding point
	// in another of the same size.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T>
	Matrix<RESULT_T, W, H> ptwise_mul(const Matrix<OTHER_T, W, H>& other) const {
		Matrix<RESULT_T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			result.data[i] = this->data[i] * other.data[i];
		}

		return std::move(result);
	}

	// Same as ptwise_mul, but for division.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T>
	Matrix<RESULT_T, W, H> ptwise_div(const Matrix<OTHER_T, W, H>& other) const {
		Matrix<RESULT_T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			result.data[i] = this->data[i] / other.data[i];
		}

		return std::move(result);
	}

	T max() const {
		T m = 0;
		for (size_t i = 0; i < data.size(); i++) {
			if (m < data[i]) m = data[i];
		}

		return m;
	}

	template<std::floating_point ST>
	size_t calc_field_width(this const Matrix<ST, W, H>& self) {
		return std::format("{:.{}f}", self.max(), str_float_precision).size() + 2;
	}

	template<std::integral ST>
	size_t calc_field_width(this const Matrix<ST, W, H>& self) {
		T m = self.max();
		if (m == 0) return 1;

		return (size_t)std::log10(std::abs(self.max())) + 1;
	}

	template<std::floating_point ST>
	static std::string fmt_item(ST val, size_t field_width) {
		return std::format("{:{}.{}f}", val, field_width, str_float_precision);
	}

	template<std::integral ST>
	static std::string fmt_item(ST val, size_t field_width) {
		return std::format("{:{}d}", val, field_width);
	}

	explicit operator std::string() const {
		size_t field_width = calc_field_width();

		std::stringstream out;

		unsigned int col = 0;
		for (T val : data) {
			out << fmt_item(val, field_width) << " ";
			
			col++;
			if (col >= W) {
				col = 0;
				out << std::endl;
			}
		}

		return out.str();
	}
};

auto calc_dct_mat() {
	Matrix<double, 8, 8> mat;

	double rt2 = std::sqrt(2);
	for (size_t y = 0; y < mat.height; y++) {
		for (size_t x = 0; x < mat.width; x++) {
			double c = (x == 0) ? (1/(2*rt2)) : 0.5;
			double cos_arg = (std::numbers::pi * x * (2*y + 1)) / 16;

			mat[x, y] = c * std::cos(cos_arg);
		}
	}

	return mat;
}

const Matrix<double, 8, 8> dct_mat = calc_dct_mat();
const Matrix<double, 8, 8> dct_mat_t = dct_mat.transpose();

// For a standard JPEG impl, we have no idea what application segments
// contain, so leave it as unparsed data.
struct AppSegment {
	unsigned int n; // APPn
	std::vector<uint8_t> data;

	explicit operator std::string() const {
		return std::format("APP{}: {}", n, string_from_byte_vec(data));
	}
};

struct QuantizationTable {
	bool set = false;
	Matrix<uint16_t, 8, 8> data;

	QuantizationTable() = default;
	QuantizationTable(const std::initializer_list<uint16_t>& ilist) : data{ilist}, set{true} {}

	uint8_t precision() const {
		if (data.max() <= 255) {
			return 0;
		} else {
			return 1;
		}
	}

	explicit operator std::string() const {
		if (!set) return "not set";
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

struct HuffmanCode {
	uint16_t code = 0; // huffman code
	uint8_t bits = 0;  // number of bits for this code
	uint8_t value = 0; // value of the code

	// Show the code as a binary string.
	std::string code_str() const {
		if (bits == 0) return "BITS=0";

		auto out = std::string(bits, '\0');
		int c = code;
		for (int b = bits - 1; b >= 0; b--) {
			out[b] = bool(c & 0x01) ? '1' : '0';
			c = c >> 1;
		}

		return std::move(out);
	}

	explicit operator std::string() const {
		return std::format("{:2X}: {}", value, code_str());
	}
};

struct HuffmanTable {
	bool set = false;

	// This vector contains the three tables HUFFSIZE, HUFFCODE and HUFFVAL,
	// as specified in section C.
	std::vector<HuffmanCode> codes;

	// Tables for increased lookup speed.
	std::array<uint16_t, 16> size_ptrs {0};
	std::array<uint8_t, 16> size_amts {0};

	// Value lookup table. A little big, but faster.
	std::array<HuffmanCode*, 255> value_ptrs;

	HuffmanTable() = default;
	HuffmanTable(const std::initializer_list<HuffmanCode> ilist) {
		codes.insert(codes.end(), ilist.begin(), ilist.end());
		populate_codes();
		set = true;
	}

	// Generate the table HUFFCODE, as specified in figure C.2
	// Requires HUFFSIZE to be populated.
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

	const HuffmanCode& lookup_value(uint8_t value) const {
		if (!set) {
			throw std::runtime_error("Cannot lookup value in unset table");
		}

		return *value_ptrs[value];
	}

	int16_t lookup_code(uint16_t lookup_code, uint8_t bits) const {
		if (!set) {
			throw std::runtime_error("Cannot lookup code in unset table");
		}

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
		if (!set) {
			return "not set";
		}

		if (HUFFMAN_EXTEND) {
			std::stringstream out;
			out << std::format("VALUES...\n");

			for (auto& code : codes) {
				out << std::string(code) << '\n';
			}

			return out.str();
		} else {
			return std::format("{} VALUES OMITTED", codes.size());
		}
	}

	// A bit hacky... but most Jpeg files I've seen in the wild all use
	// the tables suggested by the JPEG standard. To save some time, this
	// function will spit out the codes formatted as an initializer list
	// that can be used in C++ code.
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

struct ArithmeticTable {
	bool set = false;
	TableClass _class = TableClass::dc;
	uint8_t value = 0;
};

struct ScanComponentParams {
	uint8_t component_selector = 0;
	uint8_t dc_entropy_coding_selector = 0;
	uint8_t ac_entropy_coding_selector = 0;

	explicit operator std::string() const {
		return std::format(
			"ID={} | DCSEL={} | ACSEL={}",
			component_selector,
			dc_entropy_coding_selector,
			ac_entropy_coding_selector
		);
	}
};

struct Scan {
	uint8_t num_components = 0;
	std::vector<ScanComponentParams> component_params;

	uint8_t start_spectral_sel = 0;
	uint8_t end_spectral_sel = 63;
	uint8_t succ_approx_high = 0;
	uint8_t succ_approx_low = 0;

	uint16_t restart_interval = 0;

	// Entropy-coded data as it is encoded in the file, with RST markers
	// removed in favor of breaking up the intervals into sub-vectors.
	std::vector<std::vector<uint8_t>> entropy_coded_data;

	size_t total_bytes_data() const {
		size_t sum = 0;

		for (auto& segment : entropy_coded_data) {
			sum += segment.size();
		}

		return sum;
	}

	explicit operator std::string() const {
		std::stringstream out;

		out << std::format(
			"Ss={} | Se={} | Ah={} | Al={} | RSTi={}\n",
			start_spectral_sel,
			end_spectral_sel,
			succ_approx_high,
			succ_approx_low,
			restart_interval
		);

		for (auto& comp : component_params) {
			out << std::format(
				"COMPONENT {}\n",
				std::string(comp)
			);
		}

		out << std::format("{} SEGMENT(S)\n", entropy_coded_data.size());
		out << std::format("{} BYTES IMAGE DATA", total_bytes_data());

		return out.str();
	}
};

struct AcCoeff {
	uint8_t r = 0;
	uint8_t s = 0;
	int16_t value = 0;

	bool is_eob() const { return r == 0 && s == 0; }
};

struct FrameComponentParams {
	uint8_t identifier = 0;
	uint8_t horizontal_sampling_factor = 0;
	uint8_t vertical_sampling_factor = 0;
	uint8_t qtable_selector = 0;

	explicit operator std::string() const {
		return std::format(
			"ID={} | H={} | V={} | QTABLE={}",
			identifier,
			horizontal_sampling_factor,
			vertical_sampling_factor,
			qtable_selector
		);
	}
};

struct Frame {
	StartOfFrameInfo sof_info;

	uint8_t sample_precision = 8;
	uint16_t num_lines = 0;
	uint16_t samples_per_line = 0;
	uint8_t num_components = 0;

	std::vector<FrameComponentParams> component_params;
	std::vector<Scan> scans;

	int h_max() const {
		int max = 0;

		for (auto& params : component_params) {
			if (params.horizontal_sampling_factor > max) {
				max = params.horizontal_sampling_factor;
			}
		}

		return max;
	}

	int v_max() const {
		int max = 0;

		for (auto& params : component_params) {
			if (params.vertical_sampling_factor > max) {
				max = params.vertical_sampling_factor;
			}
		}

		return max;
	}

	int mcu_length() const {
		int out = 0;

		for (auto& params : component_params) {
			out += params.vertical_sampling_factor *
				params.horizontal_sampling_factor;
		}

		return out;
	}

	explicit operator std::string() const {
		std::stringstream out;
		out << std::format(
			"MARK: {}\nP={} | X={} | Y={}\n",
			std::string(sof_info),
			sample_precision,
			samples_per_line,
			num_lines
		);

		for (auto& comp : component_params) {
			out << std::format(
				"COMPONENT {}\n",
				std::string(comp)
			);
		}

		int count = 0;
		for (auto& scan : scans) {
			out << std::format(
				"SCAN{}:\n{}\n",
				count,
				std::string(scan)
			);
			count++;
		}

		return out.str();
	}
};

// All data needed to decode a JPEG image, including tables, entropy-coded
// data, etc.
class CompressedJpegData {
private:
	// Standard supports up to 4 slots for each type of table.
	static const uint8_t num_tables = 4;

	bool valid = false;
	std::string _source_file {""};

	std::array<QuantizationTable, num_tables> _q_tables;
	std::array<HuffmanTable, num_tables> _ac_huff_tables;
	std::array<HuffmanTable, num_tables> _dc_huff_tables;
	std::array<ArithmeticTable, num_tables> arith_tables;

	std::vector<AppSegment> _app_segments;
	std::vector<std::vector<uint8_t>> _comments;
	std::vector<Frame> _frames;
public:
	bool is_valid() {
		return valid;
	}

	void set_valid() {
		valid = true;
	}

	auto& app_segments() { return _app_segments; }
	auto& comments() { return _comments; }
	auto& q_tables() { return _q_tables; }
	auto& huff_tables(TableClass _class) {
		if (_class == TableClass::dc) {
			return _dc_huff_tables;
		} else {
			return _ac_huff_tables;
		}
	}
	auto& frames() { return _frames; }

	const auto& q_tables() const { return _q_tables; }
	const auto& huff_tables(TableClass _class) const {
		if (_class == TableClass::dc) {
			return _dc_huff_tables;
		} else {
			return _ac_huff_tables;
		}
	}

	const auto& frames() const { return _frames; }

	explicit operator std::string() const {
		std::stringstream out;

		for (auto& seg : _app_segments) {
			out << std::string(seg) << '\n';
		}

		for (auto& com : _comments) {
			out << "COM " << string_from_byte_vec(com) << '\n';
		}

		for (int i = 0; i < num_tables; i++) {
			out << std::format("QTABLE{}:\n", i) << std::string(_q_tables[i]) << '\n';
		}

		for (int i = 0; i < num_tables; i++) {
			out << std::format("DC HUFFTABLE{}: ", i) << std::string(_dc_huff_tables[i]) << '\n';
			out << std::format("AC HUFFTABLE{}: ", i) << std::string(_ac_huff_tables[i]) << '\n';
		}

		if (_frames.size() == 0) {
			out << "NO FRAMES\n";
		} else {
			int count = 0;
			for (auto& f : _frames) {
				out << std::format(
					"FRAME{}:\n{}\n",
					count,
					std::string(f)
				);
			}
		}

		return out.str();
	}
};

class RawComponent;

const unsigned int block_size = 8;

/*
 * Helper class for accessing one block of image data.
 */
class BlockView {
private:
	int _mcu_x;
	int _mcu_y;

	int _block_x;
	int _block_y;

	int _x_offset;
	int _y_offset;

	RawComponent* _component;

	static void coord_check(int x, int y);
	bool coord_advance(
		int& x,
		int& y,
		const int& x_bound,
		const int& y_bound
	);
	void recalc_x_offset();
	void recalc_y_offset();
public:
	BlockView() = delete;
	BlockView(
		int mcu_x,
		int mcu_y,
		int block_x,
		int block_y,
		RawComponent& comp
	);
	BlockView(RawComponent& comp) : BlockView(0, 0, 0, 0, comp) {}

	const uint8_t& operator[](int x, int y) const;
	uint8_t& operator[](int x, int y);

	// Allow read-only access to the underlying component for access
	// to dimension information.
	const RawComponent& component() const { return *_component; }

	int mcu_x() { return _mcu_x; }
	int mcu_y() { return _mcu_y; }

	bool mcu_advance();
	bool block_advance();

	explicit operator std::string() const {
		return std::format(
			"mcux={} mcuy={} blockx={} blocky={} offset=<{}, {}>",
			_mcu_x, _mcu_y, _block_x, _block_y, _x_offset, _y_offset
		);
	}
};

class RawComponent {
private:
	// Width and height of valid component data.
	int _x_pixels;
	int _y_pixels;

	// Width and height of backing storage in MCUs (Minimum Coded Units).
	int _x_mcus;
	int _y_mcus;

	// Width and height of one MCU, in blocks.
	int _mcu_width;
	int _mcu_height;

	// Maximum bounds on the backing store. This may be a bit
	// larger than _x_pixels and/or _y_pixels if those values aren't
	// integer multiples of the MCU size in pixels. This could be
	// calculated every time from other values, but it's cached here
	// since we check these bounds on every sample access.
	int _x_size;
	int _y_size;

	std::vector<uint8_t> _data;

	static int calc_n_mcus(int samples, int mcu_len_in_blocks) {
		auto [q, r] = std::div(samples, block_size*mcu_len_in_blocks);
		return q + ((r == 0) ? 0 : 1);
	}

	int calc_x_size() const {
		return _x_mcus * _mcu_width * block_size;
	}

	int calc_y_size() const {
		return _y_mcus * _mcu_height * block_size;
	}

	int data_idx_from_xy(int x, int y) const {
		if (x >= _x_size) {
			throw std::out_of_range("RawComponent: x value out of range");
		}

		if (y >= _y_size) {
			throw std::out_of_range("RawComponent: y value out of range");
		}

		return y*_x_size + x;
	}
public:
	RawComponent() = delete;

	RawComponent(
		int x_pixels,
		int y_pixels,
		int mcu_width_blocks,
		int mcu_height_blocks,
		uint8_t default_value
	) {
		// All of the various calculated sizes are cached on
		// construction, since we're going to be using them a LOT.
		_x_pixels = x_pixels;
		_y_pixels = y_pixels;

		_mcu_width = mcu_width_blocks;
		_mcu_height = mcu_height_blocks;

		_x_mcus = calc_n_mcus(_x_pixels, mcu_width_blocks);
		_y_mcus = calc_n_mcus(_y_pixels, mcu_height_blocks);

		_x_size = calc_x_size();
		_y_size = calc_y_size();

		size_t vecsize = _x_size * _y_size;
		_data = std::vector<uint8_t>(vecsize, default_value);
	}

	RawComponent(
		int x_pixels,
		int y_pixels,
		int mcu_width_blocks,
		int mcu_height_blocks
	) : RawComponent(x_pixels, y_pixels, mcu_width_blocks, mcu_height_blocks, 0) {}

	int x_pixels() const { return _x_pixels; }
	int y_pixels() const { return _y_pixels; }

	int x_mcus() const { return _x_mcus; }
	int y_mcus() const { return _y_mcus; }

	int mcu_width() const { return _mcu_width; }
	int mcu_height() const { return _mcu_height; }

	int x_size() const { return _x_size; }
	int y_size() const { return _y_size; }

	// Direct sample access. Access is allowed to go slightly beyond
	// img_width/img_height, to simplify cases where the img size and
	// block size are misaligned.
	const uint8_t& operator[](int x, int y) const { return _data[data_idx_from_xy(x, y)]; }
	uint8_t& operator[](int x, int y) { return _data[data_idx_from_xy(x, y)]; }

	explicit operator std::string() const {
		return std::format(
			"xpixels={} ypixels={} mcu_width={} mcu_height={} xmcus={} ymcus={}",
			_x_pixels, _y_pixels, _mcu_width, _mcu_height, _x_mcus, _y_mcus
		);
	}
};

BlockView::BlockView(
	int mcu_x,
	int mcu_y,
	int block_x,
	int block_y,
	RawComponent& comp
) {
	if (mcu_x >= comp.x_mcus()) {
		throw std::out_of_range("mcu x value out of range");
	}

	if (mcu_y >= comp.y_mcus()) {
		throw std::out_of_range("mcu y value out of range");
	}

	if (block_x >= comp.mcu_width()) {
		throw std::out_of_range("block x value out of range");
	}

	if (block_y >= comp.mcu_height()) {
		throw std::out_of_range("block y value out of range");
	}

	_component = &comp;

	_mcu_x = mcu_x;
	_mcu_y = mcu_y;
	_block_x = block_x;
	_block_y = block_y;

	recalc_x_offset();
	recalc_y_offset();
}

void BlockView::recalc_x_offset() {
	_x_offset = (_mcu_x * _component->mcu_width() + _block_x) * block_size;
}

void BlockView::recalc_y_offset() {
	_y_offset = (_mcu_y * _component->mcu_height() + _block_y) * block_size;
}

void BlockView::coord_check(int x, int y) {
	if (x >= block_size) {
		throw std::out_of_range("BlockView: x value out of range");
	}

	if (y >= block_size) {
		throw std::out_of_range("BlockView: y value out of range");
	}
}

const uint8_t& BlockView::operator[](int x, int y) const {
	coord_check(x, y);
	return (*_component)[_x_offset + x, _y_offset + y];
}

uint8_t& BlockView::operator[](int x, int y) {
	coord_check(x, y);
	return (*_component)[_x_offset + x, _y_offset + y];
}

bool BlockView::coord_advance(
	int& x,
	int& y,
	const int& x_bound,
	const int& y_bound
) {
	//msg::debug("DECODE: coord_advance x={} y={}", x, y);
	//msg::debug("DECODE: x_bound={} y_bound={}", x_bound, y_bound);

	if (y >= y_bound || (y == y_bound - 1 && x >= x_bound)) {
		msg::debug("bound");
		return false;
	}

	bool recalc_y = false;

	x++;

	if (x >= x_bound) {
		x = 0;
		y++;
	}

	recalc_x_offset();
	recalc_y_offset();

	//msg::debug("DECODE: coord_advance AFTER x={} y={}", x, y);

	if (y >= y_bound || (y == y_bound - 1 && x >= x_bound)) {
		return false;
	}
	return true;
}

bool BlockView::mcu_advance() {
	_block_x = 0;
	_block_y = 0;

	return coord_advance(
		_mcu_x,
		_mcu_y,
		_component->x_mcus(),
		_component->y_mcus()
	);
}

bool BlockView::block_advance() {
	return coord_advance(
		_block_x,
		_block_y,
		_component->mcu_width(),
		_component->mcu_height()
	);
}
}
