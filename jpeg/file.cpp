export module jpeg:file;

import std;
import msg;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;
using std::uint32_t;
using std::int32_t;

import :data;

export namespace jpeg {
class JpegFileError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

class JpegSegmentWriter {
private:
	std::fstream* file;
	std::vector<uint8_t> data;

	// uint16_t to ordered pair of bytes, big endian.
	std::tuple<uint8_t, uint8_t> uint16_be(uint16_t x) {
		uint8_t h = (x & 0xFF00) >> 8;
		uint8_t l = x & 0x00FF;

		return { h, l };
	}

	void put_length() {
		auto [h, l] = uint16_be(data.size() + 2);
		file->put(h);
		file->put(l);
	}
public:
	JpegSegmentWriter(std::fstream& file) {
		this->file = &file;
	}

	void commit() {
		put_length();
		for (auto& x : data) {
			file->put(x);
		}
	}

	void write_byte(uint8_t b) {
		data.push_back(b);
	}

	void write_uint16(uint16_t x) {
		auto [h, l] = uint16_be(x);
		data.push_back(h);
		data.push_back(l);
	}

	void write_2x_uint4(uint8_t h, uint8_t l) {
		l &= 0x0F;

		data.push_back((h << 4) | l);
	}

	void write_str(const std::string& s) {
		for (auto& c : s) {
			data.push_back(c);
		}

		data.push_back('\0');
	}
};

class JpegFile {
private:
	std::fstream file;
	std::ios_base::openmode mode;
	uint16_t restart_interval = 0;

	// Ensures a parameter is under param_limit, throwing an exception
	// if not. The majority of JPEG parameters are unsigned and have a
	// lower bound of 0, so it is sufficient to check only an upper bound
	// in many cases.
	template<std::integral T>
	static void verify_param_under(const std::string& param_name, T param, uint16_t param_limit) {
		if (param >= param_limit) {
			throw JpegFileError(std::format(
				"{0} parameter invalid, {0}={1}>={2}",
				param_name,
				param,
				param_limit
			));
		}
	}

	uint8_t read_byte() {
		return (uint8_t) file.get();
	}

	void read_until_byte(uint8_t b) {
		uint8_t rb = 0;
		while (rb != b) rb = read_byte();
	}

	uint16_t read_uint16() {
		// File is big-endian
		uint16_t high = read_byte();
		uint16_t low = read_byte();

		return (high << 8) | low;
	}

	void write_marker(const Marker& m) {
		file.put(0xFF).put(m);
	}

	// Shortcut, allow writing a special marker constant directly
	void write_marker(MarkerSpecial m) {
		file.put(0xFF).put((uint8_t)m);
	}

	// Read next marker, placing the cursor after the two byte marker value.
	// If for some reason no valid marker could be read, an invalid marker
	// is returned.
	Marker read_next_marker() {
		int marker_byte = 0;

		while(true) {
			read_until_byte(0xFF);
			marker_byte = read_byte();

			if (marker_byte != 0x00 && marker_byte != 0xFF) {
				return Marker(marker_byte);
			}
		}

		return Marker(0x00);
	}

	bool read_expect_soi() {
		Marker mark = read_next_marker();
		if (!mark.is_valid()) return false;
		return mark.is_soi();
	}

	// Reads segment size. In the file, this value includes the size of the
	// length field itself (two bytes). This function returns the size
	// of everything *after* the length field, and sets the cursor to
	// the byte after it.
	uint16_t read_segment_size() {
		return read_uint16() - 2;
	}

	// Read a byte from the file, that contains two four-bit values.
	std::tuple<uint8_t, uint8_t> read_2x_uint4() {
		uint8_t b = read_byte();
		return { (b & 0xF0) >> 4, (b & 0x0F) };
	}

	// Reads a segment consisting of only a length field followed by
	// arbitrary binary data.
	void read_opaque_segment_to(std::vector<uint8_t>& dest) {
		uint16_t length = read_segment_size();
		msg::debug("READ_JPEG: reading opaque segment size {}", length);

		dest.reserve(length);
		while (length > 0) {
			dest.push_back(read_byte());
			length--;
		}
	}

	void read_marker_segment_to(CompressedJpegData& j, const Marker& marker) {
		if (marker.is_appn()) {
			read_appn_to(j, marker.parse_appn());
		} else if (marker.is_sof()) {
			read_sof_to(j, marker.parse_sof());
		} else if (marker.is_rstn()) {
			read_rstn_to(j, marker.parse_rstn());
		} else if (marker.is_res() || marker.is_jpgn()) {
			msg::warn("Skipping reserved marker {}", std::string(marker));
		} else {
			read_special_marker_segment_to(j, marker);
		}
	}

	void read_restart_interval() {
		msg::debug("READ_JPEG: process DRI");

		uint16_t segment_length = read_segment_size();

		if (segment_length != 2) {
			throw JpegFileError(std::format(
				"DRI: Bad segment length {}", segment_length
			));
		}

		// Save restart interval for later. It will be applied to any
		// following parsed SOS segments.
		restart_interval = read_uint16();
	}

	void read_appn_to(CompressedJpegData& j, int n) {
		msg::debug("READ_JPEG: process APP{}", n);

		AppSegment appseg;
		appseg.n = n;
		read_opaque_segment_to(appseg.data);
		j.app_segments().push_back(std::move(appseg));
	}

	void read_comment_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process COM");

		std::vector<uint8_t> dest;
		read_opaque_segment_to(dest);
		j.comments().push_back(std::move(dest));
	}

	void read_sof_to(CompressedJpegData& j, const StartOfFrameInfo& i) {
		msg::debug("READ_JPEG: process {}", std::string(i));

		uint16_t segment_length = read_segment_size();

		Frame frame;
		frame.sof_info = i;

		frame.sample_precision = read_byte();
		frame.num_lines = read_uint16();
		frame.samples_per_line = read_uint16();
		frame.num_components = read_byte();

		segment_length -= 6;

		while (segment_length > 0) {
			FrameComponentParams params;

			params.identifier = read_byte();

			auto [h, v] = read_2x_uint4();
			params.horizontal_sampling_factor = h;
			params.vertical_sampling_factor = v;

			verify_param_under("(Hi-1)", h - 1, 4);
			verify_param_under("(Vi-1)", v - 1, 4);

			params.qtable_selector = read_byte();

			frame.component_params.push_back(params);
			segment_length -= 3;
		}

		if (segment_length != 0) {
			throw JpegFileError("unexpected final value for length while reading SOF, abort");
		}

		j.frames().push_back(std::move(frame));
	}

	void read_rstn_to(CompressedJpegData& j, int n) {
		msg::debug("READ_JPEG: process RST{}", n);
	}

	void read_qtable_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process DQT");

		uint16_t length = read_segment_size();

		while (length > 0) {
			auto [pq_precision, tq_destination] = read_2x_uint4();

			verify_param_under("Tq", tq_destination, 4);
			verify_param_under("Pq", pq_precision, 2);

			QuantizationTable& qtable = j.q_tables()[tq_destination];
			msg::debug("READ_JPEG: install qtable to dest {}", tq_destination);

			if (qtable.set) {
				msg::warn("READ_JPEG: replacing existing table!");
			}

			for (int i = 0; i < 64; i++) {
				uint16_t value = 0;

				if (pq_precision == 0) {
					//8-bit precision
					value = read_byte();
				} else {
					//16-bit precision (pq = 1)
					value = read_uint16();
				}

				auto& valref = qtable.data.zz(i);
				valref = value;
			}

			qtable.set = true;
			length -= 65 + 64*pq_precision;
		}

		if (length != 0) {
			throw JpegFileError("unexpected final value for length while reading qtables, abort");
		}
	}

	void read_hufftable_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process DHT");

		uint16_t segment_length = read_segment_size();

		while (segment_length > 0) {
			auto [tc_table_class, th_destination] = read_2x_uint4();

			verify_param_under("Th", th_destination, 4);
			verify_param_under("Tc", tc_table_class, 2);

			segment_length--;

			auto _class = TableClass(tc_table_class);
			HuffmanTable& hufftable = j.huff_tables(_class)[th_destination];
			msg::debug(
				"READ_JPEG: install {} hufftable to dest {}",
				str_from_table_class.at(_class),
				th_destination
			);

			if (hufftable.set) {
				msg::warn("READ_JPEG: replacing existing table!");
			}

			// Read row lengths first.
			std::vector<uint8_t> row_lengths;
			for (int i = 0; i < 16; i++) {
				uint8_t len = read_byte();

				row_lengths.push_back(len);
				segment_length--;
			}

			// Using the length declarations, read huffman table.
			// Specifically, this is the HUFFVAL table.
			// HUFFSIZE is also populated during this step.
			for (int i = 0; i < 16; i++) {
				auto huff_row_length = row_lengths[i];
				while (huff_row_length > 0) {
					uint8_t val = read_byte();

					HuffmanCode entry;
					entry.value = val;
					entry.bits = i+1;

					hufftable.codes.push_back(entry);
					huff_row_length--;
					segment_length--;
				}
			}

			// Populate HUFFCODE table
			hufftable.populate_codes();
			hufftable.set = true;
		}

		if (segment_length != 0) {
			throw JpegFileError(
				"unexpected final value for length while reading hufftables, abort"
			);
		}
	}

	void read_arithtable_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process DAC");
	}

	void read_scan_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process SOS");

		if (j.frames().size() == 0) {
			throw JpegFileError("encountered SOS marker before SOF");
		}

		uint16_t segment_length = read_segment_size();

		Scan scan;
		scan.restart_interval = restart_interval; // Defined earlier by DRI
		scan.num_components = read_byte();
		segment_length--;

		verify_param_under("(Ns-1)", scan.num_components - 1, 4);

		for (int i = 0; i < scan.num_components; i++) {
			ScanComponentParams params;

			params.component_selector = read_byte();

			auto [dc_sel, ac_sel] = read_2x_uint4();
			params.dc_entropy_coding_selector = dc_sel;
			params.ac_entropy_coding_selector = ac_sel;

			verify_param_under("DcSel", dc_sel, 4);
			verify_param_under("AcSel", ac_sel, 4);

			segment_length -= 2;

			scan.component_params.push_back(params);
		}

		scan.start_spectral_sel = read_byte();
		scan.end_spectral_sel = read_byte();

		auto [h, l] = read_2x_uint4();
		scan.succ_approx_high = h;
		scan.succ_approx_low = l;

		segment_length -= 3;

		if (segment_length != 0) {
			throw JpegFileError("unexpected final value for length while reading SOS, abort");
		}

		read_entropy_coded_data_to(scan);

		j.frames().back().scans.push_back(std::move(scan));
	}

	void read_entropy_coded_data_to(Scan& s) {
		msg::debug("READ_JPEG: Consume ECS sequence");

		uint8_t b = read_byte();

		std::vector<uint8_t> segment;

		while (true) {
			while (b != 0xFF) {
				segment.push_back(b);
				b = read_byte();
			}

			// Hit possible marker.
			b = read_byte();
			Marker mark { b };

			if (!mark.is_valid()) {
				segment.push_back(0xFF);
				segment.push_back(b);
				b = read_byte();
				continue;
			}

			msg::fine("READ_JPEG: EC: Hit mark {}", std::string(mark));

			// If it's a valid marker, we're at the end of an
			// ECS, so add it to parsed segments.
			s.entropy_coded_data.push_back(std::move(segment));
			segment = std::vector<uint8_t>();

			// Reset marker indicates start of new ECS.
			if (mark.is_rstn()) {
				b = read_byte();
				continue;
			}

			// Otherwise, we are at end of ECS segments.
			// Back up the marker we ate and bail.
			file.unget().unget();
			break;
		}
	}

	void read_special_marker_segment_to(CompressedJpegData& j, const Marker& marker) {
		switch (marker.parse_special()) {
		case MarkerSpecial::dqt_define_quant_table:
			read_qtable_to(j);
			break;
		case MarkerSpecial::dht_define_huff_table:
			read_hufftable_to(j);
			break;
		case MarkerSpecial::dac_define_arith_coding:
			read_arithtable_to(j);
			break;
		case MarkerSpecial::sos_start_of_scan:
			read_scan_to(j);
			break;
		case MarkerSpecial::com_comment:
			read_comment_to(j);
			break;
		case MarkerSpecial::dri_define_rst_interval:
			read_restart_interval();
			break;
		case MarkerSpecial::dhp_define_hierarchical_progression:
		case MarkerSpecial::exp_expand_reference_components:
		case MarkerSpecial::dnl_define_num_lines:
			msg::error(
				"READ_JPEG: Skipping unsupported marker {}",
				symbol_from_special_marker[marker.parse_special()]
			);
			// TODO
			break;
		case MarkerSpecial::eoi_end_of_img:
			msg::debug("READ_JPEG: Hit EOI. Stopping...");
			break;
		default:
			throw JpegFileError(std::format(
				"READ_JPEG: special marker parsed to unknown value {}",
				int(marker)
			));
		}
	}

	void write_jfif_header(const CompressedJpegData& j) {
		msg::debug("WRITE_JPEG: output JFIF header");

		JpegSegmentWriter seg { file };

		write_marker(Marker::appn(0));
		seg.write_str("JFIF");
		seg.write_uint16(0x0102); // JFIF version
		seg.write_byte(0);        // X/Y density unit (0 = no unit)

		// X, Y density (just use pixels)
		seg.write_uint16(j.frames()[0].samples_per_line);
		seg.write_uint16(j.frames()[0].num_lines);

		// X, Y thumbnail (no thumb, so 0)
		seg.write_byte(0);
		seg.write_byte(0);

		seg.commit();
	}

	void write_hufftable(const HuffmanTable& table, uint8_t dest, TableClass cls) {
		msg::debug("WRITE_JPEG: output {} Huffman table {}", str_from_table_class.at(cls), dest);
		JpegSegmentWriter seg { file };

		write_marker(MarkerSpecial::dht_define_huff_table);
		seg.write_2x_uint4((uint8_t)cls, dest);

		for (const auto& s : table.size_amts) {
			seg.write_byte(s);
		}

		for (const auto& c : table.codes) {
			seg.write_byte(c.value);
		}

		seg.commit();
	}

	void write_qtable(const QuantizationTable& table, uint8_t dest) {
		msg::debug("WRITE_JPEG: output quantization table {}", dest);
		JpegSegmentWriter seg { file };
		uint8_t precision = table.precision();

		write_marker(MarkerSpecial::dqt_define_quant_table);
		seg.write_2x_uint4(precision, dest);

		for (int i = 0; i < 64; i++) {
			if (precision == 0) {
				seg.write_byte(table.data.zz(i));
			} else {
				seg.write_uint16(table.data.zz(i));
			}
		}

		seg.commit();
	}

	void write_all_hufftables(const CompressedJpegData& j, TableClass cls) {
		auto& tbls = j.huff_tables(cls);

		for (int dest = 0; dest < tbls.size(); dest++) {
			auto& table = tbls[dest];
			if (table.set) write_hufftable(table, dest, cls);
		}
	}

	void write_all_qtables(const CompressedJpegData& j) {
		auto& tbls = j.q_tables();

		for (int dest = 0; dest < tbls.size(); dest++) {
			auto& table = tbls[dest];
			if (table.set) write_qtable(table, dest);
		}
	}

	void write_sof(const Frame& f) {
		msg::debug("WRITE_JPEG: output SOF segment");
		JpegSegmentWriter seg { file };

		write_marker(Marker::sof(f.sof_info));
		seg.write_byte(f.sample_precision);
		seg.write_uint16(f.num_lines);
		seg.write_uint16(f.samples_per_line);
		seg.write_byte(f.num_components);

		if (f.num_components != f.component_params.size()) {
			throw JpegFileError("write_sof: mismatch in num_components");
		}

		for (const auto& params : f.component_params) {
			seg.write_byte(params.identifier);
			seg.write_2x_uint4(
				params.horizontal_sampling_factor,
				params.vertical_sampling_factor
			);
			seg.write_byte(params.qtable_selector);
		}

		seg.commit();
	}

	void write_sos(const Scan& s) {
		msg::debug("WRITE_JPEG: write SOS segment");
		JpegSegmentWriter seg { file };

		write_marker(MarkerSpecial::sos_start_of_scan);

		if (s.num_components != s.component_params.size()) {
			throw JpegFileError("write_sos: mismatch in num_components");
		}

		seg.write_byte(s.num_components);

		for (const auto& params : s.component_params) {
			seg.write_byte(params.component_selector);
			seg.write_2x_uint4(
				params.dc_entropy_coding_selector,
				params.ac_entropy_coding_selector
			);
		}

		seg.write_byte(s.start_spectral_sel);
		seg.write_byte(s.end_spectral_sel);
		seg.write_2x_uint4(s.succ_approx_high, s.succ_approx_low);

		seg.commit();
	}

	void write_dri(uint16_t restart_interval) {
		msg::debug("WRITE_JPEG: write DRI segment");
		JpegSegmentWriter seg { file };
		write_marker(MarkerSpecial::dri_define_rst_interval);
		seg.write_uint16(restart_interval);
		seg.commit();
	}

	void write_entropy_coded_data(const Scan& s) {
		msg::debug("WRITE_JPEG: writing entropy coded data...");
		uint8_t rst_n = 0;
		for (int i = 0; i < s.entropy_coded_data.size(); i++) {
			const auto& interval = s.entropy_coded_data[i];

			for (const auto& b : interval) file.put(b);

			// For every interval except the last, place a RST
			// marker after.
			if (i < s.entropy_coded_data.size() - 1) {
				write_marker(Marker::rstn(rst_n));
				rst_n++;
				if (rst_n >= 8) rst_n = 0;
			}
		}
	}

	void write_scan(const Scan& s) {
		if (s.restart_interval > 0) {
			write_dri(s.restart_interval);
		}

		write_sos(s);
		write_entropy_coded_data(s);
	}
public:
	JpegFile(std::string path, std::ios_base::openmode mode) {
		msg::debug("FILE_JPEG: try open \"{}\"...", path);

		if (mode != std::ios_base::in && mode != std::ios_base::out) {
			throw JpegFileError(
				"BmpFile(): invalid mode, pick one of: ios_base::in, ios_base::out"
			);
		}

		this->mode = mode;

		// In normal operation, we expect no errors, so if an error
		// occurs throw an exception. This includes EOF; reaching EOF
		// means a malformed JPEG file - something almost never
		// directly in the user's control.
		file.exceptions(std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);
		file.open(path, mode | std::ios_base::binary);
	}

	bool is_read_mode() {
		return mode == std::ios_base::in;
	}

	bool is_write_mode() {
		return mode == std::ios_base::out;
	}

	CompressedJpegData read() {
		Marker mark;
		CompressedJpegData j;
		bool success = true;

		if (!is_read_mode()) throw JpegFileError("JpegFile: cannot read() in write mode");

		if (!read_expect_soi()) {
			throw JpegFileError("File missing SOI marker");
		}

		while (!mark.is_eoi()) {
			mark = read_next_marker();

			if (!mark.is_valid()) {
				msg::error("Hit invalid marker while reading, stop");
				return std::move(j);
			}

			// In case of EOI, this does nothing.
			read_marker_segment_to(j, mark);
		}

		j.set_valid();
		return std::move(j);
	}

	void write(const CompressedJpegData& j) {
		if (!is_write_mode()) throw JpegFileError("JpegFile: cannot write() in read mode");

		write_marker(MarkerSpecial::soi_start_of_img);
		write_jfif_header(j);
		write_all_qtables(j);
		write_sof(j.frames()[0]);
		write_all_hufftables(j, TableClass::dc);
		write_all_hufftables(j, TableClass::ac);

		// TODO: for now only support one scan
		write_scan(j.frames()[0].scans[0]);

		write_marker(MarkerSpecial::eoi_end_of_img);
	}
};

#pragma pack(push, 1)

struct BmpImageHeader {
	uint32_t _header_size = sizeof(BmpImageHeader);
	uint32_t _width = 0;
	int32_t _height = 0;
	uint16_t _planes = 1;
	uint16_t _bit_count = 24;
	uint32_t _compression = 0;
	uint32_t _image_size = 0;
	uint32_t _x_pix_per_meter = 0;
	uint32_t _y_pix_per_meter = 0;
	uint32_t _color_used = 0;
	uint32_t _color_important = 0;

	BmpImageHeader() = default;
	BmpImageHeader(uint32_t width, uint32_t height) {
		_width = width;
		_height = height;
	}

	explicit operator std::string() const {
		return std::format(
			"_header_size={}\n"
			"_width={}\n"
			"_height={}\n"
			"_bit_count={}",
			_header_size,
			_width,
			_height,
			_bit_count
		);
	}
};

struct BmpFileHeader {
	char _type[2] { 'B', 'M' };
	uint32_t _file_size = 0;
	uint32_t _reserved = 0;
	uint32_t _pixel_offset = sizeof(BmpFileHeader) + sizeof(BmpImageHeader);

	BmpFileHeader() = default;
	BmpFileHeader(uint32_t file_size) {
		_file_size = file_size;
	}

	explicit operator std::string() const {
		return std::format(
			"TYPE: {}{}, _file_size={}, _pixel_offset={}",
			_type[0], _type[1], _file_size, _pixel_offset
		);
	}
};

#pragma pack(pop)

class BmpFileError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

struct YCbCrPixel {
	uint8_t y;
	uint8_t cr;
	uint8_t cb;

	uint8_t operator[](unsigned int idx) const {
		if (idx == 0) return y;
		if (idx == 1) return cr;
		if (idx == 2) return cb;

		return 0;
	}
};

// Given some dimension or other value v, calculate how
// much we need to add to v to reach the first integer multiple of d greater
// than or equal to v.
uint32_t calc_pad(uint32_t v, uint32_t d) {
	uint32_t pad = d - v%d;
	if (pad == d) pad = 0;
	return pad;
}

uint8_t calc_bmp_padding_bytes(uint32_t width) {
	return calc_pad(width*3, 4);
}

uint8_t clamped_convert(double x) {
	if (x < 0.0) return 0;
	if (x > 255.0) return 255;

	return (uint8_t)x;
}

// Helper class for performing subsampling and color conversion on
// a read BMP file. Not required for writing.
class BmpSubsamplingWindow {
private:
	std::vector<YCbCrPixel> rows;
	std::fstream* file;
	uint32_t width;
	uint32_t height;
	uint32_t subsamp;
	bool reverse_order;
	bool first_read = true;
	uint32_t consumed_rows = 0;

	// Our own padding in `rows` to simplify padding operations
	uint32_t width_padded;

	// BMP padding bytes. Required to properly read pixel data.
	uint32_t file_padding_bytes;

	// What image coordinate does y=0 currently correspond to?
	uint32_t zero_y_map;

	unsigned int coord(unsigned int x, unsigned int y) const {
		if (y >= subsamp) throw BmpFileError("coord: bad y");
		if (x >= width_padded) throw BmpFileError("coord: bad x");

		return y*width_padded + x;
	}

	void copy_row(unsigned int from, unsigned int to) {
		unsigned int x0_from = coord(0, from);
		unsigned int x0_to = coord(0, to);

		for (unsigned int x = 0; x < width_padded; x++) {
			rows[x0_to + x] = rows[x0_from + x];
		}
	}

	void x_copy(unsigned int y) {
		// If we have horizontal subsample padding, copy the last pixel
		// to the rest of the columns for better averages
		for (int i = width; i < width_padded; i++) {
			rows[coord(i, y)] = rows[coord(width - 1, y)];
		}
	}

	// Copy row `first_filled` to all lower rows. If first_filled == 0,
	// this does nothing.
	void y_copy_lower(unsigned int first_filled) {
		for (int y = 0; y < first_filled; y++) {
			copy_row(first_filled, y);
		}
	}

	// Copy row `last_filled` to all higher rows. If last_filled
	void y_copy_upper(unsigned int last_filled) {
		for (int y = last_filled + 1; y < subsamp; y++) {
			copy_row(last_filled, y);
		}
	}

	// Consume padding bytes, as specified by BMP standard. If
	// `file_padding_bytes` == 0, this does nothing.
	void consume_padding_bytes() {
		for (int i = 0; i < file_padding_bytes; i++) {
			file->get();
		}
	}

	// Consume a row from the file, including all BMP padding.
	void consume_row(unsigned int y) {
		unsigned int x0 = coord(0, y);
		for (unsigned int x = 0; x < width; x++) {
			double r = file->get();
			double g = file->get();
			double b = file->get();

			double luma = 0.299*r + 0.587*g + 0.114*b;
			double cb = (-0.299*r - 0.587*g + 0.886*b) / 1.772 + 128;
			double cr = ( 0.701*r - 0.587*g - 0.114*b) / 1.402 + 128;

			rows[x0 + x] = YCbCrPixel(
				clamped_convert(luma),
				clamped_convert(cr),
				clamped_convert(cb)
			);
		}

		consumed_rows++;
		consume_padding_bytes();
	}

	void init_zero_y_map() {
		if (reverse_order) {
			zero_y_map = (height / subsamp) * subsamp;
		} else {
			zero_y_map = 0;
		}
	}

	void inc_zero_y_map() {
		if (reverse_order) {
			zero_y_map -= subsamp;
		} else {
			zero_y_map += subsamp;
		}
	}
public:
	BmpSubsamplingWindow() = delete;
	BmpSubsamplingWindow(
		std::fstream& file,
		uint32_t width,
		uint32_t height,
		uint32_t subsamp,
		bool reverse_order
	) {
		this->file = &file;
		this->width = width;
		this->height = height;
		this->subsamp = subsamp;
		this->reverse_order = reverse_order;
		this->file_padding_bytes = calc_bmp_padding_bytes(width);
		this->width_padded = width + calc_pad(width, subsamp);

		rows.reserve(subsamp*width_padded);
		for (int i = 0; i < subsamp*width_padded; i++) {
			rows.push_back(YCbCrPixel(0, 0, 0));
		}
	}

	uint32_t get_zero_y_map() const { return zero_y_map; }

	bool no_more_rows() const {
		return file->eof() || consumed_rows >= height;
	}

	// Advance the window by `subsamp` rows. Returns true if there is more
	// data to process, false otherwise.
	bool advance() {
		if (no_more_rows()) return false;

		int y_start = 0;
		if (first_read && reverse_order) {
			msg::warn(
				"READ_BMP: BMP window using reverse order, "
				"skipping height pad on first advance."
			);
			y_start = calc_pad(height, subsamp);
		}

		int last_filled;
		for (int y = y_start; y < subsamp; y++) {
			consume_row(y);
			x_copy(y);
			last_filled = y;
			if (no_more_rows()) break;
		}

		y_copy_upper(last_filled);

		if (first_read) {
			y_copy_lower(y_start);
			init_zero_y_map();
			first_read = false;
		} else {
			inc_zero_y_map();
		}

		return true;
	}

	uint8_t sample(unsigned int samp_x, unsigned int samp_y, unsigned int comp) const {
		// if luma, do not subsample
		if (comp == 0) {
			unsigned int internal_y = samp_y;
			if (reverse_order) internal_y = subsamp - samp_y - 1;
			return rows[coord(samp_x, internal_y)][comp];
		}

		// Here, specific Y value doesn't matter, since we'll be
		// covering all of them in the subsampling operation.
		uint16_t result = 0;
		unsigned int x_start = samp_x * subsamp;
		unsigned int x_end = x_start + subsamp;
		for (int y = 0; y < subsamp; y++) {
			for (int x = x_start; x < x_end; x++) {
				result += rows[coord(x, y)][comp];
			}
		}

		return (uint8_t)(result / (subsamp*subsamp));
	}
};

// Implementation of a BMP file read/writer, only allowing 24-bit RGB color.
// Automatically converts the BMP's RGB to YCbCr.
class BmpFile {
private:
	std::fstream file;
	std::ios_base::openmode mode;

	template<typename T>
	void write_generic(const T& x) {
		// File is little-endian, so can just write stuff out directly
		file.write(reinterpret_cast<const char*>(&x), sizeof(x));
	}

	void write_file_header(uint32_t file_size) {
		write_generic(BmpFileHeader(file_size));
	}

	void write_image_header(uint32_t width, uint32_t height) {
		write_generic(BmpImageHeader(width, height));
	}

	template<typename T>
	T read_generic() {
		T hdr;
		file.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
		return hdr;
	}
public:
	BmpFile(std::string path, std::ios_base::openmode mode) {
		msg::debug("FILE_BMP: try open BmpFile \"{}\"...", path);

		if (mode != std::ios_base::in && mode != std::ios_base::out) {
			throw BmpFileError(
				"BmpFile(): invalid mode, pick one of: ios_base::in, ios_base::out"
			);
		}

		this->mode = mode;

		file.exceptions(std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);
		file.open(path, mode | std::ios_base::binary);
	}

	bool is_read_mode() {
		return mode == std::ios_base::in;
	}

	bool is_write_mode() {
		return mode == std::ios_base::out;
	}

	std::vector<RawComponent> read(uint8_t subsamp=2) {
		if (!is_read_mode()) throw BmpFileError("BmpFile: cannot read() in write mode");

		msg::debug("READ_BMP: reading BMP file and converting to YCbCr...");

		std::vector<RawComponent> raws;

		auto file_header = read_generic<BmpFileHeader>();
		auto t = file_header._type;
		if (t[0] != 'B' || t[1] != 'M') throw BmpFileError("not a BMP file");

		auto image_header = read_generic<BmpImageHeader>();

		if (image_header._header_size < sizeof(BmpImageHeader)) {
			throw BmpFileError("BMP file too old");
		}

		if (image_header._bit_count != 24) {
			throw BmpFileError("BMP: only 24-bit color is supported");
		}

		if (image_header._compression != 0) {
			throw BmpFileError("BMP: compression not supported");
		}

		msg::debug("READ_BMP: checks pass, init components...");

		// Luma component.
		raws.push_back(std::move(RawComponent(
			image_header._width,
			std::abs(image_header._height),
			subsamp, subsamp
		)));

		// Cr, Cb
		for (int i = 0; i < 2; i++) {
			raws.push_back(std::move(RawComponent(
				image_header._width / subsamp,
				std::abs(image_header._height) / subsamp,
				1, 1
			)));
		}

		for (int j = 0; j < 3; j++) {
			msg::debug("READ_BMP: comp {}: {}", j, std::string(raws[j]));
		}

		// Subsampling plus the rows potentially being in reverse order
		// is a massive headache when reading, so use a helper class.
		auto window = BmpSubsamplingWindow(
			file,
			image_header._width,
			std::abs(image_header._height),
			subsamp,
			image_header._height >= 0
		);

		msg::debug("READ_BMP: start conversion...");
		while (window.advance()) {
			unsigned int ymap = window.get_zero_y_map();
			unsigned int ymap_color = ymap / subsamp;

			for (int y = 0; y < subsamp; y++) {
				for (int luma_x = 0; luma_x < image_header._width; luma_x++) {
					uint8_t s = window.sample(luma_x, y, 0);
					raws[0][luma_x, y + ymap] = s;
				}
			}

			for (int comp = 1; comp < 3; comp++) {
				for (int color_x = 0; color_x < image_header._width / subsamp; color_x++) {
					uint8_t s = window.sample(color_x, 0, comp);
					raws[comp][color_x, ymap_color] = s;
				}
			}
		}

		msg::debug("READ_BMP: done.");
		return std::move(raws);
	}

	// Writes the result of a JPEG decode. The components are assumed to
	// be YCbCr, possibly with some degree of chroma subsampling.
	// This will be converted to RGB during the writing process.
	void write(std::vector<RawComponent> components) {
		if (!is_write_mode()) throw BmpFileError("BmpFile: cannot write() in read mode");

		msg::debug("WRITE_BMP: converting to RGB and writing BMP file...");

		uint32_t width = components[0].x_pixels();
		uint32_t height = components[0].y_pixels();
		uint8_t padding_bytes = calc_bmp_padding_bytes(width);

		uint32_t file_size = (
			sizeof(BmpFileHeader)
			+ sizeof(BmpImageHeader)
			+ (width * 3 + padding_bytes) * height
		);

		msg::debug(
			"WRITE_BMP: BMP fsize={} padding={}",
			file_size, padding_bytes
		);

		write_file_header(file_size);
		// height is given as a negative number so the rows are
		// ordered top-to-bottom.
		write_image_header(width, -(int32_t)height);

		uint8_t cb_subsamp_x = components[0].mcu_width() / components[2].mcu_width();
		uint8_t cb_subsamp_y = components[0].mcu_height() / components[2].mcu_height();
		uint8_t cr_subsamp_x = components[0].mcu_width() / components[1].mcu_width();
		uint8_t cr_subsamp_y = components[0].mcu_height() / components[1].mcu_height();


		for (uint32_t y = 0; y < height; y++) {
			for (uint32_t x = 0; x < width; x++) {
				double luma = components[0][x, y];

				// NOTE: This works, but it is the opposite of
				// what the JFIF standard says... component[1]
				// should be Cb, not Cr.
				double cb = components[2][
					x / cb_subsamp_x,
					y / cb_subsamp_y
				];
				double cr = components[1][
					x / cr_subsamp_x,
					y / cr_subsamp_y
				];

				double r = luma + 1.402 * (cr - 128.0);
				double g = luma - 0.34414 * (cb - 128.0) - 0.71414 * (cr - 128.0);
				double b = luma + 1.772 * (cb - 128.0);

				file.put(clamped_convert(r));
				file.put(clamped_convert(g));
				file.put(clamped_convert(b));
			}

			for (int i = 0; i < padding_bytes; i++) {
				file.put(0);
			}
		}
	}
};
}
