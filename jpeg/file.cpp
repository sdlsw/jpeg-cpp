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
// Used to indicate failure to process a JPEG file.
class JpegFileError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// Helper class for JpegFile. Each segment in a JPEG file must start with
// the number of bytes in the segment (including the 2 bytes for the length
// value). To make this easier, we copy the data to a vector first, then call
// commit() to write the whole segment at once, including the length.
//
// This is a little less efficient than calculating the segment size ahead
// of time, but in practice it doesn't matter much. By far the largest amount
// of data is in the entropy coded segments, which are not written through
// JpegSegmentWriter as they aren't preceded by a length value.
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

	// Writes the length of the segment to file.
	void put_length() {
		auto [h, l] = uint16_be(data.size() + 2);
		file->put(h);
		file->put(l);
	}
public:
	JpegSegmentWriter(std::fstream& file) {
		this->file = &file;
	}

	// Writes the complete segment to file.
	void commit() {
		put_length();
		for (auto& x : data) {
			file->put(x);
		}
	}

	// Writes a byte to the segment.
	void write_byte(uint8_t b) {
		data.push_back(b);
	}

	// Writes an unsigned 16-bit integer to the segment.
	void write_uint16(uint16_t x) {
		auto [h, l] = uint16_be(x);
		data.push_back(h);
		data.push_back(l);
	}

	// Writes two 4-bit integers as one byte - common in the JPEG
	// headers. Bits: HHHHLLLL.
	void write_2x_uint4(uint8_t h, uint8_t l) {
		l &= 0x0F;

		data.push_back((h << 4) | l);
	}

	// Writes a string to the segment, including null terminator.
	void write_str(const std::string& s) {
		for (auto& c : s) {
			data.push_back(c);
		}

		data.push_back('\0');
	}
};

// Helper class for JpegFile. Reading version of JpegSegmentWriter. Reads the
// 2-byte length field first, then counts remaining bytes in the segment as
// they are read. `verify()` may then be called to ensure the correct number
// of bytes have been consumed.
class JpegSegmentReader {
private:
	std::fstream* file;
	const std::string segment_name;
	unsigned int segment_length;
	int count = 0;
public:
	// Note: On construction, the reader will automatically attempt
	// to read the length of the segment. This may throw an exception
	// if at EOF.
	JpegSegmentReader(
		std::fstream& file,
		std::string&& segment_name
	) : file{&file}, segment_name{segment_name} {
		msg::debug("READ_JPEG: process {}", segment_name);

		// Byte count includes the 2 byte length of the byte count
		// field itself, so use regular read function. This will result
		// in `count` temporarily dropping to -2, hence it is `int`.
		segment_length = read_uint16();
		count += segment_length;
	}

	bool has_bytes() {
		return count > 0;
	}

	// Verifies the segment we read. Does nothing on success, throws
	// an exception if the segment was malformed.
	void verify() {
		if (count != 0) {
			throw JpegFileError(std::format(
				"{}: Bad segment read (len={} left={})",
				segment_name, segment_length, count
			));
		}
	}

	uint8_t read_byte() {
		auto in = (uint8_t) file->get();
		--count;
		return in;
	}

	uint16_t read_uint16() {
		uint16_t h = read_byte();
		uint16_t l = read_byte();
		return (h << 8) | l;
	}

	// Reads a byte from the file that contains two four-bit values.
	std::tuple<uint8_t, uint8_t> read_2x_uint4() {
		uint8_t b = read_byte();
		return { (b & 0xF0) >> 4, (b & 0x0F) };
	}

	// Reads everything remaining in the segment to a new vector. No need
	// to verify() after this.
	std::vector<uint8_t> read_all() {
		std::vector<uint8_t> result;

		result.reserve(count);
		while (count > 0) {
			result.push_back(read_byte());
		}

		return result;
	}
};

// JPEG file reader/writer. All read functions will throw JpegFileError
// in case of malformed JPEG files, since there is no way to reasonably recover
// from that. In practice, malformed JPEG files are uncommon.
class JpegFile {
private:
	std::fstream file;
	std::ios_base::openmode mode;
	uint16_t restart_interval = 0;

	// Ensures a header parameter is under param_limit, throwing an exception
	// if not. The majority of JPEG parameters are unsigned and have a
	// lower bound of 0, so it is sufficient to check only an upper bound.
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

	// Reads one byte from the JPEG file. Automatically converts from the
	// `int` returned by `file.get` since exceptions are enabled - checking
	// for error returns is unnecessary.
	uint8_t read_byte() {
		return (uint8_t) file.get();
	}

	// Reads bytes until the value `b` is encountered. The file cursor
	// will be placed immediately after it.
	void read_until_byte(uint8_t b) {
		uint8_t rb = 0;
		while (rb != b) rb = read_byte();
	}

	// Writes a JPEG file marker, 0xFFXX.
	void write_marker(const Marker& m) {
		file.put(0xFF).put(m);
	}

	// Shortcut, allows writing a special marker constant directly
	void write_marker(MarkerSpecial m) {
		file.put(0xFF).put((uint8_t)m);
	}

	// Reads next marker, placing the cursor after the two byte marker value.
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

	// Returns `true` if the next marker in the file is SOI (start of
	// image), and `false` otherwise. May throw std::ios_base::failure if EOF is
	// reached.
	bool read_expect_soi() {
		Marker mark = read_next_marker();
		if (!mark.is_valid()) return false;
		return mark.is_soi();
	}

	// Reads a DRI (define restart interval) segment and sets this reader's
	// current restart interval. This restart interval will be used for any
	// following SOS (start of scan) segments, until new DRI segment is
	// read.
	void read_restart_interval() {
		JpegSegmentReader seg { file, "DRI" };
		restart_interval = seg.read_uint16();
		seg.verify();
	}

	// Reads an application segment and places it in the compressed data
	// package.
	void read_appn_to(CompressedJpegData& j, int n) {
		JpegSegmentReader seg { file, std::format("APP{}", n) };
		AppSegment appseg;

		appseg.n = n;
		appseg.data = seg.read_all();
		j.app_segments().push_back(std::move(appseg));
	}

	// Reads a comment segment and places it in the compressed data
	// package.
	void read_comment_to(CompressedJpegData& j) {
		JpegSegmentReader seg { file, "COM" };
		j.comments().push_back(seg.read_all());
	}

	// Reads a SOF (start of frame) segment, and starts a new frame in the
	// compressed JPEG data. Following SOS (start of scan) segments will
	// start scans in this most recent frame.
	void read_sof_to(CompressedJpegData& j, const FrameType& i) {
		JpegSegmentReader seg { file, std::string(i) };

		Frame frame;
		frame.type = i;

		frame.sample_precision = seg.read_byte();
		frame.num_lines = seg.read_uint16();
		frame.samples_per_line = seg.read_uint16();
		frame.num_components = seg.read_byte();

		for (int comp = 0; comp < frame.num_components; comp++) {
			FrameComponentParams params;

			params.identifier = seg.read_byte();

			auto [h, v] = seg.read_2x_uint4();
			params.horizontal_sampling_factor = h;
			params.vertical_sampling_factor = v;

			verify_param_under("(Hi-1)", h - 1, 4);
			verify_param_under("(Vi-1)", v - 1, 4);

			params.qtable_selector = seg.read_byte();

			frame.component_params.push_back(params);
		}

		seg.verify();

		j.frames().push_back(std::move(frame));
	}

	// Currently does nothing, just announces that we stepped over a RSTn
	// (restart N) marker.
	void read_rstn_to(CompressedJpegData& j, int n) {
		msg::debug("READ_JPEG: process RST{}", n);
	}

	// Reads a DQT (define quantization table) segment, and installs the
	// defined tables to their appropriate places in the CompressedJpegData
	// object.
	void read_qtable_to(CompressedJpegData& j) {
		JpegSegmentReader seg { file, "DQT" };

		// No way to tell how many tables left without just tracking
		// byte count (no table count parameter).
		while (seg.has_bytes()) {
			auto [pq_precision, tq_destination] = seg.read_2x_uint4();

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
					value = seg.read_byte();
				} else {
					//16-bit precision (pq = 1)
					value = seg.read_uint16();
				}

				auto& valref = qtable.data.zz(i);
				valref = value;
			}

			qtable.set = true;
		}

		seg.verify();
	}

	// Reads a DHT (define huffman table) segment, and installs the
	// defined tables to their appropriate places in the CompressedJpegData
	// object.
	void read_hufftable_to(CompressedJpegData& j) {
		JpegSegmentReader seg { file, "DHT" };

		while (seg.has_bytes()) {
			auto [tc_table_class, th_destination] = seg.read_2x_uint4();

			verify_param_under("Th", th_destination, 4);
			verify_param_under("Tc", tc_table_class, 2);

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
				uint8_t len = seg.read_byte();
				row_lengths.push_back(len);
			}

			// Using the length declarations, read huffman table.
			// Specifically, this is the HUFFVAL table.
			// HUFFSIZE is also populated during this step.
			for (int i = 0; i < 16; i++) {
				auto huff_row_length = row_lengths[i];
				while (huff_row_length > 0) {
					uint8_t val = seg.read_byte();

					HuffmanCode entry;
					entry.value = val;
					entry.bits = i+1;

					hufftable.codes.push_back(entry);
					huff_row_length--;
				}
			}

			// Populate HUFFCODE table
			hufftable.populate_codes();
			hufftable.set = true;
		}

		seg.verify();
	}

	// Currently does nothing, since arithmetic coding is not implemented.
	void read_arithtable_to(CompressedJpegData& j) {
		msg::debug("READ_JPEG: process DAC");
	}

	// Reads a SOS (start of scan) segment and adds it to the most recent
	// frame in the CompressedJpegData object. This also consumes all of
	// the following entropy coded data segments until the next non-RST
	// marker.
	void read_scan_to(CompressedJpegData& j) {
		JpegSegmentReader seg { file, "SOS" };

		if (j.frames().size() == 0) {
			throw JpegFileError("encountered SOS marker before SOF");
		}

		Scan scan;
		scan.restart_interval = restart_interval; // Defined earlier by DRI
		scan.num_components = seg.read_byte();

		verify_param_under("(Ns-1)", scan.num_components - 1, 4);

		for (int i = 0; i < scan.num_components; i++) {
			ScanComponentParams params;

			params.component_selector = seg.read_byte();

			auto [dc_sel, ac_sel] = seg.read_2x_uint4();
			params.dc_entropy_coding_selector = dc_sel;
			params.ac_entropy_coding_selector = ac_sel;

			verify_param_under("DcSel", dc_sel, 4);
			verify_param_under("AcSel", ac_sel, 4);

			scan.component_params.push_back(params);
		}

		scan.start_spectral_sel = seg.read_byte();
		scan.end_spectral_sel = seg.read_byte();

		auto [h, l] = seg.read_2x_uint4();
		scan.succ_approx_high = h;
		scan.succ_approx_low = l;

		seg.verify();

		read_entropy_coded_data_to(scan);

		j.frames().back().scans.push_back(std::move(scan));
	}

	// Reads entropy coded data segments to the given Scan, until a non-RST
	// marker is encountered.
	void read_entropy_coded_data_to(Scan& s) {
		// No JpegSegmentReader here, entropy coded data does not
		// include a length header, unlike every other segment.
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

	// Reads a "special" segment to the given CompressedJpegData object.
	// Special markers are any marker that does not contain additional
	// information in its value. So, we can call Marker::parse_special once
	// and then use `switch`.
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
				symbol_from_special_marker.at(marker.parse_special())
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

	// Reads a segment to the given CompressedJpegData object, depending
	// on the Marker that was read just before.
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

	// Writes the APP0 segment describing this JPEG as a JFIF-compliant
	// file. Nearly all JPEGs found on the internet today are JFIF. (At
	// least all the onces I've tested)
	void write_jfif_header(const CompressedJpegData& j) {
		msg::debug("WRITE_JPEG: output JFIF header");

		JpegSegmentWriter seg { file };

		write_marker(Marker::appn(0));
		seg.write_str("JFIF");
		seg.write_uint16(0x0102); // JFIF version 1.02
		seg.write_byte(0);        // X/Y density unit (0 = no unit)

		// X, Y density (just use pixels)
		seg.write_uint16(j.frames()[0].samples_per_line);
		seg.write_uint16(j.frames()[0].num_lines);

		// X, Y thumbnail (no thumb, so 0)
		seg.write_byte(0);
		seg.write_byte(0);

		seg.commit();
	}

	// Writes a DHT (define huffman table) segment to this file.
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

	// Writes a DQT (define quantization table) segment to this file.
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

	// Writes every set Huffman table in the given CompressedJpegData object
	// to this file.
	void write_all_hufftables(const CompressedJpegData& j, TableClass cls) {
		auto& tbls = j.huff_tables(cls);

		for (int dest = 0; dest < tbls.size(); dest++) {
			auto& table = tbls[dest];
			if (table.set) write_hufftable(table, dest, cls);
		}
	}

	// Writes every set quantization table in the given CompressedJpegData
	// object to this file.
	void write_all_qtables(const CompressedJpegData& j) {
		auto& tbls = j.q_tables();

		for (int dest = 0; dest < tbls.size(); dest++) {
			auto& table = tbls[dest];
			if (table.set) write_qtable(table, dest);
		}
	}

	// Writes a SOF (start of frame) segment to this file for the given
	// frame.
	void write_sof(const Frame& f) {
		msg::debug("WRITE_JPEG: output SOF segment");
		JpegSegmentWriter seg { file };

		write_marker(Marker::sof(f.type));
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

	// Writes a SOS (start of scan) segment to this file for the given
	// scan.
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

	// Writes a DRI (define restart interval) segment to this file.
	// Should be written before any scans.
	void write_dri(uint16_t restart_interval) {
		msg::debug("WRITE_JPEG: write DRI segment");
		JpegSegmentWriter seg { file };
		write_marker(MarkerSpecial::dri_define_rst_interval);
		seg.write_uint16(restart_interval);
		seg.commit();
	}

	// Writes the entropy coded data from a given scan to this file.
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

	// Writes a given scan to this file.
	void write_scan(const Scan& s) {
		// For now this is fine, since we only ever write one scan
		// with all components interleaved. In the future for
		// multi-scan JPEGs, might want to write DRI once before all
		// the scans.
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

// Represents the image header of a BMP file, specifically the Windows 3.x
// version. The defaults listed here are used to simplify writing in
// BmpFile::write().
//
// See https://www.fileformat.info/format/bmp/egff.htm
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

// Represents the initial file header of a BMP. The `_type` field must always
// be 'BM'.
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

// Used to indicate failure to process a BMP file.
class BmpFileError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// Represents a pixel in YCbCr color space. Note that we order the color
// components as JPEGs in the wild do, YCrCb (red first).
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

// Given some dimension or other value v, calculate how much we need to add to
// v to reach the first integer multiple of d greater than or equal to v.
uint32_t calc_pad(uint32_t v, uint32_t d) {
	uint32_t pad = d - v%d;
	if (pad == d) pad = 0;
	return pad;
}

// Calculates the number of padding bytes to use per row in a BMP file, given
// the width of the rows in pixels.
uint8_t calc_bmp_padding_bytes(uint32_t width) {
	return calc_pad(width*3, 4);
}

// Converts a double value to a uint8_t, clamping the value to the limits of
// an 8 bit unsigned integer.
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

	// Our own pixel padding in `rows` to simplify padding operations. Must
	// be an integer multiple of `subsamp`.
	uint32_t width_padded;

	// BMP padding bytes. Required to properly read pixel data.
	uint32_t file_padding_bytes;

	// What image row does a sampling coordinate of y=0 currently
	// correspond to?
	//
	// For these examples, we use a `subsamp` value of 3. If
	// `reverse_order` is false, then `zero_y_map` starts at zero and
	// increments by 3 as the window advances through the file - this is
	// the same order as the image (top-down).
	//
	//                   y   windows
	//                      ______...
	// zero_y_map=0 ---> 0 |      ...
	//                   1 |      ...   BMP FILE
	//                   2 |______...      |
	// zero_y_map=3 ---> 0 |      ...      |
	//                   1 |      ...      V
	//                   2 |______...
	//                   ...
	//
	// However if `reverse_order` is true, then `zero_y_map` starts at the last
	// multiple of 3 contained in the image and decrements by 3, since the
	// BMP is storing the rows bottom-up. The sampling coordinate y is also 
	// inverted.
	//
	//                     y   windows
	//                        ______...
	//                     2 |      ...
	//                     1 |      ...   BMP FILE
	// zero_y_map=L   ---> 0 |______...      |
	//                     2 |      ...      |
	//                     1 |      ...      V
	// zero_y_map=L-3 ---> 0 |______...
	//                     ...
	//
	// IMPORTANT: All private functions not explicitly mentioning
	// zero_y_map use *internal* y coordinates: y=0 always refers to the
	// first row of the `rows` vector, regardless of `reverse_order`. 
	uint32_t zero_y_map;

	// Calculates the index into the `rows` vector given (x, y)
	// coordinates.
	unsigned int coord(unsigned int x, unsigned int y) const {
		if (y >= subsamp) throw BmpFileError("coord: bad y");
		if (x >= width_padded) throw BmpFileError("coord: bad x");

		return y*width_padded + x;
	}

	// Copies row y=`from` to y=`to`.
	void copy_row(unsigned int from, unsigned int to) {
		unsigned int x0_from = coord(0, from);
		unsigned int x0_to = coord(0, to);

		for (unsigned int x = 0; x < width_padded; x++) {
			rows[x0_to + x] = rows[x0_from + x];
		}
	}

	// If we have horizontal subsample padding, copy the last pixel
	// to the rest of the columns for better averages
	void x_copy(unsigned int y) {
		for (int i = width; i < width_padded; i++) {
			rows[coord(i, y)] = rows[coord(width - 1, y)];
		}
	}

	// Copies row y=`first_filled` to all rows of lower y coordinate. If
	// first_filled == 0, this does nothing.
	void y_copy_lower(unsigned int first_filled) {
		for (int y = 0; y < first_filled; y++) {
			copy_row(first_filled, y);
		}
	}

	// Copies row `last_filled` to all rows of higher y coordinate. If
	// last_filled == subsamp - 1 this does nothing.
	void y_copy_upper(unsigned int last_filled) {
		for (int y = last_filled + 1; y < subsamp; y++) {
			copy_row(last_filled, y);
		}
	}

	// Consumes padding bytes, as specified by BMP standard. If
	// `file_padding_bytes` == 0, this does nothing.
	void consume_padding_bytes() {
		for (int i = 0; i < file_padding_bytes; i++) {
			file->get();
		}
	}

	// Consumes a row from the file, including all BMP padding.
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

	// Initializes zero_y_map. See comment on that variable.
	void init_zero_y_map() {
		if (reverse_order) {
			int m = height % subsamp;
			if (m == 0) {
				zero_y_map = height - subsamp;
			} else {
				zero_y_map = height - m;
			}
		} else {
			zero_y_map = 0;
		}
	}

	// Advances zero_y_map. See comment on that variable.
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

		// Init rows. This vector will never be resized, only written
		// over.
		rows.reserve(subsamp*width_padded);
		for (int i = 0; i < subsamp*width_padded; i++) {
			rows.push_back(YCbCrPixel(0, 0, 0));
		}
	}

	uint32_t get_zero_y_map() const { return zero_y_map; }

	bool no_more_rows() const {
		return file->eof() || consumed_rows >= height;
	}

	// Advances the window by `subsamp` rows. Returns true if we were able
	// to advance, false otherwise.
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

	// Samples the window. If comp=0 (luma) no subsampling is performed;
	// the exact luma value at (x, y) of this window will be returned.
	//
	// If comp=1,2 (chroma), then samp_y is ignored and a subsamp^2 square
	// will be averaged to produce the subsampled value, starting at `samp_x`.
	// For example, with subsamp=3, samp_x=3:
	//
	//  x 012345678...
	// y     V
	// 0  ...sss......
	// 1  ...sss......
	// 2  ...sss......
	// 
	uint8_t sample(unsigned int samp_x, unsigned int samp_y, unsigned int comp) const {
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

	// Writes a generic type byte-for-byte out to the file.
	template<typename T>
	void write_generic(const T& x) {
		// File is little-endian, so can just write stuff out directly
		file.write(reinterpret_cast<const char*>(&x), sizeof(x));
	}

	// Writes a BMP file header. The parameters other than file_size are
	// determined by the default values specified in BmpFileHeader.
	void write_file_header(uint32_t file_size) {
		write_generic(BmpFileHeader(file_size));
	}

	// Writes a BMP image header. Same deal as write_file_header.
	void write_image_header(uint32_t width, uint32_t height) {
		write_generic(BmpImageHeader(width, height));
	}

	// Reads a generic type byte-for-byte from the file. Must be manually
	// instantiated, for example `auto x = read_generic<MyType>();`.
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

	// Reads this BMP file with some degree of chroma subsampling.
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
