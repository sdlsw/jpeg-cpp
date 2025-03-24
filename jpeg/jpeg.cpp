module;

#include <cerrno>

export module jpeg;

import std;
import msg;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export import :data;
export import :decode;

export namespace jpeg {
class JpegFile {
private:
	std::fstream file;
	uint16_t reset_interval = 0;

	uint8_t read_byte() {
		return (uint8_t) file.get();
	}

	void read_until_byte(uint8_t b) {
		uint8_t rb = 0;
		while (rb != b) rb = read_byte();
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

	uint16_t read_uint16() {
		msg::fine("READ: 16-bit int");

		// File is big-endian
		uint16_t high = read_byte();
		uint16_t low = read_byte();

		msg::fine("READ: high {} low {}", high, low);

		return (high << 8) | low;
	}

	// Reads segment size. In the file, this value includes the size of the
	// length field itself (two bytes). This function returns the size
	// of everything *after* the length field, and sets the cursor to
	// the byte after it.
	uint16_t read_segment_size() {
		msg::fine("READ: read_segment_size");
		return read_uint16() - 2;
	}

	// Read a byte from the file, that contains two four-bit values.
	std::tuple<uint8_t, uint8_t> read_2x_uint4() {
		uint8_t b = read_byte();
		return { (b & 0xF0) >> 4, (b & 0x0F) };
	}

	// Reads a segment consisting of only a length field followed by
	// arbitrary binary data.
	bool read_opaque_segment_to(std::vector<uint8_t>& dest) {
		uint16_t length = read_segment_size();
		msg::debug("READ: reading opaque segment size {}", length);

		dest.reserve(length);
		while (length > 0) {
			dest.push_back(read_byte());
			length--;
		}

		return true;
	}

	bool read_marker_segment_to(CompressedJpegData& j, const Marker& marker) {
		if (marker.is_appn()) {
			return read_appn_to(j, marker.parse_appn());
		} else if (marker.is_sof()) {
			return read_sof_to(j, marker.parse_sof());
		} else if (marker.is_rstn()) {
			read_rstn_to(j, marker.parse_rstn());
		} else if (marker.is_res() || marker.is_jpgn()) {
			msg::warn("Skipping reserved marker {}", std::string(marker));
		} else {
			return read_special_marker_segment_to(j, marker);
		}

		return true;
	}

	bool read_reset_interval() {
		msg::debug("READ: process DRI");

		uint16_t segment_length = read_segment_size();

		if (segment_length != 2) {
			msg::error("READ: DRI: Bad segment length {}", segment_length);
			return false;
		}

		// Save reset interval for later. It will be applied to any
		// following parsed SOS segments.
		reset_interval = read_uint16();
		return true;
	}

	bool read_appn_to(CompressedJpegData& j, int n) {
		msg::debug("READ: process APP{}", n);

		AppSegment appseg;
		appseg.n = n;
		if (!read_opaque_segment_to(appseg.data)) return false;
		j.app_segments().push_back(std::move(appseg));

		return true;
	}

	bool read_comment_to(CompressedJpegData& j) {
		msg::debug("READ: process COM");

		std::vector<uint8_t> dest;
		if (!read_opaque_segment_to(dest)) return false;
		j.comments().push_back(std::move(dest));

		return true;
	}

	bool read_sof_to(CompressedJpegData& j, const StartOfFrameInfo& i) {
		msg::debug("READ: process {}", std::string(i));

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

			if (params.horizontal_sampling_factor - 1 > 3 ||
			    params.vertical_sampling_factor - 1 > 4) {
				msg::error("bad sampling factor(s) while reading frame header");
				return false;
			}

			params.qtable_selector = read_byte();

			frame.component_params.push_back(params);
			segment_length -= 3;
		}

		if (segment_length != 0) {
			msg::error("unexpected final value for length while reading SOF, abort");
			return false;
		}

		j.frames().push_back(std::move(frame));
		return true;
	}

	void read_rstn_to(CompressedJpegData& j, int n) {
		msg::debug("READ: process RST{}", n);
	}

	bool read_qtable_to(CompressedJpegData& j) {
		msg::debug("READ: process DQT");

		uint16_t length = read_segment_size();

		while (length > 0) {
			auto [pq_precision, tq_destination] = read_2x_uint4();

			if (tq_destination >= 4) {
				msg::error("Tq parameter invalid, Tq={}>=4", tq_destination);
				return false;
			}

			if (pq_precision >= 2) {
				msg::error("Pq parameter invalid, Pq={}>=2", pq_precision);
				return false;
			}

			QuantizationTable& qtable = j.q_tables()[tq_destination];
			msg::debug("READ: install qtable to dest {}", tq_destination);

			if (qtable.set) {
				msg::warn("READ: replacing existing table!");
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
			msg::error("unexpected final value for length while reading qtables, abort");
			return false;
		}

		return true;
	}

	bool read_hufftable_to(CompressedJpegData& j) {
		msg::debug("READ: process DHT");

		uint16_t segment_length = read_segment_size();

		while (segment_length > 0) {
			auto [tc_table_class, th_destination] = read_2x_uint4();

			if (th_destination >= 4) {
				msg::error("Th parameter invalid, Th={}>=4", th_destination);
				return false;
			}

			if (tc_table_class >= 2) {
				msg::error("Tc parameter invalid, Tc={}>=2", tc_table_class);
				return false;
			}

			segment_length--;

			auto _class = TableClass(tc_table_class);
			HuffmanTable& hufftable = j.huff_tables(_class)[th_destination];
			msg::debug(
				"READ: install {} hufftable to dest {}",
				str_from_table_class.at(_class),
				th_destination
			);

			if (hufftable.set) {
				msg::warn("READ: replacing existing table!");
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
			msg::error("unexpected final value for length while reading hufftables, abort");
			return false;
		}

		return true;
	}

	void read_arithtable_to(CompressedJpegData& j) {
		msg::debug("READ: process DAC");
	}

	bool read_scan_to(CompressedJpegData& j) {
		msg::debug("READ: process SOS");

		if (j.frames().size() == 0) {
			msg::error("encountered SOS marker before SOF");
			return false;
		}

		uint16_t segment_length = read_segment_size();

		Scan scan;
		scan.reset_interval = reset_interval;
		scan.num_components = read_byte();
		segment_length--;

		if (scan.num_components - 1 > 3) {
			msg::error("invalid Ns field while reading SOS");
			return false;
		}

		for (int i = 0; i < scan.num_components; i++) {
			ScanComponentParams params;

			params.component_selector = read_byte();

			auto [dc_sel, ac_sel] = read_2x_uint4();
			params.dc_entropy_coding_selector = dc_sel;
			params.ac_entropy_coding_selector = ac_sel;

			if (params.dc_entropy_coding_selector > 3 ||
			    params.ac_entropy_coding_selector > 3) {
				msg::error("SOS: coding selector too large for component {}", i);
				return false;
			}

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
			msg::error("unexpected final value for length while reading SOS, abort");
			return false;
		}

		if (!read_entropy_coded_data_to(scan)) return false;

		j.frames().back().scans.push_back(std::move(scan));
		return true;
	}

	bool read_entropy_coded_data_to(Scan& s) {
		msg::debug("READ: Consume ECS sequence");

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

			msg::debug("READ: EC: Hit mark {}", std::string(mark));

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
			return true;
		}

	}

	bool read_special_marker_segment_to(CompressedJpegData& j, const Marker& marker) {
		switch (marker.parse_special()) {
		case MarkerSpecial::dqt_define_quant_table:
			return read_qtable_to(j);
		case MarkerSpecial::dht_define_huff_table:
			return read_hufftable_to(j);
		case MarkerSpecial::dac_define_arith_coding:
			read_arithtable_to(j);
			break;
		case MarkerSpecial::sos_start_of_scan:
			return read_scan_to(j);
		case MarkerSpecial::com_comment:
			return read_comment_to(j);
		case MarkerSpecial::dri_define_rst_interval:
			return read_reset_interval();
		case MarkerSpecial::dhp_define_hierarchical_progression:
		case MarkerSpecial::exp_expand_reference_components:
		case MarkerSpecial::dnl_define_num_lines:
			msg::error(
				"READ: Skipping unsupported marker {}",
				symbol_from_special_marker[marker.parse_special()]
			);
			// TODO
			break;
		case MarkerSpecial::eoi_end_of_img:
			msg::debug("READ: Hit EOI. Stopping...");
			break;
		default:
			msg::error("READ: special marker parsed to unknown value {}", int(marker));
			return false;
		}

		return true;
	}

public:
	JpegFile(std::string path) {
		msg::debug("READ: try open JpegFile \"{}\"...", path);

		// In normal operation, we expect no errors, so if an error
		// occurs throw an exception.
		file.exceptions(std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);
		file.open(path, std::ios_base::in | std::ios_base::out | std::ios_base::binary);
	}

	CompressedJpegData read() {
		Marker mark;
		CompressedJpegData j;
		bool success = true;

		if (!read_expect_soi()) {
			msg::error("File missing SOI marker");
			return std::move(j);
		}

		while (!mark.is_eoi()) {
			mark = read_next_marker();

			if (!mark.is_valid()) {
				msg::error("Hit invalid marker while reading, stop");
				return std::move(j);
			}

			// In case of EOI, this does nothing.
			success = read_marker_segment_to(j, mark);
			if (!success) return std::move(j);
		}

		j.set_valid();
		return std::move(j);
	}
};
}
