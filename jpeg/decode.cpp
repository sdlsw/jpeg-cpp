export module jpeg:decode;

import std;
import msg;

import :codingbase;
import :data;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export namespace jpeg {
// Helper class that keeps track of all the indices needed to traverse
// entropy-coded data segments for decoding.
class ScanDecodeView : public DcPredictor {
private:
	const Scan* _scan;
	unsigned int rst_idx = 0;  // Current interval in scan.
	unsigned int byte_idx = 0; // Current byte in current RST interval.
	uint8_t bit_idx = 0;       // Current bit in current byte.

	uint8_t cur_byte = 0; // Current byte.

	bool is_end_of_intervals() const {
		return rst_idx >= _scan->entropy_coded_data.size();
	}

	const std::vector<uint8_t>& cur_interval() const {
		if (is_end_of_intervals()) {
			throw std::out_of_range("unexpected end of scan data");
		}

		return _scan->entropy_coded_data[rst_idx];
	}

	// Read and populate the next byte. Returns true if there is now a
	// byte available to read, and false if there is none.
	bool next_byte() {
		if (is_end_of_intervals()) return false;

		if (byte_idx >= cur_interval().size()) {
			msg::error("DECODE: next_byte: cannot get any more bytes this interval");
			return false;
		}

		// Handle byte stuffing.
		uint8_t next_byte = cur_interval()[byte_idx];
		if (cur_byte == 0xFF && next_byte == 0x00) {
			byte_idx++;
			next_byte = cur_interval()[byte_idx];
		}

		cur_byte = next_byte;
		byte_idx++;
		bit_idx = 0;

		return true;
	}
public:
	// This function must be called from outside.
	bool next_rst() {
		rst_idx++;
		if (is_end_of_intervals()) {
			msg::debug("DECODE: reached end of intervals");
			return false;
		}

		byte_idx = 0;
		if(!next_byte()) return false;

		// DC predictor is reset on each RST interval.
		dc_reset();

		return true;
	}

	static int16_t extend(int16_t v, int16_t t) {
		int16_t v_t = 0x1 << (t - 1);

		if (v < v_t) {
			// negative value
			return (0xFFFF << t) + v + 1;
		} else {
			// positive value
			return v;
		}
	}

	ScanDecodeView(const Scan& scan) : DcPredictor(scan.num_components), _scan(&scan) {
		// Prime current byte
		next_byte();
	}

	const Scan& scan() const { return *_scan; }

	int8_t next_bit() {
		if (bit_idx >= 8) {
			if (!next_byte()) return -1;
		}

		uint8_t bit = (cur_byte >> (7 - bit_idx)) & 0x1;

		bit_idx++;
		return bit;
	}

	// Decodes a huffman-coded symbol. The symbol itself will be 8 bits;
	// a special value of -1 indicates no value could be read.
	int16_t decode_symbol(const HuffmanTable& table) {
		uint16_t code = 0;
		int8_t bit = 0;
		int16_t value = 0;

		for (int bits = 1; bits <= 16; bits++) {
			bit = next_bit();
			if (bit == -1) {
				msg::error("decode_symbol: Could not read bit!");
				return -1;
			}

			code = (code << 1) | ((uint16_t) bit);
			value = table.lookup_code(code, bits);
			if (value != -1) return value;
		}

		msg::error("decode_symbol: Could not interpret next coded symbol");
		return -1;
	}

	// Pulls N bits from the bitstream.
	int16_t next_n_bits(unsigned int n) {
		if (n == 0) return 0;

		// We use the sign bit to indicate failure, so max meaningful
		// bits we can pull is 15. This is fine for JPEG decoding,
		// which will never need to pull more that 11 bits at a time.
		if (n > 15) {
			throw std::out_of_range("next_n_bits: n > 15");
		}

		// TODO do this more efficiently?
		uint16_t bits = 0;
		int8_t cur_bit = 0;
		for (unsigned int i = 0; i < n; i++) {
			cur_bit = next_bit();
			if (cur_bit == -1) return -1;

			bits = (bits << 1) | cur_bit;
		}

		return bits;
	}

	int16_t decode_dc_coeff(const HuffmanTable& table) {
		int16_t& pred = dc_pred();
		int16_t ssss = decode_symbol(table);

		if (ssss < 0) {
			throw std::runtime_error("decode_dc_coeff: failed to decode ssss symbol");
		}

		// Don't need to do anything in this case, so just end early
		if (ssss == 0) {
			return pred;
		}

		int16_t v = next_n_bits(ssss);
		if (v < 0) throw std::runtime_error("decode_dc_coeff: failed to decode diff");

		int16_t diff = extend(v, ssss);

		int16_t coeff = pred + diff;
		pred = coeff;

		return coeff;
	}

	AcCoeff decode_ac_coeff(const HuffmanTable& table) {
		int16_t rs = decode_symbol(table);
		if (rs == -1) throw std::runtime_error("decode_ac_coeff: failed to decode rs symbol");

		AcCoeff coeff;

		coeff.r = (0xF0 & rs) >> 4;
		coeff.s = 0x0F & rs;

		if (coeff.s > 0) {
			int16_t v = next_n_bits(coeff.s);
			if (v < 0) throw std::runtime_error("decode_ac_coeff: failed to decode coeff value");

			coeff.value = extend(v, coeff.s);
		}

		return coeff;
	}
};

class JpegDecoder {
private:
	std::vector<RawComponent> decode_frame(const CompressedJpegData& data, const Frame& frame) {
		std::vector<RawComponent> components;

		int hmax = frame.h_max();
		int vmax = frame.v_max();

		msg::debug(
			"DECODE: start frame: hmax={} vmax={} mculength={}",
			hmax, vmax, frame.mcu_length()
		);

		msg::debug("DECODE: init raw components...");
		int cnt = 0;
		for (const auto& comp_info : frame.component_params) {
			float h_ratio = ((float) comp_info.horizontal_sampling_factor) / hmax;
			float v_ratio = ((float) comp_info.vertical_sampling_factor) / vmax;

			RawComponent comp = RawComponent(
				std::ceil(frame.samples_per_line * h_ratio),
				std::ceil(frame.num_lines * v_ratio),
				comp_info.horizontal_sampling_factor,
				comp_info.vertical_sampling_factor,
				0
			);

			msg::debug("DECODE: comp{}: {}", cnt, std::string(comp));

			components.push_back(std::move(comp));
			cnt++;
		}

		cnt = 0;
		for (const auto& scan : frame.scans) {
			msg::debug(
				"DECODE: scan {}: RSTi={}",
				cnt, scan.restart_interval
			);
			decode_scan(data, scan, components);
			cnt++;
		}

		msg::debug("DECODE: decode_frame: done");

		return std::move(components);
	}

	void decode_scan(
		const CompressedJpegData& data,
		const Scan& scan,
		std::vector<RawComponent>& components
	) {
		std::vector<BlockView> block_views;
		auto scan_view = ScanDecodeView(scan);

		msg::debug("DECODE: decode_scan: init blockviews");

		// Initialize blockviews and DC predictors
		for (auto& comp : components) {
			block_views.push_back(BlockView(comp));
		}

		msg::debug("DECODE: decode_scan: begin decode");

		// Decode MCUs until done.
		int mcucnt = 0;
		while (decode_mcu(data, scan_view, block_views)) {
			mcucnt++;

			if (scan.restart_interval > 0 && mcucnt >= scan.restart_interval) {
				scan_view.next_rst();
				mcucnt = 0;
			}
		}
	}

	bool decode_mcu(
		const CompressedJpegData& data,
		ScanDecodeView& scan_view,
		std::vector<BlockView>& block_views
	) {
		//msg::debug("DECODE: start MCU");

		int cnt = 0;
		for (auto& block : block_views) {
			const ScanComponentParams& params = scan_view.scan().component_params[cnt];
			int acsel = params.ac_entropy_coding_selector;
			int dcsel = params.dc_entropy_coding_selector;
			int qsel = data.frames()[0].component_params[cnt].qtable_selector;

			scan_view.select_component(cnt);

			int blockcnt = 0;
			do {
				bool could_decode = false;
				try {
					could_decode = decode_block(
						data.q_tables()[qsel],
						data.huff_tables(TableClass::dc)[dcsel],
						data.huff_tables(TableClass::ac)[acsel],
						scan_view,
						block,
						cnt
					);
					//throw std::runtime_error("stop");
				} catch (const std::exception& ex) {
					msg::error("DECODE: failed to decode block: {}", ex.what());
					msg::error(
						"DECODE: ...at block_view{}: {}",
						cnt, std::string(block)
					);

					throw ex;
				}

				if (!could_decode) {
					throw std::runtime_error("could not decode mcu");
				}

				blockcnt++;
			} while (block.block_advance());

			//msg::debug("DECODE: blockcnt={}", blockcnt);

			bool could_advance = block.mcu_advance();
			if (!could_advance) {
				msg::debug("DECODE: decode_mcu: reached end of destination data");
				return false;
			}

			cnt++;
		}

		//throw std::runtime_error("stop");
		return true;
	}

	bool decode_block(
		const QuantizationTable& q_tbl,
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl,
		ScanDecodeView& scan_view,
		BlockView& block_view,
		int cnt
	) {
		Matrix<int16_t, 8, 8> tmp_block;

		// Read block into temp space
		tmp_block[0, 0] = scan_view.decode_dc_coeff(dc_tbl);

		int zz = 1;
		while (zz <= 63) {
			AcCoeff ac = scan_view.decode_ac_coeff(ac_tbl);
			if (ac.is_eob()) {
				break;
			}

			// Skip R zero coefficients
			zz += ac.r;
			if (zz >= 64) break;

			auto& valref = tmp_block.zz(zz);
			valref = ac.value;
			zz++;
		}

		// Undo quantization, then IDCT
		auto detransform = dct_mat * (tmp_block.ptwise_mul(q_tbl.data)) * dct_mat_t;
		for (size_t y = 0; y < 8; y++) {
			for (size_t x = 0; x < 8; x++) {
				double level = detransform[x, y] + 128;
				if (level < 0) level = 0;
				if (level > 255) level = 255;

				block_view[x, y] = (uint8_t) level;
			}
		}

		return true;
	}
public:
	std::vector<RawComponent> decode(const CompressedJpegData& data) {
		msg::info("Decoding compressed data");

		if (data.frames().size() == 0) {
			throw std::runtime_error("Cannot decode() data with 0 frames");
		} else if (data.frames().size() > 1) {
			throw std::runtime_error("For now, multiple frames are unsupported.");
		}

		const Frame& frame = data.frames()[0];
		return decode_frame(data, frame);
	}
};
}
