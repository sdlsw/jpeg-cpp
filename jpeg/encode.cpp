export module jpeg:encode;

import std;
import msg;

import :codingbase;
import :data;
import :tables;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export namespace jpeg {
// Helper class that keeps track of all indices needed to entropy-code
// jpeg data.
class ScanEncodeView : public DcPredictor {
private:
	Scan* _scan;
	uint8_t n_bits = 0;
	uint8_t cur_byte = 0xFF; // bytes are padded with 1 bits

	void put_byte() {
		auto& interval = _scan->entropy_coded_data.back();

		interval.push_back(cur_byte);
		if (cur_byte == 0xFF) interval.push_back(0x00);

		cur_byte = 0xFF;
		n_bits = 0;
	}

	void put_bit(uint8_t b) {
		if (n_bits == 8) put_byte();
		uint8_t mask = 0x1 << (7 - n_bits);
		cur_byte &= ~mask; // clear bit
		cur_byte |= b << (7 - n_bits);
		n_bits++;
	}

	void put_n_bits(uint16_t v, uint8_t n) {
		if (n == 0) return;
		if (n > 16) throw std::out_of_range("put_n_bits: n > 16");

		// TODO inefficient?
		for (unsigned int i = 0; i < n; i++) {
			uint16_t mask = 0x1 << (n-i-1);
			put_bit((v & mask) >> (n-i-1));
		}
	}

	static uint8_t magnitude_category(int16_t v) {
		if (v == 0) return 0;

		v = std::abs(v);
		uint16_t mask = 0xFFFF;
		uint8_t magnitude = 0;

		while((mask & v) != 0) {
			mask <<= 1;
			magnitude++;
		}

		return magnitude;
	}
public:
	ScanEncodeView(Scan& scan) : DcPredictor(scan.num_components), _scan(&scan) {
		next_rst();
	}

	const Scan& scan() const { return *_scan; }

	void next_rst() {
		_scan->entropy_coded_data.emplace(_scan->entropy_coded_data.end());
		dc_reset();
	}

	void encode_dc_coeff(const HuffmanTable& table, int16_t coeff) {
		int16_t& pred = dc_pred();
		int16_t diff = coeff - pred;
		pred = coeff;

		uint8_t ssss = magnitude_category(diff);
		const HuffmanCode& code = table.lookup_value(ssss);

		put_n_bits(code.code, code.bits);
		if (diff >= 0) {
			put_n_bits(diff, ssss);
		} else {
			put_n_bits(diff - 1, ssss);
		}
	}

	void encode_ac_coeff(const HuffmanTable& table, uint8_t run_length, int16_t coeff) {
		uint8_t s = magnitude_category(coeff);
		uint8_t rs = (run_length << 4) | s;

		const HuffmanCode& code = table.lookup_value(rs);

		put_n_bits(code.code, code.bits);
		if (coeff >= 0) {
			put_n_bits(coeff, s);
		} else {
			put_n_bits(coeff - 1, s);
		}
	}

	void encode_ac_eob(const HuffmanTable& table) {
		encode_ac_coeff(table, 0, 0);
	}
};

class EncodePolicy {
public:
	static inline const std::string phase_name = "ENCODE";

	static void code_block(
		const QuantizationTable& q_tbl,
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl,
		ScanEncodeView& scan_view,
		BlockView& block_view
	) {
		Matrix<double, 8, 8> tmp_block;

		// Read from raw component and perform level shift
		for (size_t y = 0; y < 8; y++) {
			for (size_t x = 0; x < 8; x++) {
				tmp_block[x, y] = ((double) block_view[x, y]) - 128.0;
			}
		}

		// DCT, then quantize
		auto transform = dct_mat_t * tmp_block * dct_mat;
		Matrix<int16_t, 8, 8> quantized = transform.ptwise_div<int16_t>(q_tbl.data);

		// Encode DC
		scan_view.encode_dc_coeff(dc_tbl, (int16_t) quantized[0, 0]);

		// Encode AC
		int zz = 1;
		uint8_t run_length = 0;
		while (zz <= 63) {
			auto& valref = quantized.zz(zz);

			if (valref != 0 | run_length >= 15) {
				scan_view.encode_ac_coeff(ac_tbl, run_length, valref);
				run_length = 0;
			} else {
				run_length++;
			}

			zz++;
		}

		// Explicitly code any remaining zeros as EOB
		// FIXME Recognize run lengths longer than 16, and encode only
		// one EOB if we reach the end of the block with a very long
		// run length.
		if (run_length > 0) {
			scan_view.encode_ac_eob(ac_tbl);
		}
	}
};

class JpegEncoder : CodingBase<ScanEncodeView, EncodePolicy> {
private:
	unsigned int quality = 50;

	static void add_luma_params(Frame& frame, Scan& scan) {
		FrameComponentParams frame_params;
		frame_params.identifier = 0;
		frame_params.horizontal_sampling_factor = 2;
		frame_params.vertical_sampling_factor = 2;
		frame_params.qtable_selector = 0;

		ScanComponentParams scan_params;
		scan_params.component_selector = 0;
		scan_params.dc_entropy_coding_selector = 0;
		scan_params.ac_entropy_coding_selector = 0;

		frame.component_params.push_back(frame_params);
		scan.component_params.push_back(scan_params);
	}

	static void add_chroma_params(Frame& frame, Scan& scan, uint8_t id) {
		FrameComponentParams frame_params;
		frame_params.identifier = id;
		frame_params.horizontal_sampling_factor = 1;
		frame_params.vertical_sampling_factor = 1;
		frame_params.qtable_selector = 1;

		ScanComponentParams scan_params;
		scan_params.component_selector = id;
		scan_params.dc_entropy_coding_selector = 1;
		scan_params.ac_entropy_coding_selector = 1;

		frame.component_params.push_back(frame_params);
		scan.component_params.push_back(scan_params);
	}

public:
	CompressedJpegData encode(std::vector<RawComponent>& raws) const {
		msg::debug("ENCODE: Encoding image with quality level {}", quality);

		CompressedJpegData j;

		// Init tables.
		j.q_tables()[0] = generate_qtable(qtable_base_luma, quality);
		msg::debug("ENCODE: luma qtable\n{}", std::string(j.q_tables()[0]));

		j.q_tables()[1] = generate_qtable(qtable_base_chroma, quality);
		msg::debug("ENCODE: chroma qtable\n{}", std::string(j.q_tables()[1]));

		j.huff_tables(TableClass::dc)[0] = hufftable_dc0;
		j.huff_tables(TableClass::ac)[0] = hufftable_ac0;
		j.huff_tables(TableClass::dc)[1] = hufftable_dc1;
		j.huff_tables(TableClass::ac)[1] = hufftable_ac1;

		j.frames().emplace(j.frames().end());
		Frame& frame = j.frames()[0];

		frame.num_lines = raws[0].y_pixels();
		frame.samples_per_line = raws[0].x_pixels();
		frame.sample_precision = 8;
		frame.num_components = raws.size();

		frame.scans.emplace(frame.scans.end());
		Scan& scan = frame.scans[0];

		scan.num_components = raws.size();
		scan.restart_interval = 0;

		add_luma_params(frame, scan);
		add_chroma_params(frame, scan, 1);
		add_chroma_params(frame, scan, 2);

		code_scan(j, scan, raws);

		j.set_valid();
		return std::move(j);
	}

	JpegEncoder() = default;
	JpegEncoder(unsigned int q) : quality{q} {}
};
}
