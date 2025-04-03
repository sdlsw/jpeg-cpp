export module jpeg:encode;

import std;
import msg;

import :data;
import :tables;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export namespace jpeg {
// Helper class that keeps track of all indices needed to entropy-code
// jpeg data.
class ScanEncodeView {
private:
	Scan* _scan;
	uint8_t n_bits = 0;
	uint8_t cur_byte = 0xFF; // bytes are padded with 1 bits
	size_t cur_dc_pred = 0;
	std::vector<int16_t> dc_preds;

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
		uint8_t mask = 0xFF;
		uint8_t magnitude = 0;

		while((mask & v) != 0) {
			mask <<= 1;
			magnitude++;
		}

		return magnitude;
	}
public:
	ScanEncodeView(Scan& scan) : _scan(&scan) {
		for (int i = 0; i < _scan->num_components; i++) {
			dc_preds.push_back(0);
		}

		next_rst();
	}

	const Scan& scan() const { return *_scan; }

	void next_rst() {
		_scan->entropy_coded_data.emplace(_scan->entropy_coded_data.end());
		for (auto& pred : dc_preds) {
			pred = 0;
		}
	}

	void select_component(size_t c) {
		cur_dc_pred = c;
	}

	void encode_dc_coeff(const HuffmanTable& table, int16_t coeff) {
		int16_t diff = coeff - dc_preds[cur_dc_pred];
		dc_preds[cur_dc_pred] = coeff;

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

class JpegEncoder {
private:
	static void init_luma_params(Frame& frame, Scan& scan) {
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

	static void init_chroma_params(Frame& frame, Scan& scan, uint8_t id) {
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

	static void encode_scan(
		const CompressedJpegData& j,
		Scan& scan,
		std::vector<RawComponent>& raws
	) {
		std::vector<BlockView> block_views;
		auto scan_view = ScanEncodeView(scan);

		msg::debug("ENCODE: encode_scan: init blockviews");

		for (RawComponent& raw : raws) {
			block_views.push_back(BlockView(raw));
		}

		msg::debug("ENCODE: encode_scan: begin encode");

		int mcucnt = 0;
		while(encode_mcu(j, scan_view, block_views)) {
			mcucnt++;

			if (scan.restart_interval > 0 && mcucnt >= scan.restart_interval) {
				scan_view.next_rst();
				mcucnt = 0;
			}
		}
	}

	static bool encode_mcu(
		const CompressedJpegData& data,
		ScanEncodeView& scan_view,
		std::vector<BlockView>& block_views
	) {
		int cnt = 0;
		for (auto& block : block_views) {
			const ScanComponentParams& params = scan_view.scan().component_params[cnt];
			int acsel = params.ac_entropy_coding_selector;
			int dcsel = params.dc_entropy_coding_selector;
			int qsel = data.frames()[0].component_params[cnt].qtable_selector;

			scan_view.select_component(cnt);

			int blockcnt = 0;
			do {
				encode_block(
					data.q_tables()[qsel],
					data.huff_tables(TableClass::dc)[dcsel],
					data.huff_tables(TableClass::ac)[acsel],
					scan_view,
					block,
					cnt
				);
				//throw std::runtime_error("stop");
				blockcnt++;
			} while (block.block_advance());

			bool could_advance = block.mcu_advance();
			if (!could_advance) {
				msg::debug("ENCODE: encode_mcu: reached end of source data");
				return false;
			}

			cnt++;
		}
	}

	static void encode_block(
		const QuantizationTable& q_tbl,
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl,
		ScanEncodeView& scan_view,
		BlockView& block_view,
		int cnt
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
public:
	CompressedJpegData encode(std::vector<RawComponent>& raws) const {
		CompressedJpegData j;

		// Init tables.
		j.q_tables()[0] = qtable_base_luma;
		j.q_tables()[1] = qtable_base_chroma;
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

		init_luma_params(frame, scan);
		init_chroma_params(frame, scan, 1);
		init_chroma_params(frame, scan, 2);

		// FIXME there's no const constructor for BlockView,
		// so even though we never modify the raws here we can't
		// make it const
		encode_scan(j, scan, raws);

		j.set_valid();
		return std::move(j);
	}
};
}
