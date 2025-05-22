export module jpeg:encode;

// encode.cpp:
// Implements the JPEG encoding procedure.

import std;
import msg;

import :codingbase;
import :data.jpeg;
import :data.image;
import :tables;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export namespace jpeg {
// Helper class that keeps track of all indices needed to entropy-code
// JPEG data. Data is emitted as a bitstream of variable width codes. At the
// end of the encoding process, `flush()` MUST be called.
class ScanEncodeView : public DcPredictor {
private:
	// The scan we are encoding to.
	Scan* _scan;

	// A byte functioning as a buffer for the bitstream.
	// At the end of the encoding process, any "empty" bits in the last
	// byte are filled with 1s. Therefore, always initialize the buffer to
	// all 1s.
	uint8_t cur_byte = 0xFF;

	// Number of bits written to `cur_byte`.
	uint8_t n_bits = 0;

	// Writes `cur_byte` in the current restart interval and resets it.
	void put_byte() {
		auto& interval = _scan->entropy_coded_data.back();

		interval.push_back(cur_byte);

		// To avoid false markers getting generated in the output data,
		// stuff a 0 byte after every marker indicator.
		if (cur_byte == 0xFF) interval.push_back(0x00);

		cur_byte = 0xFF;
		n_bits = 0;
	}

	// Writes a bit (LSB of b) to the bitstream.
	void put_bit(uint8_t b) {
		uint8_t mask = 0x1 << (7 - n_bits);
		cur_byte &= ~mask; // clear bit
		cur_byte |= b << (7 - n_bits);
		n_bits++;

		if (n_bits == 8) put_byte();
	}

	// Writes the N lowest bits of V to the bitstream.
	void put_n_bits(uint16_t v, uint8_t n) {
		if (n == 0) return;
		if (n > 16) throw std::out_of_range("put_n_bits: n > 16");

		// TODO inefficient?
		for (unsigned int i = 0; i < n; i++) {
			uint16_t mask = 0x1 << (n-i-1);
			put_bit((v & mask) >> (n-i-1));
		}
	}

	// Calculates the magnitude category of a coefficient value V.
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
	ScanEncodeView(
		Scan& scan,
		bool progressive_ac = false // TODO unused for now
	) : DcPredictor(scan.num_components), _scan(&scan) {
		next_rst();
	}

	const Scan& scan() const { return *_scan; }

	// Finishes encoding the current RST interval in the scan, and advances
	// to a new one.
	void next_rst() {
		_scan->entropy_coded_data.emplace(_scan->entropy_coded_data.end());
		flush();
		dc_reset();
	}

	// Encodes a DC DCT coefficient to the scan, using the given table.
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

	// Encodes an AC coefficient symbol (consists of a run of zeros
	// followed by a non-zero AC DCT coefficient).
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

	// Encodes an EOB (end of block) symbol.
	void encode_ac_eob(const HuffmanTable& table) {
		encode_ac_coeff(table, 0, 0);
	}

	// Flushes any unwritten bits to the scan. Must be called manually after
	// the encode procedure is complete.
	void flush() {
		if (n_bits > 0) put_byte();
	}
};

// Implementation of the is_coding_policy concept for the normal encoding process.
class EncodePolicy {
public:
	static inline const std::string phase_name = "ENCODE";

	// Stub out progressive mode, not yet supported for encoding.
	static unsigned int code_block_progressive_dc(
		const CodingDebugConfig<ScanEncodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanEncodeView& scan_view,
		const HuffmanTable& dc_tbl
	) { return 0; }
	static unsigned int code_block_progressive_ac(
		const CodingDebugConfig<ScanEncodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanEncodeView& scan_view,
		const HuffmanTable& ac_tbl,
		bool all_run = false
	) { return 0; }

	// Codes an arbitrarily long run of zeros to the scan, followed by a
	// non-zero coefficient value.
	static void code_run(
		const HuffmanTable& ac_tbl,
		ScanEncodeView& scan_view,
		uint8_t run_length,
		int16_t coeff
	) {
		// Max number of zeros one symbol can code is 16.
		while (run_length >= 16) {
			scan_view.encode_ac_coeff(ac_tbl, 15, 0);
			run_length -= 16;
		}

		scan_view.encode_ac_coeff(ac_tbl, run_length, coeff);
	}

	static unsigned int code_block_sequential(
		const CodingDebugConfig<ScanEncodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanEncodeView& scan_view,
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl
	) {
		scan_view.encode_dc_coeff(dc_tbl, block_view.zz(0));

		// Encode AC
		int zz = 1;
		uint8_t run_length = 0;
		while (zz <= 63) {
			auto& valref = block_view.zz(zz);

			if (valref != 0) {
				code_run(ac_tbl, scan_view, run_length, valref);
				run_length = 0;
			} else {
				run_length++;
			}

			zz++;
		}

		// Explicitly code any remaining zeros as EOB
		if (run_length > 0) {
			scan_view.encode_ac_eob(ac_tbl);
		}

		return 0;
	}

	static void flush(ScanEncodeView& scan_view) {
		scan_view.flush();
	}
};

using ScanEncoder = ScanCoder<EncodePolicy, ScanEncodeView, Scan>;

// Class encapsulating the JPEG encoding algorithm. This encoder supports only
// baseline huffman coding, which is by far the most common operation mode in
// use on the internet today.
class JpegEncoder {
private:
	// The quality setting from 0-100 with which this encoder will compress
	// images. Lower values result in smaller file size, but worse visual
	// quality.
	unsigned int quality = 50;

	// Adds luma component parameters to a frame and associated scan.
	// The component ID will always be 0.
	static void add_luma_params(Frame& frame, Scan& scan, uint8_t subsamp) {
		FrameComponentParams& frame_params = frame.new_component();
		frame_params.identifier = 0;
		frame_params.horizontal_sampling_factor = subsamp;
		frame_params.vertical_sampling_factor = subsamp;
		frame_params.qtable_selector = 0;

		ScanComponentParams scan_params;
		scan_params.component_selector = 0;
		scan_params.dc_entropy_coding_selector = 0;
		scan_params.ac_entropy_coding_selector = 0;

		scan.component_params.push_back(scan_params);
	}

	// Adds chroma component parameters to a frame and associated scan.
	static void add_chroma_params(Frame& frame, Scan& scan, uint8_t id) {
		FrameComponentParams& frame_params = frame.new_component();
		frame_params.identifier = id;
		frame_params.horizontal_sampling_factor = 1;
		frame_params.vertical_sampling_factor = 1;
		frame_params.qtable_selector = 1;

		ScanComponentParams scan_params;
		scan_params.component_selector = id;
		scan_params.dc_entropy_coding_selector = 1;
		scan_params.ac_entropy_coding_selector = 1;

		scan.component_params.push_back(scan_params);
	}

public:
	CodingDebugConfig<ScanEncodeView> debug_config;

	// Compress the data contained in a raw YCrCb image. This will emit
	// one frame, consisting of one scan with all components interleaved.
	//
	// NOTE: This function has a side effect: `buf` will be converted
	// to DCT coefficients. If this is not desirable, make a copy of the
	// buffer to pass in.
	CompressedJpegData encode(ImageBuffer& buf) const {
		// Only YCrCb is supported, so reject anything that doesn't
		// have 3 components.
		if (buf.num_components() != 3) {
			throw std::invalid_argument("encode: can only encode 3-component images");
		}

		CompressedJpegData j;

		// Determine subsamp automatically
		// FIXME For now this is ok, since `buf` is only ever
		// generated by the BMP reader, but may produce unexpected
		// results for bundles of raw components generated in other
		// ways.
		uint8_t subsamp = static_cast<uint8_t>(buf.subsamp(1).x);
		msg::debug("ENCODE: Encoding image with quality {}, subsamp {}", quality, subsamp);

		// Init quantization tables.
		QuantizationTable lumatable = generate_qtable(qtable_base_luma, quality);
		msg::debug("ENCODE: luma qtable\n{}", std::string(lumatable));

		QuantizationTable chromatable = generate_qtable(qtable_base_chroma, quality);
		msg::debug("ENCODE: chroma qtable\n{}", std::string(chromatable));

		j.q_tables().install_table(0, std::move(lumatable));
		j.q_tables().install_table(1, std::move(chromatable));

		// Set up frame
		Frame& frame = j.new_frame();
		frame.num_lines = static_cast<uint16_t>(buf[0].backing_size_pixels().height());
		frame.samples_per_line = static_cast<uint16_t>(buf[0].backing_size_pixels().width());
		frame.sample_precision = 8;
		frame.num_components = buf.num_components();

		// Install pre-calculated huffman tables.
		frame.dc_huff_tables().install_table(0, hufftable_dc0);
		frame.dc_huff_tables().install_table(1, hufftable_dc1);
		frame.ac_huff_tables().install_table(0, hufftable_ac0);
		frame.ac_huff_tables().install_table(1, hufftable_ac1);

		// Set up scan
		Scan& scan = frame.new_scan();
		scan.num_components = buf.num_components();
		scan.restart_interval = 0;

		// Configure table selectors and sampling factors for each
		// component
		add_luma_params(frame, scan, subsamp);
		add_chroma_params(frame, scan, 1);
		add_chroma_params(frame, scan, 2);

		// Encode
		msg::debug("ENCODE: encode_frame: FDCT...");
		buf.fdct(j.get_qtables(frame));
		msg::debug("ENCODE: encode_frame: done.");

		ScanEncoder encoder { debug_config, j, frame, scan, buf };
		encoder.code_scan();

		j.set_valid();
		return std::move(j);
	}

	JpegEncoder() = default;
	JpegEncoder(unsigned int q) : quality{q} {}
};
}
