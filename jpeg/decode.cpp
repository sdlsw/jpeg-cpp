export module jpeg:decode;

// decode.cpp:
// Implements the JPEG decoding procedure.

import std;
import msg;
import :codingbase;
import :data.image;
import :data.jpeg;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

export namespace jpeg {
// Helper class that keeps track of all the indices needed to traverse
// entropy-coded data segments for decoding.
class ScanDecodeView : public DcPredictor {
private:
	const Scan* _scan;         // The scan we are decoding.
	unsigned int rst_idx = 0;  // Current interval in scan.
	unsigned int byte_idx = 0; // Current byte in current RST interval.
	uint8_t bit_idx = 0;       // Current bit in current byte.
	uint8_t cur_byte = 0;      // Current byte being decoded.

	// Checks if all data in the scan has been exhausted.
	bool is_end_of_intervals() const {
		return rst_idx >= _scan->entropy_coded_data.size();
	}

	// Gets a reference to the RST interval currently being decoded.
	const std::vector<uint8_t>& cur_interval() const {
		if (is_end_of_intervals()) {
			throw std::out_of_range("unexpected end of scan data");
		}

		return _scan->entropy_coded_data[rst_idx];
	}

	// Reads the next byte in the current interval. Returns true if there is
	// now a byte available to read, and false if there is none.
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

	// Reads the next bit from the bitstream. Returns -1 if no bit could be
	// read. Otherwise, returns 0x1 for a 1 bit, 0x0 for a 0 bit.
	int8_t next_bit() {
		if (bit_idx >= 8) {
			if (!next_byte()) return -1;
		}

		uint8_t bit = (cur_byte >> (7 - bit_idx)) & 0x1;

		bit_idx++;
		return bit;
	}

	// Sign-extends a T-bit coded coefficient value V.
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

	// Decodes a huffman-coded symbol. The symbol itself will be 8 bits;
	// a special value of -1 indicates no value could be read.
	int16_t decode_symbol(const HuffmanTable& table) {
		uint16_t code = 0;
		int8_t bit = 0;
		int16_t value = 0;

		// Huffman codes are variable length, and there is no indication
		// of that length in huffman-coded data. So we need to
		// progressively read one bit at a time and check if the
		// resulting values match anything.
		//
		// The huffman codes used in JPEG compression are generated in
		// such a way that no code is ever the prefix for a longer
		// code. No need to worry about accidentally returning the
		// wrong value.
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
public:
	// Advances to the next RST interval. Returns true if successfully
	// advanced to next interval, false otherwise.
	bool next_rst() {
		rst_idx++;
		if (is_end_of_intervals()) {
			msg::debug("DECODE: reached end of intervals");
			return false;
		}

		// Read first byte in new interval
		byte_idx = 0;
		if(!next_byte()) return false;

		// DC predictor is reset on each RST interval.
		dc_reset();

		return true;
	}

	ScanDecodeView(const Scan& scan) : DcPredictor(scan.num_components), _scan(&scan) {
		// Consume first byte right away so next_bit has something to
		// read.
		next_byte();
	}

	// Allow const access to scan to enable reading dimension info.
	const Scan& scan() const { return *_scan; }

	// Decodes the next symbol as a DC DCT coefficient.
	int16_t decode_dc_coeff(const HuffmanTable& table) {
		int16_t& pred = dc_pred();
		int16_t ssss = decode_symbol(table);

		if (ssss < 0) {
			throw std::runtime_error("decode_dc_coeff: failed to decode ssss symbol");
		}

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

	// Decodes the next symbol as an AC DCT coefficient and run length.
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

// Implementation of the is_coding_policy concept for normal decoding processes.
class DecodePolicy {
public:
	static inline const std::string phase_name = "DECODE";

	static void code_block(
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl,
		ScanDecodeView& scan_view,
		BlockView<ComponentBuffer>& block_view
	) {
		auto& dcval = block_view.zz(0);
		dcval = scan_view.decode_dc_coeff(dc_tbl);

		int zz = 1;
		while (zz <= 63) {
			AcCoeff ac = scan_view.decode_ac_coeff(ac_tbl);
			if (ac.is_eob()) {
				break;
			}

			// Skip R zero coefficients
			zz += ac.r;
			if (zz >= 64) break;

			auto& valref = block_view.zz(zz);
			valref = ac.value;
			zz++;
		}
	}

	// Flush only required on encode
	static void flush(ScanDecodeView& scan_view) {}
};

// Class encapsulating the JPEG decoding algorithm. This decoder supports only
// baseline huffman coding.
class JpegDecoder : CodingBase<DecodePolicy, ScanDecodeView> {
private:
	// Decodes a single JPEG frame. CompressedJpegData is passed as well
	// for access to the tables needed for decoding.
	ImageBuffer decode_frame(
		const CompressedJpegData& data,
		const Frame& frame
	) {
		int hmax = frame.h_max();
		int vmax = frame.v_max();

		msg::debug(
			"DECODE: start frame: hmax={} vmax={}",
			hmax, vmax
		);

		msg::debug("DECODE: init DCT buffers...");

		ImageBuffer buf { SampleFormat::DctCoeff };
		for (const auto& comp_info : frame.component_params) {
			float h_ratio = ((float) comp_info.h()) / hmax;
			float v_ratio = ((float) comp_info.v()) / vmax;

			Dimensions d {
				static_cast<size_t>(std::ceil(frame.samples_per_line * h_ratio)),
				static_cast<size_t>(std::ceil(frame.num_lines * v_ratio))
			};

			const auto& comp = buf.new_component(d);

			msg::debug(
				"DECODE: comp{}: {}",
				buf.num_components() - 1,
				std::string(comp)
			);
		}

		msg::debug("DECODE: decode_frame: Decode scans...");
		int cnt = 0;
		for (const auto& scan : frame.scans) {
			msg::debug(
				"DECODE: scan {}: RSTi={}",
				cnt, scan.restart_interval
			);
			code_scan(data, frame, scan, buf);
			cnt++;
		}

		msg::debug("DECODE: decode_frame: IDCT...");
		buf.idct(data.get_qtables(frame));
		msg::debug("DECODE: decode_frame: done");

		return std::move(buf);
	}

public:
	// Decompress the data contained in a JPEG object.
	ImageBuffer decode(const CompressedJpegData& data) {
		msg::info("Decoding compressed data");

		if (data.frames().size() == 0) {
			throw std::runtime_error("Cannot decode() data with 0 frames");
		} else if (data.frames().size() > 1) {
			throw std::runtime_error("For now, multiple frames are unsupported.");
		}

		const Frame& frame = data.frames()[0];

		if (!frame.type.is_baseline()) {
			throw std::runtime_error("For now, only baseline encoding is supported.");
		}

		return decode_frame(data, frame);
	}
};
}
