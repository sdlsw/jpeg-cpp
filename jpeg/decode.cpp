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

	// In progressive mode, AC scans behave a little differently.
	// EOB codes are extended to EOBn, which can specify large runs of zero
	// blocks with only correction bits.
	bool _progressive_ac = false;

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
			if (byte_idx >= cur_interval().size()) {
				msg::error("DECODE: next_byte: ran off end of data during byte stuffing handler");
				return false;
			}
			next_byte = cur_interval()[byte_idx];
		}

		cur_byte = next_byte;
		byte_idx++;
		bit_idx = 0;

		return true;
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

	ScanDecodeView(
		const Scan& scan,
		bool progressive_ac = false
	) : DcPredictor(scan.num_components),
	    _scan(&scan),
	    _progressive_ac(progressive_ac)
	{
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
			// normal coefficient
			int16_t v = next_n_bits(coeff.s);
			if (v < 0) throw std::runtime_error("decode_ac_coeff: failed to decode coeff value");

			coeff.value = extend(v, coeff.s);
		} else if (_progressive_ac && coeff.r < 15) {
			// coeff.s == 0 as well
			// indicates EOBn, which is parsed slightly differently
			// from normal coefficients
			int16_t v = next_n_bits(coeff.r);
			if (v < 0) throw std::runtime_error("decode_ac_coeff: failed to decode runlength");

			coeff.value = v | (1 << coeff.r);
		}

		return coeff;
	}
};

// Implementation of the is_coding_policy concept for normal decoding processes.
class DecodePolicy {
public:
	static inline const std::string phase_name = "DECODE";

	static unsigned int code_block_progressive_dc(
		const CodingDebugConfig<ScanDecodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanDecodeView& scan_view,
		const HuffmanTable& dc_tbl
	) {
		const Scan& scan = scan_view.scan();
		auto& valref = block_view.zz(0);
		auto al = scan.succ_approx_low;

		if (scan.is_first_in_band()) {
			valref = scan_view.decode_dc_coeff(dc_tbl) << al;
		} else {
			// subsequent scan
			int8_t b = scan_view.next_bit();
			if (b < 0) throw std::runtime_error(
				"ran off end of data while in correction phase!"
			);
			valref |= b << al;
		}

		return 0;
	}

	// Helper function for code_block_progressive_ac. Handles zero runs
	// from normal coefficient symbols, ZRL symbols and block skips from
	// EOBn symbols.
	//
	// runlength - The length of the run. Meaning is dependent on value of
	//      `count_nonzero`.
	// zz - Mutable reference to the ZZ index from
	//      code_block_progressive_ac incremented according to how many
	//      coefficients are modified and how long the run is.
	// count_nonzero - When true, coefficients that have history (i.e.
	//      defined in a previous scan so they're nonzero) will decrement
	//      runlength.
	//
	//      If false, then coefficients with history will not
	//      count towards the length of the run, and will only accept
	//      correction bits. Note that this means that if runlength reaches
	//      zero, prog_ac_handle_run will continue accepting correction
	//      bits for nonzero coefficients until the first zero coefficient
	//      is encountered, or the block ends - whichever comes first. This
	//      behavior is consistent with the specification of symbols with
	//      R<=15.
	//
	// stop_when_runlength_zero - If true, then prog_ac_handle_run stops
	//      once runlength hits zero. This is to accomodate ZRL symbols,
	//      which ***DO NOT*** accept more correction bits once runlength
	//      reaches zero, but they DO still accept correction bits
	//      otherwise, and those corrected coefficients don't count towards
	//      the 16 runlength.
	//
	//      If false, then the time to stop is left up to the value of
	//      `count_nonzero`.
	static void prog_ac_handle_run(
		const CodingDebugConfig<ScanDecodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanDecodeView& scan_view,
		long runlength,
		unsigned int& zz,
		bool count_nonzero = true,
		bool stop_when_runlength_zero = false
	) {
		const Scan& scan = scan_view.scan();
		auto al = scan.succ_approx_low;
		auto se = scan.end_spectral_sel;
		bool print_correction = (
			debug_config.is_enabled_for_views(block_view, scan_view) &&
			debug_config.print_correction
		);

		// First scan in a band has no correction bits, no need to
		// waste time checking coefficient history.
		if (scan.is_first_in_band()) {
			zz += runlength;
			return;
		}

		auto rl0 = runlength;

		// FIXME ugly, break this into multiple functions maybe?
		while (true) {
			if (zz > se) {
				if (runlength > 0) {
					msg::warn("Ran off end of block handling run (rl0={}, rl_now={})", rl0, runlength);
					msg::debug("block: {}", std::string(block_view));
				}
				break;
			}

			if (stop_when_runlength_zero && runlength == 0) {
				break;
			}

			auto& valref = block_view.zz(zz);

			if (valref != 0) {
				// The coefficient has been defined
				// in a previous scan, so accept a correction
				// bit and add it to the coefficient.
				int16_t b = static_cast<int16_t>(scan_view.next_bit());
				if (b < 0) throw std::runtime_error(
					"ran off end of data while in correction phase!"
				);

				if (print_correction) {
					msg::debug("  correction at zz={}: {}", zz, b);
				}

				if (valref < 0) {
					valref -= b << al;
				} else {
					valref += b << al;
				}

				if (count_nonzero) {
					if (runlength <= 1) break;
					--runlength;
				}
			} else {
				if (runlength == 0) break;

				// Always decrement runlength on zero-history
				// coefficients
				--runlength;
			}

			++zz;
		}
	}

	static unsigned int code_block_progressive_ac(
		const CodingDebugConfig<ScanDecodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanDecodeView& scan_view,
		const HuffmanTable& ac_tbl,
		bool all_run = false
	) {
		const Scan& scan = scan_view.scan();
		auto ss = scan.start_spectral_sel;
		auto se = scan.end_spectral_sel;
		auto al = scan.succ_approx_low;

		bool debug = debug_config.is_enabled_for_views(block_view, scan_view);
		bool print_ac_coeff = debug && debug_config.print_ac_coeff;
		bool print_blockn = debug && debug_config.print_blockn;

		unsigned int zz = ss;

		if (print_blockn) {
			msg::debug("block: {}", block_view.blockn());
		}

		if (all_run) {
			if (print_ac_coeff) {
				msg::debug("skipping block, no zz");
			}
			prog_ac_handle_run(
				debug_config,
				block_view,
				scan_view,
				(se - ss) + 1,
				zz
			);
			return 0;
		}

		while (zz <= se) {
			AcCoeff coeff = scan_view.decode_ac_coeff(ac_tbl);

			if (print_ac_coeff) {
				msg::debug("zz={}, coeff: {}", zz, std::string(coeff));
			}

			if (coeff.is_zrl()) {
				// TODO: Would be better to count number
				// of ZRL symbols and handle the combined
				// run length all at once when a normal symbol
				// is encountered. I tried this, but it
				// doesn't work for reasons I don't understand
				// yet. Will figure it out later.
				prog_ac_handle_run(
					debug_config,
					block_view,
					scan_view,
					16,
					zz,
					false,
					true
				);
			} else if (!coeff.is_eob()) {
				// normal coeff, <=15 run length
				prog_ac_handle_run(
					debug_config,
					block_view,
					scan_view,
					coeff.r,
					zz,
					false
				);

				if (print_ac_coeff) {
					msg::debug("  note: zz after run: {}", zz);
				}

				if (!scan.is_first_in_band() && std::abs(coeff.value) != 1) {
					throw std::runtime_error("Bad approximation coefficient");
				}

				if (zz <= se) {
					auto& valref = block_view.zz(zz);
					valref = coeff.value << al;
					++zz;
				}
			} else {
				// EOBn
				prog_ac_handle_run(
					debug_config,
					block_view,
					scan_view,
					(se - zz),
					zz
				);
				return coeff.value - 1;
			}
		}

		return 0;
	}

	// Sequential algorithm is much simpler and can skip a lot of the work
	// progressive mode does, so have a separate function for it.
	static unsigned int code_block_sequential(
		const CodingDebugConfig<ScanDecodeView>& debug_config,
		BlockView<ComponentBuffer>& block_view,
		ScanDecodeView& scan_view,
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl
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

		return 0;
	}

	// Flush only required on encode
	static void flush(ScanDecodeView& scan_view) {}
};

using ScanDecoder = ScanCoder<DecodePolicy, ScanDecodeView, const Scan>;

// Class encapsulating the JPEG decoding algorithm. This decoder supports only
// baseline huffman coding.
class JpegDecoder {
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
		for (const auto& scan : frame.scans) {
			msg::debug(
				"DECODE: scan{}: {}",
				scan.index, scan.paramstr()
			);

			ScanDecoder decoder { debug_config, data, frame, scan, buf };
			decoder.code_scan();

			if (debug_config.is_enabled_for_scan(scan) && debug_config.save_dct_coeff) {
				buf.dump_to_file(std::format("scan{}.txt", scan.index));
			}
		}

		msg::debug("DECODE: decode_frame: IDCT...");
		buf.idct(data.get_qtables(frame));
		msg::debug("DECODE: decode_frame: done");

		return std::move(buf);
	}

	static bool frame_supported(const Frame& frame) {
		auto& t = frame.type;

		// baseline support is required
		if (t.is_baseline()) return true;

		// we support one progressive mode
		return (
			t.mode == OperationMode::progressive &&
			t.coding == CodingType::huffman &&
			!t.is_differential
		);
	}

public:
	CodingDebugConfig<ScanDecodeView> debug_config;

	// Decompress the data contained in a JPEG object.
	ImageBuffer decode(const CompressedJpegData& data) {
		msg::info("Decoding compressed data");

		if (data.frames().size() == 0) {
			throw std::runtime_error("Cannot decode() data with 0 frames");
		} else if (data.frames().size() > 1) {
			throw std::runtime_error("For now, multiple frames are unsupported.");
		}

		const Frame& frame = data.frames()[0];

		if (!frame_supported(frame)) {
			throw std::runtime_error("Frame coding mode unsupported.");
		}

		return decode_frame(data, frame);
	}
};
}
