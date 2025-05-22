export module jpeg:codingbase;

// codingbase.cpp
// Defines common code shared between the encode and decode processes.

import std;
import msg;
import :data.jpeg;
import :data.image;

export namespace jpeg {
// Encapsulates various debug settings for coding processes.
template<typename ScanView>
struct CodingDebugConfig {
	// Limits debug features to specific scans. If this set is empty,
	// all scans will have debug enabled (mcu domain will still be
	// considered).
	std::set<unsigned int> enable_scans;

	// Minimum and maximum MCU values to enable debug features for
	unsigned int enable_mcu_x_min = 0;
	unsigned int enable_mcu_x_max = std::numeric_limits<int>::max();
	unsigned int enable_mcu_y_min = 0;
	unsigned int enable_mcu_y_max = std::numeric_limits<int>::max();

	// Debug print toggles
	bool print_ac_coeff = false;
	bool print_correction = false;
	bool print_blockn = false;

	// File output toggles (WARNING: these files might be very large)
	bool save_dct_coeff = false;

	// Checks if debug is enabled for a given scan.
	bool is_enabled_for_scan(const Scan& scan) const {
#ifndef DEBUG
		return false;
#endif

		if (enable_scans.size() > 0) {
			auto idx = static_cast<unsigned int>(scan.index);
			if (!enable_scans.contains(idx)) {
				return false;
			}
		}

		return true;
	}

	bool is_enabled_for_scan(const ScanView& scan_view) const {
		return is_enabled_for_scan(scan_view.scan());
	}

	// Checks if debug is enabled for a pair of views (present during
	// coding processes).
	bool is_enabled_for_views(
		const BlockView<ComponentBuffer>& block_view,
		const ScanView& scan_view
	) const {
#ifndef DEBUG
		return false;
#endif
		if (block_view.mcu_coords().x < enable_mcu_x_min) return false;
		if (block_view.mcu_coords().x > enable_mcu_x_max) return false;
		if (block_view.mcu_coords().y < enable_mcu_y_min) return false;
		if (block_view.mcu_coords().y > enable_mcu_y_max) return false;
		if (!is_enabled_for_scan(scan_view)) return false;

		return true;
	}
};

// Concept describing a coding policy.
template<typename T, typename ScanView>
concept is_coding_policy = requires(
	const HuffmanTable& dc_tbl,
	const HuffmanTable& ac_tbl,
	CodingDebugConfig<ScanView> debug_config,
	ScanView& scan_view,
	BlockView<ComponentBuffer>& block_view,
	bool all_run
) {
	// En/decodes a block using the provided tables.
	// Encode processes will read from block_view and write to scan_view,
	// decode processes will read from scan_view and write to block_view.
	{ 
		T::code_block_progressive_ac(debug_config, block_view, scan_view, ac_tbl, all_run) 
	} -> std::unsigned_integral;
	{
		T::code_block_progressive_dc(debug_config, block_view, scan_view, dc_tbl)
	} -> std::unsigned_integral;
	{
		T::code_block_sequential(debug_config, block_view, scan_view, dc_tbl, ac_tbl)
	} -> std::unsigned_integral;

	// The name for the phase, to be used in debug messages and exceptions.
	// Should be all caps, something like "ENCODE" or "DECODE".
	{ T::phase_name } -> std::convertible_to<std::string>;

	// Flushes the scan view. Only used in case of encode processes.
	T::flush(scan_view);
};

// Base class for a container of DC predictors. DC coefficients are coded
// differentially; to facilitate this, a DcPredictor keeps track of the last
// coded DC coefficient in each component.
class DcPredictor {
private:
	std::vector<int16_t> dc_preds;
	size_t cur_dc_pred;

protected:
	void dc_reset() {
		for (auto& p : dc_preds) p = 0;
	}

	int16_t& dc_pred() {
		return dc_preds[cur_dc_pred];
	}

public:
	DcPredictor(size_t num_components) {
		for (int i = 0; i < num_components; i++) {
			dc_preds.push_back(0);
		}
	}

	// En/decoder doesn't know what component it's writing to, so
	// need to select it from outside. FIXME may be able to eliminate this
	void select_component(size_t c) {
		cur_dc_pred = c;
	}
};

// Helper class consisting of a collection of BlockViews into the various
// components of an image buffer.
class BlockViewBundle {
private:
	std::vector<BlockView<ComponentBuffer>> block_views;
public:
	BlockViewBundle(const Frame& frame, const Scan& scan, ImageBuffer& buf) {
		// Bundles are created based on a JPEG scan. Some scans won't
		// contain data for all components in the frame, so only create
		// BlockViews for the ones that are actually present in the
		// scan.
		for (unsigned int comp = 0; comp < scan.num_components; comp++) {
			const auto& params = frame.get_component_params_by_scan(scan, comp);
			auto& comp_buf = buf[params.index];

			Dimensions d;

			// In one-component scans, MCU is always 1x1.
			if (scan.num_components == 1) {
				d.x = 1;
				d.y = 1;
			} else {
				d.x = params.h();
				d.y = params.v();
			}

			block_views.emplace(block_views.end(), comp_buf, d);
		}
	}

	BlockView<ComponentBuffer>& operator[](size_t n) {
		return block_views[n];
	}

	size_t size() {
		return block_views.size();
	}

	// Advances every blockview in this bundle at once. Returns true if at
	// least one view could advance, false if none of them could.
	//
	// If some views advanced but others did not, this indicates that some
	// views reached the end of their components too early. This most
	// likely indicates a broken JPEG, but it's a tolerable failure, so
	// just warn about it.
	bool mcu_advance() {
		int advance_cnt = 0;
		int noadvance_cnt = 0;

		for (auto& block : block_views) {
			if (block.mcu_advance()) {
				advance_cnt++;
			} else {
				noadvance_cnt++;
			}
		}

		if (advance_cnt > 0 && noadvance_cnt > 0) {
			msg::warn("BlockViewBundle::mcu_advance: mixed result!");
		}

		return advance_cnt != 0;
	}
};

// Class encapsulating some common code needed to en/decode scans.
template<typename CodingPolicy, typename ScanView, typename ScanType>
requires is_coding_policy<CodingPolicy, ScanView>
class ScanCoder {
private:
	const CodingDebugConfig<ScanView>* debug_config;
	const CompressedJpegData* data;
	const Frame* frame;
	ScanType* scan;
	ImageBuffer* buf;

	// If >0, instructs code_mcu_noninterleaved to "skip" this number of
	// blocks. Unused by the interleaved process. "Skipping" a block means
	// that said block contains no newly defined coefficients. However, the
	// scan may contain correction bits in the successive approximation
	// case, so we still need to call CodingPolicy::code_block_progressive_ac
	// even if the block is to be skipped.
	unsigned int blocks_to_skip = 0;

	template<typename... Args>
	unsigned int code_block_with_diagnostic(
		auto blockop,
		unsigned int component_index,
		BlockView<ComponentBuffer>& block_view,
		ScanView& scan_view,
		Args&&... args
	) {
		try {
			return (*blockop)(*debug_config, block_view, scan_view, std::forward<Args>(args)...);
		} catch (const std::exception& ex) {
			msg::error(
				"{}: failed to code block: {}",
				CodingPolicy::phase_name, ex.what()
			);
			msg::error(
				"{}: ...at block_views[{}]: {}",
				CodingPolicy::phase_name,
				component_index,
				std::string(block_view)
			);

			throw;
		}
	}

	// Codes an MCU in a noninterleaved scan. For our purposes (YCbCr
	// images), a noninterleaved scan will only contain AC coefficients.
	void code_mcu_noninterleaved(
		ScanView& scan_view,
		BlockViewBundle& block_views
	) {
		const auto& ac_huff = frame->get_huff_table(*scan, 0, TableClass::ac);
		bool skip_block = blocks_to_skip > 0;

		// No need to iterate over blocks, since in noninterleaved
		// scans 1 MCU is 1 block.
		unsigned int ret_skipcount = code_block_with_diagnostic(
			&CodingPolicy::code_block_progressive_ac,
			0,
			block_views[0],
			scan_view,
			ac_huff,
			skip_block
		);

		if (skip_block) {
			--blocks_to_skip;
		} else if (ret_skipcount > 0) {
			//msg::debug("nonzero ret_skipcount ({})", ret_skipcount);
			blocks_to_skip = ret_skipcount;
		}

		// TMP
		//throw std::runtime_error("stop");
	}

	// Codes a single component in an interleaved scan.
	void code_mcu_interleaved_component(
		ScanView& scan_view,
		BlockViewBundle& block_views,
		unsigned int comp_index
	) {
		auto& block = block_views[comp_index];
		const auto& dc_huff = frame->get_huff_table(*scan, comp_index, TableClass::dc);

		std::function<void(void)> op;

		// kinda ugly but we'll see how performance is
		if (frame->type.mode == OperationMode::progressive) {
			op = [&]() {
				code_block_with_diagnostic(
					&CodingPolicy::code_block_progressive_dc,
					comp_index,
					block,
					scan_view,
					dc_huff
				);
			};
		} else {
			const auto& ac_huff = frame->get_huff_table(*scan, comp_index, TableClass::ac);
			op = [&]() {
				code_block_with_diagnostic(
					&CodingPolicy::code_block_sequential,
					comp_index,
					block,
					scan_view,
					dc_huff,
					ac_huff
				);
			};
		}

		do {
			op();
		} while (block.block_advance());
	}

	// Codes an MCU in an interleaved scan.
	void code_mcu_interleaved(
		ScanView& scan_view,
		BlockViewBundle& block_views
	) {
		for (unsigned int comp = 0; comp < block_views.size(); ++comp) {
			scan_view.select_component(comp);
			code_mcu_interleaved_component(scan_view, block_views, comp);
		}
	}

	// Selects the proper MCU coding procedure to use.
	auto select_mcu_coding_procedure() {
		bool interleaved = (
			frame->type.is_baseline() ||
			scan->start_spectral_sel == 0 && scan->end_spectral_sel == 0
		);

		if (interleaved) {
			msg::debug("{}: using interleaved process", CodingPolicy::phase_name);
			return &ScanCoder::code_mcu_interleaved;
		} else {
			msg::debug("{}: using noninterleaved process", CodingPolicy::phase_name);
			return &ScanCoder::code_mcu_noninterleaved;
		}
	}
public:
	ScanCoder(
		const CodingDebugConfig<ScanView>& _debug_config,
		const CompressedJpegData& _data,
		const Frame& _frame,
		ScanType& _scan,
		ImageBuffer& _buf
	) : debug_config{&_debug_config}, data{&_data}, frame{&_frame}, scan{&_scan}, buf{&_buf} {}

	void code_scan() {
		msg::debug("{}: code_scan: init data views", CodingPolicy::phase_name);
		bool progressive_ac = (
			scan->num_components == 1 &&
			frame->type.mode == OperationMode::progressive
		);

		ScanView scan_view { *scan, progressive_ac };
		BlockViewBundle block_views { *frame, *scan, *buf };

		msg::debug("{}: code_scan: begin coding process", CodingPolicy::phase_name);
		auto code_mcu = select_mcu_coding_procedure();

		int mcucnt = 0;
		blocks_to_skip = 0; // just in case it's left over from prev. scan
		bool could_advance = false;
		do {
			(this->*code_mcu)(scan_view, block_views);
			++mcucnt;

			if (scan->restart_interval > 0 && mcucnt >= scan->restart_interval) {
				scan_view.next_rst();
				mcucnt = 0;
			}

			could_advance = block_views.mcu_advance();
			if (!could_advance) {
				msg::debug(
					"{}: code_scan: reached end of destination data",
					CodingPolicy::phase_name
				);
			}
		} while (could_advance);

		CodingPolicy::flush(scan_view);
	}
};
}
