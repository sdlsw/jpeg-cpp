module jpeg:codingbase;

// codingbase.cpp
// Defines common code shared between the encode and decode processes.

import std;
import msg;
import :data.jpeg;
import :data.image;

namespace jpeg {
// Concept describing a coding policy.
template<typename T, typename ScanView>
concept is_coding_policy = requires(
	const HuffmanTable& dc_tbl,
	const HuffmanTable& ac_tbl,
	ScanView& scan_view,
	BlockView<ComponentBuffer>& block_view
) {
	// En/decodes a block using the provided tables.
	// Encode processes will read from block_view and write to scan_view,
	// decode processes will read from scan_view and write to block_view.
	T::code_block(dc_tbl, ac_tbl, scan_view, block_view);

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

// Base class for en/decoders. Much of the logic regarding iterating through
// the blocks is the same between the two.
template<typename CodingPolicy, typename ScanView>
requires is_coding_policy<CodingPolicy, ScanView>
class CodingBase {
protected:
	// en/decodes a scan. ScanType is templated to allow for
	// const/non-const.
	template<typename ScanType>
	static void code_scan(
		const CompressedJpegData& data,
		const Frame& frame,
		ScanType& scan,
		ImageBuffer& buf
	) {
		msg::debug("{}: code_scan: init data views", CodingPolicy::phase_name);
		ScanView scan_view { scan };
		BlockViewBundle block_views { frame, scan, buf };

		msg::debug("{}: code_scan: begin coding process", CodingPolicy::phase_name);

		// Code MCUs until done.
		int mcucnt = 0;
		bool could_advance = false;
		do {
			code_mcu(data, frame, scan_view, block_views);
			mcucnt++;

			if (scan.restart_interval > 0 && mcucnt >= scan.restart_interval) {
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

	// en/decodes an MCU. Calls CodingPolicy::code_block once per every
	// block in the MCU windows of each component.
	static void code_mcu(
		const CompressedJpegData& data,
		const Frame& frame,
		ScanView& scan_view,
		BlockViewBundle& block_views
	) {
		const Scan& scan = scan_view.scan();

		for (int comp = 0; comp < block_views.size(); comp++) {
			auto& block = block_views[comp];
			scan_view.select_component(comp);

			const auto& dc_huff = frame.get_huff_table(scan, comp, TableClass::dc);
			const auto& ac_huff = frame.get_huff_table(scan, comp, TableClass::ac);

			do {
				code_block_with_diagnostic(
					dc_huff,
					ac_huff,
					scan_view,
					block,
					comp
				);
			} while (block.block_advance());
		}
	}

	// en/decodes a block, but prints some debug information on failure
	// before propagating the exception
	static void code_block_with_diagnostic(
		const HuffmanTable& dc_tbl,
		const HuffmanTable& ac_tbl,
		ScanView& scan_view,
		BlockView<ComponentBuffer>& block_view,
		unsigned int component
	) {
		try {
			CodingPolicy::code_block(
				dc_tbl,
				ac_tbl,
				scan_view,
				block_view
			);
		} catch (const std::exception& ex) {
			msg::error(
				"{}: failed to code block: {}",
				CodingPolicy::phase_name, ex.what()
			);
			msg::error(
				"{}: ...at block_views[{}]: {}",
				CodingPolicy::phase_name,
				component,
				std::string(block_view)
			);

			throw;
		}
	}
};
}
