module jpeg:codingbase;

import std;
import msg;

import :data;

namespace jpeg {
// Concept describing a coding policy.
template<typename T, typename ScanView>
concept is_coding_policy = requires(
	const QuantizationTable& q_tbl,
	const HuffmanTable& dc_tbl,
	const HuffmanTable& ac_tbl,
	ScanView& scan_view,
	BlockView& block_view
) {
	// En/decodes a block using the provided tables.
	// Encode processes will read from block_view and write to scan_view,
	// decode processes will read from scan_view and write to block_view.
	T::code_block(q_tbl, dc_tbl, ac_tbl, scan_view, block_view);

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

// Helper class consisting of a collection of BlockViews into a collection
// of components.
class BlockViewBundle {
private:
	std::vector<BlockView> block_views;
public:
	BlockViewBundle(std::vector<RawComponent>& components) {
		for (auto& comp : components) {
			block_views.push_back(BlockView(comp));
		}
	}

	BlockView& operator[](size_t n) {
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
		ScanType& scan,
		std::vector<RawComponent>& components
	) {
		msg::debug("{}: code_scan: init data views", CodingPolicy::phase_name);
		auto scan_view = ScanView(scan);
		auto block_views = BlockViewBundle(components);

		msg::debug("{}: code_scan: begin coding process", CodingPolicy::phase_name);

		// Code MCUs until done.
		int mcucnt = 0;
		bool could_advance = false;
		do {
			code_mcu(data, scan_view, block_views);
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
		ScanView& scan_view,
		BlockViewBundle& block_views
	) {
		for (int comp = 0; comp < block_views.size(); comp++) {
			auto& block = block_views[comp];

			const ScanComponentParams& params = scan_view.scan().component_params[comp];
			int acsel = params.ac_entropy_coding_selector;
			int dcsel = params.dc_entropy_coding_selector;
			int qsel = data.frames()[0].component_params[comp].qtable_selector;

			scan_view.select_component(comp);

			do {
				code_block_with_diagnostic(
					data.q_tables()[qsel],
					data.huff_tables(TableClass::dc)[dcsel],
					data.huff_tables(TableClass::ac)[acsel],
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
		const QuantizationTable& q_tbl,
		const HuffmanTable& ac_tbl,
		const HuffmanTable& dc_tbl,
		ScanView& scan_view,
		BlockView& block_view,
		unsigned int component
	) {

		try {
			CodingPolicy::code_block(
				q_tbl,
				ac_tbl,
				dc_tbl,
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
