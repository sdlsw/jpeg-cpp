module jpeg:codingbase;

import std;
import msg;

import :data;

namespace jpeg {
// Base class for a container of DC predictors.
class DcPredictor {
private:
	// DC predictors stored in XXXScanView classes, since that's what decides
	// when to reset them.
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

	// Decoder doesn't know what component it's writing to, so
	// need to select it from outside. FIXME may be able to eliminate this
	void select_component(size_t c) {
		cur_dc_pred = c;
	}
};

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

		// return success if at least one block view could advance
		return advance_cnt != 0;
	}
};

// Base class for en/decoders. Much of the logic regarding iterating through
// the blocks is the same between the two.
template<typename ScanView, typename CodingPolicy>
class CodingBase {
protected:
	// en/decode a scan. ScanType is templated to allow for
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
	}

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
				try {
					CodingPolicy::code_block(
						data.q_tables()[qsel],
						data.huff_tables(TableClass::dc)[dcsel],
						data.huff_tables(TableClass::ac)[acsel],
						scan_view,
						block
					);
					//throw std::runtime_error("stop");
				} catch (const std::exception& ex) {
					msg::error(
						"{}: failed to code block: {}",
						CodingPolicy::phase_name, ex.what()
					);
					msg::error(
						"{}: ...at block_view{}: {}",
						CodingPolicy::phase_name, comp, std::string(block)
					);

					throw ex;
				}
			} while (block.block_advance());
		}
	}
};
}
