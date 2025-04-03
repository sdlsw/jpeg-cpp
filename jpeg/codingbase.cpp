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
}
