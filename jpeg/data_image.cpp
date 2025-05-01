export module jpeg:data.image;

// data_image.cpp:
// Defines data structures for decoded image data.

import std;

using std::uint8_t;

const unsigned int block_size = 8;

export namespace jpeg {
// Represents an uncompressed image component.
class RawComponent {
private:
	// Width and height of valid component data.
	int _x_pixels;
	int _y_pixels;

	// Width and height of backing storage in MCU (Minimum Coded Units)
	// windows.
	int _x_mcus;
	int _y_mcus;

	// Width and height of one MCU window, in blocks. Blocks are always 8x8
	// pixels.
	int _mcu_width;
	int _mcu_height;

	// Maximum bounds on the backing store. This may be a bit
	// larger than _x_pixels and/or _y_pixels if those values aren't
	// integer multiples of the MCU window size in pixels.
	// These values are calculated from other parameters and cached since
	// they are checked on every sample access.
	int _x_size;
	int _y_size;

	// Backing image storage. This is a two-dimensional array stored in
	// row-major order. This storage is allocated once on construction
	// and never resized.
	std::vector<uint8_t> _data;

	// Calculates how many MCU windows wide/tall the backing storage should
	// be.
	static int calc_n_mcus(int samples, int mcu_len_in_blocks) {
		auto [q, r] = std::div(samples, block_size*mcu_len_in_blocks);
		return q + ((r == 0) ? 0 : 1);
	}

	// Calculates the total width of the backing array.
	int calc_x_size() const {
		return _x_mcus * _mcu_width * block_size;
	}

	// Calculates the total height of the backing array.
	int calc_y_size() const {
		return _y_mcus * _mcu_height * block_size;
	}

	// Calculates the index into the backing store given (x, y)
	// coordinates. If either coordinate is out of range, throws
	// std::out_of_range.
	int data_idx_from_xy(int x, int y) const {
		if (x >= _x_size) {
			throw std::out_of_range("RawComponent: x value out of range");
		}

		if (y >= _y_size) {
			throw std::out_of_range("RawComponent: y value out of range");
		}

		return y*_x_size + x;
	}
public:
	RawComponent() = delete;

	RawComponent(
		int x_pixels,
		int y_pixels,
		int mcu_width_blocks,
		int mcu_height_blocks,
		uint8_t default_value
	) {
		// All of the calculated sizes are cached on
		// construction, since we're going to be using them a LOT.
		_x_pixels = x_pixels;
		_y_pixels = y_pixels;

		_mcu_width = mcu_width_blocks;
		_mcu_height = mcu_height_blocks;

		_x_mcus = calc_n_mcus(_x_pixels, mcu_width_blocks);
		_y_mcus = calc_n_mcus(_y_pixels, mcu_height_blocks);

		_x_size = calc_x_size();
		_y_size = calc_y_size();

		size_t vecsize = _x_size * _y_size;
		_data = std::vector<uint8_t>(vecsize, default_value);
	}

	RawComponent(
		int x_pixels,
		int y_pixels,
		int mcu_width_blocks,
		int mcu_height_blocks
	) : RawComponent(x_pixels, y_pixels, mcu_width_blocks, mcu_height_blocks, 0) {}

	int x_pixels() const { return _x_pixels; }
	int y_pixels() const { return _y_pixels; }

	int x_mcus() const { return _x_mcus; }
	int y_mcus() const { return _y_mcus; }

	int mcu_width() const { return _mcu_width; }
	int mcu_height() const { return _mcu_height; }

	int x_size() const { return _x_size; }
	int y_size() const { return _y_size; }

	// Direct sample access. Access is allowed to go slightly beyond
	// img_width/img_height, to simplify cases where the img size and
	// block size are misaligned.
	decltype(auto) operator[](this auto& self, int x, int y) {
		return self._data[self.data_idx_from_xy(x, y)];
	}

	explicit operator std::string() const {
		return std::format(
			"xpixels={} ypixels={} mcu_width={} mcu_height={} xmcus={} ymcus={}",
			_x_pixels, _y_pixels, _mcu_width, _mcu_height, _x_mcus, _y_mcus
		);
	}
};

// Helper class for accessing one block of image data according to the scan
// block ordering rules. Acts as a window into a larger RawComponent object.
class BlockView {
private:
	// MCU coordinates.
	int _mcu_x;
	int _mcu_y;

	// Block coordinates.
	int _block_x;
	int _block_y;

	// X/Y offsets into the RawComponent this BlockView is viewing.
	// Determined by the MCU and block coordinates. They are only updated
	// on block_advance() and mcu_advance().
	int _x_offset;
	int _y_offset;

	// The component being viewed.
	RawComponent* _component;

	static void coord_check(int x, int y) {
		if (x >= block_size) {
			throw std::out_of_range("BlockView: x value out of range");
		}

		if (y >= block_size) {
			throw std::out_of_range("BlockView: y value out of range");
		}
	}

	bool coord_advance(
		int& x,
		int& y,
		const int& x_bound,
		const int& y_bound
	) {
		if (y >= y_bound || (y == y_bound - 1 && x >= x_bound - 1)) {
			return false;
		}

		bool recalc_y = false;

		x++;

		if (x >= x_bound) {
			x = 0;
			y++;
		}

		recalc_x_offset();
		recalc_y_offset();

		return true;
	}

	// Recalculate the X/Y offsets into the RawComponent associated with
	// this view.
	void recalc_x_offset() {
		_x_offset = (_mcu_x * _component->mcu_width() + _block_x) * block_size;
	}
	void recalc_y_offset() {
		_y_offset = (_mcu_y * _component->mcu_height() + _block_y) * block_size;
	}
public:
	BlockView() = delete;
	BlockView(RawComponent& comp) : BlockView(0, 0, 0, 0, comp) {}

	// Construct a BlockView with starting MCU and block coordinates.
	BlockView(
		int mcu_x,
		int mcu_y,
		int block_x,
		int block_y,
		RawComponent& comp
	) {
		if (mcu_x >= comp.x_mcus()) {
			throw std::out_of_range("mcu x value out of range");
		}

		if (mcu_y >= comp.y_mcus()) {
			throw std::out_of_range("mcu y value out of range");
		}

		if (block_x >= comp.mcu_width()) {
			throw std::out_of_range("block x value out of range");
		}

		if (block_y >= comp.mcu_height()) {
			throw std::out_of_range("block y value out of range");
		}

		_component = &comp;

		_mcu_x = mcu_x;
		_mcu_y = mcu_y;
		_block_x = block_x;
		_block_y = block_y;

		recalc_x_offset();
		recalc_y_offset();
	}

	// Accesses the sample at (x, y) of the current block pointed to by
	// this view, depending on the current block and MCU coordinates.
	uint8_t& operator[](int x, int y) {
		coord_check(x, y);
		return (*_component)[_x_offset + x, _y_offset + y];
	}

	// Allow read-only access to the underlying component for access
	// to dimension information.
	const RawComponent& component() const { return *_component; }

	int mcu_x() { return _mcu_x; }
	int mcu_y() { return _mcu_y; }

	// Advances this view to the next MCU, resetting the block coordinate
	// to (0, 0). If we are able to advance, this will return true. Failure
	// to advance will return false, indicating the view is at the end of
	// the component's data.
	bool mcu_advance() {
		_block_x = 0;
		_block_y = 0;

		return coord_advance(
			_mcu_x,
			_mcu_y,
			_component->x_mcus(),
			_component->y_mcus()
		);
	}

	// Advances this window to the next block within the current MCU.
	// If we are able to advance, this will return true. Failure to advance
	// will return false, indicating we have run out of blocks in the
	// current MCU. In this case, mcu_advance() must be called.
	bool block_advance() {
		return coord_advance(
			_block_x,
			_block_y,
			_component->mcu_width(),
			_component->mcu_height()
		);
	}

	explicit operator std::string() const {
		return std::format(
			"mcux={} mcuy={} blockx={} blocky={} offset=<{}, {}>",
			_mcu_x, _mcu_y, _block_x, _block_y, _x_offset, _y_offset
		);
	}
};
}
