module;

// for max int size macros
#include <cstdint>

export module jpeg:data.image;

// data_image.cpp:
// Defines data structures for image buffers, and implements the
// forward DCT (FDCT) and inverse (IDCT) procedures.

import std;
import :concepts;
import :data.matrix;
import :data.tables;
import :data.vector;

using std::uint8_t;
using std::int16_t;

size_t rounded_div(size_t n, size_t d) {
	return std::lround(static_cast<double>(n) / static_cast<double>(d));
}

export namespace jpeg {
constexpr Dimensions block_dims {block_size, block_size};

// note: `zz_order_from_row_order` defined in `data.matrix`
constexpr Matrix<uint8_t, block_size, block_size> zzmat { zz_order_from_row_order };

constexpr auto generate_zz_coords(const Matrix<uint8_t, block_size, block_size>& m) {
	std::array<Point, block_dims.area()> pts;

	for (unsigned int y = 0; y < block_size; y++) {
		for (unsigned int x = 0; x < block_size; x++) {
			auto& pt = pts[zzmat[x, y]];
			pt.x = x;
			pt.y = y;
		}
	}

	return pts;
}

// Mapping of ZZ index to a pair of coordinates into an 8x8 block.
const std::array<Point, block_dims.area()> coords_from_zz = generate_zz_coords(zzmat);

enum class SampleFormat {
	DctCoeff, // DCT coefficient (16 bit signed integer)
	Sample8   // 8-bit unsigned image sample
};

const std::map<SampleFormat, std::string> str_from_sample_format {
	{SampleFormat::DctCoeff, "DctCoeff"},
	{SampleFormat::Sample8,  "Sample8"}
};

// A block view. Functions as a one-block window into a larger image buffer
// object.
template<typename BufferType>
class BlockView {
private:
	// Size of an MCU in blocks.
	Dimensions _mcu_size_blocks;

	// Size of the buffer in MCUs.
	Dimensions _buffer_size_mcus;

	// The current MCU window.
	Point _mcu_coords {0, 0};

	// The current block in the MCU window.
	Point _block_coords {0, 0};

	// X/Y offsets into the RawComponent this BlockView is viewing.
	// Determined by the MCU and block coordinates. They are only updated
	// on block_advance() and mcu_advance().
	Point _samp_offset {0, 0};

	// True if this BlockView is viewing a block that is in _buffer.
	// False if we're outside those bounds. Recalculated on block_advance()
	// and mcu_advance().
	bool _in_bounds = true;

	// The buffer being viewed.
	BufferType* _buffer;

	// Dummy space referenced by operator[] if the coordinates fall outside
	// valid image data.
	int16_t empty_samp = 0;

	void recalc_offset() {
		_samp_offset = (_mcu_coords * _mcu_size_blocks + _block_coords) * block_dims;
		_in_bounds = _buffer->backing_size_pixels().contains(_samp_offset);
	}

	// Advances a point through a space of given dimensions, in row-major
	// order.
	bool coords_advance(Point& coords, const Dimensions& dims) {
		bool done = (
			coords.y >= dims.height() ||
			(coords.y == dims.height() - 1 && coords.x >= dims.width() - 1)
		);

		if (done) return false;

		coords.x++;

		if (coords.x >= dims.width()) {
			coords.x = 0;
			coords.y++;
		}

		recalc_offset();

		return true;
	}
public:
	BlockView() = delete;
	BlockView(
		BufferType& buf,
		const Dimensions& mcu_size_blocks
	) : 
		_buffer{&buf},
		_mcu_size_blocks{mcu_size_blocks},
		_buffer_size_mcus{
			buf.backing_size_blocks().region_divide(mcu_size_blocks)
		}
		{}

	// If no MCU size is specified, we set the MCU size to be the whole
	// buffer. This way, we can repeatedly call block_advance() to
	// iterate over the entire buffer.
	BlockView(BufferType& buf) : BlockView(buf, buf.backing_size_blocks()) {}

	int16_t& operator[](const Point& p) {
#ifdef DEBUG
		block_dims.bounds_check(p);
#endif

		// Sample coordinates may fall outside of component buffer,
		// if it is not MCU-aligned. In this case, we return a
		// reference to a dummy value.
		if (!_in_bounds) {
			empty_samp = 0;
			return empty_samp;
		}

		return (*_buffer)[_samp_offset + p];
	}

	int16_t& operator[](size_t x, size_t y) {
		Point p { x, y };
		return (*this)[p];
	}

	// Accesses this block in zig-zag order.
	int16_t& zz(this auto& self, size_t z) {
#ifdef DEBUG
		if (z >= block_dims.area()) throw std::out_of_range("zz: out of range");
#endif
		return self[coords_from_zz[z]];
	}

	// Allow read-only access to the underlying buffer for access to
	// dimension information.
	const BufferType& buffer() const { return *_buffer; }

	// Advances this view to the next MCU, resetting the block coordinates
	// to (0, 0). If we are able to advance, this will return true. Failure
	// to advance will return false, indicating the view is at the end of
	// the component's data.
	bool mcu_advance() {
		_block_coords.x = 0;
		_block_coords.y = 0;
		return coords_advance(_mcu_coords, _buffer_size_mcus);
	}

	// Advances this window to the next block within the current MCU.
	// If we are able to advance, this will return true. Failure to advance
	// will return false, indicating we have run out of blocks in the
	// current MCU. In this case, mcu_advance() must be called.
	bool block_advance() {
		return coords_advance(_block_coords, _mcu_size_blocks);
	}

	// Copies the data at the current block to a matrix.
	void copy_to_matrix(Matrix<double>& mat) {
		for (size_t y = 0; y < block_size; y++) {
			for (size_t x = 0; x < block_size; x++) {
				mat[x, y] = (double) (*this)[x, y];
			}
		}
	}

	// Copies the data in a matrix to the current block.
	void copy_from_matrix(const Matrix<double>& mat) {
		for (size_t y = 0; y < block_size; y++) {
			for (size_t x = 0; x < block_size; x++) {
				(*this)[x, y] = (int16_t) mat[x, y];
			}
		}
	}

	// Calculates the inverse DCT of the current block. This operation is
	// in-place, and will overwrite the data currently at this block.
	void idct(const QuantizationTable& table) {
		Matrix<double> mat;
		copy_to_matrix(mat);

		mat.ptwise_mul_inplace(table.data);
		auto detransform = dct_mat * mat * dct_mat_t;

		detransform += 128.0;
		detransform.clamp(0.0, 255.0);
		copy_from_matrix(detransform);
	}

	// Calculates the forward DCT of the current block. This operation is
	// in-place.
	void fdct(const QuantizationTable& table) {
		Matrix<double> mat;
		copy_to_matrix(mat);
		mat -= 128.0;

		auto transform = dct_mat_t * mat * dct_mat;
		transform.ptwise_div_inplace(table.data);

		copy_from_matrix(transform);
	}

	explicit operator std::string() const {
		return std::format(
			"mcu={}, block={}, offset={}",
			std::string(_mcu_coords),
			std::string(_block_coords),
			std::string(_samp_offset)
		);
	}
};

// A buffer for uncompressed image data for one image component. The image may
// be stored as either DCT coefficients, or 8 bit detransformed samples.
//
// FDCT/IDCT may be performed in-place, since an 8x8 DCT block represents
// an 8x8 block of pixels.
class ComponentBuffer {
private:
	// Width and height of valid component data.
	Dimensions _valid_size;

	// Width and height of backing storage in blocks and pixels.
	Dimensions _backing_size_blocks;
	Dimensions _backing_size_pixels;

	// Format of samples in _data
	SampleFormat _sampfmt;

	// Backing image storage. This is a two-dimensional array stored in
	// row-major order. This storage is allocated once on construction
	// and never resized. int16 is enough to store any of the possible
	// sample formats.
	std::vector<int16_t> _data;

	size_t data_idx_from_pt(const Point& p) {
#ifdef DEBUG
		_backing_size_pixels.bounds_check(p);
#endif
		return p.y * _backing_size_pixels.width() + p.x;
	}

	// Performs a DCT operation on this component. A DCT operation is a member
	// function pointer from BlockView.
	template<typename DctOp>
	void do_dct_op(
		const QuantizationTable& table,
		SampleFormat required_fmt,
		SampleFormat target_fmt,
		DctOp op
	) {
		if (_sampfmt != required_fmt) {
			throw std::runtime_error(std::format(
				"_sampfmt must be {}",
				str_from_sample_format.at(required_fmt)
			));
		}

		BlockView<ComponentBuffer> view { *this };

		do {
			(view.*op)(table);
		} while (view.block_advance());

		_sampfmt = target_fmt;
	}
public:
	ComponentBuffer() = delete;
	ComponentBuffer(
		const Dimensions& image_size,
		SampleFormat sampfmt,
		int16_t default_value
	) : _sampfmt{sampfmt}, _valid_size{image_size} {
		_backing_size_blocks = image_size.region_divide(block_dims);
		_backing_size_pixels = _backing_size_blocks * block_size;

		_data.insert(_data.end(), _backing_size_pixels.area(), default_value);
	}
	ComponentBuffer(
		const Dimensions& image_size,
		SampleFormat sampfmt
	) : ComponentBuffer(image_size, sampfmt, 0) {}

	const auto& backing_size_blocks() const { return _backing_size_blocks; }
	const auto& backing_size_pixels() const { return _backing_size_pixels; }
	const auto& valid_size() const { return _valid_size; }
	auto sample_format() const { return _sampfmt; }

	decltype(auto) operator[](this auto& self, const Point& p) {
		return self._data[self.data_idx_from_pt(p)];
	}

	decltype(auto) operator[](this auto& self, unsigned int x, unsigned int y) {
		Point p { x, y };
		return self[p];
	}

	// Performs the inverse DCT on this component in-place, setting the
	// sample format to Sample8.
	void idct(const QuantizationTable& table) {
		do_dct_op(
			table,
			SampleFormat::DctCoeff,
			SampleFormat::Sample8,
			&BlockView<ComponentBuffer>::idct
		);
	}

	// Performs the forward DCT on this component in-place, setting the
	// sample format to DctCoeff.
	void fdct(const QuantizationTable& table) {
		do_dct_op(
			table,
			SampleFormat::Sample8,
			SampleFormat::DctCoeff,
			&BlockView<ComponentBuffer>::fdct
		);
	}

	explicit operator std::string() const {
		return std::format(
			"valid_size={}, backing_size_blocks={}, backing_size_pixels={}",
			std::string(_valid_size),
			std::string(_backing_size_blocks),
			std::string(_backing_size_pixels)
		);
	}
};

// A buffer for uncompressed image data for one or more image components. The
// data may be stored as either DCT coefficients, or 8 bit detransformed
// samples.
class ImageBuffer {
private:
	std::vector<ComponentBuffer> comps;
	SampleFormat _sampfmt;

	// Performs a DCT operation on this image. A DCT operation is a
	// member function pointer from ComponentBuffer.
	template<typename DctOp>
	void do_dct_op(
		const std::vector<const QuantizationTable*>& tables,
		SampleFormat target_fmt,
		DctOp op
	) {
		if (tables.size() != comps.size()) {
			throw std::invalid_argument("tables/comp size() mismatch");
		}

		// No need to check sample format here, first component
		// will do the check and throw if it's wrong.
		for (unsigned int i = 0; i < comps.size(); i++) {
			(comps[i].*op)(*(tables[i]));
		}

		_sampfmt = target_fmt;
	}
public:
	ImageBuffer() = delete;
	ImageBuffer(SampleFormat fmt) : _sampfmt{fmt} {}

	uint8_t num_components() const {
		return static_cast<uint8_t>(comps.size());
	}

	auto sample_format() { return _sampfmt; }

	ComponentBuffer& new_component(const Dimensions& dims) {
		if (comps.size() >= UINT8_MAX) {
			throw std::runtime_error("Images cannot have more than 255 components");
		}
		comps.emplace(comps.end(), dims, _sampfmt);
		return comps.back();
	}

	ComponentBuffer& operator[](size_t n) {
		return comps[n];
	}

	// Calculates the subsampling factor in the x and y direction, for
	// component n, assuming component 0 is the largest.
	Vec2<size_t> subsamp(size_t n) const {
		if (n == 0) return Vec2<size_t>(1, 1);

		const auto& a = comps[0].valid_size();
		const auto& b = comps[n].valid_size();

		return Vec2<size_t>(
			rounded_div(a.width(), b.width()),
			rounded_div(a.height(), b.height())
		);
	}

	// Runs the forward DCT procedure on this image in-place.
	void fdct(const std::vector<const QuantizationTable*>& tables) {
		do_dct_op(tables, SampleFormat::DctCoeff, &ComponentBuffer::fdct);
	}

	// Runs the inverse DCT procedure on this image in-place.
	void idct(const std::vector<const QuantizationTable*>& tables) {
		do_dct_op(tables, SampleFormat::Sample8, &ComponentBuffer::idct);
	}

	auto begin() {
		return comps.begin();
	}

	auto end() {
		return comps.end();
	}
};
}
