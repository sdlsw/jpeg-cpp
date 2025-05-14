export module jpeg:data.matrix;

import std;
import :concepts;

using std::uint8_t;

// data_matrix.cpp:
// Defines the matrix structure used for DCT transforms.

const unsigned int str_float_precision = 5;

export namespace jpeg {
constexpr unsigned int block_size = 8;
constexpr unsigned int block_area = block_size*block_size;

// For JPEG, we only encounter 8x8 blocks for DCT/quantization/etc.
// The standard document offers a table, so just hardcode it.
constexpr std::array<uint8_t, block_area> zz_order_from_row_order {
	0,   1,  5,  6, 14, 15, 27, 28,
	2,   4,  7, 13, 16, 26, 29, 42,
	3,   8, 12, 17, 25, 30, 41, 43,
	9,  11, 18, 24, 31, 40, 44, 53,
	10, 19, 23, 32, 39, 45, 52, 54,
	20, 22, 33, 38, 46, 51, 55, 60,
	21, 34, 37, 47, 50, 56, 59, 61,
	35, 36, 48, 49, 57, 58, 62, 63
};

constexpr auto invert_zzmap(const std::array<uint8_t, block_area>& a) {
	std::array<uint8_t, block_area> out;

	for (uint8_t i = 0; i < out.size(); i++) {
		out[a[i]] = i;
	}

	return out;
}

constexpr auto row_order_from_zz_order = invert_zzmap(zz_order_from_row_order);

template<arithmetic T, size_t W=block_size, size_t H=block_size>
class Matrix {
	// Every instance of Matrix must be friends with every other instance.
	template <arithmetic OTHER_T, size_t OTHER_W, size_t OTHER_H>
	friend class Matrix;
private:
	// FIXME big problem here - std::array isn't movable so all of the
	// arithmetic operations are copying entire matrices... 
	// UPDATE: In my measurements this doesn't actually seem to affect
	// performance much. TODO compare using std::vector
	std::array<T, W*H> data {0};

	constexpr void bounds_check(size_t x, size_t y) const {
		if (x >= W) throw std::out_of_range("x coordinate out of range");
		if (y >= H) throw std::out_of_range("y coordinate out of range");
	}

	template<arithmetic RESULT_T, arithmetic OTHER_T, size_t OTHER_W>
	constexpr RESULT_T calc_matmul_val(
		const Matrix<OTHER_T, OTHER_W, W>& other,
		size_t x,
		size_t y
	) const {
		// for self, start at leftmost col of row y and move
		// horizontally
		size_t self_cursor = y*W;

		// for other, start at top row of col x and move vertically
		size_t other_cursor = x;

		RESULT_T sum = 0;
		for (int i = 0; i < W; i++) {
			sum += (
				static_cast<RESULT_T>(other.data[other_cursor]) * 
				static_cast<RESULT_T>(data[self_cursor])
			);
			self_cursor++;
			other_cursor += W;
		}

		return sum;
	}

public:
	// for external access of dimensions
	static const size_t width = W;
	static const size_t height = H;

	Matrix() = default;

	constexpr Matrix(const std::initializer_list<T>& ilist) {
		if (ilist.size() != data.size()) {
			auto s = std::format(
				"Matrix: bad initializer list size {} (expected {})",
				ilist.size(),
				data.size()
			);
			throw std::invalid_argument(s);
		}

		auto iarray = ilist.begin();
		for (size_t i = 0; i < ilist.size(); i++) {
			data[i] = iarray[i];
		}
	}

	constexpr Matrix(const std::array<T, W*H>& a) {
		for (std::tuple<T&, const T&> t : std::views::zip(*this, a)) {
			std::get<0>(t) = std::get<1>(t);
		}
	}

	template<arithmetic OTHER_T>
	constexpr Matrix(const Matrix<OTHER_T, W, H>& other) {
		for (std::tuple<T&, const OTHER_T&> t : std::views::zip(*this, other)) {
			std::get<0>(t) = static_cast<T>(std::get<1>(t));
		}
	}

	constexpr auto begin() { return data.begin(); }
	constexpr auto end() { return data.end(); }

	constexpr auto begin() const { return data.cbegin(); }
	constexpr auto end() const { return data.cend(); }

	auto& operator[](this auto& self, size_t x, size_t y) {
#ifdef DEBUG
		self.bounds_check(x, y);
#endif
		return self.data[y*W + x];
	}

	// Accesses this matrix in zig-zag order. Only works for 8x8 matrices,
	// which is the size of a JPEG block.
	auto& zz(this auto& self, size_t zz) {
		if (H != block_size || W != block_size) {
			throw std::runtime_error("zz() only for 8x8 matrices");
		}

#ifdef DEBUG
		if (zz >= block_area) {
			throw std::invalid_argument("zz(): index cannot exceed 64");
		}
#endif

		return self.data[row_order_from_zz_order[zz]];
	}

	// Typical matrix multiplication. By default the number type of the
	// result will be the same as the left operand (A in A*B).
	// In order to manually specify the type of the resulting matrix, you
	// can specify with `A.operator*<desired_type>(B)`.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T, size_t OTHER_W>
	Matrix<RESULT_T, OTHER_W, H> operator*(const Matrix<OTHER_T, OTHER_W, W>& other) const {
		Matrix<RESULT_T, OTHER_W, H> result;

		for (size_t x = 0; x < OTHER_W; x++) {
			for (size_t y = 0; y < H; y++) {
				result[x, y] = calc_matmul_val<RESULT_T>(other, x, y);
			}
		}

		return std::move(result);
	}

	// Scalar matrix multiplication.
	Matrix<T, W, H> operator*=(T c) {
		for (auto& v : data) v *= c;
		return *this;
	}

	// Scalar addition. Adds c to every element of this matrix.
	Matrix<T, W, H>& operator+=(T c) {
		for (auto& v : data) v += c;
		return *this;
	}

	// Scalar subtraction. Subtracts c from every element of this matrix.
	Matrix<T, W, H>& operator-=(T c) {
		for (auto& v : data) v -= c;
		return *this;
	}

	// Creates a transposed copy of this matrix.
	Matrix<T, H, W> transposed() const {
		Matrix<T, H, W> result;

		for (size_t x = 0; x < W; x++) {
			for (size_t y = 0; y < H; y++) {
				result[y, x] = (*this)[x, y];
			}
		}

		return std::move(result);
	}

	template<arithmetic OTHER_T, typename ApplyOp>
	void ptwise_apply_inplace(
		const Matrix<OTHER_T, W, H>& other,
		ApplyOp op
	) {
		for (auto t : std::views::zip(*this, other)) {
			auto& a = std::get<0>(t);
			auto& b = std::get<1>(t);

			a = op(a, b);
		}
	}

	// Multiplies each point in this matrix by the corresponding point
	// in another of the same size.
	template<arithmetic OTHER_T>
	void ptwise_mul_inplace(const Matrix<OTHER_T, W, H>& other) {
		ptwise_apply_inplace(other, [](T a, OTHER_T b) {
			return a * static_cast<T>(b);
		});
	}

	// Same as ptwise_mul, but for division.
	template<arithmetic OTHER_T>
	void ptwise_div_inplace(const Matrix<OTHER_T, W, H>& other) {
		ptwise_apply_inplace(other, [](T a, OTHER_T b) {
			return a / static_cast<T>(b);
		});
	}

	// Forces every value in this matrix to be between `low` and `high`,
	// inclusive.
	void clamp(T low, T high) {
		for (auto& v : data) {
			if (v < low) {
				v = low;
			} else if (v > high) {
				v = high;
			}
		}
	}

	// Creates a clamped copy of this matrix, forcing every value `v`
	// to satisfy (low <= v <= high).
	Matrix<T, W, H> clamped(T low, T high) const {
		Matrix<T, W, H> result = *this;
		result.clamp(low, high);
		return result;
	}

	// Determines the maximum value in this matrix.
	T max() const {
		T m = 0;
		for (const auto& v : data) {
			if (v > m) m = v;
		}
		return m;
	}

	// Calculates how width each value in this matrix should be when
	// converting to string. Helps keep all the values lined up nicely.
	template<std::floating_point ST>
	size_t calc_field_width(this const Matrix<ST, W, H>& self) {
		return std::format("{:.{}f}", self.max(), str_float_precision).size() + 2;
	}

	// Same as above, specialized for integers.
	template<std::integral ST>
	size_t calc_field_width(this const Matrix<ST, W, H>& self) {
		T m = self.max();
		if (m == 0) return 1;

		return (size_t)std::log10(std::abs(self.max())) + 1;
	}

	// Formats a single matrix value (for float-valued matrices).
	template<std::floating_point ST>
	static std::string fmt_item(ST val, size_t field_width) {
		return std::format("{:{}.{}f}", val, field_width, str_float_precision);
	}

	// Same as above, specialized for integers.
	template<std::integral ST>
	static std::string fmt_item(ST val, size_t field_width) {
		return std::format("{:{}d}", val, field_width);
	}

	explicit operator std::string() const {
		size_t field_width = calc_field_width();

		std::stringstream out;

		unsigned int col = 0;
		for (T val : data) {
			out << fmt_item(val, field_width) << " ";
			
			col++;
			if (col >= W) {
				col = 0;
				out << std::endl;
			}
		}

		return out.str();
	}
};

// Calculates the DCT (discrete cosine transform) matrix. Used in the encoding
// and decoding processes.
auto calc_dct_mat() {
	Matrix<double> mat;

	double rt2 = std::sqrt(2);
	for (size_t y = 0; y < mat.height; y++) {
		for (size_t x = 0; x < mat.width; x++) {
			double c = (x == 0) ? (1/(2*rt2)) : 0.5;
			double cos_arg = (std::numbers::pi * x * (2*y + 1)) / 16;

			mat[x, y] = c * std::cos(cos_arg);
		}
	}

	return mat;
}

// Matrices used for the DCT transform.
const Matrix<double> dct_mat = calc_dct_mat();
const Matrix<double> dct_mat_t = dct_mat.transposed();

}
