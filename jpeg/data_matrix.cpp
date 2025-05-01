export module jpeg:data.matrix;

import std;
using std::uint8_t;

// data_matrix.cpp:
// Defines the matrix structure used for DCT transforms.

// For JPEG, we only encounter 8x8 blocks for DCT/quantization/etc.
// The standard document offers a table, so just hardcode it.
const std::array<uint8_t, 64> zigzag_order_from_row_order = {
	0,   1,  5,  6, 14, 15, 27, 28,
	2,   4,  7, 13, 16, 26, 29, 42,
	3,   8, 12, 17, 25, 30, 41, 43,
	9,  11, 18, 24, 31, 40, 44, 53,
	10, 19, 23, 32, 39, 45, 52, 54,
	20, 22, 33, 38, 46, 51, 55, 60,
	21, 34, 37, 47, 50, 56, 59, 61,
	35, 36, 48, 49, 57, 58, 62, 63
};

// However, we can use the hardcoded table to generate the inverse one:
constexpr auto generate_order_inverse(const std::array<uint8_t, 64>& order) {
	std::array<uint8_t, 64> inverse;

	for (int i = 0; i < order.size(); i++) {
		inverse[order[i]] = i;
	}

	return std::move(inverse);
}

const std::array<uint8_t, 64> row_order_from_zigzag_order = generate_order_inverse(zigzag_order_from_row_order);
const unsigned int str_float_precision = 5;

template<typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

export namespace jpeg {
template<arithmetic T, size_t W=8, size_t H=8>
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
			sum += other.data[other_cursor] * data[self_cursor];
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
	Matrix(const std::initializer_list<T>& ilist) {
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

	decltype(auto) operator[](this auto& self, size_t x, size_t y) {
		self.bounds_check(x, y);
		return self.data[y*W + x];
	}

	// Shortcut access 'operator' for zig-zag ordering. Used only for 8x8
	// matrices.
	decltype(auto) zz(this auto& self, size_t z) {
		if (W != 8 || H != 8) throw std::runtime_error("zz only valid on 8x8 matrix");
		if (z >= 64) throw std::out_of_range("zz: out of range");

		return self.data[row_order_from_zigzag_order[z]];
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

	// Creates a transposed copy of this matrix.
	Matrix<T, H, W> transpose() const {
		Matrix<T, H, W> result;

		for (size_t x = 0; x < W; x++) {
			for (size_t y = 0; y < H; y++) {
				result[y, x] = (*this)[x, y];
			}
		}

		return std::move(result);
	}

	// Multiplies each point in this matrix by the corresponding point
	// in another of the same size.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T>
	Matrix<RESULT_T, W, H> ptwise_mul(const Matrix<OTHER_T, W, H>& other) const {
		Matrix<RESULT_T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			result.data[i] = this->data[i] * other.data[i];
		}

		return std::move(result);
	}

	// Same as ptwise_mul, but for division.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T>
	Matrix<RESULT_T, W, H> ptwise_div(const Matrix<OTHER_T, W, H>& other) const {
		Matrix<RESULT_T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			result.data[i] = this->data[i] / other.data[i];
		}

		return std::move(result);
	}

	// Creates a clamped copy of this matrix, forcing every value `v`
	// to satisfy (low <= v <= high).
	Matrix<T, W, H> clamp(T low, T high) const {
		Matrix<T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			T v = this->data[i];
			if (v < low) v = low;
			if (v > high) v = high;
			result.data[i] = v;
		}

		return std::move(result);
	}

	// Multiplies every value in this matrix by a single number.
	template<arithmetic RESULT_T=T, arithmetic OTHER_T>
	Matrix<RESULT_T, W, H> scalar_mul(OTHER_T m) const {
		Matrix<RESULT_T, W, H> result;

		for (size_t i = 0; i < data.size(); i++) {
			result.data[i] = this->data[i] * m;
		}

		return std::move(result);
	}

	// Determines the maximum value in this matrix.
	T max() const {
		T m = 0;
		for (size_t i = 0; i < data.size(); i++) {
			if (m < data[i]) m = data[i];
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
	Matrix<double, 8, 8> mat;

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
const Matrix<double, 8, 8> dct_mat = calc_dct_mat();
const Matrix<double, 8, 8> dct_mat_t = dct_mat.transpose();

}
