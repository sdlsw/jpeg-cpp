export module jpeg:data.vector;

// data_vector.cpp:
// Classes representing pairs of numbers, and some derived types.

import std;
import :concepts;

export namespace jpeg {
template<arithmetic T>
struct Vec2 {
	T x = 0;
	T y = 0;

	constexpr Vec2() = default;
	constexpr Vec2(T _x, T _y) : x{_x}, y{_y} {}
	constexpr Vec2(const std::initializer_list<T>& l) {
		if (l.size() != 2) {
			throw std::invalid_argument(
				"Vec2: init list must have exactly two items"
			);
		}

		auto iter = l.begin();
		x = iter[0];
		y = iter[1];
	}

	// NOTE: Arithmetic operations on Vec2-derived types will always
	// return the same type, instead of decaying to Vec2.
	template<typename Self>
	Self& operator*=(this Self& self, const T& c) {
		self.x *= c;
		self.y *= c;
		return self;
	}

	template<typename Self>
	Self& operator*=(this Self& self, const Vec2& other) {
		self.x *= other.x;
		self.y *= other.y;
		return self;
	}

	template<typename Self, typename Rhs>
	Self operator*(this const Self& self, const Rhs& c) {
		Self v = self;
		v *= c;
		return v;
	}

	template<typename Self>
	Self& operator+=(this Self& self, const Vec2& other) {
		self.x += other.x;
		self.y += other.y;
		return self;
	}

	template<typename Self>
	Self operator+(this const Self& self, const Vec2& other) {
		Self v = self;
		v += other;
		return v;
	}

	explicit operator std::string() const {
		return std::format("({}, {})", x, y);
	}
};

// Represents a point within some Dimensions.
struct Point : public Vec2<size_t> {
	using Vec2::Vec2;
};

// Represents a 2D set of dimensions (width, height).
struct Dimensions : public Vec2<size_t> {
	using Vec2::Vec2;

	// Calculates minimum number of "segments" of size n required to
	// fit `dim`.
	static size_t calc_n_segments(size_t dim, size_t n) {
		auto [q, r] = std::lldiv(dim, n);
		return q + ((r == 0) ? 0 : 1);
	}

	constexpr size_t width() const { return x; }
	constexpr size_t height() const { return y; }

	constexpr size_t area() const {
		return width() * height();
	}

	// Divides these dimensions into regions of size given by 
	// `region_size`. The resulting dimensions will be the smallest number
	// of regions wide and tall a space needs to be to contain the region
	// described by `*this`. The overall size of the resulting dimensions
	// when multiplied by `region_size` may be larger than `*this`.
	Dimensions region_divide(const Dimensions& region_size) const {
		return {
			calc_n_segments(width(), region_size.width()),
			calc_n_segments(height(), region_size.height())
		};
	}

	// Checks if the point `p` falls within these dimensions.
	bool contains(const Point& p) const {
		return p.x < width() && p.y < height();
	}

	void bounds_check(const Point& p) const {
		if (!contains(p)) {
			throw std::out_of_range(std::format(
				"Point {} out of range of dims {}",
				std::string(p), std::string(*this)
			));
		}
	}
};
}
