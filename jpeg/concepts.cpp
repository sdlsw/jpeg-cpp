module jpeg:concepts;

import std;

// concepts.cpp:
// Internal common concepts used for templates in the jpeg module.

namespace jpeg {
	template<typename T>
	concept arithmetic = std::integral<T> or std::floating_point<T>;
};
