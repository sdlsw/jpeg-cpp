export module jpeg:file.bmp;

// file_bmp.cpp:
// Classes representing BMP files, and reader/writers for them.

import std;
import msg;
import :data.image;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;
using std::uint32_t;
using std::int32_t;

export namespace jpeg {
#pragma pack(push, 1)

// Represents the image header of a BMP file, specifically the Windows 3.x
// version. The defaults listed here are used to simplify writing in
// BmpFile::write().
//
// See https://www.fileformat.info/format/bmp/egff.htm
struct BmpImageHeader {
	uint32_t _header_size = sizeof(BmpImageHeader);
	uint32_t _width = 0;
	int32_t _height = 0;
	uint16_t _planes = 1;
	uint16_t _bit_count = 24;
	uint32_t _compression = 0;
	uint32_t _image_size = 0;
	uint32_t _x_pix_per_meter = 0;
	uint32_t _y_pix_per_meter = 0;
	uint32_t _color_used = 0;
	uint32_t _color_important = 0;

	BmpImageHeader() = default;
	BmpImageHeader(uint32_t width, uint32_t height) {
		_width = width;
		_height = height;
	}

	explicit operator std::string() const {
		return std::format(
			"_header_size={}\n"
			"_width={}\n"
			"_height={}\n"
			"_bit_count={}",
			_header_size,
			_width,
			_height,
			_bit_count
		);
	}
};

// Represents the initial file header of a BMP. The `_type` field must always
// be 'BM'.
struct BmpFileHeader {
	char _type[2] { 'B', 'M' };
	uint32_t _file_size = 0;
	uint32_t _reserved = 0;
	uint32_t _pixel_offset = sizeof(BmpFileHeader) + sizeof(BmpImageHeader);

	BmpFileHeader() = default;
	BmpFileHeader(uint32_t file_size) {
		_file_size = file_size;
	}

	explicit operator std::string() const {
		return std::format(
			"TYPE: {}{}, _file_size={}, _pixel_offset={}",
			_type[0], _type[1], _file_size, _pixel_offset
		);
	}
};

#pragma pack(pop)

// Used to indicate failure to process a BMP file.
class BmpFileError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// Represents a pixel in YCbCr color space. Note that we order the color
// components as JPEGs in the wild do, YCrCb (red first).
struct YCbCrPixel {
	uint8_t y;
	uint8_t cr;
	uint8_t cb;

	uint8_t operator[](unsigned int idx) const {
		if (idx == 0) return y;
		if (idx == 1) return cr;
		if (idx == 2) return cb;

		return 0;
	}
};

// Given some dimension or other value v, calculate how much we need to add to
// v to reach the first integer multiple of d greater than or equal to v.
size_t calc_pad(size_t v, size_t d) {
	size_t pad = d - v%d;
	if (pad == d) pad = 0;
	return pad;
}

// Calculates the number of padding bytes to use per row in a BMP file, given
// the width of the rows in pixels.
uint8_t calc_bmp_padding_bytes(size_t width) {
	return static_cast<uint8_t>(calc_pad(width*3, 4));
}

// Converts a double value to a uint8_t, clamping the value to the limits of
// an 8 bit unsigned integer.
uint8_t clamped_convert(double x) {
	if (x < 0.0) return 0;
	if (x > 255.0) return 255;

	return static_cast<uint8_t>(x);
}

// Helper class for performing subsampling and color conversion on
// a read BMP file. Not required for writing.
class BmpSubsamplingWindow {
private:
	std::vector<YCbCrPixel> rows;
	std::fstream* file;
	uint32_t width;
	uint32_t height;
	uint32_t subsamp;
	bool reverse_order;
	bool first_read = true;
	uint32_t consumed_rows = 0;

	// Our own pixel padding in `rows` to simplify padding operations. Must
	// be an integer multiple of `subsamp`.
	uint32_t width_padded;

	// BMP padding bytes. Required to properly read pixel data.
	uint32_t file_padding_bytes;

	// What image row does a sampling coordinate of y=0 currently
	// correspond to?
	//
	// For these examples, we use a `subsamp` value of 3. If
	// `reverse_order` is false, then `zero_y_map` starts at zero and
	// increments by 3 as the window advances through the file - this is
	// the same order as the image (top-down).
	//
	//                   y   windows
	//                      ______...
	// zero_y_map=0 ---> 0 |      ...
	//                   1 |      ...   BMP FILE
	//                   2 |______...      |
	// zero_y_map=3 ---> 0 |      ...      |
	//                   1 |      ...      V
	//                   2 |______...
	//                   ...
	//
	// However if `reverse_order` is true, then `zero_y_map` starts at the last
	// multiple of 3 contained in the image and decrements by 3, since the
	// BMP is storing the rows bottom-up. The sampling coordinate y is also 
	// inverted.
	//
	//                     y   windows
	//                        ______...
	//                     2 |      ...
	//                     1 |      ...   BMP FILE
	// zero_y_map=L   ---> 0 |______...      |
	//                     2 |      ...      |
	//                     1 |      ...      V
	// zero_y_map=L-3 ---> 0 |______...
	//                     ...
	//
	// IMPORTANT: All private functions not explicitly mentioning
	// zero_y_map use *internal* y coordinates: y=0 always refers to the
	// first row of the `rows` vector, regardless of `reverse_order`. 
	uint32_t zero_y_map;

	// Calculates the index into the `rows` vector given (x, y)
	// coordinates.
	unsigned int coord(unsigned int x, unsigned int y) const {
		if (y >= subsamp) throw BmpFileError("coord: bad y");
		if (x >= width_padded) throw BmpFileError("coord: bad x");

		return y*width_padded + x;
	}

	// Copies row y=`from` to y=`to`.
	void copy_row(unsigned int from, unsigned int to) {
		unsigned int x0_from = coord(0, from);
		unsigned int x0_to = coord(0, to);

		for (unsigned int x = 0; x < width_padded; x++) {
			rows[x0_to + x] = rows[x0_from + x];
		}
	}

	// If we have horizontal subsample padding, copy the last pixel
	// to the rest of the columns for better averages
	void x_copy(unsigned int y) {
		for (unsigned int i = width; i < width_padded; i++) {
			rows[coord(i, y)] = rows[coord(width - 1, y)];
		}
	}

	// Copies row y=`first_filled` to all rows of lower y coordinate. If
	// first_filled == 0, this does nothing.
	void y_copy_lower(unsigned int first_filled) {
		for (unsigned int y = 0; y < first_filled; y++) {
			copy_row(first_filled, y);
		}
	}

	// Copies row `last_filled` to all rows of higher y coordinate. If
	// last_filled == subsamp - 1 this does nothing.
	void y_copy_upper(unsigned int last_filled) {
		for (unsigned int y = last_filled + 1; y < subsamp; y++) {
			copy_row(last_filled, y);
		}
	}

	// Consumes padding bytes, as specified by BMP standard. If
	// `file_padding_bytes` == 0, this does nothing.
	void consume_padding_bytes() {
		for (unsigned int i = 0; i < file_padding_bytes; i++) {
			file->get();
		}
	}

	// Consumes a row from the file, including all BMP padding.
	void consume_row(unsigned int y) {
		unsigned int x0 = coord(0, y);
		for (unsigned int x = 0; x < width; x++) {
			double r = file->get();
			double g = file->get();
			double b = file->get();

			double luma = 0.299*r + 0.587*g + 0.114*b;
			double cb = (-0.299*r - 0.587*g + 0.886*b) / 1.772 + 128;
			double cr = ( 0.701*r - 0.587*g - 0.114*b) / 1.402 + 128;

			rows[x0 + x] = YCbCrPixel(
				clamped_convert(luma),
				clamped_convert(cr),
				clamped_convert(cb)
			);
		}

		consumed_rows++;
		consume_padding_bytes();
	}

	// Initializes zero_y_map. See comment on that variable.
	void init_zero_y_map() {
		if (reverse_order) {
			int m = height % subsamp;
			if (m == 0) {
				zero_y_map = height - subsamp;
			} else {
				zero_y_map = height - m;
			}
		} else {
			zero_y_map = 0;
		}
	}

	// Advances zero_y_map. See comment on that variable.
	void inc_zero_y_map() {
		if (reverse_order) {
			zero_y_map -= subsamp;
		} else {
			zero_y_map += subsamp;
		}
	}
public:
	BmpSubsamplingWindow() = delete;
	BmpSubsamplingWindow(
		std::fstream& file,
		uint32_t width,
		uint32_t height,
		uint32_t subsamp,
		bool reverse_order
	) {
		this->file = &file;
		this->width = width;
		this->height = height;
		this->subsamp = subsamp;
		this->reverse_order = reverse_order;
		this->file_padding_bytes = calc_bmp_padding_bytes(width);
		this->width_padded = width + static_cast<uint32_t>(calc_pad(width, subsamp));

		// Init rows. This vector will never be resized, only written
		// over.
		rows.reserve(subsamp*width_padded);
		for (unsigned int i = 0; i < subsamp*width_padded; i++) {
			rows.push_back(YCbCrPixel(0, 0, 0));
		}
	}

	uint32_t get_zero_y_map() const { return zero_y_map; }

	bool no_more_rows() const {
		return file->eof() || consumed_rows >= height;
	}

	// Advances the window by `subsamp` rows. Returns true if we were able
	// to advance, false otherwise.
	bool advance() {
		if (no_more_rows()) return false;

		unsigned int y_start = 0;
		if (first_read && reverse_order) {
			msg::warn(
				"READ_BMP: BMP window using reverse order, "
				"skipping height pad on first advance."
			);
			y_start = static_cast<unsigned int>(calc_pad(height, subsamp));
		}

		unsigned int last_filled;
		for (unsigned int y = y_start; y < subsamp; y++) {
			consume_row(y);
			x_copy(y);
			last_filled = y;
			if (no_more_rows()) break;
		}

		y_copy_upper(last_filled);

		if (first_read) {
			y_copy_lower(y_start);
			init_zero_y_map();
			first_read = false;
		} else {
			inc_zero_y_map();
		}

		return true;
	}

	// Samples the window. If comp=0 (luma) no subsampling is performed;
	// the exact luma value at (x, y) of this window will be returned.
	//
	// If comp=1,2 (chroma), then samp_y is ignored and a subsamp^2 square
	// will be averaged to produce the subsampled value, starting at `samp_x`.
	// For example, with subsamp=3, samp_x=3:
	//
	//  x 012345678...
	// y     V
	// 0  ...sss......
	// 1  ...sss......
	// 2  ...sss......
	// 
	uint8_t sample(unsigned int samp_x, unsigned int samp_y, unsigned int comp) const {
		if (comp == 0) {
			unsigned int internal_y = samp_y;
			if (reverse_order) internal_y = subsamp - samp_y - 1;
			return rows[coord(samp_x, internal_y)][comp];
		}

		// Here, specific Y value doesn't matter, since we'll be
		// covering all of them in the subsampling operation.
		uint16_t result = 0;
		unsigned int x_start = samp_x * subsamp;
		unsigned int x_end = x_start + subsamp;
		for (unsigned int y = 0; y < subsamp; y++) {
			for (unsigned int x = x_start; x < x_end; x++) {
				result += rows[coord(x, y)][comp];
			}
		}

		return (uint8_t)(result / (subsamp*subsamp));
	}
};

// Implementation of a BMP file read/writer, only allowing 24-bit RGB color.
// Automatically converts the BMP's RGB to YCbCr.
class BmpFile {
private:
	std::fstream file;
	std::ios_base::openmode mode;

	// Writes a generic type byte-for-byte out to the file.
	template<typename T>
	void write_generic(const T& x) {
		// File is little-endian, so can just write stuff out directly
		file.write(reinterpret_cast<const char*>(&x), sizeof(x));
	}

	// Writes a BMP file header. The parameters other than file_size are
	// determined by the default values specified in BmpFileHeader.
	void write_file_header(uint32_t file_size) {
		write_generic(BmpFileHeader(file_size));
	}

	// Writes a BMP image header. Same deal as write_file_header.
	void write_image_header(uint32_t width, uint32_t height) {
		write_generic(BmpImageHeader(width, height));
	}

	// Reads a generic type byte-for-byte from the file. Must be manually
	// instantiated, for example `auto x = read_generic<MyType>();`.
	template<typename T>
	T read_generic() {
		T hdr;
		file.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
		return hdr;
	}
public:
	BmpFile(std::string path, std::ios_base::openmode mode) {
		msg::debug("FILE_BMP: try open BmpFile \"{}\"...", path);

		if (mode != std::ios_base::in && mode != std::ios_base::out) {
			throw BmpFileError(
				"BmpFile(): invalid mode, pick one of: ios_base::in, ios_base::out"
			);
		}

		this->mode = mode;

		file.exceptions(std::ios_base::badbit | std::ios_base::failbit | std::ios_base::eofbit);
		file.open(path, mode | std::ios_base::binary);
	}

	bool is_read_mode() {
		return mode == std::ios_base::in;
	}

	bool is_write_mode() {
		return mode == std::ios_base::out;
	}

	// Reads this BMP file with some degree of chroma subsampling.
	ImageBuffer read(uint8_t subsamp=2) {
		if (!is_read_mode()) throw BmpFileError("BmpFile: cannot read() in write mode");

		msg::debug("READ_BMP: reading BMP file and converting to YCbCr...");

		ImageBuffer buf { SampleFormat::Sample8 };

		auto file_header = read_generic<BmpFileHeader>();
		auto t = file_header._type;
		if (t[0] != 'B' || t[1] != 'M') throw BmpFileError("not a BMP file");

		auto image_header = read_generic<BmpImageHeader>();

		if (image_header._header_size < sizeof(BmpImageHeader)) {
			throw BmpFileError("BMP file too old");
		}

		if (image_header._bit_count != 24) {
			throw BmpFileError("BMP: only 24-bit color is supported");
		}

		if (image_header._compression != 0) {
			throw BmpFileError("BMP: compression not supported");
		}

		msg::debug("READ_BMP: checks pass, init components...");

		// Luma component.
		buf.new_component(Dimensions(
			image_header._width,
			std::abs(image_header._height)
		));

		// Cr, Cb
		for (int i = 0; i < 2; i++) {
			buf.new_component(Dimensions(
				image_header._width / subsamp,
				std::abs(image_header._height) / subsamp
			));
		}

		for (int j = 0; j < 3; j++) {
			msg::debug("READ_BMP: comp {}: {}", j, std::string(buf[j]));
		}

		// Subsampling plus the rows potentially being in reverse order
		// is a massive headache when reading, so use a helper class.
		auto window = BmpSubsamplingWindow(
			file,
			image_header._width,
			std::abs(image_header._height),
			subsamp,
			image_header._height >= 0
		);

		msg::debug("READ_BMP: start conversion...");
		while (window.advance()) {
			unsigned int ymap = window.get_zero_y_map();
			unsigned int ymap_color = ymap / subsamp;

			for (unsigned int y = 0; y < subsamp; y++) {
				for (unsigned int luma_x = 0; luma_x < image_header._width; luma_x++) {
					uint8_t s = window.sample(luma_x, y, 0);
					buf[0][luma_x, y + ymap] = s;
				}
			}

			for (unsigned int comp = 1; comp < 3; comp++) {
				for (unsigned int color_x = 0; color_x < image_header._width / subsamp; color_x++) {
					uint8_t s = window.sample(color_x, 0, comp);
					buf[comp][color_x, ymap_color] = s;
				}
			}
		}

		msg::debug("READ_BMP: done.");
		return std::move(buf);
	}

	// Writes the result of a JPEG decode. The components are assumed to
	// be YCbCr, possibly with some degree of chroma subsampling.
	// This will be converted to RGB during the writing process.
	void write(ImageBuffer buf) {
		if (!is_write_mode()) throw BmpFileError("BmpFile: cannot write() in read mode");

		msg::debug("WRITE_BMP: converting to RGB and writing BMP file...");

		const auto& img_size = buf[0].valid_size();
		uint8_t padding_bytes = calc_bmp_padding_bytes(img_size.width());

		uint32_t file_size = static_cast<uint32_t>(
			sizeof(BmpFileHeader)
			+ sizeof(BmpImageHeader)
			+ (img_size.width() * 3 + padding_bytes) * img_size.height()
		);

		msg::debug(
			"WRITE_BMP: BMP fsize={} padding={}",
			file_size, padding_bytes
		);

		write_file_header(file_size);
		// height is given as a negative number so the rows are
		// ordered top-to-bottom.
		write_image_header(
			static_cast<uint32_t>(img_size.width()),
			-static_cast<int32_t>(img_size.height())
		);

		auto cb_subsamp = buf.subsamp(2);
		auto cr_subsamp = buf.subsamp(1);

		for (uint32_t y = 0; y < img_size.height(); y++) {
			for (uint32_t x = 0; x < img_size.width(); x++) {
				double luma = buf[0][x, y];

				// NOTE: This works, but it is the opposite of
				// what the JFIF standard says... component[1]
				// should be Cb, not Cr.
				double cb = buf[2][
					x / static_cast<uint32_t>(cb_subsamp.x),
					y / static_cast<uint32_t>(cb_subsamp.y)
				];
				double cr = buf[1][
					x / static_cast<uint32_t>(cr_subsamp.x),
					y / static_cast<uint32_t>(cr_subsamp.y)
				];

				double r = luma + 1.402 * (cr - 128.0);
				double g = luma - 0.34414 * (cb - 128.0) - 0.71414 * (cr - 128.0);
				double b = luma + 1.772 * (cb - 128.0);

				file.put(clamped_convert(r));
				file.put(clamped_convert(g));
				file.put(clamped_convert(b));
			}

			for (int i = 0; i < padding_bytes; i++) {
				file.put(0);
			}
		}
	}
};
}
