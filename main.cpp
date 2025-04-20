import arg;
import std;
import jpeg;
import msg;

static const std::string empty_path = "";
static const std::string default_bmp = "out.bmp";
static const std::string default_jpeg = "out.jpeg";

void usage() {
	std::cout << (
		"usage: jpeg [-h|--help] MODE [INPUT_FILE] [-o OUTPUT_FILE] [options...]\n"
		"options:\n"
		"	-h|--help      Ignores every other option and prints this message.\n"
		"	-o OUT         Specify output file. Accepted by the following modes:\n"
		"	               encode, decode, filetest-bmp, filetest-jpeg. If not \n"
		"	               specified, a default name will be chosen for you.\n"
		"	-q QUALITY     A number from 0 to 100 specifying the quality of the image.\n"
		"	               Accepted only by encode mode. Defaults to 80.\n"
		"	-s SUBSAMP     Divide dimensions of color components (Cb, Cr) by this number.\n"
		"	               Must be either 1 or 2. Defaults to 2.\n"
		"\nmain modes:\n"
		"	scan           Reads a JPEG file and dumps all of its metadata to stdout.\n"
		"	encode         Compress an input BMP file.\n"
		"	decode         Decompress an input JPEG file.\n"
		"\ntest modes:\n"
		"	filetest-bmp   Read a BMP file and immediately write it back.\n"
		"	filetest-jpeg  Read a JPEG file and immediately write it back.\n"
		"	qualitytest    Read a BMP file and encode it using multiple quality settings,\n"
		"	               from 0 to 100 in steps of 10.\n"
		"	markertest     Prints every possible marker string to stdout.\n"
	) << std::endl;
}

const std::string& in_file_arg(const arg::Args& args, const std::string& mode) {
	// 2nd positional arg is always input file
	if (args.positional.size() < 2) {
		msg::error("{} mode requires input file argument", mode);
		return empty_path;
	}

	return args.positional[1];
}

std::tuple<const std::string&, const std::string&> io_file_args(
	const arg::Args& args, const std::string& mode, const std::string& default_out
) {
	auto& in_file = in_file_arg(args, mode);

	if (args.has_opt("-o")) {
		return {in_file, args.options.at("-o")};
	} else {
		return {in_file, default_out};
	}
}

class Mode {
public:
	const std::string name;
	Mode(std::string&& name) : name{std::move(name)} {}
	virtual int run(const arg::Args& args) = 0;
};

class EncodeMode : public Mode {
public:
	EncodeMode() : Mode("encode") {}

	int run(const arg::Args& args) {
		auto [infile, outfile] = io_file_args(args, name, default_jpeg);
		if (infile == empty_path) return 1;

		unsigned int quality = args.getopt_uint("-q", 80, 0, 100);
		unsigned int subsamp = args.getopt_uint("-s", 2, 1, 2);

		jpeg::BmpFile bmp_file { infile, std::ios_base::in };
		jpeg::JpegEncoder encoder { quality };

		auto raws = bmp_file.read(subsamp);
		auto encoded = encoder.encode(raws);

		jpeg::JpegFile jpeg { outfile, std::ios_base::out };
		jpeg.write(encoded);

		return 0;
	}
};

class DecodeMode : public Mode {
public:
	DecodeMode() : Mode("decode") {}

	int run(const arg::Args& args) {
		auto [infile, outfile] = io_file_args(args, name, default_bmp);
		if (infile == empty_path) return 1;

		jpeg::JpegFile jpeg_file { infile, std::ios_base::in };
		jpeg::JpegDecoder decoder;

		auto jpeg_data = jpeg_file.read();
		auto decoded = decoder.decode(jpeg_data);

		jpeg::BmpFile bmp { outfile, std::ios_base::out };
		bmp.write(decoded);

		return 0;
	}
};

class ScanMode : public Mode {
public:
	ScanMode() : Mode("scan") {}

	int run(const arg::Args& args) {
		auto& filepath = in_file_arg(args, name);
		if (filepath == empty_path) return 1;

		jpeg::JpegFile jpeg_file { filepath, std::ios_base::in };
		jpeg::CompressedJpegData jpeg_data = jpeg_file.read();

		std::cout << "Dumping read file..." << std::endl;
		std::cout << std::string(jpeg_data) << std::endl;

		return 0;
	}
};

class QualityTestMode : public Mode {
public:
	QualityTestMode() : Mode("qualitytest") {}

	int run(const arg::Args& args) {
		auto infile = in_file_arg(args, name);
		if (infile == empty_path) return 1;

		jpeg::BmpFile in_bmp { infile, std::ios_base::in };
		auto raws = in_bmp.read();

		for (unsigned int q = 0; q <= 100; q += 10) {
			auto fname = std::format("qualitytest{:03d}.jpg", q);
			jpeg::JpegEncoder encoder { q } ;
			jpeg::JpegFile jpeg { fname, std::ios_base::out };

			jpeg.write(encoder.encode(raws));
		}

		return 0;
	}
};

// Both BmpFile and JpegFile have the same general interface, so can use a
// template for testing them.
template<typename FileType>
class FileTestMode : public Mode {
public:
	const std::string default_out_name;

	FileTestMode(
		std::string&& mode_name,
		std::string&& default_out_name
	) : Mode(std::move(mode_name)), default_out_name{default_out_name} {}

	int run(const arg::Args& args) {
		auto [infile, outfile] = io_file_args(args, name, default_out_name);
		if (infile == empty_path) return 1;

		FileType in_fobj { infile, std::ios_base::in };
		auto intermediate = in_fobj.read();

		FileType out_fobj { outfile, std::ios_base::out };
		out_fobj.write(intermediate);

		return 0;
	}
};

class MarkerTestMode : public Mode {
private:
	static void print_one_marker(unsigned char n) {
		jpeg::Marker marker { n };
		std::cout << std::string(marker) << std::endl;
	}
public:
	MarkerTestMode() : Mode("markertest") {}

	// Prints every possible JPEG file marker to make sure they're parsed
	// correctly.
	int run(const arg::Args& args) {
		for (unsigned int i=0x00; i <= 0xFF; i++) print_one_marker(i);
		return 0;
	}
};

class MatrixTestMode : public Mode {
private:
	// Convenience function for printing a matrix.
	template<typename T, size_t W, size_t H>
	static void printmat(const jpeg::Matrix<T, W, H>& mat) {
		std::cout << std::string(mat) << std::endl;
	}
public:
	MatrixTestMode() : Mode("matrixtest") {}

	// Tests matrix multiplication and the DCT transform.
	int run(const arg::Args& args) {
		jpeg::Matrix<std::uint16_t, 8, 8> imat;
		try {
			printmat(imat);
		} catch (const std::exception& ex) {
			msg::error("exc: {}", ex.what());
		}

		for (int i = 0; i < 64; i++) {
			auto& valref = imat.zz(i);
			valref = i;
		}

		printmat(imat);
		printmat(jpeg::dct_mat);
		printmat(jpeg::dct_mat_t);

		auto transformed = jpeg::dct_mat_t * imat * jpeg::dct_mat;
		printmat(transformed);

		auto detransformed = jpeg::dct_mat * transformed * jpeg::dct_mat_t;
		printmat(detransformed);

		return 0;
	}
};

// Convenience function to cut down on repetition a bit...
template<typename ModeType, typename... Args>
void push_mode(std::vector<std::unique_ptr<Mode>>& v, Args&&... args) {
	v.push_back(std::make_unique<ModeType>(std::forward<Args>(args)...));
}

// Creates a vector containing all mode objects (as unique_ptrs to base)
auto make_mode_vec() {
	std::vector<std::unique_ptr<Mode>> v;

	// Main modes
	push_mode<EncodeMode>(v);
	push_mode<DecodeMode>(v);
	push_mode<ScanMode>(v);

	// Tests
	push_mode<QualityTestMode>(v);
	push_mode<FileTestMode<jpeg::JpegFile>>(v, "filetest-jpeg", "ftest.jpeg");
	push_mode<FileTestMode<jpeg::BmpFile>>(v, "filetest-bmp", "ftest.bmp");
	push_mode<MarkerTestMode>(v);
	push_mode<MatrixTestMode>(v);

	return v;
}

int main_inner(int argc, char* argv[]) {
	arg::ArgParser parser {{
		{"--help", arg::noconsume},
		{"-h",     arg::noconsume},
		{"-o",     arg::consume},
		{"-q",     arg::consume},
		{"-s",     arg::consume}
	}};

	auto args = parser.parse(argc, argv);

	if (args.has_opt("--help") || args.has_opt("-h")) {
		usage();
		return 0;
	}

	if (args.positional.size() < 1) {
		msg::error("Requires at least one argument!");
		usage();
		return 1;
	}

	const std::string& modestr = args.positional[0];

	auto modes = make_mode_vec();
	for (const auto& mode : modes) {
		if (mode->name == modestr) return mode->run(args);
	}
	
	msg::error("unrecognized mode {}", modestr);
	usage();

	return 1;
}

int main(int argc, char* argv[]) {
	try {
		return main_inner(argc, argv);
	} catch (const std::ios_base::failure& e) {
		// FIXME Unfortunately on Windows this message is completely
		// useless... may need to replace the ios_base::failure throws
		// with something else.
		std::cerr << "FATAL (IO): " << e.code().message() << std::endl;
		return 1;
	} catch (const std::exception& e) {
		std::cerr << "FATAL: " << e.what() << std::endl;
		return 1;
	}
}
