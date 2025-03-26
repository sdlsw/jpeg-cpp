import std;
import jpeg;
import msg;

void print_one_marker(unsigned char n) {
	jpeg::Marker marker { n };
	std::cout << std::string(marker) << std::endl;
}

// Prints every possible JPEG file marker to make sure they're parsed
// correctly.
void print_all_markers() {
	for (unsigned int i=0x00; i <= 0xFF; i++) print_one_marker(i);
}

template<typename T, size_t W, size_t H>
void printmat(const jpeg::Matrix<T, W, H>& mat) {
	std::cout << std::string(mat) << std::endl;
}

void othertest() {
	std::int16_t v = 0b0000'0000'1010'0010;
	std::int16_t t = 9;

	std::int16_t extended = jpeg::ScanDecodeView::extend(v, t);

	std::cout << std::format("{:d} {:b}", extended, (std::uint16_t) extended) << std::endl;
}

void matrixtest() {
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
}

int main_inner(int argc, char* argv[]) {
	if (argc < 2) {
		std::cout << "Needs at least one argument." << std::endl;
		return 1;
	}

	std::string mode {argv[1]};
	
	if (mode == "markertest") {
		print_all_markers();
	} else if (mode == "matrixtest") {
		matrixtest();
	} else if (mode == "othertest") {
		othertest();
	} else if (mode == "scan") {
		if (argc < 3) {
			std::cout << "scan mode requires file arg" << std::endl;
			return 1;
		}

		std::string filepath {argv[2]};
		jpeg::JpegFile jpeg_file {filepath};

		jpeg::CompressedJpegData jpeg_data = jpeg_file.read();

		std::cout << "Dumping read file..." << std::endl;
		std::cout << std::string(jpeg_data) << std::endl;
	} else if (mode == "decode") {
		if (argc < 3) {
			std::cout << "decode mode requires file arg" << std::endl;
			return 1;
		}

		std::string filepath {argv[2]};
		jpeg::JpegFile jpeg_file {filepath};
		jpeg::CompressedJpegData jpeg_data = jpeg_file.read();

		jpeg::JpegDecoder decoder;
		std::vector<jpeg::RawComponent> decoded = decoder.decode(jpeg_data);

		jpeg::BmpFile bmp { "out.bmp" };
		bmp.write(decoded);
	} else {
		std::cout << "unrecognized mode" << std::endl;
		return 1;
	}

	return 0;
}

int main(int argc, char* argv[]) {
	try {
		main_inner(argc, argv);
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
