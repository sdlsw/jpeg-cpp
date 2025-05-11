export module jpeg:data.jpeg;

// data_jpeg.cpp:
// Data structures and constants common among all other files.

import std;
import msg;
import :data.matrix;
import :data.tables;

using std::uint8_t;
using std::int8_t;
using std::uint16_t;
using std::int16_t;

// Converts a byte vector to a string, inserting a space every 4 bytes:
//
// 01020304 05060708 090A0B0C...
const std::string hexchars = "0123456789ABCDEF";
std::string string_from_byte_vec(const std::vector<uint8_t>& v) {
	std::string s;
	s.reserve(
		2 * v.size() + // data chars
		v.size() / 4   // spaces
	);

	int spacercount = -1;
	for (int i = 0; i < v.size(); i++) {
		uint8_t x = v[i];
		uint8_t high = (0xF0 & x) >> 4;
		uint8_t low = 0xF & x;

		if (spacercount == 3) {
			s.push_back(' ');
			spacercount = 0;
		} else {
			spacercount++;
		}

		s.push_back(hexchars[high]);
		s.push_back(hexchars[low]);
	}

	return std::move(s);
}

export namespace jpeg {
// The standard defines four modes of operation. For now, we only support
// the baseline procedure, but doesn't hurt to recognize the others.
enum class OperationMode : uint8_t {
	baseline_seq = 0,
	extended_seq = 1,
	progressive =  2,
	lossless =     3 // always sequential
};

const std::map<OperationMode, std::string> str_from_operation_mode {
	{OperationMode::baseline_seq, "BASELINE_SEQUENTIAL"},
	{OperationMode::extended_seq, "EXTENDED_SEQUENTIAL"},
	{OperationMode::progressive,  "PROGRESSIVE"},
	{OperationMode::lossless,     "LOSSLESS"}
};

enum class CodingType : uint8_t {
	huffman =    0,
	arithmetic = 1
};

const std::map<CodingType, std::string> str_from_coding_type {
	{CodingType::huffman,    "HUFFMAN"},
	{CodingType::arithmetic, "ARITHMETIC"}
};

// Markers in JPEG files are two bytes, 0xFFnn. The first byte is always 0xFF.
const uint8_t mark_start = 0xFF;

// Bitfield marker bases and masks.
const uint8_t sof_base   = 0xC0;
const uint8_t sof_mask   = 0xF0;
const uint8_t rst_base   = 0xD0;
const uint8_t rst_mask   = 0xF8;
const uint8_t app_base   = 0xE0;
const uint8_t app_mask   = 0xF0;

// Range marker bases and masks.
const uint8_t jpgn_start = 0xF0;
const uint8_t jpgn_end   = 0xFD;
const uint8_t jpgn_mask  = 0x0F;
const uint8_t res_start  = 0x02;
const uint8_t res_end    = 0xBF;

// Invalid markers are not actually markers, so they aren't defined under
// MarkerSpecial. Just give them their own constants.
const uint8_t mark_invalid_high = 0xFF;
const uint8_t mark_invalid_low  = 0x00;

struct FrameType {
	OperationMode mode = OperationMode::baseline_seq;
	CodingType coding = CodingType::huffman;
	bool is_differential = false;

	FrameType() = default;
	FrameType(OperationMode sof_t, CodingType ct, bool diff) :
		mode{sof_t}, coding{ct}, is_differential{diff} {}

	// Constructs from marker char. In this case, the enum definitions
	// we're using line up with the values given in the standard,
	// so just use bit magic to get them.
	FrameType(uint8_t m) :
		mode{OperationMode(m & 0x03)},
		coding{CodingType((m & 0x08) >> 3)},
		is_differential{bool(m & 0x04)} {}

	// Checks if this FrameType describes the baseline process.
	// Only the mode is checked, since the baseline mode can *only* be
	// non-differential and Huffman coded.
	bool is_baseline() const {
		return mode == OperationMode::baseline_seq;
	}

	// Reconstructs the 8-bit marker value from the separate fields.
	uint8_t marker_value() const {
		return (
			sof_base |
			((uint8_t) mode) |
			(((uint8_t) coding) << 3) |
			(((uint8_t) is_differential) << 2)
		);
	}

	explicit operator std::string() const {
		return std::format(
			"SOF {} | {}DIFFERENTIAL | {}-CODED",
			str_from_operation_mode.at(mode),
			(is_differential) ? "" : "NON-",
			str_from_coding_type.at(coding)
		);
	}
};

// Special marker code assignments. "Special" markers do not
// contain any additional information in the marker code, so just hardcode
// their values.
enum class MarkerSpecial : uint8_t {
	// TEM is used only in arithmetic coding as far as I know
	tem_temporary_use =                   0x01,

	// Weird markers - These are jammed into the range occupied mostly by
	// the various SOF markers.
	dht_define_huff_table =               0xC4,
	jpg_reserved =                        0xC8,
	dac_define_arith_coding =             0xCC,

	// Other markers
	soi_start_of_img =                    0xD8,
	eoi_end_of_img =                      0xD9,
	sos_start_of_scan =                   0xDA,
	dqt_define_quant_table =              0xDB,
	dnl_define_num_lines =                0xDC,
	dri_define_rst_interval =             0xDD,
	dhp_define_hierarchical_progression = 0xDE,
	exp_expand_reference_components =     0xDF,
	com_comment =                         0xFE
};

const std::map<MarkerSpecial, std::string> symbol_from_special_marker {
	{MarkerSpecial::tem_temporary_use,                   "TEM*"},
	{MarkerSpecial::dht_define_huff_table,               "DHT"},
	{MarkerSpecial::jpg_reserved,                        "JPG"},
	{MarkerSpecial::dac_define_arith_coding,             "DAC"},
	{MarkerSpecial::soi_start_of_img,                    "SOI*"},
	{MarkerSpecial::eoi_end_of_img,                      "EOI*"},
	{MarkerSpecial::sos_start_of_scan,                   "SOS"},
	{MarkerSpecial::dqt_define_quant_table,              "DQT"},
	{MarkerSpecial::dnl_define_num_lines,                "DNL"},
	{MarkerSpecial::dri_define_rst_interval,             "DRI"},
	{MarkerSpecial::dhp_define_hierarchical_progression, "DHP"},
	{MarkerSpecial::exp_expand_reference_components,     "EXP"},
	{MarkerSpecial::com_comment,                         "COM"}
};

// Class for encapsulating all the magic numbers/bitmasks involved in
// decoding the marker values.
struct Marker {
	// Note: Technically markers are two bytes, but the first
	// byte is always 0xFF and doesn't contain any interesting
	// information, so leave it out of our representiation.
	uint8_t marker;

	bool is_sof() const {
		// SOF Markers are 0b1100xxxx, except for some special
		// markers stuffed in the middle of the SOF range.
		return (
			(marker & sof_mask) == sof_base &&
			marker != (uint8_t)(MarkerSpecial::dht_define_huff_table) &&
			marker != (uint8_t)(MarkerSpecial::jpg_reserved) &&
			marker != (uint8_t)(MarkerSpecial::dac_define_arith_coding)
		);
	}

	// Extracts frame type information from an SOF marker.
	FrameType parse_sof() const {
		// SOF markers are basically bitfields, so
		// delegate parsing to another struct...
		return FrameType(marker);
	}

	// Constructs an SOF marker.
	static Marker sof(const FrameType& type) {
		return Marker(type.marker_value());
	}

	bool is_rstn() const {
		return (marker & rst_mask) == rst_base;
	}

	// Gets the 'n' of RSTn.
	uint8_t parse_rstn() const {
		return marker & (~rst_mask);
	}

	// Constructs an RSTn marker.
	static Marker rstn(uint8_t n) {
		n &= ~rst_mask;
		return Marker(rst_base | n);
	}

	bool is_appn() const {
		return (marker & app_mask) == app_base;
	}

	// Gets the 'n' of APPn.
	uint8_t parse_appn() const {
		return marker & (~app_mask);
	}

	// Constructs an APPn marker.
	static Marker appn(uint8_t n) {
		n &= ~app_mask;
		return Marker(app_base | n);
	}

	// Checks if this marker is JPGn. These markers are specifically
	// reserved for JPEG extensions and are very uncommon. They will
	// be ignored by this implementation.
	bool is_jpgn() const {
		return marker >= jpgn_start && marker <= jpgn_end;
	}

	// Gets the 'n' of JPGn.
	uint8_t parse_jpgn() const {
		return marker & jpgn_mask;
	}

	// Checks if this marker is RES, a reserved value. These markers should
	// be ignored.
	bool is_res() const {
		return marker >= res_start && marker <= res_end;
	}

	// Checks if this marker is "special", meaning it does not fit
	// into any of the other families and does not include any additional
	// information.
	bool is_special() const {
		return !(
			is_sof() ||
			is_rstn() ||
			is_appn() ||
			is_jpgn() ||
			is_res()
		);
	}

	// Converts this marker to one of the "special" constants.
	MarkerSpecial parse_special() const {
		return MarkerSpecial { marker };
	}

	bool is_soi() const {
		return MarkerSpecial(marker) == MarkerSpecial::soi_start_of_img;
	}

	bool is_eoi() const {
		return MarkerSpecial(marker) == MarkerSpecial::eoi_end_of_img;
	}

	// Tests if a marker is a valid value. The standard specifically marks
	// two values as "invalid". They are used to encode 0xFF in the entropy
	// coded segments.
	bool is_valid() const {
		return (
			marker != mark_invalid_low &&
			marker != mark_invalid_high
		);
	}

	// Checks if a marker is "standalone". If true, then this marker does
	// not precede a marker segment.
	bool is_standalone() const {
		if (is_rstn()) return true;

		switch (marker) {
			case (uint8_t)(MarkerSpecial::soi_start_of_img):
			case (uint8_t)(MarkerSpecial::eoi_end_of_img):
			case (uint8_t)(MarkerSpecial::tem_temporary_use):
				return true;
			default:
				return false;
		}
	}

	operator int() const { return marker; };

	explicit operator std::string() const {
		if (!is_valid()) {
			return std::format("INVALID_MARKER 0x{:X}", marker);
		}

		if (is_sof()) {
			return std::string(parse_sof());
		} else if (is_rstn()) {
			return std::format("RST{}*", int(parse_rstn()));
		} else if (is_appn()) {
			return std::format("APP{}", parse_appn());
		} else if (is_jpgn()) {
			return std::format("JPG{}", parse_jpgn());
		} else if (is_res()) {
			return std::format("RESERVED 0x{:X}", marker);
		} else {
			return symbol_from_special_marker.at(parse_special());
		}
	}
};

// APPn marker segment.
struct AppSegment {
	unsigned int n; // APPn

	// For a standard JPEG impl, we have no idea what application segments
	// contain, so leave it as unparsed data.
	std::vector<uint8_t> data;

	explicit operator std::string() const {
		return std::format("APP{}: {}", n, string_from_byte_vec(data));
	}
};

// Information about one component in a scan.
struct ScanComponentParams {
	// Sets which component in the frame this scan component refers to.
	uint8_t component_selector = 0;

	// Sets which DC and AC tables should be selected while decoding this
	// data.
	uint8_t dc_entropy_coding_selector = 0;
	uint8_t ac_entropy_coding_selector = 0;

	explicit operator std::string() const {
		return std::format(
			"ID={} | DCSEL={} | ACSEL={}",
			component_selector,
			dc_entropy_coding_selector,
			ac_entropy_coding_selector
		);
	}
};

// A JPEG scan. A scan is a contiguous region of encoded data that when
// decoded, will set coefficients in the entire area of one or more components.
// That is, if a component consists of N blocks, then this scan will specify N
// blocks. Scans may contain the blocks of more than one component. In this
// case, the data will be interleaved, with the scan containing the same number
// of blocks as the sum of the number of blocks contained in each component.
// 
// The exact ordering of blocks in the scan is not trivial. See BlockView.
struct Scan {
	// The index of this scan in the scan array. Set automatically
	// by Frame::new_scan().
	unsigned int index;

	// The number of components in this scan. If >1, this scan contains
	// interleaved data.
	uint8_t num_components = 0;

	// The selection parameters for each component in this scan.
	// These parameters are stored in the same order the components will
	// occur in the coded scan data.
	std::vector<ScanComponentParams> component_params;

	// Values for operation modes other than baseline. Unused.
	uint8_t start_spectral_sel = 0;
	uint8_t end_spectral_sel = 63;
	uint8_t succ_approx_high = 0;
	uint8_t succ_approx_low = 0;

	// Number of MCUs contained in one restart interval.
	uint16_t restart_interval = 0;

	// Entropy-coded data as it is encoded in the file, with RST markers
	// removed in favor of breaking up the intervals into sub-vectors.
	std::vector<std::vector<uint8_t>> entropy_coded_data;

	// Calculates total number of entropy-coded bytes in this scan.
	size_t total_bytes_data() const {
		size_t sum = 0;

		for (const auto& segment : entropy_coded_data) {
			sum += segment.size();
		}

		return sum;
	}

	// Shorthand for getting a const reference to the nth component
	// parameter set.
	const ScanComponentParams& comp_params(uint8_t n) const {
		return component_params[n];
	}

	explicit operator std::string() const {
		std::stringstream out;

		out << std::format(
			"Ss={} | Se={} | Ah={} | Al={} | RSTi={}\n",
			start_spectral_sel,
			end_spectral_sel,
			succ_approx_high,
			succ_approx_low,
			restart_interval
		);

		for (const auto& comp : component_params) {
			out << std::format(
				"COMPONENT {}\n",
				std::string(comp)
			);
		}

		out << std::format("{} SEGMENT(S)\n", entropy_coded_data.size());
		out << std::format("{} BYTES IMAGE DATA", total_bytes_data());

		return out.str();
	}
};

// Represents a raw AC coefficient code as read from scan data.
struct AcCoeff {
	// Zero run length. If >0, then this coefficient is preceded by R
	// zeros. Max of 15.
	uint8_t r = 0;

	// Magnitude category.
	uint8_t s = 0;

	// The value of the coefficient. May be zero, in which case this
	// value codes r+1 zeros.
	int16_t value = 0;

	// Checks if this value codes EOB (end of block).
	// If this value is encountered while decoding, then the rest of the
	// block will be filled with zeros.
	bool is_eob() const { return r == 0 && s == 0; }
};

// Information about one component in a frame.
struct FrameComponentParams {
	// The ID of the component. Used in scans to identify which components
	// each scan contains.
	uint8_t identifier = 0;

	// The component index for these params.
	unsigned int index;

	// Horizontal size of a block chunk to be contained in one MCU.
	// Also determines size of this component relative to others.
	uint8_t horizontal_sampling_factor = 0;

	// Vertical size of a block chunk to be contained in one MCU.
	// Also determines size of this component relative to others.
	uint8_t vertical_sampling_factor = 0;

	// Sets the quantization table to use for this component.
	uint8_t qtable_selector = 0;

	// Shortcut functions for access
	auto h() const { return horizontal_sampling_factor; }
	auto v() const { return vertical_sampling_factor; }

	// The area (number of blocks) of an MCU region for this component
	// in interleaved scans.
	unsigned int interleave_area() const {
		return horizontal_sampling_factor * vertical_sampling_factor;
	}

	explicit operator std::string() const {
		return std::format(
			"ID={} | H={} | V={} | QTABLE={}",
			identifier,
			horizontal_sampling_factor,
			vertical_sampling_factor,
			qtable_selector
		);
	}
};

// A JPEG frame. A frame defines one image, consisting of one or more
// components (typically YCrCb). The frame structure dictates the size of each
// component, and a number of other settings, constructing the "destination"
// that compressed image data will be decoded to. Along with a frame are one
// or more "scans", which are decoded and written into the decompressed
// destination data.
struct Frame {
	// Global settings for how the frame shall be encoded.
	FrameType type;

	// The index of this frame in the frame array. Set automatically
	// by CompressedJpegData::new_frame().
	unsigned int index;

	// Number of bits on one image sample. Should always be 8 for baseline
	// operation.
	uint8_t sample_precision = 8;

	// Height (Y) of the image in pixels.
	uint16_t num_lines = 0;

	// Width (X) of the image in pixels.
	uint16_t samples_per_line = 0;

	// Number of components in this frame.
	uint8_t num_components = 0;

	// Destinations for the Huffman tables. Entropy coding tables are
	// selected on a per-scan basis, so store them under the frame.
	// FIXME may want to revisit for hierarchical, but i'm not sure anybody
	// even uses it
	TableMap<HuffmanTable> _ac_huff_tables;
	TableMap<HuffmanTable> _dc_huff_tables;

	// Parameters for each component in the scan.
	std::vector<FrameComponentParams> component_params;
	std::vector<Scan> scans;

	// Creates a new scan in this frame.
	Scan& new_scan() {
		scans.emplace(scans.end());
		_ac_huff_tables.save_dest();
		_dc_huff_tables.save_dest();

		Scan& s = scans.back();
		s.index = scans.size() - 1;

		return s;
	}

	// Adds a new component to this frame.
	FrameComponentParams& new_component() {
		component_params.emplace(component_params.end());
		FrameComponentParams& p = component_params.back();
		p.index = component_params.size() - 1;
		return p;
	}

	// Gets a set of FrameComponentParams by searching for its component ID.
	const FrameComponentParams& get_component_params_by_id(uint8_t id) const {
		for (const auto& params : component_params) {
			if (params.identifier == id) return params;
		}

		throw std::runtime_error(std::format(
			"get_component_params_by_id: could not find params by id {}",
			id
		));
	}

	// Gets a set of FrameComponentParams referred to by a component in a
	// scan.
	const FrameComponentParams& get_component_params_by_scan(
		const Scan& scan,
		uint8_t comp
	) const {
		auto comp_id = scan.comp_params(comp).component_selector;
		return get_component_params_by_id(comp_id);
	}

	template <typename Self>
	auto& ac_huff_tables(this Self& self) { return self._ac_huff_tables; }

	template <typename Self>
	auto& dc_huff_tables(this Self& self) { return self._dc_huff_tables; }

	template <typename Self>
	auto& huff_tables(this Self& self, TableClass cls) {
		if (cls == TableClass::ac) {
			return self.ac_huff_tables();
		} else {
			return self.dc_huff_tables();
		}
	}

	// Obtains a reference to the proper huffman table to use for the nth
	// component in the given scan.
	const HuffmanTable& get_huff_table(
		const Scan& scan,
		unsigned int comp_index,
		TableClass cls
	) const {
		uint8_t selector = 0;

		if (cls == TableClass::dc) {
			selector = scan.comp_params(comp_index).dc_entropy_coding_selector;
		} else {
			selector = scan.comp_params(comp_index).ac_entropy_coding_selector;
		}

		return huff_tables(cls).get(scan.index, selector);
	}

	template<typename BinaryOp>
	unsigned int accumulate_component_params(BinaryOp op) const {
		return std::accumulate(
			component_params.begin(),
			component_params.end(),
			0,
			op
		);
	}

	// Finds the maximum of a particular field among this frame's component
	// params.
	template<typename AccessOp>
	unsigned int max_component_param(AccessOp op) const {
		auto acc_op = [&op](unsigned int acc, const FrameComponentParams& p) {
			return std::max((unsigned int) (p.*op)(), acc);
		};
		return accumulate_component_params(acc_op);
	}

	// Determines the maximum horizontal sampling factor across all
	// components.
	unsigned int h_max() const {
		return max_component_param(&FrameComponentParams::h);
	}

	// Determines the maximum vertical sampling factor accross all
	// components.
	unsigned int v_max() const {
		return max_component_param(&FrameComponentParams::v);
	}

	// Determines the length in blocks of one MCU containing data from
	// every component.
	int mcu_length() const {
		auto op = [](unsigned int acc, const FrameComponentParams& p) {
			return acc + p.interleave_area();
		};
		return accumulate_component_params(op);
	}

	explicit operator std::string() const {
		std::stringstream out;
		out << std::format(
			"MARK: {}\nP={} | X={} | Y={}\n",
			std::string(type),
			sample_precision,
			samples_per_line,
			num_lines
		);

		for (const auto& comp : component_params) {
			out << std::format(
				"COMPONENT {}\n",
				std::string(comp)
			);
		}

		out << "DC HUFFMAN " << std::string(_dc_huff_tables);
		out << "AC HUFFMAN " << std::string(_ac_huff_tables);

		int count = 0;
		for (const auto& scan : scans) {
			out << std::format(
				"SCAN{}:\n{}\n",
				count,
				std::string(scan)
			);
			count++;
		}

		return out.str();
	}
};

// Encapsulates all data included in a JPEG file, including tables,
// entropy-coded data, comments and application segments, etc.
class CompressedJpegData {
private:
	bool valid = false;
	std::string _source_file {""};

	TableMap<QuantizationTable> _q_tables;

	// Application segments contained in this JPEG.
	std::vector<AppSegment> _app_segments;
	std::vector<std::vector<uint8_t>> _comments;

	// The image frames contained in this JPEG. Typically there will only
	// be one, as each frame models one image. The baseline coding
	// process only supports one.
	std::vector<Frame> _frames;
public:
	bool is_valid() {
		return valid;
	}

	void set_valid() {
		valid = true;
	}

	auto& app_segments() { return _app_segments; }
	auto& comments() { return _comments; }

	template <typename Self>
	auto& q_tables(this Self& self) { return self._q_tables; }

	auto& frames() { return _frames; }

	const auto& frames() const { return _frames; }

	// Creates a new frame in this JPEG.
	Frame& new_frame() {
		_frames.emplace(_frames.end());
		_q_tables.save_dest();

		Frame& f = _frames.back();
		f.index = _frames.size() - 1;

		return f;
	}

	// Gets a reference to the current frame.
	Frame& cur_frame() {
		if (_frames.size() == 0) throw std::runtime_error(
			"No frames, cannot get current."
		);

		return _frames.back();
	}

	// Convenience function. Creates new scan in current frame.
	Scan& new_scan() {
		return cur_frame().new_scan();
	}

	// Obtains a reference to the proper quantization table to use for the
	// nth component in the given scan.
	const QuantizationTable& get_qtable(
		const Frame& frame,
		const Scan& scan,
		unsigned int comp_index
	) const {
		auto params = frame.get_component_params_by_scan(scan, comp_index);
		return _q_tables.get(frame.index, params.qtable_selector);
	}

	// Obtains a vector of pointers to each QuantizationTable specified
	// in each component of the given frame (by index)
	auto get_qtables(const Frame& frame) const {
		std::vector<const QuantizationTable*> v;

		for (const auto& params : frame.component_params) {
			const auto& tbl = _q_tables.get(frame.index, params.qtable_selector);
			v.push_back(&tbl);
		}

		return v;
	}

	explicit operator std::string() const {
		std::stringstream out;

		for (const auto& seg : _app_segments) {
			out << std::string(seg) << '\n';
		}

		for (const auto& com : _comments) {
			out << "COM " << string_from_byte_vec(com) << '\n';
		}

		out << "QUANTIZATION " << std::string(_q_tables);

		if (_frames.size() == 0) {
			out << "NO FRAMES\n";
		} else {
			int count = 0;
			for (const auto& f : _frames) {
				out << std::format(
					"FRAME{}:\n{}\n",
					count,
					std::string(f)
				);
				count++;
			}
		}

		return out.str();
	}
};
}
