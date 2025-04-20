export module arg;

// arg.cpp:
// Basic command line argument parsing.

import std;

export namespace arg {
// Settings for whether an option is "consuming" or "non-consuming".
//
// A consuming option has an argument after it (ex. `--option 30`),
// while non-consuming options are standalone (ex. `--help`).
const bool consume = true;
const bool noconsume = false;

// Parent class of all exceptions emitted by this module.
class ArgError : public std::runtime_error {
public:
	using std::runtime_error::runtime_error;
};

// Class representing a set of command line arguments.
struct Args {
	// The string this program was invoked with (arg 0).
	std::string progname;

	// Positional arguments.
	std::vector<std::string> positional;

	// Optional arguments. Must start with either - or --. To access the
	// option's value, the map key must include the preceding dashes.
	//
	// For example, `args.options["--opt-name"]`.
	//
	// If the option was nonconsuming (no argument), then its value will be
	// "".
	std::map<std::string, std::string> options;

	// Test if a given option was specified.
	bool has_opt(const std::string& opt) const {
		return options.contains(opt);
	}

	// Gets a numeric option as an unsigned integer. If the option wasn't
	// specified, this returns `default_value` instead. In addition, if the
	// chosen value is not between `range_low` and `range_high`, an
	// exception will be thrown.
	unsigned int getopt_uint(
		const std::string& opt,
		unsigned int default_value = 0,
		unsigned int range_low = 0,
		unsigned int range_high = -1
	) const {
		if (!has_opt(opt)) return default_value;

		auto value = (unsigned int) stoi(options.at(opt));

		if (value < range_low) {
			throw ArgError(std::format(
				"{}: cannot be less than {}",
				opt, range_low
			));
		} else if (value > range_high) {
			throw ArgError(std::format(
				"{}: cannot be greater than {}",
				opt, range_high
			));
		} else {
			return value;
		}
	}
};

// Class implementing commandline argument parsing behavior.
class ArgParser {
private:
	// argv cursor. The next string consumed will be argv[n].
	unsigned int n = 0;

	// argc value currently being parsed.
	unsigned int cur_argc = 0;

	// argv array currently being parsed.
	char** cur_argv;

	// if true, the current parse has failed.
	bool fail = false;

	// Option definition map for this parser.
	std::map<std::string, bool> opt_def;

	// Validates an option definition map. Throws an exception if there's
	// something wrong, does nothing otherwise.
	static void validate_opt_def(const std::map<std::string, bool>& opt_def) {
		for (const auto& [key, _] : opt_def) {
			if (key.find('-') != 0) {
				std::string msg = std::format(
					"ArgParser: bad opt_def key '{}', "
					"keys must start with - or --",
					key
				);

				// Not ArgError for this one since this is
				// usually programmer error.
				throw std::invalid_argument(msg);
			}
		}
	}

	// Check whether an option is recognized. Throws an exception if this
	// is not the case.
	void ensure_opt_defined(const std::string& opt) const {
		if (!opt_def.contains(opt)) {
			throw ArgError(std::format(
				"unknown option '{}'", opt
			));
		}
	}

	// Consumes the next string in argv, returning a copy in a std::string
	// object.
	std::string consume_str() {
		if (n >= cur_argc) {
			fail = true;
			return "";
		}

		auto s = std::string(cur_argv[n]);
		n++;

		return std::move(s);
	}

	// Consumes an option from argv and places it in args.options.
	void consume_opt(Args& args, const std::string& opt) {
		ensure_opt_defined(opt);
		if (opt_def[opt]) {
			// consuming option
			args.options[opt] = consume_str();
			if (fail) throw ArgError(std::format(
				"opt '{}' expected argument, but none provided", opt
			));
		} else {
			// non-consuming option
			args.options[opt] = "";
		}
	}
public:
	ArgParser(const std::map<std::string, bool>& opt_def) {
		validate_opt_def(opt_def);
		this->opt_def = opt_def;
	}

	ArgParser(std::map<std::string, bool>&& opt_def) {
		validate_opt_def(opt_def);
		this->opt_def = std::move(opt_def);
	}

	// Parse an array of arguments as passed to a main() function.
	// This parser may be used multiple times.
	Args parse(int argc, char* argv[]) {
		n = 0;
		cur_argc = argc;
		cur_argv = argv;
		fail = false;

		Args args;
		std::string s;

		args.progname = consume_str();
		if (fail) throw ArgError("could not parse progname");

		while (n < cur_argc) {
			s = consume_str();

			if (s.find('-') == 0) {
				consume_opt(args, s);
			} else {
				args.positional.push_back(std::move(s));
			}
		}

		return std::move(args);
	}
};
}
