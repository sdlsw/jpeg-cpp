module msg;

import std;

static msg::LogLevel _loglevel = msg::LogLevel::debug;

namespace msg {
	LogLevel level() {
		return _loglevel;
	}

	void level(LogLevel l) {
		_loglevel = l;
	}
}

