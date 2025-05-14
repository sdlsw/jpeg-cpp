module msg;

import std;

#ifdef DEBUG
static msg::LogLevel _loglevel = msg::LogLevel::debug;
#else
static msg::LogLevel _loglevel = msg::LogLevel::info;
#endif

namespace msg {
	LogLevel level() {
		return _loglevel;
	}

	void level(LogLevel l) {
		_loglevel = l;
	}
}

