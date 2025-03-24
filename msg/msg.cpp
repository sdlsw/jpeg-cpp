export module msg;

import std;

export namespace msg {
enum class LogLevel {
	fine  = 0,
	debug = 1,
	info  = 2,
	warn  = 3,
	error = 4
};

std::map<LogLevel, std::string> str_from_loglevel {
	{LogLevel::fine,  "FINE"},
	{LogLevel::debug, "DEBUG"},
	{LogLevel::info,  "INFO"},
	{LogLevel::warn,  "WARN"},
	{LogLevel::error, "ERROR"}
};

LogLevel level();
void level(LogLevel l);

bool is_enabled(LogLevel l) {
	return int(l) >= int(level());
}

template<typename... Args>
void log_to(std::ostream& os, LogLevel l, std::format_string<Args...> format, Args&&... args) {
	if (!is_enabled(l)) return;

	os << "[" << str_from_loglevel[l] << "] " 
		<< std::format(format, std::forward<Args>(args)...) 
		<< std::endl;
}

template<typename... Args>
inline void log(LogLevel l, std::format_string<Args...> format, Args&&... args) {
	return log_to(std::cerr, l, format, std::forward<Args>(args)...);
}

template<typename... Args>
inline void fine(std::format_string<Args...> format, Args&&... args) {
	return log(LogLevel::fine, format, std::forward<Args>(args)...);
}

template<typename... Args>
inline void debug(std::format_string<Args...> format, Args&&... args) {
	return log(LogLevel::debug, format, std::forward<Args>(args)...);
}

template<typename... Args>
inline void info(std::format_string<Args...> format, Args&&... args) {
	return log(LogLevel::info, format, std::forward<Args>(args)...);
}

template<typename... Args>
inline void warn(std::format_string<Args...> format, Args&&... args) {
	return log(LogLevel::warn, format, std::forward<Args>(args)...);
}

template<typename... Args>
inline void error(std::format_string<Args...> format, Args&&... args) {
	return log(LogLevel::error, format, std::forward<Args>(args)...);
}
}
