#pragma once
//
// meshioplusplus::log — minimal, header-only logging built on C++20 std::format and
// std::source_location.
//
//   meshioplusplus::log::warn("MED: orientation for '{}' not implemented", type);
//
// * Format strings are compile-time checked (std::format_string).
// * Every message carries its call site (file:line) automatically.
// * Runtime filtering via the MESHIOPLUSPLUS_LOG_LEVEL environment variable:
//   "debug" | "info" | "warn" (default) | "error" | "off". Read once.
// * Messages go to stderr through std::osyncstream, so concurrent messages
//   (e.g. from parallel_for bodies) never interleave mid-line.
// * A filtered-out call costs one branch — no formatting, no allocation.

#include <cstdlib>
#include <format>
#include <iostream>
#include <source_location>
#include <string_view>
#include <syncstream>
#include <type_traits>
#include <utility>

namespace meshioplusplus {
namespace log {

enum class Level : int { Debug = 0, Info = 1, Warn = 2, Error = 3, Off = 4 };

inline Level threshold() {
    static const Level lvl = [] {
        const char* env = std::getenv("MESHIOPLUSPLUS_LOG_LEVEL");
        if (env == nullptr) return Level::Warn;
        std::string_view s(env);
        if (s == "debug") return Level::Debug;
        if (s == "info") return Level::Info;
        if (s == "warn" || s == "warning") return Level::Warn;
        if (s == "error") return Level::Error;
        if (s == "off" || s == "none") return Level::Off;
        return Level::Warn;
    }();
    return lvl;
}

inline bool enabled(Level lvl) {
    return static_cast<int>(lvl) >= static_cast<int>(threshold());
}

inline void write(Level lvl, std::string_view msg, const std::source_location& loc) {
    constexpr std::string_view names[] = {"debug", "info", "warning", "error"};
    std::string_view file = loc.file_name();
    if (auto p = file.find_last_of("/\\"); p != std::string_view::npos)
        file.remove_prefix(p + 1);
    std::osyncstream(std::cerr) << std::format("meshio {} [{}:{}] {}\n",
                                               names[static_cast<int>(lvl)], file,
                                               loc.line(), msg);
}

// Wraps a compile-time-checked format string together with the caller's
// source location (captured implicitly via the defaulted parameter).
template <class... Args>
struct format_with_location {
    std::format_string<Args...> fmt;
    std::source_location loc;

    template <class S>
    consteval format_with_location(  // NOLINT(google-explicit-constructor)
        const S& s, std::source_location l = std::source_location::current())
        : fmt(s), loc(l) {}
};

template <class... Args>
void debug(format_with_location<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Debug)) return;
    write(Level::Debug, std::format(f.fmt, std::forward<Args>(args)...), f.loc);
}

template <class... Args>
void info(format_with_location<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Info)) return;
    write(Level::Info, std::format(f.fmt, std::forward<Args>(args)...), f.loc);
}

template <class... Args>
void warn(format_with_location<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Warn)) return;
    write(Level::Warn, std::format(f.fmt, std::forward<Args>(args)...), f.loc);
}

template <class... Args>
void error(format_with_location<std::type_identity_t<Args>...> f, Args&&... args) {
    if (!enabled(Level::Error)) return;
    write(Level::Error, std::format(f.fmt, std::forward<Args>(args)...), f.loc);
}

}  // namespace log
}  // namespace meshioplusplus
