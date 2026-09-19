//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file classic_stream.hpp
 * @brief Streams pinned to the classic ("C") locale.
 *
 * `fast_number.hpp` makes `strtod`/`snprintf` independent of `LC_NUMERIC`,
 * which is what `setlocale` (a Qt application's `setlocale(LC_ALL, "")`, or
 * Python's `locale.setlocale(locale.LC_ALL, "")`) moves. C++ streams are driven
 * by a different knob: a stream is imbued at construction with the *global*
 * `std::locale`, and only `std::locale::global` changes that. A host that calls
 * it with a named locale then makes every stream this library constructs
 * format and parse numbers by that locale's rules, in both directions:
 *
 *  - **Reading** -- `iss >> x` on a `double` uses `std::num_get`, so `1.5` reads
 *    as `1` under a comma-decimal locale.
 *  - **Writing** -- `os << n` on an *integer* uses `std::num_put`, which applies
 *    the locale's digit grouping: `1234567` is written `1.234.567`. That is not a
 *    misread file but a corrupt one, and it needs no decimal point to happen.
 *
 * Every stream constructed under `src/cpp/` therefore goes through one of the
 * factories below, which imbue `std::locale::classic()`. In the default
 * configuration -- the global locale *is* classic, as in every test and CI run --
 * that changes no formatting decision, so output is byte-identical. The
 * `tests/python/test_no_locale_sensitive_number_io.py` guard fails a bare
 * `std::ifstream x(...)` so a new call site cannot reintroduce the exposure.
 *
 * The file-stream factories imbue *before* `open()`. `classic()`'s
 * `codecvt<char, char, mbstate_t>` is the identity, so imbuing after would be
 * harmless too, but an unopened stream sidesteps the question entirely.
 *
 * Deliberately not used by `detail/format_compat.hpp`'s `std::format` fallback:
 * that is an installed header whose function is a template, so changing it is an
 * ABI break, and it only renders log lines and exception text -- never a byte
 * of a file and never anything that is re-parsed.
 */

// System includes
#include <filesystem>
#include <fstream>
#include <ios>
#include <locale>
#include <sstream>
#include <string>
#include <utility>

namespace meshioplusplus {
namespace detail {

/**
 * Pins @p rStream to the classic locale, whatever `std::locale::global` says.
 *
 * Skips the `imbue` when the stream already carries it -- the common case, since
 * a stream takes the global locale at construction and that is classic unless a
 * host changed it. `imbue` recomputes the stream's cached facets, and several
 * readers build one stream per input line, so the comparison (a pointer test in
 * every standard library) is what keeps the default path free.
 */
inline void imbue_classic(std::ios_base& rStream) {
    if (rStream.getloc() != std::locale::classic())
        rStream.imbue(std::locale::classic());
}

/// A read stream over @p Text (taken by value: the stream copies it anyway).
inline std::istringstream make_classic_istringstream(
    std::string Text, std::ios_base::openmode Mode = std::ios_base::in) {
    std::istringstream stream(std::move(Text), Mode);
    imbue_classic(stream);
    return stream;
}

/// An empty output string stream.
inline std::ostringstream make_classic_ostringstream(
    std::ios_base::openmode Mode = std::ios_base::out) {
    std::ostringstream stream(Mode);
    imbue_classic(stream);
    return stream;
}

/// A read/write string stream initialised with @p Text.
inline std::stringstream make_classic_stringstream(
    std::string Text = std::string(),
    std::ios_base::openmode Mode = std::ios_base::in | std::ios_base::out) {
    std::stringstream stream(std::move(Text), Mode);
    imbue_classic(stream);
    return stream;
}

/// An input file stream over @p rPath; on failure `fail()` is set, as with the constructor.
inline std::ifstream make_classic_ifstream(const std::filesystem::path& rPath,
                                           std::ios_base::openmode Mode = std::ios_base::in) {
    std::ifstream stream;
    imbue_classic(stream);
    stream.open(rPath, Mode);
    return stream;
}

/// An output file stream over @p rPath; on failure `fail()` is set, as with the constructor.
inline std::ofstream make_classic_ofstream(const std::filesystem::path& rPath,
                                           std::ios_base::openmode Mode = std::ios_base::out) {
    std::ofstream stream;
    imbue_classic(stream);
    stream.open(rPath, Mode);
    return stream;
}

}  // namespace detail
}  // namespace meshioplusplus
