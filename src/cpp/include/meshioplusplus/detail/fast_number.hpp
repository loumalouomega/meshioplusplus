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
 * @file fast_number.hpp
 * @brief Locale-independent floating-point parsing and formatting.
 *
 * `std::strtod` and `std::snprintf`'s `%f`/`%e`/`%g` conversions both honour
 * the process's `LC_NUMERIC` category, which no first-party code pins: a
 * host that adopts the environment's locale (a Qt application calling
 * `setlocale(LC_ALL, "")` at startup, or Python code running
 * `locale.setlocale(locale.LC_ALL, "")`) makes every ASCII reader and writer
 * in this codebase silently misread/miswrite `1.5` as `1` (or emit `1,5`)
 * under a comma-decimal locale. `parse_double`/`snprintf_c` give every
 * caller a `strtod`/`snprintf`-compatible replacement that always uses `.`
 * as the decimal separator, regardless of the process locale.
 *
 * `std::from_chars`'s floating-point overload is deliberately NOT used as
 * the sole implementation: it is a real Emscripten/libc++ hazard (see the
 * comment in `formats/gid_read.cpp`, which predates this header and is the
 * decision this header formalizes rather than reverses) — some libc++
 * versions declare the floating-point overload without defining it, a
 * link-time rather than compile-time failure. Where it IS trustworthy
 * (`__cpp_lib_to_chars` and not Emscripten) it is used as the fast path,
 * falling back to a locale-pinned `strtod_l`, and finally to plain `strtod`
 * with a decimal-point repair for the rare platform with neither.
 *
 * `snprintf_c` deliberately does not use `std::to_chars` for writing: this
 * codebase's format strings carry width, flags and case that `to_chars`
 * does not reproduce (`"%25.16E"`, `"%.17g"`, ...), and reproducing them
 * would mean re-implementing padding — the risk the performance section of
 * the project roadmap explicitly flags ("the reference files are the gate,
 * not the claim"). Locale independence is achieved by formatting normally
 * and repairing the decimal point afterward, which is a no-op (one cheap
 * `localeconv()` call and a comparison, then nothing) in the common C-locale
 * case, so this is provably byte-identical to plain `snprintf` there.
 *
 * Integer parsing (`strtol`/`strtoll`) is deliberately out of scope: the
 * C locale's `LC_NUMERIC` affects only the radix character and the
 * (grouping-only, opt-in) thousands separator, neither of which appears in
 * the subject sequence `strtol` accepts — an integer literal has no decimal
 * point to misparse. Only floating-point parsing and formatting are
 * affected, so only those get a locale-independent replacement here.
 *
 * The `istringstream operator>>` sites elsewhere in the codebase are a
 * separate, lower-priority exposure: a default-constructed stream is
 * imbued with `std::locale::classic()` unless a host calls
 * `std::locale::global`, which is a narrower and rarer trigger than
 * `setlocale` (Qt/Python's own habit). Not addressed by this header.
 */

// System includes
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <version>

#if !defined(MESHIOPLUSPLUS_FORCE_NO_FAST_FROM_CHARS) && !defined(__EMSCRIPTEN__) && \
    defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#define MESHIOPLUSPLUS_HAS_FAST_FROM_CHARS 1
#include <charconv>
#include <system_error>
#endif

#if defined(_MSC_VER)
#include <locale.h>
// MSVC's locale.h declares `_locale_t`/`_create_locale`, not POSIX's
// `locale_t`/`newlocale` -- alias the type so the rest of this header can
// use one spelling for both.
using locale_t = _locale_t;
#define MESHIOPLUSPLUS_HAS_LOCALE_T 1
#elif defined(__EMSCRIPTEN__) || defined(__linux__) || defined(__APPLE__) || \
    defined(__FreeBSD__) || defined(__unix__)
#include <locale.h>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <xlocale.h>  // newlocale/strtod_l live here on BSD-derived libc
#endif
#define MESHIOPLUSPLUS_HAS_LOCALE_T 1
#endif

// Printf-style format-string checking, disabled on toolchains that don't
// support the GNU `format` attribute (MSVC).
#if defined(__GNUC__) || defined(__clang__)
#define MESHIOPLUSPLUS_PRINTF_LIKE(fmt_idx, first_arg_idx) \
    __attribute__((format(printf, fmt_idx, first_arg_idx)))
#else
#define MESHIOPLUSPLUS_PRINTF_LIKE(fmt_idx, first_arg_idx)
#endif

namespace meshioplusplus {
namespace detail {

#ifdef MESHIOPLUSPLUS_HAS_LOCALE_T

/**
 * @brief A `locale_t` handle fixed to the "C" numeric convention, created
 * once and never freed (a magic static; thread-safe to initialize per
 * [stmt.dcl]/4, and cheaper than creating one per call).
 * @return The C-numeric `locale_t`, or a null handle if creation failed
 * (caller must treat that as tier 2 being unavailable).
 */
inline locale_t c_numeric_locale() noexcept {
#if defined(_MSC_VER)
    static locale_t loc = _create_locale(LC_NUMERIC, "C");
#else
    static locale_t loc = newlocale(LC_NUMERIC_MASK, "C", static_cast<locale_t>(nullptr));
#endif
    return loc;
}

#endif  // MESHIOPLUSPLUS_HAS_LOCALE_T

/**
 * @brief Whether `pDecimalPoint` (from `localeconv()`) is the single
 * character `'.'`, i.e. whether the active locale's decimal point is
 * already the one this header guarantees -- the case in which every
 * function below is a pure pass-through with no repair work to do.
 */
inline bool is_c_decimal_point(const char* pDecimalPoint) noexcept {
    return pDecimalPoint != nullptr && pDecimalPoint[0] == '.' && pDecimalPoint[1] == '\0';
}

/**
 * @brief Locale-independent equivalent of `std::strtod`, cursor form.
 *
 * Same contract as `std::strtod`: parses a floating-point number starting
 * at `pFirst`, advances `rEnd` to just past the parsed text, and returns
 * `0.0` with `rEnd == pFirst` if no conversion could be performed. Always
 * parses with `.` as the decimal separator, regardless of the process
 * locale.
 *
 * @param pFirst NUL-terminated text to parse from (only the prefix up to
 * the parsed number is read; the buffer must outlive `rEnd`).
 * @param rEnd Set to the first unparsed character.
 * @return The parsed value, or `0.0` on failure (`rEnd == pFirst`).
 */
inline double parse_double(const char* pFirst, const char*& rEnd) noexcept {
#ifdef MESHIOPLUSPLUS_HAS_FAST_FROM_CHARS
    // from_chars, unlike strtod, accepts neither leading whitespace nor a
    // leading '+' -- skip both here, and fall through to the locale-pinned
    // tier on anything else from_chars declines (an otherwise malformed
    // prefix, ...), so all three tiers agree bit-for-bit.
    const char* p = pFirst;
    while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\v' || *p == '\f' || *p == '\r')
        ++p;
    if (*p == '+')
        ++p;
    // A hex float ("0x1p3") is the one case from_chars does NOT simply
    // decline: chars_format::general parses its leading "0" as a complete,
    // valid decimal literal and reports SUCCESS with ptr stopping at the
    // 'x' -- silently under-parsing relative to strtod, which greedily
    // consumes the whole hex float. Detect the "0x"/"0X" prefix (optionally
    // signed) upfront and skip straight to the locale-pinned tier for it,
    // rather than trusting a from_chars "success" that isn't one.
    const char* hex_check = (*p == '-') ? p + 1 : p;
    const bool looks_like_hex_float =
        hex_check[0] == '0' && (hex_check[1] == 'x' || hex_check[1] == 'X');
    if (!looks_like_hex_float) {
        double value = 0.0;
        const char* last = pFirst + std::strlen(pFirst);
        auto [ptr, ec] = std::from_chars(p, last, value);
        if (ec == std::errc()) {
            rEnd = ptr;
            return value;
        }
    }
#endif
#ifdef MESHIOPLUSPLUS_HAS_LOCALE_T
    locale_t loc = c_numeric_locale();
    if (loc != static_cast<locale_t>(nullptr)) {
        char* end = nullptr;
#if defined(_MSC_VER)
        double value2 = _strtod_l(pFirst, &end, loc);
#else
        double value2 = strtod_l(pFirst, &end, loc);
#endif
        rEnd = end;
        return value2;
    }
#endif
    // Tier 3: plain strtod, repairing the decimal point if the active
    // locale's isn't already '.'. Best-effort -- reached only when neither
    // from_chars nor a locale_t is available on this platform.
    const char* dp = localeconv()->decimal_point;
    if (is_c_decimal_point(dp)) {
        char* end = nullptr;
        double value3 = std::strtod(pFirst, &end);
        rEnd = end;
        return value3;
    }
    // A non-'.' single-character decimal point: copy the token into a
    // stack buffer with the separator swapped for '.', bounded to a length
    // no real numeric literal exceeds.
    char buf[64];
    std::size_t n = 0;
    const std::size_t dp_len = std::strlen(dp);
    const char* src = pFirst;
    while (*src != '\0' && n + dp_len < sizeof(buf)) {
        if (dp_len == 1 && *src == dp[0]) {
            buf[n++] = '.';
            ++src;
        } else {
            buf[n++] = *src++;
        }
    }
    buf[n] = '\0';
    char* end = nullptr;
    double value4 = std::strtod(buf, &end);
    rEnd = pFirst + (end - buf);
    return value4;
}

/**
 * @brief Locale-independent equivalent of `std::strtod(rText.c_str(), nullptr)`.
 * @param rText The text to parse.
 * @return The parsed value, or `0.0` if `rText` has no valid numeric prefix.
 */
inline double parse_double(const std::string& rText) noexcept {
    const char* end = nullptr;
    return parse_double(rText.c_str(), end);
}

/**
 * @brief Locale-independent, throwing equivalent of `std::stod(rText)`.
 *
 * For CLI option parsing and similar call sites that want `std::stod`'s
 * throw-on-failure contract rather than `parse_double`'s lenient `0.0`, so
 * they can be migrated with a plain name substitution instead of adding a
 * cursor check at every call site.
 *
 * Unlike `std::stod`, this does not distinguish "no valid prefix" from
 * "valid but out of `double` range" -- both raise the same
 * `std::invalid_argument`, never `std::out_of_range`. An out-of-range
 * literal (`"1e400"`) instead saturates to +/-`HUGE_VAL`, matching
 * `strtod`'s own C-standard behaviour. Every current call site treats any
 * exception identically (a CLI argument error), so this difference is not
 * observable in practice; it is called out here because it is the one
 * place this header's contract genuinely diverges from the standard
 * function it replaces.
 *
 * @param rText The text to parse.
 * @return The parsed value.
 * @throws std::invalid_argument if `rText` has no valid numeric prefix.
 */
inline double stod_c(const std::string& rText) {
    const char* end = nullptr;
    const double v = parse_double(rText.c_str(), end);
    if (end == rText.c_str())
        throw std::invalid_argument("stod_c: no conversion for '" + rText + "'");
    return v;
}

/**
 * @brief Locale-independent `std::vsnprintf`.
 *
 * Same contract as `std::snprintf`, with the decimal separator always `.`
 * regardless of the process locale. In the (overwhelmingly common) case
 * where the process locale's own decimal point is already `.`, this calls
 * `std::vsnprintf` and returns with no further work -- byte-identical to
 * calling `std::snprintf` directly.
 *
 * @param pBuf Destination buffer.
 * @param cap Size of `pBuf`.
 * @param pFmt A `printf`-style format string.
 * @return Same as `std::snprintf`: the number of characters that would
 * have been written (excluding the NUL), or a negative value on an
 * encoding error.
 */
MESHIOPLUSPLUS_PRINTF_LIKE(3, 4)
inline int snprintf_c(char* pBuf, std::size_t cap, const char* pFmt, ...) {
    va_list args;
    va_start(args, pFmt);
    const int n = std::vsnprintf(pBuf, cap, pFmt, args);
    va_end(args);

    const char* dp = localeconv()->decimal_point;
    if (n <= 0 || is_c_decimal_point(dp))
        return n;

    // Repair the locale's decimal point within the bytes actually written
    // (n may exceed cap on truncation -- only the written prefix exists to
    // repair). A formatted number carries at most one decimal point, so a
    // single-character locale separator (every real-world locale; a
    // multi-character one is a POSIX theoretical case with no practical
    // decimal_point that isn't length 1) is replaced in place with no
    // length change. A multi-character separator is left as a documented
    // best-effort gap rather than risking a buffer-length bug for a case
    // that does not occur on any real target this project ships for.
    const std::size_t dp_len = std::strlen(dp);
    if (dp_len != 1)
        return n;
    const std::size_t written =
        (cap == 0)
            ? 0
            : (static_cast<std::size_t>(n) < cap - 1 ? static_cast<std::size_t>(n) : cap - 1);
    for (std::size_t i = 0; i < written; ++i) {
        if (pBuf[i] == dp[0]) {
            pBuf[i] = '.';
            break;  // at most one decimal point per formatted number
        }
    }
    return n;
}

}  // namespace detail
}  // namespace meshioplusplus
