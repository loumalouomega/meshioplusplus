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
 * @file detail/row_writer.hpp
 * @brief ASCII writers' rows formatted in parallel, in fixed-size chunks, and
 * streamed in order: byte-identical to formatting them one by one.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * `write_row_chunks` hands each worker a chunk of rows and one reused buffer
 * -- not a heap string per row, the pattern `abaqus`, `ansysinp` and `lsdyna`
 * used -- and writes the chunks in row order, a bounded number at a time.
 * `CNumber` formats a number exactly as `snprintf_c` does, but resolves the
 * locale's decimal point once, outside the loop: `snprintf_c` calls
 * `localeconv()` per value, and POSIX does not require it to be thread-safe.
 * Integers go through `std::to_chars`, which is what `ostream << int` prints
 * in the classic locale. Roadmap §3, "Parallel row formatting in ASCII
 * writers".
 */

// System includes
#include <algorithm>
#include <charconv>
#include <clocale>
#include <cmath>
#include <cstring>
#include <system_error>
#include <version>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/parallel.hpp"

// Floating-point `std::to_chars` with a precision: where the standard library
// has it (the guard `fast_number.hpp` uses for `from_chars`).
#if !defined(MESHIOPLUSPLUS_FORCE_NO_FAST_FROM_CHARS) && !defined(__EMSCRIPTEN__) && \
    defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#define MESHIOPLUSPLUS_ROW_WRITER_TO_CHARS 1
#else
#define MESHIOPLUSPLUS_ROW_WRITER_TO_CHARS 0
#endif

namespace meshioplusplus {
namespace detail {

/// `snprintf_c` with the decimal point resolved once (construct it outside
/// the parallel loop): the same bytes, with no `localeconv()` per value.
class CNumber {
public:
    CNumber() : mClassic(!is_c_decimal_point(std::localeconv()->decimal_point)) {}

    /// `snprintf_c(pBuf, Cap, pFmt, ...)`: the same bytes and return value.
    int Print(char* pBuf, std::size_t Cap, const char* pFmt, ...) const {
        va_list args;
        va_start(args, pFmt);
        const int n = mClassic ? vsnprintf_classic(pBuf, Cap, pFmt, args)
                               : std::vsnprintf(pBuf, Cap, pFmt, args);
        va_end(args);
        return n;
    }

    /// Appends `v` formatted by the printf format `pFmt` (one double, any
    /// literal text around it).
    void Append(std::string& rOut, const char* pFmt, double v) const {
#if MESHIOPLUSPLUS_ROW_WRITER_TO_CHARS
        if (std::isfinite(v) && AppendToChars(rOut, pFmt, v))
            return;
#endif
        char buf[64];
        const int n = Print(buf, sizeof buf, pFmt, v);
        if (n <= 0)
            return;
        rOut.append(buf, std::min<std::size_t>(static_cast<std::size_t>(n), sizeof buf - 1));
    }

private:
#if MESHIOPLUSPLUS_ROW_WRITER_TO_CHARS
    /// The fast path for the plain precisions (roadmap 5.7.1): a format that
    /// is a literal prefix, `%.<P>e` or `%.<P>g`, and a literal suffix -- no
    /// flag, width or case -- is `std::to_chars` with that precision, which
    /// the standard specifies to produce what printf does in the "C" locale,
    /// whatever the process locale -- which is what `Print` gives. Any other
    /// format, a non-finite value (printf spells NaN signs its own way)
    /// or text the 64-byte buffer would have truncated takes `snprintf`.
    /// `test_row_writer_to_chars` compares the two over random bit patterns.
    static bool AppendToChars(std::string& rOut, const char* pFmt, double v) {
        const char* pct = std::strchr(pFmt, '%');
        if (pct == nullptr || pct[1] != '.' || pct[2] < '0' || pct[2] > '9')
            return false;
        const char* p = pct + 2;
        int precision = 0;
        while (*p >= '0' && *p <= '9') {
            precision = precision * 10 + (*p - '0');
            if (precision > 40)
                return false;
            ++p;
        }
        std::chars_format format;
        if (*p == 'e')
            format = std::chars_format::scientific;
        else if (*p == 'g')
            format = std::chars_format::general;
        else
            return false;
        const char* suffix = p + 1;
        if (std::strchr(suffix, '%') != nullptr)
            return false;
        char buf[64];
        const std::to_chars_result r = std::to_chars(buf, buf + sizeof buf, v, format, precision);
        if (r.ec != std::errc{})
            return false;
        const std::size_t prefix_len = static_cast<std::size_t>(pct - pFmt);
        const std::size_t number_len = static_cast<std::size_t>(r.ptr - buf);
        const std::size_t suffix_len = strnlen(suffix, sizeof buf);  // > 63 is refused below
        if (prefix_len + number_len + suffix_len > 63)
            return false;
        rOut.append(pFmt, prefix_len);
        rOut.append(buf, number_len);
        rOut.append(suffix, suffix_len);
        return true;
    }
#endif

    bool mClassic = false;  ///< the locale's decimal point is not '.'
};

/// Appends the decimal text of `v`, as `ostream << v` prints it in the
/// classic locale.
inline void append_int(std::string& rOut, std::int64_t v) {
    char buf[24];
    const auto r = std::to_chars(buf, buf + sizeof buf, v);
    rOut.append(buf, static_cast<std::size_t>(r.ptr - buf));
}

/**
 * @brief Formats rows `[0, NumRows)` in parallel and hands the text to
 * `rEmit(const std::string&)` in row order: `rFormat(first, last, rBuf)`
 * appends rows `[first, last)` to `rBuf`. Chunks of @p RowsPerChunk rows are
 * formatted a wave at a time, so the text held in memory stays bounded
 * however many rows there are.
 */
template <class F, class Emit>
void row_chunks(std::size_t NumRows, F&& rFormat, Emit&& rEmit, std::size_t RowsPerChunk = 2048) {
    if (NumRows == 0)
        return;
    const std::size_t chunk = std::max<std::size_t>(1, RowsPerChunk);
    const std::size_t nchunks = (NumRows + chunk - 1) / chunk;
    constexpr std::size_t kWave = 64;
    std::vector<std::string> bufs(std::min(kWave, nchunks));
    for (std::size_t first = 0; first < nchunks; first += kWave) {
        const std::size_t count = std::min(kWave, nchunks - first);
        parallel_for(
            count,
            [&](std::size_t w) {
                const std::size_t c = first + w;
                std::string& buf = bufs[w];
                buf.clear();
                rFormat(c * chunk, std::min(NumRows, (c + 1) * chunk), buf);
            },
            1);
        for (std::size_t w = 0; w < count; ++w)
            rEmit(bufs[w]);
    }
}

/// `row_chunks` written to a stream.
template <class F>
void write_row_chunks(std::ostream& rOs, std::size_t NumRows, F&& rFormat,
                      std::size_t RowsPerChunk = 2048) {
    row_chunks(
        NumRows, rFormat,
        [&](const std::string& rBuf) {
            rOs.write(rBuf.data(), static_cast<std::streamsize>(rBuf.size()));
        },
        RowsPerChunk);
}

/// `row_chunks` appended to a string.
template <class F>
void append_row_chunks(std::string& rOut, std::size_t NumRows, F&& rFormat,
                       std::size_t RowsPerChunk = 2048) {
    row_chunks(NumRows, rFormat, [&](const std::string& rBuf) { rOut += rBuf; }, RowsPerChunk);
}

}  // namespace detail
}  // namespace meshioplusplus
