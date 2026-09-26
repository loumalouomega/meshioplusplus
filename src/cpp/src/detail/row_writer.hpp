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
 * in the classic locale. Roadmap §4, "Parallel row formatting in ASCII
 * writers".
 */

// System includes
#include <algorithm>
#include <charconv>
#include <clocale>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// `snprintf_c` with the decimal point resolved once (construct it outside
/// the parallel loop): the same bytes, with no `localeconv()` per value.
class CNumber {
public:
    CNumber() {
        const char* dp = std::localeconv()->decimal_point;
        // snprintf_c repairs a single-character separator other than '.' and
        // leaves every other case alone; so does this.
        mRepair = dp != nullptr && dp[0] != '\0' && dp[1] == '\0' && dp[0] != '.';
        mDp = mRepair ? dp[0] : '.';
    }

    /// `snprintf_c(pBuf, Cap, pFmt, ...)`: the same bytes and return value.
    int Print(char* pBuf, std::size_t Cap, const char* pFmt, ...) const {
        va_list args;
        va_start(args, pFmt);
        const int n = std::vsnprintf(pBuf, Cap, pFmt, args);
        va_end(args);
        if (n <= 0 || !mRepair || Cap == 0)
            return n;
        const std::size_t written = std::min<std::size_t>(static_cast<std::size_t>(n), Cap - 1);
        for (std::size_t i = 0; i < written; ++i)
            if (pBuf[i] == mDp) {
                pBuf[i] = '.';
                break;
            }
        return n;
    }

    /// Appends `v` formatted by the printf format `pFmt` (one double).
    void Append(std::string& rOut, const char* pFmt, double v) const {
        char buf[64];
        const int n = std::snprintf(buf, sizeof buf, pFmt, v);
        if (n <= 0)
            return;
        const std::size_t len = std::min<std::size_t>(static_cast<std::size_t>(n), sizeof buf - 1);
        if (mRepair)
            for (std::size_t i = 0; i < len; ++i)
                if (buf[i] == mDp) {
                    buf[i] = '.';
                    break;
                }
        rOut.append(buf, len);
    }

private:
    bool mRepair = false;
    char mDp = '.';
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
