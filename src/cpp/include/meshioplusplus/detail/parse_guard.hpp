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
 * @file detail/parse_guard.hpp
 * @brief Checks a reader makes before trusting what a file says about itself.
 *
 * The fuzz campaign (doc/fuzzing.md) finds the same two defects in reader
 * after reader: a token indexed before the line was checked to have it, and a
 * count from a header used to size an allocation before the bytes that should
 * back it were seen. Both turn a malformed file into an out-of-bounds read or
 * an out-of-memory abort instead of a `ReadError`. These helpers make the
 * check one line at the site.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus::detail {

/**
 * @brief Throws ReadError unless @p rTokens holds at least @p Count entries.
 * @param rTokens Any container with `size()`.
 * @param Count The number of entries the caller is about to index.
 * @param pFormat Format name for the message.
 */
template <class TContainer>
void need_tokens(const TContainer& rTokens, std::size_t Count, const char* pFormat) {
    if (rTokens.size() < Count)
        throw ReadError(std::string("meshio++: ") + pFormat + ": expected " +
                        std::to_string(Count) + " fields on a line, got " +
                        std::to_string(rTokens.size()));
}

/**
 * @brief A count read from a file, checked against what the file can hold.
 *
 * @param Value The count as parsed (any integer type).
 * @param Limit The most the rest of the file could possibly describe (lines
 *        or bytes left, divided by the minimum size of one entry).
 * @param pFormat Format name for the message.
 * @param pWhat What is being counted, for the message.
 * @return @p Value as a size.
 * @throws ReadError if @p Value is negative or exceeds @p Limit.
 */
template <class TInt>
std::size_t checked_count(TInt Value, std::size_t Limit, const char* pFormat, const char* pWhat) {
    if constexpr (std::numeric_limits<TInt>::is_signed) {
        if (Value < 0)
            throw ReadError(std::string("meshio++: ") + pFormat + ": negative " + pWhat +
                            " count " + std::to_string(Value));
    }
    if (static_cast<unsigned long long>(Value) > Limit)
        throw ReadError(std::string("meshio++: ") + pFormat + ": " + pWhat + " count " +
                        std::to_string(Value) + " exceeds what the file holds (" +
                        std::to_string(Limit) + ")");
    return static_cast<std::size_t>(Value);
}

/**
 * @brief A floating-point value converted to an integer type, or ReadError.
 *
 * Casting a double outside the target's range (or NaN) is undefined
 * behaviour; formats that store integers as reals (Ansys `.rst` records,
 * `.cdb` fields written with an exponent) must go through this.
 */
template <class TInt>
TInt checked_integer(double Value, const char* pFormat) {
    // min() and max() + 1 are powers of two, so both bounds are exact doubles
    // and the comparison is exact; the negated form also rejects NaN.
    constexpr double lo = static_cast<double>(std::numeric_limits<TInt>::min());
    constexpr double hi = static_cast<double>(std::numeric_limits<TInt>::max()) + 1.0;
    if (!(Value >= lo && Value < hi))
        throw ReadError(std::string("meshio++: ") + pFormat + ": integer field out of range");
    return static_cast<TInt>(Value);
}

/**
 * @brief A file's 1-based id as a 0-based index, without overflow.
 *
 * `id - 1` is undefined for `INT64_MIN`; modular arithmetic wraps it to a
 * value the reader's own range check then rejects.
 */
inline std::int64_t zero_based(long long Id) {
    return static_cast<std::int64_t>(static_cast<std::uint64_t>(Id) - 1u);
}

/**
 * @brief The size of @p rPath in bytes, or 0 when it cannot be read.
 *
 * The bound for a text format's header counts: every entity takes at least
 * one byte of the file, usually several.
 */
inline std::size_t file_bytes(const std::string& rPath) {
    std::error_code ec;
    const auto n = std::filesystem::file_size(rPath, ec);
    return ec ? 0 : static_cast<std::size_t>(n);
}

}  // namespace meshioplusplus::detail
