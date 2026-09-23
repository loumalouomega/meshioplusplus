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
 * @file detail/zlib_inflate.hpp
 * @brief Inflate one deflate stream whose decompressed size is not known up
 * front.
 *
 * `uncompress()` needs the output size and speaks only the zlib wrapper; GiD's
 * gzipped files and FEBio's compressed `.xplt` chunks have neither the size nor,
 * in the second case, a single stream per file. This inflates exactly one
 * stream from the start of the input and reports how many input bytes it took,
 * so a caller can find the next stream right after it.
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Whether this build can inflate (`MESHIOPLUSPLUS_WITH_ZLIB`).
 */
MESHIOPLUSPLUS_API bool zlib_available();

/**
 * @brief Inflate the one deflate stream at the start of @p In.
 * @param In The compressed bytes; anything after the stream's end is left alone.
 * @param WindowBits zlib's `windowBits`: 15 for the zlib wrapper, 31 for gzip.
 * @param pConsumed When not null, receives the number of input bytes the stream
 *        occupied.
 * @param pWhat Names the caller in error messages (e.g. `"GiD"`).
 * @return The decompressed bytes.
 * @throws ReadError when the stream is corrupt or truncated, or the build has no
 *         zlib.
 */
MESHIOPLUSPLUS_API std::string zlib_inflate(std::string_view In, int WindowBits,
                                            std::size_t* pConsumed, const char* pWhat);

}  // namespace detail
}  // namespace meshioplusplus
