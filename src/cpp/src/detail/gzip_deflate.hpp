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
 * @file detail/gzip_deflate.hpp
 * @brief gzip output, byte-identical to Python's `gzip.compress(data, mtime=0)`.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. The reading
 * side is `detail::decompress_bytes` (`detail/compressed_input.hpp`).
 *
 * The layout is the one CPython writes: a 10-byte header with mtime 0, XFL 2
 * and OS 255, raw deflate at level 9, then CRC-32 and the input size, so a
 * native writer and its Python twin produce the same bytes when they share a
 * zlib. libMesh's `.gz` meshes and OpenFOAM's `writeCompression on` cases are
 * written through it.
 */

// System includes
#include <string>
#include <string_view>

namespace meshioplusplus::detail {

/// Whether this build can write gzip (it was built with zlib,
/// `MESHIOPLUSPLUS_WITH_ZLIB`).
bool gzip_compress_available();

/**
 * @brief `In` as one gzip member.
 * @param In The bytes to compress.
 * @param pWhat The format name the errors start with (`"libMesh"`).
 * @throws WriteError When zlib fails, or the build has none.
 */
std::string gzip_compress(std::string_view In, const char* pWhat);

}  // namespace meshioplusplus::detail
