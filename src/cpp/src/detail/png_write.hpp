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
 * @file detail/png_write.hpp
 * @brief A dependency-free PNG encoder: RGBA8 rows into stored (uncompressed)
 * deflate blocks, with a table-driven CRC-32 and an Adler-32 written here
 * rather than taken from zlib, so the bytes are the same on every platform and
 * in a build without zlib.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent).
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>

namespace meshioplusplus {
namespace detail {

/// CRC-32 (the IEEE polynomial, reflected, as PNG and zlib use), continuing
/// from @p Crc (start at 0).
std::uint32_t png_crc32(std::uint32_t Crc, const unsigned char* pData, std::size_t Size);

/// Adler-32, continuing from @p Adler (start at 1).
std::uint32_t png_adler32(std::uint32_t Adler, const unsigned char* pData, std::size_t Size);

/**
 * @brief A complete PNG file of an RGBA8 image, rows top to bottom.
 * @param CompressLevel 0 for stored blocks; 1-9 for zlib's deflate (throws
 *        `std::invalid_argument` in a build without zlib)
 */
std::string png_encode_rgba(const unsigned char* pRgba, std::size_t Width, std::size_t Height,
                            int CompressLevel);

}  // namespace detail
}  // namespace meshioplusplus
