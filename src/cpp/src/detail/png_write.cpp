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

// The dependency-free PNG encoder. See detail/png_write.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

// Project includes (private, not installed)
#include "png_write.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

std::array<std::uint32_t, 256> png_crc_table() {
    std::array<std::uint32_t, 256> table{};
    for (std::uint32_t n = 0; n < 256; ++n) {
        std::uint32_t c = n;
        for (int k = 0; k < 8; ++k)
            c = (c & 1U) ? (0xEDB88320U ^ (c >> 1)) : (c >> 1);
        table[n] = c;
    }
    return table;
}

void png_put_u32(std::string& rOut, std::uint32_t v) {
    rOut.push_back(static_cast<char>((v >> 24) & 0xFFU));
    rOut.push_back(static_cast<char>((v >> 16) & 0xFFU));
    rOut.push_back(static_cast<char>((v >> 8) & 0xFFU));
    rOut.push_back(static_cast<char>(v & 0xFFU));
}

void png_chunk(std::string& rOut, const char* pType, const std::string& rData) {
    png_put_u32(rOut, static_cast<std::uint32_t>(rData.size()));
    const std::size_t start = rOut.size();
    rOut.append(pType, 4);
    rOut += rData;
    const std::uint32_t crc = png_crc32(
        0, reinterpret_cast<const unsigned char*>(rOut.data() + start), rOut.size() - start);
    png_put_u32(rOut, crc);
}

// A zlib stream of stored deflate blocks (RFC 1950/1951): no compression, so
// the output depends on nothing but the input.
std::string png_zlib_stored(const std::vector<unsigned char>& rRaw) {
    std::string out;
    out.push_back(static_cast<char>(0x78));  // CM 8, a 32 KiB window
    out.push_back(static_cast<char>(0x01));  // FCHECK: (0x78 << 8 | 0x01) % 31 == 0
    const std::size_t max_block = 65535;
    std::size_t pos = 0;
    do {
        const std::size_t len = std::min(max_block, rRaw.size() - pos);
        const bool last = pos + len == rRaw.size();
        out.push_back(static_cast<char>(last ? 1 : 0));
        out.push_back(static_cast<char>(len & 0xFFU));
        out.push_back(static_cast<char>((len >> 8) & 0xFFU));
        out.push_back(static_cast<char>(~len & 0xFFU));
        out.push_back(static_cast<char>((~len >> 8) & 0xFFU));
        out.append(reinterpret_cast<const char*>(rRaw.data() + pos), len);
        pos += len;
    } while (pos < rRaw.size());
    png_put_u32(out, png_adler32(1, rRaw.data(), rRaw.size()));
    return out;
}

}  // namespace

std::uint32_t png_crc32(std::uint32_t Crc, const unsigned char* pData, std::size_t Size) {
    static const std::array<std::uint32_t, 256> table = png_crc_table();
    std::uint32_t c = Crc ^ 0xFFFFFFFFU;
    for (std::size_t i = 0; i < Size; ++i)
        c = table[(c ^ pData[i]) & 0xFFU] ^ (c >> 8);
    return c ^ 0xFFFFFFFFU;
}

std::uint32_t png_adler32(std::uint32_t Adler, const unsigned char* pData, std::size_t Size) {
    std::uint32_t a = Adler & 0xFFFFU;
    std::uint32_t b = (Adler >> 16) & 0xFFFFU;
    // 5552 is the largest run whose sums cannot overflow 32 bits before the modulo.
    std::size_t i = 0;
    while (i < Size) {
        const std::size_t end = std::min(Size, i + 5552);
        for (; i < end; ++i) {
            a += pData[i];
            b += a;
        }
        a %= 65521U;
        b %= 65521U;
    }
    return (b << 16) | a;
}

std::string png_encode_rgba(const unsigned char* pRgba, std::size_t Width, std::size_t Height,
                            int CompressLevel) {
    if (CompressLevel < 0 || CompressLevel > 9)
        throw std::invalid_argument("meshio++: png: compression level must lie in [0, 9]");
    if (Width == 0 || Height == 0 || Width > 0x7FFFFFFFU || Height > 0x7FFFFFFFU)
        throw std::invalid_argument("meshio++: png: width and height must be positive");
    // Each row is a filter byte (0, none) and its pixels.
    const std::size_t row = Width * 4;
    std::vector<unsigned char> raw((row + 1) * Height);
    for (std::size_t y = 0; y < Height; ++y) {
        raw[y * (row + 1)] = 0;
        std::copy(pRgba + y * row, pRgba + (y + 1) * row, raw.data() + y * (row + 1) + 1);
    }

    std::string idat;
    if (CompressLevel == 0) {
        idat = png_zlib_stored(raw);
    } else {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        uLongf bound = compressBound(static_cast<uLong>(raw.size()));
        idat.resize(bound);
        const int rc = compress2(reinterpret_cast<Bytef*>(idat.data()), &bound, raw.data(),
                                 static_cast<uLong>(raw.size()), CompressLevel);
        if (rc != Z_OK)
            throw std::runtime_error("meshio++: png: zlib compression failed");
        idat.resize(bound);
#else
        throw std::invalid_argument(
            "meshio++: png: compression levels 1-9 need zlib (MESHIOPLUSPLUS_WITH_ZLIB); use 0");
#endif
    }

    std::string out("\x89PNG\r\n\x1a\n", 8);
    std::string ihdr;
    png_put_u32(ihdr, static_cast<std::uint32_t>(Width));
    png_put_u32(ihdr, static_cast<std::uint32_t>(Height));
    ihdr.push_back(static_cast<char>(8));  // bit depth
    ihdr.push_back(static_cast<char>(6));  // colour type: RGBA
    ihdr.push_back(static_cast<char>(0));  // compression: deflate
    ihdr.push_back(static_cast<char>(0));  // filter method 0
    ihdr.push_back(static_cast<char>(0));  // no interlace
    png_chunk(out, "IHDR", ihdr);
    png_chunk(out, "IDAT", idat);
    png_chunk(out, "IEND", std::string());
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
