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
// System includes
#include <cstdint>
#include <string>
#include <string_view>

// Project includes
#include "gzip_deflate.hpp"
#include "meshioplusplus/exceptions.hpp"

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

namespace meshioplusplus::detail {

bool gzip_compress_available() {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    return true;
#else
    return false;
#endif
}

std::string gzip_compress(std::string_view In, const char* pWhat) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    z_stream z{};
    if (deflateInit2(&z, Z_BEST_COMPRESSION, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw WriteError(std::string(pWhat) + ": zlib deflateInit2 failed");
    std::string out("\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff", 10);
    out.resize(10 + deflateBound(&z, static_cast<uLong>(In.size())));
    z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(In.data()));
    z.avail_in = static_cast<uInt>(In.size());
    z.next_out = reinterpret_cast<Bytef*>(&out[10]);
    z.avail_out = static_cast<uInt>(out.size() - 10);
    const int rc = deflate(&z, Z_FINISH);
    const std::size_t produced = 10 + z.total_out;
    deflateEnd(&z);
    if (rc != Z_STREAM_END)
        throw WriteError(std::string(pWhat) + ": zlib deflate failed");
    out.resize(produced);
    const uLong crc =
        crc32(0L, reinterpret_cast<const Bytef*>(In.data()), static_cast<uInt>(In.size()));
    for (const std::uint32_t v :
         {static_cast<std::uint32_t>(crc), static_cast<std::uint32_t>(In.size())})
        for (int b = 0; b < 4; ++b)
            out += static_cast<char>((v >> (8 * b)) & 0xff);
    return out;
#else
    (void)In;
    throw WriteError(std::string(pWhat) +
                     ": this build has no zlib (MESHIOPLUSPLUS_WITH_ZLIB) to write .gz");
#endif
}

}  // namespace meshioplusplus::detail
