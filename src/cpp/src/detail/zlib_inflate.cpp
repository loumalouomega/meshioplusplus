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
// Streaming inflate of one deflate stream, shared by the GiD and FEBio readers.

// System includes
#include <string>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

// Project includes
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

bool zlib_available() {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    return true;
#else
    return false;
#endif
}

std::string zlib_inflate(std::string_view In, int WindowBits, std::size_t* pConsumed,
                         const char* pWhat) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    z_stream strm{};
    if (inflateInit2(&strm, WindowBits) != Z_OK)
        throw ReadError(std::string(pWhat) + ": could not initialize decompression");
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(In.data()));
    strm.avail_in = static_cast<uInt>(In.size());

    std::string out;
    std::vector<char> chunk(1 << 16);
    int rc = Z_OK;
    do {
        strm.next_out = reinterpret_cast<Bytef*>(chunk.data());
        strm.avail_out = static_cast<uInt>(chunk.size());
        rc = inflate(&strm, Z_NO_FLUSH);
        if (rc != Z_OK && rc != Z_STREAM_END) {
            inflateEnd(&strm);
            throw ReadError(std::string(pWhat) + ": decompression failed (corrupt or truncated)");
        }
        out.append(chunk.data(), chunk.size() - strm.avail_out);
    } while (rc != Z_STREAM_END);
    if (pConsumed)
        *pConsumed = In.size() - strm.avail_in;
    inflateEnd(&strm);
    return out;
#else
    (void)In;
    (void)WindowBits;
    (void)pConsumed;
    throw ReadError(std::string(pWhat) +
                    ": compressed data needs a build with -DMESHIOPLUSPLUS_WITH_ZLIB=ON");
#endif
}

}  // namespace detail
}  // namespace meshioplusplus
