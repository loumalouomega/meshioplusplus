#pragma once
//
// Byte-swap primitives backed by the single-instruction bswap intrinsics.
//
// The hand-written per-byte reversal loops previously used in the binary
// readers/writers (`for (b) buf[i*isz+b] = src[i*isz+(isz-1-b)]`) are not
// reliably recognised by compilers as a bswap, so endianness conversion was
// running at ~1 byte per cycle. These helpers compile to one instruction per
// element. Semantics are bit-identical to the manual reversal.

#include <cstdint>
#include <cstring>

#ifdef _MSC_VER
#include <stdlib.h>
#endif

namespace meshioplusplus {
namespace detail {

inline std::uint16_t bswap16(std::uint16_t v) {
#if defined(_MSC_VER)
    return _byteswap_ushort(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap16(v);
#else
    return static_cast<std::uint16_t>((v << 8) | (v >> 8));
#endif
}

inline std::uint32_t bswap32(std::uint32_t v) {
#if defined(_MSC_VER)
    return _byteswap_ulong(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap32(v);
#else
    return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8) |
           ((v & 0x00FF0000u) >> 8) | ((v & 0xFF000000u) >> 24);
#endif
}

inline std::uint64_t bswap64(std::uint64_t v) {
#if defined(_MSC_VER)
    return _byteswap_uint64(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(v);
#else
    return (static_cast<std::uint64_t>(bswap32(static_cast<std::uint32_t>(v))) << 32) |
           bswap32(static_cast<std::uint32_t>(v >> 32));
#endif
}

// Reverse `n` bytes (n in {1,2,4,8}) from src into dst. dst and src must not
// partially overlap (dst == src is fine via the temporaries).
inline void bswap_copy(char* dst, const char* src, int n) {
    switch (n) {
        case 8: {
            std::uint64_t v;
            std::memcpy(&v, src, 8);
            v = bswap64(v);
            std::memcpy(dst, &v, 8);
            break;
        }
        case 4: {
            std::uint32_t v;
            std::memcpy(&v, src, 4);
            v = bswap32(v);
            std::memcpy(dst, &v, 4);
            break;
        }
        case 2: {
            std::uint16_t v;
            std::memcpy(&v, src, 2);
            v = bswap16(v);
            std::memcpy(dst, &v, 2);
            break;
        }
        default:  // n == 1 (or unexpected): plain copy
            if (dst != src) std::memcpy(dst, src, static_cast<std::size_t>(n));
            break;
    }
}

// In-place variant.
inline void bswap_inplace(char* p, int n) { bswap_copy(p, p, n); }

}  // namespace detail
}  // namespace meshioplusplus
