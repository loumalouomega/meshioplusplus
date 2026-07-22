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
#include <algorithm>
#include <cstring>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
#include <zstd.h>
#endif
#ifdef MESHIOPLUSPLUS_HAS_LZ4
#include <lz4.h>
#endif

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

const char* b64_table() {
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
}

std::string b64encode(const unsigned char* pData, std::size_t len) {
    const char* tbl = b64_table();
    // Every 3-byte group maps to 4 output chars at a deterministic offset:
    // pre-size the output and write by index -> parallel over groups.
    const std::size_t ngroups = len / 3;  // full groups
    std::string out(((len + 2) / 3) * 4, '\0');
    parallel_for(ngroups, [&](std::size_t g) {
        const std::size_t i = g * 3;
        unsigned n =
            (unsigned(pData[i]) << 16) | (unsigned(pData[i + 1]) << 8) | unsigned(pData[i + 2]);
        char* o = out.data() + g * 4;
        o[0] = tbl[(n >> 18) & 63];
        o[1] = tbl[(n >> 12) & 63];
        o[2] = tbl[(n >> 6) & 63];
        o[3] = tbl[n & 63];
    });
    const std::size_t i = ngroups * 3;
    if (i < len) {  // trailing 1- or 2-byte group with '=' padding
        const bool two = (i + 1 < len);
        unsigned n = unsigned(pData[i]) << 16;
        if (two)
            n |= unsigned(pData[i + 1]) << 8;
        char* o = out.data() + ngroups * 4;
        o[0] = tbl[(n >> 18) & 63];
        o[1] = tbl[(n >> 12) & 63];
        o[2] = two ? tbl[(n >> 6) & 63] : '=';
        o[3] = '=';
    }
    return out;
}

std::vector<unsigned char> b64decode(const char* pS, std::size_t len) {
    static int8_t inv[256];
    static bool init = false;
    if (!init) {
        for (int i = 0; i < 256; ++i)
            inv[i] = -1;
        const char* tbl = b64_table();
        for (int i = 0; i < 64; ++i)
            inv[(unsigned char)tbl[i]] = static_cast<int8_t>(i);
        init = true;
    }
    std::vector<unsigned char> out;
    out.reserve(len / 4 * 3);
    int buf = 0, bits = 0;
    for (std::size_t i = 0; i < len; ++i) {
        char ch = pS[i];
        if (ch == '=' || ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t')
            continue;
        int v = inv[(unsigned char)ch];
        if (v < 0)
            continue;
        buf = (buf << 6) | v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<unsigned char>((buf >> bits) & 0xFF));
        }
    }
    return out;
}

const char* vtk_codec_compressor(VtkCodec codec) {
    switch (codec) {
        case VtkCodec::Zlib:
            return "vtkZLibDataCompressor";
        case VtkCodec::LZ4:
            return "vtkLZ4DataCompressor";
        case VtkCodec::ZSTD:
            return "vtkZSTDDataCompressor";
        case VtkCodec::LZMA:
            return "vtkLZMADataCompressor";
        default:
            return "";
    }
}

const char* vtk_codec_name(VtkCodec codec) {
    switch (codec) {
        case VtkCodec::Zlib:
            return "zlib";
        case VtkCodec::LZ4:
            return "lz4";
        case VtkCodec::ZSTD:
            return "zstd";
        case VtkCodec::LZMA:
            return "lzma";
        default:
            return "none";
    }
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
/**
 * @brief Compresses one block with zlib's default `compress()` (a single
 * deflate call, no streaming).
 * @throws WriteError if zlib does not return `Z_OK`.
 */
std::vector<unsigned char> zlib_compress_block(const unsigned char* pSrc, std::size_t n) {
    uLongf bound = compressBound(static_cast<uLong>(n));
    std::vector<unsigned char> out(bound);
    uLongf destLen = bound;
    int r = compress(out.data(), &destLen, pSrc, static_cast<uLong>(n));
    if (r != Z_OK)
        throw WriteError("zlib compression failed");
    out.resize(destLen);
    return out;
}

/**
 * @brief Decompresses one zlib-compressed block whose decompressed size is
 * already known.
 * @throws ReadError if zlib does not return `Z_OK`.
 */
std::vector<unsigned char> zlib_decompress(const unsigned char* pSrc, std::size_t n,
                                           std::size_t expected) {
    std::vector<unsigned char> out(expected);
    uLongf destLen = static_cast<uLongf>(expected);
    int r = uncompress(out.data(), &destLen, pSrc, static_cast<uLong>(n));
    if (r != Z_OK)
        throw ReadError("zlib decompression failed");
    out.resize(destLen);
    return out;
}
#endif  // MESHIOPLUSPLUS_HAS_ZLIB

#ifdef MESHIOPLUSPLUS_HAS_ZSTD
/**
 * @brief Compresses one block as a single raw zstd frame.
 *
 * One frame per block, matching the zlib path: VTU's own header already records
 * each block's compressed and decompressed size, so no zstd framing metadata is
 * needed on top.
 * @throws WriteError if zstd reports an error.
 */
std::vector<unsigned char> zstd_compress_block(const unsigned char* pSrc, std::size_t n) {
    const std::size_t bound = ZSTD_compressBound(n);
    std::vector<unsigned char> out(bound);
    const std::size_t written = ZSTD_compress(out.data(), bound, pSrc, n, ZSTD_CLEVEL_DEFAULT);
    if (ZSTD_isError(written))
        throw WriteError(std::string("zstd compression failed: ") + ZSTD_getErrorName(written));
    out.resize(written);
    return out;
}

/**
 * @brief Decompresses one zstd block whose decompressed size is already known.
 * @param expected Exact size from the VTU block header -- so the frame's own
 *        content-size field is never consulted.
 * @throws ReadError if zstd reports an error or the size disagrees.
 */
std::vector<unsigned char> zstd_decompress(const unsigned char* pSrc, std::size_t n,
                                           std::size_t expected) {
    std::vector<unsigned char> out(expected);
    const std::size_t written = ZSTD_decompress(out.data(), expected, pSrc, n);
    if (ZSTD_isError(written))
        throw ReadError(std::string("zstd decompression failed: ") + ZSTD_getErrorName(written));
    out.resize(written);
    return out;
}
#endif  // MESHIOPLUSPLUS_HAS_ZSTD

#ifdef MESHIOPLUSPLUS_HAS_LZ4
/**
 * @brief Compresses one block in LZ4's **raw block** format.
 *
 * Raw block, not the LZ4 *frame* format -- this is what `vtkLZ4DataCompressor`
 * emits (its `SetAccelerationLevel` knob is the tell that it calls
 * `LZ4_compress_fast`), so files written here stay readable by VTK/ParaView.
 * Acceleration 1 is LZ4's default.
 * @throws WriteError if lz4 reports an error.
 */
std::vector<unsigned char> lz4_compress_block(const unsigned char* pSrc, std::size_t n) {
    const int src_size = static_cast<int>(n);
    const int bound = LZ4_compressBound(src_size);
    if (bound <= 0)
        throw WriteError("lz4 compression failed: block too large");
    std::vector<unsigned char> out(static_cast<std::size_t>(bound));
    const int written =
        LZ4_compress_fast(reinterpret_cast<const char*>(pSrc), reinterpret_cast<char*>(out.data()),
                          src_size, bound, /*acceleration=*/1);
    if (written <= 0)
        throw WriteError("lz4 compression failed");
    out.resize(static_cast<std::size_t>(written));
    return out;
}

/**
 * @brief Decompresses one raw-block-format LZ4 block.
 * @param expected Exact decompressed size from the VTU block header. Passing it
 *        as the output capacity is what makes `LZ4_decompress_safe` bounded --
 *        a corrupt block cannot overrun the buffer.
 * @throws ReadError if lz4 reports an error.
 */
std::vector<unsigned char> lz4_decompress(const unsigned char* pSrc, std::size_t n,
                                          std::size_t expected) {
    std::vector<unsigned char> out(expected);
    const int written = LZ4_decompress_safe(reinterpret_cast<const char*>(pSrc),
                                            reinterpret_cast<char*>(out.data()),
                                            static_cast<int>(n), static_cast<int>(expected));
    if (written < 0)
        throw ReadError("lz4 decompression failed");
    out.resize(static_cast<std::size_t>(written));
    return out;
}
#endif  // MESHIOPLUSPLUS_HAS_LZ4

bool vtk_codec_available(VtkCodec codec) {
    switch (codec) {
        case VtkCodec::None:
            return true;
        case VtkCodec::Zlib:
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
            return true;
#else
            return false;
#endif
        case VtkCodec::LZ4:
#ifdef MESHIOPLUSPLUS_HAS_LZ4
            return true;
#else
            return false;
#endif
        case VtkCodec::ZSTD:
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
            return true;
#else
            return false;
#endif
        default:
            return false;  // LZMA is recognized but never implemented here
    }
}

std::string vtk_codec_build_option(VtkCodec codec) {
    switch (codec) {
        case VtkCodec::Zlib:
            return "MESHIOPLUSPLUS_WITH_ZLIB=ON";
        case VtkCodec::LZ4:
            return "MESHIOPLUSPLUS_WITH_LZ4=ON";
        case VtkCodec::ZSTD:
            return "MESHIOPLUSPLUS_WITH_ZSTD=ON";
        default:
            return "";
    }
}

std::string vtk_codec_missing_message(VtkCodec codec, bool for_write) {
    const std::string what = for_write ? "compression" : "decompression";
    if (codec == VtkCodec::LZMA)
        return "VTK XML lzma " + what + " is not implemented by the C++ core";
    return "VTK XML " + std::string(vtk_codec_name(codec)) + " " + what +
           " requires a build with -D" + vtk_codec_build_option(codec);
}

void vtk_codec_require_read(VtkCodec codec) {
    if (!vtk_codec_available(codec))
        throw ReadError(vtk_codec_missing_message(codec, /*for_write=*/false));
}

void vtk_codec_require_write(VtkCodec codec) {
    if (!vtk_codec_available(codec))
        throw WriteError(vtk_codec_missing_message(codec, /*for_write=*/true));
}

std::vector<unsigned char> vtk_codec_compress_block(VtkCodec codec, const unsigned char* pSrc,
                                                    std::size_t n) {
    switch (codec) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        case VtkCodec::Zlib:
            return zlib_compress_block(pSrc, n);
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
        case VtkCodec::ZSTD:
            return zstd_compress_block(pSrc, n);
#endif
#ifdef MESHIOPLUSPLUS_HAS_LZ4
        case VtkCodec::LZ4:
            return lz4_compress_block(pSrc, n);
#endif
        default:
            break;
    }
    (void)pSrc;
    (void)n;
    throw WriteError(vtk_codec_missing_message(codec, /*for_write=*/true));
}

std::vector<unsigned char> vtk_codec_decompress_block(VtkCodec codec, const unsigned char* pSrc,
                                                      std::size_t n, std::size_t expected) {
    switch (codec) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        case VtkCodec::Zlib:
            return zlib_decompress(pSrc, n, expected);
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
        case VtkCodec::ZSTD:
            return zstd_decompress(pSrc, n, expected);
#endif
#ifdef MESHIOPLUSPLUS_HAS_LZ4
        case VtkCodec::LZ4:
            return lz4_decompress(pSrc, n, expected);
#endif
        default:
            break;
    }
    (void)pSrc;
    (void)n;
    (void)expected;
    throw ReadError(vtk_codec_missing_message(codec, /*for_write=*/false));
}

std::uint64_t read_uint_le(const unsigned char* pP, std::size_t isz) {
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < isz; ++i)
        v |= static_cast<std::uint64_t>(pP[i]) << (8 * i);
    return v;
}

std::vector<unsigned char> vtu_decode_uncompressed(const char* pText, std::size_t len,
                                                   std::size_t hsz) {
    std::vector<unsigned char> all = b64decode(pText, len);
    if (all.size() < hsz)
        throw ReadError("VTU binary data too short");
    std::uint64_t total = read_uint_le(all.data(), hsz);
    if (all.size() < hsz + total)
        throw ReadError("VTU binary data truncated");
    return std::vector<unsigned char>(all.begin() + hsz, all.begin() + hsz + total);
}

std::vector<unsigned char> vtu_decode_blocks(const char* pText, std::size_t len, std::size_t hsz,
                                             VtkCodec codec) {
    // Every codec is checked here rather than at the per-block call so an
    // absent one is reported before any work is done -- and as a ReadError,
    // never a link error, which is what keeps the Python fallback reachable.
    vtk_codec_require_read(codec);

    std::size_t first_chars = ((hsz + 2) / 3) * 4;
    if (len < first_chars)
        throw ReadError("VTU compressed-block header too short");
    std::vector<unsigned char> hb = b64decode(pText, first_chars);
    std::uint64_t num_blocks = read_uint_le(hb.data(), hsz);

    std::size_t num_header_bytes = hsz * (3 + static_cast<std::size_t>(num_blocks));
    std::size_t num_header_chars = ((num_header_bytes + 2) / 3) * 4;
    if (len < num_header_chars)
        throw ReadError("VTU compressed-block header truncated");
    std::vector<unsigned char> header = b64decode(pText, num_header_chars);

    std::uint64_t max_block = read_uint_le(header.data() + hsz, hsz);
    std::uint64_t last_block = read_uint_le(header.data() + 2 * hsz, hsz);
    std::vector<std::uint64_t> comp_sizes(num_blocks);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        comp_sizes[k] = read_uint_le(header.data() + (3 + k) * hsz, hsz);

    std::vector<unsigned char> blockdata =
        b64decode(pText + num_header_chars, len - num_header_chars);

    // Input offsets are a (cheap, sequential) prefix sum of comp_sizes; the
    // output offset of block k is k*max_block per the VTU block scheme -> the
    // per-block inflate runs in parallel into a pre-sized buffer.
    std::vector<std::size_t> in_off(static_cast<std::size_t>(num_blocks) + 1, 0);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        in_off[static_cast<std::size_t>(k) + 1] =
            in_off[static_cast<std::size_t>(k)] + static_cast<std::size_t>(comp_sizes[k]);

    const std::size_t total = num_blocks ? static_cast<std::size_t>(num_blocks - 1) *
                                                   static_cast<std::size_t>(max_block) +
                                               static_cast<std::size_t>(last_block)
                                         : 0;
    std::vector<unsigned char> out(total);
    parallel_for(
        static_cast<std::size_t>(num_blocks),
        [&](std::size_t k) {
            std::size_t expected = (k + 1 == num_blocks) ? static_cast<std::size_t>(last_block)
                                                         : static_cast<std::size_t>(max_block);
            auto dec =
                vtk_codec_decompress_block(codec, blockdata.data() + in_off[k],
                                           static_cast<std::size_t>(comp_sizes[k]), expected);
            std::memcpy(out.data() + k * static_cast<std::size_t>(max_block), dec.data(),
                        std::min(dec.size(), expected));
        },
        /*grain=*/1);  // each block is 32 KB of inflate work
    return out;
}

std::string vtu_encode_binary(const unsigned char* pData, std::size_t nbytes, VtkCodec codec) {
    if (codec == VtkCodec::None) {
        std::vector<unsigned char> buf(4 + nbytes);
        std::uint32_t header = static_cast<std::uint32_t>(nbytes);
        std::memcpy(buf.data(), &header, 4);
        if (nbytes)
            std::memcpy(buf.data() + 4, pData, nbytes);
        return b64encode(buf.data(), buf.size());
    }

    vtk_codec_require_write(codec);
    const std::uint32_t max_block = 32768;
    std::uint32_t num_blocks = static_cast<std::uint32_t>((nbytes + max_block - 1) / max_block);
    std::uint32_t last_block_size =
        num_blocks ? static_cast<std::uint32_t>(nbytes - std::size_t(num_blocks - 1) * max_block)
                   : max_block;

    // Blocks are independent -> compress in parallel into pre-sized slots.
    std::vector<std::vector<unsigned char> > blocks(num_blocks);
    parallel_for(
        num_blocks,
        [&](std::size_t b) {
            std::size_t off = b * max_block;
            std::size_t len = std::min<std::size_t>(max_block, nbytes - off);
            blocks[b] = vtk_codec_compress_block(codec, pData + off, len);
        },
        /*grain=*/1);  // each block is 32 KB of deflate work

    std::vector<std::uint32_t> header;
    header.reserve(3 + num_blocks);
    header.push_back(num_blocks);
    header.push_back(max_block);
    header.push_back(last_block_size);
    for (const auto& b : blocks)
        header.push_back(static_cast<std::uint32_t>(b.size()));

    std::string out = b64encode(reinterpret_cast<const unsigned char*>(header.data()),
                                header.size() * sizeof(std::uint32_t));
    std::vector<unsigned char> concat;
    for (const auto& b : blocks)
        concat.insert(concat.end(), b.begin(), b.end());
    out += b64encode(concat.data(), concat.size());
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
