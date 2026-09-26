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
#include <array>
#include <cstring>
#include <limits>
#include <string>

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

// Project includes (private, not installed)
#include "vtu_decode.hpp"

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

namespace {

/** @brief The inverse alphabet: a character's 6-bit value, or -1 when it is skipped. */
const std::array<int8_t, 256>& vtub_b64_inverse() {
    // A magic static, initialised once and thread-safely ([stmt.dcl]/4): the
    // hand-rolled "static bool init" it replaces was a data race when two
    // threads decoded at once (C, Fortran, Julia or R callers).
    static const std::array<int8_t, 256> inv = [] {
        std::array<int8_t, 256> t;
        t.fill(-1);
        const char* tbl = b64_table();
        for (int i = 0; i < 64; ++i)
            t[(unsigned char)tbl[i]] = static_cast<int8_t>(i);
        return t;
    }();
    return inv;
}

/// Text bytes per decode chunk: a constant, so the chunking never depends on
/// the thread count (the output does not either -- each chunk writes its own
/// byte range -- but a constant keeps the work split reproducible).
constexpr std::size_t kVtubDecodeChunk = std::size_t(1) << 20;

// Decodes the 4-character groups [first, last) of the valid-character stream,
// reading the text from position `from`, and stores only the decoded bytes
// [Lo, Hi) of the stream, at pOut[0 .. Hi - Lo): a group wholly inside the
// window stores its bytes directly, a group across either end byte by byte.
void vtub_decode_groups_window(const std::array<int8_t, 256>& rInv, const char* pS, std::size_t len,
                               std::size_t from, std::size_t first, std::size_t last,
                               std::size_t Lo, std::size_t Hi, unsigned char* pOut) {
    std::size_t at = from;
    for (std::size_t k = first; k < last; k += 4) {
        unsigned n = 0;
        int got = 0;
        while (got < 4 && at < len) {
            const int v = rInv[(unsigned char)pS[at++]];
            if (v < 0)
                continue;
            n = (n << 6) | static_cast<unsigned>(v);
            ++got;
        }
        const std::size_t o = k / 4 * 3;
        unsigned char b[3];
        std::size_t nb = 0;
        if (got == 4) {
            b[0] = static_cast<unsigned char>(n >> 16);
            b[1] = static_cast<unsigned char>(n >> 8);
            b[2] = static_cast<unsigned char>(n);
            nb = 3;
        } else if (got == 3) {
            b[0] = static_cast<unsigned char>(n >> 10);
            b[1] = static_cast<unsigned char>(n >> 2);
            nb = 2;
        } else if (got == 2) {
            b[0] = static_cast<unsigned char>(n >> 4);
            nb = 1;
        }
        if (o >= Lo && o + nb <= Hi) {
            std::memcpy(pOut + (o - Lo), b, nb);
        } else {
            for (std::size_t i = 0; i < nb; ++i)
                if (o + i >= Lo && o + i < Hi)
                    pOut[o + i - Lo] = b[i];
        }
    }
}

}  // namespace

VtubB64::VtubB64(const char* pS, std::size_t len) : mpS(pS), mLen(len) {
    // Every character outside the alphabet -- '=' padding, whitespace, line
    // breaks, anything else -- is skipped, and the valid characters form one
    // bit stream: m of them decode to floor(6m / 8) bytes. Line-wrapped base64
    // relies on that. The text is cut into fixed chunks and each chunk's valid
    // characters counted (in parallel); the prefix sum tells every chunk which
    // 4-character groups start inside it and where their bytes go.
    const std::array<int8_t, 256>& inv = vtub_b64_inverse();
    const std::size_t nchunks = (len + kVtubDecodeChunk - 1) / kVtubDecodeChunk;
    mValid.assign(nchunks + 1, 0);
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            const std::size_t lo = c * kVtubDecodeChunk;
            const std::size_t hi = std::min(len, lo + kVtubDecodeChunk);
            std::size_t cnt = 0;
            for (std::size_t i = lo; i < hi; ++i)
                cnt += inv[(unsigned char)pS[i]] >= 0;
            mValid[c + 1] = cnt;
        },
        1);
    for (std::size_t c = 0; c < nchunks; ++c)
        mValid[c + 1] += mValid[c];
    const std::size_t m = mValid[nchunks];
    mSize = m / 4 * 3 + (m % 4 == 3 ? 2 : m % 4 == 2 ? 1 : 0);
}

void VtubB64::Decode(std::size_t Lo, std::size_t Hi, unsigned char* pOut) const {
    Hi = std::min(Hi, mSize);
    if (Lo >= Hi)
        return;
    const std::array<int8_t, 256>& inv = vtub_b64_inverse();
    const std::size_t nchunks = mValid.size() - 1;
    const std::size_t m = mValid[nchunks];
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            // The first group starting in this chunk is the first multiple of
            // four at or after its first valid character; skip the valid
            // characters before it (they finish the previous chunk's group).
            const std::size_t first = (mValid[c] + 3) / 4 * 4;
            const std::size_t last = std::min(m, mValid[c + 1]);
            if (first >= last)
                return;
            // The chunk's groups decode bytes [first/4*3, ceil(last/4)*3).
            if ((last + 3) / 4 * 3 <= Lo || first / 4 * 3 >= Hi)
                return;
            std::size_t at = c * kVtubDecodeChunk;
            for (std::size_t skip = first - mValid[c]; skip > 0; ++at)
                skip -= inv[(unsigned char)mpS[at]] >= 0;
            vtub_decode_groups_window(inv, mpS, mLen, at, first, last, Lo, Hi, pOut);
        },
        1);
}

std::vector<unsigned char> b64decode(const char* pS, std::size_t len) {
    const VtubB64 stream(pS, len);
    std::vector<unsigned char> out(stream.Size());
    stream.Decode(0, out.size(), out.data());
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
    // The expected size comes from the file's block header and sizes the
    // output buffer. No codec expands a block by more than about 2^16 (zlib's
    // ceiling is ~1032:1), so a larger claim is corruption, not an allocation.
    if (expected > (std::size_t{1} << 16) * (n + 1) + (std::size_t{1} << 20))
        throw ReadError("VTK compressed block declares an implausible decompressed size");
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
    if (hb.size() < hsz)
        throw ReadError("VTU compressed-block header too short");
    std::uint64_t num_blocks = read_uint_le(hb.data(), hsz);
    // Every block has an hsz-byte size in the header, so the count is bounded
    // by the text (and the header size below cannot wrap).
    if (num_blocks > len / hsz)
        throw ReadError("VTU compressed-block header declares more blocks than it holds");

    std::size_t num_header_bytes = hsz * (3 + static_cast<std::size_t>(num_blocks));
    std::size_t num_header_chars = ((num_header_bytes + 2) / 3) * 4;
    if (len < num_header_chars)
        throw ReadError("VTU compressed-block header truncated");
    std::vector<unsigned char> header = b64decode(pText, num_header_chars);
    if (header.size() < num_header_bytes)
        throw ReadError("VTU compressed-block header truncated");

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
    for (std::uint64_t k = 0; k < num_blocks; ++k) {
        const std::size_t at = in_off[static_cast<std::size_t>(k)];
        if (comp_sizes[k] > blockdata.size() - at)  // not at + size: that can wrap
            throw ReadError("VTU compressed blocks are larger than the data that follows");
        in_off[static_cast<std::size_t>(k) + 1] = at + static_cast<std::size_t>(comp_sizes[k]);
    }

    // No codec expands a block by more than about 2^16 (zlib's ceiling is
    // ~1032:1); a header claiming more is corrupt, and would size `out`.
    const std::uint64_t ceiling = (std::uint64_t{1} << 16) * (blockdata.size() + 1) + (1u << 20);
    if (num_blocks && (max_block > ceiling || last_block > ceiling ||
                       (num_blocks - 1) > ceiling / std::max<std::uint64_t>(max_block, 1)))
        throw ReadError("VTU compressed-block header declares an implausible size");
    const std::size_t total = num_blocks ? static_cast<std::size_t>(num_blocks - 1) *
                                                   static_cast<std::size_t>(max_block) +
                                               static_cast<std::size_t>(last_block)
                                         : 0;
    if (total > ceiling)
        throw ReadError("VTU compressed-block header declares an implausible size");
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

namespace {

// b64encode of the concatenation `pHead[0..nHead) + pData[0..n)` without
// building it: the byte stream is read in place, groups wholly inside the data
// straight from it (the uncompressed DataArray's header and payload are one
// base64 unit). The same text as encoding the concatenation.
std::string vtub_b64encode_pair(const unsigned char* pHead, std::size_t nHead,
                                const unsigned char* pData, std::size_t n) {
    const char* tbl = b64_table();
    const std::size_t len = nHead + n;
    const auto at = [&](std::size_t i) -> unsigned {
        return i < nHead ? pHead[i] : pData[i - nHead];
    };
    const std::size_t ngroups = len / 3;
    std::string out(((len + 2) / 3) * 4, '\0');
    parallel_for(ngroups, [&](std::size_t g) {
        const std::size_t i = g * 3;
        unsigned v;
        if (i >= nHead) {
            const unsigned char* d = pData + (i - nHead);
            v = (unsigned(d[0]) << 16) | (unsigned(d[1]) << 8) | unsigned(d[2]);
        } else {
            v = (at(i) << 16) | (at(i + 1) << 8) | at(i + 2);
        }
        char* o = out.data() + g * 4;
        o[0] = tbl[(v >> 18) & 63];
        o[1] = tbl[(v >> 12) & 63];
        o[2] = tbl[(v >> 6) & 63];
        o[3] = tbl[v & 63];
    });
    const std::size_t i = ngroups * 3;
    if (i < len) {  // trailing 1- or 2-byte group with '=' padding
        const bool two = (i + 1 < len);
        unsigned v = at(i) << 16;
        if (two)
            v |= at(i + 1) << 8;
        char* o = out.data() + ngroups * 4;
        o[0] = tbl[(v >> 18) & 63];
        o[1] = tbl[(v >> 12) & 63];
        o[2] = two ? tbl[(v >> 6) & 63] : '=';
        o[3] = '=';
    }
    return out;
}

}  // namespace

std::string vtu_encode_binary(const unsigned char* pData, std::size_t nbytes, VtkCodec codec) {
    return vtu_encode_binary(pData, nbytes, codec, 4);
}

std::size_t vtu_header_bytes_for(std::uint64_t maxArrayBytes) {
    return maxArrayBytes > std::numeric_limits<std::uint32_t>::max() ? 8 : 4;
}

std::string vtu_encode_binary(const unsigned char* pData, std::size_t nbytes, VtkCodec codec,
                              std::size_t hsz) {
    if (hsz != 4 && hsz != 8)
        throw WriteError("VTK XML: header_type must be 4 or 8 bytes, got " + std::to_string(hsz));
    // One little-endian header item of hsz bytes. A size that does not fit is
    // refused: a UInt32 header silently truncated a 4 GiB array's byte count.
    auto put = [hsz](std::vector<unsigned char>& rOut, std::uint64_t Value) {
        if (hsz == 4 && Value > std::numeric_limits<std::uint32_t>::max())
            throw WriteError("VTK XML: a size of " + std::to_string(Value) +
                             " does not fit a UInt32 header_type (the writer must choose UInt64)");
        for (std::size_t b = 0; b < hsz; ++b)
            rOut.push_back(static_cast<unsigned char>((Value >> (8 * b)) & 0xFF));
    };

    if (codec == VtkCodec::None) {
        std::vector<unsigned char> head;
        put(head, nbytes);  // refuses before anything is allocated
        return vtub_b64encode_pair(head.data(), head.size(), pData, nbytes);
    }

    vtk_codec_require_write(codec);
    const std::size_t max_block = 32768;
    const std::size_t num_blocks = (nbytes + max_block - 1) / max_block;
    const std::size_t last_block_size =
        num_blocks ? nbytes - (num_blocks - 1) * max_block : max_block;

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

    std::vector<unsigned char> header;
    header.reserve((3 + num_blocks) * hsz);
    put(header, num_blocks);
    put(header, max_block);
    put(header, last_block_size);
    std::size_t total = 0;
    for (const auto& b : blocks) {
        put(header, b.size());
        total += b.size();
    }

    std::string out = b64encode(header.data(), header.size());
    std::vector<unsigned char> concat;
    concat.reserve(total);
    for (const auto& b : blocks)
        concat.insert(concat.end(), b.begin(), b.end());
    out += b64encode(concat.data(), concat.size());
    return out;
}

namespace {

// Decompresses one block into `pDst` (capacity `Expected`, the size the VTU
// block header gives), with the same checks and messages as
// vtk_codec_decompress_block; returns the bytes written.
std::size_t vtub_decompress_into(VtkCodec codec, const unsigned char* pSrc, std::size_t n,
                                 unsigned char* pDst, std::size_t Expected) {
    if (Expected > (std::size_t{1} << 16) * (n + 1) + (std::size_t{1} << 20))
        throw ReadError("VTK compressed block declares an implausible decompressed size");
    switch (codec) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        case VtkCodec::Zlib: {
            uLongf dest_len = static_cast<uLongf>(Expected);
            if (uncompress(pDst, &dest_len, pSrc, static_cast<uLong>(n)) != Z_OK)
                throw ReadError("zlib decompression failed");
            return static_cast<std::size_t>(dest_len);
        }
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
        case VtkCodec::ZSTD: {
            const std::size_t written = ZSTD_decompress(pDst, Expected, pSrc, n);
            if (ZSTD_isError(written))
                throw ReadError(std::string("zstd decompression failed: ") +
                                ZSTD_getErrorName(written));
            return written;
        }
#endif
#ifdef MESHIOPLUSPLUS_HAS_LZ4
        case VtkCodec::LZ4: {
            const int written = LZ4_decompress_safe(
                reinterpret_cast<const char*>(pSrc), reinterpret_cast<char*>(pDst),
                static_cast<int>(n), static_cast<int>(Expected));
            if (written < 0)
                throw ReadError("lz4 decompression failed");
            return static_cast<std::size_t>(written);
        }
#endif
        default:
            break;
    }
    (void)pSrc;
    (void)n;
    (void)pDst;
    throw ReadError(vtk_codec_missing_message(codec, /*for_write=*/false));
}

// The payload of an uncompressed DataArray, straight into a fresh array:
// vtu_decode_uncompressed's checks, without its two intermediate buffers.
NDArray vtub_uncompressed_ndarray(const char* pText, std::size_t len, std::size_t hsz, DType dt) {
    const VtubB64 stream(pText, len);
    if (stream.Size() < hsz)
        throw ReadError("VTU binary data too short");
    unsigned char head[8] = {0};
    stream.Decode(0, hsz, head);
    const std::uint64_t total = read_uint_le(head, hsz);
    if (total > stream.Size() - hsz)
        throw ReadError("VTU binary data truncated");
    const std::size_t isz = dtype_size(dt);
    const std::size_t n = isz ? static_cast<std::size_t>(total) / isz : 0;
    NDArray a = NDArray::Uninit(dt, {n});
    if (n)
        stream.Decode(hsz, hsz + n * isz, reinterpret_cast<unsigned char*>(a.Data()));
    return a;
}

}  // namespace

NDArray vtu_decode_ndarray(const char* pText, std::size_t len, std::size_t hsz, VtkCodec codec,
                           DType dt) {
    if (codec == VtkCodec::None)
        return vtub_uncompressed_ndarray(pText, len, hsz, dt);

    // vtu_decode_blocks, decompressing each block straight into the array.
    vtk_codec_require_read(codec);
    const std::size_t first_chars = ((hsz + 2) / 3) * 4;
    if (len < first_chars)
        throw ReadError("VTU compressed-block header too short");
    const std::vector<unsigned char> hb = b64decode(pText, first_chars);
    if (hb.size() < hsz)
        throw ReadError("VTU compressed-block header too short");
    const std::uint64_t num_blocks = read_uint_le(hb.data(), hsz);
    if (num_blocks > len / hsz)
        throw ReadError("VTU compressed-block header declares more blocks than it holds");
    const std::size_t num_header_bytes = hsz * (3 + static_cast<std::size_t>(num_blocks));
    const std::size_t num_header_chars = ((num_header_bytes + 2) / 3) * 4;
    if (len < num_header_chars)
        throw ReadError("VTU compressed-block header truncated");
    const std::vector<unsigned char> header = b64decode(pText, num_header_chars);
    if (header.size() < num_header_bytes)
        throw ReadError("VTU compressed-block header truncated");
    const std::uint64_t max_block = read_uint_le(header.data() + hsz, hsz);
    const std::uint64_t last_block = read_uint_le(header.data() + 2 * hsz, hsz);
    std::vector<std::uint64_t> comp_sizes(num_blocks);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        comp_sizes[k] = read_uint_le(header.data() + (3 + k) * hsz, hsz);
    const std::vector<unsigned char> blockdata =
        b64decode(pText + num_header_chars, len - num_header_chars);
    std::vector<std::size_t> in_off(static_cast<std::size_t>(num_blocks) + 1, 0);
    for (std::uint64_t k = 0; k < num_blocks; ++k) {
        const std::size_t at = in_off[static_cast<std::size_t>(k)];
        if (comp_sizes[k] > blockdata.size() - at)
            throw ReadError("VTU compressed blocks are larger than the data that follows");
        in_off[static_cast<std::size_t>(k) + 1] = at + static_cast<std::size_t>(comp_sizes[k]);
    }
    const std::uint64_t ceiling = (std::uint64_t{1} << 16) * (blockdata.size() + 1) + (1u << 20);
    if (num_blocks && (max_block > ceiling || last_block > ceiling ||
                       (num_blocks - 1) > ceiling / std::max<std::uint64_t>(max_block, 1)))
        throw ReadError("VTU compressed-block header declares an implausible size");
    const std::size_t total = num_blocks ? static_cast<std::size_t>(num_blocks - 1) *
                                                   static_cast<std::size_t>(max_block) +
                                               static_cast<std::size_t>(last_block)
                                         : 0;
    if (total > ceiling)
        throw ReadError("VTU compressed-block header declares an implausible size");
    const std::size_t isz = dtype_size(dt);
    const std::size_t n = isz ? total / isz : 0;
    // A payload that is not a whole number of elements (a malformed file) keeps
    // the buffered path, whose last element is simply dropped.
    if (isz == 0 || total % isz != 0) {
        const std::vector<unsigned char> bytes = vtu_decode_blocks(pText, len, hsz, codec);
        NDArray a(dt, {n});
        if (n)
            std::memcpy(a.Data(), bytes.data(), n * isz);
        return a;
    }
    NDArray a = NDArray::Uninit(dt, {n});
    unsigned char* dst = reinterpret_cast<unsigned char*>(a.Data());
    parallel_for(
        static_cast<std::size_t>(num_blocks),
        [&](std::size_t k) {
            const std::size_t expected = (k + 1 == num_blocks)
                                             ? static_cast<std::size_t>(last_block)
                                             : static_cast<std::size_t>(max_block);
            unsigned char* out = dst + k * static_cast<std::size_t>(max_block);
            const std::size_t written =
                vtub_decompress_into(codec, blockdata.data() + in_off[k],
                                     static_cast<std::size_t>(comp_sizes[k]), out, expected);
            // A short block leaves the rest zero, as the zero-filled buffer did.
            if (written < expected)
                std::memset(out + written, 0, expected - written);
        },
        /*grain=*/1);
    return a;
}

}  // namespace detail
}  // namespace meshioplusplus
