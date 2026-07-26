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
 * @file vtu_binary.hpp
 * @brief Base64 and VTU "binary" `DataArray` codecs (raw and zlib/zstd/lz4
 * -compressed), shared helpers behind the VTU (VTK XML) reader/writer's
 * binary I/O.
 *
 * VTU's binary encoding wraps raw little-endian bytes (optionally compressed
 * in fixed-size blocks) as base64 text inside the XML. This header declares
 * both halves: plain base64 encode/decode (`b64encode`/`b64decode`), and the
 * VTU-specific framing on top of it — `vtu_decode_uncompressed`/
 * `vtu_encode_binary(codec=None)` for the uncompressed scheme (a
 * little-endian byte-count header followed by raw data) and
 * `vtu_decode_blocks`/`vtu_encode_binary(codec=...)` for the compressed block
 * scheme (num_blocks / max_block_size / last_block_size header, then each
 * block's compressed size, then the concatenated compressed blocks). All
 * codecs share this one block framing; only the per-block compressor differs
 * (`VtkCodec`). zlib/zstd/lz4 support is conditionally compiled on the
 * matching `MESHIOPLUSPLUS_HAS_*` macro; without it, the codec is still
 * *selectable* but decoding/encoding with it throws rather than compiling out
 * entirely — the absence is discovered at runtime and routes the caller to
 * the Python fallback. Both directions parallelize per independent unit of
 * work (base64 3-byte groups; compressed blocks) via `parallel_for`, since
 * base64/compression are genuinely compute-bound (unlike the
 * memory-bandwidth-bound gather/byteswap work elsewhere, which uses
 * `parallel_for_bw` instead).
 *
 * Every function here is called once per data array or once per ~32 KiB
 * block — both far coarser than per-scalar — so bodies (and the zlib/zstd/lz4
 * `<...>` headers they need) live in `src/cpp/src/detail/vtu_binary.cpp` rather
 * than inline here. Keeping the codec headers out of this header also keeps
 * them out of every translation unit (and amalgamation consumer) that merely
 * needs the declarations below.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief The standard base64 alphabet (RFC 4648), indexed by 6-bit value. */
MESHIOPLUSPLUS_API const char* b64_table();

/**
 * @brief Base64-encodes `len` bytes of `data`.
 *
 * Every full 3-byte input group maps to exactly 4 output characters at a
 * fixed, independently-computable offset, so the output string is
 * pre-sized once and each group is encoded directly into its slot via
 * `parallel_for` — no synchronization or intermediate buffering needed. Any
 * trailing 1- or 2-byte group is handled afterward, sequentially, with the
 * standard `'='` padding.
 * @param pData Bytes to encode.
 * @param len Number of bytes in `pData`.
 * @return The base64-encoded text, `'='`-padded to a multiple of 4 characters.
 */
MESHIOPLUSPLUS_API std::string b64encode(const unsigned char* pData, std::size_t len);

/**
 * @brief Base64-decodes `len` characters of `s`.
 *
 * Builds (and caches, in a function-local `static`) an inverse lookup table
 * from ASCII byte to 6-bit value on first call. Silently skips `'='`
 * padding and whitespace (`\n \r space \t`), and silently ignores any other
 * character outside the base64 alphabet, rather than treating either as an
 * error — VTU-embedded base64 can be split across lines.
 * @param pS Base64 text to decode (need not be NUL-terminated; length is explicit).
 * @param len Number of characters in `pS` to consider.
 * @return The decoded raw bytes.
 */
MESHIOPLUSPLUS_API std::vector<unsigned char> b64decode(const char* pS, std::size_t len);

/**
 * @brief Block-compression codec of a VTK XML "binary" DataArray.
 *
 * All codecs share VTU's block framing verbatim (num_blocks / max_block /
 * last_block / per-block compressed sizes, fixed 32 KiB blocks); only the
 * per-block compressor differs, which is what lets one framing implementation
 * serve them all.
 */
enum class VtkCodec {
    None,  ///< raw bytes behind a byte-count header
    Zlib,  ///< vtkZLibDataCompressor -- the default, and the only one always available
    LZ4,   ///< vtkLZ4DataCompressor -- a real VTK compressor
    ZSTD,  ///< vtkZSTDDataCompressor -- a meshio++ extension; VTK has no ZSTD compressor
    LZMA   ///< vtkLZMADataCompressor -- recognized but not implemented (Python fallback)
};

/** @brief The `compressor=` attribute a codec is recorded under, or "" for None. */
MESHIOPLUSPLUS_API const char* vtk_codec_compressor(VtkCodec codec);

/** @brief Short user-facing codec name (`--codec` values). */
MESHIOPLUSPLUS_API const char* vtk_codec_name(VtkCodec codec);

/**
 * @name Per-codec block dispatch
 *
 * These always exist and throw at runtime when the codec was not compiled in --
 * deliberately: a *link* error would make the absence a build failure,
 * whereas a ReadError/WriteError routes the caller to the Python fallback,
 * which is the documented contract. The message names the CMake option to
 * turn on, in the spirit of `registry_compiled_out()`.
 * @{
 */

/** @brief Whether @p codec can be decoded by this build. */
MESHIOPLUSPLUS_API bool vtk_codec_available(VtkCodec codec);

/** @brief The CMake option that would enable @p codec. */
MESHIOPLUSPLUS_API std::string vtk_codec_build_option(VtkCodec codec);

/** @brief Actionable "this build cannot do that" message. */
MESHIOPLUSPLUS_API std::string vtk_codec_missing_message(VtkCodec codec, bool for_write);

/** @throws ReadError when @p codec cannot be decoded by this build. */
MESHIOPLUSPLUS_API void vtk_codec_require_read(VtkCodec codec);

/** @throws WriteError when @p codec cannot be encoded by this build. */
MESHIOPLUSPLUS_API void vtk_codec_require_write(VtkCodec codec);

/** @brief Compress one block with @p codec. Callers must have required it. */
MESHIOPLUSPLUS_API std::vector<unsigned char> vtk_codec_compress_block(VtkCodec codec, const unsigned char* pSrc,
                                                    std::size_t n);

/** @brief Decompress one block with @p codec into @p expected bytes. */
MESHIOPLUSPLUS_API std::vector<unsigned char> vtk_codec_decompress_block(VtkCodec codec, const unsigned char* pSrc,
                                                      std::size_t n, std::size_t expected);
/** @} */

/**
 * @brief Reads a little-endian unsigned integer of `isz` bytes from `p`.
 * @param pP Buffer to read from, at least `isz` bytes.
 * @param isz Width in bytes of the integer to read (typically 4 or 8, the
 *            VTU header_type item size).
 * @return The decoded value, widened to `uint64_t`.
 */
MESHIOPLUSPLUS_API std::uint64_t read_uint_le(const unsigned char* pP, std::size_t isz);

/**
 * @brief Decodes an uncompressed VTU "binary" `DataArray`: base64 text of a
 * little-endian byte-count header followed by the raw payload.
 * @param pText Base64-encoded DataArray text.
 * @param len Length of `pText` in characters.
 * @param hsz `header_type` item size in bytes (4 for `UInt32`, 8 for `UInt64`).
 * @return The decoded raw payload bytes (header stripped).
 * @throws ReadError if the decoded data is shorter than the header, or
 *         shorter than the header declares.
 */
MESHIOPLUSPLUS_API std::vector<unsigned char> vtu_decode_uncompressed(const char* pText, std::size_t len,
                                                   std::size_t hsz);

/**
 * @brief Decodes a compressed VTU "binary" `DataArray` (the VTK block
 * compression scheme), with @p codec selecting the per-block compressor.
 *
 * The format, all base64-encoded: a header of `num_blocks`, `max_block`,
 * `last_block` (each `hsz` bytes), then `num_blocks` compressed-size
 * entries, then the concatenated compressed blocks themselves (each block
 * `max_block` bytes decompressed, except the last which is `last_block`).
 * Decoded in three passes: decode just enough base64 to learn
 * `num_blocks`, decode the rest of the header to get each block's
 * compressed size, then base64-decode the block data. Input offsets are a
 * cheap sequential prefix sum of the per-block compressed sizes; output
 * offsets are `k * max_block` by construction, so with both known up front
 * the per-block decompress calls are independent and run under
 * `parallel_for` with `grain=1` (each ~32 KiB block is a full unit of
 * inflate/decompress work, so per-block dispatch is exactly right — this is
 * compute-bound work, unlike the memory-gather use of `parallel_for_bw`
 * elsewhere).
 *
 * @param pText Base64-encoded DataArray text.
 * @param len Length of `pText` in characters.
 * @param hsz `header_type` item size in bytes (4 for `UInt32`, 8 for `UInt64`).
 * @param codec block-compression codec recorded in the file.
 * @return The decoded, decompressed raw payload bytes (all blocks concatenated).
 * @throws ReadError if @p codec was not compiled into this build, or if the
 *         header/data is truncated, or if any block fails to decompress.
 */
MESHIOPLUSPLUS_API std::vector<unsigned char> vtu_decode_blocks(const char* pText, std::size_t len, std::size_t hsz,
                                             VtkCodec codec);

/**
 * @brief Encodes raw little-endian bytes as a VTU "binary" `DataArray` text,
 * either uncompressed or compressed with @p codec (block scheme).
 *
 * Uncompressed (`codec == VtkCodec::None`): a 4-byte little-endian length
 * header followed by the raw bytes, base64-encoded as one unit.
 *
 * Compressed: splits `data` into fixed 32 KiB blocks, compresses each
 * independently under `parallel_for` with `grain=1` (each block is a full,
 * sizeable unit of compute — one whole compress call — so per-block dispatch
 * is ideal; this is compute-bound, unlike the memory-gather work that uses
 * `parallel_for_bw`), then emits the
 * `num_blocks`/`max_block`/`last_block_size`/per-block-compressed-size
 * header followed by the concatenated compressed blocks, all base64-encoded.
 *
 * @param pData Raw bytes to encode (already in the file's target byte order).
 * @param nbytes Number of bytes in `pData`.
 * @param codec block-compression codec to use, or `VtkCodec::None` for raw.
 * @return The base64-encoded VTU `DataArray` text.
 * @throws WriteError if @p codec is requested but was not compiled into this
 *         build.
 */
MESHIOPLUSPLUS_API std::string vtu_encode_binary(const unsigned char* pData, std::size_t nbytes, VtkCodec codec);

}  // namespace detail
}  // namespace meshioplusplus
