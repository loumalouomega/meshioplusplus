#pragma once
//
// base64 + VTU binary DataArray encoding/decoding, matching meshio's Python
// implementation so files roundtrip with the reference reader.
//
//   uncompressed: base64( header[UInt32 total_nbytes] + raw_data )
//   zlib:         base64( header[nblocks, blocksize, last, csizes...] )
//                 ++ base64( concat(zlib_block) )   (two separate encodings)
//
// The header dtype is UInt32 (meshio's default header_type). Host is assumed
// little-endian.

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <zlib.h>

#include "meshio/exceptions.hpp"

namespace meshio {
namespace detail {

inline const char* b64_table() {
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
}

inline std::string b64encode(const unsigned char* data, std::size_t len) {
    const char* tbl = b64_table();
    std::string out;
    out.reserve(((len + 2) / 3) * 4);
    std::size_t i = 0;
    for (; i + 2 < len; i += 3) {
        unsigned n = (unsigned(data[i]) << 16) | (unsigned(data[i + 1]) << 8) |
                     unsigned(data[i + 2]);
        out.push_back(tbl[(n >> 18) & 63]);
        out.push_back(tbl[(n >> 12) & 63]);
        out.push_back(tbl[(n >> 6) & 63]);
        out.push_back(tbl[n & 63]);
    }
    if (i < len) {
        const bool two = (i + 1 < len);
        unsigned n = unsigned(data[i]) << 16;
        if (two) n |= unsigned(data[i + 1]) << 8;
        out.push_back(tbl[(n >> 18) & 63]);
        out.push_back(tbl[(n >> 12) & 63]);
        out.push_back(two ? tbl[(n >> 6) & 63] : '=');
        out.push_back('=');
    }
    return out;
}

inline std::vector<unsigned char> b64decode(const char* s, std::size_t len) {
    static int8_t inv[256];
    static bool init = false;
    if (!init) {
        for (int i = 0; i < 256; ++i) inv[i] = -1;
        const char* tbl = b64_table();
        for (int i = 0; i < 64; ++i) inv[(unsigned char)tbl[i]] = static_cast<int8_t>(i);
        init = true;
    }
    std::vector<unsigned char> out;
    out.reserve(len / 4 * 3);
    int buf = 0, bits = 0;
    for (std::size_t i = 0; i < len; ++i) {
        char ch = s[i];
        if (ch == '=' || ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t') continue;
        int v = inv[(unsigned char)ch];
        if (v < 0) continue;
        buf = (buf << 6) | v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<unsigned char>((buf >> bits) & 0xFF));
        }
    }
    return out;
}

inline std::vector<unsigned char> zlib_compress_block(const unsigned char* src,
                                                      std::size_t n) {
    uLongf bound = compressBound(static_cast<uLong>(n));
    std::vector<unsigned char> out(bound);
    uLongf destLen = bound;
    int r = compress(out.data(), &destLen, src, static_cast<uLong>(n));
    if (r != Z_OK) throw WriteError("zlib compression failed");
    out.resize(destLen);
    return out;
}

inline std::vector<unsigned char> zlib_decompress(const unsigned char* src,
                                                  std::size_t n,
                                                  std::size_t expected) {
    std::vector<unsigned char> out(expected);
    uLongf destLen = static_cast<uLongf>(expected);
    int r = uncompress(out.data(), &destLen, src, static_cast<uLong>(n));
    if (r != Z_OK) throw ReadError("zlib decompression failed");
    out.resize(destLen);
    return out;
}

inline std::uint64_t read_uint_le(const unsigned char* p, std::size_t isz) {
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < isz; ++i)
        v |= static_cast<std::uint64_t>(p[i]) << (8 * i);
    return v;
}

// Decode an uncompressed VTU "binary" DataArray (base64 of header + raw data).
// `hsz` is the header_type item size (4 for UInt32, 8 for UInt64).
inline std::vector<unsigned char> vtu_decode_uncompressed(const char* text,
                                                          std::size_t len,
                                                          std::size_t hsz) {
    std::vector<unsigned char> all = b64decode(text, len);
    if (all.size() < hsz) throw ReadError("VTU binary data too short");
    std::uint64_t total = read_uint_le(all.data(), hsz);
    if (all.size() < hsz + total) throw ReadError("VTU binary data truncated");
    return std::vector<unsigned char>(all.begin() + hsz, all.begin() + hsz + total);
}

// Decode a zlib-compressed VTU "binary" DataArray (block scheme).
inline std::vector<unsigned char> vtu_decode_zlib(const char* text, std::size_t len,
                                                  std::size_t hsz) {
    std::size_t first_chars = ((hsz + 2) / 3) * 4;
    if (len < first_chars) throw ReadError("VTU zlib header too short");
    std::vector<unsigned char> hb = b64decode(text, first_chars);
    std::uint64_t num_blocks = read_uint_le(hb.data(), hsz);

    std::size_t num_header_bytes = hsz * (3 + static_cast<std::size_t>(num_blocks));
    std::size_t num_header_chars = ((num_header_bytes + 2) / 3) * 4;
    if (len < num_header_chars) throw ReadError("VTU zlib header truncated");
    std::vector<unsigned char> header = b64decode(text, num_header_chars);

    std::uint64_t max_block = read_uint_le(header.data() + hsz, hsz);
    std::uint64_t last_block = read_uint_le(header.data() + 2 * hsz, hsz);
    std::vector<std::uint64_t> comp_sizes(num_blocks);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        comp_sizes[k] = read_uint_le(header.data() + (3 + k) * hsz, hsz);

    std::vector<unsigned char> blockdata =
        b64decode(text + num_header_chars, len - num_header_chars);

    std::vector<unsigned char> out;
    std::size_t off = 0;
    for (std::uint64_t k = 0; k < num_blocks; ++k) {
        std::size_t expected =
            (k + 1 == num_blocks) ? static_cast<std::size_t>(last_block)
                                  : static_cast<std::size_t>(max_block);
        auto dec = zlib_decompress(blockdata.data() + off,
                                   static_cast<std::size_t>(comp_sizes[k]), expected);
        out.insert(out.end(), dec.begin(), dec.end());
        off += static_cast<std::size_t>(comp_sizes[k]);
    }
    return out;
}

// Encode raw little-endian bytes as a VTU "binary" DataArray text.
inline std::string vtu_encode_binary(const unsigned char* data, std::size_t nbytes,
                                     bool zlib_compress) {
    if (!zlib_compress) {
        std::vector<unsigned char> buf(4 + nbytes);
        std::uint32_t header = static_cast<std::uint32_t>(nbytes);
        std::memcpy(buf.data(), &header, 4);
        if (nbytes) std::memcpy(buf.data() + 4, data, nbytes);
        return b64encode(buf.data(), buf.size());
    }

    const std::uint32_t max_block = 32768;
    std::uint32_t num_blocks =
        static_cast<std::uint32_t>((nbytes + max_block - 1) / max_block);
    std::uint32_t last_block_size =
        num_blocks ? static_cast<std::uint32_t>(nbytes - std::size_t(num_blocks - 1) * max_block)
                   : max_block;

    std::vector<std::vector<unsigned char>> blocks;
    blocks.reserve(num_blocks);
    for (std::uint32_t b = 0; b < num_blocks; ++b) {
        std::size_t off = std::size_t(b) * max_block;
        std::size_t len = std::min<std::size_t>(max_block, nbytes - off);
        blocks.push_back(zlib_compress_block(data + off, len));
    }

    std::vector<std::uint32_t> header;
    header.reserve(3 + num_blocks);
    header.push_back(num_blocks);
    header.push_back(max_block);
    header.push_back(last_block_size);
    for (const auto& b : blocks) header.push_back(static_cast<std::uint32_t>(b.size()));

    std::string out = b64encode(reinterpret_cast<const unsigned char*>(header.data()),
                                header.size() * sizeof(std::uint32_t));
    std::vector<unsigned char> concat;
    for (const auto& b : blocks) concat.insert(concat.end(), b.begin(), b.end());
    out += b64encode(concat.data(), concat.size());
    return out;
}

}  // namespace detail
}  // namespace meshio
