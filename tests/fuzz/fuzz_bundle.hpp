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

/**
 * @file fuzz_bundle.hpp
 * @brief Deterministic multi-file bundle for companion-file fuzzing (doc/fuzzing.md).
 *
 * Single-file libFuzzer inputs cannot reach readers that need companions: XDMF
 * XML/heavy-data bundles (`.xdmf` plus `.h5`/`.bin` payloads) and ADIOS2 `.bp`
 * directories (`md.idx`, `md.0`, `data.0`, ...). This header packs such a set
 * into one fuzz input and unpacks it back, so findings stay ordinary file bytes
 * that replay with the unmodified replay harness.
 *
 * Format (all integers little-endian): magic `MIOB`, version byte `1`, u16
 * entry count, then per entry u16 name length, name bytes, u32 data length,
 * data bytes. Names are relative paths (`input.xdmf`, `data.h5`,
 * `input.bp/md.idx`): no absolute paths, no `..`, no empty components, at most
 * 16 entries, each name at most 256 bytes. Malformed bundles throw `ReadError`,
 * the harness's tolerated refusal, never an abort: the mutator freely explores
 * bundle structure.
 */

#pragma once

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

inline constexpr std::size_t kFuzzBundleMaxEntries = 16;
inline constexpr std::size_t kFuzzBundleMaxName = 256;

/** @brief Whether `pData[0:size]` starts with the bundle magic and version. */
inline bool fuzz_bundle_is_bundle(const std::uint8_t* pData, std::size_t size) {
    return size >= 7 && pData[0] == 'M' && pData[1] == 'I' && pData[2] == 'O' && pData[3] == 'B' &&
           pData[4] == 1;
}

/** @brief One bundle entry: relative path plus raw bytes. */
using FuzzBundleEntry = std::pair<std::string, std::vector<std::uint8_t>>;

inline std::uint16_t fuzz_bundle_read_u16(const std::uint8_t* p) {
    return static_cast<std::uint16_t>(p[0] | (static_cast<std::uint16_t>(p[1]) << 8));
}

inline std::uint32_t fuzz_bundle_read_u32(const std::uint8_t* p) {
    return static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
           (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
}

inline void fuzz_bundle_write_u16(std::vector<std::uint8_t>& rOut, std::uint16_t v) {
    rOut.push_back(static_cast<std::uint8_t>(v & 0xFF));
    rOut.push_back(static_cast<std::uint8_t>((v >> 8) & 0xFF));
}

inline void fuzz_bundle_write_u32(std::vector<std::uint8_t>& rOut, std::uint32_t v) {
    rOut.push_back(static_cast<std::uint8_t>(v & 0xFF));
    rOut.push_back(static_cast<std::uint8_t>((v >> 8) & 0xFF));
    rOut.push_back(static_cast<std::uint8_t>((v >> 16) & 0xFF));
    rOut.push_back(static_cast<std::uint8_t>((v >> 24) & 0xFF));
}

/**
 * @brief Whether a bundle entry name is safe to materialize under a scratch dir.
 * @param rName Relative path to test.
 * @return True for `a.b`, `d/f`, never for absolute paths, `..` or empty parts.
 */
inline bool fuzz_bundle_name_ok(const std::string& rName) {
    if (rName.empty() || rName.size() > kFuzzBundleMaxName || rName[0] == '/' || rName[0] == '\\')
        return false;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= rName.size(); ++i) {
        if (i == rName.size() || rName[i] == '/' || rName[i] == '\\') {
            const std::size_t len = i - start;
            if (len == 0 || len > 128)
                return false;
            if (len == 1 && rName[start] == '.')
                return false;
            if (len == 2 && rName[start] == '.' && rName[start + 1] == '.')
                return false;
            for (std::size_t j = start; j < i; ++j) {
                const unsigned char c = static_cast<unsigned char>(rName[j]);
                const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                                (c >= '0' && c <= '9') || c == '.' || c == '_' || c == '-' ||
                                c == '+';
                if (!ok)
                    return false;
            }
            start = i + 1;
        }
    }
    return true;
}

/**
 * @brief Decode bundle bytes into entries, throwing `ReadError` when malformed.
 * @param pData Raw fuzz input, which must start with the bundle magic.
 * @param size Input length in bytes.
 * @return Entries in bundle order.
 * @throws ReadError on truncated counts, too many entries or unsafe names.
 */
inline std::vector<FuzzBundleEntry> fuzz_bundle_decode(const std::uint8_t* pData,
                                                       std::size_t size) {
    if (!fuzz_bundle_is_bundle(pData, size))
        throw ReadError("fuzz bundle: missing magic");
    if (size < 7)
        throw ReadError("fuzz bundle: truncated header");
    const std::size_t count = fuzz_bundle_read_u16(pData + 5);
    if (count == 0 || count > kFuzzBundleMaxEntries)
        throw ReadError("fuzz bundle: bad entry count");
    std::vector<FuzzBundleEntry> out;
    std::size_t pos = 7;
    for (std::size_t e = 0; e < count; ++e) {
        if (pos + 2 > size)
            throw ReadError("fuzz bundle: truncated name length");
        const std::size_t name_len = fuzz_bundle_read_u16(pData + pos);
        pos += 2;
        if (name_len == 0 || name_len > kFuzzBundleMaxName || pos + name_len + 4 > size)
            throw ReadError("fuzz bundle: bad name length");
        const std::string name(reinterpret_cast<const char*>(pData + pos), name_len);
        pos += name_len;
        if (!fuzz_bundle_name_ok(name))
            throw ReadError("fuzz bundle: unsafe entry name");
        const std::size_t data_len = fuzz_bundle_read_u32(pData + pos);
        pos += 4;
        if (pos + data_len > size)
            throw ReadError("fuzz bundle: truncated entry data");
        std::vector<std::uint8_t> data(pData + pos, pData + pos + data_len);
        pos += data_len;
        out.emplace_back(name, std::move(data));
    }
    if (out.empty())
        throw ReadError("fuzz bundle: no entries");
    return out;
}

/**
 * @brief Pack entries into bundle bytes for seeds and tests.
 * @param rEntries Relative-path entries to pack (1..16, validated names).
 * @return Bundle bytes ready to feed the harness.
 * @throws ReadError on empty input or unsafe names.
 */
inline std::vector<std::uint8_t> fuzz_bundle_encode(const std::vector<FuzzBundleEntry>& rEntries) {
    if (rEntries.empty() || rEntries.size() > kFuzzBundleMaxEntries)
        throw ReadError("fuzz bundle: bad entry count");
    std::vector<std::uint8_t> out = {'M', 'I', 'O', 'B', 1};
    fuzz_bundle_write_u16(out, static_cast<std::uint16_t>(rEntries.size()));
    for (const auto& [name, data] : rEntries) {
        if (!fuzz_bundle_name_ok(name))
            throw ReadError("fuzz bundle: unsafe entry name");
        if (data.size() > 0xFFFFFFFFu)
            throw ReadError("fuzz bundle: entry too large");
        fuzz_bundle_write_u16(out, static_cast<std::uint16_t>(name.size()));
        out.insert(out.end(), name.begin(), name.end());
        fuzz_bundle_write_u32(out, static_cast<std::uint32_t>(data.size()));
        out.insert(out.end(), data.begin(), data.end());
    }
    return out;
}

}  // namespace meshioplusplus
