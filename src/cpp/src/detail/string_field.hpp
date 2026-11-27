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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file detail/string_field.hpp
 * @brief UTF-8 string and boolean-mask conventions in `field_data` (roadmap §8.2).
 *
 * - Scalar string `name`:
 *     - `field_data[name]`           = `UInt8` 1-D array of UTF-8 bytes
 *     - `field_data[name:offsets]`   = `Int64` 1-D array of shape `{1}` holding `{len}`
 * - List of `k` strings `name`:
 *     - `field_data[name]`           = `UInt8` 1-D array of concatenated UTF-8 bytes
 *     - `field_data[name:offsets]`   = `Int64` 1-D array of shape `{k + 1}` holding CSR offsets
 * - Boolean array `name` at location `loc` (`"point"`, `"cell"`, `"field"`):
 *     - stored as `UInt8` (`0`/`1`), marked by `field_data["meshio:bool:<loc>:<name>"] = UInt8[1] { 1 }`
 *
 * Decoding reads byte/bool payloads via `DType::UInt8` fast path or `read_int`
 * fallback so it works identically on MESHIO, NATIVE, and KRATOS backends
 * (where `NativeMesh::canonicalize_array` widens `UInt8` to `Int64`).
 *
 * Core-private: no installed header includes this.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus::detail {

inline constexpr std::string_view kStringOffsetsSuffix = ":offsets";
inline constexpr std::string_view kBoolMarkerPrefix = "meshio:bool:";

inline std::string string_offsets_key(std::string_view rName) {
    return std::string(rName) + std::string(kStringOffsetsSuffix);
}

inline bool is_bool_marker_key(std::string_view rName) {
    return rName.starts_with(kBoolMarkerPrefix);
}

inline std::string bool_marker_key(std::string_view rLoc, std::string_view rName) {
    return std::string(kBoolMarkerPrefix) + std::string(rLoc) + ":" + std::string(rName);
}

inline void set_bool_marker(Mesh& rMesh, std::string_view rLoc, std::string_view rName) {
    NDArray flag(DType::UInt8, {1});
    flag.As<std::uint8_t>()[0] = 1;
    rMesh.AddFieldData(bool_marker_key(rLoc, rName), std::move(flag));
}

inline bool is_bool_array(const Mesh& rMesh, std::string_view rLoc, std::string_view rName) {
    const std::string key = bool_marker_key(rLoc, rName);
    if (!rMesh.HasFieldData(key))
        return false;
    const NDArray& flag = rMesh.FieldData(key);
    return flag.Size() >= 1 && read_int(flag, 0) != 0;
}

/// Returns true if `rKey` is a companion `<base>:offsets` array for an existing
/// string field `<base>` on `rMesh`.
inline bool is_string_offsets_key(const Mesh& rMesh, std::string_view rKey) {
    if (rKey.size() <= kStringOffsetsSuffix.size() || !rKey.ends_with(kStringOffsetsSuffix))
        return false;
    const std::string base(rKey.substr(0, rKey.size() - kStringOffsetsSuffix.size()));
    return rMesh.HasFieldData(base);
}

/// Returns true if `rKey` is internal string/bool metadata (`<base>:offsets` or
/// `meshio:bool:*`) that numeric-only writers should skip.
inline bool is_string_or_bool_meta_key(const Mesh& rMesh, std::string_view rKey) {
    return is_bool_marker_key(rKey) || is_string_offsets_key(rMesh, rKey);
}

/// Store a scalar UTF-8 string in `field_data[rName]`.
inline void set_field_string(Mesh& rMesh, const std::string& rName, std::string_view rValue) {
    NDArray bytes(DType::UInt8, {rValue.size()});
    if (!rValue.empty())
        std::memcpy(bytes.Data(), rValue.data(), rValue.size());
    NDArray offs(DType::Int64, {1});
    offs.As<std::int64_t>()[0] = static_cast<std::int64_t>(rValue.size());
    rMesh.AddFieldData(rName, std::move(bytes));
    rMesh.AddFieldData(string_offsets_key(rName), std::move(offs));
}

/// Store a list of UTF-8 strings in `field_data[rName]`.
inline void set_field_strings(Mesh& rMesh, const std::string& rName,
                              const std::vector<std::string>& rValues) {
    std::size_t total = 0;
    for (const auto& s : rValues)
        total += s.size();
    NDArray bytes(DType::UInt8, {total});
    NDArray offs(DType::Int64, {rValues.size() + 1});
    auto* p_off = offs.As<std::int64_t>();
    auto* p_bytes = bytes.As<std::uint8_t>();
    std::size_t pos = 0;
    p_off[0] = 0;
    for (std::size_t i = 0; i < rValues.size(); ++i) {
        if (!rValues[i].empty())
            std::memcpy(p_bytes + pos, rValues[i].data(), rValues[i].size());
        pos += rValues[i].size();
        p_off[i + 1] = static_cast<std::int64_t>(pos);
    }
    rMesh.AddFieldData(rName, std::move(bytes));
    rMesh.AddFieldData(string_offsets_key(rName), std::move(offs));
}

/// Decode `rBytes[rBegin .. rEnd)` into a `std::string` regardless of whether
/// the backend kept `UInt8` or widened to `Int64`.
inline std::string decode_utf8_slice(const NDArray& rBytes, std::size_t Begin, std::size_t End) {
    if (End <= Begin || End > rBytes.Size())
        return {};
    const std::size_t len = End - Begin;
    std::string out(len, '\0');
    if (rBytes.Dtype() == DType::UInt8) {
        std::memcpy(out.data(), rBytes.As<std::uint8_t>() + Begin, len);
    } else {
        for (std::size_t i = 0; i < len; ++i)
            out[i] = static_cast<char>(static_cast<std::uint8_t>(read_int(rBytes, Begin + i)));
    }
    return out;
}

/// Returns true if `field_data[rName]` is encoded as a scalar UTF-8 string.
inline bool is_field_string(const Mesh& rMesh, const std::string& rName) {
    const std::string off_key = string_offsets_key(rName);
    if (!rMesh.HasFieldData(rName) || !rMesh.HasFieldData(off_key))
        return false;
    const NDArray& offs = rMesh.FieldData(off_key);
    return offs.Shape().size() == 1 && offs.Shape()[0] == 1;
}

/// Returns true if `field_data[rName]` is encoded as a list of UTF-8 strings.
inline bool is_field_string_list(const Mesh& rMesh, const std::string& rName) {
    const std::string off_key = string_offsets_key(rName);
    if (!rMesh.HasFieldData(rName) || !rMesh.HasFieldData(off_key))
        return false;
    const NDArray& offs = rMesh.FieldData(off_key);
    return offs.Shape().size() == 1 && offs.Shape()[0] >= 2;
}

/// Read a scalar UTF-8 string from `field_data[rName]`, or `std::nullopt`.
inline std::optional<std::string> get_field_string(const Mesh& rMesh, const std::string& rName) {
    if (!is_field_string(rMesh, rName))
        return std::nullopt;
    const NDArray& bytes = rMesh.FieldData(rName);
    const NDArray& offs = rMesh.FieldData(string_offsets_key(rName));
    const std::int64_t len = read_int(offs, 0);
    if (len < 0 || static_cast<std::size_t>(len) != bytes.Size())
        return std::nullopt;
    return decode_utf8_slice(bytes, 0, bytes.Size());
}

/// Read a list of UTF-8 strings from `field_data[rName]` (also accepts a scalar
/// string as a 1-element list). Returns `std::nullopt` if not a string field.
inline std::optional<std::vector<std::string>> get_field_strings(const Mesh& rMesh,
                                                                 const std::string& rName) {
    if (auto scalar = get_field_string(rMesh, rName))
        return std::vector<std::string>{std::move(*scalar)};
    if (!is_field_string_list(rMesh, rName))
        return std::nullopt;
    const NDArray& bytes = rMesh.FieldData(rName);
    const NDArray& offs = rMesh.FieldData(string_offsets_key(rName));
    const std::size_t n = offs.Size() - 1;
    std::vector<std::string> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::int64_t a = read_int(offs, i);
        const std::int64_t b = read_int(offs, i + 1);
        if (a < 0 || b < a || static_cast<std::size_t>(b) > bytes.Size())
            return std::nullopt;
        out.push_back(decode_utf8_slice(bytes, static_cast<std::size_t>(a),
                                        static_cast<std::size_t>(b)));
    }
    return out;
}

}  // namespace meshioplusplus::detail
