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
 * @file detail/crs.hpp
 * @brief Georeferencing coordinate reference system metadata (`field_data["meshio:crs"]`,
 * roadmap 9.10.2). Stores a WKT2 or `EPSG:<code>` string as UInt8 bytes.
 */

// System includes
#include <cstdint>
#include <cstring>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/// The `field_data` key storing a mesh's CRS as UTF-8 UInt8 bytes (e.g. `"EPSG:4326"` or WKT2).
inline constexpr const char* kCrsName = "meshio:crs";

/// Return the CRS string stored on @p rMesh, or empty string when absent.
inline std::string mesh_crs(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kCrsName))
        return {};
    const NDArray& a = rMesh.FieldData(kCrsName);
    if (a.Dtype() != DType::UInt8 || a.Size() == 0)
        return {};
    const auto* p = a.As<std::uint8_t>();
    std::size_t n = a.Size();
    while (n > 0 && p[n - 1] == 0)
        --n;
    return std::string(reinterpret_cast<const char*>(p), n);
}

/// Record @p rCrs (WKT2 or `EPSG:<code>`) on @p rMesh under `field_data["meshio:crs"]`.
inline void set_mesh_crs(Mesh& rMesh, const std::string& rCrs) {
    NDArray a = NDArray::Uninit(DType::UInt8, {rCrs.size()});
    if (!rCrs.empty())
        std::memcpy(a.As<std::uint8_t>(), rCrs.data(), rCrs.size());
    rMesh.AddFieldData(kCrsName, std::move(a));
}

}  // namespace detail
}  // namespace meshioplusplus
