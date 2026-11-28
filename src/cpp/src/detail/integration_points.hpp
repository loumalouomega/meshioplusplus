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
 * @file detail/integration_points.hpp
 * @brief Unified integration-point (Gauss-point and element-nodal) metadata
 *        convention in `field_data` (roadmap 8.3).
 *
 * A per-integration-point cell array is stored as a 2-D `cell_data` array of
 * shape `(num_cells, points * width)` in point-major order (all `width`
 * components of integration point `0`, then point `1`, ...). Its layout is
 * recorded in:
 *
 *   `field_data["meshio:ip:<name>"] = Int64[3] { points, width, rule_id }`
 *
 * where `rule_id` identifies the point locations in `[0, 1]` reference
 * coordinates:
 *   - `-1`: unknown / format-defined order (uniform averaging for `ip_to_cells`,
 *           mean extrapolation for `ip_to_points`).
 *   - `0`: element-nodal data (one sample per cell node in connectivity order).
 *   - `>= 1`: exactness order passed to `detail::quadrature_rule(type, rule_id)`.
 *
 * Custom reference coordinates and weights may also be attached per cell type:
 *   `field_data["meshio:ip_coords:<name>:<cell_type>"] = Float64[points, dim]`
 *   `field_data["meshio:ip_weights:<name>:<cell_type>"] = Float64[points]`
 *
 * For backward compatibility, `get_ip_layout` also recognises the legacy
 * format-specific keys `abaqus:layout:<name>`, `lsdyna_d3plot:layout:<name>`,
 * `gid:gauss_points:<name>`, and `gmsh:element_node_components:<name>`.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "quadrature.hpp"

namespace meshioplusplus::detail {

inline constexpr const char* kIpLayoutPrefix = "meshio:ip:";
inline constexpr const char* kIpCoordsPrefix = "meshio:ip_coords:";
inline constexpr const char* kIpWeightsPrefix = "meshio:ip_weights:";

/// Parsed integration-point layout for a cell_data array.
struct IpLayout {
    std::int64_t mPoints = 0;
    std::int64_t mWidth = 0;
    std::int64_t mRuleId = -1;
    bool mValid = false;
};

/// Build the `meshio:ip:<name>` field_data key.
inline std::string ip_layout_key(const std::string& rName) {
    return std::string(kIpLayoutPrefix) + rName;
}

/// Build the `meshio:ip_coords:<name>:<cell_type>` field_data key.
inline std::string ip_coords_key(const std::string& rName, const std::string& rCellType) {
    return std::string(kIpCoordsPrefix) + rName + ":" + rCellType;
}

/// Build the `meshio:ip_weights:<name>:<cell_type>` field_data key.
inline std::string ip_weights_key(const std::string& rName, const std::string& rCellType) {
    return std::string(kIpWeightsPrefix) + rName + ":" + rCellType;
}

/// Attach or overwrite `field_data["meshio:ip:<name>"] = [points, width, rule_id]`.
inline void set_ip_layout(Mesh& rMesh, const std::string& rName, std::int64_t Points,
                          std::int64_t Width, std::int64_t RuleId = -1) {
    NDArray arr(DType::Int64, {3});
    auto* p = arr.As<std::int64_t>();
    p[0] = Points;
    p[1] = Width;
    p[2] = RuleId;
    rMesh.AddFieldData(ip_layout_key(rName), std::move(arr));
}

/// Infer the smallest `quadrature_rule` exactness order (1..6) whose point count
/// equals `NumPoints` for `Type`, or `-1` if none matches.
inline std::int64_t infer_ip_rule_id(CellType Type, std::int64_t NumPoints) {
    if (NumPoints < 1)
        return -1;
    for (int order = 1; order <= 6; ++order) {
        const QuadratureRule rule = quadrature_rule(Type, order);
        if (static_cast<std::int64_t>(rule.mPoints.size()) == NumPoints)
            return static_cast<std::int64_t>(order);
    }
    return -1;
}

/// Read integration-point layout for `cell_data[rName]`, checking `meshio:ip:`
/// first and falling back to legacy per-format layout keys.
inline IpLayout get_ip_layout(const Mesh& rMesh, const std::string& rName) {
    IpLayout out;
    const std::string unified = ip_layout_key(rName);
    if (rMesh.HasFieldData(unified)) {
        const NDArray& a = rMesh.FieldData(unified);
        if (a.Size() >= 2) {
            out.mPoints = read_int(a, 0);
            out.mWidth = read_int(a, 1);
            out.mRuleId = a.Size() >= 3 ? read_int(a, 2) : -1;
            out.mValid = (out.mPoints >= 1 && out.mWidth >= 1);
            return out;
        }
    }

    // Legacy 2-element `[points, width]` keys.
    for (const char* pPrefix : {"abaqus:layout:", "lsdyna_d3plot:layout:"}) {
        const std::string key = std::string(pPrefix) + rName;
        if (rMesh.HasFieldData(key)) {
            const NDArray& a = rMesh.FieldData(key);
            if (a.Size() >= 2) {
                out.mPoints = read_int(a, 0);
                out.mWidth = read_int(a, 1);
                out.mRuleId = -1;
                out.mValid = (out.mPoints >= 1 && out.mWidth >= 1);
                return out;
            }
        }
    }

    // Helper: inspect first non-empty block column count for single-number legacy keys.
    auto first_cols = [&]() -> std::size_t {
        if (!rMesh.HasCellData(rName))
            return 0;
        const std::size_t nb = rMesh.CellDataNumBlocks(rName);
        for (std::size_t b = 0; b < nb; ++b) {
            const NDArray& arr = rMesh.CellData(rName, b);
            if (rows(arr) > 0)
                return cols(arr);
        }
        return 0;
    };

    const std::string gid_key = "gid:gauss_points:" + rName;
    if (rMesh.HasFieldData(gid_key)) {
        const NDArray& a = rMesh.FieldData(gid_key);
        if (a.Size() >= 1) {
            const std::int64_t g = static_cast<std::int64_t>(read_double(a, 0));
            const std::size_t c = first_cols();
            if (g >= 1 && c > 0 && c % static_cast<std::size_t>(g) == 0) {
                out.mPoints = g;
                out.mWidth = static_cast<std::int64_t>(c / static_cast<std::size_t>(g));
                out.mRuleId = -1;
                out.mValid = true;
                return out;
            }
        }
    }

    const std::string gmsh_key = "gmsh:element_node_components:" + rName;
    if (rMesh.HasFieldData(gmsh_key)) {
        const NDArray& a = rMesh.FieldData(gmsh_key);
        if (a.Size() >= 1) {
            const std::int64_t w = read_int(a, 0);
            const std::size_t c = first_cols();
            if (w >= 1 && c > 0 && c % static_cast<std::size_t>(w) == 0) {
                out.mPoints = static_cast<std::int64_t>(c / static_cast<std::size_t>(w));
                out.mWidth = w;
                out.mRuleId = 0;
                out.mValid = true;
                return out;
            }
        }
    }

    return out;
}

/// Resolve the quadrature rule (`mPoints` in `[0, 1]` reference coordinates and
/// `mWeights`) for `(rName, Type)` given `rLayout`. Returns true when a rule
/// matching `rLayout.mPoints` is found.
inline bool resolve_ip_quadrature(const Mesh& rMesh, const std::string& rName, CellType Type,
                                  const IpLayout& rLayout, QuadratureRule& rRule) {
    rRule = QuadratureRule{};
    if (!rLayout.mValid || rLayout.mPoints < 1)
        return false;

    const std::string ct_name = cell_type_name(Type);
    const std::string ckey = ip_coords_key(rName, ct_name);
    if (rMesh.HasFieldData(ckey)) {
        const NDArray& coords = rMesh.FieldData(ckey);
        const std::size_t npts = rows(coords);
        const std::size_t dim = cols(coords);
        if (static_cast<std::int64_t>(npts) == rLayout.mPoints && dim >= 1 && dim <= 3) {
            rRule.mPoints.resize(npts, {0.0, 0.0, 0.0});
            for (std::size_t i = 0; i < npts; ++i) {
                for (std::size_t d = 0; d < dim; ++d)
                    rRule.mPoints[i][d] = read_double(coords, i * dim + d);
            }
            const std::string wkey = ip_weights_key(rName, ct_name);
            if (rMesh.HasFieldData(wkey) && rMesh.FieldData(wkey).Size() == npts) {
                const NDArray& w = rMesh.FieldData(wkey);
                rRule.mWeights.resize(npts);
                for (std::size_t i = 0; i < npts; ++i)
                    rRule.mWeights[i] = read_double(w, i);
            } else {
                rRule.mWeights.assign(npts, 1.0 / static_cast<double>(npts));
            }
            return true;
        }
    }

    if (rLayout.mRuleId >= 1) {
        QuadratureRule q = quadrature_rule(Type, static_cast<int>(rLayout.mRuleId));
        if (static_cast<std::int64_t>(q.mPoints.size()) == rLayout.mPoints) {
            rRule = std::move(q);
            return true;
        }
        // Fallback: check if another exactness order matches `rLayout.mPoints`.
        const std::int64_t inferred = infer_ip_rule_id(Type, rLayout.mPoints);
        if (inferred >= 1) {
            rRule = quadrature_rule(Type, static_cast<int>(inferred));
            return true;
        }
    }
    return false;
}

}  // namespace meshioplusplus::detail
