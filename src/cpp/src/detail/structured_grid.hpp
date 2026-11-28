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
 * @file detail/structured_grid.hpp
 * @brief Structured-grid annotation convention (`meshio:structured:<block>`).
 *
 * Structured readers annotate each expanded structured block with
 * `field_data["meshio:structured:<block>"] = Int64[3] { ni, nj, nk }` (point
 * counts, `i` fastest). Structured writers (`vti`, `vtr`, `vts`) check this
 * annotation before falling back to geometric lattice detection, and geometry-
 * or topology-modifying operations strip the annotation so a deformed or
 * unstructured mesh is never silently miswritten as structured.
 *
 * Core-private: no installed header includes this.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "typed_view.hpp"

namespace meshioplusplus::detail {

inline constexpr std::string_view kStructuredKeyPrefix = "meshio:structured:";

inline std::string structured_key(std::size_t BlockIdx) {
    return std::string(kStructuredKeyPrefix) + std::to_string(BlockIdx);
}

inline bool is_structured_key(std::string_view rName) {
    return rName.starts_with(kStructuredKeyPrefix);
}

/// Annotate block `BlockIdx` with structured point counts `[Ni, Nj, Nk]`.
inline void set_structured_dims(Mesh& rMesh, std::size_t BlockIdx, std::int64_t Ni, std::int64_t Nj,
                                std::int64_t Nk) {
    if (Ni <= 0 || Nj <= 0 || Nk <= 0)
        return;
    NDArray arr(DType::Int64, {3});
    auto* p = arr.As<std::int64_t>();
    p[0] = Ni;
    p[1] = Nj;
    p[2] = Nk;
    rMesh.AddFieldData(structured_key(BlockIdx), std::move(arr));
}

/// Read structured point counts `[ni, nj, nk]` for block `BlockIdx`, if present.
/// Falls back to `plot3d:dims` row `BlockIdx` when `meshio:structured:<block>` is absent.
inline std::optional<std::array<std::int64_t, 3>> get_structured_dims(const Mesh& rMesh,
                                                                      std::size_t BlockIdx) {
    const std::string key = structured_key(BlockIdx);
    if (rMesh.HasFieldData(key)) {
        const NDArray& arr = rMesh.FieldData(key);
        if (arr.Size() == 3) {
            const std::int64_t ni = read_int(arr, 0);
            const std::int64_t nj = read_int(arr, 1);
            const std::int64_t nk = read_int(arr, 2);
            if (ni > 0 && nj > 0 && nk > 0)
                return std::array<std::int64_t, 3>{ni, nj, nk};
        }
    }
    if (rMesh.HasFieldData("plot3d:dims")) {
        const NDArray& arr = rMesh.FieldData("plot3d:dims");
        if (arr.Shape().size() == 2 && arr.Shape()[1] == 3 && BlockIdx < arr.Shape()[0]) {
            const std::int64_t ni = read_int(arr, BlockIdx * 3 + 0);
            const std::int64_t nj = read_int(arr, BlockIdx * 3 + 1);
            const std::int64_t nk = read_int(arr, BlockIdx * 3 + 2);
            if (ni > 0 && nj > 0 && nk > 0)
                return std::array<std::int64_t, 3>{ni, nj, nk};
        }
    }
    return std::nullopt;
}

/// Returns true if `rMesh` has any `meshio:structured:*` field_data entries.
inline bool has_structured_annotations(const Mesh& rMesh) {
    for (const std::string& name : rMesh.FieldDataNames())
        if (is_structured_key(name))
            return true;
    return false;
}

/// Strips all `meshio:structured:*` annotations from `rMesh`, returning a mesh
/// without them (no-op when none are present).
inline void strip_structured_annotations(Mesh& rMesh) {
    std::vector<std::string> keys;
    for (const std::string& name : rMesh.FieldDataNames())
        if (is_structured_key(name))
            keys.push_back(name);
    if (!keys.empty())
        rMesh = data_drop(rMesh, DataLocation::Field, keys);
}

/// Validates that a single-block mesh matches the point count, cell type and
/// cell count of a structured `[ni, nj, nk]` point grid (`i` fastest).
inline bool validate_structured_single_block(const Mesh& rMesh,
                                             const std::array<std::int64_t, 3>& rDims) {
    const std::int64_t ni = rDims[0];
    const std::int64_t nj = rDims[1];
    const std::int64_t nk = rDims[2];
    if (ni < 2 || nj < 2 || nk < 2)
        return false;
    if (rMesh.NumCellBlocks() != 1)
        return false;
    const std::size_t expected_pts =
        static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj) * static_cast<std::size_t>(nk);
    if (rMesh.NumPoints() != expected_pts)
        return false;
    const auto cb = rMesh.Cells(0);
    if (cb.IsRagged() || cb.Type() != "hexahedron")
        return false;
    const std::size_t nx = static_cast<std::size_t>(ni - 1);
    const std::size_t ny = static_cast<std::size_t>(nj - 1);
    const std::size_t nz = static_cast<std::size_t>(nk - 1);
    if (cb.NumCells() != nx * ny * nz)
        return false;

    // Spot-check first and last cell connectivity (i fastest, VTK hex order).
    const Int64View conn(cb.Conn());
    const auto idx = [&](std::size_t i, std::size_t j, std::size_t k) -> std::int64_t {
        return static_cast<std::int64_t>(
            (k * static_cast<std::size_t>(nj) + j) * static_cast<std::size_t>(ni) + i);
    };
    auto check_cell = [&](std::size_t i, std::size_t j, std::size_t k) -> bool {
        const std::size_t c = (k * ny + j) * nx + i;
        const std::int64_t want[8] = {
            idx(i, j, k),     idx(i + 1, j, k),     idx(i + 1, j + 1, k),     idx(i, j + 1, k),
            idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1), idx(i, j + 1, k + 1),
        };
        for (std::size_t m = 0; m < 8; ++m)
            if (conn[c * 8 + m] != want[m])
                return false;
        return true;
    };
    return check_cell(0, 0, 0) && check_cell(nx - 1, ny - 1, nz - 1);
}

/// Recovers a rectilinear (possibly graded, axis-aligned) 3-D structured grid
/// from `rMesh`. Uses `meshio:structured:0` when present, falling back to
/// `lattice_from_mesh` for unannotated uniform lattices.
inline bool rectilinear_from_mesh(const Mesh& rMesh, std::array<std::int64_t, 3>& rDims,
                                  std::vector<double>& rX, std::vector<double>& rY,
                                  std::vector<double>& rZ) {
    if (rMesh.PointDim() < 3)
        return false;
    auto annotated = get_structured_dims(rMesh, 0);
    std::array<std::int64_t, 3> dims{};
    if (annotated) {
        dims = *annotated;
        if (!validate_structured_single_block(rMesh, dims))
            return false;
    } else {
        LatticeSpec spec;
        if (!lattice_from_mesh(rMesh, spec))
            return false;
        dims = {spec.mDims[0] + 1, spec.mDims[1] + 1, spec.mDims[2] + 1};
    }

    const std::size_t ni = static_cast<std::size_t>(dims[0]);
    const std::size_t nj = static_cast<std::size_t>(dims[1]);
    const std::size_t nk = static_cast<std::size_t>(dims[2]);
    const std::size_t pdim = rMesh.PointDim();
    const DoubleView pts(rMesh.Points());

    std::vector<double> xs(ni), ys(nj), zs(nk);
    for (std::size_t i = 0; i < ni; ++i)
        xs[i] = pts[i * pdim + 0];
    for (std::size_t j = 0; j < nj; ++j)
        ys[j] = pts[(j * ni) * pdim + 1];
    for (std::size_t k = 0; k < nk; ++k)
        zs[k] = pts[(k * nj * ni) * pdim + 2];

    auto strictly_increasing = [](const std::vector<double>& rV) {
        for (std::size_t i = 1; i < rV.size(); ++i)
            if (!(rV[i] > rV[i - 1]))
                return false;
        return true;
    };
    if (!strictly_increasing(xs) || !strictly_increasing(ys) || !strictly_increasing(zs))
        return false;

    const double tol_x = 1e-9 * std::max(1.0, std::abs(xs.back() - xs.front()));
    const double tol_y = 1e-9 * std::max(1.0, std::abs(ys.back() - ys.front()));
    const double tol_z = 1e-9 * std::max(1.0, std::abs(zs.back() - zs.front()));
    for (std::size_t k = 0; k < nk; ++k) {
        for (std::size_t j = 0; j < nj; ++j) {
            for (std::size_t i = 0; i < ni; ++i) {
                const std::size_t p = (k * nj + j) * ni + i;
                if (std::abs(pts[p * pdim + 0] - xs[i]) > tol_x ||
                    std::abs(pts[p * pdim + 1] - ys[j]) > tol_y ||
                    std::abs(pts[p * pdim + 2] - zs[k]) > tol_z)
                    return false;
            }
        }
    }

    rDims = dims;
    rX = std::move(xs);
    rY = std::move(ys);
    rZ = std::move(zs);
    return true;
}

/// Recovers a curvilinear 3-D structured grid from `rMesh`. Accepts any point
/// positions when `meshio:structured:0` is present and matches topology;
/// otherwise falls back to `lattice_from_mesh`.
inline bool curvilinear_from_mesh(const Mesh& rMesh, std::array<std::int64_t, 3>& rDims) {
    if (auto annotated = get_structured_dims(rMesh, 0)) {
        if (!validate_structured_single_block(rMesh, *annotated))
            return false;
        rDims = *annotated;
        return true;
    }
    LatticeSpec spec;
    if (!lattice_from_mesh(rMesh, spec))
        return false;
    rDims = {spec.mDims[0] + 1, spec.mDims[1] + 1, spec.mDims[2] + 1};
    return true;
}

}  // namespace meshioplusplus::detail
