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
 * @file detail/tet_fem.hpp
 * @brief P1 tetrahedral FEM assembly (stiffness, lumped/consistent mass,
 *        per-cell gradients) over the uniform mesh API (roadmap 15.2.1).
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"

#include "fem_kernels.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus::detail {

struct TetFemCell {
    std::array<std::int64_t, 4> mNodes;
    std::array<Vec3, 4> mGrad;
    double mVolume = 0.0;
    std::int64_t mIndex = 0;
};

struct TetFem {
    std::vector<Vec3> mPoints;
    std::vector<TetFemCell> mCells;
    std::int64_t mNumInput = 0;
    std::int64_t mNumDegenerate = 0;
};

/// Returns true iff @p rMesh has at least one Tetra block and only Tetra blocks.
inline bool mesh_is_pure_tetra(const Mesh& rMesh) {
    bool has_tet = false;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.NumCells() == 0)
            continue;
        if (cell_type_from_name(cb.Type()) != CellType::Tetra)
            return false;
        has_tet = true;
    }
    return has_tet;
}

/// Build P1 tetrahedral geometry and gradients via `tet_p1_gradients`.
inline TetFem build_tet_fem(const Mesh& rMesh, const std::string& rOp) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    if (pdim < 3)
        throw std::invalid_argument("meshio++: " + rOp + ": tetrahedra require 3-D points");
    const NDArray& pts = rMesh.Points();
    TetFem out;
    out.mPoints.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        out.mPoints[i] = {read_double(pts, i * pdim + 0), read_double(pts, i * pdim + 1),
                          read_double(pts, i * pdim + 2)};
    }
    std::int64_t cell_idx = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.NumCells() == 0)
            continue;
        if (cell_type_from_name(cb.Type()) != CellType::Tetra)
            throw std::invalid_argument("meshio++: " + rOp +
                                        ": mixed cell types are not supported");
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t id = cell_idx++;
            const std::array<std::int64_t, 4> nodes = {
                read_int(conn, c * 4 + 0), read_int(conn, c * 4 + 1),
                read_int(conn, c * 4 + 2), read_int(conn, c * 4 + 3)};
            const Vec3 p[4] = {out.mPoints[static_cast<std::size_t>(nodes[0])],
                               out.mPoints[static_cast<std::size_t>(nodes[1])],
                               out.mPoints[static_cast<std::size_t>(nodes[2])],
                               out.mPoints[static_cast<std::size_t>(nodes[3])]};
            Vec3 grad[4];
            double vol = 0.0;
            if (!tet_p1_gradients(p, grad, vol)) {
                ++out.mNumDegenerate;
                continue;
            }
            TetFemCell tc;
            tc.mNodes = nodes;
            for (int k = 0; k < 4; ++k)
                tc.mGrad[static_cast<std::size_t>(k)] = grad[k];
            tc.mVolume = vol;
            tc.mIndex = id;
            out.mCells.push_back(tc);
        }
    }
    out.mNumInput = cell_idx;
    if (out.mNumInput == 0)
        throw std::invalid_argument("meshio++: " + rOp + ": mesh has no tetrahedra");
    return out;
}

/// Assemble lumped mass and P1 stiffness triplets on a tetrahedral mesh.
inline void tet_fem_assemble(const TetFem& rFem, std::vector<double>& rMass,
                             std::vector<SparseTriplet>& rStiffness) {
    rMass.assign(rFem.mPoints.size(), 0.0);
    rStiffness.reserve(rFem.mCells.size() * 16);
    for (const TetFemCell& tc : rFem.mCells) {
        const double m_lump = tc.mVolume * 0.25;
        double K[4][4];
        tet_p1_stiffness(tc.mGrad.data(), tc.mVolume, K);
        for (int a = 0; a < 4; ++a) {
            const auto na = static_cast<std::size_t>(tc.mNodes[static_cast<std::size_t>(a)]);
            rMass[na] += m_lump;
            for (int b = 0; b < 4; ++b) {
                const auto nb = static_cast<std::size_t>(tc.mNodes[static_cast<std::size_t>(b)]);
                rStiffness.push_back({na, nb, K[a][b]});
            }
        }
    }
}

}  // namespace meshioplusplus::detail
