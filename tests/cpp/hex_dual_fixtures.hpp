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
// Fixtures shared by the hexahedral dual tests (`test_hex_dual_schemes.cpp`, the full ones,
// and `test_quick_hex_dual_schemes.cpp`, their small versions that run in the pull-request
// coverage job): a paired octree refined by a predicate, and the check of its dual.

#pragma once

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/hex_dual_schemes.hpp"
#include "../../src/cpp/src/detail/hex_fit.hpp"
#include "../../src/cpp/src/detail/linear_octree.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"

namespace hdsfix {

namespace mio = meshioplusplus;
namespace d = meshioplusplus::detail;

// A paired tree: the roots split once, then every leaf whose centre (in
// units of the forest's extent) passes `Refine` split, `Rounds` times.
inline d::LinearOctree hds_tree(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel,
                                int Rounds,
                                const std::function<bool(double, double, double)>& rRefine) {
    d::LinearOctree t = d::linear_octree_roots(rRoots, MaxLevel);
    std::vector<std::uint8_t> all(t.Size(), 1);
    d::linear_octree_split(t, all, 0);
    for (int r = 0; r < Rounds; ++r) {
        std::vector<std::uint8_t> mark(t.Size(), 0);
        for (std::size_t i = 0; i < t.Size(); ++i) {
            const auto lo = t.Lo(i);
            const double s = static_cast<double>(t.Span(i));
            mark[i] =
                rRefine((static_cast<double>(lo[0]) + 0.5 * s) / static_cast<double>(t.Extent(0)),
                        (static_cast<double>(lo[1]) + 0.5 * s) / static_cast<double>(t.Extent(1)),
                        (static_cast<double>(lo[2]) + 0.5 * s) / static_cast<double>(t.Extent(2)));
        }
        d::linear_octree_split(t, mark, 0);
        d::linear_octree_pair(t, d::OctreeAdjacency::Full, 0);
    }
    d::linear_octree_pair(t, d::OctreeAdjacency::Full, 0);
    return t;
}

inline double hds_volume(const std::vector<d::Vec3>& rP, const std::array<std::int64_t, 8>& rH) {
    static constexpr int kTets[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                        {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    double vol = 0.0;
    for (const auto& t : kTets) {
        const d::Vec3& a = rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[0])])];
        const d::Vec3 b =
            d::vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[1])])], a);
        const d::Vec3 c =
            d::vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[2])])], a);
        const d::Vec3 e =
            d::vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[3])])], a);
        vol += d::vec3_dot(b, d::vec3_cross(c, e)) / 6.0;
    }
    return vol;
}

// The dual of the tree, checked: one hexahedron per interior vertex, every one
// of positive volume and a conforming, manifold mesh.
inline d::HexDualMesh hds_check(const d::LinearOctree& rT) {
    const d::HexDualMesh dual = d::hex_dual_mesh(rT);
    EXPECT_FALSE(dual.mHexes.empty());
    for (const auto& h : dual.mHexes)
        EXPECT_GT(hds_volume(dual.mPoints, h), 0.0);
    std::vector<std::vector<double>> pts;
    for (const auto& p : dual.mPoints)
        pts.push_back({p[0], p[1], p[2]});
    std::vector<std::vector<std::int64_t>> cells;
    for (const auto& h : dual.mHexes)
        cells.emplace_back(h.begin(), h.end());
    const mio::Mesh m = mt::make_mesh(pts, "hexahedron", cells);
    const mio::ConformityResult r = mio::check_conformity(m);
    EXPECT_TRUE(r.mConforming) << r.mNumHangingNodes << " hanging, " << r.mNumPartialOverlaps
                               << " overlaps, " << r.mNumUnweldedFacets << " unwelded";
    EXPECT_EQ(r.mNumNonManifoldFacets, 0);
    return dual;
}

inline std::int64_t hds_schemes(const d::HexDualMesh& rD) {
    std::int64_t n = 0;
    for (std::int64_t c : rD.mNumSchemes)
        n += c;
    return n;
}

}  // namespace hdsfix
