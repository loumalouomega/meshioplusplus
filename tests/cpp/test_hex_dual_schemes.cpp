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
// detail/hex_dual_schemes: the all-hexahedral dual of paired octrees with flat,
// concave and convex transitions and random refinements. Every interior vertex
// becomes one hexahedron of positive volume, and the hexahedra conform.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <set>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/hex_dual_schemes.hpp"
#include "../../src/cpp/src/detail/hex_fit.hpp"
#include "../../src/cpp/src/detail/linear_octree.hpp"
#include "hex_dual_fixtures.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"

namespace {

namespace mio = meshioplusplus;
namespace d = meshioplusplus::detail;
using namespace hdsfix;
}  // namespace

TEST(HexDualSchemes, TheTablesAreUnitCubes) {
    for (std::size_t t = 0; t < d::kNumHexTransitions; ++t) {
        const d::HexDualScheme& s = d::hex_dual_scheme(static_cast<d::HexTransition>(t));
        std::array<int, 3> lo{20, 20, 20}, hi{0, 0, 0};
        for (std::size_t v = 0; v < s.mNumVerts; ++v)
            for (std::size_t k = 0; k < 3; ++k) {
                lo[k] = std::min<int>(lo[k], s.mVerts[v * 3 + k]);
                hi[k] = std::max<int>(hi[k], s.mVerts[v * 3 + k]);
            }
        EXPECT_EQ(lo, (std::array<int, 3>{0, 0, 0})) << t;
        EXPECT_EQ(hi, (std::array<int, 3>{20, 20, 20})) << t;
        // Every face is used once on the cube's boundary, twice inside.
        std::vector<int> uses(s.mNumFaces, 0);
        for (std::size_t i = 0; i < s.mPolyEnds[s.mNumPolys - 1]; ++i)
            ++uses[s.mPolys[i]];
        for (std::size_t f = 0, at = 0; f < s.mNumFaces; ++f) {
            bool outer = false;
            for (std::size_t k = 0; k < 3 && !outer; ++k) {
                bool l = true, h = true;
                for (std::size_t i = at; i < s.mFaceEnds[f]; ++i) {
                    l = l && s.mVerts[s.mFaces[i] * 3 + k] == 0;
                    h = h && s.mVerts[s.mFaces[i] * 3 + k] == 20;
                }
                outer = l || h;
            }
            EXPECT_EQ(uses[f], outer ? 1 : 2) << "scheme " << t << " face " << f;
            at = s.mFaceEnds[f];
        }
    }
}

TEST(HexDualSchemes, AUniformTreeIsItsOwnDual) {
    const d::LinearOctree t =
        hds_tree({{2, 2, 2}}, 3, 0, [](double, double, double) { return false; });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_EQ(dual.mNumTransitionVertices, 0);
    EXPECT_EQ(hds_schemes(dual), 0);
    EXPECT_EQ(dual.mHexes.size(), 27u);  // 4^3 cells, 3^3 interior vertices
}

TEST(HexDualSchemes, AFlatTransition) {
    const d::LinearOctree t =
        hds_tree({{2, 2, 2}}, 4, 1, [](double X, double, double) { return X < 0.5; });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_GT(dual.mNumTransitionVertices, 0);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::Flat)], 0);
}

TEST(HexDualSchemes, AConcaveEdge) {
    const d::LinearOctree t =
        hds_tree({{2, 2, 2}}, 4, 1, [](double X, double Y, double) { return X < 0.5 || Y < 0.5; });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::Edge)], 0);
}

TEST(HexDualSchemes, AConvexEdge) {
    const d::LinearOctree t =
        hds_tree({{2, 2, 2}}, 4, 1, [](double X, double Y, double) { return X < 0.5 && Y < 0.5; });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::Convex1)], 0);
}

TEST(HexDualSchemes, AConcaveCorner) {
    const d::LinearOctree t = hds_tree({{2, 2, 2}}, 4, 1, [](double X, double Y, double Z) {
        return X < 0.5 || Y < 0.5 || Z < 0.5;
    });
    const d::HexDualMesh dual = hds_check(t);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::VertCenter)], 0);
}

TEST(HexDualSchemes, AConvexCorner) {
    const d::LinearOctree t = hds_tree({{2, 2, 2}}, 4, 1, [](double X, double Y, double Z) {
        return X < 0.5 && Y < 0.5 && Z < 0.5;
    });
    hds_check(t);
}

TEST(HexDualSchemes, ABallTwoLevelsDown) {
    const d::LinearOctree t = hds_tree({{2, 2, 2}}, 5, 2, [](double X, double Y, double Z) {
        const double r =
            std::sqrt((X - 0.5) * (X - 0.5) + (Y - 0.5) * (Y - 0.5) + (Z - 0.5) * (Z - 0.5));
        return std::fabs(r - 0.3) < 0.12;
    });
    hds_check(t);
}

TEST(HexDualSchemes, RandomTrees) {
    int failed = 0;
    std::string first;
    for (std::uint32_t seed = 1; seed <= 12; ++seed) {
        std::mt19937 rng(seed);
        const d::LinearOctree t =
            hds_tree({{2, 2, 2}}, 5, 3, [&](double, double, double) { return rng() % 6 == 0; });
        try {
            hds_check(t);
        } catch (const std::exception& e) {
            if (failed++ == 0)
                first = std::to_string(seed) + ": " + e.what();
        }
    }
    EXPECT_EQ(failed, 0) << first;
}

TEST(HexDualSchemes, TheDualOfAShellIsFairWithoutFitting) {
    // Centroid duals of the schemes: the transition hexahedra's worst corner
    // stays above a quarter (0.283 here, 1 away from the transitions).
    const d::LinearOctree t = hds_tree({{2, 2, 2}}, 5, 2, [](double X, double Y, double Z) {
        const double r =
            std::sqrt((X - 0.5) * (X - 0.5) + (Y - 0.5) * (Y - 0.5) + (Z - 0.5) * (Z - 0.5));
        return std::fabs(r - 0.3) < 0.12;
    });
    const d::HexDualMesh dual = d::hex_dual_mesh(t);
    double worst = 1.0;
    for (const auto& h : dual.mHexes)
        worst = std::min(worst, d::hex_min_scaled_jacobian(dual.mPoints, h));
    EXPECT_GT(worst, 0.25);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::Flat)], 0);
    EXPECT_GT(dual.mNumSchemes[static_cast<std::size_t>(d::HexTransition::Edge)], 0);
}
