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
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"

namespace {

namespace mio = meshioplusplus;
namespace d = meshioplusplus::detail;

// A paired tree: the roots split once, then every leaf whose centre (in
// units of the forest's extent) passes `Refine` split, `Rounds` times.
d::LinearOctree hds_tree(const std::array<std::int64_t, 3>& rRoots, std::int32_t MaxLevel,
                         int Rounds, const std::function<bool(double, double, double)>& rRefine) {
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

double hds_volume(const std::vector<d::Vec3>& rP, const std::array<std::int64_t, 8>& rH) {
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
d::HexDualMesh hds_check(const d::LinearOctree& rT) {
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

std::int64_t hds_schemes(const d::HexDualMesh& rD) {
    std::int64_t n = 0;
    for (std::int64_t c : rD.mNumSchemes)
        n += c;
    return n;
}

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
