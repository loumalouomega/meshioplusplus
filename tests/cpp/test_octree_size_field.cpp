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
// detail/octree_size_field: the frame's levels and root size, welding, the
// feature edges of a cube, the curvature size on a sphere, refinement that
// only ever adds leaves when a criterion is added, and the budget and lattice
// limits refused by name.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "../../src/cpp/src/detail/octree_size_field.hpp"

namespace {

namespace mio = meshioplusplus;
namespace det = meshioplusplus::detail;

constexpr const char* kPrefix = "meshio++: test: ";

det::SizeField sf_field(const det::OctreeSkin& rSkin, const det::TriBvh& rTree, double HMin,
                        double HMax, double Curvature, double Gap) {
    det::SizeFieldSpec spec;
    spec.mMinSize = HMin;
    spec.mMaxSize = HMax;
    spec.mCurvatureAngleDeg = Curvature;
    spec.mGapCells = Gap;
    return det::octree_size_field(rSkin, rTree, {}, spec);
}

}  // namespace

TEST(OctreeSizeField, FrameLevelsAndRootSize) {
    const det::OctreeFrame f =
        det::octree_frame({0.0, 0.0, 0.0}, {3.0, 1.0, 0.5}, 0.1, 0.8, -1, kPrefix);
    EXPECT_EQ(f.mMaxLevel, 3);
    EXPECT_DOUBLE_EQ(f.mRootSize, 0.8);
    EXPECT_DOUBLE_EQ(f.mMinSize, 0.1);
    EXPECT_EQ(f.mRoots, (std::array<std::int64_t, 3>{{4, 2, 1}}));
    // Centred on the box.
    EXPECT_DOUBLE_EQ(f.mOrigin[0] + 0.5 * 4 * 0.8, 1.5);
    EXPECT_EQ(f.LevelFor(0.8), 0);
    EXPECT_EQ(f.LevelFor(0.4), 1);
    EXPECT_EQ(f.LevelFor(0.39), 2);
    EXPECT_EQ(f.LevelFor(0.1), 3);
    EXPECT_EQ(f.LevelFor(0.01), 3);
    EXPECT_EQ(f.LevelFor(5.0), 0);
    // A cap on the level coarsens the finest leaf.
    const det::OctreeFrame c =
        det::octree_frame({0.0, 0.0, 0.0}, {3.0, 1.0, 0.5}, 0.1, 0.8, 1, kPrefix);
    EXPECT_EQ(c.mMaxLevel, 1);
    EXPECT_DOUBLE_EQ(c.mMinSize, 0.4);
    EXPECT_THROW(det::octree_frame({0, 0, 0}, {1e4, 1, 1}, 1e-3, 1e-3, -1, kPrefix),
                 std::invalid_argument);
    EXPECT_THROW(det::octree_frame({0, 0, 0}, {1, 1, 1}, 1e-9, 1.0, -1, kPrefix),
                 std::invalid_argument);
}

TEST(OctreeSizeField, SkinWeldsAndFindsTheCubesFeatures) {
    // A box surface whose blocks repeat the shared corner points.
    const mio::Mesh cube = mio::extract_surface(mio::box({{2, 2, 2}}, {{1.0, 1.0, 1.0}}));
    const det::OctreeSkin skin = det::octree_skin(cube, "", kPrefix);
    EXPECT_EQ(skin.mPoints.size(), 26u);  // 3^3 - 1
    EXPECT_EQ(skin.mTris.size(), 48u);    // 24 quads fanned
    // 12 cube edges of 2 segments each.
    EXPECT_EQ(det::octree_skin_features(skin, 30.0).size(), 24u);
    EXPECT_TRUE(det::octree_skin_features(skin, 0.0).empty());
    EXPECT_THROW(det::octree_skin(mt::line_mesh(), "", kPrefix), std::invalid_argument);
}

TEST(OctreeSizeField, CurvatureSizeOnASphere) {
    const mio::Mesh sphere = mio::icosa_sphere(2.0, 4);
    const det::OctreeSkin skin = det::octree_skin(sphere, "", kPrefix);
    const det::TriBvh tree = det::build_tri_bvh(skin.mSoup);
    const det::SizeField f = sf_field(skin, tree, 0.01, 10.0, 15.0, 0.0);
    const double want = 2.0 * std::sin(7.5 * 3.14159265358979323846 / 180.0) * 2.0;
    for (double h : f.mPointSize)
        EXPECT_NEAR(h, want, 0.05 * want);
    // The gap term on a sphere of radius 2 with 4 cells: its diameter / 4.
    const det::SizeField g = sf_field(skin, tree, 0.01, 10.0, 0.0, 4.0);
    for (double h : g.mPointSize)
        EXPECT_NEAR(h, 1.0, 0.02);
}

TEST(OctreeSizeField, AddingACriterionOnlyRefines) {
    const mio::Mesh skin_mesh = analytic::torus(1.0, 0.3, 48, 24);
    const det::OctreeSkin skin = det::octree_skin(skin_mesh, "", kPrefix);
    const det::TriBvh tree = det::build_tri_bvh(skin.mSoup);
    const det::OctreeFrame frame =
        det::octree_frame({-1.5, -1.5, -0.5}, {1.5, 1.5, 0.5}, 0.02, 0.5, -1, kPrefix);
    const det::SizeField plain = sf_field(skin, tree, 0.02, 0.5, 0.0, 0.0);
    const det::SizeField both = sf_field(skin, tree, 0.02, 0.5, 15.0, 2.0);
    const auto a = det::octree_adaptive_build(frame, plain.mPrims, -1, 0);
    const auto b = det::octree_adaptive_build(frame, both.mPrims, -1, 0);
    EXPECT_GT(b.mTree.Size(), a.mTree.Size());
    std::string why;
    EXPECT_TRUE(det::linear_octree_check(b.mTree, nullptr, &why)) << why;
    // Every leaf of `b` lies inside a leaf of `a` at least as coarse.
    for (std::size_t i = 0; i < b.mTree.Size(); ++i) {
        const auto lo = b.mTree.Lo(i);
        const std::int64_t j = a.mTree.Locate(lo[0], lo[1], lo[2]);
        ASSERT_GE(j, 0);
        EXPECT_LE(a.mTree.mLevel[static_cast<std::size_t>(j)], b.mTree.mLevel[i]);
    }
    // Cut flags are exactly the leaves with candidates.
    for (std::size_t i = 0; i < b.mTree.Size(); ++i)
        EXPECT_EQ(b.mCut[i] != 0, !b.mCand[i].empty());
    EXPECT_THROW(det::octree_adaptive_build(frame, both.mPrims, -1, 100), std::invalid_argument);
}
