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
// Tests for the geometric statistics operation.

// System includes
#include <cmath>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/stats.hpp"

namespace {

using meshioplusplus::compute_stats;
using meshioplusplus::Mesh;
using meshioplusplus::StatsReport;

Mesh unit_cube_hex() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
}

TEST(Stats, BBoxAndCentroid) {
    StatsReport s = compute_stats(unit_cube_hex());
    EXPECT_EQ(s.mNumPoints, 8);
    EXPECT_EQ(s.mNumCells, 1);
    EXPECT_DOUBLE_EQ(s.mBBoxMin[0], 0.0);
    EXPECT_DOUBLE_EQ(s.mBBoxMax[2], 1.0);
    EXPECT_DOUBLE_EQ(s.mExtent[1], 1.0);
    EXPECT_NEAR(s.mCentroid[0], 0.5, 1e-12);
    EXPECT_NEAR(s.mCentroid[2], 0.5, 1e-12);
    ASSERT_EQ(s.mCellTypeCounts.size(), 1u);
    EXPECT_EQ(s.mCellTypeCounts[0].first, "hexahedron");
    EXPECT_EQ(s.mCellTypeCounts[0].second, 1);
}

TEST(Stats, UnitCubeVolumeAndArea) {
    StatsReport s = compute_stats(unit_cube_hex());
    EXPECT_NEAR(s.mSignedVolume, 1.0, 1e-9);
    EXPECT_NEAR(s.mUnsignedVolume, 1.0, 1e-9);
    EXPECT_NEAR(s.mTotalArea, 6.0, 1e-9);
    EXPECT_EQ(s.mNumInverted, 0);
}

TEST(Stats, TriangleArea) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {2, 0, 0}, {0, 3, 0}}, "triangle", {{0, 1, 2}});
    StatsReport s = compute_stats(m);
    EXPECT_NEAR(s.mTotalArea, 3.0, 1e-12);
    EXPECT_DOUBLE_EQ(s.mSignedVolume, 0.0);
}

TEST(Stats, InvertedElementCounted) {
    // top/bottom swapped -> negative volume
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{4, 5, 6, 7, 0, 1, 2, 3}});
    StatsReport s = compute_stats(m);
    EXPECT_LT(s.mSignedVolume, 0.0);
    EXPECT_EQ(s.mNumInverted, 1);
}

}  // namespace
