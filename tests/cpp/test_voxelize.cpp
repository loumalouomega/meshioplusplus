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
// voxelize: the composition of the lattice and the distance kernel, each of
// which is tested on its own elsewhere. What is left to check here is the
// composition itself -- which cells a fill rule keeps, and that the result is
// still an ordinary, well-formed mesh.

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
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/tri_box.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::compute_stats;
using meshioplusplus::Mesh;
using meshioplusplus::VoxelFill;
using meshioplusplus::voxelize;
using meshioplusplus::VoxelOptions;
using meshioplusplus::VoxelResult;
namespace d = meshioplusplus::detail;

namespace {

// The unit cube [0,1]^3 as a closed triangle surface, outward wound.
Mesh cube_surface() {
    std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<std::vector<std::int64_t>> tris = {
        {0, 2, 1}, {0, 3, 2},  // z = 0, outward -z
        {4, 5, 6}, {4, 6, 7},  // z = 1, outward +z
        {0, 1, 5}, {0, 5, 4},  // y = 0
        {1, 2, 6}, {1, 6, 5},  // x = 1
        {2, 3, 7}, {2, 7, 6},  // y = 1
        {3, 0, 4}, {3, 4, 7},  // x = 0
    };
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

VoxelOptions res(std::int64_t N, VoxelFill Fill) {
    VoxelOptions o;
    o.mResolution = std::array<std::int64_t, 3>{{N, N, N}};
    o.mFill = Fill;
    return o;
}

}  // namespace

TEST(Voxelize, FillAllCoversTheBoundingBox) {
    const VoxelResult r = voxelize(cube_surface(), res(4, VoxelFill::All));
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 64u);
    EXPECT_EQ(r.mNumOccupied, 64);
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_DOUBLE_EQ(r.mOrigin[k], 0.0);
        EXPECT_DOUBLE_EQ(r.mSpacing[k], 0.25);
        EXPECT_EQ(r.mDims[k], 4);
    }
    // The whole point of the design: what comes out is an ordinary mesh.
    const meshioplusplus::StatsReport st = compute_stats(r.mMesh);
    EXPECT_EQ(st.mNumInverted, 0);
    EXPECT_DOUBLE_EQ(st.mSignedVolume, 1.0);
}

TEST(Voxelize, PaddingGrowsTheBoxOnEverySide) {
    VoxelOptions o = res(2, VoxelFill::All);
    o.mPadding = 0.5;
    const VoxelResult r = voxelize(cube_surface(), o);
    EXPECT_DOUBLE_EQ(r.mOrigin[0], -0.5);
    EXPECT_DOUBLE_EQ(r.mSpacing[0], 1.0);  // a 2-unit box in 2 cells
}

TEST(Voxelize, ExplicitBoundsOverrideTheMesh) {
    VoxelOptions o = res(2, VoxelFill::All);
    o.mBounds = std::array<double, 6>{{-1.0, -1.0, -1.0, 1.0, 1.0, 1.0}};
    const VoxelResult r = voxelize(cube_surface(), o);
    EXPECT_DOUBLE_EQ(r.mOrigin[0], -1.0);
    EXPECT_DOUBLE_EQ(r.mSpacing[0], 1.0);
}

TEST(Voxelize, FillSurfaceKeepsAShellAndLeavesTheInteriorEmpty) {
    // A cube at resolution 6 has a 6x6x6 lattice aligned with its own faces. The
    // surface shell is every cell except the 4x4x4 interior... except that the
    // cube's faces lie exactly ON cell boundaries, so the outermost ring of
    // cells on each side touches them. Rather than pin an exact count against
    // that reasoning, assert the structural facts: the interior is empty, the
    // boundary is not, and the kept set is a strict subset.
    const VoxelResult r = voxelize(cube_surface(), res(6, VoxelFill::Surface));
    EXPECT_GT(r.mNumOccupied, 0);
    EXPECT_LT(r.mNumOccupied, 216);

    // No kept cell may have its centre deep inside the cube.
    const Mesh& m = r.mMesh;
    const double* p = m.Points().As<double>();
    const std::int64_t* conn = m.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < m.Cells(0).NumCells(); ++c) {
        d::Vec3 centre{{0.0, 0.0, 0.0}};
        for (std::size_t i = 0; i < 8; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                centre[k] += p[static_cast<std::size_t>(conn[c * 8 + i]) * 3 + k] / 8.0;
        const double margin =
            std::min(std::min(centre[0], centre[1]),
                     std::min(centre[2], std::min(1.0 - centre[0],
                                                  std::min(1.0 - centre[1], 1.0 - centre[2]))));
        EXPECT_LT(margin, 1.0 / 6.0 + 1e-12)
            << "a cell well inside the cube was kept as surface";
    }
}

TEST(Voxelize, FillSurfaceUsesExactOverlapNotBoundingBoxes) {
    // A single thin diagonal triangle. Its bounding box covers the whole lattice,
    // so a bounding-box test would keep every cell; the exact separating-axis
    // test keeps only the ones it actually passes through.
    const Mesh diagonal =
        mt::make_mesh({{0, 0, 0}, {1, 1, 0}, {1, 1, 1}}, "triangle", {{0, 1, 2}});
    VoxelOptions o = res(8, VoxelFill::Surface);
    o.mDistance.mWatertightCheck = meshioplusplus::SdfWatertightCheck::Off;
    const VoxelResult r = voxelize(diagonal, o);
    EXPECT_GT(r.mNumOccupied, 0);
    EXPECT_LT(r.mNumOccupied, 512) << "a bounding-box test would have kept every cell";
}

TEST(Voxelize, FillInsideKeepsTheInteriorOnly) {
    // The cube's faces lie on cell boundaries at even resolutions, so use an odd
    // one and a little padding: cell centres then sit unambiguously inside or
    // outside and the count is exact.
    VoxelOptions o = res(5, VoxelFill::Inside);
    o.mBounds = std::array<double, 6>{{-0.5, -0.5, -0.5, 1.5, 1.5, 1.5}};
    const VoxelResult r = voxelize(cube_surface(), o);
    // Cell size 0.4; centres at -0.3, 0.1, 0.5, 0.9, 1.3. Three of the five lie
    // in (0, 1) on each axis, so 27 cells are inside.
    EXPECT_EQ(r.mNumOccupied, 27);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), 27u);
    const meshioplusplus::StatsReport st = compute_stats(r.mMesh);
    EXPECT_EQ(st.mNumInverted, 0);
    EXPECT_NEAR(st.mSignedVolume, 27.0 * 0.4 * 0.4 * 0.4, 1e-12);
}

TEST(Voxelize, FillInsideRefusesAnUnsignedRequest) {
    VoxelOptions o = res(4, VoxelFill::Inside);
    o.mDistance.mSign = meshioplusplus::SdfSign::Unsigned;
    try {
        voxelize(cube_surface(), o);
        FAIL() << "fill 'inside' accepted an unsigned request";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("unsigned"), std::string::npos) << e.what();
    }
}

TEST(Voxelize, ASelectiveFillPrunesTheUnreferencedPoints) {
    const VoxelResult all = voxelize(cube_surface(), res(6, VoxelFill::All));
    const VoxelResult shell = voxelize(cube_surface(), res(6, VoxelFill::Surface));
    EXPECT_LT(shell.mMesh.NumPoints(), all.mMesh.NumPoints())
        << "the subset kept points no cell references";
    // Every connectivity entry must be in range of the compacted point array.
    const std::int64_t n = static_cast<std::int64_t>(shell.mMesh.NumPoints());
    const std::int64_t* conn = shell.mMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < shell.mMesh.Cells(0).NumCells() * 8; ++i) {
        ASSERT_GE(conn[i], 0);
        ASSERT_LT(conn[i], n);
    }
}

TEST(Voxelize, OccupancyIsAttachedOnRequest) {
    VoxelOptions o = res(3, VoxelFill::All);
    o.mAttachOccupancy = true;
    const VoxelResult r = voxelize(cube_surface(), o);
    ASSERT_TRUE(r.mMesh.HasCellData("voxel:occupancy"));
    const meshioplusplus::NDArray& occ = r.mMesh.CellData("voxel:occupancy", 0);
    EXPECT_EQ(occ.Size(), 27u);
    for (std::size_t c = 0; c < occ.Size(); ++c)
        EXPECT_EQ(occ.As<std::int64_t>()[c], 1);
}

TEST(Voxelize, ResolutionAndCellSizeAreMutuallyExclusive) {
    VoxelOptions neither;
    EXPECT_THROW(voxelize(cube_surface(), neither), std::invalid_argument);

    VoxelOptions both;
    both.mResolution = std::array<std::int64_t, 3>{{2, 2, 2}};
    both.mCellSize = 0.5;
    try {
        voxelize(cube_surface(), both);
        FAIL() << "both selectors were accepted";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("exactly one"), std::string::npos) << e.what();
    }
}

TEST(Voxelize, CellSizeCoversTheBox) {
    VoxelOptions o;
    o.mCellSize = 0.3;  // does not divide 1.0
    const VoxelResult r = voxelize(cube_surface(), o);
    EXPECT_EQ(r.mDims[0], 4);  // ceil(1.0 / 0.3)
    EXPECT_DOUBLE_EQ(r.mSpacing[0], 0.3);
    EXPECT_GE(r.mOrigin[0] + 4.0 * 0.3, 1.0);
}

TEST(Voxelize, TheCellBudgetRefusesByName) {
    VoxelOptions o = res(100, VoxelFill::All);
    o.mMaxCells = 1000;
    try {
        voxelize(cube_surface(), o);
        FAIL() << "the cell budget did not refuse";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("max_cells"), std::string::npos) << e.what();
    }
}

TEST(Voxelize, FillNamesParse) {
    EXPECT_EQ(meshioplusplus::voxel_fill_from_name("all"), VoxelFill::All);
    EXPECT_EQ(meshioplusplus::voxel_fill_from_name("surface"), VoxelFill::Surface);
    EXPECT_EQ(meshioplusplus::voxel_fill_from_name("inside"), VoxelFill::Inside);
    EXPECT_THROW(meshioplusplus::voxel_fill_from_name("solid"), std::invalid_argument);
}

TEST(Voxelize, OutputIsStableAcrossRepeatedRuns) {
    // Determinism is structural here -- the lattice is index arithmetic and the
    // fill is a totally ordered search -- but the claim is cheap to check.
    const VoxelResult a = voxelize(cube_surface(), res(5, VoxelFill::Inside));
    const VoxelResult b = voxelize(cube_surface(), res(5, VoxelFill::Inside));
    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    ASSERT_EQ(a.mMesh.Cells(0).NumCells(), b.mMesh.Cells(0).NumCells());
    const double* pa = a.mMesh.Points().As<double>();
    const double* pb = b.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < a.mMesh.NumPoints() * 3; ++i)
        ASSERT_EQ(pa[i], pb[i]) << i;
    const std::int64_t* ca = a.mMesh.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* cb = b.mMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < a.mMesh.Cells(0).NumCells() * 8; ++i)
        ASSERT_EQ(ca[i], cb[i]) << i;
}

// --- the separating-axis primitive ------------------------------------------

TEST(TriBox, OverlapIsExactNotABoundingBoxTest) {
    const d::Vec3 centre{{0.0, 0.0, 0.0}};
    const d::Vec3 half{{0.5, 0.5, 0.5}};

    // Straight through the middle.
    EXPECT_TRUE(d::tri_box_overlap(centre, half, {{-2.0, 0.0, 0.0}}, {{2.0, 0.0, 0.0}},
                                   {{0.0, 2.0, 0.0}}));
    // Far away.
    EXPECT_FALSE(d::tri_box_overlap(centre, half, {{5.0, 5.0, 5.0}}, {{6.0, 5.0, 5.0}},
                                    {{5.0, 6.0, 5.0}}));
    // The case a bounding-box test gets wrong: a triangle whose bounding box
    // covers the origin but which passes well clear of it.
    EXPECT_FALSE(d::tri_box_overlap(centre, half, {{-3.0, 3.0, 0.0}}, {{3.0, 9.0, 0.0}},
                                    {{-3.0, 9.0, 0.0}}));
    // Touching the face counts as overlapping -- the documented tie rule.
    EXPECT_TRUE(d::tri_box_overlap(centre, half, {{0.5, 0.0, 0.0}}, {{2.0, 0.0, 0.0}},
                                   {{2.0, 1.0, 0.0}}));
}

TEST(TriBox, AnEdgeOnlyCrossingIsDetected) {
    // The nine edge-cross-axis tests exist for exactly this: a triangle that
    // clips a box corner without any vertex inside and without the box centre
    // being in its plane's slab.
    const d::Vec3 centre{{0.0, 0.0, 0.0}};
    const d::Vec3 half{{1.0, 1.0, 1.0}};
    EXPECT_TRUE(d::tri_box_overlap(centre, half, {{0.9, 0.9, -3.0}}, {{0.9, 0.9, 3.0}},
                                   {{3.0, 3.0, 0.0}}));
}
