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
// octree: probes on a sphere skin and a circle polyline. The uniform case must
// reproduce voxelize's lattice; a deeper tree must keep leaves at max_level only
// where the skin cuts; the inside volume must converge to the analytic one; and
// the 2:1 rule must actually hold.

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
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::octree;
using meshioplusplus::OctreeBalance;
using meshioplusplus::OctreeOptions;
using meshioplusplus::OctreeResult;
using meshioplusplus::VoxelFill;

namespace {

constexpr double kPi = 3.14159265358979323846;

// A latitude/longitude sphere of radius R, outward wound, poles welded.
Mesh uv_sphere(double R, int Nu, int Nv) {
    std::vector<std::vector<double>> pts;
    pts.push_back({0, 0, -R});
    for (int j = 1; j < Nv; ++j)
        for (int i = 0; i < Nu; ++i) {
            const double th = kPi * j / Nv, ph = 2.0 * kPi * i / Nu;
            pts.push_back({R * std::sin(th) * std::cos(ph), R * std::sin(th) * std::sin(ph),
                           -R * std::cos(th)});
        }
    pts.push_back({0, 0, R});
    const std::int64_t south = 0, north = static_cast<std::int64_t>(pts.size()) - 1;
    auto at = [&](int j, int i) { return static_cast<std::int64_t>(1 + (j - 1) * Nu + i % Nu); };
    std::vector<std::vector<std::int64_t>> tris;
    for (int i = 0; i < Nu; ++i) {
        tris.push_back({south, at(1, i + 1), at(1, i)});
        tris.push_back({north, at(Nv - 1, i), at(Nv - 1, i + 1)});
    }
    for (int j = 1; j < Nv - 1; ++j)
        for (int i = 0; i < Nu; ++i) {
            tris.push_back({at(j, i), at(j, i + 1), at(j + 1, i + 1)});
            tris.push_back({at(j, i), at(j + 1, i + 1), at(j + 1, i)});
        }
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

// A closed N-gon polyline of radius R at z = 0.
Mesh circle(double R, int N, bool Close = true) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < N; ++i)
        pts.push_back({R * std::cos(2.0 * kPi * i / N), R * std::sin(2.0 * kPi * i / N), 0.0});
    std::vector<std::vector<std::int64_t>> lines;
    for (int i = 0; i + 1 < N; ++i)
        lines.push_back({i, i + 1});
    if (Close)
        lines.push_back({N - 1, 0});
    return mt::make_mesh(std::move(pts), "line", std::move(lines));
}

const std::int64_t* cell_data_i64(const Mesh& M, const char* pName) {
    return M.CellData(pName, 0).As<std::int64_t>();
}

}  // namespace

TEST(Octree, UniformTreeIsTheVoxelizeLattice) {
    const Mesh sphere = uv_sphere(1.0, 16, 8);
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1, -1, -1, 1, 1, 1}};
    o.mMinLevel = 3;
    o.mMaxLevel = 3;
    const OctreeResult r = octree(sphere, o);
    EXPECT_EQ(r.mNumLeaves, 512);

    meshioplusplus::VoxelOptions vo;
    vo.mBounds = o.mBounds;
    vo.mResolution = std::array<std::int64_t, 3>{{8, 8, 8}};
    const meshioplusplus::VoxelResult v = meshioplusplus::voxelize(sphere, vo);
    ASSERT_EQ(r.mMesh.NumPoints(), v.mMesh.NumPoints());
    const double* a = r.mMesh.Points().As<double>();
    const double* b = v.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints() * 3; ++i)
        EXPECT_EQ(a[i], b[i]);
    EXPECT_EQ(r.mNumHanging, 0);
}

TEST(Octree, MaxLevelLeavesOnlyWhereTheSkinCuts) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 4;
    const OctreeResult r = octree(uv_sphere(1.0, 24, 12), o);
    const std::int64_t* level = cell_data_i64(r.mMesh, meshioplusplus::kOctreeLevelName);
    const std::int64_t* cut = cell_data_i64(r.mMesh, meshioplusplus::kOctreeCutName);
    for (std::int64_t n = 0; n < r.mNumLeaves; ++n)
        if (cut[n])
            EXPECT_EQ(level[n], 4);
    EXPECT_EQ(r.mFinestLevel, 4);
}

TEST(Octree, InsideVolumeConvergesToTheBall) {
    double prev_err = 1e9;
    for (std::int64_t L : {3, 4, 5}) {
        OctreeOptions o;
        o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
        o.mMaxLevel = L;
        o.mFill = VoxelFill::Inside;
        const OctreeResult r = octree(uv_sphere(1.0, 48, 24), o);
        const std::int64_t* level = cell_data_i64(r.mMesh, meshioplusplus::kOctreeLevelName);
        double vol = 0.0;
        for (std::int64_t n = 0; n < r.mNumLeaves; ++n) {
            const double e = r.mRootSize / static_cast<double>(1 << level[n]);
            vol += e * e * e;
        }
        const double err = std::fabs(vol - 4.0 / 3.0 * kPi) / (4.0 / 3.0 * kPi);
        EXPECT_LT(err, prev_err);
        prev_err = err;
    }
    EXPECT_LT(prev_err, 0.1);
}

TEST(Octree, FaceBalanceHolds) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 4;
    o.mBalance = OctreeBalance::Face;
    const OctreeResult r = octree(uv_sphere(1.0, 24, 12), o);
    const Mesh& m = r.mMesh;
    const std::int64_t* level = cell_data_i64(m, meshioplusplus::kOctreeLevelName);
    const double* pts = m.Points().As<double>();
    const std::int64_t* conn = m.Cells(0).Conn().As<std::int64_t>();
    const std::size_t n = m.Cells(0).NumCells();
    struct Box { double lo[3], hi[3]; };
    std::vector<Box> boxes(n);
    for (std::size_t c = 0; c < n; ++c) {
        for (int k = 0; k < 3; ++k) {
            boxes[c].lo[k] = pts[conn[c * 8] * 3 + k];
            boxes[c].hi[k] = pts[conn[c * 8 + 6] * 3 + k];
        }
    }
    const double eps = 1e-12;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b) {
            // Face neighbours: touching along one axis, overlapping in the other two.
            int touching = 0, overlap = 0;
            for (int k = 0; k < 3; ++k) {
                if (std::fabs(boxes[a].hi[k] - boxes[b].lo[k]) < eps ||
                    std::fabs(boxes[b].hi[k] - boxes[a].lo[k]) < eps)
                    ++touching;
                else if (boxes[a].lo[k] < boxes[b].hi[k] - eps && boxes[b].lo[k] < boxes[a].hi[k] - eps)
                    ++overlap;
            }
            if (touching == 1 && overlap == 2)
                EXPECT_LE(std::llabs(level[a] - level[b]), 1);
        }
}

TEST(Octree, QuadtreeOfACircleConvergesToItsArea) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, 0, 1.2, 1.2, 0}};
    o.mMaxLevel = 6;
    o.mFill = VoxelFill::Inside;
    const OctreeResult r = octree(circle(1.0, 256), o);
    EXPECT_EQ(r.mDimension, 2);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "quad");
    const std::int64_t* level = cell_data_i64(r.mMesh, meshioplusplus::kOctreeLevelName);
    double area = 0.0;
    for (std::int64_t n = 0; n < r.mNumLeaves; ++n) {
        const double e = r.mRootSize / static_cast<double>(1 << level[n]);
        area += e * e;
    }
    EXPECT_NEAR(area, kPi, 0.05);
    EXPECT_GT(r.mNumHanging, 0);
}

TEST(Octree, OpenPolylineIsRefusedForInside) {
    OctreeOptions o;
    o.mMaxLevel = 3;
    o.mFill = VoxelFill::Inside;
    EXPECT_THROW(octree(circle(1.0, 32, /*Close=*/false), o), std::invalid_argument);
}

TEST(Octree, BudgetIsEnforced) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 5;
    o.mMaxCells = 100;
    EXPECT_THROW(octree(uv_sphere(1.0, 24, 12), o), std::invalid_argument);
}

namespace {

// Total volume of the cells of a mesh made only of hexahedra, pyramids and
// tetrahedra (every polyhedron is skipped), by tetrahedral decomposition.
double sum_cell_volume(const Mesh& M, bool& rAllPositive) {
    const double* p = M.Points().As<double>();
    const auto tet = [&](std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t d) {
        double e[3][3];
        const std::int64_t ids[3] = {b, c, d};
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                e[j][k] = p[ids[j] * 3 + k] - p[a * 3 + k];
        return (e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
                e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0])) /
               6.0;
    };
    double total = 0.0;
    rAllPositive = true;
    for (std::size_t b = 0; b < M.NumCellBlocks(); ++b) {
        const auto cb = M.Cells(b);
        const std::string type(cb.Type());
        if (type != "hexahedron" && type != "pyramid" && type != "tetra")
            continue;
        const std::int64_t* c = cb.Conn().As<std::int64_t>();
        const std::size_t n = cb.NodesPerCell();
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            const std::int64_t* r = c + i * n;
            double v = 0.0;
            if (type == "hexahedron")
                v = 6.0 * tet(r[0], r[1], r[3], r[4]);
            else if (type == "pyramid")
                v = tet(r[0], r[1], r[2], r[4]) + tet(r[0], r[2], r[3], r[4]);
            else
                v = tet(r[0], r[1], r[2], r[3]);
            if (!(v > 0.0))
                rAllPositive = false;
            total += v;
        }
    }
    return total;
}

}  // namespace

TEST(Octree, TemplatesConserveVolumeAndStayPositive) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 4;
    o.mTransition = meshioplusplus::OctreeTransition::Templates;
    const OctreeResult r = octree(uv_sphere(1.0, 24, 12), o);
    EXPECT_EQ(r.mNumHanging, 0);
    bool positive = false;
    EXPECT_NEAR(sum_cell_volume(r.mMesh, positive), 2.4 * 2.4 * 2.4, 1e-9);
    EXPECT_TRUE(positive);
}

TEST(Octree, PolyhedraNeedABalancedTree) {
    OctreeOptions o;
    o.mMaxLevel = 3;
    o.mBalance = OctreeBalance::None;
    o.mTransition = meshioplusplus::OctreeTransition::Polyhedra;
    EXPECT_THROW(octree(uv_sphere(1.0, 24, 12), o), std::invalid_argument);
}

TEST(Octree, SnapNeedsAFillAndMovesTheOuterNodes) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 4;
    o.mSnap = true;
    EXPECT_THROW(octree(uv_sphere(1.0, 48, 24), o), std::invalid_argument);
    o.mFill = VoxelFill::Inside;
    const OctreeResult r = octree(uv_sphere(1.0, 48, 24), o);
    EXPECT_GT(r.mNumSnapped, 0);
}

TEST(Octree, SurfaceLevelStopsCutLeavesEarly) {
    OctreeOptions o;
    o.mBounds = std::array<double, 6>{{-1.2, -1.2, -1.2, 1.2, 1.2, 1.2}};
    o.mMaxLevel = 5;
    const OctreeResult full = octree(uv_sphere(1.0, 48, 24), o);
    o.mSurfaceLevel = 3;
    const OctreeResult shallow = octree(uv_sphere(1.0, 48, 24), o);
    EXPECT_LT(shallow.mNumLeaves, full.mNumLeaves);
    o.mCurvatureThreshold = 0.1;
    EXPECT_EQ(octree(uv_sphere(1.0, 48, 24), o).mNumLeaves, full.mNumLeaves);
}
