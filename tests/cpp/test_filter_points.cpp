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
// Point-cloud conditioning (roadmap 9.9): outlier removal, voxel and
// Poisson-disk thinning, what rides along, and the refusals.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/filter_points.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::filter_points;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::PointFilterOptions;

using Xyz = std::array<double, 3>;

Mesh fp_cloud(const std::vector<Xyz>& rPts) {
    std::vector<std::vector<double>> p;
    for (const auto& q : rPts)
        p.push_back({q[0], q[1], q[2]});
    return mt::make_mesh(p, "vertex", [&] {
        std::vector<std::vector<std::int64_t>> c;
        for (std::size_t i = 0; i < rPts.size(); ++i)
            c.push_back({static_cast<std::int64_t>(i)});
        return c;
    }());
}

Xyz fp_at(const Mesh& rMesh, std::size_t I) {
    const double* c = rMesh.Points().As<double>();
    const std::size_t d = rMesh.PointDim();
    return {c[I * d], c[I * d + 1], c[I * d + 2]};
}

// A tight noisy patch of the plane z = 0 plus a few far points.
std::vector<Xyz> fp_scan(std::size_t NInliers, std::size_t NOutliers) {
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::normal_distribution<double> noise(0.0, 0.002);
    std::vector<Xyz> pts;
    for (std::size_t i = 0; i < NInliers; ++i)
        pts.push_back({u(rng), u(rng), noise(rng)});
    for (std::size_t i = 0; i < NOutliers; ++i)
        pts.push_back({3.0 + u(rng), 3.0 + u(rng), 2.0 + u(rng)});
    return pts;
}

double fp_dist(const Xyz& a, const Xyz& b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                     (a[2] - b[2]) * (a[2] - b[2]));
}

}  // namespace

TEST(FilterPoints, StatisticalRemovesFarPointsAndKeepsThePatch) {
    const auto r = filter_points(fp_cloud(fp_scan(1500, 12)), [] {
        PointFilterOptions o;
        o.mStatisticalK = 12;
        o.mStdRatio = 2.0;
        return o;
    }());
    EXPECT_GE(r.mNumRemovedStatistical, 12);
    EXPECT_LE(r.mNumRemovedStatistical, 12 + 60);  // a few tail inliers may go too
    // Every output point is an input point; no outlier survives.
    for (std::size_t i = 0; i < r.mPointMap.size(); ++i)
        EXPECT_LT(r.mPointMap[i], 1500) << i;
    EXPECT_EQ(r.mNumInput, 1512);
    ASSERT_EQ(r.mMesh.NumPoints(), r.mPointMap.size());
}

TEST(FilterPoints, RadiusRemovesIsolatedPoints) {
    PointFilterOptions o;
    o.mRadius = 0.08;
    o.mMinNeighbors = 5;
    const auto r = filter_points(fp_cloud(fp_scan(2000, 8)), o);
    EXPECT_GE(r.mNumRemovedRadius, 8);
    for (const auto id : r.mPointMap)
        EXPECT_LT(id, 2000);
    // Every kept point has at least five others within the radius (in the input).
    const auto pts = fp_scan(2000, 8);
    for (const auto id : r.mPointMap) {
        int n = 0;
        for (std::size_t j = 0; j < pts.size(); ++j)
            if (static_cast<std::int64_t>(j) != id && fp_dist(pts[id], pts[j]) <= 0.08)
                ++n;
        EXPECT_GE(n, 5);
    }
}

TEST(FilterPoints, VoxelKeepsOneRealPointPerOccupiedCube) {
    // A 10 x 10 x 10 lattice of spacing 1: a voxel of side 2 holds eight points.
    std::vector<Xyz> pts;
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 10; ++j)
            for (int k = 0; k < 10; ++k)
                pts.push_back({double(i), double(j), double(k)});
    PointFilterOptions o;
    o.mVoxelSize = 2.0;
    const auto r = filter_points(fp_cloud(pts), o);
    EXPECT_EQ(r.mMesh.NumPoints(), 125u);
    EXPECT_EQ(r.mNumRemovedVoxel, 1000 - 125);
    // Output is ascending in (x, y, z) voxel key, and each survivor is a real
    // input point: the lowest-id of the eight points nearest their mean.
    std::tuple<int, int, int> prev{-1, -1, -1};
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const Xyz p = fp_at(r.mMesh, i);
        const std::tuple<int, int, int> key{int(p[0] / 2), int(p[1] / 2), int(p[2] / 2)};
        EXPECT_GT(key, prev);
        prev = key;
        EXPECT_EQ(pts[r.mPointMap[i]], p);
    }
}

TEST(FilterPoints, PoissonDiskIsSeparatedAndMaximal) {
    const auto pts = fp_scan(3000, 0);
    PointFilterOptions o;
    o.mMinDistance = 0.05;
    const auto r = filter_points(fp_cloud(pts), o);
    const std::size_t n = r.mMesh.NumPoints();
    ASSERT_GT(n, 10u);
    ASSERT_LT(n, 3000u);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            EXPECT_GE(fp_dist(fp_at(r.mMesh, i), fp_at(r.mMesh, j)), 0.05);
    // Maximal: every dropped point is within the distance of a kept one.
    const std::set<std::int64_t> kept(r.mPointMap.begin(), r.mPointMap.end());
    for (std::size_t j = 0; j < pts.size(); ++j) {
        if (kept.count(static_cast<std::int64_t>(j)))
            continue;
        bool near = false;
        for (const auto id : r.mPointMap)
            if (fp_dist(pts[j], pts[id]) < 0.05) {
                near = true;
                break;
            }
        EXPECT_TRUE(near) << j;
    }
    // Greedy in index order: point 0 always survives.
    EXPECT_EQ(r.mPointMap[0], 0);
}

TEST(FilterPoints, StepsRunInOrderAndCount) {
    PointFilterOptions o;
    o.mStatisticalK = 12;
    o.mStdRatio = 2.0;
    o.mVoxelSize = 0.1;
    o.mMinDistance = 0.12;
    const auto r = filter_points(fp_cloud(fp_scan(2000, 10)), o);
    EXPECT_EQ(r.mNumInput - r.mNumRemovedStatistical - r.mNumRemovedRadius - r.mNumRemovedVoxel -
                  r.mNumRemovedPoisson,
              static_cast<std::int64_t>(r.mMesh.NumPoints()));
    EXPECT_GT(r.mNumRemovedVoxel, 0);
    EXPECT_GT(r.mNumRemovedPoisson, 0);
}

TEST(FilterPoints, PointDataRegionsAndCellsRideAlong) {
    std::vector<Xyz> pts;
    for (int i = 0; i < 6; ++i)
        pts.push_back({double(i), 0, 0});
    pts.push_back({100, 100, 100});  // the outlier
    Mesh m = fp_cloud(pts);
    m.AddPointData("id", mt::data_array({0, 1, 2, 3, 4, 5, 6}));
    m.AddFieldData("meta", mt::data_array({1.0, 2.0}));
    NDArray ids = NDArray::Uninit(meshioplusplus::DType::Int64, {std::size_t{3}});
    ids.As<std::int64_t>()[0] = 1;
    ids.As<std::int64_t>()[1] = 6;
    ids.As<std::int64_t>()[2] = 3;
    m.AddRegion(meshioplusplus::Region("keep", meshioplusplus::RegionKind::Point, ids));
    NDArray gone = NDArray::Uninit(meshioplusplus::DType::Int64, {std::size_t{1}});
    gone.As<std::int64_t>()[0] = 6;
    m.AddRegion(meshioplusplus::Region("lost", meshioplusplus::RegionKind::Point, gone));
    NDArray cells = NDArray::Uninit(meshioplusplus::DType::Int64, {std::size_t{1}});
    cells.As<std::int64_t>()[0] = 0;
    m.AddRegion(meshioplusplus::Region("cellgroup", meshioplusplus::RegionKind::Cell, cells));

    PointFilterOptions o;
    o.mRadius = 1.5;
    o.mMinNeighbors = 1;
    const auto r = filter_points(m, o);
    EXPECT_EQ(r.mMesh.NumPoints(), 6u);
    EXPECT_EQ(r.mNumRemovedRadius, 1);
    const NDArray& id = r.mMesh.PointData("id");
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(id.As<double>()[i], static_cast<double>(i));
    EXPECT_TRUE(r.mMesh.HasFieldData("meta"));
    EXPECT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "vertex");
    // The point group follows (entry 6 is gone, the group stays); the point-only
    // group that lost its one point stays, empty; the cell group is dropped.
    const std::size_t keep = r.mMesh.FindRegion("keep", meshioplusplus::RegionKind::Point);
    ASSERT_NE(keep, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(keep).NumEntries(), 2u);
    const std::size_t lost = r.mMesh.FindRegion("lost", meshioplusplus::RegionKind::Point);
    ASSERT_NE(lost, Mesh::npos);
    EXPECT_EQ(r.mMesh.Region(lost).NumEntries(), 0u);
    EXPECT_EQ(r.mMesh.FindRegion("cellgroup", meshioplusplus::RegionKind::Cell), Mesh::npos);
    EXPECT_FALSE(r.mWarnings.empty());
}

TEST(FilterPoints, IsDeterministic) {
    PointFilterOptions o;
    o.mStatisticalK = 8;
    o.mRadius = 0.05;
    o.mMinNeighbors = 2;
    o.mVoxelSize = 0.02;
    o.mMinDistance = 0.03;
    const Mesh m = fp_cloud(fp_scan(4000, 20));
    const auto a = filter_points(m, o);
    const auto b = filter_points(m, o);
    EXPECT_EQ(a.mPointMap, b.mPointMap);
}

TEST(FilterPoints, Refusals) {
    const Mesh m = fp_cloud(fp_scan(100, 0));
    EXPECT_THROW(filter_points(m, {}), std::invalid_argument);  // no step
    PointFilterOptions bad;
    bad.mRadius = -1.0;
    EXPECT_THROW(filter_points(m, bad), std::invalid_argument);
    bad = {};
    bad.mVoxelSize = std::nan("");
    EXPECT_THROW(filter_points(m, bad), std::invalid_argument);
    bad = {};
    bad.mStatisticalK = -2;
    EXPECT_THROW(filter_points(m, bad), std::invalid_argument);
    bad = {};
    bad.mVoxelSize = 1e-30;  // far too small for the extent
    EXPECT_THROW(filter_points(m, bad), std::invalid_argument);
    Mesh nan_cloud = fp_cloud({{0, 0, 0}, {std::nan(""), 0, 0}});
    PointFilterOptions ok;
    ok.mMinDistance = 0.1;
    EXPECT_THROW(filter_points(nan_cloud, ok), std::invalid_argument);
    // An empty cloud with a step enabled is an empty result, not an error.
    const Mesh empty = fp_cloud({});
    ok.mMinDistance = 0.1;
    EXPECT_EQ(filter_points(empty, ok).mMesh.NumPoints(), 0u);
}

// ---- PointKnn::CountWithin (the radius query the radius filter uses) --------

#include "../../src/cpp/src/detail/point_knn.hpp"

TEST(PointKnnCountWithin, MatchesABruteForceCountAndHonoursTheLimit) {
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<std::array<double, 3>> pts(800);
    for (auto& p : pts)
        p = {u(rng), u(rng), u(rng)};
    const double radius = 0.3;
    const meshioplusplus::detail::PointKnn knn(pts, radius);
    for (int q = 0; q < 40; ++q) {
        const std::array<double, 3> p{u(rng), u(rng), u(rng)};
        std::size_t brute = 0;
        for (const auto& s : pts) {
            const double dx = p[0] - s[0], dy = p[1] - s[1], dz = p[2] - s[2];
            if (dx * dx + dy * dy + dz * dz <= radius * radius)
                ++brute;
        }
        EXPECT_EQ(knn.CountWithin(p, radius), brute) << q;
        EXPECT_EQ(knn.CountWithin(p, radius, 3), std::min<std::size_t>(brute, 3)) << q;
    }
    EXPECT_EQ(knn.CountWithin(pts[0], 0.0), 1u);  // the point itself, radius 0
    EXPECT_EQ(knn.CountWithin(pts[0], -1.0), 0u);
    EXPECT_EQ(meshioplusplus::detail::PointKnn().CountWithin({0, 0, 0}, 1.0), 0u);
}
