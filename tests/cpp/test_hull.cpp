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
// Convex hull and planar Delaunay triangulation of scattered points (roadmap 7.5.1).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/hull.hpp"
#include "../../src/cpp/src/detail/predicates.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

Mesh hu_cloud(const std::vector<std::array<double, 3>>& rPts) {
    Mesh m;
    NDArray p(DType::Float64, {rPts.size(), std::size_t{3}});
    for (std::size_t i = 0; i < rPts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[i * 3 + k] = rPts[i][k];
    m.AssignPoints(std::move(p));
    return m;
}

std::vector<std::array<double, 3>> hu_random(std::size_t N, unsigned Seed, bool Flat) {
    std::mt19937_64 rng(Seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<std::array<double, 3>> out(N);
    for (auto& q : out)
        q = {u(rng), u(rng), Flat ? 0.0 : u(rng)};
    return out;
}

/// Every edge of a closed, consistently wound triangle surface appears once
/// in each direction.
void hu_expect_closed(const Mesh& rM) {
    std::map<std::pair<std::int64_t, std::int64_t>, int> edges;
    const auto cb = rM.Cells(0);
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < cb.NumCells(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            ++edges[{c[3 * t + k], c[3 * t + (k + 1) % 3]}];
    for (const auto& [e, n] : edges) {
        EXPECT_EQ(n, 1);
        EXPECT_EQ(edges.count({e.second, e.first}), 1u);
    }
}

}  // namespace

TEST(Hull, CubeCornersAndInteriorPoints) {
    std::vector<std::array<double, 3>> pts;
    for (int i = 0; i < 8; ++i)
        pts.push_back({double(i & 1), double((i >> 1) & 1), double((i >> 2) & 1)});
    for (const auto& q : hu_random(200, 3, false))
        pts.push_back({0.5 + 0.4 * q[0], 0.5 + 0.4 * q[1], 0.5 + 0.4 * q[2]});
    pts.push_back({0.5, 0.5, 1.0});  // on a face: not a vertex
    const ConvexHullResult r = convex_hull(hu_cloud(pts));
    EXPECT_EQ(r.mDimension, 3);
    EXPECT_EQ(r.mNumVertices, 8);
    EXPECT_EQ(r.mNumFacets, 12);
    EXPECT_NEAR(r.mVolume, 1.0, 1e-14);
    EXPECT_NEAR(r.mArea, 6.0, 1e-14);
    hu_expect_closed(r.mMesh);
    const std::int64_t* ids = r.mMesh.PointData(kHullPointIdName).As<std::int64_t>();
    for (std::int64_t i = 0; i < 8; ++i)
        EXPECT_EQ(ids[i], i);
}

TEST(Hull, RandomCloudContainsEveryPoint) {
    const auto pts = hu_random(5000, 11, false);
    const ConvexHullResult r = convex_hull(hu_cloud(pts));
    hu_expect_closed(r.mMesh);
    EXPECT_EQ(r.mNumFacets, 2 * r.mNumVertices - 4);  // Euler, all triangles
    const auto cb = r.mMesh.Cells(0);
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    const std::int64_t* ids = r.mMesh.PointData(kHullPointIdName).As<std::int64_t>();
    for (std::size_t t = 0; t < cb.NumCells(); ++t) {
        const double* a = pts[static_cast<std::size_t>(ids[c[3 * t]])].data();
        const double* b = pts[static_cast<std::size_t>(ids[c[3 * t + 1]])].data();
        const double* d = pts[static_cast<std::size_t>(ids[c[3 * t + 2]])].data();
        for (const auto& q : pts)
            ASSERT_GE(detail::orient3d(a, b, d, q.data()), 0.0);  // nothing outside
    }
    EXPECT_GT(r.mVolume, 7.0);
    EXPECT_LT(r.mVolume, 8.0);
}

TEST(Hull, TwoDimensionalLoop) {
    std::vector<std::array<double, 3>> pts{{0, 0, 0}, {2, 0, 0},   {2, 1, 0}, {0, 1, 0},
                                           {1, 0, 0}, {1, 0.5, 0}, {0, 0, 0}};
    const ConvexHullResult r = convex_hull(hu_cloud(pts));
    EXPECT_EQ(r.mDimension, 2);
    EXPECT_EQ(r.mNumVertices, 4);  // (1, 0) is on an edge, (0, 0) repeated
    EXPECT_EQ(r.mNumFacets, 4);
    EXPECT_DOUBLE_EQ(r.mArea, 2.0);
    EXPECT_DOUBLE_EQ(r.mPerimeter, 6.0);
    EXPECT_EQ(r.mMesh.Cells(0).Type(), "line");
}

TEST(Hull, Refusals) {
    EXPECT_THROW(
        convex_hull(hu_cloud({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {3, 3, 3}}), ConvexHullOptions{3}),
        std::invalid_argument);
    EXPECT_THROW(
        convex_hull(hu_cloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}), ConvexHullOptions{3}),
        std::invalid_argument);  // coplanar
    EXPECT_THROW(convex_hull(hu_cloud({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}})),
                 std::invalid_argument);  // collinear
    EXPECT_THROW(convex_hull(hu_cloud({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}), ConvexHullOptions{4}),
                 std::invalid_argument);
    EXPECT_THROW(convex_hull(hu_cloud({{0, 0, 0}, {1, 0, 0}, {0, NAN, 0}})), std::invalid_argument);
}

TEST(Delaunay2d, EmptyCircumcircles) {
    auto pts = hu_random(2000, 5, true);
    pts.push_back(pts[7]);  // a duplicate
    for (auto& q : pts)
        q[2] = q[0] * q[0] - q[1];  // a terrain: z rides along
    const Delaunay2dResult r = delaunay_2d(hu_cloud(pts));
    EXPECT_EQ(r.mNumDuplicates, 1);
    const auto cb = r.mMesh.Cells(0);
    ASSERT_EQ(cb.Type(), "triangle");
    EXPECT_EQ(r.mMesh.NumPoints(), pts.size());
    // Euler for a triangulated convex polygon: T = 2n - h - 2.
    EXPECT_EQ(r.mNumTriangles, 2 * static_cast<std::int64_t>(pts.size() - 1) - r.mNumBoundary - 2);
    const std::int64_t* c = cb.Conn().As<std::int64_t>();
    std::vector<std::array<double, 2>> xy(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i)
        xy[i] = {pts[i][0], pts[i][1]};
    for (std::size_t t = 0; t < cb.NumCells(); ++t) {
        const double* a = xy[static_cast<std::size_t>(c[3 * t])].data();
        const double* b = xy[static_cast<std::size_t>(c[3 * t + 1])].data();
        const double* d = xy[static_cast<std::size_t>(c[3 * t + 2])].data();
        ASSERT_GT(detail::orient2d(a, b, d), 0.0);
        for (const auto& q : xy)
            ASSERT_LE(detail::incircle(a, b, d, q.data()), 0.0);
    }
    EXPECT_EQ(r.mMesh.Points().As<double>()[3 * 4 + 2], pts[4][2]);
}

TEST(Delaunay2d, GridIsCocircularButValid) {
    std::vector<std::array<double, 3>> pts;
    for (int j = 0; j < 5; ++j)
        for (int i = 0; i < 5; ++i)
            pts.push_back({double(i), double(j), 0.0});
    const Delaunay2dResult a = delaunay_2d(hu_cloud(pts));
    const Delaunay2dResult b = delaunay_2d(hu_cloud(pts));
    EXPECT_EQ(a.mNumTriangles, 32);
    EXPECT_EQ(a.mNumBoundary, 16);
    const auto ca = a.mMesh.Cells(0).Conn(), cb = b.mMesh.Cells(0).Conn();
    EXPECT_TRUE(std::equal(ca.As<std::int64_t>(), ca.As<std::int64_t>() + 96,
                           cb.As<std::int64_t>()));  // deterministic
}

TEST(Delaunay2d, Refusals) {
    EXPECT_THROW(delaunay_2d(hu_cloud({{0, 0, 0}, {1, 1, 0}, {2, 2, 0}})), std::invalid_argument);
    EXPECT_THROW(delaunay_2d(hu_cloud({{0, 0, 0}, {0, 0, 1}, {1, 0, 0}})), std::invalid_argument);
    EXPECT_THROW(delaunay_2d(hu_cloud({{0, 0, 0}, {1e300, 0, 0}, {0, 1, 0}})),
                 std::invalid_argument);
}

TEST(Hull, LatticeIsAllCoplanarFaces) {
    // Every outer lattice point lies on a face plane: only the corners are
    // extreme, whatever order the coplanar points arrive in.
    std::vector<std::array<double, 3>> pts;
    for (int k = 0; k < 6; ++k)
        for (int j = 0; j < 6; ++j)
            for (int i = 0; i < 6; ++i)
                pts.push_back({double(i), double(j), double(k)});
    std::shuffle(pts.begin(), pts.end(), std::mt19937_64(4));
    const ConvexHullResult r = convex_hull(hu_cloud(pts));
    EXPECT_EQ(r.mNumVertices, 8);
    EXPECT_EQ(r.mNumFacets, 12);
    EXPECT_DOUBLE_EQ(r.mVolume, 125.0);
    EXPECT_DOUBLE_EQ(r.mArea, 150.0);
    hu_expect_closed(r.mMesh);
}
