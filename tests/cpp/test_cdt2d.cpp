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
// detail/cdt2d: constrained Delaunay triangulations of polygons with holes.
// Every result is checked for the same invariants: the area equals the
// analytic one, every triangle is counter-clockwise with positive area (by
// the exact orient2d), every constrained boundary segment is an edge and the
// boundary edges of the triangulation are exactly those segments, the Euler
// characteristic V - E + F is 1 - holes, and every unconstrained edge passes
// the exact in-circle test (the constrained Delaunay property).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/cdt2d.hpp"
#include "../../src/cpp/src/detail/predicates.hpp"

namespace {

namespace mio = meshioplusplus;
using mio::detail::Cdt2dPolygon;
using mio::detail::Cdt2dStatus;
using Vec2 = mio::detail::Cdt2d::Vec2;
using Loops = std::vector<std::vector<Vec2>>;

double polygon_area(const std::vector<Vec2>& rLoop) {
    double a = 0.0;
    for (std::size_t i = 0; i < rLoop.size(); ++i) {
        const Vec2& p = rLoop[i];
        const Vec2& q = rLoop[(i + 1) % rLoop.size()];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

struct Checked {
    double mArea = 0.0;
    double mMinEdge = 0.0;
    double mMaxEdge = 0.0;
    double mMinAngleDeg = 0.0;
    std::size_t mEdges = 0;
};

// Check every invariant of the file comment; `Holes` for the Euler test.
Checked check_triangulation(const Cdt2dPolygon& rT, double ExpectedArea, int Holes) {
    Checked c;
    EXPECT_EQ(rT.mStatus, Cdt2dStatus::Ok) << rT.mWhy;
    const auto& pts = rT.mPoints;
    using Edge = std::pair<std::int64_t, std::int64_t>;
    std::map<Edge, std::vector<std::pair<std::size_t, int>>>
        edges;  // undirected -> (tri, opposite)
    std::set<std::int64_t> used;
    c.mMinEdge = 1e300;
    c.mMinAngleDeg = 180.0;
    double comp = 0.0;
    for (std::size_t t = 0; t < rT.mTriangles.size(); ++t) {
        const auto& tri = rT.mTriangles[t];
        const Vec2& a = pts[tri[0]];
        const Vec2& b = pts[tri[1]];
        const Vec2& cc = pts[tri[2]];
        EXPECT_GT(mio::detail::orient2d(a.data(), b.data(), cc.data()), 0.0) << "triangle " << t;
        // Neumaier's compensated sum: tens of thousands of areas stay within
        // 1e-12 of the analytic total.
        const double term = 0.5 * ((b[0] - a[0]) * (cc[1] - a[1]) - (b[1] - a[1]) * (cc[0] - a[0]));
        const double sum = c.mArea + term;
        comp +=
            std::fabs(c.mArea) >= std::fabs(term) ? (c.mArea - sum) + term : (term - sum) + c.mArea;
        c.mArea = sum;
        for (int k = 0; k < 3; ++k) {
            const std::int64_t u = tri[(k + 1) % 3];
            const std::int64_t v = tri[(k + 2) % 3];
            edges[{std::min(u, v), std::max(u, v)}].emplace_back(t, k);
            used.insert(tri[k]);
            const Vec2& p = pts[u];
            const Vec2& q = pts[v];
            const double len = std::hypot(q[0] - p[0], q[1] - p[1]);
            c.mMinEdge = std::min(c.mMinEdge, len);
            c.mMaxEdge = std::max(c.mMaxEdge, len);
            // The angle at tri[k].
            const Vec2& o = pts[tri[k]];
            const double ux = p[0] - o[0], uy = p[1] - o[1], vx = q[0] - o[0], vy = q[1] - o[1];
            const double ang = std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy);
            c.mMinAngleDeg = std::min(c.mMinAngleDeg, ang * 180.0 / M_PI);
        }
    }
    c.mArea += comp;
    c.mEdges = edges.size();
    EXPECT_NEAR(c.mArea, ExpectedArea, 1e-12 * std::fabs(ExpectedArea));
    std::set<Edge> segs;
    for (const auto& s : rT.mSegments)
        segs.insert({std::min(s[0], s[1]), std::max(s[0], s[1])});
    for (const Edge& s : segs)
        EXPECT_TRUE(edges.count(s)) << "segment " << s.first << "-" << s.second << " missing";
    for (const auto& [e, uses] : edges) {
        EXPECT_LE(uses.size(), 2u);
        if (uses.size() == 1) {
            EXPECT_TRUE(segs.count(e)) << "open edge " << e.first << "-" << e.second;
            continue;
        }
        if (segs.count(e))
            continue;
        // Unconstrained interior edge: the opposite vertex of one triangle is
        // not inside the circumcircle of the other.
        const auto& t0 = rT.mTriangles[uses[0].first];
        const std::int64_t d = rT.mTriangles[uses[1].first][uses[1].second];
        EXPECT_LE(mio::detail::incircle(pts[t0[0]].data(), pts[t0[1]].data(), pts[t0[2]].data(),
                                        pts[d].data()),
                  0.0)
            << "edge " << e.first << "-" << e.second << " is not Delaunay";
    }
    const auto v = static_cast<std::int64_t>(used.size());
    const auto e = static_cast<std::int64_t>(edges.size());
    const auto f = static_cast<std::int64_t>(rT.mTriangles.size());
    EXPECT_EQ(v - e + f, 1 - Holes);
    return c;
}

std::vector<Vec2> square(double X0, double Y0, double S, bool Ccw = true) {
    std::vector<Vec2> p = {{X0, Y0}, {X0 + S, Y0}, {X0 + S, Y0 + S}, {X0, Y0 + S}};
    if (!Ccw)
        std::reverse(p.begin(), p.end());
    return p;
}

TEST(Cdt2d, UnitSquare) {
    const auto t = mio::detail::cdt2d_triangulate_polygon({square(0, 0, 1)}, 0.0);
    check_triangulation(t, 1.0, 0);
    EXPECT_EQ(t.mTriangles.size(), 2u);
}

TEST(Cdt2d, LShape) {
    const Loops l = {{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}};
    const auto t = mio::detail::cdt2d_triangulate_polygon(l, 0.0);
    check_triangulation(t, 3.0, 0);
    EXPECT_EQ(t.mTriangles.size(), 4u);
}

TEST(Cdt2d, SquareWithSquareHoleEitherOrientation) {
    for (const bool ccw : {true, false}) {
        const auto t =
            mio::detail::cdt2d_triangulate_polygon({square(0, 0, 4), square(1, 1, 2, ccw)}, 0.0);
        check_triangulation(t, 12.0, 1);
        EXPECT_EQ(t.mTriangles.size(), 8u);
    }
}

TEST(Cdt2d, TwoHolesAndRefinement) {
    const Loops l = {square(0, 0, 10), square(1, 1, 2), square(5, 6, 3, false)};
    check_triangulation(mio::detail::cdt2d_triangulate_polygon(l, 0.0), 100.0 - 4.0 - 9.0, 2);
    const auto fine = mio::detail::cdt2d_triangulate_polygon(l, 0.5);
    const Checked c = check_triangulation(fine, 100.0 - 4.0 - 9.0, 2);
    EXPECT_LE(c.mMaxEdge, 1.5 * 0.5);
    EXPECT_GE(c.mMinEdge, 0.3 * 0.5);
    EXPECT_GT(c.mMinAngleDeg, 10.0);
}

TEST(Cdt2d, CollinearAndNearDuplicatePoints) {
    // Collinear points along every side, a near-duplicate corner and an exact
    // duplicate, in a square of side 3.
    const Loops l = {{{0, 0},
                      {1, 0},
                      {2, 0},
                      {3, 0},
                      {3, 0},
                      {3, 1},
                      {3, 2},
                      {3, 3 + 1e-14},
                      {3, 3},
                      {1.5, 3},
                      {0, 3},
                      {0, 1.5},
                      {1e-15, 1e-15}}};
    const auto t = mio::detail::cdt2d_triangulate_polygon(l, 0.0);
    // The near-duplicates merge into the first of each pair: (3, 3 + 1e-14)
    // stays and (3, 3) folds into it, (1e-15, 1e-15) into (0, 0).
    const std::vector<Vec2> merged = {{0, 0}, {1, 0},         {2, 0},   {3, 0}, {3, 1},
                                      {3, 2}, {3, 3 + 1e-14}, {1.5, 3}, {0, 3}, {0, 1.5}};
    check_triangulation(t, polygon_area(merged), 0);
    EXPECT_EQ(t.mNumBoundary, 10);
}

TEST(Cdt2d, ThinSliver) {
    for (const double w : {1e-3, 1e-6, 1e-9}) {
        std::vector<Vec2> loop;
        for (int i = 0; i <= 20; ++i)
            loop.push_back({double(i) / 20.0, 0.0});
        loop.push_back({1.0, w});
        for (int i = 19; i >= 1; --i)
            loop.push_back({double(i) / 20.0, w});
        loop.push_back({0.0, w});
        const auto t = mio::detail::cdt2d_triangulate_polygon({loop}, 0.0);
        check_triangulation(t, w, 0);
    }
}

TEST(Cdt2d, SelfIntersectingIsRefused) {
    const Loops bowtie = {{{0, 0}, {1, 1}, {1, 0}, {0, 1}}};
    EXPECT_EQ(mio::detail::cdt2d_triangulate_polygon(bowtie, 0.0).mStatus,
              Cdt2dStatus::SelfIntersecting);
    // Two loops that cross.
    const Loops crossing = {square(0, 0, 2), square(1, 1, 2)};
    EXPECT_EQ(mio::detail::cdt2d_triangulate_polygon(crossing, 0.0).mStatus,
              Cdt2dStatus::SelfIntersecting);
}

TEST(Cdt2d, DegenerateIsRefused) {
    EXPECT_EQ(mio::detail::cdt2d_triangulate_polygon({{{0, 0}, {1, 0}}}, 0.0).mStatus,
              Cdt2dStatus::Degenerate);
    EXPECT_EQ(mio::detail::cdt2d_triangulate_polygon({{{0, 0}, {1, 0}, {2, 0}}}, 0.0).mStatus,
              Cdt2dStatus::Degenerate);
}

TEST(Cdt2d, RandomStarPolygonsWithHoles) {
    std::mt19937_64 rng(29);
    std::uniform_real_distribution<double> radius(1.0, 2.0);
    for (int it = 0; it < 40; ++it) {
        const int n = 8 + it * 7;
        std::vector<Vec2> outer;
        for (int i = 0; i < n; ++i) {
            const double a = 2.0 * M_PI * double(i) / double(n);
            const double r = radius(rng);
            outer.push_back({r * std::cos(a), r * std::sin(a)});
        }
        Loops l = {outer};
        double area = polygon_area(outer);
        // A hole inside the unit disc (inside every outer radius).
        std::vector<Vec2> hole;
        for (int i = 0; i < 6; ++i) {
            const double a = 2.0 * M_PI * double(i) / 6.0;
            hole.push_back({0.5 * std::cos(a), 0.5 * std::sin(a)});
        }
        l.push_back(hole);
        area -= polygon_area(hole);
        check_triangulation(mio::detail::cdt2d_triangulate_polygon(l, 0.0), area, 1);
        check_triangulation(mio::detail::cdt2d_triangulate_polygon(l, 0.15), area, 1);
    }
}

TEST(Cdt2d, IncrementalApi) {
    // Points on a lattice (many exactly cocircular quadruples), then one
    // segment across it; the in-circle ties must not break anything.
    mio::detail::Cdt2d cdt(0, 0, 8, 8);
    std::vector<Vec2> pts;
    for (int j = 0; j <= 8; ++j)
        for (int i = 0; i <= 8; ++i)
            pts.push_back({double(i), double(j)});
    const auto ids = cdt.InsertBatch(pts);
    EXPECT_EQ(cdt.InsertBatch(pts), ids);  // exact duplicates map to the same vertices
    EXPECT_EQ(cdt.InsertSegment(ids[0], ids[80]),
              Cdt2dStatus::Ok);  // the diagonal, through 7 vertices
    for (int k = 0; k < 8; ++k)
        EXPECT_TRUE(
            cdt.IsConstrainedEdge(ids[std::size_t(k * 10)], ids[std::size_t(k * 10 + 10)]) ||
            cdt.IsConstrainedEdge(ids[std::size_t(k * 10 + 10)], ids[std::size_t(k * 10)]));
    // Through lattice vertices only, touching the diagonal at (4, 4): fine.
    EXPECT_EQ(cdt.InsertSegment(ids[8], ids[72]), Cdt2dStatus::Ok);
    // Across the constrained edge (0, 0)-(1, 1) at (0.5, 0.5): refused.
    EXPECT_EQ(cdt.InsertSegment(ids[1], ids[9]), Cdt2dStatus::SelfIntersecting);
}

TEST(Cdt2d, DanglingSegmentIsAnOpenLoop) {
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    const auto ids = cdt.InsertBatch({{0, 0}, {4, 0}, {4, 4}, {0, 4}, {1, 1}, {2, 3}});
    for (int k = 0; k < 4; ++k)
        ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.InsertSegment(ids[4], ids[5]), Cdt2dStatus::Ok);
    EXPECT_EQ(cdt.Classify(), Cdt2dStatus::OpenLoop);
}

TEST(Cdt2d, LargeRefinement) {
    // About 46 000 interior points: the walk and the flips stay fast, and
    // the lattice keeps the edges within [0.3 h, 1.5 h].
    const double h = 0.005;
    const auto t = mio::detail::cdt2d_triangulate_polygon({square(0, 0, 1)}, h);
    const Checked c = check_triangulation(t, 1.0, 0);
    EXPECT_GT(t.mPoints.size(), 40000u);
    EXPECT_LE(c.mMaxEdge, 1.5 * h);
    EXPECT_GE(c.mMinEdge, 0.3 * h);
    EXPECT_GT(c.mMinAngleDeg, 10.0);
}

}  // namespace
