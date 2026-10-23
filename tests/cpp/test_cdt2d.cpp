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
// the exact in-circle test (the constrained Delaunay property). The crossing
// finder is checked case by case and against a test of every pair, and a
// segment inserted after classification against the inside it leaves alone.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
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

// The crossings of a set of segments, in the plane: every pair that meets
// anywhere but at a shared endpoint, by exact orientation tests alone.
std::vector<std::array<std::int64_t, 2>> crossings(
    const std::vector<Vec2>& rPoints, const std::vector<std::array<std::int64_t, 2>>& rSegs) {
    return mio::detail::cdt2d_find_crossings(rPoints, rSegs);
}

TEST(Cdt2dCrossings, ClassifiesEveryWayTwoSegmentsMeet) {
    using Pairs = std::vector<std::array<std::int64_t, 2>>;
    // An X.
    EXPECT_EQ(crossings({{0, 0}, {2, 2}, {0, 2}, {2, 0}}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
    // A T: an endpoint in the other's interior.
    EXPECT_EQ(crossings({{0, 0}, {2, 0}, {1, 0}, {1, 1}}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
    // Collinear and overlapping.
    EXPECT_EQ(crossings({{0, 0}, {2, 0}, {1, 0}, {3, 0}}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
    // Collinear and apart.
    EXPECT_TRUE(crossings({{0, 0}, {1, 0}, {2, 0}, {3, 0}}, {{0, 1}, {2, 3}}).empty());
    // Two distinct points at the same place touch.
    EXPECT_EQ(crossings({{0, 0}, {1, 0}, {1, 0}, {1, 1}}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
    // Adjacent segments of a polyline meet only at their shared point...
    EXPECT_TRUE(crossings({{0, 0}, {1, 0}, {1, 1}}, {{0, 1}, {1, 2}}).empty());
    // ... and straight on is fine too, but folding back is not.
    EXPECT_TRUE(crossings({{0, 0}, {1, 0}, {2, 0}}, {{0, 1}, {1, 2}}).empty());
    EXPECT_EQ(crossings({{0, 0}, {2, 0}, {1, 0}}, {{0, 1}, {1, 2}}), (Pairs{{0, 1}}));
    // The same segment twice, either way round.
    EXPECT_EQ(crossings({{0, 0}, {1, 1}}, {{0, 1}, {1, 0}}), (Pairs{{0, 1}}));
    // A segment of zero length meets nothing.
    EXPECT_TRUE(crossings({{0, 0}, {2, 2}, {1, 1}, {1, 1}}, {{0, 1}, {2, 3}}).empty());
}

TEST(Cdt2dCrossings, NearMissesAreDecidedExactly) {
    // An endpoint one ulp off the other segment, either side: a crossing on
    // one side only. A tolerance would call both (or neither) a crossing.
    const double y = 0.1 + 0.2;  // not 0.3: the segment's line is exact anyway
    const Vec2 below{1.0, std::nextafter(y, 0.0)};
    const Vec2 above{1.0, std::nextafter(y, 1.0)};
    using Pairs = std::vector<std::array<std::int64_t, 2>>;
    // The segment (0, y)-(2, y) and one from (1, -1) up to just below/above it.
    EXPECT_TRUE(crossings({{0, y}, {2, y}, {1, -1}, below}, {{0, 1}, {2, 3}}).empty());
    EXPECT_EQ(crossings({{0, y}, {2, y}, {1, -1}, above}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
    EXPECT_EQ(crossings({{0, y}, {2, y}, {1, -1}, {1.0, y}}, {{0, 1}, {2, 3}}), (Pairs{{0, 1}}));
}

TEST(Cdt2dCrossings, AgreesWithEveryPairOnRandomPolygons) {
    // Random closed polylines (star polygons with noisy radii, some of them
    // crossing): the bucketed search finds exactly the pairs a test of every
    // pair finds.
    std::mt19937_64 rng(20261023);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    std::size_t found = 0;
    for (int trial = 0; trial < 40; ++trial) {
        std::vector<Vec2> pts;
        std::vector<std::array<std::int64_t, 2>> segs;
        const int loops = 1 + trial % 3;
        for (int l = 0; l < loops; ++l) {
            const int n = 5 + static_cast<int>(u01(rng) * 40);
            const double cx = u01(rng) * 3, cy = u01(rng) * 3;
            const auto base = static_cast<std::int64_t>(pts.size());
            for (int i = 0; i < n; ++i) {
                const double a = 2 * M_PI * (i + 0.5 * u01(rng)) / n;
                const double r = 0.5 + u01(rng);
                pts.push_back({cx + r * std::cos(a), cy + r * std::sin(a)});
            }
            for (int i = 0; i < n; ++i)
                segs.push_back({base + i, base + (i + 1) % n});
        }
        std::vector<std::array<std::int64_t, 2>> brute;
        for (std::size_t s = 0; s < segs.size(); ++s)
            for (std::size_t t = s + 1; t < segs.size(); ++t)
                if (!crossings(pts, {segs[s], segs[t]}).empty())
                    brute.push_back({static_cast<std::int64_t>(s), static_cast<std::int64_t>(t)});
        EXPECT_EQ(crossings(pts, segs), brute) << "trial " << trial;
        found += brute.size();
    }
    EXPECT_GT(found, 100u);  // the overlapping loops do cross
}

TEST(Cdt2d, WhereLocatesAPointExactly) {
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    const auto ids = cdt.InsertBatch({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
    for (int k = 0; k < 4; ++k)
        ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    using Place = mio::detail::Cdt2dPlace;
    EXPECT_EQ(cdt.Where({1, 1}).mPlace, Place::Inside);
    EXPECT_EQ(cdt.Where({5, 1}).mPlace, Place::Outside);
    EXPECT_EQ(cdt.Where({2, 0}).mPlace, Place::OnConstraint);
    EXPECT_EQ(cdt.Where({std::nextafter(4.0, 5.0), 2}).mPlace, Place::Outside);
    EXPECT_EQ(cdt.Where({std::nextafter(4.0, 0.0), 2}).mPlace, Place::Inside);
    const auto v = cdt.Where({4, 4});
    EXPECT_EQ(v.mPlace, Place::OnVertex);
    EXPECT_EQ(v.mVertex, ids[2]);
    // Where inserts nothing.
    EXPECT_EQ(cdt.NumVertices(), mio::detail::Cdt2d::kFirstVertex + 4);
}

TEST(Cdt2d, ASegmentAfterClassifyLeavesTheDomainAlone) {
    // A square with a square hole; then, after classification, a dangling
    // segment and a chord between two boundary vertices are constrained.
    // Classification counted the loops alone, so the inside is the same
    // area, and nothing is an open loop.
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    const std::vector<Vec2> outer = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    const std::vector<Vec2> hole = {{1, 1}, {1, 2}, {2, 2}, {2, 1}};
    std::vector<Vec2> all = outer;
    all.insert(all.end(), hole.begin(), hole.end());
    const auto ids = cdt.InsertBatch(all);
    for (int l = 0; l < 2; ++l)
        for (int k = 0; k < 4; ++k)
            ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(4 * l + k)],
                                        ids[std::size_t(4 * l + (k + 1) % 4)]),
                      Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    const auto area = [&] {
        double a = 0.0;
        for (const auto& t : cdt.InsideTriangles()) {
            const Vec2& p = cdt.Point(t[0]);
            const Vec2& q = cdt.Point(t[1]);
            const Vec2& r = cdt.Point(t[2]);
            a += 0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]));
        }
        return a;
    };
    EXPECT_DOUBLE_EQ(area(), 15.0);
    const std::int64_t a = cdt.Insert({2.5, 3.0});
    const std::int64_t b = cdt.Insert({3.5, 0.5});
    ASSERT_GE(a, mio::detail::Cdt2d::kFirstVertex);
    EXPECT_EQ(cdt.InsertSegment(a, b), Cdt2dStatus::Ok);            // dangling, inside
    EXPECT_EQ(cdt.InsertSegment(ids[3], ids[6]), Cdt2dStatus::Ok);  // a chord to the hole
    EXPECT_TRUE(cdt.IsConstrainedEdge(a, b) || cdt.IsConstrainedEdge(b, a));
    EXPECT_TRUE(cdt.IsConstrainedEdge(ids[3], ids[6]) || cdt.IsConstrainedEdge(ids[6], ids[3]));
    EXPECT_DOUBLE_EQ(area(), 15.0);
    for (const auto& t : cdt.InsideTriangles())
        EXPECT_GT(mio::detail::orient2d(cdt.Point(t[0]).data(), cdt.Point(t[1]).data(),
                                        cdt.Point(t[2]).data()),
                  0.0);
}

TEST(Cdt2d, CoordinatesOutsideThePredicatesRangeAreRefused) {
    // The predicates are exact while no product inside them overflows or
    // underflows; incircle takes coordinate differences to the fourth power.
    // A square of side 2^127 still triangulates; coordinates near 1e305 or
    // subnormal ones (found by the triangulate_polygon_2d fuzz target, whose
    // walk then failed with an internal error) are refused by name.
    using mio::detail::cdt2d_in_range;
    EXPECT_TRUE(cdt2d_in_range(0.0));
    EXPECT_TRUE(cdt2d_in_range(-0.0));
    EXPECT_TRUE(cdt2d_in_range(0x1p128));
    EXPECT_TRUE(cdt2d_in_range(-0x1p-128));
    EXPECT_FALSE(cdt2d_in_range(std::nextafter(0x1p128, 1e300)));
    EXPECT_FALSE(cdt2d_in_range(std::nextafter(0x1p-128, 0.0)));
    EXPECT_FALSE(cdt2d_in_range(1e-323));
    const double s = 0x1p127;
    check_triangulation(mio::detail::cdt2d_triangulate_polygon({square(-s, -s, s)}, 0.0), s * s, 0);
    for (const std::vector<Vec2>& bad :
         {std::vector<Vec2>{{5e-299, -8.9e304}, {-8.9e304, 1e-323}, {0.0, 0.0}},
          std::vector<Vec2>{{0, 0}, {1, 0}, {0.5, 1e-300}}}) {
        try {
            (void)mio::detail::cdt2d_triangulate_polygon({bad}, 0.0);
            ADD_FAILURE() << "not refused";
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("exact predicates"), std::string::npos);
        }
    }
}

}  // namespace
