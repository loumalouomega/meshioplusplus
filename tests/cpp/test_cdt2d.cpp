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
// Refinement is checked for its bound above 60-degree input angles, its
// termination on adversarial domains, the small input angles it reports and
// leaves, the order it takes triangles in, the splits it asks a caller for,
// and the guarded vertex moves smoothing uses.

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

// --- Delaunay refinement (Ruppert 1995, Chew's circumcentres) -------------

mio::detail::Cdt2dPolygon refined(const Loops& rLoops, double H, double MinAngle,
                                  std::int64_t MaxPoints = 200000) {
    mio::detail::Cdt2dPolygonQuality q;
    q.mMinAngleDeg = MinAngle;
    q.mMaxPoints = MaxPoints;
    return mio::detail::cdt2d_triangulate_polygon(rLoops, H, -1.0, q);
}

double triangle_min_angle(const mio::detail::Cdt2dPolygon& rT, std::size_t T) {
    double m = 180.0;
    const auto& tri = rT.mTriangles[T];
    for (int k = 0; k < 3; ++k) {
        const Vec2& o = rT.mPoints[tri[k]];
        const Vec2& p = rT.mPoints[tri[(k + 1) % 3]];
        const Vec2& q = rT.mPoints[tri[(k + 2) % 3]];
        const double ux = p[0] - o[0], uy = p[1] - o[1], vx = q[0] - o[0], vy = q[1] - o[1];
        m = std::min(m, std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / M_PI);
    }
    return m;
}

std::vector<Vec2> regular(double Cx, double Cy, double R, int N, double A0 = 0.0) {
    std::vector<Vec2> p;
    for (int i = 0; i < N; ++i) {
        const double a = A0 + 2.0 * M_PI * double(i) / double(N);
        p.push_back({Cx + R * std::cos(a), Cy + R * std::sin(a)});
    }
    return p;
}

TEST(Cdt2dRefine, ZeroMinAngleKeepsTheTriangulation) {
    const Loops l = {square(0, 0, 10), square(1, 1, 2), square(5, 6, 3, false)};
    const auto a = mio::detail::cdt2d_triangulate_polygon(l, 0.5);
    const auto b = refined(l, 0.5, 0.0);
    EXPECT_EQ(a.mPoints, b.mPoints);
    EXPECT_EQ(a.mTriangles, b.mTriangles);
    EXPECT_EQ(a.mSegments, b.mSegments);
}

TEST(Cdt2dRefine, InputAnglesOfSixtyDegreesOrMoreMeetTheBound) {
    // Ruppert's guarantee at B = sqrt(2) (20.7 degrees) for input angles of
    // at least 60 degrees; 25 and 30 are what it reaches in practice here.
    const std::vector<std::pair<Loops, std::pair<double, int>>> domains = {
        {{square(0, 0, 1)}, {0.1, 0}},
        {{{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}}, {0.2, 0}},
        {{square(0, 0, 10), square(1, 1, 2), square(5, 6, 3, false)}, {0.5, 2}},
        {{regular(0, 0, 3, 40), regular(0, 0, 1, 6)}, {0.4, 1}},
        {{regular(0, 0, 1, 5)}, {0.0, 0}},
        {{regular(0, 0, 1, 6), regular(0.2, 0.1, 0.3, 5)}, {0.0, 1}},
    };
    for (const double ang : {20.7, 25.0, 30.0}) {
        for (const auto& [l, hh] : domains) {
            const auto t = refined(l, hh.first, ang);
            double area = 0.0;
            for (std::size_t k = 0; k < l.size(); ++k)
                area += (k == 0 ? 1.0 : -1.0) * std::fabs(polygon_area(l[k]));
            const Checked c = check_triangulation(t, area, hh.second);
            EXPECT_GE(c.mMinAngleDeg, ang - 1e-9) << "min_angle " << ang << ", h " << hh.first;
            EXPECT_TRUE(t.mRefine.mSmallAngles.empty());
            EXPECT_EQ(t.mRefine.mNumUnresolved, 0);
            EXPECT_EQ(t.mRefine.mNumForced, 0);
            EXPECT_FALSE(t.mRefine.mBudgetReached);
        }
    }
}

TEST(Cdt2dRefine, AdversarialDomainsTerminateUnderMaxPoints) {
    // Slivers: a strip of width w along 20 boundary pieces; refinement fills
    // it to the bound, or stops at max_points with a valid triangulation.
    for (const double w : {1e-2, 1e-3, 1e-6}) {
        std::vector<Vec2> loop;
        for (int i = 0; i <= 20; ++i)
            loop.push_back({double(i) / 20.0, 0.0});
        loop.push_back({1.0, w});
        for (int i = 19; i >= 1; --i)
            loop.push_back({double(i) / 20.0, w});
        loop.push_back({0.0, w});
        const auto t = refined({loop}, 0.0, 25.0, 5000);
        const Checked c = check_triangulation(t, w, 0);
        if (w >= 1e-3) {
            EXPECT_FALSE(t.mRefine.mBudgetReached);
            EXPECT_GE(c.mMinAngleDeg, 25.0 - 1e-9);
        } else {
            EXPECT_TRUE(t.mRefine.mBudgetReached);
            EXPECT_LE(t.mPoints.size(), 42u + 5000u);
        }
    }
    // Near-collinear holes: a needle hole whose apex angles are tiny, in a
    // square; and random star polygons, whose corners are often sharp.
    {
        const Loops l = {square(0, 0, 4), {{1, 1}, {3, 1.001}, {1, 1.002}}};
        const auto t = refined(l, 0.25, 25.0);
        check_triangulation(t, 16.0 - std::fabs(polygon_area(l[1])), 1);
        EXPECT_FALSE(t.mRefine.mBudgetReached);
        EXPECT_EQ(t.mRefine.mNumUnresolved, 0);
        // The needle's tips are reflex corners of the domain, not small ones.
        EXPECT_TRUE(t.mRefine.mSmallAngles.empty());
    }
    std::mt19937_64 rng(31);
    std::uniform_real_distribution<double> radius(1.0, 2.0);
    for (int it = 0; it < 40; ++it) {
        const int n = 8 + it * 7;
        std::vector<Vec2> outer;
        for (int i = 0; i < n; ++i) {
            const double a = 2.0 * M_PI * double(i) / double(n);
            const double r = radius(rng);
            outer.push_back({r * std::cos(a), r * std::sin(a)});
        }
        const Loops l = {outer, regular(0, 0, 0.5, 6)};
        const double area = polygon_area(outer) - std::fabs(polygon_area(l[1]));
        for (const double h : {0.0, 0.15}) {
            const auto t = refined(l, h, 25.0);
            check_triangulation(t, area, 1);
            EXPECT_FALSE(t.mRefine.mBudgetReached);
            // Every triangle left below the bound is one a small input
            // angle forces.
            EXPECT_EQ(t.mRefine.mNumUnresolved, 0) << "polygon " << it << ", h " << h;
            std::size_t bad = 0;
            for (std::size_t k = 0; k < t.mTriangles.size(); ++k)
                bad += triangle_min_angle(t, k) < 25.0 - 1e-9;
            EXPECT_EQ(std::int64_t(bad), t.mRefine.mNumForced);
            if (bad > 0)
                EXPECT_FALSE(t.mRefine.mSmallAngles.empty());
        }
    }
}

TEST(Cdt2dRefine, SmallInputAnglesAreReportedAndTheirTrianglesLeft) {
    // A triangle with a sharp apex at the origin and a quadrilateral with
    // two sharp corners. Below 60 degrees the corner is reported with its
    // angle; every triangle under the bound holds the apex or has its
    // shortest edge between the apex's two segments (the angle forces it),
    // and the smallest angle those keep is reported per corner.
    for (const double phi : {45.0, 30.0, 15.0, 5.0}) {
        const double r = phi * M_PI / 180.0;
        const Loops l = {{{0, 0}, {2, 0}, {2, 2}, {1.8 * std::cos(r), 1.8 * std::sin(r)}}};
        const auto t = refined(l, 0.0, 25.0);
        check_triangulation(t, std::fabs(polygon_area(l[0])), 0);
        EXPECT_FALSE(t.mRefine.mBudgetReached);
        EXPECT_EQ(t.mRefine.mNumUnresolved, 0);
        bool found = false;
        double forced_min = 180.0;
        for (const auto& r_a : t.mRefine.mSmallAngles) {
            if (r_a.mVertex == 0) {
                found = true;
                EXPECT_NEAR(r_a.mAngleDeg, phi, 1e-9);
            }
            forced_min = std::min(forced_min, r_a.mMinForcedDeg);
        }
        EXPECT_TRUE(found);
        const Vec2 apex{0, 0};
        const Vec2 dir[2] = {{1, 0}, {std::cos(r), std::sin(r)}};
        const auto on_segment = [&](const Vec2& rP) {
            int which = -1;
            for (int k = 0; k < 2; ++k)
                if (std::fabs(rP[0] * dir[k][1] - rP[1] * dir[k][0]) <= 1e-12 &&
                    rP[0] * dir[k][0] + rP[1] * dir[k][1] > 0.0)
                    which = k;
            return which;
        };
        double seen_min = 180.0;
        for (std::size_t k = 0; k < t.mTriangles.size(); ++k) {
            const double m = triangle_min_angle(t, k);
            if (m >= 25.0 - 1e-9)
                continue;
            seen_min = std::min(seen_min, m);
            const auto& tri = t.mTriangles[k];
            bool holds = false;
            int sides = 0;
            for (int j = 0; j < 3; ++j) {
                const Vec2& p = t.mPoints[tri[j]];
                holds = holds || p == apex;
                const int w = on_segment(p);
                sides |= w >= 0 ? (1 << w) : 0;
            }
            // The other sharp corner (2, 2)-(1.8 cos, 1.8 sin)-... is not at
            // the origin: those triangles are attributed to it instead.
            if (!holds && sides != 3)
                EXPECT_EQ(t.mRefine.mSmallAngles.size(), 2u) << "phi " << phi;
        }
        if (seen_min < 180.0)
            EXPECT_NEAR(forced_min, seen_min, 1e-9);
    }
}

TEST(Cdt2dRefine, TheOrderDependsOnTheInputAlone) {
    // The worst triangle first, ties broken by sorted vertex ids, never a
    // slot: two runs give the same points and triangles.
    const Loops l = {regular(0, 0, 3, 40), regular(0.3, -0.2, 1, 7)};
    const auto a = refined(l, 0.4, 30.0);
    const auto b = refined(l, 0.4, 30.0);
    EXPECT_EQ(a.mPoints, b.mPoints);
    EXPECT_EQ(a.mTriangles, b.mTriangles);
    EXPECT_GT(a.mRefine.mNumCircumcentres, 0);
}

TEST(Cdt2dRefine, AskingForSplitsLeavesTheBoundaryToTheCaller) {
    // A square whose bottom side has one piece and an interior vertex just
    // above its middle, inside its diametral circle: asked, not split.
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    const auto ids = cdt.InsertBatch({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
    for (int k = 0; k < 4; ++k)
        ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    ASSERT_GE(cdt.InsertSteiner({2.0, 0.5}), 0);
    const auto before = cdt.ConstrainedEdges();
    mio::detail::Cdt2dRefineOptions o;
    o.mMinAngleDeg = 25.0;
    o.mSplitSegments = false;
    const auto rep = cdt.Refine(o);
    EXPECT_EQ(cdt.ConstrainedEdges(), before);
    EXPECT_EQ(rep.mNumSplits, 0);
    const auto bottom = std::minmax(ids[0], ids[1]);
    bool asked = false;
    for (const auto& r_q : rep.mRequests) {
        EXPECT_EQ(r_q.mApex, -1);
        asked = asked || (r_q.mA == bottom.first && r_q.mB == bottom.second);
    }
    EXPECT_TRUE(asked);
    // The caller's own measure: exempting every triangle refines nothing
    // (an acute needle, so no vertex encroaches a side).
    mio::detail::Cdt2d cdt2(0, 0, 4, 4);
    const auto ids2 = cdt2.InsertBatch({{0, 0}, {0.7, 0}, {0.35, 4}});
    for (int k = 0; k < 3; ++k)
        ASSERT_EQ(cdt2.InsertSegment(ids2[std::size_t(k)], ids2[std::size_t((k + 1) % 3)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt2.Classify(), Cdt2dStatus::Ok);
    o.mSplitSegments = true;
    o.mRatio = [](const mio::detail::Cdt2d::Tri&) { return -1.0; };
    const auto rep2 = cdt2.Refine(o);
    EXPECT_EQ(rep2.mNumCircumcentres + rep2.mNumSplits, 0);
    EXPECT_EQ(rep2.mNumForced, 1);
}

TEST(Cdt2dRefine, ASplitTheTriangulationRefusesDeclinesItsTriangle) {
    // Found by the fuzz target (tests/fuzz/regressions/op_triangulate_polygon_2d):
    // a sliver whose corners lie a few ulps off the line y = x / 3. The
    // rounded midpoint of the long side makes the triangles beside it turn
    // clockwise, so the split is refused; the triangle was queued again and
    // refused again forever. It is now declined and counted unresolved.
    const Loops l = {
        {{-0.0625, -0.020833333333333332}, {0.125, 0.041666666666666664}, {-5.8125, -1.9375}}};
    for (const double ang : {20.7, 25.0, 30.0}) {
        const auto t = refined(l, 0.0, ang, 100);
        ASSERT_EQ(t.mStatus, mio::detail::Cdt2dStatus::Ok) << t.mWhy;
        EXPECT_EQ(t.mPoints.size(), 3u) << ang;
        EXPECT_EQ(t.mTriangles.size(), 1u) << ang;
        EXPECT_FALSE(t.mRefine.mBudgetReached) << ang;
        EXPECT_GE(t.mRefine.mNumUnresolved + t.mRefine.mNumForced, 1) << ang;
    }
}

TEST(Cdt2dRefine, FlipEdgeNeedsAConvexQuadrilateralAndNoConstraint) {
    // A caller judging triangles by its own measure flips edges itself: only
    // an unconstrained edge between two triangles forming a strictly convex
    // quadrilateral (exact orient2d).
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    const auto ids = cdt.InsertBatch({{0, 0}, {4, 0}, {4, 4}, {0, 4}, {3, 1}});
    for (int k = 0; k < 4; ++k)
        ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    EXPECT_FALSE(cdt.FlipEdge(ids[0], ids[1]));  // constrained
    // Every interior edge of the star of (3, 1): flipping one whose
    // quadrilateral is convex keeps four triangles of total area 16.
    int flipped = 0;
    for (const std::int64_t w : cdt.Neighbours(ids[4])) {
        if (!cdt.FlipEdge(ids[4], w))
            continue;
        ++flipped;
        double area = 0.0;
        for (const auto& t : cdt.InsideTriangles()) {
            const double o = mio::detail::orient2d(cdt.Point(t[0]).data(), cdt.Point(t[1]).data(),
                                                   cdt.Point(t[2]).data());
            EXPECT_GT(o, 0.0);
            area += 0.5 * o;
        }
        EXPECT_NEAR(area, 16.0, 1e-12);
        EXPECT_EQ(cdt.InsideTriangles().size(), 4u);
        break;
    }
    EXPECT_EQ(flipped, 1);
    // (3, 1) is now joined to three corners: the edge to the fourth is gone.
    EXPECT_EQ(cdt.Neighbours(ids[4]).size(), 3u);
}

TEST(Cdt2dRefine, MoveVertexRefusesAnInversionAndKeepsDelaunay) {
    mio::detail::Cdt2d cdt(0, 0, 4, 4);
    std::vector<Vec2> pts = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    for (int j = 1; j < 4; ++j)
        for (int i = 1; i < 4; ++i)
            pts.push_back({double(i), double(j)});
    const auto ids = cdt.InsertBatch(pts);
    for (int k = 0; k < 4; ++k)
        ASSERT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    const std::int64_t centre = ids[8];  // (2, 2)
    EXPECT_FALSE(cdt.OnConstraint(centre));
    EXPECT_TRUE(cdt.OnConstraint(ids[0]));
    EXPECT_EQ(cdt.Neighbours(centre).size(), cdt.Star(centre).size());
    // Out of its star (past a neighbour): refused, nothing moves.
    EXPECT_FALSE(cdt.MoveVertex(centre, {3.5, 2.0}));
    EXPECT_EQ(cdt.Point(centre), (Vec2{2.0, 2.0}));
    // A boundary vertex never moves.
    EXPECT_FALSE(cdt.MoveVertex(ids[0], {0.1, 0.1}));
    // Within its star: moved, then flipped back to Delaunay.
    EXPECT_TRUE(cdt.MoveVertex(centre, {2.6, 2.3}));
    EXPECT_EQ(cdt.Point(centre), (Vec2{2.6, 2.3}));
    double area = 0.0;
    const auto tris = cdt.InsideTriangles();
    for (const auto& t : tris) {
        const double o = mio::detail::orient2d(cdt.Point(t[0]).data(), cdt.Point(t[1]).data(),
                                               cdt.Point(t[2]).data());
        EXPECT_GT(o, 0.0);
        area += 0.5 * o;
        for (std::int64_t v = mio::detail::Cdt2d::kFirstVertex; v < cdt.NumVertices(); ++v)
            if (v != t[0] && v != t[1] && v != t[2])
                EXPECT_LE(mio::detail::incircle(cdt.Point(t[0]).data(), cdt.Point(t[1]).data(),
                                                cdt.Point(t[2]).data(), cdt.Point(v).data()),
                          0.0);
    }
    EXPECT_NEAR(area, 16.0, 1e-12);
}

}  // namespace

namespace {

// A square of random points, its boundary constrained and classified;
// returns the inside triangles. With a metric, set
// before the first point.
std::vector<mio::detail::Cdt2d::Tri> metric_square(const mio::detail::Cdt2d::Metric& rM,
                                                   std::vector<Vec2>& rPoints) {
    mio::detail::Cdt2d cdt(0, 0, 1, 1);
    if (rM)
        cdt.SetMetric(rM);
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(0.02, 0.98);
    std::vector<Vec2> pts = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    for (int k = 0; k < 200; ++k)
        pts.push_back({u(rng), u(rng)});
    const auto ids = cdt.InsertBatch(pts);
    for (int k = 0; k < 4; ++k)
        EXPECT_EQ(cdt.InsertSegment(ids[std::size_t(k)], ids[std::size_t((k + 1) % 4)]),
                  Cdt2dStatus::Ok);
    EXPECT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    rPoints.clear();
    for (std::int64_t v = 0; v < cdt.NumVertices(); ++v)
        rPoints.push_back(cdt.Point(v));
    return cdt.InsideTriangles();
}

}  // namespace

// Roadmap 6.4.1: with a metric, Lawson's flips decide by the in-circumellipse
// test. The identity metric is the plane's test; a stretched one still
// triangulates the domain validly, its triangles rounder in that metric
// than the plane's Delaunay triangles are; an empty one restores the plane.
TEST(Cdt2d, AMetricDecidesTheFlips) {
    std::vector<Vec2> plane_pts, id_pts, aniso_pts;
    const auto plane = metric_square(nullptr, plane_pts);
    const auto identity =
        metric_square([](const Vec2&) { return std::array<double, 3>{1.0, 0.0, 1.0}; }, id_pts);
    EXPECT_EQ(identity, plane);
    EXPECT_EQ(id_pts, plane_pts);
    // Lengths along y count tenfold, more so towards x = 1.
    const auto stretch = [](const Vec2& rX) {
        const double s = 10.0 * (1.0 + rX[0]);
        return std::array<double, 3>{1.0, 0.0, s * s};
    };
    const auto aniso = metric_square(stretch, aniso_pts);
    double area = 0.0;
    for (const auto& t : aniso) {
        const double a = mio::detail::orient2d(aniso_pts[std::size_t(t[0])].data(),
                                               aniso_pts[std::size_t(t[1])].data(),
                                               aniso_pts[std::size_t(t[2])].data());
        EXPECT_GT(a, 0.0);
        area += 0.5 * a;
    }
    EXPECT_NEAR(area, 1.0, 1e-12);
    // The smallest angle, measured in the metric at the centroid, averaged.
    const auto mean_angle = [&](const std::vector<mio::detail::Cdt2d::Tri>& rT,
                                const std::vector<Vec2>& rP) {
        double sum = 0.0;
        for (const auto& t : rT) {
            const Vec2 &a = rP[std::size_t(t[0])], &b = rP[std::size_t(t[1])],
                       &c = rP[std::size_t(t[2])];
            const double s = std::sqrt(stretch({(a[0] + b[0] + c[0]) / 3, 0.0})[2]);
            const auto len = [&](const Vec2& rU, const Vec2& rV) {
                return std::hypot(rU[0] - rV[0], s * (rU[1] - rV[1]));
            };
            const double x = len(b, c), y = len(c, a), z = len(a, b);
            const auto angle = [](double P, double Q, double R) {
                return std::acos(std::clamp((Q * Q + R * R - P * P) / (2 * Q * R), -1.0, 1.0));
            };
            sum += std::min({angle(x, y, z), angle(y, z, x), angle(z, x, y)});
        }
        return sum / double(rT.size()) * 180.0 / 3.14159265358979323846;
    };
    EXPECT_GT(mean_angle(aniso, aniso_pts), mean_angle(plane, plane_pts) + 5.0);
    // Cleared, the plane's test again.
    mio::detail::Cdt2d cdt(0, 0, 1, 1);
    cdt.SetMetric(stretch);
    cdt.SetMetric(nullptr);
    std::vector<Vec2> pts = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.4, 0.45}, {0.55, 0.6}};
    mio::detail::Cdt2d ref(0, 0, 1, 1);
    cdt.InsertBatch(pts);
    ref.InsertBatch(pts);
    EXPECT_EQ(cdt.AllTriangles(), ref.AllTriangles());
}

// --- Corner lopping (Pav and Walkington 2005) ------------------------------

namespace {

struct Lopped {
    std::vector<Vec2> mPoints;
    std::vector<mio::detail::Cdt2d::Tri> mTriangles;
    mio::detail::Cdt2dRefineReport mRefine;
};

// One loop, triangulated, classified and refined with corner lopping.
Lopped lopped(const std::vector<Vec2>& rLoop, double MinAngle, bool Split = true,
              std::int64_t MaxPoints = 200000, double Beta = 0.5, double Grid = 0.0) {
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    for (const Vec2& p : rLoop)
        for (int k = 0; k < 2; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    mio::detail::Cdt2d cdt(lo[0], lo[1], hi[0], hi[1]);
    const auto ids = cdt.InsertBatch(rLoop);
    for (std::size_t k = 0; k < ids.size(); ++k)
        EXPECT_EQ(cdt.InsertSegment(ids[k], ids[(k + 1) % ids.size()]), Cdt2dStatus::Ok);
    EXPECT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
    // Interior points on a grid (those inside), some of them in the balls.
    if (Grid > 0.0)
        for (double x = lo[0] + 0.5 * Grid; x < hi[0]; x += Grid)
            for (double y = lo[1] + 0.5 * Grid; y < hi[1]; y += Grid)
                cdt.InsertSteiner({x, y});
    mio::detail::Cdt2dRefineOptions o;
    o.mMinAngleDeg = MinAngle;
    o.mMaxPoints = MaxPoints;
    o.mSplitSegments = Split;
    o.mLopCorners = true;
    o.mLopBeta = Beta;
    Lopped out;
    out.mRefine = cdt.Refine(o);
    for (std::int64_t v = 0; v < cdt.NumVertices(); ++v)
        out.mPoints.push_back(cdt.Point(v));
    out.mTriangles = cdt.InsideTriangles();
    return out;
}

double lopped_min_angle(const Lopped& rL, std::size_t T) {
    double m = 180.0;
    const auto& tri = rL.mTriangles[T];
    for (int k = 0; k < 3; ++k) {
        const Vec2& o = rL.mPoints[std::size_t(tri[k])];
        const Vec2& p = rL.mPoints[std::size_t(tri[(k + 1) % 3])];
        const Vec2& q = rL.mPoints[std::size_t(tri[(k + 2) % 3])];
        const double ux = p[0] - o[0], uy = p[1] - o[1], vx = q[0] - o[0], vy = q[1] - o[1];
        m = std::min(m, std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / M_PI);
    }
    return m;
}

}  // namespace

TEST(Cdt2dLop, TheBoundIsCorollaryThree) {
    // arcsin(1 / (2 sqrt 2)) for straight segments, falling with the
    // turning, and none past pi / 3.
    EXPECT_NEAR(mio::detail::cdt2d_lop_bound_deg(0.0),
                std::asin(1.0 / (2.0 * std::sqrt(2.0))) * 180.0 / M_PI, 1e-12);
    EXPECT_NEAR(mio::detail::cdt2d_lop_bound_deg(1e-7), mio::detail::cdt2d_lop_bound_deg(0.0),
                1e-5);
    double last = mio::detail::cdt2d_lop_bound_deg(0.0);
    for (int k = 1; k <= 10; ++k) {
        const double b = mio::detail::cdt2d_lop_bound_deg(k * M_PI / 30.0);
        EXPECT_LT(b, last) << k;
        EXPECT_GT(b, 0.0) << k;
        last = b;
    }
    EXPECT_EQ(mio::detail::cdt2d_lop_bound_deg(M_PI / 3.0 + 1e-9), 0.0);
    EXPECT_EQ(mio::detail::cdt2d_lop_bound_deg(-1.0), 0.0);
}

TEST(Cdt2dLop, SharpCornersMeetTheBoundOutsideTheirBalls) {
    // A triangle with an apex of phi at the origin and a quadrilateral with
    // two sharp corners, refined at 20 degrees (under Corollary 3's 20.7):
    // the run ends without a budget or an unresolved triangle, the acute
    // corners are lopped with their angle and lfs0, and every triangle
    // under the bound lies within twice its corner's radius (Lemma 7).
    // With interior points, some in the balls, those shrink (rule 11).
    for (const auto& [phi, grid] : std::vector<std::pair<double, double>>{{80.0, 0.0},
                                                                          {45.0, 0.0},
                                                                          {15.0, 0.0},
                                                                          {5.0, 0.0},
                                                                          {1.0, 0.0},
                                                                          {45.0, 0.11},
                                                                          {15.0, 0.07},
                                                                          {5.0, 0.05},
                                                                          {1.0, 0.03}}) {
        const double r = phi * M_PI / 180.0;
        const std::vector<Vec2> loop = {
            {0, 0}, {2, 0}, {2, 2}, {1.8 * std::cos(r), 1.8 * std::sin(r)}};
        const auto t = lopped(loop, 20.0, true, 200000, 0.5, grid);
        EXPECT_FALSE(t.mRefine.mBudgetReached) << phi;
        EXPECT_EQ(t.mRefine.mNumUnresolved, 0) << phi;
        EXPECT_TRUE(t.mRefine.mSmallAngles.empty()) << phi;  // no shells
        bool found = false;
        for (const auto& r_c : t.mRefine.mLopped) {
            EXPECT_LT(r_c.mAngleDeg, 90.0);
            EXPECT_GT(r_c.mRadius, 0.0);
            EXPECT_LE(r_c.mRadius, (std::sqrt(2.0) - 1.0) * r_c.mLfs0 * (1 + 1e-12));
            EXPECT_TRUE(r_c.mWithinTwoRadii) << phi;
            if (r_c.mVertex == 3) {
                found = true;
                EXPECT_NEAR(r_c.mAngleDeg, phi, 1e-9);
                EXPECT_NEAR(r_c.mLfs0, 1.8, 1e-12);
            }
        }
        EXPECT_TRUE(found) << phi;
        if (grid > 0.0) {
            bool shrunk = false;
            for (const auto& r_c : t.mRefine.mLopped)
                shrunk = shrunk || r_c.mRadius < 0.5 * (std::sqrt(2.0) - 1.0) * r_c.mLfs0;
            EXPECT_TRUE(shrunk) << phi;
        }
        std::int64_t bad = 0;
        for (std::size_t k = 0; k < t.mTriangles.size(); ++k) {
            if (lopped_min_angle(t, k) >= 20.0 - 1e-9)
                continue;
            ++bad;
            bool near = false;
            for (const auto& r_c : t.mRefine.mLopped) {
                const Vec2& c = t.mPoints[std::size_t(r_c.mVertex)];
                bool all = true;
                for (const auto v : t.mTriangles[k]) {
                    const Vec2& p = t.mPoints[std::size_t(v)];
                    all = all &&
                          std::hypot(p[0] - c[0], p[1] - c[1]) <= 2.0 * r_c.mRadius * (1 + 1e-9);
                }
                near = near || all;
            }
            EXPECT_TRUE(near) << "phi " << phi << ", triangle " << k;
        }
        EXPECT_EQ(bad, t.mRefine.mNumForced) << phi;
    }
}

TEST(Cdt2dLop, ASplitIsAskedAtTheBall) {
    // Asking: the boundary is left, and each subsegment at an acute corner
    // is asked to split where the ball ends, at (sqrt 2 - 1) lfs0.
    const double r = 10.0 * M_PI / 180.0;
    const std::vector<Vec2> loop = {{0, 0}, {2, 0}, {2, 2}, {1.8 * std::cos(r), 1.8 * std::sin(r)}};
    const auto t = lopped(loop, 20.0, false);
    EXPECT_EQ(t.mRefine.mNumSplits, 0);
    int at_origin = 0;
    for (const auto& r_q : t.mRefine.mRequests)
        if (r_q.mApex == 3) {
            ++at_origin;
            EXPECT_NEAR(r_q.mDistance, (std::sqrt(2.0) - 1.0) * 1.8, 1e-12);
        }
    EXPECT_EQ(at_origin, 2);
}

TEST(Cdt2dLop, TheOrderDependsOnTheInputAloneAndBetaIsChecked) {
    const std::vector<Vec2> loop = {{0, 0}, {3, 0}, {3, 0.2}, {1.5, 0.05}, {0, 0.4}};
    const auto a = lopped(loop, 20.0);
    const auto b = lopped(loop, 20.0);
    EXPECT_EQ(a.mPoints, b.mPoints);
    EXPECT_EQ(a.mTriangles, b.mTriangles);
    EXPECT_FALSE(a.mRefine.mLopped.empty());
    const auto c = lopped(loop, 20.0, true, 200000, 0.75);
    EXPECT_FALSE(c.mRefine.mBudgetReached);
    EXPECT_EQ(c.mRefine.mNumUnresolved, 0);
    EXPECT_THROW(lopped(loop, 20.0, true, 200000, 0.4), std::invalid_argument);
    EXPECT_THROW(lopped(loop, 20.0, true, 200000, 1.0), std::invalid_argument);
}

TEST(Cdt2dLop, AnAcuteKinkIsAnInputPointAndASkippedOneIsNot) {
    // The caller names only the far corner (2, 2) as an input point: the
    // acute corner at the origin (a kink of the caller's curve) is lopped
    // all the same, unless it is skipped (a pole's).
    const double r = 10.0 * M_PI / 180.0;
    const std::vector<Vec2> loop = {{0, 0}, {2, 0}, {2, 2}, {1.8 * std::cos(r), 1.8 * std::sin(r)}};
    for (const bool skip : {false, true}) {
        mio::detail::Cdt2d cdt(0, 0, 2, 2);
        const auto ids = cdt.InsertBatch(loop);
        for (std::size_t k = 0; k < ids.size(); ++k)
            ASSERT_EQ(cdt.InsertSegment(ids[k], ids[(k + 1) % ids.size()]), Cdt2dStatus::Ok);
        ASSERT_EQ(cdt.Classify(), Cdt2dStatus::Ok);
        mio::detail::Cdt2dRefineOptions o;
        o.mMinAngleDeg = 20.0;
        o.mLopCorners = true;
        o.mCorners = {ids[2]};
        if (skip)
            o.mLopSkip = {ids[0]};
        const auto rep = cdt.Refine(o);
        bool origin = false;
        for (const auto& r_c : rep.mLopped)
            origin = origin || r_c.mVertex == ids[0];
        EXPECT_EQ(origin, !skip);
    }
}
