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
// The signed-distance kernel, tested against CLOSED-FORM distance fields with no
// grid anywhere in the picture.
//
// That separation is the point. Once a distance is sampled on a generated
// lattice, a wrong distance and a wrong lattice look identical from the outside;
// here the expected answer is an expression, so a failure says which.
//
// `ASharpSpikeIsWhereTheNaiveSignFails` is the fixture that matters. It
// demonstrates the standard mistake -- taking the sign from the nearest
// TRIANGLE's normal rather than the nearest FEATURE's pseudonormal -- rather
// than asserting the correct answer and hoping.
//
// Finding the right fixture took two attempts, and the first one is worth
// recording: a *reentrant* corner does NOT expose the bug. At the L-solid's
// concave edge both incident faces give the correct sign, so the naive method
// passes there. The failure needs two incident faces whose normals are nearly
// opposite, which is a sharp spike, not a notch.

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
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::sample_distance;
using meshioplusplus::SdfSign;
using meshioplusplus::SurfaceDistanceOptions;
using meshioplusplus::SurfaceQuality;
namespace d = meshioplusplus::detail;

namespace {

// --- fixtures ---------------------------------------------------------------

// A closed prism over a counter-clockwise xy polygon, wound outward: the bottom
// cap reversed, the top cap as given, and one outward quad per footprint edge.
Mesh prism(const std::vector<std::array<double, 2>>& rFootprint, double Lo, double Hi) {
    const std::size_t n = rFootprint.size();
    std::vector<std::vector<double>> pts;
    for (std::size_t i = 0; i < n; ++i)
        pts.push_back({rFootprint[i][0], rFootprint[i][1], Lo});
    for (std::size_t i = 0; i < n; ++i)
        pts.push_back({rFootprint[i][0], rFootprint[i][1], Hi});

    std::vector<std::vector<std::int64_t>> tris;
    const auto b = [n](std::size_t i) { return static_cast<std::int64_t>(i); };
    const auto t = [n](std::size_t i) { return static_cast<std::int64_t>(n + i); };
    for (std::size_t k = 1; k + 1 < n; ++k) {
        tris.push_back({b(0), b(k + 1), b(k)});  // bottom: reversed, so -z is out
        tris.push_back({t(0), t(k), t(k + 1)});  // top: as given, so +z is out
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        tris.push_back({b(i), b(j), t(j)});
        tris.push_back({b(i), t(j), t(i)});
    }
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

Mesh box_mesh(double Half) {
    return prism({{{-Half, -Half}}, {{Half, -Half}}, {{Half, Half}}, {{-Half, Half}}}, -Half,
                 Half);
}

// The L: the square [0,2]^2 with the quadrant x>1 and y>1 removed, extruded to
// z in [0,1]. The reentrant edge runs vertically through (1, 1).
Mesh l_solid() {
    return prism({{{0.0, 0.0}}, {{2.0, 0.0}}, {{2.0, 1.0}}, {{1.0, 1.0}}, {{1.0, 2.0}},
                  {{0.0, 2.0}}},
                 0.0, 1.0);
}

bool l_solid_contains(const d::Vec3& rP) {
    if (!(rP[0] > 0.0 && rP[0] < 2.0 && rP[1] > 0.0 && rP[1] < 2.0 && rP[2] > 0.0 && rP[2] < 1.0))
        return false;
    return !(rP[0] > 1.0 && rP[1] > 1.0);
}

// A UV sphere. Coarse on purpose: the tessellation error is then large enough to
// state a bound on, rather than hiding a real error underneath it.
Mesh uv_sphere(double Radius, std::size_t Rings, std::size_t Sectors) {
    std::vector<std::vector<double>> pts;
    const double pi = 3.14159265358979323846;
    pts.push_back({0.0, 0.0, Radius});
    for (std::size_t r = 1; r < Rings; ++r) {
        const double theta = pi * static_cast<double>(r) / static_cast<double>(Rings);
        for (std::size_t s = 0; s < Sectors; ++s) {
            const double phi = 2.0 * pi * static_cast<double>(s) / static_cast<double>(Sectors);
            pts.push_back({Radius * std::sin(theta) * std::cos(phi),
                           Radius * std::sin(theta) * std::sin(phi), Radius * std::cos(theta)});
        }
    }
    pts.push_back({0.0, 0.0, -Radius});
    const std::int64_t south = static_cast<std::int64_t>(pts.size()) - 1;
    const auto id = [Sectors](std::size_t r, std::size_t s) {
        return static_cast<std::int64_t>(1 + (r - 1) * Sectors + (s % Sectors));
    };

    std::vector<std::vector<std::int64_t>> tris;
    for (std::size_t s = 0; s < Sectors; ++s)
        tris.push_back({0, id(1, s), id(1, s + 1)});
    for (std::size_t r = 1; r + 1 < Rings; ++r)
        for (std::size_t s = 0; s < Sectors; ++s) {
            tris.push_back({id(r, s), id(r + 1, s), id(r + 1, s + 1)});
            tris.push_back({id(r, s), id(r + 1, s + 1), id(r, s + 1)});
        }
    for (std::size_t s = 0; s < Sectors; ++s)
        tris.push_back({south, id(Rings - 1, s + 1), id(Rings - 1, s)});
    return mt::make_mesh(std::move(pts), "triangle", std::move(tris));
}

// --- analytic fields --------------------------------------------------------

/// The exact signed distance to an axis-aligned box of half-extent `h`.
double box_sdf(const d::Vec3& rP, double Half) {
    d::Vec3 q{{std::fabs(rP[0]) - Half, std::fabs(rP[1]) - Half, std::fabs(rP[2]) - Half}};
    const d::Vec3 outside{{q[0] > 0.0 ? q[0] : 0.0, q[1] > 0.0 ? q[1] : 0.0,
                           q[2] > 0.0 ? q[2] : 0.0}};
    const double inside = std::min(std::max(q[0], std::max(q[1], q[2])), 0.0);
    return d::vec3_norm(outside) + inside;
}

// --- helpers ----------------------------------------------------------------

NDArray as_array(const std::vector<d::Vec3>& rPts) {
    NDArray a = NDArray::Uninit(meshioplusplus::DType::Float64, {rPts.size(), std::size_t{3}});
    double* dst = a.As<double>();
    for (std::size_t i = 0; i < rPts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            dst[i * 3 + k] = rPts[i][k];
    return a;
}

// A grid of query points over a box, deliberately offset off the half-integers
// so no query lands exactly on a facet.
std::vector<d::Vec3> sample_lattice(double Lo, double Hi, std::size_t N) {
    std::vector<d::Vec3> out;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t k = 0; k < N; ++k) {
                const auto c = [&](std::size_t t) {
                    return Lo + (Hi - Lo) * (static_cast<double>(t) + 0.317) /
                                    static_cast<double>(N);
                };
                out.push_back({{c(i), c(j), c(k)}});
            }
    return out;
}

/// The classic mistake, implemented on purpose: nearest triangle, then that
/// triangle's own normal. Used only to show that the reentrant-corner test
/// distinguishes it from the real thing.
double naive_face_normal_sdf(const d::TriangleSoup& rSoup, const d::Vec3& rP) {
    double best = std::numeric_limits<double>::infinity();
    std::size_t best_t = 0;
    d::PointTriangleHit best_hit;
    for (std::size_t t = 0; t < rSoup.NumTriangles(); ++t) {
        const d::PointTriangleHit hit = d::closest_point_on_triangle(
            rP, rSoup.mCorners[t * 3], rSoup.mCorners[t * 3 + 1], rSoup.mCorners[t * 3 + 2]);
        if (hit.mDistanceSq < best) {
            best = hit.mDistanceSq;
            best_t = t;
            best_hit = hit;
        }
    }
    const d::Vec3 n =
        d::vec3_cross(d::vec3_sub(rSoup.mCorners[best_t * 3 + 1], rSoup.mCorners[best_t * 3]),
                      d::vec3_sub(rSoup.mCorners[best_t * 3 + 2], rSoup.mCorners[best_t * 3]));
    const double side = d::vec3_dot(d::vec3_sub(rP, best_hit.mPoint), n);
    return side < 0.0 ? -std::sqrt(best) : std::sqrt(best);
}

}  // namespace

// --- the closest-point primitive --------------------------------------------

TEST(PointTriangle, EveryRegionIsReachedAndCorrect) {
    const d::Vec3 a{{0.0, 0.0, 0.0}};
    const d::Vec3 b{{1.0, 0.0, 0.0}};
    const d::Vec3 c{{0.0, 1.0, 0.0}};

    struct Case {
        d::Vec3 mQuery;
        d::TriangleFeature mFeature;
        d::Vec3 mClosest;
    };
    const Case cases[] = {
        {{{-1.0, -1.0, 0.0}}, d::TriangleFeature::VertexA, {{0.0, 0.0, 0.0}}},
        {{{2.0, -1.0, 0.0}}, d::TriangleFeature::VertexB, {{1.0, 0.0, 0.0}}},
        {{{-1.0, 2.0, 0.0}}, d::TriangleFeature::VertexC, {{0.0, 1.0, 0.0}}},
        {{{0.5, -1.0, 0.0}}, d::TriangleFeature::EdgeAB, {{0.5, 0.0, 0.0}}},
        {{{1.0, 1.0, 0.0}}, d::TriangleFeature::EdgeBC, {{0.5, 0.5, 0.0}}},
        {{{-1.0, 0.5, 0.0}}, d::TriangleFeature::EdgeCA, {{0.0, 0.5, 0.0}}},
        {{{0.25, 0.25, 1.0}}, d::TriangleFeature::Face, {{0.25, 0.25, 0.0}}},
    };
    for (const Case& k : cases) {
        const d::PointTriangleHit hit = d::closest_point_on_triangle(k.mQuery, a, b, c);
        EXPECT_EQ(static_cast<int>(hit.mFeature), static_cast<int>(k.mFeature));
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_NEAR(hit.mPoint[i], k.mClosest[i], 1e-15);
        EXPECT_NEAR(hit.mDistanceSq, d::vec3_norm_sq(d::vec3_sub(k.mQuery, k.mClosest)), 1e-15);
    }
}

TEST(PointTriangle, AnObtuseTriangleDefeatsTheClampedProjection) {
    // The point of Ericson's region classification: for an obtuse triangle the
    // "project onto the plane and clamp the barycentrics" shortcut lands on the
    // wrong edge. This query is nearest to edge AB, and a clamped projection
    // would place it on BC.
    const d::Vec3 a{{0.0, 0.0, 0.0}};
    const d::Vec3 b{{4.0, 0.0, 0.0}};
    const d::Vec3 c{{3.9, 0.2, 0.0}};
    const d::Vec3 p{{2.0, -1.0, 0.0}};
    const d::PointTriangleHit hit = d::closest_point_on_triangle(p, a, b, c);
    EXPECT_EQ(static_cast<int>(hit.mFeature), static_cast<int>(d::TriangleFeature::EdgeAB));
    EXPECT_NEAR(hit.mPoint[0], 2.0, 1e-15);
    EXPECT_NEAR(hit.mPoint[1], 0.0, 1e-15);
}

TEST(PointTriangle, ADegenerateTriangleYieldsAFiniteAnswer) {
    // Zero area: the face branch would divide by zero, and a NaN then wins every
    // comparison in a min-reduction rather than losing it.
    const d::Vec3 a{{0.0, 0.0, 0.0}};
    const d::Vec3 b{{1.0, 0.0, 0.0}};
    const d::Vec3 c{{2.0, 0.0, 0.0}};  // collinear
    const d::PointTriangleHit hit = d::closest_point_on_triangle({{1.0, 3.0, 0.0}}, a, b, c);
    EXPECT_TRUE(std::isfinite(hit.mDistanceSq));
    EXPECT_NEAR(hit.mDistanceSq, 9.0, 1e-12);
}

// --- analytic distance fields -----------------------------------------------

TEST(SurfaceDistance, BoxDistanceIsExactEverywhere) {
    // A box's surface is represented exactly by triangles, so there is no
    // tessellation error to hide behind: the answer must match the closed form.
    const double half = 1.0;
    const Mesh box = box_mesh(half);
    const std::vector<d::Vec3> queries = sample_lattice(-2.0, 2.0, 9);
    const NDArray got = sample_distance(box, as_array(queries));
    const double* v = got.As<double>();
    for (std::size_t i = 0; i < queries.size(); ++i)
        EXPECT_NEAR(v[i], box_sdf(queries[i], half), 1e-12)
            << "at (" << queries[i][0] << ", " << queries[i][1] << ", " << queries[i][2] << ")";
}

TEST(SurfaceDistance, QueryClosestPointsReturnsTheActualNearestPointOnABox) {
    // The whole reason query_closest_points exists rather than reusing
    // query_distances: a caller (remesh_volume's warp step) needs the actual
    // point, not just how far away it is. A box gives a closed-form nearest
    // point exactly (a straight axis-aligned projection onto whichever face is
    // nearest), so this is checkable to machine precision, not just plausible.
    const double half = 1.0;
    const Mesh box = box_mesh(half);
    const d::TriangleSoup soup = d::build_triangle_soup(box, "");
    const SurfaceDistanceOptions opts;
    const d::DistanceQuery query = d::build_distance_query(soup, opts);

    // Outside, off the +x face: nearest point is the straight projection.
    const d::Vec3 outside{{2.0, 0.3, 0.2}};
    // Inside, off-center: nearest point is still the straight projection onto
    // whichever face is closest (here +x, since 1 - 0.2 < the other five gaps).
    const d::Vec3 inside{{0.2, -0.1, 0.05}};
    const std::vector<d::Vec3> queries = {outside, inside};
    const std::vector<d::ClosestPointHit> hits = d::query_closest_points(query, queries);
    ASSERT_EQ(hits.size(), 2u);

    ASSERT_TRUE(hits[0].mFound);
    EXPECT_NEAR(hits[0].mPoint[0], 1.0, 1e-12);
    EXPECT_NEAR(hits[0].mPoint[1], 0.3, 1e-12);
    EXPECT_NEAR(hits[0].mPoint[2], 0.2, 1e-12);
    EXPECT_NEAR(hits[0].mDistance, 1.0, 1e-12);  // |2 - 1|
    EXPECT_GE(hits[0].mSourceCell, 0);

    ASSERT_TRUE(hits[1].mFound);
    EXPECT_NEAR(hits[1].mPoint[0], 1.0, 1e-12);
    EXPECT_NEAR(hits[1].mPoint[1], -0.1, 1e-12);
    EXPECT_NEAR(hits[1].mPoint[2], 0.05, 1e-12);
    EXPECT_NEAR(hits[1].mDistance, 0.8, 1e-12);  // |1 - 0.2|
}

TEST(SurfaceDistance, QueryClosestPointsAgreesWithQueryDistancesOnWhichTriangleIsNearest) {
    // The two functions share sd_nearest_triangle, so they cannot disagree
    // about which triangle is nearest -- checked here via the source-cell id
    // and the (unsigned) distance, on a fixture with genuine tie-break
    // structure (a UV sphere: many equidistant-ish facets).
    const double r = 1.0;
    const Mesh sphere = uv_sphere(r, 12, 24);
    const d::TriangleSoup soup = d::build_triangle_soup(sphere, "");
    const SurfaceDistanceOptions opts;
    const d::DistanceQuery query = d::build_distance_query(soup, opts);
    const std::vector<d::Vec3> queries = sample_lattice(-2.0, 2.0, 5);

    const std::vector<d::DistanceHit> dist_hits = d::query_distances(query, queries, opts);
    const std::vector<d::ClosestPointHit> point_hits = d::query_closest_points(query, queries);
    ASSERT_EQ(dist_hits.size(), point_hits.size());
    for (std::size_t i = 0; i < queries.size(); ++i) {
        ASSERT_TRUE(point_hits[i].mFound) << "index " << i;
        EXPECT_EQ(point_hits[i].mSourceCell, dist_hits[i].mSourceCell) << "index " << i;
        EXPECT_NEAR(point_hits[i].mDistance, std::fabs(dist_hits[i].mSignedDistance), 1e-12)
            << "index " << i;
    }
}

TEST(SurfaceDistance, SphereDistanceMatchesTheAnalyticFieldWithinTessellationError) {
    const double r = 1.0;
    const Mesh sphere = uv_sphere(r, 24, 48);
    const std::vector<d::Vec3> queries = sample_lattice(-2.0, 2.0, 7);
    const NDArray got = sample_distance(sphere, as_array(queries));
    const double* v = got.As<double>();
    // A UV sphere is inscribed, so its facets sit slightly inside the true
    // surface; the chord sagitta bounds the error.
    const double pi = 3.14159265358979323846;
    const double tol = r * (1.0 - std::cos(pi / 24.0)) + 1e-9;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        const double want = d::vec3_norm(queries[i]) - r;
        EXPECT_NEAR(v[i], want, tol) << "at index " << i;
        if (std::fabs(want) > tol)
            EXPECT_EQ(v[i] < 0.0, want < 0.0) << "sign at index " << i;
    }
}

TEST(SurfaceDistance, ASharpSpikeIsWhereTheNaiveSignFails) {
    // The whole reason the pseudonormal exists, and it takes the right fixture
    // to show it. A *reentrant* corner does not: at the L's concave edge both
    // incident faces happen to give the correct sign, so the naive method passes
    // there and an earlier version of this test proved nothing.
    //
    // The failure needs two incident faces whose normals are nearly OPPOSITE,
    // which is what a sharp spike has. This sliver prism's tip subtends about
    // 1.7 degrees, so its two side faces point almost exactly -y and +y. For a
    // query just off the tip, whichever of the two the nearest-triangle search
    // happens to return decides the sign -- and one of them is wrong.
    const Mesh spike = prism({{{0.0, 0.0}}, {{10.0, 0.0}}, {{10.0, 0.3}}}, 0.0, 1.0);
    const d::TriangleSoup soup = d::build_triangle_soup(spike, "");

    // A fan of points just outside the tip, sweeping the direction from which
    // they approach it.
    std::vector<d::Vec3> queries;
    const double pi = 3.14159265358979323846;
    for (int a = 0; a < 24; ++a) {
        const double theta = pi * (0.5 + static_cast<double>(a) / 24.0);  // up and to the left
        queries.push_back({{0.05 * std::cos(theta), 0.05 * std::sin(theta), 0.5}});
    }

    const NDArray got = sample_distance(spike, as_array(queries));
    const double* v = got.As<double>();

    std::size_t naive_wrong = 0;
    for (std::size_t i = 0; i < queries.size(); ++i) {
        // Every one of these is outside: they sit on a circle about the tip, in
        // the half-plane the solid does not occupy.
        EXPECT_GT(v[i], 0.0) << "pseudonormal sign wrong at (" << queries[i][0] << ", "
                             << queries[i][1] << ")";
        if (naive_face_normal_sdf(soup, queries[i]) < 0.0)
            ++naive_wrong;
    }
    EXPECT_GT(naive_wrong, 0u)
        << "the fixture does not distinguish the pseudonormal from the face normal";
}

TEST(SurfaceDistance, SignIsCorrectThroughoutTheLSolid) {
    const Mesh solid = l_solid();
    const std::vector<d::Vec3> queries = sample_lattice(-0.5, 2.5, 11);
    const NDArray got = sample_distance(solid, as_array(queries));
    const double* v = got.As<double>();
    for (std::size_t i = 0; i < queries.size(); ++i)
        EXPECT_EQ(v[i] < 0.0, l_solid_contains(queries[i]))
            << "at (" << queries[i][0] << ", " << queries[i][1] << ", " << queries[i][2] << ")";
}

TEST(SurfaceDistance, WindingNumberAgreesWithThePseudonormalOnAClosedSurface) {
    const Mesh solid = l_solid();
    const std::vector<d::Vec3> queries = sample_lattice(-0.5, 2.5, 7);
    SurfaceDistanceOptions wn;
    wn.mSign = SdfSign::WindingNumber;
    const NDArray a = sample_distance(solid, as_array(queries));
    const NDArray b = sample_distance(solid, as_array(queries), wn);
    const double* va = a.As<double>();
    const double* vb = b.As<double>();
    for (std::size_t i = 0; i < queries.size(); ++i) {
        EXPECT_NEAR(std::fabs(va[i]), std::fabs(vb[i]), 1e-12) << i;
        EXPECT_EQ(va[i] < 0.0, vb[i] < 0.0) << "sign disagreement at index " << i;
    }
}

TEST(SurfaceDistance, UnsignedIsAlwaysNonNegative) {
    SurfaceDistanceOptions un;
    un.mSign = SdfSign::Unsigned;
    const NDArray got = sample_distance(l_solid(), as_array(sample_lattice(-0.5, 2.5, 5)), un);
    const double* v = got.As<double>();
    for (std::size_t i = 0; i < got.Size(); ++i)
        EXPECT_GE(v[i], 0.0);
}

// --- the accelerator must not be observable ---------------------------------

TEST(SurfaceDistance, TheBucketSizeDoesNotChangeTheAnswer) {
    // The reason every comparison is totally ordered on (distance, triangle id).
    // A cube axis-aligned with the query lattice is chosen deliberately: its
    // coplanar triangle pairs are exactly equidistant from the points above them,
    // so ties are constant rather than rare. A fixture that avoids ties would
    // test nothing here.
    const Mesh box = box_mesh(1.0);
    const std::vector<d::Vec3> queries = sample_lattice(-2.0, 2.0, 7);
    const NDArray reference = sample_distance(box, as_array(queries));

    for (double cell : {0.1, 0.37, 1.0, 4.0}) {
        SurfaceDistanceOptions o;
        o.mGridCellSize = cell;
        const NDArray got = sample_distance(box, as_array(queries), o);
        ASSERT_EQ(got.Size(), reference.Size());
        const double* a = got.As<double>();
        const double* b = reference.As<double>();
        for (std::size_t i = 0; i < got.Size(); ++i)
            ASSERT_EQ(a[i], b[i]) << "cell size " << cell << " changed the distance at " << i;
    }
}

TEST(SurfaceDistance, TheNearestCellIsAlsoInvariantToTheBucketSize) {
    // The distance can be invariant while the *identity* of the winning triangle
    // is not, which is what would break sdf:closest_cell. Ties are what expose it.
    const Mesh box = box_mesh(1.0);
    Mesh query = box_mesh(1.0);
    SurfaceDistanceOptions base;
    base.mRecordClosestCell = true;
    const auto ids = [&](double Cell) {
        SurfaceDistanceOptions o = base;
        o.mGridCellSize = Cell;
        const meshioplusplus::SurfaceDistanceResult r =
            meshioplusplus::distance_to_surface(query, box, o);
        const NDArray& a = r.mMesh.PointData("sdf:closest_cell");
        std::vector<std::int64_t> out(a.Size());
        for (std::size_t i = 0; i < a.Size(); ++i)
            out[i] = a.As<std::int64_t>()[i];
        return out;
    };
    const std::vector<std::int64_t> reference = ids(0.0);
    for (double cell : {0.13, 0.5, 3.0})
        EXPECT_EQ(ids(cell), reference) << "cell size " << cell << " changed the nearest cell";
}

// --- what is wrong with the surface -----------------------------------------

TEST(SurfaceDistance, WatertightCheckCountsTheDefects) {
    const SurfaceQuality closed = meshioplusplus::surface_watertight_check(box_mesh(1.0));
    EXPECT_TRUE(closed.mWatertight);
    EXPECT_EQ(closed.mBoundaryEdges, 0);
    EXPECT_EQ(closed.mNonManifoldEdges, 0);
    EXPECT_EQ(closed.mInconsistentPairs, 0);
    EXPECT_EQ(closed.mDegenerateTriangles, 0);

    // A single triangle is a sheet: three boundary edges, and it says so in
    // numbers rather than as a bare "not watertight".
    const Mesh sheet = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    const SurfaceQuality open = meshioplusplus::surface_watertight_check(sheet);
    EXPECT_FALSE(open.mWatertight);
    EXPECT_EQ(open.mBoundaryEdges, 3);
}

TEST(SurfaceDistance, WatertightCheckErrorModeThrowsByName) {
    const Mesh sheet = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    SurfaceDistanceOptions o;
    o.mWatertightCheck = meshioplusplus::SdfWatertightCheck::Error;
    try {
        sample_distance(sheet, as_array({{{5.0, 5.0, 5.0}}}), o);
        FAIL() << "the watertight check did not refuse";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("boundary edge"), std::string::npos) << e.what();
    }
}

TEST(SurfaceDistance, AVolumeMeshIsRefusedByName) {
    try {
        meshioplusplus::surface_watertight_check(mt::tet_mesh());
        FAIL() << "a volume mesh was accepted as a surface";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("extract_surface"), std::string::npos) << e.what();
    }
}

// --- the band ---------------------------------------------------------------

TEST(SurfaceDistance, ABandedRunAgreesExactlyInsideTheBand) {
    // The band is an optimization, so where it did not clamp it must produce the
    // identical bytes -- otherwise it is a second, quietly different answer.
    const Mesh box = box_mesh(1.0);
    // A lattice rather than another box, so the query points span a range of
    // distances and the band actually has both sides to compare.
    const Mesh query = meshioplusplus::grid({{6, 6, 6}}, {{-2.0, -2.0, -2.0}},
                                            {{4.0 / 6.0, 4.0 / 6.0, 4.0 / 6.0}});
    const meshioplusplus::SurfaceDistanceResult full =
        meshioplusplus::distance_to_surface(query, box);

    SurfaceDistanceOptions o;
    o.mBand = 0.75;
    const meshioplusplus::SurfaceDistanceResult banded =
        meshioplusplus::distance_to_surface(query, box, o);

    const NDArray& fd = full.mMesh.PointData("sdf:distance");
    const NDArray& bd = banded.mMesh.PointData("sdf:distance");
    ASSERT_TRUE(banded.mMesh.HasPointData("sdf:band"));
    const NDArray& flag = banded.mMesh.PointData("sdf:band");
    ASSERT_EQ(fd.Size(), bd.Size());
    std::size_t compared = 0;
    for (std::size_t i = 0; i < fd.Size(); ++i) {
        if (flag.As<std::int64_t>()[i] == 0)
            continue;
        ASSERT_EQ(fd.As<double>()[i], bd.As<double>()[i]) << "inside the band at " << i;
        ++compared;
    }
    EXPECT_GT(compared, 0u) << "the band clamped everything, so nothing was compared";
    EXPECT_GT(banded.mNumBanded, 0) << "the band clamped nothing, so it was not exercised";
}
