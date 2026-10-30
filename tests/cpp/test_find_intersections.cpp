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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/tri_tri.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/find_intersections.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"

using namespace meshioplusplus;
using mt::make_mesh;

namespace {

Mesh fi_sphere(double Cx, int Subdivisions = 3) {
    return icosa_sphere(1.0, Subdivisions, {{Cx, 0.0, 0.0}});
}

Mesh fi_merged(const Mesh& rA, const Mesh& rB) {
    return merge({&rA, &rB}).mMesh;
}

// Every triangle with points of its own: a seam everywhere.
Mesh fi_unwelded(const Mesh& rMesh) {
    const NDArray& pts = rMesh.Points();
    const auto cells = rMesh.Cells(0);
    const NDArray& conn = cells.Conn();
    const std::size_t nc = cells.NumCells();
    std::vector<std::vector<double>> p;
    std::vector<std::vector<std::int64_t>> c;
    for (std::size_t i = 0; i < nc; ++i) {
        std::vector<std::int64_t> row;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto v = static_cast<std::size_t>(conn.As<std::int64_t>()[3 * i + k]);
            p.push_back({pts.As<double>()[3 * v], pts.As<double>()[3 * v + 1],
                         pts.As<double>()[3 * v + 2]});
            row.push_back(static_cast<std::int64_t>(p.size() - 1));
        }
        c.push_back(row);
    }
    return make_mesh(p, "triangle", c);
}

std::int64_t fi_count_sum(const Mesh& rMesh) {
    std::int64_t s = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(kIntersectionsCountName, b);
        for (std::size_t i = 0; i < a.Size(); ++i)
            s += a.As<std::int64_t>()[i];
    }
    return s;
}

}  // namespace

TEST(TriTriKernel, SegmentOfTwoCrossingTriangles) {
    const double a[3] = {0, 0, 0}, b[3] = {2, 0, 0}, c[3] = {0, 2, 0};
    const double d[3] = {0.5, -1, -1}, e[3] = {0.5, -1, 1}, f[3] = {0.5, 3, 0};
    const detail::TriTriSegment s = detail::tri_tri_segment(a, b, c, d, e, f);
    ASSERT_EQ(s.mKind, detail::TriTriContact::Segment);
    const double lo = std::min(s.mP[1], s.mQ[1]), hi = std::max(s.mP[1], s.mQ[1]);
    EXPECT_NEAR(s.mP[0], 0.5, 1e-14);
    EXPECT_NEAR(s.mQ[0], 0.5, 1e-14);
    EXPECT_NEAR(lo, 0.0, 1e-14);
    EXPECT_NEAR(hi, 1.5, 1e-14);
    EXPECT_NEAR(s.mP[2], 0.0, 1e-14);
    EXPECT_NEAR(s.mQ[2], 0.0, 1e-14);
    // The decision is tri_tri_intersect's, on pseudo-random pairs.
    std::uint64_t state = 12345;
    const auto next = [&]() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) / 9007199254740992.0;
    };
    int met = 0;
    for (int k = 0; k < 2000; ++k) {
        double t[6][3];
        for (auto& r_p : t)
            for (double& r_c : r_p)
                r_c = next();
        const bool decided = detail::tri_tri_intersect(t[0], t[1], t[2], t[3], t[4], t[5]);
        const detail::TriTriSegment seg =
            detail::tri_tri_segment(t[0], t[1], t[2], t[3], t[4], t[5]);
        EXPECT_EQ(decided, seg.mKind != detail::TriTriContact::None);
        met += decided;
        if (seg.mKind == detail::TriTriContact::Segment)
            for (int i = 0; i < 3; ++i) {
                EXPECT_GE(seg.mP[i], -1e-12);
                EXPECT_LE(seg.mP[i], 1.0 + 1e-12);
            }
    }
    EXPECT_GT(met, 0);
    const double far[3][3] = {{5, 5, 5}, {6, 5, 5}, {5, 6, 5}};
    EXPECT_EQ(detail::tri_tri_segment(a, b, c, far[0], far[1], far[2]).mKind,
              detail::TriTriContact::None);
    const double cp[3][3] = {{0.5, 0.5, 0}, {3, 0.5, 0}, {0.5, 3, 0}};
    EXPECT_EQ(detail::tri_tri_segment(a, b, c, cp[0], cp[1], cp[2]).mKind,
              detail::TriTriContact::Coplanar);
}

TEST(TriTriKernel, SharedCornersAreNotCrossings) {
    const double o[3] = {0, 0, 0}, x[3] = {1, 0, 0}, y[3] = {0, 1, 0};
    // A shared edge, bent: no crossing.
    const double up[3] = {0.5, 0.5, 1};
    EXPECT_TRUE(detail::tri_tri_intersect(o, x, y, o, x, up));
    EXPECT_FALSE(detail::tri_tri_intersect_unshared(o, x, y, o, x, up));
    // A shared edge, folded flat onto the other: a crossing.
    const double fold[3] = {0.5, 0.25, 0};
    EXPECT_TRUE(detail::tri_tri_intersect_unshared(o, x, y, o, x, fold));
    // Flat but on the other side of the edge: no crossing.
    const double flat[3] = {0.5, -1, 0};
    EXPECT_FALSE(detail::tri_tri_intersect_unshared(o, x, y, x, o, flat));
    // A shared corner, the other triangle pierced by an edge.
    const double p[3] = {0.3, 0.3, -1}, q[3] = {0.3, 0.3, 1};
    EXPECT_TRUE(detail::tri_tri_intersect_unshared(o, x, y, o, p, q));
    // A shared corner, the other triangle away.
    const double r[3] = {-1, 0, 1}, s[3] = {-1, 1, 1};
    EXPECT_FALSE(detail::tri_tri_intersect_unshared(o, x, y, o, r, s));
    // A shared corner, an edge running into the other triangle in its plane.
    const double t[3] = {0.2, 0.2, 0}, u[3] = {0, 0, 1};
    EXPECT_TRUE(detail::tri_tri_intersect_unshared(o, x, y, o, t, u));
    // The same triangle twice.
    EXPECT_TRUE(detail::tri_tri_intersect_unshared(o, x, y, y, o, x));
}

TEST(FindIntersections, CleanSurfacesHaveNoPairs) {
    const FindIntersectionsResult r = find_intersections(fi_sphere(0.0));
    EXPECT_TRUE(r.mClean);
    EXPECT_EQ(r.mNumPairs, 0);
    EXPECT_EQ(r.mNumSegments, 0);
    EXPECT_EQ(fi_count_sum(r.mMesh), 0);
    // An unwelded seam at every edge is still no crossing.
    EXPECT_TRUE(find_intersections(fi_unwelded(fi_sphere(0.0, 2))).mClean);
    // Quads, each fanned in two.
    EXPECT_TRUE(find_intersections(extract_surface(box({{3, 2, 2}}))).mClean);
    EXPECT_TRUE(find_intersections(uv_sphere(1.0, 12, 6, {{0, 0, 0}}, "quad")).mClean);
}

TEST(FindIntersections, TwoOverlappingSpheres) {
    const Mesh a = fi_sphere(0.0), b = fi_sphere(1.0);
    const FindIntersectionsResult self = find_intersections(fi_merged(a, b));
    EXPECT_FALSE(self.mClean);
    EXPECT_GT(self.mNumPairs, 0);
    EXPECT_EQ(self.mNumCoplanar, 0);
    EXPECT_EQ(static_cast<std::int64_t>(self.mPairs.size()), self.mNumPairs);
    EXPECT_EQ(fi_count_sum(self.mMesh), 2 * self.mNumPairs);
    // Every segment lies near the circle x = 0.5, r = sqrt(0.75), and the
    // segments add up to its length.
    ASSERT_GT(self.mNumSegments, 0);
    const NDArray& p = self.mSegments.Points();
    double length = 0.0;
    for (std::int64_t s = 0; s < self.mNumSegments; ++s) {
        const double* e0 = p.As<double>() + 6 * s;
        const double* e1 = e0 + 3;
        for (const double* e : {e0, e1}) {
            EXPECT_NEAR(e[0], 0.5, 0.05);
            EXPECT_NEAR(std::hypot(e[1], e[2]), std::sqrt(0.75), 0.05);
        }
        length += std::hypot(e1[0] - e0[0], e1[1] - e0[1], e1[2] - e0[2]);
    }
    EXPECT_NEAR(length, 2.0 * M_PI * std::sqrt(0.75), 0.1);
    EXPECT_EQ(self.mSegments.CellData(kIntersectionsCellAName, 0).Size(),
              static_cast<std::size_t>(self.mNumSegments));

    // The two-mesh form finds the same pairs, one side each.
    const FindIntersectionsResult two = find_intersections(a, b);
    EXPECT_EQ(two.mNumPairs, self.mNumPairs);
    EXPECT_EQ(two.mNumSegments, self.mNumSegments);
    EXPECT_GT(two.mNumCellsA, 0);
    EXPECT_GT(two.mNumCellsB, 0);
    EXPECT_EQ(fi_count_sum(two.mMesh), two.mNumPairs);
    EXPECT_EQ(fi_count_sum(two.mMeshB), two.mNumPairs);
    const std::int64_t na = static_cast<std::int64_t>(a.Cells(0).NumCells());
    for (const auto& r_pair : self.mPairs) {
        EXPECT_LT(r_pair[0], na);
        EXPECT_GE(r_pair[1], na);
    }

    // Deterministic, and the cap keeps the counts whole.
    const FindIntersectionsResult again = find_intersections(fi_merged(a, b));
    EXPECT_EQ(again.mPairs, self.mPairs);
    FindIntersectionsOptions capped;
    capped.mMaxPairs = 3;
    capped.mSegments = false;
    const FindIntersectionsResult c = find_intersections(a, b, capped);
    EXPECT_TRUE(c.mTruncated);
    EXPECT_EQ(c.mPairs.size(), 3u);
    EXPECT_EQ(c.mNumPairs, two.mNumPairs);
    EXPECT_EQ(c.mNumSegments, 0);
    EXPECT_FALSE(c.mWarnings.empty());
}

TEST(FindIntersections, ApartAndTouching) {
    EXPECT_TRUE(find_intersections(fi_sphere(0.0), fi_sphere(3.0)).mClean);
    // A fold-over inside one surface, coplanar.
    const Mesh fold = make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0.25, 0}}, "triangle",
                                {{0, 1, 2}, {0, 1, 3}});
    const FindIntersectionsResult r = find_intersections(fold);
    EXPECT_EQ(r.mNumPairs, 1);
    EXPECT_EQ(r.mNumCoplanar, 1);
    EXPECT_EQ(r.mNumSegments, 0);
    // Two parts touching at a corner count in the two-mesh form.
    const Mesh t1 = make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    const Mesh t2 = make_mesh({{0, 0, 0}, {-1, 0, 1}, {-1, 1, 1}}, "triangle", {{0, 1, 2}});
    EXPECT_EQ(find_intersections(t1, t2).mNumPairs, 1);
    EXPECT_TRUE(find_intersections(fi_merged(t1, t2)).mClean);
}

TEST(FindIntersections, RejectsWhatItCannotTest) {
    EXPECT_THROW(find_intersections(box()), std::invalid_argument);
    FindIntersectionsOptions bad;
    bad.mMaxPairs = -1;
    EXPECT_THROW(find_intersections(fi_sphere(0.0, 1), bad), std::invalid_argument);
    bad = {};
    bad.mRegion = "nope";
    EXPECT_THROW(find_intersections(fi_sphere(0.0, 1), bad), std::invalid_argument);
    // Lines only: nothing to test, with a warning.
    const FindIntersectionsResult r = find_intersections(mt::line_mesh());
    EXPECT_TRUE(r.mClean);
    EXPECT_FALSE(r.mWarnings.empty());
}
