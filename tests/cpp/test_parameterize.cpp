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
// parameterize: a UV atlas of a triangle surface (roadmap 7.3.3).

#include <cmath>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/parameterize.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// A torus of R = 2, r = 0.6, N x M quads split in two.
Mesh prm_torus(int N = 24, int M = 12) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            const double u = 2.0 * M_PI * i / N, v = 2.0 * M_PI * j / M;
            pts.push_back({(2.0 + 0.6 * std::cos(v)) * std::cos(u),
                           (2.0 + 0.6 * std::cos(v)) * std::sin(u), 0.6 * std::sin(v)});
        }
    std::vector<std::vector<std::int64_t>> tris;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            const std::int64_t a = i * M + j, b = ((i + 1) % N) * M + j;
            const std::int64_t c = ((i + 1) % N) * M + (j + 1) % M, d = i * M + (j + 1) % M;
            tris.push_back({a, b, c});
            tris.push_back({a, c, d});
        }
    return mt::make_mesh(pts, "triangle", tris);
}

// The checks every atlas passes: uv in the unit square, every triangle
// counter-clockwise in uv, parents in range, seams only where charts meet.
void expect_valid_atlas(const Mesh& rIn, const ParameterizeResult& rR) {
    const Mesh& m = rR.mMesh;
    ASSERT_TRUE(m.HasPointData(kUvName));
    const NDArray& uv = m.PointData(kUvName);
    const NDArray& parents = m.PointData(kParameterizeParentPointName);
    ASSERT_EQ(uv.Size(), 2 * m.NumPoints());
    EXPECT_EQ(m.NumPoints(), rIn.NumPoints() + static_cast<std::size_t>(rR.mNumSeamPoints));
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        const double u = detail::read_double(uv, 2 * i), v = detail::read_double(uv, 2 * i + 1);
        EXPECT_GE(u, -1e-12);
        EXPECT_LE(u, 1.0 + 1e-12);
        EXPECT_GE(v, -1e-12);
        EXPECT_LE(v, 1.0 + 1e-12);
        const std::int64_t p = detail::read_int(parents, i);
        ASSERT_GE(p, 0);
        ASSERT_LT(static_cast<std::size_t>(p), rIn.NumPoints());
        if (i < rIn.NumPoints())
            EXPECT_EQ(p, static_cast<std::int64_t>(i));  // the first copy keeps its index
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_EQ(detail::read_double(m.Points(), 3 * i + k),
                      detail::read_double(rIn.Points(), 3 * static_cast<std::size_t>(p) + k));
    }
    const NDArray& conn = m.Cells(0).Conn();
    std::size_t flipped = 0;
    for (std::size_t c = 0; c < m.Cells(0).NumCells(); ++c) {
        double s[3], t[3];
        for (std::size_t k = 0; k < 3; ++k) {
            const std::size_t p = static_cast<std::size_t>(detail::read_int(conn, 3 * c + k));
            s[k] = detail::read_double(uv, 2 * p);
            t[k] = detail::read_double(uv, 2 * p + 1);
        }
        if (!((s[1] - s[0]) * (t[2] - t[0]) - (s[2] - s[0]) * (t[1] - t[0]) > 0.0))
            ++flipped;
    }
    EXPECT_EQ(static_cast<std::int64_t>(flipped), rR.mNumFlipped);
    EXPECT_GE(rR.mStretch, 1.0 - 1e-9);
    EXPECT_GT(rR.mPackingEfficiency, 0.0);
    EXPECT_LE(rR.mPackingEfficiency, 1.0 + 1e-12);
}

}  // namespace

TEST(Parameterize, ADiscIsOneChartWithNoFlippedTriangle) {
    const Mesh disc = disk(1.0, 6, {0.0, 0.0, 0.0}, "triangle");
    const ParameterizeResult r = parameterize(disc, {});
    EXPECT_EQ(r.mNumCharts, 1);
    EXPECT_EQ(r.mNumFlipped, 0);
    EXPECT_EQ(r.mNumSeamPoints, 0);
    expect_valid_atlas(disc, r);
    // A flat disc flattens isometrically (up to the solver).
    EXPECT_LT(r.mStretch, 1.0 + 1e-6);
}

TEST(Parameterize, ASphereAndATorusGetAtlasesWithoutFlips) {
    const Mesh sphere = icosa_sphere(1.0, 3);
    const ParameterizeResult rs = parameterize(sphere, {});
    EXPECT_GE(rs.mNumCharts, 2);
    EXPECT_EQ(rs.mNumFlipped, 0);
    EXPECT_GT(rs.mNumSeamPoints, 0);
    expect_valid_atlas(sphere, rs);

    const Mesh ring = prm_torus();
    const ParameterizeResult rt = parameterize(ring, {});
    EXPECT_GE(rt.mNumCharts, 2);
    EXPECT_EQ(rt.mNumFlipped, 0);
    expect_valid_atlas(ring, rt);
}

TEST(Parameterize, ChartsAreRecordedAndDataRidesAlong) {
    Mesh sphere = icosa_sphere(1.0, 2);
    std::vector<std::vector<double>> h;
    for (std::size_t i = 0; i < sphere.NumPoints(); ++i)
        h.push_back({detail::read_double(sphere.Points(), 3 * i + 2)});
    sphere.AddPointData("h", mt::points_from(h));
    const ParameterizeResult r = parameterize(sphere, {});
    const Mesh& m = r.mMesh;
    const NDArray& hh = m.PointData("h");
    const NDArray& parents = m.PointData(kParameterizeParentPointName);
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        EXPECT_EQ(detail::read_double(hh, i),
                  detail::read_double(sphere.PointData("h"),
                                      static_cast<std::size_t>(detail::read_int(parents, i))));
    std::set<std::int64_t> charts;
    const NDArray& ids = m.CellData(kParameterizeChartName, 0);
    for (std::size_t c = 0; c < ids.Size(); ++c)
        charts.insert(detail::read_int(ids, c));
    EXPECT_EQ(static_cast<std::int64_t>(charts.size()), r.mNumCharts);
    // Deterministic.
    const ParameterizeResult again = parameterize(sphere, {});
    EXPECT_EQ(detail::read_double(again.mMesh.PointData(kUvName), 7),
              detail::read_double(m.PointData(kUvName), 7));
}

TEST(Parameterize, PaddingSeparatesCharts) {
    const Mesh sphere = icosa_sphere(1.0, 2);
    ParameterizeOptions tight;
    tight.mPadding = 0.0;
    ParameterizeOptions loose;
    loose.mPadding = 0.2;
    EXPECT_GT(parameterize(sphere, tight).mPackingEfficiency,
              parameterize(sphere, loose).mPackingEfficiency);
}

TEST(Parameterize, RefusesWhatIsNotATriangleSurface) {
    EXPECT_THROW(parameterize(box({1, 1, 1}), {}), std::invalid_argument);
    EXPECT_THROW(parameterize(disk(1.0, 4), {}), std::invalid_argument);  // quads
    ParameterizeOptions bad;
    bad.mPadding = -1.0;
    EXPECT_THROW(parameterize(disk(1.0, 4, {0.0, 0.0, 0.0}, "triangle"), bad),
                 std::invalid_argument);
    bad.mPadding = 0.02;
    bad.mMaxDepth = -1;
    EXPECT_THROW(parameterize(disk(1.0, 4, {0.0, 0.0, 0.0}, "triangle"), bad),
                 std::invalid_argument);
}
