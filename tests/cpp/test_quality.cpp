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
// Tests for the mesh quality operation: regular elements give ideal metric
// values, deliberately-inverted/degenerate elements are flagged.

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace {

using meshioplusplus::compute_quality;
using meshioplusplus::Mesh;
using meshioplusplus::QualityReport;

double q_val(const QualityReport& rRep, const std::string& rName, std::size_t block,
             std::size_t cell) {
    for (const auto& e : rRep.mCellArrays)
        if (e.first == rName)
            return e.second[block].template As<double>()[cell];
    return std::nan("");
}

const double kSqrt2 = std::sqrt(2.0);

}  // namespace

TEST(Quality, UnitCubeHexIsIdeal) {
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:volume", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:skewness", 0, 0), 0.0, 1e-12);
    EXPECT_EQ(r.mNumInverted, 0);
    EXPECT_EQ(r.mNumDegenerate, 0);
}

TEST(Quality, RegularTetraIsIdeal) {
    const double s = kSqrt2;
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {s, 0, 0},
                            {s / 2, s * std::sqrt(3.0) / 2, 0},
                            {s / 2, s * std::sqrt(3.0) / 6, s * std::sqrt(2.0 / 3.0)}},
                           "tetra", {{0, 1, 2, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_GT(q_val(r, "quality:volume", 0, 0), 0.0);
    EXPECT_NEAR(q_val(r, "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    const double ideal_dih = std::acos(1.0 / 3.0) * 180.0 / M_PI;  // ~70.5288
    EXPECT_NEAR(q_val(r, "quality:min_dihedral", 0, 0), ideal_dih, 1e-9);
    EXPECT_NEAR(q_val(r, "quality:max_dihedral", 0, 0), ideal_dih, 1e-9);
    EXPECT_EQ(r.mNumInverted, 0);
}

TEST(Quality, EquilateralTriangle2D) {
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {0.5, std::sqrt(3.0) / 2}}, "triangle", {{0, 1, 2}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:min_angle", 0, 0), 60.0, 1e-9);
    EXPECT_NEAR(q_val(r, "quality:max_angle", 0, 0), 60.0, 1e-9);
}

TEST(Quality, UnitSquare2D) {
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "quad", {{0, 1, 2, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_NEAR(q_val(r, "quality:volume", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:skewness", 0, 0), 0.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:aspect_ratio", 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(q_val(r, "quality:warpage", 0, 0), 0.0, 1e-12);
}

TEST(Quality, ReferenceWedgeAndPyramidScaledJacobian) {
    const double h = std::sqrt(3.0) / 2;
    Mesh w = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0.5, h, 0}, {0, 0, 1}, {1, 0, 1}, {0.5, h, 1}},
                           "wedge", {{0, 1, 2, 3, 4, 5}});
    EXPECT_NEAR(q_val(compute_quality(w), "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);

    const double ph = 1.0 / kSqrt2;
    Mesh p = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, ph}}, "pyramid",
                           {{0, 1, 2, 3, 4}});
    EXPECT_NEAR(q_val(compute_quality(p), "quality:scaled_jacobian", 0, 0), 1.0, 1e-12);
}

TEST(Quality, InvertedTetraIsFlagged) {
    // A regular tetra with two nodes swapped -> negative signed volume.
    const double s = kSqrt2;
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {s, 0, 0},
                            {s / 2, s * std::sqrt(3.0) / 2, 0},
                            {s / 2, s * std::sqrt(3.0) / 6, s * std::sqrt(2.0 / 3.0)}},
                           "tetra", {{0, 2, 1, 3}});
    QualityReport r = compute_quality(m);
    EXPECT_LT(q_val(r, "quality:volume", 0, 0), 0.0);
    EXPECT_EQ(q_val(r, "quality:inverted", 0, 0), 1.0);
    EXPECT_EQ(r.mNumInverted, 1);
}

TEST(Quality, DegenerateTriangleIsFlagged) {
    // Three collinear points -> zero area.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, "triangle", {{0, 1, 2}});
    QualityReport r = compute_quality(m);
    EXPECT_EQ(q_val(r, "quality:degenerate", 0, 0), 1.0);
    EXPECT_EQ(r.mNumDegenerate, 1);
    EXPECT_TRUE(std::isnan(q_val(r, "quality:aspect_ratio", 0, 0)));
}

TEST(Quality, AttachQualityAddsCellData) {
    Mesh m = mt::tet_mesh();
    Mesh out = meshioplusplus::attach_quality(m);
    EXPECT_TRUE(out.HasCellData("quality:volume"));
    EXPECT_TRUE(out.HasCellData("quality:inverted"));
    EXPECT_EQ(out.CellDataNumBlocks("quality:volume"), m.NumCellBlocks());
}

TEST(Quality, QuadraticTetraUsesCorners) {
    // tetra10: the metrics should match the underlying linear tetra corners.
    Mesh lin = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    Mesh quad = mt::tet10_mesh();
    // Just assert tetra10 is scored (finite volume), exercising the corner path.
    EXPECT_TRUE(std::isfinite(q_val(compute_quality(quad), "quality:volume", 0, 0)));
    EXPECT_GT(q_val(compute_quality(lin), "quality:volume", 0, 0), 0.0);
}

TEST(Quality, DeterministicAcrossRuns) {
    using meshioplusplus::DType;
    using meshioplusplus::NDArray;

    // A jittered mixed hex+tetra mesh spanning several histogram chunks; the
    // whole report must come back byte-identical on every run.
    const std::size_t N = 20;  // 8000 hexes
    const std::size_t np = (N + 1) * (N + 1) * (N + 1);
    NDArray pts = NDArray::Uninit(DType::Float64, {np, 3});
    double* p = pts.As<double>();
    std::size_t idx = 0;
    for (std::size_t k = 0; k <= N; ++k)
        for (std::size_t j = 0; j <= N; ++j)
            for (std::size_t i = 0; i <= N; ++i) {
                p[idx * 3 + 0] = static_cast<double>(i) + 0.2 * std::sin(1.7 * idx);
                p[idx * 3 + 1] = static_cast<double>(j) + 0.2 * std::sin(2.3 * idx + 1.0);
                p[idx * 3 + 2] = static_cast<double>(k) + 0.2 * std::sin(3.1 * idx + 2.0);
                ++idx;
            }
    auto nid = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * (N + 1) + j) * (N + 1) + i);
    };
    NDArray hexes = NDArray::Uninit(DType::Int64, {N * N * N, 8});
    std::int64_t* h = hexes.As<std::int64_t>();
    std::size_t c = 0;
    for (std::size_t k = 0; k < N; ++k)
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < N; ++i) {
                h[c * 8 + 0] = nid(i, j, k);
                h[c * 8 + 1] = nid(i + 1, j, k);
                h[c * 8 + 2] = nid(i + 1, j + 1, k);
                h[c * 8 + 3] = nid(i, j + 1, k);
                h[c * 8 + 4] = nid(i, j, k + 1);
                h[c * 8 + 5] = nid(i + 1, j, k + 1);
                h[c * 8 + 6] = nid(i + 1, j + 1, k + 1);
                h[c * 8 + 7] = nid(i, j + 1, k + 1);
                ++c;
            }
    const std::size_t ntets = 600;
    NDArray tets = NDArray::Uninit(DType::Int64, {ntets, 4});
    std::int64_t* t = tets.As<std::int64_t>();
    for (std::size_t q = 0; q < ntets; ++q) {
        t[q * 4 + 0] = h[q * 8 + 0];
        t[q * 4 + 1] = h[q * 8 + 1];
        t[q * 4 + 2] = h[q * 8 + 3];
        t[q * 4 + 3] = h[q * 8 + 4];
    }
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("hexahedron", std::move(hexes));
    m.AddCellBlock("tetra", std::move(tets));

    const QualityReport first = compute_quality(m);
    for (int run = 0; run < 3; ++run) {
        const QualityReport again = compute_quality(m);
        ASSERT_EQ(again.mNumInverted, first.mNumInverted);
        ASSERT_EQ(again.mNumDegenerate, first.mNumDegenerate);
        ASSERT_EQ(again.mMetrics.size(), first.mMetrics.size());
        for (std::size_t mi = 0; mi < first.mMetrics.size(); ++mi) {
            const auto& a = first.mMetrics[mi].second;
            const auto& b = again.mMetrics[mi].second;
            ASSERT_EQ(b.mCount, a.mCount);
            // Bitwise, not merely close (and NaN-safe for count-0 metrics).
            ASSERT_EQ(std::memcmp(&b.mMin, &a.mMin, sizeof(double)), 0) << "run " << run;
            ASSERT_EQ(std::memcmp(&b.mMax, &a.mMax, sizeof(double)), 0) << "run " << run;
            ASSERT_EQ(std::memcmp(&b.mMean, &a.mMean, sizeof(double)), 0) << "run " << run;
            ASSERT_EQ(b.mHistogram, a.mHistogram)
                << "run " << run << ", metric " << first.mMetrics[mi].first;
        }
        ASSERT_EQ(again.mCellArrays.size(), first.mCellArrays.size());
        for (std::size_t mi = 0; mi < first.mCellArrays.size(); ++mi) {
            const auto& av = first.mCellArrays[mi].second;
            const auto& bv = again.mCellArrays[mi].second;
            ASSERT_EQ(bv.size(), av.size());
            for (std::size_t bidx = 0; bidx < av.size(); ++bidx) {
                ASSERT_EQ(bv[bidx].Nbytes(), av[bidx].Nbytes());
                ASSERT_EQ(std::memcmp(bv[bidx].Data(), av[bidx].Data(), av[bidx].Nbytes()), 0)
                    << "run " << run << ", block " << bidx;
            }
        }
    }
}
