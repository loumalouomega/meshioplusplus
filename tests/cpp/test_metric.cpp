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
// The metric toolbox: intersection, normalisation, clamping, gradation
// (roadmap 7.4.3).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/metric.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// An m x m unit square of triangles carrying 2-D metrics "a" and "b".
Mesh mt_square(int M, double A11, double A22, double B11, double B22) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= M; ++j)
        for (int i = 0; i <= M; ++i)
            pts.push_back({double(i) / M, double(j) / M, 0.0});
    std::vector<std::vector<std::int64_t>> tris;
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i) {
            const std::int64_t a = j * (M + 1) + i, b = a + 1, c = a + M + 2, d = a + M + 1;
            tris.push_back({a, b, c});
            tris.push_back({a, c, d});
        }
    Mesh mesh = mt::make_mesh(pts, "triangle", tris);
    NDArray a(DType::Float64, {pts.size(), std::size_t{3}});
    NDArray b(DType::Float64, {pts.size(), std::size_t{3}});
    for (std::size_t p = 0; p < pts.size(); ++p) {
        double* pa = a.As<double>() + 3 * p;
        double* pb = b.As<double>() + 3 * p;
        pa[0] = A11, pa[1] = 0.0, pa[2] = A22;
        pb[0] = B11, pb[1] = 0.0, pb[2] = B22;
    }
    mesh.AddPointData("a", std::move(a));
    mesh.AddPointData("b", std::move(b));
    return mesh;
}

}  // namespace

TEST(Metric, IntersectionOfCrossedMetricsTakesTheFinerEachWay) {
    PrepareMetricOptions o;
    o.mArrays = {"a", "b"};
    const PrepareMetricResult r = prepare_metric(mt_square(4, 100.0, 1.0, 1.0, 100.0), o);
    const double* m = r.mMesh.PointData("metric").As<double>();
    EXPECT_NEAR(m[0], 100.0, 1e-9);
    EXPECT_NEAR(m[1], 0.0, 1e-9);
    EXPECT_NEAR(m[2], 100.0, 1e-9);
    // A metric intersected with itself is itself.
    o.mArrays = {"a", "a"};
    const PrepareMetricResult self = prepare_metric(mt_square(2, 9.0, 4.0, 1.0, 1.0), o);
    const double* same = self.mMesh.PointData("metric").As<double>();
    EXPECT_NEAR(same[0], 9.0, 1e-12);
    EXPECT_NEAR(same[2], 4.0, 1e-12);
}

TEST(Metric, NormalisationReachesTheTargetComplexity) {
    PrepareMetricOptions o;
    o.mArrays = {"a"};
    o.mTargetVertices = 5000.0;
    for (const double p : {0.0, 1.0, 2.0}) {
        o.mNorm = p;
        const PrepareMetricResult r = prepare_metric(mt_square(8, 400.0, 25.0, 1.0, 1.0), o);
        EXPECT_NEAR(r.mComplexityBefore, std::sqrt(400.0 * 25.0), 1e-9);
        EXPECT_NEAR(r.mComplexityAfter, 5000.0, 1e-6) << p;
        // a uniform metric keeps its shape: the anisotropy ratio is unchanged
        const double* m = r.mMesh.PointData("metric").As<double>();
        EXPECT_NEAR(m[0] / m[2], 16.0, 1e-9);
    }
}

TEST(Metric, ClampAndGradation) {
    Mesh mesh = mt_square(10, 1.0, 1.0, 1.0, 1.0);
    // one point asks for a size of 0.01
    double* a = const_cast<double*>(mesh.PointData("a").As<double>());
    a[3 * 60] = a[3 * 60 + 2] = 1e4;
    PrepareMetricOptions o;
    o.mArrays = {"a"};
    o.mHmin = 0.02;
    o.mHmax = 0.5;
    const PrepareMetricResult clamp = prepare_metric(mesh, o);
    const double* c = clamp.mMesh.PointData("metric").As<double>();
    EXPECT_NEAR(c[3 * 60], 1.0 / (0.02 * 0.02), 1e-6);  // the finest size, clamped
    EXPECT_NEAR(c[0], 4.0, 1e-9);                       // 1 / hmax^2
    o.mGradation = 1.5;
    const PrepareMetricResult graded = prepare_metric(mesh, o);
    EXPECT_GT(graded.mNumGraded, 0);
    EXPECT_GT(graded.mSweeps, 0);
    const double* g = graded.mMesh.PointData("metric").As<double>();
    // the neighbour at 0.1 now asks for less than 0.5 but more than 0.02
    const double h_next = 1.0 / std::sqrt(g[3 * 61]);
    EXPECT_LT(h_next, 0.5);
    EXPECT_GT(h_next, 0.02);
}

TEST(Metric, Refusals) {
    Mesh mesh = mt_square(2, 1.0, 1.0, 1.0, 1.0);
    PrepareMetricOptions o;
    EXPECT_THROW(prepare_metric(mesh, o), std::invalid_argument);
    o.mArrays = {"missing"};
    EXPECT_THROW(prepare_metric(mesh, o), std::invalid_argument);
    o.mArrays = {"a"};
    o.mNorm = 0.5;
    EXPECT_THROW(prepare_metric(mesh, o), std::invalid_argument);
    o.mNorm = 2.0;
    o.mGradation = 1.0;
    EXPECT_THROW(prepare_metric(mesh, o), std::invalid_argument);
    o.mGradation = 0.0;
    o.mHmin = 1.0;
    o.mHmax = 0.5;
    EXPECT_THROW(prepare_metric(mesh, o), std::invalid_argument);
}
