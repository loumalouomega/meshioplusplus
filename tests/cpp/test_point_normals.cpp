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
// Point-cloud normals by local plane fits (roadmap 7.5.2).

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/point_normals.hpp"

using namespace meshioplusplus;

namespace {

Mesh pn_cloud(const std::vector<std::array<double, 3>>& rPts) {
    Mesh m;
    NDArray p(DType::Float64, {rPts.size(), std::size_t{3}});
    for (std::size_t i = 0; i < rPts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[i * 3 + k] = rPts[i][k];
    m.AssignPoints(std::move(p));
    return m;
}

std::vector<std::array<double, 3>> pn_sphere(std::size_t N, double Radius, double Cx) {
    std::mt19937_64 rng(9);
    std::normal_distribution<double> g;
    std::vector<std::array<double, 3>> out(N);
    for (auto& q : out) {
        double x = g(rng), y = g(rng), z = g(rng);
        const double r = std::sqrt(x * x + y * y + z * z);
        q = {Cx + Radius * x / r, Radius * y / r, Radius * z / r};
    }
    return out;
}

}  // namespace

TEST(PointNormals, SphereNormalsAreRadialAndOutward) {
    const auto pts = pn_sphere(3000, 2.0, 0.0);
    PointNormalsOptions o;
    o.mRecordVariation = true;
    const PointNormalsResult r = estimate_normals(pn_cloud(pts), o);
    EXPECT_EQ(r.mNumDegenerate, 0);
    EXPECT_EQ(r.mNumComponents, 1);
    const double* n = r.mMesh.PointData("normals").As<double>();
    double worst = 1.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double d =
            (n[3 * i] * pts[i][0] + n[3 * i + 1] * pts[i][1] + n[3 * i + 2] * pts[i][2]) / 2.0;
        worst = std::min(worst, d);
    }
    EXPECT_GT(worst, 0.99);  // within about 8 degrees of radial, never inward
    EXPECT_TRUE(r.mMesh.HasPointData(kNormalsVariationName));
}

TEST(PointNormals, TwoSpheresAreTwoPieces) {
    auto a = pn_sphere(800, 1.0, 0.0);
    const auto b = pn_sphere(800, 1.0, 10.0);
    a.insert(a.end(), b.begin(), b.end());
    const PointNormalsResult r = estimate_normals(pn_cloud(a));
    EXPECT_EQ(r.mNumComponents, 2);
    const double* n = r.mMesh.PointData("normals").As<double>();
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double cx = i < 800 ? 0.0 : 10.0;
        EXPECT_GT(n[3 * i] * (a[i][0] - cx) + n[3 * i + 1] * a[i][1] + n[3 * i + 2] * a[i][2], 0.9);
    }
}

TEST(PointNormals, PlaneFacesUpAndViewpointFlips) {
    std::vector<std::array<double, 3>> pts;
    for (int j = 0; j < 20; ++j)
        for (int i = 0; i < 20; ++i)
            pts.push_back({0.1 * i, 0.1 * j, 0.0});
    const PointNormalsResult up = estimate_normals(pn_cloud(pts));
    const double* n = up.mMesh.PointData("normals").As<double>();
    for (std::size_t i = 0; i < pts.size(); ++i)
        EXPECT_DOUBLE_EQ(n[3 * i + 2], 1.0);
    PointNormalsOptions o;
    o.mOrientation = NormalOrientation::Viewpoint;
    o.mViewpoint = {1.0, 1.0, -5.0};
    const PointNormalsResult down = estimate_normals(pn_cloud(pts), o);
    const double* d = down.mMesh.PointData("normals").As<double>();
    for (std::size_t i = 0; i < pts.size(); ++i)
        EXPECT_DOUBLE_EQ(d[3 * i + 2], -1.0);
}

TEST(PointNormals, LineIsDegenerateAndRefusals) {
    std::vector<std::array<double, 3>> line;
    for (int i = 0; i < 10; ++i)
        line.push_back({double(i), 2.0 * i, 0.0});
    const PointNormalsResult r = estimate_normals(pn_cloud(line));
    EXPECT_EQ(r.mNumDegenerate, 10);
    EXPECT_TRUE(std::isnan(r.mMesh.PointData("normals").As<double>()[0]));
    EXPECT_THROW(estimate_normals(pn_cloud({{0, 0, 0}, {1, 0, 0}})), std::invalid_argument);
    PointNormalsOptions o;
    o.mK = 1;
    EXPECT_THROW(estimate_normals(pn_cloud(line), o), std::invalid_argument);
    EXPECT_THROW(normal_orientation_from_name("inward"), std::invalid_argument);
}
