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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include <array>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/fem_kernels.hpp"
#include "../../src/cpp/src/detail/point_knn.hpp"
#include "../../src/cpp/src/detail/quadrature.hpp"
#include "../../src/cpp/src/detail/reference_element.hpp"

using namespace meshioplusplus;
using namespace meshioplusplus::detail;

namespace {

double integrate_rule(const QuadratureRule& rRule, const auto& rFn) {
    double s = 0.0;
    for (std::size_t i = 0; i < rRule.mPoints.size(); ++i)
        s += rRule.mWeights[i] * rFn(rRule.mPoints[i][0], rRule.mPoints[i][1], rRule.mPoints[i][2]);
    return s;
}

}  // namespace

TEST(ReferenceElement, QuadratureExactness) {
    // Line: \int_0^1 u^5 du = 1/6
    {
        const auto r = quadrature_rule(CellType::Line, 5);
        EXPECT_NEAR(integrate_rule(r, [](double u, double, double) { return std::pow(u, 5); }),
                    1.0 / 6.0, 1e-14);
    }
    // Triangle: \int_T u^2 v^2 dA = 2! 2! / (2+2+2)! = 4 / 720 = 1/180
    for (int ord : {4, 5, 6, 8}) {
        const auto r = quadrature_rule(CellType::Triangle, ord);
        EXPECT_NEAR(integrate_rule(r, [](double u, double v, double) { return u * u * v * v; }),
                    1.0 / 180.0, 1e-13)
            << "order=" << ord;
    }
    // Quad: \int_{[0,1]^2} u^3 v^2 = 1/12
    {
        const auto r = quadrature_rule(CellType::Quad, 4);
        EXPECT_NEAR(integrate_rule(r, [](double u, double v, double) { return u * u * u * v * v; }),
                    1.0 / 12.0, 1e-14);
    }
    // Tetra: \int_K u^2 v w dV = 2! 1! 1! / (2+1+1+3)! = 2 / 5040 = 1/2520
    for (int ord : {4, 6}) {
        const auto r = quadrature_rule(CellType::Tetra, ord);
        EXPECT_NEAR(integrate_rule(r, [](double u, double v, double w) { return u * u * v * w; }),
                    1.0 / 2520.0, 1e-13)
            << "order=" << ord;
    }
    // Hex: \int_{[0,1]^3} u^2 v^2 w^2 = 1/27
    {
        const auto r = quadrature_rule(CellType::Hexahedron, 4);
        EXPECT_NEAR(
            integrate_rule(r, [](double u, double v, double w) { return u * u * v * v * w * w; }),
            1.0 / 27.0, 1e-14);
    }
    // Wedge: \int u v w^2 = (1/24) * (1/3) = 1/72
    {
        const auto r = quadrature_rule(CellType::Wedge, 4);
        EXPECT_NEAR(integrate_rule(r, [](double u, double v, double w) { return u * v * w * w; }),
                    1.0 / 72.0, 1e-13);
    }
    // Pyramid: volume = 1/3, \int w dV = \int_0^1 w (1-w)^2 dw = 1/12
    {
        const auto r = quadrature_rule(CellType::Pyramid, 3);
        EXPECT_NEAR(integrate_rule(r, [](double, double, double) { return 1.0; }), 1.0 / 3.0, 1e-14);
        EXPECT_NEAR(integrate_rule(r, [](double, double, double w) { return w; }), 1.0 / 12.0,
                    1e-14);
    }
}

TEST(ReferenceElement, ShapeFunctionsAndGradients) {
    // Test Quad8 serendipity and Tetra10 quadratic reproduction + analytical derivatives.
    {
        const std::vector<Vec3> q8 = {{0, 0, 0},   {2, 0, 0},   {2, 3, 0},   {0, 3, 0},
                                      {1, 0, 0},   {2, 1.5, 0}, {1, 3, 0},   {0, 1.5, 0}};
        auto fn = [](const Vec3& p) { return 1.0 + 2.0 * p[0] - 3.0 * p[1] + p[0] * p[1] + 0.5 * p[0] * p[0]; };
        std::vector<double> vals;
        std::vector<Vec3> grads;
        double meas = 0.0;
        ASSERT_TRUE(shape_eval_and_jacobian(CellType::Quad8, q8, {0.3, 0.7, 0.0}, vals, grads, meas));
        EXPECT_NEAR(meas, 6.0, 1e-12);
        Vec3 x{0, 0, 0};
        double f = 0.0;
        Vec3 df{0, 0, 0};
        for (std::size_t i = 0; i < q8.size(); ++i) {
            const double fi = fn(q8[i]);
            f += vals[i] * fi;
            for (int d = 0; d < 3; ++d) {
                x[d] += vals[i] * q8[i][d];
                df[d] += grads[i][d] * fi;
            }
        }
        EXPECT_NEAR(f, fn(x), 1e-12);
        EXPECT_NEAR(df[0], 2.0 + x[1] + x[0], 1e-12);
        EXPECT_NEAR(df[1], -3.0 + x[0], 1e-12);
    }
    {
        const auto lat = lagrange::vtk_lattice(lagrange::Shape::Tetra, 2);
        std::vector<Vec3> t10;
        for (const auto& c : lat)
            t10.push_back({c[0] * 1.0, c[1] * 1.5, c[2] * 2.0});
        auto fn = [](const Vec3& p) {
            return 3.0 + p[0] * p[0] - 2.0 * p[1] * p[2] + 0.5 * p[2] * p[2] + p[0];
        };
        std::vector<double> vals;
        std::vector<Vec3> grads;
        double meas = 0.0;
        ASSERT_TRUE(shape_eval_and_jacobian(CellType::Tetra10, t10, {0.2, 0.3, 0.25}, vals, grads,
                                            meas));
        // Jacobian = diag(2, 3, 4), det = 24
        EXPECT_NEAR(meas, 24.0, 1e-12);
        Vec3 x{0, 0, 0};
        double f = 0.0;
        Vec3 df{0, 0, 0};
        for (std::size_t i = 0; i < t10.size(); ++i) {
            const double fi = fn(t10[i]);
            f += vals[i] * fi;
            for (int d = 0; d < 3; ++d) {
                x[d] += vals[i] * t10[i][d];
                df[d] += grads[i][d] * fi;
            }
        }
        EXPECT_NEAR(f, fn(x), 1e-12);
        EXPECT_NEAR(df[0], 2.0 * x[0] + 1.0, 1e-12);
        EXPECT_NEAR(df[1], -2.0 * x[2], 1e-12);
        EXPECT_NEAR(df[2], -2.0 * x[1] + x[2], 1e-12);
    }
}

TEST(ReferenceElement, BernsteinLagrangeRoundTrip) {
    for (lagrange::Shape s : {lagrange::Shape::Line, lagrange::Shape::Triangle,
                              lagrange::Shape::Quad, lagrange::Shape::Tetra,
                              lagrange::Shape::Hexahedron, lagrange::Shape::Wedge}) {
        for (int p : {1, 2, 3}) {
            const std::size_t n = lagrange::num_nodes(s, p);
            const auto b2l = lagrange::bezier_to_lagrange_matrix(s, p);
            const auto l2b = lagrange::lagrange_to_bezier_matrix(s, p);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < n; ++j) {
                    double dot = 0.0;
                    for (std::size_t k = 0; k < n; ++k)
                        dot += b2l[i * n + k] * l2b[k * n + j];
                    EXPECT_NEAR(dot, (i == j) ? 1.0 : 0.0, 1e-11);
                }
            }
        }
    }
}

TEST(ReferenceElement, TetP1FemKernels) {
    const Vec3 p[4] = {{0, 0, 0}, {2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
    Vec3 g[4];
    double vol = 0.0;
    ASSERT_TRUE(tet_p1_gradients(p, g, vol));
    EXPECT_NEAR(vol, 4.0, 1e-14);
    EXPECT_NEAR(g[1][0], 0.5, 1e-14);
    EXPECT_NEAR(g[2][1], 1.0 / 3.0, 1e-14);
    EXPECT_NEAR(g[3][2], 0.25, 1e-14);
    double K[4][4], M[4][4];
    tet_p1_stiffness(g, vol, K);
    tet_p1_mass_consistent(vol, M);
    double m_sum = 0.0;
    for (int i = 0; i < 4; ++i) {
        double k_row = 0.0;
        for (int j = 0; j < 4; ++j) {
            k_row += K[i][j];
            m_sum += M[i][j];
        }
        EXPECT_NEAR(k_row, 0.0, 1e-13);
    }
    EXPECT_NEAR(m_sum, vol, 1e-14);
}

TEST(ReferenceElement, PointKnnWithinRadius) {
    std::vector<Vec3> pts = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {2, 0, 0}, {0.5, 0.5, 0}};
    PointKnn knn(pts, 1.0);
    std::vector<std::pair<double, std::int64_t>> out;
    knn.WithinRadius({0, 0, 0}, 1.01, out);
    ASSERT_EQ(out.size(), 4u);
    EXPECT_EQ(out[0].second, 0);
    EXPECT_EQ(out[1].second, 1);
    EXPECT_EQ(out[2].second, 2);
    EXPECT_EQ(out[3].second, 4);
}
