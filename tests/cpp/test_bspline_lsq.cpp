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
// The OpenCASCADE-free B-spline least-squares kernel of the skin-to-CAD fits:
// basis derivatives, evaluation, the thin-plate term, the regularised fit and
// the Newton projection (parameter correction).

#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/bspline_lsq.hpp"
#include "../../src/cpp/src/detail/nurbs_eval.hpp"

using namespace meshioplusplus::detail;

namespace {

// A bicubic surface with 6 x 5 poles over [0, 2] x [-1, 1], a wavy height field.
BsplineSurface wavy() {
    BsplineSurface s = bspline_shape(3, 3, 6, 5, 0.0, 2.0, -1.0, 1.0);
    for (std::size_t j = 0; j < s.mNv; ++j)
        for (std::size_t i = 0; i < s.mNu; ++i) {
            double* p = s.mPoles.data() + 3 * (i + s.mNu * j);
            p[0] = 2.0 * static_cast<double>(i) / 5.0;
            p[1] = -1.0 + 2.0 * static_cast<double>(j) / 4.0;
            p[2] = 0.3 * std::sin(1.7 * static_cast<double>(i)) *
                   std::cos(0.9 * static_cast<double>(j));
        }
    return s;
}

}  // namespace

TEST(BsplineLsq, ClampedKnots) {
    EXPECT_EQ(bspline_clamped_knots(3, 6, 0.0, 3.0),
              (std::vector<double>{0, 0, 0, 0, 1, 2, 3, 3, 3, 3}));
    EXPECT_THROW(bspline_clamped_knots(3, 3, 0.0, 1.0), std::invalid_argument);
    EXPECT_THROW(bspline_clamped_knots(2, 4, 1.0, 1.0), std::invalid_argument);
}

TEST(BsplineLsq, BasisDerivativesMatchFiniteDifferences) {
    const std::vector<double> k = bspline_clamped_knots(3, 7, 0.0, 1.0);
    const double u = 0.37, h = 1e-6;
    const std::size_t span = nurbs_find_span(3, k, 7, u);
    std::vector<double> d(3 * 4), dp(3 * 4), dm(3 * 4), n(4);
    bspline_basis_ders(span, u, 3, 2, k, d.data());
    bspline_basis_ders(span, u + h, 3, 2, k, dp.data());
    bspline_basis_ders(span, u - h, 3, 2, k, dm.data());
    nurbs_basis(span, u, 3, k, n.data());
    double sum = 0.0;
    for (std::size_t j = 0; j < 4; ++j) {
        EXPECT_NEAR(d[j], n[j], 1e-15);
        EXPECT_NEAR(d[4 + j], (dp[j] - dm[j]) / (2 * h), 1e-6);
        EXPECT_NEAR(d[8 + j], (dp[4 + j] - dm[4 + j]) / (2 * h), 1e-5);
        sum += d[j];
    }
    EXPECT_NEAR(sum, 1.0, 1e-15);
}

TEST(BsplineLsq, ReproducesBicubicExactly) {
    const BsplineSurface ref = wavy();
    std::vector<double> us, vs, xyz;
    for (int j = 0; j <= 24; ++j)
        for (int i = 0; i <= 30; ++i) {
            const double u = 2.0 * i / 30.0, v = -1.0 + 2.0 * j / 24.0;
            BsplineEval e;
            bspline_eval(ref, u, v, 0, e);
            us.push_back(u);
            vs.push_back(v);
            xyz.insert(xyz.end(), e.mS, e.mS + 3);
        }
    BsplineSurface fit = bspline_shape(3, 3, 6, 5, 0.0, 2.0, -1.0, 1.0);
    ASSERT_TRUE(bspline_lsq_fit(fit, us, vs, xyz, nullptr, 0.0));
    for (std::size_t i = 0; i < ref.mPoles.size(); ++i)
        EXPECT_NEAR(fit.mPoles[i], ref.mPoles[i], 1e-10);
    // Too few samples and no smoothing: singular, refused.
    BsplineSurface under = bspline_shape(3, 3, 6, 5, 0.0, 2.0, -1.0, 1.0);
    const std::vector<double> u1 = {0.5}, v1 = {0.0}, x1 = {0.0, 0.0, 1.0};
    EXPECT_FALSE(bspline_lsq_fit(under, u1, v1, x1, nullptr, 0.0));
    // The thin-plate term leaves affine maps free, so smoothing needs three
    // samples off a line, and then gives the plane through them.
    const std::vector<double> u3 = {0.5, 1.5, 1.0}, v3 = {-0.5, -0.5, 0.5},
                              x3 = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 1.0, 1.0};
    EXPECT_FALSE(bspline_lsq_fit(under, u1, v1, x1, nullptr, 1e-3));
    ASSERT_TRUE(bspline_lsq_fit(under, u3, v3, x3, nullptr, 1e-3));
    BsplineEval e;
    bspline_eval(under, 0.2, 0.7, 0, e);
    EXPECT_NEAR(e.mS[2], 1.0, 1e-9);
}

TEST(BsplineLsq, SmoothingMatrixKillsAffine) {
    // Poles at the Greville abscissae of an affine map reproduce it; its second
    // derivatives vanish, so the thin-plate energy does too.
    BsplineSurface s = bspline_shape(3, 3, 7, 6, 0.0, 1.0, 0.0, 1.0);
    SkylineMatrix k = bspline_normal_envelope(s);
    bspline_add_thin_plate(s, 1.0, k);
    const auto greville = [](const std::vector<double>& rK, std::size_t I) {
        return (rK[I + 1] + rK[I + 2] + rK[I + 3]) / 3.0;
    };
    std::vector<double> c(s.NumPoles());
    for (std::size_t j = 0; j < s.mNv; ++j)
        for (std::size_t i = 0; i < s.mNu; ++i)
            c[i + s.mNu * j] = 2.0 + 3.0 * greville(s.mKnotsU, i) - 5.0 * greville(s.mKnotsV, j);
    double energy = 0.0, scale = 0.0;
    for (std::size_t a = 0; a < c.size(); ++a)
        for (std::size_t b = 0; b < c.size(); ++b)
            if (k.Contains(a, b)) {
                energy += c[a] * k.At(a, b) * c[b];
                scale += std::abs(k.At(a, b));
            }
    EXPECT_GT(scale, 0.0);
    EXPECT_NEAR(energy, 0.0, 1e-9 * scale);
    // ...and a bump has positive energy.
    std::vector<double> bump(s.NumPoles(), 0.0);
    bump[3 + s.mNu * 3] = 1.0;
    EXPECT_GT(k.At(3 + s.mNu * 3, 3 + s.mNu * 3), 0.0);
}

TEST(BsplineLsq, NewtonProjectionConverges) {
    const BsplineSurface s = wavy();
    for (int k = 0; k < 20; ++k) {
        const double u = 0.1 + 1.8 * k / 19.0, v = -0.8 + 1.6 * ((k * 7) % 20) / 19.0;
        BsplineEval e;
        bspline_eval(s, u, v, 1, e);
        double n[3] = {e.mSu[1] * e.mSv[2] - e.mSu[2] * e.mSv[1],
                       e.mSu[2] * e.mSv[0] - e.mSu[0] * e.mSv[2],
                       e.mSu[0] * e.mSv[1] - e.mSu[1] * e.mSv[0]};
        const double ln = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        // A point 0.01 off the surface along its normal: its foot is (u, v).
        const double off = 0.01;
        const double p[3] = {e.mS[0] + off * n[0] / ln, e.mS[1] + off * n[1] / ln,
                             e.mS[2] + off * n[2] / ln};
        double pu = u + 0.05, pv = v - 0.04;
        const double d = bspline_project(s, p, pu, pv);
        EXPECT_NEAR(d, off, 1e-9);
        EXPECT_NEAR(pu, u, 1e-6);
        EXPECT_NEAR(pv, v, 1e-6);
    }
    // A point beyond the domain projects onto its border.
    const double far[3] = {-5.0, 0.0, 0.0};
    double pu = 1.0, pv = 0.0;
    bspline_project(s, far, pu, pv);
    EXPECT_EQ(pu, 0.0);
}

TEST(BsplineLsq, GaussLegendreIsExact) {
    for (int n = 1; n <= 8; ++n) {
        std::vector<double> x, w;
        gauss_legendre(n, x, w);
        // ∫ t^k over [-1, 1] = 2 / (k + 1) for even k, 0 for odd, exact to 2n - 1.
        for (int k = 0; k <= 2 * n - 1; ++k) {
            double s = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i)
                s += w[i] * std::pow(x[i], k);
            EXPECT_NEAR(s, k % 2 ? 0.0 : 2.0 / (k + 1), 1e-13) << n << " " << k;
        }
    }
}
