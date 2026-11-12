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
// System includes
#include <cmath>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/bspline_curve.hpp"
#include "../../src/cpp/src/detail/kinematic_fit.hpp"

using namespace meshioplusplus::detail;

namespace {

constexpr double kKftPi = 3.14159265358979323846;

// (I + 1) x (J + 1) points through `rAt(u, v)`, u and v in [0, 1].
template <class F>
void kft_grid(int I, int J, F&& rAt, std::vector<PatchVec>& rPts, std::vector<PatchTri>& rTris) {
    rPts.clear();
    rTris.clear();
    for (int j = 0; j <= J; ++j)
        for (int i = 0; i <= I; ++i)
            rPts.push_back(rAt(static_cast<double>(i) / I, static_cast<double>(j) / J));
    auto id = [&](int i, int j) { return static_cast<std::size_t>(j * (I + 1) + i); };
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < I; ++i) {
            rTris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            rTris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
}

double kft_lobe(double A) {
    return 1.0 + 0.2 * std::cos(3.0 * A) + 0.08 * std::sin(2.0 * A);
}

double kft_wave(double Z) {
    return 1.0 + 0.15 * std::sin(2.0 * kKftPi * Z);
}

KinematicOptions kft_options(double Tol) {
    KinematicOptions o;
    o.mTolerance = Tol;
    return o;
}

double kft_axis_error(const PatchVec& rA, const PatchVec& rB) {
    const double d = std::abs(rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2]);
    return 1.0 - d;
}

}  // namespace

TEST(KinematicFit, HalfALobedPrismIsAnExtrusion) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    kft_grid(
        48, 12,
        [](double U, double V) {
            const double a = kKftPi * U;
            return PatchVec{kft_lobe(a) * std::cos(a), kft_lobe(a) * std::sin(a), 1.5 * V};
        },
        pts, tris);
    const KinematicFit f = kinematic_fit(pts, tris, kft_options(6e-3));
    ASSERT_TRUE(f.mOk) << f.mWhy;
    EXPECT_EQ(f.mKind, KinematicKind::Extrusion);
    EXPECT_LT(kft_axis_error(f.mAxis, {0.0, 0.0, 1.0}), 1e-9);
    EXPECT_LE(f.mMaxDist, 3e-3);
    // The exact surface and its B-spline agree.
    for (double t : {0.1, 0.5, 0.9})
        for (double s : {0.0, 0.3}) {
            const double tt = f.mProfile.T0() + t * (f.mProfile.T1() - f.mProfile.T0());
            const PatchVec p = f.At(tt, s);
            BsplineEval e;
            bspline_eval(f.mSurface, tt, s, 0, e);
            EXPECT_NEAR(e.mS[0], p[0], 1e-12);
            EXPECT_NEAR(e.mS[1], p[1], 1e-12);
            EXPECT_NEAR(e.mS[2], p[2], 1e-12);
        }
}

TEST(KinematicFit, HalfAWavyWallIsARevolution) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    kft_grid(
        48, 32,
        [](double U, double V) {
            const double a = kKftPi * U, z = 2.0 * V;
            return PatchVec{kft_wave(z) * std::cos(a), kft_wave(z) * std::sin(a), z};
        },
        pts, tris);
    const KinematicFit f = kinematic_fit(pts, tris, kft_options(6e-3));
    ASSERT_TRUE(f.mOk) << f.mWhy;
    EXPECT_EQ(f.mKind, KinematicKind::Revolution);
    EXPECT_LT(kft_axis_error(f.mAxis, {0.0, 0.0, 1.0}), 1e-6);
    // The axis passes through the z axis.
    EXPECT_NEAR(std::hypot(f.mOrigin[0], f.mOrigin[1]), 0.0, 2e-3);
    EXPECT_LE(f.mMaxDist, 3e-3);
}

TEST(KinematicFit, AnAdvancingScrewIsHelical) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    // A rounded profile r(s) = 1 + 0.1 s^2 (s in [-1, 1]) at z = 0.3 s,
    // screwed with a pitch of 0.8 a turn through three quarters of a turn.
    kft_grid(
        16, 64,
        [](double U, double V) {
            const double s = 2.0 * U - 1.0, a = 1.5 * kKftPi * V;
            const double r = 1.0 + 0.1 * s * s, z = 0.3 * s + 0.8 * a / (2.0 * kKftPi);
            return PatchVec{r * std::cos(a), r * std::sin(a), z};
        },
        pts, tris);
    const KinematicFit f = kinematic_fit(pts, tris, kft_options(4e-3));
    ASSERT_TRUE(f.mOk) << f.mWhy;
    EXPECT_EQ(f.mKind, KinematicKind::Helical);
    EXPECT_LT(kft_axis_error(f.mAxis, {0.0, 0.0, 1.0}), 1e-4);
}

TEST(KinematicFit, PrimitivesAreRefused) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    const auto refused = [&](const char* pName) {
        const KinematicFit f = kinematic_fit(pts, tris, kft_options(4e-3));
        EXPECT_FALSE(f.mOk) << pName << " fitted as " << static_cast<int>(f.mKind);
    };
    kft_grid(16, 16, [](double U, double V) { return PatchVec{U, V, 0.0}; }, pts, tris);
    refused("plane");
    kft_grid(
        32, 16,
        [](double U, double V) { return PatchVec{std::cos(kKftPi * U), std::sin(kKftPi * U), V}; },
        pts, tris);
    refused("cylinder");
    kft_grid(
        32, 16,
        [](double U, double V) {
            const double r = 1.0 + 0.4 * V;
            return PatchVec{r * std::cos(kKftPi * U), r * std::sin(kKftPi * U), V};
        },
        pts, tris);
    refused("cone");
    kft_grid(
        32, 16,
        [](double U, double V) {
            const double th = 0.3 + 1.2 * V, ph = kKftPi * U;
            return PatchVec{std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th)};
        },
        pts, tris);
    refused("sphere");
}

TEST(KinematicFit, BandsAreTold) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    const int n = 32, m = 4;
    for (int j = 0; j <= m; ++j)
        for (int i = 0; i < n; ++i) {
            const double a = 2.0 * kKftPi * i / n;
            pts.push_back({std::cos(a), std::sin(a), static_cast<double>(j) / m});
        }
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) {
            const std::size_t a = static_cast<std::size_t>(j * n + i),
                              b = static_cast<std::size_t>(j * n + (i + 1) % n);
            tris.push_back({a, b, b + n});
            tris.push_back({a, b + n, a + n});
        }
    EXPECT_TRUE(kinematic_is_band(tris));
    EXPECT_FALSE(kinematic_fit(pts, tris, kft_options(1e-3)).mOk);
}

TEST(KinematicFit, Reproducible) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    kft_grid(
        48, 32,
        [](double U, double V) {
            const double a = kKftPi * U, z = 2.0 * V;
            return PatchVec{kft_wave(z) * std::cos(a), kft_wave(z) * std::sin(a), z};
        },
        pts, tris);
    const KinematicFit a = kinematic_fit(pts, tris, kft_options(6e-3));
    const KinematicFit b = kinematic_fit(pts, tris, kft_options(6e-3));
    ASSERT_TRUE(a.mOk);
    EXPECT_EQ(a.mSurface.mPoles, b.mSurface.mPoles);
    EXPECT_EQ(a.mU, b.mU);
}

TEST(BsplineCurve, ProjectAndSecondDerivative) {
    std::vector<PatchVec> poly;
    for (int i = 0; i <= 40; ++i) {
        const double t = i / 40.0;
        poly.push_back({t, 0.3 * std::sin(3.0 * t), 0.1 * t * t});
    }
    CurveFitOptions co;
    co.mTolerance = 1e-5;
    const CurveFitResult fit = bspline_curve_fit(poly, co);
    ASSERT_TRUE(fit.mOk);
    const BsplineCurve& c = fit.mCurve;
    const double h = 1e-5;
    for (double f : {0.2, 0.5, 0.8}) {
        const double t = c.T0() + f * (c.T1() - c.T0());
        PatchVec d1, d2, d1p, d2p, d1m, d2m;
        bspline_curve_eval2(c, t, d1, d2);
        bspline_curve_eval2(c, t + h, d1p, d2p);
        bspline_curve_eval2(c, t - h, d1m, d2m);
        for (std::size_t k = 0; k < 3; ++k)
            EXPECT_NEAR(d2[k], (d1p[k] - d1m[k]) / (2.0 * h), 1e-4 * (1.0 + std::abs(d2[k])));
        // A point off the curve along its normal projects back to t.
        const PatchVec p = bspline_curve_eval(c, t);
        const PatchVec off{p[0] - 0.01 * d1[1], p[1] + 0.01 * d1[0], p[2]};
        double tt = 0.0;
        const double d = bspline_curve_project(c, off, tt);
        EXPECT_NEAR(tt, t, 1e-7);
        EXPECT_NEAR(d, 0.01 * std::hypot(d1[0], d1[1]), 1e-8);
    }
}
