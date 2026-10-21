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
// The OpenCASCADE-free half of fit_cad's general freeform faces: disc
// parametrisation (least-squares conformal maps, mean-value coordinates) and
// the iterative B-spline surface fit.

#include <cmath>
#include <cstddef>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/bspline_fit.hpp"
#include "../../src/cpp/src/detail/param_disc.hpp"

using namespace meshioplusplus::detail;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Grid {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
    std::vector<double> mU, mV;  // the grid's own parameters
};

// An (I + 1) x (J + 1) grid through `rAt(u, v)`, u, v in [0, 1], two triangles
// per cell; with `Wrap`, u is periodic (a closed band).
template <class F>
Grid grid(int I, int J, bool Wrap, F&& rAt) {
    Grid g;
    const int ni = Wrap ? I : I + 1;
    for (int j = 0; j <= J; ++j)
        for (int i = 0; i < ni; ++i) {
            const double u = static_cast<double>(i) / I, v = static_cast<double>(j) / J;
            g.mPts.push_back(rAt(u, v));
            g.mU.push_back(u);
            g.mV.push_back(v);
        }
    const auto id = [&](int i, int j) { return static_cast<std::size_t>(j * ni + (i % ni)); };
    for (int j = 0; j < J; ++j)
        for (int i = 0; i < I; ++i) {
            g.mTris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            g.mTris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
    return g;
}

Grid saddle(int N) {
    return grid(N, N, false, [](double u, double v) {
        const double x = 2.0 * u - 1.0, y = 2.0 * v - 1.0;
        return PatchVec{x, y, 0.3 * (x * x - y * y)};
    });
}

Grid wavy(int N) {
    return grid(N, N, false, [](double u, double v) {
        return PatchVec{u, v, 0.05 * std::sin(4.0 * kPi * u) * std::sin(3.0 * kPi * v)};
    });
}

}  // namespace

TEST(ParamDisc, RecognisesDiscs) {
    EXPECT_TRUE(param_is_disc(saddle(6).mTris));
    // A band closed around an axis has two boundary loops.
    const Grid band = grid(16, 4, true, [](double u, double v) {
        return PatchVec{std::cos(2 * kPi * u), std::sin(2 * kPi * u), v};
    });
    bool ok = true;
    EXPECT_EQ(param_boundary_loops(band.mTris, ok).size(), 2u);
    EXPECT_TRUE(ok);
    EXPECT_FALSE(param_is_disc(band.mTris));
    EXPECT_FALSE(param_disc(band.mPts, band.mTris).mOk);
    EXPECT_FALSE(param_is_disc({}));
}

TEST(ParamDisc, LscmIsAnIsometryOnAPlane) {
    const Grid g =
        grid(8, 5, false, [](double u, double v) { return PatchVec{3.0 * u, 2.0 * v, 0.0}; });
    const DiscParam p = param_lscm(g.mPts, g.mTris);
    ASSERT_TRUE(p.mOk);
    EXPECT_EQ(p.mNumFlips, 0u);
    // A flat patch maps onto itself up to a rigid motion: lengths are kept.
    for (const PatchTri& t : g.mTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)],
                              b = t[static_cast<std::size_t>((k + 1) % 3)];
            const double d3 = std::hypot(g.mPts[a][0] - g.mPts[b][0], g.mPts[a][1] - g.mPts[b][1]);
            const double d2 = std::hypot(p.mU[a] - p.mU[b], p.mV[a] - p.mV[b]);
            EXPECT_NEAR(d2, d3, 1e-8);
        }
}

TEST(ParamDisc, LscmFlattensACurvedDiscWithoutFlips) {
    const Grid g = saddle(20);
    const DiscParam p = param_lscm(g.mPts, g.mTris);
    ASSERT_TRUE(p.mOk);
    EXPECT_EQ(p.mNumFlips, 0u);
    const DiscParam d = param_disc(g.mPts, g.mTris);
    EXPECT_TRUE(d.mOk);
    // LSCM is taken first, turned onto its principal axes: the same map up
    // to a rigid motion, so every edge keeps its length.
    EXPECT_EQ(d.mNumFlips, 0u);
    for (const PatchTri& t : g.mTris) {
        const double lp = std::hypot(p.mU[t[0]] - p.mU[t[1]], p.mV[t[0]] - p.mV[t[1]]);
        const double ld = std::hypot(d.mU[t[0]] - d.mU[t[1]], d.mV[t[0]] - d.mV[t[1]]);
        EXPECT_NEAR(ld, lp, 1e-9 * (1.0 + lp));
    }
}

TEST(ParamDisc, MeanValueIsBijectiveOnACircle) {
    // A tall bump: steep sides, still a disc.
    const Grid g = grid(16, 16, false, [](double u, double v) {
        const double x = 2.0 * u - 1.0, y = 2.0 * v - 1.0;
        return PatchVec{x, y, 1.5 * std::exp(-4.0 * (x * x + y * y))};
    });
    const DiscParam p = param_mean_value(g.mPts, g.mTris);
    ASSERT_TRUE(p.mOk);
    EXPECT_EQ(p.mNumFlips, 0u);
    // The boundary lies on the circle.
    bool ok = true;
    const auto loops = param_boundary_loops(g.mTris, ok);
    ASSERT_EQ(loops.size(), 1u);
    const double r0 = std::hypot(p.mU[loops[0][0]], p.mV[loops[0][0]]);
    for (std::size_t i : loops[0])
        EXPECT_NEAR(std::hypot(p.mU[i], p.mV[i]), r0, 1e-9 * r0);
}

TEST(ParamDisc, DetectsFlips) {
    const Grid g = saddle(6);
    std::vector<double> u = g.mU, v = g.mV;
    EXPECT_EQ(param_flips(g.mTris, u, v), 0u);
    // Move an interior vertex across its neighbours.
    const std::size_t c = 3 * 7 + 3;
    u[c] += 0.4;
    v[c] += 0.4;
    EXPECT_GT(param_flips(g.mTris, u, v), 0u);
}

TEST(BsplineFit, ReproducesABicubic) {
    // Points on a polynomial surface of degree 3 in each parameter: a cubic
    // B-spline with no interior knot holds it exactly.
    const Grid g = grid(12, 12, false, [](double u, double v) {
        return PatchVec{u, v, 0.2 * u * u * u - 0.3 * u * v * v + 0.1 * v};
    });
    BsplineFitOptions o;
    o.mTolerance = 1e-9;
    o.mSmoothing = 0.0;
    o.mStartPoles = 4;
    o.mMargin = 0.0;
    o.mCentroidWeight = 0.0;
    o.mCentroidTolerance = 1.0;  // a flat facet's centroid is off any smooth surface
    const BsplineFitResult r = bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o);
    EXPECT_TRUE(r.mOk);
    EXPECT_LE(r.mMaxVertexError, 1e-9);
    EXPECT_EQ(r.mSurface.mNu, 4u);
}

TEST(BsplineFit, FitsASaddleOverItsConformalMap) {
    const Grid g = saddle(24);
    const DiscParam p = param_disc(g.mPts, g.mTris);
    ASSERT_TRUE(p.mOk);
    BsplineFitOptions o;
    o.mTolerance = 1e-4;
    const BsplineFitResult r = bspline_fit_patch(g.mPts, g.mTris, p.mU, p.mV, o);
    EXPECT_TRUE(r.mOk);
    EXPECT_LE(r.mMaxVertexError, 1e-4);
    EXPECT_LE(r.mMaxCentroidError, 2e-4);
    for (std::size_t i = 1; i < r.mHistory.size(); ++i)
        EXPECT_LE(r.mHistory[i], r.mHistory[i - 1]);
}

TEST(BsplineFit, KnotInsertionReachesTheTolerance) {
    const Grid g = wavy(40);
    BsplineFitOptions o;
    o.mTolerance = 5e-4;
    o.mStartPoles = 4;
    const BsplineFitResult r = bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o);
    EXPECT_TRUE(r.mOk);
    EXPECT_GT(r.mRefinements, 0);
    EXPECT_GT(r.mSurface.mNu, 4u);
    EXPECT_LE(r.mMaxVertexError, 5e-4);
    ASSERT_GE(r.mHistory.size(), 2u);
    for (std::size_t i = 1; i < r.mHistory.size(); ++i)
        EXPECT_LE(r.mHistory[i], r.mHistory[i - 1]);
    EXPECT_LT(r.mHistory.back(), r.mHistory.front());
}

TEST(BsplineFit, RespectsThePoleBudget) {
    const Grid g = wavy(40);
    BsplineFitOptions o;
    o.mTolerance = 1e-6;
    o.mStartPoles = 4;
    o.mMaxPoles = 6;
    const BsplineFitResult r = bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o);
    EXPECT_FALSE(r.mOk);
    EXPECT_LE(r.mSurface.mNu, 6u);
    EXPECT_LE(r.mSurface.mNv, 6u);
    EXPECT_LE(r.mIterations, o.mMaxIterations);
}

TEST(BsplineFit, IsDeterministic) {
    const Grid g = wavy(20);
    BsplineFitOptions o;
    o.mTolerance = 5e-4;
    const BsplineFitResult a = bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o);
    const BsplineFitResult b = bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o);
    EXPECT_EQ(a.mSurface.mPoles, b.mSurface.mPoles);
    EXPECT_EQ(a.mSurface.mKnotsU, b.mSurface.mKnotsU);
    EXPECT_EQ(a.mHistory, b.mHistory);
}

TEST(BsplineFit, RefusesDegenerateInput) {
    BsplineFitOptions o;
    o.mTolerance = 1e-3;
    EXPECT_FALSE(bspline_fit_patch({}, {}, {}, {}, o).mOk);
    const Grid g = saddle(4);
    std::vector<double> zero(g.mPts.size(), 0.0);
    EXPECT_FALSE(bspline_fit_patch(g.mPts, g.mTris, zero, zero, o).mOk);
    o.mTolerance = 0.0;
    EXPECT_FALSE(bspline_fit_patch(g.mPts, g.mTris, g.mU, g.mV, o).mOk);
}
