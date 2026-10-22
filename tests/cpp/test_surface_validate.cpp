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
// detail/surface_validate (fit_cad's freeform modes, v17.12.0): the noise
// floor of a patch, and a fitted surface sampled against its patch.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "analytic_surfaces.hpp"

#include "../../src/cpp/src/detail/bspline_fit.hpp"
#include "../../src/cpp/src/detail/param_disc.hpp"
#include "../../src/cpp/src/detail/surface_validate.hpp"

using namespace meshioplusplus;
using detail::PatchTri;
using detail::PatchVec;

namespace {

// An n x n grid over the unit square at height `F(x, y)`, plus noise of
// deviation `Sigma` in z.
template <class Height>
void grid(int n, Height F, double Sigma, std::vector<PatchVec>& rPts,
          std::vector<PatchTri>& rTris) {
    std::uint64_t state = 3;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const double x = static_cast<double>(i) / n, y = static_cast<double>(j) / n;
            double g = -6.0;
            for (int k = 0; k < 12; ++k)
                g += static_cast<double>(analytic::splitmix64(state) >> 11) * 0x1.0p-53;
            rPts.push_back({x, y, F(x, y) + Sigma * g});
        }
    const auto id = [n](int i, int j) { return static_cast<std::size_t>(j * (n + 1) + i); };
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            rTris.push_back({id(i, j), id(i + 1, j), id(i + 1, j + 1)});
            rTris.push_back({id(i, j), id(i + 1, j + 1), id(i, j + 1)});
        }
}

}  // namespace

TEST(SurfaceValidate, NoiseFloorOfAPlaneIsZero) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(20, [](double, double) { return 0.0; }, 0.0, pts, tris);
    EXPECT_LT(detail::patch_noise_sigma(pts, tris), 1e-12);
}

TEST(SurfaceValidate, NoiseFloorFindsTheNoise) {
    for (double sigma : {1e-4, 1e-3}) {
        std::vector<PatchVec> pts;
        std::vector<PatchTri> tris;
        // A curved, finely tessellated surface: the curvature must not read as noise.
        grid(40, [](double x, double y) { return 0.2 * (x * x - y * y); }, sigma, pts, tris);
        const double est = detail::patch_noise_sigma(pts, tris);
        EXPECT_GT(est, 0.6 * sigma) << sigma;
        EXPECT_LT(est, 1.6 * sigma) << sigma;
    }
}

// A fitted surface measured both ways: its samples stay as near the patch as
// its fit says the vertices are, and the reverse distance is the fit's own.
TEST(SurfaceValidate, FittedSurfaceIsWithinItsTolerance) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    grid(
        24,
        [](double x, double y) {
            return 0.3 * std::exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / 0.05);
        },
        0.0, pts, tris);
    const detail::DiscParam par = detail::param_ladder(pts, tris, detail::ParamLadder::Auto);
    ASSERT_TRUE(par.mOk);
    detail::BsplineFitOptions o;
    o.mTolerance = 2e-3;
    // A centroid lies its chordal sag (about 2.6e-3 at the peak) off the surface.
    o.mCentroidTolerance = 1e-2;
    const detail::BsplineFitResult fit = detail::bspline_fit_patch(pts, tris, par.mU, par.mV, o);
    ASSERT_TRUE(fit.mOk) << fit.mMaxVertexError << " " << fit.mMaxCentroidError << " "
                         << fit.mSurface.mNu << "x" << fit.mSurface.mNv << " " << fit.mIterations;
    const detail::SurfaceValidation v =
        detail::surface_validate(fit.mSurface, pts, tris, fit.mU, fit.mV);
    EXPECT_GT(v.mNumSamples, fit.mSurface.NumPoles());
    EXPECT_LE(v.mMaxReverse, fit.mMaxVertexError * 1.0001 + 1e-12);
    EXPECT_LT(v.mMaxForward, 8e-3);
    EXPECT_LE(v.mRmsForward, v.mMaxForward);
    // Deterministic.
    const detail::SurfaceValidation w =
        detail::surface_validate(fit.mSurface, pts, tris, fit.mU, fit.mV);
    EXPECT_EQ(v.mMaxForward, w.mMaxForward);
    EXPECT_EQ(v.mMaxReverse, w.mMaxReverse);
}
