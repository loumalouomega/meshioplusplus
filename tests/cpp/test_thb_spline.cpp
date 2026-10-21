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

// Truncated hierarchical B-splines (detail/thb_spline).

// System includes
#include <cmath>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/bspline_fit.hpp"
#include "../../src/cpp/src/detail/dense_chol.hpp"
#include "../../src/cpp/src/detail/thb_spline.hpp"

namespace {

using meshioplusplus::detail::BsplineFitOptions;
using meshioplusplus::detail::PatchTri;
using meshioplusplus::detail::PatchVec;
using meshioplusplus::detail::ThbBasis;

// Every finest coefficient summed over the basis: one where the basis is a
// partition of unity.
std::vector<double> coverage(const ThbBasis& rB) {
    const std::size_t n = rB.FineSize();
    std::vector<double> sum(n * n, 0.0);
    for (const auto& f : rB.mFunctions)
        for (const auto& c : f.mFine)
            sum[c.first] += c.second;
    return sum;
}

// A grid patch over [0, 1]^2 on z = Height(x, y), its parameters (x, y).
void grid_patch(std::size_t N, double (*Height)(double, double), std::vector<PatchVec>& rPts,
                std::vector<PatchTri>& rTris, std::vector<double>& rU, std::vector<double>& rV) {
    for (std::size_t j = 0; j <= N; ++j)
        for (std::size_t i = 0; i <= N; ++i) {
            const double x = static_cast<double>(i) / N, y = static_cast<double>(j) / N;
            rPts.push_back({x, y, Height(x, y)});
            rU.push_back(x);
            rV.push_back(y);
        }
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t a = i + (N + 1) * j, b = a + 1, c = a + N + 1, d = c + 1;
            rTris.push_back({a, b, d});
            rTris.push_back({a, d, c});
        }
}

double bump(double X, double Y) {
    return 0.15 * std::exp(-((X - 0.3) * (X - 0.3) + (Y - 0.35) * (Y - 0.35)) / 0.004);
}

}  // namespace

TEST(ThbSpline, UnrefinedIsTheTensorProduct) {
    ThbBasis b = meshioplusplus::detail::thb_hierarchy(3, 4, 3);
    meshioplusplus::detail::thb_build(b);
    EXPECT_EQ(b.NumLevels(), 1);
    ASSERT_EQ(b.mFunctions.size(), 49u);
    for (const auto& f : b.mFunctions) {
        ASSERT_EQ(f.mFine.size(), 1u);
        EXPECT_EQ(f.mFine[0].second, 1.0);
    }
}

TEST(ThbSpline, RefinedEverywhereIsTheFinerTensorProduct) {
    ThbBasis b = meshioplusplus::detail::thb_hierarchy(3, 3, 2);
    std::fill(b.mRefined[0].begin(), b.mRefined[0].end(), 1);
    meshioplusplus::detail::thb_build(b);
    EXPECT_EQ(b.mFunctions.size(), 81u);  // (2 * 3 + 3)^2, all of level 1
    for (const auto& f : b.mFunctions)
        EXPECT_EQ(f.mLevel, 1);
}

TEST(ThbSpline, TruncatedBasisIsAPartitionOfUnityAndIndependent) {
    // Three levels: a corner refined, and inside it a smaller corner again.
    ThbBasis b = meshioplusplus::detail::thb_hierarchy(3, 4, 3);
    for (std::size_t j = 0; j < 2; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            b.mRefined[0][i + 4 * j] = 1;
    for (std::size_t j = 0; j < 2; ++j)
        for (std::size_t i = 0; i < 3; ++i)
            b.mRefined[1][i + 8 * j] = 1;
    meshioplusplus::detail::thb_build(b);
    EXPECT_EQ(b.NumLevels(), 3);
    const std::size_t n = b.FineSize();
    EXPECT_EQ(n, 19u);
    EXPECT_GT(b.mFunctions.size(), 49u);
    EXPECT_LT(b.mFunctions.size(), n * n);
    for (double c : coverage(b))
        EXPECT_NEAR(c, 1.0, 1e-13);
    // The Gram matrix of their fine coefficients is positive definite.
    const std::size_t nd = b.mFunctions.size();
    std::vector<double> g(nd * nd, 0.0);
    std::vector<std::vector<double>> dense(nd, std::vector<double>(n * n, 0.0));
    for (std::size_t k = 0; k < nd; ++k)
        for (const auto& c : b.mFunctions[k].mFine)
            dense[k][c.first] = c.second;
    for (std::size_t a = 0; a < nd; ++a)
        for (std::size_t c = 0; c < nd; ++c) {
            double s = 0.0;
            for (std::size_t f = 0; f < n * n; ++f)
                s += dense[a][f] * dense[c][f];
            g[a * nd + c] = s;
        }
    EXPECT_TRUE(meshioplusplus::detail::chol_factor(g, nd));
}

TEST(ThbSpline, LocalBumpTakesFewerDegreesOfFreedom) {
    std::vector<PatchVec> pts;
    std::vector<PatchTri> tris;
    std::vector<double> u, v;
    grid_patch(64, bump, pts, tris, u, v);
    BsplineFitOptions o;
    o.mTolerance = 2e-3;
    o.mMaxPoles = 64;
    o.mMaxIterations = 40;
    meshioplusplus::detail::ThbFitInfo info;
    const auto thb = meshioplusplus::detail::thb_fit_patch(pts, tris, u, v, o, &info);
    const auto knots = meshioplusplus::detail::bspline_fit_patch(pts, tris, u, v, o);
    ASSERT_TRUE(thb.mOk);
    ASSERT_TRUE(knots.mOk);
    EXPECT_LE(thb.mMaxVertexError, o.mTolerance);
    EXPECT_GT(info.mNumLevels, 1);
    // The hierarchy spends its freedom near the bump, knot lines everywhere.
    EXPECT_LT(info.mNumDofs, knots.mSurface.NumPoles());
    // Flattened, it is one tensor-product B-spline over the finest knots.
    EXPECT_EQ(thb.mSurface.mNu, thb.mSurface.mNv);
    EXPECT_GE(thb.mSurface.NumPoles(), info.mNumDofs);
    // And the same bytes again.
    meshioplusplus::detail::ThbFitInfo again_info;
    const auto again = meshioplusplus::detail::thb_fit_patch(pts, tris, u, v, o, &again_info);
    EXPECT_EQ(again.mSurface.mPoles, thb.mSurface.mPoles);
    EXPECT_EQ(again_info.mNumDofs, info.mNumDofs);
}
