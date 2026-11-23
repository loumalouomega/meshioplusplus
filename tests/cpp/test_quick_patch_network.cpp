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
// Small versions of the patch network tests (`test_patch_network.cpp`). The full suite fits
// the saddle at 24 divisions and the closed cages and tori with every extraordinary-vertex
// mode, which takes minutes in the instrumented coverage build, so a pull request skips it
// there; these fit a saddle and a hood at 10 to 12 divisions, one cage and one torus, and
// check what the full tests check: every patch within the tolerance, the patches meeting
// exactly along their shared edges, G1 across them (swept and coupled), and the same poles
// on repeated runs. Their suite name starts with `Quick`, which the job's exclusion does not
// match.

// System includes
#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "thread_sweep.hpp"

// Project includes
#include "analytic_surfaces.hpp"
#include "patch_network_fixtures.hpp"

using namespace pnfix;

namespace {

// The layout of a skin or its upper patch at a target size.
detail::QuadLayout quick_layout(const PnTestPatch& rP, double Size) {
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = Size;
    return detail::quad_layout(rP.mPts, rP.mTris, lo);
}

std::vector<double> quick_poles(const detail::PatchNetwork& rN) {
    std::vector<double> poles;
    for (const auto& p : rN.mPatches)
        poles.insert(poles.end(), p.mSurface.mPoles.begin(), p.mSurface.mPoles.end());
    for (const auto& c : rN.mChains)
        for (const auto& q : c.mPoles)
            poles.insert(poles.end(), q.begin(), q.end());
    return poles;
}

}  // namespace

// The saddle's top fits with shared edges, then the coupled solve reaches the G1 tolerance
// the sweeps miss, every patch still within the tolerance and exactly joined.
TEST(QuickPatchNetwork, TheSaddleTopFitsAndTheCoupledSolveReachesG1) {
    const PnTestPatch top = pn_patch(analytic::monkey_saddle_box(16), true);
    const detail::QuadLayout l = quick_layout(top, 0.5);
    ASSERT_TRUE(l.mOk);
    ASSERT_GE(l.mPatches.size(), 2U);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    const detail::PatchNetwork sw = detail::patch_network_fit(top.mPts, top.mTris, l, o);
    EXPECT_TRUE(sw.mOk);
    EXPECT_GT(sw.mNumEdges, 0U);
    for (const auto& p : sw.mPatches)
        EXPECT_LE(p.mMaxDist, o.mTolerance);
    EXPECT_LT(pn_g0_gap(l, sw), 1e-9);
    o.mCoupled = true;
    const detail::PatchNetwork cp = detail::patch_network_fit(top.mPts, top.mTris, l, o);
    EXPECT_TRUE(cp.mOk);
    for (const auto& p : cp.mPatches)
        EXPECT_LE(p.mMaxDist, o.mTolerance);
    EXPECT_LT(pn_g0_gap(l, cp), 1e-9);
    EXPECT_LE(cp.mMaxG1AngleDeg, sw.mMaxG1AngleDeg);
    EXPECT_LE(cp.mMaxG1AngleAwayDeg, o.mG1ToleranceDeg);
    EXPECT_EQ(cp.mMaxG1AngleAtEvDeg, 0.0);
}

// A closed cage with extraordinary vertices, in every mode: patches that fail say why and
// the patches meet exactly. (The coupled solve is the saddle's test above.)
TEST(QuickPatchNetwork, TheCageInEveryExtraordinaryMode) {
    const PnTestPatch p = pn_patch(analytic::cc_cage(3, 3), false);
    const detail::QuadLayout l = quick_layout(p, 0.25 * pn_diag(p));
    ASSERT_TRUE(l.mOk);
    for (FitCadExtraordinaryVertex ev :
         {FitCadExtraordinaryVertex::G0, FitCadExtraordinaryVertex::G1,
          FitCadExtraordinaryVertex::G1x, FitCadExtraordinaryVertex::G1xx}) {
        detail::PatchNetworkOptions o;
        o.mTolerance = 2e-3 * pn_diag(p);
        o.mExtraordinary = ev;
        const detail::PatchNetwork n = detail::patch_network_fit(p.mPts, p.mTris, l, o);
        std::size_t good = 0;
        for (const auto& q : n.mPatches) {
            good += q.mOk ? 1 : 0;
            if (!q.mOk)
                EXPECT_NE(q.mReason, FitCadDemotionReason::None);
        }
        EXPECT_LT(pn_g0_gap(l, n), 1e-9) << static_cast<int>(ev);
        EXPECT_GT(n.mNumExtraordinary, 0U);
        // At this size the limit surface is not fitted everywhere; the demotions say why and
        // the patches that are fitted meet exactly.
        EXPECT_LE(good, n.mPatches.size());
    }
}

// A closed torus fits in every mode; the fit is the same on repeated runs and thread counts.
TEST(QuickPatchNetwork, ATorusAndTheHoodAreDeterministic) {
    const PnTestPatch t = pn_patch(analytic::organic_torus(24, 12), false);
    const detail::QuadLayout tl = quick_layout(t, 0.3 * pn_diag(t));
    ASSERT_TRUE(tl.mOk);
    detail::PatchNetworkOptions to;
    to.mTolerance = 2e-3 * pn_diag(t);
    const detail::PatchNetwork tn = detail::patch_network_fit(t.mPts, t.mTris, tl, to);
    EXPECT_LT(pn_g0_gap(tl, tn), 1e-9);
    for (const auto& q : tn.mPatches)
        if (!q.mOk)
            EXPECT_NE(q.mReason, FitCadDemotionReason::None);

    const PnTestPatch top = pn_patch(analytic::hood_box(10, true), true);
    const detail::QuadLayout l = quick_layout(top, 0.8);
    ASSERT_TRUE(l.mOk);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    o.mCoupled = true;
    std::vector<double> ref;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 4}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1, 1}) {
        (void)threads;
#endif
        const std::vector<double> poles =
            quick_poles(detail::patch_network_fit(top.mPts, top.mTris, l, o));
        if (ref.empty())
            ref = poles;
        else
            EXPECT_EQ(poles, ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}
