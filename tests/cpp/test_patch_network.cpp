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
// detail/bspline_curve and detail/patch_network (fit_cad's network mode,
// v17.12.0): exact knot insertion and sub-curves, shared edges that meet
// exactly, G1 across them, the extraordinary-vertex modes, determinism.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "thread_sweep.hpp"

#include "analytic_surfaces.hpp"
#include "patch_network_fixtures.hpp"

#include "../../src/cpp/src/detail/bspline_curve.hpp"
#include "../../src/cpp/src/detail/patch_network.hpp"
#include "../../src/cpp/src/detail/quad_layout.hpp"

using namespace meshioplusplus;
using namespace analytic;
using detail::PatchTri;
using detail::PatchVec;
using namespace pnfix;

TEST(BsplineCurve, KnotInsertionAndSegmentsAreExact) {
    std::vector<PatchVec> poly;
    for (int k = 0; k <= 40; ++k) {
        const double t = k / 40.0;
        poly.push_back({t, t * t, 0.3 * t * t * t});
    }
    detail::CurveFitOptions o;
    o.mTolerance = 1e-6;
    o.mInterpolate = {13, 27};
    const detail::CurveFitResult r = detail::bspline_curve_fit(poly, o);
    ASSERT_TRUE(r.mOk) << r.mMaxDist;
    // The held points are exact.
    for (std::size_t k : {std::size_t{0}, std::size_t{13}, std::size_t{27}, std::size_t{40}})
        EXPECT_LT(pn_dist(detail::bspline_curve_eval(r.mCurve, r.mParams[k]), poly[k]), 1e-12);
    // Insertion keeps the geometry.
    detail::BsplineCurve c = r.mCurve;
    detail::bspline_curve_insert_knot(c, 0.37, 2);
    detail::bspline_curve_insert_knot(c, 0.81, 4);
    for (int k = 0; k <= 50; ++k) {
        const double t = k / 50.0;
        EXPECT_LT(
            pn_dist(detail::bspline_curve_eval(c, t), detail::bspline_curve_eval(r.mCurve, t)),
            1e-13);
    }
    // A segment, forwards and backwards, is the same piece.
    const detail::BsplineCurve s = detail::bspline_curve_segment(r.mCurve, 0.2, 0.7);
    const detail::BsplineCurve b = detail::bspline_curve_segment(r.mCurve, 0.7, 0.2);
    for (int k = 0; k <= 20; ++k) {
        const double f = k / 20.0;
        const PatchVec want = detail::bspline_curve_eval(r.mCurve, 0.2 + 0.5 * f);
        EXPECT_LT(pn_dist(detail::bspline_curve_eval(s, f), want), 1e-13);
        EXPECT_LT(pn_dist(detail::bspline_curve_eval(b, 1.0 - f), want), 1e-13);
    }
    // Two pieces joined are the whole.
    const detail::BsplineCurve l = detail::bspline_curve_segment(r.mCurve, 0.0, 0.4);
    const detail::BsplineCurve h = detail::bspline_curve_segment(r.mCurve, 0.4, 1.0);
    const detail::BsplineCurve j = detail::bspline_curve_join({l, h}, {0.4, 0.6});
    for (int k = 0; k <= 20; ++k) {
        const double t = k / 20.0;
        EXPECT_LT(
            pn_dist(detail::bspline_curve_eval(j, t), detail::bspline_curve_eval(r.mCurve, t)),
            1e-12);
    }
}

TEST(PatchNetwork, SaddleTopFitsWithSharedEdges) {
    const PnTestPatch top = pn_patch(monkey_saddle_box(24), true);
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = 0.4;
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, lo);
    ASSERT_TRUE(l.mOk);
    ASSERT_GE(l.mPatches.size(), 4U);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    const detail::PatchNetwork n = detail::patch_network_fit(top.mPts, top.mTris, l, o);
    EXPECT_TRUE(n.mOk);
    for (const auto& p : n.mPatches) {
        EXPECT_TRUE(p.mOk);
        EXPECT_LE(p.mMaxDist, o.mTolerance);
    }
    EXPECT_LT(pn_g0_gap(l, n), 1e-9);
    EXPECT_LE(n.mMaxG1AngleDeg, 2.0);
    EXPECT_GT(n.mNumEdges, 0U);
    std::printf("  saddle: %zu patches, %zu edges, G1 %.3f deg, %d sweeps\n", n.mPatches.size(),
                n.mNumEdges, n.mMaxG1AngleDeg, n.mSweeps);
}

TEST(PatchNetwork, ClosedSkinsAndEveryExtraordinaryMode) {
    struct Case {
        const char* mName;
        Mesh mMesh;
    };
    std::vector<Case> cases;
    cases.push_back({"cc_cage3", cc_cage(3, 4)});
    cases.push_back({"organic_torus", organic_torus(48, 24)});
    for (const Case& c : cases)
        for (FitCadExtraordinaryVertex ev :
             {FitCadExtraordinaryVertex::G0, FitCadExtraordinaryVertex::G1,
              FitCadExtraordinaryVertex::G1x, FitCadExtraordinaryVertex::G1xx}) {
            SCOPED_TRACE(c.mName);
            const PnTestPatch p = pn_patch(c.mMesh, false);
            detail::QuadLayoutOptions lo;
            lo.mTargetSize = 0.25 * pn_diag(p);
            const detail::QuadLayout l = detail::quad_layout(p.mPts, p.mTris, lo);
            ASSERT_TRUE(l.mOk);
            detail::PatchNetworkOptions o;
            o.mTolerance = 2e-3 * pn_diag(p);
            o.mExtraordinary = ev;
            const detail::PatchNetwork n = detail::patch_network_fit(p.mPts, p.mTris, l, o);
            std::size_t good = 0;
            double worst = 0.0;
            for (const auto& q : n.mPatches) {
                good += q.mOk ? 1 : 0;
                worst = std::max(worst, q.mMaxDist);
                if (!q.mOk)
                    EXPECT_NE(q.mReason, FitCadDemotionReason::None);
            }
            EXPECT_LT(pn_g0_gap(l, n), 1e-9);
            // The cube cage's limit surface: every patch within the tolerance,
            // the double knot of G1xx aside.
            if (std::string(c.mName) == "cc_cage3" && ev != FitCadExtraordinaryVertex::G1xx) {
                EXPECT_EQ(good, n.mPatches.size());
                EXPECT_EQ(n.mNumExtraordinary, 8U);
            }
            std::printf(
                "  %-14s ev %d: %zu/%zu patches ok, worst %.2e (tol %.2e), G1 %.3f deg, %zu EVs, "
                "%d sweeps\n",
                c.mName, static_cast<int>(ev), good, n.mPatches.size(), worst, o.mTolerance,
                n.mMaxG1AngleDeg, n.mNumExtraordinary, n.mSweeps);
        }
}

TEST(PatchNetwork, DeterministicAcrossThreads) {
    const PnTestPatch top = pn_patch(hood_box(16, true), true);
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = 0.5;
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, lo);
    ASSERT_TRUE(l.mOk);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    std::vector<double> ref;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 8}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1}) {
        (void)threads;
#endif
        const detail::PatchNetwork n = detail::patch_network_fit(top.mPts, top.mTris, l, o);
        std::vector<double> poles;
        for (const auto& p : n.mPatches)
            poles.insert(poles.end(), p.mSurface.mPoles.begin(), p.mSurface.mPoles.end());
        if (ref.empty())
            ref = poles;
        else
            EXPECT_EQ(poles, ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}

// The coupled solve (doc/fit_cad_freeform_plan.md step 7): the edge curves
// move with both patches they bound, so the saddle top reaches the G1
// tolerance the sweeps miss, every patch still within the tolerance and the
// patches still meeting exactly on their shared curves.
TEST(PatchNetwork, CoupledSolveReachesG1OnTheSaddle) {
    const PnTestPatch top = pn_patch(monkey_saddle_box(24), true);
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = 0.4;
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, lo);
    ASSERT_TRUE(l.mOk);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    const detail::PatchNetwork sw = detail::patch_network_fit(top.mPts, top.mTris, l, o);
    o.mCoupled = true;
    const detail::PatchNetwork cp = detail::patch_network_fit(top.mPts, top.mTris, l, o);
    EXPECT_TRUE(cp.mOk);
    for (const auto& p : cp.mPatches) {
        EXPECT_TRUE(p.mOk);
        EXPECT_LE(p.mMaxDist, o.mTolerance);
    }
    EXPECT_LT(pn_g0_gap(l, cp), 1e-9);
    EXPECT_LE(cp.mMaxG1AngleDeg, sw.mMaxG1AngleDeg);
    EXPECT_LE(cp.mMaxG1AngleAwayDeg, o.mG1ToleranceDeg);
    // No extraordinary vertex on the saddle top: nothing at them.
    EXPECT_EQ(cp.mMaxG1AngleAtEvDeg, 0.0);
    EXPECT_EQ(cp.mMaxG1AngleAwayDeg, cp.mMaxG1AngleDeg);
    std::printf("  saddle: sweeps %.3f deg, coupled %.3f deg (kept %d)\n", sw.mMaxG1AngleDeg,
                cp.mMaxG1AngleDeg, cp.mCoupledKept ? 1 : 0);
}

// On a closed cage with extraordinary vertices: never fewer patches within
// the tolerance and no larger angle than the sweeps, G0 exact, and the angle
// at the extraordinary vertices reported apart from the one away from them.
TEST(PatchNetwork, CoupledSolveOnTheCage) {
    const PnTestPatch p = pn_patch(cc_cage(3, 4), false);
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = 0.25 * pn_diag(p);
    const detail::QuadLayout l = detail::quad_layout(p.mPts, p.mTris, lo);
    ASSERT_TRUE(l.mOk);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(p);
    const detail::PatchNetwork sw = detail::patch_network_fit(p.mPts, p.mTris, l, o);
    o.mCoupled = true;
    const detail::PatchNetwork cp = detail::patch_network_fit(p.mPts, p.mTris, l, o);
    const auto good = [](const detail::PatchNetwork& rN) {
        std::size_t g = 0;
        for (const auto& q : rN.mPatches)
            g += q.mOk ? 1 : 0;
        return g;
    };
    EXPECT_GE(good(cp), good(sw));
    EXPECT_LE(cp.mMaxG1AngleDeg, sw.mMaxG1AngleDeg);
    EXPECT_LT(pn_g0_gap(l, cp), 1e-9);
    EXPECT_EQ(cp.mNumExtraordinary, 8U);
    EXPECT_LE(cp.mMaxG1AngleAwayDeg, cp.mMaxG1AngleDeg);
    EXPECT_LE(cp.mMaxG1AngleAtEvDeg, cp.mMaxG1AngleDeg);
    EXPECT_GT(cp.mMaxG1AngleAtEvDeg, 0.0);
    std::printf(
        "  cage: sweeps %.3f deg (away %.3f, at EVs %.3f), coupled %.3f (away %.3f, at "
        "EVs %.3f), kept %d\n",
        sw.mMaxG1AngleDeg, sw.mMaxG1AngleAwayDeg, sw.mMaxG1AngleAtEvDeg, cp.mMaxG1AngleDeg,
        cp.mMaxG1AngleAwayDeg, cp.mMaxG1AngleAtEvDeg, cp.mCoupledKept ? 1 : 0);
}

TEST(PatchNetwork, CoupledSolveIsDeterministic) {
    const PnTestPatch top = pn_patch(hood_box(16, true), true);
    detail::QuadLayoutOptions lo;
    lo.mTargetSize = 0.5;
    const detail::QuadLayout l = detail::quad_layout(top.mPts, top.mTris, lo);
    ASSERT_TRUE(l.mOk);
    detail::PatchNetworkOptions o;
    o.mTolerance = 2e-3 * pn_diag(top);
    o.mCoupled = true;
    std::vector<double> ref;
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    const int saved = mio_test::max_threads();
    for (int threads : {1, 2, 8}) {
        mio_test::set_threads(threads);
#else
    for (int threads : {1, 1}) {
        (void)threads;
#endif
        const detail::PatchNetwork n = detail::patch_network_fit(top.mPts, top.mTris, l, o);
        std::vector<double> poles;
        for (const auto& p : n.mPatches)
            poles.insert(poles.end(), p.mSurface.mPoles.begin(), p.mSurface.mPoles.end());
        for (const auto& c : n.mChains)
            for (const auto& q : c.mPoles)
                poles.insert(poles.end(), q.begin(), q.end());
        if (ref.empty())
            ref = poles;
        else
            EXPECT_EQ(poles, ref);
    }
#if MESHIOPLUSPLUS_TEST_THREAD_SWEEP
    mio_test::set_threads(saved);
#endif
}
