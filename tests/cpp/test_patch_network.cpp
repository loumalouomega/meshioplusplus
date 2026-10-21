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
// detail/bspline_curve and detail/patch_network (fit_cad's network mode,
// v17.12.0): exact knot insertion and sub-curves, shared edges that meet
// exactly, G1 across them, the extraordinary-vertex modes, determinism.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
#include <omp.h>
#endif

#include "analytic_surfaces.hpp"

#include "../../src/cpp/src/detail/bspline_curve.hpp"
#include "../../src/cpp/src/detail/patch_network.hpp"
#include "../../src/cpp/src/detail/quad_layout.hpp"

using namespace meshioplusplus;
using namespace analytic;
using detail::PatchTri;
using detail::PatchVec;

namespace {

struct PnTestPatch {
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
};

PnTestPatch pn_patch(const Mesh& rMesh, bool TopOnly) {
    PnTestPatch p;
    const double* x = rMesh.Points().As<double>();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        p.mPts.push_back({x[i * 3], x[i * 3 + 1], x[i * 3 + 2]});
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t t = 0; t < rMesh.Cells(0).NumCells(); ++t) {
        const PatchTri tri{static_cast<std::size_t>(c[t * 3]),
                           static_cast<std::size_t>(c[t * 3 + 1]),
                           static_cast<std::size_t>(c[t * 3 + 2])};
        if (TopOnly) {
            const PatchVec& a = p.mPts[tri[0]];
            const PatchVec& b = p.mPts[tri[1]];
            const PatchVec& d = p.mPts[tri[2]];
            const double nx = (b[1] - a[1]) * (d[2] - a[2]) - (b[2] - a[2]) * (d[1] - a[1]);
            const double ny = (b[2] - a[2]) * (d[0] - a[0]) - (b[0] - a[0]) * (d[2] - a[2]);
            const double nz = (b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]);
            if (!(nz > 0.3 * std::sqrt(nx * nx + ny * ny + nz * nz)))
                continue;
        }
        p.mTris.push_back(tri);
    }
    return p;
}

double pn_diag(const PnTestPatch& rP) {
    PatchVec lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const PatchTri& t : rP.mTris)
        for (std::size_t v : t)
            for (std::size_t j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], rP.mPts[v][j]);
                hi[j] = std::max(hi[j], rP.mPts[v][j]);
            }
    return std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                     (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

double pn_dist(const PatchVec& rA, const PatchVec& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

PatchVec pn_eval(const detail::BsplineSurface& rS, double U, double V) {
    detail::BsplineEval e;
    detail::bspline_eval(rS, U, V, 0, e);
    return {e.mS[0], e.mS[1], e.mS[2]};
}

void pn_side_uv(std::size_t Side, double S, double& rU, double& rV) {
    const double u[4] = {S, 1.0, 1.0 - S, 0.0}, v[4] = {0.0, S, 1.0, 1.0 - S};
    rU = u[Side];
    rV = v[Side];
}

// The largest gap between the two patches along every shared stretch.
double pn_g0_gap(const detail::QuadLayout& rL, const detail::PatchNetwork& rN) {
    double worst = 0.0;
    for (std::size_t c = 0; c < rL.mChains.size(); ++c) {
        struct Use {
            std::size_t mPatch, mSide;
            double mTb, mTe, mS0, mS1;
        };
        std::vector<Use> uses;
        for (std::size_t i = 0; i < rL.mPatches.size(); ++i)
            for (std::size_t j = 0; j < 4; ++j) {
                const auto& segs = rL.mSides[i][j].mSegments;
                if (segs.size() != 1)
                    continue;
                if (segs[0].mChain != c)
                    continue;
                uses.push_back({i, j, rN.mChainParams[c][segs[0].mBegin],
                                rN.mChainParams[c][segs[0].mEnd], 0.0, 1.0});
            }
        for (std::size_t a = 0; a < uses.size(); ++a)
            for (std::size_t b = a + 1; b < uses.size(); ++b) {
                if (rN.mPatches[uses[a].mPatch].mSurface.NumPoles() == 0 ||
                    rN.mPatches[uses[b].mPatch].mSurface.NumPoles() == 0)
                    continue;
                const double lo = std::max(std::min(uses[a].mTb, uses[a].mTe),
                                           std::min(uses[b].mTb, uses[b].mTe));
                const double hi = std::min(std::max(uses[a].mTb, uses[a].mTe),
                                           std::max(uses[b].mTb, uses[b].mTe));
                for (int k = 0; k <= 10 && hi > lo; ++k) {
                    const double t = lo + (hi - lo) * k / 10.0;
                    double ua, va, ub, vb;
                    pn_side_uv(uses[a].mSide, (t - uses[a].mTb) / (uses[a].mTe - uses[a].mTb), ua,
                               va);
                    pn_side_uv(uses[b].mSide, (t - uses[b].mTb) / (uses[b].mTe - uses[b].mTb), ub,
                               vb);
                    worst = std::max(
                        worst, pn_dist(pn_eval(rN.mPatches[uses[a].mPatch].mSurface, ua, va),
                                       pn_eval(rN.mPatches[uses[b].mPatch].mSurface, ub, vb)));
                }
            }
    }
    return worst;
}

}  // namespace

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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    const int saved = omp_get_max_threads();
    for (int threads : {1, 2, 8}) {
        omp_set_num_threads(threads);
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
#ifdef MESHIOPLUSPLUS_PARALLEL_OPENMP
    omp_set_num_threads(saved);
#endif
}
