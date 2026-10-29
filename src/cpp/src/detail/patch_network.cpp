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
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>

// Project includes
#include "meshioplusplus/log.hpp"

// Project includes (private, not installed)
#include "nurbs_eval.hpp"
#include "param_disc.hpp"
#include "patch_network.hpp"
#include "sparse_solve.hpp"
#include "surface_validate.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kPnPi = 3.14159265358979323846;

PatchVec pn_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
double pn_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
PatchVec pn_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double pn_norm(const PatchVec& rA) {
    return std::sqrt(pn_dot(rA, rA));
}

// A piece of a patch side: chain parameters `[mTb, mTe]` (in the side's
// direction) on side parameters `[mS0, mS1]`.
struct PnSegment {
    std::size_t mChain = 0;
    double mTb = 0.0, mTe = 0.0;
    double mS0 = 0.0, mS1 = 1.0;
};

// A patch side: its curve on [0, 1] from corner j to j + 1, and its pieces.
struct PnSide {
    BsplineCurve mCurve;
    std::vector<PnSegment> mSegments;
};

// A patch while it is fitted.
struct PnPatch {
    std::vector<PatchTri> mTris;
    std::vector<std::size_t> mVerts;  ///< the region points it uses, ascending
    std::array<PnSide, 4> mSides;
    std::vector<double> mExtraU, mExtraV;  ///< knots added by refinement or G1xx
    BsplineSurface mSurface;
    std::vector<double> mU, mV;                           ///< per region point
    std::vector<double> mU0, mV0;                         ///< the square map's, to start again from
    std::array<bool, 4> mEv{false, false, false, false};  ///< per corner
    std::array<PatchVec, 4> mEvNormal{};
};

// The surface parameters of side `Side`'s parameter `S`.
void pn_side_uv(std::size_t Side, double S, double& rU, double& rV) {
    switch (Side) {
        case 0:
            rU = S;
            rV = 0.0;
            break;
        case 1:
            rU = 1.0;
            rV = S;
            break;
        case 2:
            rU = 1.0 - S;
            rV = 1.0;
            break;
        default:
            rU = 0.0;
            rV = 1.0 - S;
            break;
    }
}

// Rebuild a patch's surface from its sides and extra knots: the boundary
// rows are the side curves, the interior poles a bilinearly blended (Coons)
// guess when `Fresh`, else kept where the knots did not change.
bool pn_build_surface(PnPatch& rP, std::size_t MaxPerDirection, std::size_t MaxTotal) {
    BsplineCurve c0 = rP.mSides[0].mCurve;
    BsplineCurve c2 = bspline_curve_reverse(rP.mSides[2].mCurve);
    BsplineCurve c1 = rP.mSides[1].mCurve;
    BsplineCurve c3 = bspline_curve_reverse(rP.mSides[3].mCurve);
    std::vector<double> ku = bspline_knot_union(c0.mKnots, c2.mKnots);
    std::vector<double> kv = bspline_knot_union(c1.mKnots, c3.mKnots);
    if (!rP.mExtraU.empty()) {
        std::vector<double> e = {0.0, 0.0, 0.0, 0.0};
        e.insert(e.end(), rP.mExtraU.begin(), rP.mExtraU.end());
        e.insert(e.end(), {1.0, 1.0, 1.0, 1.0});
        std::sort(e.begin(), e.end());
        ku = bspline_knot_union(ku, e);
    }
    if (!rP.mExtraV.empty()) {
        std::vector<double> e = {0.0, 0.0, 0.0, 0.0};
        e.insert(e.end(), rP.mExtraV.begin(), rP.mExtraV.end());
        e.insert(e.end(), {1.0, 1.0, 1.0, 1.0});
        std::sort(e.begin(), e.end());
        kv = bspline_knot_union(kv, e);
    }
    bspline_curve_refine_to(c0, ku);
    bspline_curve_refine_to(c2, ku);
    bspline_curve_refine_to(c1, kv);
    bspline_curve_refine_to(c3, kv);
    const std::size_t nu = c0.mPoles.size(), nv = c1.mPoles.size();
    if (c2.mPoles.size() != nu || c3.mPoles.size() != nv || nu < 4 || nv < 4)
        return false;
    if (nu > MaxPerDirection || nv > MaxPerDirection || nu * nv > MaxTotal)
        return false;
    BsplineSurface s;
    s.mDegreeU = s.mDegreeV = 3;
    s.mKnotsU = c0.mKnots;
    s.mKnotsV = c1.mKnots;
    s.mNu = nu;
    s.mNv = nv;
    s.mPoles.assign(3 * nu * nv, 0.0);
    const auto put = [&](std::size_t I, std::size_t J, const PatchVec& rX) {
        for (std::size_t c = 0; c < 3; ++c)
            s.mPoles[3 * (I + nu * J) + c] = rX[c];
    };
    for (std::size_t i = 0; i < nu; ++i) {
        put(i, 0, c0.mPoles[i]);
        put(i, nv - 1, c2.mPoles[i]);
    }
    for (std::size_t j = 0; j < nv; ++j) {
        put(0, j, c3.mPoles[j]);
        put(nu - 1, j, c1.mPoles[j]);
    }
    // Corners: both curves end there; take their mean.
    const auto mean = [](const PatchVec& rA, const PatchVec& rB) {
        return PatchVec{0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1]), 0.5 * (rA[2] + rB[2])};
    };
    put(0, 0, mean(c0.mPoles.front(), c3.mPoles.front()));
    put(nu - 1, 0, mean(c0.mPoles.back(), c1.mPoles.front()));
    put(nu - 1, nv - 1, mean(c2.mPoles.back(), c1.mPoles.back()));
    put(0, nv - 1, mean(c2.mPoles.front(), c3.mPoles.back()));
    // The interior: the Coons patch of the boundary rows.
    const auto at = [&](std::size_t I, std::size_t J, std::size_t C) {
        return s.mPoles[3 * (I + nu * J) + C];
    };
    for (std::size_t j = 1; j + 1 < nv; ++j)
        for (std::size_t i = 1; i + 1 < nu; ++i) {
            const double a = static_cast<double>(i) / static_cast<double>(nu - 1);
            const double b = static_cast<double>(j) / static_cast<double>(nv - 1);
            for (std::size_t c = 0; c < 3; ++c) {
                const double lu = (1.0 - b) * at(i, 0, c) + b * at(i, nv - 1, c);
                const double lv = (1.0 - a) * at(0, j, c) + a * at(nu - 1, j, c);
                const double bl = (1.0 - a) * (1.0 - b) * at(0, 0, c) +
                                  a * (1.0 - b) * at(nu - 1, 0, c) + a * b * at(nu - 1, nv - 1, c) +
                                  (1.0 - a) * b * at(0, nv - 1, c);
                s.mPoles[3 * (i + nu * j) + c] = lu + lv - bl;
            }
        }
    rP.mSurface = std::move(s);
    return true;
}

// The unit normal of a surface at (u, v), and the derivative across a side.
PatchVec pn_normal(const BsplineSurface& rS, double U, double V) {
    BsplineEval e;
    bspline_eval(rS, U, V, 1, e);
    const PatchVec su{e.mSu[0], e.mSu[1], e.mSu[2]}, sv{e.mSv[0], e.mSv[1], e.mSv[2]};
    PatchVec n = pn_cross(su, sv);
    const double l = pn_norm(n);
    return l > 0.0 ? PatchVec{n[0] / l, n[1] / l, n[2] / l} : PatchVec{0.0, 0.0, 0.0};
}

// A G1 sample: where on this patch, which derivative, the neighbour's normal.
struct PnG1Sample {
    double mU = 0.0, mV = 0.0;
    bool mAcrossU = false;  ///< the cross derivative is S_u (sides 1, 3), else S_v
    PatchVec mNormal{};
    bool mNearEv = false;
};

// Normal equations of a patch's interior poles: the rows the three
// coordinates share (data, smoothing) in one dense scalar block, the rows
// that couple them (G1, tangent planes) as triplets.
struct PnNormal {
    std::size_t mN = 0;  ///< scalar unknowns (interior poles)
    std::vector<double> mS;
    std::vector<double> mB;  ///< 3 per pole
    std::vector<SparseTriplet> mCoupled;
    explicit PnNormal(std::size_t N) : mN(N), mS(N * N, 0.0), mB(3 * N, 0.0) {}
    void Scalar(const std::vector<std::size_t>& rCols, const std::vector<double>& rVals,
                const double* pRhs, double W) {
        for (std::size_t a = 0; a < rCols.size(); ++a) {
            for (std::size_t b = 0; b < rCols.size(); ++b)
                mS[rCols[a] * mN + rCols[b]] += W * rVals[a] * rVals[b];
            for (std::size_t c = 0; c < 3; ++c)
                mB[3 * rCols[a] + c] += W * rVals[a] * pRhs[c];
        }
    }
    void Coupled(const std::vector<std::size_t>& rCols, const std::vector<double>& rVals,
                 double Rhs, double W) {
        for (std::size_t a = 0; a < rCols.size(); ++a) {
            for (std::size_t b = 0; b < rCols.size(); ++b)
                mCoupled.push_back({rCols[a], rCols[b], W * rVals[a] * rVals[b]});
            mB[rCols[a]] += W * rVals[a] * Rhs;
        }
    }
    SparseCsr Matrix() const {
        std::vector<SparseTriplet> t;
        for (std::size_t a = 0; a < mN; ++a)
            for (std::size_t b = 0; b < mN; ++b) {
                const double v = mS[a * mN + b];
                if (v != 0.0)
                    for (std::size_t c = 0; c < 3; ++c)
                        t.push_back({3 * a + c, 3 * b + c, v});
            }
        t.insert(t.end(), mCoupled.begin(), mCoupled.end());
        return csr_from_triplets(3 * mN, 3 * mN, std::move(t));
    }
};

}  // namespace

PatchNetwork patch_network_fit(const std::vector<PatchVec>& rPts,
                               const std::vector<PatchTri>& rTris, const QuadLayout& rLayout,
                               const PatchNetworkOptions& rOptions) {
    PatchNetwork out;
    const double tol = rOptions.mTolerance;
    if (!rLayout.mOk || !(tol > 0.0))
        return out;
    const double target =
        rOptions.mFreeformTolerance > 0.0 ? rOptions.mFreeformTolerance : 0.25 * tol;
    const std::size_t np = rLayout.mPatches.size();
    const std::size_t nc = rLayout.mChains.size();

    // 1. Edge curves, through the corners and T-junction points.
    std::vector<std::vector<std::size_t>> held(nc);
    for (const auto& sides : rLayout.mSides)
        for (const QuadSide& sd : sides)
            for (const QuadSegment& sg : sd.mSegments) {
                held[sg.mChain].push_back(sg.mBegin);
                held[sg.mChain].push_back(sg.mEnd);
            }
    // Free fits: each patch on its own square map, its boundary loose, to
    // put the smooth cuts on (a zigzag of mesh edges stands for a curve on
    // the smooth surface the patches beside it fit, not on the facets).
    struct PnFree {
        std::vector<PatchTri> mTris;
        std::vector<double> mU, mV;
        BsplineSurface mSurface;
        bool mOk = false;
    };
    std::vector<PnFree> free_fit(np);
    std::vector<std::vector<std::size_t>> patches_of(rPts.size());
    for (std::size_t i = 0; i < np; ++i) {
        PnFree& f = free_fit[i];
        const QuadPatch& q = rLayout.mPatches[i];
        for (std::size_t t : q.mTris)
            f.mTris.push_back(rTris[t]);
        std::vector<std::size_t> verts;
        for (const PatchTri& t : f.mTris)
            for (std::size_t v : t)
                verts.push_back(v);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        for (std::size_t v : verts)
            patches_of[v].push_back(i);
        std::array<std::size_t, 4> corners = q.mCornerPos;
        const DiscParam par = param_square(rPts, f.mTris, q.mBoundary, corners);
        if (!par.mOk)
            continue;
        f.mU = par.mU;
        f.mV = par.mV;
        const std::size_t n = std::max<std::size_t>(
            4, std::min<std::size_t>(
                   12,
                   static_cast<std::size_t>(std::sqrt(static_cast<double>(verts.size()))) / 2 + 2));
        f.mSurface = bspline_shape(3, 3, n, n, 0.0, 1.0, 0.0, 1.0);
        std::vector<double> su, sv, xyz;
        for (int pass = 0; pass < 3; ++pass) {
            su.clear();
            sv.clear();
            xyz.clear();
            for (std::size_t v : verts) {
                su.push_back(f.mU[v]);
                sv.push_back(f.mV[v]);
                xyz.insert(xyz.end(), {rPts[v][0], rPts[v][1], rPts[v][2]});
            }
            if (!bspline_lsq_fit(f.mSurface, su, sv, xyz, nullptr, 1e-6))
                break;
            f.mOk = true;
            for (std::size_t v : verts)
                bspline_project(f.mSurface, rPts[v].data(), f.mU[v], f.mV[v]);
        }
    }

    out.mChains.resize(nc);
    out.mChainParams.resize(nc);
    for (std::size_t c = 0; c < nc; ++c) {
        const QuadChain& ch = rLayout.mChains[c];
        std::vector<PatchVec> poly;
        for (std::size_t v : ch.mVerts)
            poly.push_back(rPts[v]);
        std::sort(held[c].begin(), held[c].end());
        held[c].erase(std::unique(held[c].begin(), held[c].end()), held[c].end());
        CurveFitOptions co;
        co.mTolerance = target;
        co.mMaxPoles = rOptions.mMaxPolesPerDirection;
        co.mInterpolate = held[c];
        CurveFitResult r;
        if (ch.mSmooth && poly.size() > 2) {
            // A cut of the layout's own zigzags along mesh edges: smoothed
            // (a binomial filter, its held points kept), each point moved
            // onto the free fit of a patch it belongs to, and the curve
            // fitted through those.
            std::vector<std::uint8_t> keep(poly.size(), 0);
            keep.front() = keep.back() = 1;
            for (std::size_t h : held[c])
                keep[h] = 1;
            std::vector<PatchVec> sm = poly;
            for (int pass = 0; pass < 3; ++pass) {
                std::vector<PatchVec> nx = sm;
                for (std::size_t k = 1; k + 1 < sm.size(); ++k)
                    if (!keep[k])
                        for (std::size_t j = 0; j < 3; ++j)
                            nx[k][j] = 0.25 * sm[k - 1][j] + 0.5 * sm[k][j] + 0.25 * sm[k + 1][j];
                sm.swap(nx);
            }
            for (std::size_t k = 0; k < sm.size(); ++k) {
                if (keep[k])
                    continue;
                const std::size_t v = ch.mVerts[k];
                double best = 1e300;
                PatchVec on = sm[k];
                for (std::size_t i : patches_of[v]) {
                    const PnFree& f = free_fit[i];
                    if (!f.mOk)
                        continue;
                    double u = f.mU[v], w = f.mV[v];
                    const double d = bspline_project(f.mSurface, sm[k].data(), u, w);
                    if (d < best) {
                        best = d;
                        BsplineEval e;
                        bspline_eval(f.mSurface, u, w, 0, e);
                        on = {e.mS[0], e.mS[1], e.mS[2]};
                    }
                }
                sm[k] = on;
            }
            co.mSmoothing = 1e-6;
            r = bspline_curve_fit(sm, co);
        } else {
            r = bspline_curve_fit(poly, co);
        }
        out.mChains[c] = r.mCurve;
        out.mChainParams[c] = r.mParams;
        if (!ch.mOnBoundary)
            ++out.mNumEdges;
    }
    // Extraordinary vertices: the tangent plane all curves leaving one share.
    const bool g1_ev = rOptions.mExtraordinary != FitCadExtraordinaryVertex::G0;
    std::map<std::size_t, PatchVec> ev_normal;
    for (std::size_t v : rLayout.mExtraordinary) {
        PatchVec n{0.0, 0.0, 0.0};
        for (const PatchTri& t : rTris)
            if (t[0] == v || t[1] == v || t[2] == v) {
                const PatchVec a =
                    pn_cross(pn_sub(rPts[t[1]], rPts[t[0]]), pn_sub(rPts[t[2]], rPts[t[0]]));
                for (std::size_t j = 0; j < 3; ++j)
                    n[j] += a[j];
            }
        const double l = pn_norm(n);
        if (l > 0.0)
            ev_normal[v] = {n[0] / l, n[1] / l, n[2] / l};
    }
    out.mNumExtraordinary = ev_normal.size();
    if (g1_ev)
        for (std::size_t c = 0; c < nc; ++c) {
            const std::vector<std::size_t>& cv = rLayout.mChains[c].mVerts;
            BsplineCurve& curve = out.mChains[c];
            if (curve.mPoles.size() < 4)
                continue;
            for (int end = 0; end < 2; ++end) {
                const auto it = ev_normal.find(end == 0 ? cv.front() : cv.back());
                if (it == ev_normal.end())
                    continue;
                const PatchVec& o = end == 0 ? curve.mPoles.front() : curve.mPoles.back();
                PatchVec& q = end == 0 ? curve.mPoles[1] : curve.mPoles[curve.mPoles.size() - 2];
                const double h = pn_dot(pn_sub(q, o), it->second);
                for (std::size_t j = 0; j < 3; ++j)
                    q[j] -= h * it->second[j];
            }
        }

    // 2. Patches: their sides as pieces of the chains, and their parameters.
    std::vector<PnPatch> patches(np);
    out.mPatches.resize(np);
    // Per chain, the (patch, side, segment) stretches over it, for G1.
    struct PnUse {
        std::size_t mPatch, mSide, mSeg;
    };
    std::vector<std::vector<PnUse>> chain_uses(nc);
    for (std::size_t i = 0; i < np; ++i) {
        PnPatch& p = patches[i];
        const QuadPatch& q = rLayout.mPatches[i];
        for (std::size_t f : q.mTris)
            p.mTris.push_back(rTris[f]);
        for (const PatchTri& t : p.mTris)
            for (std::size_t v : t)
                p.mVerts.push_back(v);
        std::sort(p.mVerts.begin(), p.mVerts.end());
        p.mVerts.erase(std::unique(p.mVerts.begin(), p.mVerts.end()), p.mVerts.end());
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t corner = q.mBoundary[q.mCornerPos[j]];
            const auto it = ev_normal.find(corner);
            if (it != ev_normal.end()) {
                p.mEv[j] = true;
                p.mEvNormal[j] = it->second;
            }
            PnSide& side = p.mSides[j];
            std::vector<BsplineCurve> parts;
            std::vector<double> lengths;
            for (const QuadSegment& sg : rLayout.mSides[i][j].mSegments) {
                PnSegment seg;
                seg.mChain = sg.mChain;
                seg.mTb = out.mChainParams[sg.mChain][sg.mBegin];
                seg.mTe = out.mChainParams[sg.mChain][sg.mEnd];
                parts.push_back(bspline_curve_segment(out.mChains[sg.mChain], seg.mTb, seg.mTe));
                double len = 0.0;
                const std::vector<std::size_t>& cv = rLayout.mChains[sg.mChain].mVerts;
                const std::size_t lo = std::min(sg.mBegin, sg.mEnd),
                                  hi = std::max(sg.mBegin, sg.mEnd);
                for (std::size_t t = lo; t < hi; ++t)
                    len += pn_norm(pn_sub(rPts[cv[t + 1]], rPts[cv[t]]));
                lengths.push_back(len);
                side.mSegments.push_back(seg);
                chain_uses[sg.mChain].push_back({i, j, side.mSegments.size() - 1});
            }
            double total = 0.0;
            for (double l : lengths)
                total += l;
            double at = 0.0;
            for (std::size_t k = 0; k < side.mSegments.size(); ++k) {
                side.mSegments[k].mS0 = at;
                at = k + 1 == side.mSegments.size() ? 1.0 : at + lengths[k] / total;
                side.mSegments[k].mS1 = at;
            }
            side.mCurve = bspline_curve_join(parts, lengths);
        }
        // G1xx: a double knot in the middle of the first span towards each
        // extraordinary corner.
        if (rOptions.mExtraordinary == FitCadExtraordinaryVertex::G1xx) {
            const auto first_span = [](const BsplineCurve& rC, bool AtStart) {
                const std::size_t p3 = static_cast<std::size_t>(rC.mDegree);
                return AtStart ? 0.5 * rC.mKnots[p3 + 1]
                               : 0.5 * (1.0 + rC.mKnots[rC.mKnots.size() - p3 - 2]);
            };
            if (p.mEv[0] || p.mEv[3])
                p.mExtraU.insert(p.mExtraU.end(), 2, first_span(p.mSides[0].mCurve, true));
            if (p.mEv[1] || p.mEv[2])
                p.mExtraU.insert(p.mExtraU.end(), 2, first_span(p.mSides[0].mCurve, false));
            if (p.mEv[0] || p.mEv[1])
                p.mExtraV.insert(p.mExtraV.end(), 2, first_span(p.mSides[1].mCurve, true));
            if (p.mEv[2] || p.mEv[3])
                p.mExtraV.insert(p.mExtraV.end(), 2, first_span(p.mSides[1].mCurve, false));
        }
        NetworkPatch& r = out.mPatches[i];
        if (!pn_build_surface(p, rOptions.mMaxPolesPerDirection, rOptions.mMaxTotalPoles)) {
            r.mReason = FitCadDemotionReason::BudgetExceeded;
            continue;
        }
        std::array<std::size_t, 4> corners = q.mCornerPos;
        const DiscParam par = param_square(rPts, p.mTris, q.mBoundary, corners);
        if (!par.mOk) {
            r.mReason = FitCadDemotionReason::ParamNotInjective;
            continue;
        }
        p.mU = par.mU;
        p.mV = par.mV;
        p.mU0 = par.mU;
        p.mV0 = par.mV;
    }

    // Each patch's neighbours: the patches sharing a stretch of a chain.
    std::vector<std::vector<std::size_t>> neighbours(np);
    for (std::size_t c = 0; c < nc; ++c)
        for (const PnUse& a : chain_uses[c])
            for (const PnUse& b : chain_uses[c])
                if (a.mPatch != b.mPatch)
                    neighbours[a.mPatch].push_back(b.mPatch);
    for (auto& n : neighbours) {
        std::sort(n.begin(), n.end());
        n.erase(std::unique(n.begin(), n.end()), n.end());
    }

    // Area-weighted vertex normals, for the fold test.
    std::vector<PatchVec> vertex_normal(rPts.size(), PatchVec{0.0, 0.0, 0.0});
    for (const PatchTri& t : rTris) {
        const PatchVec a = pn_cross(pn_sub(rPts[t[1]], rPts[t[0]]), pn_sub(rPts[t[2]], rPts[t[0]]));
        for (std::size_t v : t)
            for (std::size_t j = 0; j < 3; ++j)
                vertex_normal[v][j] += a[j];
    }

    // 3-6. Fit, sweep, refine.
    const double base_g1 = tol / std::max(1e-6, std::sin(rOptions.mG1ToleranceDeg * kPnPi / 180.0));
    double g1_weight = base_g1;
    std::vector<std::uint8_t> active(np, 0);
    for (std::size_t i = 0; i < np; ++i)
        active[i] = out.mPatches[i].mReason == FitCadDemotionReason::None ? 1 : 0;

    // The G1 samples of patch `I` against its neighbours' current surfaces.
    const auto g1_samples = [&](std::size_t I) {
        std::vector<PnG1Sample> out_s;
        const PnPatch& p = patches[I];
        for (std::size_t j = 0; j < 4; ++j)
            for (const PnSegment& seg : p.mSides[j].mSegments) {
                // A crease is no G1 edge: neither measured nor pulled flat.
                if (!rLayout.mChains[seg.mChain].mSmooth)
                    continue;
                for (const PnUse& use : chain_uses[seg.mChain]) {
                    if (use.mPatch == I || !active[use.mPatch])
                        continue;
                    const PnSegment& other =
                        patches[use.mPatch].mSides[use.mSide].mSegments[use.mSeg];
                    // The overlap of the two stretches, in chain parameters.
                    const double lo =
                        std::max(std::min(seg.mTb, seg.mTe), std::min(other.mTb, other.mTe));
                    const double hi =
                        std::min(std::max(seg.mTb, seg.mTe), std::max(other.mTb, other.mTe));
                    if (!(hi > lo) || !(std::abs(seg.mTe - seg.mTb) > 0.0) ||
                        !(std::abs(other.mTe - other.mTb) > 0.0))
                        continue;
                    const std::size_t ns = 12;
                    for (std::size_t k = 1; k < ns; ++k) {
                        const double t = lo + (hi - lo) * static_cast<double>(k) / ns;
                        const double s =
                            seg.mS0 + (seg.mS1 - seg.mS0) * (t - seg.mTb) / (seg.mTe - seg.mTb);
                        const double so = other.mS0 + (other.mS1 - other.mS0) * (t - other.mTb) /
                                                          (other.mTe - other.mTb);
                        PnG1Sample g;
                        pn_side_uv(j, s, g.mU, g.mV);
                        double uo = 0.0, vo = 0.0;
                        pn_side_uv(use.mSide, so, uo, vo);
                        g.mAcrossU = j == 1 || j == 3;
                        g.mNormal = pn_normal(patches[use.mPatch].mSurface, uo, vo);
                        // Near an extraordinary corner: its first span.
                        const BsplineCurve& sc = p.mSides[j].mCurve;
                        const double first = sc.mKnots[static_cast<std::size_t>(sc.mDegree) + 1];
                        const double last =
                            sc.mKnots[sc.mKnots.size() - static_cast<std::size_t>(sc.mDegree) - 2];
                        g.mNearEv = (p.mEv[j] && s < first) || (p.mEv[(j + 1) % 4] && s > last);
                        out_s.push_back(g);
                    }
                }
            }
        return out_s;
    };
    const auto g1_angle = [&](std::size_t I, bool SkipEv) {
        double worst = 0.0;
        for (const PnG1Sample& g : g1_samples(I)) {
            if (SkipEv && g.mNearEv)
                continue;
            const PatchVec n = pn_normal(patches[I].mSurface, g.mU, g.mV);
            // The angle between the tangent planes: two patches whose
            // parameters turn opposite ways have opposite normals on one
            // smooth surface (the face's sense is set apart).
            const double c = std::min(1.0, std::fabs(pn_dot(n, g.mNormal)));
            worst = std::max(worst, std::acos(c) * 180.0 / kPnPi);
        }
        return worst;
    };

    // Fit patch `I`'s interior poles.
    const auto fit = [&](std::size_t I, bool WithG1) {
        PnPatch& p = patches[I];
        BsplineSurface& s = p.mSurface;
        const std::size_t nu = s.mNu, nv = s.mNv;
        const std::size_t mi = nu - 2;
        // The scalar unknown of an interior pole, or -1 on the boundary.
        const auto scalar = [&](std::size_t Pole) -> std::ptrdiff_t {
            const std::size_t i = Pole % nu, j = Pole / nu;
            if (i == 0 || j == 0 || i + 1 == nu || j + 1 == nv)
                return -1;
            return static_cast<std::ptrdiff_t>((i - 1) + mi * (j - 1));
        };
        const auto unknown = [&](std::size_t Pole, std::size_t C) -> std::ptrdiff_t {
            const std::ptrdiff_t k = scalar(Pole);
            return k < 0 ? -1 : 3 * k + static_cast<std::ptrdiff_t>(C);
        };
        const std::size_t nint = mi * (nv - 2);
        const std::size_t dim = 3 * nint;
        PnNormal ne(nint);
        std::vector<std::size_t> cols;
        std::vector<double> vals;
        std::vector<std::size_t> rc;
        std::vector<double> rv;
        // Data: the vertices, the facet centres and the edges' middles; per
        // interior pole, how much data holds it.
        std::size_t nsamples = 0;
        std::vector<double> support(nint, 0.0);
        const auto data = [&](double U, double V, const PatchVec& rX) {
            // A foot on the domain's boundary is a point the fixed boundary
            // curve answers for (or one across a smooth edge, which belongs
            // to the neighbour): it would only bend the interior towards it.
            if (U <= 1e-9 || U >= 1.0 - 1e-9 || V <= 1e-9 || V >= 1.0 - 1e-9)
                return;
            bspline_row(s, U, V, cols, vals);
            rc.clear();
            rv.clear();
            double rhs[3] = {rX[0], rX[1], rX[2]};
            for (std::size_t k = 0; k < cols.size(); ++k) {
                const std::ptrdiff_t x = scalar(cols[k]);
                if (x < 0) {
                    for (std::size_t c = 0; c < 3; ++c)
                        rhs[c] -= vals[k] * s.mPoles[3 * cols[k] + c];
                } else {
                    rc.push_back(static_cast<std::size_t>(x));
                    rv.push_back(vals[k]);
                }
            }
            if (!rc.empty())
                ne.Scalar(rc, rv, rhs, 1.0);
            for (std::size_t k = 0; k < rc.size(); ++k)
                support[rc[k]] += rv[k];
            ++nsamples;
        };
        for (std::size_t v : p.mVerts)
            data(p.mU[v], p.mV[v], rPts[v]);
        for (const PatchTri& t : p.mTris) {
            const PatchVec c{(rPts[t[0]][0] + rPts[t[1]][0] + rPts[t[2]][0]) / 3.0,
                             (rPts[t[0]][1] + rPts[t[1]][1] + rPts[t[2]][1]) / 3.0,
                             (rPts[t[0]][2] + rPts[t[1]][2] + rPts[t[2]][2]) / 3.0};
            data((p.mU[t[0]] + p.mU[t[1]] + p.mU[t[2]]) / 3.0,
                 (p.mV[t[0]] + p.mV[t[1]] + p.mV[t[2]]) / 3.0, c);
            for (std::size_t k = 0; k < 3; ++k) {
                const std::size_t a2 = t[k], b2 = t[(k + 1) % 3];
                if (a2 > b2)
                    continue;  // each interior edge once (boundary edges once too)
                data(0.5 * (p.mU[a2] + p.mU[b2]), 0.5 * (p.mV[a2] + p.mV[b2]),
                     PatchVec{0.5 * (rPts[a2][0] + rPts[b2][0]), 0.5 * (rPts[a2][1] + rPts[b2][1]),
                              0.5 * (rPts[a2][2] + rPts[b2][2])});
            }
        }
        // Smoothing: second differences of the pole grid, much stronger
        // where little data holds the poles (a fine boundary knot vector
        // over a coarse patch), so the surface cannot fold between samples.
        const double lambda = rOptions.mSmoothing * static_cast<double>(nsamples) /
                              static_cast<double>(std::max<std::size_t>(1, nint));
        const auto starved = [&](std::size_t Pole) {
            const std::ptrdiff_t k = scalar(Pole);
            if (k < 0)
                return 0.0;
            const double sup = support[static_cast<std::size_t>(k)];
            return sup >= 1.0 ? 0.0 : 1.0 - sup;
        };
        const auto diff = [&](const std::vector<std::pair<std::size_t, double>>& rTerms, double W) {
            rc.clear();
            rv.clear();
            double rhs[3] = {0.0, 0.0, 0.0};
            for (const auto& t : rTerms) {
                const std::ptrdiff_t x = scalar(t.first);
                if (x < 0) {
                    for (std::size_t c = 0; c < 3; ++c)
                        rhs[c] -= t.second * s.mPoles[3 * t.first + c];
                } else {
                    rc.push_back(static_cast<std::size_t>(x));
                    rv.push_back(t.second);
                }
            }
            if (!rc.empty())
                ne.Scalar(rc, rv, rhs, W);
        };
        const auto near_ev = [&](std::size_t I2, std::size_t J2) {
            const std::size_t ci[4] = {0, nu - 1, nu - 1, 0}, cj[4] = {0, 0, nv - 1, nv - 1};
            for (std::size_t k = 0; k < 4; ++k)
                if (p.mEv[k] && (I2 > ci[k] ? I2 - ci[k] : ci[k] - I2) <= 2 &&
                    (J2 > cj[k] ? J2 - cj[k] : cj[k] - J2) <= 2)
                    return true;
            return false;
        };
        const bool g1x = rOptions.mExtraordinary == FitCadExtraordinaryVertex::G1x ||
                         rOptions.mExtraordinary == FitCadExtraordinaryVertex::G1xx;
        for (std::size_t j = 0; j < nv; ++j)
            for (std::size_t i = 0; i < nu; ++i) {
                const double w = lambda * (g1x && near_ev(i, j) ? 10.0 : 1.0) +
                                 1e-2 * starved(i + nu * j) * static_cast<double>(nsamples) /
                                     static_cast<double>(std::max<std::size_t>(1, nint));
                if (i >= 1 && i + 1 < nu)
                    diff({{(i - 1) + nu * j, 1.0}, {i + nu * j, -2.0}, {(i + 1) + nu * j, 1.0}}, w);
                if (j >= 1 && j + 1 < nv)
                    diff({{i + nu * (j - 1), 1.0}, {i + nu * j, -2.0}, {i + nu * (j + 1), 1.0}}, w);
                if (i + 1 < nu && j + 1 < nv)
                    diff({{i + nu * j, 1.0},
                          {(i + 1) + nu * j, -1.0},
                          {i + nu * (j + 1), -1.0},
                          {(i + 1) + nu * (j + 1), 1.0}},
                         2.0 * w);
            }
        // Extraordinary corners: the twist pole in the common tangent plane.
        if (g1_ev) {
            const std::size_t ci[4] = {0, nu - 1, nu - 1, 0}, cj[4] = {0, 0, nv - 1, nv - 1};
            const std::size_t ti[4] = {1, nu - 2, nu - 2, 1}, tj[4] = {1, 1, nv - 2, nv - 2};
            for (std::size_t k = 0; k < 4; ++k) {
                if (!p.mEv[k])
                    continue;
                const std::size_t twist = ti[k] + nu * tj[k], corner = ci[k] + nu * cj[k];
                rc.clear();
                rv.clear();
                double rhs = 0.0;
                for (std::size_t c = 0; c < 3; ++c) {
                    rc.push_back(static_cast<std::size_t>(unknown(twist, c)));
                    rv.push_back(p.mEvNormal[k][c]);
                    rhs += p.mEvNormal[k][c] * s.mPoles[3 * corner + c];
                }
                ne.Coupled(rc, rv, rhs, 1e3 * static_cast<double>(nsamples));
            }
        }
        // G1 across the shared sides.
        if (WithG1) {
            const std::vector<PnG1Sample> gs = g1_samples(I);
            const double w = g1_weight * g1_weight * static_cast<double>(nsamples) /
                             static_cast<double>(std::max<std::size_t>(1, gs.size()));
            const std::size_t p3 = 3;
            double bu[2 * (kBsplineMaxDegree + 1)], bv[2 * (kBsplineMaxDegree + 1)];
            for (const PnG1Sample& g : gs) {
                if (!g1_ev && g.mNearEv)
                    continue;
                const std::size_t su = nurbs_find_span(3, s.mKnotsU, nu, g.mU);
                const std::size_t sv = nurbs_find_span(3, s.mKnotsV, nv, g.mV);
                bspline_basis_ders(su, g.mU, 3, 1, s.mKnotsU, bu);
                bspline_basis_ders(sv, g.mV, 3, 1, s.mKnotsV, bv);
                // The current cross derivative's length, to make the row a sine.
                BsplineEval e;
                bspline_eval(s, g.mU, g.mV, 1, e);
                const double* d = g.mAcrossU ? e.mSu : e.mSv;
                const double dl = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                if (!(dl > 0.0))
                    continue;
                rc.clear();
                rv.clear();
                double rhs = 0.0;
                for (std::size_t b = 0; b <= p3; ++b)
                    for (std::size_t a = 0; a <= p3; ++a) {
                        const std::size_t pi = su - p3 + a, pj = sv - p3 + b;
                        const double coef =
                            g.mAcrossU ? bu[p3 + 1 + a] * bv[b] : bu[a] * bv[p3 + 1 + b];
                        if (coef == 0.0)
                            continue;
                        const std::size_t pole = pi + nu * pj;
                        for (std::size_t c = 0; c < 3; ++c) {
                            const double v = coef * g.mNormal[c] / dl;
                            const std::ptrdiff_t x = unknown(pole, c);
                            if (x < 0) {
                                rhs -= v * s.mPoles[3 * pole + c];
                            } else {
                                rc.push_back(static_cast<std::size_t>(x));
                                rv.push_back(v);
                            }
                        }
                    }
                if (!rc.empty())
                    ne.Coupled(rc, rv, rhs, w);
            }
        }
        const SparseCsr a = ne.Matrix();
        std::vector<double> x(dim, 0.0);
        for (std::size_t j = 1; j + 1 < nv; ++j)
            for (std::size_t i = 1; i + 1 < nu; ++i)
                for (std::size_t c = 0; c < 3; ++c)
                    x[3 * ((i - 1) + mi * (j - 1)) + c] = s.mPoles[3 * (i + nu * j) + c];
        const SparseSolveReport rep = csr_pcg(a, ne.mB.data(), x.data(), 1e-13, 20 * dim + 200);
        for (std::size_t j = 1; j + 1 < nv; ++j)
            for (std::size_t i = 1; i + 1 < nu; ++i)
                for (std::size_t c = 0; c < 3; ++c)
                    s.mPoles[3 * (i + nu * j) + c] = x[3 * ((i - 1) + mi * (j - 1)) + c];
        out.mPatches[I].mIterations += 1;
        // Parameter correction.
        for (std::size_t v : p.mVerts)
            bspline_project(s, rPts[v].data(), p.mU[v], p.mV[v]);
        return rep.mConverged;
    };

    // The distance of a point from the network near patch `I`: its own
    // surface, or, where the foot lands on that surface's boundary (the
    // smooth edge runs across the patch's zigzag one), a neighbour's.
    const auto network_distance = [&](std::size_t I, const double* pX, double& rU, double& rV) {
        double d = bspline_project(patches[I].mSurface, pX, rU, rV);
        const bool on_edge = rU <= 1e-9 || rU >= 1.0 - 1e-9 || rV <= 1e-9 || rV >= 1.0 - 1e-9;
        if (!on_edge)
            return d;
        for (std::size_t j : neighbours[I]) {
            if (!active[j])
                continue;
            // Start from the nearest of a coarse grid on the neighbour.
            const BsplineSurface& sj = patches[j].mSurface;
            double bu = 0.5, bv = 0.5, bd = 1e300;
            for (int a2 = 0; a2 <= 8; ++a2)
                for (int b2 = 0; b2 <= 8; ++b2) {
                    BsplineEval e;
                    bspline_eval(sj, a2 / 8.0, b2 / 8.0, 0, e);
                    const double dx = e.mS[0] - pX[0], dy = e.mS[1] - pX[1], dz = e.mS[2] - pX[2];
                    const double dd = dx * dx + dy * dy + dz * dz;
                    if (dd < bd) {
                        bd = dd;
                        bu = a2 / 8.0;
                        bv = b2 / 8.0;
                    }
                }
            d = std::min(d, bspline_project(sj, pX, bu, bv));
        }
        return d;
    };
    // A patch's test: every vertex within tol / 2 and every facet centre
    // within tol of the network, and, with `Validate`, the surface's samples
    // within tol of its own and its neighbours' triangles.
    const auto accept = [&](std::size_t I, bool Validate) {
        PnPatch& p = patches[I];
        NetworkPatch& r = out.mPatches[I];
        double worst = 0.0, sum = 0.0;
        for (std::size_t v : p.mVerts) {
            const double d = network_distance(I, rPts[v].data(), p.mU[v], p.mV[v]);
            worst = std::max(worst, d);
            sum += d * d;
        }
        double worst_c = 0.0;
        for (const PatchTri& t : p.mTris) {
            const double c[3] = {(rPts[t[0]][0] + rPts[t[1]][0] + rPts[t[2]][0]) / 3.0,
                                 (rPts[t[0]][1] + rPts[t[1]][1] + rPts[t[2]][1]) / 3.0,
                                 (rPts[t[0]][2] + rPts[t[1]][2] + rPts[t[2]][2]) / 3.0};
            double u = (p.mU[t[0]] + p.mU[t[1]] + p.mU[t[2]]) / 3.0;
            double w = (p.mV[t[0]] + p.mV[t[1]] + p.mV[t[2]]) / 3.0;
            worst_c = std::max(worst_c, network_distance(I, c, u, w));
        }
        r.mMaxDist = worst;
        r.mRmsDist = p.mVerts.empty() ? 0.0 : std::sqrt(sum / static_cast<double>(p.mVerts.size()));
        r.mMaxCentroidDist = worst_c;
        r.mAllowance = tol;
        if (!(worst <= 0.5 * tol) || !(worst_c <= tol))
            return false;
        // No fold: the surface's normal on a grid agrees with the skin's at
        // the nearest vertex in the parameters (a fold back along an edge
        // stays near the skin, so the distances above do not see it). The
        // grid is denser in the two spans next to each side, where the
        // corners' constraints bend the surface most.
        const auto fold_params = [](const std::vector<double>& rK) {
            std::vector<double> t;
            for (int a2 = 0; a2 <= 8; ++a2)
                t.push_back(a2 / 8.0);
            std::vector<double> k(rK.begin(), rK.end());
            k.erase(std::unique(k.begin(), k.end()), k.end());
            const std::size_t spans = k.size() > 1 ? k.size() - 1 : 0;
            for (std::size_t sp = 0; sp < spans; ++sp)
                if (sp < 2 || sp + 2 >= spans)
                    for (double f : {0.1, 0.35, 0.65, 0.9})
                        t.push_back(k[sp] + f * (k[sp + 1] - k[sp]));
            std::sort(t.begin(), t.end());
            return t;
        };
        const std::vector<double> fu = fold_params(p.mSurface.mKnotsU);
        const std::vector<double> fv = fold_params(p.mSurface.mKnotsV);
        for (double u : fu)
            for (double w : fv) {
                std::size_t near_v = p.mVerts.front();
                double bd = 1e300;
                for (std::size_t v : p.mVerts) {
                    const double d = (p.mU[v] - u) * (p.mU[v] - u) + (p.mV[v] - w) * (p.mV[v] - w);
                    if (d < bd) {
                        bd = d;
                        near_v = v;
                    }
                }
                if (pn_dot(pn_normal(p.mSurface, u, w), vertex_normal[near_v]) < 0.0) {
                    log::debug("pn: patch {} folds near ({}, {})", I, u, w);
                    return false;
                }
            }
        if (Validate) {
            std::vector<PatchTri> near = p.mTris;
            for (std::size_t j : neighbours[I])
                near.insert(near.end(), patches[j].mTris.begin(), patches[j].mTris.end());
            const SurfaceValidation val =
                surface_validate(p.mSurface, rPts, p.mTris, p.mU, p.mV, 3, &near);
            r.mReverseMaxDist = val.mMaxForward;
            if (!(val.mMaxForward <= tol))
                return false;
        }
        return true;
    };

    // Refine patch `I` at its worst sample's spans.
    const auto refine = [&](std::size_t I) {
        PnPatch& p = patches[I];
        const BsplineSurface& s = p.mSurface;
        double worst = -1.0, wu = 0.5, wv = 0.5;
        for (std::size_t v : p.mVerts) {
            // A foot on the surface's boundary is the zigzag's, not a misfit.
            if (p.mU[v] <= 1e-9 || p.mU[v] >= 1.0 - 1e-9 || p.mV[v] <= 1e-9 ||
                p.mV[v] >= 1.0 - 1e-9)
                continue;
            BsplineEval e;
            bspline_eval(s, p.mU[v], p.mV[v], 0, e);
            const PatchVec q{e.mS[0], e.mS[1], e.mS[2]};
            const double d = pn_norm(pn_sub(q, rPts[v]));
            if (d > worst) {
                worst = d;
                wu = p.mU[v];
                wv = p.mV[v];
            }
        }
        const auto mid = [](const std::vector<double>& rK, std::size_t N, double T) {
            const std::size_t sp = nurbs_find_span(3, rK, N, T);
            return 0.5 * (rK[sp] + rK[sp + 1]);
        };
        PnPatch trial = p;
        trial.mExtraU.push_back(mid(s.mKnotsU, s.mNu, wu));
        trial.mExtraV.push_back(mid(s.mKnotsV, s.mNv, wv));
        if (!pn_build_surface(trial, rOptions.mMaxPolesPerDirection, rOptions.mMaxTotalPoles))
            return false;
        p.mExtraU = trial.mExtraU;
        p.mExtraV = trial.mExtraV;
        p.mSurface = std::move(trial.mSurface);
        ++out.mPatches[I].mRefinements;
        return true;
    };

    // A. Accuracy, patch by patch: fit from the square map, correct the
    // parameters, fit again; refine where it fails and start again.
    const auto fit_from_start = [&](std::size_t I) {
        patches[I].mU = patches[I].mU0;
        patches[I].mV = patches[I].mV0;
        fit(I, false);
        fit(I, false);
    };
    for (std::size_t i = 0; i < np; ++i)
        if (active[i])
            fit_from_start(i);
    std::vector<std::uint8_t> ok(np, 0);
    for (std::size_t i = 0; i < np; ++i) {
        if (!active[i])
            continue;
        ok[i] = accept(i, false) ? 1 : 0;
        while (!ok[i] && out.mPatches[i].mRefinements < rOptions.mMaxRefinements && refine(i)) {
            fit_from_start(i);
            ok[i] = accept(i, false) ? 1 : 0;
        }
    }
    out.mSweeps = 1;

    // B. G1 sweeps on top: a patch a sweep pushes out of the tolerance keeps
    // its accurate surface (accuracy before smoothness).
    double worst_g1 = 0.0;
    for (std::size_t i = 0; i < np; ++i)
        if (active[i])
            worst_g1 = std::max(worst_g1, g1_angle(i, !g1_ev));
    for (int sweep = 1; sweep < rOptions.mMaxSweeps && worst_g1 > rOptions.mG1ToleranceDeg;
         ++sweep) {
        out.mSweeps = sweep + 1;
        for (std::size_t i = 0; i < np; ++i) {
            if (!active[i])
                continue;
            const BsplineSurface keep = patches[i].mSurface;
            const std::vector<double> ku = patches[i].mU, kv = patches[i].mV;
            const std::int64_t kr = out.mPatches[i].mRefinements;
            fit(i, true);
            // A smooth fit that leaves the tolerance gets knots where it
            // misses before it is given up: next to an extraordinary vertex
            // G1 needs more freedom than accuracy alone asked for.
            while (ok[i] && !accept(i, false) &&
                   out.mPatches[i].mRefinements < rOptions.mMaxRefinements && refine(i))
                fit(i, true);
            if (ok[i] && !accept(i, false)) {
                patches[i].mSurface = keep;
                patches[i].mU = ku;
                patches[i].mV = kv;
                out.mPatches[i].mRefinements = kr;
            }
        }
        const double before = worst_g1;
        worst_g1 = 0.0;
        for (std::size_t i = 0; i < np; ++i)
            if (active[i])
                worst_g1 = std::max(worst_g1, g1_angle(i, !g1_ev));
        if (!(worst_g1 < before))
            g1_weight = std::min(g1_weight * 2.0, base_g1 * 8.0);
    }
    out.mMaxG1AngleDeg = worst_g1;

    // The final test, with surface sampling; the G1 measure then over the
    // accepted patches only (a rejected one is no face).
    out.mOk = true;
    out.mMaxG1AngleDeg = 0.0;
    for (std::size_t i = 0; i < np; ++i) {
        NetworkPatch& r = out.mPatches[i];
        if (!active[i]) {
            out.mOk = false;
            continue;
        }
        r.mOk = accept(i, rOptions.mValidate);
        if (!r.mOk) {
            r.mReason = out.mPatches[i].mRefinements >= rOptions.mMaxRefinements
                            ? FitCadDemotionReason::HausdorffExceeded
                            : FitCadDemotionReason::BudgetExceeded;
            out.mOk = false;
        }
        r.mSurface = patches[i].mSurface;
        r.mU = patches[i].mU;
        r.mV = patches[i].mV;
    }
    for (std::size_t i = 0; i < np; ++i)
        active[i] = out.mPatches[i].mOk ? 1 : 0;
    for (std::size_t i = 0; i < np; ++i) {
        if (!active[i])
            continue;
        NetworkPatch& r = out.mPatches[i];
        r.mMaxG1AngleDeg = g1_angle(i, !g1_ev);
        out.mMaxG1AngleDeg = std::max(out.mMaxG1AngleDeg, r.mMaxG1AngleDeg);
        out.mMaxG1AngleEvDeg = std::max(out.mMaxG1AngleEvDeg, g1_angle(i, false));
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
