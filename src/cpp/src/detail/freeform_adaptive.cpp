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
#include <cmath>

// Project includes
#include "bspline_fit.hpp"
#include "freeform_adaptive.hpp"
#include "surface_validate.hpp"
#include "thb_spline.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

FitCadParamMethod fa_method(int Method) {
    switch (Method) {
        case 1:
            return FitCadParamMethod::Lscm;
        case 2:
            return FitCadParamMethod::MeanValue;
        case 3:
            return FitCadParamMethod::Plane;
        case 4:
            return FitCadParamMethod::Cylinder;
        default:
            return FitCadParamMethod::None;
    }
}

}  // namespace

FreeformAcceptance freeform_accept(const BsplineSurface& rS, const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris, std::vector<double>& rU,
                                   std::vector<double>& rV, double Tol, double SkinSigma,
                                   bool Validate) {
    FreeformAcceptance out;
    std::vector<std::uint8_t> used(rPts.size(), 0);
    for (const PatchTri& t : rTris)
        for (std::size_t v : t)
            used[v] = 1;
    // The patch's edges, for the residual's correlation across them.
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            edges.emplace_back(std::min(a, b), std::max(a, b));
        }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    // Each vertex projected (from its fitted parameters): its distance and
    // its signed residual along the surface normal.
    double worst_vertex = 0.0, sum = 0.0;
    std::size_t nv = 0;
    std::vector<double> resid(rPts.size(), 0.0);
    for (std::size_t v = 0; v < rPts.size(); ++v) {
        if (!used[v])
            continue;
        const double d = bspline_project(rS, rPts[v].data(), rU[v], rV[v]);
        BsplineEval e;
        bspline_eval(rS, rU[v], rV[v], 1, e);
        const double n[3] = {e.mSu[1] * e.mSv[2] - e.mSu[2] * e.mSv[1],
                             e.mSu[2] * e.mSv[0] - e.mSu[0] * e.mSv[2],
                             e.mSu[0] * e.mSv[1] - e.mSu[1] * e.mSv[0]};
        const double ln = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        const double side = ln > 0.0
                                ? ((rPts[v][0] - e.mS[0]) * n[0] + (rPts[v][1] - e.mS[1]) * n[1] +
                                   (rPts[v][2] - e.mS[2]) * n[2]) /
                                      ln
                                : d;
        resid[v] = side >= 0.0 ? d : -d;
        worst_vertex = std::max(worst_vertex, d);
        sum += d * d;
        ++nv;
    }
    const double rms = nv > 0 ? std::sqrt(sum / static_cast<double>(nv)) : 0.0;
    // Noise is white: its residuals do not correlate across edges, where a
    // shape the surface misses leaves residuals that do. Only white
    // residuals earn the noise allowance.
    double cross = 0.0, self = 0.0;
    for (const auto& e : edges) {
        cross += resid[e.first] * resid[e.second];
        self += 0.5 * (resid[e.first] * resid[e.first] + resid[e.second] * resid[e.second]);
    }
    const double rho = self > 0.0 ? cross / self : 0.0;
    // And no more noise than the skin itself shows (its quadric estimate): a
    // fit limited by the data's density can leave white residuals on a clean
    // patch too.
    out.mResidualWhite = rho < 0.5 && rms > 0.0;
    out.mWhite = out.mResidualWhite && SkinSigma > 0.05 * Tol;
    const double sigma = std::min(rms, SkinSigma);
    const double noise = out.mWhite ? std::min(4.5 * sigma, 3.0 * Tol) : 0.0;
    out.mSigma = out.mWhite ? sigma : 0.0;
    out.mAllowance = std::max(Tol, noise);
    const double vertex_limit = std::max(0.5 * Tol, noise);
    double worst_centroid = 0.0;
    for (const PatchTri& t : rTris) {
        const double c[3] = {(rPts[t[0]][0] + rPts[t[1]][0] + rPts[t[2]][0]) / 3.0,
                             (rPts[t[0]][1] + rPts[t[1]][1] + rPts[t[2]][1]) / 3.0,
                             (rPts[t[0]][2] + rPts[t[1]][2] + rPts[t[2]][2]) / 3.0};
        double u = (rU[t[0]] + rU[t[1]] + rU[t[2]]) / 3.0;
        double w = (rV[t[0]] + rV[t[1]] + rV[t[2]]) / 3.0;
        worst_centroid = std::max(worst_centroid, bspline_project(rS, c, u, w));
    }
    out.mMaxDist = worst_vertex;
    out.mRmsDist = rms;
    out.mMaxCentroidDist = worst_centroid;
    if (!(worst_vertex <= vertex_limit) || !(worst_centroid <= out.mAllowance))
        return out;
    if (Validate) {
        const SurfaceValidation val = surface_validate(rS, rPts, rTris, rU, rV);
        out.mReverseMaxDist = val.mMaxForward;
        out.mValidated = true;
        if (!(val.mMaxForward <= out.mAllowance))
            return out;
    }
    out.mOk = true;
    return out;
}

FreeformAdaptiveResult freeform_adaptive_fit(const std::vector<PatchVec>& rPts,
                                             const std::vector<PatchTri>& rTris,
                                             const FreeformAdaptiveOptions& rOptions) {
    FreeformAdaptiveResult out;
    const double tol = rOptions.mTolerance;
    if (!(tol > 0.0) || rTris.size() < rOptions.mMinTriangles) {
        out.mReason = FitCadDemotionReason::TooSmall;
        return out;
    }
    const double skin_sigma = patch_noise_sigma(rPts, rTris);
    out.mSigma = skin_sigma;
    // The noise allowance, bounded: no face is ever held looser than 3 tol.
    const double base =
        rOptions.mFreeformTolerance > 0.0 ? rOptions.mFreeformTolerance : 0.25 * tol;
    // The first fit may stop at the noise the skin seems to carry; whether it
    // is noise is decided from the fit's residual (below), never assumed.
    const double loose = std::max(base, std::min(4.0 * skin_sigma, 2.5 * tol));

    struct Rung {
        ParamLadder mWhich;
        std::size_t mPoles;
        int mIterations;
    };
    const std::size_t poles2 = std::min<std::size_t>(2 * rOptions.mMaxPoles, 1024);
    std::vector<Rung> rungs = {
        {rOptions.mParameterization, rOptions.mMaxPoles, rOptions.mMaxIterations},
        {rOptions.mParameterization, poles2, 2 * rOptions.mMaxIterations}};
    if (rOptions.mParameterization == ParamLadder::Auto) {
        rungs.push_back({ParamLadder::Lscm, poles2, 2 * rOptions.mMaxIterations});
        rungs.push_back({ParamLadder::Harmonic, poles2, 2 * rOptions.mMaxIterations});
    }
    // Each rung fits to the tolerance first, with no noise allowance; only a
    // skin that shows noise, whose strict fit fails with white residuals,
    // tries the same rung again at the noise's level.
    const bool noisy_skin = skin_sigma > 0.05 * tol;
    int last_method = -1;
    for (std::size_t r = 0; r < rungs.size(); ++r) {
        if (r > 0)
            ++out.mEscalations;
        const DiscParam par = param_ladder(rPts, rTris, rungs[r].mWhich);
        if (!par.mOk) {
            out.mReason =
                par.mDisc ? FitCadDemotionReason::ParamNotInjective : FitCadDemotionReason::NonDisc;
            continue;
        }
        // A re-parametrisation rung that lands on the map already tried, with
        // the same budget, would repeat it.
        if (r >= 2 && par.mMethod == last_method)
            continue;
        last_method = par.mMethod;
        for (int pass = 0; pass < 2; ++pass) {
            const bool loose_pass = pass == 1;
            const double target = loose_pass ? loose : base;
            BsplineFitOptions o;
            o.mTolerance = target;
            o.mCentroidTolerance = std::max(tol, 2.0 * target);
            o.mSmoothing = rOptions.mSmoothing;
            o.mMaxPoles = rungs[r].mPoles;
            o.mMaxIterations = rungs[r].mIterations;
            o.mMaxTotalPoles = rOptions.mMaxTotalPoles;
            o.mMargin = rOptions.mMargin;
            // THB (experimental) where it meets the target, else knot lines,
            // so it never loses a face the knots would keep.
            ThbFitInfo thb;
            BsplineFitResult fit;
            if (rOptions.mThb)
                fit = thb_fit_patch(rPts, rTris, par.mU, par.mV, o, &thb);
            if (!fit.mOk) {
                thb.mNumDofs = 0;
                fit = bspline_fit_patch(rPts, rTris, par.mU, par.mV, o);
            }
            out.mParam = fa_method(par.mMethod);
            out.mIterations += fit.mIterations;
            out.mRefinements += fit.mRefinements;
            out.mMaxDist = fit.mMaxVertexError;
            out.mRmsDist = fit.mRmsError;
            out.mMaxCentroidDist = fit.mMaxCentroidError;
            if (fit.mSurface.NumPoles() == 0) {
                out.mReason = FitCadDemotionReason::BudgetExceeded;
                break;
            }
            std::vector<double> fu = fit.mU, fv = fit.mV;
            const FreeformAcceptance acc =
                freeform_accept(fit.mSurface, rPts, rTris, fu, fv, tol,
                                loose_pass ? skin_sigma : 0.0, rOptions.mValidate);
            out.mSigma = acc.mSigma;
            out.mAllowance = acc.mAllowance;
            out.mMaxDist = acc.mMaxDist;
            out.mRmsDist = acc.mRmsDist;
            out.mMaxCentroidDist = acc.mMaxCentroidDist;
            out.mReverseMaxDist = acc.mReverseMaxDist;
            if (acc.mOk) {
                out.mOk = true;
                out.mReason = FitCadDemotionReason::None;
                out.mSurface = fit.mSurface;
                out.mNumDofs = thb.mNumDofs;
                out.mU = std::move(fu);
                out.mV = std::move(fv);
                return out;
            }
            out.mReason = acc.mValidated || fit.mOk ? FitCadDemotionReason::HausdorffExceeded
                                                    : FitCadDemotionReason::BudgetExceeded;
            // The noise level is tried only for residuals that look like
            // noise on a skin that shows some: a structured misfit is shape.
            if (loose_pass || !noisy_skin || !acc.mResidualWhite)
                break;
            ++out.mEscalations;
        }
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
