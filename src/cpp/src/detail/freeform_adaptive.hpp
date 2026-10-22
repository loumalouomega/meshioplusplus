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
#pragma once

/**
 * @file detail/freeform_adaptive.hpp
 * @brief The adaptive freeform face of `fit_cad`'s freeform modes
 * (v17.12.0): one patch, its parametrisation and B-spline fit, the rungs of
 * the escalation ladder that stay within the patch, and the checks that
 * accept the face -- free of OpenCASCADE, so both backends use it.
 *
 *  1. **Rungs.** Each rung parametrises (`param_ladder`) and fits
 *     (`bspline_fit_patch`, a cubic least-squares B-spline with parameter
 *     correction, tangent distance steps and knot insertion within
 *     `mMaxTotalPoles`): first with the given parametrisation and budgets,
 *     then with twice the poles per direction and twice the iterations
 *     (refine knots), then, under `Auto`, with least-squares conformal maps
 *     and with mean-value coordinates explicitly (re-parametrise). The fit's
 *     target at a vertex is `freeform_tolerance` (`tol / 4` by default),
 *     with no noise allowance; only on a patch whose skin seems noisy
 *     (`patch_noise_sigma` above `tol / 20`) and whose strict fit fails with
 *     white residuals does the rung try again at `min(4 sigma, 2.5 tol)`, so
 *     a clean patch is never held looser than the tolerance.
 *  2. **Noise or shape.** The vertices' signed residuals from the fit are
 *     white when they are noise and correlate across the patch's edges when
 *     the surface misses a shape (an under-resolved bump reads as noise to
 *     any local estimate). Only white residuals (edge correlation below 0.5)
 *     on a skin that shows noise itself earn the noise allowance
 *     `a = max(tol, min(4.5 sigma, 3 tol))`, `sigma` the smaller of the
 *     residual's RMS and the skin's own estimate (a fit limited by the data's
 *     density can leave white residuals on a clean patch too) -- a
 *     scanned skin is no more accurate than its noise, as a coarse
 *     tessellation is no more accurate than its chordal sag -- and a
 *     structured misfit sends the next rungs to the strict target. A clean
 *     patch has `a = tol`.
 *  3. **Acceptance.** Every vertex within `tol / 2` on a clean patch and
 *     `max(tol / 2, min(4.5 sigma, 3 tol))` on a noisy one, every facet centre
 *     within `a` (the legacy freeform test on a clean patch), by projection
 *     from their fitted parameters; and, with
 *     `mValidate`, the surface's samples inside the trim within `a` of the
 *     patch (`surface_validate`).
 *
 * A patch no rung explains comes back with a `FitCadDemotionReason`; the
 * caller splits it or demotes it. Serial and deterministic.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_cad.hpp"

// Project includes (private, not installed)
#include "bspline_lsq.hpp"
#include "param_disc.hpp"
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// What `freeform_adaptive_fit` aims for and may spend.
struct FreeformAdaptiveOptions {
    double mTolerance = 0.0;          ///< fit_cad's tolerance, > 0
    double mFreeformTolerance = 0.0;  ///< the fit's vertex target; 0: tol / 4
    double mSmoothing = 1e-12;
    std::size_t mMaxPoles = 40;  ///< per direction
    int mMaxIterations = 30;
    std::size_t mMaxTotalPoles = 400;  ///< per face
    ParamLadder mParameterization = ParamLadder::Auto;
    bool mValidate = true;  ///< sample the surface against the patch
    /// A patch with fewer triangles is not fitted (`TooSmall`).
    std::size_t mMinTriangles = 8;
    /// The parameter box widened on every side by this share of its size, so
    /// the face runs on past its patch (the least-bending extension) and an
    /// edge on its boundary projects onto its own surface.
    double mMargin = 0.05;
    /// Refine by truncated hierarchical B-splines (`thb_fit_patch`,
    /// `local_refine="thb"`, experimental) instead of knot lines where that
    /// meets the target (the face is the flattened surface), else by knot
    /// lines as without it.
    bool mThb = false;
};

/// The face, and how it got there.
struct FreeformAdaptiveResult {
    bool mOk = false;
    BsplineSurface mSurface;
    /// Per patch vertex, its parameters on the surface.
    std::vector<double> mU, mV;
    double mSigma = 0.0;      ///< the residual's RMS when white (noise), else 0
    double mAllowance = 0.0;  ///< `max(tol, min(4.5 sigma, 3 tol))`
    FitCadParamMethod mParam = FitCadParamMethod::None;
    FitCadDemotionReason mReason = FitCadDemotionReason::None;
    std::int64_t mIterations = 0;
    std::int64_t mRefinements = 0;
    /// The rungs climbed past the first.
    std::int64_t mEscalations = 0;
    /// With `mThb`, the hierarchical degrees of freedom the surface was
    /// fitted with (its flattened poles are `mSurface.NumPoles()`); 0 for a
    /// knot-line fit.
    std::size_t mNumDofs = 0;
    double mMaxDist = 0.0;  ///< largest vertex distance from the surface
    double mRmsDist = 0.0;
    double mMaxCentroidDist = 0.0;
    double mReverseMaxDist = 0.0;  ///< largest distance of a surface sample from the patch
};

/// What `freeform_accept` measured of a fitted surface against its patch.
struct FreeformAcceptance {
    bool mOk = false;
    /// The vertices' residuals do not correlate across edges.
    bool mResidualWhite = false;
    /// They are white on a skin that shows noise (`SkinSigma > tol / 20`):
    /// the noise allowance applies. `SkinSigma` 0 grants none.
    bool mWhite = false;
    /// The surface was sampled against the patch (it passed the vertex and
    /// centroid tests).
    bool mValidated = false;
    double mSigma = 0.0;
    double mAllowance = 0.0;
    double mMaxDist = 0.0;
    double mRmsDist = 0.0;
    double mMaxCentroidDist = 0.0;
    double mReverseMaxDist = 0.0;
};

/**
 * @brief The acceptance test of a freeform face (steps 2 and 3 of the file
 * comment), shared by the adaptive faces and the patch network: every vertex
 * projected from its parameters (`rU`, `rV`, updated to the feet), the noise
 * allowance from the residuals' whiteness and `SkinSigma`, every facet centre
 * within the allowance and, with `Validate`, the surface's samples too.
 */
FreeformAcceptance freeform_accept(const BsplineSurface& rS, const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris, std::vector<double>& rU,
                                   std::vector<double>& rV, double Tol, double SkinSigma,
                                   bool Validate);

/// Fit one patch (see the file comment).
FreeformAdaptiveResult freeform_adaptive_fit(const std::vector<PatchVec>& rPts,
                                             const std::vector<PatchTri>& rTris,
                                             const FreeformAdaptiveOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
