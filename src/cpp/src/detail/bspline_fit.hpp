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
 * @file detail/bspline_fit.hpp
 * @brief Fit a tensor-product B-spline surface to a parametrised triangle
 * patch within a tolerance, iteratively: the general freeform faces of
 * `fit_cad`.
 *
 * The samples are the patch's vertices and its triangles' centroids, at the
 * parameters `param_disc.hpp` gives (a centroid at the mean of its corners').
 * The surface starts as a degree-`p` spline on uniform clamped knots over the
 * parameters' bounding box (widened by `mMargin` on every side when asked),
 * fitted by least squares with a thin-plate term (`bspline_lsq.hpp`) that is
 * near zero over the knot spans holding samples and strong over those that
 * hold none (the corners of a disc's parametrisation), so the surface runs on
 * past the patch as the least-bending extension. Optional ghost samples of a
 * coarse fit (`mGhostWeight`) hold such regions harder. Corrected parameters
 * stay inside the samples' box. Then, each iteration:
 *
 *  1. **parameter correction**: every sample is projected onto the surface by
 *     Newton's method (Hoschek's correction), and the poles are fitted again
 *     at the new parameters -- point distance minimisation;
 *  2. when that stalls (less than 10% better), **tangent distance
 *     minimisation** steps (Wang, Pottmann and Liu, ACM TOG 2006): the poles
 *     minimise the squared distances along each sample's foot-point normal
 *     plus `mu` times the point distances, `mu` the Levenberg-Marquardt damping
 *     (halved on a step that improves, multiplied by ten on one that does not);
 *  3. **knot insertion**: a knot at the middle of every span, in `u` and in
 *     `v`, whose largest sample distance is above the target.
 *
 * It stops once every vertex is within `mTolerance` and every centroid within
 * `mCentroidTolerance`, or when the pole budget (`mMaxPoles` per direction) or the
 * iteration budget is spent. A knot insertion always moves the fit on (a
 * finer spline may need several before it pays), but the result is the best
 * state seen, and the history records the best so far, so it never increases.
 *
 * A **core-private** header: nothing installed names it. Serial and
 * deterministic.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes (private, not installed)
#include "bspline_lsq.hpp"
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// What `bspline_fit_patch` aims for and may spend.
struct BsplineFitOptions {
    int mDegree = 3;
    /// The largest distance of a vertex from the surface it aims for; > 0.
    double mTolerance = 0.0;
    /// The largest distance of a triangle's centroid, which lies its chordal
    /// sag off a smooth surface; at most 0 takes twice `mTolerance`.
    double mCentroidTolerance = 0.0;
    /// The thin-plate weight over the knot spans that hold samples, relative to
    /// the sample count on the unit parameter square, so it depends on neither.
    double mSmoothing = 1e-12;
    /// The weight over spans no sample reaches (the margin, a disc's corners):
    /// strong, so the surface runs on past its patch as the least-bending
    /// extension instead of wherever unconstrained poles fall.
    double mEmptySmoothing = 1e-4;
    /// The weight of a ghost sample (see the file comment) against a vertex's;
    /// 0 (the default) adds none: they hold an extension smooth at some cost
    /// in accuracy, and the span-wise smoothing usually suffices.
    double mGhostWeight = 0.0;
    /// The most poles per direction.
    std::size_t mMaxPoles = 40;
    /// The most poles in all (`mNu * mNv`); 0 (the default) sets no such cap,
    /// as before the option existed (v17.12.0).
    std::size_t mMaxTotalPoles = 0;
    /// The poles per direction to start from (at least the degree plus one).
    std::size_t mStartPoles = 6;
    int mMaxIterations = 30;
    /// The widening of the parameter box on every side, as a share of its size
    /// (0: the surface ends at the samples' box; `fit_cad` extends faces with
    /// OpenCASCADE where an edge needs it).
    double mMargin = 0.0;
    /// The weight of a triangle's centroid against a vertex's in the fits: a
    /// flat facet's centroid lies its chordal sag off a smooth surface.
    double mCentroidWeight = 0.25;
    /// The schedule (unreleased, doc/fit_cad_freeform_plan.md step 3). 0 (the
    /// default) is the fit as before. 1, continuation: the thin-plate weight
    /// starts at `mContinuationStart` times the sample count and falls two
    /// decades at a time to `mSmoothing` whenever parameter correction stops
    /// gaining 10%, before any knot is inserted (the falls are not counted
    /// against `mMaxIterations`); a knot goes in at the
    /// residual-weighted median of a span's bad samples, not its middle; and
    /// once the targets are met, interior knots are removed, the least loaded
    /// first, while the refitted surface still meets them.
    int mSchedule = 0;
    /// The thin-plate weight a continuation starts from, relative to the
    /// sample count.
    double mContinuationStart = 1e-4;
    /// The objective (unreleased, step 4). 0 (the default) as before: point
    /// distance minimisation with parameter correction, a tangent-distance
    /// step when it stalls. 1 point distance only. 2 a tangent-distance step
    /// tried first every iteration. 3 a squared-distance step (Wang, Pottmann
    /// and Liu 2006, extended to surfaces by the principal curvatures at each
    /// foot point) tried first every iteration.
    int mObjective = 0;
};

/// The fitted surface and how it got there.
struct BsplineFitResult {
    bool mOk = false;  ///< the targets were met
    BsplineSurface mSurface;
    double mMaxVertexError = 0.0;
    double mMaxCentroidError = 0.0;
    double mRmsError = 0.0;
    int mIterations = 0;
    int mRefinements = 0;  ///< knot insertion rounds
    int mRemovals = 0;     ///< knots removed (`mSchedule` 1)
    /// The largest sample distance after each kept state, non-increasing.
    std::vector<double> mHistory;
    /// The vertices' foot-point parameters on the final surface.
    std::vector<double> mU, mV;
};

/**
 * @brief The parameter box of the vertices `rTris` uses, widened by `Margin`
 * of its extent on each side: `pLo[2]` its corner and `pExt[2]` its size,
 * which map the parameters onto the unit square. `rUsed` flags the used
 * vertices. False (the outputs then unspecified) when the box is flat.
 */
bool bspline_fit_param_box(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                           const std::vector<double>& rU, const std::vector<double>& rV,
                           double Margin, double* pLo, double* pExt,
                           std::vector<std::uint8_t>& rUsed);

/**
 * @brief Fit a B-spline surface to a patch (see the file comment).
 * @param rPts the patch's vertices; @param rTris its triangles.
 * @param rU, rV the vertices' parameters (a flip-free parametrisation).
 */
BsplineFitResult bspline_fit_patch(const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris,
                                   const std::vector<double>& rU, const std::vector<double>& rV,
                                   const BsplineFitOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
