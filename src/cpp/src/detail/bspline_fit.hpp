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
    /// The largest sample distance after each kept state, non-increasing.
    std::vector<double> mHistory;
    /// The vertices' foot-point parameters on the final surface.
    std::vector<double> mU, mV;
};

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
