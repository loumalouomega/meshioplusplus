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
 * @file detail/bspline_curve.hpp
 * @brief Polynomial B-spline curves in 3D for `fit_cad`'s patch network
 * (v17.12.0): least-squares fitting of a vertex polyline with exact
 * interpolation constraints, exact knot insertion (Boehm, The NURBS Book
 * A5.1), sub-curves, reversal, C0 joins and common knot vectors -- the shared
 * edge curves both patches beside an edge take their boundary from, so the
 * two meet exactly.
 *
 * Curves are clamped (`p + 1` copies of each end knot); a fitted curve lives
 * on `[0, 1]` by chord length. The fit minimises `sum |C(t_k) - x_k|^2 +
 * lambda sum |second differences of the poles|^2` (or, with normals, the
 * residual weighted along each point's normal and lightly across it)
 * subject to `C(t_k) = x_k` at the ends and at the points asked for, solved
 * as its KKT system by Gaussian elimination with partial pivoting; each parameter is then moved
 * to its point's foot on the curve (Newton) and the fit repeated, and the
 * poles doubled until every point is within the tolerance or the budget is
 * spent.
 *
 * A **core-private** header: nothing installed names it. Serial and
 * deterministic.
 */

// System includes
#include <cstddef>
#include <vector>

// Project includes (private, not installed)
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// A polynomial B-spline curve in 3D.
struct BsplineCurve {
    int mDegree = 3;
    std::vector<double> mKnots;  ///< full, clamped: poles + degree + 1
    std::vector<PatchVec> mPoles;
    double T0() const { return mKnots[static_cast<std::size_t>(mDegree)]; }
    double T1() const { return mKnots[mPoles.size()]; }
};

/// The point at `T` (clamped to the domain) and, with `pD`, the first derivative.
PatchVec bspline_curve_eval(const BsplineCurve& rC, double T, PatchVec* pD = nullptr);

/// Insert the knot `T` `Times` times (exact; the multiplicity may reach the
/// degree + 1, which splits the curve there).
void bspline_curve_insert_knot(BsplineCurve& rC, double T, int Times);

/// Map the curve's parameters affinely onto `[T0, T1]` (same geometry).
void bspline_curve_reparam(BsplineCurve& rC, double T0, double T1);

/// The same curve run backwards, on the same domain.
BsplineCurve bspline_curve_reverse(const BsplineCurve& rC);

/// The piece of the curve between parameters `A` and `B` (backwards when
/// `A > B`), exactly, on `[0, 1]`.
BsplineCurve bspline_curve_segment(const BsplineCurve& rC, double A, double B);

/// Curves joined end to end (each one's end on the next one's start; the
/// shared ends averaged), C0 at the joins, on `[0, 1]` with the joins at the
/// cumulative fractions of `rLengths`. All of one degree.
BsplineCurve bspline_curve_join(const std::vector<BsplineCurve>& rParts,
                                const std::vector<double>& rLengths);

/// The knots (with multiplicity) both curves need to share one knot vector:
/// the union of theirs, each multiplicity the larger. Both on one domain.
std::vector<double> bspline_knot_union(const std::vector<double>& rA,
                                       const std::vector<double>& rB);

/// Insert into the curve the knots of `rKnots` (a superset of its own, same
/// domain and degree) it lacks, so its knot vector becomes `rKnots`.
void bspline_curve_refine_to(BsplineCurve& rC, const std::vector<double>& rKnots);

/// What `bspline_curve_fit` aims for.
struct CurveFitOptions {
    int mDegree = 3;
    double mTolerance = 0.0;  ///< every point within this of the curve
    std::size_t mMaxPoles = 64;
    double mSmoothing = 1e-9;  ///< relative to the points' count
    /// Points (indices) held exactly besides the two ends, ascending.
    std::vector<std::size_t> mInterpolate;
    int mParameterCorrections = 3;
    /// The poles to start from (0: the fewest the constraints allow); with
    /// an unreachable tolerance, the count fitted.
    std::size_t mStartPoles = 0;
    /// Per point, a unit normal: when given, a point's residual counts along
    /// it in full and across it (in the tangent plane) times
    /// `mTangentWeight`, and the tolerance applies along it -- a curve on a
    /// surface through a zigzag of its vertices, free to cut across the
    /// zigzag but not to leave the surface.
    std::vector<PatchVec> mNormals;
    double mTangentWeight = 1.0;
};

/// The fitted curve, every point's parameter on it and the largest distance.
struct CurveFitResult {
    bool mOk = false;
    BsplineCurve mCurve;
    std::vector<double> mParams;
    double mMaxDist = 0.0;
};

/// Fit a curve through a polyline (see the file comment); at least two
/// distinct points.
CurveFitResult bspline_curve_fit(const std::vector<PatchVec>& rPts,
                                 const CurveFitOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
