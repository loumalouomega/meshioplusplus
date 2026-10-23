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
 * @file detail/edge_sampler.hpp
 * @brief Evaluating a B-rep face's surface and trimming curves, and sampling
 * a trimmed edge once so every face that uses it shares the same 3-D points:
 * the 1-D half of `mesh_cad_surface`.
 *
 *  - **Evaluation.** `BrepSurfaceEval` and `BrepCurveEval` wrap a
 *    `NurbsPatch` and a `BrepCurve` (a trim's curve in the face's
 *    parameters) for `detail/iga_kernels.hpp`'s rational evaluator, with
 *    first derivatives, at any degree.
 *  - **A trace** is one use of an edge: the face's surface composed with the
 *    trim's curve, `X(t) = S(c(t))`, run along the loop (`mDirection`). Its
 *    arc length is integrated with five-point Gauss-Legendre on the knot
 *    intervals of the curve, bisected until halving changes an interval's
 *    length by less than `1e-12` of it, and inverted by safeguarded Newton.
 *  - **Sampling.** `edge_sample` splits a trace into `max(1, ceil(L / h))`
 *    pieces of equal arc length (at least `mMinSegments`), then halves any
 *    piece whose sagitta (the distance of the curve at a quarter, half and
 *    three quarters of the piece from its chord) exceeds the chordal cap: the
 *    "linear deflection" of OpenCASCADE's mesher. A trace shorter than the
 *    tolerance (a pole) keeps its two ends.
 *  - **Sharing.** `edge_project` carries a master trace's samples to another
 *    use of the same edge (the neighbouring face, or the same face across a
 *    seam): each point is projected onto the other trace by Gauss-Newton
 *    from its arc-length fraction, monotonically between its neighbours. The
 *    caller keeps the master's 3-D points and takes only the parameters (so
 *    the faces share the points exactly); the largest residual is reported.
 *
 * A **core-private** header: nothing installed names it. Pure functions of
 * their inputs; evaluators are read-only after construction, so traces of
 * different edges may be sampled concurrently.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/nurbs.hpp"

// Project includes (private, not installed)
#include "iga_kernels.hpp"

namespace meshioplusplus {
namespace detail {

using EdgeVec3 = std::array<double, 3>;
using EdgeVec2 = std::array<double, 2>;

/// A two-parameter `NurbsPatch` with first derivatives; parameters outside
/// its knot range are clamped to it.
class BrepSurfaceEval {
public:
    /// @throws std::invalid_argument when the patch is not a valid surface.
    BrepSurfaceEval(const NurbsPatch& rPatch, bool Rational);
    /// `pX` (3) and, when not null, `pSu` and `pSv` (3 each).
    void Eval(double U, double V, double* pX, double* pSu, double* pSv) const;
    /// The knot range: `{u0, u1, v0, v1}`.
    const std::array<double, 4>& Box() const { return mBox; }

private:
    const NurbsPatch* pPatch = nullptr;
    IgaEvalNet mNet;
    std::array<std::size_t, 2> mCounts{0, 0};
    std::array<double, 4> mBox{0, 0, 0, 0};
};

/// A `BrepCurve` (its first two coordinates are the face's parameters for a
/// trim) with its derivative, over its active range or else its knot range.
class BrepCurveEval {
public:
    /// @throws std::invalid_argument when the curve is not valid.
    explicit BrepCurveEval(const BrepCurve& rCurve);
    void Eval(double T, double* pX, double* pD) const;
    double T0() const { return mT0; }
    double T1() const { return mT1; }
    /// The distinct knots strictly inside `(T0, T1)`, ascending.
    const std::vector<double>& Breaks() const { return mBreaks; }

private:
    std::int32_t mDegree = 1;
    std::vector<double> mKnots;
    std::vector<double> mValues;
    std::vector<double> mWeights;
    std::size_t mCount = 0;
    double mT0 = 0.0;
    double mT1 = 0.0;
    std::vector<double> mBreaks;
    IgaEvalNet mNet;
};

/// One use of an edge: `X(t) = S(c(t))` run along the loop, with its arc
/// length table.
class TrimTrace {
public:
    TrimTrace(const BrepSurfaceEval* pSurface, BrepCurveEval Curve, bool Forward);
    /// The curve parameters at the start and the end of the trace.
    double Start() const { return mForward ? mCurve.T0() : mCurve.T1(); }
    double End() const { return mForward ? mCurve.T1() : mCurve.T0(); }
    EdgeVec2 Uv(double T) const;
    EdgeVec3 Point(double T) const;
    /// `dX / dt`.
    EdgeVec3 Velocity(double T) const;
    double Length() const { return mCum.empty() ? 0.0 : mCum.back(); }
    /// The parameter at arc length `L` from the start (clamped to the trace).
    double ParamAtLength(double L) const;
    const BrepSurfaceEval& Surface() const { return *pSurface; }

private:
    const BrepSurfaceEval* pSurface = nullptr;
    BrepCurveEval mCurve;
    bool mForward = true;
    std::vector<double> mT;    ///< table parameters, in trace order
    std::vector<double> mCum;  ///< arc length from the start at each

    double Speed(double T) const;
    double Gauss(double A, double B) const;
    void Refine(double A, double B, double Whole, double La, int Depth, std::vector<double>& rT,
                std::vector<double>& rL) const;
};

struct EdgeSampleOptions {
    double mTargetSize = 0.0;              ///< h > 0
    double mMaxChordalDeviation = 0.0;     ///< 0: off
    std::int64_t mMinSegments = 1;         ///< at least this many pieces
    double mTolerance = 0.0;               ///< a trace this short is degenerate
    std::int64_t mMaxSegments = 10000000;  ///< refused (std::length_error) above
};

struct EdgeSamples {
    std::vector<double> mParams;     ///< curve parameters, start to end, ends included
    std::vector<EdgeVec3> mPoints;   ///< `X` at each
    std::vector<double> mFractions;  ///< arc-length fraction of each, 0 to 1
    double mLength = 0.0;
    bool mDegenerate = false;  ///< shorter than the tolerance: its two ends only
};

EdgeSamples edge_sample(const TrimTrace& rTrace, const EdgeSampleOptions& rOptions);

/// Whether two traces of one edge run the same way (compared at both ends and
/// at a quarter and three quarters of their lengths).
bool edge_same_direction(const TrimTrace& rA, const TrimTrace& rB);

struct EdgeProjection {
    std::vector<double> mParams;  ///< on the other trace, in its own order
    double mMaxGap = 0.0;         ///< the largest distance from a sample to the other trace
};

/// Carry a master's samples to another trace of the same edge; `SameDirection`
/// as `edge_same_direction` says. The ends map to the other trace's ends.
EdgeProjection edge_project(const TrimTrace& rOther, const EdgeSamples& rMaster,
                            bool SameDirection);

}  // namespace detail
}  // namespace meshioplusplus
