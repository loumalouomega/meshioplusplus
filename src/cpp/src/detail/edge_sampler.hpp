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
 *  - **Sampling.** `edge_sample` splits a trace at its corners (the knots
 *    where the trim is only C0, a polyline's vertices), each piece into
 *    `max(1, ceil(L / h))` parts of equal arc length (at least
 *    `mMinSegments` in all), then halves any
 *    piece whose sagitta (the distance of the curve at a quarter, half and
 *    three quarters of the piece from its chord) exceeds the chordal cap: the
 *    "linear deflection" of OpenCASCADE's mesher. A trace shorter than the
 *    tolerance (a pole) keeps its two ends. `edge_refine` halves chosen
 *    pieces by arc length afterwards (a crossing to separate).
 *  - **Sharing.** `edge_project` carries a master trace's samples to another
 *    use of the same edge (the neighbouring face, or the same face across a
 *    seam): each point is projected onto the other trace by Gauss-Newton
 *    from its arc-length fraction, monotonically between its neighbours. The
 *    caller keeps the master's 3-D points and takes only the parameters (so
 *    the faces share the points exactly); the largest residual is reported.
 *  - **Closed edges started elsewhere.** A closed edge's other use may begin
 *    at another point of the loop (a cap's circle listed from another angle
 *    than the cylinder's rim). `edge_closest_length` finds where on it the
 *    master begins, and `edge_project` and `edge_same_direction` take that
 *    arc length as the use's new start, running past its end and round from
 *    its start: the use is rotated onto the master's corner.
 *
 * A **core-private** header: nothing installed names it. Pure functions of
 * their inputs; evaluators are read-only after construction, so traces of
 * different edges may be sampled concurrently.
 */

// System includes
#include <array>
#include <cstdint>
#include <functional>
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
    /// `pX`, the first derivatives `pSu`, `pSv` and the second derivatives
    /// `pSuu`, `pSuv`, `pSvv` (3 each, none null).
    void Eval2(double U, double V, double* pX, double* pSu, double* pSv, double* pSuu, double* pSuv,
               double* pSvv) const;
    /// The largest absolute principal curvature at `(U, V)`, from the first
    /// and second fundamental forms; 0 where the tangents are parallel.
    double MaxCurvature(double U, double V) const;
    /// The knot range: `{u0, u1, v0, v1}`.
    const std::array<double, 4>& Box() const { return mBox; }
    /// The largest absolute control point coordinate: the size round-off in
    /// a point or a derivative is relative to.
    double Scale() const { return mScale; }

private:
    const NurbsPatch* pPatch = nullptr;
    IgaEvalNet mNet;
    std::array<std::size_t, 2> mCounts{0, 0};
    std::array<double, 4> mBox{0, 0, 0, 0};
    double mScale = 0.0;
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
    /// Those of them where the curve is only C0 (multiplicity at least the
    /// degree): its corners, which a sampling must keep.
    const std::vector<double>& Corners() const { return mCorners; }

private:
    std::int32_t mDegree = 1;
    std::vector<double> mKnots;
    std::vector<double> mValues;
    std::vector<double> mWeights;
    std::size_t mCount = 0;
    double mT0 = 0.0;
    double mT1 = 0.0;
    std::vector<double> mBreaks;
    std::vector<double> mCorners;
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
    /// The arc length from the start to the parameter `T` (clamped to the
    /// trace): the inverse of `ParamAtLength`.
    double LengthAt(double T) const;
    /// The arc lengths of the curve's corners from the start, ascending.
    std::vector<double> CornerLengths() const;
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

    double mNegligible = 0.0;    ///< a piece no longer than this needs no accuracy
    double mAbsTolerance = 0.0;  ///< nor one whose halves agree within this
};

struct EdgeSampleOptions {
    double mTargetSize = 0.0;              ///< h > 0
    double mMaxChordalDeviation = 0.0;     ///< 0: off
    std::int64_t mMinSegments = 1;         ///< at least this many pieces
    double mTolerance = 0.0;               ///< a trace this short is degenerate
    std::int64_t mMaxSegments = 10000000;  ///< refused (std::length_error) above
    /// The local size h(x) at a point of the edge (roadmap 6.3.4); null:
    /// `mTargetSize` everywhere. With it a stretch between corners is cut at
    /// equal increments of the integral of ds / h(s), which h(x) <= the
    /// target keeps at least as fine.
    std::function<double(const EdgeVec3&)> mSize;
};

struct EdgeSamples {
    std::vector<double> mParams;     ///< curve parameters, start to end, ends included
    std::vector<EdgeVec3> mPoints;   ///< `X` at each
    std::vector<double> mFractions;  ///< arc-length fraction of each, 0 to 1
    double mLength = 0.0;
    bool mDegenerate = false;  ///< shorter than the tolerance: its two ends only
};

EdgeSamples edge_sample(const TrimTrace& rTrace, const EdgeSampleOptions& rOptions);

/// Split pieces of a sampling at their arc-length midpoints: `rPieces` are
/// piece indices (piece k runs from sample k to k + 1), ascending and
/// distinct. Splitting by arc length rather than by parameter keeps the
/// halves equal on the edge, whatever its parameterisation.
void edge_refine(const TrimTrace& rTrace, EdgeSamples& rSamples,
                 const std::vector<std::size_t>& rPieces);

/// The same, piece `rPieces[i]` split at arc-length fraction `rAt[i]` of
/// the piece (0.5 its midpoint): Delaunay refinement's concentric shells,
/// a power-of-two distance from a sharp corner at one end.
void edge_refine(const TrimTrace& rTrace, EdgeSamples& rSamples,
                 const std::vector<std::size_t>& rPieces, const std::vector<double>& rAt);

/// The points inside each piece at arc-length fractions j / Order of the
/// piece (j = 1 .. Order - 1), piece by piece: the nodes of a curved
/// element's edge along the trimmed edge (roadmap 6.6). None for a
/// degenerate edge.
std::vector<std::array<double, 3>> edge_piece_points(const TrimTrace& rTrace,
                                                     const EdgeSamples& rSamples, int Order);

/// Whether two traces of one edge run the same way (compared at both ends and
/// at a quarter and three quarters of their lengths). With `StartLengthB` at
/// least 0, `rB` is a closed trace taken to start that far along it: only
/// the quarter points are compared, measured from there.
bool edge_same_direction(const TrimTrace& rA, const TrimTrace& rB, double StartLengthB = -1.0);

/// The arc length along a trace of its point closest to `rP`: the best of 64
/// equal steps of arc length, refined by Newton on the squared distance.
double edge_closest_length(const TrimTrace& rTrace, const EdgeVec3& rP);

struct EdgeProjection {
    std::vector<double> mParams;  ///< on the other trace, in its own order
    double mMaxGap = 0.0;         ///< the largest distance from a sample to the other trace
};

/// Carry a master's samples to another trace of the same edge; `SameDirection`
/// as `edge_same_direction` says. The ends map to the other trace's ends.
/// With `StartLength` at least 0, `rOther` is a closed trace rotated to start
/// that far along it (`edge_closest_length` of the master's start): both ends
/// map there, and the samples run to the trace's end and on from its start.
EdgeProjection edge_project(const TrimTrace& rOther, const EdgeSamples& rMaster, bool SameDirection,
                            double StartLength = -1.0);

}  // namespace detail
}  // namespace meshioplusplus
