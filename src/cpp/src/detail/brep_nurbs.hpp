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
 * @file detail/brep_nurbs.hpp
 * @brief Trimmed NURBS faces without a CAD kernel (v17.12.0): the native STEP
 * writer's faces -- planes, spheres, cylinders, cones, tori and B-splines
 * bounded by loops of edges -- as `BrepFace`s on `NurbsPatch`es with their
 * trimming curves, so `fit_cad`'s native backend returns its solid as a
 * `BrepInfo` that `write_cad_json` writes as Kratos CAD JSON.
 *
 *  - **Surfaces.** A plane is a bilinear patch over its loops' box; a
 *    sphere, cylinder, cone or torus is the exact rational surface of
 *    revolution of its profile (a line or a circular arc) over the angles its
 *    loops cover, its arcs rational quadratics of at most 90 degrees (The
 *    NURBS Book, sections 7.3 and 8.5); a B-spline is itself.
 *  - **Trims.** Each loop's points are mapped to the surface's parameters
 *    (in closed form, a B-spline's by projection), the angles unwrapped along
 *    the loop. A trim that keeps one parameter fixed (a circle of latitude, a
 *    generator) is a straight line in the parameters and so exact; a circle
 *    on a plane is its exact rational arc; any other is a polyline through
 *    the mapped points, refined until the surface between them stays within
 *    a small share of the tolerance of the edge curve when one is given.
 *  - **Seams.** A loop that winds once round a periodic direction (a band's
 *    rim) is closed through the seam: two such loops of opposite windings are
 *    joined by a seam line and its copy one period over, a single one through
 *    a pole (a sphere's) or an apex (a cone's). Seams and poles are edges of
 *    the face alone.
 *  - **Orientation.** The loops come with the face on their left seen from
 *    outside; a face whose surface normal points inwards (`Sense` false) has
 *    its loops reversed and `mSwappedNormal` set, OpenCASCADE's convention.
 *
 * A face this cannot express (a loop through a pole, windings of several
 * turns or in both directions of a torus) is refused with the reason, never
 * written wrong. A **core-private** header: nothing installed names it.
 * Serial and deterministic.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"

// Project includes (private, not installed)
#include "bspline_curve.hpp"
#include "bspline_lsq.hpp"
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// A rational quadratic arc of the unit circle from angle `mA0` to `mA1`, in
/// equal segments of at most 90 degrees; its parameter runs over the angles,
/// the knot of segment `k` at `mA0 + k * (mA1 - mA0) / mSegments`.
struct UnitArc {
    double mA0 = 0.0;
    double mA1 = 0.0;
    std::size_t mSegments = 0;
    std::vector<double> mKnots;                 ///< full, clamped, degree 2
    std::vector<std::array<double, 3>> mPoles;  ///< (x, y, weight)
};

/// The arc from `A0` to `A1 > A0` (radians; more than a turn is allowed).
UnitArc unit_arc(double A0, double A1);

/// The arc's parameter at angle `A` (in `[mA0, mA1]`), exactly: on a segment
/// of half-angle `h`, `tan(psi / 2) = tan(h / 2) (2 s - 1)`.
double unit_arc_param(const UnitArc& rArc, double A);

/// A NURBS curve's point at `T`: `(x, y, z)` of its control points.
PatchVec brep_curve_point(const BrepCurve& rC, double T);

/// A NURBS patch's point at `(U, V)` (parametric dimension 2).
PatchVec nurbs_patch_point(const NurbsPatch& rP, double U, double V);

/// A straight segment from `rA` to `rB` on `[0, 1]`.
BrepCurve brep_curve_line(const PatchVec& rA, const PatchVec& rB);

/// A polyline through `rPts` (two or more), its parameter the point index.
BrepCurve brep_curve_polyline(const std::vector<PatchVec>& rPts);

/// The arc of radius `R` about `rCentre` in the plane of the unit, orthogonal
/// `rX` and `rY`, from angle `A0` to `A1` (either way round), its parameter
/// the arc's own (`[min, max]` of the two, run backwards when `A1 < A0`).
BrepCurve brep_curve_arc(const PatchVec& rCentre, const PatchVec& rX, const PatchVec& rY, double R,
                         double A0, double A1);

/// The same curve run backwards, on the same parameter range.
BrepCurve brep_curve_reverse(const BrepCurve& rC);

/// A polynomial B-spline curve as a `BrepCurve`.
BrepCurve brep_curve_bspline(const BsplineCurve& rC);

/// What a face's surface is (`fit_cad:type` numbering).
struct BrepNurbsSurface {
    int mKind = 1;  ///< 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus, 6 B-spline
    /// The frame: the origin (a sphere's or torus's centre, a cone's apex, a
    /// point of the plane or of the cylinder's axis) and orthonormal axes, `mZ`
    /// the axis of revolution or the plane's normal.
    PatchVec mOrigin{0.0, 0.0, 0.0};
    PatchVec mX{1.0, 0.0, 0.0};
    PatchVec mY{0.0, 1.0, 0.0};
    PatchVec mZ{0.0, 0.0, 1.0};
    double mRadius = 0.0;        ///< sphere, cylinder; a torus's major radius
    double mMinor = 0.0;         ///< a torus's minor radius
    double mTanHalfAngle = 0.0;  ///< a cone's: radius per height above the apex
    const BsplineSurface* pBspline = nullptr;
};

/// One trim of a loop, along the loop.
struct BrepNurbsTrim {
    /// Its points, the first and last its ends; used when `mAt` is empty.
    std::vector<PatchVec> mPoints;
    /// The exact edge curve along the loop over `[0, 1]`, when there is one.
    std::function<PatchVec(double)> mAt;
    /// On a plane, the edge is this circle's arc (`mAt` gives its points).
    bool mCircle = false;
    PatchVec mCentre{0.0, 0.0, 0.0};
    double mRadius = 0.0;
    /// The edge it lies on, and whether the loop runs along its curve.
    std::int64_t mEdge = -1;
    bool mAlong = true;
};

/// A face built by `brep_nurbs_face`.
struct BrepNurbsFace {
    NurbsPatch mPatch;
    /// Its loops (`mLoops`, `mTrimmed`, `mRational`, `mSwappedNormal`); the
    /// ids, the patch index and the control-point ids are the caller's.
    BrepFace mFace;
    /// Per trim, in loop order: the edge it lies on (`mEdge`, -1 for one of
    /// this face's own seams or poles, then `mSeam`), and whether it runs
    /// along that edge's curve.
    struct Ref {
        std::int64_t mEdge = -1;
        std::int64_t mSeam = -1;
        bool mAlong = true;
    };
    std::vector<Ref> mRefs;
    /// This face's own edges: its seams (a 3-D curve each) and poles (none).
    std::vector<BrepCurve> mSeams;
    std::vector<bool> mSeamHasCurve;
};

/**
 * @brief Build the trimmed NURBS face of a surface and its loops (see the
 * file comment).
 * @param rLoops the loops, each its trims in order, the face on their left
 *        seen against the outward normal.
 * @param Sense whether the surface's own normal is the outward one.
 * @param Tol the writer's tolerance: every loop point must lie within it of
 *        the surface, and a trim with an edge curve follows it to a share.
 * @return false, with `rWhy`, for a face it cannot express.
 */
bool brep_nurbs_face(const BrepNurbsSurface& rSurface,
                     const std::vector<std::vector<BrepNurbsTrim>>& rLoops, bool Sense, double Tol,
                     BrepNurbsFace& rOut, std::string& rWhy);

}  // namespace detail
}  // namespace meshioplusplus
