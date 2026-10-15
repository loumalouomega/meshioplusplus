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
 * @file feature_fit.hpp
 * @brief The OpenCASCADE-free half of `fit_cad`'s fillet and thread features
 * (roadmap 6.3.1): the radius of a constant-radius rolling-ball fillet between
 * two fitted primitives, and a straight-flank thread surface (a ruled helicoid)
 * fitted to a patch; and, for what those do not explain, a rolling-ball blend
 * of varying radius between any two surfaces and a helical surface of any
 * profile, with or without a known axis. Since v17.2.0 (roadmap 6.3) also the
 * threads whose pitch varies or that are cut on a taper, the blends whose
 * cross-section is a straight line (a chamfer) or a conic other than the
 * rolling ball's circle, and the boundary of a corner patch where three or
 * more blends meet.
 *
 * Not installed: a core-private helper of `operations/fit_cad.cpp`.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

// Project includes
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The distance from a point to a fitted primitive.
 * @param Type 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus.
 * @param pParams the ten `fit:patch_params` entries (see `fit_primitives.hpp`).
 * @return the unsigned distance, or +inf for an unknown type.
 */
double primitive_distance(int Type, const double* pParams, const PatchVec& rP);

/// A rolling-ball fillet's radius, fitted to a patch between two primitives.
struct FilletFit {
    bool mOk = false;
    /// The ball radius.
    double mRadius = 0.0;
    /// +1 when the ball centres lie inside the solid (a convex, rounded edge:
    /// centre = point - radius * outward normal), -1 outside (a concave fillet).
    int mSide = 0;
    /// RMS over the interior vertices of (distance from the ball centre to
    /// each primitive - radius).
    double mRms = 0.0;
    /// The ball centres of the interior vertices: points of the fillet's spine.
    std::vector<PatchVec> mCentres;
};

/**
 * @brief Fit the radius of a constant-radius rolling-ball fillet.
 *
 * Every interior vertex of the patch (one no boundary edge touches), moved by
 * `r` against or along its outward normal, should be a ball centre: `r` from
 * both primitives. The radius minimises the mean squared misfit of both
 * distances, over a log-spaced scan of `[RMin, RMax]` refined by golden
 * section, for each side; the better side wins.
 *
 * @param rPts the patch's points.
 * @param rTris its triangles, wound outwards (their normals point out of the solid).
 * @return the fit, `mOk` false when the patch has fewer than four interior
 *         vertices or the range is empty.
 */
FilletFit fillet_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                     int Type1, const double* pParams1, int Type2, const double* pParams2,
                     double RMin, double RMax);

/**
 * @brief A straight-flank thread surface: the ruled helicoid
 * `P(rho, theta) = O + rho (cos theta E1 + sin theta E2) + (C + K theta + T rho) A`,
 * `K = pitch / (2 pi)`, `T` the flank's slope (the tangent of its angle to the
 * radial plane). A thread of varying pitch adds `Q (theta - ThetaRef)^2` to
 * the axial term: its pitch changes linearly along it. A straight flank cut
 * on a taper of half-angle `alpha` is the same surface, `K` the apparent
 * advance `(1 - T tan alpha)` times the true one.
 */
struct HelicoidFit {
    bool mOk = false;
    PatchVec mOrigin{}, mAxis{}, mE1{}, mE2{};
    double mK = 0.0;  // axial advance per radian: pitch / (2 pi), signed by the hand
    double mT = 0.0;  // axial slope per unit radius
    double mC = 0.0;  // axial offset
    /// The pitch's rate of change: the axial term's quadratic coefficient
    /// about `mThetaRef` (0 for a constant pitch).
    double mQ = 0.0, mThetaRef = 0.0;
    /// The patch's extent in radius and in the unwrapped angle.
    double mRhoMin = 0.0, mRhoMax = 0.0, mThetaMin = 0.0, mThetaMax = 0.0;
    /// RMS and largest orthogonal distance of the vertices to the surface.
    double mRms = 0.0, mMax = 0.0;
    /// The surface point at `(Rho, Theta)`.
    PatchVec At(double Rho, double Theta) const;
    /// The pitch: the axial advance per turn (at `mThetaRef`).
    double Pitch() const;
    /// The local pitch at the unwrapped angle `Theta`.
    double PitchAt(double Theta) const;
};

/**
 * @brief Fit a ruled helicoid about a known axis to a patch.
 *
 * The angle of every vertex about the axis is unwrapped continuously over the
 * patch's edges (breadth-first from its lowest vertex), which makes the model
 * linear in `(C, K, T)`: one least-squares solve.
 *
 * @param rOrigin a point on the axis; @param rAxis its direction.
 * @param VariablePitch also fit the pitch's rate of change (`mQ`, one more
 *        column in the same linear solve); false keeps the constant-pitch fit.
 * @return the fit; `mOk` false when the patch touches the axis, cannot be
 *         unwrapped consistently, spans less than half a turn, or the solve is
 *         singular (a patch with no advance or no radial extent).
 */
HelicoidFit helicoid_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                         const PatchVec& rOrigin, const PatchVec& rAxis,
                         bool VariablePitch = false);

/// The point of a surface nearest to a point; false when there is none.
using ClosestPoint = std::function<bool(const PatchVec& rX, PatchVec& rFoot)>;

/**
 * @brief A smooth function of one parameter, fitted by least squares: a
 * Chebyshev series in `T` over `[mT0, mT1]`, or with `mPeriodic` a Fourier
 * series of period `mT1 - mT0`.
 */
struct FeatureLaw {
    bool mPeriodic = false;
    double mT0 = 0.0, mT1 = 1.0;
    /// Chebyshev: `T_0 .. T_d`; Fourier: `1, cos t, sin t, .., cos d t, sin d t`.
    std::vector<double> mCoef;
    double At(double T) const;
    /// The derivative in `T`.
    double Slope(double T) const;
};

/**
 * @brief Fit a law of degree @p Degree to `(rT[i], rY[i])` by least squares.
 * @return false when there are too few points or the system is singular.
 */
bool feature_law_fit(const std::vector<double>& rT, const std::vector<double>& rY, bool Periodic,
                     double T0, double T1, int Degree, FeatureLaw& rOut);

/**
 * @brief A rolling-ball blend between two surfaces given by their closest
 * points, of constant or varying radius (roadmap 6.3.1). Its parameter is the
 * arc length along the rail, the patch boundary where it touches the first
 * surface.
 */
struct RollingBallFit {
    bool mOk = false;
    /// +1 when the ball centres lie inside the solid, -1 outside.
    int mSide = 0;
    /// The radius changes along the rail by more than the tolerance.
    bool mVariable = false;
    /// The radius along the rail.
    FeatureLaw mRadius;
    /// The rail, in order; for a closed rail its first point is not repeated.
    std::vector<PatchVec> mRail;
    bool mClosed = false;
    /// Arc length at each rail point; `mLength` is the whole rail's.
    std::vector<double> mRailT;
    double mLength = 0.0;
    /// The vertices used: their rail parameter and ball centre.
    std::vector<double> mT;
    std::vector<PatchVec> mCentres;
    /// RMS and largest deviation of the vertices' own radii from the law.
    double mRms = 0.0, mMax = 0.0;
    /// The rail point and unit tangent at `T` (wrapped on a closed rail,
    /// extended along the end tangents past an open one).
    PatchVec RailAt(double T) const;
    PatchVec RailTangent(double T) const;
};

/**
 * @brief Fit a rolling-ball blend to a patch between two surfaces.
 *
 * Each interior vertex, moved by `r` against (`mSide` +1) or along its
 * outward normal, should be a ball centre `r` from both surfaces; `r` is
 * found per vertex (a log-spaced scan of `[RMin, RMax]` refined by golden
 * section), and the radii are fitted as a law along the rail: a constant when
 * that stays within the tolerance, else a Chebyshev (open rail) or Fourier
 * (closed rail) series of rising degree.
 *
 * @param rTris the patch, wound outwards.
 * @param Tol the tolerance the radii must be explained to.
 * @return the fit; `mOk` false when the patch has fewer than eight interior
 *         vertices, a tenth of them fit no ball, no rail is found, or no law
 *         explains the radii.
 */
RollingBallFit rolling_ball_fit(const std::vector<PatchVec>& rPts,
                                const std::vector<PatchTri>& rTris, const ClosestPoint& rS1,
                                const ClosestPoint& rS2, double RMin, double RMax, double Tol);

/// One cross-section of a rolling-ball blend: the ball and where it touches.
struct RollingBallSection {
    bool mOk = false;
    PatchVec mCentre{}, mFoot1{}, mFoot2{};
    double mRadius = 0.0;
};

/**
 * @brief The ball at rail parameter `T`: Gauss-Newton, in the plane through
 * the rail point normal to @p rTangent, from the nearest vertex centre moved
 * into it, for the centre `r(T)` from both surfaces.
 */
RollingBallSection rolling_ball_section(const RollingBallFit& rFit, const ClosestPoint& rS1,
                                        const ClosestPoint& rS2, double T, const PatchVec& rTangent,
                                        double Scale);

/**
 * @brief A helical surface of any profile: the profile curve, in the plane
 * through the axis, swept by the screw motion of advance `K` per radian,
 * `P(s, theta) = O + rho (cos theta E1 + sin theta E2) + (zeta + K theta) A`,
 * where the profile gives `zeta(rho)` or, with `mRhoOfZeta`, `rho(zeta)`.
 *
 * A tapered or variable-pitch thread (roadmap 6.3.3) generalises the advance
 * to `g(theta) = K theta + Q (theta - ThetaRef)^2` and moves the profile out
 * by `Taper g(theta)` as it advances, the thread cut on a cone whose radius
 * grows by `Taper` per unit length along the axis:
 * `P = O + (rho + Taper g) (cos theta E1 + sin theta E2) + (zeta + g) A`.
 */
struct HelicalFit {
    bool mOk = false;
    PatchVec mOrigin{}, mAxis{}, mE1{}, mE2{};
    double mK = 0.0;
    /// The taper (radial growth per unit advance) and the pitch's rate of
    /// change about `mThetaRef`; both 0 for a straight, constant-pitch thread.
    double mTaper = 0.0, mQ = 0.0, mThetaRef = 0.0;
    bool mRhoOfZeta = false;
    /// The profile, over the parameter `s` (`rho`, or `zeta`).
    FeatureLaw mProfile;
    /// The patch's extent in `s` and in the unwrapped angle.
    double mS0 = 0.0, mS1 = 0.0, mThetaMin = 0.0, mThetaMax = 0.0;
    /// RMS and largest distance of the vertices to the profile, in its plane.
    double mRms = 0.0, mMax = 0.0;
    /// The profile is a straight line: a ruled helicoid, a straight flank.
    bool mStraight = false;
    PatchVec At(double S, double Theta) const;
    double Pitch() const;
    /// The local pitch at the unwrapped angle `Theta`.
    double PitchAt(double Theta) const;
    /// The advance `g(Theta)`.
    double Advance(double Theta) const;
};

/**
 * @brief Fit a helical surface to a patch.
 *
 * Without an axis, the screw motion is the one whose velocity field is most
 * nearly tangent to the patch (Pottmann and Randrup, "Rotational and helical
 * surface approximation for reverse engineering", Computing 60, 1998): the
 * smallest eigenvector of a 6 x 6 normal-weighted moment, reduced to 3 x 3. A
 * surface with more symmetry than one screw (a plane, a cylinder, a sphere)
 * or with no advance (a surface of revolution) is refused. With an axis, only
 * the advance is estimated. Both are then refined, axis and advance together,
 * by Levenberg-Marquardt on the profile's residual, and the profile is fitted
 * as a Chebyshev series of the lowest degree within the tolerance.
 *
 * @param pOrigin, pAxis an axis to start from, or both null.
 * @param General also refine a taper and a linearly varying pitch (roadmap
 *        6.3.3) with the screw; false keeps a straight, constant-pitch screw.
 * @param TaperSeed the taper to start from with `General` (the tangent of a
 *        neighbouring cone's half-angle, along the axis it opens towards).
 * @return the fit; `mOk` false as `helicoid_fit` refuses a patch, when the
 *         screw is not determined, or when no profile of degree 8 or less
 *         explains the patch within half the tolerance.
 */
HelicalFit helical_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const PatchVec* pOrigin, const PatchVec* pAxis, double Tol,
                       bool General = false, double TaperSeed = 0.0);

/// The unit normal of a surface at the point of it nearest to a point (either
/// sign); false when there is none.
using NormalAt = std::function<bool(const PatchVec& rX, PatchVec& rNormal)>;

/**
 * @brief A blend between two surfaces whose cross-section is not the rolling
 * ball's circle (roadmap 6.3.2): a chamfer, whose sections are straight, or a
 * conic blend, whose sections are conic arcs tangent to both surfaces -- the
 * rational quadratic `(R1, A, R2; w)`, `A` the apex where the two tangent
 * lines meet, `w` its weight along the rail (with equal legs, `w = sin(beta /
 * 2)`, `beta` the angle at the apex, is the circle; smaller is flatter, larger
 * fuller, and 0 the chord).
 *
 * Its rails are where the patch meets the two surfaces; a cross-section lies
 * in the plane through a point of the first rail normal to the edge the two
 * surfaces would meet along there (`n1 x n2`; the rail's own tangent where
 * the surfaces are nearly parallel), ending where that plane meets the
 * second rail -- the cross-section of a chamfer or blend of the edge.
 */
struct SectionBlendFit {
    bool mOk = false;
    /// Straight cross-sections: a chamfer (of constant or varying width).
    bool mChamfer = false;
    /// The rails, in order; a closed rail's first point is not repeated.
    std::vector<PatchVec> mRail1, mRail2;
    bool mClosed = false;
    /// Arc length at each rail point, and each rail's whole length.
    std::vector<double> mRail1T, mRail2T;
    double mLength = 0.0, mLength2 = 0.0;
    /// The conic's weight along the first rail (unused for a chamfer).
    FeatureLaw mWeight;
    /// The surfaces' normals, which place a section's apex.
    NormalAt mNormal1, mNormal2;
    /// RMS and largest distance of the patch's vertices to their sections.
    double mRms = 0.0, mMax = 0.0;
    /// The normal of the section plane at rail-1 parameter `T`.
    PatchVec Normal(double T) const;
    /// The rail-1 parameter whose section plane holds `rX`.
    double ParameterOf(const PatchVec& rX) const;
    /// The section at rail-1 parameter `T` (wrapped on a closed rail,
    /// extended along the end tangents past an open one): its ends and, for a
    /// conic, its apex. False when the plane misses the second rail or, for
    /// a conic, the tangent lines do not meet ahead of both ends.
    bool Section(double T, PatchVec& rR1, PatchVec& rR2, PatchVec& rApex) const;
    /// The section's point at `U` (0 on the first rail, 1 on the second; the
    /// conic and the chamfer's line continue past both); false as `Section`.
    bool At(double T, double U, PatchVec& rP) const;
};

/**
 * @brief Fit a chamfer or a conic blend to a patch between two surfaces.
 *
 * The rails are the patch's boundary edges along each surface, chained:
 * those @p rEdgeSide puts across it (1 or 2; anything else is neither), or
 * without it those whose ends both lie within the tolerance of it, and
 * nearer it than the other.
 * Each interior vertex is placed in the section through its nearest point of
 * the first rail: within 0.4 of the tolerance of every section's chord, the
 * patch is a chamfer; otherwise each vertex's conic weight follows from its
 * barycentric coordinates in the triangle `(R1, A, R2)`, `w = b1 / (2
 * sqrt(b0 b2))`, and a law of the lowest degree whose conics pass within 0.4
 * of the tolerance of every vertex is the blend.
 *
 * @return the fit; `mOk` false when either rail has fewer than three points,
 *         fewer than six interior vertices fit their sections, or no weight
 *         law explains the patch.
 */
SectionBlendFit section_blend_fit(
    const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris, const ClosestPoint& rS1,
    const NormalAt& rN1, const ClosestPoint& rS2, const NormalAt& rN2, double Tol,
    const std::function<int(std::size_t, std::size_t)>& rEdgeSide = {});

/**
 * @brief The boundary of a corner patch where three or more blends meet (a
 * vertex blend, roadmap 6.3.1): the patch's one closed boundary loop, split
 * into runs by the neighbour across each boundary edge.
 */
struct VertexBlendBoundary {
    bool mOk = false;
    /// The loop's vertices in order, and the neighbour across the edge from
    /// each to the next.
    std::vector<std::size_t> mLoop;
    std::vector<std::int64_t> mAcross;
    /// The distinct neighbours, in the order the loop meets them.
    std::vector<std::int64_t> mNeighbours;
};

/**
 * @brief The boundary loop of a patch and who lies across it.
 * @param rAcross the neighbour across a boundary edge `(a, b)` of local
 *        vertices, or -1 when there is none.
 * @return `mOk` false unless the boundary is one closed loop meeting at least
 *         three neighbours, each along one contiguous run.
 */
VertexBlendBoundary vertex_blend_boundary(
    const std::vector<PatchTri>& rTris,
    const std::function<std::int64_t(std::size_t, std::size_t)>& rAcross);

}  // namespace detail
}  // namespace meshioplusplus
