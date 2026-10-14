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
 * profile, with or without a known axis.
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
 * radial plane).
 */
struct HelicoidFit {
    bool mOk = false;
    PatchVec mOrigin{}, mAxis{}, mE1{}, mE2{};
    double mK = 0.0;  // axial advance per radian: pitch / (2 pi), signed by the hand
    double mT = 0.0;  // axial slope per unit radius
    double mC = 0.0;  // axial offset
    /// The patch's extent in radius and in the unwrapped angle.
    double mRhoMin = 0.0, mRhoMax = 0.0, mThetaMin = 0.0, mThetaMax = 0.0;
    /// RMS and largest orthogonal distance of the vertices to the surface.
    double mRms = 0.0, mMax = 0.0;
    /// The surface point at `(Rho, Theta)`.
    PatchVec At(double Rho, double Theta) const;
    /// The pitch: the axial advance per turn.
    double Pitch() const;
};

/**
 * @brief Fit a ruled helicoid about a known axis to a patch.
 *
 * The angle of every vertex about the axis is unwrapped continuously over the
 * patch's edges (breadth-first from its lowest vertex), which makes the model
 * linear in `(C, K, T)`: one least-squares solve.
 *
 * @param rOrigin a point on the axis; @param rAxis its direction.
 * @return the fit; `mOk` false when the patch touches the axis, cannot be
 *         unwrapped consistently, spans less than half a turn, or the solve is
 *         singular (a patch with no advance or no radial extent).
 */
HelicoidFit helicoid_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                         const PatchVec& rOrigin, const PatchVec& rAxis);

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
 */
struct HelicalFit {
    bool mOk = false;
    PatchVec mOrigin{}, mAxis{}, mE1{}, mE2{};
    double mK = 0.0;
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
 * @return the fit; `mOk` false as `helicoid_fit` refuses a patch, when the
 *         screw is not determined, or when no profile of degree 8 or less
 *         explains the patch within half the tolerance.
 */
HelicalFit helical_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const PatchVec* pOrigin, const PatchVec* pAxis, double Tol);

}  // namespace detail
}  // namespace meshioplusplus
