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
 * fitted to a patch.
 *
 * Not installed: a core-private helper of `operations/fit_cad.cpp`.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
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

}  // namespace detail
}  // namespace meshioplusplus
