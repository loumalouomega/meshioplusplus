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
 * @file detail/kinematic_fit.hpp
 * @brief A patch no primitive explains, fitted as a kinematic surface: a
 * surface of linear extrusion, of revolution or a helical surface, each the
 * sweep of a planar profile by a motion (`fit_cad`'s `kinematic_surfaces`,
 * doc/fit_cad_freeform_plan.md step 5).
 *
 *  - **Extrusion.** The patch's normals are all perpendicular to one
 *    direction `d`: the smallest eigenvector of their area-weighted second
 *    moment `sum a n n^T`. The profile is the points projected on the plane
 *    normal to `d`.
 *  - **Revolution.** The normals all meet one axis: the linear line complex
 *    of the normal lines (Pottmann and Randrup, Computing 60, 1998;
 *    `screw_estimate`) has no advance. The profile is the points in the
 *    meridian half-plane, `(rho, z)`.
 *  - **Helical.** `helical_fit` (an advancing screw, any profile).
 *
 * A profile is a cubic B-spline curve through the points ordered along it
 * (by angle about their centre or along their principal direction, the
 * shorter path), averaged into bins and checked against every point by
 * projection; a profile that is a line or a circular arc is refused, those
 * being the plane, cylinder, cone, sphere and torus `fit_primitives` fits.
 *
 * Every fit also carries a B-spline surface over its own parameters
 * (`u` along the profile, `v` along the motion): exactly the extrusion, and
 * within a twentieth of the tolerance the revolution and the helical
 * surface, which a polynomial cannot be. The face is accepted by
 * `freeform_accept` on it (vertices within half the tolerance, facet
 * centres within it, the surface sampled inside its trim), the test every
 * freeform face passes; the exact description goes to the writers.
 *
 * The patch must be a topological disc: `kinematic_is_band` tells a band
 * (an extruded or revolved wall closed round), which the caller cuts first.
 * Serial; the output depends only on the input.
 *
 * Core-private (not installed).
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes (private, not installed)
#include "bspline_curve.hpp"
#include "bspline_lsq.hpp"
#include "feature_fit.hpp"
#include "patch_flatten.hpp"

namespace meshioplusplus::detail {

/// The kinds a patch is tried as, by their `fit_cad:type` values.
enum class KinematicKind : int {
    None = 0,
    Helical = 8,
    Extrusion = 11,
    Revolution = 12,
};

/// What `kinematic_fit` aims for.
struct KinematicOptions {
    double mTolerance = 0.0;  ///< `fit_cad`'s tolerance, > 0
    bool mExtrusion = true;
    bool mRevolution = true;
    bool mHelical = true;
    /// A patch with fewer triangles is not tried.
    std::size_t mMinTriangles = 32;
    /// The surfaces run on past the patch by this share of its extent.
    double mMargin = 0.05;
};

/// A kinematic surface and how far the patch lies from it.
struct KinematicFit {
    bool mOk = false;
    KinematicKind mKind = KinematicKind::None;
    /// The frame: for an extrusion `mAxis` is the direction and the profile
    /// lies in the plane through `mOrigin` spanned by `mE1`, `mE2`; for a
    /// revolution `mAxis` is the axis through `mOrigin`, the profile is in
    /// `(rho, z)` and the angle is measured from `mE1` towards `mE2`.
    PatchVec mOrigin{0.0, 0.0, 0.0};
    PatchVec mAxis{0.0, 0.0, 1.0};
    PatchVec mE1{1.0, 0.0, 0.0};
    PatchVec mE2{0.0, 1.0, 0.0};
    /// The profile, `(x, y, 0)` in the frame's coordinates.
    BsplineCurve mProfile;
    /// The extrusion's extent along `mAxis`, the revolution's angles.
    double mV0 = 0.0, mV1 = 0.0;
    /// The helical fit, for `Helical`.
    HelicalFit mHelical;
    /// The surface over `(profile, motion)` parameters; see the file comment.
    BsplineSurface mSurface;
    /// Per patch point, its parameters on `mSurface`.
    std::vector<double> mU, mV;
    double mMaxDist = 0.0;  ///< largest vertex distance from the surface
    double mRmsDist = 0.0;
    double mMaxCentroidDist = 0.0;
    double mReverseMaxDist = 0.0;  ///< largest distance of a surface sample from the patch
    /// Why a kind was refused, for the debug log.
    std::string mWhy;

    /// A point of the exact surface at profile parameter `T` and motion
    /// parameter `S` (the extrusion's height, the revolution's angle).
    PatchVec At(double T, double S) const;
    /// The profile point `(x, y)` at `T` in 3-D, in the frame.
    PatchVec ProfilePoint(double T) const;
    /// The profile as a 3-D curve (its poles mapped as `ProfilePoint` maps
    /// its points): an extrusion's at height 0, a revolution's at angle 0.
    BsplineCurve SweptCurve() const;
};

/// Whether a patch is a band: one connected surface of genus 0 with two
/// boundary loops (a wall closed round), which `kinematic_fit` refuses.
bool kinematic_is_band(const std::vector<PatchTri>& rTris);

/// Fit a disc patch (see the file comment); `mOk` false, with `mWhy`, when
/// no kind passes.
KinematicFit kinematic_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                           const KinematicOptions& rOptions);

}  // namespace meshioplusplus::detail
