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
 * @file detail/patch_network.hpp
 * @brief B-spline patches fitted on a quad layout as one network: shared
 * edge curves, exact G0, G1 by penalty across shared edges, and the
 * extraordinary-vertex modes -- `fit_cad`'s network mode (v17.12.0), free of
 * OpenCASCADE.
 *
 *  1. **Edges.** Every layout chain is fitted once (`bspline_curve_fit`, at
 *     the freeform target `tol / 4`), through its corners and T-junction
 *     points exactly.
 *  2. **Boundaries.** Each patch's four sides are the exact pieces of those
 *     curves (`bspline_curve_segment`, joined C0 where a side runs over
 *     several chains), the opposite sides' knot vectors merged by knot
 *     insertion, so a patch's boundary rows of poles *are* the shared curves
 *     and two patches meet exactly (G0) whatever their interiors.
 *  3. **Interior.** The interior poles are fitted by least squares to the
 *     patch's vertices and facet centres, parametrised on the unit square
 *     (`param_square`) and corrected by projection after each fit, with a
 *     small thin-plate-like smoothing of the pole grid; solved by the
 *     deterministic conjugate gradients of `sparse_solve.hpp`.
 *  4. **G1.** From the second sweep on, each shared side adds rows asking the
 *     patch's cross-boundary derivative to lie in the neighbour's tangent
 *     plane (`n . S_v = 0` at samples along the side, normalised to the sine
 *     of the angle), weighted so an angle at `mG1ToleranceDeg` costs what a
 *     distance at the tolerance does: block Gauss-Seidel over the patches in
 *     index order, the weight doubled (up to eight times its start) after a
 *     sweep that does not lower the largest angle, up to `mMaxSweeps` sweeps.
 *     A patch the G1 rows push out of the tolerance gets knots where it
 *     misses (step 6) and is fitted again before its accurate surface is
 *     restored: next to an extraordinary vertex G1 needs more poles than
 *     accuracy alone.
 *  5. **Extraordinary vertices** (where other than four patches meet):
 *     `G0` leaves them alone (and the G1 measure skips their first span);
 *     `G1` puts the second pole of every curve leaving one, and every
 *     patch's corner twist pole, in a common tangent plane (the area-weighted
 *     normal), so all patches share the tangent plane at the vertex (after
 *     Peters, "Patching Catmull-Clark meshes", SIGGRAPH 2000); `G1x` smooths
 *     the poles round it more; `G1xx` also inserts a double knot in the first
 *     span towards it, for the extra freedom of Loop and Schaefer's
 *     construction ("Approximating Catmull-Clark subdivision surfaces with
 *     bicubic patches", ACM TOG 2008) -- the constructions are the papers',
 *     reimplemented, no code copied.
 *  6. **Refinement.** A patch outside the tolerance gets a knot in the middle
 *     of its worst span in each direction (its boundary curves refined
 *     exactly, so its neighbours do not change), within `mMaxTotalPoles` and
 *     `mMaxRefinements`; each patch is then held to `freeform_accept`, the
 *     adaptive faces' test.
 *
 * A **core-private** header: nothing installed names it. Serial and
 * deterministic.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_cad.hpp"

// Project includes (private, not installed)
#include "bspline_curve.hpp"
#include "bspline_lsq.hpp"
#include "quad_layout.hpp"

namespace meshioplusplus {
namespace detail {

/// What `patch_network_fit` aims for and may spend.
struct PatchNetworkOptions {
    double mTolerance = 0.0;          ///< fit_cad's tolerance, > 0
    double mFreeformTolerance = 0.0;  ///< curves' and vertices' target; 0: tol / 4
    double mSmoothing = 1e-6;         ///< the pole grid's smoothing, relative
    std::size_t mMaxPolesPerDirection = 40;
    std::size_t mMaxTotalPoles = 400;  ///< per patch
    double mG1ToleranceDeg = 1.0;
    FitCadExtraordinaryVertex mExtraordinary = FitCadExtraordinaryVertex::G1x;
    int mMaxSweeps = 6;
    int mMaxRefinements = 6;
    bool mValidate = true;
};

/// One fitted patch.
struct NetworkPatch {
    bool mOk = false;
    FitCadDemotionReason mReason = FitCadDemotionReason::None;
    BsplineSurface mSurface;
    /// Per region point, its parameters on the surface (0 for points the
    /// patch does not use).
    std::vector<double> mU, mV;
    double mMaxDist = 0.0;
    double mRmsDist = 0.0;
    double mMaxCentroidDist = 0.0;
    double mReverseMaxDist = 0.0;
    double mAllowance = 0.0;
    /// The largest angle between its normal and a neighbour's along its
    /// shared sides, in degrees.
    double mMaxG1AngleDeg = 0.0;
    std::int64_t mIterations = 0;
    std::int64_t mRefinements = 0;
};

/// A fitted network.
struct PatchNetwork {
    /// Every patch within the tolerance.
    bool mOk = false;
    std::vector<NetworkPatch> mPatches;
    /// Per layout chain, its curve on `[0, 1]` and its vertices' parameters.
    std::vector<BsplineCurve> mChains;
    std::vector<std::vector<double>> mChainParams;
    double mMaxG1AngleDeg = 0.0;
    /// The largest G1 angle over the accepted patches with the samples next
    /// to extraordinary vertices included, whatever the mode: the measure
    /// two modes are compared by.
    double mMaxG1AngleEvDeg = 0.0;
    /// Chains shared by two patches.
    std::size_t mNumEdges = 0;
    std::size_t mNumExtraordinary = 0;
    int mSweeps = 0;
};

/// Fit the network of a layout (see the file comment). `rTris` the layout's
/// region; `rLayout` must be `mOk`.
PatchNetwork patch_network_fit(const std::vector<PatchVec>& rPts,
                               const std::vector<PatchTri>& rTris, const QuadLayout& rLayout,
                               const PatchNetworkOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
