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
 * @file operations/registration.hpp
 * @brief Rigid registration by iterative closest points (roadmap 7.5.3):
 * the rotation and translation that best lay one mesh or point cloud onto
 * another, as the 4x4 matrix `transform` applies. Aligning a scan to its
 * design, and two runs' meshes to each other before comparing them.
 *
 * `register_rigid(source, target, options)` (P. J. Besl and N. D. McKay, "A
 * method for registration of 3-D shapes", IEEE PAMI 14(2), 1992) alternates
 * two steps from `mInitial`:
 *
 * - **Correspondences.** Each source point, moved by the current transform,
 *   is paired with its closest point on the target: on its surface when the
 *   target has surface cells (triangles, quads, polygons, through
 *   `detail/surface_distance`), else with its nearest target point (a cloud).
 *   With `mMaxDistance` > 0 a pair farther apart is dropped (outliers, parts
 *   the target lacks).
 * - **The best rigid motion for those pairs.** `PointToPoint`: the motion
 *   minimising the summed squared distances, in closed form by B. K. P.
 *   Horn's unit quaternions ("Closed-form solution of absolute orientation
 *   using unit quaternions", JOSA A 4(4), 1987). `PointToPlane`: the motion
 *   minimising the summed squared distances to the target's tangent planes
 *   (Y. Chen and G. Medioni 1992), linearised in the rotation (K.-L. Low,
 *   "Linear least-squares optimization for point-to-plane ICP surface
 *   registration", 2004) and applied as an exact rotation; it slides along
 *   flat regions and usually converges in far fewer iterations. The plane is
 *   the target cell's normal, or for a cloud its `normals` point data
 *   (`estimate_normals`), refused by name when absent.
 *
 * It stops when an iteration's motion is below `mTolerance` (the rotation
 * in radians, the translation relative to the source's extent), or after
 * `mMaxIterations`. ICP finds the nearest local optimum: the source must
 * start roughly aligned (`mInitial`, a pre-alignment, or a coarse manual
 * transform), and symmetric shapes can settle in a symmetric pose.
 *
 * The matrix maps source coordinates onto the target's: `x_t = R x_s + t`,
 * row-major `[[R, t], [0, 0, 0, 1]]`. Only points are read: the source's cells
 * and data play no part. The correspondence search runs in parallel, the
 * solves serially; the result is the same on every backend and thread count.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The error `register_rigid` minimises.
enum class RegistrationMethod : std::int32_t {
    PointToPoint = 0,
    PointToPlane = 1,
};

/// How `register_rigid` iterates.
struct RegistrationOptions {
    RegistrationMethod mMethod = RegistrationMethod::PointToPoint;
    /// The iteration cap; at least 1.
    std::int64_t mMaxIterations = 50;
    /// Stop when an iteration rotates by less than this (radians) and
    /// translates by less than this times the source's bounding-box
    /// diagonal.
    double mTolerance = 1e-10;
    /// Drop pairs farther apart than this; 0 keeps every pair.
    double mMaxDistance = 0.0;
    /// The starting transform, row-major 4x4 (a rigid motion).
    std::array<double, 16> mInitial{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
};

/// What `register_rigid` found.
struct RegistrationResult {
    /// Row-major 4x4: source coordinates onto the target's.
    std::array<double, 16> mMatrix{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    /// Root-mean-square pair distance before the first iteration (under
    /// `mInitial`) and after the last.
    double mInitialRms = 0.0;
    double mRms = 0.0;
    /// Iterations run.
    std::int64_t mIterations = 0;
    /// Pairs kept at the end (all source points without `mMaxDistance`).
    std::int64_t mNumPairs = 0;
    /// False when `mMaxIterations` ran out first.
    bool mConverged = false;
};

/** @brief `"point_to_point"` or `"point_to_plane"`; anything else throws naming them. */
MESHIOPLUSPLUS_API RegistrationMethod registration_method_from_name(const std::string& rName);

/**
 * @brief The rigid motion laying @p rSource onto @p rTarget (see the file
 * comment).
 * @throws std::invalid_argument on an empty source or target, a non-finite
 *         coordinate, `mMaxIterations` below 1, a negative tolerance or
 *         distance, `PointToPlane` on a cloud without `normals`, a target with
 *         volume cells (naming `extract_surface`), or when no pair survives
 *         `mMaxDistance`.
 */
MESHIOPLUSPLUS_API RegistrationResult register_rigid(const Mesh& rSource, const Mesh& rTarget,
                                                     const RegistrationOptions& rOptions = {});

}  // namespace meshioplusplus
