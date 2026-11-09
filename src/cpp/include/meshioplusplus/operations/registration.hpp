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
 * @brief Rigid registration by iterative closest points (roadmap 7.5.3, 7.3):
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
 * `RegistrationInitial::Auto` drops the need for a start (roadmap 7.3.1):
 * the principal frames of the two inputs (the area-weighted second moments
 * of their surface cells, or of their points for a cloud) are laid onto each
 * other by each of the 24 proper rotations that map one frame's signed axes
 * onto the other's (the axes' signs, and their order where spreads are
 * close, are ambiguous), point-to-point ICP runs from each on a strided
 * subsample of at most 1000 source points at a coarse tolerance, and the
 * candidate of the lowest trimmed RMS (over the closest 90% of its pairs;
 * ties to the lowest index) is refined with the options given.
 * `mSourceRegion` and `mTargetRegion` restrict the pairs to a part's true
 * surface (excluding a fixture or a turntable); `mAdaptiveTrim` rejects, at
 * each iteration, the pairs beyond a multiple of the median pair distance
 * (after D. Chetverikov et al., "The trimmed iterative closest point
 * algorithm", ICPR 2002).
 *
 * `register_landmarks` is the weighted Horn solve on picked pairs, the start
 * a viewer or a script gives `register_rigid` (7.3.2). `find_symmetry_plane`
 * (7.3.6) mirrors the input across each principal plane through its centroid
 * (and an optional seed), registers the image to the original by ICP, and
 * reads the symmetry plane off the product of the reflection and the motion
 * (an improper isometry, projected onto the nearest reflection).
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
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The error `register_rigid` minimises.
enum class RegistrationMethod : std::int32_t {
    PointToPoint = 0,
    PointToPlane = 1,
};

/// Where `register_rigid` starts.
enum class RegistrationInitial : std::int32_t {
    /// From `mInitial`.
    Matrix = 0,
    /// From the best of the 24 principal-frame candidates (see the file comment).
    Auto = 1,
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
    /// `Auto` ignores `mInitial` and starts from the principal frames.
    RegistrationInitial mInitialMode = RegistrationInitial::Matrix;
    /// With `Auto`: warn when a candidate of a pose more than 5 degrees from
    /// the chosen one scores within `mSymmetricTolerance` of it (a symmetric
    /// part, whose pose the surfaces cannot decide).
    bool mSymmetricCheck = false;
    /// The relative trimmed-RMS margin of `mSymmetricCheck`.
    double mSymmetricTolerance = 0.05;
    /// Only the points of the source's cells in this cell region are moved
    /// into pairs; empty uses every source point.
    std::string mSourceRegion;
    /// Only the target's cells in this cell region are matched; empty uses
    /// every target cell. Needs a target with surface cells.
    std::string mTargetRegion;
    /// Drop, at each iteration, the pairs farther apart than this multiple of
    /// the median pair distance (adaptive trimming); 0 keeps them. Applied
    /// together with `mMaxDistance` when both are set.
    double mAdaptiveTrim = 0.0;
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
    /// Source points left out of the pairs at the end (by `mMaxDistance`,
    /// `mAdaptiveTrim` or a missing normal).
    std::int64_t mNumRejected = 0;
    /// With `RegistrationInitial::Auto`: each candidate's trimmed RMS after
    /// its coarse fit, and the index of the one refined (-1 otherwise).
    std::vector<double> mCandidateRms;
    std::int64_t mChosenCandidate = -1;
    std::vector<std::string> mWarnings;
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

/** @brief `"matrix"` or `"auto"`; anything else throws naming them. */
MESHIOPLUSPLUS_API RegistrationInitial registration_initial_from_name(const std::string& rName);

/// What `register_landmarks` found.
struct LandmarkResult {
    /// Row-major 4x4: source coordinates onto the target's.
    std::array<double, 16> mMatrix{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    /// Per pair, the distance between the moved source point and its target.
    std::vector<double> mResiduals;
    /// The weighted root-mean-square of the residuals.
    double mRms = 0.0;
};

/**
 * @brief The rigid motion laying picked source points onto their target
 * points: the weighted Horn solve, in closed form (no iteration).
 * @param rSource the source points.
 * @param rTarget their targets, as many.
 * @param rWeights one non-negative weight per pair, or empty for all ones.
 * @throws std::invalid_argument on fewer than three pairs, mismatched sizes,
 *         a non-finite coordinate or weight, a negative weight, fewer than
 *         three pairs of positive weight, or collinear source points (the
 *         rotation about their line is then undetermined).
 */
MESHIOPLUSPLUS_API LandmarkResult register_landmarks(
    const std::vector<std::array<double, 3>>& rSource,
    const std::vector<std::array<double, 3>>& rTarget, const std::vector<double>& rWeights = {});

/// Options for `find_symmetry_plane`.
struct SymmetryPlaneOptions {
    /// A seed plane tried beside the three principal planes: a point of it and
    /// its normal; a zero normal (the default) tries no seed.
    std::array<double, 3> mSeedOrigin{{0.0, 0.0, 0.0}};
    std::array<double, 3> mSeedNormal{{0.0, 0.0, 0.0}};
    /// A sample whose mirror image lies within this distance of the surface
    /// counts as symmetric; 0 uses 0.005 of the bounding-box diagonal.
    double mTolerance = 0.0;
    /// Refuse a plane whose share of symmetric samples is below this.
    double mMinShare = 0.9;
    /// The samples of the fine fit and of the share (a stride over the points).
    std::int64_t mMaxSamples = 4000;
    /// The iteration cap of each ICP run.
    std::int64_t mMaxIterations = 50;
};

/// What `find_symmetry_plane` found.
struct SymmetryPlaneResult {
    /// The plane: the input's centroid projected onto it, and its unit normal
    /// (signed so its largest-magnitude component is positive).
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};
    std::array<double, 3> mNormal{{0.0, 0.0, 0.0}};
    /// The RMS distance of the mirrored samples to the surface, and the share
    /// within `mTolerance`.
    double mRms = 0.0;
    double mShare = 0.0;
    double mTolerance = 0.0;  ///< the tolerance applied
    /// Each candidate's RMS after its coarse fit (the principal planes, from
    /// the largest spread to the least, then the seed), and the one kept.
    std::vector<double> mCandidateRms;
    std::int64_t mChosenCandidate = 0;
    std::vector<std::string> mWarnings;
};

/**
 * @brief The plane a mesh or point cloud is most nearly symmetric about (see
 * doc/registration.md): each candidate plane's mirror image is registered to
 * the original by ICP, and the plane is recovered from the product of the
 * reflection and the rigid motion.
 * @throws std::invalid_argument on an empty or non-finite input, a bad
 *         option, a target with volume cells, or when the best plane's share
 *         is below `mMinShare` (naming both).
 */
MESHIOPLUSPLUS_API SymmetryPlaneResult
find_symmetry_plane(const Mesh& rMesh, const SymmetryPlaneOptions& rOptions = {});

}  // namespace meshioplusplus
