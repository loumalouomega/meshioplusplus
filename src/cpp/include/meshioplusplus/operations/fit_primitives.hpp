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
 * @file operations/fit_primitives.hpp
 * @brief Fit a plane, sphere, cylinder, cone or torus to each patch of a
 * surface, and report how well it fits -- the measurement half of reverse
 * engineering a triangle skin into a CAD model.
 *
 * The patches are a cell-data label (by default `segment:patch`, which
 * `segment` writes); cells labelled below 0, and cells that are not surface
 * cells, belong to none. An empty `mPatchData` fits the whole surface as one
 * patch. Each patch is fitted to its cells' corner points, deterministically:
 *
 *  - **plane**: the least-squares plane through the centroid (the smallest
 *    eigenvector of the covariance); its normal points the way the cells do;
 *  - **sphere**: an algebraic fit, refined by Levenberg-Marquardt on the
 *    orthogonal distance;
 *  - **cylinder** and **cone**: the axis direction is the direction the cell
 *    normals vary least along, the axis point the least-squares meeting point
 *    of the normal lines; a cone's apex and half-angle come from a line fit of
 *    radius against height. Both are then refined;
 *  - **torus**: the tube-centre circle is found by scanning the minor radius
 *    (cell centroids moved back along their normals must be cocircular), then
 *    refined.
 *
 * A type is **accepted** when its RMS orthogonal distance is at most the
 * tolerance and at least `mMinInlierFraction` of the points lie within the
 * tolerance. The patch takes the first accepted type in the order plane,
 * sphere, cylinder, cone, torus (the simplest that explains it); when none is
 * accepted its type is `none` (0) and the report carries the best attempt's
 * residuals, never a fit rounded to a neat primitive.
 *
 * ### Output
 *
 * The input, unchanged, plus cell data -- `fit:type` (Int64: -1 no patch,
 * 0 none, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus), `fit:rms` and
 * `fit:max` (Float64, the patch's; NaN with no patch) -- and per-patch field
 * data, one row per patch label in ascending order:
 *
 *  - `fit:patch_id` (Int64): the label;
 *  - `fit:patch_type` (Int64): as `fit:type`;
 *  - `fit:patch_params` (Float64, `(m, 10)`): `x, y, z, ax, ay, az, r1, r2,
 *    angle, 0` -- plane: a point and the unit normal; sphere: the centre and
 *    `r1` the radius; cylinder: the axis point nearest the patch centroid, the
 *    unit axis and `r1` the radius; cone: the apex, the unit axis pointing into
 *    the cone and `angle` the half-angle in degrees; torus: the centre, the
 *    unit axis, `r1` the major and `r2` the minor radius. Unused entries are 0;
 *    a patch of type none has NaN throughout. An axis is signed so its
 *    largest-magnitude component is positive (a cone's points into the cone);
 *  - `fit:patch_rms`, `fit:patch_max` (Float64): the accepted fit's residuals,
 *    or the best attempt's;
 *  - `fit:patch_fraction` (Float64): the share of points within tolerance;
 *  - `fit:patch_tolerance` (Float64): the absolute tolerance applied.
 *
 * Patches are fitted in parallel, each serially, so the output is
 * byte-identical on every backend and thread count. The numpy twin follows the
 * same steps and agrees within rounding.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell data (Int64): the fitted type of each cell's patch; see `FitPrimitiveType`.
inline constexpr const char* kFitTypeName = "fit:type";
/// Cell data (Float64): the RMS distance of the cell's patch to its fit.
inline constexpr const char* kFitRmsName = "fit:rms";
/// Cell data (Float64): the largest distance of the cell's patch to its fit.
inline constexpr const char* kFitMaxName = "fit:max";

/// The values of `fit:type` and `fit:patch_type`.
enum class FitPrimitiveType : std::int32_t {
    None = 0,
    Plane = 1,
    Sphere = 2,
    Cylinder = 3,
    Cone = 4,
    Torus = 5,
};

/// `FitPrimitivesOptions::mTypes` bits: `1 << FitPrimitiveType`.
inline constexpr std::uint32_t kFitAllTypes =
    (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4) | (1u << 5);

/// What `fit_primitives` fits and accepts.
struct FitPrimitivesOptions {
    /// The cell-data label naming each cell's patch; empty fits the whole
    /// surface as one patch.
    std::string mPatchData = "segment:patch";
    /// The acceptance tolerance on the orthogonal distance; must be > 0.
    double mTolerance = 1e-3;
    /// Whether `mTolerance` is a fraction of each patch's bounding-box diagonal
    /// (true) or a length (false).
    bool mRelative = true;
    /// The share of a patch's points that must lie within tolerance; in `[0, 1]`.
    double mMinInlierFraction = 0.9;
    /// Which types to try, as `1 << FitPrimitiveType` bits; must name at least one.
    std::uint32_t mTypes = kFitAllTypes;
    /// The cap on Levenberg-Marquardt iterations per refined fit; at least 1.
    std::int64_t mMaxIterations = 50;
};

/// The annotated mesh and how many patches took each type.
struct FitPrimitivesResult {
    /// The input with the `fit:*` cell and field data.
    Mesh mMesh;
    std::int64_t mNumPatches = 0;
    std::int64_t mNumPlane = 0;
    std::int64_t mNumSphere = 0;
    std::int64_t mNumCylinder = 0;
    std::int64_t mNumCone = 0;
    std::int64_t mNumTorus = 0;
    /// Patches no type was accepted for.
    std::int64_t mNumUnfitted = 0;
};

/**
 * @brief Fit analytic primitives to the patches of a surface.
 *
 * @param rMesh a surface mesh carrying the patch labels; lines and vertices
 *        belong to no patch.
 * @param rOptions see `FitPrimitivesOptions`.
 * @throws std::invalid_argument on an option out of range, a missing label
 *         array, or a volume cell block.
 */
MESHIOPLUSPLUS_API FitPrimitivesResult fit_primitives(const Mesh& rMesh,
                                                      const FitPrimitivesOptions& rOptions = {});

}  // namespace meshioplusplus
