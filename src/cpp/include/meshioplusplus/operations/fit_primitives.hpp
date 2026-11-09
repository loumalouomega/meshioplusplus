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
 *  - `fit:patch_tolerance` (Float64): the absolute tolerance applied;
 *  - `fit:patch_form_error` (Float64): the peak-to-peak form error, the largest
 *    positive minus the largest negative signed orthogonal distance of the
 *    patch's points (flatness, sphericity, cylindricity under the least-squares
 *    reference);
 *  - `fit:patch_diameter_min`, `fit:patch_diameter_max` (Float64): a sphere's
 *    or cylinder's smallest and largest diameter the points reach (twice their
 *    least and largest radial distance);
 *  - `fit:patch_length` (Float64): a cylinder's or cone's axial length, the
 *    points' extent along the axis;
 *  - `fit:patch_arc` (Float64): a cylinder's or cone's arc span in degrees,
 *    360 less the largest angular gap between the points about the axis (360
 *    when that gap is at most twice the median one);
 *  - `fit:patch_extent` (Float64, `(m, 2)`): a plane's two in-plane principal
 *    extents, largest first.
 *
 * A value that does not apply to the patch's type (or a patch of type none) is
 * NaN.
 *
 *  - `fit:patch_sigma` (Float64): the RMS orthogonal distance of the points
 *    the fit used;
 *  - `fit:patch_num_used`, `fit:patch_num_rejected` (Int64): the points the
 *    fit used and the points sigma trimming left out.
 *
 * ### Measurement-grade fits (roadmap 7.1.1, 7.1.2)
 *
 * **Region and seed modes.** `mRegion` fits one named cell region; `mSeeds`
 * grows a region from the facet nearest each seed point while one primitive
 * explains it (the growth `segment`'s primitive method runs from its own
 * seeds, `detail/primitive_grow.hpp`), unions the regions of several seeds
 * (a bore interrupted by a cross hole) and peels their rim; the cells are
 * written as the region `mOutputRegion`. Either way the fit is one patch, the
 * usual `fit:` rows with patch id 0.
 *
 * **Sigma trimming.** With `mUsedPoints` `Sigma3`, `Sigma2` or `Sigma1` the
 * patch is refitted on the points within k times the RMS distance of the
 * points kept, until that set stops changing (at most 20 rounds), so a scan's
 * stray points, a burr or a fixture do not drag the fit. A plane or sphere
 * fit (when that is the one type asked for) starts from a least median of
 * squares consensus (Rousseeuw 1984) over a fixed schedule of minimal subsets
 * drawn by stride from the Hilbert-ordered points, with no random sampling;
 * the other types start from the least-squares fit. The share of points
 * within tolerance is then over the points used.
 *
 * Patches are fitted in parallel, each serially, so the output is
 * byte-identical on every backend and thread count. The numpy twin follows the
 * same steps and agrees within rounding.
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

/// Cell data (Int64): the fitted type of each cell's patch; see `FitPrimitiveType`.
inline constexpr const char* kFitTypeName = "fit:type";
/// Cell data (Float64): the RMS distance of the cell's patch to its fit.
inline constexpr const char* kFitRmsName = "fit:rms";
/// Cell data (Float64): the largest distance of the cell's patch to its fit.
inline constexpr const char* kFitMaxName = "fit:max";

/// Field data (Float64): each patch's peak-to-peak form error.
inline constexpr const char* kFitPatchFormErrorName = "fit:patch_form_error";
/// Field data (Float64): a sphere's or cylinder's smallest diameter.
inline constexpr const char* kFitPatchDiameterMinName = "fit:patch_diameter_min";
/// Field data (Float64): a sphere's or cylinder's largest diameter.
inline constexpr const char* kFitPatchDiameterMaxName = "fit:patch_diameter_max";
/// Field data (Float64): a cylinder's or cone's axial length.
inline constexpr const char* kFitPatchLengthName = "fit:patch_length";
/// Field data (Float64): a cylinder's or cone's arc span in degrees.
inline constexpr const char* kFitPatchArcName = "fit:patch_arc";
/// Field data (Float64, `(m, 2)`): a plane's two in-plane extents.
inline constexpr const char* kFitPatchExtentName = "fit:patch_extent";
/// Field data (Float64): the RMS orthogonal distance of the points used.
inline constexpr const char* kFitPatchSigmaName = "fit:patch_sigma";
/// Field data (Int64): the points the fit used.
inline constexpr const char* kFitPatchNumUsedName = "fit:patch_num_used";
/// Field data (Int64): the points sigma trimming left out.
inline constexpr const char* kFitPatchNumRejectedName = "fit:patch_num_rejected";

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

/// Which points a fit uses (roadmap 7.1.2).
enum class FitUsedPoints : std::int32_t {
    /// `Sigma3` in the region and seed modes, `All` in the patch mode.
    Auto = 0,
    /// Every point, as before roadmap 7.1.2.
    All = 1,
    /// Refit until the points within 3, 2 or 1 times the RMS distance of the
    /// points kept are the points kept.
    Sigma3 = 2,
    Sigma2 = 3,
    Sigma1 = 4,
};

/// The reference a fit is measured against (roadmap 7.1.4).
enum class FitCriterion : std::int32_t {
    /// The least-squares fit, as before roadmap 7.1.4.
    LeastSquares = 0,
    /// The narrowest zone holding every point (two parallel planes, two
    /// concentric spheres, two coaxial cylinders): the Chebyshev reference
    /// GD&T form tolerances are defined against.
    MinimumZone = 1,
    /// The smallest enclosing sphere or cylinder (a shaft's envelope).
    MinCircumscribed = 2,
    /// The largest sphere or cylinder no point lies inside (a hole's envelope).
    MaxInscribed = 3,
};

/// `FitCriterion` from `"least_squares"`, `"minimum_zone"`,
/// `"min_circumscribed"` or `"max_inscribed"`.
/// @throws std::invalid_argument on any other name.
MESHIOPLUSPLUS_API FitCriterion fit_criterion_from_name(const std::string& rName);

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

    // ---- roadmap 7.1 (ABI 43) ------------------------------------------------

    /// Write the roadmap 7.1 columns: the form error and extremes
    /// (`fit:patch_form_error` to `fit:patch_extent`) and `fit:patch_sigma`,
    /// `fit:patch_num_used`, `fit:patch_num_rejected`. `fit_cad` and
    /// `refine_patches` turn it off: the tables they carry never had them.
    bool mExtendedReport = true;
    /// Region mode: fit the cells of this named `Cell` region as one patch
    /// (`mPatchData` is then unused).
    std::string mRegion;
    /// Seed mode: from the facet nearest each point, grow across smooth edges
    /// while one primitive of `mTypes` (the type to fit, or several for "auto")
    /// explains the region within the tolerance, union the regions and peel
    /// their rim; the result is fitted as one patch and written as the cell
    /// region `mOutputRegion`. The tolerance, when relative, is a fraction of
    /// the whole surface's diagonal.
    std::vector<std::array<double, 3>> mSeeds;
    /// Seed mode: edges whose faces' normals turn by more are not crossed.
    double mFeatureAngleDeg = 30.0;
    /// Seed mode: boundary faces of the grown region whose normal turns from
    /// the fitted surface's by more are peeled off, ring by ring (a fillet or
    /// a flared bore mouth); 180 keeps every one.
    double mRimAngleDeg = 15.0;
    /// Seed mode: the name of the cell region the fitted cells are written as.
    std::string mOutputRegion = "fit:feature";
    /// Which points the fit uses (roadmap 7.1.2): `Auto` is `Sigma3` in the
    /// region and seed modes and `All` in the patch mode.
    FitUsedPoints mUsedPoints = FitUsedPoints::Auto;
    /// The reference (roadmap 7.1.4): after the least-squares fit is
    /// accepted, a plane, sphere or cylinder patch is refitted under it on the
    /// points used (`detail/zone_fit.hpp`); the parameters and the form
    /// columns are then the criterion's, the residual columns the
    /// least-squares fit's. A plane takes only the minimum zone; a cone or
    /// torus keeps its least-squares fit, with a warning.
    FitCriterion mCriterion = FitCriterion::LeastSquares;
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
    /// Seed mode: the cells written as `mOutputRegion` (0 otherwise).
    std::int64_t mNumFeatureCells = 0;
    /// What the run could not do (a seed nothing grew from, ...).
    std::vector<std::string> mWarnings;
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
