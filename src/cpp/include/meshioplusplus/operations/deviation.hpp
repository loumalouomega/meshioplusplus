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
 * @file operations/deviation.hpp
 * @brief Deviation analysis of a scan against its intended geometry: the
 * signed distance of every scan point to a reference surface or to a fitted
 * primitive, and the statistics an inspection reads from it (roadmap 7.4).
 *
 * ### The distance
 *
 * Against a reference **surface** (`deviation`), each scan point's distance is
 * the signed distance to the nearest reference triangle, signed by the
 * angle-weighted pseudonormal of the feature it lands on, as
 * `distance_to_surface` signs it: positive outside the reference (excess
 * material), negative inside (missing material).
 *
 * **The far-side rule.** A scan point's nearest reference triangle may be the
 * back face of a thin wall. With `mMaxNormalAngleDeg` below 180, a reference
 * triangle whose unit normal turns away from the scan point's own normal by
 * more than that angle is skipped inside the nearest-triangle search
 * (`detail::query_surface_projections_filtered`), so the answer is the nearest
 * triangle that faces the same way, still under the (distance squared,
 * triangle id) order. A scan point's normal is its `normals` point data when
 * the scan carries one, else the angle-weighted normal of the scan's own
 * triangles at it; a point with neither (a bare cloud, an orphan) is measured
 * without the rule.
 *
 * Against a fitted **primitive** (`deviation_to_primitive`), the distance is
 * the analytic signed orthogonal distance of `detail::primitive_signed_distance`
 * (positive along the primitive's normal), with no reference mesh and no
 * alignment.
 *
 * ### What comes back
 *
 * The scan, bit-identical, with the Float64 point data `deviation:distance`
 * (NaN for a point no accepted triangle reaches), and a report over the
 * measured points: min, max, mean, RMS, the population standard deviation,
 * the 1st, 5th, 50th, 95th and 99th percentiles (linear interpolation between
 * order statistics), a histogram of `mBins` equal bins over the range, the
 * share within `+-mTolerance`, and a suggested colour range: the 95th
 * percentile of `|d|` rounded up to two significant figures.
 *
 * Deterministic: each point's answer is a total-order minimum and the
 * statistics are summed in point order, so the result is identical across
 * parallel backends and thread counts.
 *
 * @see doc/deviation.md
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

/// Point data (Float64): each scan point's signed deviation.
inline constexpr const char* kDeviationDistanceName = "deviation:distance";

/// Options of `deviation` and `deviation_to_primitive`.
struct DeviationOptions {
    /// The share within `+-mTolerance` is reported when positive (NaN otherwise).
    double mTolerance = 0.0;
    /// Number of histogram bins, at least 1.
    std::int64_t mBins = 20;
    /// The histogram range; when `mRangeMin >= mRangeMax` (the default), the
    /// symmetric range `+-max |d|`. Values outside it are counted in
    /// `mNumOutsideRange`, not in a bin.
    double mRangeMin = 0.0;
    double mRangeMax = 0.0;
    /// The far-side rule's angle, degrees in `[0, 180]`; 180 turns it off.
    double mMaxNormalAngleDeg = 90.0;
    /// Restrict the reference surface to this named `Cell` region.
    std::string mReferenceRegion;
    /// Attach `deviation:distance` to the returned scan.
    bool mAttach = true;
    /// Register the scan to the reference first (`deviation` only): the
    /// automatic start, point-to-plane ICP and adaptive trimming at three
    /// median pair distances (`register_rigid`, roadmap 7.3.1 and 7.3.3),
    /// restricted to `mReferenceRegion`. The returned scan is then the
    /// registered one, its vector and tensor point data rotated with it.
    bool mRegister = false;
};

/// The measured scan and its statistics.
struct DeviationResult {
    /// The scan, bit-identical, with `deviation:distance` when asked for.
    Mesh mMesh;
    std::int64_t mNumPoints = 0;    ///< scan points
    std::int64_t mNumMeasured = 0;  ///< of those, points with a finite distance
    double mMin = 0.0;              ///< the statistics, over the measured points
    double mMax = 0.0;
    double mMean = 0.0;
    double mRms = 0.0;
    double mSigma = 0.0;  ///< population standard deviation about the mean
    /// The 1st, 5th, 50th, 95th and 99th percentiles.
    std::array<double, 5> mPercentiles{{0.0, 0.0, 0.0, 0.0, 0.0}};
    /// The histogram: `mBins + 1` edges and `mBins` counts.
    std::vector<double> mHistogramEdges;
    std::vector<std::int64_t> mHistogramCounts;
    std::int64_t mNumOutsideRange = 0;  ///< measured points outside the range
    /// The share of measured points within `+-mTolerance` (NaN without one).
    double mShareWithin = 0.0;
    /// The 95th percentile of `|d|`, rounded up to two significant figures.
    double mSuggestedRange = 0.0;
    /// Points measured without the far-side rule (no normal to compare).
    std::int64_t mNumWithoutNormal = 0;
    std::vector<std::string> mWarnings;
    /// With `mRegister`: the registration applied (row-major 4x4, scan
    /// coordinates onto the reference's; the identity otherwise) and its
    /// final RMS pair distance.
    bool mRegistered = false;
    std::array<double, 16> mMatrix{{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}};
    double mRegistrationRms = 0.0;
};

/**
 * @brief The signed deviation of every point of @p rScan from the surface
 * @p rReference.
 * @param rScan a point cloud or a surface (its points are measured).
 * @param rReference a surface mesh (triangles, quads, polygons).
 * @throws std::invalid_argument on a bad option, a volume or higher-order
 *         reference block or an unknown region.
 */
MESHIOPLUSPLUS_API DeviationResult deviation(const Mesh& rScan, const Mesh& rReference,
                                             const DeviationOptions& rOptions = {});

/**
 * @brief The signed deviation of every point of @p rScan from a fitted
 * primitive: `Type` 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus, and its
 * ten `fit:patch_params` entries.
 * @throws std::invalid_argument on an unknown type or non-finite parameters.
 */
MESHIOPLUSPLUS_API DeviationResult deviation_to_primitive(const Mesh& rScan, std::int64_t Type,
                                                          const std::array<double, 10>& rParams,
                                                          const DeviationOptions& rOptions = {});

/// Cell data (Int64): the index of each cell's deviation region in
/// `DeviationRegionsResult::mRegions`, -1 for a cell within tolerance.
inline constexpr const char* kDeviationRegionName = "deviation:region";

/// Options of `deviation_regions`.
struct DeviationRegionsOptions {
    /// A cell beyond `+mTolerance` is excess material, beyond `-mTolerance`
    /// missing material; positive and finite.
    double mTolerance = 0.0;
    /// The point data holding the signed deviation (`deviation` writes it).
    std::string mField = kDeviationDistanceName;
    /// Region names are `<prefix>excess_<k>` and `<prefix>missing_<k>`.
    std::string mPrefix = "deviation:";
    /// Drop regions smaller than this area (0 keeps every one).
    double mMinArea = 0.0;
};

/// One connected region beyond tolerance.
struct DeviationRegionRecord {
    std::string mName;       ///< its cell region's name
    std::int64_t mSign = 0;  ///< +1 excess material, -1 missing material
    std::int64_t mNumCells = 0;
    double mArea = 0.0;
    double mExtreme = 0.0;  ///< the largest (excess) or least (missing) cell value
    double mMean = 0.0;     ///< area-weighted mean deviation
    double mVolume = 0.0;   ///< the signed volume estimate, the sum of d * area
    std::array<double, 3> mCentroid{{0.0, 0.0, 0.0}};               ///< area-weighted
    std::array<double, 6> mBounds{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};  ///< lo x, y, z, hi x, y, z
};

/// What `deviation_regions` found.
struct DeviationRegionsResult {
    /// The input with one named cell region per record and `deviation:region`.
    Mesh mMesh;
    /// Largest `|mVolume|` first (ties: the region holding the lower cell).
    std::vector<DeviationRegionRecord> mRegions;
    std::int64_t mNumExcess = 0;
    std::int64_t mNumMissing = 0;
    std::int64_t mNumDropped = 0;  ///< regions below `mMinArea`
    std::vector<std::string> mWarnings;
};

/**
 * @brief The connected regions (cells joined across shared edges) of a
 * measured surface whose deviation lies beyond `+-mTolerance`: excess and
 * missing material, each with its area, extreme and mean deviation, signed
 * volume estimate, centroid and bounding box, written as named cell regions
 * (roadmap 7.4.3).
 *
 * A cell's deviation is the mean of its corners' values of `mField` (corners
 * with NaN are left out; a cell with none finite is within tolerance).
 * @throws std::invalid_argument on a bad option, a missing or non-scalar
 *         field, or a volume block.
 */
MESHIOPLUSPLUS_API DeviationRegionsResult
deviation_regions(const Mesh& rMesh, const DeviationRegionsOptions& rOptions);

/// Options of `exaggerate`.
struct ExaggerateOptions {
    /// The signed deviation point data to exaggerate.
    std::string mField = kDeviationDistanceName;
    /// Points move along their normals by `mScale * d`.
    double mScale = 10.0;
    /// `d` is clamped to `+-mClamp` first; 0 (the default) clamps at the 95th
    /// percentile of `|d|` rounded up to two significant figures, the colour
    /// range `deviation` suggests.
    double mClamp = 0.0;
    /// `d` is then averaged over the points within this fraction of the
    /// bounding-box diagonal (weights `1 - (r / R)^2`), so scan noise does not
    /// become fur; 0 leaves it unsmoothed.
    double mSmoothFraction = 0.01;
};

/// The exaggerated mesh.
struct ExaggerateResult {
    Mesh mMesh;
    double mClamp = 0.0;         ///< the clamp applied
    double mRadius = 0.0;        ///< the smoothing radius applied
    std::int64_t mNumMoved = 0;  ///< points with a normal and a finite deviation
    std::vector<std::string> mWarnings;
};

/**
 * @brief A copy of the surface @p rMesh with each point displaced along its
 * angle-weighted unit normal by `mScale` times its clamped, smoothed
 * deviation (roadmap 7.4.4): a mesh for a figure or an animation, the
 * measured scan left untouched. Connectivity and data are unchanged.
 * @throws std::invalid_argument on a missing or non-scalar field, a
 *         non-finite scale or a negative clamp or fraction.
 */
MESHIOPLUSPLUS_API ExaggerateResult exaggerate(const Mesh& rMesh,
                                               const ExaggerateOptions& rOptions = {});

}  // namespace meshioplusplus
