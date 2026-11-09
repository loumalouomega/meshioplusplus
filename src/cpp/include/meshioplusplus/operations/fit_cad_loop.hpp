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
 * @file operations/fit_cad_loop.hpp
 * @brief Closing the reverse-engineering loop (roadmap 7.7.1, 7.7.2): how much
 * of each face of a fitted solid the scan supports, and a driver that refits
 * towards a target.
 *
 * `face_support(solid, scan)` samples each face of a `fit_cad` re-mesh (its
 * `fit_cad:face` cell data; each triangle's four sub-triangle centroids,
 * weighted by area) and measures the samples against the scan. A face is
 * `Measured` when at least `mMeasuredShare` of its area lies within the
 * tolerance of the scan, `Inferred` below `mInferredShare` (no scan near:
 * hidden interiors, the far side of a through-hole, a face the solid closed
 * by construction) and `Partial` between. Scan cells marked as fill
 * (`repair:hole` at least 0, `repair`'s provenance) count as unmeasured, so a
 * hole `repair` closed before `fit_cad` reads as inferred. The solid gains the
 * field data `fit_cad:face_support` (Int64, per face, 0 measured, 1 partial, 2
 * inferred) and `fit_cad:face_measured` (Float64, the share), and the cell data
 * `fit_cad:support` (its face's class); the report gives the counts and the
 * share of the scan's points within the tolerance of the solid.
 *
 * `fit_cad_iterate(scan)` segments the scan (`segment`'s primitive method),
 * fits (`fit_primitives`) and builds the solid (`fit_cad`) in rounds: each
 * round measures the share of the scan within tolerance of the solid (and,
 * with `mTargetExplained`, the share of the skin's area analytic faces
 * explain), stops at the target, at no improvement or after `mMaxRounds`, and
 * otherwise re-segments the patches of the worst faces (the largest face
 * deviation beyond the tolerance first, then the largest faceted faces) with a
 * finer tolerance, splicing their new labels in. `mPartClass` sets the target
 * and the tolerances (`Analytic` 0.90 and `1e-3` of the diagonal, `Cast` 0.85
 * and `2e-3`, `Freeform` 0.75 and `4e-3`). The best round's solid is kept;
 * deterministic: the same input and options give the same rounds.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"

namespace meshioplusplus {

/// Field data (Int64, per face): 0 measured, 1 partial, 2 inferred.
inline constexpr const char* kFitCadFaceSupportName = "fit_cad:face_support";
/// Field data (Float64, per face): the share of its area within tolerance of the scan.
inline constexpr const char* kFitCadFaceMeasuredName = "fit_cad:face_measured";
/// Cell data (Int64): the support class of each triangle's face.
inline constexpr const char* kFitCadSupportName = "fit_cad:support";

struct FaceSupportOptions {
    /// Within this distance of the scan a sample is measured; at most 0 takes
    /// `1e-3` of the scan's bounding-box diagonal.
    double mTolerance = 0.0;
    double mMeasuredShare = 0.9;  ///< at least this share: measured
    double mInferredShare = 0.1;  ///< below this share: inferred
    /// Scan cell data naming fill cells (value >= 0), which count as unmeasured.
    std::string mHoleData = "repair:hole";
};

struct FaceSupportResult {
    Mesh mMesh;  ///< the solid with the support arrays
    std::vector<double> mFaceShare;
    std::vector<std::int64_t> mFaceClass;
    std::int64_t mNumMeasured = 0;
    std::int64_t mNumPartial = 0;
    std::int64_t mNumInferred = 0;
    /// The share of the scan's points (on measured cells) within the tolerance of the solid.
    double mShareWithin = 0.0;
    double mTolerance = 0.0;
};

/**
 * @brief Class each face of a `fit_cad` re-mesh by the scan's support (see
 * the file comment).
 * @throws std::invalid_argument on a solid without `fit_cad:face` or a scan
 * with no surface cell outside the fills.
 */
MESHIOPLUSPLUS_API FaceSupportResult face_support(const Mesh& rSolid, const Mesh& rScan,
                                                  const FaceSupportOptions& rOptions = {});

/// A part class, setting `fit_cad_iterate`'s target and tolerances.
enum class FitCadPartClass : std::int32_t { Any = 0, Analytic = 1, Cast = 2, Freeform = 3 };

/// `"any"`, `"analytic"`, `"cast"` or `"freeform"`; anything else throws naming them.
MESHIOPLUSPLUS_API FitCadPartClass fit_cad_part_class_from_name(const std::string& rName);

struct FitCadIterateOptions {
    /// Each round's `fit_cad` (its `mFaceDeviation` is turned on).
    FitCadOptions mFitCad;
    /// The first segmentation's tolerance, relative to the scan's diagonal;
    /// at most 0 takes the part class's (`1e-3` for `Any`).
    double mSegmentTolerance = 0.0;
    /// `fit_primitives`' relative tolerance; at most 0 takes the segment one.
    double mFitTolerance = 0.0;
    /// The share of the scan within tolerance to stop at; at most 0 takes the
    /// part class's (0.9 for `Any`).
    double mTargetWithin = 0.0;
    /// The share of the skin's area analytic faces must also explain; 0 off.
    double mTargetExplained = 0.0;
    FitCadPartClass mPartClass = FitCadPartClass::Any;
    std::int64_t mMaxRounds = 4;      ///< at least 1
    std::int64_t mWorstFaces = 3;     ///< faces re-segmented per round; at least 1
    double mRefineFactor = 0.25;      ///< the segment tolerance's factor per round, in (0, 1)
    std::int64_t mMinPatchCells = 1;  ///< `segment`'s smallest patch
    FaceSupportOptions mSupport;
};

/// One round of `fit_cad_iterate`.
struct FitCadIterateRound {
    double mShareWithin = 0.0;
    double mExplained = 0.0;  ///< the solid's area on analytic or feature faces, a share
    double mHausdorffMax = 0.0;
    double mSegmentTolerance = 0.0;
    std::int64_t mNumPatches = 0;
    std::int64_t mNumFaces = 0;
    /// The patches this round re-segmented for the next.
    std::vector<std::int64_t> mResegmented;
};

struct FitCadIterateResult {
    FitCadResult mFit;           ///< the best round's solid
    FaceSupportResult mSupport;  ///< its support (`mMesh` the solid with the arrays)
    Mesh mSkin;                  ///< the best round's fitted skin
    std::vector<FitCadIterateRound> mRounds;
    std::int64_t mBestRound = 0;
    double mTargetWithin = 0.0;
    bool mTargetMet = false;
    /// `"target"`, `"no improvement"`, `"rounds"` or `"nothing to refine"`.
    std::string mStopReason;
};

/**
 * @brief Segment, fit and build in rounds towards a target (see the file
 * comment); with @p rStepPath, the best solid is written there as STEP.
 * @throws std::invalid_argument on an option out of range, and whatever
 * `segment`, `fit_primitives` and `fit_cad` throw.
 */
MESHIOPLUSPLUS_API FitCadIterateResult fit_cad_iterate(const Mesh& rScan,
                                                       const std::string& rStepPath = {},
                                                       const FitCadIterateOptions& rOptions = {});

}  // namespace meshioplusplus
