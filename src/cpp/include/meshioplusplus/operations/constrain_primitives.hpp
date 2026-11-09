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
 * @file operations/constrain_primitives.hpp
 * @brief Fit related primitives together under explicit control (roadmap
 * 7.7.3): the pass `fit_cad`'s `mConstrainPrimitives` runs (coaxial groups,
 * parallel and perpendicular planes, equal radii; `doc/fit_cad.md`), on a
 * `fit_primitives` skin, with a force and a forbid list of relations between
 * patches, angle and distance tolerances in place of its own, and a report of
 * every relation applied with its residual.
 *
 * Relations name patches by their `segment:patch` label (`fit:patch_id`). A
 * forbidden relation is never applied; a forced one is tried whatever the
 * tolerances say, and kept only when the patches still pass the test they were
 * accepted under (a forced relation that fails is listed in `mUnmet`). With
 * the default options the result is the pass `fit_cad` runs, bit for bit.
 * Run it before `fit_cad` (with `mConstrainPrimitives` off) to control the
 * relations a solid is built with.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// A relation between two fitted primitives.
enum class PrimitiveRelationKind : std::int32_t {
    Coaxial = 0,        ///< one axis (a sphere: its centre on it)
    Parallel = 1,       ///< a plane's normal along an axis or another normal
    Perpendicular = 2,  ///< a plane's normal across an axis or another normal
    EqualRadius = 3     ///< one radius (cylinders, spheres; a torus's minor radius)
};

/// `"coaxial"`, `"parallel"`, `"perpendicular"` or `"equal_radius"`; anything
/// else throws naming them.
MESHIOPLUSPLUS_API PrimitiveRelationKind primitive_relation_from_name(const std::string& rName);
/// The name `primitive_relation_from_name` reads.
MESHIOPLUSPLUS_API const char* primitive_relation_name(PrimitiveRelationKind Kind) noexcept;

/// A relation named by two `segment:patch` labels (unordered).
struct PrimitiveRelation {
    std::int64_t mPatchA = -1;
    std::int64_t mPatchB = -1;
    PrimitiveRelationKind mKind = PrimitiveRelationKind::Coaxial;
};

/// A relation the pass applied.
struct PrimitiveRelationRecord {
    std::int64_t mPatchA = -1;  ///< the reference (the group's seed, the earlier plane or axis)
    std::int64_t mPatchB = -1;  ///< the patch that took it
    PrimitiveRelationKind mKind = PrimitiveRelationKind::Coaxial;
    /// Before the change: the angle between the two directions, in degrees
    /// (for `Perpendicular`, its distance from 90; 0 for `EqualRadius`).
    double mAngleDeg = 0.0;
    /// Before the change: a coaxial member's axis point (a sphere's centre)
    /// from the seed's axis, or the two radii's difference; 0 otherwise.
    double mDistance = 0.0;
    bool mForced = false;
};

struct ConstrainPrimitivesOptions {
    /// The largest distance a face may keep from its patch before a writer
    /// demotes it; at most 0 takes three times the largest
    /// `fit:patch_tolerance` (as `fit_cad`'s default tolerance does).
    double mFaceTolerance = 0.0;
    /// The largest angle, in degrees, between two directions still related;
    /// negative keeps the pass's own (what moves a patch by twice its
    /// tolerance over its extent).
    double mAngleToleranceDeg = -1.0;
    /// The largest distance between two axes (or of a centre from an axis)
    /// still coaxial, and between two radii still equal; negative keeps the
    /// pass's own (twice the patch tolerance; the smaller tolerance).
    double mDistanceTolerance = -1.0;
    /// The share of a patch's points that must stay within its tolerance;
    /// negative takes `FitPrimitivesOptions`' default.
    double mMinFraction = -1.0;
    std::vector<PrimitiveRelation> mForce;
    std::vector<PrimitiveRelation> mForbid;
};

struct ConstrainPrimitivesResult {
    /// The skin with its `fit:patch_*` table (and `fit:rms`, `fit:max`)
    /// rewritten for the patches changed.
    Mesh mMesh;
    std::int64_t mNumConstrained = 0;  ///< patches that took a shared value
    std::int64_t mNumGroups = 0;       ///< groups, families and sets of two or more
    std::vector<PrimitiveRelationRecord> mRelations;
    std::vector<PrimitiveRelation> mUnmet;  ///< forced relations not kept
    double mFaceTolerance = 0.0;            ///< the face tolerance used
};

/**
 * @brief Constrain the primitives of a `fit_primitives` skin (see the file
 * comment).
 * @throws std::invalid_argument on a mesh without the `fit:patch_*` table or
 * `segment:patch`, an unknown patch in a relation, or a tolerance that is not
 * finite.
 */
MESHIOPLUSPLUS_API ConstrainPrimitivesResult
constrain_primitives(const Mesh& rMesh, const ConstrainPrimitivesOptions& rOptions = {});

}  // namespace meshioplusplus
