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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file operations/critical_points.hpp
 * @brief Critical points (`nabla f = 0`) of a scalar field on high-order or
 * linear cells via damped Newton iteration and Hessian eigenvalue
 * classification (roadmap 15.3.3).
 *
 * Returns a `vertex` mesh carrying:
 *  - `point_data["critical:value"]` (Float64 scalar value at the critical point),
 *  - `point_data["critical:type"]` (Int32: 0 = minimum, 1 = maximum, 2 = saddle,
 *    3 = degenerate),
 *  - `point_data["critical:cell"]` (Int64 owning global cell id).
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

inline constexpr const char* kCriticalValueName = "critical:value";
inline constexpr const char* kCriticalTypeName = "critical:type";
inline constexpr const char* kCriticalCellName = "critical:cell";

/// Classification of a scalar-field critical point by Hessian eigenvalue signs.
enum class CriticalPointType : std::int32_t {
    Minimum = 0,
    Maximum = 1,
    Saddle = 2,
    Degenerate = 3,
};

/// Options for `critical_points`.
struct CriticalPointsOptions {
    /// Name of the scalar `point_data` field to analyze (required).
    std::string mArray;
    /// Gradient-norm convergence tolerance for Newton iteration.
    double mTolerance = 1e-10;
    /// Maximum Newton iterations per seed.
    std::int64_t mMaxIterations = 30;
    /// Spatial deduplication tolerance across shared cell faces/seeds.
    double mMergeTolerance = 1e-6;
    /// Reference-lattice seed subdivisions per axis (0 = auto `order + 1`).
    std::int64_t mSeedSubdivisions = 0;
};

/// Result of `critical_points`.
struct CriticalPointsResult {
    /// `vertex` mesh of deduplicated critical points.
    Mesh mMesh;
    /// Number of local minima found.
    std::int64_t mNumMinima = 0;
    /// Number of local maxima found.
    std::int64_t mNumMaxima = 0;
    /// Number of saddle points found.
    std::int64_t mNumSaddles = 0;
    /// Number of degenerate critical points found.
    std::int64_t mNumDegenerate = 0;
};

/**
 * @brief Locate and classify critical points of a scalar point_data field.
 */
MESHIOPLUSPLUS_API CriticalPointsResult critical_points(const Mesh& rMesh,
                                                        const CriticalPointsOptions& rOptions);

}  // namespace meshioplusplus
