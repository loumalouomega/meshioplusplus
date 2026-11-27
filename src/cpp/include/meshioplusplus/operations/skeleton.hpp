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
 * @file operations/skeleton.hpp
 * @brief Curve-skeleton and centerline extraction from a closed triangle
 * surface via mean-curvature flow contraction and connectivity thinning
 * (roadmap 9.13; Au et al. 2008; Tagliasacchi et al. 2012).
 *
 * Returns a 1-D `line` mesh carrying:
 *  - `point_data["skeleton:radius"]`: local inscribed radius (distance to the
 *    original surface via `detail::DistanceQuery`),
 *  - `point_data["skeleton:degree"]`: graph vertex degree (1 = endpoint,
 *    2 = regular, >= 3 = branch),
 *  - `field_data["skeleton:surface_map"]`: `(n_input_points,)` Int64 mapping
 *    each input surface vertex to its skeleton node id.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

inline constexpr const char* kSkeletonRadiusName = "skeleton:radius";
inline constexpr const char* kSkeletonDegreeName = "skeleton:degree";
inline constexpr const char* kSkeletonSurfaceMapName = "skeleton:surface_map";

/// Options for `extract_skeleton`.
struct SkeletonOptions {
    /// Maximum mean-curvature flow contraction iterations.
    std::int64_t mMaxIterations = 10;
    /// Initial Laplacian contraction weight multiplier (scaled by mean edge^2).
    double mContractionWeight = 1.0;
    /// Medial-axis attraction weight multiplier.
    double mAttractionWeight = 1.0;
    /// Stop contracting once enclosed volume ratio `|V| / |V_0|` falls below this threshold.
    double mVolumeTolerance = 1e-4;
    /// Edge-length merge tolerance on the thinned 1-D graph (0 = auto `0.25 * mean_edge`).
    double mMergeTolerance = 0.0;
    /// Recenter skeleton vertices toward the mean of their mapped surface vertices' nearest axis.
    bool mRecenter = true;
};

/// Result of `extract_skeleton`.
struct SkeletonResult {
    /// 1-D `line` mesh representing the extracted curve skeleton.
    Mesh mMesh;
    /// Number of contraction iterations executed.
    std::int64_t mContractionIterations = 0;
    /// Final volume ratio `|V_final| / |V_0|`.
    double mVolumeRatio = 0.0;
    /// Number of skeleton nodes.
    std::int64_t mNumNodes = 0;
    /// Number of skeleton line segments.
    std::int64_t mNumSegments = 0;
    /// Number of branch nodes (degree >= 3).
    std::int64_t mNumBranches = 0;
    /// Number of endpoint nodes (degree == 1).
    std::int64_t mNumEndpoints = 0;
};

/**
 * @brief Extract a 1-D curve skeleton / centerline from a triangle surface mesh.
 */
MESHIOPLUSPLUS_API SkeletonResult extract_skeleton(const Mesh& rMesh,
                                                   const SkeletonOptions& rOptions = {});

}  // namespace meshioplusplus
