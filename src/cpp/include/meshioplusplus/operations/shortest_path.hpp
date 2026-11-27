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
 * @file operations/shortest_path.hpp
 * @brief Shortest geodesic paths (`shortest_path`) and multi-source geodesic
 * Voronoi partitioning (`geodesic_voronoi`) on meshes (roadmap 9.14).
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

inline constexpr const char* kPathDistanceName = "path:distance";
inline constexpr const char* kPathSourceCellName = "path:source_cell";
inline constexpr const char* kGeodesicNearestSourceName = "geodesic:nearest_source";

/// Algorithm for `shortest_path`.
enum class ShortestPathMethod : std::int32_t {
    /// Deterministic edge-graph Dijkstra with ascending vertex-id tie-breaking.
    Dijkstra = 0,
    /// Gradient back-trace across triangle faces from `geodesic_distance`'s heat field.
    HeatTrace = 1,
};

/// Parse `"dijkstra"` or `"heat_trace"`.
MESHIOPLUSPLUS_API ShortestPathMethod shortest_path_method_from_name(const std::string& rName);

/// Options for `shortest_path`.
struct ShortestPathOptions {
    /// Source vertex id (>= 0) when `mHasSourcePoint` is false.
    std::int64_t mSource = -1;
    /// Target vertex id (>= 0) when `mHasTargetPoint` is false.
    std::int64_t mTarget = -1;
    /// If true, snap `mSourcePoint` to the nearest mesh vertex.
    bool mHasSourcePoint = false;
    /// 3-D source coordinates when `mHasSourcePoint` is true.
    std::array<double, 3> mSourcePoint = {0.0, 0.0, 0.0};
    /// If true, snap `mTargetPoint` to the nearest mesh vertex.
    bool mHasTargetPoint = false;
    /// 3-D target coordinates when `mHasTargetPoint` is true.
    std::array<double, 3> mTargetPoint = {0.0, 0.0, 0.0};
    /// Path extraction method (`Dijkstra` or `HeatTrace`).
    ShortestPathMethod mMethod = ShortestPathMethod::Dijkstra;
    /// Heat-method time factor when `mMethod == HeatTrace`.
    double mTimeFactor = 1.0;
};

/// Result of `shortest_path`.
struct ShortestPathResult {
    /// 1-D `line` mesh carrying `point_data["path:distance"]` and `cell_data["path:source_cell"]`.
    Mesh mMesh;
    /// Total arc length of the extracted path (0 if source == target, NaN if unreachable).
    double mTotalLength = 0.0;
    /// Number of line segments in `mMesh`.
    std::int64_t mNumSegments = 0;
    /// Resolved source vertex id.
    std::int64_t mSourceVertex = -1;
    /// Resolved target vertex id.
    std::int64_t mTargetVertex = -1;
    /// True if a path connecting source and target was found.
    bool mReached = false;
};

/// Options for `geodesic_voronoi`.
struct GeodesicVoronoiOptions {
    /// Seed vertex ids (0-based).
    std::vector<std::int64_t> mSources;
    /// Optional 3-D seed coordinates (flat `3 * K`), snapped to nearest mesh vertices.
    std::vector<double> mSeedPoints;
    /// If true, also attach majority-vote `cell_data["geodesic:nearest_source"]`.
    bool mAttachCellLabels = true;
};

/// Result of `geodesic_voronoi`.
struct GeodesicVoronoiResult {
    /// Copy of input mesh with `geodesic:distance` and `geodesic:nearest_source`.
    Mesh mMesh;
    /// Number of seed sources used.
    std::int64_t mNumSources = 0;
    /// Number of vertices not reachable from any seed.
    std::int64_t mNumUnreached = 0;
    /// Maximum finite distance across reached vertices.
    double mMaxDistance = 0.0;
};

/**
 * @brief Extract a shortest path between two landmarks on a mesh as a 1-D `line` mesh.
 */
MESHIOPLUSPLUS_API ShortestPathResult shortest_path(const Mesh& rMesh,
                                                    const ShortestPathOptions& rOptions);

/**
 * @brief Partition a mesh into multi-source geodesic Voronoi regions.
 */
MESHIOPLUSPLUS_API GeodesicVoronoiResult geodesic_voronoi(const Mesh& rMesh,
                                                          const GeodesicVoronoiOptions& rOptions);

}  // namespace meshioplusplus
