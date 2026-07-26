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
 * @file operations/stats.hpp
 * @brief Dependency-free geometric statistics of a mesh (read-only), the
 * geometric complement to the topological `info`.
 *
 * `compute_stats` reports the bounding box + extents, centroid, per-cell-type
 * counts, total surface area (of 2D cells plus the boundary of 3D cells), total
 * signed & unsigned volume of 3D cells, and the count of inverted (negative
 * signed-volume) cells. The mesh is not modified.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Geometric statistics of a mesh (see `compute_stats`).
struct StatsReport {
    std::int64_t mNumPoints = 0;  ///< Number of points.
    std::int64_t mNumCells = 0;   ///< Total number of cells (all blocks).
    /// Axis-aligned bounding box lower / upper corner and extent (max-min). For a
    /// 2D mesh the z components are 0.
    double mBBoxMin[3] = {0, 0, 0};
    double mBBoxMax[3] = {0, 0, 0};
    double mExtent[3] = {0, 0, 0};
    double mCentroid[3] = {0, 0, 0};  ///< Mean of the point coordinates.
    /// Per cell type (first-seen order): the type name and its cell count.
    std::vector<std::pair<std::string, std::int64_t>> mCellTypeCounts;
    double mTotalArea = 0.0;        ///< Area of 2D cells + boundary of 3D cells.
    double mSignedVolume = 0.0;     ///< Sum of signed volumes of 3D cells.
    double mUnsignedVolume = 0.0;   ///< Sum of |volume| of 3D cells.
    std::int64_t mNumInverted = 0;  ///< 3D cells with negative signed volume.
};

/**
 * @brief Compute geometric statistics of a mesh (read-only).
 * @param rMesh the mesh to measure (unmodified).
 * @return the populated `StatsReport`.
 */
MESHIOPLUSPLUS_API StatsReport compute_stats(const Mesh& rMesh);

}  // namespace meshioplusplus
