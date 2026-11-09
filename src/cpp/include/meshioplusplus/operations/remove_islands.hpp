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
 * @file operations/remove_islands.hpp
 * @brief Drop a scan's crumbs (roadmap 7.8.1): the connected components whose
 * bounding-box diagonal, area or cell count is below a threshold.
 *
 * Components are joined by shared edges (cells sharing a facet) or by shared
 * points, numbered in the order of each one's lowest cell. A component is
 * removed when its diagonal is below `mMinDiagonal`, its area (of its 2-D
 * cells) below `mMinArea`, or its cell count below `mMinCells`; with
 * `mRelative` the diagonal and the area are fractions of the whole mesh's
 * diagonal and area. Points no kept cell uses go with them. A dry run removes
 * nothing and returns the mesh with the Int64 cell data
 * `islands:component` and `islands:removed` (1 for a cell that would go).
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

enum class IslandConnectivity : std::int32_t {
    Edge = 0,   ///< cells sharing a facet (an edge of a surface)
    Point = 1,  ///< cells sharing a point
};

struct RemoveIslandsOptions {
    IslandConnectivity mConnectivity = IslandConnectivity::Edge;
    double mMinDiagonal = 0.01;  ///< below it a component goes (relative by default)
    double mMinArea = 0.0;       ///< likewise, 0 off
    std::int64_t mMinCells = 0;  ///< likewise, 0 off
    bool mRelative = true;       ///< diagonal and area as fractions of the mesh's
    bool mDryRun = false;
};

struct RemoveIslandsResult {
    Mesh mMesh;
    std::int64_t mNumComponents = 0;
    std::int64_t mNumRemoved = 0;  ///< components removed (or that would be)
    std::int64_t mNumCellsRemoved = 0;
    std::int64_t mNumPointsRemoved = 0;
    /// Per component: its cell count, diagonal and area, and whether it goes.
    std::vector<std::int64_t> mComponentCells;
    std::vector<double> mComponentDiagonal;
    std::vector<double> mComponentArea;
    std::vector<std::uint8_t> mComponentRemoved;
};

/// IslandConnectivity from "edge" or "point". @throws std::invalid_argument otherwise.
MESHIOPLUSPLUS_API IslandConnectivity island_connectivity_from_name(const std::string& rName);

/**
 * @brief Remove the small connected components of @p rMesh (see the file comment).
 * @throws std::invalid_argument on a negative or non-finite threshold.
 */
MESHIOPLUSPLUS_API RemoveIslandsResult remove_islands(const Mesh& rMesh,
                                                      const RemoveIslandsOptions& rOptions = {});

}  // namespace meshioplusplus
