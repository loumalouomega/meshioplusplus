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
 * @file operations/side_data.hpp
 * @brief Materialize side regions as lower-dimensional cell blocks and back (roadmap §8.4).
 *
 * `sides_to_cells` turns a `Side` region (or all `Side` regions) into lower-dimensional
 * cell blocks carrying `cell_data["meshio:side_of"]` (`Int64 (n, 2)` `(global_cell, local_facet)`),
 * transferring any `field_data["exodus:sset:<region>:<var>"]` values onto the new facet
 * blocks' `cell_data`. `cells_to_sides` is the inverse, recovering `RegionKind::Side`
 * regions and `exodus:sset:<region>:<var>` field arrays from `meshio:side_of` blocks.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

inline constexpr const char* kSideOfName = "meshio:side_of";
inline constexpr const char* kEdgeOfName = "meshio:edge_of";

/// Options for `sides_to_cells`.
struct SidesToCellsOptions {
    /// Name of a specific `Side` region to materialize (empty = all `Side` regions).
    std::string mRegion;
};

/// Result of `sides_to_cells`.
struct SidesToCellsResult {
    Mesh mMesh;
    std::int64_t mNumFacets = 0;
    std::vector<std::string> mWarnings;
};

/// Options for `cells_to_sides`.
struct CellsToSidesOptions {
    /// Optional region name override when a facet block has no matching `Cell` region.
    std::string mRegionName;
    /// If true, remove the materialized `meshio:side_of` facet blocks from the output mesh.
    bool mDropFacetBlocks = false;
};

/// Result of `cells_to_sides`.
struct CellsToSidesResult {
    Mesh mMesh;
    std::int64_t mNumSides = 0;
    std::vector<std::string> mWarnings;
};

MESHIOPLUSPLUS_API SidesToCellsResult sides_to_cells(const Mesh& rMesh,
                                                     const SidesToCellsOptions& rOptions = {});

MESHIOPLUSPLUS_API CellsToSidesResult cells_to_sides(const Mesh& rMesh,
                                                     const CellsToSidesOptions& rOptions = {});

}  // namespace meshioplusplus
