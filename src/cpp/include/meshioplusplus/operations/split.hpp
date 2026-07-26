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
 * @file operations/split.hpp
 * @brief Dependency-free partitioning: split one mesh into several submeshes by
 * cell type, connected component, or an integer cell_data tag.
 *
 *  - **Type**: one piece per (first-seen) cell type; key = the meshio type name.
 *  - **Component**: flood-fill over node-sharing cell adjacency (union-find);
 *    one piece per connected component; keys `"0"`, `"1"`, ...
 *  - **Tag**: one piece per distinct value of a named integer `cell_data` array
 *    (the first integer cell_data when the name is empty); key = the value.
 *
 * Each piece is pruned to only its own points, with connectivity and all data
 * remapped, and the per-piece point/cell index maps are returned (so the Python
 * shim can remap the shim-only sets). File-writing (the `out_{key}.vtu` pattern)
 * is a caller concern — the core just returns the pieces.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The partitioning criterion for `split`.
enum class SplitBy {
    Type,       ///< One piece per cell type.
    Component,  ///< One piece per connected component (node-sharing).
    Tag,        ///< One piece per distinct integer cell_data value.
    /**
     * One piece per named `Cell` region (`"regions"`, plural -- distinct from
     * the pre-existing singular `"region"` name, which remains an alias for
     * `Tag` for backward compatibility; see `split_by_from_name`). Unlike
     * every other criterion this is not a partition: regions may overlap, so
     * a cell can appear in zero, one, or several pieces.
     */
    Region,
};

/// One output piece of a `split`.
struct SplitPiece {
    std::string mKey;  ///< The piece's key (type name / component index / tag value).
    Mesh mMesh;        ///< The pruned submesh.
    /// Int64 shape `(num_points_in,)`, input point index -> this piece's point index (-1 if
    /// absent).
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> this piece's
    /// in-block index (-1 if the cell is not in this piece).
    std::vector<NDArray> mCellMaps;
};

/// The result of `split`: the ordered pieces.
struct SplitResult {
    std::vector<SplitPiece> mPieces;
};

/**
 * @brief Parse a criterion name into a `SplitBy`.
 * @param rName one of `"type"`, `"component"`, `"region"`/`"tag"`, or
 *        `"regions"` (plural -- named `Cell` regions, `SplitBy::Region`)
 * @throws std::invalid_argument on an unknown name
 */
MESHIOPLUSPLUS_API SplitBy split_by_from_name(const std::string& rName);

/**
 * @brief Split a mesh into pieces by the given criterion.
 * @param rMesh the mesh to split.
 * @param by the partitioning criterion.
 * @param rTagName for `SplitBy::Tag`, the integer `cell_data` name to split on
 *        (empty = auto-detect the first integer cell_data); ignored otherwise.
 * @return the ordered pieces with their keys and index maps.
 * @throws std::invalid_argument for `SplitBy::Tag` when no suitable cell_data
 *         tag is found.
 */
MESHIOPLUSPLUS_API SplitResult split(const Mesh& rMesh, SplitBy by, const std::string& rTagName = "");

}  // namespace meshioplusplus
