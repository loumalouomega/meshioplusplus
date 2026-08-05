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
 * @file operations/crop.hpp
 * @brief Dependency-free subsetting: extract the part of a mesh inside an
 * axis-aligned bounding box, inside a half-space, or satisfying a comparison on
 * one of its own `cell_data` arrays.
 *
 * A point is "inside" the bbox when `lo <= p <= hi` component-wise, or inside
 * the half-space when `(p - point) . normal >= 0`. A cell is kept when ALL of
 * its nodes are inside (`CropMode::All`, the default) or when ANY node is inside
 * (`CropMode::Any`). The kept cells form a new mesh with unused points pruned
 * and connectivity + all data remapped. Optionally the original point/cell ids
 * are recorded as data arrays.
 *
 * ### The predicate crop, and why `CropMode` does not apply to it
 *
 * `crop_predicate` keeps cells whose value in a scalar `cell_data` array
 * satisfies a comparison. That is deliberately **general rather than
 * inside/outside-a-surface specific**: the inside/outside case composes as
 *
 * ```
 * crop_predicate(distance_to_surface(mesh, skin, {.mLocation = Center}).mMesh,
 *                kSdfDistanceName, RefineCompare::Less, 0.0)
 * ```
 *
 * and the same one mode also crops by `quality:*`, by a material id, by
 * `partition:part`, or by anything `data_calc` can produce. A dedicated
 * crop-by-surface would have served one of those.
 *
 * `CropMode::All|Any` is **not** a parameter here, and its absence is the honest
 * answer rather than an omission: bbox and half-space test *points* and then
 * need a rule for reducing a cell's several nodes to one verdict, whereas a
 * `cell_data` predicate is already one value per cell and has nothing to reduce.
 * A mode that meant nothing would be worse than no mode.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/refine.hpp"

namespace meshioplusplus {

/// Cell-keep policy for `crop`.
enum class CropMode {
    All,  ///< Keep a cell only if every node is inside the region.
    Any,  ///< Keep a cell if any node is inside the region.
};

/**
 * @brief The result of a crop: the pruned submesh plus the applied index maps
 * (so the Python shim can remap the shim-only sets).
 */
struct CropResult {
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> output index (-1 if pruned).
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> output (-1 if dropped).
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Crop a mesh to an axis-aligned bounding box.
 * @param rMesh the input mesh.
 * @param pLo the box lower corner (3 doubles; unused axes ignored for 2D).
 * @param pHi the box upper corner (3 doubles).
 * @param mode keep-all-nodes-inside (default) vs keep-any-node-inside.
 * @param record_ids attach Int64 `crop:original_point_id` / `crop:original_cell_id`.
 * @return the pruned submesh and index maps.
 */
MESHIOPLUSPLUS_API CropResult crop_bbox(const Mesh& rMesh, const double* pLo, const double* pHi,
                                        CropMode mode = CropMode::All, bool record_ids = false);

/**
 * @brief Crop a mesh to the half-space `(p - point) . normal >= 0`.
 * @param rMesh the input mesh.
 * @param pPoint a point on the plane (3 doubles).
 * @param pNormal the plane normal (3 doubles); the kept side is where the
 *        signed distance is non-negative.
 * @param mode keep-all-nodes-inside (default) vs keep-any-node-inside.
 * @param record_ids attach Int64 `crop:original_point_id` / `crop:original_cell_id`.
 * @return the pruned submesh and index maps.
 */
MESHIOPLUSPLUS_API CropResult crop_halfspace(const Mesh& rMesh, const double* pPoint,
                                             const double* pNormal, CropMode mode = CropMode::All,
                                             bool record_ids = false);

/**
 * @brief Crop a mesh to the cells whose value in @p rArray satisfies a comparison.
 * @param rMesh the input mesh.
 * @param rArray the name of a **scalar** `cell_data` array covering every block.
 * @param Op the comparison; the shared `RefineCompare` vocabulary, evaluated by
 *        the shared `refine_compare_value` so the two operations cannot drift.
 * @param Value the right-hand side. **A non-finite cell value never matches**,
 *        whatever the comparison -- see `refine_compare_value`.
 * @param record_ids attach Int64 `crop:original_point_id` / `crop:original_cell_id`.
 * @return the pruned submesh and index maps, exactly as the other two crops.
 * @throws std::invalid_argument when @p rArray is not a `cell_data` array (the
 *         message lists what is), does not cover every block, has the wrong row
 *         count, or is not scalar. `point_data` is refused by name rather than
 *         averaged onto cells: `data to-cell` is the explicit way to do that, and
 *         doing it implicitly would make the kept set depend on an averaging rule
 *         the caller never asked for.
 */
MESHIOPLUSPLUS_API CropResult crop_predicate(const Mesh& rMesh, const std::string& rArray,
                                             RefineCompare Op, double Value,
                                             bool record_ids = false);

}  // namespace meshioplusplus
