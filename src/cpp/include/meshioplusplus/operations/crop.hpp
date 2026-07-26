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
 * @brief Dependency-free spatial subsetting: extract the part of a mesh inside
 * an axis-aligned bounding box or a half-space.
 *
 * A point is "inside" the bbox when `lo <= p <= hi` component-wise, or inside
 * the half-space when `(p - point) . normal >= 0`. A cell is kept when ALL of
 * its nodes are inside (`CropMode::All`, the default) or when ANY node is inside
 * (`CropMode::Any`). The kept cells form a new mesh with unused points pruned
 * and connectivity + all data remapped. Optionally the original point/cell ids
 * are recorded as data arrays.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

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
MESHIOPLUSPLUS_API CropResult crop_halfspace(const Mesh& rMesh, const double* pPoint, const double* pNormal,
                          CropMode mode = CropMode::All, bool record_ids = false);

}  // namespace meshioplusplus
