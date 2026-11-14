//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████    ███████      ░░░         ░░░
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
 * @file ghost_names.hpp
 * @brief The cell/point/field data names for ghost cells, owning parts and
 * global ids, shared by the format readers and writers and the core operations.
 * See doc/ghosts.md for the convention.
 *
 * A **format-private** header (the `formats/gid_common.hpp` precedent): no
 * installed header may name it.
 */

namespace meshioplusplus {
namespace detail {

/// cell_data, Int64: 0 = owned by its piece, L >= 1 = ghost at BFS layer depth L.
inline constexpr const char* kPartitionGhost = "partition:ghost";

/// cell_data / point_data, Int64: owning part (0-based) of a ghost cell
/// (DUPLICATEPOINT owner for points).
inline constexpr const char* kPartitionOwner = "partition:owner";

/// cell_data, Int64: the piece a cell belongs to (pvtu/pvtp/pvd/vtx read; merge output).
inline constexpr const char* kPartitionPart = "partition:part";

/// point_data / cell_data, Int64: 0-based input indices (partition, split).
inline constexpr const char* kPartitionOriginalPointId = "partition:original_point_id";
inline constexpr const char* kPartitionOriginalCellId = "partition:original_cell_id";

/// point_data / cell_data, Int64: vendor/global ids (radioss, dpf, refine, cgns,
/// exodus maps when non-identity).
inline constexpr const char* kGlobalNodeId = "global:node_id";
inline constexpr const char* kGlobalCellId = "global:cell_id";

/// point_data / cell_data, UInt8: VTK-native bitfield (DUPLICATECELL = 1 cells,
/// DUPLICATEPOINT = 1 points); kept verbatim by vtu/vtp/vtx/pvtu.
inline constexpr const char* kVtkGhostType = "vtkGhostType";

/// VTK-native global ids carried by vtu/vtp/vtx/pvtu.
inline constexpr const char* kVtkOriginalPointIds = "vtkOriginalPointIds";
inline constexpr const char* kVtkOriginalCellIds = "vtkOriginalCellIds";

/// cell_data, Int64: Gmsh-native 1/0 flag from `$GhostElements`.
inline constexpr const char* kGmshGhost = "gmsh:ghost";

/// VTK bit values of `vtkGhostType`.
inline constexpr unsigned char kVtkDuplicateCell = 1;
inline constexpr unsigned char kVtkDuplicatePoint = 1;

}  // namespace detail
}  // namespace meshioplusplus
