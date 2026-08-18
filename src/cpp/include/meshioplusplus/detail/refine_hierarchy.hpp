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
 * @file detail/refine_hierarchy.hpp
 * @brief Shared reader for the `refine:cell_id`/`refine:parent_id` persistent
 * hierarchy (`RefineOptions::mRecordHierarchy`, `refine.hpp`), hoisted verbatim
 * out of `operations/refine.cpp` -- `refine_attach_hierarchy` (deciding whether
 * to maintain an existing hierarchy or start a fresh one) and `undo_green`
 * (resolving the coarse mesh's own id space) both need the identical
 * Absent/Valid/Invalid read, and a second transcription would drift silently.
 *
 * Free function in `meshioplusplus::detail`, called once per operation (not per
 * element), so its body lives in `src/cpp/src/detail/refine_hierarchy.cpp`
 * rather than inline here. Built on the uniform mesh API only.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// How an existing `refine:cell_id`/`refine:parent_id` pair on a mesh relates to
/// what a caller should do with it: `Absent` (no array at all), `Valid` (usable,
/// hence maintained) or `Invalid` (malformed or non-unique, hence
/// dropped-with-a-warning and treated like `Absent`).
enum class RefineHierarchyState { Absent, Valid, Invalid };

/**
 * @brief Read a mesh's `refine:cell_id`/`refine:parent_id`, if any.
 *
 * On success `rIds` holds one id per global (block-major) cell and `rIdBase` is
 * one past the largest id in use anywhere -- over BOTH arrays, since a cell's id
 * can outlive its own row once the cell is split (the row is gone, but the id
 * must never be reissued to a different cell).
 *
 * Uniqueness is the guard here, not staleness: these ids are *values*, not
 * indices, so `reorder`/`crop`/`clean` carry them correctly with no coordinate
 * check at all (unlike `refine:entity`). A repeated id means the mesh was
 * `merge`d with another hierarchy, or a cell-splitting operation replicated the
 * array without updating it -- either way the array is this operation's own
 * bookkeeping, not user input, so it is warned-and-dropped rather than rejected
 * outright.
 * @param rMesh the mesh to read.
 * @param rBases a `block_bases` table for @p rMesh.
 * @param rIds out: one id per global cell (only meaningful on `Valid`).
 * @param rIdBase out: one past the largest id in use (only meaningful on `Valid`).
 * @return the hierarchy's state.
 */
MESHIOPLUSPLUS_API RefineHierarchyState
refine_read_hierarchy(const Mesh& rMesh, const std::vector<std::int64_t>& rBases,
                      std::vector<std::int64_t>& rIds, std::int64_t& rIdBase);

}  // namespace detail
}  // namespace meshioplusplus
