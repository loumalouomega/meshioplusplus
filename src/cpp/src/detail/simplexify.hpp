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
 * @file detail/simplexify.hpp
 * @brief `convert_cells(Simplexify)` with the exact-tetrahedron recognition
 * `poly_to_tet` adds.
 *
 * A **core-private** header (the `detail/slot_runs.hpp` precedent): it lives
 * beside the `.cpp` files, no installed header names it, and it adds nothing to
 * the API or the ABI -- `ConvertCellsOptions` keeps its layout, and
 * `convert_cells` keeps calling the decomposition with `ExactTets = false`, so
 * its output is unchanged.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"

namespace meshioplusplus {
namespace detail {

/// How `simplexify_ex` treated the polyhedra it met.
struct SimplexifyCounts {
    std::int64_t mNumExactTets = 0;    ///< polyhedra emitted as exactly one tetra
    std::int64_t mNumFannedCells = 0;  ///< polyhedra decomposed by the fan
};

/**
 * @brief Simplexify @p rMesh (regions carried), optionally emitting a
 * polyhedron that already is a tetrahedron as that one `tetra`.
 * @param pCounts when not null, incremented per polyhedron.
 */
ConvertCellsResult simplexify_ex(const Mesh& rMesh, bool RecordParentIds, bool ExactTets,
                                 SimplexifyCounts* pCounts);

}  // namespace detail
}  // namespace meshioplusplus
