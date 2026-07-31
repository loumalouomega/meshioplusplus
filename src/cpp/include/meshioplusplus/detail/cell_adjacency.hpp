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
 * @file detail/cell_adjacency.hpp
 * @brief The single owner of **cell-to-cell adjacency by shared node**: the
 * cell -> node and node -> cell incidences, and the neighbour list derived from
 * them.
 *
 * This is a different object from `detail/node_adjacency.hpp`'s
 * `NodeAdjacency`, which is a node-to-node graph. Two operations need this one
 * and must not disagree about what "neighbour" means:
 *  - `partition` grows its ghost (halo) layers breadth-first over it;
 *  - `gradient`'s least-squares stencil fits the linear model over it.
 *
 * Rows are indexed by the **global (block-major) cell index** owned by
 * `detail/cell_index.hpp` — never re-derive that numbering here.
 *
 * The per-cell node walk goes through `detail::cell_node_ids`, which is
 * dtype-agnostic (`detail::read_int`) and drops out-of-range ids. That matters:
 * connectivity is *not* guaranteed to be Int64 (a MESHIO-backed mesh often
 * carries Int32 straight from numpy), so reading it as `As<std::int64_t>()`
 * — which performs no dtype check — silently fuses two Int32 entries into one
 * garbage id. `partition`'s ghost layers did exactly that before this header
 * existed, producing halos that were quietly too small rather than failing.
 *
 * Free functions in `meshioplusplus::detail`, each called once per operation
 * (not per element), so their bodies live in
 * `src/cpp/src/detail/cell_adjacency.cpp` rather than inline here — the
 * convention `detail/node_adjacency.hpp` and `detail/subset.hpp` follow.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief A ragged incidence relation in compressed-sparse-row form.
 *
 * Row `i` is `mEntries[mOffsets[i] .. mOffsets[i + 1])`; `mOffsets` has one
 * more entry than there are rows and starts at 0.
 */
struct CellIncidence {
    std::vector<std::int64_t> mOffsets;  ///< Row offsets, size `rows + 1`.
    std::vector<std::int64_t> mEntries;  ///< Payload, size `mOffsets.back()`.
};

/**
 * @brief Global (block-major) cell -> the node ids it references.
 *
 * Rows are in connectivity order (deliberately **not** sorted, so the relation
 * stays a faithful record of the cell's own node list); a polyhedron
 * contributes the concatenation of its face rings, which may repeat a node.
 * Out-of-range ids are dropped by `cell_node_ids`.
 * @param rMesh Mesh to read connectivity from (never modified).
 * @param NumPoints Node-id upper bound; ids outside `[0, NumPoints)` are dropped.
 * @return The incidence, with one row per global cell.
 */
MESHIOPLUSPLUS_API CellIncidence build_cell_node_incidence(const Mesh& rMesh, std::size_t NumPoints);

/**
 * @brief Inverts a cell -> node incidence into node -> cell, by counting.
 *
 * Each row lists the cells referencing that node, ascending by global cell
 * index (a consequence of the counting pass walking cells in order), with a
 * cell repeated as many times as it references the node.
 * @param rCellNodes The forward incidence from `build_cell_node_incidence`.
 * @param NumPoints Number of nodes; the result has `NumPoints` rows.
 * @param NumCells Number of global cells (the forward incidence's row count).
 * @return The inverted incidence.
 */
MESHIOPLUSPLUS_API CellIncidence invert_cell_node_incidence(const CellIncidence& rCellNodes,
                                                            std::size_t NumPoints,
                                                            std::size_t NumCells);

/**
 * @brief The cells sharing at least one node with @p Cell.
 *
 * Ascending global cell index, de-duplicated, with @p Cell itself excluded.
 * This is the definition of "neighbour" that `gradient`'s least-squares stencil
 * and `partition`'s ghost growth share; the fixed ascending order is what makes
 * a neighbour-summed quantity byte-identical across mesh backends, thread
 * counts and the C++/numpy boundary.
 * @param rCellNodes Forward incidence from `build_cell_node_incidence`.
 * @param rNodeCells Inverted incidence from `invert_cell_node_incidence`.
 * @param Cell Global cell index to query.
 * @param rOut Cleared and filled with the neighbour ids.
 */
MESHIOPLUSPLUS_API void cell_node_neighbors(const CellIncidence& rCellNodes,
                                            const CellIncidence& rNodeCells, std::size_t Cell,
                                            std::vector<std::int64_t>& rOut);

}  // namespace detail
}  // namespace meshioplusplus
