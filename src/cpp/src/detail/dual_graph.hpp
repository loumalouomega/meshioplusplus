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
 * @file detail/dual_graph.hpp
 * @brief The facet-sharing dual graph of a mesh in CSR form: cells are
 * vertices, and two cells are joined when they share a facet -- a face when the
 * mesh holds any 3-D cell, else an edge.
 *
 * A **core-private** header (the `cell_components.hpp` precedent): no installed
 * header names it, and it adds nothing to the ABI. `partition`'s KaHIP method
 * and the `adjacency` operation's `"cell"` kind share it, so the graph KaHIP
 * cuts is the graph `adjacency` returns.
 *
 * Vertices are global (block-major) cell indices. Each undirected edge is
 * stored twice, each row is sorted ascending and free of duplicates and
 * self-loops. Facet keys are built in parallel into disjoint slots and grouped
 * by sorted runs; the adjacency fill is serial, so the graph is byte-identical
 * on every backend and thread count.
 *
 * Blocks whose dimension is below the dual dimension, and blocks without a
 * facet table, are isolated vertices. Polyhedra contribute their own faces.
 * Polygon blocks (ragged rows, or uniform `polygon*` rows) are isolated unless
 * `PolygonRings` is set, when each row contributes its ring of edges in a 2-D
 * dual graph: `partition` leaves it off, so its KaHIP input is unchanged.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// A cell-to-cell graph in compressed-sparse-row form.
struct DualGraphCsr {
    std::vector<std::int64_t> mXadj;    ///< Row offsets, size `NumCells + 1`.
    std::vector<std::int64_t> mAdjncy;  ///< Neighbour cells, size `mXadj.back()`.
};

/**
 * @brief Builds the facet-sharing dual graph of @p rMesh.
 * @param rMesh The mesh (never modified).
 * @param NumCells The total number of cells over all blocks.
 * @param PolygonRings Whether polygon rows take part (as ring edges, 2-D only).
 */
DualGraphCsr build_facet_dual_graph(const Mesh& rMesh, std::size_t NumCells,
                                    bool PolygonRings = false);

}  // namespace detail
}  // namespace meshioplusplus
