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
 * @file operations/adjacency.hpp
 * @brief A mesh's node or cell graph in compressed-sparse-row form, and the
 * METIS and Scotch graph files written from it.
 *
 * `adjacency(mesh, kind)` returns `(indptr, indices)`: the neighbours of vertex
 * `i` are `indices[indptr[i] .. indptr[i + 1])`, sorted ascending, with no
 * duplicates and no self-loops, and every undirected edge stored in both
 * directions, so `indptr.back()` is twice the edge count.
 *
 * ### Kinds
 *
 * - `AdjacencyKind::Node` (`"node"`) -- nodes joined by a cell edge, the graph
 *   `edge_index` returns in coordinate form. Edges come from the linear types'
 *   edge tables, the rings of polygon rows (ragged or uniform `polygon*`) and
 *   the face rings of polyhedra; a higher-order, VTK-Lagrange or `custom`
 *   block, whose edges are not tabulated, contributes none (its nodes stay
 *   isolated unless another block joins them).
 * - `AdjacencyKind::Clique` (`"clique"`) -- nodes that share any cell: the
 *   sparsity pattern of an assembled finite-element matrix (without the
 *   diagonal), what a preallocation or a bandwidth reordering needs. Every
 *   cell type contributes.
 * - `AdjacencyKind::Cell` (`"cell"`) -- cells that share a facet: a face when
 *   the mesh holds any 3-D cell, else an edge. Vertices are global
 *   (block-major) cell indices, and lower-dimensional cells are isolated. This
 *   is the dual graph `partition`'s KaHIP method cuts, and the graph METIS,
 *   Scotch and KaHIP take; polyhedra and polygons take part.
 *
 * ### Guarantees
 *
 * Built through the uniform mesh API, so every mesh backend gives the same
 * graph; the fill passes are serial and the per-row sorts independent, so the
 * graph is byte-identical on every parallel backend and thread count.
 *
 * ### Graph files
 *
 * `write_graph` writes a graph in METIS's format (`.graph`: a header `n m`,
 * then one line per vertex of its 1-based neighbours) or Scotch's (`.grf`:
 * version `0`, then `n 2m`, then `0 000`, then one line per vertex of its
 * degree and 0-based neighbours). Both formats describe an undirected graph,
 * so the writer refuses an asymmetric one, a self-loop and an out-of-range
 * neighbour.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Which graph `adjacency` builds.
enum class AdjacencyKind : std::int32_t {
    Node = 0,    ///< Nodes joined by a cell edge.
    Clique = 1,  ///< Nodes that share a cell (FEM sparsity).
    Cell = 2,    ///< Cells that share a facet (the dual graph).
};

/// A graph in compressed-sparse-row form.
struct Adjacency {
    /// Int64 `(n + 1,)`: row offsets; `mIndptr[0] == 0`.
    NDArray mIndptr;
    /// Int64 `(mIndptr[n],)`: the neighbours of each row, ascending.
    NDArray mIndices;
};

/// The graph-file dialect `write_graph` writes.
enum class GraphFormat : std::int32_t {
    Metis = 0,   ///< METIS / ParMETIS / KaHIP `.graph` (1-based).
    Scotch = 1,  ///< Scotch `.grf` (0-based).
};

/// The size of a graph `write_graph` wrote.
struct GraphSummary {
    std::int64_t mNumVertices = 0;  ///< Rows of the graph.
    std::int64_t mNumEdges = 0;     ///< Undirected edges (`indptr.back() / 2`).
    std::int64_t mMaxDegree = 0;    ///< The largest row length.
};

/** @brief `"node"`, `"clique"` or `"cell"`; anything else throws naming them. */
MESHIOPLUSPLUS_API AdjacencyKind adjacency_kind_from_name(const std::string& rName);

/** @brief The name `adjacency_kind_from_name` accepts for @p Kind. */
MESHIOPLUSPLUS_API std::string adjacency_kind_name(AdjacencyKind Kind);

/** @brief `"metis"` or `"scotch"`; anything else throws naming them. */
MESHIOPLUSPLUS_API GraphFormat graph_format_from_name(const std::string& rName);

/**
 * @brief The dialect a path's extension names: `.graph` or `.metis` for METIS,
 * `.grf` for Scotch.
 * @throws std::invalid_argument for any other extension, naming these.
 */
MESHIOPLUSPLUS_API GraphFormat graph_format_for_path(const std::string& rPath);

/**
 * @brief The @p Kind graph of @p rMesh. A connectivity entry outside the
 * points is ignored, as the shared adjacency kernels ignore it.
 */
MESHIOPLUSPLUS_API Adjacency adjacency(const Mesh& rMesh, AdjacencyKind Kind = AdjacencyKind::Node);

/** @brief The vertex count, undirected edge count and largest degree of @p rGraph. */
MESHIOPLUSPLUS_API GraphSummary graph_summary(const Adjacency& rGraph);

/**
 * @brief Writes @p rGraph to @p rPath as a METIS or Scotch graph file.
 * @throws std::invalid_argument on a malformed CSR (decreasing offsets, an
 *         out-of-range neighbour), a self-loop or an asymmetric graph;
 *         std::runtime_error when the file cannot be written.
 */
MESHIOPLUSPLUS_API GraphSummary write_graph(const std::string& rPath, const Adjacency& rGraph,
                                            GraphFormat Format);

}  // namespace meshioplusplus
