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
 * @file cell_subdivision.hpp
 * @brief The **shared** per-cell-type edge and quad-face tables used whenever a
 * cell has to gain nodes at the midpoints of its own sub-entities: the
 * `elevate` mode of `operations/convert_cells.hpp` (one node per edge) and the
 * uniform subdivision of `operations/refine.hpp` (one node per edge, per quad
 * face, and for the hexahedron one per body).
 *
 * Unlike `cell_edges.hpp`, which deliberately covers **surface** cells only
 * (it backs the 2D branch of the boundary extractor and is gated on
 * `cell_type_dimension(ct) == 2`), this header covers every type that can gain
 * mid-entity nodes, volume types included. The 2D rows are not duplicated here
 * -- `cell_refine_edges` projects `cell_edges()` for `Triangle`/`Quad`, so that
 * table stays the single owner of the 2D edge order.
 *
 * **Row order is a contract, not a convenience.** Each table's row `k` maps to
 * a fixed slot in the corresponding meshio/VTK higher-order node layout, and
 * both consumers rely on that:
 *  - `cell_refine_edges(T)` row `k` -> node `num_corners(T) + k` of the
 *    serendipity quadratic type (`tetra10` nodes 4-9 are edges 01,12,02,03,13,23;
 *    `hexahedron20` nodes 8-19; `wedge15` nodes 6-14; `pyramid13` nodes 5-12).
 *  - `cell_refine_quad_faces(T)` row `k` -> the face-centre node of the full
 *    Lagrange type (`hexahedron27` nodes 20-25, `wedge18` nodes 15-17, `quad9`
 *    node 8).
 *
 * Face rows here are wound in that VTK centre-node order, which is *not* the
 * outward winding of `cell_faces.hpp`. That is harmless for both consumers: a
 * face contributes only a sorted node key and the mean of its corners, and both
 * are winding-independent. `tests/cpp/test_cell_subdivision.cpp` cross-checks
 * every row against `cell_faces()`'s own centre-node columns so the two tables
 * cannot drift.
 *
 * These tables are looked up once per cell **block** (not per cell or per node),
 * so the lookups live in `src/cpp/src/detail/cell_subdivision.cpp` rather than
 * inline here.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

/// The two local corner indices of one edge of a cell.
using CellEdgePair = std::array<std::uint8_t, 2>;

/// The four local corner indices of one quadrilateral face of a cell.
using CellQuadFace = std::array<std::uint8_t, 4>;

/**
 * @brief The edges of a cell type, ordered to match its quadratic mid-node
 * slots, or an empty list for types with no mid-edge layout.
 * @param Type The cell type to query.
 * @return Reference to the process-wide edge table (empty if unsupported).
 */
MESHIOPLUSPLUS_API const std::vector<CellEdgePair>& cell_refine_edges(CellType Type);

/**
 * @brief The quadrilateral faces of a cell type, ordered to match its
 * full-Lagrange face-centre node slots, or an empty list for types with no quad
 * face (`line`, `triangle`, `tetra`).
 *
 * A 2D `Quad` reports its own single face, so that a mesh mixing a
 * `hexahedron` block with a boundary `quad` block subdivides conformingly
 * across the dimension boundary.
 *
 * @param Type The cell type to query.
 * @return Reference to the process-wide quad-face table (empty if none).
 */
MESHIOPLUSPLUS_API const std::vector<CellQuadFace>& cell_refine_quad_faces(CellType Type);

}  // namespace detail
}  // namespace meshioplusplus
