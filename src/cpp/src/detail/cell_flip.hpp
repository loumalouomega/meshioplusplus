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
 * @file detail/cell_flip.hpp
 * @brief The node permutation that reverses a cell's orientation, and where it
 * sends the cell's facets (roadmap 3.2.5, `mirror`).
 *
 * A reflection turns every cell inside out: a tetrahedron's volume goes
 * negative, a triangle's normal points the other way. Permuting the nodes of
 * each cell by `cell_flip_permutation` puts it right side out again: the
 * permuted cell `new[i] = old[P[i]]` covers the same nodes with the opposite
 * orientation, its mid-edge, mid-face and centre nodes still in their slots.
 * Every permutation is an involution and keeps node 0 in place.
 *
 * The facets move with the nodes: `cell_flip_facet_map(Type)[f]` is the
 * number, in the flipped cell, of what was facet `f` (a face of a 3-D cell as
 * `cell_faces` numbers them, an edge of a 2-D one as `cell_edges` does), so a
 * `Side` region entry `(cell, f)` becomes `(cell, map[f])`. The maps are found
 * from the facet tables themselves, not written out.
 *
 * Covered: `vertex`, `line` to `line11`, `triangle`, `triangle6`, `triangle7`,
 * `quad`, `quad8`, `quad9`, `tetra`, `tetra10`, `pyramid`, `pyramid13`,
 * `pyramid14`, `wedge`, `wedge15`, `wedge18`, `hexahedron`, `hexahedron20` and
 * `hexahedron27`. A polygon reverses its ring after node 0 (edge `k` becomes
 * edge `n - 1 - k`) and a polyhedron each face ring (its faces keep their
 * numbers); neither has a fixed table, so the caller does those. Every other
 * type returns `nullptr`.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus::detail {

/// The orientation-reversing node permutation of @p Type, or `nullptr`.
const std::vector<int>* cell_flip_permutation(CellType Type);

/// Per facet of @p Type, its number after the flip; empty for a type without
/// facets (a vertex or a line), `nullptr` when the type has no permutation.
const std::vector<int>* cell_flip_facet_map(CellType Type);

}  // namespace meshioplusplus::detail
