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
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus::detail {

/// The orientation-reversing node permutation of @p Type, or `nullptr`.
const std::vector<int>* cell_flip_permutation(CellType Type);

/// Per facet of @p Type, its number after the flip; empty for a type without
/// facets (a vertex or a line), `nullptr` when the type has no permutation.
const std::vector<int>* cell_flip_facet_map(CellType Type);

/// How one block's cells are flipped: a fixed table, a polygon ring or a
/// polyhedron's face rings (`mirror` and `flip_normals`, roadmap 7.8.2).
struct BlockFlip {
    const std::vector<int>* pPerm = nullptr;
    const std::vector<int>* pFacets = nullptr;
    std::vector<int> mPolygonPerm;  ///< a rectangular polygon block's own table
    bool mRaggedPolygon = false;
    bool mPolyhedron = false;

    /// The node permutation of a rectangular block.
    const std::vector<int>& Perm() const { return mPolygonPerm.empty() ? *pPerm : mPolygonPerm; }
};

/// The flip of a block of cells.
/// @throws Unsupported for a type with no orientation flip, naming @p pOp.
template <class CellViewT>
BlockFlip block_flip(const CellViewT& rCells, const char* pOp) {
    BlockFlip f;
    const std::string type(rCells.Type());
    const CellType ct = cell_type_from_name(type);
    if (rCells.IsPolyhedron()) {
        f.mPolyhedron = true;
        return f;
    }
    if (ct == CellType::Polygon) {
        if (rCells.IsRagged()) {
            f.mRaggedPolygon = true;
            return f;
        }
        // A polygon with rectangular storage: every row has the same size.
        const int npc = static_cast<int>(rCells.NodesPerCell());
        f.mPolygonPerm.push_back(0);
        for (int i = npc - 1; i >= 1; --i)
            f.mPolygonPerm.push_back(i);
        return f;
    }
    f.pPerm = rCells.IsRagged() ? nullptr : cell_flip_permutation(ct);
    f.pFacets = cell_flip_facet_map(ct);
    if (f.pPerm == nullptr ||
        static_cast<std::size_t>(cell_type_num_nodes(ct)) != rCells.NodesPerCell())
        throw Unsupported(std::string("meshio++: ") + pOp +
                          ": no orientation flip for cell type '" + type +
                          "' (linear and quadratic types, polygons and polyhedra only)");
    return f;
}

/// The local facet a flipped cell gives facet @p F of a cell with @p N facets
/// (a polygon's N is its node count), or -1 when F is not a facet of the cell.
inline std::int64_t block_flip_facet(const BlockFlip& rFlip, std::int64_t F, std::int64_t N) {
    if (F < 0)
        return -1;
    if (rFlip.mPolyhedron)
        return F < N ? F : -1;
    if (rFlip.mRaggedPolygon || !rFlip.mPolygonPerm.empty())
        return F < N ? N - 1 - F : -1;
    if (rFlip.pFacets == nullptr || F >= static_cast<std::int64_t>(rFlip.pFacets->size()))
        return -1;
    return (*rFlip.pFacets)[static_cast<std::size_t>(F)];
}

}  // namespace meshioplusplus::detail
