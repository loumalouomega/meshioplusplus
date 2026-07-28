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
 * @file refine_templates.hpp
 * @brief The per-cell-type subdivision templates `operations/refine.hpp` uses,
 * indexed by which of the cell's edges are **split**.
 *
 * This is the sibling of `detail/cell_subdivision.hpp`: that header owns the
 * *order* of a cell's edges and quad faces, this one owns what to emit once you
 * know which of them carry a new node. It exists separately from `refine.cpp`
 * because the tables have properties worth asserting directly
 * (`tests/cpp/test_refine_templates.cpp`) rather than only sampling through a
 * refined mesh, and because the Python reference implementation is pinned
 * against them across the `_core` boundary.
 *
 * ### The mask
 *
 * A cell's state is a bitmask over `cell_refine_edges(Type)`: bit `k` is set iff
 * edge `k` carries a new mid-edge node. Local node ids in a template address the
 * same flat per-cell space `refine.cpp` already builds:
 *
 *     [0, num_corners)             the parent's own corner nodes
 *     [num_corners, +num_edges)    one per edge, cell_refine_edges() order
 *     [.., +num_quad_faces)        one per quad face, cell_refine_quad_faces() order
 *     [.., +1)                     the body centre (hexahedron only)
 *
 * which coincides with each type's own full-Lagrange numbering (`line3`,
 * `triangle6`, `quad9`, `tetra10`, `wedge18`, `hexahedron27`). A template only
 * ever references a slot whose entity actually exists for its mask -- an edge
 * slot only when that bit is set, a face-centre slot only when all four of that
 * face's edges are set, the body only on the full mask. That correspondence is
 * a *derived* rule, not a tabulated one: see `refine.hpp`.
 *
 * ### Admissibility and promotion
 *
 * Not every mask has a same-type subdivision. The **admissible** masks per type
 * are chosen so that the set is closed under intersection and contains the full
 * mask -- a Moore family -- which is exactly what makes "the smallest admissible
 * superset of `m`" well defined. `refine_promote_mask` is that closure operator:
 * monotone and idempotent, so the mesh-wide fixed point it drives is unique and
 * independent of the order cells are visited in. Determinism across mesh
 * backends and thread counts follows from the algebra rather than from a
 * convention about traversal.
 *
 *  | type         | admissible masks                                   | children       |
 *  |--------------|----------------------------------------------------|----------------|
 *  | `line`       | all 2                                              | 1, 2           |
 *  | `triangle`   | all 8                                              | 1, 2, 3, 4     |
 *  | `quad`       | none, either opposite pair, all                    | 1, 2, 2, 4     |
 *  | `tetra`      | none, 6 singles, all 15 pairs, 4 face-triples, all  | 1, 2, 3/4, 4, 8|
 *  | `wedge`      | none, the 6 triangle edges, the 3 verticals, all   | 1, 4, 2, 8     |
 *  | `hexahedron` | unions of the 3 parallel edge classes (8 of them)  | 1, 2, 2, 2, 4, 4, 4, 8 |
 *
 * The quadrilateral row is forced, not chosen: a quadrangulation of an `n`-gon
 * satisfies `4Q = B + 2I`, so an odd boundary count is impossible and a quad
 * with one or three split edges has *no* all-quad subdivision at any number of
 * interior nodes. The tetrahedron's three-edges-at-a-common-vertex mask is
 * excluded deliberately: it would put the ambiguous two-edge case on all three
 * incident faces at once.
 *
 * ### The one choice made from global ids
 *
 * A face with two *adjacent* split edges has a remnant quadrilateral needing a
 * diagonal, and the cell on the other side of that face must pick the same one.
 * The template therefore stores **both** variants plus the two local corner ids
 * `mTieA`/`mTieB` that decide between them: use `mChildrenAlt` when the global
 * node id at `mTieB` is smaller than the one at `mTieA`. Since the rule reduces
 * to "the diagonal starting at the smaller-id surviving corner", two neighbours
 * that number the face differently still agree. Both variants always have the
 * same child count, so a caller can size its output before resolving the tie.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

/// One `(cell type, split-edge mask)` subdivision.
struct RefineMaskTemplate {
    /// The children, each a list of local node ids (see the file comment).
    /// Empty exactly when the mask is not admissible.
    std::vector<std::vector<std::uint8_t>> mChildren;
    /// The other diagonal, or empty when this mask makes no such choice.
    std::vector<std::vector<std::uint8_t>> mChildrenAlt;
    /// Local **corner** ids deciding between the two variants: use
    /// `mChildrenAlt` when the global node id at `mTieB` is smaller than the one
    /// at `mTieA`. Meaningless unless `HasVariant()`.
    std::uint8_t mTieA = 0;
    std::uint8_t mTieB = 0;

    /// Whether this mask has a same-type subdivision at all.
    bool IsAdmissible() const { return !mChildren.empty(); }
    /// Whether resolving this template needs the global-id tie-break.
    bool HasVariant() const { return !mChildrenAlt.empty(); }
    /// The child count, the same for both variants.
    std::size_t NumChildren() const { return mChildren.size(); }
};

/**
 * @brief Whether `refine` has subdivision templates for a cell type.
 * @param Type The cell type to query.
 * @return `true` for `line`, `triangle`, `quad`, `tetra`, `wedge`,
 *   `hexahedron`; `false` for every higher-order type, `pyramid` and the rest.
 */
MESHIOPLUSPLUS_API bool refine_type_supported(CellType Type);

/**
 * @brief The mask with every edge of a cell type split (its uniform, "red"
 * refinement), or `0` for an unsupported type.
 * @param Type The cell type to query.
 */
MESHIOPLUSPLUS_API std::uint16_t refine_full_mask(CellType Type);

/**
 * @brief The smallest admissible superset of a split-edge mask.
 *
 * Monotone and idempotent, so iterating it over a mesh converges to a unique
 * fixed point regardless of the order cells are visited in.
 *
 * @param Type The cell type.
 * @param Mask The edges currently split, as a bitmask over
 *   `cell_refine_edges(Type)`.
 * @param Propagate When `true`, any non-empty mask promotes straight to the
 *   full split instead (`RefineClosure::Propagate`).
 * @return The promoted mask, or `0` for an unsupported type.
 */
MESHIOPLUSPLUS_API std::uint16_t refine_promote_mask(CellType Type, std::uint16_t Mask,
                                                     bool Propagate);

/**
 * @brief The subdivision for one `(type, mask)` pair.
 * @param Type The cell type.
 * @param Mask An **admissible** mask (call `refine_promote_mask` first).
 * @return Reference to the process-wide template. A non-admissible mask, or an
 *   unsupported type, yields a template whose `IsAdmissible()` is `false`.
 */
MESHIOPLUSPLUS_API const RefineMaskTemplate& refine_mask_template(CellType Type,
                                                                  std::uint16_t Mask);

}  // namespace detail
}  // namespace meshioplusplus
