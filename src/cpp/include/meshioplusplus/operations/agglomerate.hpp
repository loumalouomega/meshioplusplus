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
 * @file operations/agglomerate.hpp
 * @brief Polyhedral coarsening: merge groups of cells into single larger
 * polyhedral cells, the many-to-one counterpart to `subdivide`.
 *
 * `decimate` raises by name on a polyhedron, pointing at
 * `convert_cells(Simplexify)` — its fixed-template QEM edge collapse has no
 * analogue for merging arbitrary polyhedral cells. `agglomerate` is a
 * genuinely different algorithm: greedy seed-and-grow over the mesh's
 * shared-face dual (`detail::build_global_faces`), absorbing face-adjacent
 * neighbours into a group until it reaches a target size, then emitting one
 * polyhedron per group whose faces are exactly that group's *external*
 * boundary — every face shared by two members of the same group cancels out
 * of the result, by construction, never by a tolerance.
 *
 * ### The construction
 *
 * 1. `detail::build_global_faces(rMesh)` gives the mesh's volume cells as a
 *    compact-space face dual (`mOwner`/`mNeighbour` per face). A mesh with any
 *    non-manifold face (used by three or more cells) is **refused** — the
 *    owner/neighbour classification below is only well-defined on a
 *    2-manifold face, and guessing would silently misclassify a boundary.
 * 2. Greedy seed-and-grow, serial and deterministic: cells are seeded in
 *    ascending compact-id order; a group absorbs its unclaimed face-neighbour
 *    with the largest *accumulated* shared-face area (summed over every face
 *    the group's current members share with that neighbour), ties broken by
 *    ascending compact id, until `mTargetGroupSize` is reached or no unclaimed
 *    neighbour remains (a short group at a mesh boundary or pocket is
 *    expected, not an error).
 * 3. Each group emits **one** polyhedron cell: walk every member's faces:
 *    a face whose *other* side is also in the group is internal and dropped
 *    (this happens from both sides, so it is never emitted twice); every
 *    other face is the group's own external boundary and is kept, wound
 *    forward or reversed exactly as `detail::GlobalFaces::mCellFaces`' sign
 *    already records for that member — no new orientation logic needed, since
 *    the merged cell simply inherits each member's own local winding at its
 *    surviving faces.
 *
 * This is deliberately face-adjacency (`build_global_faces`'s `mOwner`/
 * `mNeighbour`), never `detail::cell_adjacency.hpp`'s node-adjacency (used by
 * `partition`'s ghost layers and `gradient`'s stencil): merging on shared-node
 * adjacency could fuse two cells touching only at a single pinch-point vertex,
 * producing a non-manifold union, and the two headers key their cells in
 * genuinely different index spaces (`cell_adjacency.hpp` global block-major,
 * `GlobalFaces` compact volume-cell) that would silently misindex if mixed.
 *
 * ### Output structure
 *
 * Non-volume blocks (2D/1D boundary markers, and any 3D block with no
 * `cell_faces` row) pass through unchanged, in their original relative
 * position. Every volume cell is consumed into **one** new `polyhedron`
 * block, emitted at the position the *first* original volume block occupied
 * — so a mesh whose volume cells already form one contiguous run keeps its
 * overall block order. Cells inside that block have whatever face/node count
 * their own group boundary produces; `AddPolyhedronBlock` stores genuinely
 * ragged CSR with no same-shape constraint, so there is no grouping-by-size
 * step to do (the same simplification `subdivide` already established).
 *
 * ### Points and data
 *
 * Points are **not** compacted: a group can leave interior nodes unreferenced
 * (a node shared only by the members' now-internal faces), the same
 * never-prune-or-renumber precedent `subdivide` set for its own orphan case.
 * `clean(remove_orphans=True)` is the documented follow-up for a caller who
 * wants a minimal point set. `point_data` and `field_data` therefore need no
 * remapping at all — they pass through unchanged.
 *
 * `cell_data` for a pass-through block is copied verbatim. For the merged
 * block, each group's row is its **first member's** row (ascending compact
 * id within the group) — the same keep-first convention `clean`'s point weld
 * already uses, since there is no single principled value for an array like a
 * material tag once several cells with (possibly) different values merge. An
 * array whose block count does not match the input mesh is not per-cell data
 * shaped this operation understands and is dropped with a warning rather than
 * guessed at.
 *
 * ### What this does not do (yet)
 *
 * Coplanar boundary-face merging (fusing two adjacent group-boundary faces on
 * the same plane into one larger polygon, rather than leaving the edge
 * between them) and a shape-quality (e.g. sphericity) absorption gate are
 * both deferred follow-ups, not shipped here — see the roadmap.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is
 * not in the format registry.
 */

// System includes
#include <cstddef>
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Options for `agglomerate`.
struct AgglomerateOptions {
    /// Approximate number of member cells to absorb into each output group;
    /// must be at least 1. `1` means every cell is its own group (an
    /// identity transform in everything but representation).
    std::size_t mTargetGroupSize = 8;
};

/// The result of `agglomerate`: the coarsened mesh plus the cell index map.
struct AgglomerateResult {
    /// The coarsened mesh.
    Mesh mMesh;
    /// Int64, **flat**, shape `(total input cell count,)` — input global
    /// (block-major) cell index -> output global cell index. Unlike
    /// `SubdivideResult`/`ConvertCellsResult`'s per-input-block
    /// `std::vector<NDArray>`, this is a single array: an output cell's index
    /// is a function of which group it joined, not which input block it came
    /// from.
    NDArray mCellMap;
};

/**
 * @brief Merge groups of cells into single larger polyhedral cells.
 * @param rMesh the mesh to coarsen.
 * @param rOptions the target group size (defaults to 8).
 * @return the coarsened mesh plus the flat cell index map.
 * @throws std::invalid_argument when `mTargetGroupSize == 0`, or when the
 *         mesh contains a face shared by three or more cells (non-manifold).
 */
MESHIOPLUSPLUS_API AgglomerateResult agglomerate(const Mesh& rMesh,
                                                 const AgglomerateOptions& rOptions = {});

}  // namespace meshioplusplus
