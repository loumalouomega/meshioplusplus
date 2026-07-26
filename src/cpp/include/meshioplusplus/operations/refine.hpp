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
 * @file refine.hpp
 * @brief Uniform mesh refinement: subdivide every cell into congruent children
 * of the **same** cell type, interpolating `point_data` onto the new nodes.
 *
 * This is the resolution-increasing counterpart to the resolution-preserving
 * `convert_cells` and the resolution-reducing `crop`/`clean`. One level applies
 * a fixed template per cell type:
 *
 *  - `line` -> 2 `line` (1 new mid-edge node)
 *  - `triangle` -> 4 `triangle` (3 mid-edge nodes; the standard 1-to-4 split)
 *  - `quad` -> 4 `quad` (4 mid-edge + 1 face-centre node)
 *  - `tetra` -> 8 `tetra` (6 mid-edge nodes)
 *  - `wedge` -> 8 `wedge` (9 mid-edge + 3 quad-face-centre nodes)
 *  - `hexahedron` -> 8 `hexahedron` (12 mid-edge + 6 face-centre + 1 body node)
 *
 * New nodes are placed at the midpoint of the entity that defines them and
 * carry the mean of that entity's corner values for every `point_data` array.
 * Original points keep their indices; new nodes are appended. Each parent's
 * `cell_data` row is replicated to its children.
 *
 * **Conformity.** Mid-edge nodes and quad-face-centre nodes are *shared*
 * between every cell that touches the entity, so the refined mesh has no
 * hanging nodes. Only the hexahedron body node is per-cell. (Sharing the face
 * centres is not optional: with a per-cell copy, two hexahedra meeting at a
 * face would reference distinct coincident nodes and the mesh would be
 * topologically split along every interior face.)
 *
 * Note that the tetrahedron's interior diagonal -- fixed at the opposite-edge
 * pair `(0,1)`-`(2,3)`, i.e. `tetra10` nodes 4-9 -- is chosen for determinism
 * only, **not** for conformity: it is strictly interior, so a face's
 * subdivision is fixed by that face's own mid-edge nodes whatever the
 * neighbour does. This is the opposite of `convert_cells`' hex-simplexify
 * diagonal 0-6, whose endpoints lie on the boundary and which therefore *is*
 * conformity-critical.
 *
 * **Volume and orientation.** Children inherit the parent's orientation, so a
 * well-oriented input refines to a well-oriented output (zero newly-inverted
 * cells). Volume is conserved exactly for `line`/`triangle`/`quad`/`tetra`
 * always, and for `wedge`/`hexahedron` when the parent is affine (a right
 * prism / parallelepiped). For a general trilinear hexahedron the eight
 * children's volumes do not sum to the parent's, because the parent's bilinear
 * faces are replaced by four different bilinear patches -- that is a property
 * of the geometry, not of this implementation.
 *
 * **Block structure is preserved 1:1**: the output has exactly
 * `NumCellBlocks()` blocks, in input order, which is what keeps the
 * one-array-per-block `cell_data` invariant trivially correct.
 *
 * Constructs that would break the same-type contract raise rather than guess:
 * higher-order cells (`tetra10`, ...; linearize first), `pyramid` (whose
 * uniform refinement is 6 pyramids + 4 tetrahedra; simplexify first), and
 * ragged polygon/polyhedron blocks.
 *
 * Determinism: the templates are fixed and the new-node numbering comes from a
 * **serial** dedup pass over a `parallel_for`-filled disjoint-slot buffer
 * (`src/cpp/src/operations/surface.cpp`'s phase-split idiom), never from a
 * concurrent hash insert. Output is byte-identical across mesh backends and
 * thread counts.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format -- it is not in the format
 * registry.
 */

// System includes
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Options for `refine`.
struct RefineOptions {
    /// How many times to apply the subdivision templates. `0` (or less) returns
    /// an unchanged clone; `n` multiplies the cell count of a supported block
    /// by `children_per_cell^n`.
    int mLevels = 1;
    /// Attach an Int64 `refine:parent_cell` `cell_data` array recording, per
    /// output cell, the index of the **original** input cell it descends from
    /// *within its own block* (blocks correspond 1:1). Across several levels
    /// this is the original ancestor, not the immediate parent.
    bool mRecordParentIds = false;
};

/// The result of `refine`: the refined mesh plus the index maps.
struct RefineResult {
    /// The refined mesh.
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> output point index.
    /// Refinement never prunes, so this is the identity; it is returned in full
    /// so callers need not depend on that.
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> the
    /// index of its **first** child in the corresponding output block. A cell's
    /// children are contiguous, so cell `c` owns
    /// `[map[c], c + 1 < n ? map[c + 1] : num_cells_out)`.
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Uniformly refine a mesh, subdividing every cell into same-type
 * children.
 * @param rMesh The mesh to refine (unchanged).
 * @param rOptions Level count and parent-id recording.
 * @return The refined mesh and its point/cell index maps.
 * @throws std::invalid_argument on a higher-order, `pyramid`, or ragged cell
 *   block, none of which can be subdivided into same-type children.
 */
MESHIOPLUSPLUS_API RefineResult refine(const Mesh& rMesh, const RefineOptions& rOptions = {});

}  // namespace meshioplusplus
