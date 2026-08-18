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
 * @file operations/subdivide.hpp
 * @brief Polyhedral refinement: split every 3D cell into one polyhedral child
 * per face, by connecting each face's boundary to a new interior point.
 *
 * `refine` and `decimate` both raise by name on a polyhedron, pointing at
 * `convert_cells(Simplexify)` — both are built on fixed same-type subdivision
 * templates, and an arbitrary polyhedron has none. `subdivide` is the answer
 * for polyhedral refinement specifically: rather than a per-type template
 * table, it goes through `detail::cell_rings` — the same uniform face-ring
 * abstraction `gradient` and `detail::cell_measure` already use — which treats
 * a tabulated type (`tetra`, `hexahedron`, `wedge`, `pyramid`, and their
 * quadratic variants, reduced to corners) and an existing `polyhedronN` block
 * identically. So `subdivide` needs **no per-type table at all**: the same
 * code handles every 3D cell type the mesh already supports.
 *
 * ### The construction
 *
 * For each eligible 3D cell:
 *
 * 1. `detail::cell_rings` gives the cell's faces as global node-id rings,
 *    uniformly whether the cell is tabulated or already a polyhedron.
 * 2. `detail::orient_rings` repairs winding so every face points outward;
 *    a cell whose faces are not a closed orientable surface **throws**,
 *    naming the cell — the same guard `convert_cells(Simplexify)` uses for
 *    its own polyhedron branch.
 * 3. One new point is added per cell: the plain arithmetic mean of the
 *    cell's own corner node coordinates (`rings.mNodes`) — the same ad hoc
 *    corner average `convert_cells.cpp`'s polyhedron branch already computes
 *    inline. This is deliberately **not** `poly_measure().mCentroid`, the
 *    volume centroid, a different point: the corner average is what makes
 *    the children's total volume *literally the same sum* as the parent's
 *    own `poly_measure` computation, rather than merely equal by the
 *    divergence theorem (both true; only the former is checkable by exact
 *    equality rather than a tolerance).
 * 4. For each face, one polyhedron child is emitted whose boundary is that
 *    face **unchanged** (original winding, so a neighbouring cell across it
 *    still sees the identical face) plus one new triangle per face edge,
 *    connecting that edge to the cell's new interior point.
 *
 * This makes the result **automatically conforming** — no closure, no 2:1
 * balance, no `refine:hanging`/`refine:entity` analogue is needed, because a
 * shared face between two input cells is never touched, only the interior of
 * each cell is subdivided.
 *
 * ### Output structure
 *
 * **One polyhedron output block per input 3D block**, with genuinely mixed
 * cell shapes inside — `AddPolyhedronBlock` stores cells as ragged CSR with
 * no constraint that they share a node or face count, so there is no need to
 * group children by distinct node count into `polyhedronN` blocks the way
 * some *readers* (CGNS's `NFACE_n`, OpenFOAM, EnSight) do for their own
 * format-compatibility reasons — that convention is not a requirement of the
 * uniform mesh API. Keeping one output block per input block also keeps the
 * block correspondence 1:1, which is exactly what `CellMapKind::FirstChild`
 * needs to carry regions correctly (see `mCellMaps` below).
 *
 * Non-3D blocks (2D/1D boundary markers) and 3D blocks with no
 * `detail::cell_faces` row (the 3D Lagrange family) pass through unchanged.
 *
 * ### Regions
 *
 * Point and Cell regions survive (via `CellMapKind::FirstChild`, the same
 * shape `convert_cells` already uses for its own one-to-many splits — a
 * parent's children occupy a contiguous run in the corresponding output
 * block). Named **Side** regions do not: `FirstChild` drops them
 * unconditionally, since a child's facets have no correspondence with the
 * parent's — the same limitation `convert_cells(Simplexify)` already has and
 * documents, not a new gap this operation introduces.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is
 * not in the format registry.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Options for `subdivide`.
struct SubdivideOptions {
    /// Attach an Int64 `subdivide:parent_cell` `cell_data` array recording,
    /// per output cell, the index of the input cell it came from *within its
    /// own block* (blocks correspond 1:1).
    bool mRecordParentIds = false;
};

/// The result of `subdivide`: the subdivided mesh plus the cell index map.
struct SubdivideResult {
    /// The subdivided mesh.
    Mesh mMesh;
    /// Per input block, Int64 shape `(num_cells_in_block,)`: input cell -> the
    /// index of its **first** child in the corresponding output block
    /// (`CellMapKind::FirstChild`; see the file docs). A non-3D or otherwise
    /// ineligible block passes through unchanged, so its entries are the
    /// identity. A cell whose rings could not be built at all (out-of-range
    /// connectivity) contributes no children, giving it an empty run.
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Split every eligible 3D cell into one polyhedral child per face.
 * @param rMesh the mesh to subdivide.
 * @param rOptions whether to record parent ids (defaults to false).
 * @return the subdivided mesh plus the per-block cell index map.
 * @throws std::invalid_argument when a cell's faces are not a closed
 *         orientable surface (naming the cell and block).
 */
MESHIOPLUSPLUS_API SubdivideResult subdivide(const Mesh& rMesh,
                                             const SubdivideOptions& rOptions = {});

}  // namespace meshioplusplus
