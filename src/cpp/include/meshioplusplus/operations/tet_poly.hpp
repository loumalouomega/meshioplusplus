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
 * @file operations/tet_poly.hpp
 * @brief Tetrahedral <-> polyhedral conversion: `tet_to_poly` builds the median
 * dual of a tetrahedral mesh (one polyhedron per vertex), and `poly_to_tet`
 * decomposes a polyhedral mesh back into tetrahedra.
 *
 * ### `tet_to_poly`: the median dual
 *
 * This is the tet-to-poly conversion CFD meshers perform: every vertex of the
 * tetrahedral mesh becomes one polyhedral cell, bounded by
 *
 * - one **interior face** per edge `(v, w)` leaving the vertex: the ring of
 *   the tetrahedron centroids and triangle centroids around that edge, wound
 *   so that its normal points from `v` to `w` (the cell of `w` gets the same
 *   ring reversed, so the result is conforming); where the edge reaches the
 *   boundary the ring is open and closes through the edge midpoint;
 * - **boundary faces** where the vertex lies on the boundary: each boundary
 *   triangle `(v, a, b)` contributes the planar quadrilateral
 *   `v, mid(v,a), centroid, mid(v,b)`, and the quadrilaterals around `v` are
 *   fused into one polygon per *run* of consecutive triangles whose normals
 *   differ by at most `mFeatureAngleDeg` (and that lie on the same side and in
 *   the same input triangle block). A run that goes all the way round drops
 *   `v` itself; any other run keeps it, so feature edges and corners survive.
 *
 * A dual cell is built per **(vertex, tetra block, face-connected component of
 * that block's tetrahedra around the vertex)**: a vertex on the interface
 * between two blocks yields one cell on each side (a material interface stays
 * sharp), and a pinched vertex one cell per lobe.
 *
 * With the default `mFeatureAngleDeg = 1` every boundary face is planar, and
 * the dual's total volume equals the tetrahedral mesh's up to rounding.
 *
 * #### Output structure
 *
 * - Points: the input points keep their indices (an interior vertex is left
 *   unreferenced; `clean(remove_orphans=True)` compacts them), followed by
 *   the midpoints of the edges that reach a boundary (sorted edge order), the
 *   centroid of every triangle face (global face order) and the centroid of
 *   every tetrahedron (compact cell order). `point_data` is extended to the
 *   new points by averaging their source nodes.
 * - Each `tetra` block becomes one `polyhedron` block, in place.
 * - Each `triangle` block becomes a `polygon` block holding the dual boundary
 *   faces over its triangles (from the lower-index block's side on an
 *   interface). A triangle that is not on a boundary is dropped with a
 *   warning, and so is every other block.
 * - `cell_data` of a polyhedron block is the |volume|-weighted mean over the
 *   cell's tetrahedra for floating-point arrays, and the value of the largest
 *   tetrahedron (ties: lowest index) for integer ones; polygon blocks weight
 *   by triangle area the same way. `field_data` passes through.
 * - Point regions survive (points are not renumbered); Cell and Side regions
 *   are kept as empty groups with a warning, since a dual cell belongs to a
 *   vertex, not a tetrahedron.
 *
 * ### `poly_to_tet`
 *
 * `convert_cells(Simplexify)` with one addition: a polyhedron that is really a
 * tetrahedron (four triangular faces over four distinct nodes, closed and
 * orientable) becomes that one positively oriented `tetra`, with no new
 * points. Every other cell takes simplexify's own path: the volume-conserving
 * corner-average fan for a polyhedron, the fixed templates for a hexahedron,
 * wedge or pyramid, a fan for a polygon. `poly_to_tet(tet_to_poly(m))` is a
 * valid, volume-conserving tetrahedral mesh, but not `m`.
 *
 * Both operations are standard C++ over the uniform mesh API, so they compile
 * under every mesh backend.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Options for `tet_to_poly`.
struct TetToPolyOptions {
    /// Consecutive boundary triangles around a vertex share one dual boundary
    /// face while their normals differ by at most this many degrees; in
    /// `[0, 180)`.
    double mFeatureAngleDeg = 1.0;
    /// Also copy each vertex's `point_data` onto its dual cells (and dual
    /// boundary polygons) as `cell_data` of the same name; a name already used
    /// by a `cell_data` array is skipped with a warning.
    bool mPointDataToCellData = false;
    /// Attach an Int64 `dual:vertex` `cell_data` array: the input vertex each
    /// output cell was built around.
    bool mRecordVertexIds = false;
};

/// The result of `tet_to_poly`.
struct TetToPolyResult {
    /// The dual mesh.
    Mesh mMesh;
    /// Int64, flat, shape `(total output cell count,)`: output global
    /// (block-major) cell index -> the input vertex it was built around.
    NDArray mVertexMap;
    /// Points appended after the input points.
    std::int64_t mNumNewPoints = 0;
    /// Boundary faces of the dual cells (each side of an interface counted).
    std::int64_t mNumBoundaryFaces = 0;
    /// Vertices that yielded more than one dual cell (block interfaces and
    /// pinched vertices).
    std::int64_t mNumSplitVertices = 0;
};

/**
 * @brief Build the median dual of a tetrahedral mesh: one polyhedron per vertex.
 * @param rMesh the tetrahedral mesh (never modified).
 * @param rOptions the feature angle and data options.
 * @return the dual mesh, the per-cell source vertex and the counts.
 * @throws std::invalid_argument when a volume block is not linear `tetra`, the
 *         mesh has no tetrahedra, its points are not 3D, a face is shared by
 *         three or more tetrahedra, or the feature angle is out of range.
 */
MESHIOPLUSPLUS_API TetToPolyResult tet_to_poly(const Mesh& rMesh,
                                               const TetToPolyOptions& rOptions = {});

/// Options for `poly_to_tet`.
struct PolyToTetOptions {
    /// Attach an Int64 `convert:parent_cell` `cell_data` array recording, per
    /// output cell, the index of the input cell it came from within its block.
    bool mRecordParentIds = false;
};

/// The result of `poly_to_tet`: `ConvertCellsResult`'s fields plus counts.
struct PolyToTetResult {
    /// The tetrahedral (simplex) mesh.
    Mesh mMesh;
    /// Int64 `(num_points_in,)`, input point -> output point.
    NDArray mPointMap;
    /// Per input block, Int64 `(num_cells_in_block,)`: input cell -> the index
    /// of its first child in the corresponding output block.
    std::vector<NDArray> mCellMaps;
    /// Polyhedra recognised as tetrahedra and emitted as exactly one `tetra`.
    std::int64_t mNumExactTets = 0;
    /// Polyhedra decomposed by the corner-average fan.
    std::int64_t mNumFannedCells = 0;
};

/**
 * @brief Decompose a polyhedral (or mixed) mesh into tetrahedra, recognising
 * polyhedra that already are tetrahedra.
 * @param rMesh the mesh to decompose (never modified).
 * @param rOptions whether to record parent ids.
 * @return the simplex mesh plus the point/cell index maps and the counts.
 * @throws std::invalid_argument when a polyhedron's faces are not a closed
 *         orientable surface.
 */
MESHIOPLUSPLUS_API PolyToTetResult poly_to_tet(const Mesh& rMesh,
                                               const PolyToTetOptions& rOptions = {});

}  // namespace meshioplusplus
