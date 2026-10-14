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
 * @file operations/tessellate.hpp
 * @brief Isoparametric tessellation of curved (quadratic) cells onto a
 * reference lattice, with a recorded map from every output point and cell back
 * to the source mesh -- the C++ twin of `meshioplusplus.tessellate`.
 *
 * Every cell of a curved type -- `triangle6`, `quad8`, `quad9`, `tetra10` or
 * `hexahedron27` -- is sampled on a reference lattice and mapped through its
 * own quadratic shape functions; every other cell (linear types,
 * `hexahedron20`, `wedge15`/`18`, `pyramid13`/`14`, the VTK Lagrange family,
 * polygon and polyhedron blocks) passes through unchanged.
 *
 *  - **Lattice.** `quad8`/`quad9`/`hexahedron27` take `max(levels, 1)`
 *    divisions per axis; `triangle6`/`tetra10` take `levels` recursive
 *    same-type splits (`2^levels` divisions per edge), the tetrahedron through
 *    `refine`'s fixed interior diagonal.
 *  - **Output cells.** Curved surface cells become one `triangle` block (a
 *    sub-quad splits along the diagonal from its first-created corner to the
 *    opposite one); curved volume cells become one `tetra` block (a
 *    sub-hexahedron becomes 12 tetrahedra around its centroid). The blocks
 *    are the pass-through blocks in input order, then `triangle`, then `tetra`.
 *  - **Point identity.** A lattice point is keyed by the entity it lies on --
 *    an original corner, an edge (its two corner ids and the quantised
 *    fraction), a face (its canonical corner cycle and quantised coordinates)
 *    or the cell's interior -- so two cells sharing an edge or a face share its
 *    points. Output points are sorted by that key: the kept source points
 *    first, ascending, then edge, face and interior points. A source point no
 *    cell uses is dropped.
 *  - **Data.** `tessellate:source_point` (Int64 point data: the source point,
 *    or -1), `tessellate:source_cell` (Int64 cell data: the global block-major
 *    source cell) and `tessellate:sub_index` (Int64 cell data: the child index).
 *    With `mFields`, every point data array is interpolated (Float64) and every
 *    cell data row is copied to its children; field data is copied. With
 *    `mRecordStencil`, `tessellate:stencil` and `tessellate:weights` (`(P, 27)`)
 *    are added too.
 *  - **Regions.** Point regions are carried; cell regions keep their
 *    pass-through cells; side regions are dropped. Every drop is warned about.
 *
 * The output matches the numpy twin: the same points in the same order, the
 * same connectivity, stencils and data layout, with coordinates and
 * interpolated values equal within rounding.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Int64 point data: the source point a tessellated point was kept from, or -1.
inline constexpr const char* kTessellateSourcePointName = "tessellate:source_point";
/// Int64 cell data: the global (block-major) source cell of each output cell.
inline constexpr const char* kTessellateSourceCellName = "tessellate:source_cell";
/// Int64 cell data: a cell's index among its source cell's children.
inline constexpr const char* kTessellateSubIndexName = "tessellate:sub_index";
/// Int64 `(P, 27)` point data (opt-in): each point's stencil into the source points.
inline constexpr const char* kTessellateStencilName = "tessellate:stencil";
/// Float64 `(P, 27)` point data (opt-in): the matching interpolation weights.
inline constexpr const char* kTessellateWeightsName = "tessellate:weights";
/// The stencil width: `hexahedron27`'s 27 nodes.
inline constexpr std::size_t kTessellateStencilWidth = 27;

/// How `tessellate` subdivides.
struct TessellateOptions {
    /// Divisions per axis (tensor cells) or recursive splits (simplices); >= 0.
    std::int64_t mLevels = 2;
    /// False passes every cell through unchanged.
    bool mCurved = true;
    /// Interpolate point data and copy cell data onto the output.
    bool mFields = true;
    /// Also store the stencil and weights as point data.
    bool mRecordStencil = false;
};

/// The tessellated mesh, its stencils and its bookkeeping.
struct TessellateResult {
    /// The tessellated mesh, with the `tessellate:*` arrays.
    Mesh mMesh;
    /// `(P, 27)` row-major: each output point's stencil into the source points,
    /// padded with -1.
    std::vector<std::int64_t> mStencil;
    /// `(P, 27)` row-major: the matching weights, padded with 0.
    std::vector<double> mWeights;
    std::int64_t mNumSourcePoints = 0;
    std::int64_t mNumSourceCells = 0;
    /// Source cells of a curved type that were tessellated.
    std::int64_t mNumCurvedSourceCells = 0;
    /// Source cells passed through unchanged.
    std::int64_t mNumPassThroughSourceCells = 0;
    std::int64_t mNumPoints = 0;
    std::int64_t mNumCells = 0;
    /// Over every `tetra` block of the output: triangular facets used once.
    std::int64_t mNumBoundaryFacets = 0;
    /// Facets used twice.
    std::int64_t mNumInteriorFacets = 0;
    /// Facets used three or more times: zero for a conforming tessellation.
    std::int64_t mNumFacetsWithBadCount = 0;
};

/**
 * @brief Tessellate a mesh's curved cells.
 *
 * @param rMesh Any mesh.
 * @param rOptions See `TessellateOptions`.
 * @return The tessellated mesh and its bookkeeping.
 * @throws std::invalid_argument for negative levels, or a point data array
 *         whose row count is not the number of points.
 */
MESHIOPLUSPLUS_API TessellateResult tessellate(const Mesh& rMesh,
                                               const TessellateOptions& rOptions = {});

}  // namespace meshioplusplus
