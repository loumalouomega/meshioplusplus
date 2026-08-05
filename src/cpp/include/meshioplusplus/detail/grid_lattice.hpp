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
 * @file detail/grid_lattice.hpp
 * @brief The single owner of the **regular hexahedron lattice**: an axis-aligned
 * box of `nx * ny * nz` cells over a corner point grid, as an ordinary `Mesh`.
 *
 * Three consumers need exactly this object and must not disagree about its
 * numbering:
 *  - `grid()`, the public primitive constructor;
 *  - `voxelize()`, whose dense output *is* a lattice;
 *  - a future `read_vti`, since VTK ImageData is an implicit lattice that has to
 *    be expanded into explicit cells before it can be a `Mesh`.
 *
 * ### The numbering is inherited, not invented
 *
 * Points run **x fastest, then y, then z**, and cells run in that same order:
 *
 * ```
 * vid(i, j, k) = (k * (ny + 1) + j) * (nx + 1) + i
 * cid(i, j, k) = (k *  ny      + j) *  nx      + i
 * ```
 *
 * with each cell's eight nodes in the meshio/VTK `hexahedron` winding (bottom
 * ring counter-clockwise, then the top ring), so `detail/cell_faces.hpp`'s
 * "base (0,1,2,3) normal points toward the top (4,5,6,7)" holds and every cell of
 * a positively-spaced lattice is a right parallelepiped with unit scaled
 * Jacobian. That is a free structural oracle: a transposed axis shows up as an
 * inverted cell, not as a subtly wrong picture.
 *
 * This is the ordering the repo's existing C++/Python byte-identity fixtures
 * already agree on (`tests/cpp/test_refine.cpp`'s `hex_grid`, and its twin in
 * `tests/python/test_refine.py`). A handful of older fixtures use a transposed
 * linearization; they are deliberately left alone, since changing them would
 * rewrite committed expected values for no gain.
 *
 * x-fastest is also VTK ImageData's own ordering, which is what will make the
 * `.vti` writer a straight copy and what makes `cell_data` reshape to
 * `arr[z, y, x]` — the C-order tensor layout voxel tooling expects.
 *
 * ### Two properties worth stating as contracts
 *
 * **Shared corners are deduplicated arithmetically.** Neighbouring cells reference
 * the same node because the index formula says so — there is no hash, no
 * tolerance and no welding pass, so conformity and determinism here are
 * structural rather than defended.
 *
 * **Coordinates are `origin + index * spacing`, never accumulated.** An
 * accumulating `x += h` makes the last plane's coordinates depend on `nx` in the
 * low bits, which would make the numpy twin's loop-versus-`arange` choice
 * observable.
 *
 * Free functions in `meshioplusplus::detail`, called once per operation rather
 * than per element, so their bodies live in
 * `src/cpp/src/detail/grid_lattice.cpp` — the convention
 * `detail/cell_adjacency.hpp` and `detail/subset.hpp` follow. Anonymous-namespace
 * helpers there are prefixed `lat_`; note `grid_` is already owned by
 * `detail/spatial_hash.hpp`'s `grid_quantize`, which is a different thing
 * entirely (a bucket hash, not a mesh).
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The geometry of a regular lattice: where it starts, how big a cell is,
 * and how many there are along each axis.
 *
 * A cell count of zero on any axis describes an empty lattice, which is a legal
 * (if useless) request rather than an error -- `grid(0, 0, 0)` is an empty mesh,
 * not a throw.
 */
struct LatticeSpec {
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};   ///< The lo corner.
    std::array<double, 3> mSpacing{{1.0, 1.0, 1.0}};  ///< Cell size per axis.
    std::array<std::int64_t, 3> mDims{{0, 0, 0}};     ///< Cell counts per axis.
};

/// Number of points a lattice has: `(nx + 1)(ny + 1)(nz + 1)`.
MESHIOPLUSPLUS_API std::int64_t lattice_num_points(const LatticeSpec& rSpec);

/// Number of cells a lattice has: `nx * ny * nz`.
MESHIOPLUSPLUS_API std::int64_t lattice_num_cells(const LatticeSpec& rSpec);

/**
 * @brief Build the lattice from `rSpec` bounds.
 * @param rLo the lo corner of the box to cover.
 * @param rHi the hi corner; each component must be >= the matching `rLo`.
 * @param rDims cell counts per axis.
 * @return the spec whose origin is `rLo` and whose spacing divides the box
 *         exactly into `rDims` cells. An axis with zero cells gets zero spacing.
 */
MESHIOPLUSPLUS_API LatticeSpec lattice_from_bounds(const std::array<double, 3>& rLo,
                                                   const std::array<double, 3>& rHi,
                                                   const std::array<std::int64_t, 3>& rDims);

/**
 * @brief Build the lattice covering `rLo`..`rHi` with cells of at most
 *        `rCellSize`, growing the box up rather than shrinking a cell.
 *
 * Each axis gets `ceil(extent / cell)` cells (at least one when the extent is
 * positive), so the lattice covers the box completely and its hi corner may sit
 * slightly beyond `rHi`. Covering is the property callers depend on; an exactly
 * fitting box would need a per-axis cell size that is not the one asked for.
 */
MESHIOPLUSPLUS_API LatticeSpec lattice_from_cell_size(const std::array<double, 3>& rLo,
                                                      const std::array<double, 3>& rHi,
                                                      const std::array<double, 3>& rCellSize);

/**
 * @brief Materialize a lattice as a `Mesh` with one `hexahedron` cell block.
 * @param rSpec the lattice geometry.
 * @return the mesh; points are Float64 `(n, 3)`, connectivity Int64 `(m, 8)`.
 *         An empty lattice yields a mesh with no points and no cell block at all
 *         rather than a zero-row one, so `NumCellBlocks() == 0` says "nothing
 *         here" without a caller having to look at the row count.
 */
MESHIOPLUSPLUS_API Mesh lattice_build_mesh(const LatticeSpec& rSpec);

/**
 * @brief The axis-aligned bounding box of a mesh's points.
 *
 * Hoisted verbatim out of `compute_stats`, which computes the same reduction as
 * part of a much larger report. `min`/`max` are associative and exact, so the
 * chunked parallel form and a serial one agree bit-for-bit -- unlike the
 * centroid sum in the same loop, which is order-dependent and deliberately stays
 * in `stats.cpp`.
 *
 * @param rMesh the mesh to measure.
 * @param rLo receives the lo corner; unchanged axes beyond `PointDim()` read 0.
 * @param rHi receives the hi corner.
 * @return false when the mesh has no points, in which case `rLo`/`rHi` are zeroed
 *         -- an empty bbox has no defensible value and callers must decide.
 */
MESHIOPLUSPLUS_API bool point_bbox(const Mesh& rMesh, std::array<double, 3>& rLo,
                                   std::array<double, 3>& rHi);

}  // namespace detail
}  // namespace meshioplusplus
