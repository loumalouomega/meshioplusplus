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
 * @file operations/primitives.hpp
 * @brief Parametric mesh generators: `rectangle`, `box`, `disk`, `uv_sphere`,
 * `icosa_sphere`, `cylinder` and `reference_cells`.
 *
 * Almost every operation transforms a mesh you already have; `grid`
 * (`operations/voxelize.hpp`) was the one exception. These generators build a
 * shape from parameters, with the contract `grid` already keeps:
 *
 *  - coordinates are computed from indices, never accumulated (an angle is
 *    `2*pi*k/n`, not the previous point rotated);
 *  - points and cells follow a documented order (see `doc/primitives.md`);
 *  - a `MaxCells` budget refuses a request by name instead of failing to
 *    allocate, and a zero count yields an empty mesh (points `(0, 3)`, no cell
 *    block);
 *  - surfaces are oriented outward, planar cells counter-clockwise seen from
 *    +z and volumes positively;
 *  - boundaries come back as named `Side` regions (`xmin` ... `zmax` on a box,
 *    `rim` on a disk, `side`/`top`/`bottom` on a cylinder);
 *  - `Order == 2` places the mid-edge nodes on the exact geometry (a sphere or
 *    a cylinder is genuinely curved), appended after the corners in
 *    first-appearance order over the cell-major edge list of `cell_refine_edges`.
 *
 * Seams and poles are welded points, never duplicates. The Python twin
 * (`_primitives.py`) is byte-identical for `rectangle`, `box`, `icosa_sphere`
 * and `reference_cells` (index arithmetic and `sqrt` only); `disk`, `uv_sphere`
 * and `cylinder` use `sin`/`cos`, which differ between libms, so the twins agree
 * on topology and to a few ulps on coordinates.
 *
 * @see doc/primitives.md
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief A regular 2-D grid of quadrilaterals or triangles at `z = origin[2]`.
 * @param rDims cell counts `(nx, ny)`; a zero yields an empty mesh.
 * @param rSize side lengths `(lx, ly)`, positive on every axis with cells.
 * @param rOrigin the lo corner.
 * @param rCellType `"quad"` or `"triangle"`.
 * @param rSplit for triangles: `"zigzag"` (alternating diagonals), `"up"`
 *        (diagonal from the lo corner), `"down"` (the other diagonal) or
 *        `"center"` (four triangles around an added cell-centre node).
 * @param Order 1 or 2 (`quad8` / `triangle6`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `xmin`, `xmax`, `ymin`, `ymax`.
 * @throws std::invalid_argument on a bad argument or an exceeded budget.
 */
MESHIOPLUSPLUS_API Mesh rectangle(const std::array<std::int64_t, 2>& rDims = {{1, 1}},
                                  const std::array<double, 2>& rSize = {{1.0, 1.0}},
                                  const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                                  const std::string& rCellType = "quad",
                                  const std::string& rSplit = "zigzag", int Order = 1,
                                  std::int64_t MaxCells = 20000000);

/**
 * @brief A regular 3-D grid of hexahedra, tetrahedra (six per hexahedron, the
 * Kuhn split around the main diagonal, which conforms across cells) or wedges
 * (two per hexahedron).
 * @param rDims cell counts `(nx, ny, nz)`; a zero yields an empty mesh.
 * @param rSize side lengths, positive on every axis with cells.
 * @param rOrigin the lo corner.
 * @param rCellType `"hexahedron"`, `"tetra"` or `"wedge"`.
 * @param Order 1 or 2 (`hexahedron20` / `tetra10` / `wedge15`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `xmin`, `xmax`, `ymin`, `ymax`, `zmin`,
 *         `zmax`.
 */
MESHIOPLUSPLUS_API Mesh box(const std::array<std::int64_t, 3>& rDims = {{1, 1, 1}},
                            const std::array<double, 3>& rSize = {{1.0, 1.0, 1.0}},
                            const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                            const std::string& rCellType = "hexahedron", int Order = 1,
                            std::int64_t MaxCells = 20000000);

/**
 * @brief A disk in the `z = center[2]` plane.
 * @param Radius the radius, positive.
 * @param N quads: cells across a diameter (an `N x N` grid pushed onto the disk
 *        by the squircle map). triangles: rings between the centre and the rim
 *        (`6 N^2` triangles).
 * @param rCenter the centre.
 * @param rCellType `"quad"` or `"triangle"`.
 * @param Order 1 or 2 (`quad8` / `triangle6`; the rim is curved).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with one side region, `rim`.
 */
MESHIOPLUSPLUS_API Mesh disk(double Radius = 1.0, std::int64_t N = 8,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             const std::string& rCellType = "quad", int Order = 1,
                             std::int64_t MaxCells = 20000000);

/**
 * @brief A sphere surface in latitude / longitude rings, oriented outward.
 *
 * Point 0 is the north pole `(0, 0, +r)`, then `Nv - 1` rings of `Nu` points
 * from north to south, then the south pole.
 *
 * @param Radius the radius, positive.
 * @param Nu points per ring (at least 3).
 * @param Nv latitude bands (at least 2).
 * @param rCenter the centre.
 * @param rCellType `"triangle"`, or `"quad"` for quads between the rings with
 *        triangle caps.
 * @param Order 1 or 2 (`triangle6` / `quad8`).
 * @param MaxCells the cell budget; non-positive disables it.
 */
MESHIOPLUSPLUS_API Mesh uv_sphere(double Radius = 1.0, std::int64_t Nu = 16, std::int64_t Nv = 8,
                                  const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                  const std::string& rCellType = "triangle", int Order = 1,
                                  std::int64_t MaxCells = 20000000);

/**
 * @brief A sphere surface from an icosahedron whose triangles are split in four
 * @p Subdivisions times (`20 * 4^Subdivisions` triangles), oriented outward.
 * @param Radius the radius, positive.
 * @param Subdivisions the number of 1-to-4 splits (zero is the icosahedron).
 * @param rCenter the centre.
 * @param Order 1 or 2 (`triangle6`).
 * @param MaxCells the cell budget; non-positive disables it.
 */
MESHIOPLUSPLUS_API Mesh icosa_sphere(double Radius = 1.0, std::int64_t Subdivisions = 2,
                                     const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                     int Order = 1, std::int64_t MaxCells = 20000000);

/**
 * @brief A solid cylinder along z.
 * @param Radius the radius, positive.
 * @param Height the height, positive.
 * @param Nxy cells across a diameter (see `disk`).
 * @param Nz layers along the axis.
 * @param rCenter the centre of the cylinder.
 * @param rCellType `"hexahedron"` (the squircle disk, extruded) or `"wedge"`
 *        (the triangle disk, extruded).
 * @param Order 1 or 2 (`hexahedron20` / `wedge15`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `side`, `top`, `bottom`.
 */
MESHIOPLUSPLUS_API Mesh cylinder(double Radius = 1.0, double Height = 1.0, std::int64_t Nxy = 8,
                                 std::int64_t Nz = 8,
                                 const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                 const std::string& rCellType = "hexahedron", int Order = 1,
                                 std::int64_t MaxCells = 20000000);

/**
 * @brief One cell of each requested type in its reference coordinates, laid out
 * side by side along x without touching.
 * @param rTypes meshio type names; empty means every supported type in the
 *        order `vertex`, `line`, `triangle`, `quad`, `tetra`, `hexahedron`,
 *        `wedge`, `pyramid`, `line3`, `triangle6`, `quad8`, `tetra10`,
 *        `hexahedron20`, `wedge15`, `pyramid13`, `polyhedron`.
 * @param Pairs two cells of each type sharing a facet instead of one cell
 *        (linear types only).
 * @throws std::invalid_argument on an unknown or repeated type.
 */
MESHIOPLUSPLUS_API Mesh reference_cells(const std::vector<std::string>& rTypes = {},
                                        bool Pairs = false);

}  // namespace meshioplusplus
