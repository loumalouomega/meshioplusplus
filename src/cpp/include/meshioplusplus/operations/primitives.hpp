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
 * `icosa_sphere`, `cylinder`, the planar domains `reference_triangle`, `ngon`, `annulus`,
 * `lshape`, `plate_with_hole` and `honeycomb`, and `reference_cells`.
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
 *    `rim` on a disk, `side`/`top`/`bottom` on a cylinder, `inner`/`outer` on an
 *    annulus, `hole` on a plate);
 *  - `Order == 2` places the mid-edge nodes on the exact geometry (a sphere or
 *    a cylinder is genuinely curved), appended after the corners in
 *    first-appearance order over the cell-major edge list of `cell_refine_edges`.
 *
 * Seams and poles are welded points, never duplicates. The Python twin
 * (`_primitives.py`) is byte-identical for `rectangle`, `box`, `icosa_sphere`,
 * `reference_triangle`, `lshape`, `honeycomb` and `reference_cells` (index arithmetic and
 * `sqrt` only); `disk`, `uv_sphere`, `cylinder`, `ngon`, `annulus` and
 * `plate_with_hole` use `sin`/`cos`, which differ between libms, so the twins
 * agree on topology and to a few ulps on coordinates.
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
 * @brief The reference triangle `(0,0) (1,0) (0,1)`, scaled by @p rSize, refined
 * into `N^2` triangles, at `z = origin[2]`.
 *
 * Nodes run in barycentric lattice order: row `j` (y) holds the nodes `i = 0 ..
 * N - j` (x fastest). Cells run row by row, each lattice triangle followed by the
 * inverted one to its right.
 *
 * @param N divisions per edge; zero yields an empty mesh.
 * @param rSize the leg lengths `(lx, ly)`, positive.
 * @param rOrigin the right-angle corner.
 * @param Order 1 or 2 (`triangle6`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `ymin`, `xmin`, `hypotenuse`.
 */
MESHIOPLUSPLUS_API Mesh reference_triangle(std::int64_t N = 4,
                                           const std::array<double, 2>& rSize = {{1.0, 1.0}},
                                           const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                                           int Order = 1, std::int64_t MaxCells = 20000000);

/**
 * @brief A regular polygon fanned from its centre (`Rings == 1`), or refined into
 * `Rings` rings of `Sides * k` nodes whose triangles are those of `disk`.
 * @param Sides the number of sides, at least 3; zero yields an empty mesh.
 * @param Radius the circumradius, positive.
 * @param Rings radial divisions; zero yields an empty mesh.
 * @param rCenter the centre.
 * @param Order 1 or 2 (`triangle6`; the sides are straight).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with one side region, `rim`.
 */
MESHIOPLUSPLUS_API Mesh ngon(std::int64_t Sides = 6, double Radius = 1.0, std::int64_t Rings = 1,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             int Order = 1, std::int64_t MaxCells = 20000000);

/**
 * @brief An annulus, or a sector of one, in the `z = center[2]` plane.
 *
 * Nodes run ring by ring from the inner radius out, each ring counter-clockwise
 * from angle 0. A full ring (`Angle == 2 pi`) welds its last node onto its first;
 * a sector keeps both end rays as distinct node chains, the seam pair that
 * `periodic` links.
 *
 * @param InnerRadius the inner radius, positive.
 * @param OuterRadius the outer radius, above @p InnerRadius.
 * @param Ntheta cells around; at least 3 for a full ring.
 * @param Nr cells across the radial gap.
 * @param Angle the opening angle in `(0, 2 pi]`.
 * @param rCenter the centre.
 * @param rCellType `"quad"` or `"triangle"` (each quad cut along its first
 *        diagonal).
 * @param Order 1 or 2 (`quad8` / `triangle6`; the rims are curved).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `inner`, `outer` and, for a sector, `start`
 *         and `end`.
 */
MESHIOPLUSPLUS_API Mesh annulus(double InnerRadius = 0.5, double OuterRadius = 1.0,
                                std::int64_t Ntheta = 16, std::int64_t Nr = 4,
                                double Angle = 6.283185307179586,
                                const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                const std::string& rCellType = "quad", int Order = 1,
                                std::int64_t MaxCells = 20000000);

/**
 * @brief The square `[-Size/2, Size/2]^2` around @p rCenter minus its fourth
 * quadrant (`x > 0`, `y < 0`), the re-entrant corner at the centre.
 * @param N cells along half a side; the mesh has `3 N^2` quads.
 * @param Size the side of the bounding square, positive.
 * @param rCenter the centre of the bounding square.
 * @param rCellType `"quad"` or `"triangle"` (each quad cut along its first
 *        diagonal).
 * @param Order 1 or 2 (`quad8` / `triangle6`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `xmin`, `xmax`, `ymin`, `ymax`, `notch_x`
 *         (the notch wall on `x = 0`) and `notch_y` (the one on `y = 0`).
 */
MESHIOPLUSPLUS_API Mesh lshape(std::int64_t N = 4, double Size = 2.0,
                               const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                               const std::string& rCellType = "quad", int Order = 1,
                               std::int64_t MaxCells = 20000000);

/**
 * @brief A square plate around a circular hole, in four quad blocks.
 *
 * Nodes run ring by ring from the hole out (`k`), `4 Nt` positions around each
 * (`p`, counter-clockwise from the south-east corner direction). A node blends the
 * hole circle into the square linearly in the ring index.
 *
 * @param Radius the hole radius, positive and below half of @p Size.
 * @param Size the side of the plate.
 * @param Nt cells around the hole per 90 degree block.
 * @param Nr cells from the hole to the outer square.
 * @param rCenter the centre of the hole.
 * @param rCellType `"quad"` or `"triangle"` (each quad cut along its first
 *        diagonal).
 * @param Order 1 or 2 (`quad8` / `triangle6`; the hole is curved).
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with side regions `hole`, `xmin`, `xmax`, `ymin`, `ymax`.
 */
MESHIOPLUSPLUS_API Mesh plate_with_hole(double Radius = 0.25, double Size = 1.0,
                                        std::int64_t Nt = 8, std::int64_t Nr = 4,
                                        const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                        const std::string& rCellType = "quad", int Order = 1,
                                        std::int64_t MaxCells = 20000000);

/**
 * @brief A hexagonal tiling of pointy-top hexagons as one ragged `polygon` block.
 *
 * Hexagons run row by row (x fastest), odd rows shifted half a hexagon right;
 * nodes are numbered by first appearance over the cells, each hexagon's corners
 * counter-clockwise from the bottom one.
 *
 * @param rDims hexagons per row and rows `(nx, ny)`; a zero yields an empty mesh.
 * @param Edge the side length, positive.
 * @param rOrigin the lower-left corner of the bounding box.
 * @param MaxCells the cell budget; non-positive disables it.
 * @return a mesh with one side region, `boundary` (the hexagon edges seen once).
 */
MESHIOPLUSPLUS_API Mesh honeycomb(const std::array<std::int64_t, 2>& rDims = {{4, 4}},
                                  double Edge = 1.0,
                                  const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
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
