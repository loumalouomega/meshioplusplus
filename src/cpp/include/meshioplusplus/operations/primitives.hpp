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
 * `lshape`, `plate_with_hole` and `honeycomb`, `reference_cells`, the curves and
 * point clouds `line`, `circle`, `helix`, `sphere_points` and `box_points`, the
 * surfaces `tetra_sphere`, `octa_sphere`, `cube_sphere`, `geo_sphere`, `torus`,
 * `moebius`, `tube`, `cone`, `capsule`, `box_surface` and `platonic`, and the
 * solids `ball`, `spherical_shell`, `pipe`, `fichera` and `solid_torus`.
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
 *    first-appearance order over the cell-major edge list of `cell_refine_edges`;
 *    `Complete` adds the face and body centres of the full-Lagrange types
 *    (`quad9`, `hexahedron27`, `wedge18`; `line3`, `triangle6` and `tetra10`
 *    already are), appended after every mid-edge node in first-appearance order
 *    over the cells;
 *  - `Order` 3 to 10 builds the `VTK_LAGRANGE_*` family instead, its nodes on
 *    the equispaced lattice of each cell mapped through the exact geometry and
 *    shared between neighbours, appended in first-appearance order over the
 *    cells, each cell's in VTK order.
 *
 * Seams and poles are welded points, never duplicates. The Python twin
 * (`_primitives.py`) is byte-identical for `rectangle`, `box`, `icosa_sphere`,
 * `reference_triangle`, `lshape`, `honeycomb` and `reference_cells` (index arithmetic and
 * `sqrt` only), as for `line`, `box_points`, `tetra_sphere`, `octa_sphere`,
 * `cube_sphere`, `geo_sphere`, `box_surface`, `platonic`, `ball` and `fichera`
 * (and `spherical_shell` with `Ratio == 1`); the shapes that take `sin`, `cos`
 * or `pow` agree on topology and to a few ulps on coordinates, since those
 * differ between libms.
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
 * @param Order 1 or 2 (`quad8` / `triangle6`); see the overloads with
 *        `Complete` for full-Lagrange and higher orders.
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

// ---------------------------------------------------------------------------
// The existing generators with `Complete` (full-Lagrange order 2) and `Order`
// up to 10; the shorter signatures above forward here with `Complete = false`.
// ---------------------------------------------------------------------------

MESHIOPLUSPLUS_API Mesh rectangle(const std::array<std::int64_t, 2>& rDims,
                                  const std::array<double, 2>& rSize,
                                  const std::array<double, 3>& rOrigin,
                                  const std::string& rCellType, const std::string& rSplit,
                                  int Order, std::int64_t MaxCells, bool Complete);
MESHIOPLUSPLUS_API Mesh box(const std::array<std::int64_t, 3>& rDims,
                            const std::array<double, 3>& rSize,
                            const std::array<double, 3>& rOrigin, const std::string& rCellType,
                            int Order, std::int64_t MaxCells, bool Complete);
MESHIOPLUSPLUS_API Mesh disk(double Radius, std::int64_t N, const std::array<double, 3>& rCenter,
                             const std::string& rCellType, int Order, std::int64_t MaxCells,
                             bool Complete);
MESHIOPLUSPLUS_API Mesh uv_sphere(double Radius, std::int64_t Nu, std::int64_t Nv,
                                  const std::array<double, 3>& rCenter,
                                  const std::string& rCellType, int Order, std::int64_t MaxCells,
                                  bool Complete);
MESHIOPLUSPLUS_API Mesh icosa_sphere(double Radius, std::int64_t Subdivisions,
                                     const std::array<double, 3>& rCenter, int Order,
                                     std::int64_t MaxCells, bool Complete);
MESHIOPLUSPLUS_API Mesh cylinder(double Radius, double Height, std::int64_t Nxy, std::int64_t Nz,
                                 const std::array<double, 3>& rCenter, const std::string& rCellType,
                                 int Order, std::int64_t MaxCells, bool Complete);
MESHIOPLUSPLUS_API Mesh reference_triangle(std::int64_t N, const std::array<double, 2>& rSize,
                                           const std::array<double, 3>& rOrigin, int Order,
                                           std::int64_t MaxCells, bool Complete);
MESHIOPLUSPLUS_API Mesh ngon(std::int64_t Sides, double Radius, std::int64_t Rings,
                             const std::array<double, 3>& rCenter, int Order, std::int64_t MaxCells,
                             bool Complete);
MESHIOPLUSPLUS_API Mesh annulus(double InnerRadius, double OuterRadius, std::int64_t Ntheta,
                                std::int64_t Nr, double Angle, const std::array<double, 3>& rCenter,
                                const std::string& rCellType, int Order, std::int64_t MaxCells,
                                bool Complete);
MESHIOPLUSPLUS_API Mesh lshape(std::int64_t N, double Size, const std::array<double, 3>& rCenter,
                               const std::string& rCellType, int Order, std::int64_t MaxCells,
                               bool Complete);
MESHIOPLUSPLUS_API Mesh plate_with_hole(double Radius, double Size, std::int64_t Nt,
                                        std::int64_t Nr, const std::array<double, 3>& rCenter,
                                        const std::string& rCellType, int Order,
                                        std::int64_t MaxCells, bool Complete);

// ---------------------------------------------------------------------------
// Curves and point clouds.
// ---------------------------------------------------------------------------

/**
 * @brief A straight segment from @p rStart to @p rEnd in @p N `line` cells.
 * @param N the cell count; a zero yields an empty mesh.
 * @param rStart the first point.
 * @param rEnd the last point, distinct from @p rStart.
 * @param Order 1, 2 (`line3`) or 3 to 10 (`VTK_LAGRANGE_CURVE`).
 * @param MaxCells the cell budget; non-positive disables it.
 * @param Complete accepted for symmetry (`line3` is already complete).
 * @return a mesh with point regions `start` and `end`.
 */
MESHIOPLUSPLUS_API Mesh line(std::int64_t N = 8,
                             const std::array<double, 3>& rStart = {{0.0, 0.0, 0.0}},
                             const std::array<double, 3>& rEnd = {{1.0, 0.0, 0.0}}, int Order = 1,
                             std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A closed polyline of @p N `line` cells on the circle of @p Radius in
 * the `z = center[2]` plane, counter-clockwise from +x; the seam is welded.
 * @param N the cell count, at least 3; a zero yields an empty mesh.
 * @param Order 1, 2 (`line3` on the arc) or 3 to 10.
 */
MESHIOPLUSPLUS_API Mesh circle(double Radius = 1.0, std::int64_t N = 16,
                               const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                               int Order = 1, std::int64_t MaxCells = 20000000,
                               bool Complete = false);

/**
 * @brief A helix (a spring) of @p N `line` cells around the z axis through
 * @p rCenter: `Turns` turns of @p Radius, rising @p Pitch per turn from
 * `z = center[2]`, counter-clockwise from +x.
 * @param N the cell count; more than `2 * Turns`, so no cell spans half a turn.
 * @param Order 1, 2 (`line3` on the helix) or 3 to 10.
 * @return a mesh with point regions `start` and `end`.
 */
MESHIOPLUSPLUS_API Mesh helix(double Radius = 1.0, double Pitch = 1.0, double Turns = 2.0,
                              std::int64_t N = 32,
                              const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                              int Order = 1, std::int64_t MaxCells = 20000000,
                              bool Complete = false);

/**
 * @brief @p N points on the sphere of @p Radius on a Fibonacci lattice, as
 * `vertex` cells: point `i` at height `z = 1 - (2 i + 1) / N` and azimuth
 * `i * pi * (3 - sqrt 5)`.
 * @param N the point count; a zero yields an empty mesh.
 */
MESHIOPLUSPLUS_API Mesh sphere_points(std::int64_t N = 100, double Radius = 1.0,
                                      const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                      std::int64_t MaxCells = 20000000);

/**
 * @brief @p N points drawn uniformly in the box `[origin, origin + size]`, as
 * `vertex` cells. The draws come from splitmix64 seeded with @p Seed, three per
 * point (x, y, z), each `(next() >> 11) * 2^-53`, so the cloud is the same on
 * every platform and in the Python twin.
 * @param N the point count; a zero yields an empty mesh.
 */
MESHIOPLUSPLUS_API Mesh box_points(std::int64_t N = 100,
                                   const std::array<double, 3>& rSize = {{1.0, 1.0, 1.0}},
                                   const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                                   std::uint64_t Seed = 0, std::int64_t MaxCells = 20000000);

// ---------------------------------------------------------------------------
// Surfaces.
// ---------------------------------------------------------------------------

/**
 * @brief A sphere surface from a regular tetrahedron split in four
 * @p Subdivisions times (`4 * 4^subdivisions` triangles), as `icosa_sphere`.
 */
MESHIOPLUSPLUS_API Mesh tetra_sphere(double Radius = 1.0, std::int64_t Subdivisions = 3,
                                     const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                     int Order = 1, std::int64_t MaxCells = 20000000,
                                     bool Complete = false);

/**
 * @brief A sphere surface from a regular octahedron split in four
 * @p Subdivisions times (`8 * 4^subdivisions` triangles), as `icosa_sphere`.
 */
MESHIOPLUSPLUS_API Mesh octa_sphere(double Radius = 1.0, std::int64_t Subdivisions = 2,
                                    const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                    int Order = 1, std::int64_t MaxCells = 20000000,
                                    bool Complete = false);

/**
 * @brief A sphere surface from the six faces of a cube, each an `N x N` grid,
 * pushed radially onto the sphere (`6 N^2` quads, or twice as many triangles).
 * Points are the cube's boundary lattice, x fastest, then y, then z; cells run
 * face by face (`xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`).
 * @param rCellType `"quad"` or `"triangle"` (each quad cut along its first
 *        diagonal).
 */
MESHIOPLUSPLUS_API Mesh cube_sphere(double Radius = 1.0, std::int64_t N = 4,
                                    const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                    const std::string& rCellType = "quad", int Order = 1,
                                    std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A geodesic sphere: every face of an icosahedron cut into an `N x N`
 * triangle lattice (frequency @p N, `20 N^2` triangles) pushed radially onto
 * the sphere. Unlike `icosa_sphere`, any frequency, not only powers of two.
 */
MESHIOPLUSPLUS_API Mesh geo_sphere(double Radius = 1.0, std::int64_t N = 4,
                                   const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                   int Order = 1, std::int64_t MaxCells = 20000000,
                                   bool Complete = false);

/**
 * @brief A torus around the z axis: tube centre circle of @p MajorRadius, tube
 * of @p MinorRadius; @p Nu cells around the axis and @p Nv around the tube,
 * periodic both ways and welded. Genus one: Euler characteristic zero.
 * @param rCellType `"quad"` or `"triangle"`.
 */
MESHIOPLUSPLUS_API Mesh torus(double MajorRadius = 1.0, double MinorRadius = 0.25,
                              std::int64_t Nu = 32, std::int64_t Nv = 16,
                              const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                              const std::string& rCellType = "quad", int Order = 1,
                              std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A Moebius band of centre-line @p Radius and @p Width with
 * @p NumTwists half-twists: @p Nu cells along it and @p Nv across. An odd count
 * has no consistent orientation (one boundary loop); an even count is an
 * orientable band (two loops).
 * @return a mesh with one side region, `boundary`.
 */
MESHIOPLUSPLUS_API Mesh moebius(double Radius = 1.0, double Width = 0.5, std::int64_t Nu = 32,
                                std::int64_t Nv = 4, std::int64_t NumTwists = 1,
                                const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                const std::string& rCellType = "quad", int Order = 1,
                                std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief An open cylinder surface around the z axis, centred on @p rCenter:
 * @p Ntheta cells around and @p Nz along, oriented outward.
 * @return a mesh with side regions `bottom` and `top` (its two rims).
 */
MESHIOPLUSPLUS_API Mesh tube(double Radius = 1.0, double Height = 1.0, std::int64_t Ntheta = 16,
                             std::int64_t Nz = 4,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             const std::string& rCellType = "quad", int Order = 1,
                             std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A closed cone surface: base disk of @p Radius at `z = center[2]`, apex
 * @p Height above it, @p Nz rings on the slant and on the base. The rim is a
 * sharp edge.
 * @param rCellType `"quad"` (triangle fans at the apex and the base centre) or
 *        `"triangle"`.
 * @return a mesh with cell regions `side` and `base`.
 */
MESHIOPLUSPLUS_API Mesh cone(double Radius = 1.0, double Height = 1.0, std::int64_t Ntheta = 16,
                             std::int64_t Nz = 4,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             const std::string& rCellType = "quad", int Order = 1,
                             std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A capsule: a cylinder of @p Radius and length @p Height closed by two
 * hemispheres, centred on @p rCenter along z. @p Nz rings along the cylinder
 * and @p Ncap on each cap; smooth (G1) where they meet.
 * @param rCellType `"quad"` (triangle fans at the poles) or `"triangle"`.
 */
MESHIOPLUSPLUS_API Mesh capsule(double Radius = 0.5, double Height = 1.0, std::int64_t Ntheta = 16,
                                std::int64_t Nz = 4, std::int64_t Ncap = 4,
                                const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                const std::string& rCellType = "quad", int Order = 1,
                                std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief The boundary of a box as a closed surface: each face an
 * `n_a x n_b` grid, oriented outward, sharing its edges with its neighbours.
 * Points and cells as `cube_sphere`.
 * @return a mesh with cell regions `xmin` ... `zmax`, one per face.
 */
MESHIOPLUSPLUS_API Mesh box_surface(const std::array<std::int64_t, 3>& rDims = {{1, 1, 1}},
                                    const std::array<double, 3>& rSize = {{1.0, 1.0, 1.0}},
                                    const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                                    const std::string& rCellType = "quad", int Order = 1,
                                    std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A platonic solid's surface with one cell per face, inscribed in the
 * sphere of @p Radius: `"tetrahedron"`, `"octahedron"` and `"icosahedron"` in
 * triangles, `"cube"` in quads and `"dodecahedron"` as a `polygon` block of
 * pentagons. Linear only.
 */
MESHIOPLUSPLUS_API Mesh platonic(const std::string& rName = "icosahedron", double Radius = 1.0,
                                 const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                 std::int64_t MaxCells = 20000000);

// ---------------------------------------------------------------------------
// Solids.
// ---------------------------------------------------------------------------

/**
 * @brief A ball of hexahedra by the seven-block butterfly: an inner cube of
 * half side `Radius / 2` and six blocks blending its faces into the sphere, so
 * no cell degenerates at the centre. The lattice is `(2N + 1)^3` points, x
 * fastest, `N` cells from the centre to the surface.
 * @param rCellType `"hexahedron"` or `"tetra"` (six Kuhn tetrahedra per cell).
 * @return a mesh with one side region, `surface`.
 */
MESHIOPLUSPLUS_API Mesh ball(double Radius = 1.0, std::int64_t N = 4,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             const std::string& rCellType = "hexahedron", int Order = 1,
                             std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A thick spherical shell of hexahedra between two radii: the
 * `cube_sphere` lattice of @p N cells per face edge on each of `Nr + 1` layers.
 * @param Ratio the thickness of each layer over the one inside it (1 is
 *        uniform; above 1 grows outward, below 1 refines toward the outer wall).
 * @return a mesh with side regions `inner` and `outer`.
 */
MESHIOPLUSPLUS_API Mesh spherical_shell(double InnerRadius = 0.5, double OuterRadius = 1.0,
                                        std::int64_t N = 4, std::int64_t Nr = 2, double Ratio = 1.0,
                                        const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                        int Order = 1, std::int64_t MaxCells = 20000000,
                                        bool Complete = false);

/**
 * @brief A thick-walled pipe of hexahedra around the z axis, centred on
 * @p rCenter: @p Ntheta cells around, @p Nr through the wall, @p Nz along.
 * @param Ratio the radial thickness of each layer over the one inside it.
 * @return a mesh with side regions `inner`, `outer`, `bottom` and `top`.
 */
MESHIOPLUSPLUS_API Mesh pipe(double InnerRadius = 0.5, double OuterRadius = 1.0,
                             double Height = 1.0, std::int64_t Ntheta = 16, std::int64_t Nr = 2,
                             std::int64_t Nz = 4, double Ratio = 1.0,
                             const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                             int Order = 1, std::int64_t MaxCells = 20000000,
                             bool Complete = false);

/**
 * @brief The Fichera corner: the cube `[-Size/2, Size/2]^3` around @p rCenter
 * minus its positive octant, the 3-D L-shape with a re-entrant vertex at the
 * centre. @p N cells per half side (`7 N^3` hexahedra).
 * @param rCellType `"hexahedron"` or `"tetra"`.
 * @return a mesh with side regions `xmin` ... `zmax` and the notch walls
 *         `notch_x`, `notch_y`, `notch_z` (on `x = 0`, `y = 0`, `z = 0`).
 */
MESHIOPLUSPLUS_API Mesh fichera(std::int64_t N = 2, double Size = 2.0,
                                const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                const std::string& rCellType = "hexahedron", int Order = 1,
                                std::int64_t MaxCells = 20000000, bool Complete = false);

/**
 * @brief A solid torus of hexahedra: the squircle disk of `disk` (@p N cells
 * across) as the tube's cross-section, swept around the z axis in @p Nu steps
 * and welded.
 * @return a mesh with one side region, `surface`.
 */
MESHIOPLUSPLUS_API Mesh solid_torus(double MajorRadius = 1.0, double MinorRadius = 0.25,
                                    std::int64_t N = 4, std::int64_t Nu = 32,
                                    const std::array<double, 3>& rCenter = {{0.0, 0.0, 0.0}},
                                    int Order = 1, std::int64_t MaxCells = 20000000,
                                    bool Complete = false);

}  // namespace meshioplusplus
