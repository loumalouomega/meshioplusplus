# Parametric primitives

Almost every operation transforms a mesh you already have. [`grid`](voxelize.md#grid-as-a-primitive-constructor) was the one exception; these thirty-four generators build a shape from parameters, so a test, a notebook, a demo or an MCP session needs no fixture file. [`extrude` and `revolve`](extrude.md) sweep the planar ones into solids with the same contract: extruding `rectangle` gives `box` exactly.

| Generator | Builds | Cell types | Side regions |
|---|---|---|---|
| `rectangle` | a regular 2-D grid at `z = origin[2]` | `quad`, `triangle` (`zigzag`, `up`, `down`, `center`) | `xmin`, `xmax`, `ymin`, `ymax` |
| `box` | a regular 3-D grid | `hexahedron`, `tetra` (six per cell), `wedge` (two per cell) | `xmin` … `zmax` |
| `disk` | a disk | `quad` (squircle map), `triangle` (hexagonal rings) | `rim` |
| `uv_sphere` | a latitude / longitude sphere surface | `triangle`, `quad` with triangle caps | none |
| `icosa_sphere` | an icosahedron-based sphere surface | `triangle` | none |
| `cylinder` | a solid cylinder along z | `hexahedron`, `wedge` | `side`, `top`, `bottom` |
| `reference_triangle` | a refined reference triangle at `z = origin[2]` | `triangle` | `ymin`, `xmin`, `hypotenuse` |
| `ngon` | a regular polygon, fanned from its centre or in rings | `triangle` | `rim` |
| `annulus` | an annulus, or a sector of one | `quad`, `triangle` | `inner`, `outer`, and `start`, `end` for a sector |
| `lshape` | the square minus its fourth quadrant | `quad`, `triangle` | `xmin`, `xmax`, `ymin`, `ymax`, `notch_x`, `notch_y` |
| `plate_with_hole` | a square plate around a circular hole | `quad`, `triangle` | `hole`, `xmin`, `xmax`, `ymin`, `ymax` |
| `honeycomb` | a tiling of pointy-top hexagons | one ragged `polygon` block | `boundary` |
| `reference_cells` | one cell of each type in its reference coordinates | every linear type, the quadratic ones, a polyhedron | none |
| `line` | a straight segment from `start` to `end` | `line` | point regions `start`, `end` |
| `circle` | a closed polyline on a circle | `line` | none |
| `helix` | a spring around the z axis | `line` | point regions `start`, `end` |
| `sphere_points` | points on a sphere (a Fibonacci lattice) | `vertex` | none |
| `box_points` | seeded uniform points in a box (splitmix64) | `vertex` | none |
| `tetra_sphere` | a subdivided tetrahedron pushed onto the sphere | `triangle` | none |
| `octa_sphere` | a subdivided octahedron pushed onto the sphere | `triangle` | none |
| `cube_sphere` | a cube's six gridded faces pushed onto the sphere | `quad`, `triangle` | none |
| `geo_sphere` | a geodesic sphere of any frequency | `triangle` | none |
| `torus` | a torus, periodic both ways | `quad`, `triangle` | none |
| `moebius` | a Moebius band with any number of half-twists | `quad`, `triangle` | `boundary` |
| `tube` | an open cylinder surface | `quad`, `triangle` | `bottom`, `top` |
| `cone` | a closed cone surface with a sharp rim | `quad` with triangle fans, `triangle` | cell regions `side`, `base` |
| `capsule` | a cylinder closed by two hemispheres | `quad` with triangle fans, `triangle` | none |
| `box_surface` | the closed boundary of a box | `quad`, `triangle` | cell regions `xmin` … `zmax` |
| `platonic` | a platonic solid, one cell per face | `triangle`, `quad`, or a `polygon` block of pentagons | none |
| `ball` | a solid ball by the seven-block butterfly | `hexahedron`, `tetra` | `surface` |
| `spherical_shell` | a thick shell between two radii, graded | `hexahedron` | `inner`, `outer` |
| `pipe` | a thick-walled pipe, graded | `hexahedron` | `inner`, `outer`, `bottom`, `top` |
| `fichera` | the cube minus its positive octant | `hexahedron`, `tetra` | `xmin` … `zmax`, `notch_x`, `notch_y`, `notch_z` |
| `solid_torus` | a solid torus | `hexahedron` | `surface` |

```python
import meshioplusplus as mp

plate = mp.rectangle(dims=(8, 4), size=(2.0, 1.0), cell_type="triangle", split="zigzag")
block = mp.box(dims=(4, 4, 4), cell_type="tetra", order=2)   # tetra10, six side regions
ball = mp.uv_sphere(radius=2.0, nu=24, nv=12, cell_type="quad", order=2)
pipe = mp.cylinder(radius=0.5, height=3.0, nxy=6, nz=12, cell_type="wedge")
cells = mp.reference_cells(["tetra", "hexahedron20", "polyhedron"])

# the planar domains
lshape = mp.lshape(n=8)                        # the re-entrant corner at the origin
ring = mp.annulus(0.5, 1.0, ntheta=32, nr=6, angle=3.14159)   # a sector with a seam pair
plate = mp.plate_with_hole(radius=0.25, size=1.0, nt=8, nr=6, order=2)
cells = mp.honeycomb(dims=(6, 5), edge=0.5)    # a ragged polygon block

# curves, point clouds, surfaces and solids
spring = mp.helix(radius=0.5, pitch=0.2, turns=5, n=200)
cloud = mp.box_points(n=10000, size=(1, 2, 3), seed=7)   # the same cloud everywhere
donut = mp.torus(major_radius=1.0, minor_radius=0.3, nu=48, nv=24)
band = mp.moebius(num_twists=1)                 # no consistent orientation
cone = mp.cone(ntheta=32, nz=4)                 # a sharp rim, cell regions side/base
shell = mp.spherical_shell(0.9, 1.0, n=8, nr=6, ratio=1.3)   # wall-resolved layers
corner = mp.fichera(n=4)                        # the 3-D re-entrant corner

# full-Lagrange and higher-order cells on the exact geometry
quad9 = mp.cube_sphere(n=4, order=2, complete=True)
cubic = mp.ball(n=2, order=3)                   # VTK_LAGRANGE_HEXAHEDRON
```

## The contract

Every generator keeps the contract `grid` already keeps:

- **Coordinates come from indices.** A point is computed from its own integers, never accumulated from its neighbour, so the last plane does not depend on the cell count in the low bits. An angle is `2πk/n`, not the previous point rotated.
- **A documented order.** Points and cells follow the order below, so a result can be indexed without a lookup.
- **A budget.** `max_cells` (default 20 000 000, non-positive lifts it) refuses an oversized request by name before allocating, as `grid` does.
- **A zero count is an empty mesh.** Points have shape `(0, 3)` and there is no cell block, so `NumCellBlocks() == 0` means there is nothing there. A negative count, a non-positive size or an unknown choice raises `ValueError` (`std::invalid_argument` in C++) naming the shape and the value.
- **Orientation.** Planar cells are counter-clockwise seen from +z, surfaces are oriented outward and volumes are positive.
- **Welded seams and poles.** A sphere's poles and the seam of a ring are single points, never duplicates, so [`periodic`](periodic.md), [`clean`](clean.md) and Euler-characteristic checks see no boundary that is not there.
- **Named boundaries.** The side regions in the table are [`Side` regions](regions.md): one `(cell, facet)` entry per boundary facet, `dim` 2 for a volume mesh and 1 for a planar one or an open surface's rim. Every boundary facet belongs to exactly one region. A closed surface has no boundary, so `cone` and `box_surface` name their patches as `Cell` regions instead (`dim` 2), and the open curves `line` and `helix` name their ends as `Point` regions.

## Point and cell order

- **`rectangle`**: points run x fastest, then y (`id = j * (nx + 1) + i`). Quad cells run in the same order with nodes `(i, j)`, `(i+1, j)`, `(i+1, j+1)`, `(i, j+1)`. Triangles split each quad: `up` cuts from the lo corner (`a b c`, `a c d`), `down` from the other (`a b d`, `b c d`), `zigzag` alternates by the parity of `i + j` starting with `up`, and `center` adds one node per cell after the lattice (`m = (nx + 1)(ny + 1) + cell`) and emits `a b m`, `b c m`, `c d m`, `d a m`.
- **`box`**: points run x fastest, then y, then z, and hexahedra are `grid`'s cells, so `box(dims, size)` and `grid(dims, spacing=size / dims)` agree. A tetra box emits the six Kuhn tetrahedra of each hexahedron around its main diagonal `0-6` (`0126`, `0236`, `0376`, `0746`, `0456`, `0516`); every one contains the diagonal, so neighbouring cells share whole faces and the mesh conforms. A wedge box emits two prisms per hexahedron, `(0 1 2 | 4 5 6)` and `(0 2 3 | 4 6 7)`.
- **`disk`, quads**: an `(n+1)²` lattice, x fastest, on `s, t ∈ [-1, 1]` (`s = (2i - n) / n`) pushed onto the disk by `x = r s √(1 - t²/2)`, `y = r t √(1 - s²/2)`. Every lattice-boundary point lands exactly on the circle. The four corner cells have an interior angle of 180° at the rim node, the price of a map with no seam.
- **`disk`, triangles**: a centre point, then ring `k = 1 … n` of `6k` points at radius `rk/n` and angle `2πm/6k`. Cells run ring by ring and sector by sector, each sector emitting an outward triangle then, between neighbours, an inward one: `6n²` triangles in all.
- **`uv_sphere`**: point 0 is the north pole `(0, 0, +r)`, then `nv - 1` rings of `nu` points from north to south (`1 + (j-1) nu + i`), then the south pole. Poles are exact (`x = y = 0`). Triangle spheres emit the north cap, then the body (two triangles per quad, `j` outer, `i` inner), then the south cap; `cell_type="quad"` gives a `quad` block for the body followed by a `triangle` block for both caps.
- **`icosa_sphere`**: the 12 icosahedron vertices, then each level's edge midpoints in first-appearance order over the faces, normalized onto the sphere. Each face `(v0, v1, v2)` becomes `(v0, m01, m20)`, `(v1, m12, m01)`, `(v2, m20, m12)`, `(m01, m12, m20)`. There are `20 · 4^s` triangles and `10 · 4^s + 2` points.
- **`cylinder`**: layer `k = 0 … nz` of the disk's points at height `h(2k - nz) / 2nz`, layer by layer; cells run layer by layer in the disk's cell order, the wedge cylinder extruding the triangle disk.
- **`reference_triangle`**: nodes run in barycentric lattice order: row `j` (y) holds the nodes `i = 0 … n - j` (x fastest), `id(i, j) = j (n + 1) - j (j - 1) / 2 + i`, at `origin + size * (i, j) / n`. Cells run row by row: for each `i` the lattice triangle `(i, j) (i+1, j) (i, j+1)`, then the inverted one `(i+1, j) (i+1, j+1) (i, j+1)` to its right (absent in the last column), `n²` in all.
- **`ngon`**: the centre node, then ring `k = 1 … rings` of `sides · k` nodes, sector by sector: node `t` of sector `s` is `((k - t) V_s + t V_{s+1}) / rings`, with `V_s = radius (cos 2πs/sides, sin 2πs/sides)`, so every ring lies on the straight chords between the vertices. The triangles are those of the triangle `disk` with `sides` sectors instead of six, `sides · rings²` in all; `rings = 1` is the plain fan.
- **`annulus`**: ring `k = 0 … nr` runs from the inner to the outer radius (`inner + (outer - inner) k / nr`, exactly `outer` at the last), each ring `ntheta` nodes (`ntheta + 1` for a sector) at the angle `angle · i / ntheta`, counter-clockwise from +x; node id `k · np + i`. Cells run ring by ring: quad `(k, i) (k+1, i) (k+1, i+1) (k, i+1)`, or its triangles `(a b c)` and `(a c d)`. A full ring (`angle = 2π`) welds node `ntheta` onto node 0; a sector keeps the two end rays as distinct node chains, so [`periodic`](periodic.md) can link `start` to `end` through a rotation by `angle`.
- **`lshape`**: the `(2n + 1)²` lattice of `[-size/2, size/2]²` at `center + size (i - n, j - n) / 2n`, x fastest, without the nodes strictly inside the removed quadrant `x > 0, y < 0`; `3 n²` cells, x fastest too, cut `(a b c)`, `(a c d)` for triangles. The re-entrant corner is the node at `center`. `notch_x` is the notch wall on `x = center[0]` and `notch_y` the one on `y = center[1]`.
- **`plate_with_hole`**: `nr + 1` rings of `4 nt` nodes, node id `k · 4nt + p`. Position `p` is at the angle `-π/4 + 2π p / 4nt` on the hole circle and, on the outer square, at the point of the east side `(h, h (2m/nt - 1))` (`h = size / 2`, `m = p mod nt`) turned by `p div nt` quarter turns; node `(k, p)` blends the two linearly in `k / nr`, so the outer square is exact and the four blocks (east, north, west, south) share their seam lines. Cells run ring by ring, `(k, p) (k+1, p) (k+1, p+1) (k, p+1)`, `4 nt nr` in all.
- **`honeycomb`**: pointy-top hexagons on the vertex lattice `x = (√3/2) edge · u`, `y = (edge/2) · v`; hexagon `(i, j)` is centred at `(2i + 1 + (j mod 2), 3j + 2)`, rows run x fastest and odd rows are shifted right by half a hexagon. Each hexagon lists its corners counter-clockwise from the bottom one, and nodes are numbered by first appearance over the cells. The block is a `polygon` block of rows of six, ragged in the core, so it exercises the ragged-polygon paths ([polyhedra and polygons](polyhedra.md)). The `boundary` region holds the hexagon edges seen once (facet `f` of a polygon is its edge from node `f` to node `f + 1`). `order=2` is refused by name: a hexagon has no quadratic form.
- **`line`**: node `i = 0 … n` at `start + (end - start) i / n` (exactly `end` at the last), cells `(i, i + 1)`; point regions `start` (node 0) and `end` (node `n`).
- **`circle`**: node `i = 0 … n - 1` at the angle `2π i / n`, counter-clockwise from +x; cells `(i, i + 1 mod n)`, so the seam is welded.
- **`helix`**: node `i = 0 … n` at the angle `2π turns · i / n` and height `pitch · turns · i / n` above `center`; cells `(i, i + 1)`. No cell may span half a turn (`n > 2 turns`), so a mid-edge node takes the short way round.
- **`sphere_points`**: point `i = 0 … n - 1` at `z = 1 - (2i + 1) / n`, radius `√(1 - z²)` and azimuth `i · π (3 - √5)` (the golden angle), as one `vertex` cell each.
- **`box_points`**: three draws per point (x, y, z) from splitmix64 seeded with `seed`, each `(next() >> 11) · 2⁻⁵³` in `[0, 1)`, scaled into `[origin, origin + size]`. `std::uniform_real_distribution` is implementation-defined and numpy's generators have no C++ counterpart, so both twins carry the same small generator: the cloud is the same on every platform and in every binding.
- **`tetra_sphere`, `octa_sphere`**: as `icosa_sphere` from the regular tetrahedron `(±a, ±a, ±a)` with an even number of minus signs (`a = 1/√3`) or the octahedron `±x, ±y, ±z`, giving `4 · 4^s` and `8 · 4^s` triangles.
- **`cube_sphere`, `box_surface`**: the nodes of the `nx × ny × nz` lattice that lie on the box boundary, in lattice order (x fastest, then y, then z, interior skipped). Cells run face by face (`xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`), on each face the second in-plane index outer (`(j, k)` on the x faces, `(i, k)` on the y faces, `(i, j)` on the z faces), oriented outward; triangles cut each quad `(a b c)`, `(a c d)`. `box_surface` places the lattice on `origin + size · (i, j, k) / dims` and names one cell region per face; `cube_sphere` takes the `n³` lattice of `[-1, 1]³` and pushes every point (and every mid-edge node, averaged on the flat face) radially onto the sphere.
- **`geo_sphere`**: the 12 icosahedron vertices, then `n - 1` nodes on each edge (edges in first-appearance order over the faces' `v0 v1`, `v1 v2`, `v2 v0`, counted from the lower vertex), then each face's interior nodes `(n - i - j, i, j)` over its corners (`j` outer), all pushed onto the sphere: `10 n² + 2` points and `20 n²` triangles, any frequency. Each face is cut as `reference_triangle` cuts the unit triangle. (`uv_sphere` with `cell_type="quad"` already closes its poles with triangle caps; `geo_sphere` is the geodesic sphere.)
- **`torus`**: node `(i, j) = i · nv + j` at the angle `2π i / nu` around the axis and `2π j / nv` around the tube, `x = (R + r cos v) cos u`, `y = (R + r cos v) sin u`, `z = r sin v`; cells `(i, j) (i+1, j) (i+1, j+1) (i, j+1)` modulo both counts, oriented outward.
- **`moebius`**: node `(i, j) = i (nv + 1) + j` at the angle `u = 2π i / nu` and across-band parameter `w = (2j - nv) / nv`, `x = (R + w W/2 cos(t u / 2)) cos u`, `z = w W/2 sin(t u / 2)` for `t` half-twists. The last row of cells meets the first ring upside down when `t` is odd (node `j` there is node `nv - j`), so there is no consistent orientation: [`compute_normals`](normals.md) reports the `inconsistent_pairs` and `repair(fix_orientation=True)` cannot clear them. The `boundary` region holds both rims (one loop for odd `t`, two for even).
- **`tube`**: node `(k, i) = k · ntheta + i` on ring `k = 0 … nz` at height `h (2k - nz) / 2nz` and angle `2π i / ntheta`; cells `(k, i) (k, i+1) (k+1, i+1) (k+1, i)`, outward. Side regions `bottom` and `top` are the two rims.
- **`cone`**: node 0 is the apex, `height` above the base centre `center`; then `2 nz - 1` rings of `ntheta` nodes: down the slant to the rim (ring `nz`, radius `radius`) and in across the base; then the base centre. Cells as `uv_sphere`'s between its poles: with `quad`, a `quad` block for the bands and a `triangle` block for the fans at the apex and the centre. The rim is a sharp edge for [`feature_edges`](feature_edges.md); cell regions `side` and `base`.
- **`capsule`**: the north pole, `2 ncap + nz - 1` rings, the south pole, with the profile parameterised by arc length over the radius (a quarter circle, the straight side of length `height`, a quarter circle), so mid-edge nodes land on the caps and the side alike. Laid out as `cone`.
- **`platonic`**: the vertices on the sphere of `radius` (the cube's at `(±a, ±a, ±a)`, the dodecahedron's at those plus `(0, ±a/φ, ±aφ)` and its cyclic permutations, the icosahedron's as `icosa_sphere`'s), one cell per face oriented outward; the dodecahedron is a ragged `polygon` block of pentagons, as `honeycomb`'s. Linear only.
- **`ball`**: the `(2n + 1)³` lattice of `[-1, 1]³`, x fastest, as `box`'s cells (or its Kuhn tetrahedra). A lattice point `p` of max-norm `s` stays at `radius · p` inside the inner cube (`s ≤ ½`); outside it moves along its ray from the centre, blending the inner cube's face (`radius · p / 2s`) into the sphere (`radius · p / |p|`) linearly in `2s - 1`. That is the seven-block butterfly: no cell degenerates at the centre, and the `surface` region is the sphere exactly.
- **`spherical_shell`**: layer `k = 0 … nr` of `cube_sphere`'s boundary lattice; cells run layer by layer, face by face, each hexahedron the outward surface quad on layer `k` then on layer `k + 1`. The radial parameter `t = k / nr` maps to `inner + (outer - inner) (q^(nr t) - 1) / (q^nr - 1)` for `ratio = q` (`t` itself for `q = 1`), so each layer is `q` times as thick as the one inside it.
- **`pipe`**: node `(k, m, i) = (k (nr + 1) + m) ntheta + i` at height `h (2k - nz) / 2nz`, radial layer `m` (graded as the shell) and angle `2π i / ntheta`; hexahedra `(k, m, i) (k, m+1, i) (k, m+1, i+1) (k, m, i+1)` and the same on layer `k + 1`.
- **`fichera`**: the `(2n + 1)³` lattice of `[-size/2, size/2]³` around `center`, x fastest, without the nodes strictly inside the removed octant `x, y, z > 0`; `7 n³` cells, x fastest. `notch_x`, `notch_y` and `notch_z` are the three re-entrant walls on `x = 0`, `y = 0` and `z = 0`.
- **`solid_torus`**: slice `i = 0 … nu - 1` at the angle `2π i / nu` holds `disk`'s `(n + 1)²` squircle lattice as the cross-section, node `i (n + 1)² + jj (n + 1) + ii`; hexahedra join slice `i` to slice `i + 1 mod nu`, the cross-section quad clockwise in `(s, t)` so the swept volume is positive.
- **`reference_cells`**: one block per requested type, in the order requested, offset by 2 (3 with `pairs`) along x so no two touch. Simplex and tensor cells span the unit simplex or cube, a wedge is the unit triangle times `[0, 1]`, a pyramid has a unit-square base and its apex at `(½, ½, 1)`, and a `polyhedron` is the unit cube as six outward faces (a `polyhedron8` block). Quadratic types add exact edge midpoints. `pairs=True` builds two cells that share a facet for `line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge` and `pyramid`.

## Quadratic and higher-order cells

`order=2` returns serendipity cells (`quad8`, `triangle6`, `hexahedron20`, `tetra10`, `wedge15`) whose mid-edge nodes lie on the exact geometry. A sphere's edge midpoints are on the sphere, a cylinder's and a disk's rim edges are on the circle, and the seams and poles stay welded. These are the curved cells the notebooks' `curved()` helper draws, and a fixture for [`tessellate`](tessellation.md).

Corners keep their linear numbering. Mid-edge node `corners + e` of a cell is appended on first appearance, in cell-major order over the edge list of [`convert_cells`](convert_cells.md)'s `elevate` mode (`_ELEVATE` in Python, `cell_refine_edges` in C++), and shared by every cell on that edge. Each mid-node is evaluated from the shape's own map at the averaged parameters, never as a chord midpoint: a lattice coordinate for the grids, an angle for the rings (averaged along the shorter arc, with a pole or a disk centre taking the other end's azimuth), a unit vector re-normalized for the icosahedron.

`complete=True` with `order=2` returns the full-Lagrange types instead: `quad9`, `hexahedron27` and `wedge18` (`line3`, `triangle6` and `tetra10` already are complete, so those blocks do not change). Each face and body centre is the average of its corners in the same parameter space, so a `cube_sphere` quad's centre is on the sphere too. The centres are appended after every mid-edge node, in first-appearance order over the cells and each cell's centres in local node order (`cell_faces`' `hexahedron27` and `wedge18` order), and shared by every cell on the face.

`order` 3 to 10 returns the `VTK_LAGRANGE_*` family (`VTK_LAGRANGE_CURVE`, `_TRIANGLE`, `_QUADRILATERAL`, `_TETRAHEDRON`, `_HEXAHEDRON`, `_WEDGE`) at that order. Each cell's nodes are the equispaced lattice of its reference cell in VTK's order (`vtk_lattice`, VTK 9.1's edge order on the hexahedron), each named by the corners it combines and their integer weights (`lattice_weights`), so two cells on one edge or face share its nodes; a new node is appended on first appearance over the cells, at the same weighted combination of the corners' parameters. The side regions are found on the linear cells before elevation and keep their facet numbers, the corner facets of the Lagrange cell. A shape with a `polygon` or `vertex` block (`honeycomb`, `platonic`, the point clouds) has no higher order and refuses `order` above 1 by name. The planar domains follow the same rule: `reference_triangle`, `ngon` and `lshape` have straight edges, so their map is the identity and the averaged position is the exact midpoint; `annulus` places them on the circles as `disk` does; `plate_with_hole` evaluates its blend at the averaged ring index and angular position (the angle wrapping around the hole), so the hole's mid-nodes are on the circle and the outer ones on the straight sides.

## Numerical parity of the twins

Every generator has a numpy twin that agrees with the C++ core, and `tests/python/test_primitives.py` and `test_primitives_shapes.py` pin it. `rectangle`, `box`, `icosa_sphere`, `reference_triangle`, `lshape`, `honeycomb`, `reference_cells`, `line`, `box_points`, `tetra_sphere`, `octa_sphere`, `cube_sphere`, `geo_sphere`, `box_surface`, `platonic`, `ball` and `fichera` (and `spherical_shell` with `ratio=1`) use index arithmetic and `sqrt` only, which IEEE 754 rounds correctly, so the twins agree to the byte. The rest use `sin`, `cos` or `pow`, which differ between glibc, MSVC, Emscripten's libc and numpy's SIMD kernels, so their topology and regions are exact and their coordinates agree to a few ulps. The higher orders keep the same split: a Lagrange node is the same weighted sum in both twins, accumulated in the same corner order.

## Test helpers left hand-rolled

Most of the hand-written grid helpers of the test suites now call `rectangle` and `box`. These stay as they are, each for a recorded reason: `_unit_cube` (three Python copies and `test_voxelize.py`) numbers its corners in VTK order, not `box`'s x-fastest order; `_hex_cube` and `_tet_cube` carry extra `point_data` and `cell_data`, and `np.linspace` differs from `size * i / n` in the last bit above three cells; the i-major helpers of `test_smooth` and `test_openfoam` pin a different numbering; `test_region_analysis.cpp`'s `quad_grid` is asserted against region counts the generators' side regions would change; and the jittered `ft_*` fixtures of `golden_fixtures.hpp` feed pinned digests. The notebook helpers (`hex_grid` in 03 and 09, `hex_block` in 25, `plate` in 09, the squircle `disk` of 27) move to the generators at the next edit of each notebook, since a change re-executes it. `reference_cells` could replace the test-local canonical mesh of the [format conformance matrix](conformance.md), whose coordinates are exact in binary but which keeps a deliberately inexact `0.1` vertex, so the swap needs that vertex kept.

## Not generated

Generators driven by a user callback (no callback crosses the flat ABI; numpy users build `Mesh(points, cells)` directly, and [`isosurface`](isosurface.md) on a `grid` covers implicit shapes); scanned or modelled assets (the bunny, the teapot, the bracket fixture: data rather than parameters, each under its own licence); and CSG combinations (a [non-goal](roadmap.md#non-goals-and-decisions-taken)): two-body interface and contact fixtures are primitives placed with `transform` and joined with `merge`. The catalogue follows [meshzoo](https://github.com/meshpro/meshzoo)'s names and parameterisations where the two overlap, taken from its README only.

## Surfaces

| Surface | Entry |
|---|---|
| Python | `mp.rectangle`, `mp.box`, `mp.disk`, `mp.uv_sphere`, `mp.icosa_sphere`, `mp.cylinder`, `mp.reference_triangle`, `mp.ngon`, `mp.annulus`, `mp.lshape`, `mp.plate_with_hole`, `mp.honeycomb`, `mp.reference_cells`, `mp.line`, `mp.circle`, `mp.helix`, `mp.sphere_points`, `mp.box_points`, `mp.tetra_sphere`, `mp.octa_sphere`, `mp.cube_sphere`, `mp.geo_sphere`, `mp.torus`, `mp.moebius`, `mp.tube`, `mp.cone`, `mp.capsule`, `mp.box_surface`, `mp.platonic`, `mp.ball`, `mp.spherical_shell`, `mp.pipe`, `mp.fichera`, `mp.solid_torus` (not `mp.triangle`, which is the Triangle format) |
| C++ | `operations/primitives.hpp` |
| [C](c_api.md) | `mio_rectangle` … `mio_solid_torus` over one `mio_gen_opts`, and `mio_reference_cells` |
| [Fortran](fortran.md) | `mio_rectangle` … `mio_reference_cells` |
| [Julia](julia.md) | the Python names (`line` takes `stop` for its end point, `end` being reserved) |
| [R](r.md) | `mio_rectangle` … `mio_reference_cells` |
| [WASM](wasm.md) | the Python names in camelCase: `m.rectangle`, `m.uvSphere`, `m.plateWithHole`, `m.spherePoints`, `m.boxSurface`, `m.solidTorus` … |
| [CLI](cli.md#generate) | `meshioplusplus generate SHAPE OUTFILE`, in both CLIs |
| [MCP](mcp.md) | one tool per generator, none reading an input |
| [Pipeline](pipeline.md#generate) | `{ "Op": "Generate", "Shape": … }`; a document whose first step is `Generate` needs no `Input` |

The C options struct is shared by the thirty-three shapes and each reads only the fields its name suggests (`counts[0..2]`, `size`, `origin` as the lo corner or the centre, `radius`, `height`, `radius2`, `angle`, `num_twists`, `ratio`, `seed`, `cell_type`, `split`, `order`, `complete`, `max_cells`); `mio_gen_opts_init` fills it with counts of 8, unit sizes, an annulus of inner radius 0.5 over a full turn, one half-twist, a uniform wall and `order = 1`. `mio_annulus`, `mio_spherical_shell` and `mio_pipe` read `radius` as the outer radius and `radius2` as the inner; `mio_torus` and `mio_solid_torus` read `radius2` as the minor radius and `mio_moebius` as the band's width; `mio_helix` reads `height` as the pitch and `angle` as the total sweep in radians; `mio_line` reads `origin` as its start and `size` as its end point; `mio_platonic` reads `counts[0]` as the face count (4, 6, 8, 12 or 20); `mio_lshape`, `mio_plate_with_hole`, `mio_honeycomb` and `mio_fichera` read the scalar `size[0]` as the bounding side, the plate's side, the hexagon edge and the cube's side, and `mio_plate_with_hole` needs `radius` below half of `size[0]`, so the init defaults of 1 and 1 must be changed. `radius2`, `angle`, `num_twists`, `complete`, `ratio` and `seed` come from the struct's reserved tail; a zero `radius2`, `angle` or `ratio` (a struct filled by an older caller) reads as the shape's default.
