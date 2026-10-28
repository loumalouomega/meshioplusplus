# Parametric primitives

Almost every operation transforms a mesh you already have. [`grid`](voxelize.md#grid-as-a-primitive-constructor) was the one exception; these thirteen generators build a shape from parameters, so a test, a notebook, a demo or an MCP session needs no fixture file. [`extrude` and `revolve`](extrude.md) sweep the planar ones into solids with the same contract: extruding `rectangle` gives `box` exactly.

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
```

## The contract

Every generator keeps the contract `grid` already keeps:

- **Coordinates come from indices.** A point is computed from its own integers, never accumulated from its neighbour, so the last plane does not depend on the cell count in the low bits. An angle is `2πk/n`, not the previous point rotated.
- **A documented order.** Points and cells follow the order below, so a result can be indexed without a lookup.
- **A budget.** `max_cells` (default 20 000 000, non-positive lifts it) refuses an oversized request by name before allocating, as `grid` does.
- **A zero count is an empty mesh.** Points have shape `(0, 3)` and there is no cell block, so `NumCellBlocks() == 0` means there is nothing there. A negative count, a non-positive size or an unknown choice raises `ValueError` (`std::invalid_argument` in C++) naming the shape and the value.
- **Orientation.** Planar cells are counter-clockwise seen from +z, surfaces are oriented outward and volumes are positive.
- **Welded seams and poles.** A sphere's poles and the seam of a ring are single points, never duplicates, so [`periodic`](periodic.md), [`clean`](clean.md) and Euler-characteristic checks see no boundary that is not there.
- **Named boundaries.** The side regions in the table are [`Side` regions](regions.md): one `(cell, facet)` entry per boundary facet, `dim` 2 for a volume mesh and 1 for a planar one. Every boundary facet belongs to exactly one region.

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
- **`reference_cells`**: one block per requested type, in the order requested, offset by 2 (3 with `pairs`) along x so no two touch. Simplex and tensor cells span the unit simplex or cube, a wedge is the unit triangle times `[0, 1]`, a pyramid has a unit-square base and its apex at `(½, ½, 1)`, and a `polyhedron` is the unit cube as six outward faces (a `polyhedron8` block). Quadratic types add exact edge midpoints. `pairs=True` builds two cells that share a facet for `line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge` and `pyramid`.

## Quadratic cells and `order=2`

`order=2` returns serendipity cells (`quad8`, `triangle6`, `hexahedron20`, `tetra10`, `wedge15`) whose mid-edge nodes lie on the exact geometry. A sphere's edge midpoints are on the sphere, a cylinder's and a disk's rim edges are on the circle, and the seams and poles stay welded. These are the curved cells the notebooks' `curved()` helper draws, and a fixture for [`tessellate`](tessellation.md).

Corners keep their linear numbering. Mid-edge node `corners + e` of a cell is appended on first appearance, in cell-major order over the edge list of [`convert_cells`](convert_cells.md)'s `elevate` mode (`_ELEVATE` in Python, `cell_refine_edges` in C++), and shared by every cell on that edge. Each mid-node is evaluated from the shape's own map at the averaged parameters, never as a chord midpoint: a lattice coordinate for the grids, an angle for the rings (averaged along the shorter arc, with a pole or a disk centre taking the other end's azimuth), a unit vector re-normalized for the icosahedron. Full-Lagrange types (`quad9`, `hexahedron27`, `wedge18`) are not generated. The planar domains follow the same rule: `reference_triangle`, `ngon` and `lshape` have straight edges, so their map is the identity and the averaged position is the exact midpoint; `annulus` places them on the circles as `disk` does; `plate_with_hole` evaluates its blend at the averaged ring index and angular position (the angle wrapping around the hole), so the hole's mid-nodes are on the circle and the outer ones on the straight sides.

## Numerical parity of the twins

Every generator has a numpy twin that agrees with the C++ core, and `tests/python/test_primitives.py` pins it. `rectangle`, `box`, `icosa_sphere`, `reference_triangle`, `lshape`, `honeycomb` and `reference_cells` use index arithmetic and `sqrt` only, which IEEE 754 rounds correctly, so the twins agree to the byte. `disk`, `uv_sphere`, `cylinder`, `ngon`, `annulus` and `plate_with_hole` use `sin` and `cos`, which differ between glibc, MSVC, Emscripten's libc and numpy's SIMD kernels, so their topology and regions are exact and their coordinates agree to a few ulps.

## Surfaces

| Surface | Entry |
|---|---|
| Python | `mp.rectangle`, `mp.box`, `mp.disk`, `mp.uv_sphere`, `mp.icosa_sphere`, `mp.cylinder`, `mp.reference_triangle`, `mp.ngon`, `mp.annulus`, `mp.lshape`, `mp.plate_with_hole`, `mp.honeycomb`, `mp.reference_cells` (not `mp.triangle`, which is the Triangle format) |
| C++ | `operations/primitives.hpp` |
| [C](c_api.md) | `mio_rectangle` … `mio_honeycomb` over one `mio_gen_opts`, and `mio_reference_cells` |
| [Fortran](fortran.md) | `mio_rectangle` … `mio_reference_cells` |
| [Julia](julia.md) | `rectangle`, `box`, `disk`, `uv_sphere`, `icosa_sphere`, `cylinder`, `reference_triangle`, `ngon`, `annulus`, `lshape`, `plate_with_hole`, `honeycomb`, `reference_cells` |
| [R](r.md) | `mio_rectangle` … `mio_reference_cells` |
| [WASM](wasm.md) | `m.rectangle`, `m.box`, `m.disk`, `m.uvSphere`, `m.icosaSphere`, `m.cylinder`, `m.referenceTriangle`, `m.ngon`, `m.annulus`, `m.lshape`, `m.plateWithHole`, `m.honeycomb`, `m.referenceCells` |
| [CLI](cli.md#generate) | `meshioplusplus generate SHAPE OUTFILE`, in both CLIs |
| [MCP](mcp.md) | one tool per generator, none reading an input |
| [Pipeline](pipeline.md#generate) | `{ "Op": "Generate", "Shape": … }`; a document whose first step is `Generate` needs no `Input` |

The C options struct is shared by the twelve shapes and each reads only the fields its name suggests (`counts[0..2]`, `size`, `origin` as the lo corner or the centre, `radius`, `height`, `radius2` and `angle` for an annulus, `cell_type`, `split`, `order`, `max_cells`); `mio_gen_opts_init` fills it with counts of 8, unit sizes, an annulus of inner radius 0.5 over a full turn, and `order = 1`. `mio_annulus` reads `radius` as its outer radius, `mio_lshape`, `mio_plate_with_hole` and `mio_honeycomb` read the scalar `size[0]` as the bounding side, the plate's side and the hexagon edge, and `mio_plate_with_hole` needs `radius` below half of `size[0]`, so the init defaults of 1 and 1 must be changed. `radius2` and `angle` come from the struct's reserved tail, and a zero in either (a struct filled by an older caller) reads as the default.
