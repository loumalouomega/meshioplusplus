# Regular grids — `grid` and `voxelize`

`grid` builds a regular hexahedron lattice from nothing. `voxelize` builds one around a mesh and, optionally, keeps only the cells its surface passes through or encloses.

```python
import meshioplusplus as mp

background = mp.grid([32, 32, 32])                      # a lattice, from nothing
shell = mp.voxelize(mp.read("bunny.stl"), resolution=[64, 64, 64], fill="surface")
solid = mp.voxelize(mp.read("bunny.stl"), resolution=[64, 64, 64], fill="inside")
```

![voxelize's three fill modes around the bracket's surface: all, surface, inside](/images/voxelize_fills.png)

## The output is an ordinary mesh, and that is the whole design

A voxel grid here is one `hexahedron` cell block over a shared corner lattice — not a bespoke object with its own accessors. Once it *is* a [`Mesh`](mesh_data_model.md), everything already works on it:

```python
mp.write("grid.vtu", solid)                             # every writer
mp.write("grid.svg", solid, color_by="sdf:distance")    # colouring
mp.view(solid)                                          # both viewers
mp.crop(solid, bbox=...)                                # every operation
mp.isosurface(field, "sdf:distance", [0.0])             # contouring
```

None of that needed a line of new code, and none of it would have been available from a dedicated grid type. Two alternatives were considered and rejected:

- **`custom` cells** (for a future adaptive octree). `CellType::Custom` reports −1 nodes and −1 dimension, which makes the block invisible to `stats`, `quality`, `extract_surface`, `gradient` and `refine` and unwritable by most formats — forfeiting exactly the property that motivates the choice.
- **A new `voxel` cell type.** VTK's type 11 is deliberately unmapped in this codebase; adding it would mean a row in all 42 format tables to buy an implicit node ordering nothing needs.

## Numbering

Points run **x fastest, then y, then z**, and cells run in the same order:

```
vid(i, j, k) = (k * (ny + 1) + j) * (nx + 1) + i
cid(i, j, k) = (k *  ny      + j) *  nx      + i
```

with each cell's eight nodes in the meshio/VTK `hexahedron` winding. That ordering is *inherited*, not chosen: it is what this repository's existing C++/Python byte-identity fixtures already agree on, and it is VTK ImageData's own layout, so `cell_data` reshapes to `arr[z, y, x]` — the C-order tensor voxel tooling expects.

Two consequences worth stating as contracts:

- **Shared corners are deduplicated arithmetically.** Neighbouring cells reference the same node because the index formula says so — no hash, no tolerance, no welding pass. Conformity and determinism here are *structural*.
- **Every cell is a right parallelepiped**, so `attach_quality` reports a scaled Jacobian of exactly `1.0`. That is a free correctness check: a transposed axis shows up as an inverted cell rather than as a picture that looks slightly wrong, and the test suite uses it as one.

## Fill rules

| `fill` | keeps | needs |
|---|---|---|
| `all` (default) | every cell of the bounding box | nothing — the mesh contributes only its extent |
| `surface` | cells a surface triangle passes through | a surface; works on an open sheet |
| `inside` | cells whose centre is inside the surface | a surface closed enough for the chosen sign |

`surface` uses **exact triangle/box overlap** (the separating-axis theorem), not a bounding-box test — a long diagonal triangle overlaps far more boxes than it enters. A triangle lying exactly on a cell face marks both neighbouring cells; picking one side would need a tie-break with no geometric justification, and marking both keeps the occupied set a closed cover of the surface.

A selective fill prunes the points no kept cell references, numbering the survivors in ascending original order.

## Sizing

Give **exactly one** of `resolution` and `cell_size`. Defaulting one of them would silently pick a grid you did not choose, and the cost is cubic in that choice:

| grid | cells | points + connectivity |
|---|---|---|
| 128³ | 2.1 M | 0.19 GB |
| 256³ | 16.8 M | 1.48 GB |
| 512³ | 134 M | 11.8 GB |

`max_cells` (default 20 000 000, a little above 256³) refuses by name above that rather than letting the allocation fail. `cell_size` **covers** the box: each axis gets `ceil(extent / cell)` cells, so the lattice may extend slightly past the requested bounds rather than clipping them.

`padding` and `padding_relative` grow the box on every side — the latter as a fraction of the bounding-box diagonal.

## CLI

```bash
meshioplusplus voxelize bunny.stl shell.vtu --resolution 64,64,64 --fill surface
meshioplusplus voxelize bunny.stl solid.vtu --cell-size 0.5 --fill inside
meshioplusplus voxelize bunny.stl box.vtu   --resolution 32,32,32 --padding-relative 0.1
```

A negative bound needs the `--bounds=` form, as everywhere else in the CLI.

## Pipeline

```json
{ "Op": "Voxelize", "Resolution": [64, 64, 64], "Fill": "surface" }
```

Unlike every other step this one replaces its input's geometry rather than transforming it, which is exactly what makes it useful in a chain: read a skin, voxelize it, write a grid.

## `grid` as a primitive constructor

Every other operation transforms a mesh you already have; `grid` is the first that creates one. It takes cell counts, an origin and a spacing, and it is the same lattice `voxelize` builds:

```python
mp.grid([4, 4, 4], origin=(-1.0, -1.0, -1.0), spacing=(0.5, 0.5, 0.5))
```

A zero cell count on any axis yields an empty mesh — a legal, if useless, request rather than an error.

## Determinism

Output is byte-identical across the three mesh backends, across thread counts, and across the C++/numpy boundary (`tests/python/test_voxelize.py::test_cpp_matches_python`). The lattice is index arithmetic, and the fill rules are a totally ordered search; neither has an ordering question to answer.

## See also

- [Signed distance](sdf.md) — what `fill="inside"` uses, and the distance field itself.
- [Isosurfaces](isosurface.md) — contour a field on the grid.
