# Level-set cutting — `cut_level_set`, `cut_quadrature`

`cut_level_set` returns the part of a mesh that one or several level sets select: the pieces of the cells they cross, on one side of them or on their zero sets, and the cells wholly inside. `cut_quadrature` returns quadrature rules on that part, the cut-cell quadrature of unfitted finite elements (CutFEM). Both are [CutCells](https://github.com/sclaus2/CutCells)' (S. Claus, ONERA, MIT), an optional dependency of the C++ core: without it both refuse by name, pointing at `MESHIOPLUSPLUS_WITH_CUTCELLS`.

```python
import numpy as np
import meshioplusplus as mp

box = mp.box((8, 8, 8))  # hexahedra
box.point_data["ball"] = np.linalg.norm(box.points - 0.5, axis=1) - 0.3
inside, report = mp.cut_level_set(box, "ball", return_report=True)
print(report["measure"], report["num_cut_cells"])        # the ball's volume, its cut cells
skin = mp.cut_level_set(box, "ball", selection="ball = 0")  # its surface
points = mp.cut_quadrature(box, "ball", order=3)           # point_data quadrature:weight
mp.cut_has_cutcells()                                      # False without CutCells
```

## What it gives that the other cutting operations do not

- [`slice`](/slice) and [`isosurface`](/isosurface) give only the interface, by marching tetrahedra on simplices; higher-order cells are linearised first.
- [`crop`](/crop) keeps whole cells.
- MMG's [`discretize_levelset`](/mmg) remeshes linear tetrahedra and triangles conformally to the level set and drops the data.
- `cut_level_set` cuts every cell family (`line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge`, `pyramid`) into pieces of its own family (a hexahedron cut through two opposite faces stays one hexahedron), keeps the uncut cells as they are, combines several level sets with sharp edges where they meet, reads a quadratic level set on `triangle6` and `tetra10` cells, and carries the data.

## Level sets and the selection

- **Level sets.** Each is a scalar `point_data` array (`level_sets`, one value per node) or a plane (`planes`: `(origin, normal)`, `(x - origin) . normal` with the normal normalised, named `plane0`, `plane1`, ... after the arrays). A level set is linear in each cell's corners, and quadratic on a mesh whose cut cells are all `triangle6` or all `tetra10` (all six or ten nodes; CutCells' Basix node order is converted from VTK's).
- **Selection.** CutCells' expression over the level sets' names: `phi < 0`, `phi > 0`, `phi = 0` (the zero set: an interface of one dimension less), joined with `and` and `or`, such as `"wall < 0 and hole > 0"` or `"phi1 = 0 and phi2 < 0"` (the part of one interface inside another level set). Empty means `"<first> < 0"`. Names must start with a letter and hold letters, digits and `_`; an unknown name or a syntax error is refused naming the level sets.
- **Cells.** The cells of the highest topological dimension are cut; lower-dimensional blocks (a volume mesh's boundary triangles, say) are left out with a warning. Polygons and polyhedra are refused (`convert_cells(mode="simplexify")` first). Geometry is straight-sided: a curved higher-order cell is cut through its corners, and its other nodes only carry the quadratic level set and its data.

## Backends and output

- **`lut`** cuts each cell by lookup tables into straight pieces (`triangulate=True` splits them into simplices). It is exact for a level set that is linear in each cell: a planar cut's volume, area and interface come out to rounding.
- **`quadrays`** integrates curved interfaces to high order (Bernstein bounds certify the level set's sign on boxes, which are bisected where they cannot), and `cut_level_set` gives its pieces as Lagrange cells of `degree` (`VTK_LAGRANGE_HEXAHEDRON`, `..._QUADRILATERAL`, `..._CURVE`), the curved interface included. `"auto"`, the default, takes quadrays for a quadratic level set and the lookup tables otherwise.
- **The mesh.** Pieces are welded (CutCells gives each its own vertices), grouped into one block per type, and carry `cut:parent_cell` (the input cell, global block-major) and `cut:cut` (1 on a piece of a cut cell, 0 on a whole cell). `point_data` is interpolated in each new vertex's parent with the parent's shape functions (quadratic ones on `triangle6`/`tetra10` with a quadratic level set; integer arrays take the nearest node's value), `cell_data` and cell regions follow the parent, point and side regions are dropped with a warning, and `field_data` is kept.
- **The report.** `num_cells` (cells of the cut dimension), `num_cut_cells`, `num_uncut_cells` (wholly in the part), `num_pieces` (output cells), `dimension` of the part, `measure` (its volume, area or length, by a degree-1 rule of the same backend), `selection`, `backend` and `warnings`.

## Quadrature

`cut_quadrature(mesh, ..., order=k)` returns rules of degree `k` (exact for polynomials of degree `k` on straight pieces) on the part: one `vertex` per point at its physical position, with point data `quadrature:weight` (physical: the weights of the part sum to its measure), `quadrature:reference` (the point in the parent cell's reference element, `tdim` components), `quadrature:rule` (one rule per covered cell piece) and `quadrature:parent_cell`. Uncut cells inside the part get their own full-cell rules unless `include_uncut=False`. The report gives `num_rules`, `num_points`, `sum_weights`, `dimension`, `selection`, `backend` and `warnings`.

## Measured

On a unit cube cut by the plane `x = 0.3`, every family (hexahedra, tetrahedra, wedges, the cube as six pyramids, and quads and triangles in 2-D) gives the volume 0.3 on one side, 0.7 on the other and the interface area 1, all to 1e-12; two planes give the quarter, three quarters and half faces their selections ask for. A sphere of radius 0.3 on boxes of 4, 8 and 16 tetrahedra a side converges at second order (the error at least 2.5 times smaller each halving), and the same sphere as a quadratic level set on `tetra10` cells is at least five times more accurate than the linear one on the same mesh (`tests/cpp/test_cut.cpp`).

## Building with CutCells

CutCells is bring-your-own, never vendored or downloaded by CMake. `build/build-cutcells.sh PREFIX` builds the commit `build/deps-pins.sh` pins (CutCells has no release tags), static and without OpenMP; then configure meshio++ with `-DMESHIOPLUSPLUS_WITH_CUTCELLS=ON -DCUTCELLS_ROOT=PREFIX` (or `build/configure.sh --cutcells-root PREFIX`). The PyPI wheels for Linux and macOS link it (`build/build-permissive-libs.sh`); the Windows wheels, both WebAssembly artifacts and the default source build do not, and refuse by name. `cut_has_cutcells()` (C++, Python), `_core.__has_cutcells__`, `mio_has_cutcells()` (C) and their binding twins report it.

## Bindings

- **C++:** `meshioplusplus/operations/cut.hpp`: `cut_level_set(mesh, CutOptions)`, `cut_quadrature`, `cut_has_cutcells`.
- **Python:** `cut_level_set(mesh, level_sets, planes=, selection=, include_uncut=, backend=, degree=, triangulate=, record_parent_ids=, return_report=)` and `cut_quadrature(..., order=)`.
- **CLIs:** `meshioplusplus cut IN OUT --level-set NAME --plane OX,OY,OZ,NX,NY,NZ --selection EXPR [--cut-only] [--backend B] [--degree K] [--triangulate] [--json]` and `cut-quadrature ... --order K`, in both CLIs.
- **Pipeline:** the `CutLevelSet` step (`LevelSets`, `Planes` as six numbers per plane, `Selection`, `IncludeUncut`, `Backend`, `Degree`, `Triangulate`, `RecordParentIds`), in both engines.
- **MCP:** the `cut_level_set` and `cut_quadrature` tools (planes as six numbers each), reporting `cutcells_available`.
- **C:** `mio_cut_level_set` and `mio_cut_quadrature` with `mio_cut_opts` (104 bytes, `mio_cut_opts_init`) and `mio_cut_report` (104 bytes), and `mio_has_cutcells`; [Fortran](/fortran), [Julia](/julia), [R](/r) and [WASM](/wasm) (`cutLevelSet`, `cutQuadrature`, `hasCutcells`) wrap them.

## Limits

Straight-sided geometry only (CutCells' affine cells); polygons and polyhedra are refused; quadratic level sets only on meshes of `triangle6` or `tetra10` cells; at most 64 level sets in one cut (CutCells' mask); serial. The WebAssembly artifacts do not link CutCells yet.
