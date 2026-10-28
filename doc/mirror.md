# Reflection — `mirror`

`mirror` reflects a mesh across a plane and keeps every cell right side out (roadmap 3.2.5). A reflecting [`transform`](./transform.md) — a matrix of negative determinant — moves the points and leaves the connectivity as it was, so every cell comes out inverted: a tetrahedron of negative volume, a triangle whose normal points inward, and `Side` regions naming facets by their old local numbers. `mirror` reflects the points across the plane through `origin` with normal `normal`, permutes the nodes of every cell so its orientation is valid again, and moves each `Side` region entry to its facet's new number.

```python
import meshioplusplus as mp

half = mp.read("half_model.vtu")                       # the x >= 0 half of a symmetric part
whole, report = mp.mirror(half, keep_original=True, return_report=True)
print(report["num_points_on_plane"], report["num_dropped_sides"])
```

## How it works

Each point `x` goes to `x - 2 ((x - o) · n) n`. Each cell's nodes are permuted by a fixed orientation-reversing table (`detail/cell_flip`): every table keeps node 0, is its own inverse, and keeps mid-edge, mid-face and centre nodes in the slots of the edges and faces they belong to. Covered are `vertex`, `line` to `line11`, `triangle`, `triangle6`, `triangle7`, `quad`, `quad8`, `quad9`, `tetra`, `tetra10`, `pyramid`, `pyramid13`, `pyramid14`, `wedge`, `wedge15`, `wedge18`, `hexahedron`, `hexahedron20` and `hexahedron27`; a polygon's ring and each face ring of a polyhedron are reversed after their first node. Other types (`triangle10` and the higher Lagrange orders, the VTK Lagrange and Bézier families) are refused by name. The facet each local facet number moves to is found from the same `cell_faces` and `cell_edges` tables that number `Side` regions, so a region names the same nodes after the reflection as before.

## Keeping the original

With `keep_original`, the result is the input followed by its image: each block holds its own cells, then their reflected copies, and the points are the input's, then the copies'. With `weld` (the default) a point on the plane — within `tolerance` of the points' bounding-box diagonal — is its own image, so the two halves share it and a half model becomes a whole one, conforming across the plane ([`check_conformity`](./check_conformity.md) finds nothing on it). A cell lying entirely on the plane would be its own image too, so its copy is dropped (`num_dropped_cells`), and a `Side` facet on the plane is interior once both halves meet, so its entries are dropped (`num_dropped_sides`, a warning per region): the symmetry plane's own region stays, empty, as named groups do. Without `weld` the plane's points are duplicated and nothing is dropped.

`Cell` and `Point` regions name both copies; `Side` regions both copies of their facets.

## Data

Point and cell data are carried row by row, a copy's row its source's. With `reflect_vector_data`, float arrays of three columns are reflected as vectors (`v' = H v`) and of nine as row-major tensors (`A' = H A H`), with `H = I - 2 n nᵀ`, on the reflected copies only — the convention `transform`'s `rotate_vector_data` uses. A pseudo-vector (a vorticity, a normal stored as data) reflects with the opposite sign and is not singled out. Field data and property sets are copied.

Two-dimensional points stay two-dimensional: the plane must then be a line of the xy plane (a normal whose z component is 0).

## Options and report

| option | meaning |
|---|---|
| `origin` | a point of the plane (`(0, 0, 0)`) |
| `normal` | the plane's normal, any non-zero length (`(1, 0, 0)`) |
| `keep_original` | the input followed by its image (false) |
| `weld` | with `keep_original`, share the points on the plane (true) |
| `tolerance` | on the plane: a fraction of the points' bounding-box diagonal (1e-9) |
| `reflect_vector_data` | reflect float vector and tensor data (false) |

The report has `num_flipped_cells` (reflected cells written), `num_points_on_plane`, `num_dropped_cells`, `num_dropped_sides` and `warnings`.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.mirror(mesh, origin=(0, 0, 0), normal=(1, 0, 0), *, keep_original=False, weld=True, tolerance=1e-9, reflect_vector_data=False, return_report=False)`; no numpy twin |
| C++ | `mirror(const Mesh&, const MirrorOptions&)` → `MirrorResult` (`operations/mirror.hpp`) |
| C | `mio_mirror(mesh, opts, report)` with `mio_mirror_opts` (104 bytes, through `mio_mirror_opts_init`) and `mio_mirror_report` (72 bytes), `mio_mirror_warnings` |
| Fortran | `mio_mirror(mesh, origin, normal, keep_original, weld, tolerance, reflect_vector_data, report, stat, errmsg)` |
| Julia | `mirror(mesh; origin, normal, keep_original, weld, tolerance, reflect_vector_data)` |
| R | `mio_mirror(mesh, origin, normal, keep_original, weld, tolerance, reflect_vector_data)` |
| WASM | `mirror(mesh, {origin, normal, keepOriginal, weld, tolerance, reflectVectorData})` |
| CLI | `mirror INFILE OUTFILE [--origin=x,y,z] [--normal=x,y,z] [--keep-original] [--no-weld] [--tolerance T] [--reflect-data] [--json]`, both CLIs |
| MCP | `mirror(input_path, output_path, origin, normal, keep_original, weld, tolerance, reflect_vector_data)` |
| Pipeline | `Mirror` with `Point`, `Normal`, `KeepOriginal`, `Weld`, `Tolerance`, `ReflectData`; the step reports `NumFlippedCells`, `NumPointsOnPlane` and `NumDroppedCells` |
