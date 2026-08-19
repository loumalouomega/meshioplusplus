# Geometric statistics

`meshioplusplus.compute_stats(mesh)` reports **geometric** statistics of a mesh — the geometric complement to the topological [`info`](/cli#meshioplusplus-info) verb. It computes the bounding box + extents, centroid, per-cell-type counts, total surface area (of 2D cells plus the boundary of 3D cells), total signed & unsigned volume of 3D cells, and the count of inverted (negative signed-volume) cells — without modifying the mesh. It is a mesh **operation** (like [quality metrics](/mesh_quality)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![Cell-type counts reported by `compute_stats`](/images/stats_celltypes.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("part.vtu")
s = meshioplusplus.compute_stats(mesh)

print(s["bbox_min"], s["bbox_max"], s["extent"])
print(s["centroid"])
print(s["cell_type_counts"])          # {"tetra": 1234, ...}
print(s["total_area"], s["signed_volume"], s["unsigned_volume"])
print(s["num_inverted"])              # cells with negative volume
```

`compute_stats` returns a dict with keys `num_points`, `num_cells`, `bbox_min`/`bbox_max`/`extent`/`centroid` (3-tuples), `cell_type_counts` (`{type: count}`), `total_area`, `signed_volume`, `unsigned_volume`, and `num_inverted`.

## What is measured

- **Bounding box / extent / centroid** — over all points (a parallel reduction).
- **Per-cell-type counts** — the number of cells of each type.
- **Total area** — the area of 2D cells (triangles, quads) plus the boundary surface area of 3D cells (computed via [surface extraction](/extract_surface)).
- **Signed & unsigned volume** — of 3D cells (tetra / hex / wedge / pyramid), via the divergence theorem over each cell's outward-wound faces. `signed` sums the volumes; `unsigned` sums their magnitudes.
- **Inverted cells** — 3D cells with negative signed volume (wrong orientation).

Since v9.16.0 the volume comes from the shared polyhedral kernel, so:

- **`polyhedron` blocks are measured** rather than skipped. Their faces are oriented first (winding is repaired, not required — see [Polyhedra](/polyhedra)), so they never report as inverted; a cell that is not a closed orientable surface is excluded from the totals and warned about.
- **Ragged `polygon` blocks contribute their area**, which they previously did not.
- A cell with **non-planar faces** measures slightly differently than before v9.16.0: each face is fanned about its corner average rather than its first node, so the answer no longer depends on where the face's node list happens to start. Planar faces are unaffected.

## CLI

```bash
meshioplusplus stats mesh.vtu             # human-readable table
meshioplusplus stats mesh.vtu --json      # machine-readable JSON
```

See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_stats(mesh, &report)` fills a `mio_stats_report` struct (the scalar measures; per-cell-type counts are not carried across the flat ABI). See the [C API reference](/c_api).
- **Fortran** — `report = mesh%stats()` returns a `mio_stats_report` derived type. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `stats(mesh)` returns an object with `numPoints`, `bboxMin`/…, `cellTypeCounts`, `totalArea`, `signedVolume`, etc. See the [WebAssembly reference](/wasm).
