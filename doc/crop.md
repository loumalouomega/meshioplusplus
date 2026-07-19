# Crop (bounding box / half-space)

`meshioplusplus.crop(mesh, ...)` extracts the part of a mesh inside a region — an
axis-aligned **bounding box** or a **half-space** (a plane, keeping one side) —
pruning the now-unused points and remapping connectivity and all data. It is a
mesh **operation** (like [split](/split) and [merge](/merge)), not a file
format, and uses only standard C++/numpy, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("domain.vtu")

# axis-aligned bounding box
inside = meshioplusplus.crop(mesh, bbox=[0, 0, 0, 1, 1, 1])

# half-space: keep the side where (p - point) . normal >= 0
half = meshioplusplus.crop(mesh, plane=([0.5, 0, 0], [1, 0, 0]))

# keep a cell if ANY node is inside (default is ALL)
loose = meshioplusplus.crop(mesh, bbox=[0, 0, 0, 1, 1, 1], mode="any")

meshioplusplus.write("cropped.vtu", inside)
```

## Region and mode

| argument | meaning |
|---|---|
| `bbox` | `(xmin, ymin, zmin, xmax, ymax, zmax)` or `(lo3, hi3)` |
| `plane` | `(point3, normal3)` — keeps the side where `(p − point) · normal ≥ 0` |
| `mode` | `"all"` (default) keeps a cell only if **every** node is inside; `"any"` keeps it if **any** node is inside |
| `record_ids` | attach `crop:original_point_id` point_data and `crop:original_cell_id` cell_data of the source indices |

A point is inside the bbox when `lo ≤ p ≤ hi` component-wise. The point-inside
test is parallelized; exactly one of `bbox`/`plane` must be given.

## What changes

- Only the **kept cells** remain; **unused points are pruned** and connectivity
  is remapped so it stays valid.
- **`point_data`** / **`cell_data`** are subset to the kept points/cells;
  **`field_data`** is preserved.
- **`point_sets` / `cell_sets`** are remapped, dropping entries outside the kept
  region (done in the Python layer).

The numpy fallback handles rectangular cell blocks; ragged/polyhedron blocks are
handled by the C++ core only.

## CLI

```bash
meshioplusplus crop in.vtu out.vtu --bbox 0,0,0,1,1,1
meshioplusplus crop in.vtu out.vtu --plane 0.5,0,0,1,0,0 --mode any
meshioplusplus crop in.vtu out.vtu --bbox 0,0,0,1,1,1 --record-ids
```

Values that begin with `-` (e.g. a negative bbox corner) need the
`--bbox=-1,-1,-1,1,1,1` form so argparse doesn't read them as a flag. See the
[CLI reference](/cli).

## Other languages

- **C API** — `mio_crop_bbox(mesh, lo, hi, mode, record_ids)` and
  `mio_crop_plane(mesh, point, normal, mode, record_ids)` (`mode`: 0 = all, 1 =
  any). See the [C API reference](/c_api).
- **Fortran** — `mesh%crop_bbox(lo, hi, mode=...)` and `mesh%crop_plane(point,
  normal, mode=...)`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `cropBbox(mesh, lo, hi, mode, recordIds)` and
  `cropPlane(mesh, point, normal, mode, recordIds)`. See the
  [WebAssembly reference](/wasm).
