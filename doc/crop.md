# Crop (bounding box / half-space / data predicate)

`meshioplusplus.crop(mesh, ...)` extracts part of a mesh — the part inside an axis-aligned **bounding box**, inside a **half-space** (a plane, keeping one side), or the cells whose value in one of the mesh's own **`cell_data` arrays** satisfies a comparison — pruning the now-unused points and remapping connectivity and all data. It is a mesh **operation** (like [split](/split) and [merge](/merge)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![Half of a surface mesh cropped out by a bounding box](/images/crop_bbox.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("domain.vtu")

# axis-aligned bounding box
inside = meshioplusplus.crop(mesh, bbox=[0, 0, 0, 1, 1, 1])

# half-space: keep the side where (p - point) . normal >= 0
half = meshioplusplus.crop(mesh, plane=([0.5, 0, 0], [1, 0, 0]))

# keep a cell if ANY node is inside (default is ALL)
loose = meshioplusplus.crop(mesh, bbox=[0, 0, 0, 1, 1, 1], mode="any")

# a predicate on the mesh's own cell data
bad = meshioplusplus.crop(mesh, where=("quality:scaled_jacobian", "<", 0.3))

meshioplusplus.write("cropped.vtu", inside)
```

## Region and mode

| argument | meaning |
|---|---|
| `bbox` | `(xmin, ymin, zmin, xmax, ymax, zmax)` or `(lo3, hi3)` |
| `plane` | `(point3, normal3)` — keeps the side where `(p − point) · normal ≥ 0` |
| `where` | `(array_name, comparison, value)` — keeps cells whose value in a **scalar `cell_data`** array satisfies the comparison, one of `<`, `<=`, `>`, `>=`, `==`, `!=` |
| `mode` | `"all"` (default) keeps a cell only if **every** node is inside; `"any"` keeps it if **any** node is inside. `bbox`/`plane` only — see below |
| `record_ids` | attach `crop:original_point_id` point_data and `crop:original_cell_id` cell_data of the source indices |

A point is inside the bbox when `lo ≤ p ≤ hi` component-wise. The point-inside test is parallelized; exactly one of `bbox`/`plane`/`where` must be given.

## The predicate crop

`where=` is deliberately **general rather than inside/outside-a-surface specific**. Inside/outside then composes:

```python
f = meshioplusplus.distance_to_surface(mesh, skin, location="center")
inside = meshioplusplus.crop(f, where=("sdf:distance", "<", 0.0))
```

![crop(where=('sdf:distance', '<', 0.0)) keeping the inside of the bracket](/images/crop_where_sdf.png)

and the same one mode also crops by `quality:*`, by a material id, by `partition:part`, or by anything [`data calc`](/data_operations) can produce. A dedicated crop-by-surface would have served exactly one of those.

The comparison is [`refine`](/refine)'s vocabulary, evaluated by the same function, so the two operations cannot drift on the boundary cases. In particular **a non-finite cell value never matches**, whatever the comparison — including `!=`, where IEEE says `NaN != 1.0` is true. `compute_quality` deliberately reports NaN where a metric does not apply, and predicating over `quality:*` on a mixed mesh is the headline use case.

::: warning `mode` does not apply
`bbox` and `plane` test **points**, and then need a rule for reducing a cell's several nodes to one verdict. A `cell_data` predicate is already one value per cell and has nothing to reduce, so passing `mode=` alongside `where=` is an **error** rather than being silently ignored.
:::

`point_data` is likewise refused **by name** rather than averaged onto the cells: [`data to-cell`](/data_operations) is the explicit way to move it, and doing it implicitly would make the kept set depend on an averaging rule nobody asked for.

## What changes

- Only the **kept cells** remain; **unused points are pruned** and connectivity is remapped so it stays valid.
- **`point_data`** / **`cell_data`** are subset to the kept points/cells; **`field_data`** is preserved.
- **`point_sets` / `cell_sets`** are remapped, dropping entries outside the kept region (done in the Python layer).

The numpy fallback handles rectangular cell blocks; ragged/polyhedron blocks are handled by the C++ core only.

## CLI

```bash
meshioplusplus crop in.vtu out.vtu --bbox 0,0,0,1,1,1
meshioplusplus crop in.vtu out.vtu --plane 0.5,0,0,1,0,0 --mode any
meshioplusplus crop in.vtu out.vtu --bbox 0,0,0,1,1,1 --record-ids
meshioplusplus crop field.vtu inside.vtu --where "sdf:distance < 0"
```

`--where` takes `NAME OP VALUE` as one argument, scanned longest-operator-first so `<=` is not read as `<` plus a stray `=`. Array names routinely contain a colon (`sdf:distance`, `quality:scaled_jacobian`) and never an operator, so the split is unambiguous.

Values that begin with `-` (e.g. a negative bbox corner) need the `--bbox=-1,-1,-1,1,1,1` form so argparse doesn't read them as a flag. See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_crop_bbox(mesh, lo, hi, mode, record_ids)`, `mio_crop_plane(mesh, point, normal, mode, record_ids)` (`mode`: 0 = all, 1 = any) and `mio_crop_predicate(mesh, array, compare, value, record_ids)` (`compare`: a `mio_refine_compare`). See the [C API reference](/c_api).
- **Fortran** — `mesh%crop_bbox(lo, hi, mode=...)`, `mesh%crop_plane(point, normal, mode=...)` and `mesh%crop_predicate(array, compare=, value=)`. See the [Fortran reference](/fortran).
- **Julia** — `crop_bbox`, `crop_plane` and `crop_predicate(m, array; compare, value)`. See the [Julia reference](/julia).
- **R** — `mio_crop_bbox`, `mio_crop_plane` and `mio_crop_predicate(mesh, array, compare, value)`. See the [R reference](/r).
- **WebAssembly / JavaScript** — `cropBbox(mesh, lo, hi, mode, recordIds)`, `cropPlane(mesh, point, normal, mode, recordIds)` and `cropPredicate(mesh, array, compare, value, recordIds)`. See the [WebAssembly reference](/wasm).
- **Pipeline** — the `Crop` step takes `Where`/`Compare`/`Value` alongside its existing `Bbox` and `Point`+`Normal`. See the [pipeline reference](/pipeline).
