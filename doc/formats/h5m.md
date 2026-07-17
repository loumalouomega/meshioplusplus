# MOAB H5M (`.h5m`)

The [MOAB](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) mesh format, stored in HDF5 under a `tstt` root group.

| | |
|---|---|
| **Format name** | `h5m` |
| **Extensions** | `.h5m` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.h5m")
meshioplusplus.h5m.write("out.h5m", mesh,
    add_global_ids=True,
    compression="gzip",
    compression_opts=4,
)
```

- **`add_global_ids`** — write a `GLOBAL_ID` node tag (`1..n`) if the mesh doesn't already have one.
- **`compression`** / **`compression_opts`** — HDF5 gzip filter and level.

## File structure

```
tstt/
  nodes/
    coordinates              # attr start_id=1
    tags/<name>               # per point-data key; 2D data via an HDF5 ARRAY dtype
  tags/<name>/                # global tag registry: committed datatype "type", attr class=2
  elements/<H5M_TYPE>/         # e.g. "Tet4", "Tri3", "Edge2", "Hex8", "Prism6", "Pyramid5", "Quad4"
    connectivity               # 1-based, attr start_id
    element_type                # attr, enum dtype, scalar
  elemtypes                   # committed enum datatype at the tstt level
  history                     # dataset of 3 fixed-length byte strings (module, version, timestamp)
  sets/tags/                  # empty group — "MOAB wants this"
  (attr) max_id                # running global-id counter, u8
```

2D point-data arrays are stored not as `(n,k)` datasets but as `(n,)` datasets of `k`-tuples via an HDF5 ARRAY/compound datatype — the C++ core replicates this exactly via `H5Tarray_create2`.

## Cell types

Read table (`H5M type` → meshio++ type):

| H5M | meshio++ | H5M | meshio++ |
|---|---|---|---|
| `Edge2` | `line` | `Pyramid5` | `pyramid` |
| `Tri3` | `triangle` | `Quad4` | `quad` |
| `Tet4` | `tetra` | `Hex8` | `hexahedron` |
| `Prism6` | `wedge` | | |

**Write only supports three types**: `line→Edge2`, `triangle→Tri3`, `tetra→Tet4`. Any other cell type is skipped (with a warning in Python; no warning in the equivalent C++ path).

MOAB element-type enum values (stored as the `element_type` attribute): `Edge=1, Tri=2, Quad=3, Polygon=4, Tet=5, Pyramid=6, Prism=7, Knife=8, Hex=9, Polyhedron=10`.

## Data mapping

- Arbitrary `point_data` keys → datasets under `nodes/tags/<name>`, plus a registry entry under `tstt/tags/<name>`.
- `GLOBAL_ID` — conventional auto-added point tag (`add_global_ids=True`).
- **No cell_data support end-to-end** — see quirks below.

## Quirks & limitations

- Element/cell tags and MOAB "sets" are **not read at all**, even though MOAB itself supports them — the reader ignores `elements/*/tags` and the `sets` group entirely.
- The reference Python writer's cell-data code path has a **pre-existing bug**: it iterates `mesh.cell_data.items()` while treating the last `elements` sub-group from a prior loop as if it applied to every cell type — effectively dead/incorrect code for meshio++'s actual `cell_data` schema. **The C++ writer deliberately does not attempt to replicate this bug** — it simply never writes cell data. The shim additionally only attempts the C++ write path when `mesh.cell_data` is empty, so any cell-data-bearing mesh always falls back to the (buggy) Python writer.
- Indices are 1-based in the file; `+1`/`-1` applied on write/read.
- All-1-based `start_id` attributes track a running global-id counter shared across nodes and every element block.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise through `h5py`.
- No reference fixture exists under `tests/meshes/h5m/`; `tests/test_moab.py` round-trips synthetic meshes only.
