# MOAB H5M (`.h5m`)

The [MOAB](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) mesh format,
stored in HDF5 under a `tstt` root group.

| | |
|---|---|
| **Format name** | `h5m` |
| **Extensions** | `.h5m` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.h5m")
meshio.h5m.write("out.h5m", mesh,
    add_global_ids=True,
    compression="gzip",
    compression_opts=4,
)
```

- **`add_global_ids`** — write a `GLOBAL_ID` node tag if absent.
- **`compression`** / **`compression_opts`** — HDF5 gzip filter.

## File structure

`tstt/nodes/coordinates` (+ `tags`), `tstt/elements/<Type>/connectivity`
(1-based, with an enum `element_type` attribute), committed tag datatypes, a
history dataset and `max_id`.

## Cell types

Read: `line`, `triangle`, `quad`, `tetra`, `pyramid`, `wedge`, `hexahedron`.
Write: `line`, `triangle`, `tetra` (other types are skipped with a warning,
matching the reference behaviour).

## Data mapping

- Node tags → `point_data` (including the auto `GLOBAL_ID`).

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise
  through `h5py`.
