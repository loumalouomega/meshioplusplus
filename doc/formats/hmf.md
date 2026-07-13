# HMF (`.hmf`)

An experimental HDF5 mesh container used by meshio, reusing the XDMF topology
names. **The format may change at any time.**

| | |
|---|---|
| **Format name** | `hmf` |
| **Extensions** | `.hmf` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.hmf")
meshio.hmf.write("out.hmf", mesh, compression="gzip", compression_opts=4)
```

- **`compression`** / **`compression_opts`** — HDF5 gzip filter.

## File structure

Root attributes `type="hmf"`, `version="0.1-alpha"`; `/domain/grid` with
`Topology{k}` datasets (XDMF `TopologyType` attribute), a `Geometry` dataset
(`GeometryType` = `X`/`XY`/`XYZ`), and `NodeAttributes`/`CellAttributes` groups.

## Cell types

Multi-block, using the XDMF topology names.

## Data mapping

- `NodeAttributes` → `point_data`; `CellAttributes` (raw) → `cell_data`.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise
  through `h5py`. (The C++ reader correctly round-trips multi-block cell data,
  which the reference `h5py` reader mishandles.)
