# CGNS (`.cgns`)

The [CGNS](https://cgns.github.io/) (CFD General Notation System) format, stored
in its HDF5 (ADF/HDF) container.

| | |
|---|---|
| **Format name** | `cgns` |
| **Extensions** | `.cgns` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.cgns")
meshio.cgns.write("out.cgns", mesh, compression="gzip", compression_opts=4)
```

- **`compression`** / **`compression_opts`** — HDF5 gzip filter.

## File structure

`Base/Zone1/GridCoordinates/CoordinateX|Y|Z` and
`Base/Zone1/GridElements/ElementRange` + `ElementConnectivity` (1-based),
matching meshio's minimal CGNS dialect.

## Cell types

`tetra` only.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise
  through `h5py`; behaviour is identical either way. See
  [native acceleration](../formats.md#native-acceleration-and-fallbacks).
