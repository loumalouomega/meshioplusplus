# CGNS (`.cgns`)

The [CGNS](https://cgns.github.io/) (CFD General Notation System) format, in this implementation stored as a minimal, non-standard subset within an HDF5 container — a tetrahedra-only mesh, not the full CGNS/SIDS specification.

| | |
|---|---|
| **Format name** | `cgns` |
| **Extensions** | `.cgns` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.cgns")
meshioplusplus.cgns.write("out.cgns", mesh, compression="gzip", compression_opts=4)
```

- **`compression`** / **`compression_opts`** — HDF5 gzip filter and level.

## File structure

```
Base/
  Zone1/
    GridCoordinates/
      CoordinateX/" data"     # dataset name is literally " data" (leading space)
      CoordinateY/" data"
      CoordinateZ/" data"
    GridElements/
      ElementRange/" data"           -> [1, n_cells]   (1-based inclusive)
      ElementConnectivity/" data"    -> flat 1-based tetra connectivity
```

The leading-space dataset name `" data"` in every leaf group is not part of the real CGNS/HDF5 spec — it's this implementation's own ad hoc convention, consistent between the Python and C++ writers.

## Cell types

`tetra` only — this is the only cell type the reader accepts and the only one the writer emits. Reading a connectivity array that doesn't reshape to exactly 4 columns raises `ReadError("Can only read tetrahedra.")`.

## Data mapping

None — no point_data, cell_data, or field_data is read or written by this format.

## Quirks & limitations

- Read requires `"Base"` and `"Base/Zone1"` to be present, else `ReadError('Expected "Base" in file. Malformed CGNS?')` (and similarly for `"Zone1"`).
- Indices are 1-based in the file; `-1`/`+1` conversions are applied on read/write while preserving the connectivity array's original integer dtype.
- The writer only ever emits a `"tetra"` cell block — any other cell type present in `mesh.cells` is silently ignored (not warned).
- This is the least complete format meshio++ supports: no compression-aware reading (compression is a write-only concept here — HDF5 handles decompression transparently), no field/point/cell data at all.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise through `h5py` — behavior is identical either way.
- No reference fixture exists under `tests/python/meshes/cgns/`; the only test round-trips a synthetic `tet_mesh`.
