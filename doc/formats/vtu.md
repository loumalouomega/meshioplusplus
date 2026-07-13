# VTU — VTK XML UnstructuredGrid (`.vtu`)

The [serial VTK XML](https://vtk.org/Wiki/VTK_XML_Formats) UnstructuredGrid
format: an XML container with inline (ASCII or base64 binary) or appended data
arrays, optionally compressed.

| | |
|---|---|
| **Format name** | `vtu` |
| **Extensions** | `.vtu` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.vtu")
meshio.vtu.write("out.vtu", mesh,
    binary=True,
    compression="zlib",  # "zlib", "lzma", or None
    header_type=None,    # "UInt32" or "UInt64"
)
```

- **`binary`** — base64-encoded binary (`True`) or ASCII arrays.
- **`compression`** — block compression filter for binary data.
- **`header_type`** — integer type for the binary block header.

## File structure

`<VTKFile><UnstructuredGrid><Piece>` with `<Points>`, `<Cells>` (connectivity /
offsets / types), `<PointData>` and `<CellData>` `DataArray`s. Binary arrays use
the VTK base64 header + (optionally zlib/lzma block-compressed) payload.

## Cell types

The full VTK cell set, including VTK Lagrange high-order cells
(`VTK_LAGRANGE_*`).

## Data mapping

- `<PointData>` / `<CellData>` → `point_data` / `cell_data`.
- Cell sets are supported.

## Notes

- The C++ core handles ASCII, uncompressed binary and **zlib** (when built with
  `MESHIO_WITH_ZLIB`; otherwise the Python stdlib handles it). **lzma**,
  appended/raw-binary data and `polyhedron` cells fall back to the Python
  implementation. See [native acceleration](../formats.md#native-acceleration-and-fallbacks).
