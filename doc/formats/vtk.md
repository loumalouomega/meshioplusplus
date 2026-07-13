# VTK legacy (`.vtk`)

The [VTK legacy](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
file format (UNSTRUCTURED_GRID), versions **4.2** and **5.1**, in ASCII and
big-endian binary.

| | |
|---|---|
| **Format name** | `vtk` (v5.1), `vtk42`, `vtk51` |
| **Extensions** | `.vtk` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.vtk")
meshio.vtk.write("out.vtk", mesh, binary=True)   # version via file_format
```

- **`binary`** — big-endian binary (`True`) or ASCII.
- Version: `file_format="vtk"` / `"vtk51"` writes 5.1; `"vtk42"` writes 4.2.

## File structure

A text header, `POINTS`, and a `CELLS`/`CELL_TYPES` block. Version 4.2 uses the
interleaved `[count, idx…]` cell layout; 5.1 uses separate offsets/connectivity.
`POINT_DATA`/`CELL_DATA` sections carry field data.

## Cell types

The full VTK cell set including VTK Lagrange cells. Cell/point sets supported in
v5.1.

## Notes

- The C++ core handles both versions, ASCII and big-endian binary, for
  UNSTRUCTURED_GRID. Structured grids, SCALARS/VECTORS-only files, 2-component
  vector padding and `polyhedron` cells fall back to the Python implementation.
