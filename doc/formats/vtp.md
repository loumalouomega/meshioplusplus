# VTP — VTK XML PolyData (`.vtp`)

The [VTK XML PolyData format](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html): the same XML container as [VTU](./vtu.md) with a `<PolyData>` grid holding `<Verts>/<Lines>/<Polys>/<Strips>` connectivity+offsets sections instead of `<Cells>`. Surface meshes only.

| | |
|---|---|
| **Format name** | `vtp` |
| **Extensions** | `.vtp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib compression via the optional zlib build in C++, stdlib `zlib` in Python) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("surface.vtp")
meshioplusplus.vtp.write("out.vtp", mesh, binary=True, compression="zlib")
```

- **`binary`** — base64 "binary" DataArrays (default `True`) vs ASCII.
- **`compression`** — `"zlib"` (default), `None`, or `"lzma"` (Python-only; the C++ path handles `None`/`"zlib"`).
- **`header_type`** — VTU-style header integer type; leave `None` (default `UInt32`) for the C++ fast path.

## Cell types

| PolyData section | meshio++ types |
|---|---|
| `Verts` | `vertex` (single-node rows only; poly-vertices are rejected) |
| `Lines` | `line` (two-node rows only; polylines are rejected) |
| `Polys` | `triangle` (3 nodes), `quad` (4), `polygon` (≥5) |
| `Strips` | not supported (rejected) |

Anything else — volume cells, quadratic cells, polyhedra — cannot be held by PolyData and raises `WriteError` in both the C++ and Python writers.

## Data mapping

- `point_data` / `cell_data` round-trip as in VTU. Cell data follows **VTK's canonical PolyData cell order** — Verts, then Lines, then Polys — in both directions; the writer reorders the mesh's blocks (and their cell_data) into that partition, preserving relative order within each section.
- Ragged (jagged) `polygon` blocks are written as-is; on read, polygons are grouped by node count into rectangular blocks.

## Quirks & limitations

- **PolyData has no cell-type array**, so the triangle/quad vs 3-/4-noded-polygon distinction cannot survive a round-trip: a `polygon` cell with 3 or 4 nodes reads back as `triangle`/`quad`. (Genuine ≥5-noded polygons are unaffected.)
- Triangle strips, poly-vertex/poly-line rows, multiple `<Piece>`s, `<AppendedData>`, and lzma compression are not handled by the C++ reader (the Python fallback additionally covers lzma).
- The C++ writer's output layout mirrors the VTU writer (same `<DataArray>` encoding, `%.11e` ASCII floats, 32 KiB zlib blocks).
