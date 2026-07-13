# Modulef Formatted Mesh — MFM (`.mfm`)

The Modulef Formatted Mesh is a compact ASCII format used by
[FEconv](https://github.com/victorsndvg/FEconv) (a simplified NOPO/Modulef
mesh). A file holds a **single element type** (non-hybrid): an 8-integer header,
the connectivity, per-entity reference arrays, the vertex coordinates and a
per-element subdomain array.

| | |
|---|---|
| **Format name** | `mfm` |
| **Extensions** | `.mfm` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.mfm")
meshio.mfm.write("out.mfm", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format string (default `".16e"`).

## File structure

```
nel nnod nver dim lnn lnv lne lnf   # header
<connectivity, 1-based, element-major>
<face refs (dim==3)> <edge refs (dim>=2)> <vertex refs>
<vertex coordinates, dim per vertex>
<per-element subdomain array>
```

`lnv`/`lne`/`lnf` (local vertices/edges/faces) identify the element type; only
vertex coordinates are stored.

## Cell types

Linear elements only: `line`, `triangle`, `quad`, `tetra`, `hexahedron`,
`wedge`. Because MFM stores only vertex coordinates, higher-order (P2) elements
would be straight-sided, so meshio rejects them to avoid silent data loss.

## Data mapping

- Per-element subdomain / material → `cell_data["mfm:ref"]`.

## Notes

- Single element type per file: writing a mesh with more than one cell type
  raises a `WriteError`.
- Fully handled by the C++ core (with a Python fallback for buffers / non-default
  `float_fmt`).
