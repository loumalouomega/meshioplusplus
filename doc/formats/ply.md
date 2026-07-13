# PLY (`.ply`)

The [Polygon File Format](https://en.wikipedia.org/wiki/PLY_(file_format)) (aka
Stanford triangle format): a header describing named element groups with typed
properties, followed by ASCII or (little/big-endian) binary data.

| | |
|---|---|
| **Format name** | `ply` |
| **Extensions** | `.ply` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("bunny.ply")
meshio.ply.write("out.ply", mesh, binary=True)
```

- **`binary`** — binary (`True`, default) or ASCII.

## File structure

A `ply` magic line, a `format` line, `element vertex N` / `element face M`
declarations with `property` lines, then the data. Faces are a *list* property
(vertex count + indices).

## Cell types

Faces are grouped by vertex count into `vertex`, `line`, `triangle`, `quad`, and
`polygon` blocks.

## Data mapping

- Vertex properties beyond `x`/`y`/`z` (e.g. `nx,ny,nz`, `confidence`,
  `intensity`) → `point_data`.

## Notes

- The C++ core handles ASCII and little/big-endian binary. Files with faces
  carrying extra (non-index) properties, list vertex properties, or buffer I/O
  fall back to the Python implementation.
