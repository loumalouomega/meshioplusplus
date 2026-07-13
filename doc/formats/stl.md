# STL (`.stl`)

The [stereolithography](https://en.wikipedia.org/wiki/STL_(file_format))
format: a flat list of triangle facets (each with a normal and three
vertices), in ASCII or binary. STL has no shared-vertex table — every facet
repeats its own vertex coordinates.

| | |
|---|---|
| **Format name** | `stl` |
| **Extensions** | `.stl` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("part.stl")
meshio.stl.write("out.stl", mesh, binary=False)
```

- **`binary`** — binary (`True`) or ASCII (`False`, the default — unlike
  most binary-capable formats in meshio, which default to `True`).

## File structure

**ASCII**:
```
solid [name]
 facet normal nx ny nz
  outer loop
   vertex x y z
   vertex x y z
   vertex x y z
  endloop
 endfacet
 ...
endsolid
```

**Binary**: an 80-byte free-form header, then a little-endian `uint32`
triangle count, then that many fixed 50-byte records: `float32[3]` normal,
`float32[3][3]` facet vertices, `int16` attribute byte count (conventionally
0, not enforced).

**Binary-vs-ASCII detection**: if the file is under 80 bytes, it's treated
as ASCII. Otherwise, the 80-byte header and triangle count are read and the
expected total size `84 + num_triangles*50` is compared against the actual
file size; if they match, it's binary, otherwise the reader rewinds and
parses as ASCII. This heuristic deliberately avoids the naive "does the file
start with the literal word `solid`" check, since binary STL files sometimes
also start with that word in their free-form header.

## Cell types

`triangle` only.

## Data mapping

- `cell_data["facet_normals"]` — if present, written verbatim; otherwise the
  writer computes normals from the cross product of two edge vectors.

## Quirks & limitations

- **Vertex de-duplication**: all raw (possibly duplicate) triangle vertices
  are uniquified in **first-occurrence order** to build a shared point
  table — meaning point *indices* are not preserved across a round-trip, only
  the geometry.
- **dtype**: ASCII coordinates parse as `float64`; binary coordinates parse
  as `float32` (matching the on-disk format) — this dtype difference is
  intentional, not a bug.
- **Empty STL**: a file with zero triangles produces a `Mesh` with **no**
  cell blocks at all, not an empty `"triangle"` block.
- Mixed cell types on write: a warning names the discarded types, and only
  triangles are written.
- ASCII parsing uses a fast custom line reader that only looks at the last 3
  whitespace-separated tokens per line (discarding any leading keyword like
  `vertex`), for performance over `np.loadtxt`.

## Notes

- Fully handled by the C++ core (ASCII and binary).
- No reference fixture exists under `tests/meshes/stl/`; tests round-trip
  synthetic meshes only.
