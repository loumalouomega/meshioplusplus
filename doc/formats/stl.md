# STL (`.stl`)

The [stereolithography](https://en.wikipedia.org/wiki/STL_(file_format)) format:
a flat list of triangle facets (each with a normal and three vertices), in ASCII
or binary. STL has no shared-vertex table.

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

- **`binary`** — binary (`True`) or ASCII (`False`, default).

## File structure

ASCII: `solid … facet normal … outer loop … vertex x y z … endloop endfacet …
endsolid`. Binary: an 80-byte header, a `uint32` triangle count, then 50 bytes
per triangle. On read, meshio de-duplicates vertices in first-occurrence order
to build a shared point table.

## Cell types

`triangle` only.

## Notes

- Because STL stores raw facet coordinates, point indices are not preserved
  across a round-trip (the geometry is, up to vertex de-duplication).
- Fully handled by the C++ core (ASCII and binary).
