# WKT / TIN (`.wkt`)

The [Well-Known Text](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry)
representation of a
[Triangulated Irregular Network](https://en.wikipedia.org/wiki/Triangulated_irregular_network):
`TIN (((x y z, x y z, x y z, x y z)), …)`.

| | |
|---|---|
| **Format name** | `wkt` |
| **Extensions** | `.wkt` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("surface.wkt")
meshio.wkt.write("out.wkt", mesh)
```

`write` takes no keyword arguments.

## File structure

A single `TIN (…)` expression containing one `((p0, p1, p2, p0))` closed
linestring per triangle. Points are de-duplicated by exact value in
first-occurrence order; the repeated closing vertex is dropped.

## Cell types

`triangle` only.

## Notes

- Because points are de-duplicated, point indices are not preserved across a
  round-trip (the triangle geometry is).
- Fully handled by the C++ core (whitespace-tolerant parenthesis parsing).
