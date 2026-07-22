# WKT / TIN (`.wkt`)

The [Well-Known Text](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) representation of a [Triangulated Irregular Network](https://en.wikipedia.org/wiki/Triangulated_irregular_network): `TIN (((x y z, x y z, x y z, x y z)), …)` — one closed 4-point ring per triangle (the 4th point repeats the 1st, closing the ring).

| | |
|---|---|
| **Format name** | `wkt` |
| **Extensions** | `.wkt` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("surface.wkt")
meshioplusplus.wkt.write("out.wkt", mesh)
```

`write` takes no keyword arguments.

## File structure

A single `TIN (...)` expression. Each triangle is `((p0, p1, p2, p0))` — a 4-point closed linestring, one triangle per pair of nested parentheses. Points are **de-duplicated by exact floating-point value** in first-occurrence order (two points differing even in the last bit remain distinct — no epsilon tolerance); the repeated closing point of each ring is dropped once the 3 unique corner indices are recovered. A ring whose last point doesn't equal its first raises an error ("not a closed linestring").

The C++ reader parses this by tracking **parenthesis depth** rather than matching literal substrings — the triangle's point list sits at depth 3 (`TIN` → 1, the triangle polygon → 2, its linestring → 3) — which makes it naturally tolerant of arbitrary whitespace/newlines between and inside the nested parentheses.

## Cell types

`triangle` only — WKT never produces any other cell type.

## Data mapping

None — plain `Mesh(points, [CellBlock("triangle", ...)])`, no point_data, cell_data, or field_data.

## Quirks & limitations

- Point de-duplication is by **exact** value (no tolerance), so numerically- identical-but-differently-rounded coordinates are treated as distinct points if they don't parse to the exact same float.
- The write→read round trip for arbitrary meshes is currently disabled in the test suite (marked skip, pending re-enablement) — treat WKT primarily as a well-tested **reader** for now; the writer is implemented and does work, just not exercised by an automatic round-trip test.
- Whitespace handling is deliberately permissive on read — the reference `whitespaced.wkt` fixture has newlines and irregular spacing sprinkled throughout, including a point given with mixed-precision text (`0.00`, `.1`, `0.`) that still parses to the same set of unique points.

## Notes

- `tests/python/meshes/wkt/simple.wkt` — 2 triangles sharing an edge, 4 unique points, point-sum reference `4`.
- `tests/python/meshes/wkt/whitespaced.wkt` — the same 2 triangles with irregular internal whitespace and mixed-precision coordinate text, point-sum reference `3.2`.
- Fully handled by the C++ core.
