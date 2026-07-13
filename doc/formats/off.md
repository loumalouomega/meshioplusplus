# OFF (`.off`)

The [Object File Format](https://segeval.cs.princeton.edu/public/off_format.html):
a minimal ASCII surface format — a vertex/face/edge count header, the vertex
coordinates, then one line per face (leading vertex count + indices).

| | |
|---|---|
| **Format name** | `off` |
| **Extensions** | `.off` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("surface.off")
meshio.off.write("out.off", mesh)
```

`write` takes no keyword arguments.

## File structure

```
OFF
<nverts> <nfaces> <nedges>
x y z              # nverts lines
3 i j k            # nfaces lines: leading vertex count (must be 3) + indices
```

The first line must be exactly `"OFF"` (`ReadError` otherwise). The counts
line's edge count is parsed but discarded. Every face row's leading count
must be `3` — any other value raises `ReadError("Can only read triangular
faces")`.

## Cell types

`triangle` only.

## Data mapping

None — OFF carries no point_data, cell_data, or field_data; `Mesh(points,
cells)` only.

## Quirks & limitations

- Text-mode strictness: the Python reader requires a text-mode stream
  (raises if given bytes); the writer always opens the file itself in binary
  mode regardless of the caller's context.
- No boundary/edge data is ever produced — the edge count is read-and-
  discarded on read, and always written as `0`.

## Notes

- Fully handled by the C++ core.
- No reference fixture exists under `tests/meshes/off/`; tests round-trip a
  synthetic `tri_mesh` and check `.off`/`.0.off` extension dispatch.
