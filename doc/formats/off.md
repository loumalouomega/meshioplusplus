# OFF (`.off`)

The [Object File Format](https://segeval.cs.princeton.edu/public/off_format.html): a minimal ASCII surface format — a vertex/face/edge count header, the vertex coordinates, then one line per face (leading vertex count + indices).

| | |
|---|---|
| **Format name** | `off` |
| **Extensions** | `.off` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("surface.off")
meshioplusplus.off.write("out.off", mesh)
```

`write` takes no keyword arguments.

## File structure

```
OFF
<nverts> <nfaces> <nedges>
x y z                    # nverts lines
n i0 i1 ... i(n-1)       # nfaces lines: leading vertex count + indices
```

The first line must be exactly `"OFF"` (`ReadError` otherwise). The counts line's edge count is parsed but discarded. Faces are grouped by vertex count into `triangle` (3), `quad` (4), or `polygon` (any other count ≥ 3) cell blocks — a run of same-count faces stays in one block until the count changes, mirroring the OBJ reader. A leading count below 3 raises `ReadError`.

## Cell types

`triangle`, `quad`, `polygon` (a `polygon` block written back out must be rectangular — every face in the block the same vertex count; the writer's Python reference also accepts a ragged `polygon` block, since it can iterate row-by-row). Any other cell type is skipped on write with a warning.

## Data mapping

None — OFF carries no point_data, cell_data, or field_data; `Mesh(points, cells)` only.

## Quirks & limitations

- Text-mode strictness: the Python reader requires a text-mode stream (raises if given bytes); the writer always opens the file itself in binary mode regardless of the caller's context.
- No boundary/edge data is ever produced — the edge count is read-and- discarded on read, and always written as `0`.
- Faces are grouped into blocks purely by **run** (consecutive same-count faces), not by regrouping every triangle/quad/polygon in the file into one block each — a file that interleaves counts (e.g. tri, quad, tri) produces three separate blocks in file order, same as the OBJ reader.

## Notes

- Fully handled by the C++ core.
- `tests/meshes/off/cube_example.off` (6 quad faces) / `cube_example_as_triangs.off` (the same cube pre-triangulated into 12 triangles) — a unit cube, once per representation (from [issue #35](https://github.com/loumalouomega/meshioplusplus/issues/35), which reported that non-triangular faces were rejected outright); other tests round-trip synthetic `tri_mesh`/`quad_mesh`/`polygon_mesh_one_cell` fixtures and check `.off`/`.0.off` extension dispatch.
