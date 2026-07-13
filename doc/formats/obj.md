# Wavefront OBJ (`.obj`)

The [Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file)
geometry format: a line-oriented ASCII format with `v` (vertices), `vt`/`vn`
(texture/normal coordinates), `f` (faces) and `g` (groups).

| | |
|---|---|
| **Format name** | `obj` |
| **Extensions** | `.obj` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("model.obj")
meshio.obj.write("out.obj", mesh)
```

`write` takes no keyword arguments.

## File structure

Blank/`#`-comment lines are skipped.

- `v x y z` — a point row (all tokens after `v` are parsed as floats; extra
  columns are simply appended, not rejected).
- `vn ...` / `vt ...` — vertex normal / texture-coordinate rows (arbitrary
  column count, stored raw).
- `s 1` / `s off` — smooth-shading toggle, ignored.
- `f i1[/t1[/n1]] i2/... ...` — a face; only the part before the first `/`
  of each token is used as the vertex index (1-based). All faces in a "run"
  must have the same vertex count; the run breaks (opening a new cell block)
  as soon as the count changes, or a new `g` line appears.
- `g <name>` — starts a new group and increments a running group-id counter
  (starting at -1, so faces before any `g` get group id -1). Only the
  group's running **id** is kept — the name string itself is not stored.
- Any other keyword is silently ignored.

Empty trailing groups (e.g. a `g` line with no faces after it) are dropped
after the full file is scanned.

## Cell types

Faces are grouped by vertex count into `triangle` (3), `quad` (4), or
`polygon` (else).

## Data mapping

- `point_data["obj:vn"]` — vertex normals, present only if any `vn` lines
  were seen.
- `point_data["obj:vt"]` — texture coordinates, present only if any `vt`
  lines were seen.
- `cell_data["obj:group_ids"]` — one int array per cell block, the `g`-index
  that produced that block (`-1` if before the first `g`).

## Quirks & limitations

- The `/vt`/`/vn` per-face index references (e.g. `f 1/2/3`) are discarded
  entirely — only the leading vertex index survives; there is no per-face
  UV/normal index array, only the flat per-vertex `vn`/`vt` blocks.
- A new `g` line always starts a fresh cell block, **even if** the vertex
  count is unchanged from the previous group — group boundaries and vertex-
  count changes both independently trigger a new block.

## Notes

- `tests/meshes/obj/elephav.obj` — a real triangle-mesh scan (vertex-only `v`
  lines plus `f` triangle lines).
- Fully handled by the C++ core.
