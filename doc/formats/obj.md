# Wavefront OBJ (`.obj`)

The [Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) geometry format: a line-oriented ASCII format with `v` (vertices), `vt`/`vn` (texture/normal coordinates), `f` (faces) and `g` (groups).

| | |
|---|---|
| **Format name** | `obj` |
| **Extensions** | `.obj` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

The native row parser views the live trimmed input row rather than copying it first; result arrays/group names remain owned. See the [source-copy evidence](../benchmarks.md#source-copy-removal-and-temporary-ownership-guards), including the retained unfavorable writer observation; no writer speedup is claimed.

```python
import meshioplusplus

mesh = meshioplusplus.read("model.obj")
meshioplusplus.obj.write("out.obj", mesh)
```

`write` takes no keyword arguments.

## File structure

Blank/`#`-comment lines are skipped.

- `v x y z` — a point row (all tokens after `v` are parsed as floats; extra columns are simply appended, not rejected).
- `vn ...` / `vt ...` — vertex normal / texture-coordinate rows (arbitrary column count, stored raw).
- `s 1` / `s off` — smooth-shading toggle, ignored.
- `f i1[/t1[/n1]] i2/... ...` — a face; only the part before the first `/` of each token is used as the vertex index (1-based). All faces in a "run" must have the same vertex count; the run breaks (opening a new cell block) as soon as the count changes, or a new `g` line appears.
- `g <name>` — starts a new group and increments a running group-id counter (starting at -1, so faces before any `g` get group id -1). Only the group's running **id** is kept — the name string itself is not stored.
- Any other keyword is silently ignored.

Empty trailing groups (e.g. a `g` line with no faces after it) are dropped after the full file is scanned.

## Cell types

Faces are grouped by vertex count into `triangle` (3), `quad` (4), or `polygon` (else).

## Data mapping

- `point_data["obj:vn"]` — vertex normals, present only if any `vn` lines were seen.
- `point_data["obj:vt"]` — texture coordinates, present only if any `vt` lines were seen; when `uv` is read (below), only if there is one `vt` row per point.
- `point_data["uv"]` — texture coordinates by face corner: when every corner names a `vt` (`f v/vt`, `f v/vt/vn`), each point gets the `u, v` its corners name. A point whose corners name different `vt` rows (a UV seam) is split: the first copy keeps its index, the others are appended in the order met, and `point_data["obj:parent_point"]` names the point each copies (the convention of [parameterize](../parameterize.md)). Both engines read it.
- Written: with `uv` point data (one pair per point), one `vt u v` per point and every face corner as `v/vt`, in place of the raw `obj:vt` rows; without it, faces are plain `f v` as before.
- `cell_data["obj:group_ids"]` — one int array per cell block, the `g`-index that produced that block (`-1` if before the first `g`).

## Quirks & limitations

- The `/vn` per-face index references are discarded — only the vertex and texture indices are read. A face with no texture index on some corner leaves `uv` unset for the whole file. Raw `obj:vt`/`obj:vn` rows that are not one per point are kept as the file has them, which the Python `Mesh` refuses (point data must have one row per point), except `obj:vt` once `uv` holds the coordinates.
- A new `g` line always starts a fresh cell block, **even if** the vertex count is unchanged from the previous group — group boundaries and vertex- count changes both independently trigger a new block.

## Notes

- `tests/python/meshes/obj/elephav.obj` — a real triangle-mesh scan (vertex-only `v` lines plus `f` triangle lines).
- Fully handled by the C++ core.
