# Wavefront OBJ (`.obj`)

The [Wavefront OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) geometry
format: a line-oriented ASCII format with `v` (vertices), `vt`/`vn` (texture /
normal coordinates), `f` (faces) and `g` (groups).

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

## Cell types

Faces are grouped by vertex count into `triangle`, `quad` and `polygon` blocks.

## Data mapping

- `vn` (vertex normals) → `point_data["obj:vn"]`.
- `vt` (texture coordinates) → `point_data["obj:vt"]`.
- `g` group membership → `cell_data["obj:group_ids"]`.

## Notes

- Fully handled by the C++ core.
