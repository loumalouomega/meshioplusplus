# OFF (`.off`)

The [Object File Format](https://segeval.cs.princeton.edu/public/off_format.html):
a minimal ASCII surface format — a vertex count / face count header, the vertex
coordinates, then one line per face (vertex count + indices).

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
<x y z>            # nverts lines
3 i j k            # nfaces lines: leading vertex count + indices
```

## Cell types

`triangle` only.

## Notes

- Fully handled by the C++ core.
