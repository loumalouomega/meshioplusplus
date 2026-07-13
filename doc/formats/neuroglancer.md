# Neuroglancer precomputed (no extension)

The [Neuroglancer precomputed](https://github.com/google/neuroglancer/tree/master/src/neuroglancer/datasource/precomputed#mesh-representation-of-segmented-object-surfaces)
mesh representation for segmented object surfaces: a small binary format (a
vertex count, the vertex coordinates, then the triangle indices).

| | |
|---|---|
| **Format name** | `neuroglancer` |
| **Extensions** | *(none — pass `file_format="neuroglancer"`)* |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh", file_format="neuroglancer")
meshio.neuroglancer.write("out", mesh)
```

`write` takes no keyword arguments. Because the format has no canonical file
extension, always pass `file_format="neuroglancer"`.

## File structure

A little-endian `uint32` vertex count, `3 * nverts` `float32` coordinates, then
the triangle vertex indices as `uint32`.

## Cell types

`triangle` only.

## Notes

- Implemented in pure Python (no C++ core path).
