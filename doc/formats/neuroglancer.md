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

`write` takes no keyword arguments. `register_format` is called with an
**empty extension list** (`register_format("neuroglancer", [], read,
{"neuroglancer": write})`), so this format is never auto-detected from a
filename — `file_format="neuroglancer"` must be passed explicitly on both
read and write.

## File structure

Raw little-endian binary, no header/magic beyond the leading count:

1. `<I` (`uint32`) — number of vertices, `n`.
2. `n * 3` × `<f` (`float32`) — flattened `x0 y0 z0 x1 y1 z1 ...` vertex
   coordinates.
3. The **remainder of the file** — flattened `<I` (`uint32`) triangle
   vertex indices, 3 per triangle; the remaining byte length must be a
   multiple of 12.

Read uses `np.frombuffer` rather than `np.fromfile`, specifically so the
same code works transparently when handed a `gzip.GzipFile`-like file
object, not just a plain path.

## Cell types

`triangle` only.

## Data mapping

None — no point_data, cell_data, or field_data; a bare
`Mesh(points, [CellBlock("triangle", ...)])`.

## Quirks & limitations

- **Write silently casts points to `float32`** regardless of the input
  mesh's dtype — unlike most other meshio writers, there is no warning for
  this precision loss.
- Triangle-index validation is `np.any(triangles > num_vertices)` — using
  `>` rather than `>=` — so an index **exactly equal to** `num_vertices`
  (out of range, since indices are 0-based) would not be caught.
- No C++ implementation — this format is Python-only.

## Notes

- `tests/meshes/neuroglancer/simple1` (100 bytes, no file extension) —
  checked for `ref_sum=20` and `ref_num_cells=4`.
- Implemented in pure Python (no C++ core path).
