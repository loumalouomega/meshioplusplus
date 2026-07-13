# Medit / INRIA (`.mesh`, `.meshb`)

The [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) mesh
format (also the INRIA `libMeshb` format): keyword sections in ASCII (`.mesh`)
or binary GMF (`.meshb`).

| | |
|---|---|
| **Format name** | `medit` |
| **Extensions** | `.mesh`, `.meshb` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.mesh")
meshio.medit.write("out.mesh", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format for ASCII output.

## File structure

Keyword sections: `MeshVersionFormatted`, `Dimension`, `Vertices`, then element
keywords `Edges`, `Triangles`, `Quadrilaterals`, `Tetrahedra`, `Prisms`,
`Pyramids`, `Hexahedra` — each a count then rows ending in a reference tag. The
version drives the coordinate dtype (1 → float32, 2 → float64).

## Cell types

`line`, `triangle`, `quad`, `tetra`, `wedge`, `pyramid`, `hexahedron`.

## Data mapping

- Vertex reference → `point_data["medit:ref"]`.
- Per-element reference → `cell_data["medit:ref"]`.

## Notes

- The C++ core handles the ASCII `.mesh` variant. The binary `.meshb` GMF variant
  (little- and big-endian) falls back to the Python implementation.
