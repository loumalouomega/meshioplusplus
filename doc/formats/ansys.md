# Ansys / Fluent mesh (`.msh`)

The Ansys Fluent `.msh` mesh format: hexadecimal, parenthesis-nested sections, in
ASCII or binary.

| | |
|---|---|
| **Format name** | `ansys` |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.msh", file_format="ansys")
meshio.ansys.write("out.msh", mesh, binary=True)
```

- **`binary`** — binary (`True`) or ASCII.

## File structure

Sections keyed by a hexadecimal index: `(1 header)`, `(2 dim)`, node declaration
and data `(10 …)` / `(3010 …)`, and cell declaration and data `(12 …)` /
`(2012 …)` / `(3012 …)`. Node indices are 1-based on disk.

## Cell types

The Fluent element-type map: `triangle`, `tetra`, `quad`, `hexahedron`,
`pyramid`, `wedge`. 2D and 3D meshes are supported.

## Notes

- The C++ core handles the header/dimension/node/cell sections in ASCII and
  binary. Face sections `(13 …)` with a data body fall back to the Python
  implementation; zone specs `(39)`/`(45)` are skipped.
- `.msh` is shared with [`gmsh`](./gmsh.md) and [`freefem`](./freefem.md); on
  auto-detection `ansys` is tried first. Pass `file_format` to disambiguate.
