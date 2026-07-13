# FreeFem++ mesh (`.msh`)

The [FreeFem++](https://freefem.org/) `.msh` mesh format (as handled by
[FEconv](https://github.com/victorsndvg/FEconv)). ASCII, with a volume-element
block and a boundary-element block, each entity carrying an integer
region/boundary label.

| | |
|---|---|
| **Format name** | `freefem` |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

# `.msh` is shared with ansys and gmsh — be explicit on read:
mesh = meshio.read("mesh.msh", file_format="freefem")
meshio.freefem.write("out.msh", mesh)
```

`write` takes no keyword arguments.

## File structure

- **2D** — header `nver ntri nedge`, then `nver` vertices `x y ref`, `ntri`
  triangles `a b c ref`, and `nedge` boundary edges `a b ref`.
- **3D** — header `nver ntet ntri`, then vertices `x y z ref`, tetrahedra
  `a b c d ref`, and boundary triangles `a b c ref`.

All indices are 1-based. The spatial dimension is inferred from the first vertex
line.

## Cell types

Linear only: 2D meshes use `triangle` (volume) + `line` (boundary); 3D meshes use
`tetra` (volume) + `triangle` (boundary).

## Data mapping

- Vertex labels → `point_data["freefem:ref"]`.
- Element/boundary labels → `cell_data["freefem:ref"]`.

## Notes

- The `.msh` extension is shared with [`gmsh`](./gmsh.md) and `ansys`. On
  auto-detection meshio tries the registered formats in order and `freefem` is
  attempted last; pass `file_format="freefem"` to select it explicitly.
- Fully handled by the C++ core.
