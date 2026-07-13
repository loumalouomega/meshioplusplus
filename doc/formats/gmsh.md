# Gmsh (`.msh`)

The [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) mesh format,
supporting the MSH file versions **2.2**, **4.0** and **4.1**, in both ASCII and
binary.

| | |
|---|---|
| **Format name** | `gmsh` (v4.1), `gmsh22` (v2.2) |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.msh")          # version auto-detected
meshio.gmsh.write("out.msh", mesh,
    fmt_version="4.1",  # "2.2", "4.0", or "4.1"
    binary=True,
    float_fmt=".16e",
)
```

Via the generic API, `file_format="gmsh"` writes v4.1 and `file_format="gmsh22"`
writes v2.2.

- **`fmt_version`** — output MSH version.
- **`binary`** — write binary (`True`) or ASCII.
- **`float_fmt`** — coordinate format for ASCII output.

## File structure

Section-based (`$Nodes`/`$EndNodes`, `$Elements`/`$EndElements`,
`$PhysicalNames`, `$NodeData`/`$ElementData`, `$Entities`, `$Periodic`). v4.x
adds block structure with `size_t` ids. Node tags need not be contiguous; meshio
remaps them.

## Cell types

The full Gmsh element table including second-order elements. meshio applies the
Gmsh ↔ meshio node-order permutations (e.g. `tetra10`, `hexahedron20`,
`hexahedron27`, `wedge15`, `pyramid13`).

## Data mapping

- `$PhysicalNames` → `field_data`.
- `$NodeData` / `$ElementData` → `point_data` / `cell_data`.
- `gmsh:physical`, `gmsh:geometrical`, `gmsh:dim_tags`,
  `gmsh:bounding_entities` cell/point data.
- Physical groups → `cell_sets`; `gmsh_periodic` for periodic meshes.

## Notes

- The C++ core handles the common 2.2 and 4.1 read/write paths (ASCII and
  binary). Files with a `$Entities` section, `$Periodic` blocks, `gmsh:dim_tags`
  on write, or Gmsh 4.0 transparently fall back to the Python implementation.
- `.msh` is shared with [`ansys`](./ansys.md) and [`freefem`](./freefem.md);
  pass `file_format` to disambiguate on read.
