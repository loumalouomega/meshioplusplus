# DOLFIN XML (`.xml`)

The legacy [DOLFIN/FEniCS](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html)
XML mesh format. The mesh lives in one file; each cell-data array is stored in a
sibling `<stem>_<name>.xml` file.

| | |
|---|---|
| **Format name** | `dolfin-xml` |
| **Extensions** | `.xml` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.xml")
meshio.dolfin.write("out.xml", mesh)
```

`write` takes no keyword arguments.

## File structure

```xml
<dolfin><mesh celltype="triangle|tetrahedron" dim="2|3">
  <vertices size="N"><vertex index="i" x="…" y="…" [z="…"]/></vertices>
  <cells size="M"><triangle index="i" v0="…" v1="…" v2="…"/></cells>
</mesh></dolfin>
```

Vertices and cells are placed by their `index` attribute.

## Cell types

`triangle` and `tetra` only.

## Data mapping

- Each `cell_data` array → a sibling file `<stem>_<name>.xml` holding a
  `<mesh_function>`; the reader scans the mesh file's directory for them.

## Notes

- Fully handled by the C++ core (via the vendored pugixml + `std::filesystem`).
