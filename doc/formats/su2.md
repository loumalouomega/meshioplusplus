# SU2 (`.su2`)

The [SU2](https://su2code.github.io/docs_v7/Mesh-File/) mesh format: an ASCII
format with `NDIME`, `NPOIN`, `NELEM` volume cells and `NMARK` boundary markers.

| | |
|---|---|
| **Format name** | `su2` |
| **Extensions** | `.su2` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.su2")
meshio.su2.write("out.su2", mesh)
```

`write` takes no keyword arguments.

## File structure

`NDIME= <2|3>`, `NPOIN= …` with the coordinates, `NELEM= …` with VTK-type-coded
volume cells, then `NMARK= …` boundary-marker blocks (each with a tag name and
its boundary cells).

## Cell types

Volume cells use VTK type codes; boundary markers become `line` (2D) or
`triangle`/`quad` (3D) blocks.

## Data mapping

- Boundary-marker tag → `cell_data["su2:tag"]`.

## Notes

- Fully handled by the C++ core (per-type binning by VTK code, boundary-block
  merging).
