# Tecplot (`.dat`, `.tec`)

The [Tecplot ASCII](http://paulbourke.net/dataformats/tp/) data format: a
`VARIABLES` list and one or more finite-element `ZONE`s.

| | |
|---|---|
| **Format name** | `tecplot` |
| **Extensions** | `.dat`, `.tec` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("field.dat")
meshio.tecplot.write("out.dat", mesh)
```

`write` takes no keyword arguments.

## File structure

`VARIABLES = "X" "Y" "Z" …`, then a `ZONE` with a `ZONETYPE`
(`FETRIANGLE`/`FEQUADRILATERAL`/`FETETRAHEDRON`/`FEBRICK`), a `DATAPACKING`
(`BLOCK` or `POINT`), an optional `VARLOCATION` for cell-centred data, and the
node/connectivity data.

## Cell types

A single FE zone type per file: `triangle`, `quad`, `tetra`, `hexahedron`.

## Data mapping

- Nodal variables → `point_data`; cell-centred variables (via `VARLOCATION`) →
  `cell_data`.

## Notes

- The C++ core handles single-zone FE meshes (BLOCK and POINT packing,
  `VARLOCATION`). Multiple cell types, buffers, and some adversarial
  whitespace/quoting cases fall back to the Python implementation.
