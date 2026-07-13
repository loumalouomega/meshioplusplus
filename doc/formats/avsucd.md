# AVS-UCD (`.avs`)

The [AVS Unstructured Cell Data](https://lanl.github.io/LaGriT/pages/docs/read_avs.html)
format: an ASCII header of counts, node coordinates, cells (with a per-cell
material), then node- and cell-data sections.

| | |
|---|---|
| **Format name** | `avsucd` |
| **Extensions** | `.avs` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.avs")
meshio.avsucd.write("out.avs", mesh)
```

`write` takes no keyword arguments.

## File structure

Header `nnodes ncells ndata_node ndata_cell nmodel`, then `id x y z` nodes, then
`id material type n1 n2 …` cells, then node-data and cell-data blocks (scalar or
multi-component).

## Cell types

`line`, `triangle`, `quad`, `tetra`, `pyramid`, `wedge`, `hexahedron`, with the
AVS ↔ meshio node orderings.

## Data mapping

- Per-cell material → `cell_data["avsucd:material"]`.
- Node/cell data sections → `point_data` / `cell_data`.

## Notes

- Fully handled by the C++ core.
