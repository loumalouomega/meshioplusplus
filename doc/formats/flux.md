# FLUX mesh (`.pf3`)

The [Altair FLUX](https://www.altair.com/flux/) `.pf3` mesh format (as handled by
[FEconv](https://github.com/victorsndvg/FEconv)). ASCII with French keyword
headers.

| | |
|---|---|
| **Format name** | `flux` |
| **Extensions** | `.pf3` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.pf3")
meshio.flux.write("out.pf3", mesh)
```

`write` takes no keyword arguments.

## File structure

A header of `<count> NOMBRE …` keyword lines (dimensions, element counts, point
count, region counts) is followed by:

- `DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS` — one 12-integer record per element
  (the 4th field is the region reference, the 7th the type descriptor, the 8th
  the node count) followed by its 1-based connectivity.
- `COORDONNEES DES NOEUDS` — one `index x y z` line per node.

## Cell types

The type descriptor selects the element type: `vertex`, `line`/`line3`,
`triangle`/`triangle6`, `quad`/`quad8`, `tetra`/`tetra10`, `wedge`/`wedge15`,
`hexahedron`/`hexahedron20`, `pyramid`. Hybrid meshes are supported.

## Data mapping

- Per-element region reference → `cell_data["pf3:ref"]`.

## Notes

- Fully handled by the C++ core.
- Region *names* (which FLUX may store in binary) are not read. Node orderings
  round-trip losslessly but may differ from FLUX's internal ordering for some
  element types.
