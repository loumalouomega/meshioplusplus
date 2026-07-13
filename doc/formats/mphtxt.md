# COMSOL text mesh (`.mphtxt`)

The [COMSOL Multiphysics](https://www.comsol.com) `.mphtxt` text mesh format
(as handled by [FEconv](https://github.com/victorsndvg/FEconv)). ASCII, storing
a version, tag/type name tables and one or more mesh objects; comments run from
`#` to end of line.

| | |
|---|---|
| **Format name** | `mphtxt` |
| **Extensions** | `.mphtxt` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.mphtxt")
meshio.mphtxt.write("out.mphtxt", mesh)
```

`write` takes no keyword arguments.

## File structure

A mesh object stores the space dimension, the node coordinates, and a series of
**element-type blocks** — one per COMSOL element type (`vtx`, `edg`, `tri`,
`quad`, `tet`, `prism`, `pyr`, `hex`, and their second-order `*2` variants). Each
block lists its node count, its connectivity, and a per-element *geometric entity
index*.

## Cell types

Hybrid meshes are supported. COMSOL ↔ meshio type map:

| COMSOL | meshio | COMSOL | meshio |
|---|---|---|---|
| `edg` | `line` | `edg2` | `line3` |
| `tri` | `triangle` | `tri2` | `triangle6` |
| `quad` | `quad` | `quad2` | `quad9` |
| `tet` | `tetra` | `tet2` | `tetra10` |
| `prism` | `wedge` | `hex` | `hexahedron` |
| `pyr` | `pyramid` | `hex2` | `hexahedron27` |

The `quad` and `hexahedron` connectivities are reordered between COMSOL's z-order
and meshio's counter-clockwise convention.

## Data mapping

- Per-element geometric entity index → `cell_data["mphtxt:geom"]`.

## Notes

- Fully handled by the C++ core.
- Higher-order node orderings round-trip losslessly but may differ from COMSOL's
  internal ordering for some element types.
