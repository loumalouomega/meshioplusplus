# Netgen (`.vol`, `.vol.gz`)

The [Netgen](https://github.com/ngsolve/netgen) neutral mesh format (`.vol`),
optionally gzip-compressed (`.vol.gz`).

| | |
|---|---|
| **Format name** | `netgen` |
| **Extensions** | `.vol`, `.vol.gz` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.vol")
meshio.netgen.write("out.vol", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format.

## File structure

`mesh3d`, `dimension`, `geomtype`, `points`, and the `pointelements` /
`edgesegments` / `surfaceelements` / `volumeelements` sections, each with a
per-dimension column layout. A single integer cell index is stored per element.

## Cell types

The full Netgen ↔ meshio node permutations, including `triangle6`, `quad8`,
`tetra10`, `hexahedron20`, `wedge15`, `pyramid13`.

## Data mapping

- Per-cell index → `cell_data["netgen:index"]`.
- Periodic identifications → `mesh.info["netgen:identifications"]` /
  `["netgen:identificationtypes"]`.

## Notes

- The C++ core handles the common `.vol` path (points, cells, `netgen:index`).
  The gzip `.vol.gz` container, `identifications`, material/bc/`cd*names`
  sections (driven by `field_data`), and the two-line `edgesegmentsgi2` variant
  fall back to the Python implementation.
