# FLAC3D (`.f3grid`)

The [Itasca FLAC3D](https://www.itascacg.com/software/flac3d) grid format
(`.f3grid`): ASCII or binary, with separate ZONE (3D) and FACE (2D) sections and
cell groups.

| | |
|---|---|
| **Format name** | `flac3d` |
| **Extensions** | `.f3grid` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("grid.f3grid")
meshio.flac3d.write("out.f3grid", mesh, float_fmt=".16e", binary=False)
```

- **`float_fmt`** — coordinate format (ASCII).
- **`binary`** — binary (`True`) or ASCII.

## File structure

`GRIDPOINTS`, `ZONES` (3D cells) and `FACES` (2D cells), each with the meshio ↔
FLAC3D node orderings; `ZGROUP`/`FGROUP` blocks name cell groups. Binary is
auto-detected by a null byte in the first bytes.

## Cell types

Zones: `tetra`, `pyramid`, `wedge`, `hexahedron` (and B7 → hexahedron); faces:
`triangle`, `quad`. Zone corners are reordered to a right-handed system via the
scalar triple product.

## Data mapping

- Global cell ids → `cell_data["cell_ids"]`.
- `ZGROUP`/`FGROUP` groups → `cell_sets`.

## Notes

- The C++ core handles points and zone/face cells (ASCII and binary). Cell groups
  (`cell_sets`) fall back to the Python implementation.
