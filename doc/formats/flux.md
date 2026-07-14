# FLUX mesh (`.pf3`)

The [Altair FLUX](https://www.altair.com/flux/) `.pf3` mesh format (as handled
by [FEconv](https://github.com/victorsndvg/FEconv)). ASCII with French
keyword headers.

| | |
|---|---|
| **Format name** | `flux` |
| **Extensions** | `.pf3` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.pf3")
meshioplusplus.flux.write("out.pf3", mesh)
```

`write` takes no keyword arguments.

## File structure

Header lines are located by **substring search** for their French label
(rather than by fixed line position), so header ordering is somewhat
tolerant:

- `dim` — the line containing `"NOMBRE DE DIMENSIONS"` (leading integer).
- `nel` — the line containing `"D'ELEMENTS"` but *not* any of `VOLUMIQUES`,
  `SURFACIQUES`, `LINEIQUES`, `PONCTUELS`, `MACRO` (i.e. the grand-total line,
  not a per-category count).
- `nnod` — the line containing `"NOMBRE DE POINTS"` but not `"INTEGRATION"`.

The element block starts after the line containing `"DESCRIPTEUR DE
TOPOLOGIE"`; the coordinate block starts after `"COORDONNEES DES NOEUDS"`.
Between those markers, `nel` element records are parsed as one continuous
token stream (line breaks inside a record are tolerated), each a 12-integer
header followed by its connectivity:

| field | meaning |
|---|---|
| 0 | element id (ignored on read) |
| 1 | `desc1` (category — ignored on read) |
| 2 | `desc2` (subtype — ignored on read) |
| 3 | **region reference** → `cell_data["pf3:ref"]` |
| 4 | node count (redundant with field 7) |
| 5 | unused (0) |
| 6 | **`desc3`** (type code — selects the meshio++ type) |
| 7 | node count (drives how many ids follow) |
| 8-11 | unused (0 0 0 0) |

followed by the element's 1-based node ids. After `"COORDONNEES DES
NOEUDS"`, `nnod` records of `node_index x1 x2 ... x_dim` follow; the leading
index is read and discarded (rows are assumed already in file order).

## Cell types

The `desc3` field selects the type:

| desc3 | type | desc3 | type |
|---|---|---|---|
| 2 | `vertex` | 11 | `tetra10` |
| 3 | `line` | 12 | `wedge` |
| 4 | `line3` | 13 | `wedge15` |
| 5 | `triangle` | 15 | `hexahedron` |
| 6 | `triangle6` | 16 | `hexahedron20` |
| 7 | `quad` | 17 | `pyramid` |
| 8 | `quad8` | 10 | `tetra` |

meshio++ → `(desc1, desc2, desc3)` on write:

| type | desc1, desc2, desc3 | type | desc1, desc2, desc3 |
|---|---|---|---|
| `vertex` | 1, 1, 2 | `tetra10` | 5, 15, 11 |
| `line` | 2, 2, 3 | `wedge` | 6, 207, 12 |
| `line3` | 2, 3, 4 | `wedge15` | 6, 307, 13 |
| `triangle` | 3, 7, 5 | `hexahedron` | 7, 2202, 15 |
| `triangle6` | 3, 7, 6 | `hexahedron20` | 7, 3303, 16 |
| `quad` | 4, 202, 7 | `pyramid` | 8, 4202, 17 |
| `quad8` | 4, 303, 8 | `tetra` | 5, 4, 10 |

Hybrid meshes are supported.

## Data mapping

- `cell_data["pf3:ref"]` — per-element region reference (from header field 3),
  one array per cell block.

## Quirks & limitations

- Unlike UNV/gmsh/mphtxt, **no node-order permutation table** is applied —
  node ids pass through in file order directly. This round-trips losslessly
  through meshio++ but is not guaranteed to match FLUX's own internal node
  ordering convention for every element type.
- Header-line detection is French-text substring matching; a real FLUX file
  with reworded headers (not expected in practice, but possible) would break
  parsing.
- Region *names* (which FLUX may store as binary data alongside the mesh) are
  not read at all — only the numeric per-element reference.
- Several always-placeholder header fields on write: region counts are always
  `1`/`0`/`0`/`0`/`0`/`0`, and both "max nodes per element" and "max
  integration points per element" fields are hardcoded to `20` regardless of
  actual mesh content.

## Notes

- Fully handled by the C++ core.
- No reference fixture exists under `tests/meshes/flux/`; tests round-trip
  every supported linear and second-order type.
