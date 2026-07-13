# I-DEAS Universal — UNV (`.unv`)

The [I-DEAS Universal File](https://www.ceas3.uc.edu/sdrluff/) format is an ASCII
interchange format used by SDRC I-DEAS, Salome, Code-Aster and many FE tools. A
file is a sequence of **datasets**, each delimited by a line containing `-1`,
followed by a 6-digit dataset id and its records.

| | |
|---|---|
| **Format name** | `unv` |
| **Extensions** | `.unv` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.unv")          # or meshio.unv.read("mesh.unv")
meshio.unv.write("out.unv", mesh)
```

`write` takes no keyword arguments.

## Datasets handled

- **2411** (and legacy **781**) — *nodes*: a `(4I10)` record (label, coordinate
  systems, colour) followed by a `(1P3D25.16)` coordinate record. Node labels
  are remapped to contiguous 0-based indices.
- **2412** — *elements*: a `(6I10)` record (label, **FE descriptor id**,
  physical- and material-property ids, colour, node count) followed by the
  connectivity. Beam-type elements carry an extra orientation record, which is
  written and skipped automatically.
- **2467 / 2477** — *permanent groups*: mapped to `point_sets` / `cell_sets`.

## Element types & node ordering

The FE descriptor id selects the element type:

| descriptor(s) | meshio type |
|---|---|
| 11, 21 | `line` |
| 22, 24 | `line3` |
| 41, 81, 91 | `triangle` |
| 42, 82, 92 | `triangle6` |
| 44, 84, 94, 122 | `quad` |
| 45, 85, 95 | `quad8` |
| 111 | `tetra` |
| 118 | `tetra10` |
| 112 | `wedge` |
| 115 | `hexahedron` |
| 116 | `hexahedron20` |

Parabolic (second-order) elements use the **Salome/Code-Aster mid-node
"sandwich" ordering** (`corner, mid, corner, …`), which differs from meshio's
"corners first, then edge nodes" convention. meshio applies the permutation on
read and its inverse on write for `line3`, `triangle6`, `quad8`, `tetra10` and
`hexahedron20`.

## Data mapping

- Per-element physical-property id → `cell_data["unv:pid"]`.
- Named groups (dataset 2467): node groups → `point_sets`, element groups →
  `cell_sets`.

## Notes

- The C++ core handles nodes and elements (datasets 2411/2412) with the
  parabolic reordering and property ids. Files containing groups (2467) or field
  datasets (2414/55/…) transparently fall back to the Python reader, and meshes
  carrying `point_sets`/`cell_sets` fall back to the Python writer.
- This reader/writer was implemented against the
  [FEconv](https://github.com/victorsndvg/FEconv) format documentation. Node
  orderings round-trip losslessly but may differ from a specific tool's internal
  ordering for some higher-order elements.
