# I-DEAS Universal — UNV (`.unv`)

The [I-DEAS Universal File](https://www.ceas3.uc.edu/sdrluff/) format is an
ASCII interchange format used by SDRC I-DEAS, Salome, Code-Aster and many FE
tools. A file is a sequence of **datasets**, each delimited by a line
containing only `-1`, followed by a dataset-id line and the dataset's records.

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

Both `read`/`write` take no keyword arguments.

## File structure

```
    -1
  2411
<node records>
    -1
    -1
  2412
<element records>
    -1
    -1
  2467
<group records>
    -1
```

Datasets are split by scanning for lines that strip to exactly `"-1"`; the
line right after such a delimiter is the numeric dataset id, and everything up
to the next `-1` is the dataset body.

### 2411 / 781 — nodes

Two-line records: `label CS1 CS2 color` (only `label`, field 0, is used),
followed by a coordinate line. Coordinates may use Fortran exponent notation
(`D`/`d` instead of `E`/`e`), which is normalized before `float()` parsing.
Node labels are arbitrary integers (not necessarily 1..N or contiguous); a
`label → 0-based index` map is built while reading so later datasets (elements,
groups) can resolve them.

### 2412 — elements

Record 1 is 6 integers: `label fedesc pid ... ... num_nodes` (fields 0, 1, 2,
5 are used; fields 3-4 are ignored). If `fedesc` is one of the four beam
descriptors (`11, 21, 22, 24`), an extra 3-integer orientation record follows
immediately and is **discarded on read, always rewritten as `0 0 0` on
write** — beam orientation data does not round-trip through meshio.
Node-label records follow, gathered (possibly across several lines) until
`num_nodes` labels have been collected.

### 2467 / 2477 — permanent groups

Record 1 has ≥8 integers; field 7 is the entity count `n`. The next line is
the group name. Then `4*n` integers follow, laid out as `(entity_type, tag, 0,
0)` quadruples (typically 2 quadruples — 8 integers — per line). Entity type
`8` = node → the group becomes a `point_sets` entry; entity type `7` = element
→ `cell_sets`. Other entity types are collected internally but never emitted.

## Cell types & node ordering

The **FE descriptor id** (record-1 field 1 of a 2412 dataset) selects the
element type:

| descriptor(s) | meshio type | descriptor(s) | meshio type |
|---|---|---|---|
| 11, 21 | `line` | 111 | `tetra` |
| 22, 24 | `line3` | 118 | `tetra10` |
| 41, 81, 91 | `triangle` | 112 | `wedge` |
| 42, 82, 92 | `triangle6` | 115 | `hexahedron` |
| 44, 84, 94, 122 | `quad` | 116 | `hexahedron20` |
| 45, 85, 95 | `quad8` | | |

On write, one canonical descriptor is chosen per meshio type: `line→21`
(beam), `line3→24` (beam), `triangle→91`, `triangle6→92`, `quad→94`,
`quad8→95`, `tetra→111`, `tetra10→118`, `wedge→112`, `hexahedron→115`,
`hexahedron20→116`.

Parabolic (second-order) elements use the **Salome/Code-Aster mid-node
"sandwich" ordering** — corner, mid-node, corner, mid-node, … — which differs
from meshio's "all corners, then all edge nodes" convention. The full
permutation table (`meshio_position[i] = unv_position[table[i]]` on read;
inverse applied on write) is:

| type | permutation |
|---|---|
| `line3` | `[0, 2, 1]` |
| `triangle6` | `[0, 3, 1, 4, 2, 5]` |
| `quad8` | `[0, 4, 1, 5, 2, 6, 3, 7]` |
| `tetra10` | `[0, 4, 1, 5, 2, 6, 7, 8, 9, 3]` |
| `hexahedron20` | `[0, 8, 1, 9, 2, 10, 3, 11, 12, 13, 14, 15, 4, 16, 5, 17, 6, 18, 7, 19]` |

## Data mapping

- `cell_data["unv:pid"]` — the element's physical/material property id
  (record-1 field 2 of the 2412 dataset); defaults to `1` on write if absent.
- `point_sets` / `cell_sets` — from 2467/2477 groups, keyed by group name
  (node groups → `point_sets`, element groups → `cell_sets`).

## Quirks & limitations

- Beam-type orientation data (the extra record after beam elements) is read
  and discarded; the writer always emits a placeholder `0 0 0` record — this
  is a genuinely lossy round-trip for beam-oriented meshes.
- Node/element labels need not be contiguous or 1-based; meshio maintains
  explicit label→index maps throughout so groups referencing raw labels
  resolve correctly regardless of numbering gaps.
- The C++ reader fully implements datasets **2411 and 2412** (nodes,
  elements, with the sandwich reorder and property ids). Any file containing
  a **2467/2477 group dataset**, or a **field dataset** (2414/55/56/57), makes
  the C++ reader throw and the shim falls back to the Python reader.
- The C++ writer defers to Python whenever the mesh carries `point_sets` or
  `cell_sets` (since groups have no C++ writer path).

## Notes

- Implemented against the [FEconv](https://github.com/victorsndvg/FEconv)
  format documentation and its FE-descriptor table.
- No reference fixture exists under `tests/meshes/unv/`; tests round-trip
  synthetic meshes covering every supported linear and parabolic type, plus an
  explicit `test_groups` exercising `point_sets`/`cell_sets`.
