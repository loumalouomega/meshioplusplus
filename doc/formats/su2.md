# SU2 (`.su2`)

The [SU2](https://su2code.github.io/docs_v7/Mesh-File/) mesh format: an ASCII
format with `NDIME`, `NPOIN`, `NELEM` volume cells and `NMARK` boundary
markers, using VTK-style numeric type codes.

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

Line-oriented `KEY= value` records (`%` = comment; blank lines skipped;
malformed lines without `=` just warn and are skipped).

- `NDIME= 2|3` — dimension.
- `NPOIN= n [extra columns]` — the first point row is read separately to
  auto-detect whether an extra trailing global-index column is present
  (some SU2 files append one); the remaining `n-1` rows follow with the same
  column count, then the extra column(s) are stripped from all rows.
- `NELEM= n` or `MARKER_ELEMS= n` — reads exactly `n` element lines; the
  first line is peeked to determine the node count (from the VTK type code)
  and whether an extra trailing column is present, then all `n` lines are
  parsed as one integer block and binned by VTK type code (elements of
  different types within one block are separated into distinct cell
  blocks). `NELEM` cells get `su2:tag = 0`; `MARKER_ELEMS` cells get the
  current marker's tag id.
- `NMARK= n` — expected marker count (soft-checked; a mismatch only warns).
- `MARKER_TAG= <value>` — sets the tag used for the following
  `MARKER_ELEMS` block. If the tag isn't parseable as an integer, an
  auto-incrementing integer id is assigned instead, with a warning that
  string tags aren't supported.
- Boundary cell blocks of the **same type** from separate `MARKER_ELEMS`
  sections are merged into one block per type after the full-file scan,
  concatenating their `su2:tag` values in the same order.

Write: `NDIME=` from `points.shape[1]`; `NPOIN=` + coordinates; volume cells
(`triangle`/`quad` for 2D, `tetra`/`hexahedron`/`wedge`/`pyramid` for 3D)
under `NELEM=`, each row prefixed by its type code; boundary markers grouped
by the first integer-typed `cell_data` array found (via the same
"first-int-array" convention used by several other formats), one
`MARKER_TAG=`/`MARKER_ELEMS=` pair per distinct tag value.

## Cell types

| code | nodes | meshio type |
|---|---|---|
| 3 | 2 | `line` |
| 5 | 3 | `triangle` |
| 9 | 4 | `quad` |
| 10 | 4 | `tetra` |
| 12 | 8 | `hexahedron` |
| 13 | 6 | `wedge` |
| 14 | 5 | `pyramid` |

## Data mapping

- `cell_data["su2:tag"]` — volume cells always get tag `0`; boundary
  (`MARKER_ELEMS`) cells get their marker's tag id (an auto-incrementing int
  starting from 1 for the first non-numeric string tag encountered, if any).

## Quirks & limitations

- Unsupported cell types on write trigger a warning, but a latent message
  bug means the warning text shows Python's `type` builtin rather than the
  actual cell type name — cosmetic only, doesn't affect behavior.
- Only **one** integer cell-data array can be used as the boundary-marker
  tag source on write; if others exist, they're dropped with a warning.
- Point coordinates use `np.savetxt`'s default format in the Python writer
  (`%.18e`) vs. `%.16e` in the C++ writer — a formatting difference, not a
  round-trip correctness issue.

## Notes

- `tests/meshes/su2/square.su2` (official SU2-docs example: structured 2D
  quad mesh, 9 points/8 quads/markers) — checked for 16 total cells, 9
  points, 4 unique tags summing to 20.
- `tests/meshes/su2/mixgrid.su2` (mixed-cell 3D grid) — checked for 30
  cells, 16 points, 6 unique tags summing to 62.
- Fully handled by the C++ core.
