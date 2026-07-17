# Ansys / Fluent mesh (`.msh`)

The Ansys Fluent `.msh` mesh format: fully parenthesis-nested "Scheme-like" sections, in ASCII, binary, or a mix of both within one file.

| | |
|---|---|
| **Format name** | `ansys` |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.msh", file_format="ansys")
meshioplusplus.ansys.write("out.msh", mesh, binary=True)
```

- **`binary`** — write binary (`True`) or ASCII section bodies.

## File structure

Every section is `(<index> ...)`; the file is opened and always read in binary mode, since ASCII and binary sections may be mixed in one file. The index may be a bare decimal (ASCII payload) or prefixed `20`/`30` (binary payload — `20xx` = float32/int32, `30xx` = float64/int64).

- `(0 ...)` comment, `(1 "...")` header, `(2 <dim>)` dimensionality — bracket- skipped, not otherwise interpreted.
- **Nodes** — `(<pfx>10 (zone-id first last type ND) ( <data> ))`: header integers are **hexadecimal**; `first`/`last` give the point-count range, `ND` is the spatial dimension. Body is `(last-first+1) × ND` values — ASCII one point per line, binary a raw float32/float64 block. A self-contained line (equal `(` and `)` counts, no body) is a pure declaration and is skipped.
- **Cells** — `(<pfx>12 (zone-id first last zone-type element-type) ( <data> ))`: `zone-type == 0` marks a **dead zone**, skipped with no cells produced. `element-type` selects a fixed node count (table below); `mixed` (element-type 0) zones are parsed structurally but their body is **not** decoded at all (Fluent's own mixed-cell-type encoding is unresolved in this reader).
- **Faces** — `(<pfx>13 (zone-id first last type element-type) ( <data> ))`: body rows are `n0 n1 ... cr cl` (face nodes plus 2 adjacent-cell ids, the latter discarded). ASCII `mixed` (element-type 0) faces **are** parsed, each row prefixed by its own per-row type index; **binary mixed faces are not supported** and raise `ReadError`. Faces are folded into the same flat cell list as volume cells — in the resulting `Mesh.cells`, boundary faces and interior cells are only distinguished by `cell.type`.
- `(39 ...)` and `(45 ...)` (zone specifications) are always warned-and- skipped, never parsed for content.

All cells/faces across every zone are collected into one flat list, then every connectivity array has the **first point-zone's `first` index** subtracted, so files whose node numbering doesn't start at 0 or 1 still normalize correctly.

Write emits, in order:

```
(1 "meshio++ VERSION")
(2 DIM)
(10 (0 1 N_POINTS 0))                          -- node-count declaration
(12 (0 1 TOTAL_CELLS 0))                       -- cell-count declaration
(10|3010 (1 1 N 1 DIM)( DATA ))                -- node block
(12|2012|3012 (1 FIRST LAST 1 ANSYS_TYPE)( DATA ))   -- per cell block
```

## Cell types

Volume/element-type codes (`element-type` field of a `12` section):

| code | meshio++ type | nodes |
|---|---|---|
| 0 | mixed (unhandled) | — |
| 1 | `triangle` | 3 |
| 2 | `tetra` | 4 |
| 3 | `quad` | 4 |
| 4 | `hexahedron` | 8 |
| 5 | `pyramid` | 5 |
| 6 | `wedge` | 6 |

Face-type codes (`13` sections, read-only — not used on write): `0`=mixed, `2`=`line`(2), `3`=`triangle`(3), `4`=`quad`(4).

meshio++ → Ansys type codes on write: `triangle:1, tetra:2, quad:3, hexahedron:4, pyramid:5, wedge:6` (no writer support for `mixed`/polyhedral).

## Data mapping

None — `point_data`/`cell_data`/`field_data` are always empty; this format carries geometry and zone/boundary structure only.

## Quirks & limitations

- All connectivity and zone-header integers are **hexadecimal**, in both ASCII bodies and headers — the defining quirk of this format among the ones meshio++ supports.
- Binary vs. ASCII is signalled purely by an optional `"20"`/`"30"` prefix glued onto the section-index digits (e.g. `2010` = binary float32 nodes, `3012` = binary int64 cells).
- Dead zones (`zone-type == 0`) produce no cells at all.
- `mixed` cell zones (element-type 0) are structurally skipped — Fluent's own encoding for heterogeneous cell zones is not decoded, and no cells result from them.
- The C++ reader defers **any** face section (`13`) carrying a data body to the Python fallback — meaning any Fluent `.msh` with real boundary face zones (a very common real-world case) is always parsed by Python, not C++.
- 2D/3D validity (`dim in {2,3}`) is only checked on write, not on read.

## Notes

- No reference fixture exists under `tests/meshes/ansys/`; tests round-trip synthetic meshes (`empty_mesh`, `tri_mesh`, `tri_mesh_2d`, `quad_mesh`, `tri_quad_mesh`, `tet_mesh`, `hex_mesh`, `pyramid_mesh`, `wedge_mesh`), parametrized over both ASCII and binary.
- `.msh` is shared with [`gmsh`](./gmsh.md) and [`freefem`](./freefem.md); on auto-detection `ansys` is tried first. Pass `file_format` to disambiguate.
