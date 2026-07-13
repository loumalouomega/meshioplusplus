# FLAC3D (`.f3grid`)

The [Itasca FLAC3D](https://www.itascacg.com/software/flac3d) grid format
(`.f3grid`): ASCII or binary, with separate ZONE (3D) and FACE (2D) sections
and named cell groups.

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

- **`float_fmt`** — coordinate format (ASCII only).
- **`binary`** — binary (`True`) or ASCII (`False`, default).

## File structure

**Format auto-detection**: the first 8 bytes are checked for a null byte; if
found, the file is binary, else ASCII.

**Binary** (little-endian): 8 bytes of a header pair whose meaning is
undocumented on read (`_read`'s comment: "not sure what the first bytes
represent") but whose exact values (`1375135718, 3`) **are** reproduced by
the writer as a magic number; then `uint32` node count, per-node
`(point_id: uint32, x,y,z: float64×3)`; then for `zone` and `face` in that
order: `uint32` cell count, per-cell `(cell_id: uint32, num_verts: uint32,
node_ids: uint32×num_verts)` — if `num_verts == 7` (a degenerate "B7"
hexahedron-as-7-node encoding), the last node is duplicated to make 8; then
`uint32` group count and, per group, length-prefixed name/slot strings plus
a `uint32`-counted id list.

**ASCII**: `G <id> <x> <y> <z>` per point; `Z <TypeAbbrev> <id> <n0> <n1>
...` per zone (`F` for faces), with the same `"B7"` degenerate-hex handling;
`ZGROUP "<name>" SLOT <slot>` / `FGROUP "<name>" SLOT <slot>` followed by
whitespace-separated member id lines until a line starting `*`, `ZGROUP`, or
`FGROUP`.

**Cell-block grouping**: cells are grouped into blocks of **consecutive
same-typed cells in file order** — not all cells of one type merged
globally — so a file alternating types produces alternating blocks (see the
reference-file cell list below for a concrete example).

**Right-handed zone reorder** (the format's central quirk): FLAC3D requires
each zone's first four corner nodes to form a right-handed system. For every
zone cell, meshio computes the scalar triple product of the first three edge
vectors (from the "primary" meshio→FLAC3D order's first four nodes); if
positive, the primary order is used, otherwise a "flipped" alternate order.
This check only happens **on write** — the read-side reorder is a fixed,
unconditional permutation, since a well-formed file is assumed to already
store correctly-handed zones.

## Cell types & node ordering

| meshio type | FLAC3D abbrev | primary order | flipped order (write-only) |
|---|---|---|---|
| `triangle` | `T3` | `[0,1,2]` | — |
| `quad` | `Q4` | `[0,1,2,3]` | — |
| `tetra` | `T4` | `[0,1,2,3]` | `[0,2,1,3]` |
| `pyramid` | `P5` | `[0,1,3,4,2]` | `[0,3,1,4,2]` |
| `wedge` | `W6` | `[0,1,3,2,4,5]` | `[0,2,3,1,5,4]` |
| `hexahedron` | `B8` | `[0,1,3,4,2,7,5,6]` | `[0,3,1,4,2,5,7,6]` |

Read-side order (FLAC3D → meshio): `triangle/quad/tetra` unchanged;
`pyramid: [0,1,4,2,3]`; `wedge: [0,1,3,2,4,5]`; `hexahedron:
[0,1,4,2,3,6,7,5]`.

## Data mapping

- `cell_data["cell_ids"]` — the original FLAC3D global cell id, split per
  block (faces are numbered first, then zones, in the reader's internal
  concatenation order).
- `mesh.cell_sets[label]` — from `ZGROUP`/`FGROUP`, keyed by the group name.

## Quirks & limitations

- **Read/write section-order asymmetry**: the writer emits `* ZONES` before
  `* FACES` in the ascii file, but the reader concatenates faces-then-zones
  internally when building `Mesh.cells` — so the returned cell-block order
  after a read does not match the on-disk section order. This is harmless
  (the reader doesn't depend on write order) but a genuine structural
  asymmetry worth knowing about when comparing a file's raw section layout
  to `mesh.cells`.
- `ZGROUP`/`FGROUP` group labels written out are hardcoded to `SLOT 1`
  (ascii) / a fixed `"Default"` slot string (binary) — the original slot
  name from a read file is not preserved on a subsequent write.
- The C++ reader **throws immediately** if it encounters any `ZGROUP`/
  `FGROUP`/binary-group section, deferring the whole file to Python — this
  matches the shim's write-side gate (C++ write is only attempted when
  `mesh.cell_sets` is empty).

## Notes

- `tests/meshes/flac3d/flac3d_mesh_ex.f3grid` (ascii, "FLAC3D 7.00 Release
  118") and `flac3d_mesh_ex_bin.f3grid` (its binary counterpart) — checked
  for point-sum ≈307.0 and an exact ordered list of 12 alternating cell
  blocks: `[(quad,15), (triangle,3), (hexahedron,45), (pyramid,9),
  (hexahedron,18), (wedge,9), (hexahedron,6), (wedge,3), (hexahedron,6),
  (wedge,3), (pyramid,6), (tetra,3)]` — a direct illustration of the
  "consecutive same-type runs" block-splitting rule above — plus 5 named
  cell_sets (`Brick1`, `Pyramid2`, `Tetrahedron4`, `Wedge3`, and an
  `FGROUP "bottom"`).
- The C++ core handles points and zone/face cells in both ASCII and binary
  (with the determinant-based reorder); cell groups always defer to Python.
