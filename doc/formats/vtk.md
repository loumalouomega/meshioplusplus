# VTK legacy (`.vtk`)

The [VTK legacy](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
file format (`UNSTRUCTURED_GRID`), versions **4.2** and **5.1**, in ASCII and
big-endian binary. Two independently-implemented sub-readers handle the two
versions.

| | |
|---|---|
| **Format name** | `vtk` (writes v5.1), `vtk42`, `vtk51` |
| **Extensions** | `.vtk` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.vtk")
meshioplusplus.vtk.write("out.vtk", mesh, binary=True)   # version via file_format
```

- **`binary`** — big-endian binary (`True`) or ASCII.
- Version: `file_format="vtk"`/`"vtk51"` writes 5.1 (default); `"vtk42"`
  writes 4.2.

Binary numeric data is **always big-endian** regardless of host platform —
an explicit VTK-wiki-documented convention, not a meshio++ choice.

## File structure

Line 1: `# vtk DataFile Version <ver>` (the version string's exact value —
`"5.1"` vs. anything else — is what selects the 4.2 vs. 5.1 sub-reader,
so this also works for genuinely older version strings). Line 2: title
(skipped). Line 3: `ASCII` or `BINARY`. Then `DATASET <TYPE>` where `TYPE` is
one of `UNSTRUCTURED_GRID`, `STRUCTURED_POINTS`, `STRUCTURED_GRID`,
`RECTILINEAR_GRID`.

**4.2 `CELLS` layout** — interleaved: `CELLS <num_cells> <total_ints>` then,
per cell, `[n, p0, p1, ..., p_{n-1}]` (int32, big-endian in binary mode);
followed by a separate `CELL_TYPES <num_cells>` section.

**5.1 `CELLS` layout** — no official published spec; reverse-engineered from
real files and a ParaView forum discussion:
```
CELLS <num_offsets> <num_conn_items>
OFFSETS <dtype>
<num_offsets values, first is 0, last equals len(connectivity)>
CONNECTIVITY <dtype>
<num_conn_items values>
```
`<dtype>` is a literal type token like `vtktypeint64`.

`POINT_DATA`/`CELL_DATA` sub-sections: `SCALARS <name> <type> [numComp]` +
`LOOKUP_TABLE <name>` (the lookup table itself is consumed but discarded);
`VECTORS <name> <type>` (3 components); `TENSORS <name> <type>` (3×3);
`FIELD FieldData <n>` (a list of `<name> <numComp> <numTuples> <type>`
blocks, each optionally preceded by a `METADATA` sub-block skipped up to the
next blank line).

Structured dataset types (`STRUCTURED_POINTS`/`STRUCTURED_GRID`/
`RECTILINEAR_GRID`) are converted to unstructured line/quad/hex cells,
generated in Fortran (column-major) point/cell order, matching VTK's own
convention for these types.

## Cell types

Shared with VTU — see [VTU](./vtu.md#cell-types) for the full numeric
type-code table. Only `wedge` needs a node-order permutation relative to VTK
(`[0,2,1,3,5,4]`, self-inverse); every other type uses natural order.

## Data mapping

Generic `SCALARS`/`VECTORS`/`TENSORS`/`FIELD` blocks map 1:1 to
`point_data`/`cell_data`; no reserved key names. `point_sets`/`cell_sets`
round-trip as extra data arrays (v5.1 only), same as VTU.

## Quirks & limitations

- **The C++ reader only supports `UNSTRUCTURED_GRID`** — any other
  `DATASET` type (structured points/grid, rectilinear grid) always falls
  back to Python.
- The 4.2 and 5.1 sub-readers use two genuinely different cell-reconstruction
  algorithms (4.2: list-based per-block append; 5.1: shared offset-diff/
  vectorized helper also used by VTU) — this is historical rather than
  deliberate, but means bugs in one don't necessarily affect the other.
- `_cpp_ok(mesh)` gate: the C++ path is skipped for meshes with polyhedron
  cells, or with any 2-component vector data — because the Python writer
  pads 2-component vectors to 3 components (**mutating the input mesh in
  place**), which the C++ writer deliberately does not replicate.
- `COLOR_SCALARS` sections are read and discarded (only to advance the file
  cursor correctly).
- The registered write-dict alias `vtk51` currently maps to the same
  function as `vtk42` in the format registry (both point at the 4.2 writer)
  — use the `fmt_version="5.1"` kwarg via `meshioplusplus.vtk.write` directly, or the
  default `file_format="vtk"`, to reliably get a 5.1 file.

## Notes

- `tests/meshes/vtk/00_image.vtk`/`01_image.vtk` (`STRUCTURED_POINTS`,
  generating 81/100 and 72/147 points/cells), `02_structured.vtk`
  (`STRUCTURED_GRID`), `03-05_rectilinear.vtk` (`RECTILINEAR_GRID`
  variants), `06_unstructured.vtk` (hexahedron, 12/42), `06_color_scalars.vtk`
  (5 points/2 cells, exercises `COLOR_SCALARS`), `gh-935.vtk` (triangle
  regression test), `rbc_001.vtk` (996 cells, a red-blood-cell mesh).
- The C++ core handles both versions in ASCII and big-endian binary for
  `UNSTRUCTURED_GRID` only.
