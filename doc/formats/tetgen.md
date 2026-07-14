# TetGen (`.node` / `.ele`)

The [TetGen](https://wias-berlin.de/software/tetgen/fformats.html) mesh
format: a pair of sibling files sharing a stem — `<name>.node` (points,
attributes, boundary markers) and `<name>.ele` (tetrahedra, region
attributes).

| | |
|---|---|
| **Format name** | `tetgen` |
| **Extensions** | `.node`, `.ele` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.node")   # reads the .node/.ele pair
meshioplusplus.tetgen.write("out.node", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format (default `".16e"`; the C++ fast path
  is only used at this default).

Either path (`.node` or `.ele`) selects the sibling pair.

## File structure

**`.node`**: header `npoints dim nattrs nbmarkers` (`dim` must be 3), then
rows `idx x y z attr1..attrN marker1..markerM`. The **node index base**
(0 or 1) is auto-detected from the first row's `idx` value, and all indices
are then required to be **exactly consecutive** from that base (`ReadError`
otherwise).

**`.ele`**: header `ntets 4 nattrs`, then rows `idx n0 n1 n2 n3
attr1..attrK`. Connectivity is shifted by the `.node` file's detected index
base (not hardcoded), so files using either 0- or 1-based numbering read
correctly.

Write: for `.node`, attribute/marker keys are partitioned into at most one
"ref" key (the first `point_data` key containing the substring `":ref"`, or
else the first key present) plus the remaining keys as plain attributes;
header comments record the attribute/marker names. Ref columns use plain
`{}` (`str()`) formatting in the Python writer, not the float format string —
the C++ writer instead special-cases exact-integer values to print as plain
integers, falling back to `%.16e` otherwise, a refinement that can format
float-valued refs slightly differently between the two backends. For
`.ele`, the same "first `:ref`-containing key floats to the front" rule
applies to `cell_data` keys, and **each `tetra` cell block gets its own
header line with its own index counting restarting at 0** — a mesh with
multiple `tetra` blocks would therefore produce duplicate element ids across
blocks (TetGen expects globally unique ids).

## Cell types

`tetra` only — TetGen only ever represents tetrahedra by construction.

## Data mapping

- `point_data["tetgen:attr{k}"]` — the k-th node attribute column
  (`tetgen:attr1`, `tetgen:attr2`, ...).
- `point_data["tetgen:ref"]` / `"tetgen:ref2"` / ... — boundary marker
  columns.
- `cell_data["tetgen:ref"]` / `"tetgen:ref2"` / ... — region attribute
  columns (a single-element list, since TetGen only has one cell block).

## Quirks & limitations

- The format spans two files and **cannot** be read from or written to a
  buffer.
- Node index base auto-detection plus the consecutive-numbering check means
  a `.node` file with gaps in its index sequence is rejected outright, not
  silently accepted with holes.
- Writing multiple `tetra` blocks resets the element-id counter per block —
  a genuine round-trip risk for multi-block meshes (uncommon in practice,
  since TetGen conventionally produces exactly one tetrahedron block).

## Notes

- `tests/meshes/tetgen/mesh.node` (89 points, 1 attribute column named
  `"moje_data"` + 1 boundary-marker column named `"medit:ref"`) and
  `mesh.ele` (304 tetrahedra, 1 region attribute `"medit:ref"`) — checked via
  `mesh.point_data["tetgen:ref"].sum() == 12` and
  `mesh.cell_data["tetgen:ref"][0].sum() == 373`.
- Fully handled by the C++ core (Python fallback only for buffers or a
  non-default `float_fmt`).
