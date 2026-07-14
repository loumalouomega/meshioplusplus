# Ansys MAPDL coded database (`.cdb` / `.inp`)

An autonomous reader/writer for the Ansys **MAPDL "coded database"** format —
distinct from the [Fluent `.msh` format](ansys.md) also named "ansys" in
meshio. It parses `ET`/`ETBLOCK`, `NBLOCK`, `EBLOCK`, and `CMBLOCK` blocks
directly, with no dependency on any other format module.

| | |
|---|---|
| **Format name** | `ansysInp` |
| **Extensions** | `.cdb`, `.inp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.ansysInp.read("model.cdb")
meshio.ansysInp.write("out.cdb", mesh)
```

Both `read(filename)` and `write(filename, mesh)` take no keyword arguments.

**Note on `.inp`**: this format registers **both** `.cdb` and `.inp` as
extensions, colliding with [Abaqus](abaqus.md)'s pre-existing `.inp`
registration. Since `abaqus` is imported before `ansysInp` in
`src/meshio/__init__.py`, plain extension-based dispatch
(`meshio.read("x.inp")`) still resolves to **Abaqus** by default — pass
`file_format="ansysInp"` explicitly, or call `meshio.ansysInp.read`/`write`
directly, to select this format for a `.inp` file.

## File structure

Whitespace/keyword-delimited MAPDL command blocks, each starting with a
header line and ending at a sentinel row:

```
/PREP7
ET,<slot>,<ansys element type>              -- one per cell type in use
NBLOCK,6,SOLID,<n>,<n>
(3i9,6e20.13)                               -- format spec: int width, float width
<id><cs><...>  x  y  z                      -- one row per node, fixed-width fields
N,R5.3,LOC,      -1,                        -- NBLOCK terminator
EBLOCK,19,SOLID,<ntot>,<ntot>
(19i9)
<elem header (11 ints)><node ids...>        -- continuation line if >8 nodes
        -1                                  -- EBLOCK terminator
CMBLOCK,<name>,NODE,<n>                     -- named point set
(8i10)
<ids, possibly negative range markers>
CMBLOCK,<name>,ELEM,<n>                     -- named cell set
(8i10)
<ids>
FINISH
```

- **Field widths** are read from the format-spec line following each block
  header (`_int_width`/`_real_width`, parsed via regex against patterns like
  `3i9` / `6e20.13`), not hardcoded — but the writer always emits `i9`/`e20.13`
  widths.
- **`ETBLOCK`** (an alternative to individual `ET,` lines) associates a
  numeric element-type slot with an underlying Ansys element type id; either
  form populates the same `etype_lib` mapping used to classify `EBLOCK` rows.
- **`EBLOCK`** rows: fields are `(mat, type, real, secnum, esys, birth, death,
  solkey, nodes_per_elem, ..., elem_id, node_ids...)`; only fields 1 (etype
  slot), 8 (node count) and 10 (element id) are used, plus however many
  trailing node ids follow (spilling onto a continuation line if there are
  more than 8).
- **`CMBLOCK`** (named component = point/cell set): entity `NODE` → point set,
  entity starting `ELEM` → cell set. Negative values encode a **range
  marker**: `-k` after a base value `b` expands to `range(b+1, k+1)` — an
  `ReadError` is raised if a negative value appears before any base value.
- Any line matching `_is_data_line`'s exclusion list (known keywords,
  `KEYWORD,` command syntax, `!`/`/` comments) is treated as non-data and
  stops a block's row-reading loop early.

## Cell types

Ansys element type ids are grouped into 4 families by node-count-independent
type id, then combined with the node count actually present to resolve a
meshio type:

| family | Ansys element type ids |
|---|---|
| `solid` | 5, 45, 70, 87, 90, 92, 95, 162, 185, 186, 187, 226, 227, 285 |
| `shell` | 28, 43, 63, 93, 131, 132, 181, 281 |
| `plane` | 25, 42, 77, 82, 182, 183, 223 |
| `line` | 1, 3, 4, 21, 180, 188, 189, 288, 289 |

| (family, nodes) | meshio | (family, nodes) | meshio |
|---|---|---|---|
| (solid, 4) | `tetra` | (shell/plane, 3) | `triangle` |
| (solid, 10) | `tetra10` | (shell/plane, 6) | `triangle6` |
| (solid, 8) | `hexahedron` | (shell/plane, 4) | `quad` |
| (solid, 20) | `hexahedron20` | (shell/plane, 8) | `quad8` |
| (solid, 6) | `wedge` | (line, 2) | `line` |
| (solid, 15) | `wedge15` | (line, 3) | `line3` |
| (solid, 5) | `pyramid` | | |
| (solid, 13) | `pyramid13` | | |

Write uses a fixed reverse mapping (one Ansys type id per meshio type,
regardless of which id the file was originally read with):
`tetra→285, tetra10→187, hexahedron→185, hexahedron20→186, wedge→185,
wedge15→186, pyramid→185, pyramid13→186, triangle→181, triangle6→281,
quad→181, quad8→281, line→188, line3→189`. An unmapped meshio cell type
raises `WriteError`.

## Data mapping

- `mesh.point_sets[name]` — from `CMBLOCK ...,NODE,...`, 0-based point
  indices.
- `mesh.cell_sets[name]` — from `CMBLOCK ...,ELEM,...`, one array per cell
  block (in the order blocks were first encountered), 0-based local indices.
- No point_data, cell_data, or field_data — only geometry, connectivity, and
  named sets are represented.

## Quirks & limitations

- 2D input meshes are padded to 3D with a zero z-column on write (MAPDL has
  no native 2D coordinate concept).
- The writer always emits exactly one `NBLOCK`/`EBLOCK` pair (no attempt to
  preserve an original file's exact block layout, field widths, or element
  type ids) — a read→write round trip is not byte-identical, though it is
  semantically equivalent (same points, cells, sets).
- Read and write go through the C++ core (`meshio._core.ansysinp_read`/
  `ansysinp_write`), with the Python reference as an automatic fallback for
  buffers. CMBLOCK components (`point_sets`/`cell_sets`) travel through a
  dedicated `AnsysInfo` side-channel struct rather than the Mesh conversion
  layer.

## Notes

- No dedicated `tests/meshes/` reference fixture; `tests/test_ansysInp.py`
  builds MAPDL snippets and meshes inline (including negative-range
  `CMBLOCK` expansion, multi-line `EBLOCK` continuation, and round-trips
  through `io.StringIO`/temp files) rather than shipping a binary fixture.
  The internal-helper tests import the Python module directly, so the pure-
  Python reference stays exercised alongside the C++ path.
- Ported from [Simvia's meshlane fork](https://github.com/simvia-tech/meshlane)
  (see `CHANGELOG.md`) — this format did not exist upstream before that.
