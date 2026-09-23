# Ansys MAPDL coded database (`.cdb` / `.inp`)

The coded database is the text archive Ansys MAPDL writes with `CDWRITE` and reads with `CDREAD`, and what Ansys Workbench, HyperMesh and other preprocessors export for MAPDL. meshio++ reads its mesh (nodes, elements, element types and components) and writes one MAPDL reads back. It is distinct from the [Fluent `.msh` format](ansys.md), which meshio++ calls `ansys`. Results live in MAPDL's binary [`.rst`](./ansys_rst.md) files.

| | |
|---|---|
| **Format name** | `ansysInp` |
| **Extensions** | `.cdb`, `.inp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.cdb")
mesh.regions                           # CMBLOCK components
mesh.cell_data["ansys:element"]        # element routine per cell (185, 186 ...)
meshioplusplus.write("out.cdb", mesh)

meshioplusplus.ansysInp.read("model.cdb", lenient=True)  # skip unmappable elements
```

```bash
meshioplusplus convert model.cdb model.vtu
```

**`.inp` is shared with [Abaqus](abaqus.md)**, which is imported first, so `meshioplusplus.read("x.inp")` reads Abaqus. Pass `file_format="ansysInp"` to read a MAPDL `.inp`.

## File structure

A `.cdb` file is a list of MAPDL commands. Four of them open blocks of fixed-width rows, and a Fortran format line after the block header gives the column widths:

```
ET,1,186                                  -- element type slot 1 is SOLID186
KEYOPT,1,2,0
NBLOCK,6,SOLID,     321,     321
(3i9,6e21.13e3)                           -- 3 integers of 9 columns, then reals of 21
        1        0        0 1.0000000000000E+000 0.0000000000000E+000 ...
N,R5.3,LOC,       -1,
EBLOCK,19,SOLID,      40,      40
(19i10)
         1         1         1         1         0         0         0         0        20         0         1        12 ...
        ...                               -- nodes past the eighth continue on the next line
        -1
CMBLOCK,FIXED,NODE,       2               -- a node component
(8i10)
         1       -40                      -- nodes 1 to 40
```

- **Format lines** are parsed as Fortran formats (`detail::parse_fortran_format`), so each field is cut at its own width: `(3i9,6e21.13e3)` from MAPDL, `(1i7,2i9,6e21.13)` from Workbench, `(3i8,6e16.9)` from HyperMesh, repeat counts and groups included. The same parser serves the LS-DYNA reader.
- **`ET`** (by number, `ET,1,186`, or by name, `ET,1,SOLID186`), **`ETBLOCK`** and **`KEYOPT`** (and its abbreviations, such as `KEYOP`) give each element type slot its routine number and key options.
- **`NBLOCK`** rows hold the node number and coordinates. Nodal rotation angles, when a row has them, are dropped with a warning.
- **`EBLOCK`** rows (the `SOLID` layout MAPDL writes) hold the material, type, real constant and section numbers, the node count, the element number and the nodes. A non-solid `EBLOCK` (which Workbench can write for contact elements) is skipped with a warning.
- **`CMBLOCK`** lists a component's node or element numbers, a negative value closing a run opened by the value before it. A block whose header count is too large ends at the next command.
- A `!` comment after a command is ignored. Other commands (`MP`, `R`, `SECTYPE`, `D`, `F` ...) are skipped.

## Cell types

An element becomes a cell by its routine's category, as in the open readers pymapdl-reader and [mapdl-archive](https://github.com/akaszynski/mapdl-archive):

| Category | Routines (examples) | meshio++ |
|---|---|---|
| point | MASS21, 71, 175 | `vertex` |
| line | LINK180, BEAM189, PIPE289, SURF151/153 ... | `line`, or `line3` when the third node is set |
| linear line | BEAM188, 214, 216, 217 (the extra nodes orient the beam) | `line` |
| shell or plane | SHELL181/281, PLANE182/183, SURF152/154, TARGE170, CONTA174 ... | `quad`/`quad8`; `triangle`/`triangle6` when K = L |
| brick | SOLID45, 185, 186, SOLSH190, 272, 273 ... | `hexahedron`/`hexahedron20`, or a degenerate form |
| tetrahedron | SOLID92, 187, 285 ... | `tetra`/`tetra10` |

MESH200's shape comes from its KEYOPT(1). **Degenerate bricks** are resolved by their repeated nodes: a hexahedron when O ≠ P, a wedge when K = L and O = P, a pyramid when M = N = O = P, and a tetrahedron when also K = L. A shell with K = L is a triangle. A shell of five nodes (SURF152's extra orientation node) is linear, and so is an 8-node shell whose last two midside nodes are equal (a contact face written without midsides).

**Missing midside nodes** are written as node `0`, or left off the end of a row (SOLID92 rows sometimes carry nine nodes). meshio++ creates each one at its edge midpoint, shared between the elements on that edge, with one warning. An element whose type has no category, or whose corner node is `0` (a TARGE170 pilot node), is a `ReadError` naming it; read with `lenient=True` to skip such elements with a warning.

## Data mapping

| `.cdb` | meshio++ |
|---|---|
| element routine, type slot, `MAT`, `REAL`, `SECNUM` | `cell_data["ansys:element"]`, `"ansys:type"`, `"ansys:mat"`, `"ansys:real"`, `"ansys:secnum"` |
| `CMBLOCK` `NODE` / `ELEM` component | `point` / `cell` region (so also `point_sets` / `cell_sets`) |

Cells are grouped into one block per meshio++ type, in the order the types first appear.

## Writing

The writer emits `/PREP7`, one `ET` per element type slot, an `NBLOCK` in MAPDL's `(3i9,6e21.13e3)` layout, a `SOLID` `EBLOCK` in `(19i10)`, one `CMBLOCK` per point or cell region (runs of consecutive numbers packed as `first, -last`) and `FINISH`. Nodes and elements are numbered from 1 in mesh order. Each type is written in the layout MAPDL expects:

| meshio++ | written as |
|---|---|
| `hexahedron`, `wedge`, `pyramid` | SOLID185, the wedge and pyramid as degenerate bricks |
| `hexahedron20`, `wedge15`, `pyramid13` | SOLID186, likewise |
| `tetra` / `tetra10` | SOLID285 / SOLID187 |
| `quad`, `triangle` / `quad8`, `triangle6` | SHELL181 / SHELL281, the triangle as a degenerate quad |
| `line` / `line3` / `vertex` | BEAM188 / BEAM189 / MASS21 |

`ansys:element` and `ansys:type` are kept when the cell's type fits that routine's layout, so a MAPDL model keeps its element types through a round trip. Otherwise, with a warning, the default above is written. `ansys:mat`, `ansys:real` and `ansys:secnum` are written back, defaulting to 1. Side regions have no component equivalent and are dropped with a warning and a provenance note. Polygons, polyhedra and other types without a MAPDL element are a `WriteError`. The provenance block is written as `!` comments. Both engines write the same bytes.

## Notes

- **Validation.** The test fixtures in `tests/python/meshes/ansys/` are decks from mapdl-archive (MIT) that MAPDL, Workbench and HyperMesh wrote. Both engines read them into the cells mapdl-archive builds, compared as sets of node coordinates, and the components mapdl-archive finds. The one difference: mapdl-archive puts three missing triangle midside nodes at the origin, where meshio++ puts them at their edge midpoints. MAPDL itself was not available, so a written deck is checked by reading it back, not by `CDREAD`.
- `AnsysInfo`, the C++ side channel that carried components before v16.3.0, is still filled on read and still written, for code that uses it. Regions are the primary route.
