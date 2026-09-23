# Code_Aster native mesh (`.mail`)

The `.mail` file is [Code_Aster](https://code-aster.org)'s own ASCII mesh, the one `LIRE_MAILLAGE(FORMAT='ASTER')` reads and that Salome-Meca and Code_Aster's `PRE_GMSH`/`PRE_IDEAS` converters write. It is a sequence of **blocks**, each opened by a keyword and closed by `FINSF`, and the file ends at `FIN`. Nodes and elements are identified by **name** (`N1`, `M1`, `NO123`…), not by number, and groups refer to those names.

| | |
|---|---|
| **Format name** | `code_aster` |
| **Extensions** | `.mail` (also recognised by content: a leading `TITRE` or `COOR_nD`) |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.mail")        # points, cells, GROUP_MA/GROUP_NO as regions
meshioplusplus.write("out.mail", mesh)          # or meshioplusplus.code_aster.write(...)
meshioplusplus.write("out.med", mesh)           # the same mesh as MED, for Salome
```

Neither `read` nor `write` takes extra options. Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes.

## Syntax

The reader follows Code_Aster's own reader (`lrmast.F90` and the routines it calls) rather than a looser reading of the manual:

- Tokens are separated by blanks or commas. `%` starts a comment that runs to the end of the line. Keywords are case-insensitive, and names are case-sensitive.
- **Only the first 80 columns of a line are read**, as in Code_Aster (`lirlig.F90` reads `A80`). A line with content past column 80 gets one warning, and that content is ignored.
- A record is a **token stream**. An element's node names, or a node's coordinates, may continue on the next lines. Code_Aster's own writer puts an element name and 8 nodes on a line, then 8 per continuation line.
- Block-header options (`NOM=`, `NBOBJ=`, `NUMIN=`, `AUTEUR=`, `DATE=`, …) may follow the keyword on the same line or on the next ones. `KEY = VALUE` with spaces is read as `KEY=VALUE`. Only `NOM=` is used, as a group's name.
- Coordinates are Fortran reals or integers (`1`, `1.`, `.5`, `1.E+0`, `10.D-1`). `inf`, `nan` and other spellings are errors.

| Keyword | Maps to |
|---|---|
| `COOR_1D` / `COOR_2D` / `COOR_3D` | the points, with a point dimension of 1, 2 or 3 |
| `POI1` | `vertex` |
| `SEG2` / `SEG3` / `SEG4` | `line` / `line3` / `line4` |
| `TRIA3` / `TRIA6` / `TRIA7` | `triangle` / `triangle6` / `triangle7` |
| `QUAD4` / `QUAD8` / `QUAD9` | `quad` / `quad8` / `quad9` |
| `TETRA4` / `TETRA10` | `tetra` / `tetra10` |
| `PENTA6` / `PENTA15` / `PENTA18` | `wedge` / `wedge15` / `wedge18` |
| `PYRAM5` / `PYRAM13` | `pyramid` / `pyramid13` |
| `HEXA8` / `HEXA20` / `HEXA27` | `hexahedron` / `hexahedron20` / `hexahedron27` |
| `GROUP_MA` | a `cell` region; `dim` is the highest dimension of its elements |
| `GROUP_NO` | a `point` region |
| `TITRE`, `DUMP`, `DEBUG`, `GROUP_FA`, `SYS_UNIT`, `SYS_COOR`, `MACRO_AR`, `MACRO_FA`, `MACRO_EL`, `MATERIAU` | skipped silently, as Code_Aster does |
| anything else (`SEG22`, `TR3QU4`, a typo) | skipped with one warning per keyword |

Each element block becomes one cell block, in file order. A group's name is its `NOM=` value, or else the block's first token. Groups have no number, so regions carry no `tag`. A name used for both a `GROUP_MA` and a `GROUP_NO` gives two regions. Two blocks with the same name and kind are merged, with a warning.

## Node order

**The `.mail` node order is not MED's**, although the two formats come from the same code. Code_Aster's MED reader (`lrmtyp.F90`) permutes every 3-D element between the two. meshio++ keeps the `.mail` order in the `"code_aster"` tables of the [node-ordering registry](../node_ordering.md):

- The linear cells, `TETRA10` and `PYRAM13` use meshio++'s (VTK's) order unchanged.
- `HEXA20`/`HEXA27` and `PENTA15`/`PENTA18` list the mid-edge nodes of the bottom ring, then the **vertical** mid-edges, then the top ring (VTK puts the verticals last), then `HEXA27`'s face centres (`1234`, `1265`, `2376`, `3487`, `1485`, `5678`) and body centre, or `PENTA18`'s quadrilateral face centres.

The tables come from Code_Aster's gmsh reader (`inigms.F90`) composed with meshio++'s gmsh tables. They were also run over the 565 `.mail` meshes of Code_Aster's own test suite (read locally, not redistributed). Both engines read every one identically, and more than 97% of the quadratic edges have their mid-node within 5% of the edge midpoint (the rest are genuinely curved or collapsed edges). The tables are also checked against Code_Aster's MED reader: both routes give the same element, up to a symmetry of the reference cell. The test fixtures are written in Code_Aster's numbering from the edge lists above. The tests check that every mid-edge node lands on its meshio++ edge midpoint and that every solid is positively oriented.

## Writing

- A `%` provenance block, then one `COOR_nD` block (`n` = the point dimension), one element block per cell block, a `GROUP_NO` per point region and a `GROUP_MA` per cell region, then `FIN`.
- Nodes are named `N1…` and elements `M1…`, numbered in block order. Names are limited to 8 characters, so more than 9,999,999 nodes or elements is a `WriteError` (write MED instead).
- Coordinates are written as `%.16E`, which round-trips a double exactly. Every line stays **within 80 columns**: a record that would cross column 80 continues on an indented next line.
- Group names are sanitised for Code_Aster. Anything but ASCII letters, digits and `_` becomes `_`, and the name is truncated to 24 characters. A name that then collides with an earlier one of the same kind gets `_1`, `_2`, …. Each change gets a warning.
- **Dropped, with a warning and a provenance note:** side regions (`.mail` has no facet group) and every data array (a `.mail` mesh has none). A cell type with no keyword (polygons, polyhedra, higher-order Lagrange cells) is a `WriteError`.

## Errors

These follow Code_Aster, where each is fatal. A node or element name defined twice, an element or group that names an undefined node or element, a block with no `FINSF`, a coordinate that is not a number, or a token where a keyword should be all raise `ReadError`, naming the culprit and its line. A group listing the same member twice gets a warning, and the member is kept once. A file with no `FIN` is read to the end with a warning.

## Notes

- The names in the file are not kept. A `.mail` → `.mail` round trip renames nodes and elements `N…`/`M…`, so a command file (`.comm`) that refers to individual nodes by name needs the groups instead, which are kept.
- Connectivity is read as written, never reoriented. Code_Aster accepts meshes whose volume cells are all mirrored (a few of its own test meshes are). Such a mesh reads fine, and [`attach_quality`](../mesh_quality.md) flags its cells as `quality:inverted`.
- Code_Aster drops an empty group. meshio++ keeps it as an empty region, as it does for every format.
- `TRIA7` needed a new cell type, `triangle7` (v16.0.0); see [cell types](../cell_types.md).
- **Not read:** `.comm` command files, and results. Code_Aster writes results as MED (`.rmed`), which the [MED](./med.md) reader handles.

The syntax reference is [U3.01.00, *The Code_Aster mesh file*](https://biba1632.gitlab.io/code-aster-manuals/docs/user/u3.01.00.html). The reader's and writer's rules above come from Code_Aster's source (`bibfor/modelisa/lrmast.F90`, `lirlig.F90`, `stkcoo.F90`, and `bibfor/stbtrias/` for its writer).
