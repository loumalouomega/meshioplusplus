# MSC Marc input deck (`.dat`) and post file (`.t19`)

MSC Marc's input deck (the `.dat` file Marc Mentat writes and the solver reads) and its formatted post file (`.t19`, the ASCII twin of the binary `.t16` results file). meshio++ reads both: the deck's mesh and sets, and the post file's mesh, sets and results, one increment at a time.

| | |
|---|---|
| **Format names** | `marc` (input deck), `marc_t19` (formatted post file) |
| **Extensions** | `.dat` (shared with Tecplot: a `.dat` file that opens as a Marc deck is Marc's), `.t19`; both also recognised by content |
| **Read / Write** | ✓ / ✓ for the deck (since v16.17.0, by name); `marc_t19` ✓ / —, [read-only by design](../conformance.md#marc-t19): Marc writes it, and no tool reads one written elsewhere |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("job.dat")                 # a Marc deck (a Tecplot .dat still reads as Tecplot)
post = meshioplusplus.read("job.t19", time_step=-1)   # the last increment
meshioplusplus.marc.time_values("job.t19")            # every increment's time
for time, mesh in meshioplusplus.read_sequence("job.t19"):
    ...
```

```bash
meshioplusplus convert job.dat job.vtu
meshioplusplus convert job.t19 'inc_{step}.vtu'      # one .vtu per increment
```

`read_t19` (and `read` on a `.t19`) also takes `points_only` and `arrays` (read only the named arrays). The binary `.t16` is not read: Marc writes the `.t19` beside it when the `POST` option asks for a formatted file (its fourth field set to 1), and `.t16` files need Marc's own PyPost library: the [Marc `.t16` route](../routes/marc_t16.md) (v16.13.0) exports one through it.

## Writing

```python
meshioplusplus.write("job.dat", mesh, file_format="marc")   # ".dat" alone writes Tecplot
meshioplusplus.marc.write("job.dat", mesh)
```

Since v16.17.0 both engines write an input deck, and write the same bytes, in `EXTENDED` format (10-column integers, 20-column reals, the shortest exact spelling of a real when it fits, so coordinates come back exactly): `TITLE`, a `$` comment holding the provenance record, `EXTENDED`, `SIZING` (the element and node counts), one `ELEMENTS` line per type, `END`, then `CONNECTIVITY` (at most 14 nodes on an element's first line, continued 14 to a line), `COORDINATES` (three per node), a `DEFINE NODE SET` per point region and a `DEFINE ELEMENT SET` per cell region (ascending runs of three or more as `a TO b`, six items to a line, `C` continuing), the face and edge sets read from a deck back from their `marc:face_set:`/`marc:edge_set:` field data, and `END OPTION`. `.dat` stays Tecplot's for a write by extension, whatever the file it replaces holds (see [Which `.dat` is Marc's](#which-dat-is-marc-s)): name the format.

- **Element types.** A cell keeps its `marc:type` when that type's element is the cell (or a brick of type 7, 43, 117… for a tetrahedron, pyramid or wedge read from a degenerate brick, written back that way); a Herrmann, bubble or generalized-plane-strain type, whose extra nodes a cell does not keep, falls back to the default with a warning. The defaults: in a 3-D mesh `hexahedron` 7, `hexahedron20` 21, `tetra` 134, `tetra10` 127, `wedge` 136, `pyramid` 7 as a degenerate brick, `quad` 75 and `triangle` 138 (shells), `quad8` 22, `line` 9 (truss) and `line3` 64; in a planar mesh (two coordinates, or every z 0, and no volume cells) the plane-strain `quad` 11, `triangle` 6, `quad8` 27 and `triangle6` 125. A 3-D `triangle6` has no default (Marc has no 6-node shell) and is dropped, as are cells with no type (`vertex`, `hexahedron27`…).
- **Numbers.** Elements keep `marc:element` when every number is positive and unique, else are numbered 1, 2…; nodes are numbered 1, 2… in point order.
- **Set names** read back as one token: blanks, commas and `$` become `_`; `define` and `end`, which the line after a set could start with, become `define_set` and `end_set`; a name already taken, case aside, gets `_2`, `_3`….
- **Dropped, with a warning and a provenance note:** side regions (Marc numbers an element's faces and edges its own way, which is not mapped to facets), every data array but `marc:element` and `marc:type`, and cells with no type. A mesh with no points is a `WriteError`.
- **Unchecked against Marc:** no licence was available, so `SIZING`'s fields, the one-type-per-line `ELEMENTS` and the defaults for shells and plane strain follow Volume C and Mentat's own decks, not a run of Marc; see [Awaiting a licensed run](../roadmap.md#awaiting-a-licensed-run).

## Which `.dat` is Marc's

`.dat` is Tecplot's extension too. A `.dat` file is read as a Marc deck when its first line that is not a `$` comment starts with a Marc parameter (`title`, `sizing`, `extended`, `elements`, `version`, `table`, `processor`, `alloc`, `setname`, `end` ...) with no `=` in it, and an `END`, `CONNECTIVITY` or `COORDINATES` line follows in the first 64 KB. Anything else stays Tecplot's: a Tecplot header's `TITLE = "..."` has its `=`. Both engines apply the same check: the C++ registry in its format resolution (as it does for MFEM's `.mesh`), the Python reader by refusing a file that fails it so that the Tecplot reader is tried next. `meshioplusplus.sniff_format` recognises both a deck and a post file whatever their names.

## The input deck

A deck has three sections: parameters, up to `END`; the model definition, up to `END OPTION`; and the history definition (load cases), which is not read. From the parameters only `EXTENDED` matters: it doubles every field's width. The model definition options read are:

- **`COORDINATES`**: a header (coordinates per node, node count), then a node number and its coordinates per line. Nodes need not be numbered in order; a node defined twice keeps its last definition. Decks with two coordinates per node (plane and axisymmetric models) get `z = 0`.
- **`CONNECTIVITY`**: a header, then per element its number, its Marc element type and its nodes; an element with more than 14 nodes continues on the next lines.
- **`DEFINE ELEMENT SET`** and **`DEFINE NODE SET`** (and their ordered `ELSQ`/`NDSQ` forms), **`DEFINE EDGE SET`** and **`DEFINE FACE SET`**: see [Sets](#sets).
- **`INCLUDE`** (v16.12.0): the line `INCLUDE file` (a blank or a comma before the name) is replaced by the lines of the file it names, relative to the including file, anywhere in the deck and recursively (a chain deeper than 16 is refused as a cycle). A missing file is an error.

Every other option is skipped, its data lines with it. Fields are read in any of Marc's three layouts, which a deck may mix line by line:

- **fixed**: 5-column integers and 10-column reals;
- **extended** (after `EXTENDED`): 10-column integers and 20-column reals. Marc Mentat writes reals with 16 digits and an exponent without a letter, `6.666666666666666-1`, and a negative value fills its field, so it touches the one before it (`7-1.666666666666666-1`): fields are cut by column, never split on blanks;
- **free**: comma-separated items; a lone item is followed by a comma (`1,`).

Reals may use `E` or `D` exponents or none. `$` starts a comment line.

## Element types

Every Marc element type meshio++ maps lists its nodes in meshio++'s order (corners first, the face 1-2-3(-4) numbered so that its normal points into the element, then the mid-edge nodes of the bottom face, the top face and the vertical edges), so no [node permutation](../node_ordering.md) applies. Checked against the Volume B descriptions and on real Marc Mentat decks, where every cell is positively oriented and every mid-edge node sits at its edge's midpoint.

| meshio++ cell | Marc element types |
|---|---|
| `line` | 9, 31, 52, 98; rebar 165, 166, 167 |
| `line3` | 64; rebar 168, 169, 170 (which list their middle node second) |
| `triangle` | 2, 6, 138, 158, 201; Herrmann 155, 156 (their fourth node, a centroid bubble, is dropped) |
| `triangle6` | 124, 125, 126, 200; Herrmann 128, 129 |
| `quad` | 3, 10, 11, 18, 75, 139, 140; rebar 143, 144, 145, 147; composite 151, 152; Herrmann 80, 82, 83, 118, 119 (their fifth node, the pressure node, is dropped) and 81 (also its two generalized plane strain nodes) |
| `quad8` | 22, 26, 27, 28, 30, 53, 54, 55; Herrmann 32, 33, 58, 59, 63, 66, and 34, 60 (their two generalized plane strain nodes dropped); rebar 46, 48, 142, 148, and 47 (as 34); composite 153, 154 |
| `tetra` | 134, 135; Herrmann 157 (its fifth node, a bubble, is dropped) |
| `tetra10` | 127, 133; Herrmann 130 |
| `wedge` | 136, 137 |
| `hexahedron` | 7, 43, 117, 123; rebar 146; composite 149; Herrmann 84, 120 (their ninth node, the pressure node, is dropped); with repeated nodes a wedge, pyramid or tetrahedron, as LS-DYNA's and Radioss's bricks |
| `hexahedron20` | 21, 44, 57; Herrmann 35, 61; rebar 23; composite 150 (with repeated corners, kept as is with a warning) |

The Herrmann (mixed displacement–pressure), rebar and composite types (v16.12.0) are read from Volume B (2000 to 2005 editions): the ones whose pressure unknown sits at the corners, and the rebar and composite elements, list the nodes of their plain twins; the others add a pressure node (without coordinates, so it need not be defined), a centroid bubble node or generalized plane strain nodes after the geometric nodes, which are dropped. Marc's element library up to 2005 has no interface, gasket or cohesive element and no 15-node pentahedron (its "15-node" solid is a collapsed 20-node brick).

Elements of other types are skipped with one warning listing them; their node count must be known to read a deck whose elements continue on further lines, so an element of an unknown type with more than 14 nodes is refused. `cell_data["marc:element"]` holds each cell's element number and `cell_data["marc:type"]` its Marc type.

## Sets

`DEFINE ELEMENT SET name` and `DEFINE NODE SET name` become `cell` and `point` [regions](../regions.md) named as the deck names them. A set's list is read left to right:

- element or node numbers, and ranges `a TO b` (or `THROUGH`) with an optional `BY step`, counting up or down;
- the names of sets defined before it, of the same kind;
- `AND` (the default), `EXCEPT` and `INTERSECT` combining what follows with what came before;
- a line ending in `C` (or `CONTINUE`) continues on the next.

A member element with no cell (a skipped type) drops out of the region; a member number the deck never defines is an error.

`DEFINE EDGE SET name` and `DEFINE FACE SET name` (v16.12.0) list `element:number` members. Marc numbers an element's edges and faces its own way (Volume A, which was not available), so they are not mapped to meshio++ facets: each set is `field_data["marc:edge_set:<name>"]` or `["marc:face_set:<name>"]`, one row per member, `(cell, Marc edge or face number)`, the cell −1 for an element with no cell. Integration-point and other set kinds are skipped with a warning.

## The post file

The formatted post file (post file revision 9 or later, Volume D's PLDUMP2000) is a sequence of blocks, `=beg=5xxnn (name)` ... `=end=`, in 13-column fields (`i13`, `e13.6`; a negative real touches the field before it). The model blocks come first, then one group of blocks per increment between `****` and `----`.

| `.t19` | meshio++ |
|---|---|
| element connectivities (507), nodal coordinates (508) | points and cells, with the element types above |
| sets (513): element and node sets | `cell` and `point` regions |
| sets (513): edge (12) and face (13, and the ordered 18, 19) sets | `field_data["marc:edge_set:<name>"]`, `["marc:face_set:<name>"]` as for the deck |
| an increment's time; its frequency (modal, harmonic) or buckling factor | `field_data["meshio:time"]` |
| the increment and sub-increment numbers (517) | `field_data["marc:increment"]`, `"marc:subincrement"` |
| each nodal vector (524): `Displacement`, `Reaction Force`, `Temperature` ... | point data under the file's name, `(points, components)` or `(points,)` |
| the imaginary part of a complex harmonic vector | point data `<name>@imag` |
| element post codes at the integration points (523) | cell data per code, flattened point-major as the [Abaqus `.fil`](./abaqus_fil.md) reader does so that every writer holds it: `(cells, points × components)` with `field_data["marc:layout:<name>"] = [points, components]`; `(cells, components)` or `(cells,)` when elements have one integration point, with no layout |

**Element post codes.** Each code is named by its label in the file, or else by its meaning in Volume C's Table 3-3 (`Equivalent Von Mises Stress` for 17, `Temperature` for 9 ...), or `post code <n>`. A tensor (codes 301 total strain, 311 stress, 321 plastic strain, 341 Cauchy stress, 401 elastic strain, 411 global stress ... written as six consecutive codes) becomes one six-component array, components 11 22 33 12 23 31, that is `xx yy zz xy yz zx`; the order is checked on a real Marc post file, whose von Mises stress (code 17) matches its stress tensor's. A code of layer `n` (code + 1000 n) gets the suffix `@layer<n>`. A cell whose element has no cell (a skipped type) drops its values.

An increment that remeshes the model (v16.12.0; its remeshing flag set in block 517, the model blocks 502 to 514 repeated after its block 519) brings its own mesh: its step, and every later one until the next remeshing, is built on that model, with its own sets; element post codes carry over when the new model does not repeat them. Post files of revision 8 or earlier (Marc K7 and before, Volume D's PLDUMP layout, without block markers) are refused: no such file was found to check a reader against. Contact bodies, springs, distributed loads, tyings and global variables are skipped.

## Validation

No Marc licence was available. The readers are checked against:

- **Real Marc Mentat 2020 decks**: the twelve decks of [DAMASK](https://github.com/damask-multiphysics/DAMASK)'s Marc element library tests (element types 6, 7, 11, 21, 27, 54, 57, 117, 125, 127, 134 and 136, extended format), not committed (they are AGPL). Every cell is positively oriented, every mid-edge node sits at its edge's midpoint, and the element and node sets resolve (the grain sets partition the elements).
- **A real `.t19`** written by Marc (two 8-node bricks, from the FEDES project's examples), not committed: its mesh, sets, loads and stress tensor read as the file describes them.
- **Fixtures written from the manuals**, in `tests/python/meshes/marc/` by `tools/gen_marc_fixtures.py`: the same model in fixed and extended format, a free-format deck with a collapsed brick and an unknown type, a plane-strain deck, a deck that INCLUDEs its coordinates and sets (nested, with Herrmann, rebar and composite elements and edge and face sets), a two-increment post file, and a post file whose second increment remeshes.

The C++ and Python readers agree bit for bit on all of them. The writer (v16.17.0) reads back every fixture and every DAMASK deck to the same points, cells, element numbers and sets, and the same element types except the Herrmann elements of `include_main.dat`, whose pressure nodes a cell does not keep; the C++ and Python writers give the same bytes on all of them.

## Not read

- The history definition (load cases, boundary conditions, loads), and model definition options other than the three above (materials, geometry properties, boundary conditions, tables).
- The mapping of edge and face sets to facets (their numbering is Volume A's), which are kept as field data instead.
- The binary post file `.t16` (see the [route](../routes/marc_t16.md)), and post files of revision 8 or earlier.
- Marc element types not in the table above, and element types added after 2005.
