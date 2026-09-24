# LS-DYNA keyword input (`.k`, `.key`, `.dyn`)

The [LS-DYNA keyword input format](https://lsdyna.ansys.com/manuals-download/): a keyword-driven ASCII deck of `*NODE`, `*ELEMENT_*`, `*PART` and `*SET_*` cards, split over as many `*INCLUDE`d files as the model needs. OpenRadioss reads it natively too, so it is the input format of two explicit solvers.

| | |
|---|---|
| **Format name** | `lsdyna` |
| **Extensions** | `.k`, `.key`, `.dyn` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

meshio++ reads the mesh (nodes, elements, parts and sets) and nothing else. Materials, sections, contacts, loads and every other keyword are skipped without a message: they are not geometry. The results of a run, in the binary `d3plot` family, are read by [`lsdyna_d3plot`](./lsdyna_d3plot.md).

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("crash.k")           # follows *INCLUDE
meshioplusplus.write("out.k", mesh)             # nodes, elements, *PART, *SET_*

for r in mesh.regions:
    print(r.kind, r.name, r.tag, r.dim, len(r.entries))
```

Both engines (the C++ core and the Python reference) read the whole format and write byte-identical files. Buffers are read by the Python reader.

## Card formats

A card is one line whose fields sit in fixed columns. Four layouts exist, they can be **mixed in one file**, and the reader tells them apart per card:

| Layout | How it is selected | Fields |
|---|---|---|
| Standard | the default | the widths each card declares: `*NODE` is `I8, 3E16, 2F8`, an element card is ten `I8`, a set or `*PART` card is `I10` |
| Long | `*KEYWORD LONG=Y` (or `S`, `K`) for the file, or a `+` suffix on one keyword (`*NODE+`) | **every** field, integer and real, is 20 columns |
| I10 | `*KEYWORD I10=Y` for the file, or a `%` suffix on one keyword | only the 8-column integer fields widen to 10; a different mechanism from long |
| Free | a comma anywhere on the line | comma-separated, any width |

A `-` suffix (`*NODE-`) puts one keyword back to standard inside a long file. A `*KEYWORD` card sets the default for the file it appears in and for the files that file includes. Fortran real spellings are accepted (`1.5D+01`, `2.5-1` for `2.5e-1`).

## Elements

| Card | meshio++ cell | Notes |
|---|---|---|
| `*ELEMENT_SOLID` | `hexahedron`, `wedge`, `pyramid`, `tetra`, `tetra10` | see the collapse table; one card (`eid pid n1..n8`) or two lines (`eid pid`, then up to ten nodes) |
| `*ELEMENT_SHELL` | `triangle`, `quad` | `n4 == n3` or blank is a triangle; `n5..n8` (mid-side nodes) are ignored with a warning |
| `*ELEMENT_TSHELL` | as `*ELEMENT_SOLID` | thick shells are hexahedra |
| `*ELEMENT_BEAM` | `line` | the orientation node `n3` is dropped; an optional second card is skipped |
| `*ELEMENT_DISCRETE` | `line` | |
| `*ELEMENT_MASS` | `vertex` | the mass is dropped |

`*ELEMENT_SHELL` accepts the `_THICKNESS`, `_BETA`, `_MCID` and `_OFFSET` variants (one extra card per option, skipped); `*ELEMENT_SOLID_ORTHO` skips its two extra cards. Any other variant, and the `*ELEMENT_SOLID_TET4TOTET10` / `_H8TOH20` conversion directives, is skipped with a warning.

Element ids are unique per family (solids, shells, beams and so on may reuse a number, as LS-DYNA allows) and a duplicate within a family is an error. A card that references an undefined node is an error.

### Degenerate elements

LS-DYNA has no tetra, pyramid or wedge card: they are hexahedra with repeated nodes, collapsed on read and expanded on write.

| Element | LS-DYNA nodes | Read | Written |
|---|---|---|---|
| tetra | `n1 n2 n3 n4 n4 n4 n4 n4` | `tetra (n1 n2 n3 n4)` | this form |
| tetra | `n1 n2 n3 n3 n5 n5 n5 n5` | `tetra (n1 n2 n3 n5)` | |
| pyramid | `n1 n2 n3 n4 n5 n5 n5 n5` | `pyramid` | this form |
| wedge | `n1 n2 n3 n4 n5 n5 n6 n6` | `wedge (n1 n5 n2 n4 n6 n3)` | this form |
| wedge | `n1 n2 n3 n3 n5 n6 n7 n7` | `wedge (n1 n2 n3 n5 n6 n7)` | |
| wedge | any side edge collapsed in both faces, e.g. `n1 n2 n3 n1 n5 n6 n7 n5` | `wedge (n1 n2 n3 n5 n6 n7)`, the faces' cyclic order kept | |
| shell triangle | `n1 n2 n3 n3` | `triangle` | this form |

The wedge spellings are read because they occur in real decks (the second is what most preprocessors write; the edge-collapsed form, added in v16.7.0 and shared with the [Radioss](./radioss.md) reader, is what Radioss decks use). The node order is mapped so that the reading has a positive volume in the meshio++ convention; the tests check this numerically. A degenerate pattern not in the table stays a `hexahedron`.

## Parts and sets

A `*PART` becomes a **cell region**: the title is the name (`Part <pid>` when the title is blank), the `pid` is the `tag`, and the `dim` is the highest topological dimension of its elements. An empty part is kept. The `*PART_CONTACT`, `_COMPOSITE` and `_INERTIA` variants read the same first two cards.

| Keyword | Region |
|---|---|
| `*SET_NODE` | `point` |
| `*SET_SOLID`, `_SHELL`, `_TSHELL`, `_BEAM`, `_DISCRETE` | `cell` |
| `*SET_PART` | `cell`: the union of the listed parts' elements |
| `*SET_SEGMENT` | `side`: each segment's corner nodes are matched to a (cell, facet) pair |

A set region has `dim = -1` and the set id as `tag`. Its name is the title of a `_TITLE` set, otherwise `<FAMILY> set <sid>` (`NODE set 20`). The `_LIST`, `_GENERATE`, `_LIST_GENERATE` and `_TITLE` suffixes are all read. Two sets that would get the same name, tag and dimension (a titled solid set and shell set with the same title and id) are told apart by a ` [<family>]` suffix.

A segment with `n4 == n3` (or a blank `n4`) is a triangle. A segment shared by two cells goes to the lower-numbered one, a shell element's own face is facet 0, and a segment that matches no face is dropped with one warning. Entries of a set that name an undefined node or element are dropped with one warning as well.

## Blocks

A cell block is one cell type of one element keyword, so blocks follow the sections of the deck: a `*ELEMENT_SOLID` that mixes hexahedra and wedges yields two blocks, and the next keyword extends the previous block when its first type matches. Cells are numbered block by block, which is the order a `Cell` region's indices use.

## `*INCLUDE`

`*INCLUDE`, `*INCLUDE_NO_TRANSFORM` and `*INCLUDE_TRANSFORM` name one file per line (a line ending in `+` continues on the next). A name is looked up relative to the including file, then in each `*INCLUDE_PATH` / `*INCLUDE_PATH_RELATIVE` directory, then in the working directory; a backslash is accepted as a separator. Everything is collected into one state before ids are resolved, so a set may refer to elements defined in a later file. Nesting is limited to eight levels.

- A file that cannot be found is a **warning**, not an error: decks routinely include material and contact files that are not shipped with the mesh. An element that then references a missing node is an error.
- `*INCLUDE_TRANSFORM` offsets are not applied; the file is included untransformed, with a warning.

## Skipped content

| Content | Behaviour |
|---|---|
| Any other keyword | skipped silently |
| `$` comment lines | skipped |
| `*END` | stops reading the file it appears in |
| `*PARAMETER` substitution (`&name` in a mesh card) | the card is skipped; one warning reports how many |
| PGP-encrypted blocks | skipped with a warning |
| 20-node hexahedra (a third element line) | not supported |

## Writing

The writer emits `*KEYWORD`, the provenance comment, `*NODE`, one `*ELEMENT_*` keyword per cell block, one `*PART` per part, one `*SET_*_LIST_TITLE` per remaining region, and `*END`. Coordinates are written in the 16-column real field with the shortest scientific spelling that round-trips, and as many digits as fit otherwise (11 significant digits, 10 for a negative number). A mesh whose ids do not fit an 8-column field is refused rather than written in a different format.

- **Parts.** A cell region with a dimension (`dim >= 0`) is a part when none of its cells is already in an earlier part; its tag is the `pid` when that is positive and unused, otherwise the next free number. Every other cell region is written as a set, so no region is lost. Cells that no part claims form one part per cell block, named after the cell type. Every part gets the placeholder section and material ids `1`, `1`.
- **Sets.** A cell region that spans several element families is written as one set per family, with a warning. Point regions become `*SET_NODE_LIST_TITLE`, side regions `*SET_SEGMENT_TITLE` (the facet's corner nodes; a triangle repeats its last node).
- **Cell types.** `vertex`, `line`, `triangle`, `quad`, `tetra`, `pyramid`, `wedge`, `hexahedron` and `tetra10` are written; anything else is a `WriteError` that names the type.
- A region name that starts with `*` or `$` is written with a leading space so it is not read back as a keyword or a comment.

## Quirks & limitations

- **Sections and materials are placeholders.** The output is meant for a mesher to LS-DYNA hand-off: add real `*SECTION_*` and `*MAT_*` cards before running it. It has not been loaded in LS-PrePost.
- **Part sets and parts are not distinguished on write.** A `*SET_PART` is read as the elements of its parts and written back as an element set.
- **Empty parts come back as empty sets.** A part with no elements has `dim = -1`, which the writer treats as a set.
- **Shell nodes past `n4`, beam orientation nodes, element masses, `*NODE` constraints (`tc`, `rc`) and all part properties are dropped.**
- **Duplicate nodes** are an error: LS-DYNA itself refuses them.
- Nodes are always read as three coordinates; a 2-D mesh is written with `z = 0`.
