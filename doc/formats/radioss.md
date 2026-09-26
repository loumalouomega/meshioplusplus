# OpenRadioss / Radioss starter deck (`*_0000.rad`)

The input deck of the OpenRadioss explicit solver (and of Altair Radioss) in its native "block format": the **starter** file, `<run>_0000.rad`. Engine files (`_0001.rad`, run control) are refused. OpenRadioss also reads LS-DYNA `.k` decks natively, which meshio++ reads as [`lsdyna`](./lsdyna.md).

| | |
|---|---|
| **Format name** | `radioss` |
| **Extensions** | `.rad` (also recognised by content: `#RADIOSS STARTER`, or `/BEGIN` as the first non-comment line) |
| **Read / Write** | ✓ / ✓ (the mesh, parts, subsets, groups and surfaces; since v16.17.0) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("crash_0000.rad")
```

`read` takes no options. Both engines (the C++ core and the pure-Python reference) read the same meshes.

## Writing

```python
meshioplusplus.write("crash_0000.rad", mesh)                   # the mesh and its regions
meshioplusplus.radioss.write("crash_0000.rad", mesh, stubs=True)  # also placeholder materials and properties
```

Since v16.17.0 both engines write a starter deck, and write the same bytes: `#RADIOSS STARTER`, a `#` comment holding the provenance record, `/BEGIN` (the run name is the file name without `_0000`, the input version 2022, or the `radioss:version` the mesh was read with when it is 2017 or later, and the same input and work units, `kg m s`, so the starter converts nothing), `/PART` and `/SUBSET`, `/NODE`, one element card per run of cells of one block and one part, the groups and the surfaces, and `/END`. Integers fill 10 columns and reals 20: a real is written in its shortest exact spelling (the one Python's `repr` gives) when that fits, so coordinates come back exactly, and in 15 significant digits otherwise.

| meshio++ type | Card |
|---|---|
| `hexahedron` | `/BRICK` |
| `wedge`, `pyramid` | `/BRICK` with repeated nodes (`1 3 6 4 2 2 5 5`, `1 2 3 4 5 5 5 5`), which needs no formulation of its own, unlike `/PENTA6` |
| `tetra`, `tetra10` | `/TETRA4`, `/TETRA10` |
| `hexahedron20` | `/BRIC20`, through the `"radioss"` node-order table |
| `quad`, `triangle` | `/SHELL`, `/SH3N` |
| `line` | `/TRUSS` |

- **Parts.** With `radioss:part` cell data (a mesh read from a deck), its ids are the parts; a part's title is the name of the cell region of the same tag and cells, and its property and material are `radioss:property` and `radioss:material`. Without it, each cell region whose cells are all written and could share one property (solids, shells or trusses) becomes a part, its id the region's tag when positive and free; the other cells make one part per block, titled by the cell type.
- **Subsets, groups and surfaces.** Another cell region that is exactly a union of whole parts becomes a `/SUBSET` (the smallest first, so a subset inside another is its child); any other cell region a `/GRBRIC`, `/GRSHEL`, `/GRSH3N` or `/GRTRUS` per element family it holds; a point region `/GRNOD/NODE`; a side region `/SURF/SEG`, a solid's face by its corner nodes and a shell as itself. Ids are the region's tag when positive and free in its keyword, else the next free one. A title is one line of at most 100 characters; one starting with `#`, `$` or `/` gets a space in front, so it is not read as a comment or a keyword.
- **`stubs=True`** also writes a `/MAT/LAW1` (ρ = 1, E = 1, ν = 0.3) per material and a `/PROP/SOLID`, `/PROP/SHELL` (thickness 1) or `/PROP/TRUSS` (area 1) per property, so the starter accepts the deck on its own. A `/BRIC20`-only part gets `Isolid` 16; a part mixing `/BRIC20` and `/BRICK` the default, which the starter changes to 16 with a warning. A property two parts share but need different stubs for is split, a material 0 becomes 1, and a part mixing solids, shells or trusses is a `WriteError`. The values are placeholders: replace them with the model's own.
- **Dropped, with a warning and a provenance note:** cells with no card (`vertex`, other quadratic cells, polygons and polyhedra), side entries on lines, and every data array; `/SURF/PLANE` surfaces are written back from their `radioss:surf_plane:` field data. A mesh with no points is a `WriteError`.
- **What a round trip changes.** Blocks of one type written in several parts read back as one block; a region's order and the ids of regions without a tag follow the writer's; a segment names a face by its nodes, so a face two solids share, or two coincident shells, come back once, on the first cell. A `line` is a `/TRUSS`, so a zero-length one (a joint spring read from `/SPRING`) is kept but refused by the starter.

## The deck

- `#RADIOSS STARTER`, `/BEGIN` (run name, input version, units), then `/KEYWORD/option/id` blocks up to `/END`; anything after `/END` is ignored.
- Lines starting with `#` or `$` are comments; `#include file` inlines another file (relative to the including one, at most 8 deep, ended early by `#enddata`); a missing include is a warning.
- Fields are 10 columns wide for integers and 20 for reals from input version 51 on, 8 and 16 before; a line containing a comma is split on commas instead. The version is `field_data["radioss:version"]`.
- **Version 4.x decks** (v16.12.0): a deck without `/BEGIN` whose `#RADIOSS STARTER` line carries its version (`#RADIOSS STARTER      41…`, as Radioss 4.1 and 4.4 wrote them) is read in the 8- and 16-column format, and its titles are the keywords' last parts (`/PART/1/CUIVRE`, `/GRNOD/NODE/2/FIXED`), the data starting on the next line; a `/SKEW/FIX` then has no origin line.
- A keyword's id is its first integer field: `/BOX/RECTA/3/1` is box 3 given in unit system 1, `/SURF/PART/EXT/12` is surface 12.
- **Units** (v16.11.0): `/BEGIN`'s input and work units (mass, length and time, 20 columns each) set the length scale. Radioss runs in the work units, so coordinates and box dimensions are multiplied by input length ÷ work length, and the factor is `field_data["radioss:length_scale"]` when it is not 1. A length unit is an SI prefix and `m` (`mm`, `mum`, `km`…), `in`, `ft`, `yd`, `mi` or a number of metres; an unknown one leaves lengths as written, with a warning.
- **Local units** (v16.12.0): a `/UNIT/<id>` block (a title, then mass, length and time units) defines a unit system that `/NODE/<unit>`, `/BOX/…/<id>/<unit>`, `/SKEW/FIX/<id>/<unit>`, `/SURF/PLANE` and `/SURF/ELLIPS` can name; their lengths are then multiplied by that system's length ÷ the work length instead. A unit system that is not defined, or used without `/BEGIN` units, leaves the input units, with a warning.
- **Engine decks** (v16.12.0): reading `<run>_0001.rad` (`#RADIOSS ENGINE`) gives a mesh with no points whose field data holds every keyword and the numbers of its lines, as `radioss:engine:<keyword>` (`radioss:engine:RUN/<run>/1` the end time, `ANIM/DT` the start and interval of the animation files, `TFILE/4` the time history interval; an output request such as `ANIM/ELEM/SIGX` is an empty array). Reading a starter deck `<run>_0000.rad` adds the same field data from the `<run>_0001.rad` beside it.

| Card | meshio++ type |
|---|---|
| `/NODE` | points |
| `/BRICK/part` | `hexahedron`, or the tetrahedron, pyramid or wedge a brick with repeated nodes stands for (`1 2 3 3 4 4 4 4`, `1 2 3 4 5 5 5 5`, `1 2 3 3 5 6 7 7`, `1 2 3 1 5 6 7 5`, …) |
| `/PENTA6/part` | `wedge` |
| `/TETRA4/part` / `/TETRA10/part` | `tetra` / `tetra10` |
| `/BRIC20/part` | `hexahedron20`; with a zero mid-edge node, its corners as a `hexahedron` (with a warning) |
| `/SHELL/part` | `quad`, or `triangle` when its last two nodes coincide |
| `/SH3N/part`, `/QUAD/part`, `/TRIA/part` | `triangle`, `quad`, `triangle` |
| `/BEAM/part`, `/TRUSS/part`, `/SPRING/part` | `line` (the beam's orientation node is dropped; a single-node spring is skipped) |
| `/TSHELL`, `/SHEL16`, `/SPHCEL`, `/RIVET`, `/XELEM` | skipped, with one warning naming them |

The part id is the keyword's last field. Every cell's part, and the property and material its `/PART` names, are the `radioss:part`, `radioss:property` and `radioss:material` cell data.

A `/TETRA4` or `/TETRA10` written with the mirrored winding is reoriented, with a warning: OpenRadioss's own QA deck `INT_25` has every tetrahedron that way, while gmsh writes them positive.

## Node order

`/BRIC20` lists the bottom ring of mid-edge nodes, then the vertical ones, then the top ring, so it uses the `"radioss"` table of the [node-ordering registry](../node_ordering.md) (gmsh's `getVertexRAD`). `/TETRA10` and the other cards are in meshio++'s order.

## Regions

| Card | Region |
|---|---|
| `/PART/id` | a `cell` region named by its title (else `Part <id>`), tag = part id; parts with no elements are kept, empty |
| `/SUBSET/id` | a `cell` region of the parts that name this subset or one below it |
| `/GRNOD/<sub>/id` | a `point` region |
| `/GRBRIC`, `/GRSHEL`, `/GRSH3N`, `/GRQUAD`, `/GRTRIA`, `/GRBEAM`, `/GRTRUS`, `/GRSPRI` | a `cell` region |
| `/SURF/SEG/id` | a `side` region: each segment's facet (a shell segment is the shell's own face, as in LS-DYNA; a 2-node segment in a 2-D deck is a cell edge) |
| `/SURF/PART`, `/SUBSET`, `/MAT`, `/PROP` (v16.11.0) | the shells' own faces of those parts; with `/EXT` also the solids' faces no other of those solids shares, with `/ALL` every solid face |
| `/SURF/GRBRIC/EXT`, `/FREE` (v16.11.0) | the brick group's faces no other brick of the group shares (`EXT`), or no other solid of the model (`FREE`) |
| `/SURF/GRSHEL`, `/GRSH3N`, `/GRTRIA` (v16.11.0) | the group's own faces |
| `/SURF/SURF` (v16.11.0) | the union of other surfaces; a negative id removes one |
| `/SURF/BOX`, `/SURF/BOX2` (v16.12.0) | the shells with all (`BOX`) or any (`BOX2`) node in the box; with `/EXT` also the faces of the model's solids no other solid shares, with `/ALL` every solid face, chosen by their nodes the same way |
| `/SURF/PLANE`, `/SURF/ELLIPS` (v16.12.0) | analytical surfaces with no segments: `field_data["radioss:surf_plane:<id>"]` holds M and M1 (the plane through M, normal M → M1), `radioss:surf_ellips:<id>` the degree, the centre, the semi-axes and the rows of the skew's axes (the identity without a skew) |

A group's `<sub>` may be:

- the entity itself (`/GRBRIC/BRIC`, ids; a negative id removes one), `PART` (the parts' cells of that family, or their nodes for `GRNOD`), `SUBSET`, another group of the same keyword, and for `GRNOD` an element group or a `/SURF/SEG` surface;
- `GENE` (v16.11.0): pairs of first and last ids, every entity whose id falls in one; `GEN_INCR`: triples of first, last and step;
- `BOX` and `BOX2` (v16.11.0): box ids (a negative one removed). A node is taken when it lies inside the box; an element, under `BOX` when all its nodes do, under `BOX2` when any does, as the starter decides.

Boxes are `/BOX/RECTA` (the corners, from two nodes or coordinates), `/BOX/CYLIN` (the axis ends and a diameter), `/BOX/SPHER` (the centre and a diameter), and `/BOX/BOX` (other boxes; a negative id subtracts one). A `/BOX/RECTA` with a skew (v16.12.0) has its edges along the skew's axes: a `/SKEW/FIX` (origin, X and Y axes; Z = X × Y, Y recomputed as Z × X) turns the corners and the tested point into its frame; a skew that is not a `/SKEW/FIX` leaves the box empty, with a warning. Other generators are not resolved: the region is written empty with a warning. Every region is named by its title line (else `<KEYWORD>_<id>`), tag = its id; a name already taken gets ` [KEYWORD]` appended.

Materials, properties, loads, contacts, rigid bodies and every other keyword are skipped.

## Results

The animation files Radioss writes (`<run>A001`, …) are read by [`radioss_anim`](./radioss_anim.md) (v16.11.0). The time-history files (`<run>T01`, …) are read by [`radioss_th`](./radioss_th.md) (v16.12.0).

## Validation

Surface forms other than those above (`/SURF/DSURF`, Madymo surfaces…) are skipped with a warning.

**The writer** (v16.17.0) was checked with the OpenRadioss starter (build `latest-20260728`): every fixture and every one of the 81 QA decks that has a mesh (76; five keep theirs in include files the suite does not ship) was read, written with `stubs=True` and run through the starter, which accepts 75 with no error; the other is the zero-length truss above. Read back, 70 of the 76 give the same points, cells, parts, groups and surfaces, and the other six differ only in the surface entries on faces two solids share or on coincident shells. The C++ and Python writers give the same bytes on all of them, with and without stubs. `test_radioss.py` runs the starter on the fixtures when `MESHIOPLUSPLUS_OPENRADIOSS` names an OpenRadioss install.

The fixtures under `tests/python/meshes/radioss/` are written by `tools/gen_radioss_fixtures.py` from the card layouts of OpenRadioss's `hm_cfg_files`, plus two decks gmsh wrote, which read to the same cells as gmsh's own `.msh` of the same mesh; OpenRadioss's own QA decks are CC BY-NC and are not copied. Outside the repository the reader was run on the 81 starter decks of OpenRadioss's `qa-tests/`: every one reads, both engines agree, and every cell has a positive volume. For v16.11.0 the box, generator and surface semantics were taken from the OpenRadioss starter source (read, not copied) and the QA decks that use them (`NEW_BOX_V11`) resolve with both engines alike. v16.12.0: the only Radioss 4.x deck in OpenRadioss's QA suite (`ALE2D/bar2_v41/data/BAR2V41BD02`, written by Radioss 4.1 in 2000; not copied) reads to the same points, cells and regions, in both engines, as the 2017 translation of it that the suite also keeps. `/UNIT`, skewed boxes, box and analytical surfaces and engine decks follow the Radioss 2022 Reference Guide and are checked on synthetic decks.
