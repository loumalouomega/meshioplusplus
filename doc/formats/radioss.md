# OpenRadioss / Radioss starter deck (`*_0000.rad`)

The input deck of the OpenRadioss explicit solver (and of Altair Radioss) in its native "block format": the **starter** file, `<run>_0000.rad`. Engine files (`_0001.rad`, run control) are refused. OpenRadioss also reads LS-DYNA `.k` decks natively, which meshio++ reads as [`lsdyna`](./lsdyna.md).

| | |
|---|---|
| **Format name** | `radioss` |
| **Extensions** | `.rad` (also recognised by content: `#RADIOSS STARTER`, or `/BEGIN` as the first non-comment line) |
| **Read / Write** | ✓ / — |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("crash_0000.rad")
```

`read` takes no options. Both engines (the C++ core and the pure-Python reference) read the same meshes.

## The deck

- `#RADIOSS STARTER`, `/BEGIN` (run name, input version, units), then `/KEYWORD/option/id` blocks up to `/END`; anything after `/END` is ignored.
- Lines starting with `#` or `$` are comments; `#include file` inlines another file (relative to the including one, at most 8 deep, ended early by `#enddata`); a missing include is a warning.
- Fields are 10 columns wide for integers and 20 for reals from input version 51 on, 8 and 16 before; a line containing a comma is split on commas instead. The version is `field_data["radioss:version"]`. **Units are not applied**: coordinates are read as written.

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

A group's `<sub>` may be the entity itself (`/GRBRIC/BRIC`, ids; a negative id removes one), `PART` (the parts' cells of that family, or their nodes for `GRNOD`), `SUBSET`, another group of the same keyword, and for `GRNOD` an element group or a `SURF`. Other generators (`BOX`, `GENE`, …) are not resolved: the region is written empty with a warning. Every region is named by its title line (else `<KEYWORD>_<id>`), tag = its id; a name already taken gets ` [KEYWORD]` appended.

Materials, properties, loads, contacts, rigid bodies and every other keyword are skipped.

## Results

Radioss writes animation files (`A001`, …) and time history (`T01`), which meshio++ does not read. OpenRadioss converts them itself: `anim_to_vtk` writes legacy VTK, which meshio++ reads, and `th_to_csv` writes the time history as CSV (both in the [OpenRadioss tools](https://github.com/OpenRadioss/OpenRadioss/tree/main/tools)).

## Validation

The fixtures under `tests/python/meshes/radioss/` are written by `tools/gen_radioss_fixtures.py` from the card layouts of OpenRadioss's `hm_cfg_files`, plus two decks gmsh wrote, which read to the same cells as gmsh's own `.msh` of the same mesh; OpenRadioss's own QA decks are CC BY-NC and are not copied. Outside the repository the reader was run on the 81 starter decks of OpenRadioss's `qa-tests/`: every one reads, both engines agree, and every cell has a positive volume.
