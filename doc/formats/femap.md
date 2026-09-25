# Femap neutral file (`.neu`)

The neutral file is [Femap](https://www.sw.siemens.com/en-US/simcenter/femap/)'s documented ASCII interchange format, and the only open route into or out of a Femap model: a `.modfem` model is exported to it by Femap ([Femap `.modfem` route](../routes/femap_modfem.md), v16.13.0). Solvers write their results to it too: NX Nastran through Femap, MYSTRAN and EMSolution. It is a sequence of **data blocks**. Each opens with a line holding `-1`, then the block id, and closes at the next `-1` line. Records are comma-separated.

| | |
|---|---|
| **Format name** | `femap` |
| **Extensions** | `.neu` (also recognised by content: a lone `-1` line, then `100`) |
| **Read / Write** | ✓ (mesh, groups and results) / ✓ (mesh and groups) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.neu")                   # the first output set
last = meshioplusplus.read("model.neu", time_step=-1)     # the last one
times = meshioplusplus.femap.time_values("model.neu")     # one value per output set
meshioplusplus.write("out.neu", mesh)                     # Femap 8.2 layout, mesh and results
```

`read` also takes `points_only` and `arrays`, which narrow which output vectors are read. Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes.

## Blocks

Record layouts change with the Femap version in block `100`. The reader reads every record **by position and length** rather than by version, so one reader covers Femap 4.41 to 2401 (24.1):

| Block | Read as |
|---|---|
| `100` header | the version |
| `403` nodes | points: `x, y, z` are fields 11–13 in every version |
| `404` elements | cells. A record is seven lines, then, from version 4.5, one node list per non-zero list flag (rigid and weld elements), each ending at a `-1` line. The element's property and type become the `femap:property` and `femap:type` cell data |
| `402` properties | their titles name the `property_<id>` cell regions (tag = property id). From Femap 2401 a record repeats its value count and lists as many integers (function references) before the outline counts; the reader recognises that list by its shape |
| `408` groups | a group's node list becomes a `point` region and its element list a `cell` region of the group's title (tag = group id). Its rules are not evaluated: the lists are what the group holds |
| `450` output sets | **steps**: `time_step` selects one, and its value is `field_data["meshio:time"]` (with `femap:set` holding its id). The lines that follow a set's notes change with the version (none up to 9, one in 11.0, three from 11.2, seven in 2020.1), so the reader finds the next record by its shape |
| `451` / `1051` output vectors | the selected set's vectors become `point_data` (nodal, entity type 7) or `cell_data` (elemental, 8), named by their titles, NaN where a vector has no value. `451` holds `id,value` records; `1051` holds `start,end,value…` ranges, and in Femap 11 also `id,value` records |
| anything else (materials, loads, constraints, views…) | skipped |

## Element topologies

Femap stores every element's nodes in **20 slots laid out as a degenerate brick**: corners `0–3` on the bottom and `4–7` on the top, mid-edge nodes `8–11` on the bottom edges, `12–15` on the vertical edges and `16–19` on the top edges. A tetrahedron's apex is therefore slot 4, and its mid-edge nodes are the brick's.

| Topology | meshio++ type | Slots, in meshio++ node order |
|---|---|---|
| 0 / 1 | `line` / `line3` | 0 1 / 0 1 2 |
| 2 / 3 | `triangle` / `triangle6` | 0 1 2 / 0 1 2 4 5 6 |
| 4 / 5 | `quad` / `quad8` | 0–3 / 0–7 |
| 6 / 10 | `tetra` / `tetra10` | 0 1 2 4 / 0 1 2 4 8 9 10 12 13 14 |
| 7 / 11 | `wedge` / `wedge15` | 0 1 2 4 5 6 / 0 1 2 4 5 6 8 9 10 16 17 18 12 13 14 |
| 8 / 12 | `hexahedron` / `hexahedron20` | 0–7 / 0–11 16–19 12–15 |
| 9 | `vertex` | 0 |
| 14 | `pyramid` | 0–4 |

Rigid (13, 18), multi-list (15), contact (16) and weld (17) elements are skipped with one warning per topology, as is the 13-node pyramid (19), whose slots are not documented anywhere open.

## Writing

- The **Femap 8.2 layout**, the one the fixtures verify: blocks `100` (its title carries the provenance tag), `402`, `403`, `404`, `408`, and `450`/`451` for the results (v16.11.0).
- Each element's property is `femap:property` (else 1) and its type `femap:type` (else 1 rod, 17/18 linear/parabolic plate, 25/26 linear/parabolic solid, or 27 mass, by topology). A property's title comes from the cell region the reader makes of it: same tag, same cells. Other point and cell regions become groups, and a point region and a cell region with the same name share one group.
- Coordinates are written with 17 significant digits, which round-trips a double.
- **Results** (v16.11.0): one output set (`450`) whose id is `femap:set` (else 1) and whose value is `meshio:time` (else 0), holding one `451` vector per numeric point or cell data array, and one per component of a multi-component array (`<name>_0`, `<name>_1`, …): point arrays as nodal vectors (entity 7), cell arrays as elemental ones (entity 8). `femap:` arrays are not results; NaN values are left out of a vector, as Femap leaves out entities with no value. Reading the file back gives the same arrays, component by component; femap_neutral_parser reads every written value.
- **Dropped, with a warning and a provenance note:** non-numeric or ragged data arrays and other field data, side regions, and cells with no Femap topology (`quad9`, `hexahedron27`, polygons…).

## Errors

A file with no nodes (a geometry-only or results-only file has no mesh to read), a truncated or malformed record, a node or element id defined twice, an element that names an undefined node or leaves a slot its topology needs empty, and an output set that does not exist raise `ReadError`, naming the line. A group record that does not fit the layout ends that block's reading with a warning; the groups read before it are kept.

## Verification

No Femap installation was available. The layouts follow FrontISTR's `neu2fstr` converter and were checked against real files, all MIT-licensed:

- 36 Femap 8.2 models from FrontISTR, covering every solid topology, with node groups. The tests check that every mid-edge node lands on its meshio++ edge midpoint.
- EMSolution's 4.41 files, including pyramids, with 13 output sets written both as `451` and as `1051`, which the tests read to the same arrays.
- Result blocks written by Femap 2020.1, Femap 8.2 and MYSTRAN, from femap_neutral_parser. All 2,770 values of the 8.2 and MYSTRAN files equal femap_neutral_parser's own reading.

Real Femap 9.1–11.1 exports and Femap 2401 files (public, but without a licence, so not in the tests) were also read locally (v16.11.0): the 11.0 and 11.1 laminate models with their `1051` results, and the 2401 property layout, which led to the change above. Still unverified: element and group records from Femap 10 and from 12 on (the 2401 files hold geometry and properties but no elements), and the 13-node pyramid; no public file has one. A written file has not been imported into Femap.
