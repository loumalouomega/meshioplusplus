# MSC Patran 2 neutral file (`.pat`, `.out`)

The Patran 2 neutral file is MSC Patran's legacy ASCII interchange format. Patran, Cubit and ANSA still export it, and Fluent and Feko import it. It is a sequence of **packets**. Each opens with a fixed-width header card `(I2,8I8)`, the fields `IT, ID, IV, KC, N1…N5`, where `IT` is the packet type and `KC` the number of data cards that follow.

| | |
|---|---|
| **Format name** | `patran` |
| **Extensions** | `.pat`, `.out` (also recognised by content: a first card `25` or `26` in the fixed header columns) |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.pat")      # points, cells, components as regions
meshioplusplus.write("out.pat", mesh)        # or meshioplusplus.patran.write(...)
```

Neither `read` nor `write` takes extra options. Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes.

## Packets

| Packet | Read as |
|---|---|
| `25` title | skipped (the writer puts the provenance tag here) |
| `26` summary | skipped |
| `01` node | a point: `ID` is the node id, the first data card `(3E16.9)` its coordinates. Patran always writes them in the global frame, so the frame fields of the second card are ignored |
| `02` element | a cell: `ID` is the element id and `IV` its shape; the node count on the first data card tells linear from quadratic. The property id there becomes the `patran:property` cell data |
| `21` named component | a name card, then `(type, id)` pairs. Type `5` (node) entries become a `point` region, the element types (`6` bar, `7` tri, `8` quad, `9` tet, `10` pyramid, `11` wedge, `12` hex) a `cell` region of the same name. Both are tagged with the component number |
| `99` | end of file |
| `06` distributed load | a row of the `patran:distributed_load` and `patran:distributed_load_values` field data (v16.12.0; see [Loads and boundary conditions](#loads-and-boundary-conditions)) |
| `07` node force, `08` node displacement | `patran:force:<set>`, `patran:displacement:<set>` point data (v16.12.0) |
| `10` node temperature, `11` element temperature | `patran:temperature:<set>` point data, `patran:element_temperature:<set>` cell data (v16.12.0) |
| anything else (materials, properties, coordinate frames, MPCs, …) | skipped by its card count |

| Shape (`IV`) | Nodes | meshio++ type |
|---|---|---|
| 2 bar | 2 / 3 | `line` / `line3` |
| 3 tri | 3 / 6 / 7 | `triangle` / `triangle6` / `triangle7` |
| 4 quad | 4 / 8 / 9 | `quad` / `quad8` / `quad9` |
| 5 tet | 4 / 10 | `tetra` / `tetra10` |
| 6 pyramid | 5 / 13 | `pyramid` / `pyramid13` |
| 7 wedge | 6 / 15 | `wedge` / `wedge15` |
| 8 hex | 8 / 20 | `hexahedron` / `hexahedron20` |

Any other shape and node count (a 27-node hex, a point element) is skipped with one warning per combination. Node and element ids may have gaps and are renumbered. Elements that no component names are grouped by property into `property_<pid>` cell regions (tag = pid), so a file with no components still has its parts.

## Loads and boundary conditions

Since v16.12.0 the load and boundary-condition packets are read, per load or constraint set (the header's `IV`, `<set>` below), and written back from the same arrays:

- **`07` node forces and `08` node displacements** (constraints): `patran:force:<set>` and `patran:displacement:<set>`, point data of six components (three translations, three rotations), `NaN` where the packet gives no value (its component flag is off, or the node has no packet). A node whose packet names a coordinate frame other than 0 has it in `patran:force_frame:<set>` or `patran:displacement_frame:<set>` (an integer array, only present when some frame is nonzero); the values stay in that frame.
- **`10` node and `11` element temperatures:** `patran:temperature:<set>` point data and `patran:element_temperature:<set>` cell data; `NaN` where there is no packet or its data flag `N1` says the value is a dummy.
- **`06` distributed loads** (pressures and line loads on an element's face or edge; one element may carry several): one row per packet in two field data tables. `patran:distributed_load` (integers, 20 columns) holds the element's cell index, the set, `LTYPE`, `EFLAG`, `GFLAG`, the six component flags, the eight node flags and the face or edge number; `patran:distributed_load_values` (reals, 54 columns) holds the `PDATA` values in the file's order (the components at the centroid, then at each flagged node), `NaN`-padded.

Packets that name an undefined node or element are skipped with a warning.

## Result files

`meshioplusplus.patran.read(path, results={name: file, ...})` (C++ `read_patran(path, {{name, file}, ...})`; a list of paths is named by their stems) reads Patran 2.5 result files onto the neutral file's mesh (v16.12.0), matched by node and element id:

- **Nodal** (`.nod`, and `.dis` displacements or forces): `NODID` and `NWIDTH` values per node, point data named `name`.
- **Element** (`.els`): `ID`, shape and `NWIDTH` values per element, cell data named `name`.

One column is a 1-D array, more a 2-D one; `NaN` where the file has no record; records naming unknown ids are skipped with a warning. Both the text layout (`(I8,5E13.7)` records for nodal files, `(2I8)` then `(6E13.7)` lines for element ones) and the Fortran unformatted binary one (4-byte words, 4- or 8-byte reals, either byte order, told apart by the first record's length) are read. The MCP `convert` tool takes them as `patran_results`. Beam results files, and writing result files, are not supported.

## Node order

Linear cells, `line3`, `triangle6`, `triangle7`, `quad8`, `quad9`, `tetra10` and `pyramid13` use meshio++'s (VTK's) order (the 7- and 9-node faces add their centre last). **`hexahedron20` and `wedge15` do not.** Patran lists the bottom ring of mid-edge nodes, then the **vertical** mid-edges, then the top ring, where VTK puts the vertical ones last. They use the `"patran"` tables of the [node-ordering registry](../node_ordering.md). The order is the Patran Reference Manual's Element Library. The fixtures are written from its edge lists by `tools/gen_patran_fixtures.py`, and the tests check that every mid-edge node lands on its meshio++ edge midpoint.

## Writing

- Packets `25` (the title card carries the provenance tag), `26`, `01`, `02` and `99`, the loads and boundary conditions (`06`, `07`, `08`, `10`, `11`) from the `patran:` arrays a read produces, plus one `21` per region name. A point region and a cell region with the same name share one component. The component number is the region's tag when positive and unused, else the next free number.
- Coordinates are written `E16.9`, so **they keep ten significant digits**: a round trip is exact to about `1e-9` relative, not bit for bit.
- The element property comes from `patran:property`, else 1.
- Component names are limited to 12 characters (the `A12` name card). A longer or colliding name is truncated, with a numbered suffix if needed, and a warning.
- **Dropped, with a warning and a provenance note:** cells with no Patran shape (vertices, `hexahedron27`, polygons…), side regions, and every data array other than `patran:property` and the load arrays. More than 99,999,999 nodes or elements do not fit the `I8` id fields and are a `WriteError`.

## Errors

A truncated packet, a malformed field, a node or element id defined twice, or an element that names an undefined node raise `ReadError`, naming the line. A component that names undefined entities, or entities of a type other than node and element, is read without them, with a warning. A file with no `99` packet is read to the end with a warning.

## Notes

- **Verification.** Patran, Cubit and ANSA were not available. The reader and writer follow the Patran 2 Neutral File guide and the Element Library, and the generated fixtures are written from them. Real exports were then read (v16.11.0): P3/PATRAN 3.0 and PATRAN 2.5 files from WARP3D's examples (hexahedra) and Tahoe's benchmarks (quadrilaterals with named components; a 2,409-node model with bars, 9-node quadrilaterals and seven components), and a CUBIT 13.2 export (node sets, no end packet). Every cell is positively oriented and the components land on the right entities; the 9-node quadrilaterals are why `quad9` and `triangle7` were added. Three of them are committed under `tests/python/meshes/patran/real/` with their licences. No ANSA or HyperMesh export was found.
- The pyramid shape code (`6`) and the component entity code for a pyramid (`10`) come from Patran's results template, not from the packet 02 documentation.
- Loads and result files (v16.12.0) were checked on WARP3D's own run: its packet `08` constraints, and the Patran 2.5 result files it wrote (text and binary, nodal and element; binary element stresses averaged at the nodes agree with the text nodal stresses of the same step to a correlation of 0.999). A trimmed copy is committed as `real/warp3d_ssy.*`. Packets `06`, `07`, `10` and `11` follow the Reference Manual only: no real file carrying them was found.
- **Not read:** materials, properties, coordinate frames (packet `05`; a load's frame is kept as its number) and MPCs.
