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
| anything else (materials, properties, loads, …) | skipped by its card count |

| Shape (`IV`) | Nodes | meshio++ type |
|---|---|---|
| 2 bar | 2 / 3 | `line` / `line3` |
| 3 tri | 3 / 6 | `triangle` / `triangle6` |
| 4 quad | 4 / 8 | `quad` / `quad8` |
| 5 tet | 4 / 10 | `tetra` / `tetra10` |
| 6 pyramid | 5 / 13 | `pyramid` / `pyramid13` |
| 7 wedge | 6 / 15 | `wedge` / `wedge15` |
| 8 hex | 8 / 20 | `hexahedron` / `hexahedron20` |

Any other shape and node count (a 27-node hex, a point element) is skipped with one warning per combination. Node and element ids may have gaps and are renumbered. Elements that no component names are grouped by property into `property_<pid>` cell regions (tag = pid), so a file with no components still has its parts.

## Node order

Linear cells, `line3`, `triangle6`, `quad8`, `tetra10` and `pyramid13` use meshio++'s (VTK's) order. **`hexahedron20` and `wedge15` do not.** Patran lists the bottom ring of mid-edge nodes, then the **vertical** mid-edges, then the top ring, where VTK puts the vertical ones last. They use the `"patran"` tables of the [node-ordering registry](../node_ordering.md). The order is the Patran Reference Manual's Element Library. The fixtures are written from its edge lists by `tools/gen_patran_fixtures.py`, and the tests check that every mid-edge node lands on its meshio++ edge midpoint.

## Writing

- Packets `25` (the title card carries the provenance tag), `26`, `01`, `02` and `99`, plus one `21` per region name. A point region and a cell region with the same name share one component. The component number is the region's tag when positive and unused, else the next free number.
- Coordinates are written `E16.9`, so **they keep ten significant digits**: a round trip is exact to about `1e-9` relative, not bit for bit.
- The element property comes from `patran:property`, else 1.
- Component names are limited to 12 characters (the `A12` name card). A longer or colliding name is truncated, with a numbered suffix if needed, and a warning.
- **Dropped, with a warning and a provenance note:** cells with no Patran shape (vertices, `quad9`, `hexahedron27`, polygons…), side regions, and every data array other than `patran:property`. More than 99,999,999 nodes or elements do not fit the `I8` id fields and are a `WriteError`.

## Errors

A truncated packet, a malformed field, a node or element id defined twice, or an element that names an undefined node raise `ReadError`, naming the line. A component that names undefined entities, or entities of a type other than node and element, is read without them, with a warning. A file with no `99` packet is read to the end with a warning.

## Notes

- **Verification.** Patran, Cubit and ANSA were not available. The reader and writer follow the Patran 2 Neutral File guide and the Element Library, and the fixtures are written from them, so a file from a real Cubit or Patran export has not been read yet. That check stays in the [roadmap](../roadmap.md).
- The pyramid shape code (`6`) and the component entity code for a pyramid (`10`) come from Patran's results template, not from the packet 02 documentation.
- **Not read:** loads and boundary conditions (packets 06–08, 10…), materials and properties, and the result files `.nod`/`.els`/`.dis`.
