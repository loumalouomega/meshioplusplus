# Nastran (`.bdf`, `.fem`, `.nas`)

The [MSC/NX Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) bulk-data format: fixed-width card entries (`GRID`, `CTRIA3`, `CTETRA`, `CHEXA`, …) in small-field, large-field, or free (comma-separated) layout.

| | |
|---|---|
| **Format name** | `nastran` |
| **Extensions** | `.bdf`, `.fem`, `.nas` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.bdf")
meshioplusplus.nastran.write("out.bdf", mesh,
    point_format="fixed-large",   # "fixed-small", "fixed-large", or "free"
    cell_format="fixed-small",
)
```

- **`point_format`** / **`cell_format`** — the field layout for `GRID` and element cards, respectively.

## File structure

Parsing starts at a line beginning `"BEGIN BULK"` and stops at `"ENDDATA"`. Three field layouts are auto-detected per line:

- **Small field (fixed)**: 10 fields of exactly 8 characters each.
- **Large field (fixed)**: 8 + 4×16 + 8 characters; the keyword carries a trailing `*` (e.g. `GRID*`), and the record continues onto a second physical line starting with `*` — the two 16-char halves of a continued field are re-merged before further parsing.
- **Free field**: comma-separated (detected by the presence of a `,` on the line).

Continuation is signalled either by an explicit `+`/`*` marker as the first character of the next line, or *implicitly* when the last field of one line and the first field of the next are both blank.

`GRID`/`GRID*`: `[id, ref(optional), x, y, z]`. Coordinates may use Nastran's compressed-exponent float notation (e.g. `1.5+1`, `.7E1`) — decoded by inserting an implicit `e` before a bare `+`/`-` not already following `e`/`E`.

Element cards: `[element_id, ref(optional), node_ids...]`, except `CBAR`/`CBEAM`/`CBUSH`/`CBUSH1D`/`CGAP`, which only take the first 2 node ids from the card — a 3rd field present on those cards (an orientation vector or grid id) is **discarded on read**. Second-order solids (`CTETRA`, `CPYRA`, `CPENTA`, `CHEXA`) are auto-upgraded to their 10/13/15/20-node meshio++ counterparts whenever a card lists more node ids than the linear element's base count.

## Cell types & node ordering

| Nastran | meshio++ | Nastran | meshio++ |
|---|---|---|---|
| `CBEAM`, `CBUSH`, `CBUSH1D`, `CROD`, `CGAP`, `CBAR` | `line` | `CTETRA` | `tetra` |
| `CTRIAR`, `CTRIA3` | `triangle` | `CTETRA_`* | `tetra10` |
| `CTRAX6`, `CTRIAX6`, `CTRIA6` | `triangle6` | `CPYRAM`, `CPYRA` | `pyramid` |
| `CQUADR`, `CSHEAR`, `CQUAD4` | `quad` | `CPYRA_`* | `pyramid13` |
| `CQUAD8` | `quad8` | `CPENTA` | `wedge` |
| `CQUAD9` | `quad9` | `CPENTA_`* | `wedge15` |
| `CELAS1` | `vertex` | `CHEXA` | `hexahedron` |
| | | `CHEXA_`* | `hexahedron20` |

(`*` = "fictive" type names used internally to represent the auto-detected second-order upgrade of the base solid card.)

Node-order permutations (involutions — the same array is used in both directions unless noted):

| type | permutation |
|---|---|
| `triangle6`/`CTRAX6`/`CTRIAX6` | to-VTK `[0,2,4,1,3,5]`, to-Nastran `[0,3,1,4,2,5]` |
| `hexahedron20` (`CHEXA_`) | `[0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]` |
| `wedge15` (`CPENTA_`) | `[0,1,2,3,4,5,6,7,8,12,13,14,9,10,11]` |

## Data mapping

- `point_data["nastran:ref"]` — the GRID card's optional reference field.
- `cell_data["nastran:ref"]` — the element card's optional reference field, one array per block.
- `mesh.points_id` / `mesh.cells_id` — mesh-level attributes (not data-dict entries) holding the original GRID/element ids; set only by the Python reader.

## Quirks & limitations

- The **16-character float encoding** used for large-field `GRID*` cards is the trickiest part of this format: the Python writer uses `np.format_float_scientific(precision=11)` then swaps `e→E`; the C++ writer instead searches increasing precision (0 through 11) for the shortest string that round-trips exactly via `strtod`, then trims trailing mantissa zeros. Both target the same 16-character field limit but are not guaranteed to produce byte-identical text.
- `CBAR`/`CBEAM`/`CBUSH`/`CBUSH1D`/`CGAP`'s 3rd node/orientation field is always dropped on read — a genuinely lossy round-trip for those 1D element types.
- The 2nd-order-solid upgrade is a **heuristic** ("more node ids on the card than the linear element's count ⇒ treat as quadratic"), not a version flag — any card with unexpected extra trailing tokens could be misclassified.
- **The C++ reader is sentinel-gated**: it only parses files carrying the exact literal comment string the C++ writer itself emits (`"meshioplusplus-cpp-nastran"`) as the first `$` comment line. Any real-world Nastran file — including this project's own reference `.fem` fixtures — lacks that sentinel and is therefore always parsed by the (more permissive, general-purpose) Python reader. This is the single most consequential interop rule for this format. Since v10.15.0 the C++ writer also emits the standard `Written by meshio++ v<version>` provenance line (see [`doc/formats.md`](../formats.md#provenance)) as a second `$` comment line, *after* the sentinel — the sentinel stays first, so the reader's gate is unaffected. The Python writer emits that same tag but never the sentinel, which is why nastran is the one format where the two engines' output legitimately differs beyond the tag itself.
- The C++ writer only emits the `fixed-large`/`fixed-small` point/cell format combination; the shim only attempts the C++ path for exactly that combination and only when no `nastran:ref` data is present.
- Points are force-promoted to 3D on write if given 2D, with a warning.

## Notes

- `tests/python/meshes/nastran/cylinder.fem` and `cylinder_cells_first.fem` — HyperMesh/Optistruct-generated meshes matching `tests/python/meshes/med/cylinder.med` geometrically; used for a point-sum and per-type connectivity-sum checksum (`{line:241, triangle:171, quad:721, pyramid:1180, tetra:5309}`).
