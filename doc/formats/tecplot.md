# Tecplot (`.dat`, `.tec`)

The [Tecplot ASCII](http://paulbourke.net/dataformats/tp/) data format: a
`VARIABLES` list and one or more finite-element `ZONE`s. meshio++ only reads and
writes a **single** FE zone.

| | |
|---|---|
| **Format name** | `tecplot` |
| **Extensions** | `.dat`, `.tec` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("field.dat")
meshioplusplus.tecplot.write("out.dat", mesh)
```

`write` takes no keyword arguments.

## File structure

```
VARIABLES = "X" "Y" "Z" "phi" ...
ZONE T="..." N=<nodes> E=<elements> F=FEPOINT|FEBLOCK ET=TRIANGLE|... [VARLOCATION=([a-b]=CELLCENTERED)]
<data — layout depends on F/DATAPACKING and VARLOCATION>
<connectivity, 1-based>
```

`VARIABLES` supports multi-line continuation and quoted multi-word names
(re-joined if a quoted name gets split across whitespace tokens); `X`/`x` and
`Y`/`y` must be present or `ReadError`. `ZONE` header parsing tolerates
multi-line continuation (it keeps reading as long as the next line's first
token fails to parse as a float — i.e. as long as it still looks like header
text) and handles a quoted zone title (`T="..."`) plus either `F=` (accepts
only `FEPOINT`/`FEBLOCK`) + `ET=`, or `DATAPACKING=`+`ZONETYPE=` (folded into
an equivalent `"FE" + DATAPACKING` internal representation). Only the
**first** zone in a file is read — the parser breaks immediately after it,
silently ignoring any subsequent zones.

`VARLOCATION=([a-b]=CELLCENTERED)` (1-based, inclusive ranges, comma-
separated `[i]` or `[i-j]` entries) marks which variables are cell-centered;
without it, cell-centered-ness is instead inferred from `NV=` (a node-
variable count — everything after it is cell-centered).

Data itself: `FEBLOCK` reads one variable's full array before moving to the
next; `FEPOINT` reads one full-variable-tuple row per node. `X`/`x`
(optionally `Y`/`Z`, case-insensitively) become point coordinates; everything
else becomes `point_data` or `cell_data` per the location flags.

## Cell types & node ordering

| Tecplot zone type | meshio++ |
|---|---|
| `LINESEG` / `FELINESEG` | `line` |
| `TRIANGLE` / `FETRIANGLE` | `triangle` |
| `QUADRILATERAL` / `FEQUADRILATERAL` | `quad` |
| `TETRAHEDRON` / `FETETRAHEDRON` | `tetra` |
| `BRICK` / `FEBRICK` | `hexahedron` |

On write, `pyramid`/`wedge`/`hexahedron` all degrade to `FEBRICK` (8-node
brick), padding with duplicated corner nodes as needed. Write-side node-order
tables (used when the mesh has exactly one supported cell type):

| type | order |
|---|---|
| `pyramid` | `[0,1,2,3,4,4,4,4]` |
| `wedge` | `[0,1,4,3,2,2,5,5]` |

When the mesh has **2 or more** cell types, all are degraded into a single
`FEQUADRILATERAL` (if all 2D) or `FEBRICK` (if all 3D) zone using a slightly
different padding table (`triangle→[0,1,2,2]`, `tetra→[0,1,2,2,3,3,3,3]`; the
degenerate-corner tables for `pyramid`/`wedge`/`quad`/`hexahedron` are
unchanged). If the mesh mixes 2D and 3D cell types, the 2D cells are dropped
entirely with a warning.

## Data mapping

Point/cell variable names are used verbatim as `point_data`/`cell_data` keys
(no `tecplot:` prefix); `X`/`Y`/`Z` (or lowercase) are reserved for
coordinates and excluded from the data dicts.

## Quirks & limitations

- Only the **first zone** in a multi-zone file is read; the rest are
  silently ignored, not merged or errored on.
- The multi-cell-type write path (degrading everything to one zone via the
  "order_2" tables) exists **only in the Python writer** — the C++ writer
  throws `WriteError` if more than one distinct cell type is present, which
  forces the Python fallback for any such mesh.
- Data columns are wrapped at 20 values per line on write.

## Notes

- `tests/meshes/tecplot/quad_zone_comma.tec` / `quad_zone_space.tec` /
  `quad_zone_multivar.tec` — a single quad zone (`N=4, E=1,
  ET=QUADRILATERAL`), `FEBLOCK` packing, one cell-centered variable via
  `VARLOCATION=([4]=CELLCENTERED)`; the three files vary the delimiter style
  around zone-header keys (comma vs. plain space vs. an extra variable) to
  exercise the tolerant zone-header parser. `quad_zone_space.tec` in
  particular has a zone title that is *literally the string* `"VARLOCATION"`
  with spaced `=` signs — an adversarial case the C++ reader throws on
  cleanly, letting the Python reader take over.
- The C++ core handles single-zone FE meshes (BLOCK/POINT packing,
  `VARLOCATION`).
