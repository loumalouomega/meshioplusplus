# Tecplot (`.dat`, `.tec`)

The [Tecplot ASCII](http://paulbourke.net/dataformats/tp/) data format: a `VARIABLES` list and one or more finite-element `ZONE`s. meshio++ writes a **single** FE zone; on read, `ReadOptions::mTimeStep` (since v11.3.0) selects one zone of a transient file's `SOLUTIONTIME`/`STRANDID` timeline — see [Selecting a time step](#selecting-a-time-step).

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

`VARIABLES` supports multi-line continuation and quoted multi-word names (re-joined if a quoted name gets split across whitespace tokens); `X`/`x` and `Y`/`y` must be present or `ReadError`. Every `ZONE` header is parsed (not just the first), tolerating multi-line continuation (it keeps reading as long as the next line's first token fails to parse as a float — i.e. as long as it still looks like header text) and handling a quoted zone title (`T="..."`) plus either `F=` (accepts only `FEPOINT`/`FEBLOCK`) + `ET=`, or `DATAPACKING=`+`ZONETYPE=` (folded into an equivalent `"FE" + DATAPACKING` internal representation). Which zone's *data* is actually decoded is `mTimeStep`'s job — see [Selecting a time step](#selecting-a-time-step).

`VARLOCATION=([a-b]=CELLCENTERED)` (1-based, inclusive ranges, comma- separated `[i]` or `[i-j]` entries) marks which variables are cell-centered; without it, cell-centered-ness is instead inferred from `NV=` (a node- variable count — everything after it is cell-centered).

Data itself: `FEBLOCK` reads one variable's full array before moving to the next; `FEPOINT` reads one full-variable-tuple row per node. `X`/`x` (optionally `Y`/`Z`, case-insensitively) become point coordinates; everything else becomes `point_data` or `cell_data` per the location flags.

## Cell types & node ordering

| Tecplot zone type | meshio++ |
|---|---|
| `LINESEG` / `FELINESEG` | `line` |
| `TRIANGLE` / `FETRIANGLE` | `triangle` |
| `QUADRILATERAL` / `FEQUADRILATERAL` | `quad` |
| `TETRAHEDRON` / `FETETRAHEDRON` | `tetra` |
| `BRICK` / `FEBRICK` | `hexahedron` |

On write, `pyramid`/`wedge`/`hexahedron` all degrade to `FEBRICK` (8-node brick), padding with duplicated corner nodes as needed. Write-side node-order tables (used when the mesh has exactly one supported cell type):

| type | order |
|---|---|
| `pyramid` | `[0,1,2,3,4,4,4,4]` |
| `wedge` | `[0,1,4,3,2,2,5,5]` |

When the mesh has **2 or more** cell types, all are degraded into a single `FEQUADRILATERAL` (if all 2D) or `FEBRICK` (if all 3D) zone using a slightly different padding table (`triangle→[0,1,2,2]`, `tetra→[0,1,2,2,3,3,3,3]`; the degenerate-corner tables for `pyramid`/`wedge`/`quad`/`hexahedron` are unchanged). If the mesh mixes 2D and 3D cell types, the 2D cells are dropped entirely with a warning.

## Data mapping

Point/cell variable names are used verbatim as `point_data`/`cell_data` keys (no `tecplot:` prefix); `X`/`Y`/`Z` (or lowercase) are reserved for coordinates and excluded from the data dicts.

## Selecting a time step

A transient Tecplot file marks each `ZONE`'s place in a series with `SOLUTIONTIME=<t>` and, when several distinct series share one file, `STRANDID=<id>` groups zones belonging to the same one. Since v11.3.0 (roadmap §1 tier B1) `ReadOptions::mTimeStep` (0-based, negative counts from the end — the `ResolveTimeStep` contract, see [selective reads](../selective_read.md#reading-one-time-step)) resolves against the **timeline**: every zone sharing the first zone's `STRANDID` (or, when no zone has one, every zone that carries a `SOLUTIONTIME` at all), sorted by `SOLUTIONTIME`. Only that one zone's data body is decoded — connectivity and coordinates included, since each zone is a complete, self-contained FE block.

**Without any `SOLUTIONTIME` at all**, behaviour is unchanged from before: only the first zone is read. What changed is that a file with more than one such (non-transient) zone now `log::warn`s naming the count, rather than silently discarding the rest with no diagnostic at all — several *static* zones sharing one file is not itself a timeline, and concatenating them is a documented roadmap remainder (§7), not attempted here.

`read_tecplot_metadata` is the native (no data-body decode) counterpart: every `ZONE` header is scanned the same way — including the token-budget walk that locates where one zone's data ends and the next one's header begins — without decoding the tokens themselves. `mNumPoints`/`mCellBlocks` describe the resolved timeline's first zone (transient zones typically share topology); `mTimeValues` is the timeline's `SOLUTIONTIME`s in the same sorted order `mTimeStep` indexes into, empty when no zone carries one.

## Quirks & limitations

- Several **non-transient** zones (no `SOLUTIONTIME` anywhere) are not concatenated; only the first is read, with a warning when there is more than one. A genuinely transient file reads correctly — see [Selecting a time step](#selecting-a-time-step).
- The multi-cell-type write path (degrading everything to one zone via the "order_2" tables) exists **only in the Python writer** — the C++ writer throws `WriteError` if more than one distinct cell type is present, which forces the Python fallback for any such mesh.
- Data columns are wrapped at 20 values per line on write.
- The C++ writer emits a single, non-transient zone; a transient *write* path is not implemented (`sequence_write_supports_time("tecplot")` is `false`).

## Notes

- `tests/python/meshes/tecplot/quad_zone_comma.tec` / `quad_zone_space.tec` / `quad_zone_multivar.tec` — a single quad zone (`N=4, E=1, ET=QUADRILATERAL`), `FEBLOCK` packing, one cell-centered variable via `VARLOCATION=([4]=CELLCENTERED)`; the three files vary the delimiter style around zone-header keys (comma vs. plain space vs. an extra variable) to exercise the tolerant zone-header parser. `quad_zone_space.tec` in particular has a zone title that is *literally the string* `"VARLOCATION"` with spaced `=` signs — an adversarial case the C++ reader throws on cleanly, letting the Python reader take over.
- The C++ core handles FE meshes (BLOCK/POINT packing, `VARLOCATION`), writing a single zone and reading one zone selected by `mTimeStep` out of a file that may carry several.
