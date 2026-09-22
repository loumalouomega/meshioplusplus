# Tecplot (`.dat`, `.tec`)

The [Tecplot ASCII](http://paulbourke.net/dataformats/tp/) data format: a `VARIABLES` list and one or more finite-element `ZONE`s. meshio++ writes one zone per cell block; on read, `ReadOptions::mTimeStep` (since v11.3.0) selects one **step** of the file's `SOLUTIONTIME`/`STRANDID` timeline — a step being every zone sharing one distinct `SOLUTIONTIME` — see [Selecting a time step](#selecting-a-time-step). A non-transient file (no `SOLUTIONTIME` anywhere) is a single step holding every zone in the file — see [Multiple zones](#multiple-zones).

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

Data itself: `FEBLOCK` reads one variable's full array before moving to the next; `FEPOINT` reads one full-variable-tuple row per node. `X`/`x` (optionally `Y`/`Z`, case-insensitively) become point coordinates; everything else becomes `point_data` or `cell_data` per the location flags. `VARSHARELIST=([a-b]=n)` (1-based, same range syntax) reuses variables `a..b`'s data from zone `n` instead of repeating it; `PASSIVEVARLIST=([a,b,...])` marks variables this zone carries no data for at all (read back as `NaN`); `CONNECTIVITYSHAREZONE=n` reuses zone `n`'s element connectivity outright. Both `VARSHARELIST` and `CONNECTIVITYSHAREZONE` are resolved recursively (a zone may share from a zone that itself shares), with a cycle guarded against.

## Cell types & node ordering

| Tecplot zone type | meshio++ |
|---|---|
| `LINESEG` / `FELINESEG` | `line` |
| `TRIANGLE` / `FETRIANGLE` | `triangle` |
| `QUADRILATERAL` / `FEQUADRILATERAL` | `quad` |
| `TETRAHEDRON` / `FETETRAHEDRON` | `tetra` |
| `BRICK` / `FEBRICK` | `hexahedron` |

On write, `pyramid`/`wedge` degrade to `FEBRICK` (8-node brick), padding with duplicated corner nodes. Write-side node-order tables:

| type | order |
|---|---|
| `pyramid` | `[0,1,2,3,4,4,4,4]` |
| `wedge` | `[0,1,4,3,2,2,5,5]` |

A mesh with several cell blocks — of the same type or mixed — writes one `ZONE` per block (v15.5.0, roadmap §1.1); there is no longer a single-type restriction or a 2D/3D-mixing drop. See [Multiple zones](#multiple-zones).

## Multiple zones

**Write:** each cell block becomes its own `ZONE`, in block order. The first zone carries the point coordinates and every `point_data` array (nodal fields, shared by construction — there is only one point array for the whole mesh); every later zone reuses them via `VARSHARELIST=([1-k]=1)` rather than repeating them. A `cell_data` array is written `CELLCENTERED` on the zones whose block actually carries it, and `PASSIVEVARLIST` on the others (read back there as `NaN`, not zero or an error). A zone's title (`T="..."`) is the name of a `cell` [region](../regions.md) whose entries exactly match that block's cells, else `block_<i>`.

**Read:** a **non-transient** file (no `SOLUTIONTIME` anywhere) is a single step over every zone in the file — concatenated as separate, unwelded cell blocks (points are *not* deduplicated across zones, except where `VARSHARELIST` already ties a zone's coordinates to an earlier one's, in which case that zone's points are recognized as literally the same array rather than appended again). Each zone gets its own `tecplot:zone` cell-data entry (the zone's 0-based position in the file) and its own named `cell` region (the zone's own title, de-duplicated across the step, or `zone_<i>` when the title is empty). A **transient** file's steps are formed the same way, grouped by distinct `SOLUTIONTIME` instead of by "the whole file" — see [Selecting a time step](#selecting-a-time-step).

This closes the roadmap gap where several static zones were previously not a timeline and only the first was read (with a warning); they are now read in full, as multiple parts of one step, matching how SU2's multizone files and VTKHDF's composites are modeled.

## Data mapping

Point/cell variable names are used verbatim as `point_data`/`cell_data` keys (no `tecplot:` prefix); `X`/`Y`/`Z` (or lowercase) are reserved for coordinates and excluded from the data dicts. `cell_data["tecplot:zone"]` names each cell's zone (its 0-based position within the step) for a file with more than one zone. Named [regions](../regions.md) — see [Multiple zones](#multiple-zones) and the [format matrix](../regions.md#the-format-matrix).

## Selecting a time step

A transient Tecplot file marks each `ZONE`'s place in a series with `SOLUTIONTIME=<t>` and, when several distinct series share one file, `STRANDID=<id>` groups zones belonging to the same one. Since v11.3.0 (roadmap §1 tier B1) `ReadOptions::mTimeStep` (0-based, negative counts from the end — the `ResolveTimeStep` contract, see [selective reads](../selective_read.md#reading-one-time-step)) resolves against the **timeline**: every zone sharing the first zone's `STRANDID` (or, when no zone has one, every zone that carries a `SOLUTIONTIME` at all), sorted by `SOLUTIONTIME`, with every zone sharing one distinct `SOLUTIONTIME` value forming one step. Every zone of the resolved step is decoded and concatenated — see [Multiple zones](#multiple-zones).

**Without any `SOLUTIONTIME` at all**, the whole file is the (only) step: every zone in it, concatenated the same way a transient step's zones are (v15.5.0, roadmap §1.1) — before, only the first such zone was read, with the rest silently discarded.

`read_tecplot_metadata` is the native (no data-body decode) counterpart: every `ZONE` header is scanned the same way — including the token-budget walk that locates where one zone's data ends and the next one's header begins — without decoding the tokens themselves. `mNumPoints`/`mCellBlocks` describe every zone of the resolved timeline's first step; `mTimeValues` is the timeline's `SOLUTIONTIME`s in the same sorted order `mTimeStep` indexes into (one value per step, taken from that step's first zone), empty when no zone carries one.

## Quirks & limitations

- Data columns are wrapped at 20 values per line on write.
- A transient *write* path is not implemented (`sequence_write_supports_time("tecplot")` is `false`); every write is a single, non-transient set of zones — see [Multiple zones](#multiple-zones).
- `VARSHARELIST`/`CONNECTIVITYSHAREZONE`/`PASSIVEVARLIST` are read by the C++ core only; the Python reference reader stays single-zone and geometry-only, and now warns (rather than raising `KeyError`) on any zone-header field it does not itself understand, including these three.
- `POINT`/`FEPOINT` data packing has no sharing equivalent (there is nothing to point `VARSHARELIST` at within one interleaved row), so a `POINT`-packed zone that names one anyway is read as if the field were unshared, ordinary data.

## Notes

- `tests/python/meshes/tecplot/quad_zone_comma.tec` / `quad_zone_space.tec` / `quad_zone_multivar.tec` — a single quad zone (`N=4, E=1, ET=QUADRILATERAL`), `FEBLOCK` packing, one cell-centered variable via `VARLOCATION=([4]=CELLCENTERED)`; the three files vary the delimiter style around zone-header keys (comma vs. plain space vs. an extra variable) to exercise the tolerant zone-header parser. `quad_zone_space.tec` in particular has a zone title that is *literally the string* `"VARLOCATION"` with spaced `=` signs; a quote-and-paren-aware tokenizer (v15.5.0, roadmap §1.1) lets the C++ reader handle it directly — earlier it threw and let the Python reader take over.
- The C++ core handles FE meshes (BLOCK/POINT packing, `VARLOCATION`, `VARSHARELIST`, `PASSIVEVARLIST`, `CONNECTIVITYSHAREZONE`), writing and reading every zone of a step.
