# Tecplot (`.dat`, `.tec`, `.plt`)

The [Tecplot ASCII](http://paulbourke.net/dataformats/tp/) data format, a `VARIABLES` list and one or more finite-element or ordered `ZONE`s, and its binary form `.plt` (read only, since v16.10.0; see [Binary `.plt`](#binary-plt)). meshio++ writes one zone per cell block; on read, `ReadOptions::mTimeStep` (since v11.3.0) selects one **step** of the file's `SOLUTIONTIME`/`STRANDID` timeline — a step being every zone sharing one distinct `SOLUTIONTIME` — see [Selecting a time step](#selecting-a-time-step). A non-transient file (no `SOLUTIONTIME` anywhere) is a single step holding every zone in the file — see [Multiple zones](#multiple-zones).

| | |
|---|---|
| **Format name** | `tecplot` |
| **Extensions** | `.dat`, `.tec`, `.plt` (a file starting with `#!TDV` is binary whatever its extension; a `.dat` file that opens as an MSC Marc input deck is read as [`marc`](./marc.md#which-dat-is-marc-s) instead, v16.8.0) |
| **Read / Write** | ✓ / ✓ (ASCII); ✓ / — (`.plt`) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("field.dat")
flow = meshioplusplus.read("flow.plt", time_step=-1)  # binary, last step
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

`VARIABLES` supports multi-line continuation and quoted multi-word names (re-joined if a quoted name gets split across whitespace tokens); `X`/`x` and `Y`/`y` must be present or `ReadError`. Since v16.12.0 a name with a unit suffix counts as a coordinate when no exact one exists (`X(M)`, `X [m]`, as older Tecplot files label them), and data values may be separated by commas as well as blanks. Every `ZONE` header is parsed (not just the first), tolerating multi-line continuation (it keeps reading as long as the next line's first token fails to parse as a float — i.e. as long as it still looks like header text) and handling a quoted zone title (`T="..."`) plus either the old `F=` (`FEPOINT`/`FEBLOCK` with `ET=`, or `POINT`/`BLOCK` for an ordered zone) or `DATAPACKING=`+`ZONETYPE=`. As the Data Format Guide says, a zone without `ZONETYPE` is `ORDERED` and one without `DATAPACKING` is `BLOCK` (since v16.10.0; before, `POINT` was assumed). An ordered zone's size is `I=`, `J=`, `K=` (each 1 when absent). Which zone's *data* is actually decoded is `mTimeStep`'s job — see [Selecting a time step](#selecting-a-time-step).

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
| `ORDERED`, one dimension > 1 | `line` |
| `ORDERED`, two dimensions > 1 | `quad` |
| `ORDERED`, `I`, `J`, `K` > 1 | `hexahedron` |
| `FEPOLYGON` | `polygon` (ragged) |
| `FEPOLYHEDRON` | `polyhedron<N>`, one block per distinct node count |

Ordered zones (since v16.10.0, ASCII and binary) become regular cells over their dimensions longer than one, with `i` varying fastest and corners in VTK order. Their cell-centred variables have `(I-1)(J-1)(K-1)` values, but a binary file stores them over the node dimensions with the last dimension longer than one shortened by one; the ghost values at the ends are dropped. The writer writes FE zones only, so an ordered zone read in is written back as an FE zone of the same cells.

Face-based zones (`FEPOLYGON`/`FEPOLYHEDRON`, since v16.11.0) are described in [Face-based zones](#face-based-zones).

On write, `pyramid`/`wedge` degrade to `FEBRICK` (8-node brick), padding with duplicated corner nodes. Write-side node-order tables:

| type | order |
|---|---|
| `pyramid` | `[0,1,2,3,4,4,4,4]` |
| `wedge` | `[0,1,4,3,2,2,5,5]` |

A mesh with several cell blocks — of the same type or mixed — writes one `ZONE` per block (v15.5.0, roadmap §1.1); there is no longer a single-type restriction or a 2D/3D-mixing drop. See [Multiple zones](#multiple-zones).

## Face-based zones

Since v16.11.0 `FEPOLYGON` and `FEPOLYHEDRON` zones are read (ASCII and `.plt`) and written (ASCII). Such a zone has no element connectivity: its header gives `FACES`, `TOTALNUMFACENODES` (optional for polygons, where every face is an edge of two nodes), `NUMCONNECTEDBOUNDARYFACES` and `TOTALNUMBOUNDARYCONNECTIONS`, and after its data comes a face map — the node count of each face (polyhedra only), each face's nodes, then the element on the left and on the right of each face. Its data must be `DATAPACKING=BLOCK`.

Cells are rebuilt from the face map using the Data Format Guide's right-hand rule: a polyhedral face's right-hand normal points at its right element, so the face is outward for its left element and is reversed for its right one; a polygonal face (an edge) walked from its first node to its second has its left element on the left. So every rebuilt polyhedron has outward faces and every polygon is wound counter-clockwise when the file follows the rule. A polygon's ring starts at its first edge in the file and is reversed only when most of its edges disagree with that direction. A neighbour of `0` (ASCII) or `-1` (`.plt`) means none. A negative neighbour points into the boundary-connection section, i.e. at an element of another zone: that face still bounds this zone's cell, and the cross-zone link itself is skipped.

A polygonal zone is one ragged `polygon` block. A polyhedral zone is one `polyhedron<N>` block per distinct node count, in the order the counts first occur (the naming every polyhedral reader uses), so its cells can come out reordered: cell-centred variables, `tecplot:zone` and the zone's region follow each cell into its block. A zone sharing another zone's face map through `CONNECTIVITYSHAREZONE` must match its type and sizes. Cells are not reconstructed into `tetra`/`hexahedron`: a hexahedron in a polyhedral zone stays a `polyhedron8`.

On write, a `polygon*` block becomes an `FEPOLYGON` zone and a `polyhedron*` block an `FEPOLYHEDRON` zone. Each face shared by two cells is written once, with the cell that lists it first as its left element and the second as its right, so a mesh whose faces are wound outward writes a face map that follows the right-hand rule. No boundary connections are written.

Refused with `ReadError`: a missing `FACES` (or `TOTALNUMFACENODES` for polyhedra), `POINT` packing, a face with fewer than two nodes, node or element numbers out of range, an element with fewer than four faces (polyhedra) or whose edges do not close one loop (polygons), and a face map cut short.

## Multiple zones

**Write:** each cell block becomes its own `ZONE`, in block order. The first zone carries the point coordinates and every `point_data` array (nodal fields, shared by construction — there is only one point array for the whole mesh); every later zone reuses them via `VARSHARELIST=([1-k]=1)` rather than repeating them. A `cell_data` array is written `CELLCENTERED` on the zones whose block actually carries it, and `PASSIVEVARLIST` on the others (read back there as `NaN`, not zero or an error). A zone's title (`T="..."`) is the name of a `cell` [region](../regions.md) whose entries exactly match that block's cells, else `block_<i>`.

**Read:** a **non-transient** file (no `SOLUTIONTIME` anywhere) is a single step over every zone in the file — concatenated as separate, unwelded cell blocks (points are *not* deduplicated across zones, except where `VARSHARELIST` already ties a zone's coordinates to an earlier one's, in which case that zone's points are recognized as literally the same array rather than appended again). Each zone gets its own `tecplot:zone` cell-data entry (the zone's 0-based position in the file) and its own named `cell` region (the zone's own title, de-duplicated across the step, or `zone_<i>` when the title is empty). A **transient** file's steps are formed the same way, grouped by distinct `SOLUTIONTIME` instead of by "the whole file" — see [Selecting a time step](#selecting-a-time-step).

This closes the roadmap gap where several static zones were previously not a timeline and only the first was read (with a warning); they are now read in full, as multiple parts of one step, matching how SU2's multizone files and VTKHDF's composites are modeled.

## Binary `.plt`

Since v16.10.0 the reader decodes the binary format of the Data Format Guide's appendix A, version `#!TDV112` (Tecplot 360 2009 onwards), in either byte order; since v16.12.0 the older versions as well (see [Older versions](#older-versions)). A binary file decodes into the same zone model as an ASCII one, so everything on this page (zones and regions, sharing, passive variables, time steps) applies to both, and a `.plt` reads identically to the ASCII file it was made from. In detail:

- Header: the title, variable names and zone records (FE and ordered zone types, variable locations, `STRANDID` and `SOLUTIONTIME`). Geometry, text, custom-label, user and auxiliary-data records are skipped. User-defined face-neighbour connections are skipped too.
- Data: `float`, `double`, `int32`, `int16` and `byte` values; bit-packed variables are refused. Variable sharing, passive variables and connectivity sharing are honoured. FE connectivity is zero-based in the file (1-based in Tecplot 7 files, see [Older versions](#older-versions)).
- Face-based zones: the face map is stored as face-node offsets (polyhedra only), face nodes, left and right elements (0-based, `-1` for none), then the boundary connections (whose face count in the zone header includes one extra for "no neighbour"). See [Face-based zones](#face-based-zones).
- Anything else raises `ReadError`: another version, a zone of an unknown type, or a file cut short.

### Older versions

Since v16.12.0 the older binary versions are read too: `#!TDV71` and `#!TDV75` (Tecplot 7) and `#!TDV100` to `#!TDV111`, as well as `#!TDV112` and `113`. There is no public specification for them; the layouts follow VisIt's `TecplotFile.C` and were checked against the Tecplot example files VisIt tests with (see `tests/python/meshes/tecplot/visit/`). What changes by version:

- The FileType word appears only from 111. The zone record gains its parent zone from 107 and the user face-neighbour words from 108; the point/block packing word follows the zone type until 112, so a zone before 112 can be point-packed (one record of every variable per node).
- Up to 102 the zone type follows an unused `-1`, with no strand or solution time; the data section has neither the passive-variable list nor the min/max pairs.
- From 104 a cell-centred variable of an ordered zone drops the ghost values of its slowest direction (TecIO's `tecxxx.cpp`); before 104 it keeps them.
- Tecplot 7 files keep only the zone kind (block or point, ordered or finite element), the sizes and the element type (triangle, quadrilateral, tetrahedron, brick) in the zone record. The data section has one extra word before the formats, no sharing, passive or min/max lists, and one more word before the connectivity, which is 1-based. Their text and geometry records have a 2-D anchor and no clipping word.

Verified: 71 and 75 (both byte orders, block and point packing, ordered and FE zones, text and geometry records), 106, 107 and 108 read in full, and every file with an ASCII twin reads as it does. Versions 100 to 105 and 109 to 111 had no sample: they are read as the rules above say.

`.szplt` (SZL) files are a different, undocumented format and are not read; save them as `.plt` in Tecplot. The Python reference reader (`meshioplusplus.tecplot`) reads `.plt` too, value for value like the core.

## Data mapping

Point/cell variable names are used verbatim as `point_data`/`cell_data` keys (no `tecplot:` prefix); `X`/`Y`/`Z` (or lowercase) are reserved for coordinates and excluded from the data dicts. A variable that is nodal in some zones of a step and cell-centred in others is both a `point_data` and a `cell_data` array, each `NaN` where that zone stores it the other way (since v16.10.0). A zone reuses an earlier zone's points only when all its nodal variables, coordinates included, are shared from that zone; otherwise its points are its own. `cell_data["tecplot:zone"]` names each cell's zone (its 0-based position within the step) for a file with more than one zone. Named [regions](../regions.md) — see [Multiple zones](#multiple-zones) and the [format matrix](../regions.md#the-format-matrix).

## Selecting a time step

A transient Tecplot file marks each `ZONE`'s place in a series with `SOLUTIONTIME=<t>` and, when several distinct series share one file, `STRANDID=<id>` groups zones belonging to the same one. Since v11.3.0 (roadmap §1 tier B1) `ReadOptions::mTimeStep` (0-based, negative counts from the end — the `ResolveTimeStep` contract, see [selective reads](../selective_read.md#reading-one-time-step)) resolves against the **timeline**: every zone sharing the first zone's `STRANDID` (or, when no zone has one, every zone that carries a `SOLUTIONTIME` at all), sorted by `SOLUTIONTIME`, with every zone sharing one distinct `SOLUTIONTIME` value forming one step. Every zone of the resolved step is decoded and concatenated — see [Multiple zones](#multiple-zones).

**Without any `SOLUTIONTIME` at all**, the whole file is the (only) step: every zone in it, concatenated the same way a transient step's zones are (v15.5.0, roadmap §1.1) — before, only the first such zone was read, with the rest silently discarded.

`read_tecplot_metadata` is the native (no data-body decode) counterpart: every `ZONE` header is scanned the same way — including the token-budget walk that locates where one zone's data ends and the next one's header begins — without decoding the tokens themselves. `mNumPoints`/`mCellBlocks` describe every zone of the resolved timeline's first step; `mTimeValues` is the timeline's `SOLUTIONTIME`s in the same sorted order `mTimeStep` indexes into (one value per step, taken from that step's first zone), empty when no zone carries one.

## Quirks & limitations

- Data columns are wrapped at 20 values per line on write.
- A transient *write* path is not implemented (`sequence_write_supports_time("tecplot")` is `false`); every write is a single, non-transient set of zones — see [Multiple zones](#multiple-zones).
- The Python reference reader decodes every zone, time step and sharing list as the C++ core does (since v16.10.0; it was single-zone and geometry-only before), so `time_step` no longer needs the core.
- A `VARSHARELIST` entry without a zone number shares from the previous zone, as the Data Format Guide says.
- `POINT`/`FEPOINT` data packing has no sharing equivalent (there is nothing to point `VARSHARELIST` at within one interleaved row), so a `POINT`-packed zone that names one anyway is read as if the field were unshared, ordinary data.

## Notes

- `tests/python/meshes/tecplot/quad_zone_comma.tec` / `quad_zone_space.tec` / `quad_zone_multivar.tec` — a single quad zone (`N=4, E=1, ET=QUADRILATERAL`), `FEBLOCK` packing, one cell-centered variable via `VARLOCATION=([4]=CELLCENTERED)`; the three files vary the delimiter style around zone-header keys (comma vs. plain space vs. an extra variable) to exercise the tolerant zone-header parser. `quad_zone_space.tec` in particular has a zone title that is *literally the string* `"VARLOCATION"` with spaced `=` signs; a quote-and-paren-aware tokenizer (v15.5.0, roadmap §1.1) lets the C++ reader handle it directly — earlier it threw and let the Python reader take over.
- The C++ core handles FE and ordered zones (BLOCK/POINT packing, `VARLOCATION`, `VARSHARELIST`, `PASSIVEVARLIST`, `CONNECTIVITYSHAREZONE`), writing and reading every zone of a step.
- `tests/python/meshes/tecplot/plt/` holds `.plt` files written by Tecplot's TecIO library, each beside the ASCII file with the same content (`tools/gen_tecplot_plt_fixtures.py`). `poly_2d` and `poly_3d` are the face-based fixtures (TecIO's `TECPOLYFACE142`/`TECPOLYBCONN142`), including a boundary connection between two polyhedral zones. No Tecplot licence was available, so no `#!TDV112` file written by `preplot` or Tecplot 360 itself has been read, and ParaView's Tecplot reader does not read face-based zones, so it could not serve as a second opinion on them.
- `tests/python/meshes/tecplot/visit/` holds Tecplot's own example files in versions 71, 75, 106 and 108 from VisIt's test data (BSD-3-Clause), each checked against its ASCII twin or a twin in another version (v16.12.0).
