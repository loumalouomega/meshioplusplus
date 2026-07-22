# Netgen (`.vol`, `.vol.gz`)

The [Netgen](https://github.com/ngsolve/netgen) neutral mesh format (`.vol`), plain text and keyword-block based, optionally gzip-compressed (`.vol.gz`).

| | |
|---|---|
| **Format name** | `netgen` |
| **Extensions** | `.vol`, `.vol.gz` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.vol")
meshioplusplus.netgen.write("out.vol", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format.

## File structure

The file starts with a literal `mesh3d` line, then keyword blocks in arbitrary order, each `<keyword>\n<count>\n<count lines of data>`; blank lines and `#`-comments are freely interspersed and skipped anywhere. Relevant blocks:

- `dimension` — a single int (2 or 3).
- `geomtype` — a single int; unexpected values just warn.
- `points` — `<N>` then N coordinate rows.
- `pointelements` — dim-0 (vertex) cells.
- `edgesegments`/`edgesegmentsgi`/`edgesegmentsgi2` — dim-1 (line) cells (see the two-line variant below).
- `surfaceelements`/`surfaceelementsgi`/`surfaceelementsuv` — dim-2 cells; the node count is read from a fixed column position per row (variable per row: triangle=3, quad=4, triangle6=6, quad8=8).
- `volumeelements` — dim-3 cells; node count similarly variable (tetra=4, pyramid=5, wedge=6, hexahedron=8, tetra10=10, pyramid13=13, wedge15=15, hexahedron20=20).
- `materials`/`bcnames`/`cd2names`/`cd3names` — codimension-indexed name tables (codim 0/1/2/3 relative to mesh dimension) → `field_data[name] = [idx, edim]`.
- `identifications` — periodic node-pair table, `N` rows of `(node1, node2, identification_type_id)`.
- `identificationtypes` — a single row of `N` type ids.
- `face_colours`, `singular_edge_left/right`, `singular_face_inside/outside`, `singular_points` — read and discarded.

A single integer **cell index** (region/material marker) is stored per element across every block, exposed as `cell_data["netgen:index"]`.

## Cell types & node ordering

Node-count → type, per topological dimension:

| dim | counts → types |
|---|---|
| 0 | 1 → `vertex` |
| 1 | 2 → `line` |
| 2 | 3→`triangle`, 6→`triangle6`, 4→`quad`, 8→`quad8` |
| 3 | 4→`tetra`, 5→`pyramid`, 6→`wedge`, 8→`hexahedron`, 10→`tetra10`, 13→`pyramid13`, 15→`wedge15`, 20→`hexahedron20` |

Full Netgen→meshio++ node permutation table (meshio++[i] = netgen[table[i]]); the meshio++→Netgen direction uses the exact per-entry inverse:

| type | permutation |
|---|---|
| `triangle6` | `[0,1,2,5,3,4]` |
| `quad8` | `[0,1,2,3,4,7,5,6]` |
| `tetra` | `[0,2,1,3]` |
| `tetra10` | `[0,2,1,3,5,7,4,6,9,8]` |
| `pyramid` | `[0,3,2,1,4]` |
| `pyramid13` | `[0,3,2,1,4,7,6,8,5,9,12,11,10]` |
| `wedge` | `[0,2,1,3,5,4]` |
| `wedge15` | `[0,2,1,3,5,4,7,8,6,13,14,12,9,11,10]` |
| `hexahedron` | `[0,3,2,1,4,7,6,5]` |
| `hexahedron20` | `[0,3,2,1,4,7,6,5,10,9,11,8,16,19,18,17,14,13,15,12]` |

(`line`/`triangle`/`quad`/`vertex` use natural order.)

## Data mapping

- `cell_data["netgen:index"]` — the single per-cell region/material marker; Netgen cannot store the field's *name*, so on write meshio++ prefers a `netgen:index` entry if present, else the first integer-dtype cell_data array found.
- `field_data[name] = [idx, edim]` — codimension-domain names (materials, boundary-condition names, co-dim-2/3 names).
- `mesh.info["netgen:identifications"]` (an `(N,3)` int array: node1/node2/type id) and `mesh.info["netgen:identificationtypes"]` (a `(1,N)` int array) — periodic identification data.

## Quirks & limitations

- **1-based** in file; `-1`/`+1` applied on read/write.
- Row layout is column-position-dependent: `surfaceelements` reads its node count from a fixed column and node ids start at another fixed column; `volumeelements` uses different fixed positions — these encode Netgen's distinct per-dimension record schema (`surfnr bcnr domin domout np p1 p2 ...` for surfaces; `matnr np p1 p2 ...` for volumes).
- `edgesegmentsgi2` has a **two-physical-line variant**, triggered by a header line exactly `surf1 surf2 p1 p2` — each cell's data is then split across two lines instead of one. **The C++ reader does not implement this two-line variant** (nor `identifications`, `materials`/`bcnames`/etc., `face_colours`, or `singular_*`) — any of these tokens make the C++ reader throw and defer to Python.
- Only **one** integer cell-data array can be stored per file (a Netgen format limitation, not a meshio++ choice); when reading back, it is always named `"netgen:index"` regardless of its original name.
- The `.vol.gz` gzip container is handled entirely in Python (via `gzip.open`); the C++ reader/writer explicitly refuse the `.gz` suffix.
- `identifications`/`identificationtypes` are stored in `mesh.info`, which has no C++-core representation — any mesh carrying them, or with non-empty `field_data` (materials/bc names), is routed to the Python writer.

## Notes

- `tests/python/meshes/netgen/periodic_1d.vol`, `periodic_2d.vol`, `periodic_3d.vol` — each carries `identifications` data (used to test the `netgen:identifications`/`netgen:identificationtypes` round-trip, which always forces the Python path per the rules above).
- The C++ core handles the common `.vol` path (points, cells, `netgen:index`) for both ascii read/write.
