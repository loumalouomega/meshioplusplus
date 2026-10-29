# Intersections

`find_intersections(mesh)` reports where a surface intersects itself, and `find_intersections(mesh, other)` where two surfaces interpenetrate or touch: the intersecting cell pairs, the segments along which their triangles cross as a line mesh, and a count per cell. It is the precondition check that [`repair`](/repair), [`shrinkwrap`](/shrinkwrap), [`tetrahedralize`](/tetrahedralize) and [`remesh_volume`](/remesh_volume) assume and do not make themselves.

```python
import meshioplusplus as mp

skin = mp.read("part.stl")

flagged, report = mp.find_intersections(skin, return_report=True)
report["clean"]          # False
report["num_pairs"]      # 37 intersecting cell pairs
report["pairs"][:3]      # (n, 2) int64 global cell ids
mp.write("flags.vtu", flagged)                # intersections:count per cell
mp.write("segments.vtu", report["segments"])  # where the triangles cross

# Two parts: interpenetration, touching included.
_, both = mp.find_intersections(housing, insert, return_report=True)
```

## What is tested

The surfaces are reduced to triangles the way every surface operation reduces them — triangle blocks as they are, quads and polygons fanned from their first corner — and lines and vertices are ignored. Volume, polyhedron and higher-order blocks are refused: run [`extract_surface`](/extract_surface) or `linearize` first.

Candidate pairs are the triangles whose bounding boxes meet, found over a bucket grid sized like [`check_conformity`](/check_conformity)'s (a cell as large as the triangles' mean extent, never finer than 1/256 of the diagonal). Each candidate is decided by an exact triangle–triangle test built on Shewchuk's `orient3d` and `orient2d` predicates (the core's `detail/predicates`, in the spirit of Guigue and Devillers, *Journal of Graphics Tools* 8(1), 2003). The answer is therefore the exact one for the input's double coordinates: there is no tolerance, and no class of "possible" hits to review by hand.

| Input | What counts |
|---|---|
| one surface | Two triangles of one cell (a quad's fan) are never paired. Triangles that share corners — bitwise-equal coordinates, so an unwelded seam is not mistaken for a crossing — count only when they meet beyond what they share: on a shared edge, when they fold flat onto each other; on a shared corner, when an edge of one crosses the other or runs into it from the corner. Any other contact counts, touching included. |
| two surfaces | Every contact counts, touching included: two parts resting on each other are reported. Measure the gap between parts that may touch with [`thickness`](/thickness) instead. |

## What comes back

| Key | Meaning |
|---|---|
| `pairs` | the intersecting cell pairs as an `(n, 2)` Int64 array of global (block-major) cell ids, the cell of `mesh` first (for one surface, the lower id first), sorted; `max_pairs` caps it |
| `segments` | a line mesh, one `line` per crossing segment with two points of its own (not welded), and `intersections:cell_a` / `intersections:cell_b` cell data naming the cells it came from |
| `other` | `other`, bit-identical, with its own `intersections:count` (two surfaces only) |
| `num_pairs` | cell pairs, every one (not capped) |
| `num_triangle_pairs` | triangle pairs that meet |
| `num_coplanar` | of those, pairs lying in one plane and overlapping: they meet over an area, not a segment, so they are flagged and counted but give no segment |
| `num_segments` | cells of `segments`: one per triangle pair that crosses along a segment; pairs that only touch at a point are counted and flagged but give none |
| `num_cells_a`, `num_cells_b` | cells of each surface that intersect anything |
| `truncated`, `clean` | the pair list was capped; nothing intersects |

The returned mesh is the input, bit-identical — points, connectivity, block order, data and regions — with one Int64 cell-data array, `intersections:count`: how many cells each cell intersects, 0 for a clean cell. `attach_flags=False` leaves it off and `segments=False` skips the segment mesh.

The segment end points are computed in double — where the edges that cross each triangle's plane cross it, clipped to the line both planes share — so they are as accurate as the coordinates allow, while the decision behind them is exact. Two overlapping spheres give a closed ring of segments on the circle where they meet.

Deterministic: candidates are collected per triangle and sorted, and the results are reduced in triangle order, so the output is the same on every parallel backend and thread count.

## Core only

The check runs in the compiled core and has no numpy twin: a second exact triangle–triangle test would be a second implementation of the same predicates to keep in step. Without the core, the Python entry point raises `NotImplementedError` naming it.

## CLI

```sh
meshioplusplus find-intersections IN [OUT] [--other FILE] [--segments FILE]
                                  [--region R] [--region-other R] [--max-pairs N] [--json]
```

Exit status 0 when nothing intersects, 1 when something does and 2 when the input could not be checked, as for [`check-conformity`](/check_conformity#cli). `OUT` writes the mesh with `intersections:count`, `--segments` the segment mesh; `--json` prints the counts, `pairs` and `warnings`.

## Pipeline

```json
{ "Op": "FindIntersections", "Other": "", "Region": "", "RegionOther": "", "Require": true, "AttachFlags": false }
```

A gate by default: a surface that intersects itself (or meets `Other`, a mesh file) stops the pipeline with the pair count as the error. `"Require": false` only reports (`NumPairs`, `NumTrianglePairs`, `NumCoplanar`), and `"AttachFlags": true` writes `intersections:count` onto the mesh that flows on. See [pipelines](/pipeline).

## Bindings

C: `mio_find_intersections(a, b_or_NULL, &opts, &report)` returns an opaque `mio_find_intersections_result`; take its meshes with `mio_find_intersections_result_take_mesh`, `_take_mesh_b` and `_take_segments`, borrow the flattened pairs with `mio_find_intersections_result_pairs`, and free it with `mio_find_intersections_result_free`. `mio_find_intersections_opts` (64 bytes, initialized by `mio_find_intersections_opts_init`) and `mio_find_intersections_report` (96 bytes) both have reserved tails; `mio_find_intersections_warnings` returns the warnings. Fortran `mio_find_intersections(mesh[, other, region, region_other, max_pairs, segments, other_flagged, pairs, report])`, Julia `find_intersections(m[, other]; region, region_other, segments, max_pairs)`, R `mio_find_intersections(mesh, other, ...)` (the last three with 1-based pairs), WASM `findIntersections(mesh, {other, region, regionOther, segments, maxPairs})` and the MCP tool `find_intersections` (with optional `output_path` and `segments_path`).
