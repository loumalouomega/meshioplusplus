# Conformity check

`check_conformity(mesh)` reports whether a mesh is conforming — whether every two cells meet in a whole shared face, edge or corner with the same point ids on both sides — and, where it is not, which points hang and which facets only partly overlap. It answers the question [`find_interface`](/region_adjacency) leaves to the caller (`"conforming"` or `"proximity"`?), and gives [`refine`](/refine) and [`undo_green`](/undo_green), whose conformity was only discussed in comments, a checker.

```python
import meshioplusplus as mp

mesh = mp.read("assembly.vtu")

checked, report = mp.check_conformity(mesh, return_report=True)
report["conforming"]          # False
report["num_hanging_nodes"]   # 12
report["interface_mode"]      # "proximity": pass that to find_interface
mp.write("flags.vtu", checked)  # conformity:hanging marks the hanging nodes
```

## What it reports

The mesh's volume cells are examined when it has any, else its 2-D cells; lines and vertices are ignored, and a mesh of nothing else reports dimension 0 with a warning. Only the **unmatched facets** — a face one volume cell uses, an edge one 2-D cell uses — can hold a defect: a node hanging on a face is always a corner of an unmatched facet on the finer side, and a face its neighbour matches cannot hold one.

| Defect | Rule | Counted as |
|---|---|---|
| hanging node (3-D) | a point some cell uses lies on a face of a volume cell — inside it or on one of its edges — that does not use it | `num_hanging_nodes`, flag 1 |
| T-junction (2-D) | a point some cell uses lies on the interior of an edge of a 2-D cell that does not use it | `num_t_junctions`, flag 2 |
| partial overlap | two unmatched facets lie in one plane (on one line) and overlap over a positive area (length) without being one facet | `num_partial_overlaps`, per pair |
| unwelded facets | two unmatched facets have coincident corners but different point ids | `num_unwelded_facets`, per pair |

A cell's own nodes are never flagged against it, so the mid-edge nodes of quadratic cells do not hang. Two points closer than the tolerance are one location: a point on a facet's corner is an unwelded duplicate, not a hanging node. A hanging node also shows as partial overlaps — the coarse face against each fine face on it — so a mesh with hanging nodes is never `"conforming"` in either count.

The report also has `num_boundary_facets` (the unmatched facets examined, the domain's boundary included), `num_non_manifold_facets` (facets three or more cells use), `dimension`, `conforming` (none of the four defects), `warnings` (blocks with no known face topology, polyhedra that are not closed) and `interface_mode`:

| `interface_mode` | When | What to do |
|---|---|---|
| `"conforming"` | nothing reported | `find_interface(mode="conforming")` matches the facets by their point ids |
| `"weld"` | only unwelded facets | merge the coincident points ([`clean`](/clean) or `merge(weld=True)`), then `"conforming"` |
| `"proximity"` | hanging nodes, T-junctions or partly overlapping facets | `find_interface(mode="proximity")`, with the coarse side as master |

## Output

The input, bit-identical — points, connectivity, block order, every data array and region — with one `int64` point-data array, `conformity:hanging`: 0, 1 for a hanging node, 2 for a T-junction. `attach_flags=False` leaves it off.

## Tolerance

`tolerance` (default `1e-9`) is a fraction of the bounding-box diagonal of the examined facets' corners: a distance at or below it counts as zero, in every test — a point on a facet, a point on a corner, a facet in another's plane. A midpoint computed as `(a + b) / 2` is within a few units in the last place of its edge, far inside the default; points an algorithm moved off a facet by more than the tolerance are not on it.

## Probes

The tests hold the check to meshes whose defects are known by construction: `refine(closure="balanced")` flags exactly the points `refine:hanging` marks, the `"redgreen"` and `"propagate"` closures report nothing, an [`octree`](/octree) with `transition="polyhedra"` or `"templates"` reports no hanging node and with `"hanging"` flags exactly `octree:hanging`, and `undo_green`'s restored green parents leave T-junctions beside the red children they meet. The notebook `example/python/43_conformity_and_planes.ipynb` runs the same probes on an octree of the Stanford bunny and on two parts meshed apart.

## Core only

The check runs in the compiled core and has no numpy twin: its geometric tests would be a second implementation of the same tolerances to keep in step, and the wheels always carry the core. Without the core, the Python entry point raises `NotImplementedError` naming it.

## CLI

```sh
meshioplusplus check-conformity IN [OUT] [--tolerance T] [--json]
```

Exit status 0 when the mesh is conforming, 1 when it is not and 2 when it could not be checked, as for [`check`](/cli#meshioplusplus-check). `OUT` writes the mesh with `conformity:hanging`.

## Pipeline

```json
{ "Op": "CheckConformity", "Tolerance": 1e-9, "Require": true, "AttachFlags": false }
```

A gate by default: a mesh that is not conforming stops the pipeline with the counts as the error. `"Require": false` only reports (`NumHangingNodes`, `NumTJunctions`, `NumPartialOverlaps`, `NumUnweldedFacets`), and `"AttachFlags": true` writes `conformity:hanging` onto the mesh that flows on. See [pipelines](/pipeline).

## Bindings

C: `mio_check_conformity(mesh, &opts, &report)` with `mio_conformity_opts` (48 bytes, initialized by `mio_conformity_opts_init`) and `mio_conformity_report` (104 bytes, `interface_mode` a `mio_conformity_mode`), both with reserved tails, and `mio_check_conformity_warnings`. Fortran `m%check_conformity([tolerance, attach_flags, report, warnings])`, Julia `check_conformity(m; tolerance, attach_flags)`, R `mio_check_conformity(mesh, tolerance, attach_flags)`, WASM `checkConformity(mesh, {tolerance, attachFlags})` and the MCP tool `check_conformity` (with an optional `output_path`).
