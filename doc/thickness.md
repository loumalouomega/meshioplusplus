# Thickness

`thickness(mesh)` measures the wall thickness of a surface at every facet, and `thickness(mesh, other)` the clearance from each facet to a second part, by a ray along the facet's normal to the opposite surface. Thin walls before a cast or a print, the gap between two parts of an assembly, and the facets where a boundary layer would run into the far side are the questions it answers.

```python
import meshioplusplus as mp

shell = mp.read("bracket.stl")

measured, report = mp.thickness(shell, threshold=1.5, return_report=True)
report["min"], report["min_cell"]   # the thinnest wall and where it is
report["num_below"]                 # facets thinner than 1.5
mp.write("thickness.vtu", measured) # thickness:distance per cell

# Clearance between two parts, along each facet's outward normal.
_, gap = mp.thickness(housing, insert, return_report=True)
```

## The ray

From each 2-D cell's area-weighted centroid a ray is cast along the cell's unit normal (the sum of its fan triangles' normals, so a quad or a polygon has one) through a bounding-volume hierarchy of the target surface's triangles. The kernel is the core's one ray–triangle test (Möller and Trumbore, *Journal of Graphics Tools* 2(1), 1997), shared with [`octree_mesh`](/octree_mesh) and [`boundary_layer`](/boundary_layer); a ray through a shared edge or corner does not slip between its triangles.

| `direction` | The ray | Stops at | Measures |
|---|---|---|---|
| `"inward"` (default for one surface) | against the normal | the first triangle it leaves through, one whose normal points along the ray: the far side of the wall | wall thickness |
| `"outward"` (default for two surfaces) | along the normal | the first triangle it enters, one whose normal faces it: the next wall across a gap | clearance, to the same surface or to `other` |

The surface must be consistently oriented with outward normals, as a closed shell [`repair`](/repair) has oriented is; a facet wound the other way measures the opposite side. A cell's own triangles are never hit, so a facet does not measure itself, and in a concave corner a neighbouring wall is legitimately the nearest. A ray that finds nothing within `max_distance` (default 0: the diagonal of both surfaces' bounding box) is a miss.

## Beyond one ray

A ray is exact between parallel walls and too long by 1/cos θ where the walls are not, and it needs a consistently oriented surface, which a scan often is not (roadmap 7.5).

- **`method="sphere"`**: from the ray's midpoint, the distance to the nearest target triangle outside the facet's own one-ring (the triangles sharing a vertex with it), doubled, and never more than the ray. Between parallel walls it is the ray; where the faces diverge it is the perpendicular distance to the opposite face (one extra closest-point query per facet, filtered as `deviation`'s far-side rule is). On a wedge whose faces meet at 30 degrees, a point at distance `a` from the edge reads `a sin 30°` by the sphere and `a tan 30°`, 1/cos 30° longer, by the ray. The maximal-inscribed-sphere definition collapses at sharp edges, so it is not offered. A ray that misses is a miss here too.
- **`method="cone"`**: the shortest hit of a fixed fan of `cone_rays` rays (16) within `cone_angle` degrees (15) of the normal: the central ray, then a Fibonacci spiral in the facet's frame, so the same facet always casts the same rays. It finds a chamfer's or a tapered rib's narrow spot, which the straight ray passes.
- **`two_sided=True`**: scan winding is unreliable, so rays run both ways along the normal line and either face of every triangle counts; the shorter hit wins. On a hollow box with every third facet flipped the wall reads 1 everywhere it is seen, where the one-sided rays of the flipped facets miss. In a cavity narrower than the wall it reads the cavity. With `max_distance` 0 the ray is then a fifth of the bounding-box diagonal, so a facet with no wall in reach is NaN rather than the distance to the far side of the part.
- **`max_normal_deviation`** (degrees, default 90): a triangle more edge-on to the ray than this is stepped over, so a cone ray grazing a convex edge does not escape through it; 90 accepts any triangle facing the right way, today's rule, and 60 suits a scan.

The defaults (`"ray"`, one-sided, 90 degrees) give the output they always did.

## What comes back

The input, bit-identical, with two cell-data arrays: Float64 `thickness:distance` (NaN for a miss) and Int64 `thickness:hit_cell`, the global cell of the target it reached (-1 for a miss). Cells that are not measured, such as lines, read as misses. `attach=False` leaves both off.

| Key | Meaning |
|---|---|
| `direction` | the direction the rays ran, `"inward"` or `"outward"` |
| `num_cells`, `num_hit`, `num_missed` | cells measured, rays that reached a wall, rays that did not |
| `min`, `max`, `mean` | over the hits (0 when nothing was hit) |
| `min_cell`, `min_point` | the global cell of the least distance and where its ray started |
| `num_below` | hits below `threshold`, when one is given |
| `max_distance` | the ray length used |
| `method` | `"ray"`, `"sphere"` or `"cone"` |
| `p2`, `p95` | the 2nd and 95th percentiles of the hit distances, the range a colour ramp is set from |

A ray measures along the normal, not to the nearest point: on a wall of uniform thickness the two agree, at a taper the ray reads the distance across the wall along the facet's own normal. A facet with no area has no normal and reads as a miss, with a warning.

Deterministic: each ray stops at the nearest accepted triangle with ties to the lower triangle, whatever order the tree is walked in, and the summary is accumulated in cell order, so the result is the same on every parallel backend and thread count.

## Core only

The measurement runs in the compiled core and has no numpy twin. Without the core, the Python entry point raises `NotImplementedError` naming it.

## CLI

```sh
meshioplusplus thickness IN [OUT] [--other FILE] [--direction auto|inward|outward]
                         [--max-distance D] [--min T] [--region R] [--region-other R]
                         [--method ray|sphere|cone] [--cone-angle DEG] [--cone-rays N]
                         [--two-sided] [--max-normal-deviation DEG] [--json]
```

`OUT` writes the mesh with both arrays. `--min T` counts the facets below `T` and exits with status 1 when there are any, as [`hausdorff --max`](/hausdorff#cli) does.

## Pipeline

```json
{ "Op": "Thickness", "Other": "", "Direction": "auto", "MaxDistance": 0, "Min": 0, "Region": "", "RegionOther": "", "Require": true, "Method": "ray", "ConeAngle": 15, "ConeRays": 16, "TwoSided": false, "MaxNormalDeviation": 90 }
```

The two arrays flow on with the mesh. With `Min`, a facet below it stops the pipeline unless `"Require": false`; the step reports `NumHit`, `NumMissed`, `NumBelow`, `Min`, `Max`, `Mean`, `P2` and `P95`. See [pipelines](/pipeline).

## Bindings

C: `mio_thickness(mesh, target_or_NULL, &opts, &report)` with `mio_thickness_opts` (72 bytes, initialized by `mio_thickness_opts_init`, `direction` a `mio_thickness_direction`, `method` a `mio_thickness_method`, `two_sided`, `cone_angle`, `cone_rays`, `max_normal_deviation`) and `mio_thickness_report` (144 bytes, with `method`, `p2` and `p95`), and `mio_thickness_warnings`. Fortran `mio_thickness(mesh[, other, direction, max_distance, threshold, region, region_other, report, stat, errmsg, method, cone_angle, cone_rays, two_sided, max_normal_deviation])`, Julia `thickness(m[, other]; direction, max_distance, threshold, region, region_other, method, cone_angle, cone_rays, two_sided, max_normal_deviation)`, R `mio_thickness(mesh, other, ...)` and Octave/MATLAB `mio.thickness(mesh, ...)` (all with a 1-based `min_cell`), WASM `thickness(mesh, {other, direction, maxDistance, threshold, region, regionOther, method, coneAngle, coneRays, twoSided, maxNormalDeviation})` and the MCP tool `thickness` (with an optional `output_path` and a `min_thickness` gate).
