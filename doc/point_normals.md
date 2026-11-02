# Point-cloud normals

`estimate_normals(points, k=16)` gives every point of a cloud a normal, from the points alone: a scan read from PLY, PCD or XYZ has positions and no surface, and shading it, offsetting it, registering it point-to-plane or reconstructing a surface from it all need one. The result is `point_data["normals"]`, the array the [PCD](/formats/pcd) and [XYZ](/formats/xyz) writers emit as normal columns. For a mesh with faces, [`compute_normals`](./normals.md) takes them from the faces instead. [`reconstruct_surface`](./reconstruct_surface.md) turns the oriented cloud into a triangle surface.

```python
import meshioplusplus as mp

scan = mp.read("part_scan.ply")
out = mp.estimate_normals(scan)                                    # k = 16, oriented consistently
out = mp.estimate_normals(scan, orientation="viewpoint", viewpoint=(0, 0, 2.5))   # towards the scanner
out, report = mp.estimate_normals(scan, k=24, record_variation=True, return_report=True)
report["num_degenerate"], report["num_components"]
mp.write("part_scan.pcd", out)                                     # x y z normal_x normal_y normal_z
```

## The fit

Each point and its `k` nearest neighbours (the exact search of [`proximity_graph`](./proximity_graphs.md), ties to the lower index) are fitted with a plane: the normal is the eigenvector of the neighbourhood's covariance for its least eigenvalue (Hoppe, DeRose, Duchamp, McDonald and Stuetzle, "Surface reconstruction from unorganized points", SIGGRAPH 1992). A larger `k` smooths noise and rounds sharp edges; 10 to 30 suits most scans.

`record_variation=True` attaches `normals:variation`, the surface variation `λ₀ / (λ₀ + λ₁ + λ₂)` of the neighbourhood (Pauly, Gross and Kobbelt 2002): 0 on a plane, up to 1/3 for isotropic scatter. It rises at edges, at high curvature and with noise, so it serves as a curvature and noise indicator.

A point whose neighbourhood spans no plane (its neighbours on one line, or all coincident) gets a `NaN` normal and is counted in `num_degenerate`.

## Orientation

A fitted plane has no side. `orientation` chooses it:

| `orientation` | Side |
| --- | --- |
| `"consistent"` (default) | Hoppe et al.'s propagation: over the symmetrised neighbour graph, weighted `1 − |nᵢ · nⱼ|`, a minimum spanning tree is grown from each connected piece's highest point, and each normal is flipped to agree with its parent's, so the side changes only where the normals turn least. Each piece is then turned as a whole to point away from its centroid: outward for a closed scan. Where that does not decide (a flat patch), the highest point's normal points to +z; failing that, its largest component is positive. |
| `"viewpoint"` | Towards `viewpoint`, the scanner's position: right for a single scan from one side. |
| `"none"` | The largest component positive: deterministic but not consistent across the cloud. |

The propagation needs the pieces to be sampled densely enough that neighbours lie on the same sheet: two sheets closer than the point spacing (a thin plate scanned from both sides) can be bridged, and a piece then gets one side for both. A larger `k` makes this worse; scan such parts per side and orient by `viewpoint`. `num_components` counts the pieces oriented independently.

## Results

On 4 000 random points on a sphere (`k` 12) every normal is within 4° of radial and outward; on a torus sampled with 6 000 random points (`k` 10) every one is within 10° of the exact normal and outward. 500 000 points on a sphere take about 2 s with `orientation="none"` and 5 s consistently oriented on 4 cores: the fits run in parallel, and the spanning tree serially. The result is the same on every backend and thread count.

Data-only: points, cells, regions and every other array pass through; an existing `normals` array is replaced. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent): the spanning tree's ties would have to be broken bit for bit the same way, or whole regions flip. The Python function raises `NotImplementedError` without the compiled core.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `estimate_normals(points, *, k=16, orientation="consistent", viewpoint=(0, 0, 0), record_variation=False, return_report=False)` |
| Pipeline | `{"Op": "EstimateNormals", "K": 16, "Orientation": "consistent", "Viewpoint": [0, 0, 0], "RecordVariation": false}`; counters `NumDegenerate`, `NumComponents` |
| Both CLIs | `meshioplusplus point-normals IN OUT [--k N] [--orientation consistent\|viewpoint\|none] [--viewpoint X,Y,Z] [--variation] [--json]` |
| C | `mio_point_normals_opts_init(&o); o.k = 24; mio_mesh* out = mio_estimate_normals(mesh, &o, &num_degenerate, &num_components);` |
| Fortran | `out = m%estimate_normals(k=24, orientation="viewpoint", viewpoint=[0d0, 0d0, 2.5d0], stat=ierr)` |
| Julia | `estimate_normals(mesh; k=16, orientation="consistent", viewpoint=(0, 0, 0), record_variation=false)` |
| R | `mio_estimate_normals(mesh, k = 16, orientation = "consistent", viewpoint = c(0, 0, 0), record_variation = FALSE)` |
| WASM | `m.estimateNormals(mesh, {k, orientation, viewpoint, recordVariation})` returns `{mesh, numDegenerate, numComponents}` |
| MCP | tool `estimate_normals(input_path, output_path, k, orientation, viewpoint, record_variation)` |
