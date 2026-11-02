# Rigid registration

`register_rigid(source, target)` finds the rotation and translation that best lay one mesh or point cloud onto another, by iterative closest points (Besl and McKay, "A method for registration of 3-D shapes", IEEE PAMI 14(2), 1992), and returns the 4 × 4 matrix [`transform`](./transform.md) applies. Uses: aligning a scan to its design before measuring the deviation, aligning two runs' meshes before comparing them, and stitching scans.

```python
import meshioplusplus as mp

scan = mp.read("part_scan.ply")
design = mp.read("part.stl")
matrix, report = mp.register_rigid(scan, design, method="point_to_plane")
aligned = mp.transform(scan, matrix=matrix)
report["rms"], report["iterations"], report["converged"]

deviation = mp.hausdorff_distance(aligned, design)   # after alignment, what is left is shape
```

## How it works

From `initial` (the identity unless given), each iteration:

1. **Pairs** every source point, under the current transform, with its closest point on the target: on the target's surface when it has surface cells (triangles, quads, polygons; the exact closest-point query of [`sdf`](./sdf.md)), else with its nearest target point (a point cloud). With `max_distance` > 0, pairs farther apart are dropped: outliers, and parts the target does not have.
2. **Solves for the best rigid motion** for those pairs and applies it:
   - `"point_to_point"` minimises the summed squared distances between the pairs, in closed form by Horn's unit quaternions ("Closed-form solution of absolute orientation using unit quaternions", JOSA A 4(4), 1987);
   - `"point_to_plane"` minimises the summed squared distances to the target's tangent planes at the paired points (Chen and Medioni 1992), linearised in the rotation (Low, "Linear least-squares optimization for point-to-plane ICP surface registration", 2004) and applied as an exact rotation. A direction the planes do not constrain (sliding along a plane or a cylinder) is left still. The plane is the target cell's normal, or for a cloud its `normals` point data ([`estimate_normals`](./point_normals.md)); a cloud without one is refused by name.

It stops when a step rotates by less than `tolerance` radians and moves the source's centre by less than `tolerance` times the source's bounding-box diagonal, or after `max_iterations` (`converged` false, with a warning; the last iterate is returned).

**Point-to-plane against a surface converges in a handful of iterations; point-to-point against a surface only linearly**, because its closest points keep sliding along the faces: from the same start, onto the same ellipsoid, point-to-plane reaches 1e-9 in about 5 iterations where point-to-point needs hundreds. Against a cloud of the same points, both are exact within a few iterations. Prefer point-to-plane whenever the target has faces or normals.

## What it will not do

ICP finds the **nearest** local optimum. The two must start roughly aligned (a few degrees and a fraction of the size apart): a larger offset can settle in a wrong pose, and a symmetric part (a cylinder, a sphere, a plate) in any of its symmetric poses. Give `initial` a coarse alignment (a manual transform, or the matrix from aligning a few picked landmarks) when they do not. It does not scale or shear: the matrix is a rotation and a translation.

## Results

The roadmap probe, a rotated (0.12 rad) and shifted copy of a lumpy ellipsoid, registers back within 1e-9 by point-to-plane onto the surface or onto the cloud, and by point-to-point onto the cloud (`tests/python/test_registration.py`). Only points are read: the source's cells and data play no part, and `transform` moves them afterwards (`rotate_vector_data=True` turns vector data with the points). The correspondence search runs in parallel and the solves serially; the result is the same on every backend and thread count. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent): the stopping test is rounded, so a second implementation could stop a step early or late.

## Report

`initial_rms` and `rms` (root-mean-square pair distance under `initial` and at the end), `iterations`, `num_pairs` (pairs kept at the end) and `converged`.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `register_rigid(source, target, *, method="point_to_point", max_iterations=50, tolerance=1e-10, max_distance=0.0, initial=None) -> (matrix, report)` |
| Pipeline | `{"Op": "RegisterRigid", "Target": "design.stl", "Method": "point_to_plane", "MaxIterations": 50, "Tolerance": 1e-10, "MaxDistance": 0, "RotateData": false}` moves the mesh onto `Target`; counters `InitialRms`, `Rms`, `Iterations`, `Converged` |
| Both CLIs | `meshioplusplus register SOURCE TARGET [--output OUT] [--method point_to_point\|point_to_plane] [--max-iterations N] [--tolerance T] [--max-distance D] [--rotate-data] [--json]` |
| C | `mio_registration_opts_init(&o); o.method = MIO_REGISTER_POINT_TO_PLANE; mio_register_rigid(source, target, &o, &report);` then `report.matrix` |
| Fortran | `call scan%register_rigid(design, matrix, method="point_to_plane", rms=rms, stat=ierr)` |
| Julia | `r = register_rigid(scan, design; method="point_to_plane")`; `r.matrix` (4 × 4), `r.rms`, ... |
| R | `r <- mio_register_rigid(scan, design, method = "point_to_plane")`; `r$matrix` (4 × 4), `r$rms`, ... |
| WASM | `m.registerRigid(source, target, {method, maxIterations, tolerance, maxDistance, initial})` returns `{matrix, initialRms, rms, iterations, numPairs, converged}` (`matrix` 16 numbers, row-major) |
| MCP | tool `register_rigid(source_path, target_path, output_path, method, max_iterations, tolerance, max_distance, rotate_data)` |
