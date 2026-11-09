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

## A global start

`initial="auto"` needs no start: a scanner's frame is arbitrary. The principal frames of the two inputs (the area-weighted second moments of their surface cells, exact per triangle, or of their points for a cloud; axes by decreasing spread, right-handed) are laid onto each other by each of the 24 proper rotations that map one frame's signed axes onto the other's, since the axes' signs, and their order where two spreads are close, are ambiguous. Point-to-point ICP runs from each candidate on a strided subsample of at most 1000 source points (30 iterations, tolerance 1e-6), each is scored by its trimmed RMS (over the closest 90% of its pairs), and the lowest, ties to the lowest index, is refined with the options given. The report lists each candidate's RMS (`candidate_rms`) and the one refined (`chosen_candidate`).

```python
matrix, report = mp.register_rigid(scan, design, method="point_to_plane", initial="auto", symmetric_check=True)
report["warnings"]   # a symmetric part: another candidate fits as well at another pose
```

`symmetric_check` warns when a candidate whose pose is more than 5 degrees from the chosen one scores within `symmetric_tolerance` (5%) of it: a part symmetric under a rotation (a box, a cylinder) fits in several poses, and the surfaces alone cannot tell which is meant. The probe, a copy of a non-symmetric L bracket rotated by 150 degrees about a skew axis, registers back to 1e-9 with `initial="auto"` and does not with the identity start (`tests/cpp/test_registration.cpp`).

## Regions and adaptive trimming

`source_region` restricts the pairs to the points of a cell region's surface cells of the source, and `target_region` the matches to a cell region of the target: the part's true surface, without the fixture, support or turntable the scan picked up. `adaptive_trim=k` drops, at each iteration, the pairs farther apart than `k` times the median pair distance (after Chetverikov et al., "The trimmed iterative closest point algorithm", ICPR 2002), and `max_distance="auto"` in Python is `adaptive_trim=3`; with `max_distance` set as well, both apply. `num_rejected` counts the source points left out at the end. The probe: a scan with a fixture plate the nominal lacks registers within 1e-6 with the part's region named and is pulled off by the plate without it; trimming alone rejects the plate's far pairs and lands closer.

## Landmarks

`register_landmarks(source_points, target_points, weights=None)` is the weighted Horn solve on three or more picked pairs, in closed form: the start a viewer or a script gives `register_rigid`. It returns the matrix and each pair's residual. Fewer than three pairs, fewer than three of positive weight, or collinear source points (the rotation about their line is then free) are refused by name; three exact pairs recover a known transform to 1e-12.

## Auto-align from the scan alone

`auto_align(mesh)` poses a part without a nominal: on its largest flat face, or a cylindrical part on the end of its axis, so parts scanned in arbitrary frames come out in one.

- **The histogram.** Each facet's unit normal, weighted by its area, falls in the bin of its nearest vertex of a subdivided icosahedron on the Gauss sphere (`subdivisions`, 642 bins at the default 3). A bin scores the area of every bin within `angle_tolerance` (10 degrees) of it; the best, ties to the lowest bin, is the dominant direction, refined to the area-weighted mean of the facet normals within the tolerance, and the best bin perpendicular to it (within the tolerance) the second direction.
- **An opening.** The facets' summed vector area is the vector area of the boundary loops (Stokes), so an open scan's opening joins the histogram as a pseudo-facet of that area, facing away from the sum, when it exceeds 0.1% of the total area: a part scanned on a turntable lies on its opening. This needs a consistent winding, as the faces do.
- **Two poses.** `face` agrees with the area whose normals lie within the tolerance of one of the frame's three axes, `axis` with the area whose normals are perpendicular to, or along, the least-eigenvalue direction of the area-weighted normal scatter (a cylinder's axis). `axis` is taken when it beats `face` by more than 0.05. `face` sends the dominant direction to -Z (the face down) and the second to +X; `axis` sends the axis to +Z, the end with the larger area facing along it down, and the dominant bin around it to +X.

The result is translated so its bounding box is centred on the Z axis with its lowest point at z = 0. The report gives the matrix, the `mode`, the input directions sent to -Z or +Z (`primary`) and +X (`secondary`), the agreeing `share` of the chosen pose and both poses' (`face_share`, `axis_share`, to reject a pose) and the `opening_area`. A 3 × 2 × 1 block rotated arbitrarily comes back on a 3 × 2 face, axis-aligned, with a share of 1; an open tube stands on its axis.

## The symmetry plane

`find_symmetry_plane(mesh)` finds the plane a part is most nearly symmetric about, for the mirror step of reverse engineering (model half, mirror the rest; `mirror(plane="auto")`) and for a datum. A symmetry plane's normal is an eigenvector of the second moments, so the candidates are the three principal planes through the centroid, and an optional seed plane. Each candidate's mirror image (a strided subsample of at most 500 points) is registered to the original by ICP (point-to-plane onto a surface, point-to-point onto a cloud, adaptive trimming at three medians); the motion composed with the reflection is an improper isometry, and its plane is read off the eigenvector of its linear part's -1 eigenvalue (the symmetric part's least eigenvector) and its translation. The candidate of the lowest trimmed RMS is refined three times on up to `max_samples` points, each from the plane the last gave, so the plane is its own fixed point. The plane's exact mirror image is then measured against the surface: `rms`, and the `share` of samples within `tolerance` (0.005 of the diagonal by default). A plane whose share is below `min_share` (0.9) is refused by name, so an asymmetric part gets no plane. The normal is signed so its largest-magnitude component is positive, and the reported origin is the centroid projected onto the plane. The probe, a bracket symmetric about a tilted plane with Gaussian noise of 2e-4, finds the plane within 1e-4 radians and 1e-4 of the diagonal.

## What it will not do

Without `initial="auto"`, ICP finds the **nearest** local optimum. The two must start roughly aligned (a few degrees and a fraction of the size apart): a larger offset can settle in a wrong pose, and a symmetric part (a cylinder, a sphere, a plate) in any of its symmetric poses. Give `initial` a coarse alignment (a manual transform, `register_landmarks` on a few picked pairs, or `"auto"`) when they do not. The automatic start cannot decide the pose of a symmetric part (`symmetric_check` says so). It does not scale or shear: the matrix is a rotation and a translation.

## Results

The roadmap probe, a rotated (0.12 rad) and shifted copy of a lumpy ellipsoid, registers back within 1e-9 by point-to-plane onto the surface or onto the cloud, and by point-to-point onto the cloud (`tests/python/test_registration.py`). Only points are read: the source's cells and data play no part, and `transform` moves them afterwards (`rotate_vector_data=True` turns vector data with the points). The correspondence search runs in parallel and the solves serially; the result is the same on every backend and thread count. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent): the stopping test is rounded, so a second implementation could stop a step early or late.

## Report

`initial_rms` and `rms` (root-mean-square pair distance under `initial`, or the chosen candidate's coarse fit, and at the end), `iterations`, `num_pairs` (pairs kept at the end), `converged`, `num_rejected` (source points left out at the end), `candidate_rms` and `chosen_candidate` (with `initial="auto"`; -1 otherwise) and `warnings`.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `register_rigid(source, target, *, method="point_to_point", max_iterations=50, tolerance=1e-10, max_distance=0.0, initial=None, symmetric_check=False, symmetric_tolerance=0.05, source_region=None, target_region=None, adaptive_trim=0.0) -> (matrix, report)`; `register_landmarks(source_points, target_points, weights=None) -> (matrix, report)`; `find_symmetry_plane(mesh, *, seed=None, tolerance=0.0, min_share=0.9, max_samples=4000, max_iterations=50) -> report`; `auto_align(mesh, *, angle_tolerance=10.0, subdivisions=3, region=None, rotate_vector_data=False, return_report=False)`; `mirror(mesh, plane="auto")` |
| Pipeline | `{"Op": "RegisterRigid", "Target": "design.stl", "Method": "point_to_plane", "MaxIterations": 50, "Tolerance": 1e-10, "MaxDistance": 0, "RotateData": false, "Initial": "identity", "SymmetricCheck": false, "SourceRegion": "", "TargetRegion": "", "AdaptiveTrim": 0}` moves the mesh onto `Target`; counters `NumRejected`, `InitialRms`, `Rms`, `Iterations`, `Converged`. `{"Op": "RegisterLandmarks", "SourcePoints": [x, y, z, ...], "TargetPoints": [...], "Weights": [...], "RotateData": false}` (counter `Rms`), `{"Op": "AutoAlign", "AngleTolerance": 10, "Subdivisions": 3, "Region": "", "RotateData": false}` (counters `Share`, `AxisMode`) and `{"Op": "Mirror", "Plane": "auto"}` (counter `SymmetryShare`) |
| Both CLIs | `meshioplusplus register SOURCE TARGET [--output OUT] [--method point_to_point\|point_to_plane] [--max-iterations N] [--tolerance T] [--max-distance D] [--initial identity\|auto] [--symmetric-check] [--source-region R] [--target-region R] [--adaptive-trim K] [--rotate-data] [--json]`; `register-landmarks "x,y,z;..." "x,y,z;..." [--weights w,...] [--mesh IN --output OUT] [--json]`; `auto-align IN OUT [--angle-tolerance DEG] [--subdivisions N] [--region R] [--rotate-data] [--json]`; `symmetry-plane IN [--seed-origin x,y,z --seed-normal x,y,z] [--tolerance T] [--min-share S] [--max-samples N] [--json]`; `mirror IN OUT --plane auto` |
| C | `mio_registration_opts_init(&o); o.method = MIO_REGISTER_POINT_TO_PLANE; mio_register_rigid(source, target, &o, &report);` then `report.matrix`; `mio_registration_ex_opts` and `mio_register_rigid_ex` for `initial_mode = MIO_REGISTER_INITIAL_AUTO`, the regions and `adaptive_trim` (`mio_register_rigid_candidates`, `mio_register_rigid_warnings`); `mio_register_landmarks(src, dst, weights, n, matrix, residuals, &rms)`; `mio_find_symmetry_plane(mesh, &opts, &report)`; `mio_auto_align(mesh, &opts, &report)` |
| Fortran | `call scan%register_rigid(design, matrix, method="point_to_plane", initial="auto", rms=rms, stat=ierr)`; `call mio_register_landmarks(src, dst, matrix)`; `call mio_find_symmetry_plane(mesh, report)`; `aligned = mio_auto_align(mesh, report=r)` |
| Julia | `r = register_rigid(scan, design; method="point_to_plane", initial="auto")`; `r.matrix` (4 × 4), `r.rms`, ...; `register_landmarks(src, dst)`, `find_symmetry_plane(m)`, `auto_align(m)` |
| R | `r <- mio_register_rigid(scan, design, method = "point_to_plane", initial = "auto")`; `r$matrix` (4 × 4), `r$rms`, ...; `mio_register_landmarks()`, `mio_find_symmetry_plane()`, `mio_auto_align()` |
| Octave / MATLAB | `mio.register_rigid(scan, design, 'initial', 'auto')`, `mio.register_landmarks`, `mio.find_symmetry_plane`, `mio.auto_align` |
| WASM | `m.registerRigid(source, target, {method, maxIterations, tolerance, maxDistance, initial, symmetricCheck, sourceRegion, targetRegion, adaptiveTrim})` returns `{matrix, initialRms, rms, iterations, numPairs, converged, numRejected, candidateRms, chosenCandidate, warnings}` (`matrix` 16 numbers, row-major; `initial` 16 numbers or `'auto'`); `m.registerLandmarks(src, dst, weights)`, `m.findSymmetryPlane(mesh, options)`, `m.autoAlign(mesh, options)` |
| MCP | tools `register_rigid(source_path, target_path, output_path, method, max_iterations, tolerance, max_distance, rotate_data, initial, symmetric_check, source_region, target_region, adaptive_trim)`, `register_landmarks`, `find_symmetry_plane`, `auto_align`, and `mirror(plane="auto")` |
