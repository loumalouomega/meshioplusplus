# Geodesic distance

`geodesic_distance(mesh, sources)` measures, for every point of a triangle surface, the distance **along the surface** to the nearest source point, and returns the mesh with `point_data["geodesic:distance"]`. It is the heat method of Crane, Weischedel and Wardetzky ("Geodesics in Heat", ACM Transactions on Graphics 32(5), 2013) over the surface's own finite-element operators. Uses: distance-from-a-feature fields as inputs to a physics-ML model, seeding (points spaced by surface distance rather than straight-line distance), and the distance from an inlet, a weld line or a contact patch.

```python
import meshioplusplus as mp

surface = mp.read("bracket_skin.stl")
out = mp.geodesic_distance(surface, [0, 117])          # point ids
out = mp.geodesic_distance(surface, "weld_line")       # the points of a region
d = out.point_data["geodesic:distance"]

out, report = mp.geodesic_distance(surface, 0, time_factor=4.0, return_report=True)
report["max_distance"], report["num_unreached"], report["converged"]
```

## How it works

Over the triangles' cotangent stiffness `K` and lumped (one-third area) mass `M`:

1. **Heat.** Solve `(M + t K) u = δ`, one backward-Euler step of the heat equation from `δ` (1 at the sources), for a short time `t = time_factor · h²` with `h` the mean edge length.
2. **Direction.** On every triangle, `X = −∇u / |∇u|`: the unit vector pointing away from the sources. Only the direction of the heat flow is kept, which is what makes the method insensitive to how fast the heat decays.
3. **Distance.** Solve the Poisson problem `K φ = −∇·X`, with one point of each connected piece of the surface held at 0, then shift `φ` so its least value over each piece is 0.

The heat solve is Gauss–Seidel, not conjugate gradients: far from the sources the heat falls to 10⁻³⁰ and below, which a residual-based stop resolves to noise, while sweeps that stop only once no value changes by more than the tolerance relative to itself resolve each value however small (the property of a direct factorisation that Crane et al. rely on). Gauss–Seidel converges for any symmetric positive definite system, so the negative weights of obtuse triangles do not stop it; it takes a few hundred sweeps on the meshes below (a sweep per edge-hop across the surface, roughly). It runs serially in point order, and the Poisson solve is Jacobi-preconditioned conjugate gradients, so the result is the same on every parallel backend and thread count.

## Accuracy

The heat method is approximate: it measures distance on a slightly smoothed surface, at the scale `√t ≈ √time_factor · h`. On the unit sphere (an icosahedral sphere projected onto it), from the north pole:

| Points | `h` | pole → equator error | mean error, `time_factor` 1 | mean error, `time_factor` 4 |
| --- | --- | --- | --- | --- |
| 642 | 0.151 | 0.021 | 0.021 | 0.017 |
| 2 562 | 0.075 | 0.015 | 0.014 | 0.009 |
| 10 242 | 0.038 | 0.012 | 0.011 | 0.005 |
| 40 962 | 0.019 | 0.011 | 0.010 | 0.004 |

At Crane et al.'s `time_factor = 1` most of the error is a scale, about 0.6% short, that refinement does not remove; a larger factor reduces it on smooth surfaces and smooths the distance near sharp features and the source. On the unit square, triangulated 40 × 40, the distance from the centre is within 0.006 of the Euclidean one on average (0.022 at worst, near the corners). A boundary carries the natural (Neumann) condition, which bends the distance near it slightly; Crane et al.'s average of Neumann and Dirichlet solutions is not applied.

## Scope

- **The surface is the mesh's `triangle` blocks.** A volume cell is refused naming `extract_surface`, a quad or polygon naming `convert_cells(mode="simplexify")`, a higher-order surface cell naming `linearize`. Line and vertex blocks ride along untouched.
- **Sources** are point ids (0-based in Python, C, WASM and the CLIs; 1-based in Fortran, Julia and R), and/or the points of a region: a point region's entries, a cell region's cells' points, a side region's facets' points. A Python string naming a `point_set` works too. A source outside the mesh, or in no triangle, is an error, and so is no source at all.
- A point in no triangle, and every point of a connected piece with no source, gets `NaN` (`num_unreached`). A zero-area triangle adds nothing (`num_degenerate`).
- Everything else passes through: points, connectivity, every other data array, regions.

There is no numpy twin (the [`sobolev_deform`](./sobolev_deform.md) precedent): both solves stop on a rounded test, so a second implementation could stop a sweep early or late and disagree. The Python function raises `NotImplementedError` without the compiled core.

## Report

`num_sources` (distinct source points), `num_unreached`, `num_degenerate`, `max_distance` (the largest finite distance), `time_step` (`t`), `heat_iterations` and `poisson_iterations`, and `converged` (false when a solve hit its iteration cap; the last iterate is still used, with a warning).

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `geodesic_distance(mesh, sources=None, *, source_region=None, time_factor=1.0, tolerance=1e-10, max_iterations=0, return_report=False)` |
| Pipeline | `{"Op": "GeodesicDistance", "Sources": [0], "SourceRegion": "inlet", "TimeFactor": 1, "Tolerance": 1e-10, "MaxIterations": 0}`; counters `NumSources`, `NumUnreached`, `MaxDistance`, `Converged` |
| Both CLIs | `meshioplusplus geodesic IN OUT --source 0,117 [--source-region NAME] [--time-factor F] [--tolerance T] [--max-iterations N] [--json]` |
| C | `mio_geodesic_opts_init(&o); o.sources = ids; o.num_sources = n; mio_mesh* out = mio_geodesic_distance(mesh, &o, &report);` |
| Fortran | `out = m%geodesic_distance(sources=[1_int64], max_distance=d, converged=ok, stat=ierr)` |
| Julia | `r = geodesic_distance(mesh, [1]; time_factor=1.0)`; `r.mesh`, `r.max_distance`, ... |
| R | `r <- mio_geodesic_distance(mesh, 1)`; `r$mesh`, `r$max_distance`, ... |
| WASM | `m.geodesicDistance(mesh, [0], {sourceRegion, timeFactor, tolerance, maxIterations})` |
| MCP | tool `geodesic_distance(input_path, output_path, sources, source_region, time_factor, ...)` |
