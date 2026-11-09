# Deviation

`deviation(scan, reference)` measures how far every point of a scan lies from the geometry it should have: the signed distance to a reference surface, positive outside it (excess material) and negative inside (missing material), and the statistics an inspection reads from it. Against a fitted primitive (`deviation(scan, fitted, patch=k)`, or a `(type, params)` pair) it reads the analytic distance instead, with no reference mesh and no alignment: the flatness map of a face, the bow of a printed plate.

```python
import meshioplusplus as mp

scan = mp.read("scan.ply")
nominal = mp.read("nominal.stl")

measured, report = mp.deviation(scan, nominal, tolerance=0.05, return_report=True)
report["mean"], report["sigma"], report["percentiles"]["p95"]
report["share_within"]          # the share of points within +-0.05
report["suggested_range"]       # a colour range: p95 of |d|, two significant figures
mp.write("deviation.vtu", measured)  # deviation:distance per point

# Against a fitted cylinder: patch 3 of a mesh fit_primitives wrote.
fitted = mp.fit_primitives(mp.segment(scan))
_, bore = mp.deviation(scan, fitted, patch=3, return_report=True)

# Against an element of a measured mesh's table (roadmap 7.4.2; doc/measure.md).
measured, _ = mp.measure(scan, "inspection.json")
_, face = mp.deviation(measured, measured, element="A", return_report=True)
```

## The distance

Against a surface, each scan point's distance is the distance to the nearest reference triangle, signed by the angle-weighted pseudonormal of the feature it lands on, as [`distance_to_surface`](/sdf) signs it. The reference must be consistently oriented with outward normals.

**The far-side rule.** A scan point next to a thin wall can lie nearer the wall's back face than the face it was scanned from. With `max_normal_angle` below 180 (default 90 degrees), a reference triangle whose unit normal turns away from the scan point's own normal by more than that angle is skipped inside the nearest-triangle search, so the answer is the nearest triangle that faces the same way, still under the (distance squared, triangle id) order that makes the search deterministic. A scan point's normal is its `normals` point data when the scan carries one, else the angle-weighted normal of the scan's own triangles at it; a point with neither (a bare point cloud, an orphan point) is measured without the rule and counted in `num_without_normal`, with a warning.

Against a primitive, the distance is the signed orthogonal distance `fit_primitives` takes its residuals from: along a plane's normal, away from a sphere's centre or a cylinder's, cone's or torus' axis. `max_normal_angle` and `reference_region` do not apply.

## Registering first

`register=True` lays the scan onto the reference before measuring it, so a scan in its scanner's frame needs no manual pre-alignment: [`register_rigid`](./registration.md) with `initial="auto"` (the best of the 24 principal-frame candidates), point-to-plane ICP and adaptive trimming at three median pair distances, restricted to `reference_region`. The scan returned is then the registered one (its vector and tensor point data, `normals` among them, rotated with it), and the report adds `registered`, `matrix` (4 × 4, scan onto reference) and `registration_rms`; the registration's warnings join the report's. A scan of a symmetric part may register in a symmetric pose, which a deviation report cannot tell from the right one. Against a primitive there is nothing to register to, and `register=True` is refused by name.

## What comes back

The scan, bit-identical, with the Float64 point data `deviation:distance` (NaN for a point no accepted triangle reaches; `attach=False` leaves it off), and a report over the measured points:

| Key | Meaning |
|---|---|
| `num_points`, `num_measured` | scan points, and those with a finite distance |
| `min`, `max`, `mean`, `rms`, `sigma` | `sigma` is the population standard deviation about the mean |
| `percentiles` | `p1`, `p5`, `p50`, `p95`, `p99`, by linear interpolation between order statistics |
| `histogram_edges`, `histogram_counts` | `bins` equal bins over `range` (default `+-max |d|`) |
| `num_outside_range` | measured points outside the range, counted in no bin |
| `share_within` | the share within `+-tolerance` (NaN without one) |
| `suggested_range` | the 95th percentile of `|d|`, rounded up to two significant figures: a colour range a deviation map reads well with |
| `num_without_normal` | points measured without the far-side rule |
| `warnings` | what the run could not do |
| `registered`, `matrix`, `registration_rms` | with `register=True` only |

Deterministic: each point's answer is a total-order minimum and the statistics are summed in point order, so the result is the same on every parallel backend and thread count.

## Excess and missing material

`deviation_regions(measured, tolerance)` turns the map into what an inspection report lists: the connected regions (cells joined across shared edges) whose deviation lies beyond `+tolerance` (excess material) or below `-tolerance` (missing material). A cell's deviation is the mean of its corners' `deviation:distance` (`field=` names another point array; corners with NaN are left out). Each region becomes a named cell region, `deviation:excess_<k>` or `deviation:missing_<k>` (`prefix=`), so [`split`](/split) or [`crop`](/crop) can isolate it, and the cell data `deviation:region` holds each cell's region index (-1 within tolerance).

```python
measured = mp.deviation(scan, nominal)
marked, report = mp.deviation_regions(measured, 0.1, return_report=True)
for g in report["regions"]:          # largest |volume| first
    g["name"], g["area"], g["extreme"], g["volume"], g["centroid"]
```

Each record carries `name`, `sign` (+1 excess, -1 missing), `num_cells`, `area`, `extreme` (the largest or the least cell value), the area-weighted `mean`, `volume` (the signed volume estimate, the sum of each cell's deviation times its area), the area-weighted `centroid` and the `bounds` (lo x, y, z, hi x, y, z). Regions are sorted by `|volume|`, ties by the lowest cell they hold; `min_area` drops the small ones and counts them in `num_dropped`. The surface must be 2-D cells (a volume block is refused by name).

## Exaggerated maps

`exaggerate(measured, scale=20)` makes a mesh for a figure or an animation: each point moves along its angle-weighted unit normal by `scale` times its deviation, so a 0.05 mm bow shows. The deviation is clamped to `+-clamp` first (0, the default, takes `suggested_range`, the 95th percentile of `|d|` rounded up to two significant figures), then averaged over the points within `smooth_fraction` of the bounding-box diagonal (default 0.01; weights `1 - (r / R)^2`) so scan noise does not become fur. The measured scan is left untouched; connectivity and every array are copied, and the points come back as Float64. The report gives the `clamp` and smoothing `radius` used and `num_moved`; a point without a normal or a finite deviation stays where it is, with a warning. [`colorize`](/colorize) bakes the same map into colours.

## Core only

The measurement runs in the compiled core and has no numpy twin. Without the core, the Python entry point raises `NotImplementedError` naming it.

## CLI

```sh
meshioplusplus deviation SCAN REFERENCE [OUT] [--patch ID] [--tolerance T] [--min-share S]
                         [--bins N] [--range LO,HI] [--max-normal-angle DEG]
                         [--reference-region R] [--register] [--json]
```

```sh
meshioplusplus deviation-regions IN OUT --tolerance T [--field F] [--prefix P] [--min-area A] [--json]
meshioplusplus exaggerate IN OUT [--field F] [--scale S] [--clamp C] [--smooth-fraction F] [--json]
```

`OUT` writes the scan with `deviation:distance` (with `--register`, the registered scan). With `--patch`, `REFERENCE` is a mesh [`fit-primitives`](/fit_primitives) wrote and the scan is measured against that patch label's primitive. `--min-share S` with a `--tolerance` exits with status 1 when the share within tolerance is below `S`.

## Pipeline

```json
{ "Op": "Deviation", "Reference": "nominal.stl", "Patch": -1, "Tolerance": 0, "MinShare": 0, "Require": true, "Bins": 20, "MaxNormalAngle": 90, "ReferenceRegion": "", "Register": false }
```

```json
{ "Op": "DeviationRegions", "Tolerance": 0.1, "Field": "deviation:distance", "Prefix": "deviation:", "MinArea": 0 }
{ "Op": "Exaggerate", "Field": "deviation:distance", "Scale": 10, "Clamp": 0, "SmoothFraction": 0.01 }
```

`deviation:distance` flows on with the mesh. `Patch` of 0 or more measures against that patch of a fitted reference. With `Tolerance` and `MinShare`, a share below `MinShare` stops the pipeline unless `"Require": false`; the step reports `NumMeasured`, `Min`, `Max`, `Mean`, `Rms`, `Sigma` and `ShareWithin`, and with `Register` also `RegistrationRms`. `DeviationRegions` writes the regions and reports `NumExcess`, `NumMissing` and `NumDropped`; `Exaggerate` reports `NumMoved` and `Clamp`. See [pipelines](/pipeline).

## Bindings

C: `mio_deviation(scan, reference, &opts, &report)` and `mio_deviation_primitive(scan, type, params, &opts, &report)` with `mio_deviation_opts` (88 bytes, initialized by `mio_deviation_opts_init`; `register_first`) and `mio_deviation_report` (176 bytes; `registered`, `registration_rms`), both with reserved tails; `mio_deviation_histogram(counts, edges, capacity)` copies the last call's histogram, `mio_deviation_matrix(matrix)` its registration and `mio_deviation_warnings` its warnings. Fortran `mio_deviation(scan, reference[, tolerance, bins, range_min, range_max, max_normal_angle, reference_region, report, register_first])` and `mio_deviation_primitive(scan, ptype, params[, ...])`, Julia `deviation(scan, reference; ...)` and `deviation(scan, type, params; ...)`, R `mio_deviation(scan, reference, type, params, ...)`, Octave/MATLAB `mio.deviation(scan, 'reference', ref, ...)` or `mio.deviation(scan, 'type', t, 'params', p)`, WASM `deviation(scan, {reference | primitive: {type, params}, tolerance, bins, range, maxNormalAngle, referenceRegion, register})` and the MCP tool `deviation` (with an optional `output_path`, `patch`, `register` and a `min_share` gate). Julia, R and Octave/MATLAB take `register` too. `deviation_regions`: C `mio_deviation_regions(mesh, &opts, &report)` (`mio_deviation_regions_opts`, 64 bytes; `mio_deviation_regions_report`, 72 bytes) with `mio_deviation_region_get(index, &info)` (`mio_deviation_region_info`, 136 bytes), `mio_deviation_region_name` and `mio_deviation_regions_warnings`; Fortran `mio_deviation_regions(mesh, tolerance[, field, prefix, min_area, report])` and `mio_deviation_region(index)` (1-based), Julia `deviation_regions(mesh, tolerance; ...)`, R `mio_deviation_regions(mesh, tolerance, ...)`, Octave/MATLAB `mio.deviation_regions(mesh, tolerance, ...)` (a struct array of regions), WASM `deviationRegions(mesh, {tolerance, field, prefix, minArea})` and the MCP tool `deviation_regions`. `exaggerate`: C `mio_exaggerate(mesh, &opts, &report)` (`mio_exaggerate_opts` and `mio_exaggerate_report`, 64 bytes each, `mio_exaggerate_warnings`), Fortran `mio_exaggerate(mesh[, field, scale, clamp, smooth_fraction, report])`, Julia `exaggerate(mesh; ...)`, R `mio_exaggerate(mesh, ...)`, Octave/MATLAB `mio.exaggerate(mesh, ...)`, WASM `exaggerate(mesh, {field, scale, clamp, smoothFraction})` and the MCP tool `exaggerate`.
