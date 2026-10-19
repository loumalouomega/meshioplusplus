# Fit primitives

`fit_primitives(mesh)` fits a plane, sphere, cylinder, cone or torus to each patch of a surface and reports how well each fits. A recovered radius, axis, apex or plane normal is a measurement — the radius of a drilled hole, the axis of a shaft — whether or not the surface ever goes to CAD. Patches come from [`segment`](/segment) by default; see [skin to CAD](/cad_fit) for the whole workflow.

```python
import meshioplusplus as mp

seg = mp.segment(mp.read("bracket.stl"))
fit, report = mp.fit_primitives(seg, return_report=True)
report["num_plane"], report["num_cylinder"], report["num_unfitted"]

fd = fit.field_data
fd["fit:patch_type"]       # per patch: 0 none, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus
fd["fit:patch_params"]     # per patch: (x, y, z, ax, ay, az, r1, r2, angle, 0)
fd["fit:patch_rms"]        # per patch: RMS orthogonal distance of its points
fit.cell_data["fit:type"]  # per cell: its patch's type

# a whole surface as one patch, an absolute tolerance, only round types
mp.fit_primitives(part, patch_data="", tolerance=0.01, relative=False,
                  types=["cylinder", "cone"])
```

## Fitting and accepting

Each patch is fitted to the corner points of its cells, every type deterministically (no random sampling):

| Type | Initial guess | Refined |
|---|---|---|
| plane | least-squares plane through the centroid (smallest eigenvector of the covariance) | — (already optimal) |
| sphere | algebraic fit (a linear least-squares problem) | yes |
| cylinder | axis along the direction the cell normals vary least; axis point where the normal lines meet, in the least-squares sense; radius the mean distance | yes |
| cone | the same axis; apex and half-angle from a line fit of radius against height | yes |
| torus | the minor radius scanned: cell centroids moved back along their normals by it must be cocircular; the circle gives the centre, axis and major radius | yes |

Refinement is Levenberg-Marquardt on the orthogonal distance of every point, with a forward-difference Jacobian and at most `max_iterations` steps. A type is **accepted** when its RMS distance is at most the tolerance and at least `min_inlier_fraction` of the points lie within it. The patch takes the first accepted type in the order plane, sphere, cylinder, cone, torus — the simplest that explains it, since a plane is also a cylinder of infinite radius. A sphere is the one exception: a short bore or a narrow band of a cone also fits a sphere through its rims within the tolerance, so when a sphere is accepted the later types are fitted too, and a cylinder or cone whose facet centres lie no further from it than from the sphere (a torus: half as far) is taken instead (since v17.9.0). When no type is accepted the patch's type is `none` (0), its parameters are NaN and its residuals are the best attempt's: a fit that fails the tolerance is reported as failing, never rounded to a neat primitive.

| Option | Default | Meaning |
|---|---|---|
| `patch_data` | `segment:patch` | the integer cell-data label; cells below 0 belong to no patch; `""` fits the whole surface as one patch |
| `tolerance` | `1e-3` | the acceptance tolerance on the orthogonal distance |
| `relative` | `True` | `tolerance` is a fraction of each patch's bounding-box diagonal |
| `min_inlier_fraction` | `0.9` | the share of points that must be within tolerance |
| `types` | all five | the candidates, tried in the fixed order above |
| `max_iterations` | 50 | the Levenberg-Marquardt cap per fit |

Points on a tessellated analytic surface lie on it exactly, so a tessellated cylinder recovers its radius to rounding (about 1e-15) however coarse the tessellation; noise raises the residual, not the type. A patch too small for a type (fewer points than the type has parameters) skips it.

## Output

The input, unchanged, plus cell data `fit:type` (`int64`: -1 no patch, 0 none, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus), `fit:rms` and `fit:max` (`float64`, the patch's; NaN with no patch), and one field-data row per patch label, ascending:

| Field data | Type | Content |
|---|---|---|
| `fit:patch_id` | `int64` | the patch label |
| `fit:patch_type` | `int64` | as `fit:type` |
| `fit:patch_params` | `float64`, `(m, 10)` | see below |
| `fit:patch_rms`, `fit:patch_max` | `float64` | the accepted fit's residuals, or the best attempt's |
| `fit:patch_fraction` | `float64` | the share of points within tolerance |
| `fit:patch_tolerance` | `float64` | the absolute tolerance applied |

`fit:patch_params` holds `x, y, z, ax, ay, az, r1, r2, angle, 0`:

| Type | `x, y, z` | `ax, ay, az` | `r1` | `r2` | `angle` |
|---|---|---|---|---|---|
| plane | a point on the plane (the mean of the patch's points) | the unit normal, pointing the way the cells face | | | |
| sphere | the centre | | the radius | | |
| cylinder | the axis point nearest the patch centroid | the unit axis | the radius | | |
| cone | the apex | the unit axis, pointing into the cone | | | the half-angle, degrees |
| torus | the centre | the unit axis | the major radius | the minor radius | |

Unused entries are 0. A cylinder's or torus's axis is signed so its largest-magnitude component is positive. Patches are fitted in parallel, each serially, so the output is byte-identical on every backend and thread count; the numpy twin follows the same steps and agrees within rounding (it sums in numpy's order).

## CLI

```sh
meshioplusplus fit-primitives IN OUT [--patch-data NAME] [--tolerance T] [--absolute] \
    [--min-inlier-fraction F] [--types plane,cylinder,...] [--max-iterations N] \
    [--json] [--quiet]
```

## Other languages

```c
mio_fit_primitives_opts opts;
mio_fit_primitives_opts_init(&opts);
opts.types = (1u << MIO_FIT_PLANE) | (1u << MIO_FIT_CYLINDER);
mio_fit_primitives_report report;
mio_mesh* fit = mio_fit_primitives(segmented, &opts, &report);
```

```fortran
fit = seg%fit_primitives(tolerance=1.0e-3_real64, num_plane=nplane)
```

```julia
r = fit_primitives(seg; types=("plane", "cylinder"))
r.num_cylinder
```

```r
r <- mio_fit_primitives(seg, types = c("plane", "cylinder"))
```

```js
const r = m.fitPrimitives(seg, 'segment:patch', 1e-3, true, 0.9, ['plane', 'cylinder']);
r.mesh.field_data['fit:patch_params'];
```

It is also a settings-pipeline step, `{"Op": "FitPrimitives", "PatchData": "segment:patch", "Tolerance": 0.001, "Relative": true, "MinInlierFraction": 0.9, "Types": ["plane", "cylinder"], "MaxIterations": 50}` — see [pipelines](/pipeline) — and the MCP tool `fit_primitives`, whose report lists every patch's type, parameters and residuals.
