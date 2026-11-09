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
| `region` | none | fit this named `Cell` region as one patch ([below](#a-feature-from-a-seed-or-a-region)) |
| `seeds` | none | grow and fit a feature from the facet nearest each point ([below](#a-feature-from-a-seed-or-a-region)) |
| `feature_angle` | 30 | seeds: edges whose faces turn by more are not crossed |
| `rim_angle` | 15 | seeds: rim faces turning more from the fit are peeled |
| `output_region` | `fit:feature` | seeds: the cell region the fitted cells are written as |
| `used_points` | `"auto"` | `"all"`, or sigma trimming with `"3sigma"`, `"2sigma"`, `"1sigma"` ([below](#robust-fits-sigma-trimming)); `"auto"` is `"3sigma"` with a region or seeds and `"all"` otherwise |

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
| `fit:patch_form_error` | `float64` | the peak-to-peak form error: the largest positive minus the largest negative signed orthogonal distance of the patch's points (flatness, sphericity, cylindricity, under the least-squares reference) |
| `fit:patch_diameter_min`, `fit:patch_diameter_max` | `float64` | sphere and cylinder: twice the smallest and the largest radial distance the points reach |
| `fit:patch_length` | `float64` | cylinder and cone: the points' extent along the axis |
| `fit:patch_arc` | `float64` | cylinder and cone: the arc span in degrees, 360 less the largest angular gap between the points about the axis |
| `fit:patch_extent` | `float64`, `(m, 2)` | plane: the two in-plane principal extents, largest first |
| `fit:patch_sigma` | `float64` | the RMS orthogonal distance of the points the fit used |
| `fit:patch_num_used`, `fit:patch_num_rejected` | `int64` | the points the fit used, and those sigma trimming left out |

`fit:patch_params` holds `x, y, z, ax, ay, az, r1, r2, angle, 0`:

| Type | `x, y, z` | `ax, ay, az` | `r1` | `r2` | `angle` |
|---|---|---|---|---|---|
| plane | a point on the plane (the mean of the patch's points) | the unit normal, pointing the way the cells face | | | |
| sphere | the centre | | the radius | | |
| cylinder | the axis point nearest the patch centroid | the unit axis | the radius | | |
| cone | the apex | the unit axis, pointing into the cone | | | the half-angle, degrees |
| torus | the centre | the unit axis | the major radius | the minor radius | |

Unused entries are 0. A cylinder's or torus's axis is signed so its largest-magnitude component is positive. Patches are fitted in parallel, each serially, so the output is byte-identical on every backend and thread count; the numpy twin follows the same steps and agrees within rounding (it sums in numpy's order).

## A feature from a seed or a region

A metrology report measures one feature: this bore, that face. `region="name"` fits the cells of a named cell region (as [`select`](/select) writes) as one patch. `seeds=[(x, y, z), ...]` finds it: from the facet nearest each point, a region grows across smooth edges (below `feature_angle`) while one primitive of `types` explains it within the tolerance (a single type, `types=["cylinder"]`, or several for "auto"). The growth is the one [`segment`](/segment)'s primitive method runs from its own seeds (`detail/primitive_grow.hpp`, shared): a plane through the seed facet and a curved surface through its two-ring neighbourhood, each grown and trimmed to what its own refit explains, the larger kept. It is tried exactly first (the corners within 2% of the tolerance, as a CAD tessellation lies on its surfaces, which keeps a fillet's tangent faces out), and within the whole tolerance when nothing holds exactly, as a noisy scan wants. The regions of several seeds are united, so one primitive is fitted across a bore a cross hole interrupts. Then the **rim is peeled**: boundary faces whose normal turns from the fitted surface's by more than `rim_angle` go, ring by ring, so a flared bore mouth or a fillet's first strip does not bend the fit. The cells are written as the cell region `output_region` and fitted as patch 0, the usual `fit:` rows; `num_feature_cells` in the report counts them, and a seed nothing grows from leaves a warning. A relative tolerance is then a fraction of the whole surface's diagonal, not of the feature's.

```python
# the bore at (0.6, 0, 0.25): its cells only, its radius to rounding
fit, rep = mp.fit_primitives(part, seeds=[(0.6, 0.0, 0.25)], types=["cylinder"],
                             return_report=True)
fit.field_data["fit:patch_params"][0][6]        # the radius
[r for r in fit.regions if r.name == "fit:feature"]
```

The seed mode runs in the compiled core only (a discrete branch the numpy twin does not repeat, as `segment`'s primitive method); the twin fits a region with `used_points="all"`.

## Robust fits: sigma trimming

Every point weighs the same in a least-squares fit, so a scan's stray points, a burr or a fixture drag it. With `used_points="3sigma"` (the default for a region or seeds), `"2sigma"` or `"1sigma"`, the patch is refitted on the points within k times the RMS distance of the points kept, until that set stops changing (at most 20 rounds); 5% of gross outliers ten sigma off a noisy cylinder move the recovered radius by less than one sigma, where the plain fit moves by more. A plane or sphere fit, when that is the one type asked for, starts from a least median of squares consensus (Rousseeuw 1984): minimal subsets of 3 or 4 points drawn by a fixed stride from the points in Hilbert order, the subset with the least median squared residual kept, its inliers within 2.5 robust standard deviations. There is no random sampling, so the fit is byte-identical on every backend and thread count. The other types start from the least-squares fit, whose outliers the first round of trimming drops. The share of points within tolerance (`fit:patch_fraction`) is then over the points used, and the tolerance stays the whole patch's.

Every fit reports `fit:patch_sigma`, the RMS orthogonal distance of the points it used, and `fit:patch_num_used` and `fit:patch_num_rejected` (0 with `"all"`); the form error and extremes are over the points used.

### Form error and extremes

The last five columns are what a metrology report carries beside the fit. The signed distance is positive on the side the surface's normal points to (along a plane's normal, away from a sphere's centre or a cylinder's, cone's or torus's axis), so the form error is the width of the zone between two offsets of the least-squares surface that holds every point. The arc span counts angles within 1e-9 degrees as one (the rims of a tessellated cylinder share them) and reads 360 when the largest gap is at most twice the median gap, as on a closed tessellation, whose gaps are all alike; a half cylinder reads 180. A value that does not apply to the type, or to a patch of type none, is NaN. The MCP tool reports them per patch as `form_error`, `diameter_min`, `diameter_max`, `length`, `arc_deg` and `extent`.

## CLI

```sh
meshioplusplus fit-primitives IN OUT [--patch-data NAME] [--tolerance T] [--absolute] \
    [--min-inlier-fraction F] [--types plane,cylinder,...] [--max-iterations N] \
    [--region NAME] [--seed X,Y,Z]... [--feature-angle DEG] [--rim-angle DEG] \
    [--output-region NAME] [--used-points auto|all|3sigma|2sigma|1sigma] [--json] [--quiet]
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
const bore = m.fitPrimitives(part, '', 1e-3, true, 0.9, ['cylinder'], 50,
                             { seeds: [[0.6, 0, 0.25]], usedPoints: '3sigma' });
```

The seed, region and trimming options are `mio_fit_primitives_ex(mesh, &opts, &report)` in C (`mio_fit_primitives_ex_opts`, 192 bytes, embedding the `mio_fit_primitives_opts`; `used_points` a `mio_fit_used_points`; `mio_fit_primitives_warnings`), the optional `region`, `seeds` (`3 x n`), `feature_angle`, `rim_angle`, `output_region`, `used_points` and `num_feature_cells` arguments of the Fortran method, the keywords of the Julia function, the R arguments (`seeds` an `n x 3` matrix), the Octave/MATLAB options and the WASM function's last argument.

It is also a settings-pipeline step, `{"Op": "FitPrimitives", "PatchData": "segment:patch", "Tolerance": 0.001, "Relative": true, "MinInlierFraction": 0.9, "Types": ["plane", "cylinder"], "MaxIterations": 50, "Region": "", "Seeds": [0.6, 0, 0.25], "FeatureAngle": 30, "RimAngle": 15, "OutputRegion": "fit:feature", "UsedPoints": "auto"}` (`Seeds` a flat list of x, y, z triples) — see [pipelines](/pipeline) — and the MCP tool `fit_primitives`, whose report lists every patch's type, parameters, residuals, form error and extremes, sigma and points used. `fit_cad` and `refine_patches` write the table without the roadmap 7.1 columns (`FitPrimitivesOptions::mExtendedReport` off), the table they have always carried.
