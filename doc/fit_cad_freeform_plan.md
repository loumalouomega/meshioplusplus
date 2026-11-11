# Freeform fit CAD: plan

The design pass for the next round of [`fit_cad`](./fit_cad.md)'s freeform faces: what the code does today, checked against the source, what each step adds, and how each is measured. Every step is behind an option that is off by default, so the [golden digests](./fit_cad.md#getting-opencascade) and the instruction-count rows `fit_cad_adaptive` and `fit_cad_network` do not move unless a step says why.

## Goals

- **G1, fewer facets.** Lower the faceted area share on every fixture, never raise it, and keep every output a valid closed solid.
- **G2, exact kinematic faces.** Surfaces of linear extrusion, of revolution and helical surfaces where the geometry is one, rational NURBS only where that is exact.
- **G3, accuracy per control point.** A better objective, parameter correction and knot placement.
- **G4, G1 in the fit.** Tangent continuity across network patches enforced by the solve, not only measured.
- **G5, validation.** Two-sided and normal deviation, B-rep validity and a STEP round trip, all in the report.

Not in scope: rational weight optimisation by default (an optional flag at most, step 9), T-splines, LR B-spline entities (STEP has none), GPL dependencies and any machine-learning dependency in the core. The default `freeform_mode` does not change in this work; a change is proposed at the end with the measurements behind it.

## What the code does today

**Orchestration.** `fit_cad` (`operations/fit_cad.cpp`) runs the OpenCASCADE provider (`providers/occt/fit_cad_occt.cpp`, `fc_occt_fit_cad_once`) or the native backend (`fcn_native_fit_cad`); the two share their freeform helpers through `detail/fit_cad_native.hpp`. The native freeform ladder is `fcn_freeform`: the network first in `network` mode, then `freeform_adaptive_fit` per patch, then on a miss an early re-segmentation cut along sharp edges, the topology cut (`patch_split_to_discs`) and the halving (`patch_bisect`) up to `max_freeform_splits`, then facets with the fit's reason; under `auto` what is left goes through the network, and last come `resegment_rounds` rounds of variational re-segmentation. The OpenCASCADE backend runs fillets and threads, the freeform stage (network, adaptive, general and the legacy `GeomAPI_PointsToBSplineSurface` grid), the quartering of demoted adaptive faces, re-segmentation, blends and vertex blends, then builds, checks and demotes until the solid is valid; `fc_build` names `trim_failed`, `sew_failed` and `brep_check_failed`, the loop `kernel_mesh_coarse`, `hausdorff_exceeded` and `other`. Each patch's record becomes a row of `fit_cad:face_table`.

**Acceptance.** `freeform_accept` (`detail/freeform_adaptive`) holds the vertices within half the tolerance, the facet centres within the tolerance and, with `validate_by_surface_sampling`, the Greville points and a grid per knot span inside the trims within the tolerance of the skin (`detail/surface_validate`). A skin that shows noise (`patch_noise_sigma` above a twentieth of the tolerance) whose residuals do not correlate across edges is held to `min(4.5 sigma, 3 tol)` instead. The patch network has its own test (`detail/patch_network`, `accept`).

**Numerics.** `detail/bspline_lsq` holds the clamped cubic basis, the exact thin-plate energy by Gauss-Legendre quadrature, a skyline LDLᵀ least-squares solve and Newton projection. `bspline_fit_patch` (`detail/bspline_fit`) fits the vertices and, at a quarter of their weight, the facet centroids with a fixed thin-plate weight; each iteration refits after projecting the samples (Hoschek's parameter correction, kept on a gain of 10%), else takes one damped tangent-distance step, else inserts a knot at the middle of every span whose samples miss, then refits from scratch. There is no smoothing continuation, no knot removal and no squared-distance minimisation. The patch network shares its edge curves exactly (G0 by construction) and pulls neighbours' cross-boundary derivatives into each other's tangent planes by a penalty in block Gauss-Seidel sweeps, with its boundary curves fixed. `detail/feature_fit`'s `helical_fit` estimates a screw from the patch's normals by Pottmann and Randrup's linear line complex, but refuses a rotation (no advance) and cannot express a translation; its faces are grid approximations, built only by OpenCASCADE. `detail/brep_nurbs` constructs exact rational surfaces of revolution for the four primitives inline, and the native STEP writer (`detail/step_brep`) writes polynomial `B_SPLINE_SURFACE_WITH_KNOTS` faces but no rational, revolution or extrusion surface.

**Options and ABI.** `FitCadOptions` (288 bytes) and `FitCadFaceRecord` (96) are pinned in `tests/cpp/test_abi_layout.cpp`; any new member is a [Tier A](./abi.md) change, and the ABI is 45. The C API's `mio_fit_cad_opts_ex2` has `reserved32` and seven reserved `int64` slots left, `mio_fit_cad_report_ex3` eight. The Python report's legacy keys are a fixed list, and the Python goldens digest every report key, so new keys appear only when their option is on.

**Measurement.** The golden fixtures (box, capped cylinder, turned part, filleted plate, saddle box, wavy band, threaded rod, rounded box, noisy turned part) are pinned in the legacy mode for the native backend and OpenCASCADE 7.6.3 and 7.9.3. The freeform fixtures (`fcfix::freeform_fixtures`: monkey saddle, Gaussian bump, S band, hood, organic torus, sphere with a handle, double torus, Catmull-Clark cages of valence 3, 5 and 6, G1 loaf, with noisy and warped variants) feed `FitCadModes.FacetAreaTable`, which records the faceted area share only. The topology-optimised hook is a data file (`example/data/hook.mdpa`, 18 to 29 minutes a mode with OpenCASCADE).

**Limitations stated today.** G1 next to extraordinary vertices is measured, not held; the B-rep checks still demote fitted faces (the hood's corners); regions of genus one and more are halved geometrically; the noise allowance is capped at three tolerances; helical faces are approximations and OpenCASCADE-only; no surface of revolution or extrusion beyond the five primitives; no knot removal; `extend_surfaces` has not paid off; the native backend has no feature faces.

## Decisions that hold for every step

- **Options.** Few, all off by default: `validation_report` and `validation_acceptance` (step 2), `freeform_fit_schedule` (`"legacy"` or `"continuation"`, step 3), `freeform_objective` (`"default"`, `"pdm"`, `"tdm"` or `"sdm"`, step 4), `kinematic_surfaces` (step 5), `freeform_g1_solve` (`"sweeps"` or `"coupled"`, step 7), and, if they land, `regularity_constraints` (step 8) and `freeform_rational_weights` (step 9). Each reaches Python, both CLIs, the pipeline, the MCP tool, C, Fortran, Julia, R, Octave and WASM.
- **ABI.** One bump, 45 to 46, with the first commit that grows `FitCadOptions` or `FitCadResult`; the later steps of the same release append to the same structs. `FitCadFaceRecord` and the 13 columns of `fit_cad:face_table` stay as they are; per-face measures go to a new `fit_cad:face_deviation` field data array. The C API gains `mio_fit_cad4` and `mio_fit_cad_report_ex4`, which holds every new report field of the plan with a reserved tail; the new options take the reserved slots of `mio_fit_cad_opts_ex2`, whose size does not change.
- **Types.** `fit_cad:type` 11 is a surface of linear extrusion and 12 a surface of revolution; a kinematic helical face is type 8, as a thread face is.
- **Determinism.** New code is serial or uses the chunked parallel reductions; percentiles come from `nth_element` over a fixed sample order. `tools/bench_ops.sh --hash` across backends and thread counts, and a one-against-four-thread digest test per new path, prove it.
- **Dependencies.** None new. PCL's `on_nurbs`, GoTools, QuadWild, Point2CAD and CGAL are references or external benchmarks only.

## Steps

### Step 0: environment and this page

- [x] This page.
- [x] OpenCASCADE 7.6.3 (Ubuntu's packages), a Python environment and the standalone test trees with and without the kernel.
- [x] The goldens pass (native and OpenCASCADE 7.6.3) before and after every step.

### Step 1: the KPI harness

A `fit_cad_kpi` executable beside the tests (it reuses `fit_cad_fixtures.hpp` and `analytic_surfaces.hpp`) runs the golden and freeform fixtures, three new kinematic fixtures (an extruded spline profile, a revolved spline profile, a helical sweep) and, on request, the hook, in every mode and with both backends where built. Per run it records validity, closedness, the volume ratio, the faceted area share, faces by type, control points, demotions, mesh-to-CAD and CAD-to-mesh maxima and 95th percentiles, the largest normal deviation and G1 angle, the STEP write-read-check round trip and the time. `tools/fit_cad_kpi.py` runs it, writes `benchmark/fit_cad_kpi.csv`, regenerates the table in [benchmarks](./benchmarks.md) and, with `--compare`, prints each change's deltas and flags any validity regression or rise in the faceted share. A small subset runs as a test so the harness does not rot.

- [x] Harness (`meshioplusplus_fit_cad_kpi`), script (`tools/fit_cad_kpi.py`), baseline table.

### Step 2: two-sided deviation, normal deviation, round trip

A backend-independent `detail/brep_deviation` samples each face of the solid's NURBS B-rep over its trimmed domain (a fixed grid kept by an even-odd test against its trims, evaluated with rational derivatives), measures each sample's distance and normal angle to the skin, projects the skin's vertices and facet centres onto the B-rep, and reports maxima, 95th and 99th percentiles. The native legacy mode, whose `brep` stays empty, builds one for the measurement only. A STEP written by either backend is read back, checked and its volume compared by a new kernel provider hook where OpenCASCADE is available. `validation_report` turns all of it on; `validation_acceptance` also demotes a face whose sampled deviation leaves its allowance.

- [x] Helper and tests, options, report fields, `fit_cad:face_deviation`, bindings, ABI 46.

Departures: the mesh-to-CAD side measures a skin point whose foot falls outside its face's trims against the face's trim curves on the surface, not against the sample grid: the first version charged a facet centre on a planar side face, next to its trimmed top edge, 2.4 tolerances where the true distance is nil. Acceptance by the validation measures (`validation_acceptance`) applies to the adaptive, network and general freeform faces; the legacy grid faces keep their own test.

### Step 3: the fit schedule

`freeform_fit_schedule="continuation"` changes `bspline_fit_patch`: the thin-plate weight starts high and falls by a decade a stage to `freeform_smoothing`, each stage iterating parameter correction until the parameters or the error stop moving; a knot goes in at the residual-weighted median of the span's bad samples, by exact insertion, and the fit continues from there; once the target is met, knots are removed (Piegl and Tiller, *The NURBS Book*, 2nd ed., algorithm A5.8, across every row) while the face still passes. Measured as poles at equal tolerance.

- [x] The schedule (`freeform_fit_schedule="continuation"`), with tests.

Departures: a knot is removed by refitting the samples without it, not by algorithm A5.8, since the face must follow the samples, not the surface it had; the thin-plate weight falls two decades at a time (one decade at a time spent most of the iteration budget lowering it) and those falls are not charged to the budget; a step still counts as progress from a 10% gain (a 1% threshold left the fit crawling by parameter correction instead of inserting knots). On a smooth sinusoid the continuation ends with more poles than the default schedule (14 × 13 against 11 × 11 at the same target): mid-span knots suit it.

### Step 4: squared-distance minimisation

`freeform_objective="sdm"` warms up with point-distance steps and continues with the squared-distance error term of Wang, Pottmann and Liu (ACM TOG 25(2), 2006, [doi:10.1145/1138450.1138453](https://doi.org/10.1145/1138450.1138453)) extended to surfaces: at each foot point the principal curvatures and directions weigh the tangential residuals by `d / (d - rho)` on the convex side and drop them on the other, the normal residual always counts. `"pdm"` and `"tdm"` force one objective. PCL's `on_nurbs` (BSD) is read as a reference only.

- [x] Principal curvatures helper (`bspline_principal_curvatures`), SDM step (`freeform_objective`), tests.

### Step 5: kinematic surfaces

`kinematic_surfaces` tries each patch no primitive or feature explains, before the freeform ladder, as a surface of linear extrusion, of revolution or a helical surface: the line complex of its normals (Pottmann and Randrup, Computing 60, 1998) as a full 6 × 6 problem, so a translation and a rotation are found as well as a screw, refined by Levenberg-Marquardt; the points projected to the profile plane or meridian half-plane are fitted by a cubic B-spline curve; the face is accepted by the freeform test. OpenCASCADE builds `Geom_SurfaceOfLinearExtrusion` and `Geom_SurfaceOfRevolution`; the native writer writes `SURFACE_OF_LINEAR_EXTRUSION` and `SURFACE_OF_REVOLUTION` (or, if a reader disagrees, the exact rational B-spline of step 6) and gains helical faces as B-splines swept from the exact profile.

- [ ] Kinematic fit with tests, both backends, report counts, KPI.

### Step 6: rational output only where exact

The revolution construction inside `brep_nurbs_face` becomes `brep_nurbs_revolve`, with `brep_nurbs_extrude` beside it, so the NURBS B-rep (and Kratos CAD JSON) of a kinematic face is exact; the four primitives' output stays bit for bit. Tests check each rational surface against its analytic or kinematic original to `1e-9`.

- [ ] Refactor, extrusion, tests.

### Step 7: a coupled G1 network

`freeform_g1_solve="coupled"` fits each network region as one sparse system: the interior poles of every patch and of every shared edge curve are unknowns (so both patches still take the same boundary row and G0 holds), the tangent-plane conditions on shared edges are penalty rows whose weight rises in steps, solved by the deterministic conjugate gradients and alternated with parameter correction. A region whose coupled fit passes fewer patches than the sweeps keeps the sweeps. The largest G1 angle is reported away from and at extraordinary vertices separately: tensor-product cubics cannot hold G1 exactly at a vertex of valence other than four, and the page says so rather than claiming it.

- [ ] Coupled solve, tests, KPI on the hood box and the cube cage.

### Steps 8 and 9, if time allows

Regularity constraints (tangency, offset axes, angles snapped to 0, 45 or 90 degrees, each re-checked and reported) extend [`constrain_primitives`](./constrain_primitives.md) behind `regularity_constraints`; rational weights are optimised by bounded Gauss-Newton behind `freeform_rational_weights`, kept per face only where the deviation falls at equal poles. Either is merged only if the KPI table shows a gain, and is otherwise recorded under [measured and parked](./fit_cad.md#measured-and-parked).

## Measurements

The KPI table of each stage is in [benchmarks](./benchmarks.md#fit-cad-freeform-kpis) once step 1 lands; each stage's deltas against the previous one are recorded here.
