# Refine patches

`refine_patches(mesh)` refines the patches of a segmented surface until one surface explains each: it moves patch boundaries, splits patches holding several surfaces and merges patches holding one, so [`fit_primitives`](/fit_primitives) and [`fit_cad`](/fit_cad) find more real CAD faces. It sits between [`segment`](/segment) and the fit in the [skin-to-CAD](/cad_fit) chain (since v17.11.0).

[`segment`](/segment) draws its boundaries once. Segmented by angles, a fillet tangent to a wall shares the wall's patch and no primitive fits the pair; a cylinder cut at a smooth edge stays in two pieces; a freeform patch with a cylinder inside it stays freeform. `refine_patches` starts from those `segment:patch` labels and fixes them.

```python
import meshioplusplus as mp

skin = mp.read("bracket.stl")
seg = mp.segment(skin)                                      # by 30-degree creases

ref = mp.refine_patches(seg)                                # tangent junctions separated
ref, report = mp.refine_patches(seg, return_report=True)
report["initial_explained_area_fraction"], report["explained_area_fraction"]
report["num_splits"], report["num_merges"], report["stop"]

ref = mp.refine_patches(seg, tolerance=2e-3)                # a noisy scan
ref = mp.refine_patches(seg, allow_split=False)             # move and merge only
ref = mp.refine_patches(seg, ransac_seed=7)                 # another RANSAC sequence

solid = mp.fit_cad(ref, "bracket.step")                     # reads the refined fit table
solid = mp.fit_cad(skin, "bracket.step", refine_segmentation=True)  # the same, in one call
```

## How patches are refined

The method is a hybrid variational shape approximation: the Lloyd iteration over surface proxies of Cohen-Steiner, Alliez and Desbrun ([*Variational shape approximation*, SIGGRAPH 2004](https://doi.org/10.1145/1015706.1015817)), with the plane, sphere, cylinder and other primitive proxies of Wu and Kobbelt ([*Structure recovery via hybrid variational surface approximation*, Computer Graphics Forum 2005](https://doi.org/10.1111/j.1467-8659.2005.00852.x)). With `t` the `tolerance` (times the surface's bounding-box diagonal when `relative`):

- **Proxies.** Each patch is explained by a plane, sphere, cylinder, cone or torus — [`fit_primitives`](/fit_primitives)' fit, accepted as `fit_primitives` accepts one (RMS distance within `t`, nine points in ten within it), so the refinement's verdict and the fit table it writes agree — or else by a *freeform* proxy, a quadric height field over its best-fit plane, or, for an isolated patch too small to fit anything, a *facets* proxy.
- **Lloyd step.** Each patch keeps the cell its proxy explains best as its seed, and the patches grow back from their seeds together, one cell at a time, cheapest first. A cell's cost against a proxy is its area times `d² / t² + w (1 - |n · n_p|)`, with `d` its centroid's distance to the proxy, `n` its normal, `n_p` the proxy's normal there and `w` the `normal_weight`. Growth never crosses a feature edge (a dihedral angle above `feature_angle`, as in `segment`), an open, non-manifold or inconsistently wound edge, or a cell outside every patch (`-1`), so a patch stays edge-connected and inside the smooth area it started in. The proxies are refitted, and a step that raised the energy (the sum of the costs) is undone.
- **Split** (`allow_split`). A freeform patch is searched for a primitive inside it by seeded RANSAC (Schnabel, Wahl and Klein, [*Efficient RANSAC for point-cloud shape detection*, Computer Graphics Forum 2007](https://doi.org/10.1111/j.1467-8659.2007.01016.x)): hypotheses fitted to the neighbourhoods of cells drawn by a splitmix64 generator seeded by `ransac_seed`, each scored by the area of the connected cells it explains within `t`. The best becomes its own patch, the rest of the freeform patch its connected pieces. A split is kept when the cells the primitive takes cost less under it than under the freeform proxy they leave.
- **Merge** (`allow_merge`). Two neighbouring patches one primitive explains together are merged, cheapest first (by the union's RMS), each patch in at most one merge per iteration, when the union costs no more than the two patches.
- **Absorb.** A patch with fewer than `min_patch_cells` cells joins the neighbour it shares the most smooth edges with (ties to the lowest patch), whatever its proxy: two triangles always lie on a plane, so a fragment the flooding left is not a face. A small patch bounded by creases has no smooth neighbour and stays.

Refinement stops when the partition no longer changes (`unchanged`, or `energy_rose` when the Lloyd step was undone and no split or merge applied), when an iteration adds less than `epsilon` to the share of the area primitives explain (`small_gain`), or after `max_iterations` iterations (`iteration_cap`).

| Option | Default | Meaning |
|---|---|---|
| `tolerance` | `1e-3` | how far a point may lie from its patch's proxy; > 0 |
| `relative` | `True` | `tolerance` is a fraction of the surface's bounding-box diagonal, not a length |
| `normal_weight` | `1.0` | the weight `w` of the normal term in a cell's cost; >= 0 |
| `epsilon` | `1e-3` | stop when an iteration adds less than this to the explained area share; >= 0 |
| `feature_angle` | `30` | the dihedral angle, in degrees, above which an edge is a feature no patch crosses; [0, 180] |
| `max_iterations` | `10` | the most iterations; >= 0 (`0` rewrites the outputs of the input partition) |
| `min_patch_cells` | `3` | a patch with fewer cells joins a smooth neighbour; >= 1 |
| `ransac_seed` | `0x6d657368696f2b2b` | the seed of the RANSAC split (the bytes of `"meshio++"`); in [0, 2⁶⁴) |
| `types` | all | the primitives a proxy may be, from `plane`, `sphere`, `cylinder`, `cone`, `torus`; an empty list leaves every patch freeform |
| `allow_split` | `True` | split a primitive out of a freeform patch |
| `allow_merge` | `True` | merge neighbouring patches one primitive explains |
| `prefix` | `"patch"` | the refined patches' cell regions are named `<prefix>:<k>` |

`feature_angle` should be the one the input was segmented with: a smaller one stops a patch at an edge `segment` let it cross, and a larger one lets refinement reach across an edge `segment` respected.

## Output

The input — points, connectivity, block order, data and other regions bit-identical — with:

| Output | Type | Content |
|---|---|---|
| `segment:patch` (cell data), `segment:adjacency` (field data), `<prefix>:<k>` (cell regions) | | rewritten for the refined patches, numbered in the order of each patch's lowest cell, so a partition refinement leaves alone comes back exactly as it went in |
| `fit:*` (cell and field data) | | [`fit_primitives`](/fit_primitives) run on the refined patches at the absolute tolerance `t`, so `fit_cad` reads it unchanged |
| `refine:patch_id` (field data) | `int64` | one row per patch, in label order |
| `refine:patch_proxy` | `int64` | the proxy: 0 facets, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus, 6 freeform |
| `refine:patch_iterations` | `int64` | the iterations whose step changed the patch |
| `refine:patch_rms`, `refine:patch_max` | `float64` | the RMS and largest distance of the patch's points from its proxy |
| `refine:patch_cells` | `int64` | the patch's cells |
| `refine:patch_demotion` | `int64` | why the patch has no primitive: 0 none (a primitive explains it), 1 too small to fit with no neighbour to join, 2 no primitive within the tolerance |
| `refine:explained_area_fraction` | `float64`, one value | the share of the area primitives explain |
| `refine:hausdorff` | `float64`, two values | the largest and the RMS distance of a skin point from its patch's proxy; the largest can exceed the tolerance, since a primitive is accepted with one point in ten beyond it |

With `return_report=True` the report holds `num_patches`, `num_iterations` (the Lloyd steps kept), `num_splits`, `num_merges`, `num_absorbed`, `num_primitive`, `num_freeform` and `num_facets` (patches by proxy), the absolute `tolerance` applied, `explained_area_fraction` and `initial_explained_area_fraction` (before and after), `hausdorff_max` and `hausdorff_rms`, `initial_energy` and `energy`, and `stop` (`"unchanged"`, `"small_gain"`, `"iteration_cap"` or `"energy_rose"`).

A mesh without `segment:patch` is refused with a `ValueError`, as is an option out of range; surface cells only, as `segment`.

## In fit_cad

`fit_cad(..., refine_segmentation=True)` (C++ `FitCadOptions::mRefineSegmentation`) runs `refine_patches` after any segmentation and fit it runs itself, with `max_refine_iterations`, `refine_epsilon`, `normal_weight`, `allow_merge`, `allow_split`, `ransac_seed`, `refine_tolerance` and `relative_tolerance`; its report then adds `refined`, `refine_iterations`, `num_refine_splits`, `num_refine_merges`, `initial_explained_area_fraction` and `explained_area_fraction`, and the solid carries the `refine:*` field data. The option is off by default, so a default `fit_cad` answer is what it was before it existed. See [fit CAD](/fit_cad).

## Determinism

The proxies of the patches a step changed are refitted in parallel, each serially and stored by patch; the flooding, the moves and every sum run serially in a fixed order, and the RANSAC draws come from the seeded splitmix64, never a standard-library distribution (whose output differs between standard libraries). The output is byte-identical for a seed on every backend and thread count. There is no numpy twin: without the compiled core the Python function raises `NotImplementedError`.

## Limitations

- **A genuinely freeform surface can be cut into strips.** Within a tolerance, a smooth freeform region is covered by many narrow planes, cylinders or tori, and the split takes them one by one; the explained area share rises but the faces are not what a designer drew. `allow_split=False` keeps such a region one freeform patch, which [`fit_cad`](/fit_cad#freeform-faces) turns into one B-spline face; general freeform faces, the second phase of this work, are the eventual alternative.
- **Blends and corners are not primitives.** A rolling-ball fillet of varying section, a vertex blend or a thread flank is explained by `fit_cad`'s feature fits, not by a primitive, so refinement treats it as freeform and may split it into tori and cones. `fit_cad(refine_segmentation=True)` guards against this by building the solid both ways and keeping the better one ([refining the patches first](/fit_cad#refining-the-patches-first)); a standalone `refine_patches` call does not.
- **The energy guarantee covers the Lloyd step only.** A Lloyd step never raises the energy and a merge never does; a split trades the freeform rest's cost for the cells a primitive takes, and absorbing a small patch is not an energy move at all, so the final energy can exceed the initial one when either applies.
- **Boundaries move within smooth areas only.** Growth never crosses a feature edge, so a boundary on one stays where `segment` put it, and a surface segmented at a `feature_angle` other than refinement's can keep boundaries refinement cannot move.

## CLI

```sh
meshioplusplus refine-patches IN OUT [-i FORMAT] [--tolerance T] [--absolute] \
    [--normal-weight W] [--epsilon E] [--feature-angle DEG] [--max-iterations N] \
    [--min-patch-cells N] [--ransac-seed S] [--types plane,cylinder,...] \
    [--no-split] [--no-merge] [--prefix P] [--json] [--quiet]
```

An input without `segment:patch` is segmented first with `segment`'s defaults at `--feature-angle`, which refinement uses too. `--ransac-seed` takes a decimal or `0x` hexadecimal number; `--json` prints the report. `to-cad --refine` runs the refinement inside [`fit_cad`](/fit_cad#cli).

## Other languages

```cpp
#include "meshioplusplus/operations/refine_patches.hpp"

meshioplusplus::RefinePatchesOptions opts;
opts.mTolerance = 2e-3;
opts.mAllowSplit = true;
meshioplusplus::RefinePatchesResult r = meshioplusplus::refine_patches(seg, opts);
r.mMesh; r.mExplainedAreaFraction; r.mStop;            // RefineStop::Unchanged, ...
```

```c
mio_refine_patches_opts opts;
mio_refine_patches_opts_init(&opts);
opts.tolerance = 2e-3;
mio_refine_patches_report report;
mio_mesh* ref = mio_refine_patches(seg, &opts, &report);
```

```fortran
ref = seg%refine_patches(tolerance=2e-3_real64, explained_area_fraction=share, stat=ierr)
```

```julia
r = refine_patches(seg; tolerance=2e-3)
r.mesh, r.explained_area_fraction
```

```r
r <- mio_refine_patches(seg, tolerance = 2e-3)
r$mesh
```

```js
const r = m.refinePatches(mesh);
r.mesh.cell_data['segment:patch'];
```

It is also a settings-pipeline step, `{"Op": "RefinePatches", "Tolerance": 0.001, "Relative": true, "NormalWeight": 1, "Epsilon": 0.001, "FeatureAngle": 30, "MaxIterations": 10, "MinPatchCells": 3, "RansacSeed": 7, "Types": ["plane", "cylinder"], "AllowSplit": true, "AllowMerge": true, "Prefix": "patch"}`, reporting `NumPatches`, `NumIterations`, `NumSplits`, `NumMerges`, `NumPrimitive`, `NumFreeform`, `ExplainedAreaFraction`, `HausdorffMax` and `Stop` (0 unchanged, 1 small gain, 2 iteration cap, 3 energy rose) — see [pipelines](/pipeline). `RansacSeed`, a JSON number, must be an integer in [0, 2⁵³]; left out, the default seed applies. The `FitCad` step takes `RefineSegmentation`, `MaxRefineIterations`, `RefineEpsilon`, `NormalWeight`, `AllowMerge`, `AllowSplit`, `RansacSeed`, `RefineTolerance` and `RelativeTolerance`, and reports `Refined`, `NumRefineSplits`, `NumRefineMerges` and `ExplainedAreaFraction` when it refined. The MCP server exposes it as the tool `refine_patches`, which segments an input without `segment:patch` first, and `fit_cad`'s refinement through the same arguments as the Python function.
