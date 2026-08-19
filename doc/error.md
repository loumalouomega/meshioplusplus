# Error estimation — Zienkiewicz-Zhu recovery + marking

`estimate_error(mesh, array, marking=…)` estimates the per-cell recovered-gradient error of a `point_data` field, and can optionally mark cells for refinement. It is a mesh **operation**, not a file format, uses only standard C++/numpy, and runs under every mesh backend.

meshio++ could already differentiate a field ([`gradient`](/gradient)) and selectively refine a mesh by any scalar `cell_data` predicate ([`refine`'s `--where`](/refine)). This is the piece that closes the loop: it produces the indicator `refine --where` consumes, so the whole adaptive-mesh-refinement cycle — estimate, mark, refine — is three verbs, none of them new numerical code.

![The bracket's recovered-gradient error indicator, and the mesh refined where estimate_error marked it](/images/estimate_error.png)

```python
import meshioplusplus as mp

mesh = mp.read("solution.vtu")

out = mp.estimate_error(mesh, "T")                 # cell_data["error:zz"]
out, report = mp.estimate_error(
    mesh, "T", marking="dorfler", marking_value=0.6, return_report=True,
)
report["global_error"]  # sqrt(sum of eta_K^2) over evaluable cells
report["num_marked"]    # cells marked in cell_data["error:marked"]

adapted = mp.refine(out, where="error:marked > 0.5")
```

## Where it lives, and why

Like [`gradient`](/gradient), this consumes and produces data arrays but **reads** geometry (it composes `gradient` and the point↔cell averaging round trip), so it lives in the mesh-operations layer (`operations/error.hpp`), not in the `data_*` family. It is reachable as `meshioplusplus data estimate-error` in both CLIs for the same reason `gradient` is — that is where a user looks for it.

## Composition, not a new kernel

The indicator is the standard Zienkiewicz-Zhu (ZZ) superconvergent-patch recovery estimator, built entirely from operations that already ship:

1. **Differentiate** the field once per cell — `gradient(mLocation=Cell)`, Green-Gauss — giving the *raw* (unrecovered) per-cell gradient `g_K`.
2. **Recover** a smoothed nodal gradient by averaging the raw gradient onto the points, measure-weighted (`cell_data_to_point_data(weight=Measure)`), then back onto the cells (`point_data_to_cell_data`), giving `G_K`.
3. **Indicate**: `eta_K = sqrt(|measure_K| · sum((G_K − g_K)²))` — the recovered gradient error in the energy-like (H1-seminorm) norm. `measure_K` reuses the same primitive `data_average`'s `Measure` weighting already used in step 2, so the two agree by construction.

`global_error` is `sqrt(sum_K eta_K²)` over the evaluable cells — the estimated global error in that same norm.

This is why the operation exists as composition rather than a bespoke computation: the recovery step *is* the existing averaging round trip, reused exactly as `gradient(location="point")` itself composes `cell_data_to_point_data`.

## Marking

A marking pass turns the continuous indicator into a boolean `cell_data` array, `error:marked` (one array per block, Int64 0/1), attached only when `marking != "none"` — so [`refine`](/refine)'s own `where` selector needs no change at all:

| `marking` | Meaning of `marking_value` | Marks |
| --- | --- | --- |
| `"none"` (default) | ignored | nothing — estimate only |
| `"absolute"` | an indicator threshold | every evaluable cell whose indicator exceeds it |
| `"fraction"` | a fraction in `(0, 1]` of cells | the largest-indicator fraction, `k = int(ratio · N + 0.5)` |
| `"dorfler"` | a bulk fraction theta in `(0, 1]` | the smallest indicator-descending prefix whose cumulative `eta_K²` covers `theta` of the total |

`"dorfler"` (the Dörfler / bulk-chunk criterion) is the usual choice for genuine AMR loops: it marks *enough* cells to reduce a fixed fraction of the total error, rather than a fixed *count* of cells (`"fraction"`) or an absolute, problem-scale-dependent threshold (`"absolute"`).

Ties in the `"fraction"`/`"dorfler"` descending ranking break on ascending global cell index, so the result does not depend on how the per-cell work happened to schedule.

## NaN, never a guess

A cell that cannot be evaluated — an unsupported/degenerate `gradient` cell, a recovery neighbourhood with no finite contribution, or a cell whose `measure_K` cannot be computed — reads:

- `NaN` in `error:zz` (the continuous indicator array), and
- **`0`, never `NaN`**, in `error:marked` (a `NaN` there would make every `--where` comparison false anyway, but writing it as `0` is what keeps the array's own dtype meaningfully boolean).

Such cells are counted in `num_skipped` and excluded from `global_error` and from every marking policy's cell count.

## Naming

`error:zz` (Float64) is always attached; `error:marked` (Int64) only when marking is requested. Default names are fixed and can be overridden:

```python
out = mp.estimate_error(mesh, "T", output="my_indicator")
out = mp.estimate_error(
    mesh, "T", marking="absolute", marking_value=1.0, marked="hot_cells",
)
```

Unlike `gradient`'s deliberately unusual `<input>:<operator>` convention, these are fixed constants (`error:zz`, `error:marked`) following the repo's *usual* `prefix:name` scheme (`quality:*`, `partition:*`, `refine:*`) — the estimator has exactly one continuous output and one boolean output regardless of which field produced them, so there is nothing for the name to disambiguate.

## Byte-identity and its limits

The raw (pre-recovery) gradient stays bit-exact across the C++/numpy boundary, since `gradient` itself needs no tolerance. The **recovery** step, and so the final indicator, matches the numpy twin only to within a numeric tolerance — the same, already-accepted precedent `data_average`'s own `Measure` weighting has (`test_data_location.py::test_cpp_matches_python` compares it with `np.allclose`, not exact equality). Output stays byte-identical across the three mesh backends and across thread counts.

The numpy twin raises `NotImplementedError` for:

- a mesh with a ragged/polygon/polyhedron block, and
- a mesh with a quadratic 3-D cell type (`tetra10`, `hexahedron20`/`27`, `wedge15`),

both for the same reason `_convert_cells.py`'s simplexify twin is exempt: `detail::cell_measure` repairs a closed 3-D cell's winding via `orient_rings` — a discrete branch on a sign — and the Python reference does not attempt to replicate it. A mesh of ordinary linear cells (`tetra`, `hexahedron`, `wedge`, `pyramid`, and any 1-D/2-D type) is unaffected.

## Worked composition — the adaptive loop

```python
import meshioplusplus as mp

mesh = mp.read("solution.vtu")
estimated, report = mp.estimate_error(
    mesh, "T", marking="dorfler", marking_value=0.6, return_report=True,
)
print(f"global error: {report['global_error']:.4g}, "
      f"marked {report['num_marked']} of {len(mesh.cells[0].data)} cells")

adapted = mp.refine(estimated, where="error:marked > 0.5", closure="redgreen")
mp.write("adapted.vtu", adapted)
```

```sh
meshioplusplus data estimate-error solution.vtu estimated.vtu \
    --array T --marking dorfler --marking-value 0.6
meshioplusplus refine estimated.vtu adapted.vtu \
    --where "error:marked > 0.5" --closure redgreen
```

Re-running the loop on `adapted.vtu` (once a fresh solve has produced a new `T`) picks up where the last pass left off, since `refine`'s own [`refine:cell_id`/`refine:parent_id`](/refine#refinecell_id-and-refineparent_id) (`record_hierarchy=True`) tracks parent/child relationships across passes if you need them.

## CLI

```sh
meshioplusplus data estimate-error IN OUT --array NAME \
    [--method zz] \
    [--marking none|absolute|fraction|dorfler] \
    [--marking-value V] \
    [--output NAME] [--marked NAME] [--overwrite] [--quiet]
```

Both CLIs produce the same array names and report the same counters (subject to the recovery step's own numeric tolerance, above). See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

```c
double global_error = 0.0;
int64_t skipped = 0, marked = 0;
mio_mesh* out = mio_estimate_error(mesh, "T", "zz", "absolute", 1e-6, NULL,
                                   NULL, 0, &global_error, &skipped, &marked);
```

```fortran
type(mio_mesh) :: e
real(real64) :: gerr
integer(int64) :: nskip, nmark
e = m%estimate_error('T', marking='absolute', marking_value=1.0d-6, &
                     global_error=gerr, num_skipped=nskip, num_marked=nmark)
```

```julia
e = estimate_error(mesh, "T"; marking=:dorfler, marking_value=0.6)
e.mesh, e.global_error, e.num_marked
```

```r
e <- mio_estimate_error(mesh, "T", marking = "dorfler", marking_value = 0.6)
e$mesh; e$global_error; e$num_marked
```

```js
const e = await m.estimateError(mesh, 'T', 'zz', 'dorfler', 0.6);
e.mesh.cell_data['error:marked'];
```

On the flat ABIs (C, Fortran, Julia, R, WASM) the counters are out-parameters or fields of a returned record, never an opaque handle — `gradient`/`smooth`'s shape, since there are no index maps to hand back.

This operation is also reachable as an `EstimateError` step in the [settings pipeline](/pipeline) and in the browser viewer's `convertSurfaceOps` chain (`{op: 'estimateError', array: 'T', marking: 'dorfler', markingValue: 0.6}`), since it is a pure data step: geometry is untouched.
