# Field integration

`data_integrate(mesh, arrays=…)` computes a cell-measure-weighted **total** and **mean** of one or more `cell_data` arrays — the natural companion to [`gradient`](/gradient): `gradient` differentiates a field, this integrates one. It is a mesh **operation**, not a file format, uses only standard C++/numpy, and runs under every mesh backend.

meshio++ could already differentiate ([`gradient`](/gradient)), estimate error from a field ([error estimation](/error)) and summarize a data array's raw statistics ([`data_info`](/data_info)) — but none of those gives a physically meaningful total. A density field's total mass, a heat-flux field's total power, a material tag's occupied volume: all three need the value **weighted by the cell it lives on**, which `data_info`'s plain unweighted mean of raw values cannot give.

```python
import meshioplusplus as mp

mesh = mp.read("solution.vtu")

report = mp.data_integrate(mesh, arrays=["density"])
arr = report[0]
arr["domain"]["total_per_component"]   # sum(value * |measure|) -- total mass
arr["domain"]["mean_per_component"]    # the measure-weighted mean density
arr["regions"]                         # the same, independently, per named Cell region
```

## Where it lives, and why

The five [data operations](/data_operations) are *defined* by never touching geometry. This one consumes data arrays but **reads** geometry (each cell's own length/area/volume via `detail::cell_measure`), so it lives in the mesh-operations layer (`operations/data_integrate.hpp`) alongside `gradient` and `estimate_error`, not in the `data_*` family.

It is still reachable as `meshioplusplus data integrate` in both CLIs, because that is where a user looks for it — the same deliberate inconsistency `gradient` and `estimate_error` already accept in the `data` group.

## What it computes

For each requested array and each component `k`:

```
domain_measure[k] = Σ |measure(cell)|              over cells with a computable
                                                     measure AND a finite value[k]
total[k]          = Σ value[k](cell) * |measure(cell)|   over the same cells
mean[k]           = total[k] / domain_measure[k]          (NaN if the denominator is 0)
```

`measure` is the cell's own length (1-D), area (2-D) or volume (3-D), the same primitive `gradient`'s Green-Gauss method and `compute_quality` use.

**Exclusion is symmetric between geometry and value, and deliberately does not fall back to a unit weight.** A cell whose measure cannot be computed (a ragged or unsupported-type block, or a degenerate one) is excluded from **both** the numerator and the denominator — never given a fallback weight of 1, unlike `cell_data_to_point_data`'s own `Measure` weighting, because a silent unit-weight substitution would corrupt a physical total in a way it only softens an average. Excluded cells are counted in `num_skipped`.

A **non-finite value** in one component excludes that cell from that *component's* numerator **and** denominator too — not just zeroed into the numerator, or the mean would be silently biased. Counted per component in `num_nan_per_component`. This mirrors the geometric exclusion rule one level down, component by component: a NaN in component 1 of a 3-component array does not disturb components 0 and 2.

There is deliberately **no whole-array (cross-component) total or mean**, unlike `data_info`'s collapsed statistics over a flattened numeric stream. Summing a vector field's x/y/z totals into one number has no general physical meaning, unlike a flattened mean, which needs only "these are numbers".

There is deliberately **no `nan_policy`**, like `gradient` and `data_info`: this is a reduction with nothing to exclude a value *from* — non-finite values are always excluded and always counted, with no configurable policy.

## `cell_data` only

A `point_data`-only name fails, naming `point_data_to_cell_data` (CLI `data to-cell`) as the fix — the mirror image of `gradient`'s own contract, which throws on a `cell_data` name pointing at `data to-point`. Proper `point_data` integration needs shape-function quadrature, a separate and larger job than this one.

`arrays` empty (or omitted) means every `cell_data` array, in sorted name order.

## Regions

Every named `Cell` region present on the mesh gets its **own independent** entry in `regions`, alongside the whole-mesh `domain` entry — mirroring [`split`](/split)'s own "split by region" contract: regions are **not a partition**. A cell belonging to two regions contributes fully to both; a cell in none contributes to neither.

`Point`/`Side` regions are skipped, exactly as `split` skips them for the same reason — there is no cell-measure meaning for a group of points or facets.

```python
mesh.regions.append(mp.Region(name="inlet", kind="cell", entries=[0, 1, 4]))
report = mp.data_integrate(mesh, arrays=["density"])
for region in report[0]["regions"]:
    print(region["name"], region["total_per_component"])
```

## The mesh is never modified

Like `data_info` and `compute_stats`, this returns a read-only report; the input mesh comes back untouched. There is no counterpart to `gradient`'s mesh-with-a-new-array-attached shape.

## Determinism

Output is byte-identical across the three mesh backends, across thread counts, and across the C++-core / numpy-fallback boundary (pinned by `tests/python/test_data_integrate.py::test_cpp_matches_python`). Per-cell measures are computed **once**, shared across every array and every region, via `parallel_for`; every weighted-sum reduction is chunked-`parallel_for`-then-serially-merged (`detail::accumulate_weighted`), the same idiom `data_info`'s `detail::accumulate_stats` already uses.

## CLI

```sh
meshioplusplus data integrate IN --array NAME [--array NAME2 ...] [--json]
```

Report-only, like `data info` — there is no `OUT`, since the mesh is never modified. Both CLIs produce the same report shape. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

```c
mio_data_integrate* result = mio_data_integrate_create(mesh, names, count);
mio_field_integral_info entry;
mio_data_integrate_entry(result, 0, &entry);   // whole-mesh: num_cells, num_skipped, ...
double total, mean, domain_measure; int64_t num_nan;
mio_data_integrate_component(result, 0, 0, &total, &mean, &domain_measure, &num_nan);
mio_data_integrate_free(result);
```

```fortran
type(mio_field_integral_info), allocatable :: fi(:)
real(real64), allocatable :: totals(:), means(:), domain_measures(:)
fi = m%data_integrate(arrays=['density '], totals=totals, means=means, &
                      domain_measures=domain_measures)
! Per-region breakdown of one named array is a separate call:
fi = m%data_integrate_region('density', totals=totals)
```

```julia
report = data_integrate(mesh, ["density"])
report[1].domain.components[1].total
report[1].regions   # one entry per named Cell region
```

```r
report <- mio_data_integrate(mesh, "density")
report[[1]]$domain$components["total"]
report[[1]]$regions   # a list, one per named Cell region
```

```js
const report = m.dataIntegrate(mesh, ['density']);
report[0].domain.totalPerComponent;
report[0].regions;   // one entry per named Cell region
```

Two conventions differ across these surfaces and are worth pinning down:

- On the flat ABIs (C, Fortran, Julia, R, WASM) this is a **read-only report handle/value**, never an opaque mesh-carrying result — there is no mesh to hand back, so `data_info`'s shape is the right one, extended with a region axis.
- **Fortran alone splits domain and region access into two calls** (`data_integrate` for the whole mesh, `data_integrate_region` for one named array's per-region breakdown) rather than nesting three axes (array × region × component) into one call's optional out-arguments.
