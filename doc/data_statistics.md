# Data statistics per region

`data_statistics` summarises data arrays beyond [`data_info`](/cli#data-info)'s minimum, maximum and mean: the **sum, variance and standard deviation, quantiles and a histogram**, for the whole mesh and for every named [region](/regions). It is the counterpart of VTK's `vtkDescriptiveStatistics` and `vtkComputeQuantiles`, and a **Python-only layer**: `data_info` and its structs are unchanged, so no binding or ABI moves.

```python
import meshioplusplus as mp

mesh = mp.read("flow.vtu")
for a in mp.data_statistics(mesh, ["p"], quantiles=(0.05, 0.5, 0.95), bins=20):
    print(a["name"], a["domain"]["mean"], a["domain"]["quantiles"]["0.5"])
    for r in a["regions"]:
        print("  ", r["name"], r["std"], r["histogram"][0]["counts"])
```

`data_statistics(mesh, arrays=None, *, location="point", per_region=True, quantiles=(0.05, 0.5, 0.95), bins=0, ddof=0, weighted=False)` returns one dict per numeric array, `{"name", "location", "num_components", "domain", "regions"}`. `domain` and every entry of `regions` (which also has a `name`) hold, **per component**:

| Key | Meaning |
| --- | --- |
| `count`, `num_nan`, `num_inf` | the finite values the statistics use, and the others, counted rather than dropped silently |
| `sum`, `mean`, `variance`, `std` | over the finite values; `ddof` is the variance's delta degrees of freedom (0, the population, as `numpy.var`) |
| `min`, `max` | the extremes |
| `quantiles` | a dict keyed by the quantile (`"0.5"`), linear interpolation as `numpy.quantile` |
| `histogram` | with `bins`, `{"edges", "counts"}` over the component's own range |

An array with no finite value gives `NaN` statistics and a count of 0.

## Regions

- **Point data** has a row for every point region and, through the nodes of its cells (each point once), for every cell region. A side region has no point image and is skipped.
- **Cell data** has a row for every cell region.
- Where several regions share a name and kind, the one with the least `(dim, tag)` is used, as everywhere else. `per_region=False` gives the whole mesh only.

## Weighted statistics

`weighted=True` (cell location) weights `mean`, `variance` and `std` by each cell's unsigned measure, so the mean equals the mean of [`data_integrate`](/field_integration). The quantiles, the histogram and the extremes stay unweighted, and `ddof` is ignored.

## Errors

A name that is not an array at the location, `field_data` (it has no per-point or per-cell correspondence), a quantile outside `[0, 1]`, negative `bins` and `weighted` outside the cell location are each refused with a `ValueError`.

## CLI and MCP

```bash
meshioplusplus data stats flow.vtu --arrays p,U --quantiles 0.05,0.5,0.95 --bins 20 --json
meshioplusplus data stats flow.vtu --location cell --weighted --no-regions
```

`--json` prints the list above. The verb is Python-only; the MCP server carries it as the `data_statistics` tool.
