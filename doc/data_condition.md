# Value conditioning (clamp / normalize / standardize)

`data_condition` conditions the *values* of data arrays without changing the geometry. It is a [data operation](/data_operations), not a file format.

| Mode | Transform |
| --- | --- |
| `"clamp"` | `x → min(max(x, lo), hi)` |
| `"normalize"` | affine map from the array's own `[min, max]` onto `[lo, hi]` (default `[0, 1]`) |
| `"standardize"` | `(x − mean) / stddev` — zero mean, unit standard deviation |

![A field normalized to `[0, 1]`](/images/data_condition_normalize.png)

```python
import meshioplusplus as mp

mesh = mp.read("part.vtu")

mesh = mp.data_condition(mesh, "point", ["T"], mode="clamp", lo=0.0, hi=100.0)
mesh = mp.data_condition(mesh, "cell", ["damage"], mode="normalize")
mesh = mp.data_condition(mesh, "point", ["T"], mode="standardize")

# Keep the original and write the result under a new name.
mesh = mp.data_condition(mesh, "point", ["T"], mode="normalize", suffix="_n")
```

Omit the name list to condition every array at that location.

## Scope: components vs magnitude

`scope="component"` (the default) treats each trailing component independently, so each column of a vector array is conditioned on its own statistics.

`scope="magnitude"` computes the statistics over each row's Euclidean magnitude and rescales whole rows, **preserving direction** — the natural choice for a velocity or displacement field.

```python
mp.data_condition(mesh, "point", ["velocity"], mode="normalize",
                  scope="magnitude")
```

## Statistics

Statistics come from the finite values only — see the [NaN policy](/data_operations#nan-policy). For a `cell_data` array they are computed **jointly across all cell blocks** before a single transform is applied to each: normalizing one block against its own extremes would make the pieces mutually incomparable.

Degenerate reductions are handled rather than producing `NaN`:

- `normalize` on a constant array (or one with no finite values) fills the target lower bound and warns.
- `standardize` with zero standard deviation fills `0` and warns.

## Data types

`clamp` preserves the input dtype, so clamping an `int32` array leaves it `int32` (pass `preserve_dtype=False` to force `float64`). `normalize` and `standardize` always produce `float64` — a rescaled integer field is not an integer field.

## CLI

```bash
meshioplusplus data clamp     in.vtu out.vtu --point T --min 0 --max 100
meshioplusplus data normalize in.vtu out.vtu --cell damage --to 0,1
meshioplusplus data normalize in.vtu out.vtu --point T --zero-mean
```

`--magnitude` selects the magnitude scope, `--suffix` writes to a new name instead of in place, and `--nan ignore|replace|fail` (with `--nan-value`) selects the non-finite policy. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

- **C API** — `mio_data_condition(mesh, location, names, count, mode, lo, hi, scope, nan_policy, nan_replacement, suffix)` with the `MIO_COND_*`, `MIO_SCOPE_*` and `MIO_NAN_*` enumerations. See the [C API reference](/c_api).
- **Fortran** — `m%data_condition(MIO_DATA_POINT, ["T"], mode=MIO_COND_NORMALIZE)`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataCondition(mesh, "point", ["T"], "normalize", 0, 1, "component", "ignore", 0, "")`. See the [WebAssembly reference](/wasm).
