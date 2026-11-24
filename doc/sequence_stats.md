# Statistics and derivatives over time

Two streaming operations over the steps of a [sequence](/sequences): the per-point (or per-cell) **statistics** of data arrays, and the **time derivative** of an array. Both take what a sequence takes (a glob, a list of paths or a `TimeSeries`) and hold only a few steps in memory at a time, so the length of the sequence is not limited by it.

```python
import meshioplusplus as mp

stats = mp.sequence_statistics("run_*.vtu", ["U", "p"], ("mean", "std", "min", "max"))
stats.point_data["U:mean"]          # the mean velocity at every point over the run
stats.point_data["p:std"]           # the standard deviation of the pressure

steps = mp.sequence_derivative("run_*.vtu", "p")      # a generator of (time, mesh)
mp.write_sequence("dpdt_{step}.vtu", steps)           # each mesh carries "p:dt"
```

This is the counterpart of VTK's `vtkTemporalStatistics` and a time-derivative filter. It is a **Python-only layer**: the sequence readers it builds on are Python's, and the flat bindings (C, Fortran, Julia, R, Octave) already traverse a sequence step by step through `mio_sequence_read`, which is the hook for doing the same there.

## Statistics

`sequence_statistics(source, arrays=None, stats=("mean", "std"), *, location="point", ddof=0, …)` returns **the first step's mesh with its data replaced by the statistics**: for every array `a` and requested statistic `s`, `<a>:<s>` at the location (a `cell_data` statistic has one array per cell block, like any cell array).

| Statistic | Meaning |
| --- | --- |
| `mean` | the mean over the steps, by Welford's update |
| `std` | the standard deviation, `sqrt(M2 / (n − ddof))`: **population** (`ddof=0`, as `numpy.std`) by default, `ddof=1` for the sample's; `0` for a single step |
| `min`, `max` | the least and greatest value |
| `sum` | the sum |

- **One step in memory at a time**: Welford's update keeps the mean and the sum of squared differences, so the variance is accurate for a series with a large offset where the sum-of-squares formula loses it (a test pins `1e8 + t`). Memory is about the arrays of two steps.
- Every step must have the **same points, cells and array shapes**; a step that differs is refused by name, with its index and time, never skipped. `arrays` defaults to every numeric array at the location of the first step; `field_data` is refused (no per-point correspondence).
- A **non-finite value propagates** into its point's `mean`, `std` and `sum` (and `min`/`max` follow `numpy.minimum`): nothing is filtered, so a `NaN` is visible instead of silently dropped.
- Equal weights per step. Unequal time spacing is *not* a weight: a time-weighted mean is `data_calc` over the steps' integral, not this.

`return_report=True` adds `{"num_steps", "first_time", "last_time", "arrays"}`.

## Time derivative

`sequence_derivative(source, array, *, location="point", scheme="central", output=None, …)` is a **generator of `(time, mesh)`**, the mesh being the step's own with the derivative added as `<array>:dt` (or `output`). It holds at most three steps, so it streams into [`write_sequence`](/sequences) or a loop.

| `scheme` | Formula |
| --- | --- |
| `"central"` (default) | the second-order three-point formula for unequal steps, `−h₂/(h₁(h₁+h₂)) fᵢ₋₁ + (h₂−h₁)/(h₁h₂) fᵢ + h₁/(h₂(h₁+h₂)) fᵢ₊₁` with `h₁ = tᵢ − tᵢ₋₁`, `h₂ = tᵢ₊₁ − tᵢ` (the central difference when the steps are equal); **one-sided and first order at the first and the last step** |
| `"forward"` | `(fᵢ₊₁ − fᵢ) / (tᵢ₊₁ − tᵢ)`; a backward difference at the last step |

The derivative of a field linear in time is exact on both schemes and on unequal steps; the central formula is exact on a quadratic in the interior. It needs **at least two steps**, **strictly increasing times** and one topology, each refused by name otherwise. The array is named `<array>:dt` (the `name:suffix` convention the other derived arrays follow), not `d<array>/dt`.

## CLI

```bash
meshioplusplus sequence-stats 'run_*.vtu' stats.vtu --arrays U,p --stats mean,std,max
meshioplusplus sequence-derivative 'run_*.vtu' 'dpdt_{step}.vtu' --array p --scheme central
```

`sequence-stats` takes `--arrays`, `--stats`, `--ddof`, `--location` and `--time-from`, and `--json` prints `output`, `num_steps`, `first_time`, `last_time` and `arrays`; `sequence-derivative` takes `--array`, `--scheme`, `--output-name` and the same `--location` and `--time-from`, and `--json` prints `output`, `num_steps` and `files`. Both are Python CLI verbs only. The MCP server carries them as the `sequence_statistics` and `sequence_derivative` tools.
