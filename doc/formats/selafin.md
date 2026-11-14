# Selafin — TELEMAC results (`.slf`)

The TELEMAC-MASCARET result and geometry format (SERAFIN/Selafin), roadmap 2.4.3.4. Both engines read and write it (`formats/selafin.cpp`, `meshioplusplus.selafin`).

| | |
|---|---|
| **Format name** | `selafin` |
| **Extensions** | `.slf`, `.srf` |
| **Read / Write** | ✓ / ✓ (2-D, one step) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("r2d_case.slf")                     # the first step
last = meshioplusplus.selafin.read("r2d_case.slf", time_step=-1)
times = meshioplusplus.selafin.time_values("r2d_case.slf")
```

The file is big-endian Fortran sequential records: an 80-character title (its last eight `SERAFIN ` or `SERAFIND`), `NBV1 NBV2`, a 32-character name and unit per variable, `IPARAM` (a date record follows when its tenth value is 1), `NELEM NPOIN NDP 1`, the 1-based `IKLE`, `IPOBO`, `X` and `Y`, then per step a time record and one record per variable. Single or double precision is read from the size of the `X` record.

- `NDP` 3 is a `triangle` block, 4 a `quad` block (`tetra` when `IPARAM[6]`, the plane count, exceeds 1), 6 a `wedge` block (TELEMAC-3D prisms); other counts are refused by name.
- Each variable is point data under its name (blanks trimmed), its unit `field_data["selafin:unit:<name>"]`; in a 3-D file the `ELEVATION Z` (`COTE Z`) variable of the chosen step is also the points' z.
- `IPOBO` is the Int32 point data `selafin:ipobo`, the title `field_data["selafin:title"]` (bytes).
- `time_step` picks the step (0 the first, negatives from the end), its time `field_data["meshio:time"]`. Only complete steps count, so a file a solver left half-written reads up to its last whole step.

## Writing

```python
meshioplusplus.write("out.slf", mesh)
```

A 2-D double-precision (`SERAFIND`) file with one step: one triangle or quad block, every floating-point point data array as a variable (names and units cut to 16 characters), `selafin:ipobo` when present (0 otherwise), and the time from `field_data["meshio:time"]` (0 otherwise). Any other mesh is refused by name. Both engines write the same bytes.

## Verification

Single- and double-precision, dated, multi-step, truncated and 3-D prism files are built in `tests/python/test_selafin.py` from the layout above; no TELEMAC-written `r2d_*.slf` has been read yet (TELEMAC's examples are GPL data, so they are not copied here).
