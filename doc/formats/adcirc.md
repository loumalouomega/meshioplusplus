# ADCIRC — `fort.14` grids and `fort.63`/`fort.64` results

The [ADCIRC](https://adcirc.org) coastal circulation model's grid (`fort.14`, also `.grd`) and its ASCII global output (roadmap 2.4.3.3). Both engines read and write the grid (`formats/adcirc.cpp`, `meshioplusplus.adcirc`).

| | |
|---|---|
| **Format name** | `adcirc` |
| **Extensions** | `.14`, `.grd` |
| **Read / Write** | ✓ (grid and results) / ✓ (grid) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("run/fort.14")                      # fort.63/fort.64 beside it too
last = meshioplusplus.adcirc.read("run/fort.14", time_step=-1)
grid = meshioplusplus.adcirc.read("grid.grd", results=["run/fort.63"])
times = meshioplusplus.adcirc.time_values("run/fort.14")
```

- The grid: a title, `NE NP`, the nodes (`id x y depth`) and the triangles (`id 3 n1 n2 n3`). Node ids may skip; text after the numbers on a line (`! comment`, `= NE NP`) is ignored. Points are `(x, y, 0)` with `depth` point data, one `triangle` block, and the title in `field_data["adcirc:title"]` (bytes).
- The boundaries: each open boundary is a point region `open_boundary_<k>`, each land boundary `land_boundary_<k>`, tagged with its type (`IBTYPEE`/`IBTYPE`; an open boundary with no type is 0). The barrier columns after a node of a type 3, 4 or 5 land boundary (13, 23, 24, …) are read past. A region is a sorted set, so each string's own order is also `field_data["adcirc:boundary:<name>"]`.
- Results: `fort.63` (`zeta`, elevation) and `fort.64` (`velocity`, two components) next to a file named `fort.14` are read automatically (`siblings=False` skips them), or any files through `results=` (`fort.73` is `pressure`, `fort.74` `wind_velocity`, any other `adcirc:fort.NN`). Full and sparse (`TIME IT NNONDEFAULT DEFAULT`) datasets are read. `time_step` picks the dataset (0 the first, negatives from the end), its time `field_data["meshio:time"]`; the result files must have as many datasets as each other.

## Writing

```python
meshioplusplus.write("fort.14", mesh)
```

The grid only: a triangle mesh, `depth` point data (0 when absent), and the `open_boundary_*`/`land_boundary_*` point regions as boundaries, in their number's order and each in its kept order. A land boundary whose type needs barrier columns is refused by name, since they are not kept; other cells are refused too. Both engines write the same bytes (`%.17g`).

## Verification

The layouts are built in `tests/python/test_adcirc.py` and `tests/cpp/test_coastal_formats.cpp` from the ADCIRC manual's description; no ADCIRC-written case (the quarter-annulus example) has been read yet.
