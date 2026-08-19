# Mesh quality metrics

`meshioplusplus.compute_quality(mesh)` scores every cell of every block on a fixed set of geometric quality metrics and returns an aggregate report; `meshioplusplus.attach_quality(mesh)` returns a copy of the mesh with those metrics attached as `cell_data` (ready to write out and visualise). Metrics are evaluated on each cell's **corner nodes** (quadratic variants reduce to their linear parent), using only standard C++/numpy math, so they run under every mesh backend.

![Scaled-Jacobian distribution across every cell in a mesh](/images/quality_histogram.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("part.vtu")
report = meshioplusplus.compute_quality(mesh)
print(report["num_inverted"], "inverted cells")
print(report["metrics"]["quality:scaled_jacobian"]["min"])

annotated = meshioplusplus.attach_quality(mesh)   # metrics as cell_data
meshioplusplus.write("part_quality.vtu", annotated)
```

## Metrics

| metric | triangle | quad | tetra | hexahedron | wedge | pyramid |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| `quality:volume` (area in 2D) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `quality:scaled_jacobian` | – | – | ✓ | ✓ | ✓ | ✓ |
| `quality:aspect_ratio` | ✓ | ✓ | ✓ | – | – | – |
| `quality:skewness` | – | ✓ | – | ✓ | – | – |
| `quality:min_angle` / `quality:max_angle` | ✓ | – | – | – | – | – |
| `quality:warpage` | – | ✓ | – | – | – | – |
| `quality:min_dihedral` / `quality:max_dihedral` | – | – | ✓ | – | – | – |
| `quality:inverted` / `quality:degenerate` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

A metric that does not apply to a cell type is `NaN` for that cell. Ideal (perfectly regular) elements give the ideal value: `scaled_jacobian = 1`, `aspect_ratio = 1`, `skewness = 0`, `warpage = 0`, triangle angles `= 60°`, tetra dihedral `≈ 70.53°`.

- **`quality:inverted`** is `1.0` when the signed volume/area is negative (a 2D triangle/quad embedded in 3D is unsigned, so never flagged).
- **`quality:degenerate`** is `1.0` when the cell is near-zero (collapsed volume/area, a zero-length edge, or a vanishing Jacobian); its divide-prone metrics are then `NaN`.

Quadratic and higher-order variants (`tetra10`, `hexahedron20`/`27`, `quad8`/`9`, `triangle6`, `wedge15`, `pyramid13`/`14`) are scored on their corners. **Polyhedron blocks** get a reduced set since v9.16.0 — `volume`, `inverted` and `degenerate`, computed through `detail/polyhedron.hpp`; every metric defined against a reference element the cell does not have stays `NaN`, and `inverted` there means *unorientable* rather than negative-volume (see [Polyhedra](/polyhedra)). Ragged polygon and other unsupported blocks contribute all-`NaN` arrays and are excluded from the summaries and counts.

## The report

`compute_quality` returns a dict:

- `num_cells`, `num_inverted`, `num_degenerate` — integer counts.
- `metrics` — `metric name → {min, max, mean, count, histogram}` (a 10-bin histogram over `[min, max]`; `count` is the number of cells with a finite value).
- `cell_arrays` — `metric name → list of arrays` (one `(n, 1)` array per cell block, `NaN` where the metric does not apply) — the arrays `attach_quality` writes into `cell_data`.

![A mesh's boundary coloured by attach_quality's scaled Jacobian](/images/quality_boundary.png)

## Cross-language

Available in every binding surface: Python (`compute_quality`/`attach_quality`), the C API (`mio_attach_quality` + `mio_quality_counts`), Fortran (`mesh%attach_quality()` / `mesh%quality_counts(...)`), and WASM (`attachQuality`), plus the CLI verb `meshioplusplus quality`.

## Implementation

The metrics are in `src/cpp/src/operations/quality.cpp` (over the uniform mesh API, with the per-cell loop parallelised) with reusable vector math in `src/cpp/include/meshioplusplus/detail/geometry.hpp`. A pure-numpy twin (`meshioplusplus._quality._compute_quality_py`) mirrors the formulas and is checked against the C++ core in `tests/python/test_quality.py`.
