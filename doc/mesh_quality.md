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

## Quality gate

`check_quality(mesh, require)` (v16.24.0) turns the metrics into a pass/fail answer — the check a CI job over meshes scripts. Each threshold is tested against the **per-cell** values, not the report's histograms, whose bins span each metric's own data range and so cannot say how many cells fall below a bound.

```python
import meshioplusplus as mp

r = mp.check_quality(mesh, "scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%")
r["passed"]                    # every check passed
for c in r["checks"]:
    print(c["name"], c["violations"], c["evaluated"], c["worst"], c["worst_cell"])
```

**The specification text**, which every surface reads the same way: clauses separated by `;`, `,` or newlines (`#` starts a comment, so a gate file is one clause per line), each `METRIC >= VALUE` or `METRIC <= VALUE`, optionally followed by `@ FRACTION` — the fraction of evaluated cells allowed to violate the bound, a number in `[0, 1]` or a percentage (`@1%`). `METRIC` is any metric above, with or without its `quality:` prefix. A structure works too in Python and C++: `{"metric", "min", "max", "max_fraction"}` / `QualityThreshold`.

**Which cells.** A cell where the metric does not apply (NaN — `min_angle` on a quad, say) is not evaluated; a threshold that applies to no cell at all passes vacuously and says so in a warning. Inverted and degenerate cells are gated by count: `max_inverted` and `max_degenerate` default to 0, and a negative limit disables the check.

**Each check reports** its `name`, `metric`, bounds, how many cells it `evaluated` and how many are `violations` (and the `fraction`), the `worst` value — the one closest to, or furthest beyond, the bound — with its global block-major `worst_cell`, and whether it `passed`.

```sh
meshioplusplus check part.vtu --require "scaled_jacobian >= 0.2" --require "aspect_ratio <= 5 @ 1%"
meshioplusplus check part.vtu --gate quality.gate --json
```

The `check` verb (both CLIs) exits **0** when every check passes, **1** when one fails and **2** when the check could not run (an unreadable file, a malformed specification) — so a CI job can tell "the mesh is bad" from "the gate is broken". It is also the pipeline step `{"Op": "QualityGate", "Require": ["scaled_jacobian >= 0.2"], "MaxInverted": 0, "MaxDegenerate": 0}`, which leaves the mesh untouched and stops the pipeline with the summary as its error, and the MCP tool `check_quality`. The flat bindings take the specification text and return the pass/fail and counts plus the summary text both CLIs print: C `mio_check_quality(mesh, spec, max_inverted, max_degenerate, &report, buf, len)`, Fortran `m%check_quality(spec, …)`, Julia `check_quality(m; require=…)`, R `mio_check_quality(m, require)`, WASM `checkQuality(mesh, require)` (which also returns every check).

## Cross-language

Available in every binding surface: Python (`compute_quality`/`attach_quality`), the C API (`mio_attach_quality` + `mio_quality_counts`), Fortran (`mesh%attach_quality()` / `mesh%quality_counts(...)`), and WASM (`attachQuality`), plus the CLI verb `meshioplusplus quality`.

## Implementation

The metrics are in `src/cpp/src/operations/quality.cpp` (over the uniform mesh API, with the per-cell loop parallelised) with reusable vector math in `src/cpp/include/meshioplusplus/detail/geometry.hpp`. A pure-numpy twin (`meshioplusplus._quality._compute_quality_py`) mirrors the formulas and is checked against the C++ core in `tests/python/test_quality.py`.
