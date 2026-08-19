# Mesh comparison (diff)

`meshioplusplus.diff(mesh_a, mesh_b, ...)` compares two meshes and reports whether they are equivalent **within a tolerance**, plus a structured description of any differences. `meshioplusplus.meshes_equal(a, b, ...)` is the convenient boolean wrapper for test suites. It is a mesh **operation** (like [reordering](/reorder), [surface extraction](/extract_surface) and [quality metrics](/mesh_quality)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![Per-point coordinate error from `diff`, visualised on the mesh](/images/diff_error.png)

```python
import meshioplusplus

a = meshioplusplus.read("reference.vtu")
b = meshioplusplus.read("candidate.vtu")

# quick boolean, ideal in a regression test
assert meshioplusplus.meshes_equal(a, b, atol=1e-8, rtol=1e-6)

# full structured report
report = meshioplusplus.diff(a, b, atol=1e-8, rtol=1e-6)
print(report["verdict"])                       # "identical" / "equal within tolerance" / "different"
print(report["points"]["max_abs_error"])
```

## Verdicts

| verdict | meaning |
|---|---|
| `"identical"` | structurally equal and every value bitwise-equal |
| `"equal within tolerance"` | structurally equal; values differ but stay within tolerance |
| `"different"` | a structural difference, or a value beyond tolerance |

The verdict is the most severe outcome over every compared section.

## Tolerance

Values are compared elementwise with the rule

```
abs_err <= atol + rtol * |expected|
```

where `expected` is the value from `mesh_a`. Defaults are `atol=1e-12` and `rtol=1e-9`. For a round-trip through an ASCII format, pass a looser `atol` (e.g. `1e-6`).

## What gets compared

- **Points** — same count and dimension? per-coordinate errors within tolerance (max abs / max rel error, flat index of the worst offender, count exceeding).
- **Cells** — same cell blocks (types + counts)? same connectivity per block? The first few mismatching cell indices and the total mismatch count are reported per block.
- **`point_data` / `cell_data` / `field_data`** — same set of keys? For shared keys, same shape and values within tolerance (per-array max abs/rel error and mismatch counts); keys present in only one mesh are listed.
- **`point_sets` / `cell_sets`** — added / removed / changed named sets, under a `"sets"` section. These live only on the Python `Mesh`, so they are compared in the Python layer and are **not** seen by `meshes_equal` (nor by the C / Fortran / WebAssembly bindings).

## Ordered vs. unordered

By default the comparison is **ordered**: points and cells are compared index-by-index — the fast, O(N), parallel path for the common regression-test case where only numerical values may have drifted.

With `unordered=True` the points of the two meshes are matched by **spatial proximity within tolerance** using a bucket-grid coordinate hash (no O(N²) scan), then connectivity is compared up to that node correspondence. This makes two meshes with the same geometry but a **shuffled node order** compare equal. It is conservative: if a robust one-to-one correspondence cannot be established it sets `correspondence_failed` and reports `"different"` rather than guessing.

```python
# same geometry + connectivity, nodes numbered differently -> equal
assert meshioplusplus.meshes_equal(a, b, unordered=True)
```

## Report structure

`diff` returns a nested `dict`:

- `verdict`, `unordered`, `correspondence_failed`, `messages`
- `point_count_mismatch`, `num_points_a/b`, `points` (an *array diff*)
- `block_count_mismatch`, `num_blocks_a/b`, `blocks` (a list of *block diffs*)
- `point_data`, `cell_data`, `field_data` (each a *data diff*)
- `sets` (`point_sets` / `cell_sets`, each with `only_in_a` / `only_in_b` / `changed`)

An *array diff* has `shape_mismatch`, `size_a/b`, `max_abs_error`, `max_rel_error`, `worst_index` (flat index), `num_exceeding`, `exact`. A *block diff* has `type_a/b`, `count_a/b`, `type_mismatch`, `count_mismatch`, `conn_mismatch_count`, `first_mismatches`. A *data diff* has `only_in_a`, `only_in_b`, and `shared` (a list of array diffs).

## CLI

```bash
meshioplusplus diff a.vtu b.vtu                      # human report, exit 0 if equal
meshioplusplus diff a.vtu b.vtu --atol 1e-8 --rtol 1e-6
meshioplusplus diff a.msh b.vtu --unordered          # tolerant node matching
meshioplusplus diff a.vtu b.vtu --exact              # only bitwise-identical passes
meshioplusplus diff a.vtu b.vtu --quiet              # no output, only exit code
```

The command sets a **nonzero exit code when the meshes differ** and **zero when they are equal** (identical, or — unless `--exact` — equal within tolerance), so it drops straight into CI, shell scripts and Makefiles:

```bash
meshioplusplus diff expected.vtu actual.vtu --atol 1e-8 || echo "regression!"
```

See the [CLI reference](/cli).

## Other languages

The operation is exposed across every binding surface:

- **C API** — `mio_diff(a, b, atol, rtol, unordered, &result)` fills a `mio_diff_result` (read the verdict with `mio_diff_result_verdict`, the point summary with `mio_diff_result_point_summary`, per-block flags with `mio_diff_result_block`; free with `mio_diff_result_free`); the convenience `mio_meshes_equal(a, b, atol, rtol, unordered, &equal)`. See the [C API reference](/c_api).
- **Fortran** — `mesh%equals(other, atol=..., rtol=..., unordered=...)` (logical) and `mesh%diff(other, ...)` (verdict `0`/`1`/`2`). See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `diff(meshA, meshB, atol, rtol, unordered)` returns a report object; `meshesEqual(meshA, meshB, atol, rtol, unordered)`. See the [WebAssembly reference](/wasm).
