# Data operations

meshio++ carries three kinds of data alongside a mesh's geometry: `point_data` (one row per point), `cell_data` (one array *per cell block*, one row per cell) and `field_data` (mesh-global). The **data operations** act on those arrays rather than on the geometry, which none of them ever modifies — points, connectivity, block order and block types come through bit-identical.

| Operation | What it does | Page |
| --- | --- | --- |
| `data_manage` / `data_drop` / `data_keep` / `data_rename` | Rewrite which arrays a mesh carries and under what names | [Array management](/data_manage) |
| `point_data_to_cell_data` / `cell_data_to_point_data` | Move data between locations by averaging | [Location averaging](/data_average) |
| `data_calc` | Derive a new array from an elementwise expression | [Expressions](/data_calc) |
| `data_condition` | Clamp, normalize or standardize values | [Value conditioning](/data_condition) |
| `data_info` | Read-only per-array summary | [Data summary](/data_info) |

They are mesh **operations** (like [quality metrics](/mesh_quality) or [geometric statistics](/stats)), not file formats, and use only standard C++/numpy, so they run under every mesh backend and are reachable from every binding surface. All nine CLI verbs live under the `meshioplusplus data` group — see the [CLI reference](/cli#meshioplusplus-data).

::: tip Near neighbours that are *not* one of these
[`gradient`](/gradient), its companion [`hessian`](/hessian) one order further, and [`data_integrate`](/field_integration) also consume and produce data arrays, and are reachable as `meshioplusplus data gradient` / `meshioplusplus data hessian` / `meshioplusplus data integrate` — but they **read geometry and topology** (face areas, cell volumes, cell adjacency), so they are mesh operations and none of the rules on this page (the NaN policy, the dtype table) describe them. They are grouped under `data` in the CLI because that is where a user looks for them.
:::

## Non-finite values {#nan-policy}

Every data operation follows one rule for `NaN` and `±inf`:

- They are **always excluded from every reduction** — min, max, mean, standard deviation, and the accumulators used by point/cell averaging. An array that is entirely non-finite therefore reduces to `NaN`.
- What reaches the **output** is chosen by `nan_policy`:
  - `"ignore"` (the default) passes the value through unchanged;
  - `"replace"` substitutes `nan_replacement`;
  - `"fail"` raises, naming the array and the first offending index.
- [`data_info`](/data_info) never raises — it *counts* non-finite values instead, which is the point of the report.
- Division by zero inside [`data_calc`](/data_calc) produces an IEEE infinity or `NaN` and is deliberately **not** an error.

## Locations

![Point data, cell data and field data beside the geometry they never modify, and the five data operations that act on them](/diagrams/data_locations.svg)

Everywhere a location is named, these spellings are accepted:

| Canonical | Aliases |
| --- | --- |
| `point_data` | `point`, `node` |
| `cell_data` | `cell`, `element` |
| `field_data` | `field` |

## Multiple cell blocks

`cell_data` is stored as one array per cell block. Every operation that writes `cell_data` produces exactly `len(mesh.cells)` arrays, and every operation that reads a named `cell_data` array validates that its block count matches the mesh, raising a clear error otherwise. [`data_condition`](/data_condition) computes its statistics **jointly across all blocks** before applying a single transform to each, and [`data_calc`](/data_calc) is evaluated once per block.

## Data types

| Operation | dtype rule |
| --- | --- |
| `data_manage` | Byte-identical passthrough; nothing is reinterpreted. |
| `point_data_to_cell_data` / `cell_data_to_point_data` | **Always `float64`** — the mean of an integer field is not an integer. |
| `data_calc` | Arithmetic is always performed in `double`; the result is stored as `float64`. |
| `data_condition` | `clamp` preserves the input dtype (unless `preserve_dtype=False`); `normalize` and `standardize` always produce `float64`. |
| `data_info` | Read-only; reports the dtype as stored. |

::: tip
The `NATIVE` and `KRATOS` mesh backends canonicalize integer widths on ingest (`int32` → `int64`), so tests should assert the dtype *kind* rather than an exact width. Only the Python path, which is pinned to the `MESHIO` backend, sees the original width.
:::

## Point and cell sets

`point_sets` and `cell_sets` index geometry, which these operations never touch, so they pass through unchanged. Like every other operation, they are handled in the Python shim and never reach the C++ core — so they are **not** carried by the C API, Fortran or WebAssembly surfaces.
