# Data summary (`data_info`)

`meshioplusplus.data_info(mesh)` reports a per-array summary of everything a mesh's data maps carry, without modifying it. It is the **data** view that complements the topological [`info`](/cli#meshioplusplus-info) verb and the geometric [`stats`](/stats): `info` tells you about points and cells, `stats` about the shape they occupy, and `data_info` about the fields attached to them. It is a [data operation](/data_operations), not a file format.

```python
import meshioplusplus as mp

mesh = mp.read("part.vtu")
for a in mp.data_info(mesh):
    print(a["location"], a["name"], a["dtype"],
          a["num_entries"], a["num_components"],
          a["min"], a["max"], a["mean"], a["num_nan"])
```

## Returned fields

`data_info` returns a list of dicts, ordered `point_data`, then `cell_data`, then `field_data`, each group sorted by name. Each dict carries:

| Key | Meaning |
| --- | --- |
| `location` | `"point_data"`, `"cell_data"` or `"field_data"` |
| `name` | the array name |
| `dtype` | the dtype as stored, e.g. `"f8"`, `"i4"` |
| `shape` | the shape as stored (for `cell_data`, block 0's) |
| `num_blocks` | `cell_data`: the number of cell blocks; 1 otherwise |
| `num_entries` | rows — points, cells summed over all blocks, or the field length |
| `num_components` | product of the trailing dimensions (1 for a scalar) |
| `num_values` | `num_entries × num_components` |
| `min` / `max` / `mean` | over the finite values; `NaN` when there are none |
| `min_per_component` / `max_per_component` / `mean_per_component` | the same, per component |
| `num_nan` / `num_inf` / `num_finite` | value counts |
| `inconsistent_blocks` | `True` if a `cell_data` array's blocks disagree in component count |

Unlike the other data operations, `data_info` **never raises** on non-finite values — it counts them, which is the point of the report. They are still excluded from `min`/`max`/`mean`, per the [NaN policy](/data_operations#nan-policy). An inconsistent `cell_data` array is likewise reported rather than raised, so the rest of the summary stays usable.

## CLI

```bash
meshioplusplus data info mesh.vtu           # human-readable table
meshioplusplus data info mesh.vtu --json    # machine-readable JSON
```

```
<meshio++ data summary>
  location    name                 dtype  comp  entries          min          max         mean   nan   inf
  --------------------------------------------------------------------------------------------------------
  point_data  T                    f8        1        6            0           12            6     0     0
  point_data  velocity             f8        3        6            0            2          0.5     0     0
  cell_data   mat                  f8        1        3            1            3            2     0     0
```

The Python CLI and the native binary produce the same table. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

- **C API** — `mio_data_info_create(mesh)` returns an opaque `mio_data_info*`; read it with `mio_data_info_count`, `mio_data_info_name` (caller-buffer protocol), `mio_data_info_entry` (fills a `mio_data_array_info` struct) and `mio_data_info_component`, then release it with `mio_data_info_free`. See the [C API reference](/c_api).
- **Fortran** — `arrays = m%data_info(keys=names)` returns an array of `mio_data_array_info` plus, optionally, the matching names. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataInfo(mesh)` returns an array of objects with camelCase keys (`numEntries`, `numComponents`, `minPerComponent`, …). See the [WebAssembly reference](/wasm).
