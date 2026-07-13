# Exodus II (`.e`, `.exo`, `.ex2`)

The [Exodus II](https://nschloe.github.io/meshio/exodus.pdf) format, stored in
netCDF.

| | |
|---|---|
| **Format name** | `exodus` |
| **Extensions** | `.e`, `.exo`, `.ex2` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `netCDF4` (or a C++ build with netCDF) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.exo")
meshio.exodus.write("out.exo", mesh)
```

`write` takes no keyword arguments.

## File structure

netCDF dimensions `num_nodes`/`num_dim`/`num_elem`/`num_el_blk`, a transposed
`coord` variable, `connect{k}` element blocks (1-based, with an `elem_type`
attribute), point data via `name_nod_var`/`vals_nod_var{k}`, cell data via
`vals_elem_var{i}eb{k}`, and node sets via `ns_*`.

## Cell types

A large Exodus ↔ meshio type map including second-order elements.

## Data mapping

- Point data (with `X/Y/Z` and `R/Z` component recombination) → `point_data`;
  per-block cell data → `cell_data`.
- Node sets → `point_sets`.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_NETCDF`,
  otherwise through `netCDF4`. Node sets and info/QA records fall back to the
  Python implementation. Only the first time-step is read.
