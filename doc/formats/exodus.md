# Exodus II (`.e`, `.exo`, `.ex2`)

The [Exodus II](https://nschloe.github.io/meshio/exodus.pdf) format, stored
in netCDF using its classic variable/dimension conventions.

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

Global attrs: `title`, `version=5.1f`, `api_version=5.1f`,
`floating_point_word_size=8`. Dimensions: `num_nodes`, `num_dim`, `num_elem`,
`num_el_blk`, `num_node_sets`, `len_string=33`, `len_line=81`, `four=4`,
`time_step` (unlimited). Key variables:

- `time_whole(time_step)` — a dummy single `0.0` timestep is always written.
- `coor_names(num_dim, len_string)` — single-character `"X"`, `"Y"`, `"Z"`.
- `coord(num_dim, num_nodes)` — transposed relative to meshio's `(n, dim)`
  layout — or, alternatively, separate `coordx`/`coordy`/`coordz(num_nodes)`
  variables (both styles are accepted on read).
- `eb_prop1(num_el_blk)` — arbitrary distinct per-block ids (their exact
  values don't matter, only that they differ, per a ParaView requirement
  noted in the source).
- `connect{k}(num_el_in_blk{k}, num_nod_per_el{k})` for `k = 1..num_el_blk`,
  with a text `elem_type` attribute; 1-based node indices.
- `name_nod_var`/`vals_nod_var{k}` — point-data names and values (only the
  **first** timestep is read/written; extra timesteps trigger a warning).
- `name_elem_var`/`vals_elem_var{idx}[eb{block}]` — cell data indexed by
  `(variable index, element block)`, later concatenated across blocks in
  block order and re-split by target cell-block size.
- `ns_names`/`node_ns{k}` — node (point) sets, 1-based.
- `info_records`/`qa_records` — free-text info strings.

## Cell types

A large type table; representative entries:

| Exodus | meshio | Exodus | meshio |
|---|---|---|---|
| `SPHERE` | `vertex` | `TETRA`, `TET4` | `tetra4` (note: **not** `"tetra"`) |
| `BEAM`, `BEAM2`, `BAR2` | `line` | `TETRA4` | `tetra4` |
| `BEAM3` | `line3` | `TETRA8` | `tetra8` |
| `SHELL4`, `QUAD4` | `quad` | `TETRA10` | `tetra10` |
| `SHELL8`, `QUAD8` | `quad8` | `TETRA14` | `tetra14` |
| `SHELL9`, `QUAD9` | `quad9` | `PYRAMID` | `pyramid` |
| `HEX8`, `HEXAHEDRON` | `hexahedron` | `WEDGE` | `wedge` |
| `HEX20` | `hexahedron20` | `TRI3`, `TRIANGLE` | `triangle` |
| `HEX27` | `hexahedron27` | `TRI6` | `triangle6` |

The write-side reverse map picks one canonical Exodus name per meshio type
(e.g. `hexahedron → HEX8`, `tetra → TETRA`, `tetra4 → TET4` — a distinct
entry from plain `tetra`).

## Data mapping

- `point_data` — arbitrary names, with automatic recombination: names ending
  in `X`/`Y`/`Z` are checked for sibling `Y`/`Z` names and, if found, stacked
  into a 3-component vector; names ending `_R`/`_Z` are checked for a
  sibling and stacked into a 2-component vector.
- `cell_data` — arbitrary names, split per cell block by node count.
- `point_sets` — node sets (Exodus's own "node set" concept), 1-based in
  file.
- `mesh.info` — free-text strings from `info_records`/`qa_records`.

## Quirks & limitations

- The point-data name recombination (`categorize()`) has a **deliberately
  preserved quirk**: the check for a paired variable uses Python truthiness
  on the found array index, so an index of exactly `0` is treated the same
  as "not found". This is a latent edge case in the reference implementation
  that the C++ port reproduces on purpose, rather than silently fixing —
  changing it would make the two implementations disagree on some inputs.
- Only the **first timestep** is ever read (a warning is emitted if more
  exist) — matches a known ParaView writer limitation referenced in the
  source.
- `info_records`/`qa_records`/`ns_names`/`node_ns*` (info strings and node
  sets) are **explicitly unsupported by the C++ reader**, which throws and
  routes to Python whenever any of them are present in the file.
- The C++ writer does not support `mesh.point_sets` at all; the shim only
  attempts the C++ write path when `point_sets` is empty.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_NETCDF`,
  otherwise through `netCDF4`.
- No reference fixture exists under `tests/meshes/exodus/`; tests round-trip
  synthetic meshes, including one exercising `point_sets` (which always
  forces the Python path).
