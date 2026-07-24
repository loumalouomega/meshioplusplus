# Exodus II (`.e`, `.exo`, `.ex2`)

The [Exodus II](https://nschloe.github.io/meshio/exodus.pdf) format, stored in netCDF using its classic variable/dimension conventions.

| | |
|---|---|
| **Format name** | `exodus` |
| **Extensions** | `.e`, `.exo`, `.ex2` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `netCDF4` (or a C++ build with netCDF) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.exo")
meshioplusplus.exodus.write("out.exo", mesh)

# A multi-step file: 0 (the default) is the first step, negative counts
# from the end. Out of range raises, naming the available count.
last = meshioplusplus.read("run.exo", time_step=-1)

# How many steps are there, without reading the arrays?
meshioplusplus.read_metadata("run.exo")["time_values"]  # e.g. [0.0, 0.5, 1.0]
```

`write` takes no keyword arguments. `read` takes `time_step` (see
[Time steps](#time-steps)).

## File structure

Global attrs: `title`, `version=5.1f`, `api_version=5.1f`, `floating_point_word_size=8`. Dimensions: `num_nodes`, `num_dim`, `num_elem`, `num_el_blk`, `num_node_sets`, `len_string=33`, `len_line=81`, `four=4`, `time_step` (unlimited). Key variables:

- `time_whole(time_step)` — a dummy single `0.0` timestep is always written.
- `coor_names(num_dim, len_string)` — single-character `"X"`, `"Y"`, `"Z"`.
- `coord(num_dim, num_nodes)` — transposed relative to meshio++'s `(n, dim)` layout — or, alternatively, separate `coordx`/`coordy`/`coordz(num_nodes)` variables (both styles are accepted on read).
- `eb_prop1(num_el_blk)` — arbitrary distinct per-block ids (their exact values don't matter, only that they differ, per a ParaView requirement noted in the source).
- `connect{k}(num_el_in_blk{k}, num_nod_per_el{k})` for `k = 1..num_el_blk`, with a text `elem_type` attribute; 1-based node indices.
- `name_nod_var`/`vals_nod_var{k}` — point-data names and values, sliced at the requested time step (see [Time steps](#time-steps)); the writer emits one dummy step.
- `name_elem_var`/`vals_elem_var{idx}[eb{block}]` — cell data indexed by `(variable index, element block)`, later concatenated across blocks in block order and re-split by target cell-block size.
- `eb_names(num_el_blk, len_string)` — element-block names, read into `Cell` regions.
- `ns_names`/`ns_prop1`/`node_ns{k}` — node sets, 1-based; read into `Point` regions.
- `ss_names`/`ss_prop1`/`elem_ss{k}`/`side_ss{k}` — side sets; read into `Side` regions.
- `info_records`/`qa_records` — free-text info strings, read into `mesh.info`.

## Cell types

A large type table; representative entries:

| Exodus | meshio++ | Exodus | meshio++ |
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

The write-side reverse map picks one canonical Exodus name per meshio++ type (e.g. `hexahedron → HEX8`, `tetra → TETRA`, `tetra4 → TET4` — a distinct entry from plain `tetra`).

## Data mapping

- `point_data` — arbitrary names, with automatic recombination: names ending in `X`/`Y`/`Z` are checked for sibling `Y`/`Z` names and, if found, stacked into a 3-component vector; names ending `_R`/`_Z` are checked for a sibling and stacked into a 2-component vector.
- `cell_data` — arbitrary names, split per cell block by node count.
- `regions` — element blocks, node sets and side sets (see [Named regions](#named-regions)). `point_sets` is a compat view over the `Point` regions.
- `mesh.info` — free-text strings from `info_records`/`qa_records`.

## Named regions

Since v8.6.0 the reader maps all three of Exodus's grouping concepts onto
[regions](../regions.md):

| Exodus | Region kind | Name from | `tag` from |
|---|---|---|---|
| `connect{k}` (element block) | `cell` | `eb_names`, else `"Block <id>"` | `eb_prop1` |
| `node_ns{k}` (node set) | `point` | `ns_names`, else `"Nodeset <id>"` | `ns_prop1` |
| `elem_ss{k}`+`side_ss{k}` (side set) | `side` | `ss_names`, else `"Sideset <id>"` | `ss_prop1` |

Cell entries are **global block-major** cell indices; side entries are
`(global cell, local facet)` pairs. Two element blocks of the *same* element type
stay distinguishable, since the name and tag come from per-block arrays — which
is exactly what a materials assignment depends on.

**Exodus side numbering is not meshio++ facet numbering.** Exodus orders an
element's sides its own way, so the facet column is remapped through
`exo_face_index` (the twin of `abq_face_index`) rather than stored raw; a gtest
pins every entry against `detail/cell_faces.hpp` by node set. An unmappable
`(cell type, side)` pair is skipped rather than stored pointing at the wrong
face.

**Reading only.** The writer still emits no `eb_names` and no side sets, so a
region written here would not come back — Exodus is a region *source*, not a
round-trip target. Recorded as such in `tests/python/test_region_roundtrip.py`'s
`READ_ONLY_REGIONS`.

## Time steps

Exodus stores every data array with a leading `time_step` dimension. Before
v8.6.0 the reader always took the first step and warned "Skipping some time
data"; now `ReadOptions::mTimeStep` selects one:

- `0` (the default) is the first step, so existing behaviour is unchanged.
- Negative values count from the end; `-1` is the last step, which is what
  "the final state of the simulation" means without knowing the count up front.
- Out of range is an **error naming the available count**, never a silent clamp
  — quietly returning step 0 for a request of step 7 is the failure this option
  exists to remove.

`read_metadata(...)["time_values"]` reports the recorded times (from
`time_whole`), so a step request is checkable before it is issued. Exodus is
registered as an options-aware reader, so `reader_supports_options("exodus")` is
true and the option reaches every binding: `time_step=` in Python,
`mio_read_opts.time_step` in C, `m%read(..., time_step=)` in Fortran,
`ReadOptions(time_step=)` in Julia, `mio_read(time_step=)` in R,
`readMeshSelective(path, {timeStep})` in WASM, and `--time-step=N` in both CLIs.

## Quirks & limitations

- The point-data name recombination (`categorize()`) has a **deliberately preserved quirk**: the check for a paired variable uses Python truthiness on the found array index, so an index of exactly `0` is treated the same as "not found". This is a latent edge case in the reference implementation that the C++ port reproduces on purpose, rather than silently fixing — changing it would make the two implementations disagree on some inputs.
- `qa_records`/`info_records` are **provenance strings**, and `NDArray` has no string dtype, so they cannot ride on the mesh. They travel in an `ExodusInfo` side-channel struct (the `MedInfo`/`OpenFoamInfo` pattern) that the pybind binding attaches as `mesh.info`. The **flat bindings** (C, Fortran, Julia, R, WASM) construct one and drop it — the same documented gap `MedInfo` already has.
- Before v8.6.0 the C++ reader **threw** on `qa_records`, `info_records`, `ns_names` and `node_ns*`, routing the whole file to Python. Since every file SEACAS, Cubit or Sierra writes carries `qa_records`, that made Exodus entirely unreadable from **WASM**, which has no Python fallback to defer to. Fixed; see [Time steps](#time-steps) and [Named regions](#named-regions) for what those variables now produce.
- The C++ writer does not support `mesh.point_sets` at all; the shim only attempts the C++ write path when `point_sets` is empty. It also emits no `eb_names` and no side sets, which is why regions are read-only here.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_NETCDF`, otherwise through `netCDF4`.
- Read/written through the C++ core when built with netCDF, otherwise through `netCDF4`.
- Tests round-trip synthetic meshes, **plus** a hand-authored SEACAS/Cubit-shaped fixture built at test time by `tests/python/exodus_fixture.py` (`qa_records`, two same-type element blocks, two node sets, one side set, three time steps). It is hand-authored because meshio++'s own writer emits none of those, so a round-trip test could not have caught the `qa_records` defect.
