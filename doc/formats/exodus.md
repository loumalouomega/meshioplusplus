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

`write` takes no keyword arguments. `read` takes `time_step` (see [Time steps](#time-steps)).

## File structure

Global attrs: `title`, `version=5.1f`, `api_version=5.1f`, `floating_point_word_size=8`. Dimensions: `num_nodes`, `num_dim`, `num_elem`, `num_el_blk`, `num_node_sets`, `len_string=33`, `len_line=81`, `four=4`, `time_step` (unlimited). Key variables:

- `time_whole(time_step)` — one step. The writer takes its value from `field_data["exodus:time"]` (v9.9.0) and both readers set that key to the time of the step they returned, so one frame of a transient solve keeps its label. Absent, it is `0.0`.
- `coor_names(num_dim, len_string)` — single-character `"X"`, `"Y"`, `"Z"`.
- `coord(num_dim, num_nodes)` — transposed relative to meshio++'s `(n, dim)` layout — or, alternatively, separate `coordx`/`coordy`/`coordz(num_nodes)` variables (both styles are accepted on read).
- `eb_prop1(num_el_blk)` — arbitrary distinct per-block ids (their exact values don't matter, only that they differ, per a ParaView requirement noted in the source).
- `connect{k}(num_el_in_blk{k}, num_nod_per_el{k})` for `k = 1..num_el_blk`, with a text `elem_type` attribute; 1-based node indices. The attribute is trimmed of trailing NULs and spaces before the type lookup — see [Quirks](#quirks-limitations).
- `attrib{k}(num_el_in_blk{k}, num_att_in_blk{k})` / `attrib_name{k}(num_att_in_blk{k}, len_string)` — per-element attributes, read and written as `cell_data` under the `exodus:attr:` prefix (see [Element attributes](#element-attributes)).
- `name_nod_var`/`vals_nod_var{k}` — point-data names and values, sliced at the requested time step (see [Time steps](#time-steps)); the writer emits one step.
- `name_elem_var`/`elem_var_tab`/`vals_elem_var{j}eb{k}` — ordinary (non-attribute) `cell_data`, one variable × block array each, read and written since v9.9.0 (see [Data mapping](#data-mapping)).
- `name_elem_var`/`vals_elem_var{idx}[eb{block}]` — cell data indexed by `(variable index, element block)`, later concatenated across blocks in block order and re-split by target cell-block size.
- `eb_names(num_el_blk, len_string)` — element-block names, read into `Cell` regions and, since v9.9.0, written back from them (see [Named regions](#named-regions)).
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

`SPHERE` is the one-node element particle codes use — peridynamics solvers such as [PeriLab](https://github.com/PeriHub/PeriLab.jl) write whole meshes of them — and maps to meshio++'s `vertex`. Its radius is not part of the connectivity; it lives in a per-element attribute (below).

## Data mapping

- `point_data` — arbitrary names, with automatic recombination: names ending in `X`/`Y`/`Z` are checked for sibling `Y`/`Z` names and, if found, stacked into a 3-component vector; names ending `_R`/`_Z` are checked for a sibling and stacked into a 2-component vector.
- `cell_data` — arbitrary names, split per cell block by node count. Written as element variables since v9.9.0 (`vals_elem_var{j}eb{k}`); before that, dropped.
- `cell_data["exodus:attr:<NAME>"]` — per-element attributes (see [Element attributes](#element-attributes)).
- `regions` — element blocks, node sets and side sets (see [Named regions](#named-regions)). `point_sets` is a compat view over the `Point` regions.
- `mesh.info` — free-text strings from `info_records`/`qa_records`.

## Element attributes

Exodus stores a fixed number of floating-point **attributes** per element of a block, in `attrib{k}` and named by `attrib_name{k}`. They are the standard home for a value the connectivity cannot express: a `SPHERE`/`CIRCLE` element's **radius**, a beam's cross-section area, a shell's thickness. Since v9.3.0 they read and write as ordinary `cell_data`, under the `exodus:attr:` prefix:

```python
mesh = meshioplusplus.read("particles.exo")
radius = mesh.cell_data["exodus:attr:RADIUS"]   # one array per cell block
```

The prefix is what makes the mapping unambiguous in both directions. On read it keeps an attribute from colliding with a same-named element *variable* (`name_elem_var`), which is a genuinely different concept — attributes are constant in time, element variables are per-time-step. On write it is the only signal telling the writer which `cell_data` arrays belong in `attrib{k}`; everything else in `cell_data` is left alone.

Three rules follow from `cell_data` having exactly one array per cell block while Exodus attributes are per block:

- Values are always **float64** on read, whatever the on-disk type — Exodus attributes are floating point by definition, and the NaN below needs somewhere to live.
- A block the file gives no such attribute is filled with **NaN**. There is no "absent" to report in a `cell_data` array, and the name is still information.
- On write, a block whose values are **all** non-finite is left out again. That is exactly the NaN the reader fills in, so a file where only some blocks carry an attribute round-trips instead of gaining NaN attributes; the only thing lost is a genuinely all-NaN attribute, which carries no information anyway.

An attribute is one value per element, so a multi-component array under this prefix is a `WriteError` naming it rather than a silent flatten. Unnamed attributes (a blank `attrib_name{k}`, which SEACAS writes often enough) are named by their 1-based column, `exodus:attr:attribute1` and so on, so two blocks agree on which one is "the first".

## Named regions

Since v8.6.0 the reader maps all three of Exodus's grouping concepts onto [regions](../regions.md):

| Exodus | Region kind | Name from | `tag` from |
|---|---|---|---|
| `connect{k}` (element block) | `cell` | `eb_names`, else `"Block <id>"` | `eb_prop1` |
| `node_ns{k}` (node set) | `point` | `ns_names`, else `"Nodeset <id>"` | `ns_prop1` |
| `elem_ss{k}`+`side_ss{k}` (side set) | `side` | `ss_names`, else `"Sideset <id>"` | `ss_prop1` |

Cell entries are **global block-major** cell indices; side entries are `(global cell, local facet)` pairs. Two element blocks of the *same* element type stay distinguishable, since the name and tag come from per-block arrays — which is exactly what a materials assignment depends on.

**Exodus side numbering is not meshio++ facet numbering.** Exodus orders an element's sides its own way, so the facet column is remapped through `exo_face_index` (the twin of `abq_face_index`) rather than stored raw; a gtest pins every entry against `detail/cell_faces.hpp` by node set. An unmappable `(cell type, side)` pair is skipped rather than stored pointing at the wrong face.

**Element blocks round-trip since v9.9.0**; node sets and side sets are still read-only. The writer recovers `eb_names` by inverting the read path: a `Cell` region whose canonical entries are exactly the contiguous global range `[base_k, base_k+1)` *is* block `k`'s name. Entries are sorted and de-duplicated by `AddRegion`, so "covers exactly this block" is a first/last check rather than a set comparison. `eb_names` is written only when at least one block is actually named, so a region-less mesh's bytes are unchanged — and a block with no matching region gets `""`, which is what SEACAS itself writes.

The writer still emits no node sets and no side sets, so a `Point` or `Side` region written here does not come back.

## Time steps

Exodus stores every data array with a leading `time_step` dimension. Before v8.6.0 the reader always took the first step and warned "Skipping some time data"; now `ReadOptions::mTimeStep` selects one:

- `0` (the default) is the first step, so existing behaviour is unchanged.
- Negative values count from the end; `-1` is the last step, which is what "the final state of the simulation" means without knowing the count up front.
- Out of range is an **error naming the available count**, never a silent clamp — quietly returning step 0 for a request of step 7 is the failure this option exists to remove.

`read_metadata(...)["time_values"]` reports the recorded times (from `time_whole`), so a step request is checkable before it is issued. Exodus is registered as an options-aware reader, so `reader_supports_options("exodus")` is true and the option reaches every binding: `time_step=` in Python, `mio_read_opts.time_step` in C, `m%read(..., time_step=)` in Fortran, `ReadOptions(time_step=)` in Julia, `mio_read(time_step=)` in R, `readMeshSelective(path, {timeStep})` in WASM, and `--time-step=N` in both CLIs.

## Quirks & limitations

- **A NUL-terminated `elem_type` used to fail the read.** netCDF text attributes carry an explicit length, and NetCDF.jl — which is what [PeriLab](https://github.com/PeriHub/PeriLab.jl) and other Julia solvers write Exodus with — counts the C string's terminating NUL as part of it. So a `SPHERE` block arrives as the 7 characters `"SPHERE\0"`, which matched no key in the C++ reader's type table: the read failed with `Exodus: unknown element type SPHERE`, the NUL invisible in the message because `std::runtime_error::what()` is a `const char*` that stops at it. `netCDF4` strips the NUL on the way in, so the Python reference never saw this and the shim's silent fallback hid it everywhere **except WASM**, which has no fallback — which is how it surfaced as [VSCode-MDPA-Preview#63](https://github.com/loumalouomega/VSCode-MDPA-Preview/issues/63) rather than as a Python bug. Fixed in v9.3.0: both readers now trim trailing NULs and spaces before the lookup. The same normalization covers fixed-width writers that pad with spaces.
- The point-data name recombination (`categorize()`) has a **deliberately preserved quirk**: the check for a paired variable uses Python truthiness on the found array index, so an index of exactly `0` is treated the same as "not found". This is a latent edge case in the reference implementation that the C++ port reproduces on purpose, rather than silently fixing — changing it would make the two implementations disagree on some inputs.
- `qa_records`/`info_records` are **provenance strings**, and `NDArray` has no string dtype, so they cannot ride on the mesh. They travel in an `ExodusInfo` side-channel struct (the `MedInfo`/`OpenFoamInfo` pattern) that the pybind binding attaches as `mesh.info`. The **flat bindings** (C, Fortran, Julia, R, WASM) construct one and drop it — the same documented gap `MedInfo` already has.
- Before v8.6.0 the C++ reader **threw** on `qa_records`, `info_records`, `ns_names` and `node_ns*`, routing the whole file to Python. Since every file SEACAS, Cubit or Sierra writes carries `qa_records`, that made Exodus entirely unreadable from **WASM**, which has no Python fallback to defer to. Fixed; see [Time steps](#time-steps) and [Named regions](#named-regions) for what those variables now produce.
- The C++ writer does not support `mesh.point_sets` at all; the shim only attempts the C++ write path when `point_sets` is empty. It also emits no node sets and no side sets, which is why only element-block regions round-trip.
- **Ordinary `cell_data` round-trips since v9.9.0** as element variables (`name_elem_var`, an all-ones `elem_var_tab` truth table, and one `vals_elem_var{j}eb{k}` per variable × block). Before that neither writer emitted any of it, so every `cell_data` array except the `exodus:attr:`-prefixed ones was silently dropped while `point_data` round-tripped. The prefix stays explicit for the same reason it always was: an attribute is a *different* Exodus concept (constant in time, one column per element), so "write every cell_data as an attribute" was never the right rule. Trailing dimensions become extra netCDF dimensions exactly as the nodal path already does, so a vector cell field survives — standard Exodus element variables are scalar per element, making a k>1 array a meshio++ extension of the same kind the nodal path already is.
- **Fixed in v9.9.0: a heap buffer overflow in the C++ reader.** Assembling a `cell_data` array allocated a scalar `{total}` buffer and then `memcpy`'d each block's full `Nbytes()` into it, so a multi-component element variable wrote `n*k` bytes into `n` bytes of space. Pre-existing, but unreachable until this release's writer started emitting element variables, and not reachable from any real SEACAS file — none carries a multi-component element variable, which is exactly why a format's own writer is not a sufficient test oracle for its reader.
- **One step per write.** A `Mesh` is one state, so the writer emits a single `time_step`; `field_data["exodus:time"]` labels it. A genuinely multi-step Exodus writer is a stateful object of the shape `XdmfTimeSeriesWriter` has and is a follow-up, not something a `(path, mesh)` writer can express.
- A file with `num_dim = 2` (`coordx`/`coordy`, no `coordz`) reads back with **3-component** points whose `z` column is zero, on both paths. That is upstream meshio's behaviour and is kept.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_NETCDF`, otherwise through `netCDF4`.
- Read/written through the C++ core when built with netCDF, otherwise through `netCDF4`.
- The regression file itself is committed at `tests/python/meshes/exodus/DCBmodel_PD_solid.e` (Git-LFS): a real PeriLab double-cantilever-beam run — 504 one-node `SPHERE` particles in four blocks, 2-D coordinates, nine nodal fields, ten time steps, and a damage field that goes from 0 to 0.48 across them, so a reader silently pinned to step 0 fails rather than merely differs. It is BSD-3-Clause, redistributed unmodified with its notice and licence text; see that directory's `README.md`. A hand-authored fixture can reproduce the *shape*, but only a file written by that toolchain has the NUL-terminated `elem_type` naturally, which is why both exist.
- Tests round-trip synthetic meshes, **plus** two hand-authored fixtures built at test time by `tests/python/exodus_fixture.py`: a SEACAS/Cubit-shaped one (`qa_records`, two same-type element blocks, two node sets, one side set, three time steps) and a peridynamics-shaped one (2-D `coordx`/`coordy`, NUL-terminated `elem_type`, one-node `SPHERE` blocks, a `RADIUS` attribute on one block only). They are hand-authored because meshio++'s own writer emits none of that, so a round-trip test could not have caught either the `qa_records` defect or the `elem_type` one. The NUL-terminated attribute cannot be written through `netCDF4` at all (it strips the NUL), so the fixture patches the classic-format bytes: an attribute is `[nc_type][nelems][value padded to 4 bytes]`, and since `len("SPHERE")` is not a multiple of 4 the padded size is the same for 6 and 7 characters, so bumping the recorded count rewrites four bytes in place and shifts nothing. The gtest counterpart (`tests/cpp/test_netcdf_formats.cpp`) builds its fixture through the netCDF C API instead, where the length is simply a parameter.
