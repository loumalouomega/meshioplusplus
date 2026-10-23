# Julia

## Exodus series and pipeline Version 2

`ExodusSeries(path)` (or its do-block form) owns the fixed-grid Exodus writer. Call `write_points_cells!(series, mesh)` once and `write_data!(series, time, mesh)` for each step; `flush!`, `finalize!`, `finalized`, `num_steps` and idempotent `close` follow the XDMF lifecycle. netCDF is required; named `MeshioError`s report missing support. Geometry, sets and attributes are static; field schema is fixed by the first step. See [Exodus](./formats/exodus.md#stateful-series-writing).

Existing `run_pipeline_*` and report routines accept [Version 2](./pipeline.md#version-2-spatial-multi-mesh-steps) without new parser dependencies: auxiliary path `Inputs` and terminal `Output.Pattern` with `{key}`/`{part}`. Transient sequence schemas remain Version 1.

## Pipeline reports

`run_pipeline_file_report(path)`, `run_pipeline_json_report(text)` and `run_sequence_file_report` / `run_sequence_json_report` return structured report JSON (`steps` and `warnings`) without adding package dependencies. The native owner is always released in `finally`; existing status-only APIs are unchanged. See [pipeline reports](./pipeline.md#structured-reports-on-the-flat-abi).

## MED named meshes

`med_mesh_names(path)` enumerates file link order; `read_med_named(path, name; options=ReadOptions())` returns an owning `Mesh` (close it normally); `write_med_multi(path, meshes, names; version="4.1.0")` borrows input meshes only for the call. Names must be unique/nonempty, without `/` or `@`. HDF5 is required; ordinary nodal/element profiles expand with NaN fill. See [MED](./formats/med.md#native-named-meshes-and-profiles).

Parameterized writes support `mio.write(mesh, path; encoding="binary", codec="lzf")` for PCD compressed binary (no optional dependency; ASCII/non-PCD combinations fail), and `mio.write(mesh, "run.frd"; encoding="binary")` writes CalculiX's binary result layout. `write_gltf(mesh, path; color_by="temperature", cmap="turbo", up_axis="z", scale=0.001)` exposes all native glTF options; `component` is 1-based or `nothing` for magnitude, and unset `vmin`/`vmax` selects automatic bounds. Native reads accept XDMF2/3 and absolute DataItem references and carry Netgen names/periodic arrays; Netgen periodic arrays use `field_data` and keep file node ids 1-based, while `.vol.gz` requires native zlib.

XDMF series `write_points_cells!(series, mesh)` now stores the mesh's fixed point/cell/side regions once with the shared topology. Existing step, flush and append calls retain them without new arguments. See [shared named regions](xdmf_time_series.md#shared-named-regions).

The shared native reader supports multiple pieces without welding and appended raw/base64 arrays in `vtu`, `vtp`, `vts`, `vtr` and `vti`, with UInt32/UInt64 headers and either byte order. VTP/VTS/VTR/VTI writers remain inline; optional codecs follow the C library build. Legacy `vtk` also reads structured points, structured grids and rectilinear grids in ASCII and big-endian binary through the existing read API. See [formats](formats.md).

meshio++ ships a Julia package, `MeshioPlusPlus`, layered on the [C API](/c_api) via `ccall` — the same way the [Fortran](/fortran) module is, and aimed at the same HPC audience:

```julia
import MeshioPlusPlus as mio
using MeshioPlusPlus

m = mio.read("bracket.msh")
println(m)                     # Mesh(9231 points, 42145 cells in 3 blocks)
surf = extract_surface(m)
mio.write(surf, "bracket_surface.vtu")
```

::: warning This binding is not MIT
`bindings/julia/` is released under the **GNU General Public License, version 3 (GPL-3.0)**, not the MIT licence covering the rest of meshio++. GPL-3.0 is a **copyleft** license, not a permission-required one: anyone — including a company — may use, modify or sell it commercially with **no permission needed**, but *distributing* it or a modified version of it must be under GPL-3.0 too, with source available. Purely private/internal use that is never distributed carries no obligation.

The C API and C++ core this binding calls are unaffected and stay MIT. Calling that stable, non-GPL C ABI via `ccall`/`dlopen` at runtime is the standard "linking exception" case — it does not require the C library to also be GPL. See [`bindings/julia/LICENSE`](https://github.com/loumalouomega/meshioplusplus/blob/master/bindings/julia/LICENSE).
:::

`read(path; format="gmsh")` reads non-periodic Gmsh 2.2, 4.0 and 4.1 through the native library. The 4.0 path accepts ASCII and binary with 4- or 8-byte producer counts. Periodic files use `read_with_info` / `write_with_info`; output remains 4.1 (`gmsh`) or 2.2 (`gmsh22`). Info-less reads refuse `$Periodic` rather than lose it. See [Gmsh](formats/gmsh.md).

## Building and installing

The package binds the **installed** C library; no C++ is compiled by it.

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

Then either put `/opt/meshioplusplus/lib` on the loader path, or point the package straight at the library:

```julia
ENV["MESHIOPLUSPLUS_LIB"] = "/opt/meshioplusplus/lib/libmeshioplusplus.so"

using Pkg
Pkg.develop(path="bindings/julia/MeshioPlusPlus")
```

`MESHIOPLUSPLUS_LIB` is checked first, then the standard loader path (`Libdl.find_library`). Failing both raises an error naming these two build commands.

::: warning A build with HDF5 can fail to load from Julia on Debian/Ubuntu
If the C API was built with `-DMESHIOPLUSPLUS_WITH_HDF5=ON`, loading it from a Julia process can fail with `libcurl.so.4: version 'CURL_OPENSSL_4' not found`. This is a real Debian/Ubuntu + Julia interaction, not a meshio++ bug: `libhdf5-dev` there pulls in `libhdf5_openmpi`, which links the *system* `libcurl`, while Julia bundles its own `libcurl` — and depending on the process's library-resolution order (an IJulia kernel subprocess is more likely to hit it than an interactive `julia` invocation), the two can conflict. Build with `-DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF` if this binding is all you need from the library, as the `julia` CI job does.
:::

::: tip No registration yet, and no JLL
GPL-3.0 **is** OSI-approved, so the package is eligible for Julia's General registry — but registering it is a separate follow-up step, not done yet, so `Pkg.add("MeshioPlusPlus")` will not work until then; install by path or by URL in the meantime.

There is also deliberately no [JLL](https://docs.binarybuilder.org/stable/jll/): shipping a binary artifact through BinaryBuilder/Yggdrasil is a real distribution step and belongs in a follow-up rather than being faked here.
:::

## Array layout and 1-based indexing

The C core stores points row-major as `(num_points, dim)` and connectivity as `(num_cells, nodes_per_cell)`. Because Julia is column-major, the **same memory** is naturally

```julia
points(m)                      # (dim, num_points)
connectivity(m, 1)             # (nodes_per_cell, num_cells)
point_data(m, "displacement")  # (components, num_points)
```

so nothing is ever transposed — exactly the reasoning the Fortran module documents. Node indices are **1-based** here, and the ±1 shift happens inside the **copying** accessors, where a copy is made anyway:

| accessor | copies? | node indices |
|---|---|---|
| `points(m)` | yes | n/a — coordinates carry no indices |
| `points_ptr(m)` | **no** (zero-copy borrow) | n/a |
| `connectivity(m, i)` | yes | **1-based** |
| `connectivity_ptr(m, i)` | **no** | **0-based** — the ABI's own |

A borrow cannot be shifted without copying it, which is the whole point of a borrow — hence the two names. `point_data_ptr`, `cell_data_ptr` and `field_data_ptr` follow the same pattern.

Index maps and permutations (`refine`, `convert_cells`, `subdivide`, `agglomerate`, `decimate`, `partition`, `reorder`) are copies, so they are 1-based too, and the C API's `-1` "pruned / absent" sentinel becomes **`0`** — never a valid 1-based index. That is verbatim the Fortran rule; the bindings agree deliberately.

`partition_labels` is the exception: those are part **ids**, not indices, so they stay in `0:nparts-1`.

## The borrow window

`points_ptr` and friends return a `MeshBorrow`, a real zero-copy view of the buffer inside the C++ core (`unsafe_wrap`). The C API's rule is that such a pointer stays valid until the next **mutating** `mio_mesh_*` call on that mesh, or until it is freed; read-only accessors never invalidate it.

`MeshBorrow` enforces exactly that rather than trusting the caller:

```julia
b = points_ptr(m)
b[1, 2]                                # fine
add_point_data!(m, "T", temperatures)  # a MUTATING call ends the window
b[1, 2]                                # BorrowError -- not a stale read
```

It also holds a reference to the owning `Mesh`, so the mesh cannot be garbage-collected while a view of it is alive. `parent(b)` returns the raw `Array` for hot loops; that escape hatch is **unchecked** and valid only inside the window.

## Memory management

A `Mesh` releases its handle through a **finalizer** — the one real difference from the Fortran module, where meshes are freed explicitly with `call m%free()`. `close(m)` releases one deterministically and is idempotent.

`refine` takes an optional cell selection: at most one of `cells` (global block-major, **1-based** here), `region` (a cell region selects its cells, a point region every cell with any node in it; a side region is an error) and `where_array` + `where_op` + `where_value`, plus `closure` (`"redgreen"`, local, or `"propagate"`, which reaches the whole edge-connected component) and `record_levels`. With no selector every cell is refined. `record_hierarchy` attaches `refine:cell_id`/`refine:parent_id` — the persistent parent/child hierarchy a multigrid caller resolves across the sequence of meshes it keeps, riding the raw C-side numbering unshifted like `partition_labels`' part ids — and forces `refine:entity` (the multigrid prolongation stencil) to be attached even when the closure leaves no hanging node. See [refine](/refine#refinecell_id-and-refineparent_id).

Operations producing an opaque C result (`split`, `partition`, `reorder`, `refine`, `decimate`, `convert_cells`, `subdivide`, `agglomerate`) always **transfer ownership** of the mesh out of that result rather than handing back a borrow into it, so a piece stays valid after the result is gone:

```julia
for (key, piece) in mio.split(m; by="type")
    mio.write(piece, "part_$key.vtu")   # still valid; the result is long gone
end
```

## Names that would shadow `Base`

`read`, `write`, `convert`, `merge`, `split` and `diff` are **not exported**, because they collide with `Base`. Call them qualified:

```julia
import MeshioPlusPlus as mio
m = mio.read("bracket.msh")
```

Everything else — the accessors, the setters, `regions`, and the remaining operations — is exported. `close` and `isopen` are proper `Base` method extensions.

## Error handling

Every failure raises a `MeshioError` carrying the C API's own thread-local message; a status code never reaches the caller.

```julia
julia> mio.read("/nope.vtu")
ERROR: MeshioError (read error): meshio++: cannot open file '/nope.vtu'
```

Using a borrow outside its window raises `BorrowError` instead.

## Named regions

```julia
add_region!(m, "inlet", :point, [1, 3, 5])
add_region!(m, "solid", :cell, [1, 2]; dim=3, tag=17)
add_region!(m, "wall", :side, Int64[1 2; 0 2])   # (cell, facet) pairs

for r in regions(m)
    println(r.name, " ", r.kind, " ", size(r.entries, 2), " entries")
end
```

Entries are 1-based, with one exception documented in [`doc/regions.md`](/regions): for a `:side` region the second row is a **facet ordinal within the cell type**, not a mesh index, so it is passed through unshifted — the same rule Fortran's `partition_labels` follows.

## Selective reads and time steps

`ReadOptions` narrows what a read materializes, and since v8.6.0 also picks which time step of a multi-step file to decode:

```julia
using MeshioPlusPlus

m = MeshioPlusPlus.read("big.vtu"; options = ReadOptions(points_only = true))

# 0 (the default) is the first step; negative counts from the end.
last = MeshioPlusPlus.read("run.exo"; options = ReadOptions(time_step = -1))

meta = read_metadata("run.exo")
meta.time_values          # [0.0, 0.5, 1.0] -- how many steps `time_step` may name
```

Format inference is the shared registry's, the C API's: the extension picks the format, with seven exceptions: a `.mesh` file whose first line names an MFEM mesh is read as `mfem` rather than `medit` (v16.5.0), Z88's fixed file names (`z88i1.txt`, `z88o2.txt`, …) are `z88` rather than xyz (v16.7.0), a `.dat` file that opens as an MSC Marc input deck is `marc` rather than `tecplot` (v16.8.0), a file named `d3plot`, which has no extension, is LS-DYNA's `lsdyna_d3plot` (v16.9.0), a file named like `crashA001` (a stem, `A` and three or more digits, no extension) is an OpenRadioss animation file, `radioss_anim` (v16.11.0), a file named `binout` (or `binout0000`…) is LS-DYNA's `lsdyna_binout`, and one named like `crashT01` (a stem, `T` and two digits, no extension) an OpenRadioss time-history file, `radioss_th` (both v16.12.0). The two content checks are for reading: a write goes by the extension alone, so a `.dat` is written as Tecplot and a `.mesh` as Medit whatever the file being replaced holds (v16.17.0), and a Marc deck is written by naming `marc`. Formats added to the registry (Patran, Femap and MFEM in v16.5.0; libMesh, Z88, Abaqus `.fil` and Radioss in v16.7.0; Marc's `marc` and `marc_t19` and the `ansys_rst_cyclic` full rotor in v16.8.0; `lsdyna_d3plot` and `nastran_op2` in v16.9.0; Tecplot's binary `.plt` in v16.10.0; `radioss_anim`, the libMesh writer and its `.gz`/`.bz2` names in v16.11.0; `lsdyna_binout` and `radioss_th` in v16.12.0; `vtx`, DOLFINx's ADIOS2 `.bp` directories, and `szplt`, Tecplot SZL, in v16.13.0, when the library was built with ADIOS2 and TecIO; the `radioss` and `marc` writers in v16.17.0, the `radioss_th`, `nastran_op2`, `lsdyna_binout` and `abaqus_fil` writers after them) reach this binding with no code change.

Out of range is an error naming the available count, never a silent clamp. `meta.time_values` is empty for a format with no time concept, so `length(...)` is always safe. Honoured by `exodus`; see [Selective reads](/selective_read).

## Transient (time-series) XDMF writing

`XdmfSeries` is the write half of the above, and the one writer `write` cannot express: a series is a **stateful** multi-call object, so it gets its own handle rather than a keyword. The grid goes out once and each solve appends a cheap step. See [XDMF time series](/xdmf_time_series).

```julia
using MeshioPlusPlus

s = XdmfSeries("simulation.xdmf")          # "HDF" by default
write_points_cells!(s, mesh)               # the static grid, once
for k in 0:nsteps-1
    solve!(mesh)
    write_data!(s, k * dt, mesh)           # point_data/cell_data only
end
num_steps(s)                               # 10
finalize!(s)                               # close(s) would do this too
close(s)
```

The function form closes the series afterwards even if the body throws:

```julia
XdmfSeries("simulation.xdmf"; data_format = "XML") do s
    write_points_cells!(s, mesh)
    write_data!(s, 0.0, mesh)
end
```

`data_format` is `"HDF"` (the default; needs an HDF5-enabled library), `"XML"` (everything inline in the `.xdmf`) or `"Binary"`; `gzip_level` applies to `"HDF"` datasets only and is negative (no compression) by default. An unknown format, or `"HDF"` against a library built without HDF5, throws a `MeshioError` from the constructor.

Two things worth knowing before reading the result back:

- the `.xdmf` light data is **buffered until the series is finalized**, so the file is only readable after `finalize!(s)` (or `close(s)`, which finalizes first);
- `close` swallows a write failure during that implicit finalize — call `finalize!` explicitly to see one.

Like `Mesh`, the handle is released by a GC finalizer and `close` is the deterministic, idempotent form. The name is `finalize!` rather than `finalize` because `Base.finalize` runs an object's GC finalizer and means something quite different. Reading a finished series back is the ordinary `MeshioPlusPlus.read(path; options = ReadOptions(time_step = k))`.

## Sets ↔ integer data

`sets_to_data(mesh, :point; data_name=nothing, join_char="-", order=String[])` and `data_to_sets(mesh, :cell, "material")` operate on region-backed memberships and return a new mesh. Labels remain zero-based (uncovered rows are -1), even though region accessors are 1-based. See [sets/data conversions](data_manage.md#sets--integer-data).

## Documented gaps

These are gaps in the **C ABI**, shared with the [Fortran](/fortran) and [R](/r) bindings; the Julia package invents no workaround for any of them:

- **Four Python-only formats.** `pmsh`, `zarr`, `cae` and `usd` (v10.35.0, the physics-ML data path) are registered in the Python layer only, not in the shared C++ dispatch registry, so this surface cannot read or write them. They are export targets for a training pipeline rather than interchange formats; see [formats](/formats).

- the **`frozen` pin mask** of `smooth` and `decimate`;
- **per-cell-type counts** in `stats` — use `cell_block_types` with `cell_block_info`;
- ~~ragged block connectivity~~ — **closed in v9.15.0**: `polygon_block` / `polyhedron_block` read them as nested 1-based vectors and `add_polygon_block!` / `add_polyhedron_block!` build them. `connectivity` still throws, since a ragged block has no matrix. See [Polyhedra and ragged cells](/polyhedra);
- the combined **`data_manage`** — `data_drop` / `data_keep` / `data_rename` compose to the same effect;
- **Exodus provenance strings** (`qa_records` / `info_records`): they ride the `ExodusInfo` side channel, which like `MedInfo` does not cross the flat ABI, so `mesh.info` has no counterpart here. Geometry, data, regions and time steps are unaffected.

One further limitation is a *format* one rather than an ABI one, and easy to trip over: **gmsh does not currently round-trip named regions.** The writer emits the `$PhysicalNames` entry but does not attach the physical tag to an entity in `$Elements`, so a reader finds nothing to rebuild the group from. This is pre-existing meshio++ behaviour, reproducible from Python; `abaqus` round-trips regions correctly.

## Tests

```sh
MESHIOPLUSPLUS_LIB=/opt/meshioplusplus/lib/libmeshioplusplus.so \
  julia --project=bindings/julia/MeshioPlusPlus -e 'using Pkg; Pkg.test()'
```

The suite uses the same deliberately non-square fixture as [`tests/fortran/test_fortran_api.f90`](https://github.com/loumalouomega/meshioplusplus/blob/master/tests/fortran/test_fortran_api.f90) — 5 points × 3 dims, 2 tetrahedra × 4 nodes, 3-component vector data — so a transposed mapping or a missed shift cannot cancel out and pass anyway. It pins the column-major identity, the 1-based/0-based accessor pair, the borrow window, regions, and every operation.

## Software rendering (v16.34.0)

A deterministic software rasterizer and its terminal, HTML and PNG encodings, with no display or GPU: see [terminal rendering](/tui) for what is drawn and every option. Every option is a keyword argument of the same name as the C `mio_render_opts` field.

```julia
frame = render(mesh; width=800, height=600, supersample=2,
               color_by="temperature", colorbar=true, clip=(2, 98), isolines=8)
frame.rgba        # Array{UInt8,3} (4, width, height); frame.rgba[:, x, y] is a pixel, row 1 on top
frame.cell_ids    # (width, height): the input cell at each pixel, -1 for none
frame.range       # (vmin, vmax), or nothing when no field is mapped
frame.notes       # the colour range, ticks and keys

print(render_text(mesh; cols=80, rows=24, encoding="braille", color_by="temperature"))
png = render_png(mesh; width=800, height=600)          # Vector{UInt8}
write_snapshot("part.png", mesh)                       # .png .txt .ansi .html .cast
```

**Cut-aways (v16.36.0).** `cutaway=[px, py, pz, nx, ny, nz]` (a point and the normal of the side kept), or a collection of up to two such planes, with `cutaway_tint=rgba(...)`, clips the geometry away and tints the back faces now in view: `render(mesh; view="-x", cutaway=[0.5, 0, 0, 1, 0, 0])`. A wrong length is an `ErrorException`; a zero normal a `MeshioError`.

**Streamlines (v16.37.0).** `streamlines="velocity"` follows a vector point array over the mesh's cells, with `stream_seeds` (default 40), `stream_length` (model diagonals each way, default 0.5) and `stream_color=rgba(...)`: `render(mesh; streamlines="velocity", stream_seeds=60)`. An unknown array or a seed count outside 1 to 10000 is a `MeshioError`.

`render` copies the frame out of the library, so a `Frame` owns its arrays and nothing needs closing. Colours are `0xRRGGBBAA` integers or `(r, g, b[, a])` tuples (`rgba(r, g, b, a)` builds the integer); `shading`, `edges`, `scale` and `diagnostic` take their names (`"smooth"`, `"feature"`, `"log"`, `"orientation"`, ...). `render_text(...; format="html")` returns a self-contained page, which a notebook shows with `display("text/html", ...)`. An unknown option or enum name is an `ErrorException` listing the choices; a library failure (an unknown array, a bad range) is a `MeshioError`. `detect_color_depth(; no_color, color_term, term)` returns `"truecolor"`, `"256"`, `"16"` or `"mono"`.

## v16.25.0 additions

- `agglomerate(mesh; target_group_size=8, merge_coplanar_faces=false, coplanar_angle=1.0, min_sphericity=0.0) -> (; mesh, cell_map, num_faces_merged, num_rejected)` — the coplanar-face merge and sphericity gate of [agglomerate](/agglomerate).
- `resample(seq, out_path, times; method=:linear, clamp=false, blend_points=false, format="")` — resample a `Sequence` onto new times (`:linear`, `:nearest` or `:previous`), at most two source meshes live; and `blend_steps(a, b, w; blend_points=false) -> Mesh`. See [Resampling onto new times](/sequences#resampling-onto-new-times). Both are exported.

## v15.4.0 additions

- `compute_normals(mesh; point_normals=true, cell_normals=false, weight=:angle, split_angle=nothing, record_parent_ids=false, region="") -> (; mesh, quality, num_isolated, num_undefined, num_degenerate, num_split_points, num_added_points)` — point and cell normals of a surface, written as `normals` (`(n, 3)` point data; `(cells, 3)` cell data with `cell_normals`). `split_angle=nothing` gives one smooth normal per point; a number of degrees in `[0, 180]` duplicates points at creases so every point carries exactly one normal, appending the copies after the original points while cells keep their numbering. Never reorients: check `quality.inconsistent_pairs`. See [normals](/normals).
- `check_quality(mesh; require="", max_inverted=0, max_degenerate=0) -> (; passed, num_checks, num_failed, num_cells, num_inverted, num_degenerate, summary)` — the [quality gate](/mesh_quality#quality-gate).
- `feature_edges(mesh; feature_angle=30.0, feature=true, boundary=true, non_manifold=true, inconsistent=true, region="") -> (; mesh, num_feature, num_boundary, num_non_manifold, num_inconsistent)` — the sharp, open, non-manifold and inconsistently wound edges of a surface (or of a volume mesh's skin) as a `line` mesh with `feature:kind`/`feature:angle`. See [feature edges](/feature_edges).
- `tessellate(mesh; levels=2, curved=true, fields=true, record_stencil=false) -> (; mesh, num_source_points, num_source_cells, num_curved_source_cells, num_pass_through_source_cells, num_points, num_cells, num_boundary_facets, num_interior_facets, num_facets_with_bad_count)` — isoparametric tessellation of `triangle6`/`quad8`/`quad9`/`tetra10`/`hexahedron27` cells into triangles and tetrahedra, with `tessellate:source_point` point data (0-based source point, or -1), `tessellate:source_cell`/`tessellate:sub_index` cell data and, with `record_stencil`, the `(27, num_points)` `tessellate:stencil`/`tessellate:weights`; see [tessellation](/tessellation).
- `segment(mesh; feature_angle=30.0, normal_tolerance=0.0, curvature_tolerance=0.0, min_patch_cells=1, region="", prefix="patch", method="angle", tolerance=1e-3) -> (; mesh, num_patches, num_merged, num_feature_edges, num_boundary_edges, num_non_manifold_edges)` — smooth patches bounded by the feature edges, as `segment:patch` cell data, `segment:adjacency` field data and one cell region per patch, `method="primitive"` ending a patch where its plane, cylinder, cone, sphere or torus does; see [segment](/segment).
- `fit_primitives(mesh; patch_data="segment:patch", tolerance=1e-3, relative=true, min_inlier_fraction=0.9, types=("plane", "sphere", "cylinder", "cone", "torus"), max_iterations=50) -> (; mesh, num_patches, num_plane, num_sphere, num_cylinder, num_cone, num_torus, num_unfitted)` — a primitive fitted to each patch, with `fit:*` cell data and `fit:patch_*` field data; see [fit primitives](/fit_primitives).
- `refine_patches(mesh; tolerance=1e-3, relative=true, normal_weight=1.0, epsilon=1e-3, feature_angle=30.0, max_iterations=10, min_patch_cells=3, ransac_seed=0x6d657368696f2b2b, types=nothing, allow_split=true, allow_merge=true, prefix="patch") -> (; mesh, num_patches, num_iterations, num_splits, num_merges, num_absorbed, num_primitive, num_freeform, num_facets, tolerance, explained_area_fraction, initial_explained_area_fraction, hausdorff_max, hausdorff_rms, initial_energy, energy, stop)` (v17.11.0) — refines the `segment:patch` labels by variational shape approximation, moving boundaries, splitting a freeform patch where seeded RANSAC finds a primitive inside it and merging neighbours one primitive explains, and returns the mesh with the patches, regions and `fit:*` data rewritten plus `refine:*` field data. `types=nothing` allows every primitive, `ransac_seed` is a `UInt64` (exact over its whole range, unlike R's doubles) and `stop` is `:unchanged`, `:small_gain`, `:iteration_cap` or `:energy_rose`; a mesh without `segment:patch` or an option out of range throws a `MeshioError`. See [refine patches](/refine_patches).
- `iga_refine(info; subdivisions=1, elevate=0, sample_subdivisions=0, trim_samples=16, edges=true, trim=true) -> (Mesh, FormatInfo)`, `iga_sample(info; sample_subdivisions=0, …) -> (Mesh, FormatInfo)` and `iga_evaluate(info, patch, params; derivatives=0, trim_samples=16) -> (; points, fields, derivatives, normals, inside)` (`patch` 1-based, `params` a `(d, N)` matrix; absent parts are `nothing`) act on the NURBS model of a `read_with_info` handle (MFEM NURBS, CAD JSON, STEP, IGES); the returned `FormatInfo` writes back through `write_with_info`. See [isogeometric layer](/iga).
- `mesh_cad_surface(info; target_size, max_chordal_deviation=0.0, min_points_per_edge=0, tolerance=0.0, record_face_ids=true, edges=true, match_unreferenced_trims=true, strict=false, max_points=50_000_000) -> (; mesh, num_points, num_triangles, …, faces, warnings)` triangulates the faces of a `cad_json`, `step` or `iges` handle into one conforming surface mesh; `faces` holds each face's 1-based index, `status` and `reason`. See [mesh generation from CAD](/mesh_cad).
- `fit_cad(mesh; step_path="", tolerance=0.0, exact_edges=true, types=("plane", "sphere", "cylinder", "cone", "torus"), mesh_deflection=0.0, max_demotion_rounds=8, freeform=true, freeform_degree=3, freeform_max_control=64, freeform_split=true, fillets=true, threads=true, merge_faces=true, backend="auto", refine_segmentation=false, max_refine_iterations=10, refine_epsilon=1e-3, normal_weight=1.0, allow_merge=true, allow_split=true, ransac_seed=0x6d657368696f2b2b, refine_tolerance=1e-3, relative_tolerance=true, freeform_general=false, freeform_tolerance=0.0, freeform_smoothing=1e-12, max_freeform_poles=40, max_freeform_iterations=30, max_freeform_splits=2, sew_tolerance=0.0, extend_surfaces=false, snap_junctions=false, face_deviation=false, freeform_mode="legacy", freeform_parameterization="auto", extraordinary_vertex="g1x", local_refine="knots", freeform_g1_tolerance_deg=1.0, freeform_max_control_points=400, layout_target_patch_size=0.1, layout_max_valence_deviation=2, resegment_rounds=2, layout_snap_to_features=true, validate_by_surface_sampling=true) -> (; mesh, backend, valid, closed, tolerance, volume, skin_volume, volume_ratio, max_deviation, num_faces, num_plane, num_sphere, num_cylinder, num_cone, num_torus, num_bspline, num_faceted_faces, num_patches, num_demoted, num_exact_edges, num_approx_edges, num_solids, num_split_patches, num_fillet, num_helical, num_variable_fillet, num_blend, num_rounded_thread, num_vertex_blend, num_chamfer, num_conic_blend, num_tapered_thread, num_variable_pitch_thread, refined, refinement_kept, refine_iterations, num_refine_splits, num_refine_merges, initial_explained_area_fraction, explained_area_fraction, num_freeform_general, num_freeform_refinements, avg_freeform_poles, num_extended_edges, num_snapped_junctions, hausdorff_max, hausdorff_rms, brep_checked, freeform_mode, facet_area_fraction, num_freeform_adaptive, num_network_patches, num_network_edges, num_extraordinary_vertices, num_resegmented, num_escalations, max_g1_angle_deg, max_freeform_distance, rms_freeform_distance, num_face_records)` — the fitted patches as a CAD solid with analytic faces (and, for the patches no primitive explains, rolling-ball fillet and blend faces with `fillets`, helical thread faces with `threads` and, with `freeform`, B-spline faces, a patch folding over every flattening split into pieces with `freeform_split`), written as STEP to `step_path` when given; `backend` `"occt"` needs a library built with OpenCASCADE, which `has_occt()` reports, and throws a `MeshioError` naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON` otherwise, `"native"` writes the fitted primitives as analytic STEP faces with no CAD kernel, and `"auto"` picks OpenCASCADE when present; `backend` in the result names the one that ran. Since v17.11.0, `refine_segmentation=true` first refines the patches as `refine_patches` does (through `mio_fit_cad2`), with `refine_tolerance` (relative while `relative_tolerance`), `normal_weight`, `refine_epsilon`, `max_refine_iterations`, `allow_split`, `allow_merge` and `ransac_seed`; `refined` says whether it ran, `refinement_kept` whether the solid was built from the refined patches (the OpenCASCADE backend also builds it from the patches as given and keeps the better), and the refinement's counts and explained area shares are zero when it did not. `freeform_general=true` tries a patch the plane and cylinder freeform fit leaves unexplained as a general B-spline face (`freeform_tolerance`, `freeform_smoothing`, `max_freeform_poles`, `max_freeform_iterations`, `max_freeform_splits`; `num_freeform_general`, `num_freeform_refinements` and `avg_freeform_poles` report it), and, OpenCASCADE only, `sew_tolerance` (0: the tolerance) sets the sewing tolerance, `extend_surfaces` extends freeform surfaces past their patches (`num_extended_edges`), `snap_junctions` places the corners where analytic faces meet by least squares (`num_snapped_junctions`) and `face_deviation` measures each face against the skin (`hausdorff_max` and `hausdorff_rms`, `NaN` when not measured, and the `fit_cad:face_hausdorff` / `fit_cad:face_rms` field data); `closed` says whether every shell is closed. Since v17.12.0 the [freeform modes](/fit_cad#freeform-modes) go through `mio_fit_cad3`: `freeform_mode` (`"legacy"`, `"adaptive"`, `"network"` or `"auto"`), `freeform_parameterization` (`"auto"`, `"projection"`, `"lscm"`, `"harmonic"`), `extraordinary_vertex` (`"g0"`, `"g1"`, `"g1x"`, `"g1xx"`) and `local_refine` (`"knots"`, `"thb"`) as strings, an unknown one an `ArgumentError`, with the layout and budget keywords above; the result always carries `brep_checked`, `freeform_mode`, `facet_area_fraction` (the share of the area left on facet faces), the network and escalation counts, `max_g1_angle_deg`, the freeform distances and `num_face_records`, and outside `"legacy"` the mesh carries the face records as `fit_cad:face_table` field data (one row per patch tried as a surface). See [fit CAD](/fit_cad).
- `hausdorff_distance(a, b; face_samples=0, region_a="", region_b="") -> (; distance, a_to_b, b_to_a, mean_a_to_b, rms_a_to_b, mean_b_to_a, rms_b_to_a, num_samples_a, num_samples_b, worst_point_a, worst_point_b)`. See [Hausdorff distance](/hausdorff).
- `edit_regions(mesh, op, inputs; output="", kind=nothing, dim=nothing, tag=nothing, keep_inputs=true) -> Mesh` (`op` one of `:union`, `:intersection`, `:difference`, `:rename`, `:retag`, `:delete`) and `remove_region!(mesh, index)` (1-based). See [editing regions](/regions#editing-regions).
- `region_adjacency(mesh; regions=nothing) -> Mesh` returns conforming shared facets between Cell regions, falling back to multiple cell blocks when fewer than two Cell regions exist — see [region adjacency](/region_adjacency).
- `find_interface(mesh_a, region_a, region_b; mesh_b=nothing, mode=:conforming, master=:a, gap_tolerance=0, angle_tolerance=30, overlap_tolerance=0) -> (; mesh, report)` returns the master facets and source Side entries; `contact_pairs(slave, slave_points, master_cells; master_mesh=nothing, tolerance=0, require_complete=false)` returns copied point/facet projections; `split_interface(mesh, side; add_cohesive=false)` duplicates the incident-cell fans. Point and global cell ids are 1-based; local facet/subfacet ordinals remain 0-based. See [interfaces and contact](/region_adjacency).
- `match_periodic_nodes(mesh, slave, master; translate=nothing, matrix=nothing, atol=1e-8, require_complete=true) -> (; slave, master, unmatched, num_fixed, max_residual)` — 1-based node ids; `matrix` is 16 numbers, row-major. See [periodic node pairs](/periodic).

## v10.9.0 additions

- `hessian(mesh, array; method=:green_gauss, location=:cell, output="", overwrite=false) -> (; mesh, num_skipped, num_fallback)` — the Hessian (second derivative) of a **scalar point-data** field, [`gradient`](@ref)'s companion one order further. A composition of TWO `gradient` calls, not a new numerical kernel: the field is differentiated once (point location), then that `(n, 3)` gradient is differentiated again with the default `:gradient` operator, producing `(n, 9)` — the flattened row-major 3x3 Hessian, `H[i,j]` at index `i*3+j`. `method` is forwarded to BOTH internal passes. See [`doc/hessian.md`](hessian.md).

A field that is at most LINEAR has an exactly zero Hessian everywhere — the one mesh-shape-independent guarantee. For a genuinely quadratic field the composition is exact on a structured/symmetric mesh away from its own boundary and a good, standard, but genuinely approximate curvature estimate on an irregular mesh. Input must have exactly one component — a vector field's Hessian is a separate quantity per component.

A curvature-driven refinement indicator needs no new function: `norm(...)` in [`data_calc`](@ref) on the 9-component output is exactly its Frobenius norm, ready for [`refine`](@ref)'s `where` selector.

`hessian` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v10.8.0 additions

- `data_integrate(mesh, names=String[])` — cell-measure-weighted total/mean of one or more `cell_data` arrays, [`gradient`](@ref)'s integration counterpart (`gradient` differentiates a field, this integrates one), returning a `Vector{NamedTuple}` with `name`, `num_components`, `domain` (`num_cells`, `num_skipped`, `components` — a vector of `(total, mean, domain_measure, num_nan)` NamedTuples) and `regions` (the same shape, one per named Cell region present). See [`doc/field_integration.md`](field_integration.md).

Every sum is weighted by `|measure(cell)|`; a cell whose measure is not computable, or a component whose value is non-finite, is excluded from that component's numerator **and** denominator, never given a fallback weight of 1. Regions are not a partition: a cell in two regions contributes fully to both. A `point_data`-only name throws a `MeshioError` naming `data_point_to_cell` as the fix.

`data_integrate` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v10.5.0 additions

- `undo_green(coarse, fine)` — green-element undo: restores `fine`'s transitional (closure-only) cells back to their original parent, read verbatim from `coarse` — a **lookup, not a reconstruction**, since [`refine`](@ref) never renumbers or prunes points, so a green parent's exact connectivity and cell_data are already sitting, byte-for-byte, in `coarse` at the row `fine`'s `refine:parent_id` names. Returns `(; mesh, num_groups_undone, num_cells_removed)`. See [`doc/undo_green.md`](undo_green.md).

A **two-mesh** operation, like [`interpolate`](@ref): `coarse` is the mesh a prior `refine(coarse, ...; record_hierarchy=true, record_levels=true)` call was run on, `fine` is that call's output — both flags are required, `record_hierarchy` alone does not imply `record_levels`. The six reserved `refine:*` arrays are always dropped from the output; only a single-pass (`levels=1`) hierarchy is supported, a deeper multi-level hierarchy being refused by name. Unlike `subdivide`/`agglomerate`, this operation has **no winding repair or discrete sign branch anywhere in it** — it is pure array bookkeeping, which on the Python side means it has a full numpy reference implementation rather than being C++-core-only (this binding always calls the installed C library either way).

`undo_green` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

- `conservative_interpolate(source, target; arrays=String[], default_value=0.0, on_conflict="error")` — mass-preserving field transfer, [`interpolate`](@ref)'s sibling: conserves `sum(target value * target measure) == sum(source value * source measure)` over the region the two meshes share, a property `interpolate`'s pointwise sampling does not have. Both meshes are simplexified first, accepting ragged/polyhedron blocks for free. Unlike `interpolate`, an empty `arrays` means every source point_data **and** cell_data array. Output arrays are always Float64. Like `undo_green`, a **two-mesh** operation, and like `subdivide`/`agglomerate` this is C++-core only — the 3D clip kernel's discrete branches could disagree with a second implementation near a degenerate overlap, so there is no pure-Julia fallback either; the installed C library is always called. See [`doc/conservative_interpolate.md`](conservative_interpolate.md).

## v10.4.0 additions

- `agglomerate(mesh; target_group_size=8)` — polyhedral coarsening, the many-to-one counterpart to [`subdivide`](@ref): greedy seed-and-grow over the mesh's shared-face dual, absorbing face-adjacent neighbours by accumulated shared-face area until each group reaches `target_group_size` members, then emitting one polyhedron per group whose faces are exactly its external boundary — conserving volume exactly, since internal faces are simply dropped rather than re-triangulated. Returns `(; mesh, cell_map)`. See [`doc/agglomerate.md`](agglomerate.md).

Unlike every other opaque-result operation here, `cell_map` is a **single flat**, not block-indexed, 1-based array: an agglomerated cell's output index is a function of which group it joined, not which input block it came from, so `_result_map` (the same single-array reader `split`'s node map uses) reads it rather than `_result_cell_maps`. Like `subdivide`, there is no `point_map` — points are never pruned or renumbered, so `clean(mesh; remove_orphans=true)` is the documented follow-up for a minimal point set. A non-manifold input (a face shared by three or more cells) raises a `MeshioError` naming the face rather than guessing.

`agglomerate` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v10.3.0 additions

- `subdivide(mesh; record_parent_ids=false)` — polyhedral refinement: one polyhedral child per face of every eligible 3D cell, connected to a new interior point, returning `(; mesh, cell_maps)`. Needs no per-type template table — tabulated types (reduced to corners for a quadratic variant) and existing polyhedron blocks are handled uniformly through the same `detail::cell_rings`/`orient_rings` machinery `gradient` uses — and is automatically conforming, unlike [`refine`](@ref). See [`doc/subdivide.md`](subdivide.md).

Unlike [`convert_cells`](@ref), there is no `point_map`: `subdivide` never prunes or renumbers an original point. `cell_maps[b]` is 1-based input cell → the index of its **first** child (one per face) in the corresponding output block, the same `FirstChild` shape `convert_cells` already uses for its own one-to-many splits.

`subdivide` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v10.2.0 additions

- `estimate_error(mesh, array; method=:zz, marking=:none, marking_value=0.0, output="", marked="", overwrite=false)` — the Zienkiewicz-Zhu recovery-based error indicator of a **point-data** field, plus optional marking, returning `(; mesh, global_error, num_skipped, num_marked)`. A composition of `gradient` with the point↔cell averaging round trip, not a new kernel. See [`doc/error.md`](error.md).

`error:zz` is always attached; `error:marked` too when `marking` is not `:none`, so `refine`'s own `where` selector needs no change to consume it — `refine(mesh, where="error:marked > 0.5")`. Cells that cannot be evaluated read `NaN` in `error:zz` and `0` (never `NaN`) in `error:marked`, counted in `num_skipped` and excluded from `global_error`/`num_marked`.

`estimate_error` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v9.11.0 additions

- The [settings pipeline](pipeline.md): `run_pipeline_file(settings_path)` and `run_pipeline_json(json_text)` run a whole `settings.json` (read → operation chain → write; PascalCase ops/keys) through the C++ engine, and `pipeline_has_json()` reports whether the loaded library carries the JSON parser — a build without it raises a `MeshioError` naming `-DMESHIOPLUSPLUS_WITH_JSON=ON` rather than missing a symbol. The flat ABI carries JSON text only; the structured run report is a recorded follow-up. All three are exported (none shadows `Base`).

## v9.10.0 additions

- `gradient(mesh, array; operator=:gradient, method=:green_gauss, location=:cell, output="", component=-1, overwrite=false)` — the gradient, divergence or curl of a **point-data** field, returning `(; mesh, num_skipped, num_fallback)`. See [`doc/gradient.md`](gradient.md).

Two spellings to note. `method` is given as an underscored symbol (`:green_gauss`, `:least_squares`) and translated to the C API's hyphenated name by the binding, because a Julia symbol cannot carry a hyphen. And `component` is negative for **every** component here — deliberately the opposite of `isosurface`, where negative means the row magnitude.

`gradient` shadows nothing in `Base`, so unlike `read`/`write`/`split` it is exported.

## v9.1.0 additions

- `ReadOptions(; lenient=true)` — see [`doc/selective_read.md`](selective_read.md).
- `ReadOptions(; piece=k)` — keep one piece of a partitioned file (VTKHDF partitions or composite blocks): `0` is the first, negative counts from the end, `nothing` (the default) merges every piece into one mesh with one cell region per piece; out of range fails the read naming the piece count. See [`doc/selective_read.md`](selective_read.md#picking-a-piece).
- `ReadOptions(; drop_ghosts=true)` — remove the ghost cells (halo) of a partitioned `.pvtu`, `.pvtp` or `.pvd`: every cell with a `vtkGhostType` bit set, and the points only they used. `false` (the default) keeps them; every other reader ignores it. See [`doc/selective_read.md`](selective_read.md#dropping-ghost-cells).
- XDMF series: `flush!(s)`, `finalized(s)`, `XdmfSeries(path; mode=:append, auto_flush=false)`, and `write_data!(s, t, Dict("u" => values))` for writing a step from raw arrays with no `Mesh` in between. An `n x k` matrix is transposed on the way out, since Julia is column-major and the C ABI expects `k` components per entity row-major.

`flush!` is named with a bang for the same reason `finalize!` is: `Base.flush` means "flush this IO stream". 

## Format side channels

Gmsh `read_with_info` returns a `FormatInfo`; `gmsh_info(info)` copies `bounding_entities` (signed tags per 4.1 cell block) and `periodic` links `(entity, affine, node_pairs)`. `entity` is `(dimension, slave entity tag, master entity tag)` in raw file ids, `affine` has 0 or 16 coefficients, and `node_pairs` is `(2, N)` with **1-based point rows**, retaining ordering and duplicates. `write_with_info(mesh, info, path; format="gmsh")` or `"gmsh22"` restores the records. `close(info)` releases the handle, not the copies. Operations do not remap the channel; rebuild it after point/topology edits.

MDPA's `mdpa_info(info).submodelparts` entries contain `name`, `data`, `tables`, `geometry_ids` and `constraint_ids`. The last two lists retain original file ids, ordering and duplicates without Julia's point-index shift; mesh operations do not remap them. `write_with_info` restores nested membership-only parts, with constraints themselves kept opaque in `raw_blocks`.

MFEM NURBS meshes (v16.40.0): `read_with_info(path; format="mfem")` returns the patches an `MFEM NURBS mesh` was sampled from, and `nurbs_info(info)` copies them into a `NamedTuple`: `patches` (`degrees`, `knots`, `control_points` as `(space_dimension, points)`, `weights`, `fields` as `(components, points)` matrices), `field_names`, `topology` (`nothing` when the source states none) and the cell map `cell_patch`, `cell_span` (`(3, cells)`) and `boundary_cell_patch`. Vertices, patches and knot spans are 1-based; attributes and knot-vector ids stay as the file spelled them. `write_with_info(mesh, info, path; format="mfem")` writes the patches back ([MFEM](formats/mfem.md#the-patches-themselves-the-nurbs-side-channel)).

Kratos CAD JSON (and, with OpenCASCADE, STEP and IGES): `read_with_info(path; format="cad_json")` returns the B-rep a mesh was sampled from, and `brep_info(info)` copies it into a `NamedTuple`: `patches` (`degrees`, the full `knots`, `control_points` as `(3, points)` with direction 1 fastest, `weights`), `faces` (`id`, `name`, `patch`, `swapped_normal`, `trimmed`, `rational`, `reduced_knots`, `control_point_ids`, `loops` of `(type, trims, extra)`, `embedded_edges`, `surface_extra`, `extra`; a trim is `(trim_index, direction, curve, extra)`), `edges` (`id`, `name`, `curve` or `nothing`, `topology`, `extra`), `vertices` (`id`, `name`, `topology`, `extra`), `breps` and the file's `extra`. A curve is `(degree, knots, reduced_knots, rational, control_points, weights, ids, active_range, extra)`. Indices (a face's `patch`, a brep's lists) are 1-based; ids (-1 where the file names the entity instead) and trim indices stay as the file spelled them. `write_with_info` writes the B-rep back as CAD JSON; the `MIO_BREP_*` section constants are exported for the raw accessors.

`read_with_info(path; format="", options=nothing)` returns `(mesh, info)`, `info` a `FormatInfo` for a format that has a side channel (MDPA: tables, geometries, `Mesh` blocks, sub-model-part data, text `ModelPartData`, raw blocks such as `Constraints`) and `nothing` otherwise; `write_with_info(mesh, info, path)` puts it back and `close(info)` releases it. `mdpa_info(info)` copies the MDPA side channel into a `NamedTuple` (`tables`, `geometries` with 1-based `(nodes, n)` connectivity, `mesh_blocks`, `submodelparts`, `raw_blocks`, …); file ids stay as the file spelled them. See [MDPA](formats/mdpa.md#the-blocks-the-mesh-cannot-hold-v16-26-0).

## Sequences (transient / multi-file datasets)

`Sequence` wraps the C API's ordered plan over a set of files (or the steps inside one multi-step file). Ordering is natural-numeric, so `out_9.vtu` precedes `out_10.vtu`; nothing is read until `read_step`, and the sequence caches nothing, so a 500-step dataset is traversable without materialising it.

```julia
seq = Sequence("out_*.vtu")               # or Sequence(["a.vtu", "b.vtu"])
for i in 1:length(seq)
    mesh = read_step(seq, i)              # OWNED; one mesh alive at a time
    @show MeshioPlusPlus.time(seq, i), MeshioPlusPlus.time_source(seq, i)
end
to_timeseries(seq, "series.xdmf")         # fan-in
to_timeseries(seq, "series2.xdmf"; ascii=true)  # "XML" data format, no HDF5 needed
resample(seq, "u_{index}.vtu", 0:0.1:2)   # onto new times, two meshes live
close(seq)

timeseries_to_sequence("series.xdmf", "step_{step}.vtu")   # fan-out
run_sequence_file("transient.json")                        # per-step chain
```

`Sequence` is also iterable (`for mesh in seq`), and has a `do`-block form that closes the handle even if the body throws.

**Not exported**, because they shadow `Base`: `path`, `step`, `time`, `time_source`. Reach them as `MeshioPlusPlus.time(seq, i)` — the same rule the binding already applies to `read`/`write`/`split`.

See [sequences](sequences.md) for the ordering rule, the time-value precedence and the streaming guarantee.

- `grid(dims; origin, spacing)`, `voxelize(m; resolution, fill, ...)`, `sample_distance(surface, points)`, `distance_to_surface(query, surface)` and `surface_watertight_check(m)` (v9.24.0) — regular grids and signed distance. `grid` shadows nothing in `Base`, so it is exported plainly. See [`doc/voxelize.md`](voxelize.md) and [`doc/sdf.md`](sdf.md).
- `rectangle(; dims, size, origin, cell_type, split, order)`, `box`, `disk`, `uv_sphere`, `icosa_sphere`, `cylinder`, the planar domains `reference_triangle`, `ngon`, `annulus`, `lshape`, `plate_with_hole` and `honeycomb`, and `reference_cells(types; pairs)` — parametric generators that take no input mesh, like `grid`; keyword arguments mirror the Python functions (`size` and `center` instead of numpy arrays), `cell_type` and `split` are strings, and an unknown one throws `ArgumentError`. None of the names shadows `Base`. See [parametric primitives](primitives.md).
- `octree(m; root_cells, min_level, max_level, balance, fill, transition, snap, surface_level, curvature_threshold, feature_angle, size_array, plane, ...)` — an octree (quadtree for a planar `line` skin) around `m`, fine only where the skin is; returns `(; mesh, dimension, num_leaves, num_hanging, finest_level, origin, root_size)`. See [`doc/octree.md`](octree.md).
- `octree_mesh(m; output, min_size, max_size, max_level, bounds, curvature_angle, feature_angle, gap_cells, size_array, surface_region, skin_grading, conformity, fill, keep_outside, attach_distance, attach_quality, max_cells, max_tets, ...)` — mesh the inside of a skin on a balanced octree: graded tetrahedra by default, `output=:castellated` for the hexahedral leaves. Every keyword left at `nothing` keeps the `mio_octree_mesh_opts_init` default (the Python defaults), and an unknown `output`, `conformity`, `skin_grading` or `fill` throws `ArgumentError`. Returns `(; mesh, output, guarantee_applies, ..., num_inverted, num_warnings, warnings)` with every `mio_octree_mesh_report` field and `warnings::Vector{String}`; the struct mirrors are size-checked by `_check_abi_layout`. See [`doc/octree_mesh.md`](octree_mesh.md).
- `compute_sdf(surface; structure, resolution, root_resolution, max_depth, ...)` (v9.25.0) — the grid and the field in one call, returning `(; mesh, dims, origin, spacing, max_depth, num_banded, quality)`. `structure=:octree` refines only near the surface and sizes itself from `root_resolution`/`max_depth`, so passing `resolution` or `cell_size` with it is an error; its output is 1-irregular.
- `crop_predicate(m, array; compare, value, record_ids)` (v9.25.0) — keep the cells whose scalar `cell_data` value satisfies a comparison, the same vocabulary `refine`'s `where_op` uses. There is deliberately no `mode`: a per-cell value has nothing for an all/any rule to reduce. See [`doc/crop.md`](crop.md).
