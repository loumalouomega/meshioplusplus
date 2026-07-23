<!--pytest-codeblocks:skipfile-->
# Changelog

This document records every released version of meshio++ — new formats, new operations,
notable enhancements, and breaking changes. Breaking changes are called out explicitly as
**Breaking:**; everything else is additive unless stated otherwise.

**Keep this file current: add an entry in the same change as every version bump.** See the
"Version bumps" section of `CLAUDE.md`.

## v8.5.0 (2026-07-23)

**Parallelization pass over the newer operations** — an audit of every
operation's serial loops parallelized the remaining safe ones (independent
iterations writing disjoint slots), with output staying **byte-identical**
across backends and thread counts: `detail/subset.cpp`'s connectivity
remaps/index maps (benefits crop, split and partition at once), quality's
per-cell array assembly and histogram (fixed-chunk partials merged serially in
chunk order; the min/max/sum summary stays serial so reported values are
unchanged), reorder's per-cell sort keys and connectivity rebuilds, merge's
ragged builds and dedup-filtered gathers, crop's kept-cell test (phase-split:
parallel flag pass + serial compaction) and clean's rep→final remaps. The
deliberately-serial determinism passes (first-seen dedup, FP scatter
accumulation, greedy/BFS loops, stable argsorts) are untouched. New
`Quality.DeterministicAcrossRuns` gtest.

**KOKKOS parallel backend** (`-DMESHIOPLUSPLUS_PARALLEL_BACKEND=KOKKOS`) — a
fifth backend for `meshioplusplus::parallel_for`, running on
`Kokkos::DefaultHostExecutionSpace` (host deliberately: loop bodies capture
host pointers, so device offload is served by the DLPack/CuPy handoff, not by
this backend). Bring-your-own like KaHIP (never picked by `AUTO`; point
`Kokkos_DIR` at an installed Kokkos ≥ 3.4), lazily initialized only when the
embedding application hasn't initialized Kokkos itself, `parallel_for_bw`'s
4-thread bandwidth cap preserved by range partitioning. Building Kokkos
alongside `MESHIOPLUSPLUS_BUILD_C_API=ON` needs Kokkos built with
`-DCMAKE_POSITION_INDEPENDENT_CODE=ON` — its default static archives aren't
PIC and fail to link into the shared `libmeshioplusplus.so`. CI: the new
`kokkos` job (cached PIC source build of Kokkos 4.5.01, full gtest suite, plus
a 2-thread pool re-run).

**`NDArray` buffer-allocator hook** (`meshioplusplus::set_buffer_allocator`) —
the `doc/gpu.md` Phase-2 enabler: every owning `NDArray` buffer is now
allocated through an optional process-global `BufferAllocator` (plain C
callbacks), so readers can fill e.g. CUDA pinned memory directly, removing the
staging copy of a later host→device transfer. Each buffer keeps a
`shared_ptr` reference to the allocator it was born with, so uninstalling the
hook never orphans live buffers; content, zero-copy and determinism contracts
are unchanged (views are unaffected). The Python/CuPy wiring
(`pinned_reads()`) is the recorded follow-up in `doc/gpu.md`.

## v8.4.0 (2026-07-23)

**GPU handoff at the I/O boundary** (`to_dlpack` / `to_cupy` / `from_cupy`) —
move a mesh's arrays to and from device memory through the standard exchange
protocols, with no file round-trip. Python-only
(`src/python/meshioplusplus/_gpu.py`); the C++/WASM/C/Fortran core is
untouched and stays dependency-free. Stated plainly: **a host→device move is
always a bus transfer** — what this removes is the file round-trip and every
*extra* copy on either side of that one transfer; "zero-copy" applies only to
host buffer sharing and on-device adoption, never to the transfer itself.

- `to_dlpack(mesh)` returns a payload of host numpy arrays, each natively
  speaking `__dlpack__`/`__dlpack_device__` (`kDLCPU`) — consumable by
  `np.from_dlpack`, `torch.from_dlpack`, JAX, Numba, CuPy, … It honors the
  interop layer's `zero_copy_only` **host** buffer-sharing contract.
- `to_cupy(mesh, float32=…, int32=…, pinned=…, stream=…)` transfers the
  payload to the CUDA device: points, per-block connectivity, point/cell data,
  and named regions as device **index arrays** (`Side` `(cell, facet)` pairs
  included — the one interop target that keeps them). Optional pinned-memory
  staging with async DMA on a caller-supplied stream. Deliberately no
  `zero_copy_only` parameter here — it would be a lie on a 100%-copy path.
- `from_cupy(payload)` rebuilds an ordinary `Mesh` with one deliberate
  device→host copy per array; arrays may be anything exposing DLPack or
  `__cuda_array_interface__` (CuPy, torch, Numba). Host DLPack exporters are
  adopted zero-copy without CuPy installed.
- Dtypes stay canonical float64/int64 unless explicitly downcast
  (`float32=True` / `int32=True`) — recorded in the warned `notes`, and int32
  index casts are range-checked, never wrapped.
- Deliberately **no `[gpu]` pip extra**: CuPy wheels are CUDA-version-specific
  (`cupy-cuda13x` / `cupy-cuda12x` / `cupy-cuda11x` / ROCm), so a pinned extra would break for
  most users; the install error names the wheel recipe instead. Kept out of
  `[all]` and `[interop]`. Docs: [`doc/gpu.md`](doc/gpu.md).
- CI covers the pure payload layer and the DLPack **host** round-trip only;
  the CUDA device path runs under gated tests (`importorskip("cupy")` plus a
  real device check) and is **not covered by public CI** — its lines read as
  uncovered on the Codecov patch check by design.
- Housekeeping: the `CLAUDE.md` version-bump checklist grew from six to
  **eight** files — the Julia (`bindings/julia/MeshioPlusPlus/Project.toml`)
  and R (`bindings/r/meshioplusplus/DESCRIPTION`) manifests added in v8.3.0
  carry versions too and are exactly the kind that drift unnoticed.

## v8.3.0 (2026-07-23)

**Julia and R bindings** — the next two languages of the scientific-computing
audience after Python, C, Fortran and JavaScript. Both sit on the **existing**
C API (`libmeshioplusplus`, the installed pure-C99 header), exactly as the
Fortran module does: no new C++ is written, and the C++/WebAssembly/Python core
is untouched.

- **Julia** — `bindings/julia/MeshioPlusPlus/`, docs
  [`doc/julia.md`](doc/julia.md). `ccall` into the installed shared library,
  discovered through `MESHIOPLUSPLUS_LIB` or the standard loader path. A `Mesh`
  wraps the opaque handle with a GC **finalizer** (unlike Fortran, which frees
  explicitly). Genuine **zero-copy borrows** — `points_ptr`,
  `connectivity_ptr`, `point_data_ptr`, … — are returned as a `MeshBorrow`
  that both keeps its owning mesh alive and records the mesh's mutation
  generation, so using one after a mutating call raises `BorrowError` rather
  than reading stale memory. That is the C header's rule 3 *enforced*, not
  merely documented.
- **R** — `bindings/r/meshioplusplus/`, docs [`doc/r.md`](doc/r.md). Plain
  `.Call` over R's own C API rather than Rcpp, keeping the dependency
  footprint at zero; the handle is an external pointer with a registered
  finalizer. **R is copy-only**: R vectors are R-managed, so the C API's
  zero-copy borrow cannot survive into R without ALTREP machinery that is out
  of scope, and every accessor copies. There is therefore no `_ptr` accessor
  at all — the 0-based reader is named `mio_connectivity_raw`, deliberately
  not `_ptr`, so nobody reads it as a borrow. R also has no native 64-bit
  integer, so `int64` arrays arrive as `double` (exact to 2^53, no `bit64`
  dependency), with the stored dtype reported in a `"dtype"` attribute.

Both cover the full Fortran surface: lifecycle, read/write (including the
selective-read options and file metadata), every setter and getter, named
regions, and all ~24 operations — including the ones returning an opaque C
result, which are always drained through `_take_mesh` so every mesh handed to
the caller owns its handle and no piece can dangle when its result is freed.

Both follow the Fortran module's two conventions verbatim, because Julia and R
are column-major too: points shaped `(dim, n)` and connectivity
`(nodes_per_cell, n)` are the **same memory** as the C API's row-major shapes,
so **nothing is transposed anywhere**; and connectivity is **1-based**, with
the ±1 shift applied inside the copying accessors only. Index maps and
permutations shift the same way, with the C API's `-1` "pruned / absent"
sentinel becoming `0` — never a valid 1-based index. `partition_labels`
returns part *ids* rather than indices and is deliberately left unshifted.

> **Licence exception, deliberate and the one thing to know:** the Julia
> binding in `bindings/julia/` is **not MIT**. It is released under the
> **GNU General Public License, version 3 (GPL-3.0)** — a copyleft license,
> not a permission-required one: anyone, including a company, may use,
> modify or sell it commercially with **no permission needed**; the
> condition is on *conveying* (distributing) it — a distributed copy or
> modified version must be under GPL-3.0 too, with source available.
> Purely private/internal use that is never distributed carries **no
> obligation** at all. Because GPL-3.0 **is** OSI-approved, the package is
> eligible for Julia's General registry (registration itself is a separate
> follow-up, not done yet), and it still installs by path or URL in the
> meantime; there is also no BinaryBuilder JLL yet. Everything else — the
> C++ core, the C API the binding calls, and the R binding — remains MIT.

Neither binding invents a workaround for the C ABI's documented gaps (point and
cell sets beyond regions, the `frozen` pin masks, per-cell-type counts in the
statistics report, ragged block connectivity, the combined `data_manage`); each
repeats the same list in its own README and doc page. CI gains a `julia` and an
`r` job, both mirroring the existing external-consumer smoke test: build the C
API, `cmake --install` it, then consume it exactly as a user would — `Pkg.test`
and `R CMD check --as-cran` respectively — and each also asserts that a missing
library fails with a message naming how to build one.

**Julia and R example notebooks** — `example/julia/*.ipynb` and `example/r/*.ipynb`, the same
three-notebook tour (`01_read_and_visualize`, `02_convert_and_inspect`, `03_mesh_operations`) as the
existing C++ notebooks, called through the new bindings on their own Jupyter kernels ([IJulia](https://github.com/JuliaLang/IJulia.jl)
/ [IRkernel](https://irkernel.github.io/)). Since the flat C API these bindings ride on can't drive the SVG
writer's per-call data-driven colouring either (a gap the C++ notebooks don't hit, since they call
`write_svg` directly), quality/field renders that the C++ tour shows as a coloured mesh are shown as a
small chart instead — hand-rolled SVG bar/histogram charts in Julia, plain base-R graphics in R, neither
needing a plotting-library dependency. Writing these caught two real, now-fixed defects: `smooth()`'s `mu`
default had been hardcoded `-0.53` in **both** new bindings — an invented value, never checked against the
real default (`-0.34`, matching Fortran and the Python bindings) — which is wide enough to reliably degrade
a tangled mesh rather than recover it; and R's data setters were found to always write `Float64` regardless
of the R vector's storage mode, meaning `mio_split(by = "region")` cannot be driven by a tag built fresh in
R (only one already present in a read file, or produced by the C++ core itself), now a documented `doc/r.md`
gap distinct from Julia, which has no such restriction. Building the C API with HDF5 on for Julia notebooks
also surfaced a real Debian/Ubuntu + Julia interaction — an IJulia kernel's `dlopen` of a
`libhdf5_openmpi`-linked library can fail on a `libcurl` symbol-version mismatch against Julia's own bundled
`libcurl` — documented in `doc/julia.md` and fixed in the `julia` CI job (HDF5 off, not just netCDF).

## v8.2.0 (2026-07-23)

**In-memory interoperability with the wider ecosystem, without a file
round-trip.** A new Python-only module,
`src/python/meshioplusplus/_interop.py`, converts a `Mesh` to and from the
in-memory objects of its main consumers, sharing the underlying numpy buffers
wherever the target accepts them as they are. Nothing here touches the
C++/WASM/C/Fortran core, which stays dependency-free — this is pure Python over
the numpy the readers already return.

- **PyVista** — `to_pyvista(mesh, zero_copy_only=False)` /
  `from_pyvista(grid)`. Builds the VTK 9 `connectivity`/`offsets`/`celltypes`
  triple from the mesh's blocks in order, so mixed-type meshes are the normal
  case. This exists because PyVista's own `from_meshio` targets the upstream
  `meshio` package and does not recognize `meshioplusplus`.
- **trimesh** — `to_trimesh` / `from_trimesh`. trimesh holds triangles only, so
  non-triangle input is routed through meshio++'s *existing* operations rather
  than a reimplementation: volume meshes through `extract_surface`,
  higher-order cells through `convert_cells("linearize")`, quads and polygons
  through `convert_cells("simplexify")`.
- **Apache Arrow / Parquet** — `to_arrow` / `from_arrow` / `write_parquet` /
  `read_parquet`, plus a `meshioplusplus data export IN OUT.parquet
  [--location point|cell]` sub-verb in the Python CLI's `data` group (the
  native C++ CLI has no counterpart). This is a **tabular export of data
  arrays for analytics — not a mesh format**: it moves `point_data`/`cell_data`
  into pandas/polars/DuckDB, does not round-trip geometry, and is deliberately
  *not* registered in the format registry, so
  `meshioplusplus convert mesh.vtu out.parquet` does not work. Multi-component
  arrays become Arrow `fixed_size_list` columns rather than `_0`/`_1`/`_2`
  suffix columns, which would lose the shape, and the mesh's counts, cell
  types, version, location and region names ride along in the schema metadata.

**The zero-copy contract.** Buffers are shared when the target accepts the array
as-is (contiguous, supported dtype, right shape). Every `to_*` takes
`zero_copy_only`: `False` (default) records each copy in a `notes` list surfaced
as a warning; `True` raises instead, naming the array and the reason. It governs
arrays that exist in the `Mesh` — `points`, each data array, and single-block
connectivity — but not *derived* ones: VTK's `offsets`/`celltypes` have no
meshio++ counterpart and a multi-block mesh has no single connectivity array to
share, so making those fatal would reject every mixed-type mesh. Copies are
forced by an int32→int64 connectivity widening, a `wedge` block (whose meshio
node order differs from VTK's), a 2-D mesh's point padding, and any
cross-block concatenation. Returned wrappers hold references to every shared
array, so they stay valid after the source mesh is garbage-collected.

**Regions per target.** `Point` and `Cell` regions export as int8 `region:<name>`
mask arrays, plus a JSON sidecar (`meshioplusplus:regions` in PyVista's
`field_data`, trimesh's `metadata`) carrying each region's `dim`/`tag`, which a
mask alone cannot express — so a gmsh physical group's integer tag survives a
PyVista round-trip exactly. `Side` regions are dropped with a warning naming
them: neither target has a `(cell, local facet)` concept. For trimesh, regions
follow the composed operations' own documented behaviour and get no second
policy — a pure-triangle input keeps them, anything routed through
`extract_surface` loses them.

**New extras**, all pip-friendly and all Python-only: `[pyvista]`, `[trimesh]`,
`[arrow]`, and the aggregate `[interop]`. They are deliberately **not** in
`[all]`, which means "the optional deps the *formats* need" (h5py/netCDF4);
`[viewer]`/`[kahip]`/`[codecs]` stand apart for the same reason. Every
third-party import is lazy and raises a named `pip install meshioplusplus[...]`
error when absent.

**Open3D and DOLFINx are deferred to Phase 2.** `has_open3d()`/`has_dolfinx()`
ship returning `False` and `to_open3d`/`to_dolfinx` raise `NotImplementedError`
naming the phase. The constraints are recorded now in
[`doc/interop.md`](doc/interop.md): Open3D's `Vector3dVector` typically copies
(so the zero-copy claim differs per structure) and its wheel is ~400 MB, while
DOLFINx needs a single-cell-type mesh, a `ufl`/`basix` domain, an MPI
communicator and a VTK→basix node-ordering permutation, and is conda/apt-only so
it can never be a pip extra.

Docs: [`doc/interop.md`](doc/interop.md).

## v8.1.0 (2026-07-23)

**Named groups of entities are now a first-class part of the core.** A `Region`
— a named group of *points*, *cells* or, for the first time, *cell facets* —
lives in `meshioplusplus::Region`, is carried by all three mesh backends, crosses
the pybind11 boundary natively, and is visible from the C API, Fortran,
WebAssembly and the native CLI. This is Phase 1 of unifying what every format
spells differently: gmsh physical groups, Exodus blocks and sets, Abaqus
`*NSET`/`*ELSET`/`*SURFACE`, MED families, UNV groups, Ansys components,
OpenFOAM patches, Kratos SubModelParts.

What round-trips now:

- **Gmsh** (`.msh` 2.2) — physical groups map to `cell` regions carrying their
  dimension and their integer tag, derived from the `gmsh:physical` cell_data
  and `$PhysicalNames` on read and synthesized into both on write. A mesh whose
  groups came from another format now writes real physical groups.
- **Abaqus** (`.inp`) — `*NSET` → point regions, `*ELSET` → cell regions and
  `*SURFACE` → **side** regions, in both the C++ core and the Python reference.
  The C++ reader gained full parity with the Python one along the way
  (`GENERATE`, set-of-set references, `ELSET=` on an `*ELEMENT` line,
  `*INCLUDE`), so a set-carrying `.inp` no longer falls back.
- **Side sets have no precedent anywhere in meshio++** — they are only reachable
  through `.regions`, not through `point_sets`/`cell_sets`.

`point_sets` and `cell_sets` keep working exactly as they always have. They are
now views over the mesh's regions: reading one materializes the historical shape,
writing one creates or replaces the matching regions. Two accommodations make
that lossless rather than nearly so — a `cell_sets` entry that is not cell-index
data (gmsh's `gmsh:bounding_entities`, whose entity tags are signed) is kept
verbatim through a passthrough, and a `None` block becomes an empty array, which
every consumer already treated identically.

Also in this release:

- **`detail/region_remap.hpp`** — one shared remapper replaces nine near-identical
  per-operation Python implementations. `crop`, `split`, `merge` (with source-id
  namespacing), `reorder`, `clean`, `partition`, `convert_cells`, `refine` and
  `decimate` remap regions in C++; `transform`, `smooth`, `interpolate` and the
  data operations pass them through; `slice`, `isosurface` and surface extraction
  drop them with a warning naming the operation, never silently.
- **`detail/cell_index.hpp`** — the single owner of the global (block-major) cell
  index that `cell` and `side` regions are defined against.
- **KRATOS backend**: regions materialize as SubModelParts, taking precedence
  over the names the integer-tag inference would have claimed.
- **C API**: `mio_regions_create` / `_count` / `_name` / `_info` / `_entries` /
  `_free` plus `mio_mesh_add_region`. **Fortran**: `m%regions(...)` and
  `m%add_region(...)`. **WebAssembly**: regions travel on the mesh object, so
  `readMesh` / `writeMesh` / `convert` carry them with no new call.
- **Native CLI**: `info` prints point/cell/side sets, closing its documented
  omission, and `diff` compares regions in the core rather than only in the
  Python shim.

Deferred to Phase 2, explicitly: Exodus blocks/node sets/side sets, MED families
and groups (absorbing `MedInfo`), UNV and Ansys (absorbing `UnvInfo`/`AnsysInfo`),
OpenFOAM boundary patches, XDMF Sets, and VTU/VTP — which have no native set
concept, so a convention has to be chosen rather than invented silently. A
`regions` CLI group and a region-aware `split --by region` are deferred with them.

Two behaviour changes worth calling out. A remapped set now comes back **sorted**
rather than in the permutation's own order, because region entries are canonical
(sorted, de-duplicated) so that membership comparison is exact. And an `*ELSET`
declared inside an `*INCLUDE`d Abaqus file is now carried through the merge; it
used to be dropped under a TODO in `_abaqus.py`.

## v8.0.0 (2026-07-23)

**The WebAssembly build now ships every format the C++ core has.** `cgns`,
`h5m`, `hmf`, `med` and `exodus` — the five that need HDF5 or netCDF — are
readable and writable from `@meshioplusplus/wasm`, as is XDMF's `Format="HDF"`
data path. There is no longer a WASM-specific format gap: 41 formats, 40
readable, 41 writable. `availableFormats()` reports them, and every wasm-written
file was verified to read back correctly through the native h5py/netCDF4-backed
Python package.

- New `build/build-wasm-deps.sh` source-builds a wasm32 **HDF5 1.14.6** and
  **netcdf-c 4.9.3** (pinned, SHA256-checked) into a self-contained prefix.
  CMake still never downloads anything — it only *finds* the result, via
  `CMAKE_FIND_ROOT_PATH`. `build/configure-wasm.sh` runs it automatically the
  first time and gained `--with/--without-hdf5`, `--with/--without-netcdf` and
  `--deps-prefix`. `--without-hdf5` reproduces the old, smaller artifact.
- **Breaking:** writing `.xdmf` from the WASM build now emits an HDF companion
  `<base>.h5` beside the XML instead of inlining the heavy data, because the
  registry's XDMF writer default follows the build — the same rule every native
  build already obeyed. A JS caller must pull **two** files out of the virtual
  filesystem, not one. Reading all three XDMF data formats is unaffected.
- **Breaking:** the published `.wasm` grows from ~2.3 MB to ~5.5 MB (statically
  linked libhdf5 + libnetcdf). Build with `--without-hdf5` for the small
  artifact; the JS API is identical either way.
- MED **cannot write named fields** in the WASM build: the C++ MED writer defers
  a mesh carrying `point_data`/`cell_data` to the Python reference writer, and
  this build has no Python to defer to, so it throws by name. MED geometry,
  `point_tags`/`cell_tags` and families write normally.
- Fixed a latent, silent **stack overflow** in the WASM build, found while
  adding the above: HDF5's and netCDF-4's frames overrun Emscripten's default
  64 KiB stack, which grows down into the static data segment. One Exodus write
  clobbered libc++'s locale facets, after which every `istream >> number` in the
  module — i.e. every ASCII reader, `gmsh`/`obj`/`off`/`vtk` included — trapped.
  The wasm target now links with `-sSTACK_SIZE=4MB` and
  `-sSTACK_OVERFLOW_CHECK=1`, so a recurrence aborts loudly instead of
  corrupting unrelated state.
- CI: `wasm.yml` now also builds and smoke-tests on pull requests touching the
  wasm surface or `src/cpp/` (it was tag-only), with the dependency prefix
  cached.

## v7.16.0 (2026-07-22)

New **`decimate`** operation — reduce a surface mesh's face count by greedy
quadric-error-metric (Garland–Heckbert) edge collapse, preserving shape,
boundaries and features: the resolution-*reducing* inverse of `refine`,
completing the pair. Surface meshes only (`quad`/`polygon` blocks are
triangulated first, so the output is all-triangle with the block structure
kept 1:1); a volume mesh raises by name pointing at `extract_surface`.

- Stopping criteria (exactly one): `ratio` (fraction of faces to keep),
  `target_faces` (absolute, within one collapse), or `max_error` (collapse
  while the cheapest quadric error is below it).
- Placement `optimal` (quadric minimizer, midpoint fallback when the 3x3
  system is ill-conditioned) / `midpoint` / `endpoint`. Float `point_data` is
  blended along the collapsed edge (clamped parameter); integer arrays keep
  the survivor's value.
- Boundary vertices (once-used-edge test) and feature vertices (face normals
  differing by more than `feature_angle`, default 30°) are pinned by default;
  an optional `frozen` mask pins more. The link condition and a normal-flip
  guard reject any collapse that would change topology, create a non-manifold
  edge, or fold the surface — rejections are counted in the report.
- Deterministic: parallel setup with pinned FP order, serial greedy loop with
  a total heap ordering. Output is byte-identical across the three mesh
  backends, thread counts, and the C++/numpy-fallback boundary
  (`test_cpp_matches_python`).
- Exposed everywhere: Python `meshioplusplus.decimate` (with `return_report`),
  C API `mio_decimate` + opaque `mio_decimate_result` (maps + counters;
  `frozen` is a documented flat-ABI gap), Fortran `m%decimate(...)`, WASM
  `decimate(...)` (also a `convertSurfaceOps` pipeline op, so the browser
  viewer's worker can reach it), and a `decimate` verb in both CLIs.
- Docs: `doc/decimate.md`, README "Decimation" section, notebook demo cell.

Also: fixed two path strings in `CLAUDE.md` corrupted by the repository
restructure's `cpp/` → `src/cpp/` rewrite, and its "five version files"
sentence (there are six).

## v7.15.0 (2026-07-22)

New **`isosurface`** operation — the level set of a scalar field: the locus where
a `point_data` array equals a given isovalue, as a mesh one topological dimension
below the cut cells (a 3D volume mesh yields a `triangle`/`quad` surface, a 2D
surface mesh a `line` contour). This is the **data-driven sibling of `slice`**:
slice cuts where `dot(x - origin, normal) = 0`, isosurface where
`f(x) - isovalue = 0`.

- `isosurface(mesh, array, isovalues, component=, record_parent_ids=)`: the field
  must be `point_data` — `cell_data` is piecewise constant, so there is no
  crossing to locate and no level set to draw; naming one raises, pointing at
  `cell_data_to_point_data` (`meshioplusplus data to-point`) as the fix. A
  multi-component array reduces to `component`, or to the row magnitude when that
  is unset.
- **Several isovalues land in one mesh**, cut in ascending order (sorted, exact
  duplicates dropped) and concatenated with the section blocks merged by cell
  type, so a single-isovalue call has exactly slice's block structure. Each
  contour cell carries a Float64 `iso:value` (the level) and an Int64
  `iso:index` (its ordinal) — the latter because `split`'s tag criterion needs an
  integer array, which makes `split --by region --tag iso:index` the
  one-mesh-per-contour recipe.
- **The contoured field reads back as exactly the isovalue** on the cut points;
  every other `point_data` array is interpolated at the crossing (Float64, exact
  for a linear field). The exception is a magnitude-reduced multi-component
  array, where `|lerp(v)| != lerp(|v|)` mathematically and the value stays
  approximate.
- Degeneracy rule, uniform with slice: a node whose value is exactly the isovalue
  counts as being on the **positive** side, so a plateau lying at the isovalue
  emits its boundary once, not twice. Contours are watertight (crossings on
  shared edges dedupe to one node) and wound toward increasing field. An isovalue
  outside the field's range is an empty contour, not an error.
- `record_parent_ids` attaches an Int64 `iso:parent_cell`; each contour cell
  inherits its parent's `cell_data`. Contour points are all new, so
  `point_sets`/`cell_sets` are not carried.
- Output is byte-identical across the three mesh backends, thread counts and the
  C++/numpy boundary. Exposed on every binding surface (pybind `_core`, C API
  `mio_isosurface`, Fortran `m%isosurface`, WASM `isosurface`), as the CLI verb
  `isosurface IN OUT --array NAME --values v1,v2 [--component I]
  [--record-parent-ids]` in both CLIs, and as an operation chip in the browser
  viewer. Docs: [`doc/isosurface.md`](doc/isosurface.md).

Fixed: `src/viewer/package-lock.json` recorded `@meshioplusplus/wasm@7.10.0`
while the package itself was at 7.14.0, which makes `npm ci` hard-fail
(`does not satisfy`) and so broke the viewer jobs in `ci.yml` and `docs.yml`.
The lock's `"../wasm"` version is now part of the version bump — see the
"Version bumps" section of `CLAUDE.md`.

Internal, **no behaviour change**: slice's marching-tetrahedra cutter — the
simplexify, the sign-mask case table, the watertight edge dedup, the winding, the
degeneracy rule and the `point_data`/`cell_data` carry — was hoisted verbatim out
of `operations/slice.cpp` into the shared `detail/marching.hpp` (and out of
`_slice.py` into `_marching.py`), the way `spatial_hash.hpp` was hoisted from
merge in v7.13.0 and `space_filling.hpp` from reorder in v7.6.0. Slice's output is
byte-identical across the hoist and its test suites are unchanged; the only
parameterized rule is the winding, which slice resolves against its fixed plane
normal and isosurface against the local field gradient.

## v7.14.0 (2026-07-22)

New **`slice`** operation — the planar cross-section of a mesh: the actual
intersection of the mesh with a plane, one topological dimension below the cut
cells (a 3D volume mesh yields a `triangle`/`quad` surface, a 2D surface mesh a
`line` mesh). Unlike `crop` (plane mode), which keeps whole cells on one side,
`slice` computes the intersection and lowers the dimension.

- `slice(mesh, origin=, normal=, record_parent_ids=)`: robust marching
  tetrahedra on a simplexified input (every 3D cell becomes a tetra, every 2D
  cell a triangle), so each cell's cross-section is a well-defined convex
  primitive — a hex/wedge section is therefore the union of its simplices'
  sections. The signed-distance crossing `t = d_i/(d_i - d_j)` is computed from
  the sorted edge endpoints and deduped by that edge key, so shared edges yield
  a single output node (the section is watertight). Degeneracy rule: a node on
  the plane (`d == 0`) is classified positive, making the sign mask total, so a
  plane grazing a shared face is emitted exactly once (no double emission);
  collapsed primitives are dropped. Section faces are wound so their Newell
  normal points toward the `+normal` side. Each section cell inherits its
  parent's `cell_data`; `record_parent_ids` attaches an Int64
  `slice:parent_cell`. `point_data` is interpolated at the cut (Float64 output);
  `point_sets`/`cell_sets` are not carried (the section is new topology). Output
  is byte-identical across the three mesh backends, thread counts, and the
  C++-core/numpy-fallback boundary.
- Exposed on every surface: pybind `_core` + the numpy fallback (`slice` shadows
  the built-in only as a module attribute), C API `mio_slice`, Fortran
  type-bound `m%slice`, WASM `slice`, and the `slice IN OUT --origin --normal`
  verb in both CLIs.
- The browser viewer's planar "section" now routes through `slice` (the true
  cross-section) instead of the previous crop-half-space + re-skin.
- Docs: new `doc/slice.md`, CLI reference entry, README "Slicing /
  cross-sections" section, and a notebook demo in both
  `example/python/03_mesh_operations.ipynb` and the C++ mirror.

## v7.13.0 (2026-07-21)

New **`interpolate`** operation — cross-mesh field transfer, the first two-mesh
operation that moves data (diff compares, merge concatenates; neither
resamples): sample a source mesh's data arrays onto a target mesh, returning a
copy of the target with its geometry, connectivity, own data and sets preserved
exactly.

- `interpolate(source, target, method=, arrays=, extrapolate=, default_value=,
  on_conflict=)`: source `point_data` sampled at the target's points,
  `cell_data` by nearest source-cell centroid (always, whatever the method).
  `method="nearest"` (default) copies the nearest source point's value
  bit-for-bit (dtype-preserving); `method="barycentric"` simplexifies the
  source first and interpolates linearly — exact on a linear field, Float64
  output, with `default_value`/`extrapolate` covering target points outside the
  source domain. `on_conflict` is `error`/`overwrite`/`suffix` (`name +
  "_interp"`). Output is byte-identical across the three mesh backends, thread
  counts, and the C++-core/numpy-fallback boundary.
- Exposed on every surface: pybind `_core` + the numpy fallback, C API
  `mio_interpolate` (arrays as `char**` + count, `NULL`/`<= 0` = all
  point_data), Fortran module-level `mio_interpolate`, WASM `interpolate`, and
  the `interpolate SOURCE TARGET OUT` verb in both CLIs.
- New shared `detail/spatial_hash.hpp`: merge's weld bucket grid hoisted
  verbatim (merge's output stays byte-identical) and extended with the
  expanding-shell / box-insert queries interpolate needs.
- Docs: new `doc/interpolate.md`, CLI reference entry, README "Field transfer"
  section, and a notebook demo (`example/03_mesh_operations.ipynb`).

## v7.12.0 (2026-07-21)

The OFF reader/writer (C++ core and Python reference) gain **quad and polygon
face support** — [issue #35](https://github.com/loumalouomega/meshioplusplus/issues/35)
reported that a valid OFF file using quad faces was rejected outright with
"Can only read triangular faces". OFF's own spec allows faces of any vertex
count; only this implementation (and, as it turns out, upstream meshio too)
had hard-coded the triangle-only assumption.

- `read_off` now groups faces by vertex count into `triangle` (3), `quad` (4),
  or `polygon` (else) cell blocks, exactly like the sibling OBJ reader in the
  same file: a run of same-count faces stays in one block until the count
  changes. A leading count below 3 remains a hard `ReadError`.
- `write_off` now writes every `triangle`/`quad`/`polygon` cell block (in mesh
  order); any other cell type is skipped with a warning instead of silently
  dropped. A `polygon` block written by the C++ path must be rectangular; the
  Python reference writer also accepts a ragged `polygon` block.
- New fixtures `tests/meshes/off/cube_example.off` (6 quad faces) and
  `cube_example_as_triangs.off` (the same cube pre-triangulated) back the
  regression tests.

## v7.11.0 (2026-07-21)

The SVG and TikZ writers gain **data-driven colouring**: a `color_by` scalar
field mapped through a built-in colormap to per-face fills, with an optional
colorbar. They already rendered 3D meshes by projecting the extracted skin;
until now every face got the same flat fill, so the vector figures could show
shape but never a field. This is the vector complement to `screenshot()`
(v7.8.0), which covers the raster side.

- **`color_by` / `component` / `cmap` / `vmin` / `vmax` / `nan_color` /
  `colorbar`** on `write_svg` / `write_tikz`, the pybind bindings, the Python
  shims and both CLIs' `convert`. Appended after the camera arguments, in that
  order, everywhere.
  - **Point data** colours a face by the mean of its corner values; **cell
    data** by its owning cell's value. For a projected volume mesh the owner is
    found through the `"surface:parent_cell"` provenance, so a per-cell material
    or quality metric lands on the right skin facet.
  - Multi-component arrays reduce to `component` or to their row magnitude.
  - The range defaults to the finite range of the **drawn faces** — so the
    visible figure spans the whole colorbar. This differs from ParaView, which
    ranges over the whole array. Non-finite values are excluded from the range
    and drawn in `nan_color`.
  - `--color-by NAME [--component I] [--cmap …] [--vmin V] [--vmax V]
    [--nan-color C] [--colorbar]` on `meshioplusplus convert`, in **both** the
    Python and the native CLI, rejected for any output but `.svg`/`.tikz`.
- **Built-in colormaps** (`detail/colormap.{hpp,cpp}` + its `_colormap.py`
  twin, generated by `tools/gen_colormaps.py`): viridis, coolwarm and turbo as
  256-entry uint8 LUTs. **No new dependency** — matplotlib is needed only to
  *regenerate* the tables, never to use them. Storing full 256-entry tables
  rather than interpolated control points is what keeps the C++ and Python
  writers byte-identical: it removes every floating-point interpolation from the
  colour path, leaving a single index expression. The C++ side is
  declaration/definition-split like `projection.hpp`: the header declares the
  API only, and the table data plus function bodies live in `colormap.cpp`.
- **Colouring is a documented gap on the flat bindings.** The C API, Fortran and
  WebAssembly surfaces reach SVG/TikZ through the shared registry, whose
  `(path, mesh)` writer entries structurally cannot carry parameters, so they
  keep emitting the fixed default styling — as with the point/cell-set gaps in
  `diff`/`merge`/`split`.
- **`tools/gen_doc_images.py`** regenerates the committed doc/README figures
  from the bundled sample mesh: the coloured SVGs via this feature, the shaded
  screenshots via `screenshot()`. It complements the PyVista notebooks rather
  than replacing them — those still own the executable demonstrations.
- **Behaviour change (cosmetic):** the SVG writer's *default* stroke width now
  prints as `1` rather than `1.0`. The C++ core always emitted `%g`; the Python
  reference used `str()`, and that single token was the only thing keeping the
  two from being byte-identical. Fixing it lets `tests/test_svg.py` assert full
  byte equality, as `tests/test_tikz.py` already did — which is what now guards
  every fill, the `<style>` block and the colorbar. Rendered output is
  unaffected. Pass `stroke_width=` explicitly to pin an exact value.

With `color_by` unset, SVG and TikZ output is byte-identical to v7.10.0, and the
flat 2D path is physically untouched.

## v7.10.0 (2026-07-21)

The native CLI gets `view` and `screenshot`, restoring verb parity with the
Python one — which v7.9.0 had broken without saying so. The browser viewer
gains click-to-inspect.

- **`view` and `screenshot` in the C++ binary**, backed by
  [Polyscope](https://polyscope.run) (MIT), vendored as a git submodule.
  Volume rendering and slice planes: the things a surface renderer, and so the
  browser viewer, structurally cannot do.
  - **Optional and off by default** (`-DMESHIOPLUSPLUS_WITH_POLYSCOPE=ON`,
    `build/configure.sh --with-polyscope`). The prebuilt release binaries do
    **not** include it and are unchanged — their whole point is being
    dependency-free single files, and Polyscope needs OpenGL, GLFW and X11.
  - It attaches to the **CLI target only, never to the core**. Unlike KaHIP —
    which lives in the core because partitioning is a core operation — viewing
    is not, so `_core`, the C API, the Fortran module and the WebAssembly build
    cannot acquire an OpenGL dependency through it.
  - The verbs exist in **every** build and report the flag when it is off,
    rather than silently not existing.
  - Polyscope vendors its own submodules, so enabling it needs
    `git submodule update --init --recursive`.
- **Click-to-inspect in the browser viewer.** An **Inspect** toggle reports a
  clicked cell's id and type, every cell and point array's value there (all
  components), the nearest vertex, and the originating volume cell for a solid
  — with the picked cell outlined. On click only, never on hover, and disabled
  above two million cells rather than made slow.
- **Fixed:** `CLAUDE.md`'s "documented gaps vs the Python CLI" list omitted the
  two verbs v7.9.0 added, and its verb list was a release out of date.

Internal: `gather_cell_data_onto_surface` moves from the WebAssembly binding
into `operations/surface.cpp`, since the CLI now needs it too. The mesh →
Polyscope mapping is deliberately free of Polyscope headers, so it compiles and
is tested in the default build with no GL.

## v7.9.0 (2026-07-21)

The browser viewer stops being a generic vtk.js app and starts running meshio++
itself. It previously used 6 of the 40 functions the WASM package exposes.

- **Mesh operations in the browser.** Quality, clean, smooth, refine, partition
  and a sectioning cut, applied to the mesh you opened and re-rendered in
  place, with no server and no upload. They compose, each shows as a chip you
  can remove, and **undo is exact**: the worker keeps the original file bytes
  and replays the remaining pipeline, so nothing needs an inverse and nothing
  accumulates rounding.
- **WASM: `convertSurfaceOps`.** One binding applies an operation pipeline and
  writes the renderable surface, all inside C++. Chaining the individual
  operation bindings would route the mesh through the flat JS representation on
  every step and destroy every multi-component array — the thing
  `convertSurface` exists to prevent. An empty pipeline is byte-identical to
  `convertSurface`, so the plain and post-operation display paths cannot drift.
- **Viewer polish**: a DOM colour legend with an editable range replacing
  `vtkScalarBarActor`, an orientation cube whose faces snap the camera,
  surface/wireframe/points, an opacity slider, and Fit/PNG buttons.
- **`viewer/` is now TypeScript**, with `tsc --noEmit` in CI. The shared worker
  protocol means a mismatch between what the client sends and what the worker
  handles is a compile error rather than a runtime surprise.
- **A TikZ icon set** under `icons/`, built with the same `pdflatex` +
  `dvisvgm` pair `logo/build.sh` uses, generated into a typed module so every
  icon reference is checked.
- **The offline page carries results.** `view(backend="browser")` gained
  `quality=True` to bake in per-cell metrics, always embeds the volume mesh's
  geometric statistics (the page renders only the boundary, so it cannot derive
  them), and `color_by` now works instead of warning that it does not.

Fixed:

- Point-data colouring used only half the colormap. `interpolateScalarsBefore
  Mapping` renders a range spanning zero entirely in the warm half in vtk.js
  32.9.0 — measured at 0 blue pixels against 275k red ones on a symmetric
  field — so it is now off, with the evidence recorded beside the setting.
- Surface edges z-fought into faint dashes: the polygon offset was on the body
  mapper, but a negative offset moves *toward* the viewer, so the surface was
  drawn in front of its own wireframe.
- The `kahip` CI job failed at its Python step with `ImportError: libkahip.so`.
  scikit-build-core strips the RPATH from the installed extension so the wheel
  stays relocatable, so the loader needs the prefix on its path;
  [`doc/partition.md`](doc/partition.md) now warns about this for users too.

## v7.8.0 (2026-07-21)

meshio++ can now *show* you a mesh. One entry point, `view()`, with two
backends — a native desktop window and a browser — plus a hosted demo that
doubles as a client-side format converter.

- **`view(mesh, backend=...)` — interactive visualization.** `"polyscope"`
  opens a native window; `"browser"` renders with vtk.js, inline in a notebook
  or in your default browser; `"auto"` picks polyscope when it is installed and
  a display is available. Also `screenshot(mesh, path)`, which renders
  headlessly and so works from CI and a docs build, and `has_viewer()`. New CLI
  verbs `meshioplusplus view` and `meshioplusplus screenshot`.
  - **Polyscope is a Python-only optional dependency**, behind a new
    `[viewer]` extra (`pip install meshioplusplus[viewer]`). It never reaches
    the C++/WASM/C/Fortran core, and the browser backend needs nothing from it.
    A missing install raises naming the command that fixes it.
  - The mesh → renderer mapping is pure and separately tested: no renderer
    import, no display, no mutation of the input. Volume meshes route through
    `convert_cells(simplexify)` only where they must — polyscope holds
    tetrahedra and hexahedra directly, so a hexahedral mesh keeps its
    hexahedra. Every lossy step is reported rather than done quietly.
- **A browser viewer at `viewer/`**, deployed to GitHub Pages alongside the
  docs, that consumes the published `@meshioplusplus/wasm` package exactly as
  an external user would. Drag in any of the ~36 formats the WASM build
  supports, colour by point or cell data with a scalar bar, and convert and
  download to any writable format. Everything runs client-side: no server, no
  upload. The same bundle, built without the WASM, ships in the wheel as the
  offline render path for `view(backend="browser")`.
  - It renders **VTP**, not VTU, and shows a volume mesh by its boundary:
    vtk.js has no unstructured-grid model at all.
- **WASM: `availableFormats()`** returns the reader and writer names this build
  actually supports, so a consumer no longer has to hardcode a table that
  drifts from the build. **`convertSurface()`** reads, extracts the boundary,
  linearizes and writes in one call without materializing a JS mesh — which is
  what keeps multi-component (vector/tensor) arrays, since the flat JS mesh
  representation cannot carry them. Boundary facets also now inherit their
  owning cell's data.
  - Fixed: `index.d.ts` referenced a `MeshMetadata` type it never defined,
    a TS2304 for any consumer without `skipLibCheck`.

## v7.7.0 (2026-07-21)

A new dependency-free mesh operation that improves element *shape* in place — the
counterpart to `quality`, which only measures it, and the complement to `refine`, which
changes resolution without changing shape.

- **`smooth` — relax point coordinates toward their edge-neighbour centroids.** A pure
  coordinate move: connectivity, `cell_data`, `field_data` and `point_data` values all pass
  through unchanged, and the points array keeps its input dtype.
  - Two operators. **Laplacian** (`x <- x + lambda*L(x)`) smooths strongly but shrinks;
    **Taubin** alternates a `+lambda` pass with a larger `-mu` pass and does not, which is
    why it is the default. On a jittered 8x8 quad grid over 40 iterations, Laplacian
    contracts the bounding box by 57% where Taubin contracts it by 3.6%.
  - The neighbour graph is **edge** adjacency, not the element clique, so a structured
    hexahedron block is a fixed point rather than being bevelled toward a sphere.
  - Boundary nodes are pinned by default (`fix_boundary`), as are geometric corners and
    creases (`preserve_features`, `feature_angle_deg`, default 30°), nodes named in an
    optional `frozen` mask, and nodes of blocks whose edge topology is unknown — the
    higher-order family, the VTK-Lagrange types and `custom` — since an unknown
    neighbourhood gives no defined smoothing target.
  - An **inversion guard** (on by default) rejects any move that would turn a valid cell
    inverted, counting the rejections. It is "do no harm", not "preserve the sign": a cell
    that arrives already inverted imposes no constraint, so smoothing can still repair a
    tangled region rather than locking the tangle in.
  - Deterministic by construction: Jacobi updates from the previous pass's positions,
    neighbour sums in ascending node id, no hashing or sorting in the update loop, and the
    boundary set built with `surface.cpp`'s serial-dedup phase split. Output is
    byte-identical across the MESHIO/NATIVE/KRATOS backends and across thread counts.
  - On every surface: Python `smooth`, C API `mio_smooth` (plain mesh plus nullable counter
    out-params, like `mio_clean`), Fortran `m%smooth`, WASM `smooth`, and a `smooth` verb in
    both CLIs.
- **Internal:** the node-adjacency CSR moves to `detail/node_adjacency.hpp` with a
  `Clique | Edge` kind and is now shared with `reorder`, which keeps `Clique` and whose
  output is unchanged.

## v7.6.0 (2026-07-20)

A new mesh operation for domain decomposition, plus the repo's first optional
partitioning dependency.

- **`partition` — decompose a mesh into N balanced parts.** The count-driven
  complement to the criterion-driven `split`. Two methods:
  - **SFC** (always available, the default fallback): cells are ordered along a
    Hilbert space-filling curve of their centroids (the same key transforms
    `reorder` uses, now shared via `detail/space_filling.hpp`) and cut into
    `nparts` contiguous ranges — equal-weight part sizes differ by at most one
    cell, and with `weights=<cell_data>` the cut follows the weight prefix sum.
    Deterministic and byte-identical across mesh backends, thread counts, and
    the C++-core/numpy-fallback boundary (pinned by a test).
  - **KaHIP** (optional, the quality path): the shared-face dual graph is handed
    to KaHIP's serial `kaffpa()` with configurable `imbalance` (default 0.03),
    `mode` (`fast`/`eco`/`strong`, default `eco` — eco/strong carry the edge-cut
    wins) and `seed`. Ported from the Kratos KaHIPApplication
    (KratosMultiphysics/Kratos#14453).
  - `partition_labels` returns just the block-aligned Int64 assignment
    (`partition:part`); `record_ids` attaches `partition:original_point_id`/
    `partition:original_cell_id` to each piece. Pieces keep the input block
    structure 1:1 (unlike `split`), so they recombine into the input.
    `ghost_layers` is reserved (raises) in v1.
  - On every surface: Python `partition`/`partition_labels`, C API
    `mio_partition`/`mio_partition_labels` (opaque result handle with zero-copy
    map getters), Fortran `m%partition`/`m%partition_labels`, WASM
    `partition`/`partitionLabels`, and a `partition` verb in both CLIs
    (`OUT_{part}.vtu` expansion, `--labels-only`).
- **New optional dependency `MESHIOPLUSPLUS_WITH_KAHIP`** (OFF by default; MIT
  like meshio++, so no licensing implication). Located via the new
  `cmake/FindKaHIP.cmake` (`KAHIP_ROOT` prefix / pkg-config) — never vendored or
  auto-downloaded, and only the serial `kaffpa` interface is linked (no
  ParHIP/MPI). The installed library's 32/64-bit index width is detected at
  runtime via `kahip_sizeof_idx()`. Conan `with_kahip`, vcpkg feature `kahip`,
  pip extra `meshioplusplus[kahip]` (the MIT `kahip` wheel, which also gives
  pip-only installs the quality path). Requesting `method="kahip"` without any
  KaHIP backend fails with an error naming the option — never a silent
  downgrade to SFC.

## v7.5.0 (2026-07-20)

A new dependency-free mesh operation that *increases* resolution — the counterpart to
`convert_cells`, which preserves it, and `crop`/`clean`, which reduce it.

- **`refine` — uniform mesh refinement.** Subdivides every cell into congruent children of
  the **same** cell type, one fixed template per type: `line` → 2, `triangle` → 4,
  `quad` → 4, `tetra` → 8, `wedge` → 8, `hexahedron` → 8. `levels=n` applies the templates
  `n` times.
  - New nodes sit at the midpoints of the parent's edges, quad faces and (hexahedron only)
    body, and carry the mean of that entity's corner values for every `point_data` array, so
    a linear field is interpolated exactly. Each parent's `cell_data` row is replicated to
    its children, and block structure is preserved 1:1.
  - **The refined mesh has no hanging nodes.** Mid-edge *and* quad-face-centre nodes are
    shared between every cell touching the entity — only the hexahedron body node is
    per-cell. (Per-cell face centres would leave two adjacent hexahedra referencing distinct
    coincident nodes, splitting the mesh topologically along every interior face.)
  - Children inherit the parent's orientation, so a well-oriented input refines to zero
    newly-inverted cells. Volume is conserved exactly for `line`/`triangle`/`quad`/`tetra`
    always, and for `wedge`/`hexahedron` when the parent is affine — for a general trilinear
    hexahedron the children's volumes do not sum to the parent's, which is a property of the
    geometry rather than of the implementation.
  - The tetrahedron's interior diagonal is fixed at the opposite-edge pair `(0,1)`–`(2,3)`
    for determinism only; being strictly interior, it does not affect conformity — unlike
    `convert_cells`' hex-simplexify diagonal, whose endpoints lie on the boundary.
  - Higher-order cells (linearize first), `pyramid` (whose uniform refinement is 6 pyramids
    + 4 tetrahedra, breaking the same-type contract), and ragged polygon/polyhedron blocks
    raise by name rather than being silently passed through, which would produce hanging
    nodes next to refined neighbours.
  - Output is byte-identical across the MESHIO/NATIVE/KRATOS backends and any thread count:
    the templates are fixed and the new-node numbering comes from a serial dedup pass over a
    parallel-filled buffer, never a concurrent hash insert.

  Exposed as Python `meshioplusplus.refine(mesh, levels=1, record_parent_ids=False)`, C
  `mio_refine` (+ the `mio_refine_result` handle), Fortran `m%refine(levels, ...)`, WASM
  `refine(...)`, and a `refine` verb on both the Python and native CLIs.

- **`wedge18` is now skinnable.** `cell_faces` gained the `wedge18` row (completing the
  family alongside the existing `hexahedron27` and `pyramid14` entries), so
  `extract_surface`, `extract_skin` and `compute_quality` now handle `wedge18` meshes
  instead of warning and skipping them. Mirrored in the pure-Python `_skin.py` twin.

- Internal: the per-type edge tables that `convert_cells`' `elevate` mode and `refine` both
  need are now owned by one shared `detail/cell_subdivision.hpp` (which delegates the 2D
  rows to `detail/cell_edges.hpp`) rather than being transcribed independently. No
  behaviour change.

## v7.4.0 (2026-07-20)

A new dependency-free mesh operation, plus a WebAssembly fix that makes every
*existing* geometry operation reachable from the published package for the first time.

- **`convert_cells` — convert a mesh's element representation.** A mesh operation (not a
  file format), dependency-free and available on every binding surface, with three modes:
  - `linearize` — every higher-order cell becomes its linear base (`tetra10` → `tetra`,
    `hexahedron27` → `hexahedron`, ...), keeping the corner connectivity verbatim and
    pruning the nodes that become unreferenced (connectivity, `point_data` and
    `point_sets` are remapped). Cell count is unchanged, so `cell_data` passes through.
  - `simplexify` — every cell is decomposed into simplices of the same topological
    dimension: quad → 2 triangles, polygon(n) → (n−2)-triangle fan, hexahedron → 6 tetra
    (a canonical Freudenthal fan around the main diagonal 0–6), wedge → 3 tetra,
    pyramid → 2 tetra. The children reuse the parent's own corner nodes, so no points are
    added, and each parent's `cell_data` row is replicated to its children. Higher-order
    input is linearized first. Every emitted simplex is positively oriented for a
    well-oriented input, and volume is conserved — both pinned by tests.
  - `elevate` — linear cells are promoted to their serendipity quadratic counterpart
    (`triangle` → `triangle6`, `hexahedron` → `hexahedron20`, ...), creating one new node
    per unique edge at the edge midpoint with `point_data` set to the endpoint mean. The
    full-Lagrange targets that need face/body centres (`quad9`, `hexahedron27`) are an
    explicit non-goal of this version and raise by name.

  All three modes are idempotent on cells they do not apply to, so they are safe on a
  mixed-order mesh, and output is byte-identical across the MESHIO/NATIVE/KRATOS backends
  and any thread count (the mid-edge numbering is assigned by a serial pass over a
  parallel-filled buffer, never a concurrent hash insert).

  Exposed as Python `meshioplusplus.convert_cells(mesh, mode=..., record_parent_ids=...)`,
  C `mio_convert_cells` (+ the `mio_convert_cells_result` handle), Fortran
  `m%convert_cells(mode, ...)`, WASM `convertCells(...)`, and a `convert-cells` verb on
  both the Python and native CLIs.

- **WebAssembly: every geometry operation is now reachable from `loadMeshioPlusPlus()`.**
  `wasm/src/index.mjs` previously forwarded only file I/O and the five `data_*` operations,
  so `extractSurface`, `extractSkin`, `attachQuality`, `sniffFormat`, `reorder`,
  `computeBandwidth`, `diff`, `meshesEqual`, `merge`, `transform`, `clean`, `cropBbox`,
  `cropPlane`, `split`, `stats` and `meshBackend` were bound in `js_bindings.cpp` but
  unreachable through the package's own API — the same class of bug fixed for the data
  operations in v7.2.1. All of them (plus the new `convertCells`) are now forwarded by the
  wrapper, declared in `wasm/index.d.ts`, and exercised **through the wrapper** by
  `wasm/test/smoke.mjs`, which is what would have caught the original breakage.

## v7.3.0 (2026-07-20)

Three additive I/O performance features. **Default behaviour is unchanged**: `read()`,
`write()` and every existing file are byte-for-byte as before, and the default build gains
no new dependency.

- **Selective / partial reads and `read_metadata()`** — read only what you need.
  - `read(path, points_only=True)` returns geometry (points *and* connectivity) with no data
    arrays; `read(path, arrays=["u", "v"])` returns only the named ones. `arrays=None` means
    every array and `arrays=[]` means none — a deliberate distinction, preserved all the way
    down to the C ABI. Names absent from a file are ignored, not an error.
  - `read_metadata(path)` summarizes a file — point/cell counts, per-block cell types,
    data-array names — without loading the heavy arrays.
  - **VTU, VTP, XDMF and Gmsh** skip the unwanted array bodies outright, and all four plus
    **Gmsh 4.1** have native header-only metadata paths (Gmsh 2.2 declines and falls back, since
    it stores a type per element). Every other format is
    read in full and filtered, which is correct but not faster; `read_metadata`'s
    `fell_back_to_full_read` says which happened, so a summary never implies a saving that did
    not occur. XDMF is the cheapest case (every `<DataItem>` declares its shape, so counts are
    exact without touching any payload, and on the HDF path without opening the `.h5` at all).
  - Honest scope note: for VTU/VTP the file is still read and XML-parsed — pugixml always
    materializes PCDATA. What is skipped is base64 decoding, decompression, allocation and
    byte-swapping. That is a large constant factor, not an asymptotic change.
  - Exposed everywhere: Python `read`/`read_metadata`, pybind `points_only`/`arrays` kwargs,
    C `mio_read_ex` + `mio_read_opts` + the opaque `mio_read_metadata` handle (`mio_read` is
    unchanged), Fortran `mesh%read(..., points_only=, arrays=)` and `mio_read_metadata`, WASM
    `readMeshSelective`/`readMetadata`, and both CLIs (`info --fast`,
    `convert --points-only|--arrays a,b`).
- **Memory-mapped reading** — `detail::FileSource` maps whole files where that pays (POSIX
  `mmap`, Windows `MapViewOfFile`, always buffered under Emscripten), removing a full-file copy
  and the peak-RSS doubling it causes. `Auto` maps regular files at or above 16 MiB
  (`MESHIOPLUSPLUS_MMAP_THRESHOLD` overrides); anything unmappable falls back silently, so
  mapping is advisory and never fails a read. This is a memory-footprint feature more than a
  throughput one. All five whole-file readers use it — gmsh, vtk, ensight, ugrid's ASCII branch
  and openfoam, the last gaining most since it previously paid for two extra full-file copies.
- **Optional zstd and lz4 codecs** for VTK XML block compression — `MESHIOPLUSPLUS_WITH_ZSTD`
  / `_LZ4` (both **off** by default) → `_core.__has_zstd__` / `__has_lz4__`, plus Conan
  `with_zstd`/`with_lz4` and vcpkg `zstd`/`lz4` features. **zlib remains the default
  everywhere.**
  - `lz4` writes `vtkLZ4DataCompressor` in LZ4's raw block format — a real VTK compressor, so
    such files stay readable by VTK and ParaView. Verified both directions against VTK 9.6.
  - `zstd` writes `vtkZSTDDataCompressor`, which is a **meshio++ extension**: VTK ships no ZSTD
    compressor, so ParaView will report an unknown compressor rather than misread the file.
  - A build without a codec reports an error naming the CMake option to enable, and the
    pure-Python reference supports both via the new `codecs` extra
    (`pip install "meshioplusplus[codecs]"`). A file needing a codec available in *neither* is
    genuinely unreadable — a new failure class, and it fails by name.
  - `--codec zlib|lz4|zstd` on both CLIs' `compress`; **rejected** for formats with no block
    codec rather than silently ignored.

## v7.2.1 (2026-07-20)

- Fix: the eight v7.2.0 data-operation bindings were compiled into the WASM module but never forwarded by `wasm/src/index.mjs`'s ergonomic wrapper, so `dataInfo`/`dataCalc`/etc. were unreachable from `loadMeshioPlusPlus()`'s return value (`m.dataInfo is not a function`); also updates `wasm/index.d.ts`'s ambient TypeScript declarations, which had likewise never been extended. No other bindings affected.

## v7.2.0 (2026-07-19)

- **New data operations** — five dependency-free operations acting on a mesh's
  `point_data` / `cell_data` / `field_data` arrays rather than on its geometry, which none of
  them ever modifies:
  - **`data_manage`** (`data_drop` / `data_keep` / `data_rename`): rewrite which arrays a mesh
    carries and under what names. Values, dtypes and shapes are copied verbatim; an unknown key
    raises listing every available key. Phases apply in the order keep → drop → rename.
  - **`point_data_to_cell_data` / `cell_data_to_point_data`**: move data between locations by
    averaging. Point→cell is the mean over each cell's own nodes; cell→point is the mean over
    the incident cells, optionally weighted by each cell's |measure| (area/volume). Output is
    always `float64`, since a mean is not an integer.
  - **`data_calc`**: derive a new array from an elementwise expression over existing arrays at
    the same location. The evaluator is a hand-written tokenizer plus recursive-descent parser
    supporting `+ - * /`, unary minus, parentheses, numeric literals, array names, and
    `abs`/`sqrt`/`min`/`max`/`norm` — **no external parser library and no evaluation of
    arbitrary code**. Identifiers may contain `:` (so `gmsh:physical` works) and backtick
    quoting handles names with spaces.
  - **`data_condition`**: clamp to `[lo, hi]`, normalize onto a target range (default `[0, 1]`),
    or standardize to zero mean / unit standard deviation — per component or by row magnitude.
    For `cell_data` the statistics are computed jointly across all cell blocks.
  - **`data_info`**: a read-only per-array summary (location, dtype, shape, components, entry
    count, min/max/mean whole-array and per component, NaN/inf counts) — the data-side
    complement to the topological `info` and the geometric `stats`.
- **Documented NaN/inf policy**, shared by all five: non-finite values are always excluded from
  every reduction, and `nan_policy` (`ignore` / `replace` / `fail`) decides only what reaches
  the output. `data_info` never raises — it counts them.
- **New nested CLI group `meshioplusplus data <verb>`** with nine verbs (`info`, `rename`,
  `drop`, `keep`, `to-cell`, `to-point`, `calc`, `clamp`, `normalize`), in both the Python CLI
  and the Python-free native binary. This is the project's first two-level subcommand.
- Exposed across every binding surface: pybind (`meshioplusplus.data_*`), the C API
  (`mio_data_*` plus the opaque `mio_data_info` handle), the Fortran module (type-bound
  `data_*` procedures), and WASM (`dataCalc`, `dataInfo`, …). Fortran additionally exports
  `STRBUF_LEN`, which consumers need in order to declare the `keys` out-argument of `split`
  and `data_info`.
- Documented flat-ABI gap: the combined `data_manage` (keep + drop + rename in one call) is not
  exposed over the C ABI; the three primitives compose to the same effect.
- Documented at `doc/data_manage.md`, `doc/data_average.md`, `doc/data_calc.md`,
  `doc/data_condition.md` and `doc/data_info.md`. Not breaking.

## v7.1.0 (2026-07-19)

- **New CLI verbs for the editing and statistics operations**: `meshioplusplus transform`,
  `clean`, `crop`, `split` and `stats` in both the Python CLI and the native binary, with the
  matching README/`doc/cli.md` documentation.

## v7.0.0 (2026-07-19)

- **Breaking: the default branch moved from `main` to `master`** — CI workflow references and
  the README badges were updated accordingly; consumers pinning the branch in a URL need to
  follow.
- **Breaking: `find_package(meshioplusplus)` consumers must require version 7.0**; the packaged
  CMake config version was bumped with the release.
- **New mesh operations**: `transform` (affine transform of point coordinates, with
  translate/scale/rotate/matrix/units builders and optional vector/tensor rotation), `clean`
  (one-pass weld / drop-degenerate / drop-duplicate / remove-orphans), `crop` (subset by
  bounding box or half-space), `split` (partition by cell type, connected component, or integer
  tag), and `stats` (bounding box, centroid, per-cell-type counts, area, signed/unsigned volume,
  inverted-cell count). All are exposed across Python, the C API, Fortran, WASM and both CLIs,
  and share `detail/subset.hpp` for the prune-and-remap step.
- **arm64 support across the release artifacts**: Linux arm64 wheels, native CLI binaries and
  Conan packages are now built natively on GitHub's hosted arm64 runners (no QEMU).
- **Static runtime linking** (`MESHIOPLUSPLUS_STATIC_RUNTIME`, hoisted to a top-level option) so
  the prebuilt CLI binaries carry no `libstdc++`/`libgcc_s` dependency, including the MSVC
  static CRT on Windows.

## v6.9.0 (2026-07-19)

- **New `diff` operation and CLI verb**: compare two meshes with absolute/relative tolerances,
  reporting per-section differences (points, cells, each data map) with max abs/rel error and
  the worst index. An `unordered` mode matches points by proximity via a bucket-grid hash.
  `meshioplusplus diff A B` **exits non-zero when the meshes differ**, for use in CI and
  Makefiles. Named `point_sets`/`cell_sets` are compared in the Python shim only.
- **New `merge` operation and CLI verb**: combine two or more meshes, either concatenating or
  welding coincident points within a tolerance (spatial-hash bucket grid, no O(N²)), with
  `source_mesh_id` tagging, intersection/fill data policies, and optional duplicate-cell
  dropping.

## v6.8.0 (2026-07-19)

- **New `reorder` operation and CLI verb**: renumber nodes and elements by reverse Cuthill–McKee,
  Morton order, or Hilbert order, as a pure permutation that preserves all geometry and data and
  returns the applied permutations. `compute_bandwidth` reports the connectivity bandwidth.
- **New standalone, Python-free native CLI binary** (`meshioplusplus_cli`, installed as
  `meshioplusplus`), built over the shared format registry and the operations layer. Prebuilt
  statically-linked binaries for Linux x86_64/arm64, macOS universal and Windows x86_64 are
  attached to every `v*` tag's GitHub Release.
- **Doxygen API reference** for the C++/C headers, generated in CI and published alongside the
  VitePress site at `/api/`.

## v6.7.0 (2026-07-19)

- **New `extract_surface` operation**: the dimension-aware generalization of `extract_skin` —
  boundary faces of a 3D volume mesh, or boundary edges of a 2D surface mesh — sharing one
  implementation with the skin extractor. Optional `record_parent_ids` attaches the owning input
  cell index as `surface:parent_cell`.
- **New `sniff_format` operation**: content-based format detection from a file's leading bytes,
  wired in as a **read-only** fallback when the extension yields nothing. It returns a format
  only on a confident signature match, never guessing at ambiguous magics.

## v6.6.0 (2026-07-18)

- **New skin extraction** (`meshioplusplus.extract_skin(mesh, linearize=False)`): derives the boundary surface of a 3D volume mesh (tetra/hexahedron/wedge/pyramid + tetra10/hexahedron20/27/wedge15/pyramid13/14 → triangle/quad and their quadratic variants), implementing the face-hashing algorithm of Kratos Multiphysics' `SkinDetectionProcess` (credit: Kratos, BSD — algorithm reimplemented, no code copied). C++ core (uniform mesh API, all three backends) with a byte-equivalent numpy fallback. Points are compacted; `point_data` follows; `cell_data`/`field_data`/sets are dropped. Documented at `doc/extract_skin.md`.
- **Breaking: STL and PLY write the extracted skin of volume meshes by default.** `stl.write`/`ply.write` (and the C++ `write_stl`/`write_ply`, the flat registry, and therefore the WASM/C-API/Fortran surfaces) gain a `skin=True` parameter: a mesh containing supported volume cells now writes its boundary skin (STL triangulates quads; PLY compacts the vertex table) instead of silently dropping the volume cells and writing an empty/vertex-only file. Pre-existing surface blocks are dropped with a warning in that mode. Pass `skin=False` for the previous behavior.
- **SVG and TikZ render 3D meshes**: genuinely non-flat input no longer raises — the boundary skin (or the surface cells of a 3D shell mesh) is projected through an orthographic camera (`azimuth`/`elevation`/`roll` in degrees, default the classic CAD isometric view) and painted back-to-front (painter's algorithm). Flat 2D output is byte-identical to previous releases; the TikZ C++/Python byte-identity guarantee extends to the 3D path.
- **New logo**: the Stanford Bunny (`example/Bunny.stl`, "Stanford Bunny — Digitized!" by MakerBot, CC-BY, thingiverse thing:88208), decimated and rendered through meshio++'s own TikZ 3D machinery with the blue→teal palette (`logo/gen_logo_tikz.py`).

## v6.5.0 (2026-07-18)

- **New `ensight` format** (`.case`/`.geo`): EnSight Gold geometry, read **and** write, ASCII and C-binary (foreign byte order auto-detected on read; the writer ports Kratos's `EnSightOutput` Gold write logic onto the meshio++ mesh API). Multi-part files concatenate into one point array with the owning part tagged as `cell_data["ensight:part"]`; `nsided`/`nfaced` sections read into polygon/polyhedron blocks (write of ragged blocks raises). Backed by the C++ core with a full-fidelity pure-Python fallback.
- **New `vtp` format** (`.vtp`, VTK XML PolyData): read and write of surface meshes (`vertex`/`line`/`triangle`/`quad`/`polygon`), reusing the VTU base64/zlib stack (`binary`/`compression`/`header_type` parameters mirror `vtu`; lzma is Python-only). The shared VTK-XML `<DataArray>` helpers moved into `detail/vtk_xml.hpp` (VTU output is unchanged).
- **New `triangle` format** (`.node`/`.ele`/`.poly`): Shewchuk's Triangle, the 2D analogue of tetgen — `.node`/`.ele` pairs (`triangle`/`triangle6`) plus the `.poly` PSLG (segments as `line` cells; holes/regions skipped). `.node`/`.ele` still default to tetgen; the reader dispatcher falls through to `triangle` for 2D pairs, and writes need `file_format="triangle"` (only `.poly` defaults to triangle).
- All three formats are registered in the shared dispatch registry and therefore reachable from the WASM, C API, and Fortran flat bindings — WASM now ships 36 formats (35 readable / 36 writable). Not breaking.

## v6.4.0 (2026-07-18)

- **New `tikz` format** (`.tikz`): a write-only, 2D TikZ/PGF (LaTeX) writer, the LaTeX counterpart to `svg`. Emits a directly `pdflatex`-compilable `standalone` document by default (`standalone=False` for a bare `tikzpicture` snippet). **Both `svg` and `tikz` are now backed by the C++ core** (`write_svg`/`write_tikz`) with the pure-Python reference kept as fallback, and are registered in the shared dispatch registry, so they are additionally reachable (write-only, fixed default styling) from the WASM, C API, and Fortran flat bindings — WASM now ships 33 writable formats. No change to existing APIs; documented at `doc/formats/tikz.md`/`svg.md`. Not breaking.

## v6.3.2 (2026-07-17)

- **Coverage extended**: the `coverage` job now instruments the C API (`bindings_c/c_api.cpp` + its gtest suite, previously dark) and drops structurally-unreachable code (the non-MESHIO mesh-backend headers, covered by the separate `cpp-tests` matrix, and the generated `single_include/`) from the denominator. New tests lift the darkest paths — a `ply` C++ suite (the one format that had none), UGRID binary/endian flavours, malformed-input `ReadError` cases across ply/ugrid/su2/tetgen/vtk/xdmf/med, and Python public-API error paths + CLI edge cases. Tests/CI only; no API change. Not breaking.

## v6.3.1 (2026-07-17)

- **Coverage CI properly wired up**: the combined Python + C++ `coverage` job now runs (it was gated behind a red `lint` check and had never executed), uploads to Codecov under the `python`/`cpp` flags, and gates PRs (project/patch statuses flipped from informational to blocking). `pyproject.toml` gains `[tool.coverage.run]` (`relative_files`) so `coverage.xml` paths match Codecov's flag filters. Tooling/CI only; no API change. Not breaking.

## v6.3.0 (2026-07-17)

- **Single-header, header-only C++ distribution**: `single_include/meshioplusplus/meshioplusplus.hpp`, generated by `tools/amalgamate.sh` and kept up to date by CI. Declarations are always visible; `#define MESHIOPLUSPLUS_IMPLEMENTATION` before including it in one translation unit pulls in the implementation, pugixml bundled and no external dependencies by default (HDF5/netCDF/zlib/Eigen stay opt-in behind their existing `MESHIOPLUSPLUS_HAS_*` macros). No API change; documented at `doc/single_header.md`. Not breaking.

## v6.2.0 (2026-07-17)

- **New C API and Fortran interface** for HPC consumers: an installable
  `libmeshioplusplus` shared library with a stable pure-C99 header
  (`mio_*` functions; pkg-config + `find_package(meshioplusplus)` support)
  and a modern OO Fortran 2008 module (`type(mio_mesh)` with type-bound
  procedures) on top of it via ISO_C_BINDING. Off by default
  (`MESHIOPLUSPLUS_BUILD_C_API` / `MESHIOPLUSPLUS_BUILD_FORTRAN`, or
  `build/configure.sh --c-api` / `--fortran`); Python/WASM artifacts are
  unaffected. The WASM binding's format-dispatch tables moved into a shared
  core registry (`meshioplusplus/registry.hpp`) used by both flat bindings —
  no JS API change. Not breaking; documented at `doc/c_api.md` and
  `doc/fortran.md`.

## v6.0.0 (2026-07-14)

- **Default C++ parallel backend is now `AUTO`** (prefers OpenMP, then
  STL+TBB, then sequential) instead of `STL` — the old default silently ran
  sequentially on libstdc++ without TBB. Published wheels are now parallel.
  `meshioplusplus._core.__parallel_backend__` reports the active backend. The
  binary-format read/write paths were also optimised (bulk-buffered I/O); output
  is unchanged (byte-identical). **Source builds** should now run
  `git submodule update --init` to fetch the vendored **Eigen** (used for the
  MED transpose); it is optional — builds without it fall back to a plain loop.
- **Project renamed to meshio++** (machine identifier `meshioplusplus`, used
  wherever a literal `+` isn't valid: the Python package/import name, PyPI
  distribution, CLI entry point, C++ namespace, CMake project/targets, and
  build macros). This is a clean break with no compatibility shim — `import
  meshio` / `pip install meshio` no longer refer to this project; use `import
  meshioplusplus` / `pip install meshioplusplus` going forward. The public API
  surface, file formats, and behavior are otherwise unchanged from v5.x.
- Added two new formats, `ansysInp` (Ansys/APDL coded database, `.cdb`/`.inp`) and
  `openfoam` (OpenFOAM polyMesh, read-only), and significantly extended MED/Salome
  support: multi-mesh files (`meshioplusplus.med.read_med_multi`/`write_med_multi`), ragged
  polygon/Voronoi cell blocks, MED 4.1 bitmask metadata, node-orientation fixes,
  quadratic `triangle7`/`quad9`/polygon type support, mesh-level metadata
  (`mesh_name`/`description`/`unit_time`/`unit_coords`), and preserving Gmsh physical
  groups as MED families on write. `meshioplusplus.med.read`/`write` now always use the Python
  implementation (the C++ `meshioplusplus._core.med_read`/`med_write` bindings remain directly
  callable for the narrower/faster behavior). This work originates from
  [Simvia's `meshlane` fork](https://github.com/simvia-tech/meshlane) of meshio,
  contributed by Mariam Kesba, Fatima-Zahra Noussi, and Lucas Sovre, and has been
  brought back into this repository.
- The C++ core is now C++20 (previously C++17) with `std::format`-based logging
  (`MESHIOPLUSPLUS_LOG_LEVEL` env var) and a compile-time-selectable parallel
  backend for hot loops (`-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ|STL|OPENMP|TBB`,
  default STL).

## v5.1.0 (Dec 11, 2021)

- CellBlocks are no longer tuples, but classes. You can no longer iterate over them like
  ```python
  for cell_type, cell_data in cells:
      pass
  ```
  Instead, use
  ```python
  for cell_block in cells:
      cell_block.type
      cell_block.data
  ```

## v5.0.0 (Aug 06, 2021)

- meshio now only provides one command-line tool, `meshio`, with subcommands like
  `info`, `convert`, etc. This replaces the former `meshio-info`, `meshio-convert` etc.

## v4.4.0 (Apr 29, 2021)

- Polygons are now stored as `"polygon"` cell blocks, not `"polygonN"` (where `N` is the
  number of nodes per polygon). One can simply retrieve the number of points via
  `cellblock.data.shape[1]`.

## v4.0.0 (Feb 18, 2020)

- `mesh.cells` used to be a dictionary of the form

  ```python
  {
    "triangle": [[0, 1, 2], [0, 2, 3]],
    "quad": [[0, 7, 1, 10], ...]
  }
  ```

  From 4.0.0 on, `mesh.cells` is a list of tuples,

  ```python
  [
    ("triangle", [[0, 1, 2], [0, 2, 3]]),
    ("quad", [[0, 7, 1, 10], ...])
  ]
  ```

  This has the advantage that multiple blocks of the same cell type can be accounted
  for. Also, cell ordering can be preserved.

  You can now use the method `mesh.get_cells_type("triangle")` to get all cells of
  `"triangle"` type, or use `mesh.cells_dict` to build the old dictionary structure.

- `mesh.cell_data` used to be a dictionary of the form

  ```python
  {
    "triangle": {"a": [0.5, 1.3], "b": [2.17, 41.3]},
    "quad": {"a": [1.1, -0.3, ...], "b": [3.14, 1.61, ...]},
  }
  ```

  From 4.0.0 on, `mesh.cell_data` is a dictionary of lists,

  ```python
  {
    "a": [[0.5, 1.3], [1.1, -0.3, ...]],
    "b": [[2.17, 41.3], [3.14, 1.61, ...]],
  }
  ```

  Each data list, e.g., `mesh.cell_data["a"]`, can be `zip`ped with `mesh.cells`.

  An old-style `cell_data` dictionary can be retrieved via `mesh.cell_data_dict`.
