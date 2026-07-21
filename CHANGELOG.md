<!--pytest-codeblocks:skipfile-->
# Changelog

This document records every released version of meshio++ — new formats, new operations,
notable enhancements, and breaking changes. Breaking changes are called out explicitly as
**Breaking:**; everything else is additive unless stated otherwise.

**Keep this file current: add an entry in the same change as every version bump.** See the
"Version bumps" section of `CLAUDE.md`.

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
