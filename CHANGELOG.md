<!--pytest-codeblocks:skipfile-->
# Changelog

This document records every released version of meshio++ — new formats, new operations,
notable enhancements, and breaking changes. Breaking changes are called out explicitly as
**Breaking:**; everything else is additive unless stated otherwise.

**Keep this file current: add an entry in the same change as every version bump.** See the
"Version bumps" section of `CLAUDE.md`.

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
