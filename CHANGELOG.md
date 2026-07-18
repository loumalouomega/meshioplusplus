<!--pytest-codeblocks:skipfile-->
# Changelog

This document only describes _breaking_ changes in meshio++. If you are interested in bug
fixes, enhancements etc., best follow [the meshio++ project on
GitHub](https://github.com/loumalouomega/meshioplusplus).

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
