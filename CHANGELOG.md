<!--pytest-codeblocks:skipfile-->
# Changelog

This document only describes _breaking_ changes in meshio++. If you are interested in bug
fixes, enhancements etc., best follow [the meshio++ project on
GitHub](https://github.com/loumalouomega/meshioplusplus).

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
