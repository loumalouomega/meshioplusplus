# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

meshio is a Python library for reading and writing many mesh file formats used in scientific computing and FEM (Finite Element Method). It provides a unified `Mesh` data structure that all format readers/writers convert to and from.

This branch (`port-to-c++`) adds a **C++ core** (`meshio._core`, pybind11) built with **scikit-build-core + CMake**, replacing the old setuptools build. Most formats read/write through the C++ core with a pure-Python fallback; the public Python API is unchanged (a drop-in replacement).

## Keeping docs in sync (required)

**Every time a feature or format is added or changed, update all three of `CLAUDE.md`, `README.md`, and the `doc/` site in the same change.** At minimum: the format table and format-specific notes in [`doc/formats.md`](doc/formats.md) (plus any other relevant `doc/` page), the supported-formats list in `README.md`, and the architecture/"Adding a new format" guidance here. A feature is not "done" until CLAUDE.md, README.md, and the docs all reflect it.

## Commands

**Build/install for development** (native paths on). PEP 668 systems: use the in-repo uv venv (`uv venv --python 3.12 .venv`), then:
```bash
CMAKE_ARGS="-DMESHIO_WITH_HDF5=ON -DMESHIO_WITH_NETCDF=ON -DMESHIO_WITH_ZLIB=ON" \
  uv pip install --python .venv --no-build-isolation -e .
```
Plain environments: `pip install -e ".[all]"`. The optional deps (HDF5, netCDF, zlib) are auto-detected; when absent the C++ paths compile out and Python fallbacks (h5py/netCDF4/stdlib zlib) take over.

**Run tests:** `pytest tests/` (or `.venv/bin/python -m pytest tests/`).
Single file/test: `pytest tests/test_gmsh.py::test_gmsh22`.

**Lint / format:** `pre-commit run -a` (isort, black, flake8). Or `black --check . && flake8 . && isort --check .`.

**CLI usage:**
```bash
meshio convert input.msh output.vtk
meshio info mesh.msh
meshio ascii mesh.msh         # convert to ASCII
meshio binary mesh.msh        # convert to binary
```

**Docs (VitePress):** `cd doc && npm install && npm run docs:build` (dev: `npm run docs:dev`).

## Architecture

**Core data model** (`src/meshio/_mesh.py`, pure Python, unchanged):
- `Mesh`: holds `points`, `cells` (list of `CellBlock`), `point_data`, `cell_data`, `field_data`, `point_sets`, `cell_sets`, optionally `gmsh_periodic`/`info`
- `CellBlock`: a cell type string (e.g. `"triangle"`, `"tetra10"`) + numpy node-index array
- `_common.py`: `num_nodes_per_cell`; `topological_dimension` in `_mesh.py`

**C++ core** (`cpp/`, `bindings/`):
- `cpp/include/meshio/` headers (`ndarray.hpp`, `mesh.hpp`, `types.hpp`, `detail/*`, `formats/*`), `cpp/src/formats/*.cpp`, compiled with `bindings/_core.cpp` into `meshio._core`.
- **Zero-copy at the I/O boundary**: readers return capsule-backed writeable numpy; writers view numpy memory (`bindings/np_conversions.hpp`). The conversion layer carries points/cells/point_data/cell_data/field_data, but **not** `mesh.info`, `cell_sets`, or `point_sets` — formats that need those defer to Python.
- **Optional deps** (`CMakeLists.txt`): `MESHIO_WITH_HDF5`/`_NETCDF`/`_ZLIB` options → `MESHIO_HAS_*` compile definitions and `_core.__has_hdf5__`/`__has_netcdf__` flags. `#ifdef`-guarded sources become empty TUs when off. Shared HDF5 helpers in `cpp/include/meshio/detail/hdf5_util.hpp`.

**Format registration** (`src/meshio/_helpers.py`): each format's `__init__.py` calls `register_format(name, extensions, reader, {name: writer})`. `read`/`write` auto-detect from extension. Modules are imported (and thus registered) in `src/meshio/__init__.py`.

**Format module layout** (the shim pattern, e.g. `src/meshio/su2/`):
- `__init__.py` — imports `_core`, defines `read`/`write` that **try the C++ function and fall back** to the Python reference on any exception, then `register_format(...)`.
- `_<format>.py` — the pure-Python reference reader/writer.
- Multi-version formats (Gmsh, VTK) dispatch to version submodules.

**Cell type naming**: meshio uses its own names (`triangle`, `tetra10`, `hexahedron20`, …). Each format maps between its native names and meshio's (e.g. `gmsh/common.py`).

**Tests** (`tests/`): `helpers.py` has `Mesh` fixtures + `write_read()` (the round-trip pattern). Per-format `test_<format>.py` parametrize over meshes. `tests/meshes/` and `tests/input/` hold read-only reference files (Git-LFS). HDF5/netCDF suites use `pytest.importorskip`.

**CLI** (`src/meshio/_cli/`): `convert`, `info`, `ascii`, `binary`, `compress`, `decompress`, registered in `_main.py`.

**CI** (`.github/workflows/`): `ci.yml` (3-OS test matrix; Linux/macOS build native paths, Windows builds them off and uses Python fallbacks), `wheels.yml` (cibuildwheel with native paths off + PyPI trusted publishing on `v*` tags), `docs.yml` (VitePress → GitHub Pages).

## Adding a new format

1. Create `src/meshio/<format>/_<format>.py` (Python reference `read`/`write`) and `__init__.py` (shim + `register_format`).
2. Add the C++ implementation: `cpp/src/formats/<format>.cpp` + `cpp/include/meshio/formats/<format>.hpp`, bind it in `bindings/_core.cpp`, and switch `__init__.py` to the try-C++/Python-fallback shim. New `.cpp` files are auto-globbed by CMake.
3. Import the module in `src/meshio/__init__.py` (both the import tuple and `__all__`).
4. Add `tests/test_<format>.py` using `helpers.write_read()`. Verify a direct C++ roundtrip and cross-compat with the Python reference.
5. **Update `README.md`, `doc/formats.md`, and this file** (see "Keeping docs in sync").

Do not copy test data from GPL sources (e.g. FEconv examples) into this MIT repo — generate reference files via round-trip instead.
