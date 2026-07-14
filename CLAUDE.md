# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

meshio++ (Python package/import name `meshioplusplus`) is a Python library for reading and writing many mesh file formats used in scientific computing and FEM (Finite Element Method). It provides a unified `Mesh` data structure that all format readers/writers convert to and from.

This branch (`port-to-c++`) adds a **C++ core** (`meshioplusplus._core`, pybind11) built with **scikit-build-core + CMake**, replacing the old setuptools build. Most formats read/write through the C++ core with a pure-Python fallback.

## Keeping docs in sync (required)

**Every time a feature or format is added or changed, update all three of `CLAUDE.md`, `README.md`, and the `doc/` site in the same change.** At minimum: the format table and format-specific notes in [`doc/formats.md`](doc/formats.md) (plus any other relevant `doc/` page), the supported-formats list in `README.md`, and the architecture/"Adding a new format" guidance here. A feature is not "done" until CLAUDE.md, README.md, and the docs all reflect it.

## Commands

**Build/install for development** (native paths on). PEP 668 systems: use the in-repo uv venv (`uv venv --python 3.12 .venv`), then:
```bash
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON" \
  uv pip install --python .venv --no-build-isolation -e .
```
Plain environments: `pip install -e ".[all]"`. The optional deps (HDF5, netCDF, zlib) are auto-detected; when absent the C++ paths compile out and Python fallbacks (h5py/netCDF4/stdlib zlib) take over.

**Standalone C++ build**: `build/configure.sh` (Linux/macOS) and `build/configure.bat` (Windows) configure a CMake tree next to themselves (`build/cpp-<type>`); flags: `--backend SEQ|STL|OPENMP|TBB`, `--tests`, `--build`, `--with/-out-hdf5/netcdf/zlib`, `--tbb-dir <path>`. They print the matching `CMAKE_ARGS … pip install` line.

**Run tests:** `pytest tests/` (or `.venv/bin/python -m pytest tests/`).
Single file/test: `pytest tests/test_gmsh.py::test_gmsh22`.

**Lint / format:** `pre-commit run -a` (isort, black, flake8). Or `black --check . && flake8 . && isort --check .`.

**CLI usage:**
```bash
meshioplusplus convert input.msh output.vtk
meshioplusplus info mesh.msh
meshioplusplus ascii mesh.msh         # convert to ASCII
meshioplusplus binary mesh.msh        # convert to binary
```

**Docs (VitePress):** `cd doc && npm install && npm run docs:build` (dev: `npm run docs:dev`).

**Logo:** built with TikZ. `logo/build.sh` runs `logo/gen_logo_tikz.py`
(numpy+matplotlib triangulation of an "FE blob" → `logo/_mesh_icon.tex`),
compiles `logo/logo.tex`/`logo-icon.tex` with `pdflatex`, and converts to SVG
via `dvisvgm` (PNG via PyMuPDF). Committed assets: `logo/logo-with-text.svg`
(README banner + `doc/public/logo.svg`), `logo/logo-icon.svg` (favicon/nav).
The old pygmsh generator `logo/logo.py` is superseded (see `logo/README.md`).

**Example notebooks** (`example/*.ipynb`): read the bundled `example/example.msh`
with meshio++ and render/convert it (PyVista off-screen via a VTU round-trip,
matplotlib fallback). Committed **with outputs**; re-execute with
`PYVISTA_OFF_SCREEN=true jupyter nbconvert --to notebook --execute --inplace example/*.ipynb`.

**Benchmarks** (`benchmark/`): `bench.py` times read/write on formats both
libraries support — meshio++ vs the legacy pure-Python `meshio` (imported from
`/home/vicente/src/meshio_legacy/src` via `sys.path`, no build). `inputs.py`
provides the real `example/example.msh` bracket (the headline input) + a
synthetic numpy tet grid (also used for a size-scaling sweep). `01_benchmark.ipynb`
runs it, writes `results.csv`, and regenerates the plots in
`doc/public/benchmarks/` (`benchmark_times`/`_speedup` = the bracket,
`benchmark_scaling` = speedup vs mesh size) shown on `doc/benchmarks.md`. The needed extras
(`pyvista matplotlib jupyter nbconvert ipykernel`) are installed into `.venv`
with `uv pip install`.

## Architecture

**Core data model** (`src/meshioplusplus/_mesh.py`, pure Python, unchanged):
- `Mesh`: holds `points`, `cells` (list of `CellBlock`), `point_data`, `cell_data`, `field_data`, `point_sets`, `cell_sets`, optionally `gmsh_periodic`/`info`
- `CellBlock`: a cell type string (e.g. `"triangle"`, `"tetra10"`) + numpy node-index array
- `_common.py`: `num_nodes_per_cell`; `topological_dimension` in `_mesh.py`

**C++ core** (`cpp/`, `bindings/`):
- `cpp/include/meshioplusplus/` headers (`ndarray.hpp`, `mesh.hpp`, `types.hpp`, `detail/*`, `formats/*`), `cpp/src/formats/*.cpp`, compiled with `bindings/_core.cpp` into `meshioplusplus._core`.
- **Zero-copy at the I/O boundary**: readers return capsule-backed writeable numpy; writers view numpy memory (`bindings/np_conversions.hpp`). The conversion layer carries points/cells/point_data/cell_data/field_data, but **not** `mesh.info`, `cell_sets`, or `point_sets` — formats that need those either defer to Python or carry them out-of-band via a **side-channel struct** the binding `setattr`s onto the Python `Mesh` (e.g. `MedInfo`, `AnsysInfo` for `point_sets`/`cell_sets`, `OpenFoamInfo` for `cell_tags`).
- **Ragged cell blocks**: `meshioplusplus::CellBlock` also holds optional ragged data — `polygon_rows` (1-level: jagged polygons) and `polyhedron_rows` (2-level: list of faces per cell) — for types that can't fit a rectangular `NDArray`. Zero-copy applies only to rectangular blocks; ragged blocks are **copied** across the boundary. `py_to_mesh`'s `allow_ragged` flag is **off by default** (so every rectangular-only writer keeps safely rejecting ragged meshes → Python fallback) and only the ragged-aware bindings (MED write) opt in.
- **Optional deps** (`CMakeLists.txt`): `MESHIOPLUSPLUS_WITH_HDF5`/`_NETCDF`/`_ZLIB`/`_EIGEN` options → `MESHIOPLUSPLUS_HAS_*` compile definitions and `_core.__has_hdf5__`/`__has_netcdf__` flags. `#ifdef`-guarded sources become empty TUs (or fall back to a plain loop) when off. Shared HDF5 helpers in `cpp/include/meshioplusplus/detail/hdf5_util.hpp`. **Eigen** is vendored as a git submodule at `cpp/third_party/eigen` (header-only; run `git submodule update --init` for source builds) and used for the MED Fortran↔C transpose (`med.cpp`, guarded by `MESHIOPLUSPLUS_HAS_EIGEN`); when the submodule is absent the code uses a hand-written transpose, so sdists without submodules still build.
- **C++ standard: C++20** (pinned in `CMakeLists.txt` twice: `CMAKE_CXX_STANDARD` + `target_compile_features`). macOS needs deployment target ≥ 13.3 for `std::format` (set in `ci.yml`/`wheels.yml`).
- **Logging** (`cpp/include/meshioplusplus/log.hpp`): `meshioplusplus::log::debug/info/warn/error("fmt {}", …)` — `std::format` (compile-time-checked) + `std::source_location`, thread-safe via `std::osyncstream`, filtered at runtime by the `MESHIOPLUSPLUS_LOG_LEVEL` env var (`debug|info|warn|error|off`, default `warn`). No printf/cerr elsewhere; errors remain exceptions.
- **Parallelism** (`cpp/include/meshioplusplus/parallel.hpp`): `meshioplusplus::parallel_for(n, f, grain=2048, max_threads=0)` with a compile-time backend selected by `MESHIOPLUSPLUS_PARALLEL_BACKEND` = `AUTO|SEQ|STL|OPENMP|TBB`. **Default is `AUTO`**, which prefers **OpenMP** (portable: libgomp on manylinux, MSVC built-in, libomp on macOS; needs no TBB), else the STL(+TBB) probe, else SEQ. `_core.__parallel_backend__` reports the active backend from Python (STL-without-TBB is effectively sequential — that's why AUTO avoids it). Adding a backend (Kokkos, …) = one CMake branch + one `#elif` in `parallel.hpp`. Iterations must be independent; the first exception is captured and rethrown after the join; `n <= grain` runs sequentially. **Two flavors:** `parallel_for` (all cores — for compute-bound loops: zlib/base64 in `detail/vtu_binary.hpp`, ASCII formatting) and **`parallel_for_bw`** (thread-capped, `parallel_bandwidth_threads=4` — for memory-bandwidth-bound loops: byte-swap/transpose/gather, which *regress* past ~4 threads because they saturate memory bandwidth). Use `parallel_for_bw` for byte I/O and index gathers; `parallel_for` only when the per-element work is real compute. Hoist per-element dtype switches with `detail::dispatch_dtype(dt, []<class T>(){…})` (`value_io.hpp`) before parallelizing. **Binary I/O rule:** never write per-element to a `std::ostream` (byte-at-a-time `os.put`/tiny `os.write` was the pre-optimization VTK/Gmsh bottleneck) — build one pre-sized buffer (fused gather+byteswap where possible) and emit it with a single `os.write`; slurp input files with a bulk seek+read, not `istreambuf_iterator`.

**Format registration** (`src/meshioplusplus/_helpers.py`): each format's `__init__.py` calls `register_format(name, extensions, reader, {name: writer})`. `read`/`write` auto-detect from extension. Modules are imported (and thus registered) in `src/meshioplusplus/__init__.py`.

**Format module layout** (the shim pattern, e.g. `src/meshioplusplus/su2/`):
- `__init__.py` — imports `_core`, defines `read`/`write` that **try the C++ function and fall back** to the Python reference on any exception, then `register_format(...)`.
- `_<format>.py` — the pure-Python reference reader/writer.
- Multi-version formats (Gmsh, VTK) dispatch to version submodules.

**Cell type naming**: meshio++ uses its own names (`triangle`, `tetra10`, `hexahedron20`, …). Each format maps between its native names and meshio++'s (e.g. `gmsh/common.py`).

**Tests** (`tests/`): `helpers.py` has `Mesh` fixtures + `write_read()` (the round-trip pattern). Per-format `test_<format>.py` parametrize over meshes. `tests/meshes/` and `tests/input/` hold read-only reference files (Git-LFS). HDF5/netCDF suites use `pytest.importorskip`.

**CLI** (`src/meshioplusplus/_cli/`): `convert`, `info`, `ascii`, `binary`, `compress`, `decompress`, registered in `_main.py`.

**CI** (`.github/workflows/`): `ci.yml` (3-OS test matrix; Linux/macOS build native paths, Windows builds them off and uses Python fallbacks), `wheels.yml` (cibuildwheel with native paths off + PyPI trusted publishing on `v*` tags), `docs.yml` (VitePress → GitHub Pages).

## Adding a new format

1. Create `src/meshioplusplus/<format>/_<format>.py` (Python reference `read`/`write`) and `__init__.py` (shim + `register_format`). `register_format` may also live directly in `_<format>.py` for a Python-only format with no shim to wrap (e.g. `mdpa`, `neuroglancer`, `svg`, `ansysInp`, `openfoam`).
2. Add the C++ implementation: `cpp/src/formats/<format>.cpp` + `cpp/include/meshioplusplus/formats/<format>.hpp`, bind it in `bindings/_core.cpp`, and switch `__init__.py` to the try-C++/Python-fallback shim. New `.cpp` files are auto-globbed by CMake. **This step is optional** — a format can stay Python-only indefinitely (precedent: `mdpa`, `neuroglancer`, `svg`, and the Simvia-contributed `ansysInp`/`openfoam`); don't force a C++ port just to satisfy this checklist.
3. Import the module in `src/meshioplusplus/__init__.py` (both the import tuple and `__all__`).
4. Add `tests/test_<format>.py` using `helpers.write_read()`. Verify a direct C++ roundtrip and cross-compat with the Python reference (if a C++ path exists).
5. **Update `README.md`, `doc/formats.md`, and this file** (see "Keeping docs in sync").

Do not copy test data from GPL sources (e.g. FEconv examples) into this MIT repo — generate reference files via round-trip instead. Code/fixtures from other MIT-licensed forks in the same lineage (e.g. the [meshlane](https://github.com/simvia-tech/meshlane) fork) may be copied directly; credit the source in `CITATION.cff`/`CHANGELOG.md`.

**Note on `med`:** the C++ core (built with HDF5) handles the mesh-representation part of MED exactly — points, point/cell tags, families with `GRO` group names, mesh-level metadata (`mesh_name`/`description`/`unit_time`/`unit_coords`/`point_tag_groups`/`cell_tag_groups`, via `MedInfo`), node-orientation permutations, and `POG`/`POG2` ragged polygons — and `meshioplusplus.med.read`/`write` use it by default. It **raises** (→ Python fallback) for the constructs it does not replicate byte-for-byte: `CHA` fields (MED-4.1 bitmask / units / step metadata), the `gmsh:physical`→family bridging, non-default profiles/ELGA, and multi-mesh files (`read_med_multi`/`write_med_multi` are always Python). The C++ MED reader iterates `MAI` cell blocks in HDF5 **creation order** (`group_links_crt`, matching h5py `track_order`) since block order aligns cell_data/cell_sets; the shim reconstructs `point_sets`/`cell_sets` from families via the shared Python helpers. See [`doc/formats/med.md`](doc/formats/med.md#quirks-limitations).
