# Repository guidance

meshio++ (`meshioplusplus`) reads, writes, and operates on scientific meshes. Python exposes a unified `Mesh`; the C++20 core (`meshioplusplus._core`, pybind11) builds with scikit-build-core and CMake. Many formats have Python reference implementations and C++ acceleration.

## Keeping this guide lean

Keep this file focused on working rules, commands, and architectural constraints. Put feature behavior and format caveats in the relevant `doc/` page, release history in `CHANGELOG.md`, and local setup details outside shared guidance. Link to those sources instead of repeating their contents. Update this guide when a change affects how contributors should work; do not append an entry for every feature or release.

## Change checklist

- Update `README.md` and relevant `doc/` pages with user-facing changes. Format changes belong in [doc/formats.md](doc/formats.md) and their `doc/formats/<format>.md` page. Registry or binding changes also require the affected [WASM](doc/wasm.md), [C](doc/c_api.md), [Fortran](doc/fortran.md), [Julia](doc/julia.md), and [R](doc/r.md) docs.
- Keep the MCP server in sync with the Python public API: update both `src/python/meshioplusplus/mcp/_tools.py` and `_server.py`, `tests/python/test_mcp.py`, and [doc/mcp.md](doc/mcp.md) in the same change. This includes signatures, defaults, report fields, and format capabilities. The `test_every_operation_has_a_tool` parity guard checks `__all__` against tool `wraps`; use `_NOT_TOOLS` exemptions only for APIs that cannot be expressed through file paths, such as in-memory interop and interactive viewing.
- When implementing a gap in [doc/roadmap.md](doc/roadmap.md), remove the completed item or narrow it to the remaining work.
- Demonstrate new user-facing capabilities in `example/python/*.ipynb` with graphical output: PyVista off-screen rendering embedded through matplotlib, with a matplotlib fallback when GL is unavailable. Follow `03_mesh_operations.ipynb`; re-execute and commit outputs. Mirror notable additions in the C++ notebooks where practical.
- Keep prose in `doc/**/*.md` on one line per paragraph or list item; do not hard-wrap it. Code fences are exempt. Add new pages to the appropriate sidebar group in `doc/.vitepress/config.mts`.
- Run relevant tests and lint checks before opening or updating a PR. Do not wait for CI to catch formatting errors.

## Commands

Run from the repository root unless specified otherwise.

| Task | Command / notes |
| --- | --- |
| Development environment | `uv venv --python 3.12 .venv` if needed; install build requirements from `pyproject.toml` before using `--no-build-isolation` |
| Editable native build | `CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON" uv pip install --python .venv --no-build-isolation -e .` |
| Standard editable install | `pip install -e ".[all]"` in a virtual environment |
| Python tests | `.venv/bin/python -m pytest tests/python/`; scope to affected files while developing |
| Python lint | `isort --check <files>`, `black --check <files>`, `flake8 <files>`; full hooks: `pre-commit run -a` |
| Standalone C++ tests | `build/configure.sh --tests --build`, then `ctest --test-dir <configured-tree> --output-on-failure` |
| Native CLI | `build/configure.sh --cli --build` |
| C / Fortran libraries | `build/configure.sh --c-api` or `--fortran`, then build and install the configured tree |
| Installable C++ API | `build/configure.sh --install-cpp --cpp-backends MESHIO,NATIVE`; see [doc/cpp_api.md](doc/cpp_api.md) |
| WASM | `build/configure-wasm.sh --build`; requires Emscripten; see [doc/wasm.md](doc/wasm.md) |
| Single header | `tools/amalgamate.sh` to regenerate; `--check` to check freshness; `--smoke` to regenerate and compile |
| C++ includes | `tools/include-cleanup.sh --check` with HDF5/netCDF installed; unused includes fail, missing-include suggestions are advisory |
| Docs (VitePress) | In `doc/`: `npm install`, then `npm run docs:build` (or `npm run docs:dev`) |
| Doxygen reference | In `doc/`: `doxygen Doxyfile`; after the VitePress build, copy `doxygen-build/html` to `.vitepress/dist/api` |
| Python notebook outputs | `PYVISTA_OFF_SCREEN=true jupyter nbconvert --to notebook --execute --inplace example/python/*.ipynb` |

`build/configure.bat` is the Windows counterpart. Configure scripts support `--backend` and `--mesh-backend`; non-MESHIO mesh backends disable Python. Quote multi-entry `--cpp-backends` lists on Windows. Optional native libraries are feature-gated; check the configured capabilities rather than assuming they are available.

## Architecture and implementation rules

### Source map

| Area | Location / reference |
| --- | --- |
| Python mesh model and dispatch | `src/python/meshioplusplus/_mesh.py`, `_helpers.py`, `__init__.py` |
| Native interfaces and implementations | `src/cpp/include/meshioplusplus/`, `src/cpp/src/` |
| Python conversion and bindings | `bindings/python/np_conversions.hpp`, `_core.cpp` |
| Shared native format registry | `src/cpp/src/registry.cpp`; consumed by the native CLI, C API, and WASM |
| Language bindings | `bindings/{c,fortran,julia,r,wasm}/` |
| CLIs | `src/python/meshioplusplus/_cli/`, `src/cpp/cli/` |
| Browser viewer / npm library | `src/viewer/` / `src/wasm/`; separate packages and dependencies |
| Tests | `tests/python/`, `tests/cpp/`; per-format round trips use `tests/python/helpers.py` |
| CI | `.github/workflows/`; consult the affected workflow for its dependency/backend matrix |

See [architecture](doc/architecture.md), [mesh data model](doc/mesh_data_model.md), and [mesh backends](doc/cpp_backends.md) for details.

### C++ core

- Use the uniform mesh API (`mesh_api.hpp`) in formats, operations, shared helpers, and flat bindings. Use `AssignPoints`, `AddCellBlock`, data accessors, and `CellView` methods rather than backend members. Code must work with MESHIO, NATIVE, and KRATOS; KRATOS meshes are not copy-constructible, so clone through the uniform API. Hoist mesh/data lookups out of hot loops.
- `bindings/python/np_conversions.hpp` is the backend-specific exception: Python requires MESHIO. Rectangular arrays use capsule-backed, writable NumPy storage on reads and views on writes. Do not mutate caller-owned buffers. Copy regions before canonicalizing; carry format-specific metadata through explicit side channels when the common conversion cannot represent it.
- Support ragged polygons/polyhedra through `Row`/`Face` accessors and the shared [polyhedron kernel](doc/polyhedra.md). Python ragged conversions copy and require `allow_ragged` opt-in. A polygon can have rectangular storage: use its cell type, not just `IsRagged()`, to decide support.
- Reuse `detail/region_remap.*` for [regions](doc/regions.md). Select the correct map shape (`Direct`, `FirstChild`, or `Global`); retain empty named groups and warn when dropping regions. Region cell indices are global block-major; compatibility `cell_sets` entries are local to each block. Property sets are keyed by ID, not cell index; preserve them according to the operation's documented contract.
- Preserve buffer ownership when cloning arrays, including rank-0 scalars: an empty shape with a buffer is one element; a default absent array is zero elements. Use owned copies where results outlive inputs.
- Data-only operations must leave points, connectivity, and block order bit-identical. Reuse `detail/data_ops.*`, topology tables, subdivision helpers, and shared geometry kernels instead of introducing parallel implementations.
- Use `registry_write_ex` / `WriteOptions` for parameterized native writes. The basic registry writer has fixed defaults. Flat bindings have no Python fallback; test their native paths directly.
- Use `parallel_for` for independent hot loops; preserve deterministic ordering and reduction behavior. Default parallel backend selection is `AUTO`; do not assume it means a particular runtime.
- Prefer unordered containers for order-independent lookup. Preserve sorted data-name iteration and deterministic output. Direct mutations of MESHIO data maps must call `InvalidateNameCaches()`.
- Follow `.clang-format`: K&R braces, four spaces. Types and methods use `PascalCase`, members `mPascalCase`, reference/pointer parameters `rName`/`pName` (excluding forwarding/move references), and free functions/local variables/files `snake_case`. Keep Python/JS API string names stable. Do not format third-party code.
- Use `meshioplusplus::log::{debug,info,warn,error}` from `log.hpp`; errors remain exceptions. Avoid ad hoc `printf`/`cerr` logging.
- Keep number I/O locale-independent. Floating-point parsing and formatting go through `detail/fast_number.hpp` (`parse_double`, `snprintf_c`), never bare `strtod`/`stod`/`snprintf("%f")`; every stream is built with `detail/classic_stream.hpp` (`auto in = detail::make_classic_ifstream(rPath)`), never `std::ifstream in(...)`, because a global locale otherwise groups integers (`1.234.567`) and misreads `1.5`. `tests/python/test_no_locale_sensitive_number_io.py` fails either mistake.
- Prefix anonymous-namespace helpers uniquely per `.cpp`: amalgamation combines them into one translation unit. Regenerate and commit `src/single_include/meshioplusplus/meshioplusplus.hpp` after any `src/cpp/` change; edit the source, not the generated header.

### ABI and bindings

Follow [doc/abi.md](doc/abi.md) for every installed C++ header change. Layout changes and changes to existing inline/template bodies or default arguments require a bump to `MESHIOPLUSPLUS_ABI_VERSION` in `abi_version.hpp`. Pure additions use a release-specific [ABI review](doc/abi_reviews.md) naming each changed header. Keep `tests/cpp/test_abi_layout.cpp` current; layout tests cannot detect inline-body changes.

C entry points must catch exceptions through the existing guards and return status/errors. Setters copy; getter borrows expire on mutation. Preserve C option-struct layouts and reserved tails. Update Fortran, Julia, R, and WASM wrappers with native API changes, including layout checks, exports (especially R's `NAMESPACE`), defaults, and reports. See the language-specific docs for indexing and ownership conventions.

### Adding a new format or operation

1. For a format, add the Python reference and registration under `src/python/meshioplusplus/<format>/`, then import/export it in the package `__init__.py`. A native implementation is optional. If added, put it in `src/cpp/{include/meshioplusplus,src}/formats/`, bind it in `_core.cpp`, register readers/writers and extension defaults in `registry.cpp`, and add the C++/Python fallback shim. In the shim, ask `core_declined(exc, "<fmt>", "read"|"write", filename)` from `_fallback.py` inside `except Exception as exc:` and re-raise when it returns `False`; never write a bare `except Exception: pass` (`tests/python/test_no_broad_core_swallow.py` fails it). Set `MESHIOPLUSPLUS_STRICT_CORE=1` to make every fallback an error and prove the core handles a file.
2. For a native operation, use `operations/` and expose it through Python, C, Fortran, Julia, R, WASM, and both CLIs as appropriate. Add pipeline support where the operation fits the pipeline contract. Python-only integrations can remain Python-only.
3. Test direct C++ behavior as well as Python dispatch: broad fallback handlers can hide a broken native implementation. Check cross-compatibility with the Python reference, optional dependencies, relevant mesh backends, and meaningful edge cases. Preserve documented numerical parity; do not silently substitute a different algorithm where parity is unavailable.
4. Follow the change checklist above, including MCP, docs, examples, and generated outputs.

Reference fixtures in `tests/python/meshes/` and `tests/python/input/` use Git LFS. Do not copy GPL fixtures into this MIT repository; generate them or use compatible licensed sources with attribution in `CITATION.cff`/`CHANGELOG.md`.

## Version bumps

Update all ten files to the same release version:

| File | Field |
| --- | --- |
| `pyproject.toml` | `version` |
| `CMakeLists.txt` | `project(... VERSION ...)` |
| `conanfile.py` | `version` |
| `packages/vcpkg/meshioplusplus/vcpkg.json` | `version` |
| `src/wasm/package.json` | `version` |
| `src/viewer/package-lock.json` | Local `@meshioplusplus/wasm` dependency; regenerate with `npm install` in `src/viewer/` using npm ≥10, do not hand-edit |
| `bindings/julia/MeshioPlusPlus/Project.toml` | `version` |
| `bindings/r/meshioplusplus/DESCRIPTION` | `Version:` |
| `src/cpp/include/meshioplusplus/version.hpp` | `MESHIOPLUSPLUS_VERSION_MAJOR/MINOR/PATCH` and `MESHIOPLUSPLUS_VERSION_STRING` |
| `bindings/c/include/meshioplusplus/meshioplusplus.h` | `MIO_VERSION_MAJOR/MINOR/PATCH` |

Also, in the same change:

- Update the `find_package(meshioplusplus X.Y.Z EXACT CONFIG REQUIRED COMPONENTS CXX)` pins in `README.md` and `doc/cpp_api.md`.
- Refresh the four `BASELINE_HASHES` in `tests/python/test_io_baseline.py`. First substitute the old version into newly written bytes and verify the old hashes: a version bump must not conceal an unrelated output regression.
- Add a dated `CHANGELOG.md` entry; prefix breaking items with **Breaking:**. A bugfix release can use a short entry.
- Regenerate the single header. Assess C++ ABI changes separately using the policy above; a release bump alone does not imply an ABI bump.

## Topic references

Use [doc/formats.md](doc/formats.md) and `doc/formats/<format>.md` for supported formats, limitations, and ordering rules. Operation details live in their named `doc/` pages.

- [Provenance](doc/provenance.md), [pipelines](doc/pipeline.md), [sequences](doc/sequences.md), [selective reads](doc/selective_read.md), [codecs](doc/codecs.md).
- [Viewer](doc/viewer.md), [interop](doc/interop.md), [GPU handoff](doc/gpu.md), [Blender](doc/blender.md), [PhysicsNeMo](doc/physicsnemo.md), [datasets](doc/datasets.md).
- Notebook setup: [C++](example/cpp/README.md), [Julia](example/julia/README.md), [R](example/r/README.md). Use each language's binding and rendering helpers.
- Generated visuals: [diagrams](doc/diagrams/README.md), [icons](doc/icons/README.md), [logo](doc/logo/README.md), and `tools/gen_doc_images.py`. Edit sources and regenerate committed assets.
