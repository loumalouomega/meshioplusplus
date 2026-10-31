# Split Python packages: plan

The design pass for roadmap 4.7: one MIT base wheel, `meshioplusplus`, and add-on wheels that carry the copyleft libraries, `meshioplusplus-occt`, `meshioplusplus-netgen` and `meshioplusplus-mmg`, each a small extension linked against the base wheel's shared core, in the way KratosMultiphysics ships its applications beside `KratosMultiphysics` with an everything extra on top. v18.0.0 built every phase except Windows for the permissive libraries and for plugins; [What shipped in v18.0.0](#what-shipped-in-v18-0-0) records how the build differs from the design below, and the roadmap holds what is left (4.7).

## The problem

The PyPI wheel is one extension, `meshioplusplus._core`, with OpenCASCADE 7.6.3 linked into it statically (`wheels.yml`, `-DMESHIOPLUSPLUS_REQUIRE_OCCT=ON`). Every wheel is therefore `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`, including for a user who never reads a STEP file, and every redistributor inherits the LGPL's relinking obligation, which the project meets only by documenting a rebuild from source ([third-party notices](https://github.com/loumalouomega/meshioplusplus/blob/main/THIRD_PARTY_NOTICES.md)). Netgen and MMG are never in a wheel at all, so `tetrahedralize` keeps no surface and `adapt`, `discretize_levelset`, `mesh_domain` and `move` raise for every pip user. Meanwhile libraries with no licence conflict (HDF5, netCDF, CGNS, zstd, lz4, bzip2, KaHIP) are left out too, and the formats they serve run on Python fallbacks.

## The target

| Wheel | Carries | Licence | Depends on |
| --- | --- | --- | --- |
| `meshioplusplus` | `_core`, the shared core `libmeshioplusplus_core`, its headers and CMake package; zlib, pugixml, nlohmann/json, gidpost, then HDF5, netCDF, CGNS, zstd, lz4, bzip2, KaHIP and Eigen (phase 5) | MIT (MIT AND MPL-2.0 once Eigen is in) | numpy, rich |
| `meshioplusplus-occt` | `meshioplusplus_occt._occt`: STEP and IGES reading, `fit_cad`'s OpenCASCADE backend, OpenCASCADE as shared libraries | MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0) | `meshioplusplus==X.Y.Z` |
| `meshioplusplus-netgen` | `meshioplusplus_netgen._netgen`: `tetrahedralize`'s Netgen backend, nglib and ngcore as shared libraries | MIT AND LGPL-2.1-only | `meshioplusplus==X.Y.Z` |
| `meshioplusplus-mmg` | `meshioplusplus_mmg._mmg`: `adapt`, `discretize_levelset`, `mesh_domain`, `move` (without ELAS), libmmg as a shared library | MIT AND LGPL-3.0-or-later | `meshioplusplus==X.Y.Z` |

The base package's extras then select add-ons: `meshioplusplus[occt]`, `[netgen]`, `[mmg]`, and `[addons]` installing all three. Because each add-on bundles its LGPL library as a shared object next to its own extension, a recipient can swap that library without rebuilding anything, which is the LGPL's preferred form of compliance, and a user who installs only the base wheel installs no copyleft code at all.

Every wheel ships the same version, and an add-on pins the base exactly (`meshioplusplus==X.Y.Z`), so the C++ interface between them only has to be stable *within* one release. That is what Kratos relies on, and it is why the project's ABI policy ([C++ ABI](./abi.md), 36 bumps since v9.0.0) does not need to become stricter for this.

## How the pieces fit

### The provider registry (core)

The core stops testing `MESHIOPLUSPLUS_HAS_OCCT`, `HAS_NETGEN` and `HAS_MMG` and asks a registry instead. A new installed header, `meshioplusplus/providers.hpp`, holds one typed slot per backend:

```cpp
struct CadKernelProvider {
    std::string mName, mVersion;                       // "occt", OCC_VERSION_COMPLETE
    std::function<Mesh(const std::string&, const CadReadOptions&)> mReadStep, mReadIges;
    std::function<FitCadResult(const Mesh&, const std::string&, const FitCadOptions&)> mFitCad;
};
struct TetMesherProvider { std::string mName, mVersion; std::function<TetrahedralizeResult(...)> mMesh; };
struct RemesherProvider  { std::string mName, mVersion; bool mElas; int mIndexBits;
                           std::function<...> mAdapt, mDiscretizeLevelset, mMeshDomain, mMove; };

MESHIOPLUSPLUS_API void register_cad_kernel(CadKernelProvider);
MESHIOPLUSPLUS_API const CadKernelProvider* cad_kernel() noexcept;   // nullptr: none
// ... the same pair for the tetrahedral mesher and the remesher
```

Rules for the registry:

- **The core's control flow keeps its shape.** `cad_read_has_occt()`, `fit_cad_has_occt()`, `tetrahedralize_has_netgen()`, `mmg_available()` and `mmg_has_elas()` keep their signatures and read the registry. `read_step` dispatches to `cad_kernel()->mReadStep` or throws the `Unsupported` it throws today (naming the add-on instead of a CMake flag). `fit_cad` and `tetrahedralize` fall back to their native backends exactly as now.
- **One registration path for every build.** The monolithic builds (source builds, the native CLI, Conan, vcpkg, Spack, both WebAssembly artifacts) compile the providers into the core and register them through an explicit `register_builtin_providers()`. The core calls it once, lazily, from each capability query: no static initializers, and so no initialization-order fragility. A split build compiles none of them in, and the add-on registers its own when it is imported.
- **Formats register at runtime too.** `registry.cpp`'s tables become mutable singletons behind a lock, with `registry_register_reader(name, ReadFn, ReadExFn, MetadataFn)` and the writer and extension counterparts. The existing `static const` initializer lists become the first registrations, so the iteration order, and every test that pins it, stays the same. `step` and `iges` stay registered as stubs that throw `Unsupported`, as they do now, and a provider replaces them.
- **The pipeline gains no hook.** `FitCad`, `Tetrahedralize`, `Adapt` and the other steps already call the C++ functions, which reach the provider, so `pipe_op_table()` and the dispatch chain stay as they are.
- **What moves out of the core.** `detail/occt_brep.*`, the OpenCASCADE branch of `formats/cad_read.cpp`, the OpenCASCADE backend of `operations/fit_cad.cpp` (`fc_occt_fit_cad`, about 4,400 lines), the Netgen branch of `operations/tetrahedralize.cpp` and the MMG implementation in `operations/mmg.cpp` move to `src/cpp/providers/{occt,netgen,mmg}/`. The core keeps the option structs, the native backends and the dispatch.

### The shared core (base wheel)

`_core` links a shared `libmeshioplusplus_core` instead of embedding `meshioplusplus_core_obj`. The library already exists: `MESHIOPLUSPLUS_INSTALL_CPP` builds `meshioplusplus_core_<backend>` with hidden visibility, `MESHIOPLUSPLUS_API` exports and `SOVERSION ${MESHIOPLUSPLUS_ABI_VERSION}`. The wheel build turns on the MESHIO variant of it and installs, under `meshioplusplus/`:

- `lib/` (`.so.36`, `.dylib`, or `.dll` with its import library);
- `include/meshioplusplus/`;
- `lib/cmake/meshioplusplus/`, so an add-on finds the core by `find_package(meshioplusplus CONFIG)` against `meshioplusplus.get_cmake_dir()`, as pybind11's `get_cmake_dir()` works.

Two consequences need care:

- **One C++ runtime.** Two shared objects now pass `std::string`, `std::function` and exceptions between them, so neither may carry a private copy of libstdc++: `MESHIOPLUSPLUS_STATIC_RUNTIME`'s `-static-libstdc++` and `--exclude-libs,ALL` go away on Linux in favour of the manylinux system libstdc++. On Windows the CRT is the dynamic `/MD`, which it already is for `_core`. The exception classes the add-ons throw or catch (`ReadError`, `WriteError`, `Unsupported`) are exported with `MESHIOPLUSPLUS_API` so that their type identity crosses the boundary.
- **The core-internal surface.** The provider sources call `detail::` helpers that no installed header declares today (fit CAD's patch network, the B-spline fitters, `face_mesh`). They are exported under a `MESHIOPLUSPLUS_INTERNAL_API` macro and documented as valid only within one release, which the exact version pin guarantees. Compiling private copies of them into each add-on instead would break the one-definition rule across shared objects.

### The add-on wheels

Each add-on is its own project under `packages/python/meshioplusplus-<name>/`, with a `pyproject.toml` whose scikit-build-core build points `cmake.source-dir` at the repository root with `-DMESHIOPLUSPLUS_ADDON=<name>`. Its build requirements are `meshioplusplus==X.Y.Z`, pybind11 and the library itself (OpenCASCADE from `build/build-occt.sh --shared`, Netgen and MMG from new `build/build-netgen.sh` and the existing `build/build-mmg.sh`). The extension module is about a hundred lines:

```cpp
PYBIND11_MODULE(_occt, m) {
    meshioplusplus::providers::register_occt();     // the CAD kernel and the step/iges readers
    m.attr("__version__") = MESHIOPLUSPLUS_VERSION_STRING;
    m.attr("__occt_version__") = OCC_VERSION_COMPLETE;
}
```

Meshes never cross a Python boundary: the base `_core` converts the arguments, calls `fit_cad`, and the core reaches the provider in-process. The add-on's Python package is just its `__init__.py`, which imports `_occt`, and an entry point:

```toml
[project.entry-points."meshioplusplus.plugins"]
occt = "meshioplusplus_occt"
```

Linking against the core, and repairing each platform's wheel:

| Platform | The add-on's extension finds the core through | The repair tool |
| --- | --- | --- |
| Linux | an RPATH of `$ORIGIN/../meshioplusplus/lib` | `auditwheel repair --exclude libmeshioplusplus_core.so.36` bundles OpenCASCADE but not the core |
| macOS | an RPATH of `@loader_path/../meshioplusplus/lib` | `delocate-wheel --exclude libmeshioplusplus_core` |
| Windows | the base package's `__init__.py` calling `os.add_dll_directory()` on its `lib` directory before importing `_core` | `delvewheel repair --exclude meshioplusplus_core.dll` |

### Discovery (Python)

A new `meshioplusplus/_plugins.py` loads the entry points of the group `meshioplusplus.plugins` once, lazily. It is called from the capability queries (`has_occt()`, a new public `has_netgen()`, `has_mmg()`) and from the paths that need a provider (`step._read_cad`, `_fit_cad`, `_tetrahedralize`, `_mmg`), and checks each add-on's `__version__` against the base, refusing a mismatch by name. `MESHIOPLUSPLUS_NO_PLUGINS=1` skips discovery, so a test run can prove the base wheel alone.

The `_core.__has_occt__`, `__has_netgen__` and `__has_mmg__` attributes stay for compatibility. Read after discovery, they report the registry rather than a compile-time constant, through a module-level `__getattr__` on the pybind module's Python wrapper.

## Phases

Each phase lands on its own, with the monolithic builds green and unchanged in behaviour throughout.

0. **Spikes (S).** Each spike answers one question by experiment before the code it would change is written:
   - Exceptions and `std::function` crossing two shared objects built in the manylinux_2_34 image with gcc-toolset-13 (whose libstdc++ links part of itself statically into each object), and on macOS and Windows.
   - Whether `netgen-mesher`'s PyPI wheels expose a linkable nglib, which would let `meshioplusplus-netgen` depend on them instead of bundling.
   - A scikit-build-core project building against an installed base wheel's CMake package.
1. **The provider registry (M).**
   - `providers.hpp`, `register_builtin_providers()` and the runtime format registration.
   - The provider sources moved under `src/cpp/providers/`, still compiled into the monolithic core.
   - Exit: every existing test passes in every CI configuration (OpenCASCADE, Netgen and MMG jobs included), `fit_cad`'s goldens are unchanged, and a new C++ test registers a fake provider of each kind and sees it used.
2. **The shared core in the base wheel (M).**
   - `_core` linked against `libmeshioplusplus_core`, with headers, the CMake package and `get_cmake_dir()` in the wheel, built with OpenCASCADE OFF.
   - The repair configuration for all three platforms.
   - Exit: the base wheel's licence is `MIT`, `nm`/`dumpbin` finds no OpenCASCADE symbol in it, and the Python suite passes with the OpenCASCADE tests skipping by name.
   - **Breaking:** pip users lose STEP/IGES reading and the OpenCASCADE `fit_cad` until they install the add-on. This phase and phase 3 therefore ship together in one major release (v18.0.0).
3. **`meshioplusplus-occt` (L).** OpenCASCADE built shared by `build-occt.sh`, the add-on project, its entry point and discovery, and a `wheels.yml` job building it on the same runner as the base and testing the pair. Exit: `tests/wheel_occt_smoke.py` passes with base plus add-on installed, and the suite's OpenCASCADE tests run instead of skipping.
4. **`meshioplusplus-netgen` and `meshioplusplus-mmg` (M each).** The same template. MMG ships without ELAS (ELAS makes a build GPL-3.0), so `move` raises naming it; Netgen is 6.2.2401, the version the native CI job and the full WebAssembly artifact use. Exit: the `tetrahedralize` and MMG suites run against the installed wheels.
5. **The permissive libraries in the base wheel (M).** HDF5, netCDF and CGNS built static with hidden symbols (h5py brings its own HDF5 into the same process, so none of ours may be exported), zstd, lz4, bzip2, KaHIP and Eigen. The native HDF5-backed formats then run in the wheel instead of their Python fallbacks. Exit: the conformance matrix is regenerated from a wheel, and h5py and our core read the same file in one process.
6. **Extras, docs and release (S).**
   - The extras, and the version-bump table in `AGENTS.md` gaining the add-on `pyproject.toml` files.
   - `LICENSE-BINARIES.md` and `THIRD_PARTY_NOTICES.md` per wheel; [installation](./installation.md) rewritten around the split, with a migration note.
   - The PyPI trusted publishers for the three new projects.
7. **Native consumers, optional (M).** The release CLI binaries and the C API could load the same providers as shared plugins (`MESHIOPLUSPLUS_PLUGIN_PATH`, an exported `extern "C" meshioplusplus_register_plugin()`), so the CLI archives could split too. Until then they keep OpenCASCADE linked statically, and their licence is unchanged.

Altogether **L to XL**: about a month for phases 0 to 3, which deliver the MIT base wheel and OpenCASCADE as an add-on, and as much again for the rest.

## Risks

- **Two copies of a C++ runtime.** If spike 0.1 shows the toolset's partly static libstdc++ cannot pass exceptions between the core and an add-on, the fallback is to build the add-ons in the same CI job as the base, from the same toolchain, and to require exactly matching wheel tags. Kratos does the same.
- **h5py and our HDF5 in one process (phase 5).** This is safe only if no HDF5 symbol of ours is visible. The phase's exit test exists for this reason.
- **Users who never pin.** `pip install -U meshioplusplus` without the add-on leaves an old `meshioplusplus-occt` behind. Discovery refuses the mismatch by name rather than crashing, and the add-on's exact pin makes pip upgrade it alongside whenever the add-on is installed through the extra.
- **Windows DLL resolution.** `os.add_dll_directory` must run before any extension import. It goes in the first lines of `meshioplusplus/__init__.py`, and an import test runs on a clean Windows runner.

## What stays the same

- **Source and monolithic builds.** `build/configure.sh`, Conan, vcpkg and Spack build one core with whatever they find, as now.
- **Both WebAssembly artifacts.** They link their providers statically; the lean/full split is already the browser's version of this design.
- **The Python API.** No function is added or removed by the split. A capability that is missing raises the same `NotImplementedError` or `Unsupported` it raises in a source build without the library, now naming the add-on to install.

## Decided here

- **The everything option is a separate `[addons]` extra, not `[all]`.** The first draft put the add-ons in `[all]`. It cannot work: an add-on pins the base wheel's exact version, so `pip install -e ".[all]"` from a source checkout (every contributor, and the benchmark and fuzz workflows) cannot resolve them until they are on PyPI at that version. `[all]` stays h5py, netCDF4 and polyscope, and `[addons]` installs the three LGPL wheels. A separate `meshioplusplus-all` metapackage would add a release artifact without adding a choice.
- **The add-on module names are `meshioplusplus_occt`, `meshioplusplus_netgen` and `meshioplusplus_mmg`.** They are separate top-level packages rather than `meshioplusplus.occt` namespace subpackages, because a namespace package cannot share the base's regular package directory.

## What shipped in v18.0.0

Phases 0 to 3 landed together, as planned. Where the build differs from the design above:

- **The registry holds the CAD kernel only.** `providers.hpp` has `CadKernelProvider` with `register_cad_kernel()`, `unregister_cad_kernel()`, `cad_kernel()` and `register_builtin_providers()`. The tetrahedral-mesher and remesher slots, and the runtime format registration, come with the Netgen and MMG add-ons: `step` and `iges` were already registered stubs that dispatch to `read_step` and `read_iges`, so the OpenCASCADE add-on needs no new registry entry. `cad_kernel()` returns a `shared_ptr` and the registry retains every provider it ever held, so the `const char*` of `fit_cad_occt_version()` stays valid across a replacement.
- **The provider sources moved** to `src/cpp/providers/occt/` (`occt_brep.*`, `occt_provider.*`, `fit_cad_occt.cpp`, and the add-on's module in `addon/`). The OpenCASCADE backend of `fit_cad` reaches the native backend's helpers (the freeform networks, the refinement, the face table, the report) through `src/cpp/src/detail/fit_cad_native.hpp`, which exports eight functions rather than the whole `detail::` surface.
- **The shared core exports every symbol** (default visibility, `WINDOWS_EXPORT_ALL_SYMBOLS`) instead of tagging each helper with a `MESHIOPLUSPLUS_INTERNAL_API` macro. The provider sources call dozens of `detail::` helpers that no installed header declares, and the exact version pin is what makes that safe. A Windows export table has a 65,535-symbol limit, and the first Windows wheel build confirmed it (LNK1189): on Windows `cmake/gen_core_def.py` writes the export table from the compiled objects, keeping the `meshioplusplus` namespace.
- **The C++ runtime stays static on Linux, and exported.** Spike 0.1 measured the design's first risk. Two shared objects that each link libstdc++ statically abort on the first exception that crosses between them when `--exclude-libs,ALL` hides the runtime's symbols, which is what the single-extension wheels used; with the symbols exported, the first-loaded copy serves every object and exceptions, `std::function` and `std::string` cross correctly. The shared-core mode therefore keeps `-static-libstdc++` (manylinux's libstdc++ is older than gcc-toolset-13's) and drops `--exclude-libs`. The spike was a toy program on one toolchain; the wheel build on the CI runners is its real test.
- **Discovery** is `meshioplusplus/_plugins.py`. `_core.__has_occt__` is read live through a module `__getattr__`, so an add-on registered after `_core` is imported is seen.
- **`get_cmake_dir()`** is in `meshioplusplus/__init__.py` and is not in `__all__`: it describes the install, not a mesh operation.
- **Add-on build.** `MESHIOPLUSPLUS_ADDON=occt` makes the root `CMakeLists.txt` include `cmake/addon_occt.cmake` and stop; the add-on's `pyproject.toml` sets it, with `cmake.source-dir` at the repository root. `build/build-occt.sh PREFIX --shared` builds the shared libraries it links.
- **Netgen and MMG are providers too** (phase 4). `providers.hpp` gained `TetMesherProvider` and `RemesherProvider` beside `CadKernelProvider`. `tetrahedralize` keeps its validation, region bookkeeping, native octree route and statistics, and calls the mesher once per closed region with that region's points and triangles (`TetMeshRegionRequest`), taking back its points and tetrahedra; the nglib code, its mutex and the boundary-preservation check of the points and facets live in `src/cpp/providers/netgen`. The MMG implementation moved to `src/cpp/providers/mmg` whole, with `adapt`, `discretize_levelset`, `mesh_domain` and `move` behind `RemesherProvider`; `metric_from_hessian` and the signal guard stay in the core. `cmake/MmgProbe.cmake` (library file, 64-bit indices, ELAS) is shared by the root build, the plugins and the add-on. The builtin registration is per slot, so unregistering one kind never suppresses another's. `_core.__has_netgen__`, `__has_mmg__`, `__has_mmg_elas__` and `__mmg_index_bits__` are live like `__has_occt__`, and `has_netgen()` joins `has_mmg()` and `has_occt()`.
- **The Netgen and MMG wheels** follow the OpenCASCADE template (`cmake/addon_{netgen,mmg}.cmake`, `packages/python/meshioplusplus-{netgen,mmg}`, `wheels.yml` jobs). `build/build-netgen.sh` builds Netgen 6.2.2401 shared (the pin is `build/deps-pins.sh`, shared with the WebAssembly build); the MMG wheel uses `build/build-mmg.sh` without `--with-elas`, so `move` raises naming it and the wheel stays LGPL-3.0 rather than GPL-3.0.
- **The permissive libraries** (phase 5): `build/build-permissive-libs.sh` builds bzip2, zstd, lz4, HDF5, netCDF, cgnslib and KaHIP as static, position-independent, hidden-visibility archives; `MESHIOPLUSPLUS_PREFER_STATIC_DEPS` makes CMake look the archives up first, and on Linux the shared core links them with `--exclude-libs` naming each archive (`MESHIOPLUSPLUS_HIDDEN_ARCHIVES`, never `ALL`, which would hide libstdc++). macOS needs no hiding: its two-level namespace binds each symbol to its own image. Eigen is a header-only submodule and makes the base wheel `MIT AND MPL-2.0`. The wheel smoke test checks that no dependency symbol is exported and that h5py and netCDF4 read what the core wrote, in one process. Windows is the one platform left on the Python fallbacks.
- **Plugins** (phase 7): `MESHIOPLUSPLUS_PLUGIN_PATH` and `load_plugin`/`load_plugins`/`mio_load_plugins`, as [the C++ API page](./cpp_api.md#providers-and-plugins) describes. The scan runs once from `register_builtin_providers()`, after the compiled-in providers and outside the registry lock; a plugin fills only an empty slot. Plugins are built without linking the core, so the host exports it: the CLI and the test binary are linked `-rdynamic` (`ENABLE_EXPORTS`), and the C API library promotes itself to global scope with `dladdr` and `dlopen(RTLD_NOLOAD | RTLD_GLOBAL)`. Windows is not supported.
