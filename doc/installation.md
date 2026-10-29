# Installation

## Basic install

```
pip install meshioplusplus
```

The base install requires NumPy and Rich (the CLI's output library). Most text-based formats work without any additional dependencies.

## The `[all]` extra

```
pip install meshioplusplus[all]
```

This pulls in:

| Package | Required for |
|---------|-------------|
| `h5py` | CGNS, H5M, HMF, MED, XDMF (HDF data format) |
| `netCDF4` | Exodus |
| `polyscope` | the desktop viewer backend of `view()` and `screenshot()` |

The PyPI wheels are built without HDF5 and netCDF, so these two extras are how the HDF5 and netCDF formats work there. Features that only the native HDF5 build has, such as polyhedral CGNS sections, need a source build with HDF5 (see [CGNS](formats/cgns.md)).

## What the PyPI wheels contain

The wheels link these libraries statically into `meshioplusplus._core`:

| Component | Licence |
|---|---|
| OpenCASCADE 7.6.3 (`fit_cad`, STEP and IGES reading) | LGPL-2.1 with the OCCT exception |
| zlib 1.3.1 (vendored; VTU compression and the GiD writer) | Zlib |
| gidpost 2.14 (vendored; the GiD writer) | BSD-2-Clause-Views |
| pugixml 1.14 (vendored; XML parsing) | MIT |
| nlohmann/json 3.12.0 (the JSON pipeline parser) | MIT |

They leave out HDF5 and netCDF (the `[all]` extra provides h5py and netCDF4), zstd, lz4 and bzip2 (the `[codecs]` extra, or a source build), ADIOS2, TecIO, KaHIP and MMG (a source build; KaHIP also has the `[kahip]` extra), Polyscope (the `[viewer]` extra) and Eigen, which is not compiled without HDF5.

Linux wheels are built with the SEQ parallel backend. macOS and Windows wheels use OpenMP; the macOS build uses Homebrew's libomp, which delocate is expected to bundle into the wheel (to be confirmed in a wheel listing).

The wheels are licensed `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`, as [`LICENSE-BINARIES.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/LICENSE-BINARIES.md) states. Each wheel ships `LICENSE`, `LICENSE-BINARIES.md`, `THIRD_PARTY_NOTICES.md` and the OpenCASCADE licence texts in its `.dist-info/licenses` directory.

## Spack

meshio++ is also packaged for [Spack](https://spack.io), upstream in [`spack/spack-packages`](https://github.com/spack/spack-packages) ([PR #5624](https://github.com/spack/spack-packages/pull/5624)) — no local checkout of this repo needed:

```
spack install py-meshioplusplus +hdf5 +netcdf +zlib
```

This is `py-meshioplusplus` (a `PythonPackage`), the Python bindings — the Spack counterpart of `pip install meshioplusplus[all]`, built via scikit-build-core against a Spack-provided `hdf5`/`netcdf-c`/`zlib-api`, with the matching `py-h5py`/`py-netcdf4`/`py-zstandard`/`py-lz4` fallbacks pulled in per variant. A separate `meshioplusplus` package (`CMakePackage`) builds the standalone C API / Fortran module instead, for HPC codes with no Python involved — see [C API → Package managers](c_api.md#package-managers-conan-vcpkg-spack). Both take the same `+hdf5`/`+netcdf`/`+zlib`/`+zstd`/`+lz4`/`+bzip2`/`+kahip`/`+cgnslib`/`+adios2` variants (all named identically to the Conan options / vcpkg features), and both `conflicts` below GCC 10 (meshio++ needs a C++20 toolchain). The recipes live upstream rather than in this repo, so a new meshio++ release is spack-installable only once someone bumps them: that bump is part of the release checklist in [AGENTS.md](https://github.com/loumalouomega/meshioplusplus/blob/master/AGENTS.md#version-bumps) (after the tag, `spack checksum` and one `version(...)` line per recipe). Until it lands, pin the last version the recipes carry — they trail the tags by design, not by accident.

## Development install

```
git clone https://github.com/loumalouomega/meshioplusplus.git
cd meshioplusplus
pip install -e ".[all]"
```

Run the test suite with:

```
pytest tests/python/
```

or via tox (tests against Python 3.9 and 3.12):

```
tox
```

## Building from source (C++ core)

meshio++'s core is C++20, built through scikit-build-core + CMake when you `pip install` from source. The optional native paths (HDF5, netCDF, and zlib, which the wheels vendor instead) are auto-detected; the CMake options can be passed through `CMAKE_ARGS`:

```
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON" \
  pip install --no-build-isolation -e .
```

Off by default and never auto-enabled: `-DMESHIOPLUSPLUS_WITH_ZSTD=ON` / `_LZ4=ON` add the VTK XML codecs of [compression codecs](codecs.md), and `-DMESHIOPLUSPLUS_WITH_BZIP2=ON` (v16.12.0; libbz2 through CMake's `FindBZip2`) lets the core read and write [libMesh](formats/libmesh.md)'s `.bz2` meshes, which otherwise go through Python's `bz2` (`_core.__has_bzip2__` reports it; Conan `with_bzip2`, vcpkg feature `bzip2`).

Two read-only formats need a library meshio++ never builds by default (v16.13.0). `-DMESHIOPLUSPLUS_WITH_ADIOS2=ON` links ADIOS2 (its own CMake config, `adios2::cxx` — `adios2::cxx11` before ADIOS2 2.9; set `ADIOS2_DIR` if it is not found) so the core reads [DOLFINx VTX](formats/vtx.md) `.bp` directories; without it they go through the `adios2` Python package (`pip install meshioplusplus[adios2]`); `_core.__has_adios2__` reports it (Conan `with_adios2`, bring-your-own; vcpkg feature `adios2`). `-DMESHIOPLUSPLUS_WITH_TECIO=ON -DTECIO_ROOT=<dir>` links a TecIO you already have (`<dir>` holding `TECIO.h` and `libtecio`) so the core reads [Tecplot SZL](formats/szplt.md) `.szplt` files; without it they go through a shared TecIO named by `MESHIOPLUSPLUS_TECIO_LIBRARY`; `_core.__has_tecio__` reports it (Conan `with_tecio`, bring-your-own). `build/configure.sh` takes `--with-adios2`, `--with-tecio <dir>` and `--with-bzip2`. The release wheels carry neither library.

OpenCASCADE is on by default (`MESHIOPLUSPLUS_WITH_OCCT=ON`, auto-detected like HDF5): a build links an OpenCASCADE 7.6 or later it finds through its own CMake config (`OpenCASCADE_ROOT`, `OpenCASCADE_DIR`), so [`fit_cad`](fit_cad.md) builds trimmed, sewn and checked CAD faces and STEP/IGES read. On Ubuntu `apt-get install libocct-foundation-dev libocct-modeling-data-dev libocct-modeling-algorithms-dev libocct-data-exchange-dev libtbb-dev` gives a shared one; `build/build-occt.sh PREFIX` builds the pinned 7.6.3 static, as the PyPI wheels and the release CLI binaries link it. `-DMESHIOPLUSPLUS_REQUIRE_OCCT=ON` makes a missing kernel a configure error, and `--without-occt` (`-DMESHIOPLUSPLUS_WITH_OCCT=OFF`) leaves it out; without it `fit_cad` falls back to its kernel-free native backend with a warning and STEP/IGES reading raises naming the flag. OpenCASCADE is LGPL-2.1 with the OCCT exception: the published binaries carry it statically under [`THIRD_PARTY_NOTICES.md`](https://github.com/loumalouomega/meshioplusplus/blob/master/THIRD_PARTY_NOTICES.md). `_core.__has_occt__` reports it (Conan `with_occt`, on by default and auto-detected; vcpkg feature `occt`, a default feature).


The [MMG backend](mmg.md) is off by default and never in the wheels: `-DMESHIOPLUSPLUS_WITH_MMG=ON -Dmmg_DIR=<prefix>/lib/cmake/mmg` (`build/configure.sh --with-mmg --mmg-root <prefix>`) links a shared MMG install, which `build/build-mmg.sh <prefix>` builds (MMG 5.8.0; `--with-elas` adds ELAS for `move`, which makes the build GPL-3.0). MMG is LGPL-3.0; see its [licence section](mmg.md#licence). Put `<prefix>/lib` on `LD_LIBRARY_PATH` when importing a Python build, whose RPATH scikit-build-core strips.
For development, `-DMESHIOPLUSPLUS_SANITIZE="address;undefined"` (`build/configure.sh --sanitize address,undefined`, GCC or Clang) builds the core, the C API and the tests with the sanitizers, and `-DMESHIOPLUSPLUS_BUILD_FUZZERS=ON` (`--fuzzers`, Clang) the libFuzzer reader target; both refuse the Python extension and configure a separate `-san` tree ([fuzzing and sanitizers](fuzzing.md)).

### Standalone C++ build

For using the C++ library directly (without Python), two configure scripts live in `build/`:

```
./build/configure.sh --backend OPENMP --tests --build     # Linux/macOS
build\configure.bat --backend STL --tests --build         # Windows
```

They create a CMake tree under `build/cpp-<build-type>` and print the follow-up build/ctest commands.

### Mesh backends

Standalone C++ builds can swap the in-memory mesh structure itself with `--mesh-backend` (CMake: `MESHIOPLUSPLUS_MESH_BACKEND`):

```
./build/configure.sh --mesh-backend NATIVE --tests --build   # fastest pure-C++ structure
./build/configure.sh --mesh-backend KRATOS --tests --build   # Kratos-style ModelPart
```

- `MESHIO` (default) — mirrors the Python `meshio.Mesh`; **required** when the pybind11 extension is built (PyPI wheels always use it).
- `NATIVE` — canonical Float64/Int64 storage, `CellType` enum, CSR ragged blocks; the WebAssembly build uses it.
- `KRATOS` — a Kratos-Multiphysics-style `ModelPart` behind the same API, with a header-only bridge to the real `Kratos::ModelPart`.

All formats work identically under every backend. See [C++ mesh backends](cpp_backends.md) for the full story.

### Parallelism

The C++ core parallelizes its hot loops through a compile-time-selected backend (`meshioplusplus::parallel_for`):

```
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=AUTO     # default: OpenMP, else STL(+TBB), else SEQ
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=STL      # C++17 parallel algorithms
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=TBB
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=KOKKOS   # Kokkos host execution space (bring-your-own)
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ      # sequential
```

Notes:

- `AUTO` (the default) prefers OpenMP because it is the portable choice — libgomp on manylinux (the Linux wheels are built SEQ instead), built into MSVC, libomp on macOS — and needs no TBB. It falls back to the STL backend, then to sequential.
- With GCC/libstdc++ the STL backend requires TBB (`apt install libtbb-dev`); when TBB is unusable, CMake warns and falls back to the sequential backend. This is why `AUTO` does not pick STL first: without TBB it runs sequentially.
- `_core.__parallel_backend__` reports the backend actually compiled in.
- MSVC's STL backend needs nothing extra; Apple's libc++ has no parallel STL (use OpenMP via `brew install libomp`, or SEQ).
- `KOKKOS` is bring-your-own, like KaHIP: it is never picked by `AUTO`, and CMake locates an installed [Kokkos](https://kokkos.org) via `find_package(Kokkos CONFIG)` — point `Kokkos_DIR` (or `CMAKE_PREFIX_PATH`) at the install prefix. It runs on Kokkos's **host** execution space deliberately: meshio++'s loop bodies work on host memory, so device (CUDA/HIP/SYCL) dispatch is not meaningful here — GPU data movement is what the [DLPack/CuPy handoff](gpu.md) is for. meshio++ initializes Kokkos lazily only if the embedding application hasn't already done so (apps wanting full lifecycle control should call `Kokkos::initialize()` before their first meshio++ call). Not available under Emscripten.
- The design is open to new backends (HPX, …): one CMake branch plus one `#elif` block in `src/cpp/include/meshioplusplus/parallel.hpp`.

### Logging

The C++ core logs through `std::format`-based helpers with source locations. Control verbosity with the `MESHIOPLUSPLUS_LOG_LEVEL` environment variable: `debug`, `info`, `warn` (default), `error`, or `off`.

The Python layer logs through the standard `logging` module under the logger name `meshioplusplus`. It never prints to stdout, so a CLI pipeline's output stays clean. At `DEBUG` it reports every ambiguous-extension candidate that declined a file and every C++ fast path that fell back to its Python twin; at `WARNING` it reports a fast path that failed unexpectedly. Set `MESHIOPLUSPLUS_STRICT_CORE=1` to turn every such fallback into an error instead; see [when the native path declines](./formats.md#when-the-native-path-declines).

### JavaScript / WebAssembly

The same C++ core also compiles to WebAssembly for use in the browser or Node.js, published as [`@meshioplusplus/wasm`](https://www.npmjs.com/package/@meshioplusplus/wasm) (`npm install @meshioplusplus/wasm`). Building it from source needs the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html) instead of a native compiler:

```sh
./build/configure-wasm.sh --build
node tests/wasm/smoke.mjs
```

See [WebAssembly / JavaScript](./wasm.md) for the full usage guide, the supported-format list (79 readable and 74 writable formats, the HDF5- and netCDF-backed ones included since v8.0.0), and known limitations.
