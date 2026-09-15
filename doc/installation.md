# Installation

## Basic install

```
pip install meshioplusplus
```

The base install only requires NumPy. Most text-based formats work without any additional dependencies.

## Full install (all optional dependencies)

```
pip install meshioplusplus[all]
```

This pulls in:

| Package | Required for |
|---------|-------------|
| `h5py` | CGNS, H5M, HMF, MED, XDMF (HDF data format) |
| `netCDF4` | Exodus |

## Spack

meshio++ is also packaged for [Spack](https://spack.io), upstream in [`spack/spack-packages`](https://github.com/spack/spack-packages) ([PR #5624](https://github.com/spack/spack-packages/pull/5624)) — no local checkout of this repo needed:

```
spack install py-meshioplusplus +hdf5 +netcdf +zlib
```

This is `py-meshioplusplus` (a `PythonPackage`), the Python bindings — the Spack counterpart of `pip install meshioplusplus[all]`, built via scikit-build-core against a Spack-provided `hdf5`/`netcdf-c`/`zlib-api`, with the matching `py-h5py`/`py-netcdf4`/`py-zstandard`/`py-lz4` fallbacks pulled in per variant. A separate `meshioplusplus` package (`CMakePackage`) builds the standalone C API / Fortran module instead, for HPC codes with no Python involved — see [C API → Package managers](c_api.md#package-managers-conan-vcpkg-spack). Both take the same `+hdf5`/`+netcdf`/`+zlib`/`+zstd`/`+lz4`/`+kahip` variants (all named identically to the Conan options / vcpkg features), and both `conflicts` below GCC 10 (meshio++ needs a C++20 toolchain). Since the recipes live upstream rather than in this repo, a new meshio++ release needs no action here to become spack-installable — a `version(...)` line is added to the upstream recipe after each tag.

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

meshio++'s core is C++20, built through scikit-build-core + CMake when you `pip install` from source. The optional native paths (HDF5, netCDF, zlib) are auto-detected; the CMake options can be passed through `CMAKE_ARGS`:

```
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON" \
  pip install --no-build-isolation -e .
```

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

- `AUTO` (the default) prefers OpenMP because it is the portable choice — libgomp on manylinux, built into MSVC, libomp on macOS — and needs no TBB. It falls back to the STL backend, then to sequential.
- With GCC/libstdc++ the STL backend requires TBB (`apt install libtbb-dev`); when TBB is unusable, CMake warns and falls back to the sequential backend. This is why `AUTO` does not pick STL first: without TBB it runs sequentially.
- `_core.__parallel_backend__` reports the backend actually compiled in.
- MSVC's STL backend needs nothing extra; Apple's libc++ has no parallel STL (use OpenMP via `brew install libomp`, or SEQ).
- `KOKKOS` is bring-your-own, like KaHIP: it is never picked by `AUTO`, and CMake locates an installed [Kokkos](https://kokkos.org) via `find_package(Kokkos CONFIG)` — point `Kokkos_DIR` (or `CMAKE_PREFIX_PATH`) at the install prefix. It runs on Kokkos's **host** execution space deliberately: meshio++'s loop bodies work on host memory, so device (CUDA/HIP/SYCL) dispatch is not meaningful here — GPU data movement is what the [DLPack/CuPy handoff](gpu.md) is for. meshio++ initializes Kokkos lazily only if the embedding application hasn't already done so (apps wanting full lifecycle control should call `Kokkos::initialize()` before their first meshio++ call). Not available under Emscripten.
- The design is open to new backends (HPX, …): one CMake branch plus one `#elif` block in `src/cpp/include/meshioplusplus/parallel.hpp`.

### Logging

The C++ core logs through `std::format`-based helpers with source locations. Control verbosity with the `MESHIOPLUSPLUS_LOG_LEVEL` environment variable: `debug`, `info`, `warn` (default), `error`, or `off`.

### JavaScript / WebAssembly

The same C++ core also compiles to WebAssembly for use in the browser or Node.js, published as [`@meshioplusplus/wasm`](https://www.npmjs.com/package/@meshioplusplus/wasm) (`npm install @meshioplusplus/wasm`). Building it from source needs the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html) instead of a native compiler:

```sh
./build/configure-wasm.sh --build
node tests/wasm/smoke.mjs
```

See [WebAssembly / JavaScript](./wasm.md) for the full usage guide, the supported-format list (all 43 formats, the HDF5- and netCDF-backed ones included since v8.0.0), and known limitations.
