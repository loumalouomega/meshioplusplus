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

## Conda

```
conda install -c conda-forge meshioplusplus
```

## Development install

```
git clone https://github.com/<org>/meshioplusplus.git
cd meshioplusplus
pip install -e ".[all]"
```

Run the test suite with:

```
pytest tests/
```

or via tox (tests against Python 3.8 and 3.12):

```
tox
```

## Building from source (C++ core)

meshio++'s core is C++20, built through scikit-build-core + CMake when you
`pip install` from source. The optional native paths (HDF5, netCDF, zlib) are
auto-detected; the CMake options can be passed through `CMAKE_ARGS`:

```
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_HDF5=ON -DMESHIOPLUSPLUS_WITH_NETCDF=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON" \
  pip install --no-build-isolation -e .
```

### Standalone C++ build

For using the C++ library directly (without Python), two configure scripts live
in `build/`:

```
./build/configure.sh --backend OPENMP --tests --build     # Linux/macOS
build\configure.bat --backend STL --tests --build         # Windows
```

They create a CMake tree under `build/cpp-<build-type>` and print the follow-up
build/ctest commands.

### Parallelism

The C++ core parallelizes its hot loops through a compile-time-selected
backend (`meshioplusplus::parallel_for`):

```
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=STL      # default: C++17 parallel algorithms
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=OPENMP
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=TBB
-DMESHIOPLUSPLUS_PARALLEL_BACKEND=SEQ      # sequential
```

Notes:

- With GCC/libstdc++ the STL backend requires TBB (`apt install libtbb-dev`);
  when TBB is unusable, CMake warns and falls back to the sequential backend.
- MSVC's STL backend needs nothing extra; Apple's libc++ has no parallel STL
  (use OpenMP via `brew install libomp`, or SEQ).
- The design is open to new backends (Kokkos, …): one CMake branch plus one
  `#elif` block in `cpp/include/meshioplusplus/parallel.hpp`.

### Logging

The C++ core logs through `std::format`-based helpers with source locations.
Control verbosity with the `MESHIOPLUSPLUS_LOG_LEVEL` environment variable:
`debug`, `info`, `warn` (default), `error`, or `off`.
