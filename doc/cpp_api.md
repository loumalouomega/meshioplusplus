# C++ API

The full C++ core — `Mesh`, the [format registry](/formats), all [mesh operations](/operations) and [data operations](/data_operations), and the header-only [Kratos bridge](/cpp_backends#kratos) — installs as a normal CMake package with exported targets. It is the right entry point for a C++ consumer that wants more than the flat [C API](/c_api) can express: real `Mesh` objects, the operations layer, and `meshioplusplus::ModelPart`, none of which cross a C ABI.

::: tip Which API do I want?
Use the [C API](/c_api) if you are writing C or Fortran, want a stable ABI, or want the smallest possible dependency surface. Use this one if you are writing C++ and want the library's real types and operations.
:::

## Building and installing

Off by default, so `pip install .` and the wheels pay nothing for it. From the repo, the configure helper is the short spelling:

```sh
build/configure.sh --install-cpp --c-api --build    # --c-api optional: both APIs in one prefix
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

(`--cpp-backends MESHIO,KRATOS` trims the backend set; `configure.bat` takes the same flags on Windows, with a multi-entry list quoted). Or drive CMake directly:

```sh
cmake -S . -B build \
  -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
  -DMESHIOPLUSPLUS_INSTALL_CPP=ON \
  -DMESHIOPLUSPLUS_BUILD_C_API=ON      # optional: both APIs in one prefix
cmake --build build
cmake --install build --prefix /opt/meshioplusplus
```

which lays out:

```
include/meshioplusplus/           # the whole public header tree
  mesh.hpp  registry.hpp  kratos_bridge.hpp  export.hpp  ...
  backends/  detail/  formats/  operations/
  meshioplusplus.h                # the C header, when BUILD_C_API=ON
lib/libmeshioplusplus_core_meshio.{a,so}
lib/libmeshioplusplus_core_native.{a,so}
lib/libmeshioplusplus_core_kratos.{a,so}
lib/libmeshioplusplus.so          # the C API, when BUILD_C_API=ON
lib/cmake/meshioplusplus/         # find_package(meshioplusplus)
lib/pkgconfig/meshioplusplus-cxx.pc
bin/meshioplusplus                # the CLI, when BUILD_CLI=ON
```

`detail/` is installed on purpose: the public headers include it transitively, so an install without it does not compile.

## Consuming it

```cmake
find_package(meshioplusplus 9.0 CONFIG REQUIRED COMPONENTS CXX)
target_link_libraries(my_solver PRIVATE meshioplusplus::core)
```

```cpp
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/operations/partition.hpp"

auto mesh  = meshioplusplus::registry_read("bracket.msh", "", {});
auto parts = meshioplusplus::partition(mesh, {.mNParts = 8});
```

The two components are independent, and one `find_package` serves both:

| Component | Targets | Needs |
| --- | --- | --- |
| `CXX` | `meshioplusplus::core`, `::core_meshio`, `::core_native`, `::core_kratos` | `MESHIOPLUSPLUS_INSTALL_CPP=ON` |
| `C` | `meshioplusplus::meshioplusplus` | `MESHIOPLUSPLUS_BUILD_C_API=ON` |
| `Fortran` | `meshioplusplus::meshioplusplus_fortran` | `MESHIOPLUSPLUS_BUILD_FORTRAN=ON` |

Asking for a component the install does not carry fails at `find_package` time rather than at the link, so a misconfigured prefix is caught immediately.

## Mesh backends in one prefix

`meshioplusplus::Mesh` is [selected at compile time](/cpp_backends), so a library built for one backend cannot be linked by a translation unit compiled for another — the two disagree about the layout of nearly every argument. Rather than force each consumer to build its own private install, one prefix ships **all three backends side by side**:

```cmake
target_link_libraries(my_solver PRIVATE meshioplusplus::core_kratos)   # Kratos-shaped ModelPart
target_link_libraries(my_tool   PRIVATE meshioplusplus::core_native)   # fastest pure-C++ storage
```

`meshioplusplus::core` is an alias for whichever backend the build treated as its default (`MESHIOPLUSPLUS_MESH_BACKEND`, `MESHIO` unless set). Trim the set to save build time — each backend is a full, independent compile of the core:

```sh
cmake -S . -B build -DMESHIOPLUSPLUS_INSTALL_CPP=ON \
  -DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS="MESHIO;KRATOS"
```

Each variant carries its own `MESHIOPLUSPLUS_MESH_BACKEND_*` macro in `INTERFACE_COMPILE_DEFINITIONS`, so linking `meshioplusplus::core_kratos` compiles *your* sources with `MESHIOPLUSPLUS_MESH_BACKEND_KRATOS` whether you ask for it or not. You cannot get this wrong through CMake.

Two guards catch the cases CMake cannot:

* defining more than one backend macro is a **compile** error (`mesh.hpp`);
* defining none — so `mesh.hpp` silently assumes `MESHIO` — while linking a `NATIVE`/`KRATOS` build is a **link** error naming the backend it expected (`detail/mesh_backend_check.hpp`). This is what protects a pkg-config or hand-written-makefile consumer, which gets the include path but not the definitions. Define `MESHIOPLUSPLUS_NO_BACKEND_LINK_CHECK` to opt out.

## The Kratos bridge

`kratos_bridge.hpp` is header-only and has **no Kratos dependency at all** — it is templated on the consumer's own model-part type via a `bridge_traits` customization point, defaulting to meshio++'s own `meshioplusplus::ModelPart`. It is unreachable through the C ABI (which cannot hand out a `ModelPart`), which is a large part of why this install exists.

```cpp
#include "meshioplusplus/kratos_bridge.hpp"

Kratos::ModelPart& dest = model.CreateModelPart("Imported");
meshioplusplus::to_model_part(source, dest, [&](std::size_t id) -> auto& {
    return dest.GetProperties(id);
});
```

See [mesh backends](/cpp_backends) for `KratosMesh` and the `ModelPart` materialization rules.

## Static, shared and symbol visibility

`BUILD_SHARED_LIBS` picks the library kind, as usual. The C++ libraries are built `-fvisibility=hidden` (`VISIBILITY_INLINES_HIDDEN` too) and the public surface is annotated with `MESHIOPLUSPLUS_API` (`export.hpp`), which also drives `__declspec(dllexport/dllimport)` on Windows — so a shared build exports its documented API and nothing else, on every platform.

`SOVERSION` is `0`: like the C API's, the C++ ABI is **declared unstable pre-1.0**. Pin an exact version if you ship binaries against it.

## Dependencies

Unlike the C API — whose HDF5/netCDF/zlib are private to its shared object and resolved by the dynamic linker — the C++ targets propagate their dependencies, because you compile the real headers and, in a static build, link the real dependency graph. The generated config issues the matching `find_dependency()` calls (`ZLIB`, `zstd`, `lz4`, `netCDF`, `KaHIP`, `OpenMP`/`TBB`/`Kokkos`) and installs `FindKaHIP.cmake` beside itself so a KaHIP build resolves for a consumer with no copy of its own.

HDF5 and MPI reach you as the imported targets `HDF5::HDF5` / `MPI::MPI_C`, so no absolute library path is baked into `meshioplusplusTargets.cmake` and the prefix can be moved to a machine whose HDF5 lives somewhere else. The price is that both Find modules need the **C language enabled**, which a C++-only project otherwise has no reason to do:

```cmake
project(my_solver LANGUAGES CXX C)   # or enable_language(C) before find_package
```

An HDF5-enabled install says exactly that, by name, if you forget.

The config also exports what the build actually compiled in, for consumers that want to branch:

```cmake
if(MESHIOPLUSPLUS_WITH_HDF5)  # ...also _NETCDF, _ZLIB, _ZSTD, _LZ4, _KAHIP, _EIGEN
  # XDMF/MED/CGNS are available
endif()
message(STATUS "backends: ${MESHIOPLUSPLUS_MESH_BACKENDS}")
```

## pkg-config

For build systems that do not use CMake:

```sh
g++ my_solver.cpp $(pkg-config --cflags --libs meshioplusplus-cxx) -o my_solver
```

`meshioplusplus-cxx.pc` is deliberately **separate** from the C API's `meshioplusplus.pc`, whose `Cflags` must stay valid for a plain-C compile and so can never carry `-std=c++20` or the backend macro. It covers the default backend only; an install with several exposes the others through CMake alone.

## Package managers

Both self-hosted packaging paths expose the option:

```sh
conan create . -o meshioplusplus/*:with_cxx_api=True \
               -o meshioplusplus/*:cxx_api_backends="MESHIO,KRATOS"

vcpkg install "meshioplusplus[cxx-api,cxx-api-kratos]" \
  --overlay-ports=packages/vcpkg
```

See the [C API's package-manager notes](/c_api#package-managers-conan-vcpkg) for the caveats that apply to both (the `libaec` overlay, and neither package being on ConanCenter or the vcpkg registry yet).

## Limitations

* The C++ **ABI is unstable** (`SOVERSION 0`). Rebuild consumers against a matching version.
* Building several backends multiplies compile time — each is a full, independent compile of the core. Trim `MESHIOPLUSPLUS_INSTALL_CPP_BACKENDS` to the ones you need.
* An HDF5-enabled install requires **`C` among your project's languages** (above).
* Names under `detail/` are installed (the public headers need them) but are **not** a stable API.
* `KratosMesh::InvalidateBlocks()`'s rebuild cannot recover a SubModelPart's *nesting* (regions are flat) or a region's `mDim`/`mTag`, which a SubModelPart has nowhere to store — see [mesh backends](/cpp_backends).

## Versioning: what to pin

The two installed components make **different** compatibility promises, and the
CMake package's `SameMajorVersion` mode describes only the first.

- **`COMPONENTS C`** (`libmeshioplusplus`, the flat C ABI) — `SOVERSION 0`, and
  the option structs (`mio_read_opts`, `mio_write_opts`, `mio_xdmf_series_opts`)
  grow only into their `reserved` tails, so a binary compiled against 9.0
  headers keeps running against a 9.1 `.so`. Pin the major:

  ```cmake
  find_package(meshioplusplus 9 CONFIG REQUIRED COMPONENTS C)
  ```

- **`COMPONENTS CXX`** — **no ABI promise at all.** `Mesh`, `ModelPart` and
  `GeometricalEntity` are header-defined types whose layout changes with the
  headers; v9.1.0 added a member to `GeometricalEntity`, for instance. The
  library and every consumer translation unit must be compiled from the *same*
  meshio++ version, and a consumer must be rebuilt whenever meshio++ is:

  ```cmake
  find_package(meshioplusplus 9.1 EXACT CONFIG REQUIRED COMPONENTS CXX)
  ```

  For distribution packaging that means an exact `= 9.1.x` dependency, not a
  range. The `SOVERSION 0` on the C++ variants exists so the files install
  cleanly, not as a compatibility claim.

The mesh-backend macro rides in `INTERFACE_COMPILE_DEFINITIONS`, so a CMake
consumer cannot disagree about the backend by accident; a non-CMake consumer
that defines none gets an undefined `mesh_backend_is_<backend>` at link time
instead (`MESHIOPLUSPLUS_NO_BACKEND_LINK_CHECK` opts out).

## `MESHIOPLUSPLUS_NO_STD_SPAN`

Boost's uBLAS and MSVC's `<span>` both use an internal macro named
`_BACKUP_ITERATOR_DEBUG_LEVEL`, so any MSVC translation unit including both
fails inside `<span>` itself ([boostorg/ublas#77](https://github.com/boostorg/ublas/issues/77)).
Kratos uses uBLAS, so this is a real collision. Three facts about the escape
hatch, all of which CI now gates:

- It guards **only** the `<span>` include and the inline `NativeMesh::ConnSpan()`
  accessor. No member variable is involved, so it is **ABI-neutral**.
- The guard is `#ifndef`, so **a consumer can define it themselves** —
  `target_compile_definitions(app PRIVATE MESHIOPLUSPLUS_NO_STD_SPAN)` — against
  *any* prefix, including a distro or Conan build that left the CMake option
  off. A private meshio++ build is not required.
- A build that sets the CMake option exports it (through
  `INTERFACE_COMPILE_DEFINITIONS` and `meshioplusplus-cxx.pc`), so consumers of
  such a prefix inherit the same choice by default.

## Parallelism: meshio++ is serial

There is **no MPI in the library**: no distributed reader or writer, no
communicator anywhere in the API, and none planned. `FindMPI` appearing in the
generated package config is HDF5's transitive requirement (a parallel HDF5 needs
`mpi.h` even for purely serial use of its API), not meshio++'s.

The intended distributed workflow is to decompose and then let each rank handle
its own piece: `partition(mesh, {nparts, ghost_layers})` produces exactly the
shared-node halo an MPI assembly needs when `mGhostLayers > 0`. All parallelism
inside meshio++ is intra-process (`MESHIOPLUSPLUS_PARALLEL_BACKEND`).
