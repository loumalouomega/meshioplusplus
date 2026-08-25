# Vendored gidpost 2.14

This directory is a hardcopy of a subset of [gidpost 2.14](https://www.gidsimulation.com/downloads/gidpost-2-14-library-to-write-postprocess-results-for-gid-in-ascii-binary-or-hdf5-format/),
CIMNE's C library for writing [GiD](https://www.gidsimulation.com/) postprocess
files in ASCII, compressed binary and HDF5. It backs meshio++'s `gid` writer
(`src/cpp/src/formats/gid.cpp`).

- **Upstream version**: 2.14 (`GP_VERSION_MAJOR`/`GP_VERSION_MINOR` in `source/gidpost.h`).
- **Licence**: BSD-2-Clause-Views, Copyright (c) 2015-2024 CIMNE. See `LICENSE`
  in this directory (copied verbatim from upstream) and the notice requirements
  it carries — permissive, no obligation on the rest of this MIT repository
  beyond keeping the notice with this code.
- **Vendored near-unmodified**: line endings were normalized (CRLF -> LF), and
  `source/gidpostHDF5.c` carries **one two-part upstream bug fix** (below).
  Everything else meshio++ needs (`NDEBUG`, `ENABLE_HDF5`) is supplied as a
  compile definition from `CMakeLists.txt`, never a source edit.

## The one source fix, and why (re-apply it on any version bump)

Both parts are in `new_CurrentHdf5WriteData()` / `delete_CurrentHdf5WriteData()`
in `source/gidpostHDF5.c`, and both are marked inline with `meshio++ fix`
comments. Together they made gidpost's HDF5 flavour **single-use per process,
and deadlock-prone afterwards**:

1. `G_num_HDF5_files` was incremented on open but **never decremented** on
   close, so it counted files *ever* opened rather than files *currently* open.
   The second HDF5 file opened in a process -- even strictly after the first had
   been closed -- therefore tripped the "more than one HDF5 file" thread-safety
   check and was refused. It is now a file-scope static that
   `delete_CurrentHdf5WriteData()` decrements.
2. That refusal path did `return NULL` **without `_UNLOCK_`**, leaving gidpost's
   process-global `pthread_mutex_t` held forever, so every subsequent gidpost
   call in the process deadlocked -- including ones for the unrelated ASCII and
   binary flavours.

Found by meshio++'s own test suite, which writes several HDF5 files in one
process; a single-file program never reaches either bug. Report upstream if
gidpost gains an issue tracker.

## What was excluded, and why

| Excluded | Reason |
|---|---|
| `cfortran/`, `source/gidpostforAPI.c` | gidpost's own licence names `cfortran.h` as third-party code under its own distribution policy; `gidpostforAPI.c` is the only file that includes it. |
| `source/gidpostfor.c`, `source/gidpostfor.h` | gidpost's Fortran binding. meshio++ has its own Fortran module over its own C API and needs none of this layer. |
| `binaries/` (~96 MB prebuilt libraries), `doc/` (the ~7 MB `gidpost.pdf` reference manual), `examples/`, `gidpost-swig/`, `win/`, `fortran_module/` | Not vendored assets — reference material and build artefacts, not source this library compiles. |
| Upstream `CMakeLists.txt`, `source/CMakeLists.txt`, `gidpost_config.h.in`, `set_nv_compilers.sh`, upstream `README.md` | meshio++ supplies its own build integration (`CMakeLists.txt`'s `MESHIOPLUSPLUS_WITH_GIDPOST` option); these files are not used. |

## Layout

The `source/` subdirectory is kept flat and unmodified relative to upstream so
its own `#include "gidpostInt.h"`-style includes resolve without change.

- 9 `.c` files: `gidpost`, `gidpostFILES`, `gidpostHDF5`, `gidpostHash`,
  `gidpostInt`, `hashtab`, `hdf5c`, `lookupa`, `recycle`.
- 14 headers: the above plus `gidpostMutex.h`, `gidpost_types.h`,
  `gidpost_functions.h`, `gidpost_functions_deprecated.h`, `standard.h`.

`gidpostHDF5.c` and `hdf5c.c` are only compiled when
`MESHIOPLUSPLUS_WITH_HDF5=ON` (see `CMakeLists.txt`) — `gidpostHDF5.c` does not
self-gate on `ENABLE_HDF5` (its own `#ifdef ENABLE_HDF5` is commented out), so
the CMake source list conditionally excludes it rather than relying on the
macro alone.

## Building against it

The vendored sources are compiled directly into meshio++'s core source list
(see `MESHIOPLUSPLUS_WITH_GIDPOST` in the root `CMakeLists.txt`) rather than a
separate library target, so they inherit the same compiler flags the C++ core
needs (in particular Emscripten's `-sUSE_ZLIB=1`, required for `<zlib.h>` to
resolve under WASM). `gidpostInt.h` includes `<zlib.h>` unconditionally, so
this library cannot be built without zlib.
