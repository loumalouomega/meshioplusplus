# C API

The C++ core also ships as an installable shared library, `libmeshioplusplus`, with a stable pure-C99 header — the natural entry point for HPC codes written in C (and the foundation of the [Fortran interface](/fortran)). Like the [WebAssembly binding](/wasm), it is a flat, whole-mesh API over the same C++ core, built exclusively on the uniform mesh API, so it works identically under every [mesh backend](/cpp_backends).

## Building and installing

The C API is off by default (a plain `pip install` never builds it). From the repo:

```sh
build/configure.sh --c-api --build          # add --fortran for the Fortran module
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

(or pass `-DMESHIOPLUSPLUS_BUILD_C_API=ON` to a direct CMake configure). The install lays out:

```
include/meshioplusplus/meshioplusplus.h     # the only installed header
lib/libmeshioplusplus.so[.0.6.1.0]
lib/cmake/meshioplusplus/                   # find_package(meshioplusplus)
lib/pkgconfig/meshioplusplus.pc             # pkg-config
```

Consume it with pkg-config:

```sh
gcc my_solver.c $(pkg-config --cflags --libs meshioplusplus) -o my_solver
```

or CMake:

```cmake
find_package(meshioplusplus 6.1 REQUIRED)
target_link_libraries(my_solver PRIVATE meshioplusplus::meshioplusplus)
```

HDF5/netCDF/zlib are detected at configure time exactly as for the Python build; they are private dependencies of the shared library (consumers never link them directly).

## Package managers (Conan & vcpkg)

The C API also ships as a **Conan** recipe (root [`conanfile.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/conanfile.py)) and a **vcpkg** overlay port ([`ports/meshioplusplus/`](https://github.com/loumalouomega/meshioplusplus/tree/main/ports/meshioplusplus)). Both are self-hosted in the repo and drive the same `-DMESHIOPLUSPLUS_BUILD_C_API=ON -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF` install path, so consumers get the identical `meshioplusplus::meshioplusplus` target:

```sh
# Conan (options: with_hdf5 / with_netcdf / with_zlib / with_eigen / fortran)
conan create . -o meshioplusplus/*:with_hdf5=True -o meshioplusplus/*:with_netcdf=True

# vcpkg (features: hdf5 / netcdf / zlib -- all on by default)
vcpkg install meshioplusplus --overlay-ports=ports
```

Both are validated in CI on every PR and on `v*` release tags (`.github/workflows/packages.yml`). Two caveats: the shared library is **shared-only** (no static build yet), and the vendored **Eigen** submodule is off in both recipes (the MED transpose falls back to a hand-written loop), since it is absent from a source tarball. Neither is submitted to Conan Center / the upstream vcpkg registry today.

## Example

The complete, CI-tested example lives at [`doc/examples/c_api_example.c`](https://github.com/loumalouomega/meshioplusplus/blob/main/doc/examples/c_api_example.c):

```c
#include <meshioplusplus/meshioplusplus.h>

/* Build a mesh from raw arrays and write it (setters copy your buffers). */
mio_mesh* m = mio_mesh_create();
mio_mesh_set_points(m, MIO_FLOAT64, num_points, 3, xyz);            /* (n, 3) row-major */
mio_mesh_add_cell_block(m, "tetra", num_cells, 4, MIO_INT64, conn); /* 0-based indices */
mio_mesh_add_point_data(m, "temperature", MIO_FLOAT64, 1, (int64_t[]){num_points}, temp);
if (mio_write("out.vtu", m, NULL) != MIO_OK)          /* NULL: infer format from extension */
    fprintf(stderr, "%s\n", mio_last_error());
mio_mesh_free(m);

/* Read one back (getters are zero-copy borrows into mesh-owned memory). */
mio_mesh* r = mio_read("in.msh", NULL);               /* .msh defaults to gmsh */
const void* pts; mio_dtype dt;
mio_mesh_get_points(r, &pts, &dt);                    /* valid until mutation/free */
mio_mesh_free(r);

/* Or convert file-to-file without touching the data. */
mio_convert("in.msh", NULL, "out.vtk", NULL);
```

## The contract in five rules

1. **Errors**: every fallible function returns a `mio_status` (`MIO_OK == 0`) or `NULL`/`-1`; the message is retrievable via `mio_last_error()` (thread-local, valid until the next `mio_*` call on the same thread). No C++ exception ever crosses the ABI.
2. **Setters copy** — your buffers can be freed as soon as the call returns.
3. **Getters are zero-copy** — returned data pointers alias mesh-owned memory and stay valid until the next *mutating* call on that mesh or `mio_mesh_free()`. Read-only accessors never invalidate them.
4. **Arrays are row-major** (C order): points `(num_points, dim)`, connectivity `(num_cells, nodes_per_cell)` with **0-based** node indices.
5. **String getters** use the snprintf convention: copy at most `buflen - 1` bytes plus a NUL, return the full length (excluding the NUL), `-1` on error.

## API reference

| Group | Functions |
| --- | --- |
| Introspection | `mio_version`, `mio_mesh_backend`, `mio_format_readable`, `mio_format_writable`, `mio_last_error` |
| Cell-type metadata | `mio_cell_type_name`, `mio_cell_type_from_name`, `mio_cell_type_num_nodes`, `mio_cell_type_dimension` (the `mio_cell_type` enum mirrors the C++ table; strings like `"tetra10"` are the primary representation) |
| Lifecycle & I/O | `mio_mesh_create`, `mio_mesh_free`, `mio_read`, `mio_write`, `mio_convert` (format `NULL`/`""` = infer from extension; `.msh` → gmsh, `.inp` → abaqus) |
| Building | `mio_mesh_set_points`, `mio_mesh_add_cell_block` (int32 connectivity is widened to the core's int64), `mio_mesh_add_point_data`, `mio_mesh_append_cell_data` (one call per cell block, in block order), `mio_mesh_add_field_data` |
| Points/cells | `mio_mesh_num_points`, `mio_mesh_point_dim`, `mio_mesh_get_points`, `mio_mesh_num_cell_blocks`, `mio_mesh_cell_block_info`, `mio_mesh_cell_block_type`, `mio_mesh_cell_block_conn` |
| Named data | `mio_mesh_num_{point,cell,field}_data`, `mio_mesh_{point,cell,field}_data_name` (names in sorted order — identical on every backend), `mio_mesh_get_{point,cell,field}_data`, `mio_mesh_cell_data_num_blocks` |

Every function is documented in the installed header, [`bindings_c/include/meshioplusplus/meshioplusplus.h`](https://github.com/loumalouomega/meshioplusplus/blob/main/bindings_c/include/meshioplusplus/meshioplusplus.h).

## Format support

All formats with a C++ implementation are available — the same set as the [WASM binding](/wasm#format-support) **plus**, when the build found the dependencies, the HDF5-backed formats (`cgns`, `h5m`, `hmf`, `med`, XDMF's HDF heavy-data path) and the netCDF-backed `exodus`. Probe at runtime with `mio_format_readable()`/`mio_format_writable()`; requesting a compiled-out format fails with a message naming the missing dependency. Formats that only exist in Python (`mdpa`, `svg`, `neuroglancer`, …) are not reachable from C.

Parameterized writers use each format's Python-reference default (VTU: binary+zlib, gmsh: 4.1 binary, STL: ASCII, XDMF: HDF when built with HDF5 else XML, …); per-call writer options are a possible future addition.

## Limitations (v1)

- **Ragged cell blocks** (polygons/polyhedra of varying size) cannot be built through the C API, and on meshes read from files their connectivity is not accessible (`mio_mesh_cell_block_conn` returns `MIO_ERR_UNSUPPORTED`; counts, type and `is_ragged` still work).
- **Side-channel metadata** is dropped: `point_sets`/`cell_sets` (`ansysinp`, `unv`), MED families/groups, OpenFOAM cell tags. Use the Python API when you need those.
- A `mio_mesh` handle is not thread-safe; distinct handles may be used from distinct threads freely.
