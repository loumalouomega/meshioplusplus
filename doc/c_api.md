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
find_package(meshioplusplus REQUIRED)  # pin a minimum version if you need one, e.g. `meshioplusplus 7.0`
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

`ports/` also carries a second overlay, `ports/libaec/`, alongside `ports/meshioplusplus/` — it replaces vcpkg registry's own `libaec` port (a transitive dependency of `hdf5`'s `szip` feature) with an otherwise-identical one that sources the same release from `libaec`'s official GitHub mirror instead of its `gitlab.dkrz.de` host, which had been intermittently rate-limiting this project's release CI. `--overlay-ports=ports` picks up both automatically.

Both are validated in CI on every PR and on `v*` release tags (`.github/workflows/packages.yml`). Two caveats: the shared library is **shared-only** (no static build yet), and the vendored **Eigen** submodule is off in both recipes (the MED transpose falls back to a hand-written loop), since it is absent from a source tarball. Neither is submitted to Conan Center / the upstream vcpkg registry today, so nothing resolves `meshioplusplus` as a plain requirement out of the box; as a stopgap, every `v*` release attaches ready-to-use Linux artifacts, for both x86_64 and arm64, to its [GitHub Release](https://github.com/loumalouomega/meshioplusplus/releases) that supply the missing recipe:

```sh
# Conan: restore the release archive (adds the meshioplusplus recipe + a matching
# prebuilt binary to your local cache), then install normally -- Conan Center
# still supplies the transitive hdf5/netcdf/zlib/openssl binaries over the
# network as usual; --build=missing rebuilds meshioplusplus itself from the
# archive's recipe if your profile doesn't match the prebuilt one.
conan cache restore meshioplusplus-conan-full-linux-x86_64.tgz   # or ...-linux-arm64.tgz
conan install --requires=meshioplusplus/<version> --build=missing \
  -o meshioplusplus/*:with_hdf5=True -o meshioplusplus/*:with_netcdf=True -o meshioplusplus/*:with_zlib=True

# vcpkg, from the release archive (unzip first):
vcpkg install meshioplusplus --overlay-ports=meshioplusplus-vcpkg-overlay-port-<version>
```

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
| Lifecycle & I/O | `mio_mesh_create`, `mio_mesh_free`, `mio_read`, `mio_write`, `mio_convert` (format `NULL`/`""` = infer from extension; `.msh` → gmsh, `.inp` → abaqus, `.node`/`.ele` → tetgen — pass `"triangle"` for 2D Triangle pairs) |
| Building | `mio_mesh_set_points`, `mio_mesh_add_cell_block` (int32 connectivity is widened to the core's int64), `mio_mesh_add_point_data`, `mio_mesh_append_cell_data` (one call per cell block, in block order), `mio_mesh_add_field_data` |
| Points/cells | `mio_mesh_num_points`, `mio_mesh_point_dim`, `mio_mesh_get_points`, `mio_mesh_num_cell_blocks`, `mio_mesh_cell_block_info`, `mio_mesh_cell_block_type`, `mio_mesh_cell_block_conn` |
| Named data | `mio_mesh_num_{point,cell,field}_data`, `mio_mesh_{point,cell,field}_data_name` (names in sorted order — identical on every backend), `mio_mesh_get_{point,cell,field}_data`, `mio_mesh_cell_data_num_blocks` |
| Operations | `mio_extract_surface`, `mio_extract_skin`, `mio_attach_quality`, `mio_quality_counts`, `mio_sniff_format`, `mio_compute_bandwidth`, and reordering: `mio_reorder` (method `"rcm"`/`"morton"`/`"hilbert"`) → `mio_reorder_result` (`mio_reorder_result_mesh` borrow, `mio_reorder_result_node_perm` / `mio_reorder_result_num_cell_perms` / `mio_reorder_result_cell_perm` zero-copy borrows, `mio_reorder_result_take_mesh` ownership transfer, `mio_reorder_result_free`), and comparison: `mio_meshes_equal(a, b, atol, rtol, unordered, &equal)` plus `mio_diff(a, b, atol, rtol, unordered, &result)` → `mio_diff_result` (`mio_diff_result_verdict`, `mio_diff_result_point_summary`, `mio_diff_result_num_block_diffs` / `mio_diff_result_block`, `mio_diff_result_free`), and merging: `mio_merge(meshes, count, weld, atol, source_tag, data_policy, drop_duplicate_cells)` (array of mesh handles → a new `mio_mesh*`; `data_policy` `0`=intersection / `1`=fill), and the editing/stats bundle: `mio_transform(mesh, matrix[16], rotate_vector_data)` → `mio_mesh*`; `mio_clean(mesh, weld, atol, remove_orphans, drop_degenerate, drop_duplicate_cells, &n_welded, &n_orphan, &n_degen, &n_dup)` → `mio_mesh*`; `mio_crop_bbox(mesh, lo, hi, mode, record_ids)` / `mio_crop_plane(mesh, point, normal, mode, record_ids)` (`mode` `0`=all / `1`=any) → `mio_mesh*`; `mio_split(mesh, by, tag_name)` → `mio_split_result` (`mio_split_result_count` / `_key` / `_mesh` borrow / `_take_mesh` / `_free`); `mio_stats(mesh, &report)` → fills a `mio_stats_report` struct (scalar geometric measures; per-cell-type counts are not carried) |

Every function is documented in the installed header, [`bindings_c/include/meshioplusplus/meshioplusplus.h`](https://github.com/loumalouomega/meshioplusplus/blob/main/bindings_c/include/meshioplusplus/meshioplusplus.h).

## Format support

All formats with a C++ implementation are available — the same set as the [WASM binding](/wasm#format-support) **plus**, when the build found the dependencies, the HDF5-backed formats (`cgns`, `h5m`, `hmf`, `med`, XDMF's HDF heavy-data path) and the netCDF-backed `exodus`. This includes the write-only visualization formats `svg` and `tikz` (writable, not readable; emitted with the fixed default styling — 3D meshes render their projected boundary skin with the default isometric camera). Probe at runtime with `mio_format_readable()`/`mio_format_writable()`; requesting a compiled-out format fails with a message naming the missing dependency. Formats that only exist in Python (`mdpa`, `neuroglancer`, …) are not reachable from C.

Parameterized writers use each format's Python-reference default (VTU: binary+zlib, gmsh: 4.1 binary, STL: ASCII, XDMF: HDF when built with HDF5 else XML, …). In particular `stl`/`ply` extract and write the boundary **skin** of a 3D volume mesh, matching the Python default (see [Skin extraction](/extract_skin)); a standalone skin-extraction entry point is not part of the C API yet (documented follow-up). Per-call writer options are a possible future addition.

## Native command-line binary

The core also ships a **Python-free CLI** — the same verbs as the Python command-line tool, built directly on the C++ core, so it needs neither a Python interpreter nor the pybind11 extension.

**Prebuilt binaries:** every [GitHub Release](https://github.com/loumalouomega/meshioplusplus/releases) tag (`v*`) ships a ready-to-run, statically-linked `meshioplusplus` executable for Linux (x86_64 and arm64), macOS (universal x86_64+arm64), and Windows (x86_64) — no install step, no runtime dependency to satisfy (no `libstdc++`/vcredist). These release binaries build with `MESHIOPLUSPLUS_WITH_HDF5/NETCDF/ZLIB=OFF` and the sequential parallel backend, trading those for a single dependency-free file (see `.github/workflows/cli.yml`); build from source (below) for the full HDF5/netCDF/multi-threaded feature set.

**Build from source:** off by default; build it alongside a standalone C++ tree:

```sh
build/configure.sh --cli --build          # -DMESHIOPLUSPLUS_BUILD_CLI=ON
build/cpp-release/meshioplusplus --help
```

Pass `-DMESHIOPLUSPLUS_STATIC_RUNTIME=ON` (alongside `--cli`, via a direct CMake configure) for the same statically-linked-runtime behavior as the release binaries: static `libgcc`/`libstdc++` on GNU, static MSVC CRT (`/MT`) on Windows.

The installed executable is named `meshioplusplus` and mirrors the Python verbs:

```sh
meshioplusplus convert in.msh out.vtu
meshioplusplus info out.vtu
meshioplusplus ascii out.vtu             # in-place; also: binary / compress / decompress
meshioplusplus quality mesh.vtu
meshioplusplus extract-surface vol.vtu surf.vtu
meshioplusplus reorder in.vtu out.vtu --method rcm --report
meshioplusplus diff a.vtu b.vtu          # nonzero exit if different
meshioplusplus merge a.vtu b.vtu out.vtu
```

Format dispatch reuses the shared registry (with a content-sniff fallback on read); ASCII/binary/compress variants call the per-format writers directly. Because it links only the C++ core it inherits the same flat-surface limitations as this C API: point/cell **sets** and the `convert -s/-d` sets↔data conversions are unavailable (they live only in the Python `Mesh`), so `info`/`diff` omit sets. There is also no Python fallback, so a file whose C++ reader raises (e.g. a Gmsh file using `$Entities`) is not readable by the native CLI — use the Python CLI for those. It works under any [mesh backend](/cpp_backends) and any optional-dependency configuration; verbs that touch a compiled-out HDF5/netCDF format report the missing dependency by name.

## Limitations (v1)

- **Ragged cell blocks** (polygons/polyhedra of varying size) cannot be built through the C API, and on meshes read from files their connectivity is not accessible (`mio_mesh_cell_block_conn` returns `MIO_ERR_UNSUPPORTED`; counts, type and `is_ragged` still work).
- **Side-channel metadata** is dropped: `point_sets`/`cell_sets` (`ansysinp`, `unv`), MED families/groups, OpenFOAM cell tags. Use the Python API when you need those.
- A `mio_mesh` handle is not thread-safe; distinct handles may be used from distinct threads freely.
