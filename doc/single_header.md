# Single-header C++ (header-only)

The whole meshio++ C++ core is also available as **one self-contained header**,
[`single_include/meshioplusplus/meshioplusplus.hpp`](https://github.com/loumalouomega/meshioplusplus/blob/main/single_include/meshioplusplus/meshioplusplus.hpp).
Drop it into any project — no CMake, no submodules, nothing to link (unless you
opt into the optional formats below). It follows the well-known
[STB](https://github.com/nothings/stb) single-file-library convention:

```cpp
// in exactly ONE .cpp of your project — pulls in all the implementations:
#define MESHIOPLUSPLUS_IMPLEMENTATION
#include "meshioplusplus/meshioplusplus.hpp"

// in every other translation unit — declarations only:
#include "meshioplusplus/meshioplusplus.hpp"
```

The header bundles [pugixml](https://pugixml.org/) (MIT), so the XML-based
formats (VTU, XDMF-XML, DOLFIN) work out of the box.

## Example

```cpp
#define MESHIOPLUSPLUS_IMPLEMENTATION
#include "meshioplusplus/meshioplusplus.hpp"
#include <cstring>
using namespace meshioplusplus;

int main() {
    Mesh m;
    NDArray points(DType::Float64, {3, 3});
    double xyz[9] = {0,0,0, 1,0,0, 0,1,0};
    std::memcpy(points.Data(), xyz, sizeof(xyz));
    m.AssignPoints(std::move(points));

    NDArray conn(DType::Int64, {1, 3});
    long long tri[3] = {0, 1, 2};
    std::memcpy(conn.Data(), tri, sizeof(tri));
    m.AddCellBlock("triangle", std::move(conn));

    registry_writers().at("stl")("mesh.stl", m);       // write
    Mesh back = registry_readers().at("stl")("mesh.stl");   // read back
    return 0;
}
```

```sh
g++ -std=c++20 -I single_include main.cpp -o main
```

The public surface is the same uniform mesh API and format
[registry](https://github.com/loumalouomega/meshioplusplus/blob/main/cpp/include/meshioplusplus/registry.hpp)
the rest of the C++ core uses: `registry_readers()` / `registry_writers()`
(`name -> function`), `resolve_format(path, "")`, plus the `Mesh` type and the
`ReadError` / `WriteError` exceptions.

## Configuration macros

Everything defaults to a zero-dependency build. Define these **before** including
the header to change that:

| Macro | Effect |
| ----- | ------ |
| *(none)* | Mesh backend `MESHIO`, sequential parallel backend, no external libraries. |
| `MESHIOPLUSPLUS_MESH_BACKEND_NATIVE` / `_KRATOS` | Select a different [mesh backend](cpp_backends.md). |
| `MESHIOPLUSPLUS_PARALLEL_OPENMP` / `_TBB` / `_STL` | Enable a parallel backend (and link/compile with the matching flags). |
| `MESHIOPLUSPLUS_HAS_ZLIB` | VTU zlib compression — also link `-lz`. |
| `MESHIOPLUSPLUS_HAS_HDF5` | CGNS / HMF / H5M / MED / XDMF-HDF — also link `-lhdf5`. |
| `MESHIOPLUSPLUS_HAS_NETCDF` | Exodus — also link `-lnetcdf`. |
| `MESHIOPLUSPLUS_HAS_EIGEN` | Faster MED transpose — add Eigen to the include path. |

Optional-dependency code stays behind its `MESHIOPLUSPLUS_HAS_*` guard, so the
default include compiles with no third-party libraries at all — the same format
set as the [WebAssembly build](wasm.md).

## How it is generated

The single header is **generated** from the sources under `cpp/` by
[`tools/amalgamate.sh`](https://github.com/loumalouomega/meshioplusplus/blob/main/tools/amalgamate.sh)
(which drives `tools/amalgamate/amalgamate.py`) and **committed** to the repo.
CI regenerates it on every push and fails if the committed copy is stale, then
smoke-compiles it (declarations-only, `MESHIOPLUSPLUS_IMPLEMENTATION`, and a
two-TU link). **Do not edit the generated file by hand** — edit the sources
under `cpp/` and run:

```sh
./tools/amalgamate.sh          # regenerate single_include/…/meshioplusplus.hpp
./tools/amalgamate.sh --smoke  # regenerate + smoke-compile
```

The generator emits every header once, at file scope, in dependency order (so a
header can never be trapped inside another's `#ifdef`), then the `.cpp` bodies
under a single `MESHIOPLUSPLUS_IMPLEMENTATION` guard.
