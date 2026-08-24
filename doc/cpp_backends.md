# C++ mesh backends

The C++ core has three interchangeable **in-memory mesh backends**, selected at build time with the CMake option `MESHIOPLUSPLUS_MESH_BACKEND` (exactly one is compiled per build, like the [parallel backend](installation.md#parallelism)):

| Backend  | Structure | Use it for |
| -------- | --------- | ---------- |
| `MESHIO` (default) | `Mesh`/`CellBlock` over dtype-erased `NDArray`s, mirroring the Python `meshio.Mesh` | The Python extension (**required** for it — PyPI wheels always use MESHIO) |
| `NATIVE` | Canonical statically-typed storage: Float64 points, Int64 connectivity, `CellType` enum, CSR ragged blocks | The fastest pure-C++ consumer surface; the [WebAssembly build](wasm.md) uses it |
| `KRATOS` | A [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos)-style `ModelPart` (Nodes / Elements / Conditions / SubModelParts) | Near-costless exchange with Kratos (or CoSimIO) via the header-only bridge |

```sh
# standalone C++ build with a non-default backend (implies no Python extension)
./build/configure.sh --mesh-backend NATIVE --tests --build
./build/configure.sh --mesh-backend KRATOS --tests --build
```

Every format reader/writer is written against a **uniform mesh API** (`src/cpp/include/meshioplusplus/mesh_api.hpp`), so all ~42 formats compile and round-trip identically under every backend — the full GoogleTest suite runs per backend in CI. Selecting `NATIVE`/`KRATOS` together with `MESHIOPLUSPLUS_BUILD_PYTHON=ON` is a CMake configure error: the zero-copy numpy boundary is written against MESHIO's exact struct layout.

## The uniform mesh API

`meshioplusplus::Mesh` is a compile-time alias (`src/cpp/include/meshioplusplus/mesh.hpp`) for the selected backend type. All backends implement:

**Ingestion** (what readers use — `NDArray` is the universal staging type, handed over by move):

```cpp
mesh.AssignPoints(NDArray points);                    // (n, dim), float dtype
mesh.AddCellBlock("tetra", NDArray conn);             // (n, npc), integer dtype
mesh.AddPolygonBlock("polygon", rows);                // 1-level ragged
mesh.AddPolyhedronBlock("polyhedron", cells);         // 2-level ragged
mesh.AddPointData("temperature", NDArray data);
mesh.AddCellData("gmsh:physical", std::vector<NDArray> perBlock);
mesh.AppendCellData("medit:ref", NDArray oneBlock);   // incremental variant
mesh.AddFieldData("group", NDArray data);
```

**Accessors** (what writers use):

```cpp
mesh.NumPoints();  mesh.PointDim();  mesh.Points();   // const NDArray&
mesh.NumCellBlocks();
for (const auto cb : mesh.CellRange()) {              // Mesh::CellView values
    cb.Type();          // "tetra", ...
    cb.NumCells();  cb.NodesPerCell();  cb.Conn();    // const NDArray&
    cb.IsRagged();  cb.Row(i);  cb.Face(i, f);        // ragged access
}
mesh.PointDataNames();                                // always sorted
mesh.PointData("temperature");                        // const NDArray&
mesh.CellData("gmsh:physical", blockIndex);
```

Dtype rules: MESHIO stores arrays exactly as received; NATIVE and KRATOS canonicalize **within kind** (floats → Float64, ints → Int64 — an integer tag array never becomes float, so "first integer cell_data is the tag" format conventions survive). Owning arrays that are already canonical are *moved*, not copied, and readers produce canonical dtypes almost everywhere, so ingest is near-free. One observable consequence: under NATIVE/KRATOS, a file with Float32 points is re-written as Float64.

## The NATIVE backend

`meshioplusplus::NativeMesh` (`backends/native_mesh.hpp`) adds a fast-consumer surface on top of the uniform API:

```cpp
const double* xyz = mesh.PointsData();                    // contiguous Float64
std::span<const std::int64_t> conn = mesh.ConnSpan(0);    // per-block Int64
meshioplusplus::CellType t = mesh.BlockType(0);           // enum, not string
const auto& csr = mesh.GlobalConnectivity();              // whole-mesh CSR
// csr.mOffsets (ncells+1), csr.mConn (flat), csr.mTypes (one per cell)
```

Ragged blocks are stored CSR-style (flat node buffer + offset arrays) rather than nested vectors. `GlobalConnectivity()` is built lazily and cached.

## The KRATOS backend

`meshioplusplus::KratosMesh` (`backends/kratos_mesh.hpp`) puts a Kratos-style `ModelPart` behind the same API:

```cpp
meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh("part.msh");
meshioplusplus::ModelPart& mp = mesh.GetModelPart();   // materialized lazily
mp.NumberOfNodes();       // Ids are 1-based (node Id = point index + 1)
mp.NumberOfElements();    // cell blocks of the mesh's max topological dim
mp.NumberOfConditions();  // lower-dimension blocks
mp.GetSubModelPart("gmsh_physical_1");  // built from integer tag arrays
```

- **Elements vs Conditions**: blocks whose topological dimension equals the mesh's maximum become Elements, lower-dimension blocks Conditions (the Kratos convention, matching the [mdpa](formats/mdpa.md) reader/writer), each kind Id-numbered 1..N in block order with default Kratos names (`Element3D4N`, `SurfaceCondition3D3N`, ... — `backends/kratos_names.hpp`).
- **Tags → SubModelParts**: integer cell-data under well-known names (`gmsh:physical`, `su2:tag`, `medit:ref`, `cell_tags`, ...) automatically become SubModelParts named `<key>_<value>` containing the tagged entities and their nodes. Disable with `mesh.SetBuildSubModelPartsFromTags(false)` before the first `GetModelPart()` call. The tag arrays also stay as elemental/conditional data, so writer round-trips are byte-identical.

Since v9.2.0 the pass can also be narrowed **per key**, which matters for `.mdpa`: a Kratos properties id is read as a `gmsh:physical` cell tag, so a deck with one sub model part comes back with `gmsh_physical_0` beside it — material assignment surfacing as a group. For a genuine gmsh file that same inference is wanted, so the key stays in the table and the consumer says which meaning applies:

  ```cpp
  mesh.ExcludeTagSubModelPartKey("gmsh:physical");   // everything else still applies
  mesh.SetTagSubModelPartKeys({"cell_tags"});        // ...or an explicit allow-list
  mesh.TagSubModelPartKeys();                        // the effective set
  ```

An empty allow-list restores the default. Both setters re-materialize, so they may be called after `GetModelPart()`.

- **Property sets → `Properties`** (v9.2.0): `Begin Properties` bodies read from a `.mdpa` ride on the mesh (`AddPropertySet`/`GetPropertySet`) and materialize as `ModelPart::Properties` blocks with **real values**, not just ids — which is what makes `to_model_part`'s "apply property" overload transfer material data. They are keyed by id rather than by entity index, so no operation ever has to remap them; shape-preserving operations (`clean`, `smooth`, `transform`, `attach_quality`, the data ops) carry them through, while restructuring and multi-input ones (`merge`, `crop`, `split`, `partition`) do not.
- **point/cell data** become simplified per-entity variables (`mp.GetNodalData("temperature")`, `mp.GetElementalData(...)`, `mp.GetNodalValue("temperature", nodeId)`).
- **Lazy and write-transparent**: a plain read → write conversion never builds the ModelPart at all; writer accessors serve from the canonical staging storage, so output matches the NATIVE backend byte-for-byte.
- **Mutation**: after changing the ModelPart directly (`CreateNewNode`, ...), call `mesh.InvalidateBlocks()`; the block view is then rebuilt from the ModelPart (consecutive same-type Elements group into blocks, then Conditions). Ragged pass-through blocks (polygon/polyhedron — Kratos has no such geometry) and SubModelPart structure are not representable back and are dropped by that rebuild.

### The Kratos bridge (works from any backend)

`src/cpp/include/meshioplusplus/kratos_bridge.hpp` is header-only, templated, and independent of the selected mesh backend — `meshioplusplus::ModelPart` and the bridge compile in every build. `to_model_part` populates **any** Kratos-like class through the narrow creation API only (`CreateNewNode/CreateNewElement(name, id, nodeIds, properties)/...`), so it works with a real `Kratos::ModelPart` without meshio++ ever linking Kratos:

```cpp
#include "meshioplusplus/kratos_bridge.hpp"

// Real Kratos: map properties ids to Properties::Pointer.
Kratos::ModelPart& dest = model.CreateModelPart("FromMeshio");
meshioplusplus::to_model_part(mesh.GetModelPart(), dest, [&](auto pid) {
    return dest.HasProperties(pid) ? dest.pGetProperties(pid)
                                   : dest.CreateNewProperties(pid);
});

// And back (duck-typed via bridge_traits; specialize it for classes whose
// accessors differ from meshioplusplus::ModelPart's shape):
meshioplusplus::ModelPart mine = meshioplusplus::from_model_part(source);
```

Sub model parts (including nested ones) are copied when the destination supports `CreateSubModelPart`/`AddNodes`/`AddElements`/`AddConditions`. CoSimIO's `ModelPart` (whose `CreateNewElement` takes an `ElementType` enum) is populated with a thin loop instead — CI compile-checks that pattern against the real CoSimIO headers on every run. The conversion cost is one O(n) bulk-create pass — the same cost Kratos's own CoSimIO conversion utilities pay, because Kratos's pointer-based entity storage cannot be aliased from outside.

## Benchmarks between backends

`benchmark/bench_backends.sh` builds one benchmark binary per backend (`src/cpp/benchmark/bench_backends.cpp`, CMake option `MESHIOPLUSPLUS_BUILD_BENCHMARKS=ON`) and collates a CSV comparing ingest, accessor traversal, ModelPart materialization (KRATOS only), and full file round-trips on a synthetic tet cube. See [Benchmarks](benchmarks.md#mesh-backend-benchmarks) for results and method.

## Adding a backend

One CMake branch defining `MESHIOPLUSPLUS_MESH_BACKEND_<NAME>`, one `#elif` in `src/cpp/include/meshioplusplus/mesh.hpp`, and a `backends/<name>_mesh.hpp` implementing the uniform API (`mesh_api.hpp` documents the exact contract; `tests/cpp/test_mesh_api.cpp` is its executable form and must pass).

## Kratos entity names, properties and nesting (v9.1.0)

Three things a Kratos consumer needs that used to be lost on the way through.

**Entity names.** `GeometricalEntity` carries an optional Kratos registration name (`Name()`/`HasName()`); empty means "derive it from the cell type", which is what every caller got before. It matters because deriving is lossy in one direction only: `SmallDisplacementElement3D4N` resolves to `Tetrahedra4`, but `Tetrahedra4` only ever derives back to the canonical `Element3D4N`.

Names are stored as interned `const std::string*` from a root-owned `detail::NamePool`, not owned strings: there is one distinct name per *block* but one entity per *cell*, so an owned string would add ~320 MB to a 10 M-element model part. `std::unordered_set` is node-based, so the pointers survive both rehashing and a `ModelPart` move.

On the mesh backend they ride per block:

```cpp
meshioplusplus::MdpaInfo info;
Mesh m = meshioplusplus::read_mdpa(path, info);
for (std::size_t b = 0; b < m.NumCellBlocks(); ++b)
    m.SetBlockEntityName(b, info.mEntityNames[b].mName);
ModelPart& r_mp = m.GetModelPart();   // entities now carry the real names
```

`SetBlockEntityName`/`BlockEntityName` are KRATOS-only extras, not part of the uniform mesh API — no other backend has a ModelPart to spell names for. This is the same "fast-consumer surface" shape as NATIVE's `PointsData()`/`ConnSpan()`.

**Properties.** `ModelPart` has a real Properties store (`CreateNewProperties`/`GetProperties`/`Properties()`), and `Materialize()` reads each entity's properties id from the `gmsh:physical` cell_data instead of leaving every entity on properties 0. The values are `PropertyValue` key/value pairs (`meshioplusplus/properties.hpp`), shared with MDPA so there is one representation. Turning them into typed Kratos `Variable<T>`s is the consumer's job — see the applier overload in [`doc/cpp_api.md`](cpp_api.md).

**Nesting.** SubModelPart nesting now survives the staging round trip, flattened into region names with `/` — see [`doc/regions.md`](regions.md#nested-groups-the--convention).
