//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file mesh_api.hpp
 * @brief The uniform format-facing mesh API: the compile-time contract every
 * mesh backend implements, plus `mesh_backend_name()`.
 *
 * meshio++ has three interchangeable in-memory mesh backends, selected at
 * build time by `MESHIOPLUSPLUS_MESH_BACKEND` (exactly one is compiled, like
 * the parallel backend — see `mesh.hpp` for the dispatch):
 *
 *  - **MESHIO** (default; required when the pybind11 extension is built):
 *    the meshio-mirroring `Mesh`/`CellBlock` over dtype-erased `NDArray`s
 *    (`backends/meshio_mesh.hpp`).
 *  - **NATIVE**: canonical statically-typed storage — Float64 points, Int64
 *    connectivity, `CellType` enum, CSR-shaped ragged blocks
 *    (`backends/native_mesh.hpp`). The fastest pure-C++/WASM consumer.
 *  - **KRATOS**: a Kratos-Multiphysics-style `ModelPart`
 *    (Nodes/Elements/Conditions/SubModelParts) behind the same API
 *    (`backends/kratos_mesh.hpp`).
 *
 * Format readers/writers (and `bindings/wasm/`) MUST use only the methods
 * below — never backend-specific members — so every format compiles
 * unchanged under all backends. (`bindings/python/np_conversions.hpp` is the
 * one sanctioned exception: the Python build is pinned to MESHIO.)
 *
 * ## The contract (duck-typed; each backend implements these members)
 *
 * Reader-side ingestion — `NDArray` is the universal staging type; readers
 * build owning arrays locally (keeping the `NDArray::Uninit` fill hot loops)
 * and hand them over **by move**. MESHIO stores arrays as received; NATIVE
 * and KRATOS canonicalize *within kind* (floats → Float64, ints → Int64 —
 * never int → float, so "first integer cell_data is the tag" conventions
 * survive), moving instead of copying when the dtype already matches:
 *
 *  - `void AssignPoints(NDArray points)` — float dtype, shape `(n, dim)`.
 *  - `void AddCellBlock(std::string type, NDArray conn)` — integer dtype,
 *    shape `(n, nodes_per_cell)`.
 *  - `void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows)`
 *  - `void AddPolyhedronBlock(std::string type, std::vector<std::vector<std::vector<std::int64_t>>>
 * cells)`
 *  - `void AddPointData(std::string name, NDArray data)` /
 *    `AddFieldData(std::string name, NDArray data)` — insert-or-assign.
 *  - `void AddCellData(std::string name, std::vector<NDArray> blocks)` — one
 *    array per cell block, in block order.
 *  - `void AppendCellData(std::string name, NDArray block)` — per-block
 *    incremental variant (the medit/stl pattern).
 *
 * Writer-side accessors — cheap, but `Points()`/`Conn()`/data lookups should
 * be hoisted out of per-element hot loops (under KRATOS they may gather into
 * a lazily-built cache on first call):
 *
 *  - `std::size_t NumPoints() const`, `std::size_t PointDim() const`
 *  - `const NDArray& Points() const`
 *  - `std::size_t NumCellBlocks() const`
 *  - `CellView Cells(std::size_t i) const` — a cheap value type with
 *    `Type()` (meshio name), `NumCells()`, `NodesPerCell()` (0 if ragged),
 *    `IsRagged()`, `IsPolyhedron()`, `Conn()` (`const NDArray&`,
 *    rectangular blocks only), and ragged access `RowSize(cell)` /
 *    `Row(cell)` (polygon) and `NumFaces(cell)` / `Face(cell, face)`
 *    (polyhedron, returning `{ptr, size}`).
 *  - `detail::CellBlockRange<Mesh> CellRange() const` — iteration sugar:
 *    `for (const auto cb : rMesh.CellRange())`.
 *  - Data maps: `PointDataNames()` / `CellDataNames()` / `FieldDataNames()`
 *    return names **in sorted order** — this bakes the former
 *    `detail::sorted_keys` guarantee into the API so on-disk field order
 *    stays byte-identical across backends; `NumPointData()` /
 *    `NumCellData()` / `NumFieldData()`; `HasPointData(name)` /
 *    `HasCellData(name)` / `HasFieldData(name)`; `PointData(name)` /
 *    `FieldData(name)` (`const NDArray&`), `CellData(name, block)`
 *    (`const NDArray&`, one per cell block).
 */

// System includes
#include <cstddef>

namespace meshioplusplus {

/**
 * @brief Name of the compiled-in mesh backend, mirroring
 * `parallel_backend_name()`.
 * @return `"meshio"`, `"native"`, or `"kratos"`.
 */
constexpr const char* mesh_backend_name() {
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
    return "native";
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
    return "kratos";
#else
    return "meshio";
#endif
}

namespace detail {

/**
 * @brief Index-based range over a mesh's cell blocks, yielding
 * `TMesh::CellView` values.
 *
 * Backend-agnostic: it only requires `NumCellBlocks()` and `Cells(i)`, so a
 * single template serves every backend. Obtain one via `Mesh::CellRange()`.
 * @tparam TMesh The mesh backend type.
 */
template <class TMesh>
class CellBlockRange {
public:
    explicit CellBlockRange(const TMesh& rMesh) : mpMesh(&rMesh) {}

    class Iterator {
    public:
        Iterator(const TMesh* pMesh, std::size_t index) : mpMesh(pMesh), mIndex(index) {}
        auto operator*() const { return mpMesh->Cells(mIndex); }
        Iterator& operator++() {
            ++mIndex;
            return *this;
        }
        bool operator!=(const Iterator& rOther) const { return mIndex != rOther.mIndex; }

    private:
        const TMesh* mpMesh;
        std::size_t mIndex;
    };

    Iterator begin() const { return Iterator(mpMesh, 0); }
    Iterator end() const { return Iterator(mpMesh, mpMesh->NumCellBlocks()); }

private:
    const TMesh* mpMesh;
};

}  // namespace detail
}  // namespace meshioplusplus
