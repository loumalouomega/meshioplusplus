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
 * @file meshio_mesh.hpp
 * @brief The MESHIO mesh backend: `meshioplusplus::Mesh` and
 * `meshioplusplus::CellBlock`, the meshio-mirroring in-memory representation.
 *
 * This is the default mesh backend (see `mesh.hpp` for the compile-time
 * dispatch and `mesh_api.hpp` for the uniform format-facing API it
 * implements), and the only one compatible with the pybind11 extension —
 * `bindings/python/np_conversions.hpp` is written against these exact members.
 *
 * It mirrors the fields of the pure-Python `meshio.Mesh`: it is the type
 * every C++ format reader produces and every C++ format writer consumes.
 * The pybind11 binding layer (`bindings/python/np_conversions.hpp`) converts between
 * this type and the pure-Python `meshio.Mesh` at the I/O boundary, following
 * a "zero-copy at the boundary" strategy: `py_to_mesh` builds non-owning
 * `NDArray` *views* over the caller's numpy buffers (write path, no input
 * copy) and `mesh_to_py` moves each `NDArray`'s owned buffer into a capsule
 * backing a writeable numpy array (read path, no output copy).
 *
 * The conversion layer carries `points`, `cells`, `point_data`, `cell_data`,
 * `field_data` **and `regions`** (`region.hpp` — named groups of points, cells
 * or cell facets, which the Python `Mesh` exposes both as `.regions` and,
 * compatibly, as `point_sets`/`cell_sets`). It deliberately does **not** carry
 * `mesh.info` or `gmsh_periodic`, which remain Python-only attributes; the
 * formats that still need other out-of-band data carry it via a side-channel
 * struct the binding `setattr`s onto the Python `Mesh` after conversion (e.g.
 * `OpenFoamInfo` for `cell_tags`).
 */

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/map_order.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/properties.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

/**
 * @brief One homogeneous block of cells of a single meshio cell type.
 *
 * Mirrors a Python `meshio.CellBlock`. Most blocks are *rectangular*: `mData`
 * is a `(num_cells, nodes_per_cell)` integer `NDArray` of node indices.
 * Some formats, however, produce cells that cannot be described by a fixed
 * nodes-per-cell count, so `CellBlock` also carries two optional *ragged*
 * (jagged) representations:
 *
 *  - `mPolygonRows` — 1-level ragged: a `"polygon"` block whose cells have
 *    varying node counts (e.g. MED POG Voronoi meshes). Row `i` is the list
 *    of node ids for cell `i`.
 *  - `mPolyhedronRows` — 2-level ragged: a `"polyhedron"` block. Cell `i` is
 *    a list of faces, each face itself a list of node ids.
 *
 * Exactly one of `mData`, `mPolygonRows`, `mPolyhedronRows` is populated per
 * block (see `IsRagged()`); the unused members are left empty, so ordinary
 * rectangular blocks (the overwhelming majority) are unaffected. Zero-copy
 * numpy conversion at the binding boundary only applies to the rectangular
 * `mData` case — ragged blocks are always *copied* across the boundary, and
 * `py_to_mesh`'s `allow_ragged` flag is off by default so a rectangular-only
 * writer given a ragged mesh safely throws and triggers the Python fallback;
 * only ragged-aware bindings (e.g. MED write) opt in.
 */
struct CellBlock {
    std::string mType;  // meshio cell type, e.g. "triangle"
    NDArray mData;      // (num_cells, nodes_per_cell), integer dtype
    std::vector<std::string> mTags;

    // Ragged (jagged) representations, used only for cell types whose rows do
    // not fit a rectangular buffer. Exactly one of `mData` / `mPolygonRows` /
    // `mPolyhedronRows` is populated per block; the two ragged members are
    // empty for every rectangular block (all rectangular formats unaffected).
    //
    //  * mPolygonRows    — 1-level ragged: a "polygon" block whose cells have
    //                      varying node counts (e.g. MED POG Voronoi meshes).
    //                      Row i = mPolygonRows[i] = node ids of cell i.
    //  * mPolyhedronRows — 2-level ragged: a "polyhedron" block. Cell i is a
    //                      list of faces; each face is a list of node ids.
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronRows;

    CellBlock() = default;
    CellBlock(std::string t, NDArray d) : mType(std::move(t)), mData(std::move(d)) {}

    /**
     * @brief Whether this block uses one of the ragged representations.
     * @return `true` iff `mPolygonRows` or `mPolyhedronRows` is non-empty.
     */
    bool IsRagged() const { return !mPolygonRows.empty() || !mPolyhedronRows.empty(); }

    /**
     * @brief Number of cells in this block, whichever representation is active.
     * @return `mPolygonRows.size()`, else `mPolyhedronRows.size()`, else the
     *         first dimension of `mData` (0 if `mData` has no shape).
     */
    std::size_t NumCells() const {
        if (!mPolygonRows.empty())
            return mPolygonRows.size();
        if (!mPolyhedronRows.empty())
            return mPolyhedronRows.size();
        return mData.Shape().empty() ? 0 : mData.Shape()[0];
    }
};

/**
 * @brief The C++ in-memory mesh: points, cell blocks, and field data.
 *
 * Produced by every C++ format reader and consumed by every C++ format
 * writer; see the file-level comment for how this maps to/from the
 * pure-Python `meshio.Mesh` at the pybind11 boundary. Note what is
 * deliberately absent from this struct: `mesh.info` and `gmsh_periodic` are
 * Python-only attributes not represented here. Named groups of entities are
 * *not* in that list any more — they live in `mRegions` (see `region.hpp`).
 */
struct Mesh {
    NDArray mPoints;  // (num_points, dim)
    std::vector<CellBlock> mCells;

    // Field data. mCellData holds one NDArray per cell block, in mCells order.
    // These are unordered_map for O(1) name lookup; where key *order* is
    // observable (Python dict order, on-disk field order) call
    // detail::sorted_keys (map_order.hpp) at the consumption site.
    std::unordered_map<std::string, NDArray> mPointData;
    std::unordered_map<std::string, std::vector<NDArray>> mCellData;
    std::unordered_map<std::string, NDArray> mFieldData;

    /**
     * @brief Invalidate the sorted-name caches behind `PointDataNames()` and
     * friends.
     *
     * Unlike the NATIVE/KRATOS backends -- whose `detail::NamedItems` owns its
     * storage and can memoize safely -- this backend's three maps are PUBLIC
     * struct members, because `bindings/python/np_conversions.hpp` inserts into
     * them directly (the one sanctioned exception to the uniform-API rule).
     * **Anything that writes those maps without going through AddPointData /
     * AddCellData / AppendCellData / AddFieldData must call this.**
     */
    void InvalidateNameCaches() const { mNamesDirty = true; }

    // Named groups of points / cells / cell facets (region.hpp). Kept sorted by
    // Region::Key() with canonical (sorted, de-duplicated) entries, so region
    // order and content are identical on every backend.
    detail::RegionList mRegions;
    detail::PropertySetList mPropertySets;

    /**
     * @brief Number of points in the mesh.
     * @return The first dimension of `mPoints.Shape()`, or 0 if unset.
     */
    std::size_t NumPoints() const { return mPoints.Shape().empty() ? 0 : mPoints.Shape()[0]; }

    // -----------------------------------------------------------------
    // Uniform format-facing API (the compile-time contract shared by all
    // mesh backends — see mesh_api.hpp). Format code must go through these
    // methods, never the public members above; on this backend every method
    // is a trivial inline forward, so there is zero cost over direct access.
    // -----------------------------------------------------------------

    /**
     * @brief Cheap, copyable view over one cell block (see `mesh_api.hpp`).
     *
     * On this backend it simply wraps a `const CellBlock*`; it stays valid
     * only while the underlying `Mesh` is alive and its `mCells` vector is
     * not resized.
     */
    class CellView {
    public:
        explicit CellView(const CellBlock& rBlock) : mpBlock(&rBlock) {}
        /** @brief The meshio cell-type name (e.g. `"triangle"`). */
        const std::string& Type() const { return mpBlock->mType; }
        /** @brief Number of cells in the block (any representation). */
        std::size_t NumCells() const { return mpBlock->NumCells(); }
        /** @brief Nodes per cell for rectangular blocks; 0 for ragged ones. */
        std::size_t NodesPerCell() const {
            return mpBlock->mData.Ndim() >= 2 ? mpBlock->mData.Shape()[1] : 0;
        }
        /** @brief Whether the block uses a ragged representation. */
        bool IsRagged() const { return mpBlock->IsRagged(); }
        /** @brief Whether the block is 2-level ragged (list of faces per cell). */
        bool IsPolyhedron() const { return !mpBlock->mPolyhedronRows.empty(); }
        /** @brief Rectangular `(num_cells, nodes_per_cell)` connectivity (empty if ragged). */
        const NDArray& Conn() const { return mpBlock->mData; }
        /** @brief Node count of polygon cell @p cell (1-level ragged blocks). */
        std::size_t RowSize(std::size_t cell) const { return mpBlock->mPolygonRows[cell].size(); }
        /** @brief Node ids of polygon cell @p cell (1-level ragged blocks). */
        const std::int64_t* Row(std::size_t cell) const {
            return mpBlock->mPolygonRows[cell].data();
        }
        /** @brief Face count of polyhedron cell @p cell (2-level ragged blocks). */
        std::size_t NumFaces(std::size_t cell) const {
            return mpBlock->mPolyhedronRows[cell].size();
        }
        /** @brief `{node ids, count}` of face @p face of polyhedron cell @p cell. */
        std::pair<const std::int64_t*, std::size_t> Face(std::size_t cell, std::size_t face) const {
            const auto& r_face = mpBlock->mPolyhedronRows[cell][face];
            return {r_face.data(), r_face.size()};
        }

    private:
        const CellBlock* mpBlock;
    };

    // --- reader-side ingestion ---

    /** @brief Takes ownership of the point array (float dtype, shape `(n, dim)`). */
    void AssignPoints(NDArray points) { mPoints = std::move(points); }
    /** @brief Appends a rectangular cell block (integer dtype, shape `(n, npc)`). */
    void AddCellBlock(std::string type, NDArray conn) {
        mCells.emplace_back(std::move(type), std::move(conn));
    }
    /** @brief Appends a 1-level ragged (polygon) cell block. */
    void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows) {
        CellBlock cb;
        cb.mType = std::move(type);
        cb.mPolygonRows = std::move(rows);
        mCells.push_back(std::move(cb));
    }
    /** @brief Appends a 2-level ragged (polyhedron) cell block. */
    void AddPolyhedronBlock(std::string type,
                            std::vector<std::vector<std::vector<std::int64_t>>> cells) {
        CellBlock cb;
        cb.mType = std::move(type);
        cb.mPolyhedronRows = std::move(cells);
        mCells.push_back(std::move(cb));
    }
    /** @brief Inserts or replaces a named per-point data array. */
    void AddPointData(std::string name, NDArray data) {
        mNamesDirty = true;
        mPointData[std::move(name)] = std::move(data);
    }
    /** @brief Inserts or replaces a named per-cell data array list (one per block). */
    void AddCellData(std::string name, std::vector<NDArray> blocks) {
        mNamesDirty = true;
        mCellData[std::move(name)] = std::move(blocks);
    }
    /** @brief Appends one block's array to a named cell-data list (creating it if new). */
    void AppendCellData(const std::string& rName, NDArray block) {
        mNamesDirty = true;
        mCellData[rName].push_back(std::move(block));
    }
    /** @brief Inserts or replaces a named field-data array. */
    void AddFieldData(std::string name, NDArray data) {
        mNamesDirty = true;
        mFieldData[std::move(name)] = std::move(data);
    }

    /** @brief Fills the sorted-name caches when `mNamesDirty`. */
    void RefreshNameCaches() const {
        if (!mNamesDirty)
            return;
        mPointNames = detail::sorted_keys(mPointData);
        mCellNames = detail::sorted_keys(mCellData);
        mFieldNames = detail::sorted_keys(mFieldData);
        mNamesDirty = false;
    }

    // Memoized sorted name lists. `mutable` so a const accessor can fill them;
    // see InvalidateNameCaches() above for the one discipline they require.
    mutable std::vector<std::string> mPointNames;
    mutable std::vector<std::string> mCellNames;
    mutable std::vector<std::string> mFieldNames;
    mutable bool mNamesDirty = true;

    // --- writer-side accessors ---

    /** @brief Spatial dimension of the points (second shape entry), or 0 if unset. */
    std::size_t PointDim() const { return mPoints.Ndim() >= 2 ? mPoints.Shape()[1] : 0; }
    /** @brief The `(num_points, dim)` point array. */
    const NDArray& Points() const { return mPoints; }
    /** @brief Number of cell blocks. */
    std::size_t NumCellBlocks() const { return mCells.size(); }
    /** @brief View over cell block @p i (in insertion order). */
    CellView Cells(std::size_t i) const { return CellView(mCells[i]); }
    /** @brief Range over all cell blocks: `for (const auto cb : mesh.CellRange())`. */
    detail::CellBlockRange<Mesh> CellRange() const { return detail::CellBlockRange<Mesh>(*this); }

    /** @brief Point-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> PointDataNames() const {
        RefreshNameCaches();
        return mPointNames;
    }
    /** @brief Number of named point-data arrays. */
    std::size_t NumPointData() const { return mPointData.size(); }
    /** @brief Whether a point-data array named @p rName exists. */
    bool HasPointData(const std::string& rName) const { return mPointData.count(rName) > 0; }
    /** @brief The point-data array named @p rName (throws if absent). */
    const NDArray& PointData(const std::string& rName) const { return mPointData.at(rName); }

    /** @brief Cell-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> CellDataNames() const {
        RefreshNameCaches();
        return mCellNames;
    }
    /** @brief Number of named cell-data array lists. */
    std::size_t NumCellData() const { return mCellData.size(); }
    /** @brief Whether a cell-data list named @p rName exists. */
    bool HasCellData(const std::string& rName) const { return mCellData.count(rName) > 0; }
    /** @brief Block @p block of the cell-data list named @p rName (throws if absent). */
    const NDArray& CellData(const std::string& rName, std::size_t block) const {
        return mCellData.at(rName)[block];
    }
    /** @brief Number of blocks in the cell-data list named @p rName (throws if absent). */
    std::size_t CellDataNumBlocks(const std::string& rName) const {
        return mCellData.at(rName).size();
    }

    /** @brief Field-data names in sorted order (drives on-disk field order). */
    std::vector<std::string> FieldDataNames() const {
        RefreshNameCaches();
        return mFieldNames;
    }
    /** @brief Number of named field-data arrays. */
    std::size_t NumFieldData() const { return mFieldData.size(); }
    /** @brief Whether a field-data array named @p rName exists. */
    bool HasFieldData(const std::string& rName) const { return mFieldData.count(rName) > 0; }
    /** @brief The field-data array named @p rName (throws if absent). */
    const NDArray& FieldData(const std::string& rName) const { return mFieldData.at(rName); }

    // --- named regions (region.hpp) ---
    // NOTE: `Region(i)` hides the type name `meshioplusplus::Region` for the
    // rest of this class body — every declaration below must qualify it.

    /// Sentinel returned by `FindRegion` when no region matches.
    static constexpr std::size_t npos = detail::RegionList::npos;

    /** @brief Adds a region, replacing one with the same (kind, name, dim, tag). */
    void AddRegion(meshioplusplus::Region region) { mRegions.Add(std::move(region)); }

    // --- uniform API: property sets (properties.hpp) -----------------------
    //
    // Keyed by id, never by entity index, so no operation has to remap them.
    // Named GetPropertySet rather than PropertySet(i) on purpose: an accessor
    // called PropertySet would hide the TYPE of that name for the rest of this
    // class body, the trap mesh_api.hpp documents for Region(i).

    /** @brief Store a properties block, replacing one with the same `mId`. */
    void AddPropertySet(PropertySet propertySet) { mPropertySets.Add(std::move(propertySet)); }
    /** @brief Number of properties blocks. */
    std::size_t NumPropertySets() const { return mPropertySets.Size(); }
    /** @brief Properties block @p Index, in ascending-`mId` order. */
    const PropertySet& GetPropertySet(std::size_t Index) const { return mPropertySets.At(Index); }
    /** @brief Whether a properties block with this id exists. */
    bool HasPropertySet(std::int64_t Id) const { return mPropertySets.Has(Id); }
    /** @brief Index of the block with this id, or `npos`. */
    std::size_t FindPropertySet(std::int64_t Id) const { return mPropertySets.Find(Id); }

    /** @brief Number of named regions. */
    std::size_t NumRegions() const { return mRegions.Size(); }
    /** @brief Region @p i in `(kind, name, dim, tag)` order. */
    const meshioplusplus::Region& Region(std::size_t i) const { return mRegions.At(i); }
    /** @brief Region names in sorted order, de-duplicated across kinds. */
    std::vector<std::string> RegionNames() const { return mRegions.Names(); }
    /** @brief Whether any region carries this name, regardless of kind. */
    bool HasRegion(const std::string& rName) const { return mRegions.Has(rName); }
    /** @brief Whether a region of this name and kind exists. */
    bool HasRegion(const std::string& rName, RegionKind kind) const {
        return mRegions.Has(rName, kind);
    }
    /** @brief Index of the region with this name and kind, or `npos`. */
    std::size_t FindRegion(const std::string& rName, RegionKind kind) const {
        return mRegions.Find(rName, kind);
    }
};

}  // namespace meshioplusplus
