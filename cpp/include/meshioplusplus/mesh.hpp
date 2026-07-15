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
 * @file mesh.hpp
 * @brief `meshioplusplus::Mesh` and `meshioplusplus::CellBlock`: the C++ in-memory
 * mesh representation.
 *
 * This mirrors the fields of the pure-Python `meshio.Mesh`: it is the type
 * every C++ format reader produces and every C++ format writer consumes.
 * The pybind11 binding layer (`bindings/np_conversions.hpp`) converts between
 * this type and the pure-Python `meshio.Mesh` at the I/O boundary, following
 * a "zero-copy at the boundary" strategy: `py_to_mesh` builds non-owning
 * `NDArray` *views* over the caller's numpy buffers (write path, no input
 * copy) and `mesh_to_py` moves each `NDArray`'s owned buffer into a capsule
 * backing a writeable numpy array (read path, no output copy).
 *
 * The conversion layer carries `points`, `cells`, `point_data`, `cell_data`,
 * and `field_data` — but deliberately **not** `mesh.info`, `cell_sets`, or
 * `point_sets`, which are custom attributes that live only on the Python
 * `Mesh`. Formats that need those either defer entirely to the Python
 * fallback or carry the extra data out-of-band via a side-channel struct
 * that the binding `setattr`s onto the Python `Mesh` object after
 * conversion (e.g. `MedInfo`/`AnsysInfo` for `point_sets`/`cell_sets`,
 * `OpenFoamInfo` for `cell_tags`).
 */

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "ndarray.hpp"

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
 * deliberately absent from this struct: `mesh.info`, `point_sets`, and
 * `cell_sets` are Python-only attributes not represented here (they travel,
 * when needed, through a per-format side-channel struct instead).
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
     * @brief Number of points in the mesh.
     * @return The first dimension of `mPoints.Shape()`, or 0 if unset.
     */
    std::size_t NumPoints() const { return mPoints.Shape().empty() ? 0 : mPoints.Shape()[0]; }
};

}  // namespace meshioplusplus
