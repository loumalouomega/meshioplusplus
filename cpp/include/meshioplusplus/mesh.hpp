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
#include <map>
#include <string>
#include <vector>

// Project includes
#include "ndarray.hpp"

namespace meshioplusplus {

/**
 * @brief One homogeneous block of cells of a single meshio cell type.
 *
 * Mirrors a Python `meshio.CellBlock`. Most blocks are *rectangular*: `data`
 * is a `(num_cells, nodes_per_cell)` integer `NDArray` of node indices.
 * Some formats, however, produce cells that cannot be described by a fixed
 * nodes-per-cell count, so `CellBlock` also carries two optional *ragged*
 * (jagged) representations:
 *
 *  - `polygon_rows` — 1-level ragged: a `"polygon"` block whose cells have
 *    varying node counts (e.g. MED POG Voronoi meshes). Row `i` is the list
 *    of node ids for cell `i`.
 *  - `polyhedron_rows` — 2-level ragged: a `"polyhedron"` block. Cell `i` is
 *    a list of faces, each face itself a list of node ids.
 *
 * Exactly one of `data`, `polygon_rows`, `polyhedron_rows` is populated per
 * block (see `is_ragged()`); the unused members are left empty, so ordinary
 * rectangular blocks (the overwhelming majority) are unaffected. Zero-copy
 * numpy conversion at the binding boundary only applies to the rectangular
 * `data` case — ragged blocks are always *copied* across the boundary, and
 * `py_to_mesh`'s `allow_ragged` flag is off by default so a rectangular-only
 * writer given a ragged mesh safely throws and triggers the Python fallback;
 * only ragged-aware bindings (e.g. MED write) opt in.
 */
struct CellBlock {
    std::string type;                  // meshio cell type, e.g. "triangle"
    NDArray data;                      // (num_cells, nodes_per_cell), integer dtype
    std::vector<std::string> tags;

    // Ragged (jagged) representations, used only for cell types whose rows do
    // not fit a rectangular buffer. Exactly one of `data` / `polygon_rows` /
    // `polyhedron_rows` is populated per block; the two ragged members are
    // empty for every rectangular block (all rectangular formats unaffected).
    //
    //  * polygon_rows    — 1-level ragged: a "polygon" block whose cells have
    //                      varying node counts (e.g. MED POG Voronoi meshes).
    //                      Row i = polygon_rows[i] = node ids of cell i.
    //  * polyhedron_rows — 2-level ragged: a "polyhedron" block. Cell i is a
    //                      list of faces; each face is a list of node ids.
    std::vector<std::vector<std::int64_t>> polygon_rows;
    std::vector<std::vector<std::vector<std::int64_t>>> polyhedron_rows;

    CellBlock() = default;
    CellBlock(std::string t, NDArray d) : type(std::move(t)), data(std::move(d)) {}

    /**
     * @brief Whether this block uses one of the ragged representations.
     * @return `true` iff `polygon_rows` or `polyhedron_rows` is non-empty.
     */
    bool is_ragged() const {
        return !polygon_rows.empty() || !polyhedron_rows.empty();
    }

    /**
     * @brief Number of cells in this block, whichever representation is active.
     * @return `polygon_rows.size()`, else `polyhedron_rows.size()`, else the
     *         first dimension of `data` (0 if `data` has no shape).
     */
    std::size_t num_cells() const {
        if (!polygon_rows.empty()) return polygon_rows.size();
        if (!polyhedron_rows.empty()) return polyhedron_rows.size();
        return data.shape().empty() ? 0 : data.shape()[0];
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
    NDArray points;                                         // (num_points, dim)
    std::vector<CellBlock> cells;

    // Field data. cell_data holds one NDArray per cell block, in cells order.
    std::map<std::string, NDArray> point_data;
    std::map<std::string, std::vector<NDArray>> cell_data;
    std::map<std::string, NDArray> field_data;

    /**
     * @brief Number of points in the mesh.
     * @return The first dimension of `points.shape()`, or 0 if unset.
     */
    std::size_t num_points() const {
        return points.shape().empty() ? 0 : points.shape()[0];
    }
};

}  // namespace meshioplusplus
