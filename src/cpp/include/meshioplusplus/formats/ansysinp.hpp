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
 * @file ansysinp.hpp
 * @brief Ansys MAPDL coded database (`.cdb`, `CDWRITE` output) C++ reader/writer.
 *
 * Unrelated to the Fluent `.msh` format in ansys.hpp. Mirrors
 * `src/python/meshioplusplus/ansysInp/_ansysInp.py`.
 *
 *  - Blocks are sliced by their own Fortran format lines (`(3i9,6e21.13e3)`,
 *    `(1i7,2i9,6e21.13)`, `(19i10)` ...), parsed by `detail::parse_fortran_format`.
 *    `ET`/`ETBLOCK`/`KEYOPT` give each element type slot its routine number
 *    (`ET,1,186` or `ET,1,SOLID186`); `NBLOCK` the nodes; `EBLOCK` (the SOLID
 *    layout MAPDL writes) the elements; `CMBLOCK` the components.
 *  - An element becomes a cell by its routine's category (point, line, shell or
 *    plane, degenerate brick, native tetrahedron; MESH200 by KEYOPT(1)): a brick
 *    with repeated nodes is a wedge, pyramid or tetrahedron, a shell with K == L
 *    a triangle. A missing midside node (node 0) is placed at its edge midpoint.
 *  - Each cell carries `ansys:element` (routine), `ansys:type` (ET slot),
 *    `ansys:mat`, `ansys:real` and `ansys:secnum` cell data; `NODE`/`ELEM`
 *    components become point/cell regions (and fill @ref AnsysInfo too).
 *
 * The writer emits every type in the layout MAPDL expects (degenerate forms
 * expanded to the full brick or quad), keeping `ansys:*` data where it fits the
 * cell's shape, and writes point/cell regions as components.
 * See doc/formats/ansysinp.md.
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel carrying CMBLOCK-derived point/cell sets across the
 *        Mesh conversion boundary (the `point_sets`/`cell_sets` Python Mesh
 *        attributes are not part of the C++ Mesh/NDArray conversion layer).
 */
struct AnsysInfo {
    /** Component name -> node indices (0-based), from `CMBLOCK ...,NODE`. */
    std::map<std::string, std::vector<std::int64_t>> mPointSets;
    /**
     * Component name -> per-cell-block lists of local cell indices
     * (0-based), one inner list per mesh cell block in block order (the
     * order blocks were first encountered while reading `EBLOCK`), from
     * `CMBLOCK ...,ELEM`.
     */
    std::map<std::string, std::vector<std::vector<std::int64_t>>> mCellSets;
};

/**
 * @brief Read an Ansys MAPDL coded database (`.cdb`).
 *
 * @param rPath filesystem path to read
 * @param[out] rInfo receives the `CMBLOCK` components as point/cell sets
 *        (0-based indices; cell sets per cell block), keyed by name; the same
 *        components are also point/cell regions of the returned mesh
 * @return the mesh, with `ansys:element`/`type`/`mat`/`real`/`secnum` cell data
 * @throws ReadError if the file can't be opened, holds no block, a format line
 *         does not parse, an element uses an undefined element type or node, a
 *         `CMBLOCK` opens with a range end, or an element type has no meshio++
 *         cell type
 */
MESHIOPLUSPLUS_API Mesh read_ansysinp(const std::string& rPath, AnsysInfo& rInfo);

/**
 * @brief `read_ansysinp` with read options: `mLenient` skips elements whose type
 *        has no meshio++ cell (contact, target and other special elements)
 *        instead of failing.
 */
MESHIOPLUSPLUS_API Mesh read_ansysinp(const std::string& rPath, const ReadOptions& rOptions,
                                      AnsysInfo& rInfo);

/**
 * @brief Write `rMesh` as an Ansys MAPDL coded database.
 *
 * One `NBLOCK` (`(3i9,6e21.13e3)`) and one SOLID `EBLOCK` (`(19i10)`), nodes and
 * elements numbered from 1, 2-D points padded with z = 0. Each cell keeps its
 * `ansys:element`/`ansys:type` when the element's layout fits the cell's shape;
 * otherwise the default type for its meshio++ type is used (tetra 285, tetra10
 * 187, hexahedron/wedge/pyramid 185, their quadratic forms 186, triangle/quad
 * 181, triangle6/quad8 281, line 188, line3 189, vertex 21), with wedges,
 * pyramids and tetrahedra under 185/186 and triangles under 181/281 written in
 * their degenerate forms. `ansys:mat`/`real`/`secnum` are written when present
 * (else 1). Point and cell regions, then any `rInfo` set not named by a region,
 * become `CMBLOCK` components, runs of ids packed as ranges; side regions are
 * dropped with a warning.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param rInfo extra point/cell sets to emit as components
 * @throws WriteError for a cell type with no element type (polygons,
 *         polyhedra, higher-order Lagrange cells)
 */
MESHIOPLUSPLUS_API void write_ansysinp(const std::string& rPath, const Mesh& rMesh,
                                       const AnsysInfo& rInfo);

}  // namespace meshioplusplus
