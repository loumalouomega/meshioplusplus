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
 * @file unv.hpp
 * @brief I-DEAS Universal File (`.unv` / `.uff`) C++ reader/writer.
 *
 * A universal file is a sequence of datasets, each opened by a `-1` line and a
 * dataset-number line and closed by another `-1` line. Read:
 *
 *  - **Nodes** 2411 (double precision, `D` exponents), 781 and the legacy 15. A
 *    node defined in a Cartesian 2420 coordinate system is moved into the global
 *    one (`x = M x_local + origin`, as Salome reads it); cylindrical/spherical
 *    systems are warned about and left local.
 *  - **Elements** 2412 and the legacy 780. The FE descriptor id and node count
 *    select the type (rods/beams 11, 21-25; plane, plate, membrane, axisymmetric
 *    and thin-shell triangles/quads 41-96; solids 111-119, 312); the extra
 *    orientation record of beam elements is skipped. Parabolic types are permuted
 *    from the UNV "sandwich" order (bottom ring, vertical mid-edges, top ring).
 *    Physical and material property ids become the integer cell data `unv:pid`
 *    and `unv:mid`.
 *  - **Permanent groups** 2467/2477/2452/2435 (quadruples) and 2417/2429/2430/2432
 *    (pairs); entity type 7 = node, 8 = element. Each group becomes a Point region
 *    (its nodes) and/or a Cell region (its elements) whose tag is the group number.
 *  - **Units** 164 -> field data `unv:units` and `unv:unit_factors` (not applied).
 *  - **Results** 2414 (data at nodes / on elements), 55 (nodes), 56 (elements) and
 *    58/58b (functions at nodal DOF). Every result block and every 58 abscissa
 *    sample belongs to a **step** keyed by its analysis type, its step/mode number
 *    and its time/frequency; `ReadOptions::mTimeStep` selects one. A step's value,
 *    analysis type and number are the field data `meshio:time`, `unv:analysis` and
 *    `unv:step` (only when the file has several steps or a non-zero analysis type).
 *    Symmetric tensors are reordered to meshio++'s `xx yy zz xy yz zx`; complex data
 *    becomes `<name>_real`/`<name>_imag`; entities without a value are NaN.
 *
 * See doc/formats/unv.md for the mapping tables.
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
 * @brief Side-channel carrying permanent-group-derived point/cell sets across
 *        the Mesh conversion boundary (the `point_sets`/`cell_sets` Python
 *        Mesh attributes are not part of the C++ Mesh/NDArray layer).
 *
 * Kept for compatibility: groups also arrive as the mesh's regions. On read, node
 * members (UNV entity type 7) fill `mPointSets` (0-based node indices) and element
 * members (entity type 8) fill `mCellSets` (per-cell-block lists of 0-based local
 * cell indices). On write, a set here replaces the mesh region of the same name
 * and kind.
 */
struct UnvInfo {
    std::map<std::string, std::vector<std::int64_t>> mPointSets;
    std::map<std::string, std::vector<std::vector<std::int64_t>>> mCellSets;
};

/**
 * @brief Write a mesh as a UNV file.
 *
 * Emits 164 (when the mesh carries `unv:units` and `unv:unit_factors`), 2411 (or
 * 781) nodes, 2412 elements with one FE descriptor per type (line 21, line3 24,
 * triangle 91, triangle6 92, quad 94, quad8/quad9 95, tetra 111, tetra10 118, wedge
 * 112, wedge15 113, hexahedron 115, hexahedron20 116, pyramid 312, pyramid13 114),
 * the inverse sandwich permutation, `unv:pid`/`unv:mid` as the property ids, the
 * mesh's point and cell regions as 2467 groups (a point and a cell region sharing a
 * name are one group; side regions are dropped with a warning), and `point_data` /
 * `cell_data` as results of one step described by the field data `unv:analysis`,
 * `unv:step` and `meshio:time`.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param code_aster emit the legacy 55 (nodes) / 56 (elements) datasets in single
 *        precision instead of 2414
 * @param node_dataset node dataset id to emit: `2411` (default) or `781`
 * @note unsupported cell types are warned about and skipped.
 */
MESHIOPLUSPLUS_API void write_unv(const std::string& rPath, const Mesh& rMesh,
                                  bool code_aster = false, int node_dataset = 2411);

/**
 * @brief `write_unv` with extra groups: `rInfo`'s sets are written as 2467 groups,
 *        replacing a mesh region of the same name and kind.
 */
MESHIOPLUSPLUS_API void write_unv(const std::string& rPath, const Mesh& rMesh, const UnvInfo& rInfo,
                                  bool code_aster = false, int node_dataset = 2411);

/**
 * @brief Read a UNV file, its first step of results included.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if the file cannot be read, a record is malformed or an element
 *         references an undefined node
 */
MESHIOPLUSPLUS_API Mesh read_unv(const std::string& rPath);

/** @brief `read_unv`, additionally filling the `UnvInfo` compatibility sets. */
MESHIOPLUSPLUS_API Mesh read_unv(const std::string& rPath, UnvInfo& rInfo);

/**
 * @brief Read a UNV file with read options.
 *
 * @param rOpts `mTimeStep` selects the step (`ResolveTimeStep`); `mPointsOnly` and
 *        `mDataArrays` narrow the result arrays
 */
MESHIOPLUSPLUS_API Mesh read_unv(const std::string& rPath, const ReadOptions& rOpts);

/** @brief `read_unv` with read options, additionally filling the `UnvInfo` sets. */
MESHIOPLUSPLUS_API Mesh read_unv(const std::string& rPath, UnvInfo& rInfo,
                                 const ReadOptions& rOpts);

/**
 * @brief Summarize a UNV file: the mesh, step 0's result names, and every step's
 *        value as `mTimeValues`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_unv_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
