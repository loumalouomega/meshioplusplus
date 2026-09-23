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
 * @file femap.hpp
 * @brief Femap neutral file (`.neu`) C++ reader and mesh writer.
 *
 * A neutral file is a sequence of data blocks, each opened by a line holding
 * `-1`, then the block id, and closed by the next `-1` line; records are
 * comma-separated. Record layouts change with the Femap version in block 100,
 * so every record is read by position and length rather than by version:
 *
 *  - `403` nodes (x, y, z are fields 11-13 in every version);
 *  - `404` elements: seven lines, then (from 4.5) one node list per non-zero
 *    list flag, each closed by a `-1` line. The topology code picks the cell
 *    type, and the nodes sit in a 20-slot "degenerate brick" layout (a
 *    tetrahedron's apex is slot 4, mid-edge nodes are the brick's). The element
 *    property and type become the `femap:property` and `femap:type` cell data;
 *  - `402` properties name the `property_<id>` cell regions by their titles;
 *  - `408` groups become point and cell regions (their node and element lists);
 *  - `450` output sets are steps (`ReadOptions::mTimeStep`), and the `451` and
 *    `1051` output vectors of the selected set become point data (nodal) or cell
 *    data (elemental), NaN where a vector has no value.
 *
 * Rigid, contact, weld and multi-list elements are skipped with a warning. See
 * doc/formats/femap.md for what is verified against which Femap versions.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read a Femap neutral file.
 * @param rPath filesystem path to read
 * @param rOpts `mTimeStep` selects the output set (0 = first, negative counts
 *        from the end); `mPointsOnly`/`mDataArrays` narrow the output vectors read
 * @return the mesh, with the selected output set's vectors as data and its value
 *         as `field_data["meshio:time"]`
 * @throws ReadError if the file can't be read, holds no nodes, a record is
 *         truncated or malformed, an id is defined twice, an element names an
 *         undefined node, or the requested output set does not exist
 */
MESHIOPLUSPLUS_API Mesh read_femap(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The mesh summary plus one time value per output set.
 * @param rPath filesystem path to read
 * @param rOpts ignored
 * @return the metadata, with `mTimeValues` holding each output set's value
 */
MESHIOPLUSPLUS_API MeshMetadata read_femap_metadata(const std::string& rPath,
                                                    const ReadOptions& rOpts);

/**
 * @brief The value of every output set, in file order.
 * @param rPath filesystem path to read
 * @return one value per `450` output set
 */
MESHIOPLUSPLUS_API std::vector<double> femap_time_values(const std::string& rPath);

/**
 * @brief Write `rMesh` as a Femap 8.2 neutral file: blocks 100, 402, 403, 404
 * and 408.
 *
 * The property of each element is `femap:property` (else 1) and its type
 * `femap:type` (else one derived from the cell type); properties take their
 * titles from the cell regions the reader makes of them. Other point and cell
 * regions become groups; results are not written. Cell types without a Femap
 * topology, side regions and data arrays are dropped with a warning.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for points of dimension above 3
 */
MESHIOPLUSPLUS_API void write_femap(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
