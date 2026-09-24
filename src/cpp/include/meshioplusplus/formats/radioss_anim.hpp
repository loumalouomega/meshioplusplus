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
 * @file radioss_anim.hpp
 * @brief OpenRadioss / Radioss animation file (`<run>A001`, `A002`...) reader.
 *
 * Each animation file is one state of an explicit run, big-endian with the
 * magic `0x542C`: its time, three titles and ten flags, then the 2-D section
 * (the nodes and the 4-node facets, with nodal scalars, nodal vectors and 2-D
 * element tensors) and, as the flags say, the 3-D (8-node bricks), 1-D (2-node
 * elements) and SPH sections, masses, node and element ids, the part, subset,
 * material and property hierarchy and the time-history lists. The layout
 * follows OpenRadioss's `anim_to_vtk` converter (MIT); no code is copied.
 *
 * The mesh:
 *  - the points are the state's coordinates; `field_data["meshio:time"]` is its
 *    time, so a glob over a run's files is a transient sequence;
 *  - cells in the file's family order (1-D `line`, 2-D `quad`/`triangle`, 3-D
 *    bricks collapsed to the tetra, pyramid or wedge they stand for, SPH
 *    `vertex`), one block per type in order of first appearance;
 *  - nodal scalars and vectors as point data under their names, and
 *    `radioss:node_id`, `radioss:mass` when the file has them;
 *  - element scalars, tensors (xx yy zz xy yz zx; a 2-D tensor's out-of-plane
 *    components NaN) and 1-D force/moment sets (nine components) as cell data,
 *    NaN on the families that do not have them; `radioss:part`,
 *    `radioss:alive` (0 once an element is deleted), and when the file has
 *    them `radioss:element_id`, `radioss:mass`, `radioss:material`,
 *    `radioss:property`;
 *  - parts as `Cell` regions (tag = part id).
 * See doc/formats/radioss_anim.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Whether @p rPath's basename is an animation file's: a stem, then `A`
 *        and three or more digits (`crashA001`).
 */
MESHIOPLUSPLUS_API bool is_radioss_anim_filename(const std::string& rPath);

/**
 * @brief Read one OpenRadioss animation file.
 * @param rPath filesystem path to read
 * @return the state's mesh and data
 * @throws ReadError if the file can't be read, is truncated, has another magic
 *         number or names a node out of range
 */
MESHIOPLUSPLUS_API Mesh read_radioss_anim(const std::string& rPath);

/**
 * @brief Summarize one animation file: its one time step (`mTimeValues`), so a
 *        sequence takes each file's own time rather than its name's number.
 *
 * The file is read whole (`mFellBackToFullRead`).
 * @param rPath filesystem path to read
 * @param rOpts read options (unused)
 */
MESHIOPLUSPLUS_API MeshMetadata read_radioss_anim_metadata(const std::string& rPath,
                                                           const ReadOptions& rOpts = {});

}  // namespace meshioplusplus
