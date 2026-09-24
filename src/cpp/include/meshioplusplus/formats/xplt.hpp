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
 * @file xplt.hpp
 * @brief FEBio plot file (`.xplt`) C++ reader: the results FEBio writes.
 *
 * An `.xplt` file is a tree of tagged chunks (`u32 id, u32 size, payload`)
 * after a 4-byte magic: a root holding the header and the dictionary of plot
 * variables, the mesh, then one state per saved time. With compression on,
 * every chunk after the first mesh is its own zlib stream. Plot versions
 * 0x0030 and later (FEBio 3 and 4) are read, either byte order.
 *
 *  - Each domain is a cell block and a `Cell` region (named after the domain,
 *    else the element set holding exactly its elements, else its part), tagged
 *    with its part id. hex27 goes through the `"febio"` node-order table.
 *  - Node sets, element sets and surfaces become `Point`, `Cell` and `Side`
 *    regions, as in `.feb`. Each data surface is also a block of facet cells
 *    per facet type, after the domains, with a `Cell` region named after it
 *    and `cell_data["xplt:surface"]` (its 1-based id, 0 on the domains).
 *  - A remeshed run writes a new mesh before the states that use it; each
 *    state is read on its own mesh. Metadata describes the first.
 *  - `ReadOptions::mTimeStep` picks the state. Its time is
 *    `field_data["meshio:time"]`, its index `"xplt:step"`, its status
 *    `"xplt:status"`. Nodal variables are point data; per-element variables
 *    (`FMT_ITEM`) cell data, NaN on domains without them; per-region ones
 *    (`FMT_REGION`) cell data repeated over the domain; per-element-node ones
 *    (`FMT_NODE`, `FMT_MULT`) point data, averaged over the elements sharing a
 *    node. Global variables are field data. Surface variables live on the facet
 *    cells the same way (per facet, per surface, or averaged per facet node);
 *    one on no surface of the mesh gives no array. Edge and material-point
 *    variables are not read (a warning names them).
 *  - Symmetric tensors keep FEBio's 6 components (xx, yy, zz, xy, yz, xz).
 *
 * See doc/formats/xplt.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one state of an FEBio plot file.
 * @param rPath filesystem path to read
 * @param rOptions `mTimeStep` picks the state (negative counts from the end);
 *        `mPointsOnly`/`mDataArrays` narrow the data read; `mLenient`
 *        downgrades tet5/tet15 domains
 * @return the mesh with the chosen state's data
 * @throws ReadError for a bad magic, an unsupported plot version, a
 *         truncated mesh, or an out-of-range time step
 */
MESHIOPLUSPLUS_API Mesh read_xplt(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarise an FEBio plot file, its state times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_xplt_metadata(const std::string& rPath,
                                                   const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
