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
 * @file patran.hpp
 * @brief MSC Patran 2 neutral file (`.pat`/`.out`) C++ reader/writer.
 *
 * The neutral file is a sequence of *packets*. Each opens with a fixed-width
 * header card `(I2,8I8)` -- `IT, ID, IV, KC, N1..N5` -- followed by `KC` data
 * cards. This reader uses:
 *
 *  - `25` title and `26` summary (skipped);
 *  - `01` node: `ID` is the node id, card 1 the coordinates `(3E16.9)` (Patran
 *    always writes them in the global frame);
 *  - `02` element: `ID` is the element id and `IV` its shape (2 bar, 3 tri,
 *    4 quad, 5 tet, 6 pyramid, 7 wedge, 8 hex); the node count on card 1 tells
 *    linear from quadratic. Card 1 also carries the property id, kept as the
 *    `patran:property` cell data. Hex20 and wedge15 list their vertical mid-edge
 *    nodes before the top ring (the `"patran"` tables of
 *    `detail/node_order.hpp`); every other shape is in meshio++'s order;
 *  - `21` named component: a name card, then `(type, id)` pairs. Type 5 (node)
 *    becomes a `Point` region, the element types (6 bar ... 12 hex) a `Cell`
 *    region of the same name, tagged with the component number;
 *  - `99` end of file.
 *
 * Every other packet (materials, properties, loads, ...) is skipped by its
 * `KC`. Elements no component names are grouped by property into
 * `property_<pid>` cell regions (tag = pid). See doc/formats/patran.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Read a Patran 2 neutral file.
 * @param rPath filesystem path to read
 * @return the mesh, with named components and property groups as regions
 * @throws ReadError if the file can't be read, a packet is truncated, a field
 *         is malformed, a node or element id is defined twice, or an element
 *         names an undefined node
 */
MESHIOPLUSPLUS_API Mesh read_patran(const std::string& rPath);

/**
 * @brief Write `rMesh` as a Patran 2 neutral file.
 *
 * Emits packets 25 (the title carries the provenance line), 26, 01, 02 and 99,
 * plus one packet 21 per region name: a point region and a cell region of the
 * same name share one component. Coordinates are written `E16.9`, so they keep
 * ten significant digits. The element property comes from the
 * `patran:property` cell data (else 1). Cell types without a Patran shape,
 * side regions and other data arrays are dropped with a warning and a
 * provenance note; component names longer than 12 characters are truncated.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for points of dimension above 3 or more than 99,999,999
 *         nodes or elements (the `I8` id fields)
 */
MESHIOPLUSPLUS_API void write_patran(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
