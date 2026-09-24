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
 * @file radioss.hpp
 * @brief OpenRadioss / Altair Radioss starter deck (`*_0000.rad`) reader.
 *
 * A Radioss starter deck ("block format"): `#RADIOSS STARTER`, `/BEGIN` (run
 * name, input version `Invers`, units), then `/KEYWORD/option/id` blocks up to
 * `/END`. Lines starting with `#` or `$` are comments; `#include file` inlines
 * another file (relative to the including one, at most 8 deep). Fields are 10
 * columns wide (integers) and 20 (reals) from format 51 on, 8 and 16 before; a
 * line with a comma is split on commas.
 *
 * The mesh comes from `/NODE` and the element blocks `/BRICK`, `/PENTA6`,
 * `/TETRA4`, `/TETRA10`, `/BRIC20`, `/SHELL`, `/SH3N`, `/QUAD`, `/TRIA`,
 * `/BEAM`, `/TRUSS` and `/SPRING` (the part id is the keyword's last field).
 * A `/BRICK` with repeated nodes is the tetrahedron, pyramid or wedge it
 * stands for (`detail/degenerate_solid.hpp`); a `/SHELL` whose last two nodes
 * coincide is a triangle; a `/BRIC20` with a zero mid-edge node keeps its
 * corners. `/BRIC20` lists its vertical mid-edges before the top ring (the
 * `"radioss"` table of `detail/node_order.hpp`); the other types are in
 * meshio++'s order. Each cell's part, property and material are the
 * `radioss:part`/`radioss:property`/`radioss:material` cell data.
 *
 * Regions: every `/PART` is a `Cell` region named by its title (tag = part
 * id); `/SUBSET` a `Cell` region of the parts that name it or its children;
 * `/GRNOD` a `Point` region and `/GRBRIC`, `/GRSHEL`, `/GRSH3N`, `/GRQUAD`,
 * `/GRTRIA`, `/GRBEAM`, `/GRTRUS`, `/GRSPRI` `Cell` regions, from entity ids
 * (a negative id removes one), part ids or other groups of the same keyword;
 * `/SURF/SEG` a `Side` region (a shell segment is the shell's own face). Other
 * group generators (`BOX`, `GENE`, ...) and every non-mesh keyword (materials,
 * properties, loads, contacts) are skipped. Units are not applied; the input
 * version is `field_data["radioss:version"]`. See doc/formats/radioss.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Read an OpenRadioss / Radioss starter deck.
 * @param rPath filesystem path to read
 * @return the mesh, with parts, subsets, groups and surfaces as regions
 * @throws ReadError if the file can't be read, is an engine deck, a field is
 *         malformed, an element is cut short or names an undefined node, or a
 *         node or element id is defined twice
 */
MESHIOPLUSPLUS_API Mesh read_radioss(const std::string& rPath);

}  // namespace meshioplusplus
