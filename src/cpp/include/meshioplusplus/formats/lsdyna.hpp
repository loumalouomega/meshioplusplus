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
#pragma once

/**
 * @file lsdyna.hpp
 * @brief LS-DYNA keyword input deck (`.k` / `.key` / `.dyn`) C++ reader/writer.
 *
 * The reader collects `*NODE`, `*ELEMENT_SOLID`/`_SHELL`/`_TSHELL`/`_BEAM`/
 * `_DISCRETE`/`_MASS`, `*PART` and `*SET_*` cards from the deck and from the files
 * it `*INCLUDE`s into one state and resolves ids only once everything is read, so a
 * set may refer to elements defined in a later file.
 *
 *  - `*PART` becomes a `Cell` region: the title is the name, `pid` the tag, and the
 *    dimension is the highest of its elements'.
 *  - `*SET_NODE` becomes a `Point` region, `*SET_SOLID`/`_SHELL`/`_TSHELL`/`_BEAM`/
 *    `_DISCRETE`/`_PART` a `Cell` region and `*SET_SEGMENT` a `Side` region (each
 *    segment is matched to the cell face it names). Set regions have dimension -1
 *    and the set id as tag.
 *  - LS-DYNA has no tetra, pyramid or wedge card: they are hexahedra with repeated
 *    nodes, collapsed on read and expanded on write. A triangular shell is a quad
 *    with `n4 == n3`.
 *  - Standard, `LONG=`, `I10=` and comma-separated card formats are all read, per
 *    card (see `detail/keyword_card.hpp`).
 *
 * Materials, sections, contacts, loads and every other keyword are not parsed. See
 * doc/formats/lsdyna.md for the full mapping and its limits.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `rMesh` as an LS-DYNA keyword deck.
 *
 * Emits `*NODE`, one `*ELEMENT_*` keyword per cell block, one `*PART` per part and
 * one `*SET_*_LIST_TITLE` per remaining region, with placeholder section and
 * material ids. Cell regions with a dimension become parts; a cell that no part
 * claims goes to one part per cell block.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for a cell type with no LS-DYNA element card, or when an id no
 *         longer fits its 8-column field
 */
MESHIOPLUSPLUS_API void write_lsdyna(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an LS-DYNA keyword deck.
 * @param rPath filesystem path to read
 * @return the mesh, with parts and sets as regions
 * @throws ReadError if the file can't be read, a field is malformed, an element
 *         references an undefined node, an id is duplicated, or `*INCLUDE` nests
 *         deeper than eight levels
 */
MESHIOPLUSPLUS_API Mesh read_lsdyna(const std::string& rPath);

}  // namespace meshioplusplus
