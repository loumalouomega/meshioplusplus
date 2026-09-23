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
 * @file code_aster.hpp
 * @brief Code_Aster native mesh (`.mail`) C++ reader/writer.
 *
 * The `.mail` format is Code_Aster's own ASCII mesh (`LIRE_MAILLAGE(FORMAT=
 * 'ASTER')`, document U3.01.00). It is a sequence of blocks, each opened by a
 * keyword and closed by `FINSF`; the file ends at `FIN`. Tokens are separated by
 * spaces or commas, `%` starts a comment, and Code_Aster reads only the first 80
 * columns of a line, so this reader does too.
 *
 *  - `COOR_1D`/`COOR_2D`/`COOR_3D` hold named nodes; the point dimension follows
 *    the keyword.
 *  - Each element block (`POI1`, `SEG2/3/4`, `TRIA3/6/7`, `QUAD4/8/9`,
 *    `TETRA4/10`, `PENTA6/15/18`, `PYRAM5/13`, `HEXA8/20/27`) becomes one cell
 *    block. Records are read as a token stream, so they may wrap across lines.
 *    Node order goes through the `"code_aster"` tables of
 *    `detail/node_order.hpp`; it is **not** MED's order.
 *  - `GROUP_MA` becomes a `Cell` region and `GROUP_NO` a `Point` region, named by
 *    `NOM=` or by the block's first token. Neither carries a tag.
 *  - `TITRE`, `DUMP`, `DEBUG` and the keywords Code_Aster itself ignores are
 *    skipped; any other keyword is skipped with a warning.
 *
 * Nodes and elements are identified by name in the file. The names are not kept:
 * the writer names them `N1…` and `M1…` again. See doc/formats/code_aster.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Read a Code_Aster `.mail` mesh.
 * @param rPath filesystem path to read
 * @return the mesh, with `GROUP_MA`/`GROUP_NO` as regions
 * @throws ReadError if the file can't be read, a block is not closed, a value is
 *         malformed, a node or element name is defined twice, or an element or a
 *         group names an undefined node or element
 */
MESHIOPLUSPLUS_API Mesh read_code_aster(const std::string& rPath);

/**
 * @brief Write `rMesh` as a Code_Aster `.mail` mesh.
 *
 * Emits one coordinate block, one element block per cell block, a `GROUP_NO` per
 * point region and a `GROUP_MA` per cell region, with every line within 80
 * columns. Group names are sanitised to letters, digits and `_` and truncated to
 * 24 characters (collisions get a numbered suffix), each change with a warning.
 * Side regions and data arrays have no `.mail` equivalent and are dropped with a
 * warning.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for a cell type with no `.mail` keyword, a point dimension
 *         above 3, or more than 9,999,999 nodes or elements (names are limited to
 *         8 characters)
 */
MESHIOPLUSPLUS_API void write_code_aster(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
