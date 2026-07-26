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
 * @file dex.hpp
 * @brief FLUX field file (.dex) C++ reader/writer.
 *
 * A DEX file stores a single nodal field: a two-line `#`-delimited header
 * (`NAME`/`FORMULA` and `NB_REAL`/`NB_COMP`/`NB_POINT`), then one row per
 * point holding the point coordinates (x y z) followed by its NB_COMP field
 * values. Read here as a geometry-less Mesh (no cells) whose `points` come
 * from the coordinates and whose `point_data[<field>]` holds the values.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read a FLUX field file (.dex) into a geometry-less Mesh. */
MESHIOPLUSPLUS_API Mesh read_dex(const std::string& rPath);

/** @brief Write a mesh's first nodal field as a FLUX field file (.dex). */
MESHIOPLUSPLUS_API void write_dex(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
