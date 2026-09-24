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
 * @file detail/abaqus_types.hpp
 * @brief Abaqus element type names and the meshio++ cells they map to, shared
 *        by the `.inp` reader/writer and the `.fil` results reader.
 *
 * `abaqus_type_table` is the `.inp` table, in source order (the writer's
 * meshio++ -> Abaqus inverse keeps the last entry per cell type). The Python
 * twin is `abaqus_to_meshio_type` in `abaqus/_abaqus.py`.
 *
 * `abaqus_cell_type` also recognises element families the table does not list
 * (`CPE8R`, `DC3D20`, `SC8R`, `M3D4`, ...): a results file names the exact
 * element, so the `.fil` reader maps it by family and node count -- the node
 * lists of solid, plane, axisymmetric, shell and membrane elements are in
 * meshio++'s order for every shape it maps.
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// The `(Abaqus element type, meshio++ cell type)` rows of the `.inp` format.
MESHIOPLUSPLUS_API const std::vector<std::pair<std::string, std::string>>& abaqus_type_table();

/**
 * @brief The meshio++ cell type of an Abaqus element with @p NodeCount nodes.
 *
 * Looks the name up in `abaqus_type_table`; otherwise classifies it by prefix
 * (3-D continuum `C3D`/`DC3D`/`AC3D`/`COH3D`/`SC`, planar and axisymmetric
 * continuum `CPE`/`CPS`/`CAX`/`CGAX`/`DC2D`/`DCAX`/`COH2D`, shell `S<digit>`/
 * `STRI`, membrane and rigid `M3D`/`R3D`/`SFM3D`, line `T2D`/`T3D`/`B2`/`B3`/
 * `PIPE`/`R2D`/`RB2D`/`RB3D`/`DC1D`) and picks the shape from @p NodeCount.
 * @param Name The element type, upper case, blanks trimmed (`C3D8R`).
 * @param NodeCount How many nodes the element lists.
 * @return the cell type, or an empty string when the element has no meshio++
 *         equivalent.
 */
MESHIOPLUSPLUS_API std::string abaqus_cell_type(std::string_view Name, std::size_t NodeCount);

}  // namespace detail
}  // namespace meshioplusplus
