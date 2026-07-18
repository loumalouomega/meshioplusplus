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
 * @file svg.hpp
 * @brief SVG (Scalable Vector Graphics) 2D mesh writer (write-only).
 *
 * Draws the mesh's `line`/`triangle`/`quad` cells as `<path>` elements in a
 * single `<svg>` document — a flat-2D visualization format with no reader.
 * Points must be 2D or flat 3D (all z ~ 0); a genuinely non-flat mesh raises
 * `WriteError`. The y-axis is flipped (`max_y + min_y - y`) to convert the
 * mesh/math convention (y-up) to SVG's screen convention (y-down). Any cell
 * type other than line/triangle/quad is silently skipped. No
 * point_data/cell_data/field_data is emitted.
 */

// System includes
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh's `line`/`triangle`/`quad` cells as an SVG document.
 *
 * @param rPath        filesystem path to write
 * @param rMesh        the mesh to write (only line/triangle/quad contribute)
 * @param rFloatFmt    printf-style float format for coordinates without the
 *                     leading '%' (e.g. `".3f"`)
 * @param rStrokeWidth explicit stroke width; `std::nullopt` auto-computes it as
 *                     1% of the on-canvas width
 * @param rImageWidth  output width in user units; `std::nullopt` keeps the
 *                     mesh's own width (no scaling)
 * @param rFill        cell fill colour
 * @param rStroke      edge stroke colour
 * @throws WriteError on an unopenable output path or a non-flat 3D mesh
 */
void write_svg(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt = ".3f",
               const std::optional<std::string>& rStrokeWidth = std::nullopt,
               const std::optional<double>& rImageWidth = 100.0,
               const std::string& rFill = "#c8c5bd", const std::string& rStroke = "#000080");

}  // namespace meshioplusplus
