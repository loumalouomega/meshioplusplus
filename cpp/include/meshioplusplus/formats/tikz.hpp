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
 * @file tikz.hpp
 * @brief TikZ/PGF (LaTeX) 2D mesh writer (write-only).
 *
 * Draws the mesh's `line`/`triangle`/`quad` cells as `\draw` commands inside a
 * `tikzpicture` environment. By default it emits a full, directly
 * `pdflatex`-compilable `standalone` document; with `standalone=false` it emits
 * only the bare `tikzpicture` snippet for `\input` into a larger document. It is
 * the LaTeX counterpart to the SVG writer; unlike SVG there is no y-flip (TikZ
 * uses the math convention, y-up). Points must be 2D or flat 3D (all z ~ 0); a
 * non-flat mesh raises `WriteError`. Non-line/triangle/quad cells are silently
 * skipped, and no point_data/cell_data/field_data is emitted.
 */

// System includes
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh's `line`/`triangle`/`quad` cells as a TikZ figure.
 *
 * @param rPath       filesystem path to write
 * @param rMesh       the mesh to write (only line/triangle/quad contribute)
 * @param rFloatFmt   printf-style float format for coordinates without the
 *                    leading '%' (e.g. `".6f"`)
 * @param Standalone  when true, wrap the picture in a compilable
 *                    `\documentclass{standalone}` document; otherwise emit only
 *                    the `tikzpicture` environment
 * @param rLineWidth  TikZ line width (e.g. `"0.4pt"`); `std::nullopt` uses TikZ's
 *                    default
 * @param rFill       xcolor fill spec for the filled faces
 * @param rDraw       xcolor spec for the edge stroke
 * @param rScale      optional `\begin{tikzpicture}[scale=...]` factor;
 *                    `std::nullopt` emits no scale key
 * @throws WriteError on an unopenable output path or a non-flat 3D mesh
 */
void write_tikz(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt = ".6f",
                bool Standalone = true, const std::optional<std::string>& rLineWidth = std::nullopt,
                const std::string& rFill = "gray!30", const std::string& rDraw = "black",
                const std::optional<double>& rScale = std::nullopt);

}  // namespace meshioplusplus
