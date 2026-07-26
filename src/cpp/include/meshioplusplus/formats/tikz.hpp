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
 * uses the math convention, y-up). 2D or flat-3D input (all z ~ 0) takes the
 * classic flat path, byte-identical to previous releases. A genuinely
 * non-flat 3D mesh is rendered: the boundary skin of any supported volume
 * cells is extracted first (`extract_skin`, linearized), then the surface
 * faces are projected with an orthographic camera (azimuth/elevation/roll,
 * default the classic CAD isometric view) and drawn back-to-front
 * (painter's algorithm — see `detail/projection.hpp`). Non-drawable cells
 * are silently skipped, and no point_data/cell_data/field_data is emitted.
 */

// System includes
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
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
 * @param azimuth     camera azimuth in degrees (3D input only)
 * @param elevation   camera elevation in degrees (3D input only); the default
 *                    pair (45, atan(1/sqrt(2))) is the classic CAD isometric
 *                    view
 * @param roll        in-screen camera roll in degrees (3D input only)
 * @param rColorBy    name of a `point_data` or `cell_data` array to colour the
 *                    faces by; empty (the default) keeps the flat `rFill` and
 *                    leaves the output byte-identical to previous releases.
 *                    Point data colours a face by the mean of its corner
 *                    values, cell data by its owning cell's value -- found
 *                    through `"surface:parent_cell"` for a projected volume
 *                    mesh, whose skin is extracted with provenance in that case
 * @param rComponent  component of a multi-component array; `std::nullopt` uses
 *                    the row magnitude
 * @param rCmap       built-in colormap name (see `detail/colormap.hpp`)
 * @param rVMin       low end of the mapped range; `std::nullopt` uses the
 *                    smallest finite value among the *drawn* faces
 * @param rVMax       high end of the mapped range; `std::nullopt` likewise
 * @param rNanColor   xcolor spec for faces whose value is NaN or infinite
 * @param Colorbar    when true, append a gradient bar with min/max labels to
 *                    the right of the figure
 * @throws WriteError on an unopenable output path
 * @throws std::invalid_argument for an unknown array name or colormap, an
 *         out-of-range component, or `vmin > vmax`
 */
MESHIOPLUSPLUS_API void write_tikz(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt = ".6f",
                bool Standalone = true, const std::optional<std::string>& rLineWidth = std::nullopt,
                const std::string& rFill = "gray!30", const std::string& rDraw = "black",
                const std::optional<double>& rScale = std::nullopt, double azimuth = 45.0,
                double elevation = 35.264389682754654, double roll = 0.0,
                const std::string& rColorBy = "",
                const std::optional<int>& rComponent = std::nullopt,
                const std::string& rCmap = "viridis",
                const std::optional<double>& rVMin = std::nullopt,
                const std::optional<double>& rVMax = std::nullopt,
                const std::string& rNanColor = "gray", bool Colorbar = false);

}  // namespace meshioplusplus
