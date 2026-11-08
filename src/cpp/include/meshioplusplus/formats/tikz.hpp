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
 * @brief TikZ/PGF (LaTeX) 2D mesh writer and drawing reader.
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
 * `read_tikz` reads back only the structure the writer emits (a drawing is not
 * a mesh): see its documentation below.
 */

// System includes
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/formats/drawing_options.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/// Options of `write_tikz`. Since ABI 42 (before, positional parameters).
struct TikzWriteOptions {
    /// printf-style float format for coordinates without the leading '%'
    std::string mFloatFmt = ".6f";
    /// wrap the picture in a compilable `\documentclass{standalone}` document;
    /// otherwise emit only the `tikzpicture` environment
    bool mStandalone = true;
    /// TikZ line width (e.g. `"0.4pt"`); `std::nullopt` uses TikZ's default
    std::optional<std::string> mLineWidth;
    std::string mFill = "gray!30";  ///< xcolor fill spec for the filled faces
    std::string mDraw = "black";    ///< xcolor spec for the edge stroke
    /// optional `\begin{tikzpicture}[scale=...]` factor; `std::nullopt` emits no
    /// scale key
    std::optional<double> mScale;
    DrawingView mView;  ///< the camera, for 3D input
    /// colouring by data; `mNanColor` is an xcolor spec
    DrawingColoring mColoring{{},           std::nullopt, "viridis", std::nullopt,
                              std::nullopt, "gray",       false};
};

/**
 * @brief Write a mesh's `line`/`triangle`/`quad` cells as a TikZ figure.
 *
 * @param rPath    filesystem path to write
 * @param rMesh    the mesh to write (only line/triangle/quad contribute)
 * @param rOptions the figure's style, camera and colouring (`TikzWriteOptions`)
 * @throws WriteError on an unopenable output path
 * @throws std::invalid_argument for an unknown array name or colormap, an
 *         out-of-range component, or `vmin > vmax`
 */
MESHIOPLUSPLUS_API void write_tikz(const std::string& rPath, const Mesh& rMesh,
                                   const TikzWriteOptions& rOptions = {});

/**
 * @brief Read a TikZ figure written by `write_tikz` back as a planar mesh.
 *
 * Inverts exactly what the writer emits and refuses everything else with a
 * `ReadError` naming the construct. One `\draw[...] (x,y) -- ... -- cycle;` is a
 * face (3 vertices a `triangle`, 4 a `quad`, more a `polygon`) and one
 * `\draw[...] (x,y) -- (x,y);` a `line`; consecutive cells of one type and size
 * form one block. Points are welded on the exact text of their coordinates, in
 * order of first appearance, and are taken verbatim (there is no y-flip and the
 * picture's `scale=` is not applied), so planar input comes back exactly in x and
 * y. The colour bar (`\fill ... rectangle` and `\node`) is skipped; a
 * `fill={rgb,255:red,R;green,G;blue,B}` comes back as the `int32` `(n, 3)` cell
 * array `tikz:fill` (`-1 -1 -1` where a cell has none). `z` is always 0, and the
 * read warns once that depth is lost.
 *
 * `\path`, `plot`, `circle`, `.. controls ..`, named, polar and relative (`++`)
 * coordinates, `\foreach` and every other construct are refused.
 *
 * @param rPath    the file to read
 * @param rOptions `mLenient` is `strict=False`: without it a file that lacks the
 *                 meshio++ marker is refused; `mPointsOnly` and `mDataArrays`
 *                 narrow `tikz:fill` away
 * @throws ReadError for anything but a drawing meshio++ wrote
 */
MESHIOPLUSPLUS_API Mesh read_tikz(const std::string& rPath, const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
