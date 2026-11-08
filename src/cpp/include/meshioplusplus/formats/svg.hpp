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
 * @brief SVG (Scalable Vector Graphics) mesh writer and drawing reader.
 *
 * Draws the mesh's `line`/`triangle`/`quad` cells as `<path>` elements in a
 * single `<svg>` document. `read_svg` reads back only what the writer emits (a
 * drawing is not a mesh): see its documentation below. For 2D or
 * flat-3D input (all z ~ 0) the classic flat path is used, byte-identical to
 * previous releases. A genuinely non-flat 3D mesh is *rendered*: if it
 * contains supported volume cells its boundary skin is extracted first
 * (`extract_skin`, linearized), then the surface faces are projected with an
 * orthographic camera (azimuth/elevation/roll, default the classic CAD
 * isometric view) and painted back-to-front (painter's algorithm on face
 * centroid depth — see `detail/projection.hpp`). The y-axis is flipped
 * (`max_y + min_y - y`) to convert the mesh/math convention (y-up) to SVG's
 * screen convention (y-down). Any non-drawable cell type is silently
 * skipped. No point_data/cell_data/field_data is emitted.
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

/// Options of `write_svg`. Since ABI 42 (before, positional parameters).
struct SvgWriteOptions {
    /// printf-style float format for coordinates without the leading '%'
    std::string mFloatFmt = ".3f";
    /// explicit stroke width; `std::nullopt` auto-computes it as 1% of the
    /// on-canvas width
    std::optional<std::string> mStrokeWidth;
    /// output width in user units; `std::nullopt` keeps the mesh's own width
    /// (no scaling)
    std::optional<double> mImageWidth = 100.0;
    std::string mFill = "#c8c5bd";    ///< cell fill colour
    std::string mStroke = "#000080";  ///< edge stroke colour
    DrawingView mView;                ///< the camera, for 3D input
    /// colouring by data. When `mColorBy` is set, each `<path>` carries its own
    /// `fill` attribute, overriding the document-level rule. The colour bar
    /// widens the `viewBox` (and only the viewBox): the scaling factor, the
    /// stroke width and every mesh coordinate are unaffected
    DrawingColoring mColoring{{},           std::nullopt, "viridis", std::nullopt,
                              std::nullopt, "#808080",    false};
};

/**
 * @brief Write a mesh's `line`/`triangle`/`quad` cells as an SVG document.
 *
 * @param rPath    filesystem path to write
 * @param rMesh    the mesh to write (only line/triangle/quad contribute)
 * @param rOptions the drawing's style, camera and colouring (`SvgWriteOptions`)
 * @throws WriteError on an unopenable output path
 * @throws std::invalid_argument for an unknown array name or colormap, an
 *         out-of-range component, or `vmin > vmax`
 */
MESHIOPLUSPLUS_API void write_svg(const std::string& rPath, const Mesh& rMesh,
                                  const SvgWriteOptions& rOptions = {});

/**
 * @brief Read an SVG drawing written by `write_svg` back as a planar mesh.
 *
 * Inverts exactly what the writer emits and refuses everything else with a
 * `ReadError` naming the construct. A closed `<path d="M x y L x y ... Z">` of
 * 3 vertices is a `triangle`, of 4 a `quad`, of more a `polygon`; an open
 * 2-vertex path is a `line`; consecutive cells of one type and size form one
 * block. Points are welded on the exact text of their coordinates, in order of
 * first appearance. The y-flip is undone and `image_width`'s scale restored from
 * the meshio++ marker comment (or, without one, from the `viewBox`); the colour
 * bar (`<rect>`/`<text>`) is skipped; a path's `fill="#rrggbb"` comes back as
 * the `int32` `(n, 3)` cell array `svg:fill` (`-1 -1 -1` where a cell has none),
 * never as the data it was coloured from. `z` is always 0, and the read warns
 * once that depth is lost.
 *
 * A `<polygon>`, a `transform`, a curve (`C`/`S`/`Q`/`T`/`A`), a relative or
 * `H`/`V` command, several subpaths and any other element are refused.
 *
 * @param rPath    the file to read
 * @param rOptions `mLenient` is `strict=False`: without it a file that lacks the
 *                 meshio++ marker is refused; `mPointsOnly` and `mDataArrays`
 *                 narrow `svg:fill` away
 * @throws ReadError for anything but a drawing meshio++ wrote
 */
MESHIOPLUSPLUS_API Mesh read_svg(const std::string& rPath, const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
