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
 * @file operations/parameterize.hpp
 * @brief A UV atlas of a triangle surface: cut into charts, each flattened
 * without flipping a triangle, packed into the unit square.
 *
 * `parameterize(mesh)` returns the surface with `point_data["uv"]`, two
 * texture coordinates per point in `[0, 1]`. Three steps:
 *
 * 1. **Charts.** Each edge-connected piece that is a topological disc is one
 *    chart. A closed or holed genus-0 piece is cut into discs along a cut
 *    graph (`fit_cad`'s); anything else, or a disc that does not flatten
 *    without folding, is halved along the direction its normals spread most
 *    (again `fit_cad`'s split), up to `mMaxDepth` times; so is a chart that
 *    folds nothing but whose L2 stretch exceeds `mMaxStretch`. A piece still folding
 *    at the cap keeps its best map and counts its flipped triangles; one with
 *    no map at all becomes one chart per triangle.
 * 2. **Flattening.** Least-squares conformal maps (Lévy et al. 2002), else
 *    mean-value coordinates on a circle (Floater 2003), each accepted only with
 *    no flipped triangle. A chart is scaled so its area in the plane equals
 *    its area on the surface, so every chart has one texel density, mirrored
 *    if needed so its triangles keep their winding, and turned to its
 *    principal axes.
 * 3. **Packing.** The charts' bounding boxes go onto shelves, tallest first,
 *    `mPadding` times the square root of the surface area apart, and the
 *    whole atlas is scaled uniformly into the unit square.
 *
 * A point on a seam belongs to several charts and is split: the first copy
 * keeps its index (the chart of the first triangle that uses it) and the
 * others are appended, `point_data["parameterize:parent_point"]` naming the
 * input point of each, so point regions, lines and vertex cells stay valid.
 * Point data is copied to the copies, cell data and regions ride along, and
 * `cell_data["parameterize:chart"]` holds each triangle's chart (-1 on other
 * blocks). A point in no triangle gets `NaN` coordinates.
 *
 * **Scope.** The surface is the mesh's `triangle` blocks; `line` and `vertex`
 * blocks ride along. A volume cell, a quad, a polygon or a higher-order
 * surface cell is refused naming the operation that makes triangles.
 * Everything is serial and in input order, so the result is the same on every
 * backend. There is no numpy twin: the conformal solve is iterative.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data: Float64 `(n, 2)`, the texture coordinates in `[0, 1]`.
inline constexpr const char* kUvName = "uv";
/// Point data: Int64 `(n,)`, the input point each output point copies.
inline constexpr const char* kParameterizeParentPointName = "parameterize:parent_point";
/// Cell data: Int64, the chart of each triangle (-1 on other blocks).
inline constexpr const char* kParameterizeChartName = "parameterize:chart";

/// How `parameterize` cuts and packs.
struct ParameterizeOptions {
    /// The gap between charts, as a fraction of the square root of the
    /// surface's area (before the atlas is scaled into the unit square).
    double mPadding = 0.02;
    /// How many times a piece that does not flatten may be halved.
    std::int64_t mMaxDepth = 8;
    /// A chart whose L2 stretch exceeds this is halved again (while the depth
    /// allows), even though nothing in it folds; infinity keeps every chart
    /// that does not fold.
    double mMaxStretch = 2.0;
};

/// What `parameterize` did.
struct ParameterizeResult {
    /// The surface, split at the seams, with `uv`, the parent points and the charts.
    Mesh mMesh;
    std::int64_t mNumCharts = 0;
    /// Triangles whose orientation in the plane is not their chart's (or of no
    /// area): 0 unless a piece still folded at the depth cap.
    std::int64_t mNumFlipped = 0;
    /// Sander et al.'s area-weighted L2 stretch of the charts at their own
    /// scale (1 for a map that keeps lengths, larger with distortion).
    double mStretch = 0.0;
    /// The fraction of the unit square the charts' triangles cover.
    double mPackingEfficiency = 0.0;
    /// Points added on seams.
    std::int64_t mNumSeamPoints = 0;
};

/**
 * @brief A UV atlas of @p rMesh's triangles (see the file comment).
 * @throws std::invalid_argument on a cell type outside the scope, no triangle,
 *         a triangle naming a point outside the mesh, a negative `mPadding` or
 *         `mMaxDepth`, or an `mMaxStretch` below 1.
 */
MESHIOPLUSPLUS_API ParameterizeResult parameterize(const Mesh& rMesh,
                                                   const ParameterizeOptions& rOptions);

}  // namespace meshioplusplus
