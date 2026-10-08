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
 * @file detail/raster.hpp
 * @brief The tile rasterizer behind `render`: triangles, lines and discs into
 * a colour, depth and id buffer.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * ### The determinism contract
 *
 * Every pixel's result is a pure function of the primitive lists, whatever the
 * parallel backend and thread count:
 *
 *  - Screen positions are rounded once to fixed point with `kRasterSubBits`
 *    (8) fractional bits and clamped to `kRasterLimitPx` pixels, so every edge
 *    function is an exact 64-bit integer product.
 *  - A pixel is covered when its centre lies inside all three edges, with the
 *    top-left rule for a centre exactly on an edge: a pixel on the shared edge
 *    of two triangles belongs to exactly one of them.
 *  - Depth and shading are interpolated in double precision from those same
 *    integers, in a fixed expression order (`-ffp-contract=off` is set for the
 *    whole core, so no fused multiply-add can change a rounding).
 *  - The depth test is strict: a primitive replaces a pixel only when it is
 *    closer, so on a tie the one drawn first keeps it. Triangles are drawn in
 *    list order, then lines, then points.
 *  - The frame is cut into `kRasterTile`-pixel square tiles; one worker owns a
 *    tile, and each tile's primitive lists are built by one serial pass in list
 *    order. No two workers write the same pixel.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/operations/render.hpp"

namespace meshioplusplus {
namespace detail {

inline constexpr int kRasterSubBits = 8;
inline constexpr int kRasterTile = 32;
inline constexpr double kRasterLimitPx = 262144.0;

/// One screen-space vertex: pixel coordinates (y down, pixel centres at
/// `+0.5`) and a closeness key -- larger is nearer the camera, and the key is
/// linear in screen space for either camera.
struct RasterVertex {
    double mX = 0.0;
    double mY = 0.0;
    double mDepth = 0.0;
};

/// A triangle with a base colour scaled by an intensity interpolated from its
/// three corners (equal for flat shading).
struct RasterTri {
    std::int64_t mV[3] = {0, 0, 0};
    double mIntensity[3] = {1.0, 1.0, 1.0};
    std::array<std::uint8_t, 4> mColor = {0, 0, 0, 255};
    std::int64_t mId = -1;
};

/// A line segment `Width` pixels wide whose depth is biased toward the camera.
struct RasterLine {
    std::int64_t mA = 0;
    std::int64_t mB = 0;
    std::array<std::uint8_t, 4> mColor = {0, 0, 0, 255};
    std::int64_t mId = -1;
};

/// A disc around one vertex.
struct RasterPoint {
    std::int64_t mV = 0;
    std::array<std::uint8_t, 4> mColor = {0, 0, 0, 255};
    std::int64_t mId = -1;
};

/// What a raster pass draws.
struct RasterScene {
    std::vector<RasterVertex> mVertices;
    std::vector<RasterTri> mTris;
    std::vector<RasterLine> mLines;
    std::vector<RasterPoint> mPoints;
    /// Line width in pixels (a positive integer: the supersampling factor).
    int mLineWidth = 1;
    /// Disc radius in pixels.
    double mPointRadius = 1.0;
    /// Added to a line's or point's depth before its test, so an edge drawn on
    /// its own face wins.
    double mDepthBias = 0.0;
};

/// The buffers a raster pass fills: `mWidth * mHeight` pixels of RGBA, depth
/// and id. Pixels no primitive covers keep `mBackground` and id -1.
struct RasterTarget {
    int mWidth = 0;
    int mHeight = 0;
    std::vector<std::uint8_t> mRgba;
    std::vector<float> mDepth;
    std::vector<std::int64_t> mIds;
};

/**
 * @brief Draw @p rScene into a new target of the given size.
 * @param Background the colour uncovered pixels keep
 */
RasterTarget rasterize(const RasterScene& rScene, int Width, int Height,
                       const std::array<std::uint8_t, 4>& rBackground);

/**
 * @brief The fixed-order box filter: average each `Factor` x `Factor` block of
 * @p rTarget into one pixel, weighting colour by alpha, with integer rounding.
 * The id of an output pixel is the id at its block's `(Factor/2, Factor/2)`
 * sample.
 */
void raster_downsample(const RasterTarget& rTarget, int Factor, std::vector<std::uint8_t>& rRgba,
                       std::vector<std::int64_t>& rIds);

/**
 * @brief Resolve `RenderOptions::mView` into an azimuth and elevation, or copy
 * the options' own when it is empty. Defined in `operations/render.cpp`.
 * @throws std::invalid_argument for an unknown view name
 */
void render_view_angles(const RenderOptions& rOpt, double& rAzimuth, double& rElevation);

}  // namespace detail
}  // namespace meshioplusplus
