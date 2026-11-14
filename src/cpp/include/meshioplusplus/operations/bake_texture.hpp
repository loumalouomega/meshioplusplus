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
 * @file operations/bake_texture.hpp
 * @brief Bake a point or cell array into a texture image through the mesh's `uv`
 * (roadmap 2.4.8), so a field travels as a picture in glTF and OBJ.
 *
 * Each 2-D cell (a triangle; a quad split along its 0-2 diagonal; a polygon fanned
 * from its first node) is drawn in texture space, `x = u * width` and
 * `y = (1 - v) * height` (the image's top row is `v = 1`). A pixel whose centre lies
 * in a triangle (all three barycentric weights non-negative; a later cell draws over
 * an earlier one) takes the array's value there: the weighted value of a point
 * array, the cell's own value of a cell array. A multi-component array reduces to
 * `mComponent`, else to its magnitude. The value maps through the colormap
 * (`detail/colormap.hpp`) over `[mVMin, mVMax]` (the finite range by default);
 * a pixel no cell covers, or whose value is not finite, stays transparent.
 *
 * The result is the RGBA image and the mesh with it as `field_data["texture:png"]`
 * (PNG, zlib level 9 when the build has zlib, stored otherwise). The Python twin is
 * `meshioplusplus.bake_texture`, which draws the same pixels.
 * @see doc/bake_texture.md
 */

// System includes
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief What to bake and how. */
struct BakeTextureOptions {
    std::string mArray;          ///< a point or cell data array
    std::int64_t mWidth = 512;   ///< image width in pixels
    std::int64_t mHeight = 512;  ///< image height in pixels
    std::string mCmap = "viridis";
    std::optional<double> mVMin;   ///< low end of the colour range; the finite minimum when unset
    std::optional<double> mVMax;   ///< high end; the finite maximum when unset
    std::int64_t mComponent = -1;  ///< the component of a vector array; -1 is the magnitude
};

/** @brief The baked image and the mesh carrying it. */
struct BakeTextureResult {
    Mesh mMesh;          ///< the input with `field_data["texture:png"]`
    NDArray mImage;      ///< `(height, width, 4)` uint8, top row first
    double mVMin = 0.0;  ///< the range used
    double mVMax = 0.0;
    std::int64_t mCoveredPixels = 0;  ///< pixels some cell drew
    std::vector<std::string> mWarnings;
};

/**
 * @brief Bake @p rOptions.mArray into a texture through `point_data["uv"]`.
 * @throws std::invalid_argument when the mesh has no `uv`, the array is missing or
 *         has the wrong length, the size is not positive, or the colormap is unknown.
 */
MESHIOPLUSPLUS_API BakeTextureResult bake_texture(const Mesh& rMesh,
                                                  const BakeTextureOptions& rOptions);

}  // namespace meshioplusplus
