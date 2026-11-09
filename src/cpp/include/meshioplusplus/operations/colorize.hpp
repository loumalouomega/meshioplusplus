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
 * @file operations/colorize.hpp
 * @brief Bake a scalar field into 8-bit RGB colours through a colormap
 * (roadmap 7.4.4), so a deviation map (or any field) travels as the vertex
 * or face colours PLY viewers, MeshLab and glTF (`color_by="rgb"`) read.
 *
 * The field's values `v` map to `t = (v - vmin) / (vmax - vmin)`, clamped to
 * `[0, 1]` (`detail::color_param`: a degenerate range maps to the middle),
 * then through the colormap's 256-entry table (`detail::colormap_lookup`).
 * With `mSteps` of 2 or more the colours are **banded**: `t` is replaced by
 * the centre of the equal band it falls in, so the map reads as `mSteps`
 * flat colours at fixed steps (the last band holds `t = 1`). A non-finite
 * value takes `mNanColor`.
 *
 * The colours are written as three UInt8 arrays named `red`, `green` and
 * `blue` beside the field, on its location (point data for a point field,
 * cell data for a cell field): the names PLY's vertex and face colour
 * properties take, which the PLY writer writes as `uchar` properties under
 * those names. Points, cells and every other array are unchanged.
 *
 * @see doc/colorize.md
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The UInt8 colour channels `colorize` writes.
inline constexpr const char* kColorRedName = "red";
inline constexpr const char* kColorGreenName = "green";
inline constexpr const char* kColorBlueName = "blue";

/// Options of `colorize`.
struct ColorizeOptions {
    /// The point or cell data to colour by (required). A vector field is
    /// reduced to its magnitude unless `mComponent` picks one entry.
    std::string mField;
    /// A colormap name (`detail::colormap_names()`: viridis, coolwarm, turbo, ...).
    std::string mCmap = "viridis";
    /// The range; when `mVMin >= mVMax` (the default) the field's finite
    /// minimum and maximum (a symmetric `+-max |v|` with `mSymmetric`).
    double mVMin = 0.0;
    double mVMax = 0.0;
    /// With the automatic range: `+-max |v|`, so zero sits mid-map (a
    /// deviation map's natural range).
    bool mSymmetric = false;
    /// Bands: 0 or 1 is continuous, 2 or more bands at fixed steps.
    std::int64_t mSteps = 0;
    /// The colour of a non-finite value.
    std::string mNanColor = "#808080";
    /// Colour by this component of a vector field (-1: the magnitude).
    std::int64_t mComponent = -1;
};

/// The coloured mesh and the range used.
struct ColorizeResult {
    Mesh mMesh;
    double mVMin = 0.0;
    double mVMax = 0.0;
    /// "point" or "cell": where the colours were written.
    std::string mLocation;
};

/**
 * @brief Bake @p rOptions.mField into `red`, `green` and `blue` UInt8 data.
 * @throws std::invalid_argument on a missing field, an unknown colormap, a
 *         bad colour, component or step count.
 */
MESHIOPLUSPLUS_API ColorizeResult colorize(const Mesh& rMesh, const ColorizeOptions& rOptions);

}  // namespace meshioplusplus
