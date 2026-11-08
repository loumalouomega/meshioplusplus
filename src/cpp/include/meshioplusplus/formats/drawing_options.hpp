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
 * @file formats/drawing_options.hpp
 * @brief The camera and colouring options the two drawing writers share
 * (`write_svg`, `write_tikz`): one spelling, embedded in each writer's own
 * options struct (`SvgWriteOptions`, `TikzWriteOptions`).
 */

// System includes
#include <optional>
#include <string>

namespace meshioplusplus {

/**
 * @brief The orthographic camera of a drawing writer, used for genuinely 3D
 * input only (flat input keeps the classic 2D path). The default pair
 * (45, atan(1/sqrt(2))) is the classic CAD isometric view.
 */
struct DrawingView {
    double mAzimuth = 45.0;                  ///< camera azimuth, degrees
    double mElevation = 35.264389682754654;  ///< camera elevation, degrees
    double mRoll = 0.0;                      ///< in-screen camera roll, degrees
};

/**
 * @brief Colouring a drawing's faces by a data array.
 *
 * Empty `mColorBy` (the default) keeps the writer's flat fill and leaves the
 * output byte-identical to the uncoloured drawing. Point data colours a face by
 * the mean of its corner values, cell data by its owning cell's value -- found
 * through `"surface:parent_cell"` for a projected volume mesh, whose skin is
 * extracted with provenance in that case.
 */
struct DrawingColoring {
    /// name of a `point_data` or `cell_data` array to colour the faces by
    std::string mColorBy;
    /// component of a multi-component array; `std::nullopt` uses the row magnitude
    std::optional<int> mComponent;
    /// built-in colormap name (see `detail/colormap.hpp`)
    std::string mCmap = "viridis";
    /// low end of the mapped range; `std::nullopt` uses the smallest finite value
    /// among the *drawn* faces
    std::optional<double> mVMin;
    /// high end of the mapped range; `std::nullopt` likewise
    std::optional<double> mVMax;
    /// fill for faces whose value is NaN or infinite, in the writer's colour
    /// syntax; each writer's options give its own default
    std::string mNanColor;
    /// append a gradient bar with min/max labels to the right of the figure
    bool mColorbar = false;
};

}  // namespace meshioplusplus
