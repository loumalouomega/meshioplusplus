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
 * @file operations/offset_surface.hpp
 * @brief Offset a surface, thicken it into a shell, hollow a solid (roadmap
 * 7.8.4).
 *
 * Each point of a triangle and quad surface moves along its offset direction
 * (`detail/offset_normals.hpp`, the directions `boundary_layer` grows along:
 * the most normal normal of Aubry and Löhner), by the distance times its
 * stretch, one over the least cosine between the direction and a face around
 * it, at most 4: both faces of a sharp edge keep the distance (a cube offset
 * by `d` is the cube grown by `d`, its corners moved by `sqrt(3) * d`). A
 * positive distance moves along the faces' normals (their winding). A point
 * no direction sees from all its faces within a cosine of 0.05 moves by the
 * distance alone and is counted; a face whose normal turns over is counted as
 * folded, with a warning, as the offset of a concave region deeper than its
 * radius folds.
 *
 * `offset_surface` keeps the connectivity, the data and the regions, and moves
 * the points. `thicken` turns a surface into a closed shell: the surface
 * offset to the two sides of the wall (`Outward`: the surface and its offset
 * by `+thickness`; `Inward`: by `-thickness`; `Both`: by `+-thickness / 2`),
 * the lower one reversed, and each boundary edge joined across by a quad; the
 * points are the upper layer's then the lower's (the point data duplicated),
 * and the Int64 cell data `thicken:side` is 0 on the upper layer, 1 on the
 * lower and 2 on the rim. `hollow` gives a closed surface an inner wall: the
 * surface and its offset by `-wall`, reversed, with no rim (an open surface is
 * refused). Their cells are one `triangle` and one `quad` block; the input's
 * cell data and regions are not carried (a warning names them).
 */

// System includes
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

struct OffsetSurfaceOptions {
    double mDistance = 0.0;     ///< along the faces' normals; negative: against them
    int mSmoothIterations = 0;  ///< smoothing passes over the directions
};

struct OffsetSurfaceResult {
    Mesh mMesh;
    std::int64_t mNumNotVisible = 0;  ///< points moved by the distance alone
    std::int64_t mNumClamped = 0;     ///< points whose stretch hit its bound
    std::int64_t mNumFolded = 0;      ///< faces whose normal turned over
    double mMaxMove = 0.0;            ///< the longest point move
    std::int64_t mNumRimQuads = 0;    ///< thicken: the quads joining the layers
    std::vector<std::string> mWarnings;
};

enum class ThickenSide : std::int32_t { Inward = 0, Outward = 1, Both = 2 };

struct ThickenOptions {
    double mThickness = 0.0;
    ThickenSide mSide = ThickenSide::Inward;
    int mSmoothIterations = 0;
};

/// @throws std::invalid_argument naming the accepted sides.
MESHIOPLUSPLUS_API ThickenSide thicken_side_from_name(std::string_view Name);

/**
 * @brief Move each point of @p rMesh along its offset direction by the
 * distance times its stretch (see the file comment).
 * @throws std::invalid_argument on a non-finite distance or a cell that is not
 * a triangle or a quad.
 */
MESHIOPLUSPLUS_API OffsetSurfaceResult offset_surface(const Mesh& rMesh,
                                                      const OffsetSurfaceOptions& rOptions = {});

/**
 * @brief A closed shell of the given thickness from a surface (see the file
 * comment).
 * @throws std::invalid_argument on a thickness that is not positive and finite.
 */
MESHIOPLUSPLUS_API OffsetSurfaceResult thicken(const Mesh& rMesh, const ThickenOptions& rOptions);

/**
 * @brief A closed surface with an inner, reversed wall @p Wall inside it.
 * @throws std::invalid_argument on a wall that is not positive and finite or
 * a surface with boundary edges.
 */
MESHIOPLUSPLUS_API OffsetSurfaceResult hollow(const Mesh& rMesh, double Wall,
                                              int SmoothIterations = 0);

}  // namespace meshioplusplus
