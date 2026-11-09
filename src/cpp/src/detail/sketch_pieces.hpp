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
 * @file detail/sketch_pieces.hpp
 * @brief A `fit_curves` sketch read back from a mesh's `sketch:*` field data
 * as pieces in world coordinates (roadmap 7.6.4, 7.6.5), for the writers that
 * draw it: DXF and STEP.
 *
 * A piece's 2-D start, end and centre (`sketch:params`) are mapped through
 * its frame (`sketch:frames`, row-major 4x4, columns e1, e2, normal, origin)
 * as `origin + u e1 + v e2`, in that order.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

struct SketchWorldPiece {
    std::int64_t mType = 0;  ///< 1 line, 2 arc, 3 circle (`SketchPieceType`)
    std::array<double, 3> mStart{}, mEnd{}, mCentre{};
    std::array<double, 3> mE1{}, mE2{}, mNormal{}, mOrigin{};  ///< its frame
    double mRadius = 0.0;
    double mSweepDeg = 0.0;  ///< counter-clockwise about the normal positive
};

/// Whether @p rMesh carries the four `sketch:*` field data arrays.
bool has_sketch(const Mesh& rMesh);

/**
 * @brief The sketch's pieces in world coordinates.
 * @param pPrefix the caller's error prefix (`"meshio++: dxf: "`).
 * @throws WriteError when the arrays disagree or a piece names a missing frame.
 */
std::vector<SketchWorldPiece> sketch_world_pieces(const Mesh& rMesh, const char* pPrefix);

}  // namespace detail
}  // namespace meshioplusplus
