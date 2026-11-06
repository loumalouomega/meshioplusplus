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
 * @file detail/tet_checks.hpp
 * @brief The volume tolerance and the face-overlap count shared by
 * `check_tetrahedra` and the volume fix of `octree_mesh`.
 *
 * Private (not installed): a change here is not an ABI event.
 */

// System includes
#include <array>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus::detail {

/// A tetrahedron as four point indices.
using TetNodes = std::array<std::int64_t, 4>;

/// The signed volume of the tetrahedron @p rTet over @p rPoints, and the
/// tolerance under which it counts as zero: @p ZeroVolumeFactor times the cube
/// of its longest edge, plus the effect of the six significant digits of
/// written coordinates.
std::pair<double, double> tet_volume_and_tolerance(const std::vector<Vec3>& rPoints,
                                                   const TetNodes& rTet, double ZeroVolumeFactor);

/// The pairs of tetrahedra that share a face (and only two use it) whose
/// opposite nodes lie on the same side of it. A tetrahedron with a non-zero
/// entry in @p pSkip (one per tetrahedron, or NULL) is left out.
std::int64_t count_tet_overlaps(const std::vector<Vec3>& rPoints,
                                const std::vector<TetNodes>& rTets,
                                const std::vector<std::uint8_t>* pSkip = nullptr);

}  // namespace meshioplusplus::detail
