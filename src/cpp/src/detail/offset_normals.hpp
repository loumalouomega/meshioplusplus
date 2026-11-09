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
 * @file detail/offset_normals.hpp
 * @brief Each point's offset direction over a triangle and quad surface: the
 * most normal normal of Aubry and Löhner (2009), optionally smoothed, and the
 * stretch that keeps every face around a point at the offset distance.
 *
 * A **core-private** header: no installed header names it. `boundary_layer`
 * grows its layers along these directions and `offset_surface`, `thicken` and
 * `hollow` move points along them (roadmap 6.9, 7.8.4); one function computes
 * both, so they agree bit for bit.
 *
 * For each point, the angle-weighted mean of its faces' unit normals is
 * stepped towards the face that sees it least (100 steps, the best kept);
 * `mStretch` is one over the least cosine between the direction and a face
 * normal around the point, at most `kOffsetMaxStretch`, so moving a point by
 * `distance * stretch` keeps each face's plane `distance` away (a cube's
 * corner moves by `sqrt(3) * distance` along its diagonal). A point that no
 * direction sees from every face within `kOffsetVisible` is `mNotVisible`,
 * its stretch 0. Smoothing passes (Jacobi, so order-free) move a direction
 * towards its neighbours' mean while every face still sees it at more than
 * half its first least cosine; the stretch is then recomputed.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace meshioplusplus::detail {

/// The least cosine at which a point still sees all its faces.
inline constexpr double kOffsetVisible = 0.05;
/// Offsets are divided by the least cosine, at most this many times.
inline constexpr double kOffsetMaxStretch = 4.0;

struct OffsetNormals {
    std::vector<std::array<double, 3>> mFaceNormals;  ///< unit, times the sign
    /// Each point's faces, with the corner angle there.
    std::vector<std::vector<std::pair<std::size_t, double>>> mAround;
    std::vector<std::vector<std::int64_t>> mNeighbours;  ///< sorted, unique
    std::vector<std::array<double, 3>> mDirections;
    std::vector<double> mStretch;  ///< 0 where not visible or unused
    std::vector<double> mSeen;     ///< the first least cosine
    std::vector<char> mUsed;
    std::vector<char> mNotVisible;
};

/**
 * @brief The offset directions of the points @p rX (xyz, packed) over the
 * faces @p rFaces (a triangle's fourth entry is -1); @p Sign is 1 or -1 (the
 * faces' own normals, or the opposite side).
 */
OffsetNormals offset_normals(const std::vector<double>& rX,
                             const std::vector<std::array<std::int64_t, 4>>& rFaces, double Sign,
                             int SmoothIterations);

}  // namespace meshioplusplus::detail
