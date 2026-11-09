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
 * @file detail/polyline_chain.hpp
 * @brief Segments joined into polylines (roadmap 7.6.1): the chains of a
 * line mesh, each walked end to end in a canonical order, so a section is
 * a list of curves rather than a bag of segments.
 *
 * A point's degree is the number of segments touching it. A chain runs
 * between two points whose degree is not 2 (an end, or a branch point) or,
 * with every degree 2, round a closed loop. Canonical order, so both engines
 * and every run agree: open chains first, each started at the lowest-id
 * point whose degree is not 2 that still has an unwalked segment, then
 * loops, each started at its lowest point id; a walk always takes the
 * unwalked segment to the lowest-id neighbour, the segment of lowest index
 * on a tie. Degenerate segments (both ends one point) are left out.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// One chain: its points in walking order (a loop's last point is its first),
/// its segments in that order and whether each was walked against its own
/// orientation.
struct PolylineChain {
    std::vector<std::int64_t> mPoints;
    std::vector<std::int64_t> mSegments;
    std::vector<std::uint8_t> mReversed;
    bool mClosed = false;
};

/**
 * @brief The chains of @p rSegments (point-id pairs) over @p NumPoints
 * points, in the canonical order of the file comment.
 */
std::vector<PolylineChain> polyline_chains(
    const std::vector<std::array<std::int64_t, 2>>& rSegments, std::size_t NumPoints);

}  // namespace detail
}  // namespace meshioplusplus
