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
 * @file detail/cut_graph.hpp
 * @brief The topology of a triangle patch (Euler characteristic, boundary
 * loops, genus) and the cuts that turn a genus-0 patch that is not a disc
 * into discs, for `fit_cad`'s freeform modes (v17.12.0).
 *
 * A freeform face needs a parametrisation, and a parametrisation needs a
 * disc. `fit_cad`'s legacy stages halve such a patch geometrically until the
 * pieces happen to be discs; here the cut follows the topology instead:
 *
 *  - an annulus (two boundary loops) is cut by two shortest paths between its
 *    loops, the second starting half way round the first loop from the
 *    first, into two discs;
 *  - a sphere with `b >= 3` holes is cut by a cycle of shortest paths
 *    `L0 -> L1 -> ... -> L(b-1) -> L0` through every loop, which separates
 *    it into two discs (each path leaves a loop half way round from where
 *    the previous one arrived);
 *  - a closed genus-0 patch (no boundary) is cut along the level set of the
 *    geodesic distance difference to its two most distant vertices (the
 *    area-weighted median), into two caps.
 *
 * A patch of genus 1 or more is not cut here: `fit_cad` leaves it to the quad
 * layout of `network` mode, or to the legacy halving. Shortest paths are
 * Dijkstra's on edge lengths with ties broken by vertex index; pieces are
 * numbered in the order of their lowest triangle. Serial and deterministic.
 *
 * A **core-private** header (the `patch_flatten.hpp` precedent): nothing
 * installed names it.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "patch_flatten.hpp"

namespace meshioplusplus {
namespace detail {

/// What `patch_topology` measured.
struct PatchTopology {
    std::size_t mNumVertices = 0;  ///< vertices the triangles use
    std::size_t mNumEdges = 0;
    std::size_t mNumTriangles = 0;
    std::size_t mNumBoundaryLoops = 0;
    /// `V - E + F`.
    std::int64_t mEuler = 0;
    /// `(2 - chi - b) / 2`; meaningful only for a manifold patch.
    std::int64_t mGenus = 0;
    /// Every edge borders one or two triangles, and every boundary vertex
    /// starts one boundary edge.
    bool mManifold = true;
    bool IsDisc() const { return mManifold && mNumBoundaryLoops == 1 && mEuler == 1; }
};

/// The topology of a triangle patch (`rTris` index any points).
PatchTopology patch_topology(const std::vector<PatchTri>& rTris);

/**
 * @brief The shortest edge path, by Euclidean edge length, from any vertex in
 * `rSources` to any vertex in `rTargets`, through no vertex in `rForbidden`
 * (sources and targets themselves excepted).
 * @return the vertices from source to target; empty when none connects them.
 */
std::vector<std::size_t> patch_shortest_path(const std::vector<PatchVec>& rPts,
                                             const std::vector<PatchTri>& rTris,
                                             const std::vector<std::size_t>& rSources,
                                             const std::vector<std::size_t>& rTargets,
                                             const std::vector<std::uint8_t>& rForbidden);

/**
 * @brief Split a genus-0 patch that is not a disc into discs (see the file
 * comment).
 * @param MinTriangles a piece with fewer triangles fails the split.
 * @param rNumPieces the pieces, numbered from 0; 0 when the patch is a disc
 *        already, has genus above 0, is not manifold, or no cut gives discs
 *        of at least `MinTriangles` triangles.
 * @return per triangle its piece; empty when `rNumPieces` is 0.
 */
std::vector<std::int64_t> patch_split_to_discs(const std::vector<PatchVec>& rPts,
                                               const std::vector<PatchTri>& rTris,
                                               std::size_t MinTriangles, std::size_t& rNumPieces);

}  // namespace detail
}  // namespace meshioplusplus
