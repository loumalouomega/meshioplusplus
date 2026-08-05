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
 * @file detail/surface_distance.hpp
 * @brief The signed-distance kernel: a triangle soup, its pseudonormal tables,
 * a bucket-grid accelerator, and the nearest-point query over them.
 *
 * `operations/sdf.cpp` is the driver; everything geometric lives here so that it
 * can be tested against closed-form distance fields with no grid, no mesh
 * generation and no binding in the way.
 *
 * ### The accelerator is deliberately not a BVH
 *
 * `detail/spatial_hash.hpp` already provides `InsertBox` (which exists precisely
 * to register a simplex's quantized bounding box) and `ForEachInShell` with a
 * documented stopping rule, and `interpolate.cpp` already runs an expanding-shell
 * nearest search over it. The workload here -- many queries against a
 * roughly uniform triangle soup -- is the uniform grid's best case and the BVH's
 * worst.
 *
 * The decisive reason, though, is not speed: **a BVH would make the accelerator
 * observable in the output.** Its build order and split heuristic determine the
 * order candidates are visited in, which determines which of two equidistant
 * triangles wins. The grid can be made provably unobservable instead, by giving
 * every comparison a total order -- `(distance squared, triangle id)` -- so that
 * the answer cannot depend on visit order at all. `SurfaceDistanceOptions`
 * exposes `mGridCellSize` mainly so a test can *prove* that: the same query at
 * several bucket sizes must return byte-identical distances and nearest cells.
 *
 * ### The sign
 *
 * `closest_point_on_triangle` reports which feature of the triangle the nearest
 * point lies on, and that selects the normal:
 *
 *  - face   -> the triangle's own normal;
 *  - edge   -> the sum of the two incident triangles' unit normals;
 *  - vertex -> the angle-weighted sum of the incident triangles' unit normals.
 *
 * Using the nearest *triangle's* normal in all three cases is the classic bug:
 * exactly right on convex geometry, exactly wrong on the concave side of every
 * crease, and it renders beautifully either way.
 *
 * Building those tables is the one place with a floating-point ordering hazard,
 * because summing unit normals in a non-deterministic order changes the last
 * bits and a last-bit change can flip a sign at a near-tangent query point. The
 * tables are therefore accumulated by a **serial** pass in ascending
 * (triangle, corner) order -- the phase-split idiom `surface.cpp` uses.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief A surface reduced to triangles, with each triangle remembering which
 * input cell it came from.
 *
 * The provenance matters because `sdf:closest_cell` must name a cell of the
 * mesh the caller passed in, not of a triangulation they never saw.
 */
struct TriangleSoup {
    /// Corner coordinates, three per triangle, in `(a, b, c)` order.
    std::vector<Vec3> mCorners;
    /// Per triangle, the global (block-major) index of the input cell it came from.
    std::vector<std::int64_t> mSourceCell;
    /// Per triangle, the three vertex ids in the *welded* numbering below.
    std::vector<std::array<std::int64_t, 3>> mVertices;
    /// Distinct vertex positions, indexed by the ids in `mVertices`.
    std::vector<Vec3> mPoints;

    std::size_t NumTriangles() const { return mSourceCell.size(); }
};

/**
 * @brief Reduce a surface mesh to a triangle soup.
 * @param rSurface the mesh; `triangle` blocks are taken as they are, `quad` and
 *        rectangular `polygon` blocks are fanned exactly as
 *        `convert_cells(Simplexify)` fans them so the two cannot disagree.
 * @param rRegion restrict to this named `Cell` region; empty takes everything.
 * @throws std::invalid_argument on a 3-D or polyhedron block (pointing at
 *         `extract_surface`), on a higher-order block (pointing at `linearize`),
 *         and on an unknown region name.
 */
MESHIOPLUSPLUS_API TriangleSoup build_triangle_soup(const Mesh& rSurface,
                                                    const std::string& rRegion);

/// The four edge defect counts of a soup, and the resulting verdict.
MESHIOPLUSPLUS_API SurfaceQuality soup_quality(const TriangleSoup& rSoup);

/// An undirected edge, as the sorted pair of its endpoints' vertex ids.
using SurfaceEdgeKey = std::array<std::int64_t, 2>;

/// Hash for `SurfaceEdgeKey`, the `GridKeyHash` mixing constant and shape.
struct SurfaceEdgeKeyHash {
    std::size_t operator()(const SurfaceEdgeKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

/**
 * @brief A soup prepared for querying: the accelerator plus the normal tables.
 *
 * Held by value and moved; building one is O(triangles) and is paid once per
 * call rather than once per query, which is why the public API is a batch one.
 */
struct DistanceQuery {
    const TriangleSoup* mpSoup = nullptr;
    /// The bucket grid, holding triangle ids by quantized bounding box.
    SpatialGrid mGrid{1.0};
    /// Per triangle, its unnormalized normal (`cross(ab, ac)`).
    std::vector<Vec3> mFaceNormal;
    /// Per welded vertex, the weighted sum of incident unit face normals.
    std::vector<Vec3> mVertexNormal;
    /// Per edge (sorted endpoint ids), the sum of incident unit face normals.
    std::unordered_map<SurfaceEdgeKey, Vec3, SurfaceEdgeKeyHash> mEdgeNormal;
    /// The bucket size actually used, after the auto rule or the caller's override.
    double mCellSize = 1.0;
};

/// Build the accelerator and normal tables for @p rSoup.
MESHIOPLUSPLUS_API DistanceQuery build_distance_query(const TriangleSoup& rSoup,
                                                      const SurfaceDistanceOptions& rOptions);

/// What a single query point resolved to.
struct DistanceHit {
    double mSignedDistance = 0.0;   ///< Negative inside, per the usual convention.
    std::int64_t mSourceCell = -1;  ///< Nearest input cell, or -1 outside the band.
    bool mInBand = true;            ///< False when the value was clamped to the band.
};

/**
 * @brief Distances from @p rPoints to the soup @p rQuery was built from.
 *
 * Parallel over query points, which are independent; each point's own search is
 * serial and totally ordered, so the result does not depend on thread count.
 */
MESHIOPLUSPLUS_API std::vector<DistanceHit> query_distances(
    const DistanceQuery& rQuery, const std::vector<Vec3>& rPoints,
    const SurfaceDistanceOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus
