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
 * @file operations/hull.hpp
 * @brief The convex hull and the planar Delaunay triangulation of a point set
 * (roadmap 7.5.1): scattered points from an `xyz`, `pcd` or `ply` scan turned
 * into a footprint, a closed hull, or a terrain.
 *
 * Both read only the mesh's points (its cells are ignored), and every
 * decision is one of the exact predicates of `detail/predicates.hpp`, so the
 * result depends on the coordinates alone and is the same on every platform.
 *
 * - `convex_hull(mesh, {dimension})`: in 3-D (quickhull, C. B. Barber, D. P.
 *   Dobkin and H. Huhdanpaa, ACM TOMS 22(4), 1996, with exact `orient3d`
 *   visibility), a closed `triangle` surface wound outward; in 2-D (Andrew's
 *   monotone chain over the `x, y` coordinates with exact `orient2d`), a
 *   counter-clockwise closed `line` loop. Only the extreme points are
 *   vertices: a point on a facet or an edge is not. The output holds the hull
 *   vertices alone, in ascending input order, with their point data gathered
 *   and `hull:point_id` (Int64) naming the input point each came from; field
 *   data rides along. `dimension` 0 picks 2-D for 2-D points and for 3-D
 *   points whose `z` are all equal, else 3-D. Coplanar 3-D points and
 *   collinear 2-D ones span no hull and are refused by name.
 * - `delaunay_2d(mesh)`: the Delaunay triangulation of the points' `x, y`
 *   (the constrained triangulation of `detail/cdt2d.hpp` with the convex hull's
 *   edges as its only constraints, which are Delaunay edges, so the result is
 *   the Delaunay triangulation; among cocircular points the choice is
 *   deterministic). The output keeps every input point, its coordinates
 *   (a terrain keeps its `z`) and its point data, with one counter-clockwise
 *   `triangle` block; a point that coincides in `x, y` with an earlier one is
 *   in no triangle (`mNumDuplicates`). Fewer than three distinct points, or
 *   collinear ones, are refused by name.
 *
 * Coordinates must be finite; `delaunay_2d` also needs each within the
 * exact range of `detail::cdt2d_in_range` (zero, or a magnitude in
 * [2^-128, 2^128]). Sequential and deterministic.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data on a hull: Int64, the input point each hull vertex is.
inline constexpr const char* kHullPointIdName = "hull:point_id";

/// What `convex_hull` builds.
struct ConvexHullOptions {
    /// 2 (the hull of the `x, y` coordinates), 3, or 0 to decide from the
    /// points (see the file comment).
    std::int32_t mDimension = 0;
};

/// What `convex_hull` built.
struct ConvexHullResult {
    /// The hull: a `triangle` surface (3-D) or a `line` loop (2-D).
    Mesh mMesh;
    /// The dimension used, 2 or 3.
    std::int32_t mDimension = 0;
    /// Hull vertices (extreme points).
    std::int64_t mNumVertices = 0;
    /// Triangles (3-D) or edges (2-D).
    std::int64_t mNumFacets = 0;
    /// 2-D: the enclosed area; 3-D: the surface area.
    double mArea = 0.0;
    /// 3-D: the enclosed volume; 2-D: 0.
    double mVolume = 0.0;
    /// 2-D: the perimeter; 3-D: 0.
    double mPerimeter = 0.0;
};

/**
 * @brief The convex hull of a mesh's points (see the file comment).
 * @throws std::invalid_argument on a non-finite coordinate, a dimension other
 *         than 0, 2 or 3, 1-D points, too few distinct points, or points that
 *         span no hull (collinear in 2-D, coplanar in 3-D).
 */
MESHIOPLUSPLUS_API ConvexHullResult convex_hull(const Mesh& rMesh,
                                                const ConvexHullOptions& rOptions = {});

/// What `delaunay_2d` built.
struct Delaunay2dResult {
    /// Every input point, one `triangle` block.
    Mesh mMesh;
    std::int64_t mNumTriangles = 0;
    /// Points coinciding in `x, y` with an earlier point (in no triangle).
    std::int64_t mNumDuplicates = 0;
    /// Points on the convex hull's boundary (its vertices and the points
    /// on its edges).
    std::int64_t mNumBoundary = 0;
};

/**
 * @brief The Delaunay triangulation of a mesh's points in the `x, y` plane
 * (see the file comment).
 * @throws std::invalid_argument on a non-finite or out-of-range coordinate,
 *         1-D points, fewer than three distinct points, or collinear points.
 */
MESHIOPLUSPLUS_API Delaunay2dResult delaunay_2d(const Mesh& rMesh);

}  // namespace meshioplusplus
