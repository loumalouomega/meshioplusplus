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
 * @file operations/locate_points.hpp
 * @brief Which cell of a mesh contains each point, and where in it (roadmap 9.3).
 *
 * `locate_points(mesh, points, options)` is the search `interpolate`'s
 * barycentric mode runs, as a query of its own. The mesh is first reduced to
 * simplices (`convert_cells(Simplexify)`: a hexahedron, a prism, a pyramid, a
 * polyhedron or a quadrilateral becomes tetrahedra or triangles; a
 * higher-order cell is linearised), and each point is looked up in a grid of
 * the simplices' bounding boxes: the first containing simplex in ascending
 * order wins, and a point on a shared face belongs to the lowest-numbered one.
 * A simplex contains a point when every barycentric weight is at least
 * `-mTolerance`.
 *
 * Per point the result gives:
 *
 * - `mCell`: the global (block-major) index of the **input cell** the
 *   containing simplex came from; -1 when the point is outside every cell;
 * - `mNodes`: the four point ids of that simplex (the input mesh's own point
 *   ids; a triangle has three and its fourth entry is -1);
 * - `mBarycentric`: the four weights with respect to those nodes (a triangle's
 *   fourth weight is 0), so `sum(w[k] * f[nodes[k]])` is `interpolate`'s value
 *   of a point array `f` bit for bit.
 *
 * Volume meshes (tetrahedra after simplexification) and *planar* triangle
 * meshes are located. A triangle's containment test is made in the xy-plane
 * (its `z` is ignored there), as `interpolate`'s is, but the grid lookup that
 * precedes it is three-dimensional, so a point must lie in the mesh's plane,
 * within a grid cell. The same lookup means a point is only compared with the
 * simplices whose bounding boxes share its grid cell: `mTolerance` widens the
 * barycentric test, not the lookup, so a point outside the mesh by more than
 * the rounding that keeps it in the grid cell is reported outside whatever the
 * tolerance. A surface in space is not located by containment:
 * `distance_to_surface` projects onto it and `cast_rays` intersects it. A mesh
 * with neither tetrahedra nor triangles throws.
 *
 * The points are searched in parallel; the result is the same on every
 * backend and thread count. The mesh is never modified.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How strictly a point must lie inside a simplex.
struct LocatePointsOptions {
    /// A weight at least `-mTolerance` counts as inside (dimensionless: the
    /// weights are volume ratios); the default is `interpolate`'s.
    double mTolerance = 1e-12;
};

/// One entry per point, in the order the points were given.
struct LocatePointsResult {
    std::vector<std::int64_t> mCell;   ///< Global input cell; -1 outside.
    std::vector<std::int64_t> mNodes;  ///< Four per point; -1 outside or padding.
    std::vector<double> mBarycentric;  ///< Four per point; 0 outside or padding.
    std::int64_t mNumLocated = 0;
};

/**
 * @brief Locate points in the cells of @p rMesh (see the file comment).
 * @param rPoints three numbers per point (a planar mesh's points carry z = 0).
 * @throws std::invalid_argument on a size that is not a multiple of three, a
 *         negative or non-finite tolerance, a non-finite coordinate, or a mesh
 *         with no tetrahedron or triangle after simplexification.
 */
MESHIOPLUSPLUS_API LocatePointsResult locate_points(const Mesh& rMesh,
                                                    const std::vector<double>& rPoints,
                                                    const LocatePointsOptions& rOptions = {});

}  // namespace meshioplusplus
