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
 * @file operations/point_normals.hpp
 * @brief Normals of a point cloud (roadmap 7.5.2, first part): scans read
 * from PLY, PCD or XYZ carry points and no surface, and shading, offsetting,
 * point-to-plane registration and reconstruction all need a normal per point.
 *
 * `estimate_normals(mesh, {k, orientation})` fits a plane to each point and
 * its `k` nearest neighbours (`operations/neighbors.hpp`, exact k-nearest
 * search, ties to the lower index): the normal is the eigenvector of the
 * neighbourhood's covariance for its least eigenvalue (H. Hoppe, T. DeRose, T.
 * Duchamp, J. McDonald and W. Stuetzle, "Surface reconstruction from
 * unorganized points", SIGGRAPH 1992). A fit gives a direction, not a side;
 * the side comes from `mOrientation`:
 *
 * - `Consistent` (the default): Hoppe et al.'s propagation. Over the
 *   symmetrised k-nearest graph, weighted `1 - |n_i . n_j|`, a minimum
 *   spanning tree is grown from each connected piece's highest point (Prim,
 *   ties to the lower index), and each normal is flipped to agree with its
 *   parent's. Then each piece is turned as a whole to point away from its
 *   centroid (outward, for a closed scan); where that is undecided (a flat
 *   patch) its highest point's normal points to +z, and failing that the
 *   normal's largest component is positive.
 * - `Viewpoint`: each normal points towards `mViewpoint` (the scanner).
 * - `None`: each normal's largest component is positive (deterministic, not
 *   consistent).
 *
 * The result is `point_data["normals"]` (Float64 `(n, 3)`, the name the PCD
 * and XYZ writers emit as normal columns; an existing array is replaced) and,
 * with `mRecordVariation`, `normals:variation`, the surface variation
 * `l0 / (l0 + l1 + l2)` of the neighbourhood's eigenvalues (0 on a plane, 1/3
 * for isotropic scatter; a curvature and noise indicator, M. Pauly, M. Gross
 * and L. Kobbelt, "Efficient simplification of point-sampled surfaces", 2002).
 * A point whose neighbourhood spans no plane (all on one line, or coincident)
 * gets a `NaN` normal and takes no part in the propagation
 * (`mNumDegenerate`).
 *
 * Data-only: points, cells, regions and the other arrays pass through. The
 * fits run in parallel, the propagation serially; the result is the same on
 * every backend and thread count.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data: Float64 `(n,)`, the neighbourhood's surface variation.
inline constexpr const char* kNormalsVariationName = "normals:variation";

/// How `estimate_normals` chooses each normal's side.
enum class NormalOrientation : std::int32_t {
    Consistent = 0,  ///< minimum spanning tree propagation, then outward
    Viewpoint = 1,   ///< towards `mViewpoint`
    None = 2,        ///< largest component positive
};

/// What `estimate_normals` fits, and how it orients.
struct PointNormalsOptions {
    /// Neighbours per fit (the point itself is added); at least 2, clamped
    /// to the point count minus one.
    std::int64_t mK = 16;
    NormalOrientation mOrientation = NormalOrientation::Consistent;
    /// The point normals face, for `Viewpoint`.
    std::array<double, 3> mViewpoint{0.0, 0.0, 0.0};
    /// Also attach `normals:variation`.
    bool mRecordVariation = false;
};

/// What `estimate_normals` did.
struct PointNormalsResult {
    /// The mesh with `normals` (and `normals:variation`).
    Mesh mMesh;
    /// Points whose neighbourhood spans no plane (`NaN` normals).
    std::int64_t mNumDegenerate = 0;
    /// `Consistent`: connected pieces of the neighbour graph oriented
    /// independently; 0 otherwise.
    std::int64_t mNumComponents = 0;
};

/** @brief `"consistent"`, `"viewpoint"` or `"none"`; anything else throws naming them. */
MESHIOPLUSPLUS_API NormalOrientation normal_orientation_from_name(const std::string& rName);

/**
 * @brief Normals of a mesh's points by local plane fits (see the file comment).
 * @throws std::invalid_argument on fewer than three points, `mK` below 2, a
 *         non-finite coordinate, or points with fewer than two coordinates.
 */
MESHIOPLUSPLUS_API PointNormalsResult estimate_normals(const Mesh& rMesh,
                                                       const PointNormalsOptions& rOptions = {});

}  // namespace meshioplusplus
