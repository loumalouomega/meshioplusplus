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
 * @file operations/geodesic.hpp
 * @brief Geodesic distance on a triangle surface by the heat method (Crane,
 * Weischedel and Wardetzky, "Geodesics in Heat", ACM TOG 32(5), 2013).
 *
 * `geodesic_distance(mesh, {sources})` returns the mesh with
 * `point_data["geodesic:distance"]`: every point's distance, along the
 * surface, to the nearest source point. Three steps, over the cotangent
 * stiffness `K` and the lumped (one-third area) mass `M` of the triangles:
 *
 * 1. Heat: solve `(M + t K) u = delta` for a short time `t = mTimeFactor h^2`,
 *    `h` the mean edge length, `delta` 1 at the sources and 0 elsewhere.
 * 2. Direction: per triangle, `X = -grad u / |grad u|`, the unit vector
 *    pointing away from the sources.
 * 3. Distance: solve the Poisson problem `K phi = -div X` (`div` the
 *    integrated divergence at each vertex), then shift `phi` so its least value
 *    over each connected piece is 0.
 *
 * The heat solve is serial Gauss-Seidel stopped on every value's own relative
 * change, which keeps the far field's tiny values accurate where a
 * residual-based stop would not; the Poisson solve (one vertex of each
 * connected piece held at 0, which removes its constant null space) is
 * Jacobi-preconditioned conjugate gradients. Both are deterministic: the result
 * is the same on every backend and thread count. The heat method is first-order
 * accurate: the error is of the order of the mesh size, and halves when the
 * mesh is refined. A boundary carries the natural (Neumann) condition, which
 * bends the distance near it slightly; Crane et al.'s Dirichlet average is not
 * applied.
 *
 * **Scope.** The surface is the mesh's `triangle` blocks. A volume cell, a
 * quad, a polygon or a higher-order surface cell is refused naming the
 * operation that makes triangles (`extract_surface`, `convert_cells`); line
 * and vertex blocks ride along. A point in no triangle, and every point of a
 * connected piece with no source, gets `NaN`. A zero-area triangle adds
 * nothing and is counted. There is no numpy twin (the `sobolev_deform`
 * precedent): the iterates depend on a rounded stopping test.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data: Float64 `(n,)`, the surface distance to the nearest source.
inline constexpr const char* kGeodesicDistanceName = "geodesic:distance";

/// What `geodesic_distance` measures from, and how.
struct GeodesicOptions {
    /// Source point ids, 0-based.
    std::vector<std::int64_t> mSources;
    /// A region whose points are sources too (a point region's entries, a cell
    /// region's cells' points, a side region's facets' points); empty for none.
    std::string mSourceRegion;
    /// The heat time is `mTimeFactor * h^2`, `h` the mean edge length (Crane
    /// et al. recommend 1; larger smooths the distance, smaller sharpens it).
    double mTimeFactor = 1.0;
    /// Relative residual tolerance of both solves.
    double mTolerance = 1e-10;
    /// Iteration cap of each solve; 0 takes twice the point count plus 10, at
    /// least 1000 (Gauss-Seidel on a tiny surface can need more than 2 n + 10).
    std::int64_t mMaxIterations = 0;
};

/// What `geodesic_distance` did.
struct GeodesicResult {
    /// The input mesh with `point_data["geodesic:distance"]`.
    Mesh mMesh;
    /// The distinct source points used.
    std::int64_t mNumSources = 0;
    /// Points that got `NaN`: in no triangle, or in a piece with no source.
    std::int64_t mNumUnreached = 0;
    /// Zero-area triangles left out.
    std::int64_t mNumDegenerate = 0;
    /// The largest finite distance.
    double mMaxDistance = 0.0;
    /// The heat time `t` used.
    double mTimeStep = 0.0;
    /// Conjugate-gradient iterations of the heat and the Poisson solve.
    std::int64_t mHeatIterations = 0;
    std::int64_t mPoissonIterations = 0;
    /// False when either solve hit its cap; the last iterate is still used.
    bool mConverged = true;
};

/**
 * @brief The surface distance from @p rOptions' sources to every point.
 * @throws std::invalid_argument without a source, with a source out of range
 *         or in no triangle, a non-positive `mTimeFactor`, no triangle, or a
 *         cell type outside the scope (see the file comment).
 */
MESHIOPLUSPLUS_API GeodesicResult geodesic_distance(const Mesh& rMesh,
                                                    const GeodesicOptions& rOptions);

}  // namespace meshioplusplus
