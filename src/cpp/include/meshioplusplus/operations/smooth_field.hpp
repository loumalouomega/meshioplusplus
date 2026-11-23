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
#pragma once

/**
 * @file operations/smooth_field.hpp
 * @brief Screened-Poisson smoothing of a point field on a triangle surface
 * (roadmap 15.2.2): the field `g` solving `(M + l^2 K) g = M f` over the
 * lumped mass `M` and the cotangent stiffness `K` of `mesh_operator`, `l` a
 * length. A field varying on a scale much shorter than `l` is flattened, one
 * varying on a longer scale passes nearly unchanged; the mean (the `M`-weighted
 * average) is preserved exactly, and `l = 0` returns the field. The points do
 * not move (that is `sobolev_deform`'s job).
 *
 * `smooth_field(mesh, {array, length})` adds `point_data["<array>:smooth"]`
 * (or `mOutput`), each component solved on its own by Jacobi-preconditioned
 * conjugate gradients, serially and in input order: the same on every backend.
 * A point in no triangle keeps its value. The surface is the mesh's `triangle`
 * blocks, with `mesh_operator`'s scope and refusals.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The suffix of the smoothed array's default name (`<array>:smooth`).
inline constexpr const char* kSmoothFieldSuffix = ":smooth";

/// What `smooth_field` smooths, and how much.
struct SmoothFieldOptions {
    /// The point array to smooth: Float64-convertible, `(n,)` or `(n, k)`.
    std::string mArray;
    /// The smoothing length, in the mesh's units; 0 returns the field.
    double mLength = 0.0;
    /// The output name; empty takes `<mArray>:smooth`.
    std::string mOutput;
    /// Relative residual tolerance of the solve.
    double mTolerance = 1e-12;
    /// The iteration cap of each component's solve; 0 takes the point count plus 100.
    std::int64_t mMaxIterations = 0;
};

/// What `smooth_field` did.
struct SmoothFieldResult {
    /// The input mesh with the smoothed array.
    Mesh mMesh;
    /// Zero-area triangles left out.
    std::int64_t mNumDegenerate = 0;
    /// Conjugate-gradient iterations, summed over the components.
    std::int64_t mIterations = 0;
    /// False when a solve hit its cap; the last iterate is still used.
    bool mConverged = true;
};

/**
 * @brief Smooth @p rOptions.mArray over the triangle surface of @p rMesh.
 * @throws std::invalid_argument on an unknown or non-numeric array, one whose
 *         rows are not the point count, a negative length or tolerance, or a
 *         mesh outside `mesh_operator`'s scope.
 */
MESHIOPLUSPLUS_API SmoothFieldResult smooth_field(const Mesh& rMesh,
                                                  const SmoothFieldOptions& rOptions);

}  // namespace meshioplusplus
