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
 * @file operations/eigenmodes.hpp
 * @brief The lowest eigenpairs of a triangle surface's Laplace-Beltrami
 * operator (roadmap 7.3.2): the surface's own Fourier basis. The modes are
 * spectral positional encodings for graph networks, shape descriptors
 * (the eigenvalues are invariant to rigid motion and parameterisation), and
 * the first modes of a membrane's vibration.
 *
 * `laplacian_eigenmodes(mesh, {num_modes})` solves the generalised problem
 * `K phi = lambda M phi` over the surface's P1 elements (`detail/surface_fem`,
 * the operators `geodesic_distance` uses): `K` the cotangent stiffness, `M` the
 * lumped (one-third area) mass, with the natural (Neumann) condition on a
 * boundary. The `num_modes` smallest eigenvalues come by block inverse
 * iteration with a shift (`(K + s M)^-1 M`, `s` one over the surface's area,
 * each solve Jacobi-preconditioned conjugate gradients) and Rayleigh-Ritz
 * over twice as many vectors (at least `num_modes + 8`), so repeated
 * eigenvalues (a sphere's `l (l + 1)`, `2 l + 1` times) come out whole. It
 * stops when every wanted pair's residual `|K phi - lambda M phi|` (in the
 * `M^-1` norm, relative to `lambda + s`) is below `mTolerance`.
 *
 * The result is `point_data["laplacian:modes"]`, Float64 `(n, num_modes)`,
 * the eigenvectors `M`-orthonormal (`phi_i^T M phi_j = delta_ij`), each signed
 * so its largest-magnitude entry is positive; and
 * `field_data["laplacian:eigenvalues"]`, ascending. A connected piece has a
 * zero eigenvalue with a constant mode. A point in no triangle gets `NaN`.
 * Within a repeated eigenvalue the basis is one of many, the same on every
 * run. Everything is the same on every backend and thread count.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data: Float64 `(n, num_modes)`, the eigenvectors.
inline constexpr const char* kLaplacianModesName = "laplacian:modes";
/// Field data: Float64 `(num_modes,)`, the eigenvalues, ascending.
inline constexpr const char* kLaplacianEigenvaluesName = "laplacian:eigenvalues";

/// What `laplacian_eigenmodes` computes.
struct EigenmodeOptions {
    /// Eigenpairs wanted, the smallest; at least 1, at most the surface's points.
    std::int64_t mNumModes = 10;
    /// The relative residual every wanted pair must reach.
    double mTolerance = 1e-8;
    /// Block iterations at most; 0 takes 300.
    std::int64_t mMaxIterations = 0;
};

/// What `laplacian_eigenmodes` found.
struct EigenmodeResult {
    /// The mesh with `laplacian:modes` and `laplacian:eigenvalues`.
    Mesh mMesh;
    /// The eigenvalues, ascending.
    std::vector<double> mEigenvalues;
    /// Block iterations run.
    std::int64_t mIterations = 0;
    /// The largest relative residual among the wanted pairs.
    double mMaxResidual = 0.0;
    /// False when `mMaxIterations` ran out first (the last iterate is kept).
    bool mConverged = false;
    /// Triangles left out for (near) zero area.
    std::int64_t mNumDegenerate = 0;
};

/**
 * @brief The lowest eigenpairs of the surface Laplacian (see the file comment).
 * @throws std::invalid_argument on a volume, quad, polygon or higher-order
 *         cell (naming the operation that makes linear triangles), no
 *         triangle, `mNumModes` below 1 or above the surface's points, or a
 *         non-positive tolerance.
 */
MESHIOPLUSPLUS_API EigenmodeResult laplacian_eigenmodes(const Mesh& rMesh,
                                                        const EigenmodeOptions& rOptions = {});

}  // namespace meshioplusplus
