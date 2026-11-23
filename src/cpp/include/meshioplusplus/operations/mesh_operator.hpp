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
 * @file operations/mesh_operator.hpp
 * @brief The mesh's own finite-element operators exported as sparse matrices:
 * the cotangent Laplacian, the lumped and the consistent mass, and the
 * per-triangle gradient of a triangle surface (roadmap 15.2.1).
 *
 * The same P1 operators `geodesic_distance`, `laplacian_eigenmodes` and
 * `sobolev_deform` solve with, returned as data so an application can build its
 * own systems on them. Nothing is solved here and no physics is assembled.
 *
 * | Kind | Matrix | Shape |
 * | --- | --- | --- |
 * | `Laplacian` | `K_ij = ∫ ∇φ_i · ∇φ_j dA`, the cotangent stiffness: symmetric, positive
 * semidefinite, rows summing to 0, the natural (Neumann) boundary | `n × n` | | `MassLumped` |
 * `M_ii = Σ A_t / 3` over the triangles at `i`, diagonal | `n × n` | | `MassConsistent` | `M_ij = ∫
 * φ_i φ_j dA`, `A_t / 12 (1 + δ_ij)` per triangle | `n × n` | | `Gradient` | row `3 t + k` is
 * component `k` of the gradient of the P1 field on triangle `t`, so `G f` is the per-triangle
 * gradient, laid out `(t, k)` | `3 T × n` |
 *
 * `n` is the point count and `T` the number of triangles of the mesh's
 * `triangle` blocks in input order (the `Gradient` rows of a zero-area
 * triangle are empty). The matrices are compressed rows with ascending,
 * unique columns, built serially, so they are the same on every backend.
 *
 * **Scope.** The surface is the mesh's `triangle` blocks, as for
 * `geodesic_distance`: a volume cell, a quad, a polygon or a higher-order
 * surface cell is refused naming the operation that fixes it. Tetrahedra wait
 * on the shared simplex kernels (roadmap 3.4.4).
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Which operator `mesh_operator` returns.
enum class MeshOperatorKind : std::int32_t {
    Laplacian = 0,
    MassLumped = 1,
    MassConsistent = 2,
    Gradient = 3,
};

/// What `mesh_operator` builds.
struct MeshOperatorOptions {
    MeshOperatorKind mKind = MeshOperatorKind::Laplacian;
};

/// A compressed-row sparse matrix.
struct MeshOperatorResult {
    std::int64_t mRows = 0;
    std::int64_t mCols = 0;
    /// `mRows + 1` row starts into `mIndices` and `mValues`.
    std::vector<std::int64_t> mIndptr;
    std::vector<std::int64_t> mIndices;
    std::vector<double> mValues;
    /// Zero-area triangles left out.
    std::int64_t mNumDegenerate = 0;
};

/**
 * @brief The operator of @p rOptions over the triangle surface of @p rMesh.
 * @throws std::invalid_argument on a cell type outside the scope, a triangle
 *         naming a point outside the mesh, no triangle or every triangle
 *         degenerate.
 */
MESHIOPLUSPLUS_API MeshOperatorResult mesh_operator(const Mesh& rMesh,
                                                    const MeshOperatorOptions& rOptions = {});

}  // namespace meshioplusplus
