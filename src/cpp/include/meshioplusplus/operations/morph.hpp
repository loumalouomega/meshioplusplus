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
 * @file operations/morph.hpp
 * @brief Boundary-driven mesh morphing (roadmap 7.4.1): move a mesh's
 * boundary (or any chosen points) by a prescribed displacement and carry the
 * motion smoothly into the rest of the mesh, keeping its connectivity. Shape
 * optimisation moves a design surface; fluid-structure interaction moves the
 * fluid mesh with the structure.
 *
 * `morph(mesh, {array})` reads a displacement per point from point data
 * `array` (`(n, dim)` or `(n, 3)`). The **constrained** points take it as
 * given: the mesh's boundary points (`mFixBoundary`; the points of facets
 * used by one top-dimensional cell), the points `mConstrained` or the point
 * data `mConstrainedArray` mark (nonzero), and every point no top-dimensional
 * cell uses. Every other point's displacement is the **harmonic extension**:
 * per component, `K u = 0` at the free points with `u = d` at the constrained
 * ones, `K` the P1 stiffness of the top-dimensional cells (each decomposed
 * into simplices as `convert_cells(mode="simplexify")` does, so hexahedra,
 * wedges, pyramids, quads and polygons work) solved by Jacobi-preconditioned
 * conjugate gradients. Each simplex's stiffness is scaled by
 * `(mean measure / its measure)^mStiffening` (Jacobian-based stiffening:
 * T. E. Tezduyar et al., 1992; K. Stein, T. Tezduyar and R. Benney, "Mesh
 * moving techniques for fluid-structure interactions with large
 * displacements", J. Appl. Mech. 70, 2003), so small cells move more rigidly
 * and keep their shape. With `mSteps` > 1 the displacement is applied in
 * equal increments, the stiffness reassembled on the moved mesh each time.
 *
 * **Inversion guard.** Afterwards every linear volume cell's signed volume
 * (`tetra`, `hexahedron`, `wedge`, `pyramid`) and every surface cell's
 * orientation (its Newell normal against the original's) is compared with
 * the input's: a cell that turned inside out is counted (`mNumInverted`)
 * and, unless `mAllowInverted`, the call is refused naming `steps` and
 * `stiffening`, the remedies. Points move; connectivity, data and regions
 * pass through unchanged.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// What `morph` moves, and how.
struct MorphOptions {
    /// Point data: the displacement, `(n, dim)` or `(n, 3)` (required).
    std::string mArray;
    /// Constrain the boundary points.
    bool mFixBoundary = true;
    /// Points to constrain too: empty, or one entry per point (nonzero pins).
    std::vector<std::uint8_t> mConstrained;
    /// Point data whose nonzero entries constrain their points too.
    std::string mConstrainedArray;
    /// Stiffness scaled by `(mean measure / measure)^mStiffening`; 0 for plain
    /// Laplace. Not negative.
    double mStiffening = 1.0;
    /// Increments the displacement is applied in; at least 1.
    std::int64_t mSteps = 1;
    /// Return a mesh with inverted cells instead of refusing it.
    bool mAllowInverted = false;
    /// Relative residual of each solve.
    double mTolerance = 1e-10;
    /// Iteration cap of each solve; 0 takes `2 n + 10` (at least 1000).
    std::int64_t mMaxIterations = 0;
};

/// What `morph` did.
struct MorphResult {
    /// The moved mesh.
    Mesh mMesh;
    /// Points that took the prescribed displacement, and points solved for.
    std::int64_t mNumConstrained = 0;
    std::int64_t mNumFree = 0;
    /// Cells that turned inside out (0 unless `mAllowInverted`).
    std::int64_t mNumInverted = 0;
    /// Conjugate-gradient iterations, all solves together.
    std::int64_t mIterations = 0;
    /// The largest displacement applied to a point.
    double mMaxDisplacement = 0.0;
    /// False when a solve hit its iteration cap.
    bool mConverged = true;
};

/**
 * @brief Move a mesh by a prescribed displacement of its boundary (or chosen
 * points), harmonically extended into the rest (see the file comment).
 * @throws std::invalid_argument on a missing or misshapen array, a negative
 *         stiffening, `mSteps` below 1, a mesh with no 2-D or 3-D cell, or,
 *         unless `mAllowInverted`, a result with inverted cells.
 */
MESHIOPLUSPLUS_API MorphResult morph(const Mesh& rMesh, const MorphOptions& rOptions);

}  // namespace meshioplusplus
