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
 * @file operations/metric.hpp
 * @brief A metric toolbox that needs no MMG (roadmap 7.4.3): intersect the
 * anisotropic metrics of several fields, normalise the result to a target
 * vertex count, clamp its sizes and grade it, ready for `adapt` (or any
 * metric-driven mesher).
 *
 * `prepare_metric(mesh, {arrays})` reads one or more metric fields
 * (`point_data`, six components `m11 m21 m22 m31 m32 m33` in 3-D or three
 * `m11 m21 m22` in 2-D, the order `metric_from_hessian` writes and a Medit
 * `.sol` holds) and, per point, in this order:
 *
 * 1. **Intersection.** The fields are folded by simultaneous reduction
 *    (P. J. Frey and P.-L. George, "Mesh Generation", 2008, §10.6): with
 *    `C = A^(-1/2) B A^(-1/2) = Q diag(mu) Q^T`, `A cap B = A^(1/2) Q
 *    diag(max(1, mu)) Q^T A^(1/2)`, the smallest metric that asks for at
 *    least what each does in every direction (its unit ball is the largest
 *    ellipsoid inside both): one mesh that resolves every field's features.
 * 2. **Normalisation** (`mTargetVertices` > 0): the continuous-mesh L^p
 *    normalisation of A. Loseille and F. Alauzet ("Continuous mesh framework
 *    part II", SIAM J. Numer. Anal. 49, 2011): `M <- N^(2/d) (int
 *    det(M)^(p/(2p+d)))^(-2/d) det(M)^(-1/(2p+d)) M`, `d` the metric's
 *    dimension (or 2 for a 3-D metric on a planar mesh with no volume cell,
 *    measured by its x-y block), the integral over the cells by lumped point
 *    measures, so the
 *    metric's complexity `int sqrt(det M)`, about the vertex count a mesher
 *    produces, becomes `N` (before clamping). `mNorm` is `p`; 0 is
 *    `p = infinity`, a uniform scaling that keeps the metric's shape.
 * 3. **Clamping** (`mHmin`, `mHmax` > 0): eigenvalues to `[1/hmax^2,
 *    1/hmin^2]`, so every size asked lies in `[hmin, hmax]`.
 * 4. **Gradation** (`mGradation` > 1): F. Alauzet's metric-space gradation
 *    ("Size gradation control of anisotropic meshes", Finite Elem. Anal.
 *    Des. 46, 2010): along every edge `pq`, of length `l` in `M_p`,
 *    `M_q <- M_q cap (1 + l ln(beta))^(-2) M_p` and the same the other way,
 *    swept over the edges in order until no metric changes by more than
 *    1e-3 (or `mMaxSweeps`), so sizes grow by at most a factor `beta` per
 *    unit length; then clamped again.
 *
 * The result is `point_data[mOutput]`, the width of the inputs. Everything
 * else passes through. Serial, deterministic.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// What `prepare_metric` does.
struct PrepareMetricOptions {
    /// Point data: the metric fields, all six (3-D) or all three (2-D) wide.
    std::vector<std::string> mArrays;
    /// The `point_data` name of the result.
    std::string mOutput = "metric";
    /// The vertex count to normalise to; 0 skips the normalisation.
    double mTargetVertices = 0.0;
    /// The norm `p` of the normalisation (at least 1); 0 is infinity.
    double mNorm = 2.0;
    /// The smallest and largest size asked; 0 leaves that side unclamped.
    double mHmin = 0.0;
    double mHmax = 0.0;
    /// The gradation `beta` (> 1); 0 skips the gradation.
    double mGradation = 0.0;
    /// Gradation sweeps at most.
    std::int64_t mMaxSweeps = 100;
};

/// What `prepare_metric` did.
struct PrepareMetricResult {
    /// The mesh with `point_data[mOutput]`.
    Mesh mMesh;
    /// The complexity `int sqrt(det M)` of the intersection, and of the result.
    double mComplexityBefore = 0.0;
    double mComplexityAfter = 0.0;
    /// Points whose metric the gradation changed, and its sweeps.
    std::int64_t mNumGraded = 0;
    std::int64_t mSweeps = 0;
};

/**
 * @brief Intersect, normalise, clamp and grade metric fields (see the file
 * comment).
 * @throws std::invalid_argument on no array, a missing one, fields of mixed or
 *         wrong width, a negative target, a norm in (0, 1), a gradation in
 *         (0, 1], `mHmin` > `mHmax`, or a normalisation of a mesh with no
 *         cell of the metric's dimension.
 */
MESHIOPLUSPLUS_API PrepareMetricResult prepare_metric(const Mesh& rMesh,
                                                      const PrepareMetricOptions& rOptions);

}  // namespace meshioplusplus
