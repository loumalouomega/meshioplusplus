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
 * @file detail/brep_deviation.hpp
 * @brief The two-sided distance and the normal deviation between a trimmed
 * NURBS B-rep (`BrepInfo`) and the triangle skin it was fitted to, for
 * `fit_cad`'s validation report (doc/fit_cad_freeform_plan.md, step 2).
 *
 *  - **CAD to mesh.** Each face is sampled on a grid of `N + 1` nodes a side
 *    over the bounding box of its trims in its parameter space (its knot
 *    domain when untrimmed), a node kept when the even-odd rule puts it inside
 *    the trims (`GeometryUtils::PointInPolygon2D`, each trimming curve
 *    evaluated as a polyline), and each kept node's unsigned distance to the
 *    skin is measured.
 *  - **Mesh to CAD.** Every skin vertex and facet centre is measured to the
 *    face its facet was built into (`pTriFace`, a vertex to each of its
 *    facets' faces): Gauss-Newton from the face's nearest sample node, the
 *    foot kept when it lies inside the trims, else the nearest point of the
 *    face's trim curves on the surface. Without `pTriFace` the face is the
 *    one of the nearest triangle of the faces' trimmed domains triangulated
 *    in their parameters, the next ones tried when the foot falls outside.
 *  - **Normals.** At each facet centre's foot point, the unsigned angle (in
 *    `[0, 90]` degrees) between the facet's normal and the surface's: a
 *    face's orientation flag does not enter, so the measure is the same for
 *    every backend.
 *
 * Every pass is serial or runs on the deterministic distance queries, in a
 * fixed sample order, so nothing depends on the thread count. Maxima are
 * exact over the samples; quantiles are `coverage_quantile`'s.
 *
 * Core-private (not installed).
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"

namespace meshioplusplus::detail {

/// How densely a B-rep is sampled.
struct BrepDeviationOptions {
    /// Grid cells per parameter direction of each face (nodes: one more).
    std::size_t mSamplesPerDirection = 24;
    /// Polyline points per knot span of a trimming curve of degree above 1.
    std::size_t mTrimSamplesPerSpan = 8;
    /// Gauss-Newton steps of a mesh-to-CAD foot point.
    std::size_t mMaxNewtonIterations = 12;
};

/// One face's measures (NaN-free; 0 where a side has no sample).
struct BrepFaceDeviation {
    double mCadToMeshMax = 0.0;
    double mCadToMeshP95 = 0.0;
    double mMeshToCadMax = 0.0;
    double mNormalMaxDeg = 0.0;
    std::int64_t mNumSamples = 0;  ///< its CAD-to-mesh samples
};

/// The B-rep against the skin.
struct BrepDeviation {
    /// False when there was nothing to measure (no face, or no skin
    /// triangle); `mWhy` then says which.
    bool mOk = false;
    std::string mWhy;
    double mCadToMeshMax = 0.0;
    double mCadToMeshP95 = 0.0;
    double mCadToMeshP99 = 0.0;
    double mMeshToCadMax = 0.0;
    double mMeshToCadP95 = 0.0;
    double mMeshToCadP99 = 0.0;
    double mNormalMaxDeg = 0.0;
    double mNormalP95Deg = 0.0;
    /// Per face of `BrepInfo::mFaces`, in order.
    std::vector<BrepFaceDeviation> mFaces;
};

/**
 * @brief Measure @p rBrep against the skin of points @p rXyz (three per
 * point) and triangles @p rTri (three point ids each).
 */
BrepDeviation brep_deviation(const BrepInfo& rBrep, const std::vector<double>& rXyz,
                             const std::vector<std::int64_t>& rTri,
                             const BrepDeviationOptions& rOptions = {},
                             const std::vector<std::int64_t>* pTriFace = nullptr);

}  // namespace meshioplusplus::detail
