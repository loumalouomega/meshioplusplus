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
 * @file detail/primitive_geom.hpp
 * @brief Point queries against a fitted primitive -- the distance, the
 * closest point and the surface normal of a plane, sphere, cylinder, cone or
 * torus given by its ten `fit:patch_params` entries -- and the per-patch fit
 * `fit_primitives` runs, for the operations that grow and build on patches
 * (`segment`'s primitive method, the kernel-free B-rep `step` writer and
 * `fit_cad`).
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. Every query is
 * closed-form and serial.
 */

// System includes
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// The `fit:patch_params` row length.
inline constexpr int kPrimitiveParams = 10;

/**
 * @brief The unsigned distance from @p pP to a primitive.
 * @param Type 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus.
 * @return the distance, or +inf for any other type.
 */
double primitive_point_distance(std::int64_t Type, const double* pParams, const double* pP);

/// The point of the primitive nearest @p pP (a cone's apex for a point behind
/// it); false, leaving @p pOut as @p pP, for an unknown type or a point on the
/// axis of a surface of revolution, where the nearest point is not unique.
bool primitive_project(std::int64_t Type, const double* pParams, const double* pP, double* pOut);

/// The unit normal at the point of the primitive nearest @p pP, pointing away
/// from the centre or axis (a plane's: its own normal); false where it is not
/// defined (see `primitive_project`).
bool primitive_normal(std::int64_t Type, const double* pParams, const double* pP, double* pOut);

/// One patch fitted the way `fit_primitives` fits it.
struct PrimitiveFit {
    std::int64_t mType = 0;  ///< 0 when no type is accepted
    double mParams[kPrimitiveParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double mRms = 0.0;
    double mMax = 0.0;
};

/**
 * @brief Fit a primitive to a patch's samples, exactly as `fit_primitives`
 * fits one patch: the first type of @p Types (bits `1 << type`) in the order
 * plane, sphere, cylinder, cone, torus whose RMS distance is within @p Tol and
 * whose share of points within @p Tol is at least @p MinInlierFraction (a
 * sphere giving way to a later type whose normals agree with the facets').
 *
 * @param rPts the patch's points, three per point.
 * @param rNormals its facets' unit normals, three per facet.
 * @param rCentres its facets' centroids, three per facet.
 * @param Tol the absolute tolerance.
 * (Implemented in `operations/fit_primitives.cpp`.)
 */
PrimitiveFit fit_primitive_samples(const std::vector<double>& rPts,
                                   const std::vector<double>& rNormals,
                                   const std::vector<double>& rCentres, double Tol,
                                   std::uint32_t Types, double MinInlierFraction, int MaxIt);

}  // namespace detail
}  // namespace meshioplusplus
