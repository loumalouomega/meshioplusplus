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
 * @file detail/primitive_constraints.hpp
 * @brief Fit related primitives together (`FitCadOptions::mConstrainPrimitives`):
 * the cylinders, cones, tori and spheres of a turned part share one axis
 * exactly, planes take an axis or another plane's normal where they are
 * parallel or perpendicular to it, and cylinders, spheres and tori of nearly
 * equal radius take one radius.
 *
 * Each primitive is fitted on its own by `fit_primitives`, so the axes of a
 * noisy turned part differ in their last digits and the B-rep built from them
 * has sliver gaps where the faces should meet on one axis. This pass runs on
 * the `fit:patch_*` table `fit_cad` builds on, before either backend reads it:
 *
 *  1. **Coaxial groups.** The axial primitives (cylinder, cone, torus) are
 *     grouped in row order: a primitive joins the group of the first earlier
 *     one whose axis it is parallel to and lies on, within what moves its patch
 *     by twice its tolerance. A sphere joins a group whose axis passes through
 *     its centre. A group of two or more is refitted jointly by
 *     Levenberg-Marquardt: one axis, and each member's own parameters along it
 *     (a cylinder's radius; a cone's apex height and half-angle; a torus's
 *     centre height and radii; a sphere's centre height and radius).
 *  2. **Planes.** In row order, a plane parallel to a reference direction (the
 *     axes first, then the planes already seen) takes that direction as its
 *     normal; one perpendicular to one or two references is made exactly
 *     perpendicular to them. Its offset is refitted.
 *  3. **Equal radii.** Cylinders, spheres and torus minor radii within the
 *     tolerance of the first of their kind take their point-weighted mean.
 *
 * A change is **kept** only when the patch still passes the test
 * `fit_primitives` accepted it under (RMS distance and share of points within
 * the patch's tolerance) and its largest distance stays within the face
 * tolerance the writers demote at, or no larger than it was; otherwise the
 * primitive leaves the group and keeps its own fit. So the pass never demotes
 * a face. Everything is serial and walked in row order, so the result does not
 * depend on the thread count.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 */

// System includes
#include <cstdint>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// What `constrain_primitives` changed.
struct PrimitiveConstraintReport {
    /// Patches whose primitive took a shared axis, normal or radius.
    std::int64_t mNumConstrained = 0;
    /// Coaxial groups, parallel-normal families and equal-radius sets of two
    /// or more primitives.
    std::int64_t mNumGroups = 0;
};

/**
 * @brief Constrain the primitives of an annotated skin in place.
 *
 * @param rMesh a skin with `segment:patch` cell data and the `fit:patch_*`
 *        field data of `fit_primitives`; its `fit:patch_params`,
 *        `fit:patch_rms`, `fit:patch_max` and `fit:patch_fraction` rows and its
 *        `fit:rms`/`fit:max` cell data are rewritten for the patches changed.
 *        Without the table it is left as it is.
 * @param FaceTolerance the largest distance a face may keep from its patch
 *        before a writer demotes it; > 0.
 */
PrimitiveConstraintReport constrain_primitives(Mesh& rMesh, double FaceTolerance);

}  // namespace detail
}  // namespace meshioplusplus
