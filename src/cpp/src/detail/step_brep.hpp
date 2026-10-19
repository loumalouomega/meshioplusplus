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
 * @file detail/step_brep.hpp
 * @brief The kernel-free STEP B-rep writer behind `StepFaces::Planar` and
 * `StepFaces::Fitted`: a closed skin as one face per surface instead of one
 * per facet, with no CAD kernel in the loop.
 *
 *  - **Faces.** With fitted patches (`segment:patch` and the `fit:patch_*`
 *    table of `fit_primitives`), every accepted plane, sphere, cylinder, cone
 *    or torus patch whose points all lie within the tolerance of its surface
 *    is one `ADVANCED_FACE` on a `PLANE`, `SPHERICAL_SURFACE`,
 *    `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE` or `TOROIDAL_SURFACE`. Every
 *    other facet joins an edge-connected coplanar group (largest facet
 *    first, every corner within `1e-6` of the diagonal of the seed's plane),
 *    one `PLANE` face each, holes included; a group whose points all lie
 *    within the tolerance of an analytic neighbour's surface joins that face.
 *  - **Edges.** The mesh edges between two faces form chains from corner to
 *    corner (a point three faces meet at, or where a boundary branches). A
 *    chain between two planes is split at its corners (Douglas-Peucker) into
 *    `LINE`s; one beside a curved analytic face is a `CIRCLE` when its points
 *    lie within the tolerance of a circle about a shared axis (a bore's rim),
 *    a `LINE` when they lie on a segment, and otherwise a degree-1
 *    `B_SPLINE_CURVE_WITH_KNOTS` through them; one between a curved face and
 *    a coplanar group keeps the mesh's straight edges. A closed chain with no
 *    corner is one closed edge.
 *  - **Vertices.** A point among analytic faces only is moved onto all their
 *    surfaces (alternating projections), when that moves it no further than
 *    the tolerance; any other keeps the mesh position.
 *  - **Loops.** Each face's boundary is walked through the facet fan around
 *    every point, so a loop touching itself stays two loops; a planar face's
 *    loop of largest area is its outer bound. A periodic face has no seam
 *    edge: STEP allows its bounds without one, and CAD readers add it.
 *
 * A patch whose loops do not close is written as coplanar groups instead,
 * and counted demoted. Deterministic and serial.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/mesh.hpp"

// Project includes (private, not installed)
#include "closed_facets.hpp"

namespace meshioplusplus {
namespace detail {

/// How the B-rep writer builds its faces.
struct StepBrepOptions {
    /// Use the fitted patches of `rSurface` (else coplanar groups only).
    bool mFitted = false;
    /// The largest distance, as a length, a point may lie from its face's
    /// surface. At most 0 takes three times the largest accepted
    /// `fit:patch_tolerance`, else `1e-3` of the bounding-box diagonal.
    double mTolerance = 0.0;
};

/// What the B-rep writer wrote.
struct StepBrepResult {
    std::vector<std::int64_t> mFaceOfFacet;   ///< per facet: its face (writing order)
    std::vector<std::int64_t> mTypeOfFacet;   ///< per facet: 0 a coplanar group, else 1..5
    std::vector<std::int64_t> mPatchOfFacet;  ///< per facet: its `segment:patch`, or -1
    double mTolerance = 0.0;
    double mMaxDeviation = 0.0;  ///< the largest distance of a point from its face's surface
    double mVolume = 0.0;        ///< the skin's enclosed volume
    bool mValid = false;         ///< every face's loops closed, every edge used twice
    std::int64_t mNumFaces = 0;
    std::int64_t mNumPlane = 0;
    std::int64_t mNumSphere = 0;
    std::int64_t mNumCylinder = 0;
    std::int64_t mNumCone = 0;
    std::int64_t mNumTorus = 0;
    std::int64_t mNumFacetedFaces = 0;  ///< coplanar-group faces
    std::int64_t mNumPatches = 0;       ///< rows of the fit table
    std::int64_t mNumDemoted = 0;       ///< accepted patches written as coplanar groups
    std::int64_t mNumSolids = 0;
    std::int64_t mNumLineEdges = 0;
    std::int64_t mNumCircleEdges = 0;
    std::int64_t mNumSplineEdges = 0;
};

/**
 * @brief Write @p rShell (the facets of @p rSurface) as a STEP B-rep.
 * @param rSurface the surface the facets came from; its `segment:patch` and
 *        `fit:patch_*` data are read when `mFitted`.
 * @throws WriteError on a skin whose faces cannot be closed into loops, or
 *         when the file cannot be written.
 */
StepBrepResult step_brep_write(const std::string& rPath, const Mesh& rSurface,
                               const ClosedFacets& rShell, const StepBrepOptions& rOptions);

/// `write_step` with options, also handing back the B-rep writer's report
/// (left as constructed for `StepFaces::Facet`). Implemented in
/// `formats/step.cpp`; `fit_cad`'s native backend reads the report.
void write_step_report(const std::string& rPath, const Mesh& rMesh,
                       const StepWriteOptions& rOptions, StepBrepResult* pReport);

}  // namespace detail
}  // namespace meshioplusplus
