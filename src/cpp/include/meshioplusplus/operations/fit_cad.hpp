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
 * @file operations/fit_cad.hpp
 * @brief Turn the fitted patches of a closed triangle skin into a CAD solid --
 * plane, sphere, cylinder, cone and torus faces trimmed by the patch
 * boundaries, sewn, checked and optionally written as STEP -- through the
 * optional OpenCASCADE backend.
 *
 * The last step of the skin-to-CAD chain `segment` -> `fit_primitives` ->
 * `fit_cad`. The input is the closed, consistently wound skin those two
 * annotated: `segment:patch` per cell and the `fit:patch_*` field data. When
 * either is missing, both run first with their default options (a curvature
 * test or a tolerance of your own means running them yourself).
 *
 * ### How the solid is built
 *
 *  - **Patches.** A patch whose fit was accepted, and whose type is in
 *    `mTypes`, becomes one analytic face on its fitted surface. Every other
 *    surface cell (an unfitted patch, a type left out, a cell with no patch)
 *    becomes planar faces, one per facet, as the faceted `step` writer
 *    writes them.
 *  - **Edges.** Where two patches meet, the mesh edges between them form
 *    chains from one junction (a point three or more patches meet at) to the
 *    next; a closed chain is split in two. A chain between two analytic faces
 *    is, when `mExactEdges` is set, the branch of the two surfaces'
 *    intersection that lies within the tolerance of the chain, and otherwise
 *    a B-spline through the chain, smoothed and moved onto both surfaces,
 *    kept only within the tolerance of both. A chain beside a facet keeps the
 *    mesh's straight edges. Junctions between analytic faces
 *    only are moved onto those surfaces.
 *  - **Solid.** Each face gets parameter-space curves continuous around its
 *    loops (a boundary along a cylinder's seam or a torus's equator would
 *    otherwise fall apart); OpenCASCADE adds the seams, repairs and sews the
 *    faces and closes each shell into one solid per edge-connected piece of
 *    the skin. A patch with no boundary (a whole sphere or torus) is its
 *    surface's natural face.
 *  - **Checks, and demotion.** A face that the B-rep check rejects, or that
 *    lies further than the tolerance from the skin (or leaves part of its
 *    patch further than that from the solid), is not kept: its patch is
 *    demoted to facets and the solid is rebuilt, for at most
 *    `mMaxDemotionRounds` rounds, after which every patch is faceted. So a
 *    feature the fitter cannot explain comes back faceted, never silently
 *    wrong.
 *
 * ### Output
 *
 * The solid re-meshed by OpenCASCADE at `mMeshDeflection`, as one `triangle`
 * block (points shared between faces are welded), with cell data
 * `fit_cad:face` (Int64, the B-rep face), `fit_cad:type` (Int64, the face's
 * type: 0 a facet, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus, as
 * `fit:type`) and `fit_cad:patch` (Int64, the `segment:patch` label it was
 * built from, -1 for a facet of a cell with no patch), and the report below.
 * With a non-empty path, the solid is also written as AP214 STEP in
 * millimetres, with the header timestamp pinned so the bytes are
 * reproducible.
 *
 * Built only with `-DMESHIOPLUSPLUS_WITH_OCCT=ON` (OpenCASCADE 7.6 or later,
 * LGPL-2.1 with the OCCT exception, linked dynamically); elsewhere `fit_cad`
 * throws `Unsupported` naming that flag, and `fit_cad_has_occt` is false.
 * The work is serial, so the output does not depend on the thread count.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"

namespace meshioplusplus {

/// Cell data (Int64) of the re-meshed solid: the B-rep face each triangle lies on.
inline constexpr const char* kFitCadFaceName = "fit_cad:face";
/// Cell data (Int64): the type of that face; 0 a facet, else a `FitPrimitiveType`.
inline constexpr const char* kFitCadTypeName = "fit_cad:type";
/// Cell data (Int64): the `segment:patch` label the face was built from, or -1.
inline constexpr const char* kFitCadPatchName = "fit_cad:patch";

/// What `fit_cad` builds and how strictly it checks it.
struct FitCadOptions {
    /// The largest distance, as a length, a face may lie from the skin. At most
    /// 0 takes three times the largest `fit:patch_tolerance` of an accepted
    /// patch, or `1e-3` of the bounding-box diagonal when no patch was accepted.
    double mTolerance = 0.0;
    /// Whether an edge between two analytic faces is first tried as the exact
    /// intersection of their surfaces (true) or always approximated.
    bool mExactEdges = true;
    /// Which accepted types become analytic faces, as `1 << FitPrimitiveType`
    /// bits; must name nothing but those five. 0 facets everything.
    std::uint32_t mTypes = kFitAllTypes;
    /// The chordal deflection of the re-meshed solid; at most 0 takes a quarter
    /// of the tolerance.
    double mMeshDeflection = 0.0;
    /// How many times a failed face may be demoted and the solid rebuilt
    /// before every patch is faceted; at least 0.
    std::int64_t mMaxDemotionRounds = 8;
};

/// The re-meshed solid and what was built.
struct FitCadResult {
    /// The solid re-meshed, with the `fit_cad:*` cell data.
    Mesh mMesh;
    /// Whether the solid passed OpenCASCADE's B-rep check.
    bool mValid = false;
    /// The tolerance applied.
    double mTolerance = 0.0;
    /// The solid's volume.
    double mVolume = 0.0;
    /// The volume the skin encloses.
    double mSkinVolume = 0.0;
    /// `mVolume / mSkinVolume`.
    double mVolumeRatio = 0.0;
    /// The symmetric Hausdorff distance between the skin and the re-meshed solid.
    double mMaxDeviation = 0.0;
    /// B-rep faces in all, and by type.
    std::int64_t mNumFaces = 0;
    std::int64_t mNumPlane = 0;
    std::int64_t mNumSphere = 0;
    std::int64_t mNumCylinder = 0;
    std::int64_t mNumCone = 0;
    std::int64_t mNumTorus = 0;
    /// Planar facet faces.
    std::int64_t mNumFacetedFaces = 0;
    /// Patches in the skin.
    std::int64_t mNumPatches = 0;
    /// Accepted patches demoted to facets by the checks: the features dropped.
    std::int64_t mNumDemoted = 0;
    /// Edges between analytic faces built as exact intersections, and as B-splines.
    std::int64_t mNumExactEdges = 0;
    std::int64_t mNumApproxEdges = 0;
    /// Solids written (one per edge-connected piece of the skin).
    std::int64_t mNumSolids = 0;
};

/**
 * @brief Build a CAD solid from the fitted patches of a closed skin.
 *
 * @param rMesh a closed, consistently wound surface mesh, normally the output
 *        of `fit_primitives`; lines and vertices are ignored.
 * @param rStepPath where to write the solid as STEP; empty writes nothing.
 * @param rOptions see `FitCadOptions`.
 * @throws Unsupported in a build without OpenCASCADE.
 * @throws std::invalid_argument on an option out of range, a volume block, or
 *         a surface that is not closed and consistently wound.
 * @throws WriteError when the STEP file cannot be written.
 */
MESHIOPLUSPLUS_API FitCadResult fit_cad(const Mesh& rMesh, const std::string& rStepPath = {},
                                        const FitCadOptions& rOptions = {});

/// Whether this build carries the OpenCASCADE backend `fit_cad` needs.
MESHIOPLUSPLUS_API bool fit_cad_has_occt() noexcept;

}  // namespace meshioplusplus
