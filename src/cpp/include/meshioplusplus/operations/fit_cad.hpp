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
 * optional OpenCASCADE backend, or, without it, through the kernel-free
 * native backend.
 *
 * The last step of the skin-to-CAD chain `segment` -> `fit_primitives` ->
 * `fit_cad`. The input is the closed, consistently wound skin those two
 * annotated: `segment:patch` per cell and the `fit:patch_*` field data.
 * Without `segment:patch`, `segment` runs first with `SegmentMethod::Primitive`
 * (so a fillet comes apart from the planes it is tangent to) and
 * `fit_primitives` at the same absolute tolerance, `1e-3` of the
 * bounding-box diagonal; with the patches but without the fit table,
 * `fit_primitives` runs with its default options (a curvature test or a
 * tolerance of your own means running them yourself).
 *
 * ### How the solid is built
 *
 *  - **Patches.** A patch whose fit was accepted, and whose type is in
 *    `mTypes`, becomes one analytic face on its fitted surface. A patch no
 *    primitive explains becomes, with `mFreeform`, a B-spline face (roadmap
 *    6.3.2, below). Every other surface cell (a patch neither explains, a
 *    type left out, a cell with no patch) becomes planar faces: each group of
 *    edge-connected coplanar facets one face, its holes inner loops.
 *  - **Freeform faces.** The patch is flattened over its best-fit plane, or,
 *    when that folds it, unrolled around the axis its normals are most
 *    perpendicular to (a cylinder of the circle its points fit). A patch both
 *    fold is, with `mFreeformSplit`, split into pieces one of them does
 *    parametrise -- halved at the area-weighted median of `n . e`, `e` the
 *    direction its normals spread along most, up to four times -- and each
 *    piece becomes its own freeform face under a new patch label, numbered
 *    after the largest existing one; a piece under 16 triangles, or still
 *    folding, is left to facets. Its offset from that plane or cylinder is
 *    resampled on a grid (inside the patch by interpolation over its facets,
 *    outside it from the nearest boundary point), approximated by a C2
 *    B-spline of degree `mFreeformDegree` (`GeomAPI_PointsToBSplineSurface`)
 *    and kept once every patch vertex lies within the tolerance of it, the
 *    grid doubling up to `mFreeformMaxControl` points per direction. The
 *    face is then trimmed, sewn and checked like an analytic one.
 *  - **Fillets and threads** (roadmap 6.3.1). Before the freeform fit, a
 *    patch between two analytic faces is tried as a rolling-ball fillet
 *    (`mFillets`; a pipe of one radius, else a radius varying along it), and
 *    a patch as a helical thread face (`mThreads`; about a neighbouring
 *    cylinder's or cone's axis, else the screw the patch fits, with any
 *    profile). After it, a freeform patch between two fitted faces, one of
 *    them freeform, is tried as a rolling-ball blend.
 *  - **Chamfers, conic blends, tapered and variable-pitch threads, vertex
 *    blends** (roadmap 6.3, v17.2.0). A patch between two analytic faces that
 *    no rolling ball or thread explains is tried as a blend whose
 *    cross-section, in the plane normal to its rail on the first face, is a
 *    straight chord (a chamfer, `fit_cad:type` 10) or a conic arc tangent to
 *    both faces (a conic blend, type 7); a thread flank whose straight-flank
 *    or profile fit fails is refitted with a pitch that varies linearly along
 *    it and, for a profile, cut on a taper -- a flank beside a cone is
 *    reported tapered. Last, a patch that is one closed loop around three or
 *    more fitted neighbours (a corner where blends meet, a setback) is tried
 *    as a vertex blend (type 9): a plate surface through its vertices,
 *    tangent to every neighbour along the loop. These run only where the
 *    earlier fits failed, so a skin they do not concern builds the same
 *    solid as before. The first three ride on `mFillets` and `mThreads`, the
 *    vertex blend on `mFillets`.
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
 *    demoted to facets and the solid is rebuilt. A failure is blamed on the
 *    faces the B-rep check names, or that a free edge after sewing touches,
 *    never on every patch at once; past `mMaxDemotionRounds` rounds each
 *    round demotes those patches and their neighbours, so the loop ends with
 *    the patches nothing blames still analytic. A face's distance allowance
 *    is the tolerance or, when larger, the skin's own chordal sag over its
 *    patch. So a feature the fitter cannot explain comes back faceted, never
 *    silently wrong.
 *  - **Merging.** With `mMergeFaces`, faces on one surface are merged and
 *    collinear edges joined (`ShapeUpgrade_UnifySameDomain`), kept only when
 *    the solid stays valid with the same volume, so the faceted remainder
 *    comes out as polygons rather than triangles.
 *
 * ### Output
 *
 * The solid re-meshed by OpenCASCADE at `mMeshDeflection`, as one `triangle`
 * block (points shared between faces are welded), with cell data
 * `fit_cad:face` (Int64, the B-rep face), `fit_cad:type` (Int64, the face's
 * type: 0 a facet, 1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus, as
 * `fit:type`, 6 a freeform B-spline, 7 a rolling-ball fillet, blend or conic
 * blend, 8 a helical thread face, 9 a vertex blend, 10 a chamfer) and
 * `fit_cad:patch` (Int64, the
 * `segment:patch` label it was built from, -1 for a facet of a cell with no
 * patch), the report below, and the solid as NURBS in `mBrep` (face `i` of
 * `fit_cad:face` is its face `i`, id `i + 1`), which `write_cad_json` writes
 * as Kratos CAD JSON.
 * With a non-empty path, the solid is also written as AP214 STEP in
 * millimetres, with the header timestamp pinned so the bytes are
 * reproducible.
 *
 * The OpenCASCADE backend is built only with `-DMESHIOPLUSPLUS_WITH_OCCT=ON`
 * (OpenCASCADE 7.6 or later, LGPL-2.1 with the OCCT exception, linked
 * dynamically); elsewhere `FitCadBackend::Occt` throws `Unsupported` naming
 * that flag, `fit_cad_has_occt` is false, and `Auto` runs the native backend.
 *
 * ### The native backend
 *
 * `FitCadBackend::Native` needs no CAD kernel: every accepted plane, sphere,
 * cylinder, cone and torus patch whose type is in `mTypes` becomes one
 * analytic STEP face, bounded by line, circle and polyline edges, and the
 * rest coplanar polygons -- the kernel-free `step` writer with
 * `StepFaces::Fitted`. A patch whose corners cannot be moved onto its
 * surface within the tolerance is demoted to polygons. There are no
 * freeform, fillet or thread faces and no re-meshing: the returned mesh is
 * the skin itself with the `fit_cad:*` cell data, `mBrep` is empty,
 * `mVolume` and `mVolumeRatio` are NaN (`mSkinVolume` is the skin's), and
 * `mValid` is the writer's topological check (a closed shell, every edge
 * used twice in opposite senses).
 *
 * The work is serial, so the output does not depend on the thread count.
 */

// System includes
#include <cstdint>
#include <limits>
#include <string>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"

namespace meshioplusplus {

/// Cell data (Int64) of the re-meshed solid: the B-rep face each triangle lies on.
inline constexpr const char* kFitCadFaceName = "fit_cad:face";
/// Cell data (Int64): the type of that face; 0 a facet, else a `FitPrimitiveType`.
inline constexpr const char* kFitCadTypeName = "fit_cad:type";
/// Cell data (Int64): the `segment:patch` label the face was built from, or -1.
inline constexpr const char* kFitCadPatchName = "fit_cad:patch";
/// Field data, with `FitCadOptions::mFaceDeviation`: per face (`fit_cad:face`
/// order) the largest distance from the face's re-mesh to the skin and from
/// the skin it stands for to the solid, an `(n, 2)` float64 array (v17.11.0).
inline constexpr const char* kFitCadFaceHausdorffName = "fit_cad:face_hausdorff";
/// Field data, with `FitCadOptions::mFaceDeviation`: per face, the root mean
/// square of both sides' distances.
inline constexpr const char* kFitCadFaceRmsName = "fit_cad:face_rms";
/// `fit_cad:type` of a freeform B-spline face.
inline constexpr std::int64_t kFitCadBsplineType = 6;
/// `fit_cad:type` of a rolling-ball fillet or blend face, of constant or
/// varying radius, between analytic or freeform faces.
inline constexpr std::int64_t kFitCadFilletType = 7;
/// `fit_cad:type` of a helical thread face: a straight flank (a ruled
/// helicoid), or a profile of any shape (a rounded root or crest) swept by
/// the thread's screw motion; tapered or of varying pitch since v17.2.0.
inline constexpr std::int64_t kFitCadHelicalType = 8;
/// `fit_cad:type` of a vertex blend: the corner patch where three or more
/// blends meet, a plate surface tangent to its neighbours (v17.2.0).
inline constexpr std::int64_t kFitCadVertexBlendType = 9;
/// `fit_cad:type` of a chamfer: a blend of straight cross-sections, of
/// constant or varying width, that is no plane or cone (v17.2.0).
inline constexpr std::int64_t kFitCadChamferType = 10;

/// Which engine `fit_cad` builds the solid with (since v17.9.0).
enum class FitCadBackend : std::int32_t {
    /// OpenCASCADE when the build carries it, else `Native`.
    Auto = 0,
    /// OpenCASCADE: analytic, freeform and feature faces, sewn and checked;
    /// throws `Unsupported` in a build without it.
    Occt = 1,
    /// No CAD kernel: the patches `fit_primitives` explains become analytic
    /// faces and the rest coplanar polygons, written by the kernel-free B-rep
    /// `step` writer (`StepFaces::Fitted`); no freeform or feature faces.
    Native = 2,
};

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
    /// Whether a patch no primitive explains is tried as a B-spline face.
    bool mFreeform = true;
    /// The freeform B-spline's highest degree per direction, 3 to 8.
    std::int32_t mFreeformDegree = 3;
    /// The most grid points per direction a freeform fit may sample, 4 to 1024.
    std::int64_t mFreeformMaxControl = 64;
    /// Whether a freeform patch that neither its best-fit plane nor an unrolled
    /// cylinder parametrises (it folds over both) is split into pieces that one
    /// does, each fitted as its own B-spline face under a new patch label.
    bool mFreeformSplit = true;
    /// Whether a patch no primitive explains, running between two fitted
    /// faces, is tried as a rolling-ball fillet: of constant or varying
    /// radius between two analytic faces, and as a blend beside a freeform one;
    /// then as a chamfer or a conic blend between two analytic faces, and a
    /// corner patch among three or more fitted faces as a vertex blend.
    bool mFillets = true;
    /// Whether a patch no primitive explains is tried as a helical thread
    /// face: about the axis of a cylinder or cone beside it, or of the screw
    /// the patch itself fits, with a straight or a rounded profile, of
    /// constant or varying pitch, and cut on a taper.
    bool mThreads = true;
    /// Whether faces on one surface are merged and collinear edges joined
    /// once the solid is built (`ShapeUpgrade_UnifySameDomain`), kept only
    /// when the solid stays valid with the same volume (since v17.9.0).
    bool mMergeFaces = true;
    /// The engine (since v17.9.0).
    FitCadBackend mBackend = FitCadBackend::Auto;
    /// Whether the patches are first refined by `refine_patches` (variational
    /// shape approximation), with `mRefine`, after any segmentation and fit
    /// `fit_cad` runs itself (v17.11.0). Off by default, so the output is
    /// what it was before the option existed.
    bool mRefineSegmentation = false;
    /// The refinement's options, used when `mRefineSegmentation` is set.
    RefinePatchesOptions mRefine;
    /// Whether a patch the plane and unrolled-cylinder freeform fit leaves
    /// unexplained is tried as a general freeform B-spline face: a disc
    /// parametrisation (least-squares conformal or mean-value) and an
    /// iterative least-squares fit with parameter correction and knot
    /// insertion, the patch halved up to `mMaxFreeformSplits` times when that
    /// fails (v17.11.0). Off by default.
    bool mFreeformGeneral = false;
    /// The general fit's target at the vertices (a facet centre is held to the
    /// larger of the tolerance and twice it); at most 0 takes a quarter of the
    /// tolerance, room for the face's edges and its re-mesh. The face is then
    /// accepted by the same test as any freeform face, half the tolerance at a
    /// vertex and the tolerance at a facet centre, so a larger target never
    /// loosens it.
    double mFreeformTolerance = 0.0;
    /// The thin-plate weight of the general fit, relative to the sample count
    /// on a unit parameter square; >= 0. It only holds poles no sample reaches.
    double mFreeformSmoothing = 1e-12;
    /// The most poles per direction of a general freeform face, 4 to 1024.
    std::int64_t mMaxFreeformPoles = 40;
    /// The general fit's iteration budget, >= 1.
    std::int64_t mMaxFreeformIterations = 30;
    /// How many times a patch the general fit cannot explain is halved, >= 0.
    std::int64_t mMaxFreeformSplits = 2;
    /// The tolerance the faces are sewn at; 0 sews at the tolerance itself,
    /// as before the option existed (v17.11.0). OpenCASCADE backend only.
    double mSewTolerance = 0.0;
    /// Whether a freeform B-spline face's surface, extended past its patch on
    /// every open side (`GeomLib::ExtendSurfByLength`, by a tenth of the
    /// patch's size), is intersected with a neighbour where the surfaces
    /// themselves give no branch along their shared chain, so an edge that runs
    /// just outside the fitted area is still exact. The face keeps its own
    /// surface. OpenCASCADE backend only; off by default (v17.11.0).
    bool mExtendSurfaces = false;
    /// Whether the B-rep vertices where analytic faces meet are placed by a
    /// least-squares solve on the surfaces' tangent planes (Gauss-Newton,
    /// OpenCASCADE-free) before the averaged alternating projection; the same
    /// acceptance applies. OpenCASCADE backend only; off by default (v17.11.0).
    bool mSnapJunctions = false;
    /// Whether the two-sided distance between each face and the skin is
    /// measured, as the `fit_cad:face_hausdorff` and `fit_cad:face_rms` field
    /// data and `FitCadResult::mHausdorffMax`/`mHausdorffRms` (v17.11.0).
    /// OpenCASCADE backend only; off by default.
    bool mFaceDeviation = false;
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
    /// Freeform B-spline faces.
    std::int64_t mNumBspline = 0;
    /// Pieces that freeform patches folding over every flattening were split into.
    std::int64_t mNumSplitPatches = 0;
    /// Rolling-ball fillet faces.
    std::int64_t mNumFillet = 0;
    /// Thread faces (helical).
    std::int64_t mNumHelical = 0;
    /// Of the fillet faces, those whose radius varies along them.
    std::int64_t mNumVariableFillet = 0;
    /// Of the fillet faces, the blends beside a freeform B-spline face.
    std::int64_t mNumBlend = 0;
    /// Of the thread faces, those whose profile is not a straight line (a
    /// rounded root or crest).
    std::int64_t mNumRoundedThread = 0;
    /// Vertex blends: corner patches where three or more blends meet (v17.2.0).
    std::int64_t mNumVertexBlend = 0;
    /// Chamfers: blends of straight cross-sections (v17.2.0).
    std::int64_t mNumChamfer = 0;
    /// Of the fillet faces, the conic blends: cross-sections that are conic
    /// arcs other than the rolling ball's circle (v17.2.0).
    std::int64_t mNumConicBlend = 0;
    /// Of the thread faces, those cut on a taper (v17.2.0).
    std::int64_t mNumTaperedThread = 0;
    /// Of the thread faces, those whose pitch varies along them (v17.2.0).
    std::int64_t mNumVariablePitchThread = 0;
    /// The solid as NURBS faces, loops and edges (empty from `Native`).
    BrepInfo mBrep;
    /// The engine that built the solid: `Occt` or `Native` (since v17.9.0).
    FitCadBackend mBackend = FitCadBackend::Occt;
    /// Whether `mRefineSegmentation` refined the patches (v17.11.0); the
    /// counts below are zero when it did not.
    bool mRefined = false;
    /// Whether the solid was built from the refined patches (v17.11.0). The
    /// OpenCASCADE backend also builds it from the patches as given and keeps
    /// that one when it is better (valid where the other is not, then fewer
    /// facet faces, demotions and faces); the refinement's counts still
    /// describe what it did. False without `mRefineSegmentation`.
    bool mRefinementKept = false;
    /// The refinement's Lloyd steps kept, splits and merges.
    std::int64_t mRefineIterations = 0;
    std::int64_t mNumRefineSplits = 0;
    std::int64_t mNumRefineMerges = 0;
    /// The area share primitives explain before and after the refinement.
    double mInitialExplainedAreaFraction = 0.0;
    double mExplainedAreaFraction = 0.0;
    /// Of the B-spline faces, those the general freeform fit built (v17.11.0).
    std::int64_t mNumFreeformGeneral = 0;
    /// Their knot insertion rounds, summed.
    std::int64_t mNumFreeformRefinements = 0;
    /// Their mean pole count (0 without any).
    double mAvgFreeformPoles = 0.0;
    /// Exact edges found by intersecting the surfaces `mExtendSurfaces`
    /// extended, where the surfaces themselves gave no branch (v17.11.0).
    std::int64_t mNumExtendedEdges = 0;
    /// B-rep vertices `mSnapJunctions` placed by its least-squares solve.
    std::int64_t mNumSnappedJunctions = 0;
    /// With `mFaceDeviation`: the two-sided Hausdorff distance between the
    /// re-meshed solid and the skin, over the re-mesh's vertices and triangle
    /// centres and the skin's vertices and facet centres, and the root mean
    /// square of those distances; NaN when it was not measured.
    double mHausdorffMax = std::numeric_limits<double>::quiet_NaN();
    double mHausdorffRms = std::numeric_limits<double>::quiet_NaN();
    /// Whether every shell is closed: every edge borders two faces
    /// (v17.11.0). The native backend reports its writer's closed-loop check.
    bool mClosed = false;
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

/**
 * @brief The OpenCASCADE version this build links (`"7.6.3"`, its
 * `OCC_VERSION_COMPLETE`), or `""` in a build without it. The backend's STEP
 * bytes depend on it, so golden digests are keyed by it (v17.11.0).
 */
MESHIOPLUSPLUS_API const char* fit_cad_occt_version() noexcept;

}  // namespace meshioplusplus
