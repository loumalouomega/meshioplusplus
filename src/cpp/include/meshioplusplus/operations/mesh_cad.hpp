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
//  Provenance:      built from the published descriptions only: Shewchuk
//                   (Discrete Comput. Geom. 1997) for the robust predicates,
//                   Lawson (1977) and Anglada (Computers & Graphics 1997) for
//                   the constrained Delaunay triangulation; Ruppert (1995),
//                   Chew (1989, 1993), Shewchuk (2002), Miller, Pav and
//                   Walkington (2003), Boissonnat and Oudot (2005), Field
//                   (1988) and Chen and Xu (2004) for refinement and
//                   smoothing. No code from an
//                   implementation with an incompatible licence (Triangle,
//                   TetGen, Gmsh, CGAL) or from the CDT library is read or
//                   vendored. See doc/mesh_cad.md.
//
#pragma once

/**
 * @file operations/mesh_cad.hpp
 * @brief Mesh generation from CAD: a conforming triangulation of a B-rep's
 * faces at a target edge size (`mesh_cad_surface`), and the planar
 * constrained Delaunay triangulation it is built on (`triangulate_polygon_2d`).
 *
 * `mesh_cad_surface` takes the `BrepInfo` a CAD JSON, STEP or IGES read
 * returns (Python's `mesh.info["brep"]`) and needs no CAD kernel:
 *
 *  - each B-rep edge is sampled **once**, by arc length into pieces of about
 *    the target size (and under the optional chordal cap), and every face that
 *    uses it, or the same face across a seam, reuses those 3-D points: the
 *    faces conform by construction;
 *  - corners come from the topology (a trim's end is the next trim's start,
 *    an edge's ends are the same at each use), not from comparing
 *    coordinates; a closed edge whose uses start at different points is
 *    rotated onto one corner where a use is a loop of its own;
 *  - boundary pieces whose samples cross in a face's parameters at this size
 *    are halved by arc length, in every face that uses their edge, until
 *    they do not;
 *  - each face is triangulated in its own parameters by a constrained
 *    Delaunay triangulation with exact predicates, the inside an equilateral
 *    lattice scaled by the surface's first derivatives, and lifted onto the
 *    surface; triangles follow the surface normal, reversed for a face with a
 *    swapped normal. A face's embedded edges are constrained segments of its
 *    triangulation, so the mesh carries them as edges.
 *
 * A face that cannot be meshed (its loops meet where the curves themselves
 * meet, a loop does not close, its surface is invalid, an edge cannot be
 * carried to it) is left out and reported by index and reason, never meshed
 * wrong; an embedded edge that cannot be inserted is left out with a warning; `mStrict` turns that
 * into an exception.
 *
 * With `mMinAngleDeg` > 0 each face is refined by Ruppert's algorithm with
 * Chew's circumcentres in its plane (J. Ruppert, J. Algorithms 1995; L. P.
 * Chew, 1989 and SoCG 1993), every triangle judged on the surface, the
 * point found in the plane and lifted (Chew 1993; J.-D. Boissonnat and S.
 * Oudot, Graphical Models 2005); a boundary piece a face needs split is
 * split by arc length in every face that uses its edge, and those faces are
 * meshed again, until none asks. Bad triangles are taken worst first, ties
 * broken by vertex ids (J. R. Shewchuk, Comput. Geom. 2002). Corners below
 * 60 degrees split on concentric shells and the triangles they force are
 * left and reported with the corner (Shewchuk 2002; G. L. Miller, S. E. Pav
 * and N. J. Walkington, IMR 2003), as are poles without `mPoleRings`.
 * `mSmoothIterations` moves interior points toward their star's centroid
 * (Laplacian) or its optimal Delaunay place (L. Chen and J. Xu, J. Comput.
 * Math. 2004), each move kept only when no triangle inverts and the star's
 * smallest angle on the surface does not drop (D. A. Field, Commun. Appl.
 * Numer. Methods 1988). `mMinAngleDeg` = 0 with no smoothing is the first
 * milestone's mesh byte for byte. `mOrder` 2 or 3 writes curved `triangle6`
 * or `triangle10` elements whose nodes lie on the surface and the edges; one
 * whose Jacobian is not positive everywhere (its Bezier bound; A. Johnen,
 * J.-F. Remacle and C. Geuzaine, J. Comput. Phys. 2013) is written
 * straight-sided and counted. `mCurvatureDeviation` and `mCellsAcrossGaps`
 * size the mesh by one gradation-limited field over the whole model (P.-O.
 * Persson, Eng. Comput. 2006), measured on a first uniform mesh, so edges and
 * the faces that share them agree on it; `mPoleRings` grades the triangles
 * around a pole. There is no anisotropic metric yet: its options are in the
 * layout and refused until they land. The output is byte-identical across
 * parallel backends and thread counts. See doc/mesh_cad.md.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell data: the index of the B-rep face a triangle meshes (-1 on edge cells).
inline constexpr const char* kMeshCadFaceName = "cad:face";
/// Cell data: the index of the B-rep edge a line cell discretises (-1 on
/// triangles, and on an edge no `BrepEdge` names: an unreferenced trim or
/// embedded edge, or the boundary of an untrimmed face).
inline constexpr const char* kMeshCadEdgeName = "cad:edge";
/// The prefix of the per-face cell regions: `face:<k>`.
inline constexpr const char* kMeshCadFaceRegionPrefix = "face:";

/// Why a face was left out; explicit values (mirrored by the C enum).
enum class MeshCadFaceStatus : std::int32_t {
    Ok = 0,
    Empty = 1,             ///< its loops enclose nothing
    OpenLoop = 2,          ///< its loops do not bound a domain consistently
    SelfIntersecting = 3,  ///< its boundary meets itself in the parameters, even after splitting
    Degenerate = 4,        ///< distinct boundary samples coincide in the parameters, or a
                           ///< lifted triangle folds over (its normal turns more than 120
                           ///< degrees from the surface's)
    ProjectionFailed = 5,  ///< an edge could not be carried to this face's trim
    BudgetExceeded = 6,    ///< more points than `mMaxPoints`
    Unsupported = 7,       ///< an invalid surface or trim
    Other = 8,
};

/// `"ok"`, `"empty"`, `"open_loop"`, `"self_intersecting"`, `"degenerate"`,
/// `"projection_failed"`, `"budget_exceeded"`, `"unsupported"`, `"other"`.
MESHIOPLUSPLUS_API const char* mesh_cad_face_status_name(MeshCadFaceStatus Status);

/// How smoothing moves an interior point (mirrored by the C enum).
enum class MeshCadSmoothing : std::int32_t {
    Laplacian = 0,  ///< toward the centroid of its neighbours
    Odt = 1,        ///< to its star's area-weighted circumcentre (Chen and Xu 2004)
};

/// Where refinement puts a bad triangle's point (roadmap 6.2.9: only the
/// circumcentre yet; `Offcentre` is refused).
enum class MeshCadInsertion : std::int32_t {
    Circumcentre = 0,
    Offcentre = 1,
};

/// The metric a face is meshed in: `Mean`, one affine map of its
/// parameters per face; `Varying` (roadmap 6.4.1), the surface's own metric
/// wherever Lawson's flips and refinement's encroachment are decided.
enum class MeshCadMetric : std::int32_t {
    Mean = 0,
    Varying = 1,
};

/// `"laplacian"`, `"odt"`.
MESHIOPLUSPLUS_API const char* mesh_cad_smoothing_name(MeshCadSmoothing Method);

struct MeshCadSurfaceOptions {
    /// The target edge length `h` (required, > 0).
    double mTargetSize = 0.0;
    /// The largest distance a boundary segment may sit from its edge (its
    /// sagitta; OpenCASCADE's linear deflection); 0 is off.
    double mMaxChordalDeviation = 0.0;
    /// At least this many points inside every edge.
    std::int64_t mMinPointsPerEdge = 0;
    /// Shorter edges are degenerate (a pole) and unreferenced trims closer
    /// than this are paired; 0 is `1e-7` times the diagonal of the control
    /// nets' box.
    double mGeometricTolerance = 0.0;
    /// Record `cad:face` cell data and a `face:<k>` cell region per face.
    bool mRecordFaceIds = true;
    /// Add a `line` block of the discretised edges, with `cad:edge`.
    bool mEdges = true;
    /// Pair trims no edge references whose 3-D images coincide (ends and
    /// middle within the tolerance), for models without edge topology.
    bool mMatchUnreferencedTrims = true;
    /// Throw (`std::runtime_error`) when any face fails, instead of reporting it.
    bool mStrict = false;
    /// The most points in the result; a face or edge that would exceed it fails.
    std::int64_t mMaxPoints = 50000000;

    // Quality (roadmap 6.2), the fields 6.3 and 6.4 will use, and curved
    // meshes' order (6.6), laid out together so they took one ABI bump; a
    // field whose item has not landed is refused (std::invalid_argument)
    // unless left at its default.

    /// Refine until every triangle's smallest angle on the surface is at
    /// least this (degrees, under 60); 0 keeps the unrefined mesh. Below
    /// 20.7 degrees it is guaranteed away from corners under 60 degrees;
    /// above, a target whose shortfall is reported.
    double mMinAngleDeg = 0.0;
    /// Passes of guarded smoothing over the interior points; 0 is none.
    std::int64_t mSmoothIterations = 0;
    MeshCadSmoothing mSmoothMethod = MeshCadSmoothing::Laplacian;
    MeshCadInsertion mInsertion = MeshCadInsertion::Circumcentre;  ///< 6.2.9
    /// How fast the size may grow, at most `grading - 1` per unit length
    /// (6.3.3); 0 is 1.5, used only when a size source below is on.
    double mGrading = 0.0;
    /// Size curved regions so a triangle sits at most this from the surface
    /// (its chordal error, 6.3.2); 0 is off.
    double mCurvatureDeviation = 0.0;
    /// Put at least this many elements across a thin wall or a gap (6.3.7);
    /// 0 is off.
    double mCellsAcrossGaps = 0.0;
    MeshCadMetric mMetric = MeshCadMetric::Mean;  ///< 6.4.1
    /// Rings of points around each pole, about `h` apart, so its triangles
    /// are graded instead of a thin fan (6.4.2; at most 8); 0 is the fan.
    std::int32_t mPoleRings = 0;
    /// Mesh a face whose metric's condition number over its coarse mesh
    /// exceeds this in a flattening of that mesh instead of its parameters
    /// (6.4.4; 0 is off, otherwise above 1; order 1 only).
    double mReparamConditionNumber = 0.0;
    /// 1: linear triangles; 2 and 3: curved `triangle6` and `triangle10`
    /// (and `line3`, `line4`), their nodes on the surface and edges (6.6).
    std::int32_t mOrder = 1;
};

/// A corner that keeps triangles under `mMinAngleDeg`: two trims meeting
/// at less than 60 degrees in the face's plane, or a pole.
struct MeshCadCorner {
    std::int64_t mFace = -1;
    std::int64_t mPoint = -1;  ///< the mesh point at it
    double mU = 0.0;           ///< its parameters on the face
    double mV = 0.0;
    std::array<double, 3> mX{};
    double mAngleDeg = 0.0;  ///< the angle in the face's plane (0 at a pole)
    /// The smallest angle, on the surface, of the triangles it forces.
    double mMinForcedDeg = 180.0;
    /// A pole's reach on the surface: h, or farther where the face's metric
    /// departs more than twofold from its plane's; the triangles with a
    /// point within it are left to it. 0 for a small angle.
    double mReach = 0.0;
    bool mPole = false;
};

/// One face's outcome.
struct MeshCadFaceRecord {
    std::int64_t mFace = -1;  ///< index into `BrepInfo::mFaces`
    std::int64_t mId = -1;    ///< its `BrepFace::mId`
    MeshCadFaceStatus mStatus = MeshCadFaceStatus::Ok;
    std::string mReason;  ///< empty when meshed
    std::int64_t mNumTriangles = 0;
    std::int64_t mNumInteriorPoints = 0;
    std::int64_t mNumDegenerateDropped = 0;  ///< triangles a pole collapsed
    double mMinAngleDeg = 0.0;
    double mMaxAngleDeg = 0.0;
    double mMinEdge = 0.0;
    double mMaxEdge = 0.0;
    std::int64_t mNumRefinementPoints = 0;  ///< inside it, from refinement
    std::int64_t mNumSmoothed = 0;          ///< moves smoothing kept
    /// Its smallest corner under 60 degrees in its plane; 180 for none.
    double mMinInputAngleDeg = 180.0;
    /// Some triangle is under `mMinAngleDeg` (refinement asked for).
    bool mBelowBound = false;
    /// The smallest angle of its triangles that no corner forces: at least
    /// `mMinAngleDeg` when refinement reached its bound there.
    double mBoundDeg = 0.0;
    std::vector<MeshCadCorner> mSmallAngleCorners;
};

/// The angle histogram's bins: 10 degrees each, 0 to 180.
inline constexpr std::size_t kMeshCadAngleBins = 18;

struct MeshCadSurfaceResult {
    /// Points, one `triangle` block (and a `line` block of the edges).
    Mesh mMesh;
    std::int64_t mNumPoints = 0;
    std::int64_t mNumTriangles = 0;
    std::int64_t mNumEdgeSegments = 0;
    std::int64_t mNumFaces = 0;
    std::int64_t mNumFailed = 0;
    std::int64_t mNumEdges = 0;  ///< edges discretised (named and implicit)
    /// Triangle edges used once (a hole, or the model's own boundary), more
    /// than twice, and twice the same way round.
    std::int64_t mNumOpenEdges = 0;
    std::int64_t mNumNonManifoldEdges = 0;
    std::int64_t mNumMisorientedEdges = 0;
    std::int64_t mNumDegenerateDropped = 0;
    double mTolerance = 0.0;  ///< the geometric tolerance used
    double mMinAngleDeg = 0.0;
    double mMaxAngleDeg = 0.0;
    double mMinEdge = 0.0;
    double mMaxEdge = 0.0;
    /// The largest distance between two uses of one edge or corner: how far
    /// the faces' own trims disagree (the mesh shares the points regardless).
    double mMaxEdgeGap = 0.0;
    std::vector<MeshCadFaceRecord> mFaces;
    std::vector<std::string> mWarnings;
    /// Points refinement added (inside faces and on edges).
    std::int64_t mNumRefinementPoints = 0;
    std::int64_t mNumFacesBelowBound = 0;
    std::int64_t mNumSmoothed = 0;
    /// Every angle of every triangle on the surface, in 10-degree bins.
    std::array<std::int64_t, kMeshCadAngleBins> mAngleHistogram{};
    /// Curved elements whose Jacobian was not positive everywhere, written
    /// straight-sided instead (6.6, 6.6.1).
    std::int64_t mNumCurvedStraightened = 0;
};

/**
 * @brief Triangulate a B-rep's faces into one conforming surface mesh.
 * @throws std::invalid_argument for invalid options (`mTargetSize` not
 *         positive, a negative deviation, tolerance or count) or a B-rep with
 *         no faces; std::runtime_error naming the first failed face when
 *         `mStrict` is set.
 */
MESHIOPLUSPLUS_API MeshCadSurfaceResult
mesh_cad_surface(const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions = {});

/// Sizes a user gives (roadmap 6.3.1), beside the options: each joins the
/// size field the curvature and gap sizes make, the smallest at a point
/// winning and the field graded by `mGrading`. They only refine: the field
/// never exceeds `mTargetSize`.
struct MeshCadSizeSources {
    /// (index into `BrepInfo::mFaces`, size): every point of the face.
    std::vector<std::pair<std::int64_t, double>> mFaceSizes;
    /// (index into `BrepInfo::mEdges`, size): every point of the edge.
    std::vector<std::pair<std::int64_t, double>> mEdgeSizes;
    /// A background mesh (surface cells; a volume mesh is skinned) whose
    /// one-component point_data `mBackgroundArray` is the size at its
    /// points, read at the nearest point of the background. Borrowed.
    const Mesh* pBackground = nullptr;
    std::string mBackgroundArray = "size";
    /// Virtual edges (roadmap 6.5.2): suppress a corner between a sliver of
    /// an edge and its smooth continuation, so no face's loop and no line
    /// cell goes through it, and report the chains (doc/mesh_cad.md, item
    /// 17). Order 1 only.
    bool mVirtualEdges = false;
};

/// What the sources decided, filled when passed (roadmap 6.3.8, 6.5.2).
struct MeshCadSizeReport {
    /// Per face: the `BrepEdge` indices whose sampling is finer than the
    /// face's given size allows (its pieces under the size over the field's
    /// slack), so the face came out finer near them. Empty for a face given
    /// no size.
    std::vector<std::vector<std::int64_t>> mLimitedByEdges;
    /// Per virtual edge: the `BrepEdge` indices of its chain, in order.
    std::vector<std::vector<std::int64_t>> mVirtualEdges;
};

/**
 * @brief `mesh_cad_surface` with sizes given per face, per edge or by a
 * background mesh (roadmap 6.3.1). Empty sources reproduce the two-argument
 * overload byte for byte.
 * @throws std::invalid_argument as that overload, and for a face or edge
 *         index out of range, a size that is not finite and positive, or a
 *         background without a one-component `mBackgroundArray`.
 */
MESHIOPLUSPLUS_API MeshCadSurfaceResult mesh_cad_surface(const BrepInfo& rBrep,
                                                         const MeshCadSurfaceOptions& rOptions,
                                                         const MeshCadSizeSources& rSources,
                                                         MeshCadSizeReport* pReport = nullptr);

/// Inputs beside the options and the sources (v17.20.0), whose layouts have
/// no room left: each off by default, and with all off the mesh is the
/// four-argument overload's byte for byte.
struct MeshCadSurfaceExtras {
    /// Accept the gaps the model itself declares (roadmap 6.5.1): an edge's
    /// `"tolerance"` and its ends' `"end_tolerances"` in its `mExtra` JSON
    /// (a STEP or IGES read writes OpenCASCADE's) widen, for that edge and
    /// its corners, where uses are taken to start at one point and how far
    /// apart they may lie.
    bool mModelTolerances = false;
    /// Merge a sliver face, narrower than `h` and with no interior point, into
    /// the neighbour it meets tangentially along the most edge, through their
    /// common discrete parameterisation (roadmap 6.5.2). Order 1 only.
    bool mVirtualFaces = false;
    /// Corner lopping (roadmap 6.2.3; S. E. Pav and N. J. Walkington,
    /// "Delaunay refinement by corner lopping", IMR 2005) in place of
    /// concentric shells: each acute corner of a face gets a ball of radius
    /// (sqrt 2 - 1) times its distance to the nearest other corner, its
    /// edges split where they leave it, shrunk by half while a vertex lies
    /// inside, and refinement leaves the triangles whose point falls in it.
    /// Only with refinement (`mMinAngleDeg` or a size field).
    bool mLopCorners = false;
};

/// A face whose metric's condition number passed
/// `MeshCadSurfaceOptions::mReparamConditionNumber` (roadmap 6.4.4).
struct MeshCadReparamFace {
    std::int64_t mFace = -1;
    /// Over its coarse mesh: the largest stretch of the surface against the
    /// face's plane over the smallest.
    double mConditionNumber = 0.0;
    /// The same in the flattening's plane (0 when none was built).
    double mConditionNumberDisc = 0.0;
    /// Meshed in the flattening; otherwise `mReason` says why not.
    bool mTaken = false;
    /// 1 least-squares conformal, 2 mean-value on a circle, 0 none.
    std::int32_t mMethod = 0;
    std::string mReason;
};

/// A sliver face found (roadmap 6.5.2) and what became of it.
struct MeshCadMergedFace {
    std::int64_t mSliver = -1;  ///< index into `BrepInfo::mFaces`
    std::int64_t mHost = -1;    ///< the face it is meshed with; -1 for none
    double mWidth = 0.0;        ///< twice its area over its perimeter
    /// The `BrepEdge` indices the two share, left out of the triangulation
    /// and of the line block.
    std::vector<std::int64_t> mEdges;
    bool mMerged = false;
    std::string mReason;  ///< why not, when not merged
};

/// A place a model's own tolerance accepted what the global one refused
/// (roadmap 6.5.1).
struct MeshCadToleranceUse {
    std::int64_t mEdge = -1;  ///< `BrepEdge` index
    double mDeclared = 0.0;   ///< the tolerance the model declares there
    double mGap = 0.0;        ///< the gap it accepted
    /// `"closed_start"`, `"corner"` or `"projection"`.
    std::string mWhere;
};

/// A corner lopped (`MeshCadSurfaceExtras::mLopCorners`, roadmap 6.2.3) and
/// what Pav and Walkington (2005) state of it.
struct MeshCadLoppedCorner {
    std::int64_t mFace = -1;
    std::int64_t mPoint = -1;  ///< the corner's mesh point
    /// Its smallest acute wedge in the face's plane (degrees).
    double mAngleDeg = 0.0;
    /// The distance (in 3-D) to the nearest other corner of the face.
    double mLfs0 = 0.0;
    /// Its ball's radius at the end.
    double mLopRadius = 0.0;
    /// The most any boundary piece of the face turns (degrees).
    double mTurningDeg = 0.0;
    /// Corollary 3's angle bound for that turning, on a planar face; 0 on a
    /// curved one.
    double mStatedBoundDeg = 0.0;
    /// The triangles under `mMinAngleDeg` left in its ball, and the
    /// smallest angle among them on the surface (180 for none).
    std::int64_t mNumForced = 0;
    double mMinForcedDeg = 180.0;
    /// Every one of them within twice the radius of the corner (Lemma 7).
    bool mWithinTwoRadii = true;
    /// The face is planar, `mMinAngleDeg` is under the stated bound,
    /// refinement finished and Lemma 7 held; otherwise `mReason` says why not.
    bool mGuaranteed = false;
    std::string mReason;
};

/// What the extras decided, filled when passed.
struct MeshCadSurfaceExtrasReport {
    std::vector<MeshCadReparamFace> mReparametrised;
    std::vector<MeshCadMergedFace> mMergedFaces;
    std::vector<MeshCadToleranceUse> mToleranceUses;
    std::vector<MeshCadLoppedCorner> mLoppedCorners;
};

/**
 * @brief `mesh_cad_surface` with the extras of v17.20.0 (roadmap 6.5.1,
 * 6.5.2, 6.2.3) beside the sizes, and the reparametrised faces (6.4.4) in its
 * report.
 * @throws std::invalid_argument as the four-argument overload, and for
 *         `mVirtualFaces` with `mOrder` > 1.
 */
MESHIOPLUSPLUS_API MeshCadSurfaceResult mesh_cad_surface(
    const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions,
    const MeshCadSizeSources& rSources, const MeshCadSurfaceExtras& rExtras,
    MeshCadSizeReport* pReport = nullptr, MeshCadSurfaceExtrasReport* pExtrasReport = nullptr);

/**
 * @brief The constrained Delaunay triangulation of a polygon with holes, as a
 * `triangle` mesh in the plane z = 0.
 * @param rLoops the outer boundary first, then the holes; each a closed
 *        polyline whose last point is not repeated, in either orientation
 *        (the even-odd rule decides the domain).
 * @param TargetSize 0 keeps the boundary as given and adds no interior point;
 *        `h > 0` splits boundary segments longer than `h` and fills the inside
 *        with an equilateral lattice of spacing `h`.
 * @throws std::invalid_argument for a negative or non-finite size, a
 *         non-finite point, or a boundary that is degenerate, open or crosses
 *         itself.
 */
MESHIOPLUSPLUS_API Mesh triangulate_polygon_2d(
    const std::vector<std::vector<std::array<double, 2>>>& rLoops, double TargetSize = 0.0);

struct TriangulatePolygon2dOptions {
    double mTargetSize = 0.0;  ///< as `TargetSize` above
    /// Refine (Ruppert, splitting the boundary) until every angle is at
    /// least this (degrees, under 60); 0 keeps the unrefined triangulation.
    double mMinAngleDeg = 0.0;
    std::int64_t mMaxPoints = 1000000;  ///< the most points refinement adds
};

/**
 * @brief `triangulate_polygon_2d` with a quality bound: the same, refined by
 * Ruppert's algorithm when `mMinAngleDeg` > 0. A corner under 60 degrees
 * keeps the triangles it forces.
 * @throws std::invalid_argument as above, or for `mMinAngleDeg` outside
 *         [0, 60) or `mMaxPoints` < 0.
 */
MESHIOPLUSPLUS_API Mesh
triangulate_polygon_2d(const std::vector<std::vector<std::array<double, 2>>>& rLoops,
                       const TriangulatePolygon2dOptions& rOptions);

}  // namespace meshioplusplus
