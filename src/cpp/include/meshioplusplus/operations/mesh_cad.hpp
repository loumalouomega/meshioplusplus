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
//                   the constrained Delaunay triangulation. No code from an
//                   implementation with an incompatible licence (Triangle,
//                   TetGen, Gmsh, CGAL) or from the CDT library is read or
//                   vendored. See doc/mesh_cad.md.
//
#pragma once

/**
 * @file operations/mesh_cad.hpp
 * @brief Mesh generation from CAD: a conforming triangulation of a B-rep's
 * faces at one global target edge size (`mesh_cad_surface`), and the planar
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
 * into an exception. There are no quality guarantees (no Delaunay
 * refinement), no size field beyond the chordal cap and no anisotropic metric
 * in this first milestone. The output is byte-identical across parallel
 * backends and thread counts. See doc/mesh_cad.md.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
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
    Degenerate = 4,        ///< distinct boundary samples coincide in the parameters
    ProjectionFailed = 5,  ///< an edge could not be carried to this face's trim
    BudgetExceeded = 6,    ///< more points than `mMaxPoints`
    Unsupported = 7,       ///< an invalid surface or trim
    Other = 8,
};

/// `"ok"`, `"empty"`, `"open_loop"`, `"self_intersecting"`, `"degenerate"`,
/// `"projection_failed"`, `"budget_exceeded"`, `"unsupported"`, `"other"`.
MESHIOPLUSPLUS_API const char* mesh_cad_face_status_name(MeshCadFaceStatus Status);

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
};

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

}  // namespace meshioplusplus
