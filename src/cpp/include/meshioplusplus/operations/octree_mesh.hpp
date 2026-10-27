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
//  Provenance:      built from the published descriptions only: Sundar,
//                   Sampath and Biros (SIAM J. Sci. Comput. 2008) for the
//                   linear octree, Barill et al. (SIGGRAPH 2018) for the fast
//                   winding number, Labelle and Shewchuk (SIGGRAPH 2007) for
//                   isosurface stuffing, Bronson, Levine and Whitaker (IEEE
//                   TVCG 2014) for lattice cleaving, Livesu, Pitzalis and
//                   Cherchi (ACM TOG 2022) for the octree dual. No code from
//                   an implementation with an incompatible licence (Stellar,
//                   TetGen, snappyHexMesh) is read or vendored. See
//                   doc/octree_mesh.md.
//
#pragma once

/**
 * @file operations/octree_mesh.hpp
 * @brief A volume mesh of the inside of a skin on a balanced octree: castellated
 * hexahedra, conforming graded tetrahedra fitted to the skin (the default), or
 * all-hexahedra.
 *
 * The skin is any surface (`triangle`, `quad`, `polygon` blocks), closed or
 * slightly dirty: inside and outside come from the fast generalised winding
 * number, so holes, self-intersections and non-manifold edges are tolerated.
 * A size field (`h_min` ≤ h ≤ `h_max`, driven by curvature, sharp features,
 * the local thickness and an optional array) refines a linear octree, which
 * is 2:1-balanced across faces, edges and corners. Then:
 *
 *  - `Castellated`: the leaves inside the skin as `hexahedron` cells, with
 *    hanging nodes flagged, the `inside` and `boundary_cells` cell regions,
 *    the `surrogate_boundary` side region (the faces between inside and
 *    outside leaves) and `sdf:distance` / `sdf:closest_point` at every node,
 *    for embedded and Shifted Boundary Method solvers;
 *  - `Tetra`: a conforming graded BCC tetrahedralisation of the octree, cut
 *    to the skin by isosurface stuffing (or cleaving), sharp features
 *    snapped, then fitted, untangled and optimised;
 *  - `Hexa`: all-hexahedra by the octree dual, projected, pillowed and
 *    untangled.
 *
 * `octree` (operations/octree.hpp) is a separate operation and is not
 * affected. Every output is byte-identical across parallel backends and
 * thread counts. See doc/octree_mesh.md.
 */

// System includes
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace meshioplusplus {

/// Point data: the signed distance to the skin (negative inside).
inline constexpr const char* kOctreeMeshDistanceName = "sdf:distance";
/// Point data (3 components): the nearest point on the skin.
inline constexpr const char* kOctreeMeshClosestPointName = "sdf:closest_point";
/// Cell data (Int64): the lattice shape class of a tetrahedron (1-6), 0 for a
/// cell the stuffing or the optimiser created.
inline constexpr const char* kOctreeMeshShapeName = "octree_mesh:shape";
/// Cell data (Int64): the material of a tetrahedron (0 outside the skin, 1
/// inside, higher for nested shells under cleaving), when there is more than one.
inline constexpr const char* kOctreeMeshMaterialName = "octree_mesh:material";
/// Cell data (Int64): 1 for a hexahedron of the pillow layer along the
/// boundary of a hexahedral output, 0 for the grid's own.
inline constexpr const char* kOctreeMeshPillowName = "octree_mesh:pillow";
/// Cell region: the cells inside the skin.
inline constexpr const char* kOctreeMeshInsideRegion = "inside";
/// Cell region: the cells outside the skin (`mKeepOutside`).
inline constexpr const char* kOctreeMeshOutsideRegion = "outside";
/// Cell region: the castellated leaves the skin passes through.
inline constexpr const char* kOctreeMeshBoundaryCellsRegion = "boundary_cells";
/// Side region: the faces between inside and outside castellated leaves.
inline constexpr const char* kOctreeMeshSurrogateRegion = "surrogate_boundary";
/// Side region: the boundary faces of a tetrahedral or hexahedral output.
inline constexpr const char* kOctreeMeshBoundaryRegion = "boundary";
/// Point region: the output nodes on the skin's feature corners.
inline constexpr const char* kOctreeMeshCornersRegion = "feature_corners";
/// Point region: the output nodes on the skin's feature edges.
inline constexpr const char* kOctreeMeshFeatureEdgesRegion = "feature_edges";

/// What `octree_mesh` produces.
enum class OctreeMeshOutput {
    Castellated = 0,  ///< the leaves as hexahedra, with surrogate-boundary tags
    Tetra = 1,        ///< conforming tetrahedra fitted to the skin
    Hexa = 2,         ///< all-hexahedra by the octree dual
};

/// How a tetrahedral output conforms to the skin.
enum class OctreeMeshConformity {
    Stuffing = 0,  ///< isosurface stuffing (Labelle and Shewchuk 2007)
    Cleaving = 1,  ///< lattice cleaving of the winding-number regions (Bronson et al. 2014)
};

/// How the size field may vary along the skin.
enum class OctreeMeshSkinGrading {
    /// The size field grades along the skin too; the tetrahedra cut where the
    /// level changes are transition shapes, and the minimum dihedral is reported.
    Graded = 0,
    /// Every leaf the skin cuts is at one level, with two rings of leaves at
    /// that level around them: every tetrahedron stuffing touches is a BCC
    /// tetrahedron and the 10.7°-164.8° dihedral bound applies.
    Uniform = 1,
};

/// How a hexahedral output grades.
enum class OctreeMeshHexGrading {
    /// One grid at the skin's size everywhere.
    Uniform = 0,
    /// The octree's own grading, made conforming by the dual transition
    /// schemes of Livesu, Pitzalis and Cherchi (2022).
    Graded = 1,
};

/// Which leaves a castellated output keeps.
enum class OctreeMeshFill {
    All = 0,     ///< every leaf
    Inside = 1,  ///< the leaves whose centre is inside
    Cover = 2,   ///< the inside leaves and every leaf the skin cuts
};

/// Parse an `OctreeMeshOutput` name (`castellated`, `tet`, `hex`).
MESHIOPLUSPLUS_API OctreeMeshOutput octree_mesh_output_from_name(const std::string& rName);
/// Parse an `OctreeMeshConformity` name (`stuffing`, `cleaving`).
MESHIOPLUSPLUS_API OctreeMeshConformity octree_mesh_conformity_from_name(const std::string& rName);
/// Parse an `OctreeMeshSkinGrading` name (`graded`, `uniform`).
MESHIOPLUSPLUS_API OctreeMeshSkinGrading
octree_mesh_skin_grading_from_name(const std::string& rName);
/// Parse an `OctreeMeshHexGrading` name (`uniform`, `graded`).
MESHIOPLUSPLUS_API OctreeMeshHexGrading octree_mesh_hex_grading_from_name(const std::string& rName);
/// Parse an `OctreeMeshFill` name (`all`, `inside`, `cover`).
MESHIOPLUSPLUS_API OctreeMeshFill octree_mesh_fill_from_name(const std::string& rName);
/// The name of an `OctreeMeshOutput` (`castellated`, `tet`, `hex`).
MESHIOPLUSPLUS_API const char* octree_mesh_output_name(OctreeMeshOutput Output);

/// Options for `octree_mesh`.
struct OctreeMeshOptions {
    /// What to produce.
    OctreeMeshOutput mOutput = OctreeMeshOutput::Tetra;
    /// The smallest leaf edge; 0: `mMaxSize / 32`.
    double mMinSize = 0.0;
    /// The largest leaf edge; 0: an eighth of the skin's bounding-box diagonal.
    double mMaxSize = 0.0;
    /// Cap on the finest level (at most 19); -1: from the sizes.
    std::int64_t mMaxLevel = -1;
    /// Explicit bounds `{lo[3], hi[3]}`; unset uses the skin's bounding box.
    std::optional<std::array<double, 6>> mBounds;
    /// Padding added to every side, in world units.
    double mPadding = 0.0;
    /// Padding added to every side, as a fraction of the bounding-box diagonal.
    double mPaddingRelative = 0.1;
    /// The angle (degrees) a leaf edge may subtend on a curved skin:
    /// `h <= 2 sin(θ/2) / κ_max`; 0 disables the curvature term.
    double mCurvatureAngle = 15.0;
    /// Edges sharper than this dihedral (degrees) are features; 0 disables
    /// features.
    double mFeatureAngle = 30.0;
    /// A feature vertex where the feature turns by more than this (degrees) is
    /// a corner, as is one where the feature branches or ends.
    double mCornerAngle = 30.0;
    /// The size along feature edges; 0: `mMinSize`.
    double mFeatureSize = 0.0;
    /// Leaves across the local thickness of the skin (an inward ray cast); 0
    /// disables the thickness term.
    double mGapCells = 2.0;
    /// Also measure the gaps outside the skin (between separate parts).
    bool mGapOutside = false;
    /// Name of a `point_data` array of target sizes on the skin; empty: none.
    std::string mSizeArray;
    /// Restrict the skin to this named cell region; empty: every surface block.
    std::string mSurfaceRegion;
    /// Rings of leaves around every cut leaf kept at its level.
    std::int64_t mGradingLayers = 1;
    /// How the size field may vary along the skin.
    OctreeMeshSkinGrading mSkinGrading = OctreeMeshSkinGrading::Graded;
    /// The skin size under `Uniform` grading; 0: the smallest size the field
    /// asks for on the skin.
    double mSurfaceSize = 0.0;
    /// A point is inside where the winding number exceeds this.
    double mWindingThreshold = 0.5;
    /// Accuracy of the fast winding number (Barill et al.'s β); 0 sums every
    /// triangle exactly.
    double mWindingBeta = 2.0;
    /// How a tetrahedral output conforms to the skin.
    OctreeMeshConformity mConformity = OctreeMeshConformity::Stuffing;
    /// Stuffing's violation threshold on long lattice edges.
    double mAlphaLong = 0.24799;
    /// Stuffing's violation threshold on short lattice edges.
    double mAlphaShort = 0.40173;
    /// Snap the boundary to the skin's feature corners and edges.
    bool mPreserveFeatures = true;
    /// Rounds of fitting, untangling and optimisation; 0 disables them.
    std::int64_t mFitIterations = 5;
    /// Undo any optimisation move that lowers a node's minimum dihedral below
    /// the mesh's minimum before the pass.
    bool mDihedralFloorGuard = true;
    /// Which leaves a castellated output keeps.
    OctreeMeshFill mFill = OctreeMeshFill::Inside;
    /// Also mesh the space between the skin and the bounds (tetrahedral
    /// output), as the `outside` region.
    bool mKeepOutside = false;
    /// Attach `sdf:distance` (and, castellated, `sdf:closest_point`).
    bool mAttachDistance = true;
    /// Attach the `quality:*` cell data of `attach_quality`.
    bool mAttachQuality = false;
    /// Hexahedral output: insert a layer of hexahedra along the boundary.
    bool mPillow = true;
    /// Hexahedral output: one grid at the skin's size, or the octree's grading.
    OctreeMeshHexGrading mHexGrading = OctreeMeshHexGrading::Uniform;
    /// Hexahedral output: untangling iterations.
    std::int64_t mUntangleIterations = 50;
    /// Refuse an octree with more leaves than this.
    std::int64_t mMaxCells = 20000000;
    /// Refuse an output with more cells than this.
    std::int64_t mMaxTets = 100000000;
};

/// The result of `octree_mesh`.
struct OctreeMeshResult {
    Mesh mMesh;                                          ///< The volume mesh.
    OctreeMeshOutput mOutput = OctreeMeshOutput::Tetra;  ///< What it is.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};      ///< Lo corner of the root lattice.
    double mRootSize = 0.0;                              ///< Edge of a root cell.
    std::array<std::int64_t, 3> mRootCells{{1, 1, 1}};   ///< Root cells per axis.
    std::int64_t mMaxLevel = 0;                          ///< The finest level allowed.
    std::int64_t mFinestLevel = 0;                       ///< The deepest level present.
    std::int64_t mSurfaceLevel = -1;  ///< The skin level under `Uniform` grading, else -1.
    double mMinSize = 0.0;            ///< The resolved smallest leaf edge.
    double mMaxSize = 0.0;            ///< The resolved largest leaf edge.
    std::int64_t mNumLeaves = 0;      ///< Leaves of the balanced octree.
    std::vector<std::int64_t> mLeavesPerLevel;  ///< Leaves at each level.
    std::int64_t mNumPoints = 0;                ///< Output points.
    std::int64_t mNumCells = 0;                 ///< Output cells.
    std::int64_t mNumHexahedra = 0;             ///< Output hexahedra.
    std::int64_t mNumTetrahedra = 0;            ///< Output tetrahedra.
    std::int64_t mNumInside = 0;                ///< Cells in the `inside` region.
    std::int64_t mNumBoundaryCells = 0;         ///< Castellated cells the skin cuts.
    std::int64_t mNumSurrogateFaces = 0;        ///< Faces in the `surrogate_boundary` region.
    std::int64_t mNumHanging = 0;               ///< Hanging nodes of a castellated output.
    std::int64_t mNumCutPoints = 0;             ///< Points stuffing placed on lattice edges.
    std::int64_t mNumWarped = 0;                ///< Lattice vertices stuffing moved onto the skin.
    std::int64_t mNumWarpRefused = 0;           ///< Warps refused to keep a tetrahedron valid.
    double mMinDihedral = 0.0;                  ///< Smallest dihedral angle (tetrahedra), degrees.
    double mMaxDihedral = 0.0;                  ///< Largest dihedral angle (tetrahedra), degrees.
    double mMinScaledJacobian = 0.0;            ///< Smallest scaled Jacobian over the cells.
    std::int64_t mNumInverted = 0;              ///< Cells with a non-positive volume.
    /// The two-sided distance between the output's boundary and the skin; -1
    /// when not measured (castellated output).
    double mHausdorffToSkin = -1.0;
    std::int64_t mNumBoundaryEdgesOpen = 0;         ///< Output boundary edges used once.
    std::int64_t mNumNonManifoldBoundaryEdges = 0;  ///< Output boundary edges used 3+ times.
    std::int64_t mNumCorners = 0;                   ///< Feature corners of the skin.
    std::int64_t mNumCornersSnapped = 0;            ///< Corners an output node was snapped to.
    std::int64_t mNumFeatureChains = 0;             ///< Feature chains of the skin.
    std::int64_t mNumFeatureChainsMissed = 0;       ///< Chains not recovered by the boundary.
    std::int64_t mNumFeatureMoves = 0;              ///< Boundary nodes moved onto features.
    std::int64_t mNumOptimizedMoves = 0;            ///< Nodes the optimisation moved.
    std::int64_t mNumFlips = 0;                     ///< Flips the optimisation kept.
    /// The Labelle-Shewchuk dihedral bound applies to every tetrahedron:
    /// `Uniform` grading, nothing snapped to features, nothing refused.
    bool mGuaranteeApplies = false;
    std::int64_t mNumWindingRechecks = 0;  ///< Winding numbers re-evaluated near 0.5.
    /// Winding numbers not evaluated: the point lies deep in a region the
    /// closed skin does not cut, whose integer winding number was flooded in.
    std::int64_t mNumWindingSkipped = 0;
    QualityReport mQuality;  ///< `compute_quality` of the output.
    /// Seconds spent in each stage, in order.
    std::vector<std::pair<std::string, double>> mTimings;
    /// Everything that did not meet its target, in order.
    std::vector<std::string> mWarnings;
};

/**
 * @brief Mesh the inside of the skin @p rSkin on a balanced octree.
 * @throws std::invalid_argument on a skin with no surface cells, sizes that are
 *         not positive or `mMinSize > mMaxSize`, a lattice too fine for the
 *         63-bit keys (`root_cells << max_level >= 2^20` on an axis), more
 *         leaves than `mMaxCells`, more output cells than `mMaxTets`, an
 *         unknown `mSizeArray` or `mSurfaceRegion`, and out-of-range options,
 *         each naming the option.
 */
MESHIOPLUSPLUS_API OctreeMeshResult octree_mesh(const Mesh& rSkin,
                                                const OctreeMeshOptions& rOptions = {});

}  // namespace meshioplusplus
