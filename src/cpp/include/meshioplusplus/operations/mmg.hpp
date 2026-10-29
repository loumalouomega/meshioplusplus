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
 * @file mmg.hpp
 * @brief The optional MMG backend (roadmap 3.2): what it is built with.
 *
 * [MMG](https://www.mmgtools.org/) -- `mmg2d` (planar triangles), `mmgs`
 * (surface triangles in 3-D) and `mmg3d` (tetrahedra) -- adapts a simplicial
 * mesh to a size or an anisotropic metric, embeds a level set in it, meshes a
 * planar domain from its boundary and moves a mesh by a displacement. It is
 * LGPL-3.0: built with `-DMESHIOPLUSPLUS_WITH_MMG=ON` against an install
 * (`build/build-mmg.sh`), linked dynamically only, never vendored, and absent
 * from the release wheels and CLI binaries. Without it every MMG operation
 * throws `Unsupported` naming the flag; there is no second implementation,
 * since MMG's result cannot be reproduced by one (the `remesh_volume`
 * precedent). See `doc/mmg.md` for the build, the licence and the embedding
 * hazards (signal handlers, index width, threads) the core handles.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/quality.hpp"

namespace meshioplusplus {

/// Whether this build has the MMG backend.
MESHIOPLUSPLUS_API bool mmg_available() noexcept;

/**
 * @brief Whether the MMG this build links was built with ELAS, the
 *        linear-elasticity library `move` needs (`MMG3D_mmg3dmov`,
 *        `MMG2D_mmg2dmov`). ELAS is GPL-3.0; see `doc/mmg.md`.
 */
MESHIOPLUSPLUS_API bool mmg_has_elas() noexcept;

/// The width in bits of the linked MMG's `MMG5_int` (32 or 64), or 0 without MMG.
MESHIOPLUSPLUS_API int mmg_index_bits() noexcept;

/// How an MMG call ended when it returned a mesh: `MMG5_SUCCESS`, or
/// `MMG5_LOWFAILURE` (a valid mesh MMG could not finish; a warning, not an
/// error). `MMG5_STRONGFAILURE` throws instead.
enum class MmgStatus : std::int32_t { Success = 0, LowFailure = 1 };

/// Sizes for the entities of one named region (`Set_localParameter`); a
/// negative value keeps the global one.
struct MmgLocalParameter {
    std::string mRegion;
    double mHmin = -1.0;
    double mHmax = -1.0;
    double mHausd = -1.0;
};

/**
 * @brief The parameters every MMG operation shares. A negative size, distance
 *        or gradation leaves MMG's own default.
 */
struct MmgParameters {
    double mHmin = -1.0;      ///< minimal edge size
    double mHmax = -1.0;      ///< maximal edge size
    double mHausd = -1.0;     ///< how far the boundary may move (Hausdorff distance)
    double mHgrad = -1.0;     ///< gradation; 0 or less than 1 is rejected by MMG
    double mHgradreq = -1.0;  ///< gradation from required entities
    /// Ridge detection angle in degrees (MMG's default, 45); zero or negative
    /// turns the detection off.
    double mAngleDetection = 45.0;
    bool mOptim = false;     ///< keep the mesh's size, only improve quality
    bool mNoInsert = false;  ///< no point insertion or removal
    bool mNoSwap = false;    ///< no edge or face flip
    bool mNoMove = false;    ///< no point relocation
    bool mNoSurf = false;    ///< leave the boundary untouched (not mmgs)
    bool mNreg = false;      ///< regularise the normals
    bool mOpnbdy = false;    ///< keep internal open boundaries (mmg3d, mmg2d)
    std::int64_t mMem = -1;  ///< MMG's memory budget in MB
    int mVerbose = -1;       ///< MMG's verbosity; -1 is silent
    /// One flag per point (empty: none); a nonzero point is required: MMG
    /// neither moves nor removes it.
    std::vector<std::uint8_t> mFrozen;
    /// Regions whose entities are required: points, cells or sides.
    std::vector<std::string> mRequired;
    /// Per-region sizes.
    std::vector<MmgLocalParameter> mLocal;
};

/// The options of `adapt`.
struct AdaptOptions {
    MmgParameters mParameters;
    /// A constant target edge size (`hsiz`), when no size field or metric is given.
    double mHsiz = -1.0;
    /// A one-component `point_data` array: the target edge size at each point.
    std::string mSize;
    /// A symmetric-tensor `point_data` array, the anisotropic metric: six
    /// components (three in 2-D) in the lower-triangle order of a Medit `.sol`,
    /// `m11 m21 m22 m31 m32 m33` -- the order `metric_from_hessian` writes.
    std::string mMetric;
};

/// What an MMG operation reports beside its mesh.
struct MmgReport {
    MmgStatus mStatus = MmgStatus::Success;
    std::string mLibrary;  ///< "mmg2d", "mmgs" or "mmg3d"
    std::int64_t mNumPointsIn = 0;
    std::int64_t mNumPointsOut = 0;
    std::int64_t mNumElementsIn = 0;  ///< triangles or tetrahedra (and prisms, quads)
    std::int64_t mNumElementsOut = 0;
    /// The `quality:aspect_ratio` summary of the result's elements (1 is ideal).
    QualityMetricSummary mQuality;
    std::vector<std::string> mWarnings;
};

/// The result of `adapt`.
struct AdaptResult {
    Mesh mMesh;
    MmgReport mReport;
};

/**
 * @brief Remesh @p rMesh to a size or a metric with MMG (roadmap 3.2.1).
 *
 * A mesh with tetrahedra goes to `MMG3D_mmg3dlib` (wedges pass as MMG prisms,
 * triangles and quads as boundary faces, lines as ridges), a planar triangle
 * mesh (2-D points, or 3-D ones with one z) to `MMG2D_mmg2dlib` (quads pass,
 * lines are boundary edges) and any other triangle surface to `MMGS_mmgslib`
 * (lines are ridges). The size is `mHsiz`, the `mSize` field or the `mMetric`
 * field, else MMG's own from the mesh.
 *
 * Entities: cell and side regions, and a `medit:ref` cell array, become MMG
 * references and come back under their names (a region is never lost to a
 * renumbering; an empty one stays empty). Point regions survive when their
 * points are required (`mFrozen`, `mRequired`, `medit:RequiredVertices`,
 * `medit:Corners`), matched by their unchanged coordinates; others are dropped
 * with a warning. `medit:Corners` and `medit:Ridges` become MMG corners and
 * ridges. Boundary faces (mmg3d) or edges come back when the input had some.
 * Other point and cell data are dropped with a warning (compose with
 * `interpolate`); `field_data` is copied.
 *
 * @throws Unsupported without the MMG backend, or for a mesh beyond MMG's
 *         `MMG5_int`; std::invalid_argument for hexahedra, pyramids,
 *         higher-order or ragged blocks (`convert_cells` first), a mesh with
 *         no triangle or tetrahedron, a missing or misshapen field or region,
 *         or a frozen mask of the wrong length; std::runtime_error when MMG
 *         fails (`MMG5_STRONGFAILURE`).
 */
MESHIOPLUSPLUS_API AdaptResult adapt(const Mesh& rMesh, const AdaptOptions& rOptions);

/// One material of `discretize_levelset` (MMG's `Set_multiMat`): the cells
/// of `medit:ref` `mRef` are cut by the level set into `mRefInside` (where it
/// is below the isovalue) and `mRefOutside`, or kept whole when not `mSplit`.
struct MmgMaterial {
    std::int64_t mRef = 0;
    bool mSplit = true;
    std::int64_t mRefInside = 3;
    std::int64_t mRefOutside = 2;
};

/// The options of `discretize_levelset`.
struct DiscretizeLevelsetOptions {
    MmgParameters mParameters;
    /// A one-component `point_data` array: the level set, cut at `mIsovalue`.
    std::string mLevelSet;
    double mIsovalue = 0.0;
    /// Remove the parasitic components of the cut smaller than this fraction
    /// of the mesh's volume (MMG's `rmc`); zero or negative keeps them.
    double mRemoveSmallComponents = -1.0;
    /// Several materials (empty: the whole mesh is one, cut into 3 and 2).
    std::vector<MmgMaterial> mMaterials;
};

/// The result of `discretize_levelset`.
struct DiscretizeLevelsetResult {
    Mesh mMesh;
    MmgReport mReport;
};

/**
 * @brief Embed the isovalue of a level set in @p rMesh conformally and remesh
 *        both sides with MMG (roadmap 3.2.2).
 *
 * `MMG3D_mmg3dls` for tetrahedra, `MMG2D_mmg2dls` for a planar triangle mesh,
 * `MMGS_mmgsls` for a surface. The result's elements carry MMG's references
 * as `cell_data["medit:ref"]` -- 3 inside (the level set below the isovalue),
 * 2 outside, or the materials' -- and the cell regions `inside` and
 * `outside`; the interface is the side region `interface`, on the inside
 * elements. Side regions of the input come back as for `adapt`; its cell
 * regions do not (the cut gives the elements new references: name them with
 * materials), and point and cell data are dropped.
 *
 * @throws Unsupported without the MMG backend; std::invalid_argument as
 *         `adapt`, or for a missing or misshapen level set;
 *         std::runtime_error when MMG fails.
 */
MESHIOPLUSPLUS_API DiscretizeLevelsetResult
discretize_levelset(const Mesh& rMesh, const DiscretizeLevelsetOptions& rOptions);

/// The options of `mesh_domain`.
struct MeshDomainOptions {
    MmgParameters mParameters;
    double mHsiz = -1.0;  ///< a constant target edge size
    std::string mSize;    ///< a one-component `point_data` size at the boundary's points
    /// One point inside each hole, `x, y` pairs: the subdomains holding one
    /// are not meshed.
    std::vector<double> mHoles;
};

/// The result of `mesh_domain`.
struct MeshDomainResult {
    Mesh mMesh;
    MmgReport mReport;
    std::int64_t mNumSubdomains = 0;  ///< meshed, holes excluded
    std::int64_t mNumHoles = 0;       ///< subdomains left out
};

/**
 * @brief Triangulate the planar domain bounded by the `line` cells of
 *        @p rBoundary with MMG (`MMG2D_mmg2dmesh`, roadmap 3.2.3).
 *
 * MMG meshes every region the lines enclose and gives each its own
 * reference; it has no notion of a hole, so `mHoles` names the regions to
 * leave out by a point inside each. The result has the boundary as `line`
 * cells, each enclosed region as a cell region `subdomain:<k>` (`k` MMG's
 * reference), and the cell regions of the lines also as side regions of the
 * triangles. The size is `mHsiz`, the `mSize` field or MMG's own from the
 * boundary.
 *
 * @throws Unsupported without the MMG backend; std::invalid_argument for a
 *         boundary with other cells than lines (and vertices), non-planar
 *         points, or a hole point outside the domain; std::runtime_error when
 *         MMG fails.
 */
MESHIOPLUSPLUS_API MeshDomainResult mesh_domain(const Mesh& rBoundary,
                                                const MeshDomainOptions& rOptions);

/// The options of `move`.
struct MoveOptions {
    MmgParameters mParameters;
    /// A `point_data` vector (dim components; 3 or 2 with mmg2d): the displacement.
    std::string mDisplacement;
    /// The cell or side regions whose boundary facets move by it; the rest of
    /// the boundary stays fixed.
    std::vector<std::string> mMoving;
    /// MMG's `lag`: 0 moves only, 1 also swaps edges and faces, 2 also
    /// inserts and collapses points.
    int mLagrangian = 1;
};

/**
 * @brief Move @p rMesh by a displacement of part of its boundary and remesh
 *        as it deforms, with MMG's Lagrangian mode (`MMG3D_mmg3dmov`,
 *        `MMG2D_mmg2dmov`; roadmap 3.2.4).
 *
 * The facets of the `mMoving` regions follow `mDisplacement`; MMG extends it
 * into a layer of elements around them by linear elasticity (ELAS) and keeps
 * the rest of the boundary in place. Tetrahedra or planar triangles only.
 * Regions come back as for `adapt`.
 *
 * @throws Unsupported without the MMG backend, or with an MMG built without
 *         ELAS; std::invalid_argument for a surface mesh, a missing or
 *         misshapen displacement, or no facet in the moving regions;
 *         std::runtime_error when MMG fails.
 */
MESHIOPLUSPLUS_API AdaptResult move(const Mesh& rMesh, const MoveOptions& rOptions);

/// The options of `metric_from_hessian`.
struct MetricFromHessianOptions {
    std::string mHessian;            ///< the `point_data` Hessian: (n, 9) or (n, 4), row-major
    double mError = 1e-2;            ///< the target interpolation error
    double mHmin = 1e-3;             ///< the smallest edge the metric may ask for
    double mHmax = 1.0;              ///< the largest
    std::string mOutput = "metric";  ///< the `point_data` name of the result
};

/**
 * @brief The anisotropic metric of a Hessian, for `adapt` (roadmap 3.2.1.1).
 *
 * Per point, `M = c / error * |H|`: the Hessian symmetrised, its eigenvalues
 * taken in absolute value, scaled by `c = (1/2) (d/(d+1))^2` for P1
 * interpolation in `d` dimensions (Frey and Alauzet 2005) and clamped to
 * `[1/hmax^2, 1/hmin^2]`, so the edge sizes it asks for lie in `[hmin, hmax]`.
 * The result is `point_data[mOutput]`, six components (three for a 2-D
 * Hessian) in the lower-triangle order of a Medit `.sol`. Needs no MMG; feed
 * it [`hessian`](hessian.hpp)'s output. The mesh is otherwise unchanged.
 *
 * @throws std::invalid_argument for a missing Hessian, one of another width,
 *         or a non-positive error or size.
 */
MESHIOPLUSPLUS_API Mesh metric_from_hessian(const Mesh& rMesh,
                                            const MetricFromHessianOptions& rOptions);

}  // namespace meshioplusplus
