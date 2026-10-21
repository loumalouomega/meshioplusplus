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
 * @file operations/refine_patches.hpp
 * @brief Refine a segmentation by variational shape approximation: move patch
 * boundaries, split and merge patches until each is explained by one surface,
 * so `fit_primitives` and `fit_cad` find more real CAD faces.
 *
 * `segment` draws its boundaries once; a fillet strip claimed by the plane
 * beside it, a cylinder cut in two at a smooth edge, or a freeform patch with
 * a cylinder inside it stays that way. `refine_patches` starts from the
 * `segment:patch` labels and iterates a hybrid variational shape
 * approximation (Cohen-Steiner, Alliez and Desbrun, SIGGRAPH 2004; with the
 * primitive proxies of Wu and Kobbelt, Computer Graphics Forum 2005):
 *
 *  - **Proxies.** Each patch is explained by a plane, sphere, cylinder, cone
 *    or torus -- `fit_primitives`' fit with every point within the tolerance
 *    `t` -- or else by a *freeform* proxy, a quadric height field over its
 *    best-fit plane (no primitive explains it) or, for an isolated patch too
 *    small to fit, a *facets* proxy.
 *  - **Lloyd step.** Each patch keeps the cell its proxy explains best as its
 *    seed; the patches then grow back from their seeds together, one cell at
 *    a time, cheapest first, through a single priority queue ordered by
 *    `(cost, patch, cell)`. A cell's cost against a proxy is its area times
 *    `d² / t² + w (1 - |n · n_p|)`, `d` its centroid's distance to the proxy,
 *    `n` its normal, `n_p` the proxy's normal there and `w` `mNormalWeight`.
 *    Growth never crosses a feature edge (dihedral angle above
 *    `mFeatureAngleDeg`, as `segment`), an open, non-manifold or
 *    inconsistently wound edge, or a cell outside every patch (`-1`), so a
 *    patch stays edge-connected and inside the smooth area it started in. The
 *    proxies are refitted; if the energy (the sum of the costs) rose, the
 *    step is undone, so a Lloyd step never raises the energy.
 *  - **Split** (`mAllowSplit`). A freeform patch is searched for a primitive
 *    inside it by seeded RANSAC (Schnabel, Wahl and Klein, Computer Graphics
 *    Forum 2007): hypotheses fitted to the neighbourhood of cells drawn by a
 *    splitmix64 generator seeded by `mRansacSeed`, each scored by the area of
 *    the connected cells it explains within `t`. The best one becomes its own
 *    patch and the rest of the freeform patch its connected pieces.
 *  - **Merge** (`mAllowMerge`). Two neighbouring patches one primitive
 *    explains together are merged, cheapest first (the union's RMS), each
 *    patch in at most one merge per iteration.
 *  - **Absorb.** A freeform or facets patch smaller than `mMinPatchCells`
 *    cells joins the neighbour it shares the most smooth edges with (ties to
 *    the lowest patch). A small patch a primitive explains stays.
 *
 * A split is kept when the cells the primitive takes cost less under it than
 * under the freeform proxy they leave, and a merge when the union costs no
 * more than the two patches; absorbing is not an energy move.
 *
 * Refinement stops when the partition no longer changes (`EnergyRose` when the
 * Lloyd step was undone and no move applied), when an iteration adds less
 * than `mEpsilon` to the share of the area primitives explain, or after
 * `mMaxIterations` iterations.
 *
 * ### Output
 *
 * The input -- points, connectivity, block order, data and regions
 * bit-identical -- with:
 *
 *  - `segment:patch`, `segment:adjacency` and the `<prefix>:<k>` cell regions
 *    rewritten for the refined patches, numbered in the order of each patch's
 *    lowest cell (so a partition refinement leaves alone comes back exactly as
 *    it went in);
 *  - the `fit:*` cell and field data of `fit_primitives` run on the refined
 *    patches at the absolute tolerance `t`, so `fit_cad` reads it unchanged;
 *  - field data, one row per patch in label order: `refine:patch_id`
 *    (Int64), `refine:patch_proxy` (Int64, a `RefineProxy`),
 *    `refine:patch_iterations` (Int64, the iterations whose step changed the
 *    patch), `refine:patch_rms` and `refine:patch_max` (Float64, the
 *    distances of its points from its proxy), `refine:patch_cells` (Int64)
 *    and `refine:patch_demotion` (Int64, a `RefineDemotion`: why it has no
 *    primitive);
 *  - field data `refine:explained_area_fraction` (Float64, one value) and
 *    `refine:hausdorff` (Float64, two values: the largest and the RMS distance
 *    of a skin point from its patch's proxy).
 *
 * ### Determinism
 *
 * The proxies of the patches a step changed are refitted in parallel, each
 * serially and stored by patch; the flooding, the moves and every sum run
 * serially in a fixed order, and the RANSAC draws come from the seeded
 * splitmix64 (never a standard-library distribution, whose output differs
 * between standard libraries). The output is byte-identical on every
 * backend and thread count. There is no numpy twin.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"

namespace meshioplusplus {

/// The proxy that explains a refined patch (`refine:patch_proxy`).
enum class RefineProxy : std::int32_t {
    Facets = 0,  ///< an isolated patch too small to fit anything
    Plane = 1,
    Sphere = 2,
    Cylinder = 3,
    Cone = 4,
    Torus = 5,
    Freeform = 6,  ///< no primitive within the tolerance: a quadric height field
};

/// Why a refined patch has no primitive (`refine:patch_demotion`).
enum class RefineDemotion : std::int32_t {
    None = 0,                        ///< a primitive explains it
    TooSmall = 1,                    ///< fewer cells than a fit needs, no neighbour to join
    NoPrimitiveWithinTolerance = 2,  ///< every primitive leaves a point beyond the tolerance
};

/// Why refinement stopped.
enum class RefineStop : std::int32_t {
    Unchanged = 0,     ///< the partition did not change
    SmallGain = 1,     ///< the explained area grew by less than `mEpsilon`
    IterationCap = 2,  ///< `mMaxIterations` iterations ran
    EnergyRose = 3,    ///< the Lloyd step would have raised the energy (undone), no move applied
};

/// Field data names `refine_patches` adds.
inline constexpr const char* kRefinePatchIdName = "refine:patch_id";
inline constexpr const char* kRefinePatchProxyName = "refine:patch_proxy";
inline constexpr const char* kRefinePatchIterationsName = "refine:patch_iterations";
inline constexpr const char* kRefinePatchRmsName = "refine:patch_rms";
inline constexpr const char* kRefinePatchMaxName = "refine:patch_max";
inline constexpr const char* kRefinePatchCellsName = "refine:patch_cells";
inline constexpr const char* kRefinePatchDemotionName = "refine:patch_demotion";
inline constexpr const char* kRefineExplainedAreaName = "refine:explained_area_fraction";
inline constexpr const char* kRefineHausdorffName = "refine:hausdorff";

/// The `mRansacSeed` default: the bytes of "meshio++".
inline constexpr std::uint64_t kRefineDefaultSeed = 0x6d657368696f2b2bull;

/// Options for `refine_patches`.
struct RefinePatchesOptions {
    /// The tolerance `t` a proxy must explain every point of its patch within;
    /// a fraction of the surface's bounding-box diagonal when `mRelative`.
    double mTolerance = 1e-3;
    /// The weight `w` of the normal term in a cell's cost; >= 0.
    double mNormalWeight = 1.0;
    /// Stop when an iteration adds less than this to the explained area share; >= 0.
    double mEpsilon = 1e-3;
    /// The dihedral angle in degrees above which an edge is a feature edge no
    /// patch crosses (`segment`'s default); in [0, 180].
    double mFeatureAngleDeg = 30.0;
    /// The most iterations; >= 0 (0 rewrites the outputs of the input partition).
    std::int64_t mMaxIterations = 10;
    /// A freeform or facets patch with fewer cells joins a neighbour; >= 1.
    std::int64_t mMinPatchCells = 3;
    /// The seed of the RANSAC split's generator. Bindings whose numbers are
    /// doubles (R, JavaScript, pipeline JSON) take it exactly below 2^53.
    std::uint64_t mRansacSeed = kRefineDefaultSeed;
    /// The primitives a proxy may be, as `1 << FitPrimitiveType` bits.
    std::uint32_t mTypes = kFitAllTypes;
    /// `mTolerance` is a fraction of the bounding-box diagonal.
    bool mRelative = true;
    bool mAllowSplit = true;
    bool mAllowMerge = true;
    /// The prefix of the `<prefix>:<k>` cell regions, as `segment`'s.
    std::string mPrefix = "patch";
};

/// What `refine_patches` did.
struct RefinePatchesResult {
    Mesh mMesh;
    std::int64_t mNumPatches = 0;
    std::int64_t mNumIterations = 0;  ///< Lloyd steps kept
    std::int64_t mNumSplits = 0;
    std::int64_t mNumMerges = 0;
    std::int64_t mNumAbsorbed = 0;
    std::int64_t mNumPrimitive = 0;  ///< patches a primitive explains
    std::int64_t mNumFreeform = 0;
    std::int64_t mNumFacets = 0;
    double mTolerance = 0.0;              ///< the absolute tolerance applied
    double mExplainedAreaFraction = 0.0;  ///< the area share primitives explain
    double mInitialExplainedAreaFraction = 0.0;
    double mHausdorffMax = 0.0;  ///< the largest distance of a skin point from its proxy
    double mHausdorffRms = 0.0;
    double mInitialEnergy = 0.0;
    double mEnergy = 0.0;
    RefineStop mStop = RefineStop::Unchanged;
};

/**
 * @brief Refine the `segment:patch` labels of a surface (see the file comment).
 * @param rMesh a surface mesh carrying `segment:patch` (from `segment`).
 * @throws std::invalid_argument without `segment:patch`, on an option out of
 *         range, or on a volume block.
 */
MESHIOPLUSPLUS_API RefinePatchesResult refine_patches(const Mesh& rMesh,
                                                      const RefinePatchesOptions& rOptions = {});

}  // namespace meshioplusplus
