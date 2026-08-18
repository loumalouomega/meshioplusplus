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
 * @file decimate_volume.hpp
 * @brief Volume decimation by quadric-error-metric tet-edge collapse: reduce a
 * tetrahedral mesh's cell count while preserving its shape and boundary.
 *
 * The resolution-*reducing* counterpart of `refine` on 3D cells, and the
 * volume-mesh sibling of surface `decimate` -- a **separate operation**, not a
 * mode on it: `DecimateOptions`/`mio_decimate` are untouched by this file, and
 * `decimate()` keeps throwing by name on any 3D volume block, pointing here.
 *
 * **Tet-only in v1.** Any hexahedron/wedge/pyramid/polyhedron/ragged/
 * higher-order 3D block, or any non-3D block mixed in, throws by name pointing
 * at `convert_cells(Simplexify)`.
 *
 * **Boundary vertices participate**, with a real quadric-error objective,
 * rather than being pinned as `decimate`'s own `mPreserveBoundary` default
 * does -- `mPreserveBoundary` defaults to `false` here; set it to reproduce
 * the simpler pinned-boundary behaviour.
 *
 * ### Objective: one priority key, two regimes
 *
 * Every vertex accumulates a Garland-Heckbert plane quadric from its incident
 * **boundary triangles only** (the mesh's own outer skin, from
 * `detail::build_global_faces`'s faces with no neighbour) -- a purely interior
 * vertex, by construction, accumulates nothing and its quadric is the exact
 * zero matrix. For a candidate edge `(a, b)` the combined quadric `Q(a)+Q(b)`
 * feeds the *same* placement solve `decimate` uses (`detail::decim_place`,
 * hoisted into `detail/decimate_common.hpp`) whether the edge touches the
 * boundary or not:
 *  - at least one endpoint touches the boundary: `Q_ab` is non-degenerate (or
 *    the solve's own ill-conditioning bound already routes it to the
 *    midpoint) -- identical to `decimate`'s own code path;
 *  - both endpoints are purely interior: `Q_ab` is exactly zero, which that
 *    *same* `|det| <= 1e-12 * scale^3` bound classifies as degenerate and
 *    falls back to the midpoint -- no interior-specific placement code is
 *    needed at all.
 *
 * Scoring, however, needs an explicit split: `decim_quadric_error` against an
 * exact-zero quadric is identically zero for every interior-interior edge,
 * which would tie the whole interior of the mesh and fall through to id
 * order. The priority key is therefore `(regime, score, ids...)`, `regime = 0`
 * (`score` = quadric error) whenever `Q_ab` is non-degenerate, else `regime =
 * 1` (`score` = squared edge length) -- every boundary-touching collapse is
 * considered ahead of purely-interior ones. `mMaxError` is only really
 * meaningful for `regime == 0` entries (real quadric-error units); comparing
 * it directly against `regime == 1`'s squared length is a documented rough
 * tool for mixed meshes, not a claimed exact criterion. `mTargetRatio`/
 * `mTargetCells` are regime-agnostic and remain exactly well-defined.
 *
 * ### Validity guards
 *
 * Tet-only, so simpler than the general 3D case: a tet's *other* incident
 * cells are exactly `T(ab) = vtets[a] & vtets[b]` (its shared tets). Each
 * guard *rejects* the individual collapse (counted in `mCollapsesRejected`)
 * rather than aborting:
 *  - the **vertex-link condition**: the exact SET `Vlink(a) & Vlink(b)` (node
 *    adjacency via each vertex's own incident tets) must equal `Vlink(ab)`
 *    (the two "opposite" corners of each tet in `T(ab)`) -- pure integer set
 *    equality, no floating point;
 *  - the **duplicate-tet guard**: no surviving tet incident to `a` alone and
 *    no surviving tet incident to `b` alone may share the same "opposite
 *    triangle" (the collapse would otherwise produce two tets with identical
 *    corners);
 *  - **tet-inversion rejection** (`smooth`'s "do no harm"): a surviving tet
 *    whose signed volume is non-zero must not change sign under the candidate
 *    placement (`detail::cell_volume_from_corners`); an already-degenerate
 *    tet imposes no constraint.
 *
 * **Boundary-touching collapses additionally run the existing 2D guards**
 * (`decimate`'s own ring/shared-face link condition and normal-flip check,
 * hoisted unchanged) over the mesh's boundary skin, so the outer surface
 * cannot tear, pinch or fold independently of the interior guards above.
 *
 * **Placement stays on the boundary surface.** A boundary vertex's quadric is
 * built from boundary-triangle planes only, so its minimizer naturally lies
 * near the surface; a pinned endpoint (feature/frozen/`mPreserveBoundary`)
 * forces its own position exactly as `decimate` does.
 *
 * **Data.** Same rules as `decimate`: float `point_data` is blended at the
 * clamped edge-projection parameter; integer `point_data` keeps the
 * survivor's row; each surviving tet keeps its own `cell_data` row;
 * `field_data` passes through verbatim.
 *
 * **Block structure is preserved 1:1**: the output has exactly
 * `NumCellBlocks()` tetra blocks in input order (possibly with zero rows),
 * exactly `decimate`'s own rule.
 *
 * **Determinism.** Setup is parallel with fixed FP order (matching
 * `decimate`'s); the greedy loop is serial, driven by a priority queue with
 * lazy version-stamped deletion. Every floating-point expression is
 * transcribed token for token into `_decimate_volume.py` (the arithmetic is
 * either reused verbatim from `decimate`/`_decimate.py`, already twinned, or
 * `detail::cell_volume_from_corners`'s corner-average fan, already twinned
 * elsewhere), so output is byte-identical across the three mesh backends,
 * thread counts, and the C++/numpy-fallback boundary.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/decimate.hpp"

namespace meshioplusplus {

/// Options for `decimate_volume`. Exactly one of the three stopping criteria
/// (`mTargetRatio`, `mTargetCells`, `mMaxError`) must be set (non-negative).
struct DecimateVolumeOptions {
    /// Fraction of the tets to KEEP, in `(0, 1]`. Negative means unset.
    double mTargetRatio = -1.0;

    /// Absolute number of tets to stop at (>= 0). Negative means unset. A
    /// collapse removes `T(ab).size()` tets, so the result lands within one
    /// collapse of the target.
    std::int64_t mTargetCells = -1;

    /// Collapse only while the cheapest candidate's score is at most this.
    /// Only strictly meaningful for `regime == 0` (boundary-touching, real
    /// quadric-error units) entries -- see the file doc comment. Negative
    /// means unset.
    double mMaxError = -1.0;

    /// Where the surviving vertex goes; reuses `decimate`'s enum unchanged.
    /// Overridden to the pinned endpoint's own position whenever one endpoint
    /// of the edge is pinned.
    DecimatePlacement mPlacement = DecimatePlacement::Optimal;

    /// Pin every boundary vertex outright (the once-used-face test on the
    /// mesh's own outer skin), reproducing `decimate`'s own default instead
    /// of letting boundary vertices participate.
    bool mPreserveBoundary = false;

    /// Pin boundary vertices whose incident boundary-triangle normals
    /// pairwise differ by more than `mFeatureAngleDeg`, so corners and
    /// creases of the outer surface survive.
    bool mPreserveFeatures = true;

    /// Angle **between two incident boundary-triangle normals**, in degrees,
    /// above which their shared vertex is a feature and pinned. Only read
    /// when `mPreserveFeatures` is set.
    double mFeatureAngleDeg = 30.0;

    /// Optional caller-supplied pin mask: either empty (no extra pins) or of
    /// length `NumPoints()`, where a non-zero entry pins that vertex. Unioned
    /// with the boundary/feature pins, never subtracted from them.
    std::vector<std::uint8_t> mFrozen;
};

/// The result of `decimate_volume`: the decimated mesh, the index maps, and
/// what the run actually did.
struct DecimateVolumeResult {
    /// The decimated mesh: all-tetra blocks, 1:1 with the input blocks.
    Mesh mMesh;

    /// Int64 shape `(num_points_in,)`, input point index -> output point
    /// index. A collapsed point maps to its **survivor's** output index, and
    /// is `-1` only when the surviving vertex itself ended up unreferenced
    /// and was pruned.
    NDArray mPointMap;

    /// Per **input** tet block, Int64 shape `(num_cells_in_block,)`, input
    /// cell -> its own output index, or `-1` when it did not survive.
    std::vector<NDArray> mCellMaps;

    /// Tets removed.
    std::int64_t mTetsRemoved = 0;

    /// `num_points_in - num_points_out` (collapsed plus pruned-unreferenced).
    std::int64_t mPointsRemoved = 0;

    /// Guard-rejection **events** (`decimate`'s own `mCollapsesRejected`
    /// convention: a re-scored, re-rejected edge counts again).
    std::int64_t mCollapsesRejected = 0;

    /// The largest `regime == 0` (quadric-error) score among the committed
    /// collapses (`0.0` when nothing collapsed, or when every collapse was
    /// purely interior).
    double mMaxErrorApplied = 0.0;
};

/**
 * @brief Decimates a tetrahedral mesh by greedy quadric-error-metric tet-edge
 * collapse.
 * @param rMesh The tet mesh to decimate (never modified).
 * @param rOptions Stopping criterion, placement, pinning policy, frozen mask.
 * @return The decimated all-tetra mesh, the point/cell index maps, and the
 *         collapse/rejection summary.
 * @throws std::invalid_argument when not exactly one stopping criterion is
 *         set, on a criterion out of range, on a mis-sized `mFrozen` mask, on
 *         a non-manifold boundary face, and on any block outside the tet-only
 *         scope: non-tetra 3D cells, ragged/polyhedron blocks, higher-order
 *         tets, and any non-3D block.
 */
MESHIOPLUSPLUS_API DecimateVolumeResult decimate_volume(const Mesh& rMesh,
                                                        const DecimateVolumeOptions& rOptions = {});

}  // namespace meshioplusplus
