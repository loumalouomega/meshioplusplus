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
 * @file decimate.hpp
 * @brief Surface decimation by quadric-error-metric edge collapse: reduce a
 * surface mesh's face count while preserving its shape, boundaries and
 * features.
 *
 * This is the resolution-*reducing* inverse of `refine`, completing the pair.
 * The algorithm is the classic Garland–Heckbert construction (quadric error
 * metrics, SIGGRAPH '97): every vertex accumulates a 4x4 symmetric quadric —
 * the sum of squared-distance-to-plane forms of its incident triangles,
 * area-weighted — and every candidate edge collapse is scored by the summed
 * endpoint quadric evaluated at the collapsed vertex's position. The cheapest
 * edge is collapsed greedily, its quadrics are summed onto the survivor, the
 * survivor's incident edges are re-scored, and the loop repeats until the
 * stopping criterion is met.
 *
 * **Surface only, triangles out.** The operating type is `triangle`; `quad`
 * and rectangular `polygon` blocks are triangulated first via
 * `convert_cells(Simplexify)`, so the output is all-triangle even when the
 * input was not. A mesh containing 3D volume cells throws by name pointing at
 * `extract_surface` — QEM decimation of a *volume* mesh (tetra-collapse
 * validity, no boundary-shape objective) is a different and much harder
 * problem, deliberately out of scope rather than silently skinned. Higher-order
 * cells (`triangle6`, ...; linearize first), ragged polygon/polyhedron blocks,
 * and `line`/`vertex` blocks mixed in with the surface likewise throw by name.
 *
 * **What is pinned** (mirroring `smooth`'s vocabulary and defaults). A pinned
 * vertex never moves and is never removed; an edge between two pinned vertices
 * never enters the queue, and a collapse toward a pinned vertex keeps that
 * vertex's own position regardless of `mPlacement`. Pinned are:
 *  - boundary vertices under `mPreserveBoundary` (the classic once-used-edge
 *    test: an edge used by exactly one triangle is boundary) — this is what
 *    keeps a planar patch's outline exactly intact;
 *  - feature vertices under `mPreserveFeatures`: vertices where two incident
 *    face **normals** differ by more than `mFeatureAngleDeg` (the
 *    `vtkFeatureEdges`/`smooth` convention), which is what keeps a cube's
 *    corners and creases sharp. Conservative v1: crease vertices are fully
 *    pinned rather than allowed to slide along the crease;
 *  - the caller's `mFrozen` mask.
 *
 * **Validity guards.** Each guard *rejects* the individual collapse (counted in
 * `mCollapsesRejected`) rather than aborting:
 *  - the **link condition**: the vertices adjacent to both endpoints must be
 *    exactly the opposite vertices of the faces on the edge, and an edge used
 *    by more than two faces is non-manifold and never collapsible — so the
 *    output cannot become non-manifold or change topology;
 *  - **normal-flip rejection**: a surviving incident face whose unit normal
 *    would turn by 90 degrees or more (`n_before . n_after <= 0`) rejects the
 *    collapse. Like `smooth`'s inversion guard this is strictly "do no harm":
 *    a face that is already degenerate (zero normal) imposes no constraint;
 *  - **degenerate survivor rejection**: a surviving face collapsing to zero
 *    area yields a zero `n_after`, hence `n_before . n_after = 0`, and is
 *    rejected by the same test.
 *
 * **Placement** (`mPlacement`): `Optimal` places the survivor at the minimizer
 * of the summed quadric (a 3x3 solve), falling back to the edge midpoint when
 * the system is ill-conditioned — `|det| <= 1e-12 * max|A_ij|^3`, a
 * scale-invariant Hadamard-style test. Note an exactly planar patch is
 * singular *by construction* (the quadric is flat along the plane), so on flat
 * geometry `Optimal` degenerates to midpoint placement; that is expected, not
 * a bug. `Midpoint` always uses the edge midpoint; `Endpoint` keeps the
 * endpoint with the lower quadric error (tie -> the lower vertex id).
 *
 * **Data.** Float-kind `point_data` at the surviving vertex is blended between
 * the two endpoints at the parameter `t` obtained by projecting the placed
 * point onto the edge, **clamped to [0, 1]**, and is carried in Float64 during
 * the run with a single cast back to the source dtype at the end. That `t` is
 * exact for `Midpoint`/`Endpoint` and an approximation for `Optimal` (the
 * optimal point need not lie on the edge). **Integer** `point_data` keeps the
 * survivor's own row — a blended material id is meaningless. Each surviving
 * face keeps its own `cell_data` row (a parent's row was already replicated to
 * its triangles by the simplexify step); `field_data` passes through verbatim.
 * `mesh.info`, `point_sets`/`cell_sets` and `gmsh_periodic` never reach the
 * C++ core (the Python shim remaps the sets through the returned maps).
 *
 * **Block structure is preserved 1:1**: the output has exactly
 * `NumCellBlocks()` triangle blocks in input order (possibly with zero rows),
 * which keeps the one-array-per-block `cell_data` invariant trivially correct
 * and makes `mCellMaps` input-block-indexed.
 *
 * **Determinism.** Setup is parallel with fixed FP order: the per-face plane
 * pass fills disjoint slots, and each vertex sums its quadric **in ascending
 * incident-face index** (FP addition is not associative; the order is the
 * pin). The greedy loop is **serial** — QEM is inherently sequential — driven
 * by a priority queue with the total order (error, lower endpoint id, higher
 * endpoint id, versions); stale entries are discarded lazily via per-vertex
 * version counters. The pinned expressions (plane, quadric, 3x3 solve, error,
 * blend parameter) are transcribed token for token into `_decimate.py`, so
 * output is byte-identical across the three mesh backends, thread counts, and
 * the C++/numpy-fallback boundary
 * (`tests/python/test_decimate.py::test_cpp_matches_python`).
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format — it is deliberately not in
 * the format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// How the surviving vertex of a collapsed edge is placed.
enum class DecimatePlacement {
    /// The minimizer of the summed quadric (3x3 solve); the edge midpoint when
    /// the system is ill-conditioned (`|det| <= 1e-12 * max|A_ij|^3`).
    Optimal,
    /// The edge midpoint.
    Midpoint,
    /// The endpoint with the lower quadric error (tie -> the lower vertex id).
    Endpoint,
};

/**
 * @brief Parses a placement name.
 * @param rName One of `"optimal"`, `"midpoint"`, `"endpoint"` (case-sensitive,
 *        as elsewhere in the operations layer).
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API DecimatePlacement decimate_placement_from_name(const std::string& rName);

/// Options for `decimate`. Exactly one of the three stopping criteria
/// (`mTargetRatio`, `mTargetFaces`, `mMaxError`) must be set (non-negative).
struct DecimateOptions {
    /// Fraction of the (triangulated) faces to KEEP, in `(0, 1]`. Negative
    /// means unset.
    double mTargetRatio = -1.0;

    /// Absolute number of faces to stop at (>= 0). Negative means unset. A
    /// collapse removes one or two faces, so the result lands within one
    /// collapse of the target (see `decimate`).
    std::int64_t mTargetFaces = -1;

    /// Collapse only while the cheapest candidate's quadric error is at most
    /// this (squared mesh units). Negative means unset.
    double mMaxError = -1.0;

    /// Where the surviving vertex goes. Overridden to the pinned endpoint's own
    /// position whenever one endpoint of the edge is pinned.
    DecimatePlacement mPlacement = DecimatePlacement::Optimal;

    /// Pin boundary vertices (once-used-edge test on the triangulated mesh):
    /// they never move and are never removed, and boundary–boundary edges never
    /// enter the queue, so the outline of an open patch is preserved exactly.
    bool mPreserveBoundary = true;

    /// Pin vertices whose incident face normals pairwise differ by more than
    /// `mFeatureAngleDeg`, so corners and creases survive.
    bool mPreserveFeatures = true;

    /// Angle **between two incident face normals**, in degrees, above which
    /// their shared vertex is a feature and pinned. `0` = coplanar, `90` = the
    /// edge of a box. Only read when `mPreserveFeatures` is set.
    double mFeatureAngleDeg = 30.0;

    /// Optional caller-supplied pin mask: either empty (no extra pins) or of
    /// length `NumPoints()`, where a non-zero entry pins that vertex. Unioned
    /// with the boundary/feature pins, never subtracted from them.
    std::vector<std::uint8_t> mFrozen;
};

/// The result of `decimate`: the decimated mesh, the index maps, and what the
/// run actually did.
struct DecimateResult {
    /// The decimated mesh: all-triangle blocks, 1:1 with the input blocks.
    Mesh mMesh;

    /// Int64 shape `(num_points_in,)`, input point index -> output point index.
    /// A collapsed point maps to its **survivor's** output index — which is
    /// what makes the map usable for remapping external per-point arrays — and
    /// is `-1` only when the surviving vertex itself ended up unreferenced and
    /// was pruned.
    NDArray mPointMap;

    /// Per **input** block, Int64 shape `(num_cells_in_block,)`, input cell ->
    /// the output index (within the corresponding output block) of its first
    /// surviving triangle, or `-1` when none survived. For a `quad`/`polygon`
    /// input cell this is composed through the simplexify step's child map.
    std::vector<NDArray> mCellMaps;

    /// Triangles removed, counted on the **triangulated** mesh.
    std::int64_t mFacesRemoved = 0;

    /// `num_points_in - num_points_out` (collapsed plus pruned-unreferenced).
    std::int64_t mPointsRemoved = 0;

    /// Guard-rejection **events** — an edge re-scored after a neighbouring
    /// collapse and rejected again counts again (`smooth`'s
    /// `mNumSkippedInversion` convention).
    std::int64_t mCollapsesRejected = 0;

    /// The largest quadric error among the committed collapses (`0.0` when
    /// nothing collapsed).
    double mMaxErrorApplied = 0.0;
};

/**
 * @brief Decimates a surface mesh by greedy quadric-error-metric edge collapse.
 * @param rMesh The surface mesh to decimate (never modified).
 * @param rOptions Stopping criterion, placement, pinning policy, frozen mask.
 * @return The decimated all-triangle mesh, the point/cell index maps, and the
 *         collapse/rejection summary.
 * @throws std::invalid_argument when not exactly one stopping criterion is
 *         set, on a criterion out of range, on a mis-sized `mFrozen` mask, and
 *         on any block outside the surface scope: 3D volume cells (use
 *         `extract_surface` first), higher-order cells (linearize first),
 *         ragged polygon/polyhedron blocks, and `line`/`vertex` blocks.
 */
MESHIOPLUSPLUS_API DecimateResult decimate(const Mesh& rMesh, const DecimateOptions& rOptions = {});

}  // namespace meshioplusplus
