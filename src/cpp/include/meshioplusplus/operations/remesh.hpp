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
//  Attribution:     The clustering algorithm is derived from pyacvd
//                   (https://github.com/pyvista/pyacvd), MIT License,
//                   Copyright (c) 2017-2024 The PyVista Developers, which is
//                   itself an independent implementation of the research of
//                   S. Valette and J.-M. Chassery. See doc/remesh.md.
//
#pragma once

/**
 * @file remesh.hpp
 * @brief Isotropic and feature-preserving surface remeshing by approximated
 * centroidal Voronoi diagram (ACVD) clustering: replace a surface mesh's
 * triangulation with a new, near-uniformly-sized, well-shaped one at a
 * caller-chosen vertex count.
 *
 * This is the one resolution operation that does not work on the input's own
 * triangulation. `refine` subdivides the input's cells, `decimate` collapses
 * its edges, `subdivide`/`agglomerate` restructure it into polyhedra, `smooth`
 * moves its points -- all of them inherit the input's element shapes. QEM
 * decimation in particular can only *remove* elements, so a badly-shaped input
 * stays badly-shaped at every target count. `remesh` instead partitions the
 * surface into `mNumClusters` compact regions of near-equal area (or, under
 * `RemeshMetric::Quadric`, near-equal quadric error) and builds the **dual**
 * of that partition, so output element quality is a property of the
 * clustering rather than of the input.
 *
 * **Attribution.** The clustering engine (seeding, energy minimisation, the
 * cluster-repair passes, the dual triangulation) is derived from
 * `pyvista/pyacvd`'s `src/clustering.cpp` and `src/pyacvd/clustering.py`
 * (nanobind, MIT, (c) 2017-2024 The PyVista Developers) -- an independent
 * implementation of Valette & Chassery's published *research*, not of
 * ACVD's own CeCILL-B source, which this project never reads code from and
 * never vendors (see "Feature-preserving metric" below and doc/remesh.md for
 * the licence reasoning). `RemeshMetric::Quadric` has no pyacvd counterpart --
 * pyacvd implements only the isotropic variant -- and is a fresh synthesis
 * over meshio++'s own pre-existing Garland-Heckbert quadric machinery
 * (`detail/decimate_common.hpp`, shared with `decimate`/`decimate_volume`).
 *
 * **The isotropic algorithm** (Valette & Chassery, CGF 2004). Each input
 * vertex is an "item" carrying a weight (one third of its incident triangle
 * area, so the weights sum to the surface area) and a weighted position. A
 * cluster holds only two accumulators -- `sgamma = sum(w * x)` and
 * `srho = sum(w)` -- which is what makes testing a candidate move O(1) and the
 * whole method fast:
 *
 *   - **Seeding** grows `mNumClusters` regions breadth-first from the lowest
 *     unassigned vertex, each absorbing neighbours until it reaches its share
 *     of the remaining area. Deliberately *deterministic*: unlike pyacvd this
 *     uses no RNG at all, so a given mesh always yields the same clustering.
 *   - **Energy minimisation** sweeps the edges whose endpoints lie in
 *     different clusters and, for each, compares the current badness against
 *     moving either endpoint across. Under `Isotropic` the objective
 *     `sum(|sgamma_c|^2 / srho_c)` is *maximised* (equivalently, a "badness"
 *     `-sum(...)` is minimised, the uniform framing both metrics share); this
 *     is equivalent to minimising the centroidal Voronoi energy but avoids
 *     recomputing centroids. A cluster flagged unmodified in the previous
 *     sweep is skipped, so late sweeps are cheap.
 *   - **Repair** is not optional: seeding can leave vertices unassigned (a
 *     pocket the fronts never reached) and minimisation can split a cluster
 *     into disconnected components. Unassigned vertices are grown into from
 *     their neighbours; a split cluster keeps its largest component and the
 *     rest are re-grown, alternating with further minimisation sweeps up to
 *     `mMaxRepairPasses`.
 *   - **The dual** places one output vertex at each cluster's representative
 *     point (area-weighted centroid under `Isotropic`, see below under
 *     `Quadric`) and emits one triangle for every input triangle whose three
 *     vertices lie in three *distinct* clusters, deduplicated. Winding is
 *     fixed by comparing each output triangle's normal against the mean input
 *     normal of the clusters it joins.
 *
 * **Feature-preserving metric** (`RemeshMetric::Quadric`, "ACVDQ" in the
 * literature's naming, roadmap-named "the same clustering driven by a quadric
 * error metric"). This project has no access to and does not read ACVD's own
 * `vtkQEMetricForClustering.h` (CeCILL-B, VTK-bound) -- the design below is
 * built entirely from the general "metric-dependent discrete Voronoi diagram"
 * concept (Valette, Chassery & Prost, IEEE TVCG 2008) plus meshio++'s own,
 * already-shipped Garland-Heckbert quadric primitives
 * (`detail::decim_face_planes`/`decim_accumulate_quadrics`/
 * `decim_quadric_error`/`decim_quadric_optimal_point`, the last hoisted out of
 * `decimate`'s edge-collapse placement solve specifically so this operation
 * can share it rather than re-derive it -- "a new metric, not a second QEM
 * implementation"). Each cluster carries a 10-entry quadric *in addition to*
 * `sgamma`/`srho` (kept in both metrics -- cheap, and `Quadric` needs the
 * isotropic centroid as its solve's ill-conditioning fallback), accumulated
 * with the same O(1)-per-move update shape. The per-move badness becomes the
 * summed quadric error of the two affected clusters, each evaluated at its own
 * `decim_quadric_optimal_point`; the dual's vertex placement likewise moves
 * from the area-weighted centroid to that same optimal point. A flat region's
 * quadric is near-singular in-plane and degenerates to the centroid (the same
 * ill-conditioning fallback `decimate` already relies on); a crease or corner's
 * quadric is not, and pins the representative point onto the feature instead
 * of smearing it toward a centroid -- which is the entire quality difference
 * between the two metrics.
 *
 * **Anisotropic metric** (`RemeshMetric::Anisotropic`, `mMaxAnisotropy`).
 * Also built from the Valette/Chassery/Prost TVCG 2008 concept above -- the
 * same paper, not a second attribution. Each vertex's curvature fit
 * (`remesh_vertex_curvature`, shared with `mGradation`) already computes the
 * principal curvatures `kappa1`/`kappa2` and their tangent-plane directions;
 * this metric is the one consumer that keeps the directions instead of
 * collapsing them to `max(|kappa1|, |kappa2|)`. A per-vertex SPD tensor `M`
 * is built with in-plane eigenvalues proportional to curvature magnitude
 * (short target edges across a sharp bend, long ones along a flat one),
 * ratio clamped to `mMaxAnisotropy`, normal eigenvalue pinned to the sharper
 * in-plane one (keeps clusters on the surface), then rescaled to `det(M) = 1`
 * -- the metric only reshapes clusters, it never rescales the objective, so
 * it composes cleanly with `mGradation` (density) and stays comparable
 * vertex-to-vertex regardless of the surface's absolute curvature magnitude.
 * `M = I` exactly wherever the surface is locally flat or `mMaxAnisotropy ==
 * 1.0`. Algebraically, `E(v) = sum(w * (x-v)^T M (x-v))` expands to the same
 * `x^T A x + 2 b.x + c` form the quadric metric already accumulates and
 * minimises via `decim_quadric_optimal_point` -- so this packs into the
 * *same* 10-entry accumulator, needs no new per-cluster storage, and the
 * per-move cost is identical to `Quadric`'s.
 *
 * **Subdivision matters more than iteration count.** The clustering is a
 * *discrete* approximation of a centroidal Voronoi diagram whose resolution is
 * the number of items per cluster, so a target near the input vertex count
 * gives poor results until the input is refined. `mSubdivide` applies that many
 * uniform `refine` passes first (each multiplying the triangle count by four);
 * `remesh_suggest_subdivide` reports how many are needed to reach
 * `mSubsampleRatio` items per cluster.
 *
 * **Surface only, triangles out.** The operating type is `triangle`; `quad` and
 * rectangular `polygon` blocks are triangulated first via
 * `convert_cells(Simplexify)`. A mesh containing 3D volume cells throws by name
 * pointing at `extract_surface` (the volumetric CVD/ODT problem is genuinely
 * different, not a silent skin); higher-order cells throw pointing at
 * `convert_cells(Linearize)`; ragged polygon/polyhedron and `line`/`vertex`
 * blocks throw by name. This mirrors `decimate`'s scope exactly.
 *
 * **What does NOT survive.** The output is a new mesh with new vertices and new
 * connectivity, so -- unlike every other operation in this layer -- there is no
 * meaningful point or cell map, and `point_data`, `cell_data` and named regions
 * are **dropped** (`warn_regions_dropped` is called, never silent). Transfer a
 * field onto the result with `interpolate` or `conservative_interpolate`, which
 * is exactly the composition those operations exist for. `field_data` passes
 * through verbatim.
 *
 * **Curvature gradation** (`mGradation`, the roadmap's `area * kappa^gamma`
 * density weight). Each item's clustering weight is normally its plain
 * triangle-fan area; when `mGradation` (the exponent `gamma`) is non-zero, the
 * weight becomes `area * max(kappa, kappa_floor)^gamma`, where `kappa` is a
 * per-vertex curvature magnitude from a local osculating-paraboloid fit over
 * the vertex's 1-ring (`remesh_vertex_curvature`, built on the same
 * `detail::NodeAdjacency` graph the edge sweep already uses -- no separate
 * neighbourhood machinery). `kappa_floor` is a small internal fraction of the
 * mesh's own maximum curvature, guarding against an exact-zero weight
 * degenerating a flat vertex's cluster membership. This is the **single**
 * weight threaded everywhere `sgamma`/`srho`/the badness/the final placement
 * already use, under both metrics unchanged -- gradation is orthogonal to the
 * `Isotropic`/`Quadric` choice, not a third metric. `mGradation == 0.0` (the
 * default) means `kappa^0 == 1` identically, so curvature is never even
 * computed and every existing weight is reproduced byte-for-byte.
 *
 * **Boundaries** (`mPreserveBoundary`, default **true** -- a no-op, and so
 * free, on the closed meshes most callers and every pre-existing test use).
 * An open surface's boundary vertices are detected once
 * (`remesh_boundary_info`, a local edge-use-count pass: a triangle edge used
 * by exactly one face is a boundary edge) and seeded **before** the general
 * interior BFS, walking each boundary chain/loop first so a cluster's
 * boundary segment is anchored contiguously rather than being absorbed
 * piecemeal by whichever interior cluster reaches it first -- the "pinned"
 * half of the roadmap's "boundary items pinned and extra dual elements
 * inserted along the boundary." The second half is dual-side: a genuine
 * boundary edge whose two endpoints land in different clusters emits a `line`
 * cell between those two clusters' representative points (deduplicated by
 * cluster-id pair, the same idiom the triangle dual already uses), so the
 * output gains a second, optional `line` cell block carrying the boundary
 * polyline -- without it, a boundary-adjacent triangle simply has no "third
 * neighbour" across the missing side and is silently dropped from the
 * triangle dual, leaving no coherent output boundary at all.
 *
 * This is a clean-room design achieving the roadmap's stated boundary goals,
 * built from this project's own conventions (the pinning idiom `decimate`/
 * `smooth` already use, the dedup idiom the triangle dual already uses) --
 * not a reproduction of ACVD's own boundary-fixing algorithm, whose source
 * this project does not read (see the licence reasoning above).
 *
 * **Topology is not preserved.** Two surface sheets closer together than a
 * cluster can merge, and the genus can change. This is inherent to dualising a
 * discrete Voronoi partition, not an implementation limit -- a target vertex
 * count far below the feature scale will visibly simplify topology.
 *
 * **Output manifoldness is best-effort.** The dual of a discrete Voronoi
 * partition need not be 2-manifold, for two distinct reasons, reported
 * separately. Disconnected clusters (`mNumIsolatedClusters`) are the
 * `remesh_split_disconnected` repair loop's own concern, unchanged from
 * v10.10.0. Non-manifold **output vertices** (`mNumNonManifoldVertices`) are
 * a second, independent cause -- a vertex whose incident dual-triangle fan
 * does not form a single loop (interior) or open chain (boundary), i.e. a
 * "bowtie" -- detected by `remesh_detect_nonmanifold_vertices` and folded
 * into the **same** repair loop already driving disconnected-cluster repair
 * (its implicated clusters are unassigned and regrown/reminimised alongside
 * whatever `remesh_split_disconnected` already found, up to
 * `mMaxRepairPasses`), rather than a second, separate loop. Either counter
 * non-zero means a pathological input still produced non-manifold output;
 * check both rather than assuming.
 *
 * **Determinism, and no numpy twin.** Seeding is RNG-free, the energy sweep
 * visits edges in a fixed order and is **serial** (the objective is inherently
 * sequential -- every accepted move changes the accumulators the next test
 * reads), and the setup passes fill disjoint slots in parallel with a fixed FP
 * order. Output is therefore byte-identical across the three mesh backends and
 * thread counts, but there is deliberately **no numpy twin**: a single near-tie
 * decided differently by an independent implementation diverges into a
 * different clustering, not a last-ulp difference, the same reasoning that
 * makes `subdivide`/`agglomerate`/`decimate_volume`/`conservative_interpolate`
 * C++-core-only. `_remesh.py` raises `NotImplementedError` by name when `_core`
 * is unavailable, for any input.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format -- it is deliberately not in
 * the format registry.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The clustering objective driving `remesh`.
enum class RemeshMetric : std::uint8_t {
    /// Area-weighted centroidal distance (Valette & Chassery, CGF 2004). Fast,
    /// robust, and rounds sharp features.
    Isotropic,
    /// Garland-Heckbert quadric error ("ACVDQ" in the literature). Preserves
    /// sharp edges and corners at the cost of one extra 10-entry accumulator
    /// per cluster and a 3x3 solve per candidate move.
    Quadric,
    /// Clusters shaped by a local curvature tensor rather than isotropic
    /// distance (Valette, Chassery & Prost, IEEE TVCG 2008 -- the same paper
    /// `Quadric` cites): elongated features get elongated elements. Packs
    /// into the same 10-entry accumulator `Quadric` uses, so it costs the
    /// same per-move solve; see `RemeshOptions::mMaxAnisotropy`.
    Anisotropic,
};

/// A fixed underlying type, checked below, so appending a future enumerator
/// stays a Tier C addition (doc/abi.md) rather than accidentally widening
/// the type an already-compiled consumer stored.
static_assert(sizeof(RemeshMetric) == 1,
              "meshio++ ABI: RemeshMetric's underlying type changed width. "
              "Appending enumerators is safe and expected; widening "
              "`: std::uint8_t` is a Tier A break (doc/abi.md).");

/**
 * @brief Parses a metric name.
 * @param rName One of `"isotropic"`, `"quadric"`, `"anisotropic"`
 *        (case-sensitive, as elsewhere in the operations layer).
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API RemeshMetric remesh_metric_from_name(const std::string& rName);

/// Default for `RemeshOptions::mMaxAnisotropy` -- meaningful only under
/// `RemeshMetric::Anisotropic`, but a single named constant so the header's
/// default and the "you set this under the wrong metric" guard in
/// `remesh.cpp` cannot drift apart. Chosen by measurement, not a round
/// number: on a near-umbilic input (the icosahedron test fixture, whose
/// TRUE surface is a sphere -- equal principal curvatures everywhere, so
/// the metric has no genuine anisotropy signal to read, only discretisation
/// noise), a swept `max_anisotropy` of 1/2/4/10 produced 0/0/1/18 clusters
/// short of a 120-cluster target -- energy minimisation can fully absorb a
/// seeded cluster into its neighbours (not a disconnection `mMaxRepairPasses`
/// helps with; `RemeshResult::mNumClusters` is already documented as "never
/// higher", and more repair passes measurably did not change the count
/// here). `4.0` keeps that near-worst-case shortfall to a single cluster
/// while still giving genuinely anisotropic surfaces (a cylinder, a fillet)
/// real elongation headroom; a caller expecting to hit an exact count on
/// nearly-umbilic geometry should still prefer something closer to `2.0`.
inline constexpr double kRemeshDefaultMaxAnisotropy = 4.0;

/// Options for `remesh`.
struct RemeshOptions {
    /// Number of clusters, i.e. the number of vertices the output aims for.
    /// Must be at least 4 (fewer cannot bound a surface) and no more than the
    /// number of vertices of the subdivided input.
    std::int64_t mNumClusters = 0;

    /// Uniform `refine` passes applied to the input before clustering, each
    /// multiplying the triangle count by four. `-1` (the default) picks the
    /// smallest count reaching `mSubsampleRatio` items per cluster, capped at
    /// `mMaxSubdivide`; `0` disables subdivision entirely.
    int mSubdivide = -1;

    /// Items per cluster targeted by automatic subdivision. Higher is a finer
    /// discrete approximation and a slower run.
    double mSubsampleRatio = 10.0;

    /// Ceiling on automatically chosen subdivision, so a tiny input and a large
    /// cluster count cannot silently allocate an enormous intermediate mesh.
    /// Ignored when `mSubdivide` is set explicitly.
    int mMaxSubdivide = 4;

    /// Maximum energy-minimisation sweeps per pass. The loop also stops early
    /// as soon as a sweep changes nothing.
    int mMaxIterations = 100;

    /// How many times to alternate "split disconnected clusters, minimise
    /// again" after the first minimisation. `0` skips repair of disconnected
    /// clusters entirely (the fast, lower-quality mode).
    int mMaxRepairPasses = 10;

    /// The clustering objective. See `RemeshMetric`.
    RemeshMetric mMetric = RemeshMetric::Isotropic;

    /// Curvature-gradation exponent `gamma` in the item weight
    /// `area * kappa^gamma`. `0.0` (the default) disables gradation entirely
    /// (curvature is not even computed) and reproduces plain area weighting.
    /// Positive values concentrate clusters where the surface bends more
    /// sharply.
    double mGradation = 0.0;

    /// Detect the input's open boundary (if any) and seed it before the
    /// interior, then emit a `line` dual cell along every boundary edge whose
    /// endpoints land in different clusters. A no-op, and so free, whenever
    /// the input has no boundary edges (e.g. every mesh in this repo's own
    /// closed-surface tests).
    bool mPreserveBoundary = true;

    /// Under `RemeshMetric::Anisotropic`, the maximum ratio between the two
    /// in-plane target edge lengths (long axis / short axis) a per-vertex
    /// curvature tensor is allowed to request; `1.0` recovers the isotropic
    /// shape exactly. Must be `>= 1.0`. An error to set away from
    /// `kRemeshDefaultMaxAnisotropy` under any other metric -- there is
    /// nothing for it to shape, so a silently-ignored value would be a
    /// caller mistake with no diagnostic.
    double mMaxAnisotropy = kRemeshDefaultMaxAnisotropy;
};

/// The result of `remesh`.
struct RemeshResult {
    /// The remeshed surface: exactly one `triangle` block, new points and new
    /// connectivity. `point_data`, `cell_data` and regions are dropped;
    /// `field_data` is carried.
    Mesh mMesh;

    /// Clusters actually produced, i.e. `mMesh.NumPoints()`. It can be lower
    /// than `RemeshOptions::mNumClusters` when the surface cannot be split that
    /// finely (small or disconnected input); it is never higher.
    std::int64_t mNumClusters = 0;

    /// Energy-minimisation sweeps performed, summed over every pass. Reaching
    /// `mMaxIterations` on the first pass means the run was truncated rather
    /// than converged.
    std::int64_t mNumIterations = 0;

    /// Subdivision passes actually applied (the resolved value of
    /// `RemeshOptions::mSubdivide`).
    int mSubdivideApplied = 0;

    /// Clusters still disconnected after `mMaxRepairPasses`. Non-zero means the
    /// output may be non-manifold near those clusters.
    std::int64_t mNumIsolatedClusters = 0;

    /// Output vertices whose incident dual-triangle fan is still not a single
    /// loop/chain (a "bowtie") after `mMaxRepairPasses`. Non-zero means the
    /// output may be non-manifold at those vertices -- a distinct cause from
    /// `mNumIsolatedClusters`, reported separately.
    std::int64_t mNumNonManifoldVertices = 0;
};

/**
 * @brief Suggests a subdivision count for a target cluster count.
 * @param NumPoints Vertices of the input surface.
 * @param NumClusters Target cluster count.
 * @param SubsampleRatio Items per cluster to aim for.
 * @param MaxSubdivide Ceiling on the returned value.
 * @return The smallest count in `[0, MaxSubdivide]` reaching the ratio.
 */
MESHIOPLUSPLUS_API int remesh_suggest_subdivide(std::int64_t NumPoints, std::int64_t NumClusters,
                                                double SubsampleRatio, int MaxSubdivide);

/**
 * @brief Remeshes a surface uniformly by ACVD clustering.
 * @param rMesh The surface mesh to remesh (never modified).
 * @param rOptions Cluster count, subdivision, iteration limits and metric.
 * @return The remeshed all-triangle mesh and the run summary.
 * @throws std::invalid_argument on a cluster count below 4 or above the
 *         subdivided input's vertex count, on a non-positive limit, on
 *         `mMaxAnisotropy < 1.0`, on `mMaxAnisotropy` set away from its
 *         default under a non-`Anisotropic` metric, and on any block outside
 *         the surface scope: 3D volume cells (use `extract_surface` first),
 *         higher-order cells (linearize first), ragged polygon/polyhedron
 *         blocks, and `line`/`vertex` blocks.
 */
MESHIOPLUSPLUS_API RemeshResult remesh(const Mesh& rMesh, const RemeshOptions& rOptions);

}  // namespace meshioplusplus
