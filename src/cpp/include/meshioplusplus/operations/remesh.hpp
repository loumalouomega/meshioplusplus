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
 * **Boundaries are not specially protected in v1.** An open surface's boundary
 * vertices are ordinary items with no extra pinning or dual-edge insertion (the
 * `refine`/`decimate`-style "boundary preservation" the roadmap's own
 * "Boundaries and output manifoldness" bullet still records as open); this
 * operation runs on an open mesh without error, matching ACVD's own `-b 0`
 * default, but the outline near an open edge is not guaranteed preserved.
 *
 * **Topology is not preserved.** Two surface sheets closer together than a
 * cluster can merge, and the genus can change. This is inherent to dualising a
 * discrete Voronoi partition, not an implementation limit -- a target vertex
 * count far below the feature scale will visibly simplify topology.
 *
 * **Output manifoldness is best-effort.** The dual of a discrete Voronoi
 * partition need not be 2-manifold. Repair removes the common cause
 * (disconnected clusters), and `mNumIsolatedClusters` reports what could not be
 * fixed, but a pathological input can still produce non-manifold output; check
 * `mNumIsolatedClusters` rather than assuming.
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
enum class RemeshMetric {
    /// Area-weighted centroidal distance (Valette & Chassery, CGF 2004). Fast,
    /// robust, and rounds sharp features.
    Isotropic,
    /// Garland-Heckbert quadric error ("ACVDQ" in the literature). Preserves
    /// sharp edges and corners at the cost of one extra 10-entry accumulator
    /// per cluster and a 3x3 solve per candidate move.
    Quadric,
};

/**
 * @brief Parses a metric name.
 * @param rName One of `"isotropic"`, `"quadric"` (case-sensitive, as
 *        elsewhere in the operations layer).
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API RemeshMetric remesh_metric_from_name(const std::string& rName);

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
 *         subdivided input's vertex count, on a non-positive limit, and on any
 *         block outside the surface scope: 3D volume cells (use
 *         `extract_surface` first), higher-order cells (linearize first),
 *         ragged polygon/polyhedron blocks, and `line`/`vertex` blocks.
 */
MESHIOPLUSPLUS_API RemeshResult remesh(const Mesh& rMesh, const RemeshOptions& rOptions);

}  // namespace meshioplusplus
