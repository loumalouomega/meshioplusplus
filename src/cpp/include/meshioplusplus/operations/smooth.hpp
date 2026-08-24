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
 * @file smooth.hpp
 * @brief Geometric mesh smoothing: relax point coordinates toward their
 * edge-neighbour centroids without touching connectivity or data.
 *
 * This is the quality-*improving* member of the operations layer, and the
 * natural companion to `operations/quality.hpp`, which measures exactly what
 * this changes. Where `refine` increases resolution and `crop`/`clean` reduce
 * it, `smooth` leaves the mesh's size and structure alone and only moves nodes:
 * the output describes the same object with the same cells and the same field
 * values, but with better-shaped elements. It is the usual post-pass for a mesh
 * produced by `convert_cells`' simplexify or by a marching-cubes-style
 * extractor, both of which are topologically correct and geometrically ragged.
 *
 * Both operators are driven by the same edge-neighbour centroid displacement
 * `L(i) = mean(x[j] : j adjacent to i) - x[i]`:
 *
 *  - **Laplacian** — `x <- x + lambda * L(x)`, applied `mIterations` times.
 *    Strongly smoothing per pass and unconditionally **shrinking**: it is a
 *    low-pass filter with no pass band, so a closed surface run to convergence
 *    collapses toward a point. Use it when a few passes suffice and volume loss
 *    does not matter.
 *  - **Taubin** (the default) — each iteration is a `+lambda` pass followed by a
 *    `-mu` pass with `|mu| > lambda`. The second pass deliberately *un*-shrinks:
 *    the pair has a pass band below `k = 1/lambda + 1/mu` in which features are
 *    preserved, and a stop band above it in which noise is attenuated. With the
 *    defaults (0.33 / -0.34) hundreds of iterations leave the enclosed volume
 *    essentially unchanged, which is why it is the default.
 *
 * **The neighbour graph is EDGE adjacency, not the element clique**
 * (`detail/node_adjacency.hpp`, `NodeAdjacencyKind::Edge`). Under the clique
 * graph that `reorder` uses, a hexahedron would pull each corner toward its
 * three face diagonals and its body diagonal as hard as toward its three real
 * edge neighbours, bevelling cubes toward spheres; the edge graph leaves a
 * structured hex block a fixed point.
 *
 * **What is pinned.** A node never moves if it is
 *  - on the boundary and `mFixBoundary` is set,
 *  - a feature node under `mPreserveFeatures`,
 *  - set in the caller's `mFrozen` mask, or
 *  - referenced only by blocks whose edge topology is unknown (the higher-order
 *    family, the VTK-Lagrange types, `custom`) — an unknown neighbourhood means
 *    an undefined target, so the node holds still rather than being guessed at.
 *
 * Boundary is the classic once-used-facet test: a face used by exactly one 3D
 * cell, or — in a pure surface mesh — an edge used by exactly one 2D cell.
 * Feature nodes are boundary nodes where two incident boundary facet **normals**
 * differ by more than `mFeatureAngleDeg` (0 = coplanar, 90 = the edge of a box;
 * the `vtkFeatureEdges` convention), which is what keeps a cube's edges sharp
 * instead of rounding them off.
 *
 * **Inversion guard.** With `mGuardInversion`, a node's candidate position is
 * committed only if no incident cell's signed measure changes sign; otherwise
 * that node holds still *for that pass* and the event is counted. The guard
 * covers the cells it can measure — `tetra`/`hexahedron`/`wedge`/`pyramid` by
 * the outward face fan, 2D cells by the signed shoelace area in a 2D mesh and
 * by a normal-flip test in a 3D one. Note that last mode has no counterpart in
 * `compute_quality`, which correctly declines to call a triangle floating in 3D
 * "inverted" for want of a global orientation convention; smoothing needs the
 * weaker, purely relative question "did this facet just fold over?", so it asks
 * that one instead. A cell with no signed measure at all (`line`, polyhedron,
 * `custom`) is **skipped, not pinned**: a missing safety check is not a missing
 * target, and pinning there would freeze an entire beam-element curve for no
 * geometric reason.
 *
 * **Geometry only.** Connectivity, `cell_data`, `field_data` and `point_data`
 * *values* are carried through unchanged; only the point coordinates move, and
 * the points array keeps its input dtype. Iteration runs in `double` regardless
 * of that dtype, with a single cast on write-back, so a Float32 mesh does not
 * accumulate one rounding per pass.
 *
 * **Determinism.** Every pass is **Jacobi**: all new positions are computed from
 * the previous pass's positions and committed together, so no node ever observes
 * a half-updated neighbour and the result cannot depend on traversal order. The
 * update loop contains no hashing, no sorting and no floating-point reduction;
 * neighbour sums run in ascending neighbour id (the adjacency rows are sorted),
 * the boundary set comes from the **serial** dedup pass over a
 * `parallel_for`-filled disjoint-slot buffer (`operations/surface.cpp`'s
 * phase-split idiom), and the summary counters are folded serially out of
 * per-node slots. Output is byte-identical across mesh backends and thread
 * counts.
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

namespace meshioplusplus {

/// Which smoothing operator to apply.
///
/// Fixed `: std::uint8_t` underlying type (the `RemeshMetric`/`CellType`
/// precedent) so a width change is a mechanical ABI check, not a silent
/// layout shift; every dispatch site in `smooth.cpp` is an exhaustive
/// `switch (SmoothMethod)` with **no `default:`**, so `-Wswitch` catches the
/// next enumerator at compile time rather than a reader having to -- the
/// hazard `remesh.cpp`'s own metric branches were converted away from.
enum class SmoothMethod : std::uint8_t {
    Laplacian,  ///< `x <- x + lambda * L(x)`. Simple and strong; shrinks volume.
    Taubin,     ///< A `+lambda` then `-mu` pass pair per iteration. Shrink-free.
    /// Optimal Delaunay Triangulation: each free interior vertex moves to the
    /// volume-weighted circumcenter combination that minimises the local ODT
    /// energy (Alliez et al. 2005; Chen & Xu) -- **tet-only**, throwing by
    /// name on any other cell type, unlike `Laplacian`/`Taubin` which work on
    /// any mesh. `mMu` is silently ignored, matching `Laplacian`.
    Odt,
};

/// `sizeof(SmoothMethod)` mirror of the primary guard in `test_abi_layout.cpp`
/// -- the `RemeshMetric` precedent, so a width change is caught here too.
static_assert(sizeof(SmoothMethod) == 1,
              "meshio++ ABI: SmoothMethod's underlying type changed width. "
              "Appending enumerators is safe and expected; widening "
              "`: std::uint8_t` is a Tier A break (doc/abi.md).");

/**
 * @brief Parses a smoothing method name.
 * @param rName One of `"laplacian"`, `"taubin"`, `"odt"` (case-sensitive, as
 *        elsewhere in the operations layer).
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API SmoothMethod smooth_method_from_name(const std::string& rName);

/// Options for `smooth`.
struct SmoothOptions {
    /// The smoothing operator.
    SmoothMethod mMethod = SmoothMethod::Taubin;

    /// How many iterations to run. For `Taubin` one iteration is **two** passes
    /// (`+mLambda` then `mMu`); for `Laplacian`/`Odt` one iteration is one pass.
    /// Zero or less returns an unchanged clone.
    int mIterations = 10;

    /// Relaxation factor of the smoothing pass, which must lie in `(0, 1)`.
    /// **Negative means "this method's own default"**: `0.5` for `Laplacian`,
    /// `0.33` for `Taubin`, `0.9` for `Odt` (a near-full step toward the local
    /// ODT energy minimiser -- damped rather than exactly `1.0` so an
    /// overshoot past a concave boundary is a rejected inversion-guard move
    /// rather than a hard failure). The sensible default genuinely differs
    /// per method, and a factor of zero or less is never a meaningful
    /// request, so a negative sentinel is total; separate fields would
    /// instead let a caller set the ones the chosen method ignores.
    double mLambda = -1.0;

    /// `Taubin` only: the un-shrinking factor, which must satisfy
    /// `mMu < -mLambda < 0`. Silently ignored by `Laplacian` **and** `Odt` --
    /// deliberately not a named-error option like `RemeshOptions::mMaxAnisotropy`,
    /// since making only the newest method strict about an ignored field
    /// while its sibling stays silent would be a worse inconsistency than
    /// either rule applied uniformly, and changing `Laplacian`'s established
    /// behaviour here would be breaking.
    double mMu = -0.34;

    /// Pin every node lying on a boundary facet. Almost always what you want:
    /// with this off, a closed surface contracts and an open domain's walls
    /// creep inward.
    bool mFixBoundary = true;

    /// Additionally pin boundary nodes where incident boundary facets meet at
    /// more than `mFeatureAngleDeg`, so geometric corners and creases survive.
    /// Has no effect when `mFixBoundary` is false — every boundary node is then
    /// already free by explicit request.
    bool mPreserveFeatures = true;

    /// Angle **between two boundary facet normals**, in degrees, above which
    /// their shared nodes are treated as a feature and pinned. `0` = coplanar,
    /// `90` = the edge of a box. Only read when `mPreserveFeatures` is set.
    double mFeatureAngleDeg = 30.0;

    /// Reject a node's move when it would flip the signed measure of any
    /// incident cell, and count the rejection. Costs a second CSR (node ->
    /// incident cell) plus a corner table, so both are built only when set.
    bool mGuardInversion = true;

    /// Optional caller-supplied pin mask: either empty (no extra pins) or of
    /// length `NumPoints()`, where a non-zero entry pins that node. Unioned with
    /// the boundary/feature/unknown-topology pins, never subtracted from them.
    std::vector<std::uint8_t> mFrozen;

    /// Net displacement below `mMoveTolerance * bbox_diagonal` does not count
    /// toward `SmoothResult::mNumNodesMoved`, so a node already sitting at its
    /// own centroid is not reported as moved by its own round-off. Relative,
    /// hence scale-invariant — an absolute threshold would report everything as
    /// moved on a mesh in metres and nothing on the same mesh in kilometres.
    double mMoveTolerance = 1e-12;
};

/// The result of `smooth`: the relaxed mesh plus what the run actually did.
struct SmoothResult {
    /// The smoothed mesh: same cells, same data, moved points.
    Mesh mMesh;

    /// Number of points whose **net** displacement over the whole run exceeded
    /// `mMoveTolerance * bbox_diagonal`. Pinned nodes have a net displacement of
    /// exactly zero and never count.
    std::int64_t mNumNodesMoved = 0;

    /// The largest net displacement `max_i |x_out[i] - x_in[i]|`, in mesh units.
    /// Measured against the input rather than the previous pass, so it answers
    /// "how far did this operation move the mesh", not "did the last pass
    /// converge".
    double mMaxDisplacement = 0.0;

    /// Number of `(node, pass)` **events** in which the inversion guard rejected
    /// a move — not a node count, since one node may be rejected in several
    /// passes. Always zero when `mGuardInversion` is false.
    std::int64_t mNumSkippedInversion = 0;
};

/**
 * @brief Smooths a mesh's point coordinates, leaving topology and data intact.
 * @param rMesh The mesh to smooth (never modified).
 * @param rOptions Method, iteration count, pinning policy and guard settings.
 * @return The smoothed mesh and the run's displacement/rejection summary.
 * @throws std::invalid_argument on a mis-sized `mFrozen` mask, an `mLambda`
 *         outside `(0, 1)`, or a `Taubin` `mMu` that does not satisfy
 *         `mMu < -mLambda < 0` (which would make the pass pair amplify rather
 *         than filter).
 */
MESHIOPLUSPLUS_API SmoothResult smooth(const Mesh& rMesh, const SmoothOptions& rOptions = {});

}  // namespace meshioplusplus
