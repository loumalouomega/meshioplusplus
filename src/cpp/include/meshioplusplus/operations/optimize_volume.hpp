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
//  Attribution:     The topological-flip acceptance rule (apply a local
//                   transformation only when it strictly improves the worst
//                   incident element quality) is the published "local mesh
//                   improvement by swapping" criterion of Freitag &
//                   Ollivier-Gooch, "Tetrahedral mesh improvement using
//                   swapping and smoothing" (Int. J. Numer. Methods Eng.,
//                   1997) -- implemented from the description only. No external
//                   mesh library, and in particular no in-sphere / orientation
//                   predicate kernel, is read or vendored. See
//                   doc/optimize_volume.md.
//
#pragma once

/**
 * @file optimize_volume.hpp
 * @brief Improve a tetrahedral mesh's element quality by ODT *remeshing* --
 * relocating vertices AND changing connectivity -- without any Delaunay /
 * in-sphere predicate.
 *
 * This is the genuine "ODT remeshing" the roadmap's volumetric bullet named.
 * Its two siblings each do half the job and are deliberately kept separate:
 *
 *  - `SmoothMethod::Odt` (`operations/smooth.hpp`) is ODT *smoothing* -- it
 *    moves each free interior tet vertex to the volume-weighted circumcenter
 *    average of its incident tets, on the mesh's **fixed** connectivity. It
 *    cannot fix a badly *connected* tetrahedralization (a sliver that no
 *    vertex motion removes), only a badly *placed* one.
 *  - `remesh_volume` (`operations/remesh_volume.hpp`) *generates* a fresh tet
 *    mesh from a signed-distance lattice (surface-in, volume-out); it discards
 *    the input's tets entirely rather than improving them.
 *
 * `optimize_volume` is the missing third member: it takes an existing tet mesh
 * and raises its worst-element quality by alternating the ODT vertex
 * relocation above with quality-improving topological **flips** (2-3 and 3-2),
 * so both the vertex positions and the connectivity change. That is what makes
 * it *remeshing* rather than *smoothing*, and it is the resolution-preserving
 * quality optimiser `refine`/`decimate_volume` are not (`refine` subdivides
 * the input's own possibly-bad cells; `decimate_volume` can only remove them).
 *
 * ### Why predicate-free (the in-posture design)
 *
 * The roadmap rejects the literal Delaunay/ODT-remeshing method because a
 * robust 3D Delaunay kernel needs in-sphere/orientation predicates, exactly
 * where this project's dependency-free posture stops paying (Geogram is a
 * large dependency; TetGen/CGAL's `Mesh_3` are AGPL/GPL). This operation needs
 * **none** of that. A flip is applied iff, using only the *signed volume* of
 * the candidate tets (`detail::cell_volume_from_corners`/`detail::det3`):
 *
 *   1. every new tet is non-degenerate and the local configuration is convex
 *      (a pure signed-volume test -- NO in-sphere test), and
 *   2. the **minimum** quality (scaled Jacobian) over the new tets strictly
 *      exceeds the minimum over the tets it replaces, by `mMinImprovement`.
 *
 * Criterion 2 is Freitag & Ollivier-Gooch's improvement rule: because every
 * accepted flip strictly raises a bounded quantity (the worst incident
 * quality), the process is **monotone in worst quality** and therefore
 * terminates -- no Delaunay optimality argument, and no predicate, is needed.
 *
 * ### Boundary is invariant by construction
 *
 * A 2-3 flip acts only on an *interior* triangular face (shared by two tets);
 * a 3-2 flip only on an *interior* edge (all its incident faces interior).
 * Neither ever touches a boundary face. Combined with `mPreserveBoundary`
 * pinning boundary vertices during relocation, the output's boundary surface
 * is **byte-identical** to the input's: watertight in => watertight out, with
 * none of the coincident-edge risk `remesh_volume`'s surface warp carries. It
 * is both a design guarantee and a test oracle.
 *
 * ### What survives
 *
 * The **point set is invariant** -- relocation moves points, flips only
 * reconnect them, and no point is ever added or removed -- so `point_data`,
 * `field_data` and named **Point** regions carry through unchanged. `cell_data`
 * has no correspondence across a flip (a 2-3 flip replaces two cells with
 * three) and is **dropped with a warning**, as are named Cell and Side
 * regions; the output is a single `tetra` block. This is more generous than
 * `remesh_volume` (which drops point data too, having a genuinely new point
 * set), and it is honest: what is preserved is exactly what a preserved point
 * set can preserve.
 *
 * ### Determinism
 *
 * Per-tet quality and the face/edge adjacency of each sweep are built in
 * `parallel_for` into disjoint slots; the flip-application loop is **serial**,
 * in ascending (cell, local face/edge) order, because a flip mutates shared
 * incidence (`decimate_volume`'s greedy-loop reasoning). Output is
 * byte-identical across the three mesh backends and across thread counts.
 *
 * **C++-core only, no numpy fallback.** The flip acceptance is a discrete
 * branch on the sign of a volume and on a near-tie quality comparison; a second
 * independent implementation could land on the other side of such a tie and
 * diverge into a different connectivity, not a last-ulp difference -- the same
 * reasoning `subdivide`/`agglomerate`/`decimate_volume`/`remesh_volume` state.
 * `_optimize_volume.py` raises `NotImplementedError` by name when the compiled
 * core is unavailable, for any input.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. An operation, not a file format -- deliberately not in the registry.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `optimize_volume`.
struct OptimizeVolumeOptions {
    /// How many optimisation sweeps to run; each sweep is one ODT relocation
    /// pass followed by one flip pass. The loop stops early once a sweep both
    /// moves no vertex past its tolerance and accepts no flip (a fixed point).
    int mMaxIterations = 10;

    /// Run the ODT vertex-relocation half of each sweep. With this off and
    /// `mFlip` on, the operation is a pure connectivity optimiser (flips only);
    /// with both off it is a no-op returning an equivalent mesh.
    bool mRelocate = true;

    /// Run the topological-flip half of each sweep. With this off and
    /// `mRelocate` on, the operation reduces to ODT *smoothing*
    /// (`SmoothMethod::Odt`) -- kept as a knob so the two halves are
    /// separately testable, not as a recommended mode.
    bool mFlip = true;

    /// Pin boundary vertices during relocation (passed straight to
    /// `smooth`'s `mFixBoundary`). Default true: the flips never touch the
    /// boundary, so with this on the boundary surface is exactly preserved.
    /// Setting it false lets boundary vertices drift off the surface, the same
    /// caveat `smooth` documents -- there is no surface re-projection here.
    bool mPreserveBoundary = true;

    /// The strict quality gain, in scaled-Jacobian units, a flip must deliver
    /// to be accepted: `min(quality of new tets) > min(quality of replaced
    /// tets) + mMinImprovement`. A small positive value both guarantees
    /// progress (monotone worst quality) and prevents cycling on ties.
    double mMinImprovement = 1e-6;

    /// Optional caller-supplied pin mask for the relocation half: either empty
    /// (no extra pins) or of length `NumPoints()`, a non-zero entry pinning
    /// that vertex. Forwarded to `smooth`'s `mFrozen`. Deliberately not exposed
    /// on the flat C ABI, a documented gap like `smooth`'s / `decimate`'s own
    /// `frozen`.
    std::vector<std::uint8_t> mFrozen;
};

/// The result of `optimize_volume`: the improved mesh plus what the run did.
struct OptimizeVolumeResult {
    /// The optimised mesh: a single `tetra` block, the same point set moved,
    /// `point_data`/`field_data`/Point regions carried, `cell_data`/Cell/Side
    /// regions dropped.
    Mesh mMesh;

    /// Total flips accepted (`mNum23Flips + mNum32Flips`).
    std::int64_t mNumFlips = 0;
    /// 2-3 flips accepted (an interior face's two tets became three).
    std::int64_t mNum23Flips = 0;
    /// 3-2 flips accepted (an interior edge's three tets became two).
    std::int64_t mNum32Flips = 0;

    /// Vertices whose final position differs from their input position by more
    /// than a bbox-relative tolerance (from the relocation half).
    std::int64_t mNumVerticesMoved = 0;

    /// Tets in the output mesh.
    std::int64_t mNumTets = 0;

    /// Minimum scaled Jacobian over the input tets (orientation-normalised: an
    /// input tet is measured on its positively-oriented reordering, so a
    /// negatively-stored but well-shaped tet reads as good, and a degenerate
    /// one as 0).
    double mMinQualityBefore = 0.0;
    /// Minimum scaled Jacobian over the output tets, same convention. Never
    /// less than `mMinQualityBefore`.
    double mMinQualityAfter = 0.0;
};

/**
 * @brief ODT-remesh a tetrahedral mesh: relocate vertices and flip
 *        connectivity to raise the worst element quality.
 * @param rMesh a tetra-only volume mesh (one or more `tetra` blocks).
 * @param rOptions the sweep count, which halves to run, boundary policy and
 *        the flip-acceptance threshold.
 * @return the improved mesh (single `tetra` block) plus flip/move counters and
 *         the before/after worst quality.
 * @throws std::invalid_argument when the mesh contains a non-`tetra` 3D block
 *         (run `convert_cells(mode='simplexify')` first), a ragged/polyhedron
 *         block, a non-3D block alongside the tets (drop it via `split`), or no
 *         tetra block at all -- the same tet-only scope `smooth`'s ODT method
 *         and `decimate_volume` enforce.
 */
MESHIOPLUSPLUS_API OptimizeVolumeResult optimize_volume(const Mesh& rMesh,
                                                        const OptimizeVolumeOptions& rOptions = {});

}  // namespace meshioplusplus
