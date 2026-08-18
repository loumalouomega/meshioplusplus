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
 * @file detail/decimate_common.hpp
 * @brief Quadric-error-metric machinery shared between surface `decimate`
 * (`operations/decimate.cpp`) and volume `decimate_volume`
 * (`operations/decimate_volume.cpp`), hoisted verbatim out of the former --
 * the placement solve, the error form, the normal-flip check, and the
 * vertex-ring link-condition primitives.
 *
 * The hoist exists so the boundary half of volume decimation reuses
 * already-tested floating-point code instead of a second, silently-divergent
 * transcription: the 3x3 cofactor solve and its ill-conditioning fallback in
 * particular are exactly the kind of tricky arithmetic that drifts silently
 * between two hand-written copies. `decimate.cpp`'s own 16-test suite is the
 * regression guard that this relocation is behaviour-preserving -- every name
 * here is identical to its pre-hoist counterpart, only the namespace and
 * linkage changed.
 *
 * Also carries the generic flat-triangle-mesh plumbing (`DecimFaces`,
 * `DecimCsr`, the vertex/face CSR builder, per-face plane quadrics, per-vertex
 * quadric accumulation, feature-vertex marking): volume decimation's boundary
 * vertices are, topologically, an ordinary triangle mesh (a tet mesh's own
 * skin), so its quadric/feature setup reuses this unchanged rather than
 * re-deriving it against the volume mesh's boundary faces. `decim_build_faces`
 * (ties a `DecimFaces` to a *simplexified `Mesh`*'s own blocks) and
 * `decim_build_edges` (the open-boundary-of-a-surface-PATCH edge test) stay
 * private to `decimate.cpp` -- neither concept applies to a solid's closed
 * outer skin, which is what `decimate_volume.cpp` builds its `DecimFaces`
 * from instead (via `detail::build_global_faces`).
 *
 * Every floating-point expression here is transcribed token for token into
 * `_decimate.py`/`_decimate_volume.py` (per each operation's own docs), and
 * none of it may be contracted into FMAs across the C++/numpy boundary.
 *
 * Free functions in `meshioplusplus::detail`, built on plain flat buffers (no
 * Mesh/uniform-API dependency), so both callers can build their own face/edge
 * tables around it.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/operations/decimate.hpp"

namespace meshioplusplus {
namespace detail {

/// A flat triangle mesh's connectivity: 3 corners per face, block-major, plus
/// each "block"'s first global face id (a single-entry `{0}` when the caller
/// has no natural block subdivision, e.g. a volume mesh's boundary skin).
struct DecimFaces {
    std::vector<std::int64_t> mCorners;    // 3 * num_faces
    std::vector<std::int64_t> mBlockBase;  // per block, global id of face 0
    std::size_t mNumFaces = 0;
};

/// A vertex -> incident-face CSR (counting sort; rows ascend in face index).
struct DecimCsr {
    std::vector<std::int64_t> mXadj;
    std::vector<std::int64_t> mAdj;
};

/// Builds `rFaces`' vertex -> face incidence over `n` points.
MESHIOPLUSPLUS_API DecimCsr decim_vertex_faces_csr(const DecimFaces& rFaces, std::size_t n);

/**
 * @brief Per-face area-weighted plane quadrics and unit normals.
 *
 * `rQuadK` gets 10 entries per face, in the fixed order
 * `[aa,ab,ac,ad,bb,bc,bd,cc,cd,dd]`; `rNormals` gets 3 (zero for a degenerate
 * face). Area weighting (`w = |n|/2`) makes the accumulated per-vertex
 * quadric independent of how a planar polygon was fanned into triangles.
 */
MESHIOPLUSPLUS_API void decim_face_planes(const DecimFaces& rFaces, const std::vector<double>& rXyz,
                                          std::vector<double>& rQuadK,
                                          std::vector<double>& rNormals);

/**
 * @brief Per-vertex quadric: the sum of the incident faces' quadrics in
 * ascending face index (the CSR rows are face-ascending by construction).
 *
 * FP addition is not associative, so this fixed order is the pin across mesh
 * backends, thread counts and the numpy fallback.
 */
MESHIOPLUSPLUS_API std::vector<double> decim_accumulate_quadrics(const DecimCsr& rCsr,
                                                                 std::size_t n,
                                                                 const std::vector<double>& rQuadK);

/**
 * @brief Pins vertices whose incident face unit normals pairwise differ by
 * more than the feature angle (`CosThreshold = cos(angle)`).
 *
 * Every face participates, not just boundary facets: the creases of a closed
 * surface are interior. O(d^2) in the valence; each iteration writes only its
 * own slot, so this is safe to call under `parallel_for`.
 */
MESHIOPLUSPLUS_API void decim_mark_features(const DecimCsr& rCsr, std::size_t n,
                                            const std::vector<double>& rNormals,
                                            double CosThreshold,
                                            std::vector<std::uint8_t>& rPinned);

/// x^T Q x for the homogeneous point (x, y, z, 1), `q` the 10-entry
/// `[aa,ab,ac,ad,bb,bc,bd,cc,cd,dd]` quadric. The literal parenthesization is
/// the parity contract with the numpy twins.
MESHIOPLUSPLUS_API double decim_quadric_error(const double* q, double x, double y, double z);

/// A placement result: the surviving vertex's position and quadric error.
struct DecimPlaced {
    double mX[3] = {0.0, 0.0, 0.0};
    double mErr = 0.0;
};

/// Read-only state a placement solve needs; shared by a parallel seeding pass
/// and a serial greedy loop.
struct DecimPlaceCtx {
    const std::vector<double>* mpXyz = nullptr;
    const std::vector<double>* mpQ = nullptr;
    const std::vector<std::uint8_t>* mpPinned = nullptr;
    DecimatePlacement mPlacement = DecimatePlacement::Optimal;
};

/**
 * @brief Position and quadric error of collapsing edge (a, b).
 *
 * A pinned endpoint forces its own position; otherwise the requested
 * placement applies, with `Optimal` falling back to the midpoint when the 3x3
 * system is ill-conditioned (`|det| <= 1e-12 * max|A_ij|^3`). A vertex whose
 * quadric is exactly zero (no incident plane contributed anything -- the
 * combined quadric of two purely-interior volume-decimation vertices, for
 * instance) is a degenerate case of the SAME bound: it naturally routes to
 * the midpoint fallback with no special-casing needed.
 */
MESHIOPLUSPLUS_API DecimPlaced decim_place(const DecimPlaceCtx& rCtx, std::int64_t a,
                                           std::int64_t b);

/**
 * @brief Unnormalized triangle normal, substituting `pSub` for any corner
 * equal to `A` or `B` when given (the candidate survivor position).
 *
 * The sign of `n_before . n_after` is what a normal-flip guard branches on,
 * so normalization would be pure extra rounding.
 */
MESHIOPLUSPLUS_API void decim_face_normal(const std::vector<double>& rXyz,
                                          const std::int64_t* pCorners, std::int64_t A,
                                          std::int64_t B, const double* pSub, double pOut[3]);

/// Sorted unique vertex ring of `v`: the corners of its alive incident faces
/// (from `rVFaces`, indices into the flat 3-per-face `rCorners`), minus `v`
/// itself. Pure integer work, hence trivially deterministic.
MESHIOPLUSPLUS_API std::vector<std::int64_t> decim_vertex_ring(
    const std::vector<std::int64_t>& rVFaces, const std::vector<std::int64_t>& rCorners,
    std::int64_t v);

/// How many values two sorted vectors share (two-pointer sweep).
MESHIOPLUSPLUS_API std::size_t decim_count_common(const std::vector<std::int64_t>& rA,
                                                  const std::vector<std::int64_t>& rB);

/// Erase one value from a sorted vector, if present.
MESHIOPLUSPLUS_API void decim_sorted_erase(std::vector<std::int64_t>& rVec, std::int64_t Value);

/// Insert one value into a sorted vector, keeping it sorted.
MESHIOPLUSPLUS_API void decim_sorted_insert(std::vector<std::int64_t>& rVec, std::int64_t Value);

}  // namespace detail
}  // namespace meshioplusplus
