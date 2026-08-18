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
 * @file conservative_interpolate.hpp
 * @brief Mass-preserving cross-mesh field transfer: sample a SOURCE mesh's
 * `cell_data`/`point_data` onto a TARGET mesh so that, over the region the two
 * meshes share, `sum(target value * target measure)` equals
 * `sum(source value * source measure)` — the property `interpolate`'s
 * `Barycentric` mode does not have, since a pointwise sampler has no notion of
 * "how much of the source region a target sample stands for". CFD/FEM
 * solver-coupling and successive-remeshing workflows need exactly this
 * guarantee; `interpolate`'s modes remain the right tool for genuinely
 * pointwise sampling (probing a field at a set of coordinates, visualisation
 * resampling, ...).
 *
 * `conservative_interpolate(source, target, options)` returns a **new mesh
 * that is a copy of the target** — geometry, connectivity, and the target's
 * own data arrays preserved exactly — with the requested source arrays
 * transferred onto it. Unlike `interpolate`, empty `mArrays` means *every*
 * source `point_data` **and** `cell_data` array: there is one algorithm
 * regardless of location, so `interpolate`'s "cell_data only when named"
 * special case (a consequence of cell_data always being the coarser,
 * nearest-centroid path there) does not apply here.
 *
 * ### Algorithm
 *
 * Both meshes are first reduced to simplices via
 * `convert_cells(Simplexify, record_parent_ids=true)` — the same call
 * `interpolate`'s `Barycentric` mode already makes, which is what lets this
 * operation accept ragged and polyhedron blocks for free (Simplexify already
 * fans them into triangles/tetrahedra via a shipped, tested path) rather than
 * needing a general polygon/polyhedron clipper. Only the maximum
 * topological-dimension simplices participate: triangles in 2D, tetrahedra in
 * 3D — `interpolate`'s own rule. `source` and `target` must share the same
 * maximum topological dimension; there is no cross-dimension remap.
 *
 * For every pair of overlapping simplices (found via a bucket-grid spatial
 * hash over bounding boxes, `detail/spatial_hash.hpp`, exactly as
 * `interpolate`'s barycentric candidate search already does), the exact
 * overlap area (2D, a Sutherland-Hodgman convex polygon clip) or volume (3D, a
 * clip of the target tetrahedron against the source tetrahedron's four
 * half-spaces, since both operands are always tetrahedra) is computed and used
 * as a weight. A target cell's value is the overlap-measure-weighted mean of
 * every overlapping source cell's value — the discrete form of the
 * measure-weighted overlap integral a conservative remap computes. A target
 * cell whose covered fraction (against its own unclipped measure) falls below
 * a small relative tolerance is filled with `mDefaultValue`; the count is
 * reported through one aggregated `log::warn` per call. Output arrays are
 * always Float64 — a measure-weighted mean is not integral, the same
 * convention `interpolate`'s `Barycentric` mode and `data_average`'s
 * averaging operations already use.
 *
 * `cell_data` is transferred by this algorithm directly. `point_data` is
 * transferred by **composition**: `point_data_to_cell_data` lumps the source
 * array onto a cell proxy (the unweighted mean of a cell's own corner
 * values), the same overlap-measure algorithm remaps that proxy onto the
 * target's cells exactly, and `cell_data_to_point_data(weight=Measure)`
 * distributes the result back onto the target's points. Conservation is
 * **exact only for the middle step** — the two lumping steps that sandwich it
 * are each already-documented approximations of their own (an unweighted
 * corner mean, and a measure-weighted point mean), so the overall
 * `point_data` path is a layered approximation, not exact nodal/FEM
 * conservation. No dual-cell/control-volume machinery is built for points.
 *
 * This operation deliberately reports **no** integral/conservation
 * diagnostic of its own — measuring how well conservation held on a given
 * mesh (and by how much a partial-coverage or point_data-composed result
 * drifted) is the explicit job of the separate "field integration"
 * companion (`data integrate`/`total`/`mean` reductions, `doc/roadmap.md`),
 * not duplicated here.
 *
 * There is no `mExtrapolate` flag (unlike `interpolate`): a silent
 * nearest-source-cell fallback for uncovered cells would break the exact
 * conservation guarantee for precisely the boundary cells a caller is most
 * likely to reach for it, on the one operation whose entire purpose is that
 * guarantee.
 *
 * Determinism follows the same recipe as `interpolate`: independent
 * per-target-simplex work runs in `parallel_for` with each simplex's own
 * candidate accumulation in a fixed (deduplicated, ascending source-simplex
 * index) order, while the many-to-few scatter from simplices back onto the
 * target's original cells runs **serially** in ascending target-simplex
 * index — floating-point addition is not associative, so this scatter cannot
 * be parallelised without making the result thread-count-dependent (the same
 * rule `operations/data_average.cpp`'s cell-to-point averaging documents).
 * Output is therefore byte-identical across the three mesh backends and
 * across thread counts.
 *
 * There is **no pure-Python fallback**: the 3D clip kernel is a
 * discrete-branch geometric algorithm (half-space in/out classification,
 * cutting-plane chord deduplication, angle-sorted cap triangulation) of
 * exactly the class `subdivide`/`agglomerate`/`decimate_volume` already
 * document as unsafe to give a second, independently-written implementation —
 * near a degenerate or tangent overlap the two could silently disagree.
 * `meshioplusplus.conservative_interpolate` raises `NotImplementedError` by
 * name when `_core` is unavailable.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// What to do when a transferred array name already exists on the target.
/// A fresh enum, not a reuse of `interpolate.hpp`'s `InterpolateConflict` —
/// identical in shape and semantics, kept independent so the two operations'
/// headers/ABIs cannot be coupled by a shared type (the `decimate_volume`
/// precedent of not reaching into `decimate.hpp`).
enum class ConservativeInterpolateConflict {
    Error,      ///< Throw `std::invalid_argument` (the default).
    Overwrite,  ///< Replace the target's array.
    Suffix,     ///< Write to `name + "_interp"` instead (throws if taken too).
};

/**
 * @brief Parses a conflict-policy name.
 * @param rName One of `"error"`, `"overwrite"`, `"suffix"`.
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API ConservativeInterpolateConflict
conservative_interpolate_conflict_from_name(const std::string& rName);

/// Options for `conservative_interpolate`.
struct ConservativeInterpolateOptions {
    /// Source array names to transfer. Empty = every source point_data AND
    /// cell_data array (sorted name order); a name found in neither source
    /// location throws.
    std::vector<std::string> mArrays;

    /// The value written (to every component) for a target cell whose covered
    /// fraction of source overlap is below the coverage tolerance.
    double mDefaultValue = 0.0;

    /// What to do when a transferred name already exists on the target.
    ConservativeInterpolateConflict mOnConflict = ConservativeInterpolateConflict::Error;
};

/**
 * @brief Conservatively (measure-weighted) samples data arrays from
 * @p rSource onto @p rTarget.
 * @param rSource The mesh whose data is sampled (never modified).
 * @param rTarget The mesh receiving the samples (never modified).
 * @param rOptions Array selection, default value and conflict policy.
 * @return A copy of the target with the requested source arrays attached.
 * @throws std::invalid_argument on an empty source, mismatched maximum
 *         topological dimensions between the two meshes, an unknown array
 *         name, a name conflict under `Error` (or a taken suffix under
 *         `Suffix`), or when neither mesh has any triangle/tetrahedron
 *         simplex after simplexification.
 */
MESHIOPLUSPLUS_API Mesh conservative_interpolate(
    const Mesh& rSource, const Mesh& rTarget, const ConservativeInterpolateOptions& rOptions = {});

}  // namespace meshioplusplus
