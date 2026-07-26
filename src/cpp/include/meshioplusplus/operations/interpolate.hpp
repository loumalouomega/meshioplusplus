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
 * @file interpolate.hpp
 * @brief Cross-mesh field transfer: sample a SOURCE mesh's data arrays onto a
 * TARGET mesh. This is an operation, not a file format — deliberately not in
 * the registry — and the repo's first two-mesh op that transfers data (diff
 * compares, merge concatenates; neither resamples).
 *
 * `interpolate(source, target, options)` returns a **new mesh that is a copy
 * of the target** — geometry, connectivity, and the target's own data arrays
 * preserved exactly — with the requested source arrays sampled onto it:
 *
 * - source `point_data` is sampled **at the target's points**;
 * - source `cell_data` is sampled **at the target's cell centroids**, and is
 *   always transferred by nearest source-cell centroid regardless of the
 *   method — there is no meaningful barycentric interpolation of a piecewise
 *   constant field. Cross-location transfer (point -> cell and back) is out
 *   of scope; compose with the `data to-cell` / `to-point` operations.
 *
 * Methods:
 *
 * - `Nearest` (default; robust, any cell type — needs only source points):
 *   each target sample takes the value of the nearest source point, found via
 *   an expanding-shell search over a bucket-grid spatial hash
 *   (`detail/spatial_hash.hpp`, shared with merge's weld) — never O(N^2).
 *   Deterministic tie-break: smallest squared distance, then lowest source
 *   point index. Output arrays keep the source array's dtype (values are
 *   copied bit-for-bit, never round-tripped through double).
 * - `Barycentric` (linear; exact on a linear field): the source is first
 *   simplexified via `convert_cells` (so every cell is a triangle or
 *   tetrahedron and containment + weights are exact and uniform — on a
 *   quad/hex source this is therefore the **simplex-linear** interpolant, not
 *   true bi-/trilinear). Candidate simplices come from a bucket grid over the
 *   simplices' bounding boxes; containment is tested by the signs of the
 *   barycentric coordinates (sub-simplex volume ratios, tolerance `-1e-12` —
 *   dimensionless, hence scale-invariant), and the **first containing simplex
 *   in ascending simplex index** wins. Only the source's maximum topological
 *   dimension is sampled: tetrahedra when any exist, else triangles — whose
 *   barycentric coordinates are computed **in the xy-plane** (a planar
 *   assumption; for a curved surface embedded in 3D use `Nearest`). Output is
 *   always Float64 (a weighted mean is not integral). A target point in no
 *   simplex receives `mDefaultValue` (all components) — with one `log::warn`
 *   naming the count — unless `mExtrapolate` is set, in which case it falls
 *   back to the nearest source point's value.
 *
 * Array selection (`mArrays`): empty means every source `point_data` array, in
 * sorted name order; cell_data transfers only when named explicitly. A named
 * array present in both locations transfers both; present in neither throws.
 * If a transferred name already exists on the target, `mOnConflict` decides:
 * `Error` (default) throws, `Overwrite` replaces, `Suffix` writes to
 * `name + "_interp"` (throwing if that name is taken too).
 *
 * The target's `point_sets`/`cell_sets` ride through untouched (nothing is
 * renumbered); source sets are **not** transferred; `mesh.info` and
 * `gmsh_periodic` are not carried (they never reach the C++ core).
 *
 * Determinism: every target value is computed independently from the immutable
 * source (no cross-target dependence), the grids are built by a `parallel_for`
 * key fill into disjoint slots followed by a **serial** ascending insert
 * (surface.cpp's phase-split idiom), and every reduction — nearest scan,
 * barycentric weight accumulation (corner by corner in connectivity order) —
 * runs in a fixed index order. The grid cell size is derived without `cbrt`
 * (the smallest integer `R` with `R*R*R >= n`, then one division of the
 * largest bbox extent), so the pure-numpy fallback reproduces it bit-for-bit;
 * output is byte-identical across the three mesh backends, across thread
 * counts, and across the C++-core / numpy-fallback boundary
 * (`tests/python/test_interpolate.py::test_cpp_matches_python`).
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How a target sample draws its value from the source.
enum class InterpolateMethod {
    Nearest,      ///< Nearest source point's value (piecewise constant).
    Barycentric,  ///< Linear interpolation in the containing source simplex.
};

/**
 * @brief Parses an interpolation method name.
 * @param rName One of `"nearest"`, `"barycentric"` (case-sensitive, as
 *        elsewhere in the operations layer).
 * @return The matching enumerator.
 * @throws std::invalid_argument if the name is not recognised.
 */
MESHIOPLUSPLUS_API InterpolateMethod interpolate_method_from_name(const std::string& rName);

/// What to do when a transferred array name already exists on the target.
enum class InterpolateConflict {
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
MESHIOPLUSPLUS_API InterpolateConflict interpolate_conflict_from_name(const std::string& rName);

/// Options for `interpolate`.
struct InterpolateOptions {
    /// How target samples draw their values (cell_data is always
    /// nearest-source-cell-centroid regardless — see the file doc block).
    InterpolateMethod mMethod = InterpolateMethod::Nearest;

    /// Source array names to transfer. Empty = every source point_data array
    /// (sorted name order); cell_data only transfers when named explicitly. A
    /// name found in neither source location throws.
    std::vector<std::string> mArrays;

    /// `Barycentric` only: a target point outside the source domain falls back
    /// to the nearest source point's value instead of `mDefaultValue`.
    bool mExtrapolate = false;

    /// `Barycentric` only: the value written (to every component) for a target
    /// point outside the source domain when `mExtrapolate` is off.
    double mDefaultValue = 0.0;

    /// What to do when a transferred name already exists on the target.
    InterpolateConflict mOnConflict = InterpolateConflict::Error;
};

/**
 * @brief Samples data arrays from @p rSource onto @p rTarget.
 * @param rSource The mesh whose data is sampled (never modified).
 * @param rTarget The mesh receiving the samples (never modified).
 * @param rOptions Method, array selection, outside-domain and conflict policy.
 * @return A copy of the target with the requested source arrays attached.
 * @throws std::invalid_argument on an empty source, an unknown array name, a
 *         name conflict under `Error` (or a taken suffix under `Suffix`), a
 *         `Barycentric` source with no triangle/tetrahedron cells after
 *         simplexification, or inconsistent per-block components in a
 *         requested cell_data array.
 */
MESHIOPLUSPLUS_API Mesh interpolate(const Mesh& rSource, const Mesh& rTarget, const InterpolateOptions& rOptions = {});

}  // namespace meshioplusplus
