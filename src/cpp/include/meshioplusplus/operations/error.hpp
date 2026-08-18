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
 * @file operations/error.hpp
 * @brief A posteriori error estimation: a recovery-based indicator plus a
 * marking pass — the piece that closes the adaptive loop. `gradient` already
 * differentiates a field and selective `refine`'s `--where` already consumes
 * any scalar `cell_data` predicate; this is what produces one.
 *
 * `estimate_error` implements the Zienkiewicz-Zhu (ZZ) superconvergent-patch
 * recovery estimator: differentiate the field once per cell (`gradient`,
 * Green-Gauss, `Cell` location), recover a smoothed nodal gradient by
 * averaging onto the points and back onto the cells
 * (`cell_data_to_point_data`/`point_data_to_cell_data`, `Measure`-weighted),
 * and take the recovered-minus-raw difference as the local error in the
 * gradient. **This is composition, not a new numerical kernel** — the
 * recovery operator *is* the existing averaging round trip, reused rather
 * than reimplemented, exactly as `gradient(mLocation=Point)` itself composes
 * `cell_data_to_point_data`.
 *
 * ### The indicator
 *
 * For cell `K` with raw (unrecovered) gradient `g_K` and recovered gradient
 * `G_K` (both possibly multi-component, over all `3 * nc` derivative
 * components of the differentiated field):
 * `eta_K = sqrt(|measure_K| * sum((G_K - g_K)^2))`. That is the standard ZZ
 * local indicator in the energy-like (H1-seminorm) norm — a genuine per-cell
 * error unit, not a bare gradient difference, which is what makes indicators
 * from cells of very different size comparable. `measure_K` reuses
 * `detail::cell_measure` (`detail/data_ops.hpp`) rather than recomputing a
 * volume — the same primitive `data_average`'s `CellPointWeight::Measure`
 * already uses for the recovery weighting, so the two agree by construction.
 * `mGlobalError` is `sqrt(sum_K eta_K^2)` over the evaluable cells — the
 * estimated global error in that same norm.
 *
 * ### Marking
 *
 * A marking pass turns the continuous indicator into a boolean `cell_data`
 * array — `error:marked`, one array per block, Int64 0/1 — so that **`refine`
 * needs no change at all**: the intended use is
 * `refine --where "error:marked > 0.5"`. Three policies:
 * - `Absolute`: mark every evaluable cell whose indicator exceeds
 *   `mMarkingValue`.
 * - `Fraction`: mark the `mMarkingValue` fraction (0, 1] of evaluable cells
 *   with the largest indicator, `k = size_t(mMarkingValue * numEvaluable +
 *   0.5)` truncated — the same `int64(ratio * F + 0.5)` convention `decimate`
 *   uses for its own target-count rounding, never `llround`, whose tie rule a
 *   numpy twin cannot match.
 * - `Dorfler`: the Dörfler / bulk-chunk criterion — mark the smallest
 *   descending-indicator prefix whose cumulative `eta_K^2` covers a
 *   `mMarkingValue` fraction (0, 1] of the total.
 *
 * Ties in the `Fraction`/`Dorfler` descending sort break on ascending global
 * cell index (`detail::block_bases`' numbering), for a result independent of
 * how the per-cell work happened to schedule — a strict weak order, not a
 * convenience.
 *
 * ### Contracts
 *
 * - **Input must be `point_data`** — the same restriction `gradient` states,
 *   for the same reason: a `cell_data` field is piecewise constant and has no
 *   derivative to recover.
 * - **NaN, never a guess.** A cell `gradient` could not differentiate (an
 *   unsupported type, wrong dimension, degenerate geometry), a cell whose
 *   recovery neighbourhood contributed no finite value, or a cell whose
 *   `measure_K` cannot be computed reads `NaN` in `error:zz` and **0** (never
 *   NaN) in `error:marked`. Such cells are counted in `mNumSkipped` and
 *   excluded from `mGlobalError` and from every marking policy's cell count.
 * - **Output.** `error:zz` (Float64, one array per block) is always
 *   attached; `error:marked` (Int64, one array per block) is attached only
 *   when `mMarking != None`. Default names are fixed
 *   (`kErrorZzName`/`kErrorMarkedName`); `mOutputName`/`mMarkedName`
 *   override them.
 * - Geometry, connectivity, regions, property sets and every existing data
 *   array pass through bit-identically.
 *
 * Output is byte-identical across the three mesh backends and across thread
 * counts. Across the C++/numpy boundary it matches the numpy twin in
 * `_error.py` to within a numeric tolerance rather than bit-for-bit — the
 * same, already-accepted precedent `data_average`'s own `Measure` weighting
 * has (`test_data_location.py::test_cpp_matches_python` compares it with
 * `np.allclose`, not exact equality), inherited here because the recovery
 * step composes that exact weighting. The raw (pre-recovery) gradient stays
 * bit-exact, since `gradient` itself needs no such tolerance. The twin raises
 * `NotImplementedError` when the mesh has a ragged or polyhedral block,
 * exactly as `_convert_cells.py`'s simplexify twin does, because
 * `detail::cell_measure` runs `orient_rings`' winding repair (a discrete
 * branch on a sign) for a closed 3D cell, tabulated or not, and the Python
 * twin's `_cell_measures` does not attempt to replicate that repair — and,
 * separately, when the mesh has a quadratic 3D cell type (`tetra10`,
 * `hexahedron20`/`27`, `wedge15`), because `_cell_measures`' face table
 * (shared with `data_average`'s own recovery weighting) only covers the four
 * linear 3D types; `detail::cell_measure` handles the quadratic ones fine
 * (via the same corner reduction `gradient` itself uses), so this is a gap in
 * the *reference*, not the core, and is scoped as narrowly as that.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is
 * not in the format registry.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Which error-estimator family to use. Only `Zz` exists today; the enum
/// leaves room for a jump/Kelly estimator without an ABI-breaking rename.
enum class ErrorMethod {
    Zz,  ///< Zienkiewicz-Zhu superconvergent-patch gradient recovery.
};

/// How the continuous indicator is turned into a boolean `error:marked` array.
enum class ErrorMarking {
    None,      ///< Estimate only; `error:marked` is not attached.
    Absolute,  ///< Mark cells whose indicator exceeds `mMarkingValue`.
    Fraction,  ///< Mark the largest `mMarkingValue` fraction (0, 1] of cells.
    Dorfler,   ///< Mark the smallest prefix covering `mMarkingValue` (0, 1]
               ///< of the total squared error.
};

/// Default output array names.
inline constexpr const char* kErrorZzName = "error:zz";
inline constexpr const char* kErrorMarkedName = "error:marked";

/// Options for `estimate_error`.
struct ErrorOptions {
    /// Name of the `point_data` array to estimate the recovered-gradient
    /// error of. Required.
    std::string mArrayName;
    /// Which estimator family to use. Only `Zz` exists today.
    ErrorMethod mMethod = ErrorMethod::Zz;
    /// How to turn the indicator into `error:marked`. `None` skips marking.
    ErrorMarking mMarking = ErrorMarking::None;
    /// Meaning depends on `mMarking`: ignored for `None`; an absolute
    /// indicator threshold for `Absolute`; a fraction in (0, 1] of cells for
    /// `Fraction`; the Dörfler bulk fraction theta in (0, 1] for `Dorfler`.
    double mMarkingValue = 0.0;
    /// Indicator output name; empty selects `error:zz`.
    std::string mOutputName;
    /// Marking output name; empty selects `error:marked`. Ignored when
    /// `mMarking == None`.
    std::string mMarkedName;
    /// When false, an output name that already exists is an error rather
    /// than being overwritten.
    bool mOverwrite = false;
};

/// The mesh with the produced arrays, plus what the estimate covers.
struct ErrorResult {
    /// The input mesh with `error:zz` (and, if marking, `error:marked`)
    /// added; everything else unchanged.
    Mesh mMesh;
    /// `sqrt(sum of eta_K^2 over evaluable cells)` — the estimated global
    /// error in the recovered-gradient norm. `0` when no cell is evaluable.
    double mGlobalError = 0.0;
    /// Cells that could not be evaluated: an unsupported/degenerate
    /// `gradient` cell, a recovery neighbourhood with no finite
    /// contribution, or an unmeasurable cell. `error:zz` reads NaN there.
    std::int64_t mNumSkipped = 0;
    /// Cells marked in `error:marked`. Always 0 when `mMarking == None`.
    std::int64_t mNumMarked = 0;
};

/**
 * @brief Estimates the per-cell recovered-gradient error of a `point_data`
 * field and optionally marks cells for refinement.
 * @param rMesh the source mesh (unmodified).
 * @param rOptions which array, method, marking policy and output names.
 * @return the mesh with the produced array(s), the global error estimate and
 *         the skip/mark counters.
 * @throws std::invalid_argument on an empty or unknown array name, on a
 *         `cell_data` name (message names the fix), on a mesh with no cell of
 *         topological dimension 1 or more, on `mMarkingValue` outside (0, 1]
 *         for `Fraction`/`Dorfler`, or on an output-name collision when
 *         `mOverwrite` is false.
 */
MESHIOPLUSPLUS_API ErrorResult estimate_error(const Mesh& rMesh, const ErrorOptions& rOptions);

/**
 * @brief Parses an estimator-method name.
 * @param rName `zz` (or empty).
 * @return the matching `ErrorMethod`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API ErrorMethod error_method_from_name(const std::string& rName);

/**
 * @brief Parses a marking-policy name.
 * @param rName one of `none` (or empty), `absolute`/`abs`, `fraction`/`frac`,
 *        or `dorfler`/`bulk`.
 * @return the matching `ErrorMarking`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API ErrorMarking error_marking_from_name(const std::string& rName);

}  // namespace meshioplusplus
