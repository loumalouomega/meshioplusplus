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
 * @file operations/gradient.hpp
 * @brief Field differential operators: the gradient, divergence and curl of a
 * `point_data` field on a mesh.
 *
 * This is the operation that lets a derived quantity be contoured (`isosurface`
 * on `|grad T|`, on vorticity) and that produces the error indicator the
 * selective `refine` consumes (`refine --where "grad_mag > ..."`).
 *
 * **Why it lives in the mesh-operations layer rather than the `data_*` bundle.**
 * The five `data_*` operations are *defined* by never touching geometry
 * (`operations/data_common.hpp`, `doc/data_operations.md`); this one consumes
 * and produces data arrays but **reads** geometry and topology — face areas,
 * cell volumes, cell-to-cell adjacency — so it belongs here, with the
 * `mPascalCase` option fields the geometry operations use. It is still
 * *reachable* as `meshioplusplus data gradient` in both CLIs, because that is
 * where a user looks for it.
 *
 * ### The two methods
 *
 * **Green-Gauss** (the default) applies the divergence theorem over the cell:
 * `grad f = (1/V) * sum_faces integral_face(f n dA)`. Each face is fanned into
 * triangles about the arithmetic mean of its corners, and each sub-triangle
 * contributes `A_j * f(centroid_j)`. That is **exact for a linear field on any
 * cell** — planar faces or not — because two faces sharing an edge contribute
 * oppositely-wound triangles on that edge, so the fan surface is closed, and
 * `∮ f n dA = V grad f` is a purely algebraic identity on a closed oriented
 * piecewise-linear surface. `V` is the exact signed volume *of that same fan
 * surface*, so numerator and denominator flip together and an inverted cell
 * yields the same gradient as its positively-wound twin.
 *
 * The corner-average fan apex is **forced, not merely convenient**: the
 * quadrature needs `f` at the apex, and for linear `f` the corner average is the
 * only point whose value is known exactly without invoking shape functions
 * (`mean(f(p_i)) == f(mean(p_i))`). Do not "simplify" it to the area centroid —
 * for a trapezoidal quad face the two differ and exactness is lost. A plain
 * corner-average face value (no fan at all) is wrong by 12.5% on the trapezoid
 * `(0,0),(4,0),(2,1),(0,1)` under `f = y`.
 *
 * On a 2D cell the same theorem runs over the corner ring with the in-plane
 * outward normal `t x N`; that is exact on a **planar** cell and first-order on
 * a warped quad in 3D, which has no well-defined area or normal to begin with.
 * The result is invariant under both reversal and cyclic rotation of the ring.
 *
 * **Least-squares** fits `f(x) ~ f_c + g . (x - x_c)` over the cells sharing at
 * least one node (`detail/cell_adjacency.hpp`, the same neighbour definition
 * `partition`'s ghost layers use), with `x_c`/`f_c` the arithmetic means of the
 * cell's corners — so `(x_c, f_c)` lies exactly on a linear field and the fit is
 * exact for one, under any positive weights. Weights are `1/|d|^2`, which makes
 * the normal matrix `sum d_hat d_hat^T`: dimensionless, immune to mesh grading,
 * and trivially scale-invariant. On a 2D mesh the offsets are projected into the
 * cell's plane, so exactness holds on a **planar** mesh and degrades to
 * first-order on a curved surface, where the projection discards a real part of
 * the offset. A degenerate neighbourhood (an isolated cell, a collinear strip)
 * falls back to Green-Gauss for that cell and is counted in
 * `GradientResult::mNumFallback` — never silently wrong, never NaN when a usable
 * answer exists.
 *
 * ### Contracts
 *
 * - **Input must be `point_data`.** A `cell_data` field is piecewise constant
 *   and has no derivative; naming one throws by name and points at
 *   `cell_data_to_point_data` (CLI `data to-point`), exactly as `isosurface`
 *   does for the same reason.
 * - **Shapes.** An `nc`-component input yields a gradient of `3 * nc`
 *   components, **flat** and row-major as `[component i][derivative j]` at index
 *   `i * 3 + j`: a scalar gives `(n, 3)`, a 3-vector gives `(n, 9)`.
 *   `mComponent` selects a single component and yields 3.
 *   `Divergence` yields 1 component, `Curl` always yields 3.
 * - **Divergence and curl** need 2 or 3 components; a 2-component field reads as
 *   `(u, v, 0)`, the same padding convention `detail::read_point` applies to 2D
 *   points. On a 2D mesh `d/dz` is unrecoverable and reads 0, so a 3-component
 *   field's divergence is its *in-plane* divergence and its curl is
 *   `(d_y w, -d_x w, d_x v - d_y u)`; only a **2**-component field yields a curl
 *   with a z-component alone.
 * - **Output is always `Float64`** — a derivative is not an integer — named
 *   `<input>:gradient` / `:divergence` / `:curl` unless `mOutputName` overrides.
 *   That is deliberately `name:suffix`, not the repo's usual `prefix:name`, so
 *   that everything derived from one field sorts next to it.
 * - **Location.** `Cell` (the natural home of a per-cell derivative, the
 *   default) or `Point`, the latter obtained by composing the existing
 *   `cell_data_to_point_data` averaging rather than reimplementing it. Point
 *   values are therefore exact to within rounding, not bit-exact, for a linear
 *   field.
 * - **Unsupported cells are NaN and counted**, never approximated: blocks below
 *   the mesh's own topological dimension, 3D types with no `detail::cell_faces`
 *   row (the 3D Lagrange family), ragged polygon and polyhedron blocks, and
 *   cells whose volume/area is degenerate relative to their own size. This is
 *   `compute_quality`'s NaN-where-not-applicable convention.
 * - **Non-finite input.** Unlike every `data_*` operation there is no
 *   `NanPolicy`: Green-Gauss has no reduction to exclude a value *from*, so a
 *   single non-finite corner poisons its whole cell. Stated here rather than
 *   silently diverging from `data_common.hpp`'s policy.
 * - Geometry, connectivity, regions, property sets and every existing data
 *   array pass through bit-identically.
 *
 * Output is byte-identical across the three mesh backends, across thread counts
 * and across the C++/numpy boundary: per-cell work is independent, and every
 * accumulation inside a cell runs in a fixed face -> fan-triangle -> neighbour
 * order because floating-point addition is not associative.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles under
 * every mesh backend. This is an operation, not a file format — it is not in the
 * format registry.
 */

// System includes
#include <cstdint>
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// Which differential operator to apply.
enum class GradientOperator {
    Gradient,    ///< `grad f`; `3 * nc` components, row-major `[i][j]` at `i * 3 + j`.
    Divergence,  ///< `div u`; one component. Needs a 2- or 3-component input.
    Curl         ///< `curl u`; three components. Needs a 2- or 3-component input.
};

/// How the per-cell derivative is reconstructed.
enum class GradientMethod {
    GreenGauss,   ///< Divergence theorem over the cell's own faces/edges.
    LeastSquares  ///< Linear fit over the node-sharing neighbour cells.
};

/// Default suffix of the produced array for each operator.
inline constexpr const char* kGradientSuffix = ":gradient";
inline constexpr const char* kDivergenceSuffix = ":divergence";
inline constexpr const char* kCurlSuffix = ":curl";

/// Options for `gradient`.
struct GradientOptions {
    /// Name of the `point_data` array to differentiate. Required.
    std::string mArrayName;
    /// Which operator to apply.
    GradientOperator mOperator = GradientOperator::Gradient;
    /// How to reconstruct the per-cell derivative.
    GradientMethod mMethod = GradientMethod::GreenGauss;
    /// Where the result is stored: `Point` or `Cell`. `Field` is an error.
    DataLocation mLocation = DataLocation::Cell;
    /// Output array name; empty selects the operator's default suffix.
    std::string mOutputName;
    /// Differentiate only this component of a multi-component input.
    ///
    /// Unset means **every** component — note this is the opposite of
    /// `IsosurfaceOptions::mComponent`, where unset means the row magnitude.
    /// Setting it together with `Divergence` or `Curl` is an error. At every
    /// binding boundary the "unset" sentinel is a negative integer.
    std::optional<int> mComponent;
    /// When false, an output name that already exists at the target location is
    /// an error rather than being overwritten.
    bool mOverwrite = false;
};

/// The differentiated mesh plus what could not be computed exactly.
struct GradientResult {
    /// The input mesh with the produced array added; everything else unchanged.
    Mesh mMesh;
    /// Cells that produced a NaN row: unsupported type, wrong dimension, or a
    /// degenerate measure. Counts cells, not (cell, component) pairs.
    std::int64_t mNumSkipped = 0;
    /// Cells where `LeastSquares` hit a degenerate neighbourhood and fell back
    /// to Green-Gauss. Always 0 for `GreenGauss`.
    std::int64_t mNumFallback = 0;
};

/**
 * @brief Differentiates a `point_data` field over the mesh.
 * @param rMesh the source mesh (unmodified).
 * @param rOptions which array, operator, method, location and output name.
 * @return the mesh with the produced array, plus the skip/fallback counters.
 * @throws std::invalid_argument on an empty or unknown array name, on a
 *         `cell_data` name (the message names the fix), on a row count that
 *         disagrees with the mesh, on an out-of-range or misapplied
 *         `mComponent`, on a component count `Divergence`/`Curl` cannot use, on
 *         `DataLocation::Field`, on an output-name collision when `mOverwrite`
 *         is false, or on a mesh of topological dimension 0.
 */
MESHIOPLUSPLUS_API GradientResult gradient(const Mesh& rMesh, const GradientOptions& rOptions);

/**
 * @brief Parses an operator name.
 * @param rName one of `gradient`/`grad`, `divergence`/`div`, or `curl`/`rot`.
 * @return the matching `GradientOperator`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API GradientOperator gradient_operator_from_name(const std::string& rName);

/**
 * @brief Parses a method name.
 * @param rName one of `green-gauss`/`green_gauss`/`gg` (or empty), or
 *        `least-squares`/`least_squares`/`lsq`.
 * @return the matching `GradientMethod`.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API GradientMethod gradient_method_from_name(const std::string& rName);

}  // namespace meshioplusplus
