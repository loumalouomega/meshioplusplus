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
 * @file operations/hessian.hpp
 * @brief The second derivative (Hessian matrix) of a scalar `point_data`
 * field, for curvature-based adaptive refinement -- `gradient`'s companion
 * one order further: `gradient` differentiates a field once, this
 * differentiates it twice.
 *
 * ### Composition, not a new numerical kernel
 *
 * `hessian` is built entirely out of two calls to `gradient()`, exactly the
 * precedent `estimate_error` already sets: (1) `gradient(mArrayName,
 * location=Point)` gives the field's gradient as a genuine `point_data`
 * array (`(n, 3)`); (2) `gradient()` again on THAT array, with the default
 * `Gradient` operator, differentiates a 3-component field and so produces
 * `(n, 9)` -- the flattened row-major 3x3 Hessian, `H[i][j] = d^2f/dxi dxj`
 * at index `i*3+j`. The first call MUST use `Point` location: the second
 * call's `point_data`-only validation would otherwise reject a `cell_data`
 * intermediate. No new differentiation math is written here at all --
 * `mMethod` simply forwards to both internal `gradient()` calls, so
 * `LeastSquares` is available for free, exactly as it already is for
 * `gradient` itself.
 *
 * ### Exactness -- stated honestly, not oversold
 *
 * Green-Gauss is exact for a linear field on any cell (see `gradient.hpp`).
 * A field that is at most LINEAR therefore has an EXACTLY ZERO Hessian
 * everywhere: its gradient is a constant, and Green-Gauss of an exactly
 * zero/constant field is trivially exact. This is the one rigorous,
 * mesh-shape-independent guarantee, and it is what `test_hessian.cpp` pins
 * -- the same reasoning `test_gradient.cpp`'s own
 * `CurlOfAGradientAndDivergenceOfACurlVanish` already relies on, for the
 * identical reason.
 *
 * For a genuinely QUADRATIC field, each Green-Gauss pass alone would be
 * exact (the gradient of a quadratic field is linear), but the mandatory
 * intermediate step between the two passes --
 * `cell_data_to_point_data(weight=Uniform)`, a plain arithmetic mean of
 * incident-cell values -- is only exact for a field that is constant over
 * the averaged neighbourhood. For a genuinely-varying linear field (i.e.
 * the gradient of a true quadratic), averaging several cells'
 * individually-exact-but-different local values reproduces the true nodal
 * value only when those cells are symmetric about the shared node.
 * Measured, not merely argued: on a regular axis-aligned hex grid, every
 * INTERIOR node's 8-cell neighbourhood is exactly that symmetric, and the
 * composed Hessian comes back exact to machine precision there; the same
 * mesh's own BOUNDARY cells (a one-sided, asymmetric neighbourhood) show
 * real, bounded error (`test_hessian.cpp`'s
 * `QuadraticFieldOnAStructuredGridIsExactAwayFromTheBoundary` and
 * `...HasBoundedBoundaryError` pin both halves of this on the identical
 * mesh). So the composed Hessian is EXACT away from a mesh's own boundary
 * on a structured/symmetric mesh, and a good, standard, but genuinely
 * APPROXIMATE curvature estimate on an irregular one -- the same honesty
 * `gradient.hpp` itself uses ("first-order on a warped quad in 3D") rather
 * than a blanket exactness claim. A full quadratic-least-squares Hessian
 * recovery would restore mesh-shape independence but is a substantially
 * larger undertaking and out of scope here.
 *
 * ### Scope: scalar fields only
 *
 * The input must have exactly one component. A vector field's Hessian is a
 * separate quantity per component (`(n, 3*9)`, one 3x3 matrix per vector
 * component) -- a real but different need, refused by name rather than
 * guessed at: call `hessian` once per component.
 *
 * ### Curvature-driven refinement needs no new marking step
 *
 * `data_calc`'s `norm(...)` is a plain sum-of-squares-then-sqrt over
 * however many components its argument has, so `norm(<array>:hessian)` on
 * the 9-component output is exactly its Frobenius norm -- a scalar
 * curvature indicator with zero new code, ready for `refine --where`. See
 * `doc/hessian.md`'s worked composition (`hessian` -> `data_calc norm` ->
 * `refine --where`), which is deliberately NOT a bespoke marking pass
 * (unlike `estimate_error`'s ZZ-specific one): the composable pieces
 * already exist and are more general.
 *
 * ### Contracts
 *
 * - Input must be `point_data` with exactly one component.
 * - Private working arrays are never returned: like `estimate_error`, the
 *   actual result mesh is built from a fresh `detail::clone_mesh(rMesh)`,
 *   so a name collision with the caller's own data cannot lose anything.
 * - Geometry, connectivity, regions, property sets and every existing data
 *   array pass through bit-identically.
 *
 * Output is byte-identical across the three mesh backends and thread
 * counts. Across the C++/numpy boundary it matches the numpy twin to
 * within a numeric tolerance rather than bit-for-bit -- inherited from the
 * `Uniform`-weighted averaging step the composition relies on, the same
 * already-accepted precedent `estimate_error`'s own `Measure`-weighted
 * recovery step has.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format -- it
 * is not in the format registry.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/gradient.hpp"

namespace meshioplusplus {

/// Default suffix of the produced array.
inline constexpr const char* kHessianSuffix = ":hessian";

/// Options for `hessian`.
struct HessianOptions {
    /// Name of the scalar `point_data` array to differentiate twice. Required.
    std::string mArrayName;
    /// How each of the two internal `gradient()` passes reconstructs its
    /// derivative. Forwarded unchanged to both.
    GradientMethod mMethod = GradientMethod::GreenGauss;
    /// Where the result is stored: `Point` or `Cell`. `Field` is an error.
    DataLocation mLocation = DataLocation::Cell;
    /// Output array name; empty selects `"<array>:hessian"`.
    std::string mOutputName;
    /// When false, an output name that already exists is an error rather
    /// than being overwritten.
    bool mOverwrite = false;
};

/// The mesh with the produced `(n, 9)` array, plus what could not be computed.
struct HessianResult {
    /// The input mesh with the Hessian array added; everything else unchanged.
    Mesh mMesh;
    /// Cells where the Hessian could not be computed (the second internal
    /// `gradient()` pass's own count -- structural skip reasons are
    /// identical between the two passes, since both run over the same mesh
    /// topology). The Hessian reads NaN there.
    std::int64_t mNumSkipped = 0;
    /// Cells where either internal `gradient()` pass fell back from
    /// least-squares to Green-Gauss (summed over both passes; always 0 for
    /// `mMethod == GreenGauss`, which has no fallback concept).
    std::int64_t mNumFallback = 0;
};

/**
 * @brief The Hessian (second derivative) of a scalar `point_data` field.
 * @param rMesh the source mesh (unmodified).
 * @param rOptions which array, method, location and output name.
 * @return the mesh with the produced `(n, 9)` array and the skip/fallback
 *         counters.
 * @throws std::invalid_argument on an empty or unknown array name, on a
 *         `cell_data` name (message names the fix), on an array with more
 *         than one component (message names the per-component workaround),
 *         on a mesh with no cell of topological dimension 1 or more, or on
 *         an output-name collision when `mOverwrite` is false.
 */
MESHIOPLUSPLUS_API HessianResult hessian(const Mesh& rMesh, const HessianOptions& rOptions);

}  // namespace meshioplusplus
