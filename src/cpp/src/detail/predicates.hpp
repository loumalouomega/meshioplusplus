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
 * @file detail/predicates.hpp
 * @brief Robust geometric predicates: the exact sign of `orient2d`,
 * `incircle`, `orient3d` and `insphere` for any finite double input, used by
 * the constrained Delaunay triangulation of `detail/cdt2d` (`mesh_cad_surface`).
 *
 * The method is J. R. Shewchuk's, "Adaptive Precision Floating-Point
 * Arithmetic and Fast Robust Geometric Predicates", Discrete & Computational
 * Geometry 18(3), 1997, whose reference code is in the public domain:
 *
 *  - **A filter.** Each determinant is first evaluated in plain floating
 *    point together with its permanent (the same expression on absolute
 *    values); when `|det|` exceeds the paper's a-priori bound (stage A:
 *    `(3 + 16ε)ε`, `(10 + 96ε)ε`, `(7 + 56ε)ε` and `(16 + 224ε)ε` times the
 *    permanent) its sign is certain and is returned at once.
 *  - **Exact evaluation.** Otherwise the coordinate differences are formed
 *    exactly as two-term expansions (`two_diff`) and the determinant is
 *    evaluated exactly in expansion arithmetic (`two_sum`, `two_product` by
 *    Dekker's split, zero-eliminating expansion sums and scalings, sections
 *    2.4-2.8 of the paper). Its largest component has the sign of the exact
 *    value. The paper's intermediate adaptive stages (B, C) are folded into
 *    this exact stage: the answer is the same, at some cost on inputs that
 *    are near-degenerate.
 *
 * This is a re-implementation from the paper, not a copy of `predicates.c`.
 * It needs IEEE-754 binary64 with round-to-nearest-even (checked at compile
 * time). Floating-point contraction is switched off inside `predicates.cpp`
 * (GCC's `optimize` pragma, clang's `fp contract(off)`; the core is built
 * with `-ffp-contract=off` anyway), so a single-header consumer compiled with
 * `-mfma` gets the arithmetic the paper analyses. The design does not depend
 * on it: the exact stage multiplies only half-width operands, whose products
 * are exact fused or not, and a fused evaluation of a filter has fewer
 * roundings than the bound allows for; `test_predicates.cpp` passes built
 * with `-mfma -ffp-contract=fast`, which clang applies over the pragma.
 * Overflow and underflow of the intermediate products are outside the
 * guarantee, as in the paper.
 *
 * Each predicate returns a value whose **sign** is exact and whose magnitude
 * approximates the determinant:
 *
 *  - `orient2d(a, b, c)` > 0 when `a, b, c` turn counter-clockwise;
 *  - `incircle(a, b, c, d)` > 0 when `d` lies inside the circle through the
 *    counter-clockwise `a, b, c`;
 *  - `orient3d(a, b, c, d)` > 0 when `d` lies below the plane of `a, b, c`
 *    (those counter-clockwise seen from above), Shewchuk's convention;
 *  - `insphere(a, b, c, d, e)` > 0 when `e` lies inside the sphere through
 *    `a, b, c, d` with `orient3d(a, b, c, d)` > 0;
 *  - `diametral(a, b, p)` < 0 when `p` lies strictly inside the circle whose
 *    diameter is `ab` (the angle `apb` is obtuse): the sign of
 *    `(a - p) . (b - p)`, Ruppert's encroachment test;
 *  - `lens_angle(a, b, p)` > 0 when the angle `apb` is below 60 degrees or
 *    above 120: the sign of `3 ((a - p) . (b - p))^2 - ((a - p) x (b - p))^2`.
 *    With `diametral(a, b, p)` > 0 it says the corner `apb` is sharper than
 *    60 degrees, Ruppert's small input angle; with `diametral` < 0, that
 *    `p` lies strictly inside the diametral lens of `ab`, where `ab`
 *    subtends more than 120 degrees (Chew 1989; Shewchuk 2002).
 *
 * The last two are not in Shewchuk's paper: they are built the same way
 * (a filter with an a-priori bound on the permanent, then the exact
 * expansion), for the encroachment and small-angle decisions of Delaunay
 * refinement (`Cdt2d::Refine`), so those are exact too.
 *
 * A **core-private** header: nothing installed names it. Pure functions with
 * no state, so safe from any thread.
 */

namespace meshioplusplus {
namespace detail {

double orient2d(const double* pA, const double* pB, const double* pC);
double incircle(const double* pA, const double* pB, const double* pC, const double* pD);
double orient3d(const double* pA, const double* pB, const double* pC, const double* pD);
double insphere(const double* pA, const double* pB, const double* pC, const double* pD,
                const double* pE);
double diametral(const double* pA, const double* pB, const double* pP);
double lens_angle(const double* pA, const double* pB, const double* pP);

/// The exact-stage arithmetic, exposed for its own tests: `x + y == a * b`
/// exactly, `x` the rounded product.
void predicates_two_product(double A, double B, double& rX, double& rY);

}  // namespace detail
}  // namespace meshioplusplus
