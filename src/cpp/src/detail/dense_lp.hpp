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
 * @file detail/dense_lp.hpp
 * @brief A small dense linear-programme solver (roadmap 7.1.4): the two-phase
 * simplex method on a dense tableau, Bland's rule throughout, so it neither
 * cycles nor depends on anything but its input.
 *
 *  - `dense_lp_equality` solves `min c.y` subject to `A y = b`, `y >= 0`:
 *    phase one minimises the sum of one artificial variable per row (each
 *    row first scaled so its right-hand side is not negative), phase two the
 *    objective, the artificials kept out of the basis. It returns the
 *    optimum, its value and the duals of the equality rows (from the
 *    artificial columns' reduced costs).
 *  - `dense_lp_inequality` solves `min c.x` subject to `A x <= b` with `x`
 *    free through its dual, `max -b.y` subject to `A^T y = -c`, `y >= 0`:
 *    the tableau has one row per variable of `x` and one column per
 *    constraint, which suits the fits it serves (a handful of variables, a
 *    constraint per measured point). `x` is the dual's duals.
 *
 * A pivot element must exceed `1e-12` times its column's largest entry, and
 * a reduced cost below `-1e-12` times the objective's scale enters; on a
 * tie the lowest index wins. Sequential and deterministic.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

enum class DenseLpStatus : int {
    Optimal = 0,
    Infeasible = 1,
    Unbounded = 2,
    IterationLimit = 3,
};

struct DenseLpResult {
    DenseLpStatus mStatus = DenseLpStatus::Infeasible;
    std::vector<double> mX;      ///< the optimum
    std::vector<double> mDuals;  ///< `dense_lp_equality`: one per equality row
    double mObjective = 0.0;
    std::int64_t mPivots = 0;
};

/**
 * @brief `min c.y` subject to `A y = b`, `y >= 0`.
 * @param rA row-major, `b.size()` rows of `c.size()` entries.
 */
DenseLpResult dense_lp_equality(const std::vector<double>& rC, const std::vector<double>& rA,
                                const std::vector<double>& rB, std::int64_t MaxPivots = 1000000);

/**
 * @brief `min c.x` subject to `A x <= b`, `x` free.
 * @param rA row-major, `b.size()` rows of `c.size()` entries.
 */
DenseLpResult dense_lp_inequality(const std::vector<double>& rC, const std::vector<double>& rA,
                                  const std::vector<double>& rB, std::int64_t MaxPivots = 1000000);

}  // namespace detail
}  // namespace meshioplusplus
