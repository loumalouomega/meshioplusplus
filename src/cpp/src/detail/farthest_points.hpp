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
 * @file detail/farthest_points.hpp
 * @brief Greedy farthest-point sampling: the kernel under the Python
 * `select_points(method="farthest")` (roadmap 5.10.2).
 *
 * A **core-private** header (the `detail/slot_runs.hpp` precedent), bound to
 * Python alone: the Python layer that calls it stays Python. The arithmetic
 * is numpy's, step for step, so the selection is the twin's bit for bit:
 * squared distances as `(dx * dx + dy * dy) + dz * dz` (the core builds with
 * `-ffp-contract=off`), `np.minimum`'s NaN propagation, and `np.argmax`'s
 * first maximum (a NaN counts as the maximum, the first NaN winning).
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// What `farthest_point_order` returns.
struct FarthestPoints {
    std::vector<std::int64_t> mOrder;  ///< the selected indices, `mOrder[0]` the start
    /// `mOrder.size()` when the budget was met; otherwise the step at which
    /// every remaining distance was zero (the cloud's distinct positions ran
    /// out), and `mOrder` holds that many indices.
    std::size_t mSelected = 0;
};

/**
 * @brief Select `Count` of the `N` points at `pXyz` (row-major `(N, 3)`
 * doubles), starting from `First`, each next one the farthest from those
 * already chosen. `Count` must be in `[1, N]` and `First` in `[0, N)`.
 */
FarthestPoints farthest_point_order(const double* pXyz, std::size_t N, std::size_t Count,
                                    std::size_t First);

}  // namespace detail
}  // namespace meshioplusplus
