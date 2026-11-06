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
 * @file detail/bspline_ndu.hpp
 * @brief The Cox-de Boor triangle of The NURBS Book's A2.3 (the `ndu` table of
 * basis functions and knot differences), shared by `bspline_lsq.cpp` and
 * `iga_kernels.hpp`. Each caller keeps its own scratch storage and its own
 * derivative pass; only this loop, with its arithmetic order, is shared.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstddef>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// Fill the `(Degree + 1) x (Degree + 1)` table `rNdu` at knot span `Span`
/// and parameter `U` (row-major, entry `(R, C)` at `R * (Degree + 1) + C`).
/// `rLeft` and `rRight` take `Degree + 1` values of scratch; any indexable
/// storage works.
template <class TNdu, class TLeft, class TRight>
inline void bspline_ndu(std::size_t Span, double U, std::size_t Degree,
                        const std::vector<double>& rKnots, TNdu& rNdu, TLeft& rLeft,
                        TRight& rRight) {
    const std::size_t p = Degree;
    const auto at = [p](std::size_t R, std::size_t C) { return R * (p + 1) + C; };
    rNdu[at(0, 0)] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        rLeft[j] = U - rKnots[Span + 1 - j];
        rRight[j] = rKnots[Span + j] - U;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            rNdu[at(j, r)] = rRight[r + 1] + rLeft[j - r];
            const double temp = rNdu[at(r, j - 1)] / rNdu[at(j, r)];
            rNdu[at(r, j)] = saved + rRight[r + 1] * temp;
            saved = rLeft[j - r] * temp;
        }
        rNdu[at(j, j)] = saved;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
