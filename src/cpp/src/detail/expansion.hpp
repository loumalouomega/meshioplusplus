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
 * @file detail/expansion.hpp
 * @brief Shewchuk's exact expansion arithmetic, as the predicates of
 * `detail/predicates.hpp` use it, for exact constructions (roadmap 7.8, the
 * surface-boolean spike: points on an edge and a plane as exact homogeneous
 * coordinates).
 *
 * An expansion is a vector of nonoverlapping doubles in increasing
 * magnitude, zeros eliminated; its value is their exact sum, and empty is
 * zero. The functions are defined in `predicates.cpp`, where floating-point
 * contraction is switched off, and are not inline on purpose: whoever calls
 * them gets the paper's arithmetic whatever its own compilation flags.
 *
 * Core-private: nothing installed names it.
 */

// System includes
#include <vector>

namespace meshioplusplus {
namespace detail {

using Expansion = std::vector<double>;

/// `A - B` exactly.
Expansion exp_diff(double A, double B);
/// `A` itself (empty for zero).
Expansion exp_of(double A);
Expansion exp_sum(const Expansion& rE, const Expansion& rF);
Expansion exp_sub(const Expansion& rE, const Expansion& rF);
Expansion exp_scale(const Expansion& rE, double B);
Expansion exp_mul(const Expansion& rE, const Expansion& rF);
/// The value rounded to a double; its sign is the exact sign.
double exp_estimate(const Expansion& rE);
/// The exact sign: -1, 0 or 1.
int exp_sign(const Expansion& rE);

}  // namespace detail
}  // namespace meshioplusplus
