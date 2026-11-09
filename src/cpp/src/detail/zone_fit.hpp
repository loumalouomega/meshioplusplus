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
 * @file detail/zone_fit.hpp
 * @brief Minimum-zone and envelope fits of a plane, sphere or cylinder
 * (roadmap 7.1.4): the Chebyshev references GD&T form tolerances are defined
 * against (ISO 1101; ISO 12781 flatness, 12181 roundness, 12180
 * cylindricity) and the envelope references of a shaft and a hole.
 *
 *  - **Minimum zone**: the narrowest zone holding every point, between two
 *    parallel planes, two concentric spheres or two coaxial cylinders.
 *  - **Minimum circumscribed**: the smallest sphere or cylinder enclosing
 *    every point (a shaft's envelope).
 *  - **Maximum inscribed**: the largest sphere or cylinder no point lies
 *    inside (a hole's envelope).
 *
 * From a start (the least-squares fit), each round linearises the distances
 * in the surface's parameters (a plane's normal tilt and offset; a sphere's
 * centre; a cylinder's axis point and tilt, in the frame `lsq_frame` gives
 * its axis) and solves the linear programme of the criterion within a trust
 * region by `dense_lp_inequality` (Carr and Ferreira 1995; Shunmugam 1986: a
 * sequence of linear programmes); a step is kept when the criterion's exact
 * value improves, else the region shrinks, until it is below `1e-14` of the
 * points' extent or 200 rounds pass. Deterministic and serial.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

namespace meshioplusplus {
namespace detail {

enum class ZoneCriterion : int {
    MinimumZone = 1,
    MinCircumscribed = 2,
    MaxInscribed = 3,
};

struct ZoneFitResult {
    bool mOk = false;
    std::string mWhy;  ///< why not, when not
    /// The fit table's ten parameters (`fit:patch_params` layout): a plane's
    /// mid-zone point and unit normal; a sphere's centre and radius; a
    /// cylinder's axis point nearest the points' mean, unit axis and radius.
    /// The radius is the zone's middle, or the envelope's.
    std::array<double, 10> mParams{};
    /// The zone's width (minimum zone), or the envelope's radius.
    double mValue = 0.0;
    std::int64_t mRounds = 0;
};

/**
 * @brief Refit a plane (type 1, minimum zone only), sphere (2) or cylinder
 * (3) to @p N points (`pPts`, 3 each) under @p Criterion, from the fit-table
 * parameters @p pStart.
 */
ZoneFitResult zone_fit(std::int64_t Type, const double* pStart, const double* pPts, std::size_t N,
                       ZoneCriterion Criterion);

}  // namespace detail
}  // namespace meshioplusplus
