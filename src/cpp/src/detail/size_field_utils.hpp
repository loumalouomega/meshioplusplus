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

#pragma once

/**
 * @file detail/size_field_utils.hpp
 * @brief An isotropic size field on mesh vertices: metric edge lengths, the
 * metric midpoint of an edge, size interpolation and gradation.
 *
 * Sizes vary geometrically along an edge, `h(t) = ha (hb / ha)^t`, so the
 * metric length of an edge is `|ab| (hb - ha) / (ha hb ln(hb / ha))` and
 * `|ab| / ha` when the sizes agree (Dapogny, Dobrzynski and Frey 2014). An
 * edge is unit when its metric length is 1; `adapt_volume` splits above 4/3
 * and collapses below 4/5 (Botsch and Kobbelt 2004).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// Static helpers over an isotropic size field (Kratos-style utility class).
struct SizeFieldUtils {
    /// The metric length of the edge (a, b) with sizes @p Ha and @p Hb.
    static double EdgeLength(const Vec3& rA, const Vec3& rB, double Ha, double Hb);

    /// The metric length of an edge of Euclidean length @p Length.
    static double MetricLength(double Length, double Ha, double Hb);

    /// The fraction `t` of the way from a to b at which the metric length is
    /// half the edge's: 1/2 when the sizes agree.
    static double MidpointParameter(double Ha, double Hb);

    /// The size at fraction @p T of the way from a to b, `ha (hb / ha)^t`.
    static double Interpolate(double Ha, double Hb, double T);

    /**
     * @brief Grade @p rSize in place over the undirected @p rEdges, so that
     * `h(v) <= h(u) + (Hgrad - 1) |uv|` along every edge: a Dijkstra sweep from
     * the smallest sizes, ties broken by vertex id (deterministic).
     * @return how many sizes were lowered.
     */
    static std::int64_t Grade(const std::vector<double>& rXyz,
                              const std::vector<std::array<std::int64_t, 2>>& rEdges, double Hgrad,
                              std::vector<double>& rSize);
};

}  // namespace detail
}  // namespace meshioplusplus
