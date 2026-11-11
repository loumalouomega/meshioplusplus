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

// System includes
#include <cmath>
#include <cstddef>
#include <functional>
#include <queue>
#include <utility>

// Project includes
#include "size_field_utils.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Relative size difference below which an edge is measured as uniform: the
// closed forms divide by ln(hb / ha), which loses its digits as hb -> ha.
constexpr double kSfuUniform = 1e-6;

}  // namespace

double SizeFieldUtils::MetricLength(double Length, double Ha, double Hb) {
    if (!(Ha > 0.0) || !(Hb > 0.0))
        return Length;  // no size: Euclidean
    const double d = Hb - Ha;
    if (std::abs(d) <= kSfuUniform * Ha)
        return Length / (0.5 * (Ha + Hb));
    return Length * d / (Ha * Hb * std::log(Hb / Ha));
}

double SizeFieldUtils::EdgeLength(const Vec3& rA, const Vec3& rB, double Ha, double Hb) {
    return MetricLength(vec3_norm(vec3_sub(rB, rA)), Ha, Hb);
}

double SizeFieldUtils::MidpointParameter(double Ha, double Hb) {
    if (!(Ha > 0.0) || !(Hb > 0.0) || std::abs(Hb - Ha) <= kSfuUniform * Ha)
        return 0.5;
    // int_0^t ds / (ha r^s) = (1 - r^-t) / (ha ln r); half of the total at
    // r^-t = (1 + 1/r) / 2.
    const double r = Hb / Ha;
    return -std::log(0.5 * (1.0 + 1.0 / r)) / std::log(r);
}

double SizeFieldUtils::Interpolate(double Ha, double Hb, double T) {
    if (!(Ha > 0.0) || !(Hb > 0.0))
        return Ha + (Hb - Ha) * T;
    return Ha * std::pow(Hb / Ha, T);
}

std::int64_t SizeFieldUtils::Grade(const std::vector<double>& rXyz,
                                   const std::vector<std::array<std::int64_t, 2>>& rEdges,
                                   double Hgrad, std::vector<double>& rSize) {
    const std::size_t n = rSize.size();
    const double slope = Hgrad - 1.0;
    if (!(slope >= 0.0) || n == 0)
        return 0;
    // CSR adjacency, neighbours in edge order (the edges come sorted).
    std::vector<std::size_t> start(n + 1, 0);
    for (const auto& e : rEdges) {
        ++start[static_cast<std::size_t>(e[0]) + 1];
        ++start[static_cast<std::size_t>(e[1]) + 1];
    }
    for (std::size_t i = 0; i < n; ++i)
        start[i + 1] += start[i];
    std::vector<std::int64_t> adj(start[n]);
    std::vector<std::size_t> fill(start.begin(), start.end() - 1);
    for (const auto& e : rEdges) {
        adj[fill[static_cast<std::size_t>(e[0])]++] = e[1];
        adj[fill[static_cast<std::size_t>(e[1])]++] = e[0];
    }
    auto dist = [&](std::int64_t A, std::int64_t B) {
        const std::size_t a = static_cast<std::size_t>(A) * 3, b = static_cast<std::size_t>(B) * 3;
        const double dx = rXyz[a] - rXyz[b], dy = rXyz[a + 1] - rXyz[b + 1],
                     dz = rXyz[a + 2] - rXyz[b + 2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    };
    using Entry = std::pair<double, std::int64_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> heap;
    for (std::size_t v = 0; v < n; ++v)
        if (start[v + 1] > start[v] && rSize[v] > 0.0)
            heap.push({rSize[v], static_cast<std::int64_t>(v)});
    std::vector<std::uint8_t> lowered(n, 0);
    while (!heap.empty()) {
        const auto [h, v] = heap.top();
        heap.pop();
        if (h != rSize[static_cast<std::size_t>(v)])
            continue;  // stale
        for (std::size_t k = start[static_cast<std::size_t>(v)];
             k < start[static_cast<std::size_t>(v) + 1]; ++k) {
            const std::int64_t u = adj[k];
            const double cand = h + slope * dist(v, u);
            if (cand < rSize[static_cast<std::size_t>(u)]) {
                rSize[static_cast<std::size_t>(u)] = cand;
                lowered[static_cast<std::size_t>(u)] = 1;
                heap.push({cand, u});
            }
        }
    }
    std::int64_t count = 0;
    for (const std::uint8_t l : lowered)
        count += l;
    return count;
}

}  // namespace detail
}  // namespace meshioplusplus
