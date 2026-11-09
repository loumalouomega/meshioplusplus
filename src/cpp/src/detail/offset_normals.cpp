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
// System includes
#include <algorithm>
#include <cmath>
#include <limits>

// Project includes
#include "math_utils.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"
#include "offset_normals.hpp"

namespace meshioplusplus::detail {

using OnV3 = std::array<double, 3>;

OffsetNormals offset_normals(const std::vector<double>& rX,
                             const std::vector<std::array<std::int64_t, 4>>& rFaces, double Sign,
                             int SmoothIterations) {
    const std::size_t np = rX.size() / 3;
    const std::size_t nf = rFaces.size();
    const auto& faces = rFaces;
    const auto point = [&](std::int64_t P) {
        const auto i = static_cast<std::size_t>(P);
        return OnV3{rX[3 * i], rX[3 * i + 1], rX[3 * i + 2]};
    };
    const auto corners = [&](std::size_t F) { return faces[F][3] >= 0 ? 4 : 3; };
    // Unit face normals along the growth, and each point's faces with its
    // corner angle there.
    const double sign = Sign;
    std::vector<OnV3> fn(nf);
    std::vector<std::vector<std::pair<std::size_t, double>>> around(np);
    std::vector<std::vector<std::int64_t>> nbrs(np);
    for (std::size_t f = 0; f < nf; ++f) {
        const int k = corners(f);
        OnV3 n{0.0, 0.0, 0.0};
        for (int i = 0; i < k; ++i) {
            const OnV3 a = point(faces[f][static_cast<std::size_t>(i)]);
            const OnV3 b = point(faces[f][static_cast<std::size_t>((i + 1) % k)]);
            n = {n[0] + (a[1] - b[1]) * (a[2] + b[2]), n[1] + (a[2] - b[2]) * (a[0] + b[0]),
                 n[2] + (a[0] - b[0]) * (a[1] + b[1])};
        }
        n = detail::MathUtils::UnitByDivision(n);
        fn[f] = {sign * n[0], sign * n[1], sign * n[2]};
        for (int i = 0; i < k; ++i) {
            const std::int64_t p = faces[f][static_cast<std::size_t>(i)];
            const std::int64_t q = faces[f][static_cast<std::size_t>((i + 1) % k)];
            const std::int64_t r = faces[f][static_cast<std::size_t>((i + k - 1) % k)];
            const OnV3 u = detail::MathUtils::UnitByDivision(detail::vec3_sub(point(q), point(p)));
            const OnV3 v = detail::MathUtils::UnitByDivision(detail::vec3_sub(point(r), point(p)));
            const double angle = std::acos(std::clamp(detail::vec3_dot(u, v), -1.0, 1.0));
            around[static_cast<std::size_t>(p)].push_back({f, angle});
            nbrs[static_cast<std::size_t>(p)].push_back(q);
            nbrs[static_cast<std::size_t>(p)].push_back(r);
        }
    }
    for (auto& r_n : nbrs) {
        std::sort(r_n.begin(), r_n.end());
        r_n.erase(std::unique(r_n.begin(), r_n.end()), r_n.end());
    }
    // Each point's most normal normal (Aubry and Löhner 2009): from the
    // angle-weighted mean, steps toward the face that sees it least.
    std::vector<OnV3> dir(np, OnV3{0.0, 0.0, 0.0});
    std::vector<double> stretch(np, 0.0), seen(np, 0.0);
    std::vector<char> used(np, 0), not_visible(np, 0);
    parallel_for(
        np,
        [&](std::size_t P) {
            const auto& r_a = around[P];
            if (r_a.empty())
                return;
            used[P] = 1;
            OnV3 m{0.0, 0.0, 0.0};
            for (const auto& [f, w] : r_a)
                m = {m[0] + w * fn[f][0], m[1] + w * fn[f][1], m[2] + w * fn[f][2]};
            OnV3 n = detail::MathUtils::UnitByDivision(m);
            const auto least = [&](const OnV3& rN, std::size_t& rWho) {
                double lo = std::numeric_limits<double>::infinity();
                for (const auto& [f, w] : r_a) {
                    (void)w;
                    const double d = detail::vec3_dot(rN, fn[f]);
                    if (d < lo) {
                        lo = d;
                        rWho = f;
                    }
                }
                return lo;
            };
            std::size_t who = 0;
            OnV3 best = n;
            double best_lo = least(n, who);
            for (int it = 0; it < 100 && r_a.size() > 1; ++it) {
                const double eta = 1.0 / double(it + 2);
                n = detail::MathUtils::UnitByDivision(
                    {n[0] + eta * fn[who][0], n[1] + eta * fn[who][1], n[2] + eta * fn[who][2]});
                const double lo = least(n, who);
                if (lo > best_lo) {
                    best_lo = lo;
                    best = n;
                }
            }
            dir[P] = best;
            if (!(best_lo > kOffsetVisible)) {
                not_visible[P] = 1;
                return;
            }
            stretch[P] = 1.0 / std::max(best_lo, 1.0 / kOffsetMaxStretch);
            seen[P] = best_lo;
        },
        64);
    // Smoothed directions (6.9.1): each moved toward its neighbours' mean,
    // kept while every face around it still sees it at more than half the
    // least cosine its most normal normal had (Jacobi passes, so the result
    // does not depend on the order).
    const auto least_cos = [&](std::size_t P, const OnV3& rN) {
        double lo = std::numeric_limits<double>::infinity();
        for (const auto& [f, w] : around[P]) {
            (void)w;
            lo = std::min(lo, detail::vec3_dot(rN, fn[f]));
        }
        return lo;
    };
    for (int it = 0; it < SmoothIterations; ++it) {
        std::vector<OnV3> next = dir;
        parallel_for(
            np,
            [&](std::size_t P) {
                if (!(stretch[P] > 0.0) || nbrs[P].empty())
                    return;
                OnV3 m = dir[P];
                for (const std::int64_t q : nbrs[P]) {
                    const auto& r_q = dir[static_cast<std::size_t>(q)];
                    if (stretch[static_cast<std::size_t>(q)] > 0.0)
                        m = {m[0] + r_q[0], m[1] + r_q[1], m[2] + r_q[2]};
                }
                const OnV3 n = detail::MathUtils::UnitByDivision(m);
                const double lo = least_cos(P, n);
                if (lo > 0.5 * seen[P] && lo > kOffsetVisible)
                    next[P] = n;
            },
            64);
        dir = std::move(next);
    }
    if (SmoothIterations > 0)
        for (std::size_t p = 0; p < np; ++p)
            if (stretch[p] > 0.0)
                stretch[p] = 1.0 / std::max(least_cos(p, dir[p]), 1.0 / kOffsetMaxStretch);
    OffsetNormals out;
    out.mFaceNormals = std::move(fn);
    out.mAround = std::move(around);
    out.mNeighbours = std::move(nbrs);
    out.mDirections = std::move(dir);
    out.mStretch = std::move(stretch);
    out.mSeen = std::move(seen);
    out.mUsed = std::move(used);
    out.mNotVisible = std::move(not_visible);
    return out;
}

}  // namespace meshioplusplus::detail
