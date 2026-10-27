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
#include <cstddef>
#include <limits>
#include <set>
#include <unordered_map>
#include <vector>

// Project includes
#include "curved_quad.hpp"
#include "curved_untangle.hpp"

namespace meshioplusplus::detail {

namespace {

constexpr double kCuPi = 3.14159265358979323846;

double cu_dist(const CurvedVec3& rA, const CurvedVec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

}  // namespace

CurvedUntangleReport curved_untangle(
    int Order, const std::vector<CurvedUntangleElement>& rElems, double* pX,
    const std::vector<std::int64_t>& rFree, std::vector<std::array<double, 2>>& rUv,
    const std::function<CurvedVec3(const std::array<double, 2>&)>& rLift, int MaxSweeps) {
    return curved_untangle(
        Order, rElems, pX, rFree, rUv,
        [&](std::size_t, const std::array<double, 2>& rP) { return rLift(rP); }, nullptr,
        MaxSweeps);
}

CurvedUntangleReport curved_untangle(int Order, const std::vector<CurvedUntangleElement>& rElems,
                                     double* pX, const std::vector<std::int64_t>& rFree,
                                     std::vector<std::array<double, 2>>& rUv,
                                     const CurvedUntangleLift& rLift,
                                     const std::vector<std::array<double, 4>>* pBounds,
                                     int MaxSweeps) {
    CurvedUntangleReport out;
    const auto nn = [&](std::size_t E) {
        return rElems[E].mQuad ? static_cast<int>(curved_quad_num_nodes(Order))
                               : curved_triangle_num_nodes(Order);
    };
    const auto node = [&](std::int64_t Id) {
        const auto i = static_cast<std::size_t>(Id);
        return CurvedVec3{pX[3 * i], pX[3 * i + 1], pX[3 * i + 2]};
    };
    const auto nodes_of = [&](std::size_t E) {
        std::array<CurvedVec3, 16> x{};
        for (int i = 0; i < nn(E); ++i)
            x[static_cast<std::size_t>(i)] = node(rElems[E].mNodes[static_cast<std::size_t>(i)]);
        return x;
    };
    const auto valid = [&](std::size_t E) {
        const auto x = nodes_of(E);
        return rElems[E].mQuad ? curved_quad_check(Order, x.data(), rElems[E].mNormal).mValid
                               : curved_triangle_check(Order, x.data(), rElems[E].mNormal).mValid;
    };
    const auto bound = [&](std::size_t E) {
        const auto x = nodes_of(E);
        return rElems[E].mQuad ? curved_quad_scaled_bound(Order, x.data(), rElems[E].mNormal)
                               : curved_triangle_scaled_bound(Order, x.data(), rElems[E].mNormal);
    };
    // The free nodes by id, and the elements around each.
    std::unordered_map<std::int64_t, std::size_t> free_of;
    for (std::size_t i = 0; i < rFree.size(); ++i)
        free_of.emplace(rFree[i], i);
    std::vector<std::vector<std::size_t>> around(rFree.size());
    for (std::size_t e = 0; e < rElems.size(); ++e)
        for (int i = 0; i < nn(e); ++i) {
            const auto it = free_of.find(rElems[e].mNodes[static_cast<std::size_t>(i)]);
            if (it != free_of.end())
                around[it->second].push_back(e);
        }
    std::vector<char> bad(rElems.size(), 0);
    for (std::size_t e = 0; e < rElems.size(); ++e)
        if (!valid(e)) {
            bad[e] = 1;
            ++out.mInvalidBefore;
        }
    out.mInvalidAfter = out.mInvalidBefore;
    if (out.mInvalidBefore == 0 || rFree.empty())
        return out;
    // A node's step in its parameters: a quarter of the shortest corner edge
    // of its elements, over the surface's speed along each parameter.
    std::vector<std::array<double, 2>> step(rFree.size(), {0.0, 0.0});
    const auto init_step = [&](std::size_t F) {
        double h = std::numeric_limits<double>::infinity();
        for (const std::size_t e : around[F]) {
            const std::size_t corners = rElems[e].mQuad ? 4 : 3;
            for (std::size_t k = 0; k < corners; ++k)
                h = std::min(h, cu_dist(node(rElems[e].mNodes[k]),
                                        node(rElems[e].mNodes[(k + 1) % corners])));
        }
        const auto& uv = rUv[F];
        for (std::size_t d = 0; d < 2; ++d) {
            const double du = 1e-6 * (1.0 + std::fabs(uv[d]));
            std::array<double, 2> a = uv, b = uv;
            a[d] -= du;
            b[d] += du;
            const double speed = cu_dist(rLift(F, a), rLift(F, b)) / (2.0 * du);
            step[F][d] = speed > 0.0 && std::isfinite(h) ? 0.25 * h / speed : 0.0;
        }
    };
    std::set<std::size_t> active;
    std::vector<char> stepped(rFree.size(), 0);
    const auto activate = [&](std::size_t E) {
        for (int i = 0; i < nn(E); ++i) {
            const auto it = free_of.find(rElems[E].mNodes[static_cast<std::size_t>(i)]);
            if (it == free_of.end())
                continue;
            if (!stepped[it->second]) {
                init_step(it->second);
                stepped[it->second] = 1;
            }
            active.insert(it->second);
        }
    };
    for (std::size_t e = 0; e < rElems.size(); ++e)
        if (bad[e])
            activate(e);
    const auto star_bound = [&](std::size_t F) {
        double b = std::numeric_limits<double>::infinity();
        for (const std::size_t e : around[F])
            b = std::min(b, bound(e));
        return b;
    };
    for (int sweep = 0; sweep < MaxSweeps && !active.empty(); ++sweep) {
        ++out.mSweeps;
        bool any_step = false;
        std::set<std::size_t> touched;
        for (const std::size_t f : active) {
            if (!(step[f][0] > 0.0 || step[f][1] > 0.0))
                continue;
            any_step = true;
            const auto id = static_cast<std::size_t>(rFree[f]);
            const std::array<double, 2> uv0 = rUv[f];
            const CurvedVec3 x0 = node(rFree[f]);
            const double b0 = star_bound(f);
            double best = b0;
            std::array<double, 2> best_uv = uv0;
            CurvedVec3 best_x = x0;
            // The elements around it proved valid now must stay so: a move
            // that raises the worst bound cannot push a neighbour over.
            std::vector<char> was_valid(around[f].size(), 0);
            for (std::size_t i = 0; i < around[f].size(); ++i)
                was_valid[i] = bad[around[f][i]] ? 0 : 1;
            for (int k = 0; k < 8; ++k) {
                const double a = kCuPi * double(k) / 4.0;
                std::array<double, 2> uv{uv0[0] + step[f][0] * std::cos(a),
                                         uv0[1] + step[f][1] * std::sin(a)};
                if (pBounds != nullptr) {
                    const auto& r_b = (*pBounds)[f];
                    uv = {std::clamp(uv[0], r_b[0], r_b[1]), std::clamp(uv[1], r_b[2], r_b[3])};
                    if (uv == uv0)
                        continue;
                }
                const CurvedVec3 x = rLift(f, uv);
                for (std::size_t c = 0; c < 3; ++c)
                    pX[3 * id + c] = x[c];
                const double b = star_bound(f);
                if (!(b > best + 1e-12))
                    continue;
                bool keeps = true;
                for (std::size_t i = 0; i < around[f].size() && keeps; ++i)
                    if (was_valid[i] && !(bound(around[f][i]) > 0.0))
                        keeps = valid(around[f][i]);
                if (keeps) {
                    best = b;
                    best_uv = uv;
                    best_x = x;
                }
            }
            for (std::size_t c = 0; c < 3; ++c)
                pX[3 * id + c] = best_x[c];
            if (best > b0 + 1e-12) {
                rUv[f] = best_uv;
                ++out.mMoves;
                for (const std::size_t e : around[f])
                    touched.insert(e);
            } else {
                step[f][0] *= 0.5;
                step[f][1] *= 0.5;
                if (step[f][0] < 1e-12 * (1.0 + std::fabs(uv0[0])) &&
                    step[f][1] < 1e-12 * (1.0 + std::fabs(uv0[1])))
                    step[f] = {0.0, 0.0};
            }
        }
        // The elements a move changed are judged again: a newly invalid one
        // brings its own free nodes in.
        for (const std::size_t e : touched) {
            const bool now_bad = !valid(e);
            if (now_bad && !bad[e])
                activate(e);
            bad[e] = now_bad ? 1 : 0;
        }
        bool left = false;
        for (std::size_t e = 0; e < rElems.size() && !left; ++e)
            left = bad[e] != 0;
        if (!left || !any_step)
            break;
    }
    out.mInvalidAfter = 0;
    for (std::size_t e = 0; e < rElems.size(); ++e)
        out.mInvalidAfter += valid(e) ? 0 : 1;
    return out;
}

}  // namespace meshioplusplus::detail
