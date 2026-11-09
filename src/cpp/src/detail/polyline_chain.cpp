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
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "polyline_chain.hpp"

namespace meshioplusplus {
namespace detail {

std::vector<PolylineChain> polyline_chains(
    const std::vector<std::array<std::int64_t, 2>>& rSegments, std::size_t NumPoints) {
    // Incidence lists, each sorted by (neighbour, segment).
    std::vector<std::vector<std::array<std::int64_t, 2>>> adj(NumPoints);
    for (std::size_t s = 0; s < rSegments.size(); ++s) {
        const std::int64_t a = rSegments[s][0], b = rSegments[s][1];
        if (a == b || a < 0 || b < 0 || static_cast<std::size_t>(a) >= NumPoints ||
            static_cast<std::size_t>(b) >= NumPoints)
            continue;
        adj[static_cast<std::size_t>(a)].push_back({b, static_cast<std::int64_t>(s)});
        adj[static_cast<std::size_t>(b)].push_back({a, static_cast<std::int64_t>(s)});
    }
    for (auto& r_list : adj)
        std::sort(r_list.begin(), r_list.end());
    std::vector<std::uint8_t> used(rSegments.size(), 0);
    const auto next_edge = [&](std::int64_t P) -> const std::array<std::int64_t, 2>* {
        for (const auto& r_e : adj[static_cast<std::size_t>(P)])
            if (used[static_cast<std::size_t>(r_e[1])] == 0)
                return &r_e;
        return nullptr;
    };
    const auto walk = [&](std::int64_t Start, bool Loop) {
        PolylineChain chain;
        chain.mClosed = Loop;
        chain.mPoints.push_back(Start);
        std::int64_t p = Start;
        while (true) {
            const std::array<std::int64_t, 2>* p_e = next_edge(p);
            if (p_e == nullptr)
                break;
            const std::int64_t q = (*p_e)[0], s = (*p_e)[1];
            used[static_cast<std::size_t>(s)] = 1;
            chain.mSegments.push_back(s);
            chain.mReversed.push_back(rSegments[static_cast<std::size_t>(s)][0] == p ? 0 : 1);
            chain.mPoints.push_back(q);
            p = q;
            if (Loop ? p == Start : adj[static_cast<std::size_t>(p)].size() != 2)
                break;
        }
        return chain;
    };
    std::vector<PolylineChain> out;
    for (std::size_t p = 0; p < NumPoints; ++p)
        if (adj[p].size() != 2 && !adj[p].empty())
            while (next_edge(static_cast<std::int64_t>(p)) != nullptr)
                out.push_back(walk(static_cast<std::int64_t>(p), false));
    for (std::size_t p = 0; p < NumPoints; ++p)
        while (next_edge(static_cast<std::int64_t>(p)) != nullptr) {
            PolylineChain chain = walk(static_cast<std::int64_t>(p), true);
            chain.mClosed = chain.mPoints.back() == chain.mPoints.front();
            out.push_back(std::move(chain));
        }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
