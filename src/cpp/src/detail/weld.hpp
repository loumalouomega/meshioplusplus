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
 * @file detail/weld.hpp
 * @brief Keep-first point welding, shared by `clean` and `merge`.
 *
 * A **core-private** header (the `formats/gid_common.hpp` precedent): no
 * installed header names it.
 *
 * The rule both operations have always used: points are visited in index
 * order; a point within `atol` of an existing *representative* in its own or
 * one of the 26 neighbouring grid cells (cell size `atol`, cells visited in
 * ascending dz -> dy -> dx order, representatives in creation order) joins the
 * first such representative, otherwise it becomes a new one. The result is
 * order-dependent by design (chains A~B, B~C, A!~C), so the decision loop
 * stays serial; what this moves out of it (roadmap §3, "Welding") is
 * everything else: the cell keys and each cell's non-empty neighbour cells are
 * computed in parallel up front (cell ids come from one serial pass through a
 * flat open-addressing table), and the serial loop walks per-cell
 * representative lists by index -- no node-based hash map, no dtype switch --
 * with exactly the old visiting order, so the welded ids are unchanged.
 */

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// What `weld_keep_first` returns.
struct WeldMap {
    std::vector<std::int64_t> mRepOf;      ///< per point: its representative's id
    std::vector<std::int64_t> mRepSource;  ///< per representative: the point that created it
};

/**
 * @brief Weld @p N points keep-first within @p Atol.
 * @param rXyz Coordinates, 3 per point (z-padded for 1-D and 2-D points).
 * @param Ddim How many leading coordinates the distance uses (1..3).
 */
inline WeldMap weld_keep_first(const std::vector<double>& rXyz, std::size_t N, std::size_t Ddim,
                               double Atol) {
    WeldMap out;
    out.mRepOf.assign(N, -1);
    if (N == 0)
        return out;
    const double atol2 = Atol * Atol;

    // Cell key per point (parallel), then dense cell ids in first-seen order
    // through a flat open-addressing table (serial, O(N) -- no allocation per
    // entry, unlike an unordered_map), and each cell's non-empty neighbour
    // cells, in the dz -> dy -> dx order the grid scan used, by read-only
    // lookups in that table (parallel).
    std::vector<GridKey> keys(N);
    parallel_for(N, [&](std::size_t g) {
        keys[g] = GridKey{grid_quantize(rXyz[g * 3], Atol), grid_quantize(rXyz[g * 3 + 1], Atol),
                          grid_quantize(rXyz[g * 3 + 2], Atol)};
    });
    std::size_t cap = 16;
    while (cap < 2 * N)
        cap <<= 1;
    const std::size_t mask = cap - 1;
    const auto slot_of = [mask](const GridKey& rK) {
        std::uint64_t h = static_cast<std::uint64_t>(rK.x) * 0x9E3779B97F4A7C15ull;
        h ^= static_cast<std::uint64_t>(rK.y) * 0xC2B2AE3D27D4EB4Full + (h >> 29);
        h ^= static_cast<std::uint64_t>(rK.z) * 0x165667B19E3779F9ull + (h >> 31);
        h ^= h >> 33;
        return static_cast<std::size_t>(h) & mask;
    };
    std::vector<std::int64_t> table(cap, -1);  // slot -> cell id
    std::vector<GridKey> cell_key;
    std::vector<std::int64_t> cell_of(N);
    for (std::size_t g = 0; g < N; ++g) {
        std::size_t s = slot_of(keys[g]);
        while (table[s] >= 0 && !(cell_key[static_cast<std::size_t>(table[s])] == keys[g]))
            s = (s + 1) & mask;
        if (table[s] < 0) {
            table[s] = static_cast<std::int64_t>(cell_key.size());
            cell_key.push_back(keys[g]);
        }
        cell_of[g] = table[s];
    }
    const std::size_t ncells = cell_key.size();
    const auto find_cell = [&](const GridKey& rK) -> std::int64_t {
        for (std::size_t s = slot_of(rK);; s = (s + 1) & mask) {
            const std::int64_t c = table[s];
            if (c < 0)
                return -1;
            if (cell_key[static_cast<std::size_t>(c)] == rK)
                return c;
        }
    };
    const auto for_neighbours = [&](std::size_t c, auto&& fn) {
        const GridKey k = cell_key[c];
        for (std::int64_t dz = -1; dz <= 1; ++dz)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dx = -1; dx <= 1; ++dx) {
                    const std::int64_t nb = find_cell(GridKey{k.x + dx, k.y + dy, k.z + dz});
                    if (nb >= 0)
                        fn(nb);
                }
    };
    std::vector<std::int64_t> nbr_start(ncells + 1, 0);
    parallel_for(ncells, [&](std::size_t c) {
        std::int64_t count = 0;
        for_neighbours(c, [&](std::int64_t) { ++count; });
        nbr_start[c + 1] = count;
    });
    for (std::size_t c = 0; c < ncells; ++c)
        nbr_start[c + 1] += nbr_start[c];
    std::vector<std::int64_t> nbrs(static_cast<std::size_t>(nbr_start.back()));
    parallel_for(ncells, [&](std::size_t c) {
        std::int64_t at = nbr_start[c];
        for_neighbours(c, [&](std::int64_t nb) { nbrs[static_cast<std::size_t>(at++)] = nb; });
    });

    // The keep-first decisions, serial in point order. A cell's representatives
    // form a list linked in creation order.
    std::vector<std::int64_t> head(ncells, -1), tail(ncells, -1), next;
    for (std::size_t g = 0; g < N; ++g) {
        const std::size_t c = static_cast<std::size_t>(cell_of[g]);
        std::int64_t found = -1;
        for (std::int64_t j = nbr_start[c]; j < nbr_start[c + 1] && found < 0; ++j) {
            for (std::int64_t r = head[static_cast<std::size_t>(nbrs[static_cast<std::size_t>(j)])];
                 r >= 0; r = next[static_cast<std::size_t>(r)]) {
                const std::size_t h =
                    static_cast<std::size_t>(out.mRepSource[static_cast<std::size_t>(r)]);
                double d2 = 0.0;
                for (std::size_t d = 0; d < Ddim; ++d) {
                    const double delta = rXyz[g * 3 + d] - rXyz[h * 3 + d];
                    d2 += delta * delta;
                }
                if (d2 <= atol2) {
                    found = r;
                    break;
                }
            }
        }
        if (found >= 0) {
            out.mRepOf[g] = found;
            continue;
        }
        const std::int64_t r = static_cast<std::int64_t>(out.mRepSource.size());
        out.mRepOf[g] = r;
        out.mRepSource.push_back(static_cast<std::int64_t>(g));
        next.push_back(-1);
        if (head[c] < 0)
            head[c] = r;
        else
            next[static_cast<std::size_t>(tail[c])] = r;
        tail[c] = r;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
