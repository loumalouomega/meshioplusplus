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
 * @file detail/vertex_cell_csr.hpp
 * @brief The cells around each vertex of a fixed-corner cell list, as CSR:
 * shared by the hexahedral fit and the isosurface stuffing.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// The cells of `rCells` (`K` vertex ids each) around each of `NumPoints`
/// vertices: those of vertex `v` are `rAt[rStart[v] .. rStart[v + 1])`,
/// ascending, a cell once per corner it has there.
template <std::size_t K>
inline void vertex_cell_csr(const std::vector<std::array<std::int64_t, K>>& rCells,
                            std::size_t NumPoints, std::vector<std::int64_t>& rStart,
                            std::vector<std::int64_t>& rAt) {
    const std::size_t np = NumPoints;
    rStart.assign(np + 1, 0);
    for (const auto& c : rCells)
        for (std::int64_t v : c)
            ++rStart[static_cast<std::size_t>(v) + 1];
    for (std::size_t v = 0; v < np; ++v)
        rStart[v + 1] += rStart[v];
    rAt.resize(static_cast<std::size_t>(rStart[np]));
    {
        std::vector<std::int64_t> fill(rStart.begin(), rStart.end() - 1);
        for (std::size_t c = 0; c < rCells.size(); ++c)
            for (std::int64_t v : rCells[c])
                rAt[static_cast<std::size_t>(fill[static_cast<std::size_t>(v)]++)] =
                    static_cast<std::int64_t>(c);
    }
}

}  // namespace detail
}  // namespace meshioplusplus
