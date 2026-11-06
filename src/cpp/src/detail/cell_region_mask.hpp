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
 * @file detail/cell_region_mask.hpp
 * @brief The per-global-cell selection of a named Cell region, as the surface
 * operations take their `region` option.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

/// Per global cell, whether the region `rRegion` selects it (all cells when
/// unnamed). An unknown name throws `std::invalid_argument`, its message
/// prefixed with `pPrefix` and listing the available regions.
inline std::vector<std::uint8_t> cell_region_mask(const Mesh& rMesh, const std::string& rRegion,
                                                  std::size_t NumCells, const char* pPrefix) {
    if (rRegion.empty())
        return std::vector<std::uint8_t>(NumCells, 1);
    const std::size_t r = rMesh.FindRegion(rRegion, RegionKind::Cell);
    if (r == Mesh::npos) {
        std::string names;
        for (const std::string& n : rMesh.RegionNames())
            names += (names.empty() ? "" : ", ") + n;
        throw std::invalid_argument(std::string(pPrefix) + "no cell region named '" + rRegion +
                                    "' (available: " + (names.empty() ? "none" : names) + ")");
    }
    std::vector<std::uint8_t> selected(NumCells, 0);
    const NDArray& entries = rMesh.Region(r).mEntries;
    for (std::size_t e = 0; e < entries.Size(); ++e) {
        const std::int64_t c = detail::read_int(entries, e);
        if (c >= 0 && static_cast<std::size_t>(c) < NumCells)
            selected[static_cast<std::size_t>(c)] = 1;
    }
    return selected;
}

}  // namespace detail
}  // namespace meshioplusplus
