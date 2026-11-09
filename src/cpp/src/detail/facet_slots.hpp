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
 * @file detail/facet_slots.hpp
 * @brief `FacetSlotUtils`: the slot bookkeeping of the facet tables that
 * `extract_surface`/`extract_skin` (`surface.cpp`) and `smooth`'s boundary
 * pass (`smooth.cpp`) each built (roadmap §7).
 *
 * Both number every facet of every block they scan as one global slot range,
 * block after block, then count each facet's key with `detail::group_facets`.
 * A table block's slots are `cells x facets-per-cell`; a polyhedron block's
 * are its cells' faces, laid out by a per-cell prefix sum (`mFaceStart`) and
 * inverted by a per-slot owner (`mCellOfFace`). What each caller keys and
 * keeps stays in its own `size_of`/`fill` lambdas, so the hot loop is the
 * caller's; this header holds the shared plumbing around it.
 *
 * A **core-private** header (not installed): no API or ABI surface.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus::detail {

/// Static helpers for a global facet-slot numbering over several blocks.
struct FacetSlotUtils {
    /**
     * @brief A polyhedron block's face-slot prefix sum: `rFaceStart[c]` is cell
     * `c`'s first slot within the block, `rFaceStart[NumCells]` the block's
     * slot count, which is returned. `rBlock.NumFaces(c)` gives each count.
     */
    template <class TCellView>
    static std::size_t BuildFaceStarts(const TCellView& rBlock, std::size_t NumCells,
                                       std::vector<std::size_t>& rFaceStart) {
        rFaceStart.reserve(NumCells + 1);
        std::size_t at = 0;
        rFaceStart.push_back(0);
        for (std::size_t c = 0; c < NumCells; ++c) {
            at += rBlock.NumFaces(c);
            rFaceStart.push_back(at);
        }
        return at;
    }

    /// The inverse of `BuildFaceStarts`: the owning cell of each face slot.
    static void FillCellOfFace(const std::vector<std::size_t>& rFaceStart, std::size_t NumCells,
                               std::vector<std::uint32_t>& rCellOfFace) {
        rCellOfFace.resize(rFaceStart.back());
        parallel_for(NumCells, [&](std::size_t c) {
            for (std::size_t f = rFaceStart[c]; f < rFaceStart[c + 1]; ++f)
                rCellOfFace[f] = static_cast<std::uint32_t>(c);
        });
    }

    /**
     * @brief The block holding global slot @p Slot, given each block's first
     * slot in ascending order (one block needs no search).
     */
    static std::size_t SlotOwner(const std::vector<std::size_t>& rFirstSlots, std::size_t Slot) {
        if (rFirstSlots.size() == 1)
            return 0;
        return static_cast<std::size_t>(
                   std::upper_bound(rFirstSlots.begin(), rFirstSlots.end(), Slot) -
                   rFirstSlots.begin()) -
               1;
    }
};

}  // namespace meshioplusplus::detail
