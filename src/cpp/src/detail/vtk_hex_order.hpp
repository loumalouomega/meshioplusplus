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
 * @file detail/vtk_hex_order.hpp
 * @brief The edge order of VTK's higher-order hexahedra across file versions.
 *
 * A **core-private** header (the `typed_view.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * VTK 9.1 swapped the nodes of edges 10 and 11 of `VTK_LAGRANGE_HEXAHEDRON`
 * and `VTK_BEZIER_HEXAHEDRON` (both edges run along the third direction), and
 * marks files in the new order by version: VTK XML 2.2 and later, legacy 5.1.
 * Its readers convert an older file on read and its legacy writer converts on
 * writing 4.2. meshio++ holds the current order (`lagrange::vtk_lattice`), so a
 * reader converts a file older than those versions, and a writer of an older
 * version converts on the way out. The swap is its own inverse.
 */

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

/// Whether a VTK cell type id is a higher-order hexahedron (72 Lagrange, 79 Bezier).
inline bool vtk_hex_order_applies(std::int64_t VtkType) {
    return VtkType == 72 || VtkType == 79;
}

/**
 * @brief Swap the nodes of edges 10 and 11 of one higher-order hexahedron.
 * @param pNodes the cell's @p NumNodes nodes, swapped in place
 * @param NumNodes the cell's node count
 * @param pDegrees its degrees per direction (VTK's `HigherOrderDegrees`), or
 *        `nullptr` for an isotropic cell (its degree follows from the count)
 * @return false when the degrees do not account for the node count (nothing done)
 */
inline bool vtk_hex_swap_edges(std::int64_t* pNodes, std::size_t NumNodes,
                               const std::int64_t* pDegrees) {
    std::int64_t p[3];
    if (pDegrees) {
        std::copy(pDegrees, pDegrees + 3, p);
    } else {
        const auto q =
            static_cast<std::int64_t>(std::llround(std::cbrt(static_cast<double>(NumNodes))));
        p[0] = p[1] = p[2] = q - 1;
    }
    if (p[0] < 1 || p[1] < 1 || p[2] < 1 ||
        (p[0] + 1) * (p[1] + 1) * (p[2] + 1) != static_cast<std::int64_t>(NumNodes))
        return false;
    const std::int64_t m = p[2] - 1;  // nodes on each edge along the third direction
    const std::int64_t s = 8 + 4 * (p[0] - 1) + 4 * (p[1] - 1) + 2 * m;
    std::swap_ranges(pNodes + s, pNodes + s + m, pNodes + s + m);
    return true;
}

/**
 * @brief Swap edges 10 and 11 of every higher-order hexahedron of a VTK cell
 * stream (`rTypes` VTK ids, `rOffsets` the end offsets into @p pConn).
 *
 * @p pDegrees is the `HigherOrderDegrees` cell array (one row of three per
 * cell), or `nullptr`; a cell whose degrees cannot be told is left alone, with
 * one warning naming @p pWhat.
 */
inline void vtk_hex_swap_stream(std::int64_t* pConn, const std::vector<std::int64_t>& rOffsets,
                                const std::vector<std::int64_t>& rTypes, const NDArray* pDegrees,
                                const char* pWhat) {
    const bool degrees = pDegrees && pDegrees->Size() == 3 * rTypes.size() && !rTypes.empty();
    bool warned = false;
    for (std::size_t c = 0; c < rTypes.size() && c < rOffsets.size(); ++c) {
        if (!vtk_hex_order_applies(rTypes[c]))
            continue;
        const std::int64_t start = c == 0 ? 0 : rOffsets[c - 1];
        const std::int64_t n = rOffsets[c] - start;
        std::int64_t deg[3];
        if (degrees)
            for (std::size_t k = 0; k < 3; ++k)
                deg[k] = read_int(*pDegrees, 3 * c + k);
        if (n > 0 &&
            !vtk_hex_swap_edges(pConn + start, static_cast<std::size_t>(n),
                                degrees ? deg : nullptr) &&
            !warned) {
            log::warn(
                "{}: a higher-order hexahedron of {} nodes in a file older than VTK 9.1 has no "
                "degrees that explain its size; its edge order is left as written",
                pWhat, n);
            warned = true;
        }
    }
}

}  // namespace detail
}  // namespace meshioplusplus
