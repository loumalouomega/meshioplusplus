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

// Project includes
#include "meshioplusplus/detail/cell_adjacency.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"

namespace meshioplusplus {
namespace detail {

CellIncidence build_cell_node_incidence(const Mesh& rMesh, std::size_t NumPoints) {
    CellIncidence csr;
    csr.mOffsets.push_back(0);
    // The node walk goes through cell_node_ids, which reads connectivity
    // dtype-agnostically. Reading it as As<std::int64_t>() would be wrong on a
    // MESHIO-backed mesh carrying Int32 connectivity from numpy -- see the
    // header comment.
    std::vector<std::int64_t> row;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        for (std::size_t c = 0; c < ncells; ++c) {
            cell_node_ids(cb, c, NumPoints, row);
            csr.mEntries.insert(csr.mEntries.end(), row.begin(), row.end());
            csr.mOffsets.push_back(static_cast<std::int64_t>(csr.mEntries.size()));
        }
    }
    return csr;
}

CellIncidence invert_cell_node_incidence(const CellIncidence& rCellNodes, std::size_t NumPoints,
                                         std::size_t NumCells) {
    CellIncidence csr;
    csr.mOffsets.assign(NumPoints + 1, 0);
    for (std::int64_t nid : rCellNodes.mEntries)
        if (nid >= 0 && static_cast<std::size_t>(nid) < NumPoints)
            ++csr.mOffsets[static_cast<std::size_t>(nid) + 1];
    for (std::size_t i = 0; i < NumPoints; ++i)
        csr.mOffsets[i + 1] += csr.mOffsets[i];
    csr.mEntries.assign(static_cast<std::size_t>(csr.mOffsets.back()), 0);
    // Walking cells in ascending order leaves every row ascending by cell id.
    std::vector<std::int64_t> cursor(csr.mOffsets.begin(), csr.mOffsets.end() - 1);
    for (std::size_t c = 0; c < NumCells; ++c)
        for (std::int64_t k = rCellNodes.mOffsets[c]; k < rCellNodes.mOffsets[c + 1]; ++k) {
            const std::int64_t nid = rCellNodes.mEntries[static_cast<std::size_t>(k)];
            if (nid >= 0 && static_cast<std::size_t>(nid) < NumPoints)
                csr.mEntries[static_cast<std::size_t>(cursor[static_cast<std::size_t>(nid)]++)] =
                    static_cast<std::int64_t>(c);
        }
    return csr;
}

void cell_node_neighbors(const CellIncidence& rCellNodes, const CellIncidence& rNodeCells,
                         std::size_t Cell, std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (Cell + 1 >= rCellNodes.mOffsets.size())
        return;
    const std::int64_t self = static_cast<std::int64_t>(Cell);
    for (std::int64_t k = rCellNodes.mOffsets[Cell]; k < rCellNodes.mOffsets[Cell + 1]; ++k) {
        const std::int64_t nid = rCellNodes.mEntries[static_cast<std::size_t>(k)];
        if (nid < 0 || static_cast<std::size_t>(nid) + 1 >= rNodeCells.mOffsets.size())
            continue;
        for (std::int64_t j = rNodeCells.mOffsets[static_cast<std::size_t>(nid)];
             j < rNodeCells.mOffsets[static_cast<std::size_t>(nid) + 1]; ++j) {
            const std::int64_t other = rNodeCells.mEntries[static_cast<std::size_t>(j)];
            if (other != self)
                rOut.push_back(other);
        }
    }
    // Ascending and de-duplicated: the neighbour order is part of this helper's
    // contract, since a sum over neighbours must not depend on how the walk
    // happened to encounter them.
    std::sort(rOut.begin(), rOut.end());
    rOut.erase(std::unique(rOut.begin(), rOut.end()), rOut.end());
}

}  // namespace detail
}  // namespace meshioplusplus
