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
 * @file detail/staged_block.hpp
 * @brief Copying an unchanged rectangular or polygon cell block into an
 * operation's staged output block, and adding such a block to the result
 * mesh, shared by the operations that stage their output (agglomerate,
 * subdivide).
 *
 * A **core-private** header: no installed header names it, and it adds
 * nothing to the API or the ABI. Built on the uniform mesh API only.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

#include "always_inline.hpp"

namespace meshioplusplus {
namespace detail {

/// Stage a non-polyhedron block unchanged into @p rOut: a ragged (polygon)
/// block as one row per cell in `mPolygonRows`, a rectangular one as its flat
/// connectivity in `mConn` with `mNodesPerCell`. The caller handles
/// polyhedron blocks itself, since each staging struct holds them its own way.
template <class TBlock>
MIO_ALWAYS_INLINE void stage_flat_or_polygon(const Mesh::CellView& rBlock, TBlock& rOut) {
    if (rBlock.IsRagged()) {
        rOut.mIsRagged = true;
        rOut.mPolygonRows.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            const std::int64_t* row = rBlock.Row(c);
            rOut.mPolygonRows[c].assign(row, row + rBlock.RowSize(c));
        }
    } else {
        const NDArray& conn = rBlock.Conn();
        const std::size_t npc = rBlock.NodesPerCell();
        rOut.mNodesPerCell = npc;
        rOut.mConn.resize(rBlock.NumCells() * npc);
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c)
            for (std::size_t k = 0; k < npc; ++k)
                rOut.mConn[c * npc + k] = read_int(conn, c * npc + k);
    }
}

/// Add a block staged by stage_flat_or_polygon (or built in the same shape)
/// to @p rOut: `mPolygonRows` as a polygon block when `mIsRagged`, otherwise
/// `mConn` as an (n, mNodesPerCell) Int64 block. Moves the rows out.
template <class TBlock>
MIO_ALWAYS_INLINE void emit_flat_or_polygon(Mesh& rOut, TBlock& rBlock) {
    if (rBlock.mIsRagged) {
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mPolygonRows));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? 0 : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rBlock.mConn.size(); ++i)
            dst[i] = rBlock.mConn[i];
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

}  // namespace detail
}  // namespace meshioplusplus
