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
 * @file nastran_model.hpp
 * @brief The Nastran model half shared by the result readers (`nastran_h5`,
 *        `nastran_op2`): element cards to cell blocks, element ids to cells,
 *        property ids to regions.
 *
 * Both readers get their model as rows per card (`EID`, `PID`, the `G`
 * connectivity) and a GRID id -> point index map, and build the same mesh from
 * them: the linear or quadratic cell type of each element (a quadratic one only
 * when every mid-side node is given), `cell_data["nastran:eid"]` and
 * `["nastran:pid"]`, and one cell region per property id named
 * `<PTYPE>_<pid>` (`PID_<pid>` when the property card is unknown). An element
 * on a scalar point is dropped; one on an undefined GRID is an error.
 *
 * Python twin: `src/python/meshioplusplus/nastran/_model.py`.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief One element card with a cell type: its linear and (optional) quadratic shape. */
struct NastranCardSpec {
    const char* mCard;
    const char* mLinear;
    std::size_t mLinearNodes;
    const char* mQuadratic;  ///< nullptr: no quadratic variant
    std::size_t mQuadraticNodes;
    const int* mPermutation;  ///< quadratic connectivity: conn[k] = G[perm[k]]; nullptr = identity
};

/** @brief The spec of an element card, or nullptr when it has no cell type. */
const NastranCardSpec* nastran_card_spec(std::string_view Card);

/** @brief The elements of one card: rows of `mWidth` grid ids each. */
struct NastranCardRows {
    std::string mCard;
    std::vector<std::int64_t> mEid;
    std::vector<std::int64_t> mPid;  ///< -1 where the card has none (CONM2)
    std::vector<std::int64_t> mNodes;
    std::size_t mWidth = 0;
};

/** @brief How elements map to the cells a Nastran model's blocks hold. */
struct NastranCells {
    /// EID -> global cell. A CONM2 has none: MSC accepts one sharing its id
    /// with a structural element, so it stays out.
    std::unordered_map<std::int64_t, std::size_t> mCellIndex;
    std::vector<std::size_t> mOffsets;  ///< first global cell of each block
    std::vector<std::size_t> mSizes;    ///< cells per block
    std::size_t mNumCells = 0;
};

/**
 * @brief Add the cell blocks, `nastran:eid`/`nastran:pid` and property regions
 *        of a model to @p rMesh, whose points are the GRIDs.
 * @param rCards the cards with a spec, in the order their blocks should take
 * @param rGridIndex GRID id -> point index
 * @param rScalarPoints SPOINT/EPOINT ids (an element on one is dropped)
 * @param rPtype property id -> property card name, for the region names
 * @param rWho the reader's name, prefixed to its warnings and errors
 * @throws ReadError when an element references an undefined GRID or no node
 */
NastranCells nastran_add_cells(Mesh& rMesh, const std::vector<NastranCardRows>& rCards,
                               const std::unordered_map<std::int64_t, std::size_t>& rGridIndex,
                               const std::unordered_set<std::int64_t>& rScalarPoints,
                               const std::map<std::int64_t, std::string>& rPtype,
                               const std::string& rWho);

/**
 * @brief Keep a GRID output or definition frame as point data when any is set:
 *        `nastran:cp` (coordinates stay in the local system) or `nastran:cd`
 *        (results are in the local output system), with a warning.
 * @param rFrame "CP" or "CD"
 */
void nastran_add_frame(Mesh& rMesh, const std::string& rFrame,
                       const std::vector<std::int64_t>& rValues, const std::string& rWho);

/**
 * @brief Read a bulk-data deck as `read_nastran` does, also returning the GRID
 *        ids in point order and the element ids in global cell order (the OP2
 *        reader's sibling-deck route).
 */
Mesh nastran_read_deck(const std::string& rPath, std::vector<std::int64_t>& rGridIds,
                       std::vector<std::int64_t>& rCellIds);

}  // namespace detail
}  // namespace meshioplusplus
