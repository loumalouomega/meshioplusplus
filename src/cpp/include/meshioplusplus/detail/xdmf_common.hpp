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
 * @file xdmf_common.hpp
 * @brief XDMF cell-type-name maps and cell-data raw<->blocks conversion,
 * shared between the XDMF format implementation and HMF (which reuses
 * XDMF's topology names and raw cell-data layout).
 *
 * Ported from `src/meshio/xdmf/common.py` and the `raw_from_cell_data` /
 * `cell_data_from_raw` helpers in `src/meshio/_common.py`. XDMF (and HMF)
 * store per-cell-type-block data as one array *per cell type name string* in
 * the XML/XDMF Topology, and store cell_data for a mixed mesh as one
 * concatenated raw array per data name (all cell blocks laid end-to-end)
 * rather than one array per block — `concat_cell_data`/`split_raw_cell_data`
 * are what let this header's callers go between meshio's per-block
 * `cell_data` representation and that concatenated-raw representation.
 *
 * Each function here is called once per cell-type-name lookup or once per
 * data array (not per element), so bodies live in
 * `src/cpp/src/detail/xdmf_common.cpp` rather than inline here.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace xdmfcommon {

/**
 * @brief Maps a meshio cell-type name to its XDMF topology type name.
 * @param t meshio cell-type name (e.g. `"triangle"`, `"tetra10"`).
 * @return The corresponding XDMF `TopologyType` string (e.g. `"Triangle"`).
 * @throws WriteError if `t` has no XDMF equivalent.
 */
const char* meshio_to_xdmf(const std::string& rT);

/**
 * @brief Maps an XDMF topology type name to a meshio cell-type name.
 *
 * Accepts both the canonical XDMF spelling and common abbreviations some
 * writers emit (e.g. both `"Hexahedron_20"` and `"Hex_20"` map to
 * `"hexahedron20"`).
 * @param t XDMF `TopologyType` string as found in the file.
 * @return The corresponding meshio cell-type name.
 * @throws ReadError if `t` is not a recognized topology type.
 */
std::string xdmf_to_meshio(const std::string& rT);

/**
 * @brief Concatenates one cell-data name's per-block arrays along axis 0
 * into a single raw array, matching Python's `raw_from_cell_data`.
 *
 * Used when writing XDMF/HMF cell data for a mixed-cell-type mesh: XDMF
 * stores cell data as one flat array per data name (all blocks' rows
 * back-to-back) rather than one array per block.
 * @param rMesh The mesh whose cell data to concatenate.
 * @param rName The cell-data name; must have at least one block, all blocks
 *              sharing dtype/trailing shape.
 * @return A new array with the same trailing shape as the first block and
 *         first dimension equal to the sum of each block's row count.
 */
NDArray concat_cell_data(const Mesh& rMesh, const std::string& rName);

/**
 * @brief Splits a raw, whole-mesh cell-data array (as read from XDMF/HMF)
 * back into one `NDArray` per cell block, matching Python's
 * `cell_data_from_raw`.
 *
 * Inverse of `concat_cell_data`.
 * @param raw The concatenated array covering every cell block's rows,
 *            in cell-block order.
 * @param sizes Row count of each cell block, in the same order the blocks
 *              appear in `raw`; must sum to `raw`'s row count.
 * @return One `NDArray` per entry in `sizes`, each holding that block's slice.
 */
std::vector<NDArray> split_raw_cell_data(const NDArray& rRaw,
                                         const std::vector<std::size_t>& rSizes);

}  // namespace xdmfcommon
}  // namespace meshioplusplus
