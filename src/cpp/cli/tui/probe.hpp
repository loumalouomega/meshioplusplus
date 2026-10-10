// SPDX-License-Identifier: MIT
/// @file probe.hpp
/// @brief What the interactive viewer says about the cell under the pointer.
#ifndef MESHIOPLUSPLUS_CLI_TUI_PROBE_HPP
#define MESHIOPLUSPLUS_CLI_TUI_PROBE_HPP

// The id buffer of a frame names, per pixel, the input cell drawn there (global,
// block-major). A probe turns that number into facts: which block and cell type,
// the node ids, every cell data array's value there, the point data across its
// nodes, and the named cell regions that hold it.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::cli::tui {

/// A probed cell: text lines for the screen and the numbers behind them (for
/// comparing two pinned probes).
struct Probe {
    std::int64_t mCell = -1;
    std::vector<std::string> mLines;
    /// Scalar cell data by name, and the mean of each scalar point array over the
    /// cell's nodes. Vector arrays contribute their magnitude.
    std::map<std::string, double> mValues;
};

/// Describe the input cell `Cell` (a global, block-major index). A negative or
/// out-of-range index gives an empty probe (`mCell` -1).
Probe probe_cell(const Mesh& rMesh, std::int64_t Cell);

/// `B - A` over the names both probes have, as `name +1.5, other -2` (empty when
/// they share none).
std::string probe_difference(const Probe& rA, const Probe& rB);

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_PROBE_HPP
