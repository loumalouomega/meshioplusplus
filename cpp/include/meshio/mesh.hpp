#pragma once
//
// meshio::Mesh and meshio::CellBlock: the in-memory mesh representation,
// mirroring the fields of the Python meshio.Mesh. This is the type the C++
// readers produce and the writers consume; the binding layer converts it
// to/from the pure-Python meshio.Mesh at the boundary.

#include <map>
#include <string>
#include <vector>

#include "ndarray.hpp"

namespace meshio {

struct CellBlock {
    std::string type;                  // meshio cell type, e.g. "triangle"
    NDArray data;                      // (num_cells, nodes_per_cell), integer dtype
    std::vector<std::string> tags;

    CellBlock() = default;
    CellBlock(std::string t, NDArray d) : type(std::move(t)), data(std::move(d)) {}

    std::size_t num_cells() const {
        return data.shape().empty() ? 0 : data.shape()[0];
    }
};

struct Mesh {
    NDArray points;                                         // (num_points, dim)
    std::vector<CellBlock> cells;

    // Field data. cell_data holds one NDArray per cell block, in cells order.
    std::map<std::string, NDArray> point_data;
    std::map<std::string, std::vector<NDArray>> cell_data;
    std::map<std::string, NDArray> field_data;

    std::size_t num_points() const {
        return points.shape().empty() ? 0 : points.shape()[0];
    }
};

}  // namespace meshio
