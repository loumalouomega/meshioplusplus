#pragma once
//
// meshioplusplus::Mesh and meshioplusplus::CellBlock: the in-memory mesh representation,
// mirroring the fields of the Python meshio.Mesh. This is the type the C++
// readers produce and the writers consume; the binding layer converts it
// to/from the pure-Python meshio.Mesh at the boundary.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "ndarray.hpp"

namespace meshioplusplus {

struct CellBlock {
    std::string type;                  // meshio cell type, e.g. "triangle"
    NDArray data;                      // (num_cells, nodes_per_cell), integer dtype
    std::vector<std::string> tags;

    // Ragged (jagged) representations, used only for cell types whose rows do
    // not fit a rectangular buffer. Exactly one of `data` / `polygon_rows` /
    // `polyhedron_rows` is populated per block; the two ragged members are
    // empty for every rectangular block (all rectangular formats unaffected).
    //
    //  * polygon_rows    — 1-level ragged: a "polygon" block whose cells have
    //                      varying node counts (e.g. MED POG Voronoi meshes).
    //                      Row i = polygon_rows[i] = node ids of cell i.
    //  * polyhedron_rows — 2-level ragged: a "polyhedron" block. Cell i is a
    //                      list of faces; each face is a list of node ids.
    std::vector<std::vector<std::int64_t>> polygon_rows;
    std::vector<std::vector<std::vector<std::int64_t>>> polyhedron_rows;

    CellBlock() = default;
    CellBlock(std::string t, NDArray d) : type(std::move(t)), data(std::move(d)) {}

    bool is_ragged() const {
        return !polygon_rows.empty() || !polyhedron_rows.empty();
    }

    std::size_t num_cells() const {
        if (!polygon_rows.empty()) return polygon_rows.size();
        if (!polyhedron_rows.empty()) return polyhedron_rows.size();
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

}  // namespace meshioplusplus
