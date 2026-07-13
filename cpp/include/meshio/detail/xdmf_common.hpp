#pragma once
//
// XDMF type maps and cell-data raw<->blocks conversion, shared between the
// XDMF format (xdmf.cpp) and the HMF format (hmf.cpp reuses XDMF's topology
// names and raw cell-data layout). Ported from src/meshio/xdmf/common.py and
// src/meshio/_common.py (raw_from_cell_data / cell_data_from_raw).

#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

#include "meshio/exceptions.hpp"
#include "meshio/ndarray.hpp"

namespace meshio {
namespace xdmfcommon {

inline const char* meshio_to_xdmf(const std::string& t) {
    static const std::unordered_map<std::string, const char*> m = {
        {"vertex", "Polyvertex"}, {"line", "Polyline"}, {"line3", "Edge_3"},
        {"quad", "Quadrilateral"}, {"quad8", "Quadrilateral_8"},
        {"quad9", "Quadrilateral_9"}, {"pyramid", "Pyramid"},
        {"pyramid13", "Pyramid_13"}, {"tetra", "Tetrahedron"},
        {"triangle", "Triangle"}, {"triangle6", "Triangle_6"},
        {"tetra10", "Tetrahedron_10"}, {"wedge", "Wedge"}, {"wedge15", "Wedge_15"},
        {"wedge18", "Wedge_18"}, {"hexahedron", "Hexahedron"},
        {"hexahedron20", "Hexahedron_20"}, {"hexahedron24", "Hexahedron_24"},
        {"hexahedron27", "Hexahedron_27"}};
    auto it = m.find(t);
    if (it == m.end()) throw WriteError("XDMF: unsupported cell type " + t);
    return it->second;
}

inline std::string xdmf_to_meshio(const std::string& t) {
    static const std::unordered_map<std::string, std::string> m = {
        {"Polyvertex", "vertex"}, {"Polyline", "line"}, {"Edge_3", "line3"},
        {"Quadrilateral", "quad"}, {"Quadrilateral_8", "quad8"}, {"Quad_8", "quad8"},
        {"Quadrilateral_9", "quad9"}, {"Quad_9", "quad9"}, {"Pyramid", "pyramid"},
        {"Pyramid_13", "pyramid13"}, {"Tetrahedron", "tetra"}, {"Triangle", "triangle"},
        {"Triangle_6", "triangle6"}, {"Tri_6", "triangle6"}, {"Tetrahedron_10", "tetra10"},
        {"Tet_10", "tetra10"}, {"Wedge", "wedge"}, {"Wedge_15", "wedge15"},
        {"Wedge_18", "wedge18"}, {"Hexahedron", "hexahedron"},
        {"Hexahedron_20", "hexahedron20"}, {"Hex_20", "hexahedron20"},
        {"Hexahedron_24", "hexahedron24"}, {"Hex_24", "hexahedron24"},
        {"Hexahedron_27", "hexahedron27"}, {"Hex_27", "hexahedron27"}};
    auto it = m.find(t);
    if (it == m.end()) throw ReadError("XDMF: unsupported topology type " + t);
    return it->second;
}

// raw_from_cell_data: concatenate one name's per-block arrays along axis 0.
inline NDArray concat_cell_data(const std::vector<NDArray>& blocks) {
    std::size_t total_rows = 0;
    std::vector<std::size_t> shape = blocks.front().shape();
    for (const auto& b : blocks) total_rows += b.shape().empty() ? 0 : b.shape()[0];
    shape[0] = total_rows;
    NDArray out(blocks.front().dtype(), shape);
    std::size_t off = 0;
    for (const auto& b : blocks) {
        std::memcpy(out.data() + off, b.data(), b.nbytes());
        off += b.nbytes();
    }
    return out;
}

// cell_data_from_raw: split a concatenated array by the per-block cell counts.
inline std::vector<NDArray> split_raw_cell_data(const NDArray& raw,
                                                const std::vector<std::size_t>& sizes) {
    std::size_t ncols = raw.ndim() >= 2 ? raw.shape()[1] : 1;
    std::size_t off = 0;
    std::vector<NDArray> blocks;
    for (std::size_t bs : sizes) {
        std::vector<std::size_t> bshape = raw.shape();
        if (!bshape.empty()) bshape[0] = bs;
        NDArray b(raw.dtype(), bshape);
        std::size_t elems = bs * ncols;
        std::memcpy(b.data(), raw.data() + off * ncols * dtype_size(raw.dtype()),
                    elems * dtype_size(raw.dtype()));
        off += bs;
        blocks.push_back(std::move(b));
    }
    return blocks;
}

}  // namespace xdmfcommon
}  // namespace meshio
