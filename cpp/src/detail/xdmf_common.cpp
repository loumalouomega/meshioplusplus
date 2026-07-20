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
#include <cstring>
#include <unordered_map>

// Project includes
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace xdmfcommon {

const char* meshio_to_xdmf(const std::string& rT) {
    static const std::unordered_map<std::string, const char*> m = {
        {"vertex", "Polyvertex"},
        {"line", "Polyline"},
        {"line3", "Edge_3"},
        {"quad", "Quadrilateral"},
        {"quad8", "Quadrilateral_8"},
        {"quad9", "Quadrilateral_9"},
        {"pyramid", "Pyramid"},
        {"pyramid13", "Pyramid_13"},
        {"tetra", "Tetrahedron"},
        {"triangle", "Triangle"},
        {"triangle6", "Triangle_6"},
        {"tetra10", "Tetrahedron_10"},
        {"wedge", "Wedge"},
        {"wedge15", "Wedge_15"},
        {"wedge18", "Wedge_18"},
        {"hexahedron", "Hexahedron"},
        {"hexahedron20", "Hexahedron_20"},
        {"hexahedron24", "Hexahedron_24"},
        {"hexahedron27", "Hexahedron_27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw WriteError("XDMF: unsupported cell type " + rT);
    return it->second;
}

std::string xdmf_to_meshio(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"Polyvertex", "vertex"},
        {"Polyline", "line"},
        {"Edge_3", "line3"},
        {"Quadrilateral", "quad"},
        {"Quadrilateral_8", "quad8"},
        {"Quad_8", "quad8"},
        {"Quadrilateral_9", "quad9"},
        {"Quad_9", "quad9"},
        {"Pyramid", "pyramid"},
        {"Pyramid_13", "pyramid13"},
        {"Tetrahedron", "tetra"},
        {"Triangle", "triangle"},
        {"Triangle_6", "triangle6"},
        {"Tri_6", "triangle6"},
        {"Tetrahedron_10", "tetra10"},
        {"Tet_10", "tetra10"},
        {"Wedge", "wedge"},
        {"Wedge_15", "wedge15"},
        {"Wedge_18", "wedge18"},
        {"Hexahedron", "hexahedron"},
        {"Hexahedron_20", "hexahedron20"},
        {"Hex_20", "hexahedron20"},
        {"Hexahedron_24", "hexahedron24"},
        {"Hex_24", "hexahedron24"},
        {"Hexahedron_27", "hexahedron27"},
        {"Hex_27", "hexahedron27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw ReadError("XDMF: unsupported topology type " + rT);
    return it->second;
}

NDArray concat_cell_data(const Mesh& rMesh, const std::string& rName) {
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::size_t total_rows = 0;
    std::vector<std::size_t> shape = rMesh.CellData(rName, 0).Shape();
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto& bshape = rMesh.CellData(rName, b).Shape();
        total_rows += bshape.empty() ? 0 : bshape[0];
    }
    shape[0] = total_rows;
    NDArray out(rMesh.CellData(rName, 0).Dtype(), shape);
    std::size_t off = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& blk = rMesh.CellData(rName, b);
        std::memcpy(out.Data() + off, blk.Data(), blk.Nbytes());
        off += blk.Nbytes();
    }
    return out;
}

std::vector<NDArray> split_raw_cell_data(const NDArray& rRaw,
                                         const std::vector<std::size_t>& rSizes) {
    std::size_t ncols = rRaw.Ndim() >= 2 ? rRaw.Shape()[1] : 1;
    std::size_t off = 0;
    std::vector<NDArray> blocks;
    for (std::size_t bs : rSizes) {
        std::vector<std::size_t> bshape = rRaw.Shape();
        if (!bshape.empty())
            bshape[0] = bs;
        NDArray b(rRaw.Dtype(), bshape);
        std::size_t elems = bs * ncols;
        std::memcpy(b.Data(), rRaw.Data() + off * ncols * dtype_size(rRaw.Dtype()),
                    elems * dtype_size(rRaw.Dtype()));
        off += bs;
        blocks.push_back(std::move(b));
    }
    return blocks;
}

}  // namespace xdmfcommon
}  // namespace meshioplusplus
