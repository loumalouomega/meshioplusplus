#pragma once
//
// Shared reconstruction of meshio cell blocks from the VTK/VTU
// connectivity + end-offsets + types representation. Used by the VTU reader and
// the VTK 5.1 reader (which share this exact layout). Ported from
// vtk_cells_from_data in _vtk_common.py.

#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "meshio/exceptions.hpp"
#include "meshio/mesh.hpp"
#include "meshio/types.hpp"
#include "meshio/vtk_common.hpp"

namespace meshio {
namespace detail {

inline NDArray slice_rows(const NDArray& a, std::size_t r0, std::size_t r1) {
    std::size_t nc = a.shape().size() >= 2 ? a.shape()[1] : 1;
    std::size_t isz = dtype_size(a.dtype());
    std::size_t rowbytes = nc * isz;
    std::vector<std::size_t> shape = a.shape();
    if (shape.empty()) shape = {0};
    shape[0] = r1 - r0;
    NDArray out(a.dtype(), shape);
    if (r1 > r0)
        std::memcpy(out.data(), a.data() + r0 * rowbytes, (r1 - r0) * rowbytes);
    return out;
}

// `offsets` are end offsets (one per cell): offsets[i] is the index in
// `connectivity` just past cell i's last node.
inline void reconstruct_cells(
    const std::vector<std::int64_t>& conn, const std::vector<std::int64_t>& offsets,
    const std::vector<std::int64_t>& types,
    const std::map<std::string, NDArray>& cell_data_raw,
    std::vector<CellBlock>& out_cells,
    std::map<std::string, std::vector<NDArray>>& out_cell_data) {
    const auto& vmap = vtk_to_meshio_type();
    const std::size_t ncells = types.size();

    auto add_cd = [&](std::size_t start, std::size_t end) {
        for (const auto& kv : cell_data_raw)
            out_cell_data[kv.first].push_back(slice_rows(kv.second, start, end));
    };

    std::size_t start = 0;
    while (start < ncells) {
        std::size_t end = start + 1;
        while (end < ncells && types[end] == types[start]) ++end;

        int vtk_type = static_cast<int>(types[start]);
        if (vtk_type == 42)
            throw ReadError("polyhedron cells are not supported by the C++ reader");
        auto it = vmap.find(vtk_type);
        if (it == vmap.end())
            throw ReadError("VTK cell type " + std::to_string(vtk_type) +
                            " not supported by the C++ reader");
        const std::string& meshio_type = it->second;

        if (is_special_cell(meshio_type)) {
            std::int64_t first_node = (start == 0) ? 0 : offsets[start - 1];
            std::vector<std::int64_t> start_cn;
            start_cn.reserve(end - start + 1);
            start_cn.push_back(first_node);
            for (std::size_t i = start; i < end; ++i) start_cn.push_back(offsets[i]);
            std::vector<std::int64_t> sizes(end - start);
            for (std::size_t i = 0; i < sizes.size(); ++i)
                sizes[i] = start_cn[i + 1] - start_cn[i];

            std::size_t i = 0;
            while (i < sizes.size()) {
                std::size_t j = i;
                while (j < sizes.size() && sizes[j] == sizes[i]) ++j;
                std::int64_t sz = sizes[i];
                std::size_t m = j - i;
                NDArray data(DType::Int64, {m, static_cast<std::size_t>(sz)});
                std::int64_t* out = data.as<std::int64_t>();
                for (std::size_t r = 0; r < m; ++r) {
                    std::int64_t endoff = offsets[start + i + r];
                    std::int64_t base = endoff - sz;
                    for (std::int64_t c = 0; c < sz; ++c) out[r * sz + c] = conn[base + c];
                }
                out_cells.emplace_back(meshio_type, std::move(data));
                add_cd(start + i, start + j);
                i = j;
            }
        } else {
            auto nit = num_nodes_per_cell().find(meshio_type);
            if (nit == num_nodes_per_cell().end())
                throw ReadError("Unknown node count for cell type " + meshio_type);
            int n = nit->second;
            std::vector<int> order = vtk_to_meshio_order(vtk_type);
            std::size_t m = end - start;
            NDArray data(DType::Int64, {m, static_cast<std::size_t>(n)});
            std::int64_t* out = data.as<std::int64_t>();
            for (std::size_t r = 0; r < m; ++r) {
                std::int64_t endoff = offsets[start + r];
                std::int64_t base = endoff - n;
                for (int j = 0; j < n; ++j) {
                    int col = order.empty() ? j : order[j];
                    out[r * n + j] = conn[base + col];
                }
            }
            out_cells.emplace_back(meshio_type, std::move(data));
            add_cd(start, end);
        }
        start = end;
    }
}

}  // namespace detail
}  // namespace meshio
