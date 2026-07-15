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

#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {
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
// `connectivity` just past cell i's last node. `conn` is a raw pointer so the
// caller can pass an int64 NDArray buffer directly (VTK 5.1 connectivity is
// already vtktypeint64) without an intermediate to_int64 copy.
inline void reconstruct_cells(
    const std::int64_t* conn, const std::vector<std::int64_t>& offsets,
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
                const std::size_t ii = i;
                // Contiguous uniform-size sub-run -> block memcpy.
                const std::int64_t sub_first = start_cn[ii];
                bool sub_regular = true;
                for (std::size_t r = 0; sub_regular && r < m; ++r)
                    if (offsets[start + ii + r] !=
                        sub_first + static_cast<std::int64_t>(r + 1) * sz)
                        sub_regular = false;
                if (sub_regular) {
                    std::memcpy(out, conn + sub_first,
                                m * static_cast<std::size_t>(sz) * sizeof(std::int64_t));
                } else {
                    parallel_for_bw(m, [&](std::size_t r) {
                        std::int64_t endoff = offsets[start + ii + r];
                        std::int64_t base = endoff - sz;
                        for (std::int64_t c = 0; c < sz; ++c)
                            out[r * sz + c] = conn[base + c];
                    });
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
            const int* ord = order.empty() ? nullptr : order.data();
            const std::size_t ss = start;
            // Regular run (offsets advance by exactly n per cell) with identity
            // node order -> the run's connectivity is one contiguous slice:
            // block memcpy instead of a per-row gather.
            const std::int64_t first = (ss == 0) ? 0 : offsets[ss - 1];
            bool regular = true;
            for (std::size_t r = 0; regular && r < m; ++r)
                if (offsets[ss + r] != first + static_cast<std::int64_t>((r + 1) *
                                                                         static_cast<std::size_t>(n)))
                    regular = false;
            if (!ord && regular) {
                std::memcpy(out, conn + first,
                            m * static_cast<std::size_t>(n) * sizeof(std::int64_t));
            } else {
                parallel_for_bw(m, [&](std::size_t r) {
                    std::int64_t endoff = offsets[ss + r];
                    std::int64_t base = endoff - n;
                    for (int j = 0; j < n; ++j) {
                        int col = ord ? ord[j] : j;
                        out[r * n + j] = conn[base + col];
                    }
                });
            }
            out_cells.emplace_back(meshio_type, std::move(data));
            add_cd(start, end);
        }
        start = end;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
