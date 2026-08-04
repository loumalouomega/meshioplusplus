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
#include <algorithm>
#include <map>
#include <cstring>

// Project includes
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {
namespace detail {

void parallel_copy_i64(std::int64_t* pDst, const std::int64_t* pSrc, std::size_t n) {
    constexpr std::size_t kChunk = 1u << 19;  // 512Ki elements (4 MiB) per task
    const std::size_t nchunks = (n + kChunk - 1) / kChunk;
    if (nchunks <= 1) {
        std::memcpy(pDst, pSrc, n * sizeof(std::int64_t));
        return;
    }
    // grain=1: each chunk is already coarse (4 MiB), so dispatch per chunk —
    // otherwise the default grain (2048) would run these few chunks serially.
    parallel_for_bw(
        nchunks,
        [&](std::size_t c) {
            const std::size_t off = c * kChunk;
            const std::size_t len = std::min(kChunk, n - off);
            std::memcpy(pDst + off, pSrc + off, len * sizeof(std::int64_t));
        },
        1);
}

NDArray slice_rows(const NDArray& rA, std::size_t r0, std::size_t r1) {
    std::size_t nc = rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
    std::size_t isz = dtype_size(rA.Dtype());
    std::size_t rowbytes = nc * isz;
    std::vector<std::size_t> shape = rA.Shape();
    if (shape.empty())
        shape = {0};
    shape[0] = r1 - r0;
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);  // fully overwritten below
    if (r1 > r0)
        std::memcpy(out.Data(), rA.Data() + r0 * rowbytes, (r1 - r0) * rowbytes);
    return out;
}

std::vector<CellBlockInfo> summarize_cells(const std::vector<std::int64_t>& rOffsets,
                                           const std::vector<std::int64_t>& rTypes) {
    const auto& vmap = vtk_to_meshio_type();
    const std::size_t ncells = rTypes.size();
    std::vector<CellBlockInfo> blocks;

    std::size_t start = 0;
    while (start < ncells) {
        std::size_t end = start + 1;
        while (end < ncells && rTypes[end] == rTypes[start])
            ++end;

        const int vtk_type = static_cast<int>(rTypes[start]);
        if (vtk_type == 42)
            throw ReadError("polyhedron cells are not supported by the C++ reader");
        auto it = vmap.find(vtk_type);
        if (it == vmap.end())
            throw ReadError("VTK cell type " + std::to_string(vtk_type) +
                            " not supported by the C++ reader");
        const std::string& meshio_type = it->second;

        if (is_special_cell(meshio_type)) {
            if (rOffsets.size() < end)
                throw ReadError("VTU summary needs 'offsets' for variable-size cell type " +
                                meshio_type);
            // Split the run further wherever the per-cell node count changes.
            std::size_t i = start;
            while (i < end) {
                const std::int64_t prev = (i == 0) ? 0 : rOffsets[i - 1];
                const std::int64_t sz = rOffsets[i] - prev;
                std::size_t j = i + 1;
                while (j < end && rOffsets[j] - rOffsets[j - 1] == sz)
                    ++j;
                CellBlockInfo info;
                info.mType = meshio_type;
                info.mNumCells = j - i;
                info.mNodesPerCell = static_cast<std::size_t>(sz);
                blocks.push_back(std::move(info));
                i = j;
            }
        } else {
            auto nit = num_nodes_per_cell().find(meshio_type);
            if (nit == num_nodes_per_cell().end())
                throw ReadError("Unknown node count for cell type " + meshio_type);
            CellBlockInfo info;
            info.mType = meshio_type;
            info.mNumCells = end - start;
            info.mNodesPerCell = static_cast<std::size_t>(nit->second);
            blocks.push_back(std::move(info));
        }
        start = end;
    }
    return blocks;
}

bool cells_need_offsets(const std::vector<std::int64_t>& rTypes) {
    const auto& vmap = vtk_to_meshio_type();
    for (std::int64_t t : rTypes) {
        auto it = vmap.find(static_cast<int>(t));
        if (it != vmap.end() && is_special_cell(it->second))
            return true;
    }
    return false;
}

void reconstruct_cells(const std::int64_t* pConn, const std::vector<std::int64_t>& rOffsets,
                       const std::vector<std::int64_t>& rTypes,
                       const std::unordered_map<std::string, NDArray>& rCellDataRaw, Mesh& rMesh) {
    // The historical four-argument form: no faces stream, so type 42 still
    // refuses by name (see the header -- this must stay a distinct symbol).
    static const std::vector<std::int64_t> kNoFaceOffsets;
    reconstruct_cells(pConn, rOffsets, rTypes, rCellDataRaw, nullptr, kNoFaceOffsets, rMesh);
}

void reconstruct_cells(const std::int64_t* pConn, const std::vector<std::int64_t>& rOffsets,
                       const std::vector<std::int64_t>& rTypes,
                       const std::unordered_map<std::string, NDArray>& rCellDataRaw,
                       const std::vector<std::int64_t>* pFaces,
                       const std::vector<std::int64_t>& rFaceOffsets, Mesh& rMesh) {
    const auto& vmap = vtk_to_meshio_type();
    const std::size_t ncells = rTypes.size();

    auto add_cd = [&](std::size_t start, std::size_t end) {
        for (const auto& kv : rCellDataRaw)
            rMesh.AppendCellData(kv.first, slice_rows(kv.second, start, end));
    };

    std::size_t start = 0;
    while (start < ncells) {
        std::size_t end = start + 1;
        while (end < ncells && rTypes[end] == rTypes[start])
            ++end;

        int vtk_type = static_cast<int>(rTypes[start]);
        if (vtk_type == 42) {
            if (pFaces == nullptr || rFaceOffsets.size() != ncells)
                throw ReadError(
                    "VTU: a cell has VTK type 42 (polyhedron) but the file carries no usable "
                    "'faces'/'faceoffsets' arrays");
            // Decode this run of polyhedra, then bucket by unique node count
            // into polyhedron<N> -- the convention the OpenFOAM, EnSight, MED
            // and CGNS readers all use.
            std::vector<std::vector<std::vector<std::int64_t>>> cells;
            std::vector<std::size_t> node_counts;
            for (std::size_t c = start; c < end; ++c) {
                // faceoffsets are END offsets, so this cell's stream begins
                // where the previous polyhedral cell's ended. A non-polyhedral
                // cell carries -1 and contributes nothing.
                std::int64_t begin_at = 0;
                for (std::size_t q = 0; q < c; ++q)
                    if (rFaceOffsets[q] >= 0)
                        begin_at = rFaceOffsets[q];
                const std::int64_t end_at = rFaceOffsets[c];
                if (end_at < 0 || begin_at > end_at ||
                    static_cast<std::size_t>(end_at) > pFaces->size())
                    throw ReadError("VTU: 'faceoffsets' entry is out of range for a polyhedron");
                std::size_t at = static_cast<std::size_t>(begin_at);
                const std::int64_t nfaces = (*pFaces)[at++];
                std::vector<std::vector<std::int64_t>> faces;
                std::vector<std::int64_t> uniq;
                for (std::int64_t f = 0; f < nfaces; ++f) {
                    const std::int64_t nn = (*pFaces)[at++];
                    std::vector<std::int64_t> ring;
                    ring.reserve(static_cast<std::size_t>(nn));
                    for (std::int64_t k = 0; k < nn; ++k)
                        ring.push_back((*pFaces)[at++]);
                    uniq.insert(uniq.end(), ring.begin(), ring.end());
                    faces.push_back(std::move(ring));
                }
                std::sort(uniq.begin(), uniq.end());
                uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                node_counts.push_back(uniq.size());
                cells.push_back(std::move(faces));
            }
            std::vector<std::size_t> order;
            std::map<std::size_t, std::vector<std::size_t>> groups;
            std::size_t at_row = start;
            for (std::size_t i = 0; i < cells.size(); ++i) {
                if (groups.find(node_counts[i]) == groups.end())
                    order.push_back(node_counts[i]);
                groups[node_counts[i]].push_back(i);
            }
            for (std::size_t n : order) {
                std::vector<std::vector<std::vector<std::int64_t>>> group;
                for (std::size_t i : groups[n])
                    group.push_back(std::move(cells[i]));
                const std::size_t m = group.size();
                rMesh.AddPolyhedronBlock("polyhedron" + std::to_string(n), std::move(group));
                for (const auto& kv : rCellDataRaw)
                    rMesh.AppendCellData(kv.first, slice_rows(kv.second, at_row, at_row + m));
                at_row += m;
            }
            start = end;
            continue;
        }
        auto it = vmap.find(vtk_type);
        if (it == vmap.end())
            throw ReadError("VTK cell type " + std::to_string(vtk_type) +
                            " not supported by the C++ reader");
        const std::string& meshio_type = it->second;

        if (is_special_cell(meshio_type)) {
            std::int64_t first_node = (start == 0) ? 0 : rOffsets[start - 1];
            std::vector<std::int64_t> start_cn;
            start_cn.reserve(end - start + 1);
            start_cn.push_back(first_node);
            for (std::size_t i = start; i < end; ++i)
                start_cn.push_back(rOffsets[i]);
            std::vector<std::int64_t> sizes(end - start);
            for (std::size_t i = 0; i < sizes.size(); ++i)
                sizes[i] = start_cn[i + 1] - start_cn[i];

            std::size_t i = 0;
            while (i < sizes.size()) {
                std::size_t j = i;
                while (j < sizes.size() && sizes[j] == sizes[i])
                    ++j;
                std::int64_t sz = sizes[i];
                std::size_t m = j - i;
                NDArray data = NDArray::Uninit(DType::Int64, {m, static_cast<std::size_t>(sz)});
                std::int64_t* out = data.As<std::int64_t>();
                const std::size_t ii = i;
                // Contiguous uniform-size sub-run -> block memcpy.
                const std::int64_t sub_first = start_cn[ii];
                bool sub_regular = true;
                for (std::size_t r = 0; sub_regular && r < m; ++r)
                    if (rOffsets[start + ii + r] !=
                        sub_first + static_cast<std::int64_t>(r + 1) * sz)
                        sub_regular = false;
                if (sub_regular) {
                    parallel_copy_i64(out, pConn + sub_first, m * static_cast<std::size_t>(sz));
                } else {
                    parallel_for_bw(m, [&](std::size_t r) {
                        std::int64_t endoff = rOffsets[start + ii + r];
                        std::int64_t base = endoff - sz;
                        for (std::int64_t c = 0; c < sz; ++c)
                            out[r * sz + c] = pConn[base + c];
                    });
                }
                rMesh.AddCellBlock(meshio_type, std::move(data));
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
            NDArray data = NDArray::Uninit(DType::Int64, {m, static_cast<std::size_t>(n)});
            std::int64_t* out = data.As<std::int64_t>();
            const int* ord = order.empty() ? nullptr : order.data();
            const std::size_t ss = start;
            // Regular run (offsets advance by exactly n per cell) with identity
            // node order -> the run's connectivity is one contiguous slice:
            // block memcpy instead of a per-row gather.
            const std::int64_t first = (ss == 0) ? 0 : rOffsets[ss - 1];
            bool regular = true;
            for (std::size_t r = 0; regular && r < m; ++r)
                if (rOffsets[ss + r] !=
                    first + static_cast<std::int64_t>((r + 1) * static_cast<std::size_t>(n)))
                    regular = false;
            if (!ord && regular) {
                // Contiguous slice -> parallel block copy (fault-bound).
                parallel_copy_i64(out, pConn + first, m * static_cast<std::size_t>(n));
            } else {
                parallel_for_bw(m, [&](std::size_t r) {
                    std::int64_t endoff = rOffsets[ss + r];
                    std::int64_t base = endoff - n;
                    for (int j = 0; j < n; ++j) {
                        int col = ord ? ord[j] : j;
                        out[r * n + j] = pConn[base + col];
                    }
                });
            }
            rMesh.AddCellBlock(meshio_type, std::move(data));
            add_cd(start, end);
        }
        start = end;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
