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
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {
namespace detail {

void parallel_copy_i64(std::int64_t* pDst, const std::int64_t* pSrc, std::size_t n) {
    if (n == 0)
        return;  // an empty block may have no buffer: memcpy(nullptr, ...) is undefined
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

namespace {

/// Rows `rIdx` of `rA` (any rank >= 1), as an owning array: what a polyhedron bucket needs,
/// because its members are not contiguous in the file.
NDArray vtkcells_gather_rows(const NDArray& rA, const std::vector<std::size_t>& rIdx) {
    const std::size_t nc = rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
    const std::size_t rowbytes = nc * dtype_size(rA.Dtype());
    std::vector<std::size_t> shape = rA.Shape();
    if (shape.empty())
        shape = {0};
    shape[0] = rIdx.size();
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);  // fully overwritten below
    for (std::size_t i = 0; i < rIdx.size(); ++i)
        std::memcpy(out.Data() + i * rowbytes, rA.Data() + rIdx[i] * rowbytes, rowbytes);
    return out;
}

}  // namespace

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

void check_vtk_cell_arrays(std::size_t ConnSize, const std::vector<std::int64_t>& rOffsets,
                           const std::vector<std::int64_t>& rTypes,
                           const std::unordered_map<std::string, NDArray>& rCellDataRaw) {
    if (rOffsets.size() != rTypes.size())
        throw ReadError("VTK: " + std::to_string(rOffsets.size()) + " cell offsets but " +
                        std::to_string(rTypes.size()) + " cell types");
    std::int64_t prev = 0;
    for (const std::int64_t o : rOffsets) {
        if (o < prev || static_cast<std::uint64_t>(o) > ConnSize)
            throw ReadError("VTK: cell offsets decrease or run past the connectivity");
        prev = o;
    }
    for (const auto& [name, arr] : rCellDataRaw) {
        const std::size_t rows = arr.Shape().empty() ? arr.Size() : arr.Shape()[0];
        if (rows < rTypes.size())
            throw ReadError("VTK: cell data '" + name + "' has " + std::to_string(rows) +
                            " rows for " + std::to_string(rTypes.size()) + " cells");
    }
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

    // The last polyhedral cell's END offset into *pFaces, carried across the
    // whole loop below (not reset per run): faceoffsets are END offsets, so a
    // polyhedron's face stream begins where the *previous polyhedral* cell's
    // ended, and a polyhedral run can be preceded by an intervening
    // non-polyhedral run. Only cells actually decoded in the vtk_type == 42
    // branch below ever advance it -- a non-polyhedral cell's -1 entry is
    // simply never visited here, which is what leaves it untouched.
    std::int64_t last_face_end = 0;

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
                const std::int64_t begin_at = last_face_end;
                const std::int64_t end_at = rFaceOffsets[c];
                if (end_at < 0 || begin_at > end_at ||
                    static_cast<std::size_t>(end_at) > pFaces->size())
                    throw ReadError("VTU: 'faceoffsets' entry is out of range for a polyhedron");
                last_face_end = end_at;
                std::size_t at = static_cast<std::size_t>(begin_at);
                const auto stop = static_cast<std::size_t>(end_at);
                // Every read stays inside this cell's slice of the stream.
                const auto take = [&]() {
                    if (at >= stop)
                        throw ReadError("VTU: a polyhedron's face stream overruns its entry");
                    return (*pFaces)[at++];
                };
                const std::int64_t nfaces = take();
                std::vector<std::vector<std::int64_t>> faces;
                std::vector<std::int64_t> uniq;
                for (std::int64_t f = 0; f < nfaces; ++f) {
                    const std::int64_t nn = take();
                    if (nn < 0 || static_cast<std::uint64_t>(nn) > stop - at)
                        throw ReadError("VTU: a polyhedron's face stream overruns its entry");
                    std::vector<std::int64_t> ring;
                    ring.reserve(static_cast<std::size_t>(nn));
                    for (std::int64_t k = 0; k < nn; ++k)
                        ring.push_back(take());
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
            for (std::size_t i = 0; i < cells.size(); ++i) {
                if (groups.find(node_counts[i]) == groups.end())
                    order.push_back(node_counts[i]);
                groups[node_counts[i]].push_back(i);
            }
            for (std::size_t n : order) {
                std::vector<std::vector<std::vector<std::int64_t>>> group;
                // The bucket's members sit at these FILE rows, which are not contiguous
                // when the run mixes node counts: slicing [at_row, at_row + m) would hand
                // each block another block's cell_data.
                std::vector<std::size_t> file_rows;
                for (std::size_t i : groups[n]) {
                    group.push_back(std::move(cells[i]));
                    file_rows.push_back(start + i);
                }
                rMesh.AddPolyhedronBlock("polyhedron" + std::to_string(n), std::move(group));
                for (const auto& kv : rCellDataRaw)
                    rMesh.AppendCellData(kv.first, vtkcells_gather_rows(kv.second, file_rows));
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
            // Each cell of a fixed-size type must span exactly n entries.
            for (std::size_t c = start; c < end; ++c)
                if (rOffsets[c] - (c == 0 ? 0 : rOffsets[c - 1]) != n)
                    throw ReadError("VTK: a '" + meshio_type + "' cell does not have " +
                                    std::to_string(n) + " nodes");
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

std::size_t vtk_xml_header_bytes(const Mesh& rMesh) {
    std::uint64_t items = 3 * static_cast<std::uint64_t>(rMesh.NumPoints());
    std::uint64_t conn = 0;
    std::uint64_t ncells = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t nc = cb.NumCells();
        ncells += nc;
        if (cb.IsPolyhedron()) {
            // The face stream [nfaces, [n, nodes...] per face] bounds both it
            // and the cell's (sorted unique) connectivity.
            for (std::size_t r = 0; r < nc; ++r) {
                conn += 1 + cb.NumFaces(r);
                for (std::size_t f = 0; f < cb.NumFaces(r); ++f)
                    conn += 1 + cb.Face(r, f).second;
            }
        } else if (cb.IsRagged()) {
            for (std::size_t r = 0; r < nc; ++r)
                conn += cb.RowSize(r);
        } else {
            conn += cb.Conn().Size();
        }
    }
    items = std::max({items, conn, ncells});
    for (const auto& name : rMesh.PointDataNames())
        items = std::max<std::uint64_t>(items, rMesh.PointData(name).Size());
    for (const auto& name : rMesh.CellDataNames()) {
        std::uint64_t n = 0;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            n += rMesh.CellData(name, b).Size();
        items = std::max(items, n);
    }
    for (const auto& name : rMesh.FieldDataNames())
        items = std::max<std::uint64_t>(items, rMesh.FieldData(name).Size());
    return vtu_header_bytes_for(8 * items);
}

}  // namespace detail
}  // namespace meshioplusplus
