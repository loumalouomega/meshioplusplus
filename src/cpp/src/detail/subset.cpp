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

// Project includes
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

NDArray subset_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rIdx) {
    const std::vector<std::size_t>& shp = rSrc.Shape();
    std::size_t cols = 1;
    for (std::size_t d = 1; d < shp.size(); ++d)
        cols *= shp[d];
    std::vector<std::size_t> out_shape = shp;
    if (out_shape.empty())
        out_shape = {rIdx.size()};
    else
        out_shape[0] = rIdx.size();
    NDArray out = NDArray::Uninit(rSrc.Dtype(), out_shape);
    const std::size_t rb = cols * dtype_size(rSrc.Dtype());
    if (rb == 0)
        return out;
    const std::byte* s = rSrc.Data();
    std::byte* d = out.Data();
    parallel_for_bw(rIdx.size(), [&](std::size_t j) {
        std::memcpy(d + j * rb, s + static_cast<std::size_t>(rIdx[j]) * rb, rb);
    });
    return out;
}

SubsetResult build_cell_subset(const Mesh& rMesh,
                               const std::vector<std::vector<std::int64_t>>& rKeptCellsPerBlock,
                               const std::string& rPointIdName, const std::string& rCellIdName,
                               bool drop_empty_blocks) {
    SubsetResult res;
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    // 1) mark used points across all kept cells.
    std::vector<char> used(n, 0);
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::vector<std::int64_t>& kept = rKeptCellsPerBlock[b];
        if (cb.IsPolyhedron()) {
            for (std::int64_t c : kept)
                for (std::size_t f = 0; f < cb.NumFaces(static_cast<std::size_t>(c)); ++f) {
                    auto face = cb.Face(static_cast<std::size_t>(c), f);
                    for (std::size_t k = 0; k < face.second; ++k)
                        used[static_cast<std::size_t>(face.first[k])] = 1;
                }
        } else if (cb.IsRagged()) {
            for (std::int64_t c : kept)
                for (std::size_t k = 0; k < cb.RowSize(static_cast<std::size_t>(c)); ++k)
                    used[static_cast<std::size_t>(cb.Row(static_cast<std::size_t>(c))[k])] = 1;
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            for (std::int64_t c : kept)
                for (std::size_t k = 0; k < npc; ++k)
                    used[static_cast<std::size_t>(
                        read_int(conn, static_cast<std::size_t>(c) * npc + k))] = 1;
        }
        ++b;
    }

    // 2) compact used points -> new index; build gather list.
    std::vector<std::int64_t> new_point(n, -1);
    std::vector<std::int64_t> point_src;
    for (std::size_t g = 0; g < n; ++g)
        if (used[g]) {
            new_point[g] = static_cast<std::int64_t>(point_src.size());
            point_src.push_back(static_cast<std::int64_t>(g));
        }

    Mesh& out = res.mMesh;

    // 3) points.
    {
        NDArray newpts = NDArray::Uninit(points.Dtype(), {point_src.size(), dim});
        if (!point_src.empty() && dim > 0) {
            const std::size_t rb = dim * dtype_size(points.Dtype());
            const std::byte* s = points.Data();
            std::byte* d = newpts.Data();
            parallel_for_bw(point_src.size(), [&](std::size_t j) {
                std::memcpy(d + j * rb, s + static_cast<std::size_t>(point_src[j]) * rb, rb);
            });
        }
        out.AssignPoints(std::move(newpts));
    }

    // 4) point_data (gather by point_src when full-length).
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const std::size_t rows = a.Shape().empty() ? 0 : a.Shape()[0];
        if (rows == n)
            out.AddPointData(name, subset_gather_rows(a, point_src));
        else {
            NDArray c = a;
            c.MakeOwned();
            out.AddPointData(name, std::move(c));
        }
    }
    if (!rPointIdName.empty()) {
        NDArray ids = NDArray::Uninit(DType::Int64, {point_src.size()});
        std::int64_t* p = ids.As<std::int64_t>();
        for (std::size_t j = 0; j < point_src.size(); ++j)
            p[j] = point_src[j];
        out.AddPointData(rPointIdName, std::move(ids));
    }

    // 5) cells (remap connectivity to new point indices).
    b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::vector<std::int64_t>& kept = rKeptCellsPerBlock[b];
        if (drop_empty_blocks && kept.empty()) {
            ++b;
            continue;
        }
        const std::string type(cb.Type());
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(kept.size());
            for (std::size_t i = 0; i < kept.size(); ++i) {
                const std::size_t c = static_cast<std::size_t>(kept[i]);
                cells[i].resize(cb.NumFaces(c));
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    cells[i][f].reserve(face.second);
                    for (std::size_t k = 0; k < face.second; ++k)
                        cells[i][f].push_back(new_point[static_cast<std::size_t>(face.first[k])]);
                }
            }
            out.AddPolyhedronBlock(type, std::move(cells));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(kept.size());
            for (std::size_t i = 0; i < kept.size(); ++i) {
                const std::size_t c = static_cast<std::size_t>(kept[i]);
                rows[i].reserve(cb.RowSize(c));
                for (std::size_t k = 0; k < cb.RowSize(c); ++k)
                    rows[i].push_back(new_point[static_cast<std::size_t>(cb.Row(c)[k])]);
            }
            out.AddPolygonBlock(type, std::move(rows));
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            NDArray outconn = NDArray::Uninit(DType::Int64, {kept.size(), npc});
            std::int64_t* cd = outconn.As<std::int64_t>();
            for (std::size_t i = 0; i < kept.size(); ++i) {
                const std::size_t c = static_cast<std::size_t>(kept[i]);
                for (std::size_t k = 0; k < npc; ++k)
                    cd[i * npc + k] =
                        new_point[static_cast<std::size_t>(read_int(conn, c * npc + k))];
            }
            out.AddCellBlock(type, std::move(outconn));
        }
        ++b;
    }

    // 6) cell_data (gather kept cells per block) + optional original cell ids.
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> outblocks;
        b = 0;
        for (const auto cb : rMesh.CellRange()) {
            if (drop_empty_blocks && rKeptCellsPerBlock[b].empty()) {
                ++b;
                continue;
            }
            const NDArray& src = rMesh.CellData(name, b);
            const std::size_t rows = src.Shape().empty() ? 0 : src.Shape()[0];
            if (rows == cb.NumCells())
                outblocks.push_back(subset_gather_rows(src, rKeptCellsPerBlock[b]));
            else {
                NDArray c = src;
                c.MakeOwned();
                outblocks.push_back(std::move(c));
            }
            ++b;
        }
        out.AddCellData(name, std::move(outblocks));
    }
    if (!rCellIdName.empty()) {
        std::vector<NDArray> idblocks;
        for (const std::vector<std::int64_t>& kept : rKeptCellsPerBlock) {
            if (drop_empty_blocks && kept.empty())
                continue;
            NDArray ids = NDArray::Uninit(DType::Int64, {kept.size()});
            std::int64_t* p = ids.As<std::int64_t>();
            for (std::size_t i = 0; i < kept.size(); ++i)
                p[i] = kept[i];
            idblocks.push_back(std::move(ids));
        }
        out.AddCellData(rCellIdName, std::move(idblocks));
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        NDArray c = rMesh.FieldData(name);
        c.MakeOwned();
        out.AddFieldData(name, std::move(c));
    }

    // 7) index maps.
    res.mPointMap = NDArray::Uninit(DType::Int64, {n});
    std::int64_t* pm = res.mPointMap.As<std::int64_t>();
    for (std::size_t g = 0; g < n; ++g)
        pm[g] = new_point[g];

    b = 0;
    for (const auto cb : rMesh.CellRange()) {
        NDArray cmap = NDArray::Uninit(DType::Int64, {cb.NumCells()});
        std::int64_t* cm = cmap.As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            cm[c] = -1;
        const std::vector<std::int64_t>& kept = rKeptCellsPerBlock[b];
        for (std::size_t i = 0; i < kept.size(); ++i)
            cm[static_cast<std::size_t>(kept[i])] = static_cast<std::int64_t>(i);
        res.mCellMaps.push_back(std::move(cmap));
        ++b;
    }

    return res;
}

}  // namespace detail
}  // namespace meshioplusplus
