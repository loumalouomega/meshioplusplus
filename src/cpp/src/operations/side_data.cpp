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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/side_data.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

bool side_extract_facet_nodes(const Mesh& rMesh, const std::vector<std::int64_t>& rBases,
                              std::int64_t GlobalCell, std::int64_t LocalFacet,
                              CellType& rFacetType, std::vector<std::int64_t>& rNodes) {
    if (GlobalCell < 0 || LocalFacet < 0)
        return false;
    auto it = std::upper_bound(rBases.begin(), rBases.end(), GlobalCell);
    if (it == rBases.begin())
        return false;
    const std::size_t bi = static_cast<std::size_t>(std::distance(rBases.begin(), it) - 1);
    if (bi >= rMesh.NumCellBlocks())
        return false;
    const auto cb = rMesh.Cells(bi);
    const std::size_t local_cell = static_cast<std::size_t>(GlobalCell - rBases[bi]);
    if (local_cell >= cb.NumCells())
        return false;
    if (cb.IsRagged())
        return false;

    const CellType ctype = cell_type_from_name(cb.Type());
    const int dim = cell_type_dimension(ctype);
    const detail::Int64View conn(cb.Conn());
    const std::size_t npc = cb.NodesPerCell();

    if (dim == 3) {
        const auto faces = detail::cell_faces(ctype);
        if (static_cast<std::size_t>(LocalFacet) >= faces.size())
            return false;
        const auto& fd = faces[static_cast<std::size_t>(LocalFacet)];
        rFacetType = fd.mFaceType;
        rNodes.resize(fd.mNumNodes);
        for (std::size_t k = 0; k < fd.mNumNodes; ++k)
            rNodes[k] = conn[local_cell * npc + fd.mNodes[k]];
        return true;
    }
    if (dim == 2) {
        const auto edges = detail::cell_edges(ctype);
        if (static_cast<std::size_t>(LocalFacet) >= edges.size())
            return false;
        const auto& ed = edges[static_cast<std::size_t>(LocalFacet)];
        rFacetType = ed.mEdgeType;
        rNodes.resize(ed.mNumNodes);
        for (std::size_t k = 0; k < ed.mNumNodes; ++k)
            rNodes[k] = conn[local_cell * npc + ed.mNodes[k]];
        return true;
    }
    return false;
}

Mesh side_clone_base(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged()) {
            if (cb.IsPolyhedron()) {
                std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    cells[r].resize(cb.NumFaces(r));
                    for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                        const auto face = cb.Face(r, f);
                        cells[r][f].assign(face.first, face.first + face.second);
                    }
                }
                out.AddPolyhedronBlock(cb.Type(), cells);
            } else {
                std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    const auto* pRow = cb.Row(r);
                    rows[r].assign(pRow, pRow + cb.RowSize(r));
                }
                out.AddPolygonBlock(cb.Type(), rows);
            }
        } else {
            out.AddCellBlock(cb.Type(), detail::data_owned_copy(cb.Conn()));
        }
    }
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        out.AddRegion(rMesh.Region(r));
    return out;
}

}  // namespace

SidesToCellsResult sides_to_cells(const Mesh& rMesh, const SidesToCellsOptions& rOptions) {
    SidesToCellsResult res;
    res.mMesh = side_clone_base(rMesh);

    const std::size_t orig_blocks = rMesh.NumCellBlocks();
    std::vector<std::int64_t> bases(orig_blocks + 1, 0);
    for (std::size_t b = 0; b < orig_blocks; ++b)
        bases[b + 1] = bases[b] + static_cast<std::int64_t>(rMesh.Cells(b).NumCells());

    // Select target Side regions.
    std::vector<std::size_t> target_regions;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind != RegionKind::Side)
            continue;
        if (!rOptions.mRegion.empty() && reg.mName != rOptions.mRegion)
            continue;
        target_regions.push_back(r);
    }
    if (!rOptions.mRegion.empty() && target_regions.empty())
        throw std::invalid_argument("meshio++: sides_to_cells: no side region named '" +
                                    rOptions.mRegion + "'");

    struct FacetBlockBuild {
        CellType mType;
        std::string mRegionName;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mSideOf;  // (global_cell, local_facet) pairs
        std::vector<std::size_t> mEntryIndices;
        std::size_t mNumNodes = 0;
    };

    std::vector<FacetBlockBuild> new_blocks;
    for (std::size_t ri : target_regions) {
        const Region& reg = rMesh.Region(ri);
        const std::size_t n = reg.NumEntries();
        const std::int64_t* e = reg.Entries();

        // Group facets within this side region by their cell type (preserving order).
        std::map<CellType, std::size_t> by_type;
        std::vector<std::int64_t> nodes;
        for (std::size_t i = 0; i < n; ++i) {
            CellType ftype = CellType::Triangle;
            if (!side_extract_facet_nodes(rMesh, bases, e[2 * i], e[2 * i + 1], ftype, nodes)) {
                res.mWarnings.push_back("sides_to_cells: skipped unsupported or out-of-range side "
                                        "entry in region '" +
                                        reg.mName + "'");
                continue;
            }
            auto it = by_type.find(ftype);
            if (it == by_type.end()) {
                const std::size_t idx = new_blocks.size();
                by_type.emplace(ftype, idx);
                FacetBlockBuild fb;
                fb.mType = ftype;
                fb.mRegionName = reg.mName;
                fb.mNumNodes = nodes.size();
                new_blocks.push_back(std::move(fb));
                it = by_type.find(ftype);
            }
            FacetBlockBuild& fb = new_blocks[it->second];
            fb.mConn.insert(fb.mConn.end(), nodes.begin(), nodes.end());
            fb.mSideOf.push_back(e[2 * i]);
            fb.mSideOf.push_back(e[2 * i + 1]);
            fb.mEntryIndices.push_back(i);
            ++res.mNumFacets;
        }
    }

    // Add new facet cell blocks and track their global cell index ranges.
    std::int64_t next_global = bases.back();
    std::map<std::string, std::vector<std::int64_t>> region_new_cells;
    for (const auto& fb : new_blocks) {
        const std::size_t nf = fb.mEntryIndices.size();
        NDArray conn(DType::Int64, {nf, fb.mNumNodes});
        std::copy(fb.mConn.begin(), fb.mConn.end(), conn.As<std::int64_t>());
        res.mMesh.AddCellBlock(cell_type_name(fb.mType), std::move(conn));
        auto& ids = region_new_cells[fb.mRegionName];
        for (std::size_t i = 0; i < nf; ++i)
            ids.push_back(next_global++);
    }
    for (auto& [rname, ids] : region_new_cells) {
        NDArray r_arr(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), r_arr.As<std::int64_t>());
        res.mMesh.AddRegion(Region(rname, RegionKind::Cell, std::move(r_arr)));
    }

    const std::size_t total_blocks = res.mMesh.NumCellBlocks();

    // Carry existing cell_data (padding new facet blocks with NaN / 0).
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        blocks.reserve(total_blocks);
        const std::size_t norig = rMesh.CellDataNumBlocks(name);
        const DType dt = norig > 0 ? rMesh.CellData(name, 0).Dtype() : DType::Float64;
        const std::size_t nc = norig > 0 ? data_num_components(rMesh.CellData(name, 0)) : 1;
        for (std::size_t b = 0; b < orig_blocks; ++b)
            blocks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
        for (std::size_t b = orig_blocks; b < total_blocks; ++b) {
            const std::size_t nf = res.mMesh.Cells(b).NumCells();
            NDArray pad(dt, nc == 1 ? std::vector<std::size_t>{nf}
                                    : std::vector<std::size_t>{nf, nc});
            if (dt == DType::Float64)
                std::fill(pad.As<double>(), pad.As<double>() + pad.Size(),
                          std::numeric_limits<double>::quiet_NaN());
            else if (dt == DType::Float32)
                std::fill(pad.As<float>(), pad.As<float>() + pad.Size(),
                          std::numeric_limits<float>::quiet_NaN());
            blocks.push_back(std::move(pad));
        }
        res.mMesh.AddCellData(name, std::move(blocks));
    }

    // Add `meshio:side_of` cell_data: (-1, -1) on original blocks, (cell, facet) on facet blocks.
    if (!new_blocks.empty()) {
        std::vector<NDArray> side_of_blocks;
        side_of_blocks.reserve(total_blocks);
        for (std::size_t b = 0; b < orig_blocks; ++b) {
            const std::size_t nc_blk = res.mMesh.Cells(b).NumCells();
            NDArray pad(DType::Int64, {nc_blk, 2});
            std::fill(pad.As<std::int64_t>(), pad.As<std::int64_t>() + pad.Size(), -1);
            side_of_blocks.push_back(std::move(pad));
        }
        for (const auto& fb : new_blocks) {
            const std::size_t nf = fb.mEntryIndices.size();
            NDArray arr(DType::Int64, {nf, 2});
            std::copy(fb.mSideOf.begin(), fb.mSideOf.end(), arr.As<std::int64_t>());
            side_of_blocks.push_back(std::move(arr));
        }
        res.mMesh.AddCellData(kSideOfName, std::move(side_of_blocks));

        // Transfer `exodus:sset:<region>:<var>` field_data onto the materialized facet blocks.
        const std::string sset_prefix = "exodus:sset:";
        std::map<std::string, std::map<std::string, const NDArray*>> var_by_region;
        for (const std::string& fn : rMesh.FieldDataNames()) {
            if (fn.rfind(sset_prefix, 0) != 0)
                continue;
            const std::string rest = fn.substr(sset_prefix.size());
            const std::size_t sep = rest.rfind(':');
            if (sep == std::string::npos || sep == 0 || sep + 1 >= rest.size())
                continue;
            const std::string rname = rest.substr(0, sep);
            const std::string vname = rest.substr(sep + 1);
            var_by_region[vname][rname] = &rMesh.FieldData(fn);
        }
        for (const auto& [vname, reg_map] : var_by_region) {
            std::vector<NDArray> vblocks;
            vblocks.reserve(total_blocks);
            for (std::size_t b = 0; b < orig_blocks; ++b) {
                const std::size_t ncells = res.mMesh.Cells(b).NumCells();
                NDArray pad(DType::Float64, {ncells});
                std::fill(pad.As<double>(), pad.As<double>() + ncells,
                          std::numeric_limits<double>::quiet_NaN());
                vblocks.push_back(std::move(pad));
            }
            for (const auto& fb : new_blocks) {
                const std::size_t nf = fb.mEntryIndices.size();
                NDArray arr(DType::Float64, {nf});
                double* dst = arr.As<double>();
                auto it = reg_map.find(fb.mRegionName);
                if (it == reg_map.end()) {
                    std::fill(dst, dst + nf, std::numeric_limits<double>::quiet_NaN());
                } else {
                    const detail::DoubleView src(*it->second);
                    const std::size_t src_n = it->second->Size();
                    for (std::size_t i = 0; i < nf; ++i) {
                        const std::size_t ei = fb.mEntryIndices[i];
                        dst[i] = ei < src_n ? src[ei] : std::numeric_limits<double>::quiet_NaN();
                    }
                }
                vblocks.push_back(std::move(arr));
            }
            res.mMesh.AddCellData(vname, std::move(vblocks));
        }
    }

    return res;
}

CellsToSidesResult cells_to_sides(const Mesh& rMesh, const CellsToSidesOptions& rOptions) {
    CellsToSidesResult res;
    if (!rMesh.HasCellData(kSideOfName)) {
        res.mMesh = side_clone_base(rMesh);
        for (const std::string& name : rMesh.CellDataNames()) {
            std::vector<NDArray> blks;
            for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
                blks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
            res.mMesh.AddCellData(name, std::move(blks));
        }
        return res;
    }

    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<std::int64_t> bases(nblocks + 1, 0);
    for (std::size_t b = 0; b < nblocks; ++b)
        bases[b + 1] = bases[b] + static_cast<std::int64_t>(rMesh.Cells(b).NumCells());

    // Identify which blocks are materialized side blocks (all rows have side_of >= 0).
    std::vector<bool> is_side_block(nblocks, false);
    for (std::size_t b = 0; b < nblocks && b < rMesh.CellDataNumBlocks(kSideOfName); ++b) {
        const NDArray& so = rMesh.CellData(kSideOfName, b);
        if (so.Shape().size() == 2 && so.Shape()[1] == 2 && so.Shape()[0] > 0) {
            const detail::Int64View v(so);
            if (v[0] >= 0 && v[1] >= 0)
                is_side_block[b] = true;
        }
    }

    // Map each side block to its region name (from a matching Cell region, or option, or default).
    std::vector<std::string> block_region(nblocks);
    for (std::size_t b = 0; b < nblocks; ++b) {
        if (!is_side_block[b])
            continue;
        std::string found;
        const std::int64_t bstart = bases[b];
        const std::int64_t bend = bases[b + 1];
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
            const Region& reg = rMesh.Region(r);
            if (reg.mKind != RegionKind::Cell || reg.NumEntries() == 0)
                continue;
            const std::int64_t* e = reg.Entries();
            if (e[0] >= bstart && e[0] < bend) {
                found = reg.mName;
                break;
            }
        }
        if (found.empty())
            found = rOptions.mRegionName.empty() ? ("sides_" + std::to_string(b))
                                                 : rOptions.mRegionName;
        block_region[b] = found;
    }

    // Build output mesh (optionally dropping facet blocks).
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    std::vector<std::size_t> kept_blocks;
    for (std::size_t b = 0; b < nblocks; ++b) {
        if (rOptions.mDropFacetBlocks && is_side_block[b])
            continue;
        kept_blocks.push_back(b);
        const auto cb = rMesh.Cells(b);
        if (!cb.IsRagged())
            out.AddCellBlock(cb.Type(), detail::data_owned_copy(cb.Conn()));
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (rOptions.mDropFacetBlocks && name == kSideOfName)
            continue;
        std::vector<NDArray> blks;
        for (std::size_t b : kept_blocks)
            if (b < rMesh.CellDataNumBlocks(name))
                blks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
        if (!blks.empty())
            out.AddCellData(name, std::move(blks));
    }

    // Copy non-dropped regions.
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (rOptions.mDropFacetBlocks && reg.mKind == RegionKind::Cell && reg.NumEntries() > 0) {
            bool drop = false;
            for (std::size_t b = 0; b < nblocks; ++b)
                if (is_side_block[b] && reg.Entries()[0] >= bases[b] && reg.Entries()[0] < bases[b + 1])
                    drop = true;
            if (drop)
                continue;
        }
        out.AddRegion(reg);
    }

    // Emit Side regions and `exodus:sset:<region>:<var>` arrays from side blocks.
    std::map<std::string, std::vector<std::pair<std::int64_t, std::int64_t>>> side_entries;
    std::map<std::string, std::map<std::string, std::vector<double>>> side_vars;
    for (std::size_t b = 0; b < nblocks; ++b) {
        if (!is_side_block[b])
            continue;
        const std::string& rname = block_region[b];
        const NDArray& so = rMesh.CellData(kSideOfName, b);
        const detail::Int64View v(so);
        const std::size_t nf = so.Shape()[0];
        auto& pairs = side_entries[rname];
        for (std::size_t i = 0; i < nf; ++i) {
            pairs.emplace_back(v[2 * i], v[2 * i + 1]);
            ++res.mNumSides;
        }
        for (const std::string& cname : rMesh.CellDataNames()) {
            if (cname == kSideOfName || b >= rMesh.CellDataNumBlocks(cname))
                continue;
            const NDArray& arr = rMesh.CellData(cname, b);
            if (data_num_components(arr) != 1 || arr.Size() != nf)
                continue;
            const detail::DoubleView dv(arr);
            bool any_finite = false;
            for (std::size_t i = 0; i < nf; ++i)
                if (std::isfinite(dv[i]))
                    any_finite = true;
            if (!any_finite)
                continue;
            auto& dst = side_vars[rname][cname];
            for (std::size_t i = 0; i < nf; ++i)
                dst.push_back(dv[i]);
        }
    }

    for (const auto& [rname, pairs] : side_entries) {
        // Sort pairs in canonical order, keeping associated side_vars aligned.
        std::vector<std::size_t> order(pairs.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return pairs[a] < pairs[b]; });
        NDArray arr(DType::Int64, {order.size(), 2});
        auto* p = arr.As<std::int64_t>();
        for (std::size_t i = 0; i < order.size(); ++i) {
            p[2 * i] = pairs[order[i]].first;
            p[2 * i + 1] = pairs[order[i]].second;
        }
        out.AddRegion(Region(rname, RegionKind::Side, std::move(arr)));

        auto vit = side_vars.find(rname);
        if (vit != side_vars.end()) {
            for (const auto& [vname, vals] : vit->second) {
                NDArray fv(DType::Float64, {order.size()});
                double* fp = fv.As<double>();
                for (std::size_t i = 0; i < order.size(); ++i)
                    fp[i] = order[i] < vals.size() ? vals[order[i]] : 0.0;
                out.AddFieldData("exodus:sset:" + rname + ":" + vname, std::move(fv));
            }
        }
    }

    res.mMesh = std::move(out);
    return res;
}

}  // namespace meshioplusplus
