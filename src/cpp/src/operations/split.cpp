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
// Partition a mesh into submeshes by cell type, connected component (union-find
// over node-sharing cell adjacency), or integer cell_data tag. Each piece is
// pruned via the shared detail::build_cell_subset. Built entirely through the
// uniform mesh API so it compiles under every mesh backend. See
// operations/split.hpp for the contract.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// A named list of per-block kept-cell lists (one output piece, pre-prune).
struct SplitGroup {
    std::string key;
    std::vector<std::vector<std::int64_t>> kept;  // per input block
};

// Ensure every group has one (possibly empty) kept list per input block.
void split_init_group(SplitGroup& rG, std::size_t nblocks) {
    rG.kept.assign(nblocks, {});
}

// Union-find over global cell indices (flood-fill for connected components).
std::int64_t split_uf_find(std::vector<std::int64_t>& rParent, std::int64_t x) {
    while (rParent[static_cast<std::size_t>(x)] != x) {
        rParent[static_cast<std::size_t>(x)] =
            rParent[static_cast<std::size_t>(rParent[static_cast<std::size_t>(x)])];
        x = rParent[static_cast<std::size_t>(x)];
    }
    return x;
}

void split_uf_union(std::vector<std::int64_t>& rParent, std::int64_t a, std::int64_t b) {
    std::int64_t ra = split_uf_find(rParent, a), rb = split_uf_find(rParent, b);
    if (ra != rb)
        rParent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
}

// Visit each node index of one cell (rectangular / polygon / polyhedron).
template <class F>
void split_visit_nodes(const Mesh::CellView& rCb, std::size_t c, F&& f) {
    if (rCb.IsPolyhedron()) {
        for (std::size_t face = 0; face < rCb.NumFaces(c); ++face) {
            auto fc = rCb.Face(c, face);
            for (std::size_t k = 0; k < fc.second; ++k)
                f(fc.first[k]);
        }
    } else if (rCb.IsRagged()) {
        for (std::size_t k = 0; k < rCb.RowSize(c); ++k)
            f(rCb.Row(c)[k]);
    } else {
        const NDArray& conn = rCb.Conn();
        const std::size_t npc = rCb.NodesPerCell();
        for (std::size_t k = 0; k < npc; ++k)
            f(detail::read_int(conn, c * npc + k));
    }
}

// The same visit over a rectangular block whose connectivity was read once
// (`pConn`, `npc` ids per cell).
template <class F>
void split_visit_rect(const std::int64_t* pConn, std::size_t npc, std::size_t c, F&& f) {
    for (std::size_t k = 0; k < npc; ++k)
        f(pConn[c * npc + k]);
}

std::vector<SplitGroup> split_by_type(const Mesh& rMesh) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<SplitGroup> groups;
    std::unordered_map<std::string, std::size_t> index;
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        auto it = index.find(type);
        std::size_t gi;
        if (it == index.end()) {
            gi = groups.size();
            index[type] = gi;
            groups.push_back({});
            groups.back().key = type;
            split_init_group(groups.back(), nblocks);
        } else {
            gi = it->second;
        }
        std::vector<std::int64_t>& dst = groups[gi].kept[b];
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            dst.push_back(static_cast<std::int64_t>(c));
        ++b;
    }
    return groups;
}

std::vector<SplitGroup> split_by_component(const Mesh& rMesh) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t n = rMesh.NumPoints();

    // Global cell numbering (blocks concatenated in order).
    std::vector<std::size_t> block_base(nblocks + 1, 0);
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        block_base[bi + 1] = block_base[bi] + cb.NumCells();
        ++bi;
    }
    const std::size_t total = block_base[nblocks];

    std::vector<std::int64_t> parent(total);
    parallel_for_bw(total, [&](std::size_t i) { parent[i] = static_cast<std::int64_t>(i); });

    // Union cells that share a node (first cell per node as the anchor). The
    // sweep is order-dependent and stays serial; a rectangular block's
    // connectivity is read once, not through a dtype switch per id.
    std::vector<std::int64_t> node_first(n, -1);
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const auto visit = [&](std::size_t c, std::int64_t node) {
            const std::int64_t gc = static_cast<std::int64_t>(block_base[bi] + c);
            std::int64_t& first = node_first[static_cast<std::size_t>(node)];
            if (first < 0)
                first = gc;
            else
                split_uf_union(parent, first, gc);
        };
        if (!cb.IsPolyhedron() && !cb.IsRagged()) {
            const detail::Int64View conn(cb.Conn());
            const std::size_t npc = cb.NodesPerCell();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                split_visit_rect(conn.Data(), npc, c, [&](std::int64_t node) { visit(c, node); });
        } else {
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                split_visit_nodes(cb, c, [&](std::int64_t node) { visit(c, node); });
        }
        ++bi;
    }

    // Every root is its component's smallest cell (a union points the larger
    // root at the smaller) and parent[x] <= x, so one ascending pass resolves
    // each cell's root. Pieces are numbered in first-seen (ascending gc)
    // order: a component's first cell is its root, so the piece of a cell is
    // the number of roots before its own -- a scan of the root flags.
    for (std::size_t gc = 0; gc < total; ++gc)
        parent[gc] = parent[static_cast<std::size_t>(parent[gc])];
    std::vector<std::uint64_t> is_root(total);
    parallel_for_bw(total, [&](std::size_t gc) {
        is_root[gc] = parent[gc] == static_cast<std::int64_t>(gc) ? 1 : 0;
    });
    std::vector<std::uint64_t> piece_at(total);
    const std::uint64_t npieces =
        parallel_exclusive_scan(is_root.data(), total, piece_at.data(), std::uint64_t{0});

    std::vector<SplitGroup> groups(npieces);
    parallel_for(npieces, [&](std::size_t gi) {
        groups[gi].key = std::to_string(gi);
        split_init_group(groups[gi], nblocks);
    });
    for (std::size_t b = 0; b < nblocks; ++b)
        for (std::size_t c = 0; c < block_base[b + 1] - block_base[b]; ++c) {
            const std::size_t gc = block_base[b] + c;
            groups[piece_at[static_cast<std::size_t>(parent[gc])]].kept[b].push_back(
                static_cast<std::int64_t>(c));
        }

    return groups;
}

// One piece per named Cell region. Unlike every other criterion, this is not a
// partition: a cell belonging to two regions appears in two pieces, and a cell
// in none appears in none -- each named region is processed independently, so
// "split by region" means exactly what it says rather than reinterpreting
// overlapping groups as a forced partition.
//
// Point/Side regions produce no piece: split's contract is whole submeshes,
// and there is no sound way to turn "these facets" or "these points" into one
// without first deciding which surrounding cells to pull in -- a policy
// decision this operation does not make silently.
std::vector<SplitGroup> split_by_region(const Mesh& rMesh) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<SplitGroup> groups;
    for (const std::string& name : rMesh.RegionNames()) {
        const std::size_t idx = rMesh.FindRegion(name, RegionKind::Cell);
        if (idx == Mesh::npos)
            continue;
        const Region& region = rMesh.Region(idx);
        SplitGroup group;
        group.key = name;
        split_init_group(group, nblocks);
        const std::int64_t* entries = region.Entries();
        const std::size_t n = region.NumEntries();
        for (std::size_t e = 0; e < n; ++e) {
            const auto [b, row] = detail::global_to_block_row(bases, entries[e]);
            if (b != static_cast<std::size_t>(-1))
                group.kept[b].push_back(row);
        }
        groups.push_back(std::move(group));
    }
    return groups;
}

bool split_is_integer(DType dt) {
    switch (dt) {
        case DType::Int8:
        case DType::Int16:
        case DType::Int32:
        case DType::Int64:
        case DType::UInt8:
        case DType::UInt16:
        case DType::UInt32:
        case DType::UInt64:
            return true;
        default:
            return false;
    }
}

// The integer cell_data name to split on: `rTagName` if usable, else the first
// integer cell_data covering every block. Throws if none is found.
std::string split_resolve_tag(const Mesh& rMesh, const std::string& rTagName) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    auto usable = [&](const std::string& name) {
        if (!rMesh.HasCellData(name) || rMesh.CellDataNumBlocks(name) != nblocks)
            return false;
        std::size_t b = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& a = rMesh.CellData(name, b);
            const std::size_t rows = a.Shape().empty() ? 0 : a.Shape()[0];
            if (!split_is_integer(a.Dtype()) || rows != cb.NumCells())
                return false;
            ++b;
        }
        return true;
    };
    if (!rTagName.empty()) {
        if (!usable(rTagName))
            throw std::invalid_argument("split: cell_data '" + rTagName +
                                        "' is not an integer per-cell tag on every block");
        return rTagName;
    }
    for (const std::string& name : rMesh.CellDataNames())
        if (usable(name))
            return name;
    throw std::invalid_argument("split: no integer cell_data tag found (pass an explicit name)");
}

std::vector<SplitGroup> split_by_tag(const Mesh& rMesh, const std::string& rTagName) {
    const std::string name = split_resolve_tag(rMesh, rTagName);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<SplitGroup> groups;
    std::unordered_map<std::int64_t, std::size_t> index;
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& tag = rMesh.CellData(name, b);
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t v = detail::read_int(tag, c);
            auto it = index.find(v);
            std::size_t gi;
            if (it == index.end()) {
                gi = groups.size();
                index[v] = gi;
                groups.push_back({});
                groups.back().key = std::to_string(v);
                split_init_group(groups.back(), nblocks);
            } else {
                gi = it->second;
            }
            groups[gi].kept[b].push_back(static_cast<std::int64_t>(c));
        }
        ++b;
    }
    return groups;
}

}  // namespace

SplitBy split_by_from_name(const std::string& rName) {
    if (rName == "type")
        return SplitBy::Type;
    if (rName == "component")
        return SplitBy::Component;
    if (rName == "region" || rName == "tag")
        return SplitBy::Tag;
    if (rName == "regions")
        return SplitBy::Region;
    throw std::invalid_argument("split: unknown criterion '" + rName +
                                "' (expected 'type', 'region'/'tag', 'regions', or 'component')");
}

SplitResult split(const Mesh& rMesh, SplitBy by, const std::string& rTagName) {
    std::vector<SplitGroup> groups;
    switch (by) {
        case SplitBy::Type:
            groups = split_by_type(rMesh);
            break;
        case SplitBy::Component:
            groups = split_by_component(rMesh);
            break;
        case SplitBy::Tag:
            groups = split_by_tag(rMesh, rTagName);
            break;
        case SplitBy::Region:
            groups = split_by_region(rMesh);
            break;
    }

    SplitResult res;
    res.mPieces.resize(groups.size());
    const auto extract = [&](std::size_t i) {
        detail::SubsetResult sub = detail::build_cell_subset(rMesh, groups[i].kept, "", "",
                                                             /*drop_empty_blocks=*/true, "split");
        SplitPiece& piece = res.mPieces[i];
        piece.mKey = std::move(groups[i].key);
        piece.mMesh = std::move(sub.mMesh);
        piece.mPointMap = std::move(sub.mPointMap);
        piece.mCellMaps = std::move(sub.mCellMaps);
    };
    // Serial on purpose: every piece reads the input mesh's data names and
    // arrays, and the backends fill caches lazily behind const accessors (the
    // MESHIO name lists, the NATIVE global CSR, the KRATOS ModelPart), so one
    // Mesh is not safe to read from several threads through them.
    for (std::size_t i = 0; i < groups.size(); ++i)
        extract(i);
    return res;
}

}  // namespace meshioplusplus
