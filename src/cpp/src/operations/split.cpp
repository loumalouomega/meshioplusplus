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
#include "meshioplusplus/region.hpp"

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
    for (std::size_t i = 0; i < total; ++i)
        parent[i] = static_cast<std::int64_t>(i);

    // Union cells that share a node (first cell per node as the anchor).
    std::vector<std::int64_t> node_first(n, -1);
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t gc = static_cast<std::int64_t>(block_base[bi] + c);
            split_visit_nodes(cb, c, [&](std::int64_t node) {
                std::int64_t& first = node_first[static_cast<std::size_t>(node)];
                if (first < 0)
                    first = gc;
                else
                    split_uf_union(parent, first, gc);
            });
        }
        ++bi;
    }

    // Map component roots to piece indices in first-seen (ascending gc) order.
    std::vector<SplitGroup> groups;
    std::unordered_map<std::int64_t, std::size_t> root_index;
    std::vector<std::size_t> gc_to_group(total, 0);
    for (std::size_t gc = 0; gc < total; ++gc) {
        const std::int64_t root = split_uf_find(parent, static_cast<std::int64_t>(gc));
        auto it = root_index.find(root);
        std::size_t gi;
        if (it == root_index.end()) {
            gi = groups.size();
            root_index[root] = gi;
            groups.push_back({});
            groups.back().key = std::to_string(gi);
            split_init_group(groups.back(), nblocks);
        } else {
            gi = it->second;
        }
        gc_to_group[gc] = gi;
    }

    for (std::size_t b = 0; b < nblocks; ++b)
        for (std::size_t c = 0; c < block_base[b + 1] - block_base[b]; ++c)
            groups[gc_to_group[block_base[b] + c]].kept[b].push_back(static_cast<std::int64_t>(c));

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
    res.mPieces.reserve(groups.size());
    for (SplitGroup& g : groups) {
        detail::SubsetResult sub =
            detail::build_cell_subset(rMesh, g.kept, "", "", /*drop_empty_blocks=*/true, "split");
        SplitPiece piece;
        piece.mKey = std::move(g.key);
        piece.mMesh = std::move(sub.mMesh);
        piece.mPointMap = std::move(sub.mPointMap);
        piece.mCellMaps = std::move(sub.mCellMaps);
        res.mPieces.push_back(std::move(piece));
    }
    return res;
}

}  // namespace meshioplusplus
