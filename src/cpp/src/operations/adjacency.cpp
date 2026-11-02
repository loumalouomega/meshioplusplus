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
// The node and cell graphs of a mesh in CSR form, and the METIS/Scotch graph
// files written from them. See operations/adjacency.hpp for the contract.

// System includes
#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/adjacency.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/dual_graph.hpp"

namespace meshioplusplus {

namespace {

NDArray adjgraph_array(const std::vector<std::int64_t>& rValues) {
    NDArray a = NDArray::Uninit(DType::Int64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), a.As<std::int64_t>());
    return a;
}

// A uniform (rectangular) polygon block: `cell_refine_edges` has no row for it,
// so the shared node-adjacency kernel skips it under Kind::Edge (smooth pins
// its nodes on purpose). Its rows are rings, which is their true edge set.
bool adjgraph_uniform_polygon(const Mesh::CellView& rBlock) {
    return !rBlock.IsRagged() && !rBlock.IsPolyhedron() && rBlock.Type().rfind("polygon", 0) == 0;
}

// Adds the ring edges of the uniform polygon blocks to an edge graph: a row's
// neighbours become the sorted union of its own and the rings'.
void adjgraph_add_polygon_rings(const Mesh& rMesh, std::size_t NumPoints,
                                detail::NodeAdjacency& rCsr) {
    std::vector<std::size_t> start(NumPoints + 1, 0);
    const auto visit = [&](const auto& rEmit) {
        for (const auto cb : rMesh.CellRange()) {
            if (!adjgraph_uniform_polygon(cb))
                continue;
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            if (npc < 2)
                continue;
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                for (std::size_t k = 0; k < npc; ++k) {
                    const std::int64_t u = detail::read_int(conn, c * npc + k);
                    const std::int64_t v = detail::read_int(conn, c * npc + (k + 1) % npc);
                    if (u < 0 || v < 0 || u == v || static_cast<std::size_t>(u) >= NumPoints ||
                        static_cast<std::size_t>(v) >= NumPoints)
                        continue;
                    rEmit(u, v);
                    rEmit(v, u);
                }
            }
        }
    };
    visit([&](std::int64_t u, std::int64_t) { ++start[static_cast<std::size_t>(u) + 1]; });
    if (std::all_of(start.begin(), start.end(), [](std::size_t s) { return s == 0; }))
        return;  // no uniform polygon block: the kernel's graph is complete
    for (std::size_t i = 0; i < NumPoints; ++i)
        start[i + 1] += start[i];
    std::vector<std::int64_t> extra(start[NumPoints]);
    {
        std::vector<std::size_t> cursor(start.begin(), start.end() - 1);
        visit([&](std::int64_t u, std::int64_t v) {
            extra[cursor[static_cast<std::size_t>(u)]++] = v;
        });
    }
    // Per row: sort the ring entries, then the union with the kernel's row.
    std::vector<std::vector<std::int64_t>> rows(NumPoints);
    parallel_for(NumPoints, [&](std::size_t i) {
        auto lo = extra.begin() + static_cast<std::ptrdiff_t>(start[i]);
        auto hi = extra.begin() + static_cast<std::ptrdiff_t>(start[i + 1]);
        std::sort(lo, hi);
        hi = std::unique(lo, hi);
        const auto base_lo = rCsr.mAdj.begin() + rCsr.mXadj[i];
        const auto base_hi = rCsr.mAdj.begin() + rCsr.mXadj[i + 1];
        std::set_union(base_lo, base_hi, lo, hi, std::back_inserter(rows[i]));
    });
    detail::NodeAdjacency merged;
    merged.mXadj.resize(NumPoints + 1);
    merged.mXadj[0] = 0;
    for (std::size_t i = 0; i < NumPoints; ++i)
        merged.mXadj[i + 1] = merged.mXadj[i] + static_cast<std::int64_t>(rows[i].size());
    merged.mAdj.reserve(static_cast<std::size_t>(merged.mXadj[NumPoints]));
    for (const auto& row : rows)
        merged.mAdj.insert(merged.mAdj.end(), row.begin(), row.end());
    rCsr = std::move(merged);
}

void adjgraph_put(std::string& rOut, std::int64_t Value) {
    char buf[24];
    const auto res = std::to_chars(buf, buf + sizeof(buf), Value);
    rOut.append(buf, res.ptr);
}

// Validates a CSR as an undirected simple graph; returns its summary.
GraphSummary adjgraph_check(const Adjacency& rGraph) {
    if (rGraph.mIndptr.Size() == 0)
        throw std::invalid_argument("write_graph: indptr is empty (it holds n + 1 offsets)");
    const std::size_t n = rGraph.mIndptr.Size() - 1;
    const std::size_t nnz = rGraph.mIndices.Size();
    std::vector<std::int64_t> ptr(n + 1);
    for (std::size_t i = 0; i <= n; ++i)
        ptr[i] = detail::read_int(rGraph.mIndptr, i);
    if (ptr[0] != 0 || static_cast<std::size_t>(ptr[n]) != nnz)
        throw std::invalid_argument(
            "write_graph: indptr must start at 0 and end at the number of indices (" +
            std::to_string(nnz) + ")");
    GraphSummary summary;
    summary.mNumVertices = static_cast<std::int64_t>(n);
    std::vector<std::int64_t> idx(nnz);
    for (std::size_t k = 0; k < nnz; ++k)
        idx[k] = detail::read_int(rGraph.mIndices, k);
    for (std::size_t i = 0; i < n; ++i) {
        if (ptr[i + 1] < ptr[i])
            throw std::invalid_argument("write_graph: indptr decreases at row " +
                                        std::to_string(i));
        summary.mMaxDegree = std::max(summary.mMaxDegree, ptr[i + 1] - ptr[i]);
        for (std::int64_t k = ptr[i]; k < ptr[i + 1]; ++k) {
            const std::int64_t j = idx[static_cast<std::size_t>(k)];
            if (j < 0 || static_cast<std::size_t>(j) >= n)
                throw std::invalid_argument("write_graph: row " + std::to_string(i) +
                                            " names vertex " + std::to_string(j) +
                                            ", outside the " + std::to_string(n) + " vertices");
            if (static_cast<std::size_t>(j) == i)
                throw std::invalid_argument("write_graph: row " + std::to_string(i) +
                                            " has a self-loop, which a graph file cannot hold");
        }
    }
    // Symmetry: every (i, j) has its (j, i). Rows need not be sorted, so
    // compare the sorted (row, col) and (col, row) pair lists.
    std::vector<std::pair<std::int64_t, std::int64_t>> fwd, rev;
    fwd.reserve(nnz);
    rev.reserve(nnz);
    for (std::size_t i = 0; i < n; ++i)
        for (std::int64_t k = ptr[i]; k < ptr[i + 1]; ++k) {
            fwd.emplace_back(static_cast<std::int64_t>(i), idx[static_cast<std::size_t>(k)]);
            rev.emplace_back(idx[static_cast<std::size_t>(k)], static_cast<std::int64_t>(i));
        }
    std::sort(fwd.begin(), fwd.end());
    std::sort(rev.begin(), rev.end());
    if (std::adjacent_find(fwd.begin(), fwd.end()) != fwd.end())
        throw std::invalid_argument("write_graph: a row lists a neighbour twice");
    if (fwd != rev)
        throw std::invalid_argument(
            "write_graph: the graph is not symmetric (METIS and Scotch files are undirected)");
    summary.mNumEdges = static_cast<std::int64_t>(nnz / 2);
    return summary;
}

}  // namespace

AdjacencyKind adjacency_kind_from_name(const std::string& rName) {
    if (rName == "node")
        return AdjacencyKind::Node;
    if (rName == "clique")
        return AdjacencyKind::Clique;
    if (rName == "cell")
        return AdjacencyKind::Cell;
    throw std::invalid_argument("adjacency: unknown kind '" + rName +
                                "' (expected 'node', 'clique' or 'cell')");
}

std::string adjacency_kind_name(AdjacencyKind Kind) {
    switch (Kind) {
        case AdjacencyKind::Node:
            return "node";
        case AdjacencyKind::Clique:
            return "clique";
        case AdjacencyKind::Cell:
            return "cell";
    }
    throw std::invalid_argument("adjacency: unknown kind " +
                                std::to_string(static_cast<int>(Kind)));
}

GraphFormat graph_format_from_name(const std::string& rName) {
    if (rName == "metis")
        return GraphFormat::Metis;
    if (rName == "scotch")
        return GraphFormat::Scotch;
    throw std::invalid_argument("write_graph: unknown format '" + rName +
                                "' (expected 'metis' or 'scotch')");
}

GraphFormat graph_format_for_path(const std::string& rPath) {
    std::string ext = std::filesystem::path(rPath).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(c >= 'A' && c <= 'Z' ? c + 32 : c);
    });
    if (ext == ".graph" || ext == ".metis")
        return GraphFormat::Metis;
    if (ext == ".grf")
        return GraphFormat::Scotch;
    throw std::invalid_argument("write_graph: cannot infer the format of '" + rPath +
                                "' (use .graph or .metis for METIS, .grf for Scotch, or "
                                "name the format)");
}

Adjacency adjacency(const Mesh& rMesh, AdjacencyKind Kind) {
    Adjacency out;
    if (Kind == AdjacencyKind::Cell) {
        std::size_t total = 0;
        for (const auto cb : rMesh.CellRange())
            total += cb.NumCells();
        const detail::DualGraphCsr csr =
            detail::build_facet_dual_graph(rMesh, total, /*PolygonRings=*/true);
        out.mIndptr = adjgraph_array(csr.mXadj);
        out.mIndices = adjgraph_array(csr.mAdjncy);
        return out;
    }
    const std::size_t n = rMesh.NumPoints();
    detail::NodeAdjacency csr = detail::build_node_adjacency(rMesh, n,
                                                             Kind == AdjacencyKind::Clique
                                                                 ? detail::NodeAdjacencyKind::Clique
                                                                 : detail::NodeAdjacencyKind::Edge);
    if (Kind == AdjacencyKind::Node)
        adjgraph_add_polygon_rings(rMesh, n, csr);
    out.mIndptr = adjgraph_array(csr.mXadj);
    out.mIndices = adjgraph_array(csr.mAdj);
    return out;
}

GraphSummary graph_summary(const Adjacency& rGraph) {
    GraphSummary s;
    const std::size_t rows = rGraph.mIndptr.Size();
    if (rows == 0)
        return s;
    s.mNumVertices = static_cast<std::int64_t>(rows - 1);
    for (std::size_t i = 0; i + 1 < rows; ++i)
        s.mMaxDegree = std::max(s.mMaxDegree, detail::read_int(rGraph.mIndptr, i + 1) -
                                                  detail::read_int(rGraph.mIndptr, i));
    s.mNumEdges = static_cast<std::int64_t>(rGraph.mIndices.Size() / 2);
    return s;
}

GraphSummary write_graph(const std::string& rPath, const Adjacency& rGraph, GraphFormat Format) {
    const GraphSummary summary = adjgraph_check(rGraph);
    const std::size_t n = static_cast<std::size_t>(summary.mNumVertices);
    std::string text;
    text.reserve(16 + n * 4 + rGraph.mIndices.Size() * 8);
    if (Format == GraphFormat::Metis) {
        adjgraph_put(text, summary.mNumVertices);
        text += ' ';
        adjgraph_put(text, summary.mNumEdges);
        text += '\n';
    } else {
        text += "0\n";
        adjgraph_put(text, summary.mNumVertices);
        text += ' ';
        adjgraph_put(text, static_cast<std::int64_t>(rGraph.mIndices.Size()));
        text += "\n0 000\n";
    }
    for (std::size_t i = 0; i < n; ++i) {
        const std::int64_t lo = detail::read_int(rGraph.mIndptr, i);
        const std::int64_t hi = detail::read_int(rGraph.mIndptr, i + 1);
        bool first = true;
        if (Format == GraphFormat::Scotch) {
            adjgraph_put(text, hi - lo);
            first = false;
        }
        for (std::int64_t k = lo; k < hi; ++k) {
            if (!first)
                text += ' ';
            first = false;
            const std::int64_t j = detail::read_int(rGraph.mIndices, static_cast<std::size_t>(k));
            adjgraph_put(text, Format == GraphFormat::Metis ? j + 1 : j);
        }
        text += '\n';
    }
    auto out = detail::make_classic_ofstream(std::filesystem::path(rPath),
                                             std::ios_base::out | std::ios_base::binary);
    if (!out)
        throw std::runtime_error("write_graph: cannot open '" + rPath + "' for writing");
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!out)
        throw std::runtime_error("write_graph: failed writing '" + rPath + "'");
    return summary;
}

}  // namespace meshioplusplus
