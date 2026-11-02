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
// The adjacency operation (roadmap 7.2.2): the node, clique and cell graphs in
// CSR form, and the METIS and Scotch graph files.

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/adjacency.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

std::vector<std::int64_t> adj_vec(const NDArray& rA) {
    const std::int64_t* p = rA.As<std::int64_t>();
    return std::vector<std::int64_t>(p, p + rA.Size());
}

std::vector<std::int64_t> adj_row(const Adjacency& rG, std::size_t Row) {
    const auto ptr = adj_vec(rG.mIndptr);
    const auto idx = adj_vec(rG.mIndices);
    return std::vector<std::int64_t>(idx.begin() + ptr[Row], idx.begin() + ptr[Row + 1]);
}

void adj_expect_simple_symmetric(const Adjacency& rG) {
    const auto ptr = adj_vec(rG.mIndptr);
    const auto idx = adj_vec(rG.mIndices);
    ASSERT_FALSE(ptr.empty());
    EXPECT_EQ(ptr.front(), 0);
    EXPECT_EQ(static_cast<std::size_t>(ptr.back()), idx.size());
    std::set<std::pair<std::int64_t, std::int64_t>> pairs;
    for (std::size_t i = 0; i + 1 < ptr.size(); ++i)
        for (std::int64_t k = ptr[i]; k < ptr[i + 1]; ++k) {
            EXPECT_NE(idx[k], static_cast<std::int64_t>(i));
            if (k > ptr[i])
                EXPECT_LT(idx[k - 1], idx[k]);  // sorted, no duplicates
            pairs.emplace(static_cast<std::int64_t>(i), idx[k]);
        }
    for (const auto& [a, b] : pairs)
        EXPECT_TRUE(pairs.count({b, a})) << a << "->" << b << " has no reverse";
}

std::string adj_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// A 2x2x2 grid of hexahedra split into 6 tetrahedra each would be heavy; one
// unit cube split into 5 tetrahedra is enough: 8 nodes, 18 edges (12 cube edges
// plus 6 face diagonals), the central tetrahedron facing the other four.
Mesh adj_five_tets() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "tetra", {{0, 1, 3, 4}, {1, 2, 3, 6}, {1, 4, 5, 6}, {3, 4, 6, 7}, {1, 3, 4, 6}});
}

}  // namespace

TEST(Adjacency, NodeGraphOfFiveTetCube) {
    const Adjacency g = adjacency(adj_five_tets(), AdjacencyKind::Node);
    adj_expect_simple_symmetric(g);
    EXPECT_EQ(g.mIndptr.Size(), 9u);
    // The probe: indptr[-1] is twice the edge count.
    EXPECT_EQ(adj_vec(g.mIndptr).back(), 2 * 18);
    EXPECT_EQ(adj_row(g, 0), (std::vector<std::int64_t>{1, 3, 4}));
    EXPECT_EQ(adj_row(g, 1), (std::vector<std::int64_t>{0, 2, 3, 4, 5, 6}));
}

TEST(Adjacency, CliqueContainsNodeGraph) {
    const Mesh m = mt::hex_mesh();
    const Adjacency node = adjacency(m, AdjacencyKind::Node);
    const Adjacency clique = adjacency(m, AdjacencyKind::Clique);
    adj_expect_simple_symmetric(clique);
    // A hexahedron's 8 nodes all share the cell: the complete graph K8.
    EXPECT_EQ(adj_vec(clique.mIndptr).back(), 8 * 7);
    // Its edge graph is the cube's 12 edges, every one of which is in K8.
    EXPECT_EQ(adj_vec(node.mIndptr).back(), 2 * 12);
    EXPECT_EQ(adj_row(node, 0), (std::vector<std::int64_t>{1, 3, 4}));
}

TEST(Adjacency, CellGraphSharesFaces) {
    const Adjacency g = adjacency(adj_five_tets(), AdjacencyKind::Cell);
    adj_expect_simple_symmetric(g);
    EXPECT_EQ(g.mIndptr.Size(), 6u);
    // The central tetrahedron (cell 4) shares a face with each corner one,
    // which touch each other only along edges.
    EXPECT_EQ(adj_row(g, 4), (std::vector<std::int64_t>{0, 1, 2, 3}));
    for (std::size_t c = 0; c < 4; ++c)
        EXPECT_EQ(adj_row(g, c), (std::vector<std::int64_t>{4}));
}

TEST(Adjacency, CellGraphJoinsPolyhedraAndHexahedra) {
    // A hexahedron and a polyhedral cube stacked on it share one face.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {0, 0, 2},
                                    {1, 0, 2},
                                    {1, 1, 2},
                                    {0, 1, 2}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    m.AddPolyhedronBlock("polyhedron8", {{{
                                            {4, 7, 6, 5},
                                            {8, 9, 10, 11},
                                            {4, 5, 9, 8},
                                            {6, 7, 11, 10},
                                            {5, 6, 10, 9},
                                            {4, 8, 11, 7},
                                        }}});
    const Adjacency cells = adjacency(m, AdjacencyKind::Cell);
    EXPECT_EQ(adj_vec(cells.mIndptr), (std::vector<std::int64_t>{0, 1, 2}));
    EXPECT_EQ(adj_vec(cells.mIndices), (std::vector<std::int64_t>{1, 0}));
    // The polyhedron's face rings are its edges: 12 + 12 - 4 shared = 20.
    const Adjacency nodes = adjacency(m, AdjacencyKind::Node);
    adj_expect_simple_symmetric(nodes);
    EXPECT_EQ(adj_vec(nodes.mIndptr).back(), 2 * 20);
}

TEST(Adjacency, PolygonBlocksRaggedAndUniform) {
    // Two pentagons sharing the edge 1-2, once as a ragged block and once as a
    // uniform "polygon5" block: the same edge and dual graphs.
    const std::vector<std::vector<double>> pts = {{0, 0, 0},     {1, 0, 0},      {1, 1, 0},
                                                  {0, 1, 0},     {-0.5, 0.5, 0}, {2, 0, 0},
                                                  {2.5, 0.5, 0}, {2, 1, 0}};
    const std::vector<std::vector<std::int64_t>> rows = {{0, 1, 2, 3, 4}, {1, 5, 6, 7, 2}};
    Mesh ragged;
    ragged.AssignPoints(mt::points_from(pts));
    ragged.AddPolygonBlock("polygon", rows);
    const Mesh uniform = mt::make_mesh(pts, "polygon5", rows);
    for (const Mesh* m : std::array<const Mesh*, 2>{&ragged, &uniform}) {
        const Adjacency nodes = adjacency(*m, AdjacencyKind::Node);
        adj_expect_simple_symmetric(nodes);
        EXPECT_EQ(adj_vec(nodes.mIndptr).back(), 2 * 9);
        EXPECT_EQ(adj_row(nodes, 1), (std::vector<std::int64_t>{0, 2, 5}));
        const Adjacency cells = adjacency(*m, AdjacencyKind::Cell);
        EXPECT_EQ(adj_vec(cells.mIndices), (std::vector<std::int64_t>{1, 0}));
    }
}

TEST(Adjacency, MixedDimensionsIsolateLowerCells) {
    // A triangle pair plus a line block: in the 2-D dual graph the lines are
    // isolated vertices, while the node graph joins every block.
    Mesh m = mt::tri_mesh();
    m.AddCellBlock("line", mt::conn_from({{0, 1}, {1, 2}}));
    const Adjacency cells = adjacency(m, AdjacencyKind::Cell);
    EXPECT_EQ(adj_vec(cells.mIndptr), (std::vector<std::int64_t>{0, 1, 2, 2, 2}));
    const Adjacency nodes = adjacency(m, AdjacencyKind::Node);
    EXPECT_EQ(adj_vec(nodes.mIndptr).back(), 2 * 5);
}

TEST(Adjacency, KindAndFormatNames) {
    EXPECT_EQ(adjacency_kind_from_name("node"), AdjacencyKind::Node);
    EXPECT_EQ(adjacency_kind_from_name("clique"), AdjacencyKind::Clique);
    EXPECT_EQ(adjacency_kind_from_name("cell"), AdjacencyKind::Cell);
    EXPECT_EQ(adjacency_kind_name(AdjacencyKind::Clique), "clique");
    EXPECT_THROW(adjacency_kind_from_name("edge"), std::invalid_argument);
    EXPECT_EQ(graph_format_from_name("scotch"), GraphFormat::Scotch);
    EXPECT_THROW(graph_format_from_name("chaco"), std::invalid_argument);
    EXPECT_EQ(graph_format_for_path("a/b.graph"), GraphFormat::Metis);
    EXPECT_EQ(graph_format_for_path("b.METIS"), GraphFormat::Metis);
    EXPECT_EQ(graph_format_for_path("b.grf"), GraphFormat::Scotch);
    EXPECT_THROW(graph_format_for_path("b.txt"), std::invalid_argument);
}

TEST(Adjacency, WritesMetisAndScotch) {
    // Two triangles sharing the diagonal 0-2: 4 nodes, 5 edges.
    const Adjacency g = adjacency(mt::tri_mesh(), AdjacencyKind::Node);
    const std::string metis = mt::temp_path(".graph");
    const GraphSummary s = write_graph(metis, g, GraphFormat::Metis);
    EXPECT_EQ(s.mNumVertices, 4);
    EXPECT_EQ(s.mNumEdges, 5);
    EXPECT_EQ(s.mMaxDegree, 3);
    EXPECT_EQ(adj_slurp(metis), "4 5\n2 3 4\n1 3\n1 2 4\n1 3\n");
    const std::string scotch = mt::temp_path(".grf");
    write_graph(scotch, g, GraphFormat::Scotch);
    EXPECT_EQ(adj_slurp(scotch), "0\n4 10\n0 000\n3 1 2 3\n2 0 2\n3 0 1 3\n2 0 2\n");
    // An isolated vertex is an empty METIS line, degree 0 in Scotch.
    Mesh lonely = mt::tri_mesh();
    NDArray pts = NDArray::Uninit(DType::Float64, {5, 3});
    for (std::size_t k = 0; k < 12; ++k)
        pts.As<double>()[k] = lonely.Points().As<double>()[k];
    for (std::size_t k = 12; k < 15; ++k)
        pts.As<double>()[k] = 5.0;
    lonely.AssignPoints(std::move(pts));
    write_graph(metis, adjacency(lonely), GraphFormat::Metis);
    EXPECT_EQ(adj_slurp(metis), "5 5\n2 3 4\n1 3\n1 2 4\n1 3\n\n");
    write_graph(scotch, adjacency(lonely), GraphFormat::Scotch);
    EXPECT_EQ(adj_slurp(scotch), "0\n5 10\n0 000\n3 1 2 3\n2 0 2\n3 0 1 3\n2 0 2\n0\n");
    std::filesystem::remove(metis);
    std::filesystem::remove(scotch);
}

TEST(Adjacency, WriterRefusesMalformedGraphs) {
    const std::string path = mt::temp_path(".graph");
    const auto csr = [](std::vector<std::int64_t> Ptr, std::vector<std::int64_t> Idx) {
        Adjacency g;
        g.mIndptr = NDArray::Uninit(DType::Int64, {Ptr.size()});
        std::copy(Ptr.begin(), Ptr.end(), g.mIndptr.As<std::int64_t>());
        g.mIndices = NDArray::Uninit(DType::Int64, {Idx.size()});
        std::copy(Idx.begin(), Idx.end(), g.mIndices.As<std::int64_t>());
        return g;
    };
    EXPECT_THROW(write_graph(path, csr({0, 1, 1}, {1}), GraphFormat::Metis),
                 std::invalid_argument);  // asymmetric
    EXPECT_THROW(write_graph(path, csr({0, 1}, {0}), GraphFormat::Metis),
                 std::invalid_argument);  // self-loop
    EXPECT_THROW(write_graph(path, csr({0, 1, 2}, {1, 5}), GraphFormat::Metis),
                 std::invalid_argument);  // out of range
    EXPECT_THROW(write_graph(path, csr({0, 2, 1}, {1, 1}), GraphFormat::Metis),
                 std::invalid_argument);  // decreasing offsets
    EXPECT_THROW(write_graph(path, csr({}, {}), GraphFormat::Metis), std::invalid_argument);
    std::filesystem::remove(path);
}
