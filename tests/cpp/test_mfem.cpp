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
/**
 * @file test_mfem.cpp
 * @brief MFEM `.mesh`/`.gf` reader/writer: a hand-written order-2 mesh whose
 *        nodes follow MFEM's numbering, grid functions, the `.mesh` hand-over in
 *        the registry, a round trip, and the error paths.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
namespace detail = meshioplusplus::detail;

std::string write_file(const std::string& rBody, const std::string& rSuffix = ".mesh") {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

// Two unit squares side by side, order-2 nodes. MFEM numbers the dofs: the six
// vertices, then the seven edges in order of first appearance -- (0,1) (1,4)
// (4,3) (3,0) of the first quad, then (1,2) (2,5) (5,4) of the second -- then
// the two interiors. Each dof is placed at its entity's centre, the first
// quad's (0,1) edge bowed down to y = -0.25.
const char* kTwoQuads =
    "MFEM mesh v1.0\n"
    "# a comment\n"
    "dimension\n2\n\n"
    "elements\n2\n3 3 0 1 4 3\n7 3 1 2 5 4\n\n"
    "boundary\n2\n11 1 0 1\n12 1 2 5\n\n"
    "vertices\n6\n\n"
    "nodes\nFiniteElementSpace\nFiniteElementCollection: H1_2D_P2\nVDim: 2\nOrdering: 1\n\n"
    "0 0\n1 0\n2 0\n0 1\n1 1\n2 1\n"
    "0.5 -0.25\n1 0.5\n0.5 1\n0 0.5\n1.5 0\n2 0.5\n1.5 1\n"
    "0.5 0.5\n1.5 0.5\n";

double coord(const Mesh& rMesh, std::int64_t P, std::size_t C) {
    return detail::read_double(rMesh.Points(), static_cast<std::size_t>(P) * 2 + C);
}

}  // namespace

TEST(Mfem, ReadsOrderTwoNodesInMfemNumbering) {
    const Mesh mesh = meshioplusplus::read_mfem(write_file(kTwoQuads));
    ASSERT_EQ(mesh.NumPoints(), 15u);
    ASSERT_EQ(mesh.PointDim(), 2u);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    ASSERT_EQ(mesh.Cells(0).Type(), "quad9");
    ASSERT_EQ(mesh.Cells(1).Type(), "line3");
    const auto& q = mesh.Cells(0).Conn();
    // First quad: corners 0 1 4 3, edges (0,1)=6 (1,4)=7 (4,3)=8 (3,0)=9, centre 13.
    const std::int64_t want0[9] = {0, 1, 4, 3, 6, 7, 8, 9, 13};
    const std::int64_t want1[9] = {1, 2, 5, 4, 10, 11, 12, 7, 14};
    for (std::size_t j = 0; j < 9; ++j) {
        EXPECT_EQ(detail::read_int(q, j), want0[j]) << j;
        EXPECT_EQ(detail::read_int(q, 9 + j), want1[j]) << j;
    }
    EXPECT_DOUBLE_EQ(coord(mesh, 6, 1), -0.25);
    // The boundary edges take the element edges' dofs.
    const auto& l = mesh.Cells(1).Conn();
    EXPECT_EQ(detail::read_int(l, 2), 6);
    EXPECT_EQ(detail::read_int(l, 5), 11);
    // Attributes: cell data and regions.
    EXPECT_EQ(detail::read_int(mesh.CellData("mfem:attribute", 0), 1), 7);
    EXPECT_NE(mesh.FindRegion("attribute_3", RegionKind::Cell), Mesh::npos);
    EXPECT_NE(mesh.FindRegion("boundary_12", RegionKind::Cell), Mesh::npos);
}

TEST(Mfem, LinearGridFunctionOnACurvedMesh) {
    const std::string mesh_path = write_file(kTwoQuads);
    const std::string gf = write_file(
        "FiniteElementSpace\nFiniteElementCollection: H1_2D_P1\nVDim: 1\nOrdering: 0\n\n"
        "0\n1\n2\n3\n4\n5\n",
        ".gf");
    const Mesh mesh = meshioplusplus::read_mfem(mesh_path, {{"t", gf}});
    const auto& t = mesh.PointData("t");
    // Vertex values, then edge and centre values interpolated linearly.
    EXPECT_DOUBLE_EQ(detail::read_double(t, 4), 4.0);
    EXPECT_DOUBLE_EQ(detail::read_double(t, 6), 0.5);   // edge (0,1)
    EXPECT_DOUBLE_EQ(detail::read_double(t, 13), 2.0);  // centre of 0 1 4 3
}

TEST(Mfem, OrderTwoGridFunctionPromotesALinearMesh) {
    const std::string mesh_path = write_file(
        "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 2 0 1 2\nboundary\n0\n"
        "vertices\n3\n2\n0 0\n1 0\n0 1\n");
    const std::string gf = write_file(
        "FiniteElementSpace\nFiniteElementCollection: H1_2D_P2\nVDim: 1\nOrdering: 0\n\n"
        "0\n1\n2\n10\n11\n12\n",
        ".gf");
    const Mesh mesh = meshioplusplus::read_mfem(mesh_path, {{"u", gf}});
    ASSERT_EQ(mesh.Cells(0).Type(), "triangle6");
    ASSERT_EQ(mesh.NumPoints(), 6u);
    EXPECT_DOUBLE_EQ(coord(mesh, 3, 0), 0.5);  // edge (0,1) midpoint
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.PointData("u"), 5), 12.0);
}

TEST(Mfem, RegistryHandsMfemMeshFilesOver) {
    const std::string mfem = write_file(kTwoQuads);
    EXPECT_EQ(meshioplusplus::resolve_format(mfem, ""), "mfem");
    const std::string medit = write_file("MeshVersionFormatted 2\nDimension 3\nEnd\n");
    EXPECT_EQ(meshioplusplus::resolve_format(medit, ""), "medit");
    EXPECT_EQ(meshioplusplus::resolve_format("/nonexistent/x.mesh", ""), "medit");
    const Mesh mesh = meshioplusplus::registry_read(mfem, meshioplusplus::resolve_format(mfem, ""),
                                                    meshioplusplus::ReadOptions{});
    EXPECT_EQ(mesh.Cells(0).Type(), "quad9");
}

TEST(Mfem, RoundTripsThroughItsOwnWriter) {
    const Mesh mesh = meshioplusplus::read_mfem(write_file(kTwoQuads));
    const std::string out = mt::temp_path(".mesh");
    meshioplusplus::write_mfem(out, mesh);
    const Mesh back = meshioplusplus::read_mfem(out);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    for (std::size_t k = 0; k < mesh.NumPoints() * 2; ++k)
        EXPECT_DOUBLE_EQ(detail::read_double(back.Points(), k),
                         detail::read_double(mesh.Points(), k));
    for (std::size_t b = 0; b < 2; ++b) {
        const auto& x = mesh.Cells(b).Conn();
        const auto& y = back.Cells(b).Conn();
        ASSERT_EQ(x.Size(), y.Size());
        for (std::size_t k = 0; k < x.Size(); ++k)
            EXPECT_EQ(detail::read_int(x, k), detail::read_int(y, k));
    }
}

TEST(Mfem, ErrorsNameTheCulprit) {
    auto expect_error = [](const std::string& rBody, const std::string& rNeedle) {
        try {
            meshioplusplus::read_mfem(write_file(rBody));
            ADD_FAILURE() << "no error for " << rNeedle;
        } catch (const ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(rNeedle), std::string::npos) << e.what();
        }
    };
    expect_error("MFEM NC mesh v1.0\n", "no dimension");
    expect_error("MFEM NC mesh v1.0\ndimension\n2\nelements\n0\nboundary\n0\n",
                 "no top-level coordinates");
    expect_error("MFEM NC mesh v1.0\ndimension\n2\nelements\n0\nnodes\n", "curved non-conforming");
    expect_error(
        "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 3 0 1 2 9\nboundary\n0\n"
        "vertices\n3\n2\n0 0\n1 0\n0 1\n",
        "out of range");
    expect_error("MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 2 0 1\n", "ends");
    expect_error("not a mesh\n", "not an MFEM mesh");
}

namespace {

// One triangle, (0,0) (2,0) (0,1), with order-3 Gauss-Lobatto nodes placed by
// the affine map F(x, y) = (2x, y): edges in first-met order, each's two dofs
// from its lower vertex, then the interior dof. `u` = x + 2y on the same dofs.
std::string affine_p3_triangle(std::string& rGf) {
    const double a = (1.0 - 1.0 / std::sqrt(5.0)) / 2.0;
    const double t[2] = {a, 1.0 - a};
    std::vector<std::array<double, 2>> ref = {{0, 0}, {1, 0}, {0, 1}};
    for (double v : t)
        ref.push_back({v, 0});  // edge (0,1)
    for (double v : t)
        ref.push_back({1 - v, v});  // edge (1,2)
    for (double v : t)
        ref.push_back({0, v});  // edge (0,2), from vertex 0
    ref.push_back({1.0 / 3, 1.0 / 3});
    std::ostringstream mesh, gf;
    mesh.precision(17);
    gf.precision(17);
    mesh << "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 2 0 1 2\nboundary\n0\nvertices\n3\n\n"
            "nodes\nFiniteElementSpace\nFiniteElementCollection: H1_2D_P3\nVDim: 2\n"
            "Ordering: 1\n\n";
    gf << "FiniteElementSpace\nFiniteElementCollection: H1_2D_P3\nVDim: 1\nOrdering: 0\n\n";
    for (const auto& r : ref) {
        mesh << 2 * r[0] << " " << r[1] << "\n";
        gf << 2 * r[0] + 2 * r[1] << "\n";
    }
    rGf = gf.str();
    return mesh.str();
}

}  // namespace

TEST(Mfem, OrderThreeNodesBecomeVtkLagrangeCells) {
    std::string gf_text;
    const std::string path = write_file(affine_p3_triangle(gf_text));
    const std::string gf = write_file(gf_text, ".gf");
    const Mesh mesh = meshioplusplus::read_mfem(path, {{"u", gf}});
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "VTK_LAGRANGE_TRIANGLE");
    ASSERT_EQ(mesh.Cells(0).NodesPerCell(), 10u);
    // VTK's order: corners, edges 0-1, 1-2, 2-0, the interior; equispaced.
    const double lattice[10][2] = {{0, 0},
                                   {1, 0},
                                   {0, 1},
                                   {1.0 / 3, 0},
                                   {2.0 / 3, 0},
                                   {2.0 / 3, 1.0 / 3},
                                   {1.0 / 3, 2.0 / 3},
                                   {0, 2.0 / 3},
                                   {0, 1.0 / 3},
                                   {1.0 / 3, 1.0 / 3}};
    const meshioplusplus::NDArray& conn = mesh.Cells(0).Conn();
    const meshioplusplus::NDArray& u = mesh.PointData("u");
    for (std::size_t k = 0; k < 10; ++k) {
        const auto p = static_cast<std::size_t>(detail::read_int(conn, k));
        EXPECT_NEAR(detail::read_double(mesh.Points(), 2 * p), 2 * lattice[k][0], 1e-14) << k;
        EXPECT_NEAR(detail::read_double(mesh.Points(), 2 * p + 1), lattice[k][1], 1e-14) << k;
        EXPECT_NEAR(detail::read_double(u, p), 2 * lattice[k][0] + 2 * lattice[k][1], 1e-14) << k;
    }
    // Written back as order-3 nodes, it reads the same.
    const std::string out = mt::temp_path(".mesh");
    meshioplusplus::write_mfem(out, mesh);
    std::ifstream in(out);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("FiniteElementCollection: H1_2D_P3"), std::string::npos);
    const Mesh back = meshioplusplus::read_mfem(out);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    for (std::size_t k = 0; k < mesh.NumPoints() * 2; ++k)
        EXPECT_NEAR(detail::read_double(back.Points(), k), detail::read_double(mesh.Points(), k),
                    1e-14);
}

TEST(Mfem, NonConformingMeshReadsAsItsLeaves) {
    // A square refined once: four leaves, the edge midpoints and the centre
    // placed between their vertex parents.
    const std::string path = write_file(
        "MFEM NC mesh v1.0\ndimension\n2\nelements\n5\n0 1 3 3 1 2 3 4\n0 1 3 0 0 4 8 7\n"
        "0 1 3 0 4 1 5 8\n0 1 3 0 8 5 2 6\n0 1 3 0 7 8 6 3\nboundary\n0\n"
        "vertex_parents\n5\n4 0 1\n5 1 2\n6 2 3\n7 0 3\n8 0 2\n"
        "coordinates\n4\n2\n0 0\n2 0\n2 2\n0 2\n");
    const Mesh mesh = meshioplusplus::read_mfem(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad");
    EXPECT_EQ(mesh.Cells(0).NumCells(), 4u);
    ASSERT_EQ(mesh.NumPoints(), 9u);
    // MFEM's vertex numbers: the top-level vertices 0-3, then the others as the
    // leaves (in Hilbert order) meet them: 4, 8, 7 in the first leaf, 5, then 6.
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), 2 * 5), 1.0);  // the centre, node 8
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), 2 * 5 + 1), 1.0);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), 2 * 7), 2.0);  // node 5, midpoint of 1-2
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), 2 * 7 + 1), 1.0);
}

namespace {

// Two unit squares side by side over two ranks, each rank its own serial mesh
// (local vertices), written as <dir>/m.000000 and m.000001: with communication
// groups (ParPrint) or without (ParMesh::Save, which lists the interface edge
// as boundary on both ranks, attribute 3 on rank 0 and 4 on rank 1).
std::string two_rank_mesh(bool Groups) {
    const std::string dir = mt::temp_path("_pmesh");
    std::filesystem::create_directories(dir);
    const char* rank_body[2] = {
        "elements\n1\n1 3 0 1 3 2\n\nboundary\n%B\n\nvertices\n4\n2\n0 0\n1 0\n0 1\n1 1\n",
        "elements\n1\n1 3 0 1 2 3\n\nboundary\n%B\n\nvertices\n4\n2\n1 0\n2 0\n2 1\n1 1\n"};
    const char* boundary[2][2] = {{"1\n1 1 0 1\n", "2\n1 1 0 1\n3 1 1 3\n"},
                                  {"1\n1 1 0 1\n", "2\n1 1 0 1\n4 1 0 3\n"}};
    const char* groups[2] = {
        "mfem_serial_mesh_end\n\ncommunication_groups\nnumber_of_groups 2\n\n"
        "# number of entities in each group, followed by ranks in group\n1 0\n2 0 1\n\n"
        "total_shared_vertices 2\ntotal_shared_edges 1\n\n# group 1\nshared_vertices 2\n1\n3\n\n"
        "shared_edges 1\n1 3\n\nmfem_mesh_end\n",
        "mfem_serial_mesh_end\n\ncommunication_groups\nnumber_of_groups 2\n\n1 1\n2 0 1\n\n"
        "total_shared_vertices 2\ntotal_shared_edges 1\n\n# group 1\nshared_vertices 2\n0\n3\n\n"
        "shared_edges 1\n0 3\n\nmfem_mesh_end\n"};
    for (int r = 0; r < 2; ++r) {
        std::string body = rank_body[r];
        body.replace(body.find("%B"), 2, boundary[r][Groups ? 0 : 1]);
        std::ofstream(dir + "/m.00000" + std::to_string(r)) << "MFEM mesh v1.0\n\ndimension\n2\n\n"
                                                            << body << (Groups ? groups[r] : "");
    }
    return dir + "/m.000001";
}

}  // namespace

TEST(Mfem, ParallelRanksMergeIntoOneMesh) {
    for (const bool groups : {true, false}) {
        const std::string path = two_rank_mesh(groups);
        const Mesh mesh = meshioplusplus::read_mfem(path);
        EXPECT_EQ(mesh.NumPoints(), 6u) << "groups " << groups;
        ASSERT_EQ(mesh.NumCellBlocks(), 2u);
        EXPECT_EQ(mesh.Cells(0).Type(), "quad");
        EXPECT_EQ(mesh.Cells(0).NumCells(), 2u);
        // The two bottom edges; the interface edge the ranks both list is gone.
        EXPECT_EQ(mesh.Cells(1).NumCells(), 2u) << "groups " << groups;
        const auto& parts = mesh.CellData("partition:part", 0);
        EXPECT_EQ(detail::read_int(parts, 0), 0);
        EXPECT_EQ(detail::read_int(parts, 1), 1);
        // The shared corner (1, 1) is one point used by both quads.
        const auto& conn = mesh.Cells(0).Conn();
        std::set<std::int64_t> a, b;
        for (std::size_t k = 0; k < 4; ++k) {
            a.insert(detail::read_int(conn, k));
            b.insert(detail::read_int(conn, 4 + k));
        }
        std::vector<std::int64_t> common;
        std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(common));
        EXPECT_EQ(common.size(), 2u);
        // One rank alone.
        meshioplusplus::ReadOptions one;
        one.mPiece = 1;
        one.mPieceSet = true;
        const Mesh piece = meshioplusplus::read_mfem(path, {}, one);
        EXPECT_EQ(piece.Cells(0).NumCells(), 1u);
        EXPECT_EQ(detail::read_int(piece.CellData("partition:part", 0), 0), 1);
        one.mPiece = 2;
        EXPECT_THROW(meshioplusplus::read_mfem(path, {}, one), ReadError);
    }
}

// A quarter annulus, radii 1 and 2, as one order-2 x order-1 NURBS patch in
// MFEM's global form: the four corners, then the one interior control point
// of each curved edge (weight 1/sqrt 2) -- edges numbered as listed.
TEST(Mfem, NurbsPatchSampledAtItsLattice) {
    const std::string w = "0.70710678118654757";
    const std::string mesh_path = write_file(
        "MFEM NURBS mesh v1.0\n"
        "dimension\n2\n"
        "elements\n1\n5 3 0 1 2 3\n"
        "boundary\n0\n"
        "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\n"
        "vertices\n4\n"
        "knotvectors\n2\n2 3 0 0 0 1 1 1\n1 2 0 0 1 1\n"
        "weights\n1\n1\n1\n1\n" +
        w + "\n" + w +
        "\n"
        "FiniteElementSpace\nFiniteElementCollection: NURBS2\nVDim: 2\nOrdering: 1\n"
        "1 0\n0 1\n0 2\n2 0\n1 1\n2 2\n");
    // u = the control points' x, as a NURBS field on the same space
    const std::string gf = write_file(
        "MFEM FiniteElementSpace v1.0\nFiniteElementCollection: NURBS2\nVDim: 1\nOrdering: 0\n"
        "End: MFEM FiniteElementSpace v1.0\n\n1\n0\n0\n2\n1\n2\n",
        ".gf");
    const Mesh mesh = meshioplusplus::read_mfem(mesh_path, {{"u", gf}});
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    ASSERT_EQ(mesh.Cells(0).Type(), "quad9");
    EXPECT_EQ(mesh.Cells(1).Type(), "line3");
    EXPECT_EQ(mesh.Cells(1).NumCells(), 4u);  // the boundary MFEM builds
    EXPECT_EQ(mesh.NumPoints(), 9u);
    const auto& q = mesh.Cells(0).Conn();
    // the mid-node of edge (0,1) lies on the unit circle at 45 degrees, the
    // centre halfway out at 45 degrees; u is x there (the same rational map)
    const auto n01 = detail::read_int(q, 4);
    const auto centre = detail::read_int(q, 8);
    const double s = std::sqrt(0.5);
    EXPECT_NEAR(coord(mesh, n01, 0), s, 1e-15);
    EXPECT_NEAR(coord(mesh, n01, 1), s, 1e-15);
    EXPECT_NEAR(coord(mesh, centre, 0), 1.5 * s, 1e-15);
    EXPECT_NEAR(coord(mesh, centre, 1), 1.5 * s, 1e-15);
    EXPECT_NEAR(detail::read_double(mesh.PointData("u"), static_cast<std::size_t>(centre)), 1.5 * s,
                1e-15);
    EXPECT_EQ(detail::read_int(mesh.CellData("mfem:attribute", 0), 0), 5);
    EXPECT_THROW(meshioplusplus::read_mfem(write_file("MFEM NURBS NC-patch mesh v1.0\n")),
                 ReadError);
}

// One unit square with order-2 nodes whose bottom edge is pulled down, in
// MFEM's Bernstein (H1Pos) and serendipity (H1Ser) spaces: their degrees of
// freedom are coefficients, so the mid-edge and centre nodes are what the
// basis makes of them, not the stored values.
TEST(Mfem, BernsteinAndSerendipityNodesAreCoefficients) {
    const std::string head =
        "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 3 0 1 2 3\nboundary\n0\n"
        "vertices\n4\n\nnodes\nFiniteElementSpace\nFiniteElementCollection: ";
    const std::string corners = "0 0\n1 0\n1 1\n0 1\n";
    // Bernstein: the bottom edge coefficient (0.5, -0.2) weighs 1/2 at the
    // edge midpoint and 1/8 at the centre; the rest reproduce the square.
    const Mesh pos =
        meshioplusplus::read_mfem(write_file(head + "H1Pos_2D_P2\nVDim: 2\nOrdering: 1\n" +
                                             corners + "0.5 -0.2\n1 0.5\n0.5 1\n0 0.5\n0.5 0.5\n"));
    ASSERT_EQ(pos.Cells(0).Type(), "quad9");
    const auto& q = pos.Cells(0).Conn();
    EXPECT_NEAR(coord(pos, detail::read_int(q, 4), 1), -0.1, 1e-15);
    EXPECT_NEAR(coord(pos, detail::read_int(q, 8), 1), 0.475, 1e-15);
    // Serendipity: the edge value is nodal; its shape function is 1/2 at the
    // centre. No interior dof at order 2.
    const Mesh ser =
        meshioplusplus::read_mfem(write_file(head + "H1Ser_2D_P2\nVDim: 2\nOrdering: 1\n" +
                                             corners + "0.5 -0.1\n1 0.5\n0.5 1\n0 0.5\n"));
    const auto& s = ser.Cells(0).Conn();
    EXPECT_NEAR(coord(ser, detail::read_int(s, 4), 1), -0.1, 1e-15);
    EXPECT_NEAR(coord(ser, detail::read_int(s, 8), 1), 0.45, 1e-15);
}

// Two unit squares over two ranks as ParMesh::ParPrint writes a non-conforming
// mesh: each rank's refinement tree, its own leaf and the other's as a ghost,
// every boundary edge listed by both. The leaves merge at their shared edge.
TEST(Mfem, ParallelNonConformingRanksMerge) {
    const std::string dir = mt::temp_path("_pncmesh");
    std::filesystem::create_directories(dir);
    for (int r = 0; r < 2; ++r)
        std::ofstream(dir + "/m.00000" + std::to_string(r))
            << "MFEM NC mesh v1.0\ndimension\n2\nrank\n"
            << r
            << "\nelements\n2\n0 1 3 0 0 1 4 3\n1 2 3 0 1 2 5 4\n"
               "boundary\n6\n1 1 0 1\n1 1 1 2\n2 1 2 5\n1 1 5 4\n1 1 4 3\n3 1 3 0\n"
               "coordinates\n6\n2\n0 0\n1 0\n2 0\n0 1\n1 1\n2 1\n";
    const Mesh mesh = meshioplusplus::read_mfem(dir + "/m.000000");
    EXPECT_EQ(mesh.NumPoints(), 6u);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).NumCells(), 2u);
    EXPECT_EQ(mesh.Cells(1).NumCells(), 6u);  // each rank's own boundary, once
    const auto& parts = mesh.CellData("partition:part", 0);
    EXPECT_EQ(detail::read_int(parts, 0), 0);
    EXPECT_EQ(detail::read_int(parts, 1), 1);
    std::filesystem::remove_all(dir);
}
