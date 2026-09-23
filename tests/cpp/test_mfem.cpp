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
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mfem.hpp"
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
    expect_error("MFEM NC mesh v1.0\n", "non-conforming");
    expect_error(
        "MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 3 0 1 2 9\nboundary\n0\n"
        "vertices\n3\n2\n0 0\n1 0\n0 1\n",
        "out of range");
    expect_error("MFEM mesh v1.0\ndimension\n2\nelements\n1\n1 2 0 1\n", "ends");
    expect_error("not a mesh\n", "not an MFEM mesh");
}
