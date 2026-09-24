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
 * @file test_feconv_quirks.cpp
 * @brief The native readers on the quirks FEconv's samples exposed: FLUX's
 *        mirrored solids and truncated files, indented Gmsh files, Medit files
 *        without `Dimension`, raw and BigEndian VTU appended data and
 *        multi-piece VTU files. The Python tests read the generated fixtures
 *        (tools/gen_feconv_quirk_fixtures.py) through both engines.
 */

#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/vtu.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
namespace detail = meshioplusplus::detail;

std::string fq_write_file(const std::string& rBody, const char* pExt) {
    const std::string path = mt::temp_path(pExt);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::vector<std::int64_t> fq_row(const Mesh& rMesh, std::size_t Block, std::size_t Row) {
    const auto cb = rMesh.Cells(Block);
    const std::size_t k = detail::cols(cb.Conn());
    std::vector<std::int64_t> out(k);
    for (std::size_t j = 0; j < k; ++j)
        out[j] = detail::read_int(cb.Conn(), Row * k + j);
    return out;
}

std::string fq_flux_header(int Nel, int Nnod) {
    std::string s = " test\n       3           NOMBRE DE DIMENSIONS DU DECOUPAGE\n";
    s += "      " + std::to_string(Nel) + "           NOMBRE  D'ELEMENTS\n";
    s += "      " + std::to_string(Nel) + "           NOMBRE  D'ELEMENTS VOLUMIQUES\n";
    s += "      " + std::to_string(Nnod) + "           NOMBRE DE POINTS\n";
    return s;
}

}  // namespace

TEST(FeconvQuirks, FluxTetraIsMirrored) {
    // Base (0,1,2) clockwise seen from the apex: FLUX's layout.
    const std::string body = fq_flux_header(1, 4) +
                             " DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n"
                             "  1  5  4  7  4  0 10  4  0  0  0  0\n"
                             "  1  2  3  4\n"
                             " COORDONNEES DES NOEUDS\n"
                             "  1 0 0 0\n  2 0 1 0\n  3 1 0 0\n  4 0 0 1\n"
                             " ==== DECOUPAGE  TERMINE\n";
    const Mesh m = meshioplusplus::read_flux(fq_write_file(body, ".pf3"));
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{0, 2, 1, 3}));

    // Written back in FLUX's layout.
    const std::string out = mt::temp_path(".pf3");
    meshioplusplus::write_flux(out, m);
    EXPECT_EQ(fq_row(meshioplusplus::read_flux(out), 0, 0),
              (std::vector<std::int64_t>{0, 2, 1, 3}));
}

TEST(FeconvQuirks, FluxTruncatedSparseIds) {
    const std::string body = fq_flux_header(40, 900) +
                             " DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n"
                             "  3050  4  303  4  8  0  8  8  0  0  0  0\n"
                             "  3 121 749 112 5521 5525 5329 5522\n"
                             " COORDONNEES DES NOEUDS\n"
                             "  3 0 1 0\n  121 1 1 0\n  749 1 0 0\n  112 0 0 0\n"
                             "  5521 0.5 1 0\n  5525 1 0.5 0\n  5329 0.5 0 0\n  5522 0 0.5 0\n"
                             " ==== DECOUPAGE  TERMINE\n";
    const Mesh m = meshioplusplus::read_flux(fq_write_file(body, ".pf3"));
    EXPECT_EQ(m.NumPoints(), 8u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "quad8");
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{0, 1, 2, 3, 4, 5, 6, 7}));

    std::string bad = body;
    bad.replace(bad.find("5329"), 4, "9999");
    EXPECT_THROW(meshioplusplus::read_flux(fq_write_file(bad, ".pf3")), ReadError);
}

TEST(FeconvQuirks, GmshIndentedSections) {
    const char* lines[] = {
        "$MeshFormat", "2.2 0 8", "$EndMeshFormat",    "$Nodes",      "4",
        "1 0 0 0",     "2 1 0 0", "3 1 1 0",           "4 0 1 0",     "$EndNodes",
        "$Elements",   "1",       "1 3 2 1 1 1 2 3 4", "$EndElements"};
    std::string body;
    for (const char* pLine : lines)
        body += std::string("     ") + pLine + "   \n";
    const Mesh m = meshioplusplus::read_gmsh(fq_write_file(body, ".msh"));
    EXPECT_EQ(m.NumPoints(), 4u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "quad");
}

TEST(FeconvQuirks, MeditWithoutDimension) {
    const std::string body =
        "MeshVersionFormatted 1\n\nVertices\n4\n"
        "0 0 0 0\n1 0 0 0\n0 1 0 0\n0 0 1 0\n\nTetrahedra\n1\n1 2 3 4 5\n\nTriangles\n0\n\nEnd\n";
    const Mesh m = meshioplusplus::read_medit_ascii(fq_write_file(body, ".mesh"));
    EXPECT_EQ(m.NumPoints(), 4u);
    EXPECT_EQ(m.PointDim(), 3u);
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{0, 1, 2, 3}));

    // A 2-D file: three values per vertex row.
    const std::string flat =
        "MeshVersionFormatted 2\nVertices\n3\n0 0 1\n1 0 1\n0 1 1\nTriangles\n1\n1 2 3 0\nEnd\n";
    EXPECT_EQ(meshioplusplus::read_medit_ascii(fq_write_file(flat, ".mesh")).PointDim(), 2u);
}

namespace {

void fq_put_be(std::string& rOut, std::uint64_t V, int Size) {
    for (int i = Size - 1; i >= 0; --i)
        rOut.push_back(static_cast<char>((V >> (8 * i)) & 0xff));
}

}  // namespace

TEST(FeconvQuirks, VtuRawAppendedBigEndianPaddedOffsets) {
    // One triangle; each array a 4-byte BigEndian byte count, then the body.
    std::string payload;
    std::vector<std::size_t> offsets;
    auto add_f64 = [&](const std::vector<double>& rV) {
        offsets.push_back(payload.size());
        fq_put_be(payload, rV.size() * 8, 4);
        for (double d : rV) {
            std::uint64_t bits;
            std::memcpy(&bits, &d, 8);
            fq_put_be(payload, bits, 8);
        }
    };
    auto add_int = [&](const std::vector<std::int64_t>& rV, int Size) {
        offsets.push_back(payload.size());
        fq_put_be(payload, rV.size() * static_cast<std::size_t>(Size), 4);
        for (std::int64_t v : rV)
            fq_put_be(payload, static_cast<std::uint64_t>(v), Size);
    };
    add_f64({0, 0, 0, 1, 0, 0, 0, 1, 0});
    add_int({0, 1, 2}, 4);
    add_int({3}, 4);
    add_int({5}, 1);
    // '<' and '>' bytes inside the payload: the XML parser must never see them.
    add_int({0x3C3E3C3E}, 4);
    auto da = [&](const char* pType, const char* pName, std::size_t Off, const char* pExtra) {
        return std::string("<DataArray type=\"") + pType + "\" Name=\"" + pName + "\"" + pExtra +
               " format=\"appended\" offset=\"    " + std::to_string(Off) + "\"/>\n";
    };
    const std::string xml =
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        "byte_order=\"BigEndian\">\n<UnstructuredGrid>\n<Piece NumberOfPoints=\"   3\" "
        "NumberOfCells=\"   1\">\n<Points>" +
        da("Float64", "Points", offsets[0], " NumberOfComponents=\"3\"") + "</Points>\n<Cells>" +
        da("Int32", "connectivity", offsets[1], "") + da("Int32", "offsets", offsets[2], "") +
        da("UInt8", "types", offsets[3], "") + "</Cells>\n<CellData>" +
        da("Int32", "tag", offsets[4], "") +
        "</CellData>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_" + payload +
        "\n</AppendedData>\n</VTKFile>\n";
    const Mesh m = meshioplusplus::read_vtu(fq_write_file(xml, ".vtu"));
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_DOUBLE_EQ(detail::read_double(m.Points(), 3), 1.0);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "triangle");
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{0, 1, 2}));
    ASSERT_TRUE(m.HasCellData("tag"));
    EXPECT_EQ(detail::read_int(m.CellData("tag", 0), 0), 0x3C3E3C3E);
}

TEST(FeconvQuirks, VtuTwoPiecesMerge) {
    auto piece = [](const char* pPointData) {
        return std::string(
                   "<Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">\n"
                   "<Points><DataArray type=\"Float32\" NumberOfComponents=\"3\" "
                   "format=\"ascii\">0 0 0 1 0 0 0 1 0</DataArray></Points>\n<Cells>\n"
                   "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">0 1 "
                   "2</DataArray>\n"
                   "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">3</DataArray>\n"
                   "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">5</DataArray>\n"
                   "</Cells>\n<CellData><DataArray type=\"Int32\" Name=\"tag\" "
                   "format=\"ascii\">4</DataArray></CellData>\n<PointData>") +
               pPointData + "</PointData>\n</Piece>\n";
    };
    const std::string xml =
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
        "<UnstructuredGrid>\n" +
        piece("<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">1 2 3</DataArray>") +
        piece("") + "</UnstructuredGrid>\n</VTKFile>\n";
    const Mesh m = meshioplusplus::read_vtu(fq_write_file(xml, ".vtu"));
    EXPECT_EQ(m.NumPoints(), 6u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
    EXPECT_EQ(fq_row(m, 0, 1), (std::vector<std::int64_t>{3, 4, 5}));
    EXPECT_TRUE(m.HasCellData("tag"));
    EXPECT_FALSE(m.HasPointData("u"));  // not in every piece
}

TEST(FeconvQuirks, FluentCellsFromFaces3d) {
    // One tetrahedron known only through its faces: each face's right-hand
    // normal points into c0 (cell 1), c1 = 0 outside.
    const std::string body =
        "(2 3)\n(10 (0 1 4 0 3))\n(10 (1 1 4 1 3)(\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n))\n"
        "(13 (3 1 4 3 3)(\n1 2 3 1 0\n1 4 2 1 0\n2 4 3 1 0\n3 4 1 1 0\n))\n"
        "(12 (0 1 1 0))\n(12 (4 1 1 1 2))\n(45 (4 fluid block)())\n(45 (3 wall skin)())\n";
    const Mesh m = meshioplusplus::read_ansys(fq_write_file(body, ".msh"));
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(m.Cells(1).Type(), "triangle");
    EXPECT_EQ(m.Cells(1).NumCells(), 4u);
    // Positive volume.
    const auto t = fq_row(m, 0, 0);
    auto p = [&](std::int64_t i, int c) {
        return detail::read_double(m.Points(), static_cast<std::size_t>(i) * 3 + c);
    };
    double a[3], b[3], c[3];
    for (int k = 0; k < 3; ++k) {
        a[k] = p(t[1], k) - p(t[0], k);
        b[k] = p(t[2], k) - p(t[0], k);
        c[k] = p(t[3], k) - p(t[0], k);
    }
    const double vol = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                       a[2] * (b[0] * c[1] - b[1] * c[0]);
    EXPECT_GT(vol, 0.0);
    ASSERT_TRUE(m.HasCellData("ansys:zone"));
    EXPECT_EQ(detail::read_int(m.CellData("ansys:zone", 0), 0), 4);
    EXPECT_EQ(detail::read_int(m.CellData("ansys:zone", 1), 0), 3);
    ASSERT_EQ(m.NumRegions(), 2u);
}

TEST(FeconvQuirks, FluentTwoDimensionalFromEdges) {
    // A unit square as one quad; node zones out of order, mixed face rows.
    const std::string body =
        "(2 2)\n(10 (6 3 4 1 2) (\n1 1\n0 1\n))\n(10 (7 1 2 1 2) (\n0 0\n1 0\n))\n"
        "(13(9 1 4 3 0)(\n2 1 2 1 0\n2 2 3 1 0\n2 3 4 1 0\n2 4 1 1 0\n))\n"
        "(12 (a 1 1 1 3))\n";
    const Mesh m = meshioplusplus::read_ansys(fq_write_file(body, ".msh"));
    EXPECT_EQ(m.PointDim(), 2u);
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "quad");
    // Counter-clockwise from the first edge's start (reversed: the cell is c0).
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{1, 2, 3, 0}));
    EXPECT_EQ(m.Cells(1).Type(), "line");
}

TEST(FeconvQuirks, FluentLegacyConnectivityStillReads) {
    // meshio's old layout (a cell section with a connectivity body), which
    // meshio++ no longer writes: read as cells only.
    const std::string body =
        "(2 3)\n(10 (0 1 4 0))\n(10 (1 1 4 1 3)(\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n))\n"
        "(12 (0 1 1 0))\n(12 (1 1 1 1 2)(\n1 2 3 4\n))\n";
    const Mesh m = meshioplusplus::read_ansys(fq_write_file(body, ".msh"));
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).Type(), "tetra");
    EXPECT_EQ(fq_row(m, 0, 0), (std::vector<std::int64_t>{0, 1, 2, 3}));
    EXPECT_FALSE(m.HasCellData("ansys:zone"));
}

TEST(FeconvQuirks, FluentFaceWriterRoundTrips) {
    // The writer emits faces with c0/c1; the reader rebuilds the cells from
    // them and keeps the boundary faces as a wall zone.
    const Mesh src = mt::tet_mesh();
    for (bool binary : {false, true}) {
        const std::string path = mt::temp_path(".msh");
        meshioplusplus::write_ansys(path, src, binary);
        const Mesh m = meshioplusplus::read_ansys(path);
        ASSERT_GE(m.NumCellBlocks(), 2u);
        EXPECT_EQ(m.Cells(0).Type(), "tetra");
        ASSERT_EQ(m.Cells(0).NumCells(), src.Cells(0).NumCells());
        for (std::size_t c = 0; c < src.Cells(0).NumCells(); ++c) {
            auto a = fq_row(m, 0, c), b = fq_row(src, 0, c);
            std::sort(a.begin(), a.end());
            std::sort(b.begin(), b.end());
            EXPECT_EQ(a, b);
        }
        EXPECT_EQ(m.Cells(1).Type(), "triangle");
        EXPECT_TRUE(m.HasCellData("ansys:zone"));
        EXPECT_EQ(m.NumRegions(), 2u);  // the fluid zone and the default wall
    }
}
