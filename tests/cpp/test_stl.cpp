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

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/stl.hpp"

// System includes
#include <fstream>

namespace {

// STL stores raw triangle coordinates and re-derives point indices, so the
// point *order* is not preserved. Compare the set of triangle coordinate
// triples instead.
void stl_roundtrip(const mt::Mesh& mesh, bool binary) {
    std::string path = mt::temp_path(binary ? "_bin.stl" : "_asc.stl");
    meshioplusplus::write_stl(path, mesh, binary);
    mt::Mesh out = meshioplusplus::read_stl(path);

    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    EXPECT_EQ(out.Cells(0).NumCells(), mesh.Cells(0).NumCells());

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// Writes `rText` verbatim (no newline translation) and reads it back.
mt::Mesh read_stl_text(const std::string& rText) {
    const std::string path = mt::temp_path("_text.stl");
    {
        std::ofstream os(path, std::ios::binary);
        os << rText;
    }
    struct Cleanup {
        const std::string& rPath;
        ~Cleanup() {
            std::error_code ec;
            std::filesystem::remove(rPath, ec);
        }
    } cleanup{path};
    return meshioplusplus::read_stl(path);
}

// Two facets sharing an edge: four unique points, the second normal tilted.
const char kTwoFacets[] =
    "solid pin\n"
    "  facet normal 0 0 1\n"
    "    outer loop\n"
    "      vertex 0 0 0\n"
    "      vertex 1 0 0\n"
    "      vertex 0 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "  facet normal 1.5e-1 -2.5E+0 +3\n"
    "    outer loop\n"
    "      vertex 1 0 0\n"
    "      vertex 1 1 0\n"
    "      vertex 0 1 0\n"
    "    endloop\n"
    "  endfacet\n"
    "endsolid pin\n";

std::string with_crlf(const std::string& rText) {
    std::string out;
    for (const char c : rText) {
        if (c == '\n')
            out += '\r';
        out += c;
    }
    return out;
}

// The facets of kTwoFacets, as read: the point table, the connectivity and the
// per-facet normals (an ASCII reader returns Float64 points).
void expect_two_facets(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 4u);
    const double* p = rMesh.Points().As<double>();
    const double expected[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t d = 0; d < 3; ++d)
            EXPECT_EQ(p[3 * i + d], expected[i][d]) << "point " << i << " axis " << d;
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    ASSERT_EQ(rMesh.Cells(0).NumCells(), 2u);
    const std::int64_t* c = rMesh.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t expected_cells[6] = {0, 1, 2, 1, 3, 2};
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(c[i], expected_cells[i]) << "connectivity " << i;
    ASSERT_TRUE(rMesh.HasCellData("facet_normals"));
    const double* n = rMesh.CellData("facet_normals", 0).As<double>();
    const double expected_normals[6] = {0, 0, 1, 0.15, -2.5, 3};
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_DOUBLE_EQ(n[i], expected_normals[i]) << "normal " << i;
}

}  // namespace

TEST(Stl, TriMeshAscii) {
    stl_roundtrip(mt::tri_mesh(), false);
}
TEST(Stl, TriMeshBinary) {
    stl_roundtrip(mt::tri_mesh(), true);
}

TEST(Stl, GeometryPreserved) {
    // Verify the triangle vertex coordinates survive a round-trip.
    mt::Mesh in = mt::tri_mesh();
    std::string path = mt::temp_path(".stl");
    meshioplusplus::write_stl(path, in, false);
    mt::Mesh out = meshioplusplus::read_stl(path);
    // both meshes describe the same 2 triangles over the unit square
    EXPECT_EQ(out.Cells(0).NumCells(), 2u);
    EXPECT_GE(out.NumPoints(), 3u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Stl, SkinTetMesh) {
    // Default skin=true: a tetra mesh writes its boundary skin (6 triangles).
    std::string path = mt::temp_path("_skin.stl");
    meshioplusplus::write_stl(path, mt::tet_mesh(), /*binary=*/false);
    mt::Mesh out = meshioplusplus::read_stl(path);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    EXPECT_EQ(out.Cells(0).NumCells(), 6u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Stl, SkinHexMeshTriangulatesQuads) {
    // A hexahedron's skin is 6 quads -> 12 triangles after triangulation.
    std::string path = mt::temp_path("_skin_hex.stl");
    meshioplusplus::write_stl(path, mt::hex_mesh(), /*binary=*/true);
    mt::Mesh out = meshioplusplus::read_stl(path);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).NumCells(), 12u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Stl, SkinFalseLegacyDropsVolumeCells) {
    // skin=false keeps the legacy behavior: no triangle cells -> empty STL.
    std::string path = mt::temp_path("_legacy.stl");
    meshioplusplus::write_stl(path, mt::tet_mesh(), /*binary=*/false, /*skin=*/false);
    mt::Mesh out = meshioplusplus::read_stl(path);
    EXPECT_EQ(out.NumCellBlocks(), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The ASCII reader takes the last three tokens of every line that does not
// open with a keyword (`facet normal ...` and `vertex ...`) as one row.
TEST(Stl, AsciiTextExponentsAndSigns) {
    expect_two_facets(read_stl_text(kTwoFacets));
}

TEST(Stl, AsciiCrlfLineEndings) {
    expect_two_facets(read_stl_text(with_crlf(kTwoFacets)));
}

TEST(Stl, AsciiTabsBlankLinesAndIndentation) {
    std::string text = kTwoFacets;
    // Tabs and a form feed for the blanks, blank lines between records, and a
    // vertical tab before a keyword: all leading whitespace is skipped.
    for (std::size_t pos = 0; (pos = text.find("vertex ", pos)) != std::string::npos; pos += 8)
        text.replace(pos + 6, 1, "\t");
    text.insert(text.find("endloop"), "\n\v   \n");
    text.insert(text.find("  facet normal 1.5"), "\f");
    expect_two_facets(read_stl_text(text));
}

TEST(Stl, AsciiShortLinesAreIgnored) {
    // A line with fewer than three tokens carries no row, whatever it is.
    std::string text = kTwoFacets;
    text.insert(text.find("  facet normal 1.5"), "stray token\n42\n");
    expect_two_facets(read_stl_text(text));
}

TEST(Stl, AsciiTruncatedFacetIsMalformed) {
    std::string text = kTwoFacets;
    const std::size_t last = text.rfind("      vertex 0 1 0\n");
    text.erase(last, std::string("      vertex 0 1 0\n").size());
    EXPECT_THROW(read_stl_text(text), meshioplusplus::ReadError);
}

TEST(Stl, AsciiSmallFileKeepsItsFirstLine) {
    // Under 80 bytes the file goes straight to the ASCII reader, with no
    // header line to skip: the first facet line is a row.
    const std::string text =
        "facet normal 0 0 1\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\n";
    ASSERT_LT(text.size(), 80u);
    const mt::Mesh m = read_stl_text(text);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).NumCells(), 1u);
    EXPECT_EQ(m.NumPoints(), 3u);
}

TEST(Stl, AsciiEmptySolidHasNoCells) {
    const mt::Mesh m = read_stl_text("solid empty\nendsolid empty\n");
    EXPECT_EQ(m.NumCellBlocks(), 0u);
    EXPECT_EQ(m.NumPoints(), 0u);
}

TEST(Stl, FileOf80To83BytesIsAsciiWithoutATriangleCount) {
    // Too short for a binary header's triangle count, so ASCII, whose first
    // line is the header: nothing is left to read.
    for (const std::size_t size : {80u, 81u, 83u}) {
        const std::string text = std::string(size - 1, 'x') + "\n";
        const mt::Mesh m = read_stl_text(text);
        EXPECT_EQ(m.NumCellBlocks(), 0u) << size;
        EXPECT_EQ(m.NumPoints(), 0u) << size;
    }
}

TEST(Stl, BinaryKeepsEveryVertexInSinglePrecision) {
    const mt::Mesh in = mt::tri_mesh();
    const std::string path = mt::temp_path("_binvalues.stl");
    meshioplusplus::write_stl(path, in, /*binary=*/true);
    const mt::Mesh out = meshioplusplus::read_stl(path);
    std::error_code ec;
    std::filesystem::remove(path, ec);

    ASSERT_EQ(out.Points().Dtype(), meshioplusplus::DType::Float32);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    ASSERT_EQ(out.Cells(0).NumCells(), in.Cells(0).NumCells());
    // Every triangle corner of the result is a corner of the input (the unit
    // square's coordinates are exact in float32).
    const float* p = out.Points().As<float>();
    const std::int64_t* c = out.Cells(0).Conn().As<std::int64_t>();
    for (std::size_t i = 0; i < 3 * out.Cells(0).NumCells(); ++i) {
        const float* v = p + 3 * static_cast<std::size_t>(c[i]);
        bool found = false;
        for (std::size_t k = 0; k < in.NumPoints() && !found; ++k) {
            const double* q = in.Points().As<double>() + 3 * k;
            found = v[0] == static_cast<float>(q[0]) && v[1] == static_cast<float>(q[1]) &&
                    v[2] == static_cast<float>(q[2]);
        }
        EXPECT_TRUE(found) << "corner " << i;
    }
}
