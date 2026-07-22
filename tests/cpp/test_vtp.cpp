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
#include <filesystem>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"

namespace {

void rt(const mt::Mesh& mesh, bool binary, bool zlib) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_vtp(p, m, binary, zlib); },
                  [](const std::string& p) { return meshioplusplus::read_vtp(p); }, mesh, ".vtp");
}

mt::Mesh vertex_line_tri_mesh() {
    mt::Mesh mesh = mt::tri_mesh();
    mesh.AddCellBlock("vertex", mt::conn_from({{0}, {2}}));
    mesh.AddCellBlock("line", mt::conn_from({{0, 1}, {1, 2}}));
    return mesh;
}

// Rectangular (uniform row size) polygon block — mt::roundtrip's comparison
// requires Conn(), so the jagged case gets its own manual test below.
mt::Mesh polygon_mesh() {
    mt::Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1.5, 1, 0}, {0.5, 1.7, 0}, {-0.5, 1, 0}},
                                  "triangle", {{0, 1, 2}});
    mesh.AddCellBlock("polygon", mt::conn_from({{0, 1, 2, 3, 4}}));
    return mesh;
}

}  // namespace

TEST(Vtp, AsciiSurface) {
    rt(mt::tri_mesh(), false, false);
    rt(mt::quad_mesh(), false, false);
    rt(mt::tri_quad_mesh(), false, false);
}
TEST(Vtp, Ascii2D) {
    rt(mt::tri_mesh_2d(), false, false);
}
TEST(Vtp, AsciiVertsLines) {
    rt(vertex_line_tri_mesh(), false, false);
}
TEST(Vtp, AsciiPolygon) {
    rt(polygon_mesh(), false, false);
}
TEST(Vtp, BinaryUncompressed) {
    rt(mt::tri_mesh(), true, false);
    rt(polygon_mesh(), true, false);
}
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
TEST(Vtp, BinaryZlib) {
    rt(mt::tri_mesh(), true, true);
    rt(mt::tri_quad_mesh(), true, true);
}
#endif

TEST(Vtp, RaggedPolygonWrite) {
    // A jagged polygon block writes as Polys rows; on read the 4-noded row
    // legitimately comes back as a quad (PolyData has no cell-type array).
    mt::Mesh mesh = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1.5, 1, 0}, {0.5, 1.7, 0}, {-0.5, 1, 0}},
                                  "triangle", {{0, 1, 2}});
    mesh.AddPolygonBlock("polygon", {{0, 1, 2, 3, 4}, {0, 1, 2, 3}});
    const std::string path = mt::temp_path(".vtp");
    meshioplusplus::write_vtp(path, mesh, /*binary=*/false, /*zlib=*/false);
    mt::Mesh out = meshioplusplus::read_vtp(path);
    ASSERT_EQ(out.NumCellBlocks(), 3u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    EXPECT_EQ(out.Cells(1).Type(), "polygon");
    EXPECT_EQ(out.Cells(1).NumCells(), 1u);
    EXPECT_EQ(out.Cells(1).NodesPerCell(), 5u);
    EXPECT_EQ(out.Cells(2).Type(), "quad");
    EXPECT_EQ(out.Cells(2).NumCells(), 1u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Vtp, VolumeCellsRejected) {
    EXPECT_THROW(meshioplusplus::write_vtp(mt::temp_path(".vtp"), mt::tet_mesh(), false, false),
                 meshioplusplus::WriteError);
    EXPECT_THROW(
        meshioplusplus::write_vtp(mt::temp_path(".vtp"), mt::triangle6_mesh(), false, false),
        meshioplusplus::WriteError);
}
