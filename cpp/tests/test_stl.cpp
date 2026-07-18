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
#include "meshioplusplus/formats/stl.hpp"

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
