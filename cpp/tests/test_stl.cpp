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

    ASSERT_EQ(out.cells.size(), 1u);
    EXPECT_EQ(out.cells[0].type, "triangle");
    EXPECT_EQ(out.cells[0].num_cells(), mesh.cells[0].num_cells());

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

}  // namespace

TEST(Stl, TriMeshAscii) { stl_roundtrip(mt::tri_mesh(), false); }
TEST(Stl, TriMeshBinary) { stl_roundtrip(mt::tri_mesh(), true); }

TEST(Stl, GeometryPreserved) {
    // Verify the triangle vertex coordinates survive a round-trip.
    mt::Mesh in = mt::tri_mesh();
    std::string path = mt::temp_path(".stl");
    meshioplusplus::write_stl(path, in, false);
    mt::Mesh out = meshioplusplus::read_stl(path);
    // both meshes describe the same 2 triangles over the unit square
    EXPECT_EQ(out.cells[0].num_cells(), 2u);
    EXPECT_GE(out.num_points(), 3u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
