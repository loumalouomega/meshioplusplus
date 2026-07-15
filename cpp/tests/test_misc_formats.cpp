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
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/wkt.hpp"

TEST(Su2, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_su2(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_su2(p); };
    mt::roundtrip(w, r, mt::tri_mesh_2d(), ".su2");
    mt::roundtrip(w, r, mt::tet_mesh(), ".su2");
    mt::roundtrip(w, r, mt::hex_mesh(), ".su2");
}

TEST(Flac3d, AsciiAndBinary) {
    for (bool binary : {false, true}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_flac3d(p, m, ".16e", binary);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_flac3d(p); };
        mt::roundtrip(w, r, mt::tet_mesh(), ".f3grid");
        mt::roundtrip(w, r, mt::hex_mesh(), ".f3grid");
    }
}

TEST(Ansys, AsciiAndBinary) {
    for (bool binary : {false, true}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_ansys(p, m, binary);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_ansys(p); };
        mt::roundtrip(w, r, mt::tri_mesh_2d(), ".msh");
        mt::roundtrip(w, r, mt::tet_mesh(), ".msh");
        mt::roundtrip(w, r, mt::hex_mesh(), ".msh");
        mt::roundtrip(w, r, mt::tri_quad_mesh(), ".msh");
    }
}

TEST(Dolfin, TriangleTetra) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_dolfin(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_dolfin(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".xml");
    mt::roundtrip(w, r, mt::tri_mesh_2d(), ".xml");
    mt::roundtrip(w, r, mt::tet_mesh(), ".xml");
}

TEST(Wkt, TriangleGeometry) {
    // WKT (TIN) de-duplicates points, so point order is not preserved; check
    // that the triangle count round-trips.
    mt::Mesh in = mt::tri_mesh();
    std::string path = mt::temp_path(".wkt");
    meshioplusplus::write_wkt(path, in);
    mt::Mesh out = meshioplusplus::read_wkt(path);
    ASSERT_EQ(out.cells.size(), 1u);
    EXPECT_EQ(out.cells[0].type, "triangle");
    EXPECT_EQ(out.cells[0].num_cells(), 2u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tetgen, TetraPair) {
    // TetGen writes a .node/.ele pair sharing a stem.
    mt::Mesh in = mt::tet_mesh();
    std::string node = mt::temp_path(".node");
    meshioplusplus::write_tetgen(node, in);
    mt::Mesh out = meshioplusplus::read_tetgen(node);
    mt::expect_mesh_eq(in, out);
    std::string ele = node.substr(0, node.size() - 5) + ".ele";
    std::error_code ec;
    std::filesystem::remove(node, ec);
    std::filesystem::remove(ele, ec);
}
