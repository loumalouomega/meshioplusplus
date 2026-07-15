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

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/med.hpp"

TEST(Cgns, TetraCompressed) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_cgns(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_cgns(p); };
        mt::roundtrip(w, r, mt::tet_mesh(), ".cgns");
    }
}

TEST(H5m, LineTriangleTetra) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_h5m(p, m, true, gzip);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_h5m(p); };
        mt::roundtrip(w, r, mt::line_mesh(), ".h5m");
        mt::roundtrip(w, r, mt::tri_mesh(), ".h5m");
        mt::roundtrip(w, r, mt::tet_mesh(), ".h5m");
    }
}

TEST(Hmf, Basic) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_hmf(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_hmf(p); };
        mt::roundtrip(w, r, mt::tri_mesh(), ".hmf");
        mt::roundtrip(w, r, mt::tet_mesh(), ".hmf");
        mt::roundtrip(w, r, mt::hex_mesh(), ".hmf");
    }
}

TEST(Med, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    };
    auto r = [](const std::string& p) {
        meshioplusplus::MedInfo info;
        return meshioplusplus::read_med(p, info);
    };
    mt::roundtrip(w, r, mt::tri_mesh(), ".med");
    mt::roundtrip(w, r, mt::tet_mesh(), ".med");   // exercises node perm
    mt::roundtrip(w, r, mt::hex_mesh(), ".med");   // exercises node perm
}

TEST(Med, MetadataAndFamilies) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::MedInfo win;
    win.mesh_name = "mymesh";
    win.description = "hello";
    win.unit_coords = "mm";
    win.cell_tags[-1] = {"top"};
    win.cell_tag_groups[-1] = "FAM_-1_top";

    meshioplusplus::Mesh m = mt::tri_mesh();
    // one cell_tags block matching the single triangle block
    meshioplusplus::NDArray tag(meshioplusplus::DType::Int64, {m.cells[0].num_cells()});
    for (std::size_t i = 0; i < m.cells[0].num_cells(); ++i)
        tag.as<std::int64_t>()[i] = -1;
    m.cell_data["cell_tags"] = {tag};

    meshioplusplus::write_med(p, m, win);
    meshioplusplus::MedInfo rout;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, rout);
    EXPECT_EQ(rout.mesh_name, "mymesh");
    EXPECT_EQ(rout.description, "hello");
    EXPECT_EQ(rout.unit_coords, "mm");
    ASSERT_TRUE(rout.cell_tags.count(-1));
    EXPECT_EQ(rout.cell_tags[-1], (std::vector<std::string>{"top"}));
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, RaggedPolygons) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.points = mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}});
    meshioplusplus::CellBlock cb;
    cb.type = "polygon";
    cb.polygon_rows = {{0, 1, 2}, {1, 3, 4, 2, 5}};  // a tri and a 5-gon
    m.cells.push_back(std::move(cb));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    ASSERT_EQ(out.cells.size(), 1u);
    EXPECT_EQ(out.cells[0].type, "polygon");
    ASSERT_EQ(out.cells[0].polygon_rows.size(), 2u);
    EXPECT_EQ(out.cells[0].polygon_rows[0], (std::vector<std::int64_t>{0, 1, 2}));
    EXPECT_EQ(out.cells[0].polygon_rows[1],
              (std::vector<std::int64_t>{1, 3, 4, 2, 5}));
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
