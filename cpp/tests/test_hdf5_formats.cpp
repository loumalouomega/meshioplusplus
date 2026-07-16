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
    mt::roundtrip(w, r, mt::tet_mesh(), ".med");  // exercises node perm
    mt::roundtrip(w, r, mt::hex_mesh(), ".med");  // exercises node perm
}

TEST(Med, MetadataAndFamilies) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::MedInfo win;
    win.mMeshName = "mymesh";
    win.mDescription = "hello";
    win.mUnitCoords = "mm";
    win.mCellTags[-1] = {"top"};
    win.mCellTagGroups[-1] = "FAM_-1_top";

    meshioplusplus::Mesh m = mt::tri_mesh();
    // one cell_tags block matching the single triangle block
    const std::size_t ntri = m.Cells(0).NumCells();
    meshioplusplus::NDArray tag(meshioplusplus::DType::Int64, {ntri});
    for (std::size_t i = 0; i < ntri; ++i)
        tag.As<std::int64_t>()[i] = -1;
    m.AddCellData("cell_tags", {std::move(tag)});

    meshioplusplus::write_med(p, m, win);
    meshioplusplus::MedInfo rout;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, rout);
    EXPECT_EQ(rout.mMeshName, "mymesh");
    EXPECT_EQ(rout.mDescription, "hello");
    EXPECT_EQ(rout.mUnitCoords, "mm");
    ASSERT_TRUE(rout.mCellTags.count(-1));
    EXPECT_EQ(rout.mCellTags[-1], (std::vector<std::string>{"top"}));
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, RaggedPolygons) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2}, {1, 3, 4, 2, 5}});  // a tri and a 5-gon

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    const auto cb = out.Cells(0);
    EXPECT_EQ(cb.Type(), "polygon");
    ASSERT_TRUE(cb.IsRagged());
    ASSERT_EQ(cb.NumCells(), 2u);
    EXPECT_EQ(std::vector<std::int64_t>(cb.Row(0), cb.Row(0) + cb.RowSize(0)),
              (std::vector<std::int64_t>{0, 1, 2}));
    EXPECT_EQ(std::vector<std::int64_t>(cb.Row(1), cb.Row(1) + cb.RowSize(1)),
              (std::vector<std::int64_t>{1, 3, 4, 2, 5}));
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
