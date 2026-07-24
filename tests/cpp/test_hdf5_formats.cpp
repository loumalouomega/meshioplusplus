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

#include <fstream>

#include <hdf5.h>

#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/med.hpp"

using meshioplusplus::detail::read_double;  // NOLINT
namespace h5 = meshioplusplus::h5;

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

TEST(Med, FieldsNodalAndCellRoundTrip) {
    // The single-timestep common case: ordinary point_data/cell_data arrays,
    // no units, no component names, no multiple timesteps. Two DISTINCT cell
    // types (MED rejects two blocks of the same type), so the "field values
    // follow their cell type" property below is actually exercised.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 4}}));
    m.AddCellBlock("quad", mt::conn_from({{1, 2, 5, 4}}));

    meshioplusplus::NDArray temperature(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        temperature.As<double>()[i] = static_cast<double>(i) + 0.5;
    m.AddPointData("temperature", std::move(temperature));

    std::vector<meshioplusplus::NDArray> stress;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        meshioplusplus::NDArray blk(meshioplusplus::DType::Float64, {m.Cells(b).NumCells()});
        for (std::size_t c = 0; c < m.Cells(b).NumCells(); ++c)
            blk.As<double>()[c] = 100.0 * static_cast<double>(b + 1) + static_cast<double>(c);
        stress.push_back(std::move(blk));
    }
    m.AddCellData("stress", std::move(stress));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasPointData("temperature"));
    const meshioplusplus::NDArray& t = out.PointData("temperature");
    ASSERT_EQ(t.Size(), m.NumPoints());
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        EXPECT_DOUBLE_EQ(read_double(t, i), static_cast<double>(i) + 0.5);

    ASSERT_TRUE(out.HasCellData("stress"));
    // Field values must follow their CELL TYPE, not a fixed block index --
    // MED reorders blocks alphabetically by type code on read, so a test that
    // assumed positional stability would pass by accident, not by contract.
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const auto cb = out.Cells(b);
        const meshioplusplus::NDArray& d = out.CellData("stress", b);
        ASSERT_EQ(d.Size(), cb.NumCells());
        const double base = cb.Type() == "triangle" ? 100.0 : 200.0;
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            EXPECT_DOUBLE_EQ(read_double(d, c), base + static_cast<double>(c));
    }
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, FieldsDeferToPythonWhenUnitsPresent) {
    // A hand-set UNI/UNT attribute is real information this reader does not
    // carry -- it must decline (throw) rather than silently drop it.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray temperature(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        temperature.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("temperature", std::move(temperature));
    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid field = h5::open_group(f, "CHA/temperature");
        H5Adelete(field, "UNI");  // the writer already created it empty
        h5::write_attr_string(field, "UNI", "K");
    }

    meshioplusplus::MedInfo info;
    EXPECT_THROW(meshioplusplus::read_med(p, info), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, FieldsDeferToPythonWhenTimestepMetadataIsNonDefault) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray temperature(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        temperature.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("temperature", std::move(temperature));
    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid field = h5::open_group(f, "CHA/temperature");
        std::vector<std::string> steps = h5::group_links(field);
        ASSERT_EQ(steps.size(), 1u);
        h5::Hid ts = h5::open_group(field, steps[0]);
        h5::Hid a(H5Aopen(ts, "NDT", H5P_DEFAULT), H5Aclose);
        std::int64_t seven = 7;
        H5Awrite(a, H5T_NATIVE_INT64, &seven);
    }

    meshioplusplus::MedInfo info;
    EXPECT_THROW(meshioplusplus::read_med(p, info), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, ReadRejectsNonHdf5) {
    // A file that is not HDF5 must surface as a ReadError (the shared hdf5_util
    // open helper), not an uncaught HDF5 abort.
    std::string p = mt::temp_path(".med");
    {
        std::ofstream f(p);
        f << "not an HDF5 file\n";
    }
    meshioplusplus::MedInfo info;
    EXPECT_THROW(meshioplusplus::read_med(p, info), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
