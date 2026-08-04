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
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::detail::read_double;  // NOLINT
namespace h5 = meshioplusplus::h5;

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
    // types, so the "field values follow their cell type" property below is
    // actually exercised (same-type blocks are consolidated into one section
    // since v9.8.0 -- see Med.SameTypeBlocksAreConsolidated in this file).
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

namespace {

/// A `tri_mesh` with one nodal field, written by the C++ writer, ready for a
/// test to patch into an "enhanced" file the strict reader declines.
std::string med_field_fixture() {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray t(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        t.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("temperature", std::move(t));
    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    return p;
}

meshioplusplus::ReadOptions med_lenient() {
    meshioplusplus::ReadOptions o;
    o.mLenient = true;
    return o;
}

}  // namespace

TEST(Med, LenientReadsUnitsIntoMedInfo) {
    // Strict declines (FieldsDeferToPythonWhenUnitsPresent above); lenient must
    // read the field AND report the units, so nothing is lost even though the
    // C++ Mesh has nowhere to put a string.
    std::string p = med_field_fixture();
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid field = h5::open_group(f, "CHA/temperature");
        H5Adelete(field, "UNI");
        H5Adelete(field, "UNT");
        h5::write_attr_string(field, "UNI", "K");
        h5::write_attr_string(field, "UNT", "s");
    }

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info, med_lenient());

    ASSERT_TRUE(out.HasPointData("temperature"));
    for (std::size_t i = 0; i < out.NumPoints(); ++i)
        EXPECT_DOUBLE_EQ(read_double(out.PointData("temperature"), i), static_cast<double>(i));
    ASSERT_EQ(info.mFieldUnits.count("temperature"), 1u);
    EXPECT_EQ(info.mFieldUnits["temperature"].first, "K");
    EXPECT_EQ(info.mFieldUnits["temperature"].second, "s");

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, LenientReadsNonDefaultStepMetadataIntoMedInfo) {
    std::string p = med_field_fixture();
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
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info, med_lenient());
    ASSERT_TRUE(out.HasPointData("temperature"));
    ASSERT_EQ(info.mStepMeta.count("temperature"), 1u);
    EXPECT_EQ(std::get<0>(info.mStepMeta["temperature"]), 7);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, LenientSkipsAnElgaSupportButKeepsTheMesh) {
    // An ELNO/ELGA support is one value per node-within-cell or Gauss point --
    // a 3-D shape the uniform mesh API's (n,)/(n,k) cell_data cannot hold at
    // all. Strict declines the file; lenient drops that one field and still
    // returns the geometry, which is the whole point for a Python-less caller.
    std::string p = med_field_fixture();
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid field = h5::open_group(f, "CHA/temperature");
        std::vector<std::string> steps = h5::group_links(field);
        h5::Hid ts = h5::open_group(field, steps[0]);
        // Rename the NOE support to an ELGA-style one the reader must refuse.
        H5Lmove(ts, "NOE", ts, "NOE.TR3", H5P_DEFAULT, H5P_DEFAULT);
    }

    meshioplusplus::MedInfo strict_info;
    EXPECT_THROW(meshioplusplus::read_med(p, strict_info), meshioplusplus::ReadError);

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info, med_lenient());
    EXPECT_FALSE(out.HasPointData("temperature"));
    EXPECT_EQ(out.NumPoints(), mt::tri_mesh().NumPoints());
    EXPECT_EQ(out.NumCellBlocks(), 1u);
    ASSERT_EQ(info.mSkippedConstructs.size(), 1u);
    EXPECT_NE(info.mSkippedConstructs[0].find("ELNO/ELGA"), std::string::npos);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, LenientSkipsANamedProfile) {
    std::string p = med_field_fixture();
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid noe = h5::open_group(f, "CHA/temperature");
        std::vector<std::string> steps = h5::group_links(noe);
        h5::Hid ts = h5::open_group(noe, steps[0]);
        h5::Hid supp = h5::open_group(ts, "NOE");
        H5Adelete(supp, "PFL");
        h5::write_attr_string(supp, "PFL", "MY_PROFILE");
    }

    meshioplusplus::MedInfo strict_info;
    EXPECT_THROW(meshioplusplus::read_med(p, strict_info), meshioplusplus::ReadError);

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info, med_lenient());
    EXPECT_FALSE(out.HasPointData("temperature"));
    EXPECT_EQ(out.NumCellBlocks(), 1u);
    ASSERT_EQ(info.mSkippedConstructs.size(), 1u);
    EXPECT_NE(info.mSkippedConstructs[0].find("named profile"), std::string::npos);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, TimeStepSelectsOneStepOfAMultiStepField) {
    // A genuinely multi-step field: strict declines, `mTimeStep` picks one
    // WITHOUT needing mLenient (an explicit request, not a fallback), and
    // MedInfo reports every available step's PDT either way.
    std::string p = med_field_fixture();
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid field = h5::open_group(f, "CHA/temperature");
        std::vector<std::string> steps = h5::group_links(field);
        ASSERT_EQ(steps.size(), 1u);
        // A second step, with distinguishable values and its own PDT.
        char key[64];
        std::snprintf(key, sizeof(key), "%020lld%020lld", 2LL, -1LL);
        h5::Hid ts2 = h5::create_group(field, key);
        h5::write_attr_int(ts2, "NDT", 2);
        h5::write_attr_int(ts2, "NOR", -1);
        {
            h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
            h5::Hid a(H5Acreate2(ts2, "PDT", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT),
                      H5Aclose);
            double v = 2.5;
            H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
        }
        h5::write_attr_int(ts2, "RDT", -1);
        h5::write_attr_int(ts2, "ROR", -1);
        h5::Hid supp = h5::create_group(ts2, "NOE");
        h5::write_attr_string(supp, "GAU", "");
        h5::write_attr_string(supp, "PFL", "MED_NO_PROFILE_INTERNAL");
        h5::Hid prof = h5::create_group(supp, "MED_NO_PROFILE_INTERNAL");
        const std::size_t np = mt::tri_mesh().NumPoints();
        h5::write_attr_int(prof, "NBR", static_cast<std::int64_t>(np));
        h5::write_attr_int(prof, "NGA", 1);
        h5::write_attr_string(prof, "GAU", "");
        meshioplusplus::NDArray vals(meshioplusplus::DType::Float64, {np});
        for (std::size_t i = 0; i < np; ++i)
            vals.As<double>()[i] = 100.0 + static_cast<double>(i);
        h5::write_dataset(prof, "CO", vals);
    }

    // Strict, default step: still a decline (the Python shim keeps falling back).
    meshioplusplus::MedInfo strict_info;
    EXPECT_THROW(meshioplusplus::read_med(p, strict_info), meshioplusplus::ReadError);

    // Explicit mTimeStep = 1 (ResolveTimeStep is 0-based: 0 = first) picks the
    // SECOND step, with no mLenient needed.
    meshioplusplus::ReadOptions second;
    second.mTimeStep = 1;
    meshioplusplus::MedInfo info2;
    meshioplusplus::Mesh out2 = meshioplusplus::read_med(p, info2, second);
    ASSERT_TRUE(out2.HasPointData("temperature"));
    EXPECT_DOUBLE_EQ(read_double(out2.PointData("temperature"), 0), 100.0);
    // Every step's time is reported so a caller can choose before asking.
    ASSERT_EQ(info2.mFieldTimeValues.count("temperature"), 1u);
    ASSERT_EQ(info2.mFieldTimeValues["temperature"].size(), 2u);
    EXPECT_DOUBLE_EQ(info2.mFieldTimeValues["temperature"][1], 2.5);

    // -1 counts from the end, which here is the same second step.
    meshioplusplus::ReadOptions last;
    last.mTimeStep = -1;
    meshioplusplus::MedInfo info_last;
    meshioplusplus::Mesh out_last = meshioplusplus::read_med(p, info_last, last);
    EXPECT_DOUBLE_EQ(read_double(out_last.PointData("temperature"), 0), 100.0);

    // Lenient with the default step takes the first, warning rather than
    // throwing.
    meshioplusplus::MedInfo info1;
    meshioplusplus::Mesh out1 = meshioplusplus::read_med(p, info1, med_lenient());
    ASSERT_TRUE(out1.HasPointData("temperature"));
    EXPECT_DOUBLE_EQ(read_double(out1.PointData("temperature"), 0), 0.0);

    // An out-of-range step is a named error, never a clamp.
    meshioplusplus::ReadOptions too_far;
    too_far.mTimeStep = 99;
    meshioplusplus::MedInfo bad;
    EXPECT_THROW(meshioplusplus::read_med(p, bad, too_far), meshioplusplus::ReadError);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, PartialCellTagsBecomeFamilyZero) {
    // A block the file left FAM off of belongs to no family, and MED spells
    // that as id 0 -- so the reader reports zeros instead of throwing
    // ("partial cell tags handled by Python fallback"), which used to make an
    // ordinary Salome file unreadable wherever there is no Python.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddCellBlock("quad", mt::conn_from({{1, 4, 5, 2}}));
    std::vector<meshioplusplus::NDArray> tags;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        meshioplusplus::NDArray t(meshioplusplus::DType::Int64, {m.Cells(b).NumCells()});
        for (std::size_t c = 0; c < m.Cells(b).NumCells(); ++c)
            t.As<std::int64_t>()[c] = -1;
        tags.push_back(std::move(t));
    }
    m.AddCellData("cell_tags", std::move(tags));
    meshioplusplus::MedInfo win;
    win.mCellTags[-1] = {"grp"};
    win.mCellTagGroups[-1] = "FAM_-1_grp";
    meshioplusplus::write_med(p, m, win);

    // Delete one block's FAM so the file carries it on only one of two.
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid mai =
            h5::open_group(f, "ENS_MAA/mesh/-0000000000000000001-0000000000000000001/MAI");
        h5::Hid qu4 = h5::open_group(mai, "QU4");
        H5Ldelete(qu4, "FAM", H5P_DEFAULT);
    }

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);  // strict: must NOT throw
    ASSERT_TRUE(out.HasCellData("cell_tags"));
    // One array per block, always -- the invariant the old short-list Python
    // behaviour violated.
    ASSERT_EQ(out.CellDataNumBlocks("cell_tags"), out.NumCellBlocks());
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const auto cb = out.Cells(b);
        const meshioplusplus::NDArray& t = out.CellData("cell_tags", b);
        ASSERT_EQ(t.Size(), cb.NumCells()) << cb.Type();
        const std::int64_t expect = cb.Type() == "quad" ? 0 : -1;
        for (std::size_t c = 0; c < t.Size(); ++c)
            EXPECT_EQ(meshioplusplus::detail::read_int(t, c), expect) << cb.Type();
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, Regions) {
    // Point and Cell regions attach on read (med_attach_point_regions/
    // med_attach_cell_regions) and, when the mesh carries no native
    // point_tags/cell_tags, synthesize families on write
    // (med_point_regions_to_tags/med_cell_regions_to_tags).
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tet_mesh();  // 2 tetra, 5 points

    meshioplusplus::NDArray pts(meshioplusplus::DType::Int64, {2});
    pts.As<std::int64_t>()[0] = 0;
    pts.As<std::int64_t>()[1] = 3;
    m.AddRegion(
        meshioplusplus::Region("clamped", meshioplusplus::RegionKind::Point, std::move(pts)));

    meshioplusplus::NDArray cells(meshioplusplus::DType::Int64, {1});
    cells.As<std::int64_t>()[0] = 0;
    m.AddRegion(
        meshioplusplus::Region("solid", meshioplusplus::RegionKind::Cell, std::move(cells)));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasRegion("clamped", meshioplusplus::RegionKind::Point));
    const meshioplusplus::Region& rp =
        out.Region(out.FindRegion("clamped", meshioplusplus::RegionKind::Point));
    ASSERT_EQ(rp.NumEntries(), 2u);
    EXPECT_EQ(rp.Entries()[0], 0);
    EXPECT_EQ(rp.Entries()[1], 3);

    ASSERT_TRUE(out.HasRegion("solid", meshioplusplus::RegionKind::Cell));
    const meshioplusplus::Region& rc =
        out.Region(out.FindRegion("solid", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(rc.NumEntries(), 1u);
    EXPECT_EQ(rc.Entries()[0], 0);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, RegionsSideDroppedWithWarning) {
    // A Side region has no MED equivalent; it must be dropped without
    // preventing Point/Cell regions in the same mesh from being written.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tet_mesh();

    meshioplusplus::NDArray cells(meshioplusplus::DType::Int64, {2});
    cells.As<std::int64_t>()[0] = 0;
    cells.As<std::int64_t>()[1] = 1;
    m.AddRegion(
        meshioplusplus::Region("solid", meshioplusplus::RegionKind::Cell, std::move(cells)));

    meshioplusplus::NDArray sides(meshioplusplus::DType::Int64, {1, 2});
    sides.As<std::int64_t>()[0] = 0;
    sides.As<std::int64_t>()[1] = 1;
    m.AddRegion(meshioplusplus::Region("wall", meshioplusplus::RegionKind::Side, std::move(sides)));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    EXPECT_TRUE(out.HasRegion("solid", meshioplusplus::RegionKind::Cell));
    EXPECT_FALSE(out.HasRegion("wall", meshioplusplus::RegionKind::Side));
    EXPECT_EQ(out.NumRegions(), 1u);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, RegionsPrecedenceNativeTagsWin) {
    // A mesh that already carries point_tags/cell_tags (the MED-read shape)
    // must write exactly that, ignoring any Region -- native data wins.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    const std::size_t ntri = m.Cells(0).NumCells();
    meshioplusplus::NDArray tag(meshioplusplus::DType::Int64, {ntri});
    for (std::size_t i = 0; i < ntri; ++i)
        tag.As<std::int64_t>()[i] = -1;
    m.AddCellData("cell_tags", {std::move(tag)});

    meshioplusplus::NDArray region_cells(meshioplusplus::DType::Int64, {1});
    region_cells.As<std::int64_t>()[0] = 0;
    m.AddRegion(meshioplusplus::Region("should_be_ignored", meshioplusplus::RegionKind::Cell,
                                       std::move(region_cells)));

    meshioplusplus::MedInfo win;
    win.mCellTags[-1] = {"native"};
    win.mCellTagGroups[-1] = "FAM_-1_native";
    meshioplusplus::write_med(p, m, win);

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    EXPECT_TRUE(out.HasRegion("native", meshioplusplus::RegionKind::Cell));
    EXPECT_FALSE(out.HasRegion("should_be_ignored", meshioplusplus::RegionKind::Cell));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, GmshPhysicalGroupsBecomeFamilies) {
    // Writing a mesh through the C++ core used to throw unconditionally on
    // any `gmsh:physical` cell_data ("handled by Python fallback") -- it must
    // now synthesize families exactly as `_ensure_med_families` does in the
    // Python reference: named via `field_data` when available, else
    // "group_<id>", with id 0 meaning "no group".
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}, {1, 4, 5, 2}}));

    // Cell 0 tagged 7 (named "surf" via field_data), cell 1 tagged 0 (no group).
    meshioplusplus::NDArray phys(meshioplusplus::DType::Int64, {2});
    phys.As<std::int64_t>()[0] = 7;
    phys.As<std::int64_t>()[1] = 0;
    m.AddCellData("gmsh:physical", {std::move(phys)});
    meshioplusplus::NDArray tag_dim(meshioplusplus::DType::Int64, {2});
    tag_dim.As<std::int64_t>()[0] = 7;
    tag_dim.As<std::int64_t>()[1] = 2;
    m.AddFieldData("surf", std::move(tag_dim));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasRegion("surf", meshioplusplus::RegionKind::Cell));
    const meshioplusplus::Region& r =
        out.Region(out.FindRegion("surf", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(r.NumEntries(), 1u);
    EXPECT_EQ(r.Entries()[0], 0);
    // Cell 1 (tag 0) must not have earned a group.
    EXPECT_FALSE(out.HasRegion("group_0", meshioplusplus::RegionKind::Cell));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, GmshPhysicalUnnamedIdBecomesGroupId) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();  // 2 triangles
    const std::size_t ntri = m.Cells(0).NumCells();
    meshioplusplus::NDArray phys(meshioplusplus::DType::Int64, {ntri});
    for (std::size_t i = 0; i < ntri; ++i)
        phys.As<std::int64_t>()[i] = 42;  // no field_data name for id 42
    m.AddCellData("gmsh:physical", {std::move(phys)});

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasRegion("group_42", meshioplusplus::RegionKind::Cell));
    const meshioplusplus::Region& r =
        out.Region(out.FindRegion("group_42", meshioplusplus::RegionKind::Cell));
    EXPECT_EQ(r.NumEntries(), ntri);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, GmshPhysicalSkipsIdAlreadyANamedCellRegion) {
    // Mirrors `_ensure_med_families`' "already captured as a named cell_set"
    // rule: a physical id resolving (via field_data) to a name that is
    // *already* a Cell region name is not also added via the raw numeric
    // path, even for a cell the existing region does not itself cover.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}, {1, 4, 5, 2}}));

    // "surf" already exists as a Cell region, covering only cell 1.
    meshioplusplus::NDArray region_cells(meshioplusplus::DType::Int64, {1});
    region_cells.As<std::int64_t>()[0] = 1;
    m.AddRegion(
        meshioplusplus::Region("surf", meshioplusplus::RegionKind::Cell, std::move(region_cells)));

    // Cell 0 is tagged 7, and field_data maps 7 -> "surf" too.
    meshioplusplus::NDArray phys(meshioplusplus::DType::Int64, {2});
    phys.As<std::int64_t>()[0] = 7;
    phys.As<std::int64_t>()[1] = 0;
    m.AddCellData("gmsh:physical", {std::move(phys)});
    meshioplusplus::NDArray tag_dim(meshioplusplus::DType::Int64, {2});
    tag_dim.As<std::int64_t>()[0] = 7;
    tag_dim.As<std::int64_t>()[1] = 2;
    m.AddFieldData("surf", std::move(tag_dim));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasRegion("surf", meshioplusplus::RegionKind::Cell));
    const meshioplusplus::Region& r =
        out.Region(out.FindRegion("surf", meshioplusplus::RegionKind::Cell));
    // Cell 0 must NOT have been added -- only the original member (cell 1).
    ASSERT_EQ(r.NumEntries(), 1u);
    EXPECT_EQ(r.Entries()[0], 1);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, SameTypeBlocksAreConsolidated) {
    // MSH 4.1's canonical structure is one cell block per entity -- two
    // `triangle` blocks must consolidate into ONE `TR3` MED section rather
    // than throwing (the pre-v9.8.0 "MED files cannot have two sections of
    // the same cell type" rejection).
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_quad_mesh();  // triangle(2), quad(1), triangle(1)

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        h5::Hid mai =
            h5::open_group(f, "ENS_MAA/mesh/-0000000000000000001-0000000000000000001/MAI");
        std::vector<std::string> sections = h5::group_links(mai);
        std::sort(sections.begin(), sections.end());
        EXPECT_EQ(sections, (std::vector<std::string>{"QU4", "TR3"}));
    }

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    ASSERT_EQ(out.NumCellBlocks(), 2u);
    std::size_t n_tri = 0, n_quad = 0;
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        if (out.Cells(b).Type() == "triangle")
            n_tri = out.Cells(b).NumCells();
        if (out.Cells(b).Type() == "quad")
            n_quad = out.Cells(b).NumCells();
    }
    EXPECT_EQ(n_tri, 3u);
    EXPECT_EQ(n_quad, 1u);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, VectorFieldRoundTrip) {
    // A genuine (n, 3) point field and an (n, 6) cell field. Before v9.9.0 a
    // caller that lost its component count (notably the WASM boundary) handed
    // in a flattened (3n,) array, which wrote NCO=1/NBR=3n and produced a file
    // this reader rejects -- see Med.MisShapedFieldIsRejectedOnWrite.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    const std::size_t np = m.NumPoints();
    const std::size_t nc = m.Cells(0).NumCells();

    meshioplusplus::NDArray vel(meshioplusplus::DType::Float64, {np, 3});
    for (std::size_t i = 0; i < np * 3; ++i)
        vel.As<double>()[i] = static_cast<double>(i) + 0.25;
    m.AddPointData("velocity", std::move(vel));

    meshioplusplus::NDArray stress(meshioplusplus::DType::Float64, {nc, 6});
    for (std::size_t i = 0; i < nc * 6; ++i)
        stress.As<double>()[i] = 100.0 + static_cast<double>(i);
    m.AddCellData("stress", {std::move(stress)});

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasPointData("velocity"));
    const meshioplusplus::NDArray& v = out.PointData("velocity");
    ASSERT_EQ(v.Shape().size(), 2u);
    EXPECT_EQ(v.Shape()[0], np);
    EXPECT_EQ(v.Shape()[1], 3u);
    for (std::size_t i = 0; i < np * 3; ++i)
        EXPECT_DOUBLE_EQ(read_double(v, i), static_cast<double>(i) + 0.25);

    ASSERT_TRUE(out.HasCellData("stress"));
    const meshioplusplus::NDArray& s = out.CellData("stress", 0);
    ASSERT_EQ(s.Shape().size(), 2u);
    EXPECT_EQ(s.Shape()[0], nc);
    EXPECT_EQ(s.Shape()[1], 6u);
    for (std::size_t i = 0; i < nc * 6; ++i)
        EXPECT_DOUBLE_EQ(read_double(s, i), 100.0 + static_cast<double>(i));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, MisShapedFieldIsRejectedOnWrite) {
    // The write-side guard added in v9.9.0. A field's row count IS its entity
    // count; without this the writer emitted NBR=3n against n points and the
    // failure only surfaced on the way back in, from a different function.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray flat(meshioplusplus::DType::Float64, {m.NumPoints() * 3});
    for (std::size_t i = 0; i < flat.Size(); ++i)
        flat.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("velocity", std::move(flat));

    EXPECT_THROW(meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{}),
                 meshioplusplus::WriteError);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, VectorFieldWritesOneNomSlotPerComponent) {
    // MED's NOM is 16 characters PER COMPONENT. A k>1 field with no names of
    // its own gets MED's default `V1..Vk` spelling; a scalar keeps the single
    // blank 16-char slot it always had.
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray vel(meshioplusplus::DType::Float64, {m.NumPoints(), 3});
    for (std::size_t i = 0; i < vel.Size(); ++i)
        vel.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("velocity", std::move(vel));
    meshioplusplus::NDArray temp(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t i = 0; i < temp.Size(); ++i)
        temp.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("temperature", std::move(temp));

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        h5::Hid vfield = h5::open_group(f, "CHA/velocity");
        // read_attr_string trims trailing blanks, so compare on the tail name.
        const std::string nom = h5::read_attr_string(vfield, "NOM");
        EXPECT_EQ(nom.substr(0, 2), "V1");
        EXPECT_NE(nom.find("V2"), std::string::npos);
        EXPECT_NE(nom.find("V3"), std::string::npos);
        h5::Hid sfield = h5::open_group(f, "CHA/temperature");
        EXPECT_TRUE(h5::read_attr_string(sfield, "NOM").empty());  // one blank slot
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, GlobalNumbers) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::Mesh m = mt::tri_mesh();

    meshioplusplus::NDArray point_num(meshioplusplus::DType::Int32, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        point_num.As<std::int32_t>()[i] = static_cast<std::int32_t>(10 + i);
    m.AddPointData("med:num", std::move(point_num));

    const std::size_t ntri = m.Cells(0).NumCells();
    meshioplusplus::NDArray cell_num(meshioplusplus::DType::Int32, {ntri});
    for (std::size_t i = 0; i < ntri; ++i)
        cell_num.As<std::int32_t>()[i] = static_cast<std::int32_t>(100 + i);
    m.AddCellData("med:num", {std::move(cell_num)});

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);

    ASSERT_TRUE(out.HasPointData("med:num"));
    const meshioplusplus::NDArray& pn = out.PointData("med:num");
    ASSERT_EQ(pn.Size(), m.NumPoints());
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(pn, i), 10 + static_cast<std::int64_t>(i));

    ASSERT_TRUE(out.HasCellData("med:num"));
    const meshioplusplus::NDArray& cn = out.CellData("med:num", 0);
    ASSERT_EQ(cn.Size(), ntri);
    for (std::size_t i = 0; i < ntri; ++i)
        EXPECT_EQ(meshioplusplus::detail::read_int(cn, i), 100 + static_cast<std::int64_t>(i));

    // med:num must not have leaked into CHA as an ordinary field.
    h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    EXPECT_FALSE(h5::exists(f, "CHA"));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, NoGlobalNumbersWritesNoNumDataset) {
    std::string p = mt::temp_path(".med");
    meshioplusplus::write_med(p, mt::tri_mesh(), meshioplusplus::MedInfo{});

    h5::SilenceErrors silence;
    h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    h5::Hid noe = h5::open_group(f, "ENS_MAA/mesh/-0000000000000000001-0000000000000000001/NOE");
    EXPECT_FALSE(h5::exists(noe, "NUM"));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, NewerVersionRejected) {
    // A file whose INFOS_GENERALES declares a MED major version newer than
    // 4 must be rejected with a named diagnosis rather than an obscure
    // structural error further down.
    std::string p = mt::temp_path(".med");
    meshioplusplus::write_med(p, mt::tri_mesh(), meshioplusplus::MedInfo{});

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid infos = h5::open_group(f, "INFOS_GENERALES");
        h5::Hid a(H5Aopen(infos, "MAJ", H5P_DEFAULT), H5Aclose);
        std::int64_t five = 5;
        H5Awrite(a, H5T_NATIVE_INT64, &five);
    }

    meshioplusplus::MedInfo info;
    try {
        meshioplusplus::read_med(p, info);
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("newer than the MED 4.1"), std::string::npos);
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, OlderVersionStillReads) {
    // Only a *newer* major version is rejected -- MED 3.x files (like the
    // Python test suite's cylinder.med fixture) must keep reading.
    std::string p = mt::temp_path(".med");
    meshioplusplus::write_med(p, mt::tri_mesh(), meshioplusplus::MedInfo{});

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid infos = h5::open_group(f, "INFOS_GENERALES");
        h5::Hid a(H5Aopen(infos, "MAJ", H5P_DEFAULT), H5Aclose);
        std::int64_t three = 3;
        H5Awrite(a, H5T_NATIVE_INT64, &three);
    }

    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    EXPECT_EQ(out.NumPoints(), mt::tri_mesh().NumPoints());

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

TEST(Med, IsAnOptionsAwareReader) {
    // Registering MED in registry_readers_ex() is what lets ANY flat binding
    // pass `lenient`/`time_step` at all; before it this was false and there was
    // nowhere for either option to go -- which, with no Python fallback in a
    // WASM/C-API/Fortran build, meant a real Salome/Code_Aster file was simply
    // unreadable there. Exodus's twin assertion is in test_netcdf_formats.cpp.
    EXPECT_TRUE(meshioplusplus::registry_reader_supports_options("med"));
}

TEST(Med, PolyhedronPoeRoundTrip) {
    // MED_POLYHEDRON (POE) needs THREE 1-based arrays where a polygon needs
    // two: NOD, INN (face -> NOD) and IND (cell -> face). Before v9.19.0
    // neither the C++ nor the Python path had a POE entry at all -- contrary to
    // what doc/formats/med.md claimed, which said polyhedra were "Python-only".
    std::string p = mt::temp_path("_poe.med");
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}}});

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    const auto cb = out.Cells(0);
    EXPECT_TRUE(cb.IsPolyhedron());
    EXPECT_EQ(cb.Type(), "polyhedron8");
    ASSERT_EQ(cb.NumCells(), 1u);
    EXPECT_EQ(cb.NumFaces(0), 6u);
    // Geometry, not just arity: six faces of the right size prove nothing
    // about whether the node ids landed where they belong.
    EXPECT_NEAR(meshioplusplus::compute_stats(out).mUnsignedVolume, 1.0, 1e-12);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Med, MixedPolyhedronNodeCountsShareOnePoeSection) {
    // MED holds ONE section per type inside a MAI group, so polyhedron4 and
    // polyhedron5 must consolidate into the SAME POE. Grouping on the exact
    // meshio type string would try to create POE twice and fail at group
    // creation -- which is why the writer canonicalises every polyhedron<N>.
    std::string p = mt::temp_path("_poe_mixed.med");
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    m.AddPolyhedronBlock("polyhedron4", {{{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}});
    m.AddPolyhedronBlock("polyhedron5", {{{1, 2, 3}, {1, 3, 4}, {2, 1, 4}, {3, 2, 4}}});

    meshioplusplus::write_med(p, m, meshioplusplus::MedInfo{});
    meshioplusplus::MedInfo info;
    meshioplusplus::Mesh out = meshioplusplus::read_med(p, info);
    std::size_t total = 0;
    for (const auto cb : out.CellRange()) {
        EXPECT_TRUE(cb.IsPolyhedron());
        total += cb.NumCells();
    }
    EXPECT_EQ(total, 2u);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
