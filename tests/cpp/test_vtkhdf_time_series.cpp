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

#include <cstdint>
#include <filesystem>

#include <hdf5.h>

#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/formats/vtkhdf_time_series.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace h5 = meshioplusplus::h5;
using meshioplusplus::VtkhdfSeriesMode;
using meshioplusplus::VtkhdfTimeSeriesWriter;
using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

namespace {

using I64Vec = std::vector<std::int64_t>;

struct SeriesTempFile {
    std::string mPath = mt::temp_path(".vtkhdf");
    ~SeriesTempFile() {
        std::error_code ec;
        std::filesystem::remove(mPath, ec);
    }
};

// Two tetra sharing a face; the arrays below carry the step value so a read-back can tell steps
// apart.
mt::Mesh series_grid() {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    mt::NDArray conn = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {2, 4});
    const std::int64_t rows[8] = {0, 1, 2, 3, 1, 2, 3, 4};
    std::copy(rows, rows + 8, conn.As<std::int64_t>());
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

mt::Mesh series_step(double V, bool WithField = false) {
    mt::Mesh m = series_grid();
    m.AddPointData("u", mt::data_array({V, V, V, V, V}));
    m.AddPointData("vel", mt::data_array({V, 0, 0, V, 0, 0, V, 0, 0, V, 0, 0, V, 0, 0}, 3));
    m.AddCellData("p", {mt::data_array({10 * V, 10 * V + 1})});
    if (WithField)
        m.AddFieldData("gain", mt::data_array({V, 2 * V, 3 * V}));
    return m;
}

std::vector<double> doubles(const mt::NDArray& rA) {
    std::vector<double> out(rA.Size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = read_double(rA, i);
    return out;
}

meshioplusplus::ReadOptions step(int K) {
    meshioplusplus::ReadOptions o;
    o.mTimeStep = K;
    return o;
}

}  // namespace

TEST(VtkhdfTimeSeries, WritesGeometryOnceAndOneStepPerCall) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        EXPECT_TRUE(w.AutoFlush());  // unlike XDMF's: a flush here is cheap
        w.WritePointsCells(series_grid());
        for (int k = 0; k < 4; ++k) {
            w.WriteData(0.5 * k, series_step(k, true));
            EXPECT_EQ(w.NumSteps(), static_cast<std::size_t>(k + 1));
        }
        EXPECT_FALSE(w.Finalized());
        w.Finalize();
        EXPECT_TRUE(w.Finalized());
        w.Finalize();  // idempotent
    }
    for (int k = 0; k < 4; ++k) {
        mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, step(k));
        EXPECT_EQ(doubles(m.PointData("u")), std::vector<double>(5, k)) << k;
        EXPECT_EQ(doubles(m.CellData("p", 0)), (std::vector<double>{10.0 * k, 10.0 * k + 1})) << k;
        EXPECT_EQ(doubles(m.FieldData("gain")), (std::vector<double>{1.0 * k, 2.0 * k, 3.0 * k}))
            << k;
        EXPECT_EQ(doubles(m.FieldData(meshioplusplus::kSequenceTimeKey)),
                  std::vector<double>{0.5 * k})
            << k;
        EXPECT_EQ(m.NumPoints(), 5u);
    }
    EXPECT_EQ(doubles(meshioplusplus::read_vtkhdf(f.mPath, step(-1)).PointData("u")),
              std::vector<double>(5, 3));
}

TEST(VtkhdfTimeSeries, OnDiskShapeIsGeometryOnceAndChunkedFields) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath, 4);
        w.WritePointsCells(series_grid());
        for (int k = 0; k < 3; ++k)
            w.WriteData(k, series_step(k, true));
    }
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    EXPECT_EQ(h5::read_attr_string(g, "Type"), "UnstructuredGrid");
    EXPECT_EQ(h5::read_attr_int_array(g, "Version"), (I64Vec{2, 0}));
    // geometry x 1
    EXPECT_EQ(h5::dataset_num_rows(g, "Points"), 5u);
    EXPECT_EQ(h5::dataset_num_rows(g, "NumberOfPoints"), 1u);
    // fields x N, appended, and therefore chunked (an unlimited dimension requires it)
    EXPECT_EQ(h5::dataset_num_rows(g, "PointData/u"), 15u);
    EXPECT_EQ(h5::dataset_num_rows(g, "PointData/vel"), 15u);
    EXPECT_EQ(h5::dataset_num_rows(g, "CellData/p"), 6u);
    for (const char* name : {"PointData/u", "CellData/p", "Steps/Values", "FieldData/gain"}) {
        h5::Hid d(H5Dopen2(g, name, H5P_DEFAULT), H5Dclose);
        h5::Hid dcpl(H5Dget_create_plist(d), H5Pclose);
        EXPECT_EQ(H5Pget_layout(dcpl), H5D_CHUNKED) << name;
    }
    h5::Hid steps = h5::open_group(g, "Steps");
    EXPECT_EQ(h5::read_attr_int(steps, "NSteps"), 3);
    EXPECT_EQ(doubles(h5::read_dataset(steps, "Values")), (std::vector<double>{0, 1, 2}));
    for (const char* name : {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"})
        EXPECT_EQ(h5::read_dataset(steps, name).Size(), 3u) << name;
    const mt::NDArray po = h5::read_dataset(steps, "PointDataOffsets/u");
    EXPECT_EQ(read_int(po, 0) + read_int(po, 1) + read_int(po, 2), 0 + 5 + 10);
    EXPECT_EQ(h5::dataset_shape(steps, "FieldDataSizes/gain"), (std::vector<std::size_t>{3, 2}));
}

TEST(VtkhdfTimeSeries, EveryStepIsDurableWithoutFinalize) {
    SeriesTempFile f;
    VtkhdfTimeSeriesWriter w(f.mPath);
    w.WritePointsCells(series_grid());
    w.WriteData(0.0, series_step(1));
    w.WriteData(0.5, series_step(2));
    // A crash here would leave a readable file covering both completed steps.
    mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, step(-1));
    EXPECT_EQ(doubles(m.PointData("u")), std::vector<double>(5, 2));
}

TEST(VtkhdfTimeSeries, NamedArrayOverload) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        w.WritePointsCells(series_grid());
        for (int k = 0; k < 2; ++k)
            w.WriteData(
                k, {{"u", 1, std::vector<double>(5, k)}, {"vel", 3, std::vector<double>(15, k)}},
                {{"p", 1, std::vector<double>{1.0 * k, 2.0 * k}}});
    }
    mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, step(1));
    EXPECT_EQ(doubles(m.PointData("u")), std::vector<double>(5, 1));
    EXPECT_EQ(m.PointData("vel").Shape(), (std::vector<std::size_t>{5, 3}));
    EXPECT_EQ(doubles(m.CellData("p", 0)), (std::vector<double>{1, 2}));
}

TEST(VtkhdfTimeSeries, NamedArrayLengthMismatchIsRefused) {
    SeriesTempFile f;
    VtkhdfTimeSeriesWriter w(f.mPath);
    w.WritePointsCells(series_grid());
    EXPECT_THROW(w.WriteData(0.0, {{"u", 1, std::vector<double>(4, 0)}}),
                 meshioplusplus::WriteError);
    EXPECT_THROW(w.WriteData(0.0, {}, {{"p", 1, std::vector<double>(3, 0)}}),
                 meshioplusplus::WriteError);
}

TEST(VtkhdfTimeSeries, TheNameSetIsFixedAtTheFirstStep) {
    SeriesTempFile f;
    VtkhdfTimeSeriesWriter w(f.mPath);
    w.WritePointsCells(series_grid());
    w.WriteData(0.0, series_step(0));
    mt::Mesh extra = series_step(1);
    extra.AddPointData("brand_new", mt::data_array({1, 2, 3, 4, 5}));
    try {
        w.WriteData(1.0, extra);
        FAIL() << "a new array name must be refused";
    } catch (const meshioplusplus::WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("brand_new"), std::string::npos);
    }
    mt::Mesh fewer = series_grid();
    fewer.AddPointData("u", mt::data_array({1, 1, 1, 1, 1}));
    try {
        w.WriteData(1.0, fewer);
        FAIL() << "a dropped array must be refused";
    } catch (const meshioplusplus::WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("vel"), std::string::npos);
    }
    EXPECT_EQ(w.NumSteps(), 1u);  // neither refused step landed
}

TEST(VtkhdfTimeSeries, StepShapeMustMatchTheGrid) {
    SeriesTempFile f;
    VtkhdfTimeSeriesWriter w(f.mPath);
    w.WritePointsCells(series_grid());
    mt::Mesh bad = series_grid();
    bad.AddPointData("u", mt::data_array({1, 2, 3}));
    EXPECT_THROW(w.WriteData(0.0, bad), meshioplusplus::WriteError);
}

TEST(VtkhdfTimeSeries, OrderOfCallsIsEnforced) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        EXPECT_THROW(w.WriteData(0.0, series_step(0)), meshioplusplus::WriteError);  // no grid yet
        w.WritePointsCells(series_grid());
        EXPECT_THROW(w.WritePointsCells(series_grid()), meshioplusplus::WriteError);  // twice
        w.WriteData(0.0, series_step(0));
        w.Finalize();
        EXPECT_THROW(w.WriteData(1.0, series_step(1)), meshioplusplus::WriteError);  // finalized
    }
}

TEST(VtkhdfTimeSeries, MovedFromWriterAnswersAsFinishedAndRefusesWrites) {
    SeriesTempFile f;
    VtkhdfTimeSeriesWriter a(f.mPath);
    a.WritePointsCells(series_grid());
    a.WriteData(0.0, series_step(0));
    VtkhdfTimeSeriesWriter b(std::move(a));
    EXPECT_TRUE(a.Finalized());  // NOLINT(bugprone-use-after-move): the documented contract
    EXPECT_EQ(a.NumSteps(), 0u);
    EXPECT_FALSE(a.AutoFlush());
    EXPECT_NO_THROW(a.Flush());
    EXPECT_NO_THROW(a.Finalize());
    EXPECT_NO_THROW(a.SetAutoFlush(true));
    EXPECT_THROW(a.WritePointsCells(series_grid()), meshioplusplus::WriteError);
    EXPECT_THROW(a.WriteData(1.0, series_step(1)), meshioplusplus::WriteError);
    EXPECT_THROW(a.WriteData(1.0, {}, {}), meshioplusplus::WriteError);
    EXPECT_EQ(b.NumSteps(), 1u);
    b.WriteData(1.0, series_step(1));  // the moved-to writer carries on
    EXPECT_EQ(b.NumSteps(), 2u);
}

TEST(VtkhdfTimeSeries, AppendContinuesAnExistingSeries) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        w.WritePointsCells(series_grid());
        w.WriteData(0.0, series_step(0, true));
        w.WriteData(1.0, series_step(1, true));
    }
    {
        VtkhdfTimeSeriesWriter w(f.mPath, -1, VtkhdfSeriesMode::Append);
        EXPECT_EQ(w.NumSteps(), 2u);
        w.WritePointsCells(
            series_grid());  // a driver may call this unconditionally: it only checks the grid
        w.WriteData(2.0, series_step(2, true));
        EXPECT_EQ(w.NumSteps(), 3u);
        mt::Mesh wrong = mt::tri_mesh();
        EXPECT_THROW(w.WritePointsCells(wrong), meshioplusplus::WriteError);
    }
    for (int k = 0; k < 3; ++k) {
        mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, step(k));
        EXPECT_EQ(doubles(m.PointData("u")), std::vector<double>(5, k)) << k;
        EXPECT_EQ(doubles(m.FieldData("gain")), (std::vector<double>{1.0 * k, 2.0 * k, 3.0 * k}))
            << k;
    }
    // appending to a path that does not exist is exactly Truncate
    SeriesTempFile fresh;
    VtkhdfTimeSeriesWriter w(fresh.mPath, -1, VtkhdfSeriesMode::Append);
    EXPECT_EQ(w.NumSteps(), 0u);
    w.WritePointsCells(series_grid());
    w.WriteData(0.0, series_step(0));
    EXPECT_EQ(w.NumSteps(), 1u);
}

TEST(VtkhdfTimeSeries, AppendRefusesAFileItCannotContinue) {
    SeriesTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, series_step(0));  // a static file: no Steps group
    EXPECT_THROW(VtkhdfTimeSeriesWriter(f.mPath, -1, VtkhdfSeriesMode::Append),
                 meshioplusplus::WriteError);
}

TEST(VtkhdfTimeSeries, PolyhedralGridsCarryTheExtraStepTables) {
    // Built afresh per use: a KRATOS mesh cannot be copied.
    auto make_grid = [] {
        mt::Mesh grid;
        grid.AssignPoints(mt::points_from({{0, 0, 0},
                                           {1, 0, 0},
                                           {1, 1, 0},
                                           {0, 1, 0},
                                           {0, 0, 1},
                                           {1, 0, 1},
                                           {1, 1, 1},
                                           {0, 1, 1}}));
        grid.AddPolyhedronBlock(
            "polyhedron8",
            {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}}});
        return grid;
    };
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        w.WritePointsCells(make_grid());
        for (int k = 0; k < 2; ++k) {
            mt::Mesh m = make_grid();
            m.AddPointData("u", mt::data_array(std::vector<double>(8, k)));
            w.WriteData(k, m);
        }
    }
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    EXPECT_EQ(h5::read_attr_int_array(g, "Version"), (I64Vec{2, 5}));
    h5::Hid steps = h5::open_group(g, "Steps");
    for (const char* name :
         {"FaceConnectivityOffsets", "FaceOffsetsOffsets", "PolyhedronToFaceIdOffsets"})
        EXPECT_EQ(h5::read_dataset(steps, name).Size(), 2u) << name;
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath, step(1));
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "polyhedron8");
    EXPECT_EQ(doubles(back.PointData("u")), std::vector<double>(8, 1));
}

TEST(VtkhdfTimeSeries, ReaderReportsEveryTimeValueWithoutAFullRead) {
    SeriesTempFile f;
    {
        VtkhdfTimeSeriesWriter w(f.mPath);
        w.WritePointsCells(series_grid());
        for (double t : {0.0, 0.25, 0.75})
            w.WriteData(t, series_step(t));
    }
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_vtkhdf_metadata(f.mPath, {});
    EXPECT_FALSE(meta.mFellBackToFullRead);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.0, 0.25, 0.75}));
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
