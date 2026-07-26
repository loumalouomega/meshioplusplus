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

// System includes
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/formats/xdmf_time_series.hpp"
#include "meshioplusplus/read_options.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadOptions;
using meshioplusplus::XdmfTimeSeriesWriter;

/** @brief Deep copy, so a step mesh never aliases the fixture's storage. */
NDArray xts_copy(const NDArray& rA) {
    NDArray out(rA.Dtype(), rA.Shape());
    std::memcpy(out.Data(), rA.Data(), rA.Nbytes());
    return out;
}

// A step's data: `T` scales with the step index so a reader taking the wrong
// step is visible in the values, not merely in a count.
Mesh xts_step_mesh(const Mesh& rBase, int Step) {
    Mesh m;
    m.AssignPoints(xts_copy(rBase.Points()));
    for (const auto cb : rBase.CellRange())
        m.AddCellBlock(cb.Type(), xts_copy(cb.Conn()));

    NDArray pd(meshioplusplus::DType::Float64, {rBase.NumPoints()});
    for (std::size_t i = 0; i < rBase.NumPoints(); ++i)
        pd.As<double>()[i] = 100.0 * Step + static_cast<double>(i);
    m.AddPointData("T", std::move(pd));

    std::vector<NDArray> cd;
    for (const auto cb : rBase.CellRange()) {
        NDArray a(meshioplusplus::DType::Float64, {cb.NumCells()});
        for (std::size_t c = 0; c < cb.NumCells(); ++c)
            a.As<double>()[c] = 10.0 * Step + static_cast<double>(c);
        cd.push_back(std::move(a));
    }
    m.AddCellData("rho", std::move(cd));
    return m;
}

// Write `nsteps` steps at t = 0.5*k and hand back the .xdmf path.
std::string xts_write_series(const Mesh& rBase, int NSteps, const std::string& rDataFormat) {
    const std::string path = mt::temp_path(".xdmf");
    XdmfTimeSeriesWriter w(path, rDataFormat);
    w.WritePointsCells(rBase);
    for (int k = 0; k < NSteps; ++k)
        w.WriteData(0.5 * k, xts_step_mesh(rBase, k));
    EXPECT_EQ(w.NumSteps(), static_cast<std::size_t>(NSteps));
    w.Finalize();
    EXPECT_TRUE(w.Finalized());
    return path;
}

void xts_cleanup(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
    std::string base = rPath.substr(0, rPath.find_last_of('.'));
    std::filesystem::remove(base + ".h5", ec);
    for (int i = 0; i < 64; ++i)
        std::filesystem::remove(base + std::to_string(i) + ".bin", ec);
}

// Every step of the series must read back with exactly the values it was
// written with, and the geometry must be the shared static grid every time.
void xts_check_series(const std::string& rPath, const Mesh& rBase, int NSteps) {
    for (int k = 0; k < NSteps; ++k) {
        ReadOptions opts;
        opts.mTimeStep = k;
        Mesh m = meshioplusplus::read_xdmf(rPath, opts);
        mt::expect_mesh_eq(rBase, m);

        ASSERT_TRUE(m.HasPointData("T")) << "step " << k;
        const NDArray& pd = m.PointData("T");
        ASSERT_EQ(pd.Size(), rBase.NumPoints());
        for (std::size_t i = 0; i < pd.Size(); ++i)
            EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pd, i),
                             100.0 * k + static_cast<double>(i))
                << "step " << k << " point " << i;

        ASSERT_EQ(m.CellDataNumBlocks("rho"), m.NumCellBlocks()) << "step " << k;
        for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
            const NDArray& cd = m.CellData("rho", b);
            for (std::size_t c = 0; c < cd.Size(); ++c)
                EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(cd, c),
                                 10.0 * k + static_cast<double>(c))
                    << "step " << k << " block " << b << " cell " << c;
        }
    }
}

}  // namespace

TEST(XdmfTimeSeries, XmlRoundTrip) {
    const Mesh base = mt::tet_mesh();
    const std::string path = xts_write_series(base, 3, "XML");
    xts_check_series(path, base, 3);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, BinaryRoundTrip) {
    const Mesh base = mt::tri_mesh();
    const std::string path = xts_write_series(base, 3, "Binary");
    xts_check_series(path, base, 3);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, MixedTopologyRoundTrip) {
    // A Mixed topology in a series: the static grid packs every block once and
    // cell data is concatenated across blocks per step.
    const Mesh base = mt::tri_quad_mesh();
    const std::string path = xts_write_series(base, 2, "XML");
    xts_check_series(path, base, 2);
    xts_cleanup(path);
}

// `mTimeValues` is what makes a step request checkable before it is issued.
TEST(XdmfTimeSeries, MetadataReportsTimeValues) {
    const Mesh base = mt::tet_mesh();
    const std::string path = xts_write_series(base, 4, "XML");

    meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    ASSERT_EQ(meta.mTimeValues.size(), 4u);
    for (int k = 0; k < 4; ++k)
        EXPECT_DOUBLE_EQ(meta.mTimeValues[static_cast<std::size_t>(k)], 0.5 * k);

    // The static grid's shape still comes from the .xdmf attributes alone.
    EXPECT_EQ(meta.mNumPoints, base.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "tetra");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 2u);
    EXPECT_FALSE(meta.mHasBBox);

    // Attribute names come off the steps, not the static grid.
    ASSERT_EQ(meta.mPointDataNames.size(), 1u);
    EXPECT_EQ(meta.mPointDataNames[0], "T");
    ASSERT_EQ(meta.mCellDataNames.size(), 1u);
    EXPECT_EQ(meta.mCellDataNames[0], "rho");

    xts_cleanup(path);
}

// Negative indices count from the end, and out of range names the count rather
// than silently clamping.
TEST(XdmfTimeSeries, NegativeAndOutOfRangeSteps) {
    const Mesh base = mt::tet_mesh();
    const std::string path = xts_write_series(base, 3, "XML");

    ReadOptions last;
    last.mTimeStep = -1;
    Mesh m = meshioplusplus::read_xdmf(path, last);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.PointData("T"), 0), 200.0);

    ReadOptions bad;
    bad.mTimeStep = 7;
    EXPECT_THROW(meshioplusplus::read_xdmf(path, bad), meshioplusplus::ReadError);

    xts_cleanup(path);
}

// A default-constructed ReadOptions still means "the first step".
TEST(XdmfTimeSeries, DefaultReadTakesFirstStep) {
    const Mesh base = mt::tet_mesh();
    const std::string path = xts_write_series(base, 3, "XML");
    Mesh m = meshioplusplus::read_xdmf(path);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.PointData("T"), 0), 0.0);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, WriteDataBeforeMeshThrows) {
    const std::string path = mt::temp_path(".xdmf");
    XdmfTimeSeriesWriter w(path, "XML");
    EXPECT_THROW(w.WriteData(0.0, mt::tet_mesh()), meshioplusplus::WriteError);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, WritePointsCellsTwiceThrows) {
    const std::string path = mt::temp_path(".xdmf");
    XdmfTimeSeriesWriter w(path, "XML");
    w.WritePointsCells(mt::tet_mesh());
    EXPECT_THROW(w.WritePointsCells(mt::tet_mesh()), meshioplusplus::WriteError);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, UnknownDataFormatThrows) {
    EXPECT_THROW(XdmfTimeSeriesWriter(mt::temp_path(".xdmf"), "Zarr"), meshioplusplus::WriteError);
}

// The destructor finalizes, so a series written in a scope is complete on exit.
TEST(XdmfTimeSeries, DestructorFinalizes) {
    const Mesh base = mt::tet_mesh();
    const std::string path = mt::temp_path(".xdmf");
    {
        XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        w.WriteData(1.25, xts_step_mesh(base, 1));
    }
    ASSERT_TRUE(std::filesystem::exists(path));
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    ASSERT_EQ(meta.mTimeValues.size(), 1u);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[0], 1.25);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, FinalizeIsIdempotent) {
    const std::string path = mt::temp_path(".xdmf");
    XdmfTimeSeriesWriter w(path, "XML");
    w.WritePointsCells(mt::tet_mesh());
    w.Finalize();
    EXPECT_NO_THROW(w.Finalize());
    // Writing after finalizing is an error, not a silent no-op.
    EXPECT_THROW(w.WriteData(0.0, mt::tet_mesh()), meshioplusplus::WriteError);
    xts_cleanup(path);
}

// A series carrying no data at all still has one step per WriteData call.
TEST(XdmfTimeSeries, StepsWithoutData) {
    const Mesh base = mt::tet_mesh();
    const std::string path = mt::temp_path(".xdmf");
    {
        XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        w.WriteData(0.0, base);
        w.WriteData(1.0, base);
    }
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    EXPECT_EQ(meta.mTimeValues.size(), 2u);
    EXPECT_TRUE(meta.mPointDataNames.empty());
    Mesh m = meshioplusplus::read_xdmf(path);
    mt::expect_mesh_eq(base, m);
    xts_cleanup(path);
}

// A plain single-grid .xdmf must keep reporting no time values at all -- the
// temporal support must not invent a series where the file has none.
TEST(XdmfTimeSeries, SingleGridHasNoTimeValues) {
    const std::string path = mt::temp_path(".xdmf");
    meshioplusplus::write_xdmf(path, mt::tet_mesh(), "XML", -1);
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    EXPECT_TRUE(meta.mTimeValues.empty());
    xts_cleanup(path);
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(XdmfTimeSeries, HdfRoundTrip) {
    const Mesh base = mt::tet_mesh();
    const std::string path = xts_write_series(base, 3, "HDF");
    // One companion .h5 beside the .xdmf holding every step's data.
    EXPECT_TRUE(std::filesystem::exists(path.substr(0, path.find_last_of('.')) + ".h5"));
    xts_check_series(path, base, 3);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, HdfGzipRoundTrip) {
    const Mesh base = mt::hex_mesh();
    const std::string path = mt::temp_path(".xdmf");
    {
        XdmfTimeSeriesWriter w(path, "HDF", 4);
        w.WritePointsCells(base);
        for (int k = 0; k < 3; ++k)
            w.WriteData(0.5 * k, xts_step_mesh(base, k));
    }
    xts_check_series(path, base, 3);
    xts_cleanup(path);
}
#else
// Compiled-out contract: the symbol exists and throws by name rather than
// being absent from the library.
TEST(XdmfTimeSeries, HdfThrowsWithoutHdf5) {
    EXPECT_THROW(XdmfTimeSeriesWriter(mt::temp_path(".xdmf"), "HDF"), meshioplusplus::WriteError);
}
#endif
