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
#include <fstream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#ifndef _WIN32
#include <sys/wait.h>
#include <unistd.h>
#endif

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

// ---------------------------------------------------------------------------
// v9.1.0: Flush, append, and the NamedArray overload
// ---------------------------------------------------------------------------

TEST(XdmfTimeSeries, FlushMakesAPartialSeriesReadable) {
    for (const char* p_format : {"XML", "Binary", "HDF"}) {
        const std::string path = mt::temp_path(".xdmf");
        {
            meshioplusplus::XdmfTimeSeriesWriter w(path, p_format);
            w.WritePointsCells(mt::tri_mesh());
            for (int k = 0; k < 3; ++k) {
                w.WriteData(0.5 * k, xts_step_mesh(mt::tri_mesh(), k));
                w.Flush();
                // The file must be readable after EVERY step, not only at the end.
                const meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
                EXPECT_EQ(meta.mTimeValues.size(), static_cast<std::size_t>(k + 1))
                    << "format " << p_format << " after step " << k;
            }
            w.Finalize();
        }
        xts_cleanup(path);
    }
}

TEST(XdmfTimeSeries, AutoFlushIsOffByDefaultAndCanBeEnabled) {
    const std::string path = mt::temp_path(".xdmf");
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        EXPECT_FALSE(w.AutoFlush());
        w.WritePointsCells(mt::tri_mesh());
        w.WriteData(0.0, xts_step_mesh(mt::tri_mesh(), 0));
        // Default: nothing on disk yet -- the historical behaviour.
        EXPECT_FALSE(std::filesystem::exists(path));
        w.SetAutoFlush(true);
        EXPECT_TRUE(w.AutoFlush());
        w.WriteData(1.0, xts_step_mesh(mt::tri_mesh(), 1));
        EXPECT_TRUE(std::filesystem::exists(path));
        w.Finalize();
    }
    xts_cleanup(path);
}

#ifndef _WIN32
TEST(XdmfTimeSeries, KilledRunLeavesAReadableSeriesAndAppendContinuesIt) {
    // The headline case: a solve that dies mid-series must still leave an
    // openable .xdmf, and a restart must continue that same collection rather
    // than overwrite it.
    for (const char* p_format : {"XML", "Binary", "HDF"}) {
        const std::string path = mt::temp_path(".xdmf");
        const pid_t pid = fork();
        ASSERT_NE(pid, -1) << "fork failed";
        if (pid == 0) {
            // Child: write three steps, flushing each, then die hard. _exit and
            // SIGKILL both skip every destructor, so nothing is finalized.
            meshioplusplus::XdmfTimeSeriesWriter* p_w =
                new meshioplusplus::XdmfTimeSeriesWriter(path, p_format);
            p_w->WritePointsCells(mt::tri_mesh());
            for (int k = 0; k < 3; ++k) {
                p_w->WriteData(0.5 * k, xts_step_mesh(mt::tri_mesh(), k));
                p_w->Flush();
            }
            std::_Exit(9);  // no unwinding, no Finalize()
        }
        int status = 0;
        ASSERT_EQ(waitpid(pid, &status, 0), pid);
        ASSERT_TRUE(WIFEXITED(status) && WEXITSTATUS(status) == 9)
            << "child did not die as expected (format " << p_format << ")";

        // The killed run's output is readable and complete up to the last flush.
        const meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
        ASSERT_EQ(meta.mTimeValues.size(), 3u) << "format " << p_format;
        EXPECT_DOUBLE_EQ(meta.mTimeValues[2], 1.0);
        const Mesh recovered = meshioplusplus::read_xdmf(path);
        EXPECT_EQ(recovered.NumPoints(), mt::tri_mesh().NumPoints());

        // Restart: append two more steps to the same file.
        {
            meshioplusplus::XdmfTimeSeriesWriter w(path, p_format, -1,
                                                   meshioplusplus::XdmfSeriesMode::Append);
            EXPECT_EQ(w.NumSteps(), 3u);
            for (int k = 3; k < 5; ++k)
                w.WriteData(0.5 * k, xts_step_mesh(mt::tri_mesh(), k));
            w.Finalize();
            EXPECT_EQ(w.NumSteps(), 5u);
        }
        const meshioplusplus::MeshMetadata after = meshioplusplus::read_xdmf_metadata(path);
        ASSERT_EQ(after.mTimeValues.size(), 5u) << "format " << p_format;
        EXPECT_DOUBLE_EQ(after.mTimeValues[4], 2.0);
        // Nothing the first run wrote was overwritten: every step still reads,
        // which is what a mis-resumed heavy-data counter would break.
        for (int k = 0; k < 5; ++k) {
            meshioplusplus::ReadOptions opts;
            opts.mTimeStep = k;
            const Mesh step = meshioplusplus::read_xdmf(path, opts);
            ASSERT_TRUE(step.HasPointData("T")) << "step " << k;
            EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(step.PointData("T"), 0), 100.0 * k)
                << "format " << p_format << " step " << k;
        }
        xts_cleanup(path);
    }
}
#endif  // _WIN32

TEST(XdmfTimeSeries, AppendToAMissingPathIsJustAFreshSeries) {
    // So a restartable solver can pass Append unconditionally.
    const std::string path = mt::temp_path(".xdmf");
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML", -1,
                                               meshioplusplus::XdmfSeriesMode::Append);
        EXPECT_EQ(w.NumSteps(), 0u);
        w.WritePointsCells(mt::tri_mesh());
        w.WriteData(0.0, xts_step_mesh(mt::tri_mesh(), 0));
        w.Finalize();
    }
    EXPECT_EQ(meshioplusplus::read_xdmf_metadata(path).mTimeValues.size(), 1u);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, NamedArrayOverloadWritesTheSameShapeAsTheMeshOne) {
    const std::string path = mt::temp_path(".xdmf");
    const Mesh base = mt::tri_mesh();
    const std::size_t np = base.NumPoints();
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        for (int k = 0; k < 2; ++k) {
            meshioplusplus::XdmfTimeSeriesWriter::NamedArray t;
            t.mName = "T";
            t.mNumComponents = 1;
            t.mValues.assign(np, static_cast<double>(k));
            w.WriteData(0.5 * k, {t});
        }
        w.Finalize();
    }
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    EXPECT_EQ(meta.mTimeValues.size(), 2u);
    meshioplusplus::ReadOptions opts;
    opts.mTimeStep = 1;
    const Mesh back = meshioplusplus::read_xdmf(path, opts);
    ASSERT_TRUE(back.HasPointData("T"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.PointData("T"), 0), 1.0);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, NamedArrayOverloadRejectsAMisSizedArray) {
    const std::string path = mt::temp_path(".xdmf");
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(mt::tri_mesh());
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray bad;
        bad.mName = "u";
        bad.mNumComponents = 3;
        bad.mValues.assign(7, 0.0);  // not NumPoints * 3
        EXPECT_THROW(w.WriteData(0.0, {bad}), meshioplusplus::WriteError);
    }
    xts_cleanup(path);
}

// --- v9.2.0: Append + the NamedArray overload ------------------------------
//
// Through v9.1.0 these two features -- added in the same release, for the same
// restartable-solver consumer -- did not work together. Impl::mNumPoints /
// mNumCells were set only by WritePointsCells, which an appending writer cannot
// call (mHasMesh is already set, so it throws "called more than once"), so the
// NamedArray overload validated every array against 0 and rejected all of them
// with "expected 0 (0 x 1)". Nothing combined the two, which is why it shipped.

TEST(XdmfTimeSeries, AppendThenNamedArrayWriteData) {
    const std::string path = mt::temp_path(".xdmf");
    const Mesh base = mt::tri_mesh();
    const std::size_t np = base.NumPoints();
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray u;
        u.mName = "u";
        u.mNumComponents = 1;
        u.mValues.assign(np, 1.5);
        w.WriteData(0.0, {u});
        w.Finalize();
    }
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML", -1,
                                               meshioplusplus::XdmfSeriesMode::Append);
        ASSERT_EQ(w.NumSteps(), 1u);
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray u;
        u.mName = "u";
        u.mNumComponents = 1;
        u.mValues.assign(np, 2.5);
        // The regression: this threw before the counts were recovered.
        w.WriteData(1.0, {u});
        w.Finalize();
        EXPECT_EQ(w.NumSteps(), 2u);
    }
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_xdmf_metadata(path);
    ASSERT_EQ(meta.mTimeValues.size(), 2u);
    meshioplusplus::ReadOptions opts;
    opts.mTimeStep = 1;
    const Mesh back = meshioplusplus::read_xdmf(path, opts);
    ASSERT_TRUE(back.HasPointData("u"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.PointData("u"), 0), 2.5);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, AppendRejectsAMisSizedNamedArrayJustLikeAFreshSeries) {
    // The recovered counts must actually be USED, not merely be non-zero:
    // deriving rows from the array would make this write happily.
    const std::string path = mt::temp_path(".xdmf");
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(mt::tri_mesh());
        w.WriteData(0.0, xts_step_mesh(mt::tri_mesh(), 0));
        w.Finalize();
    }
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML", -1,
                                               meshioplusplus::XdmfSeriesMode::Append);
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray bad;
        bad.mName = "u";
        bad.mNumComponents = 3;
        bad.mValues.assign(7, 0.0);  // not NumPoints * 3
        EXPECT_THROW(w.WriteData(1.0, {bad}), meshioplusplus::WriteError);
    }
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, NamedArrayOverloadWritesCellArrays) {
    // Cell arrays were untested in either mode; only the point half had cover.
    const std::string path = mt::temp_path(".xdmf");
    const Mesh base = mt::tri_mesh();
    std::size_t ncells = 0;
    for (std::size_t b = 0; b < base.NumCellBlocks(); ++b)
        ncells += base.Cells(b).NumCells();
    ASSERT_GT(ncells, 0u);
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray rho;
        rho.mName = "rho";
        rho.mNumComponents = 1;
        rho.mValues.assign(ncells, 7.0);
        w.WriteData(0.0, {}, {rho});
        w.Finalize();
    }
    const Mesh back = meshioplusplus::read_xdmf(path);
    ASSERT_TRUE(back.HasCellData("rho"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.CellData("rho", 0), 0), 7.0);
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, AppendToAMixedTopologySeriesRecoversTheCellCount) {
    // A multi-block mesh writes TopologyType="Mixed", whose DataItem Dimensions
    // is the flat packed length, NOT the cell count -- so the recovery has to
    // read NumberOfElements. read_xdmf_metadata declines Mixed outright, which
    // is exactly why xdmf_grid_counts exists rather than reusing it.
    const std::string path = mt::temp_path(".xdmf");
    const Mesh base = mt::tri_quad_mesh();
    std::size_t ncells = 0;
    for (std::size_t b = 0; b < base.NumCellBlocks(); ++b)
        ncells += base.Cells(b).NumCells();
    ASSERT_GT(base.NumCellBlocks(), 1u);
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(base);
        w.Finalize();
    }
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML", -1,
                                               meshioplusplus::XdmfSeriesMode::Append);
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray rho;
        rho.mName = "rho";
        rho.mNumComponents = 1;
        rho.mValues.assign(ncells, 3.0);
        w.WriteData(0.0, {}, {rho});  // correct length must pass...
        meshioplusplus::XdmfTimeSeriesWriter::NamedArray bad = rho;
        bad.mValues.assign(ncells + 1, 3.0);
        EXPECT_THROW(w.WriteData(1.0, {}, {bad}), meshioplusplus::WriteError);  // ...wrong must not
        w.Finalize();
    }
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, AppendUsesTheSameStructuralResolutionAsTheReader) {
    // The static grid is found by GridType="Uniform", not by Name=="mesh". The
    // pre-9.2.0 transcription matched the literal name, so appending to another
    // producer's series set mHasMesh=false and added a SECOND static grid.
    const std::string path = mt::temp_path(".xdmf");
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML");
        w.WritePointsCells(mt::tri_mesh());
        w.WriteData(0.0, xts_step_mesh(mt::tri_mesh(), 0));
        w.Finalize();
    }
    // Rename the static grid, as a foreign producer would. Done as text: this
    // test TU cannot include pugixml (a build-only vendored PRIVATE dependency).
    {
        std::ifstream in(path);
        ASSERT_TRUE(in.good());
        std::string xml((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
        const std::string from = "Name=\"mesh\" GridType=\"Uniform\"";
        const std::string to = "Name=\"not_called_mesh\" GridType=\"Uniform\"";
        const std::size_t at = xml.find(from);
        ASSERT_NE(at, std::string::npos);
        // The step xpointers still say mesh, which is the point: this reader
        // never follows them and resolves structurally instead.
        xml.replace(at, from.size(), to);
        std::ofstream out(path);
        out << xml;
    }
    {
        meshioplusplus::XdmfTimeSeriesWriter w(path, "XML", -1,
                                               meshioplusplus::XdmfSeriesMode::Append);
        EXPECT_EQ(w.NumSteps(), 1u);
        // The mesh grid was found, so a second one must be refused.
        EXPECT_THROW(w.WritePointsCells(mt::tri_mesh()), meshioplusplus::WriteError);
    }
    xts_cleanup(path);
}

TEST(XdmfTimeSeries, MovedFromWriterIsDiagnosableNotUndefined) {
    const std::string path = mt::temp_path(".xdmf");
    meshioplusplus::XdmfTimeSeriesWriter src(path, "XML");
    src.WritePointsCells(mt::tri_mesh());
    meshioplusplus::XdmfTimeSeriesWriter dst(std::move(src));

    // Observers and idempotent operations answer as for a finished series.
    EXPECT_EQ(src.NumSteps(), 0u);            // NOLINT(bugprone-use-after-move)
    EXPECT_TRUE(src.Finalized());             // NOLINT(bugprone-use-after-move)
    EXPECT_FALSE(src.AutoFlush());            // NOLINT(bugprone-use-after-move)
    EXPECT_NO_THROW(src.SetAutoFlush(true));  // NOLINT(bugprone-use-after-move)
    EXPECT_NO_THROW(src.Flush());             // NOLINT(bugprone-use-after-move)
    EXPECT_NO_THROW(src.Finalize());          // NOLINT(bugprone-use-after-move)

    // The three that cannot silently do nothing must say so.
    meshioplusplus::XdmfTimeSeriesWriter::NamedArray u;
    u.mName = "u";
    u.mValues.assign(mt::tri_mesh().NumPoints(), 0.0);
    EXPECT_THROW(src.WritePointsCells(mt::tri_mesh()), meshioplusplus::WriteError);
    EXPECT_THROW(src.WriteData(0.0, mt::tri_mesh()), meshioplusplus::WriteError);
    EXPECT_THROW(src.WriteData(0.0, {u}), meshioplusplus::WriteError);

    dst.Finalize();
    xts_cleanup(path);
}
