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
// ParaView collection `.pvd` (roadmap §1.1, v14.1.0): a time-indexed list of VTK
// XML files. `timestep` selects the step, `part` the piece within it.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pvd.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"

using meshioplusplus::DType;
using meshioplusplus::GhostPolicy;
using meshioplusplus::kSequenceTimeKey;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::PvdSeriesWriter;
using meshioplusplus::PvtuReadOptions;
using meshioplusplus::read_pvd;
using meshioplusplus::read_pvd_metadata;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::SequenceInput;
using meshioplusplus::SequenceOutput;
using meshioplusplus::write_pvd;
using meshioplusplus::WriteError;
using meshioplusplus::detail::VtkCodec;

namespace {

namespace fs = std::filesystem;

// A directory for a collection and its pieces; removed on destruction.
struct PvdDir {
    fs::path mDir;
    PvdDir() {
        mDir = fs::path(mt::temp_path(".d"));
        fs::create_directories(mDir);
    }
    ~PvdDir() {
        std::error_code ec;
        fs::remove_all(mDir, ec);
    }
    std::string operator/(const std::string& rName) const { return (mDir / rName).string(); }
};

std::string pvd_slurp(const std::string& rPath) {
    std::ifstream in(rPath);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::size_t pvd_count(const std::string& rText, const std::string& rNeedle) {
    std::size_t n = 0;
    for (std::size_t at = rText.find(rNeedle); at != std::string::npos;
         at = rText.find(rNeedle, at + rNeedle.size()))
        ++n;
    return n;
}

std::size_t pvd_num_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += cb.NumCells();
    return n;
}

// One triangle shifted by K along x, carrying u = K.
Mesh pvd_step(double K) {
    NDArray pts(DType::Float64, {3, 3});
    double* p = pts.As<double>();
    p[0] = K;
    p[3] = 1.0 + K;
    p[6] = K;
    p[7] = 1.0;
    NDArray conn(DType::Int64, {1, 3});
    conn.As<std::int64_t>()[1] = 1;
    conn.As<std::int64_t>()[2] = 2;
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("triangle", std::move(conn));
    NDArray u(DType::Float64, {3});
    for (std::size_t i = 0; i < 3; ++i)
        u.As<double>()[i] = K;
    m.AddPointData("u", std::move(u));
    return m;
}

double pvd_x0(const Mesh& rMesh) {
    return rMesh.Points().As<double>()[0];
}

double pvd_time(const Mesh& rMesh) {
    return *rMesh.FieldData(kSequenceTimeKey).As<double>();
}

void pvd_write_index(const std::string& rPath, const std::vector<std::string>& rEntries) {
    std::ofstream out(rPath);
    out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n<Collection>\n";
    for (const std::string& e : rEntries)
        out << "<DataSet " << e << "/>\n";
    out << "</Collection>\n</VTKFile>\n";
}

// Three one-triangle steps at t = 0, 0.5, 2 in `<dir>/t.pvd`.
std::string pvd_series(const PvdDir& rDir) {
    const std::string path = rDir / "t.pvd";
    PvdSeriesWriter writer(path, /*binary=*/false, VtkCodec::None);
    writer.Write(0.0, pvd_step(0));
    writer.Write(0.5, pvd_step(1));
    writer.Write(2.0, pvd_step(2));
    writer.Finalize();
    return path;
}

}  // namespace

TEST(Pvd, APlainWriteIsAOneStepCollectionTimedByMeshioTime) {
    PvdDir d;
    Mesh m = pvd_step(1);
    NDArray t(DType::Float64, {1});
    *t.As<double>() = 0.25;
    m.AddFieldData(kSequenceTimeKey, std::move(t));
    write_pvd(d / "c.pvd", m, /*binary=*/true, /*zlib=*/false);

    const std::string text = pvd_slurp(d / "c.pvd");
    EXPECT_NE(text.find("<VTKFile type=\"Collection\""), std::string::npos);
    EXPECT_EQ(pvd_count(text, "<DataSet "), 1u);
    EXPECT_NE(text.find("timestep=\"0.25\""), std::string::npos) << text;
    EXPECT_NE(text.find("file=\"c/c_0000.vtu\""), std::string::npos) << text;
    EXPECT_TRUE(fs::exists(d.mDir / "c" / "c_0000.vtu"));

    const Mesh back = read_pvd(d / "c.pvd");
    EXPECT_EQ(pvd_time(back), 0.25);
    EXPECT_EQ(pvd_x0(back), 1.0);
    EXPECT_TRUE(back.HasPointData("u"));

    // no time recorded: step zero
    write_pvd(d / "z.pvd", pvd_step(0), false, false);
    EXPECT_NE(pvd_slurp(d / "z.pvd").find("timestep=\"0\""), std::string::npos);
}

TEST(Pvd, TheSeriesWriterRewritesTheIndexAfterEveryStep) {
    PvdDir d;
    PvdSeriesWriter writer(d / "k.pvd", false, VtkCodec::None);
    EXPECT_EQ(writer.NumSteps(), 0u);
    writer.Write(0.0, pvd_step(0));
    EXPECT_EQ(pvd_count(pvd_slurp(d / "k.pvd"), "<DataSet "), 1u);
    writer.Write(1.0, pvd_step(1));
    // no Finalize: a run killed here still leaves a collection covering both steps
    EXPECT_EQ(pvd_count(pvd_slurp(d / "k.pvd"), "<DataSet "), 2u);
    EXPECT_EQ(read_pvd_metadata(d / "k.pvd").mTimeValues, (std::vector<double>{0.0, 1.0}));
    EXPECT_EQ(writer.NumSteps(), 2u);
}

TEST(Pvd, TheSeriesWriterRefusesWhatCannotBeRepresented) {
    PvdDir d;
    {
        PvdSeriesWriter writer(d / "n.pvd", false, VtkCodec::None);
        EXPECT_THROW(writer.Write(std::numeric_limits<double>::quiet_NaN(), pvd_step(0)),
                     WriteError);
        EXPECT_THROW(writer.Write(std::numeric_limits<double>::infinity(), pvd_step(0)),
                     WriteError);
        EXPECT_THROW(writer.Finalize(), WriteError);  // no step was written
        writer.Write(1.0, pvd_step(0));
        writer.Finalize();
        writer.Finalize();                                         // idempotent
        EXPECT_THROW(writer.Write(2.0, pvd_step(0)), WriteError);  // after Finalize
    }
    PvdSeriesWriter a(d / "m.pvd", false, VtkCodec::None);
    PvdSeriesWriter b(std::move(a));
    EXPECT_THROW(a.Write(0.0, pvd_step(0)), WriteError);  // moved-from
    EXPECT_EQ(a.NumSteps(), 0u);
    EXPECT_NO_THROW(a.Finalize());
    b.Write(0.0, pvd_step(0));
    EXPECT_EQ(b.NumSteps(), 1u);
}

TEST(Pvd, TimesAreWrittenAsTheShortestExactSpelling) {
    PvdDir d;
    PvdSeriesWriter w(d / "s.pvd", false, VtkCodec::None);
    for (const double t : {0.1, 1.0 / 3.0, 1e-9, 123456.789})
        w.Write(t, pvd_step(0));
    w.Finalize();
    const std::string text = pvd_slurp(d / "s.pvd");
    EXPECT_NE(text.find("timestep=\"0.1\""), std::string::npos) << text;  // not 0.10000000000000001
    const std::vector<double> back = read_pvd_metadata(d / "s.pvd").mTimeValues;
    ASSERT_EQ(back.size(), 4u);
    EXPECT_EQ(back[0], 1e-9);  // sorted ascending, and every value exact
    EXPECT_EQ(back[1], 0.1);
    EXPECT_EQ(back[2], 1.0 / 3.0);
    EXPECT_EQ(back[3], 123456.789);
}

TEST(Pvd, ReadSelectsAStepByTimeStepAndNegativeCountsFromTheEnd) {
    PvdDir d;
    const std::string path = pvd_series(d);
    EXPECT_EQ(pvd_x0(read_pvd(path)), 0.0);
    ReadOptions o;
    o.mTimeStep = 1;
    EXPECT_EQ(pvd_x0(read_pvd(path, o)), 1.0);
    o.mTimeStep = -1;
    const Mesh last = read_pvd(path, o);
    EXPECT_EQ(pvd_x0(last), 2.0);
    EXPECT_EQ(pvd_time(last), 2.0);
    for (const int bad : {3, -4}) {
        o.mTimeStep = bad;
        try {
            read_pvd(path, o);
            FAIL() << "expected a ReadError for time_step " << bad;
        } catch (const ReadError& e) {
            const std::string msg = e.what();
            EXPECT_NE(msg.find("is out of range"), std::string::npos) << msg;
            EXPECT_NE(msg.find("3 steps"), std::string::npos) << msg;
        }
    }
}

TEST(Pvd, StepsAreTheDistinctTimesAscendingWhateverTheirOrderInTheFile) {
    PvdDir d;
    for (int i = 0; i < 3; ++i)
        meshioplusplus::write_vtu(d / ("p" + std::to_string(i) + ".vtu"), pvd_step(i), false,
                                  false);
    pvd_write_index(d / "u.pvd",
                    {"timestep=\"2.5\" file=\"p0.vtu\"", "timestep=\"0.1\" file=\"p1.vtu\"",
                     "timestep=\"1\" file=\"p2.vtu\""});
    ReadOptions o;
    o.mTimeStep = 0;
    EXPECT_EQ(pvd_x0(read_pvd(d / "u.pvd", o)), 1.0);  // t = 0.1
    o.mTimeStep = 1;
    EXPECT_EQ(pvd_x0(read_pvd(d / "u.pvd", o)), 2.0);  // t = 1
    o.mTimeStep = 2;
    EXPECT_EQ(pvd_x0(read_pvd(d / "u.pvd", o)), 0.0);  // t = 2.5
    EXPECT_EQ(read_pvd_metadata(d / "u.pvd").mTimeValues, (std::vector<double>{0.1, 1.0, 2.5}));
}

TEST(Pvd, AMissingTimestepAndPartAreZero) {
    PvdDir d;
    meshioplusplus::write_vtu(d / "a.vtu", pvd_step(0), false, false);
    pvd_write_index(d / "n.pvd", {"file=\"a.vtu\""});
    const Mesh back = read_pvd(d / "n.pvd");
    EXPECT_EQ(pvd_time(back), 0.0);
    EXPECT_EQ(pvd_num_cells(back), 1u);
}

namespace {

// Two steps x three parts, each piece a triangle shifted by 10*step + part.
std::string pvd_two_by_three(const PvdDir& rDir) {
    std::vector<std::string> entries;
    for (int t = 0; t < 2; ++t)
        for (int p = 0; p < 3; ++p) {
            const std::string name = "s" + std::to_string(t) + "p" + std::to_string(p) + ".vtu";
            meshioplusplus::write_vtu(rDir / name, pvd_step(10 * t + p), false, false);
            entries.push_back("timestep=\"" + std::to_string(t) + "\" part=\"" + std::to_string(p) +
                              "\" file=\"" + name + "\"");
        }
    pvd_write_index(rDir / "tp.pvd", entries);
    return rDir / "tp.pvd";
}

}  // namespace

TEST(Pvd, PartAndTimestepAreOrthogonal) {
    PvdDir d;
    const std::string path = pvd_two_by_three(d);
    const Mesh step0 = read_pvd(path);
    EXPECT_EQ(pvd_num_cells(step0), 3u);
    EXPECT_EQ(step0.NumPoints(), 9u);
    ASSERT_EQ(step0.NumRegions(), 3u);
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_EQ(step0.Region(i).mName, "part_" + std::to_string(i));

    ReadOptions o;
    o.mTimeStep = 1;
    EXPECT_EQ(pvd_x0(read_pvd(path, o)), 10.0);

    o.mPiece = 2;
    o.mPieceSet = true;
    const Mesh one = read_pvd(path, o);
    EXPECT_EQ(pvd_x0(one), 12.0);
    EXPECT_EQ(pvd_num_cells(one), 1u);
    EXPECT_EQ(one.NumRegions(), 0u);
    o.mPiece = -1;
    EXPECT_EQ(pvd_x0(read_pvd(path, o)), 12.0);
    o.mPiece = 3;
    try {
        read_pvd(path, o);
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("3 pieces"), std::string::npos) << e.what();
    }
}

TEST(Pvd, PartsAreOrderedByPartThenDocumentOrderAndNamedFromNameThenGroup) {
    PvdDir d;
    for (int i = 0; i < 3; ++i)
        meshioplusplus::write_vtu(d / ("p" + std::to_string(i) + ".vtu"), pvd_step(i), false,
                                  false);
    pvd_write_index(d / "g.pvd", {"timestep=\"0\" part=\"2\" file=\"p2.vtu\"",
                                  "timestep=\"0\" part=\"1\" group=\"fluid\" file=\"p1.vtu\"",
                                  "timestep=\"0\" part=\"0\" name=\"wing\" file=\"p0.vtu\""});
    const Mesh back = read_pvd(d / "g.pvd");
    ASSERT_EQ(back.NumRegions(), 3u);
    // regions are stored in name order; the entries were merged in part order
    EXPECT_TRUE(back.HasRegion("wing"));
    EXPECT_TRUE(back.HasRegion("fluid/part_1"));
    EXPECT_TRUE(back.HasRegion("part_2"));
    EXPECT_EQ(pvd_x0(back), 0.0);  // part 0 first
    ReadOptions o;
    o.mPiece = 0;
    o.mPieceSet = true;
    EXPECT_EQ(pvd_x0(read_pvd(d / "g.pvd", o)), 0.0);
}

TEST(Pvd, MetadataReportsEveryTimeWithoutOpeningALaterStep) {
    PvdDir d;
    const std::string path = pvd_series(d);
    fs::remove(d.mDir / "t" / "t_0001.vtu");
    fs::remove(d.mDir / "t" / "t_0002.vtu");
    const meshioplusplus::MeshMetadata meta = read_pvd_metadata(path);
    EXPECT_EQ(meta.mFormat, "pvd");
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.0, 0.5, 2.0}));
    EXPECT_EQ(meta.mNumPoints, 3u);  // step 0's pieces
    EXPECT_EQ(meshioplusplus::sequence_num_steps(path, ""), 3u);
}

namespace {

// `.pvd` -> `.pvtu` -> `.vtu`: a partitioned transient run, one .pvtu per step.
std::string pvd_partitioned_run(const PvdDir& rDir, int NSteps, int Layers, std::size_t& rCells) {
    const int N = 5;
    const std::size_t np = static_cast<std::size_t>((N + 1) * (N + 1));
    std::vector<std::string> entries;
    for (int t = 0; t < NSteps; ++t) {
        NDArray pts(DType::Float64, {np, 3});
        NDArray u(DType::Float64, {np});
        for (int j = 0; j <= N; ++j)
            for (int i = 0; i <= N; ++i) {
                const std::size_t k = static_cast<std::size_t>(j * (N + 1) + i);
                double* p = pts.As<double>() + 3 * k;
                p[0] = i;
                p[1] = j;
                p[2] = t;
                u.As<double>()[k] = t;
            }
        NDArray conn(DType::Int64, {static_cast<std::size_t>(2 * N * N), 3});
        std::size_t at = 0;
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const std::int64_t a = j * (N + 1) + i, b = a + 1, dd = a + N + 1, e = dd + 1;
                const std::int64_t tris[6] = {a, b, e, a, e, dd};
                std::copy(tris, tris + 6, conn.As<std::int64_t>() + 3 * at);
                at += 2;
            }
        Mesh m;
        m.AssignPoints(std::move(pts));
        m.AddCellBlock("triangle", std::move(conn));
        m.AddPointData("u", std::move(u));
        rCells = static_cast<std::size_t>(2 * N * N);

        meshioplusplus::PartitionOptions po;
        po.mNParts = 3;
        po.mMethod = meshioplusplus::PartitionMethod::SFC;
        po.mGhostLayers = Layers;
        const meshioplusplus::PartitionResult result = meshioplusplus::partition(m, po);
        std::vector<const Mesh*> ptrs;
        for (const auto& piece : result.mPieces)
            ptrs.push_back(&piece.mMesh);
        const std::string name = "step" + std::to_string(t) + ".pvtu";
        meshioplusplus::write_pvtu_pieces_codec(rDir / name, ptrs, false, VtkCodec::None);
        entries.push_back("timestep=\"" + std::to_string(0.5 * t) + "\" part=\"0\" file=\"" + name +
                          "\"");
    }
    pvd_write_index(rDir / "run.pvd", entries);
    return rDir / "run.pvd";
}

}  // namespace

TEST(Pvd, APvdOfPvtuOfVtuNestsWithoutSpecialCases) {
    PvdDir d;
    std::size_t cells = 0;
    const std::string path = pvd_partitioned_run(d, 3, 0, cells);
    for (int t = 0; t < 3; ++t) {
        ReadOptions o;
        o.mTimeStep = t;
        const Mesh back = read_pvd(path, o);
        // a step that is one .pvtu reads as that file does: its own piece regions
        ASSERT_EQ(back.NumRegions(), 3u);
        EXPECT_EQ(back.Region(0).mName, "piece_0");
        EXPECT_EQ(pvd_num_cells(back), cells);
        EXPECT_EQ(pvd_time(back), 0.5 * t);
    }
    EXPECT_EQ(read_pvd_metadata(path).mTimeValues, (std::vector<double>{0.0, 0.5, 1.0}));
}

TEST(Pvd, TheGhostPolicyReachesEveryChild) {
    PvdDir d;
    std::size_t cells = 0;
    const std::string path = pvd_partitioned_run(d, 2, 1, cells);
    const Mesh kept = read_pvd(path);
    EXPECT_GT(pvd_num_cells(kept), cells);
    EXPECT_TRUE(kept.HasCellData("vtkGhostType"));
    PvtuReadOptions drop;
    drop.mGhosts = GhostPolicy::Drop;
    const Mesh dropped = read_pvd(path, {}, drop);
    EXPECT_EQ(pvd_num_cells(dropped), cells);
    EXPECT_FALSE(dropped.HasCellData("vtkGhostType"));

    // ... and a plain .vtu child carrying a ghost flag
    Mesh piece = pvd_step(0);
    NDArray flag(DType::UInt8, {1});
    flag.As<std::uint8_t>()[0] = 1;
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(flag));
    piece.AddCellData("vtkGhostType", std::move(blocks));
    meshioplusplus::write_vtu(d / "g.vtu", piece, false, false);
    pvd_write_index(d / "g.pvd", {"timestep=\"0\" file=\"g.vtu\""});
    EXPECT_EQ(pvd_num_cells(read_pvd(d / "g.pvd")), 1u);
    EXPECT_EQ(pvd_num_cells(read_pvd(d / "g.pvd", {}, drop)), 0u);
}

TEST(Pvd, RefusesWhatItCannotRead) {
    PvdDir d;
    {  // an entry that is not one of the XML family: never legacy .vtk
        std::ofstream(d.mDir / "old.vtk") << "# vtk DataFile Version 3.0\n";
        pvd_write_index(d / "l.pvd", {"timestep=\"0\" file=\"old.vtk\""});
        try {
            read_pvd(d / "l.pvd");
            FAIL() << "expected a ReadError";
        } catch (const ReadError& e) {
            const std::string msg = e.what();
            EXPECT_NE(msg.find("unsupported piece"), std::string::npos) << msg;
            EXPECT_NE(msg.find(".vtu/.vtp/.vtm/.pvtu/.pvtp"), std::string::npos) << msg;
        }
    }
    pvd_write_index(d / "m.pvd", {"timestep=\"0\" file=\"gone/x.vtu\""});
    try {
        read_pvd(d / "m.pvd");
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("(file=)"), std::string::npos) << e.what();
    }
    pvd_write_index(d / "f.pvd", {"timestep=\"0\""});
    EXPECT_THROW(read_pvd(d / "f.pvd"), ReadError);
    pvd_write_index(d / "b.pvd", {"timestep=\"soon\" file=\"a.vtu\""});
    EXPECT_THROW(read_pvd(d / "b.pvd"), ReadError);
    std::ofstream(d.mDir / "w.pvd")
        << "<VTKFile type=\"vtkMultiBlockDataSet\"><Collection/></VTKFile>";
    EXPECT_THROW(read_pvd(d / "w.pvd"), ReadError);
    std::ofstream(d.mDir / "x.pvd") << "not xml";
    EXPECT_THROW(read_pvd(d / "x.pvd"), ReadError);

    pvd_write_index(d / "e.pvd", {});
    const Mesh empty = read_pvd(d / "e.pvd");
    EXPECT_EQ(empty.NumPoints(), 0u);
    EXPECT_EQ(empty.NumCellBlocks(), 0u);
}

// --- through the sequence engine ------------------------------------------------

TEST(PvdSequence, TheWriteCapabilityNamesPvd) {
    std::string why;
    EXPECT_TRUE(meshioplusplus::sequence_write_supports_time("pvd", why));
    EXPECT_TRUE(why.empty());
    EXPECT_FALSE(meshioplusplus::sequence_write_supports_time("vtu", why));
    EXPECT_NE(why.find("'pvd'"), std::string::npos) << why;
    EXPECT_TRUE(meshioplusplus::seq_format_may_have_steps("pvd"));
    // pvtu/pvtp are parts, not steps
    EXPECT_FALSE(meshioplusplus::seq_format_may_have_steps("pvtu"));
    EXPECT_FALSE(meshioplusplus::sequence_write_supports_time("pvtu", why));
}

TEST(PvdSequence, AFanInWritesACollectionWithTheGivenTimes) {
    PvdDir d;
    SequenceInput in;
    for (int i = 0; i < 3; ++i) {
        const std::string p = d / ("in_" + std::to_string(i) + ".vtu");
        meshioplusplus::write_vtu(p, pvd_step(i), false, false);
        in.mPaths.push_back(p);
    }
    in.mTimes = {0.0, 0.5, 2.0};
    SequenceOutput out;
    out.mPath = d / "series.pvd";
    meshioplusplus::sequence_to_timeseries(in, out);

    EXPECT_EQ(pvd_count(pvd_slurp(out.mPath), "<DataSet "), 3u);
    EXPECT_EQ(read_pvd_metadata(out.mPath).mTimeValues, (std::vector<double>{0.0, 0.5, 2.0}));
    ReadOptions o;
    o.mTimeStep = 2;
    EXPECT_EQ(pvd_x0(read_pvd(out.mPath, o)), 2.0);
    EXPECT_TRUE(fs::exists(d.mDir / "series" / "series_0002.vtu"));
}

TEST(PvdSequence, AFanOutWritesOneFilePerStep) {
    PvdDir d;
    const std::string path = pvd_series(d);
    SequenceOutput out;
    out.mPath = d / "out_{step}.vtu";
    meshioplusplus::timeseries_to_sequence(path, "", ReadOptions{}, out);
    for (int i = 0; i < 3; ++i) {
        const std::string p = d / ("out_000" + std::to_string(i) + ".vtu");
        ASSERT_TRUE(fs::exists(p)) << p;
        EXPECT_EQ(pvd_x0(meshioplusplus::read_vtu(p)), static_cast<double>(i));
    }
}

TEST(PvdSequence, TheTransientWriterHonoursEncodingAndRejectsWhatItCannotHonour) {
    PvdDir d;
    SequenceInput in;
    for (int i = 0; i < 2; ++i) {
        const std::string p = d / ("in_" + std::to_string(i) + ".vtu");
        meshioplusplus::write_vtu(p, pvd_step(i), false, false);
        in.mPaths.push_back(p);
    }
    SequenceOutput out;
    out.mPath = d / "e.pvd";
    out.mOptions.mEncoding = meshioplusplus::WriteEncoding::Ascii;
    EXPECT_NO_THROW(meshioplusplus::sequence_to_timeseries(in, out));
    const std::string piece = pvd_slurp((d.mDir / "e" / "e_0000.vtu").string());
    EXPECT_GT(pvd_count(piece, "format=\"ascii\""), 0u);
    EXPECT_EQ(pvd_count(piece, "format=\"binary\""), 0u);

    SequenceOutput bad;
    bad.mPath = d / "c.pvd";
    bad.mOptions.mCodec = VtkCodec::Zlib;
    bad.mOptions.mCodecSet = true;
    try {
        meshioplusplus::sequence_to_timeseries(in, bad);
        FAIL() << "expected a WriteError";
    } catch (const WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("PVD writer does not support Codec"),
                  std::string::npos)
            << e.what();
    }
    SequenceOutput fmt;
    fmt.mPath = d / "f.pvd";
    fmt.mOptions.mFloatFormat = ".6e";
    EXPECT_THROW(meshioplusplus::sequence_to_timeseries(in, fmt), WriteError);
}

TEST(PvdRegistry, IsReadableWritableAndDispatchedByExtension) {
    PvdDir d;
    meshioplusplus::registry_write_ex(d / "r.pvd", pvd_step(1), "", meshioplusplus::WriteOptions{});
    EXPECT_EQ(meshioplusplus::registry_extension_defaults().at(".pvd"), "pvd");
    EXPECT_EQ(meshioplusplus::registry_extension_defaults().at(".pvtu"), "pvtu");
    EXPECT_EQ(meshioplusplus::registry_extension_defaults().at(".pvtp"), "pvtp");
    EXPECT_EQ(meshioplusplus::resolve_format(d / "r.pvd", ""), "pvd");
    EXPECT_EQ(meshioplusplus::resolve_format(d / "r.pvtu", ""), "pvtu");
    EXPECT_EQ(meshioplusplus::resolve_format(d / "r.pvtp", ""), "pvtp");
    EXPECT_EQ(pvd_x0(meshioplusplus::registry_read(d / "r.pvd", "pvd", ReadOptions{})), 1.0);
    EXPECT_EQ(meshioplusplus::sniff_format(d / "r.pvd"), "pvd");
    EXPECT_EQ(meshioplusplus::registry_metadata_readers().count("pvd"), 1u);
}
