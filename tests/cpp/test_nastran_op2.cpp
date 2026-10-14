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
/**
 * @file test_nastran_op2.cpp
 * @brief Nastran OP2 reader: a one-quad model (built here block by block) with
 *        GEOM1/GEOM2/EPT tables, a displacement table and a plate stress table,
 *        in 32-bit words of both byte orders and in 64-bit words; the sibling
 *        deck route, steps, metadata, sniffing and refusals. And the writer: the
 *        bytes of the same one-quad model framed by hand, every card round
 *        tripped through the reader, series, and the writer's refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <bit>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::WriteError;

// An OP2 under construction: Fortran blocks (4-byte markers) of Ws-byte words.
struct Op2 {
    int mWs;
    bool mBig;
    std::string mBytes;

    static void Raw(std::string& rOut, std::uint64_t Bits, int N, bool Big) {
        for (int k = 0; k < N; ++k) {
            const int shift = Big ? 8 * (N - 1 - k) : 8 * k;
            rOut.push_back(static_cast<char>((Bits >> shift) & 0xff));
        }
    }
    void Block(const std::string& rPayload) {
        Raw(mBytes, static_cast<std::uint32_t>(rPayload.size()), 4, mBig);
        mBytes += rPayload;
        Raw(mBytes, static_cast<std::uint32_t>(rPayload.size()), 4, mBig);
    }
    std::string Int(std::int64_t V) const {
        std::string s;
        Raw(s,
            mWs == 4 ? static_cast<std::uint32_t>(static_cast<std::int32_t>(V))
                     : static_cast<std::uint64_t>(V),
            mWs, mBig);
        return s;
    }
    std::string Real(double V) const {
        std::string s;
        if (mWs == 4)
            Raw(s, std::bit_cast<std::uint32_t>(static_cast<float>(V)), 4, mBig);
        else
            Raw(s, std::bit_cast<std::uint64_t>(V), 8, mBig);
        return s;
    }
    std::string Text(const std::string& rS) const {  // 8 characters per 2 words
        std::string s = rS;
        s.resize(static_cast<std::size_t>(2 * mWs), ' ');
        return s;
    }
    void Marker(std::int64_t V) { Block(Int(V)); }
    void Record(const std::string& rPayload) {
        Marker(static_cast<std::int64_t>(rPayload.size() / static_cast<std::size_t>(mWs)));
        Block(rPayload);
    }
    // name, -1, 7-word header, -2 1 0, subtable name, then the records with
    // their marker triples, and the closing 0.
    void Table(const std::string& rName, const std::vector<std::string>& rRecords) {
        Record(Text(rName));
        Marker(-1);
        std::string h;
        for (int i = 0; i < 7; ++i)
            h += Int(i == 0 ? 101 : 0);
        Record(h);
        std::int64_t k = -2;
        Marker(k);
        Marker(1);
        Marker(0);
        Record(Text(rName));
        for (const std::string& r : rRecords) {
            --k;
            Marker(k);
            Marker(1);
            Marker(0);
            Record(r);
        }
        --k;
        Marker(k);
        Marker(1);
        Marker(0);
        Marker(0);
    }
    std::string Key(int A, int B, int C) const { return Int(A) + Int(B) + Int(C); }
    // A 146-word result header.
    std::string Header(int Approach, int TCode, int EType, int NumWide, int SCode, int W5 = 1,
                       int Format = 1) const {
        std::string h;
        for (int i = 0; i < 146; ++i) {
            std::int64_t v = 0;
            if (i == 0)
                v = Approach;
            if (i == 1)
                v = TCode;
            if (i == 2)
                v = EType;
            if (i == 3)
                v = 1;  // subcase
            if (i == 4)
                v = W5;  // load set, or a SORT2 table's id * 10 + device
            if (i == 8)
                v = Format;  // 1 real, 2 real/imaginary, 3 magnitude/phase
            if (i == 9)
                v = NumWide;
            if (i == 10)
                v = SCode;
            h += Int(v);
        }
        return h;
    }
};

std::string op2(int Ws, bool Big, bool Geometry) {
    Op2 o{Ws, Big, {}};
    o.Marker(3);
    o.Block(o.Int(9) + o.Int(24) + o.Int(26));  // date
    o.Marker(7);
    o.Block(std::string("NASTRAN FORT TAPE ID CODE - ")
                .substr(0, 28)
                .append(static_cast<std::size_t>(7 * Ws) - 28, ' '));
    o.Record(o.Text("NX2019.2"));
    o.Marker(-1);
    o.Marker(0);
    if (Geometry) {
        const double xy[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
        std::string grids = o.Key(4501, 45, 1);
        for (int n = 0; n < 4; ++n)
            grids += o.Int(n + 1) + o.Int(0) + o.Real(xy[n][0]) + o.Real(xy[n][1]) + o.Real(0.0) +
                     o.Int(0) + o.Int(0) + o.Int(0);
        o.Table("GEOM1S", {grids, o.Key(65535, 65535, 65535)});
        std::string quad = o.Key(2958, 51, 177) + o.Int(10) + o.Int(7);
        for (int n = 1; n <= 4; ++n)
            quad += o.Int(n);
        for (int i = 0; i < 8; ++i)
            quad += o.Int(0);  // theta, zoffs, blank, tflag, t1..t4 (14 words, NX)
        o.Table("GEOM2S", {quad, o.Key(65535, 65535, 65535)});
        std::string pshell = o.Key(2302, 23, 283) + o.Int(7);
        for (int i = 0; i < 10; ++i)
            pshell += o.Int(i == 0 ? 1 : 0);
        o.Table("EPTS", {pshell, o.Key(65535, 65535, 65535)});
    }
    std::string disp;
    for (int n = 1; n <= 4; ++n) {
        disp += o.Int(10 * n + 1) + o.Int(1);  // id*10+device, point type
        for (int c = 0; c < 6; ++c)
            disp += o.Real(n + 0.1 * c);
    }
    o.Table("OUGV1", {o.Header(11, 1, 0, 8, 0), disp});
    std::string stress = o.Int(10 * 10 + 1);
    for (int k = 1; k < 17; ++k)
        stress += o.Real(100.0 + k);
    o.Table("OES1X1", {o.Header(11, 5, 33, 17, 1), stress});
    return o.mBytes;
}

// A static step with a spring (CELAS1 20 between GRIDs 1 and 2) and its force,
// an energy table (word 5 = 0), a frequency response (magnitude/phase) and a
// transient SORT2 table; the header's POST,-1 block optional.
std::string op2_more(bool Header) {
    Op2 o{4, false, {}};
    if (Header) {
        o.Marker(3);
        o.Block(o.Int(9) + o.Int(24) + o.Int(26));
        o.Marker(7);
        o.Block(std::string("NASTRAN FORT TAPE ID CODE - "));
        o.Record(o.Text("NX2019.2"));
        o.Marker(-1);
        o.Marker(0);
    }
    std::string grids = o.Key(4501, 45, 1);
    for (int n = 0; n < 2; ++n)
        grids += o.Int(n + 1) + o.Int(0) + o.Real(n) + o.Real(0.0) + o.Real(0.0) + o.Int(0) +
                 o.Int(0) + o.Int(0);
    o.Table("GEOM1S", {grids, o.Key(65535, 65535, 65535)});
    const std::string spring =
        o.Key(601, 6, 73) + o.Int(20) + o.Int(3) + o.Int(1) + o.Int(2) + o.Int(1) + o.Int(1);
    o.Table("GEOM2S", {spring, o.Key(65535, 65535, 65535)});
    o.Table("OEF1X", {o.Header(11, 4, 11, 2, 0), o.Int(201) + o.Real(5.5)});
    o.Table("ONRGY1",
            {o.Header(11, 18, 0, 4, 0, 0), o.Int(201) + o.Real(2.0) + o.Real(100.0) + o.Real(0.5)});
    // Frequency response at 10 Hz: T1 = 2 at 90 degrees (real 0, imaginary 2).
    std::string cplx = o.Int(11) + o.Int(1) + o.Real(2.0);
    for (int c = 1; c < 6; ++c)
        cplx += o.Real(0.0);
    cplx += o.Real(90.0);
    for (int c = 1; c < 6; ++c)
        cplx += o.Real(0.0);
    const int ten = static_cast<int>(std::bit_cast<std::uint32_t>(10.0f));
    o.Table("OUGV1", {o.Header(51, 1001, 0, 14, 0, ten, 3), cplx});
    // SORT2: GRID 2 over two times, rows led by the time.
    std::string sort2;
    for (float t : {0.5f, 1.5f}) {
        sort2 += o.Real(t) + o.Int(1) + o.Real(3.0 * t);
        for (int c = 1; c < 6; ++c)
            sort2 += o.Real(0.0);
    }
    o.Table("OUGV2", {o.Header(61, 2001, 0, 8, 0, 21), sort2});
    return o.mBytes;
}

std::string write_file(const std::string& rBody, const std::string& rName = "model.op2") {
    const std::filesystem::path dir(
        mt::temp_path("_op2_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    std::ofstream((dir / rName).string(), std::ios::binary) << rBody;
    return (dir / rName).string();
}

void expect_model(const Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 4u);
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "quad");
    EXPECT_EQ(rMesh.CellData("nastran:eid", 0).As<std::int64_t>()[0], 10);
    const double* d = rMesh.PointData("DISPLACEMENT").As<double>();
    const double* r = rMesh.PointData("DISPLACEMENT_ROT").As<double>();
    for (int n = 0; n < 4; ++n)
        for (int c = 0; c < 3; ++c) {
            EXPECT_NEAR(d[3 * n + c], n + 1 + 0.1 * c, 1e-6);
            EXPECT_NEAR(r[3 * n + c], n + 1 + 0.1 * (c + 3), 1e-6);
        }
    // fibre 1: FD1 X1 Y1 TXY1 ANGLE1 MAJOR1 MINOR1 VON_MISES1 = 101..108
    EXPECT_NEAR(rMesh.CellData("STRESS:X1", 0).As<double>()[0], 102.0, 1e-4);
    EXPECT_NEAR(rMesh.CellData("STRESS:VON_MISES2", 0).As<double>()[0], 116.0, 1e-4);
    EXPECT_EQ(rMesh.FieldData("nastran:subcase").As<std::int64_t>()[0], 1);
}

}  // namespace

TEST(NastranOp2, WordSizesAndByteOrdersReadTheSame) {
    for (const auto& [ws, big] : {std::pair{4, false}, std::pair{4, true}, std::pair{8, false}}) {
        const std::string path = write_file(op2(ws, big, true));
        const Mesh mesh = meshioplusplus::read_nastran_op2(path);
        expect_model(mesh);
        ASSERT_EQ(mesh.NumRegions(), 1u);
        EXPECT_EQ(mesh.Region(0).mName, "PSHELL_7");
        EXPECT_EQ(meshioplusplus::nastran_op2_time_values(path), std::vector<double>{0.0});
    }
}

TEST(NastranOp2, SiblingDeckGivesTheMesh) {
    const std::string path = write_file(op2(4, false, false));
    EXPECT_THROW(meshioplusplus::read_nastran_op2(path), ReadError);
    std::ofstream(std::filesystem::path(path).replace_extension(".bdf").string())
        << "BEGIN BULK\n"
           "GRID,1,,0.,0.,0.\n"
           "GRID,2,,1.,0.,0.\n"
           "GRID,3,,1.,1.,0.\n"
           "GRID,4,,0.,1.,0.\n"
           "CQUAD4,10,7,1,2,3,4\n"
           "ENDDATA\n";
    expect_model(meshioplusplus::read_nastran_op2(path));
}

TEST(NastranOp2, ResolvedSniffedAndNarrowed) {
    const std::string path = write_file(op2(4, false, true));
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "nastran_op2");
    const std::string renamed = write_file(op2(8, true, true), "results.bin");
    EXPECT_EQ(meshioplusplus::sniff_format(renamed), "nastran_op2");
    const auto meta = meshioplusplus::registry_read_metadata(path, "nastran_op2", ReadOptions{});
    EXPECT_EQ(meta.mTimeValues, std::vector<double>{0.0});
    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"STRESS:X1"};
    const Mesh narrowed = meshioplusplus::read_nastran_op2(path, only);
    EXPECT_TRUE(narrowed.HasCellData("STRESS:X1"));
    EXPECT_FALSE(narrowed.HasPointData("DISPLACEMENT"));
    ReadOptions late;
    late.mTimeStep = 1;
    EXPECT_THROW(meshioplusplus::read_nastran_op2(path, late), ReadError);
    EXPECT_THROW(meshioplusplus::read_nastran_op2(write_file(op2(4, false, true).substr(0, 400))),
                 ReadError);
}

TEST(NastranOp2, SpringsEnergiesComplexAndSort2) {
    const std::string path = write_file(op2_more(true));
    const auto times = meshioplusplus::nastran_op2_time_values(path);
    ASSERT_EQ(times, (std::vector<double>{0.0, 10.0, 0.5, 1.5}));
    const Mesh s = meshioplusplus::read_nastran_op2(path);
    ASSERT_EQ(s.NumCellBlocks(), 1u);
    EXPECT_EQ(s.Cells(0).Type(), "line");  // the spring
    EXPECT_NEAR(s.CellData("ELEMENT_FORCE:F", 0).As<double>()[0], 5.5, 1e-6);
    EXPECT_NEAR(s.CellData("ENERGY:ENERGY", 0).As<double>()[0], 2.0, 1e-6);  // joined the step
    ReadOptions freq;
    freq.mTimeStep = 1;
    const Mesh f = meshioplusplus::read_nastran_op2(path, freq);
    EXPECT_NEAR(f.PointData("DISPLACEMENT_real").As<double>()[0], 0.0, 1e-6);
    EXPECT_NEAR(f.PointData("DISPLACEMENT_imag").As<double>()[0], 2.0, 1e-6);
    ReadOptions last;
    last.mTimeStep = -1;
    const Mesh t = meshioplusplus::read_nastran_op2(path, last);
    EXPECT_NEAR(t.PointData("DISPLACEMENT").As<double>()[3], 4.5, 1e-6);  // GRID 2, T1
    EXPECT_TRUE(std::isnan(t.PointData("DISPLACEMENT").As<double>()[0]));
    // PARAM,POST,-2: no header, sniffed by the first table's name.
    const std::string post2 = write_file(op2_more(false), "post2.bin");
    EXPECT_EQ(meshioplusplus::sniff_format(post2), "nastran_op2");
    EXPECT_EQ(meshioplusplus::nastran_op2_time_values(post2), times);
}

// ---------------------------------------------------------------------------------------
// The writer
// ---------------------------------------------------------------------------------------

namespace {

using meshioplusplus::NDArray;

NDArray doubles(std::initializer_list<double> Values, std::vector<std::size_t> Shape) {
    NDArray a(meshioplusplus::DType::Float64, std::move(Shape));
    std::copy(Values.begin(), Values.end(), a.As<double>());
    return a;
}

NDArray ints(std::initializer_list<std::int64_t> Values) {
    NDArray a(meshioplusplus::DType::Int64, {Values.size()});
    std::copy(Values.begin(), Values.end(), a.As<std::int64_t>());
    return a;
}

std::string read_all(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string temp_op2(const std::string& rName = "written.op2") {
    const std::filesystem::path dir(
        mt::temp_path("_op2_w_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    return (dir / rName).string();
}

// The one-quad model of `oracle_mesh`: GEOM1/GEOM2/EPT/MPT, the displacement
// table (a header of blanks after word 50) and the plate stress table.
std::string oracle_bytes() {
    Op2 o{4, false, {}};
    o.Marker(3);
    o.Block(o.Int(1) + o.Int(1) + o.Int(0));
    o.Marker(7);
    o.Block("NASTRAN FORT TAPE ID CODE - ");
    o.Record(o.Text("XXXXXXXX"));
    o.Marker(-1);
    o.Marker(0);
    const double xy[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    std::string grids = o.Key(4501, 45, 1);
    for (int n = 0; n < 4; ++n)
        grids += o.Int(n + 1) + o.Int(0) + o.Real(xy[n][0]) + o.Real(xy[n][1]) + o.Real(0.0) +
                 o.Int(0) + o.Int(0) + o.Int(0);
    const std::string end = o.Key(65535, 65535, 65535);
    o.Table("GEOM1", {grids, end});
    std::string quad = o.Key(2958, 51, 177) + o.Int(10) + o.Int(7);
    for (int n = 1; n <= 4; ++n)
        quad += o.Int(n);
    for (int i = 0; i < 8; ++i)
        quad += o.Int(0);
    o.Table("GEOM2", {quad, end});
    std::string pshell = o.Key(2302, 23, 283) + o.Int(7) + o.Int(1) + o.Real(1.0);
    for (int i = 0; i < 8; ++i)
        pshell += o.Int(0);
    o.Table("EPT", {pshell, end});
    std::string mat1 = o.Key(103, 1, 77) + o.Int(1) + o.Real(1.0) + o.Real(0.0) + o.Real(0.3);
    for (int i = 0; i < 8; ++i)
        mat1 += o.Int(0);
    o.Table("MPT", {mat1, end});
    // A 146-word header whose words from 50 on are blanks.
    const auto header = [&](int TCode, int EType, int NumWide, int SCode) {
        std::string h = o.Header(11, TCode, EType, NumWide, SCode, 1);
        std::string out = h.substr(0, 50 * 4);
        for (int i = 50; i < 146; ++i)
            out += o.Int(0x20202020);
        return out;
    };
    std::string disp;
    for (int n = 1; n <= 4; ++n) {
        disp += o.Int(10 * n + 1) + o.Int(1);
        for (int c = 0; c < 6; ++c)
            disp += o.Real(n + 0.5 * c);
    }
    o.Table("OUGV1", {header(1, 0, 8, 0), disp});
    std::string stress = o.Int(10 * 10 + 1);
    for (int k = 0; k < 16; ++k)
        stress += o.Real(100.0 + k);
    o.Table("OES1X1", {header(5, 33, 17, 1), stress});
    o.Marker(0);
    return o.mBytes;
}

Mesh oracle_mesh() {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});
    std::vector<double> disp, rot;
    for (int n = 1; n <= 4; ++n)
        for (int c = 0; c < 3; ++c) {
            disp.push_back(n + 0.5 * c);
            rot.push_back(n + 0.5 * (c + 3));
        }
    NDArray d(meshioplusplus::DType::Float64, {std::size_t{4}, std::size_t{3}});
    std::copy(disp.begin(), disp.end(), d.As<double>());
    NDArray r(meshioplusplus::DType::Float64, {std::size_t{4}, std::size_t{3}});
    std::copy(rot.begin(), rot.end(), r.As<double>());
    m.AddPointData("DISPLACEMENT", std::move(d));
    m.AddPointData("DISPLACEMENT_ROT", std::move(r));
    m.AddCellData("nastran:eid", {ints({10})});
    m.AddCellData("nastran:pid", {ints({7})});
    static const char* const names[] = {
        "FD1", "X1", "Y1", "TXY1", "ANGLE1", "MAJOR1", "MINOR1", "VON_MISES1",
        "FD2", "X2", "Y2", "TXY2", "ANGLE2", "MAJOR2", "MINOR2", "VON_MISES2"};
    for (int k = 0; k < 16; ++k)
        m.AddCellData(std::string("STRESS:") + names[k], {doubles({100.0 + k}, {1})});
    return m;
}

// Every card: one cell of each type over 30 points.
Mesh every_card_mesh() {
    Mesh m;
    NDArray pts(meshioplusplus::DType::Float64, {std::size_t{30}, std::size_t{3}});
    for (int i = 0; i < 90; ++i)
        pts.As<double>()[i] = std::fmod(0.37 * i + 0.11 * (i % 7), 1.0);
    m.AssignPoints(std::move(pts));
    const std::vector<std::pair<std::string, int>> types = {
        {"vertex", 1},     {"line", 2},  {"triangle", 3}, {"triangle6", 6},  {"quad", 4},
        {"quad8", 8},      {"quad9", 9}, {"tetra", 4},    {"tetra10", 10},   {"pyramid", 5},
        {"pyramid13", 13}, {"wedge", 6}, {"wedge15", 15}, {"hexahedron", 8}, {"hexahedron20", 20}};
    std::vector<NDArray> pids;
    for (const auto& [type, nodes] : types) {
        std::vector<std::vector<std::int64_t>> row(1);
        for (int k = 0; k < nodes; ++k)
            row[0].push_back((k + 2) % 30);
        m.AddCellBlock(type, mt::conn_from(row));
        const bool shell = type.rfind("tri", 0) == 0 || type.rfind("quad", 0) == 0;
        pids.push_back(ints({type == "vertex" ? -1 : type == "line" ? 4 : shell ? 1 : 2}));
    }
    m.AddCellData("nastran:pid", std::move(pids));
    return m;
}

Mesh step_mesh(double Scale, std::int64_t Analysis, std::int64_t Subcase = 1) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}},
                           "line", {{4, 5}});
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 3, 4}}));
    m.AddCellData("nastran:pid", {ints({4}), ints({1}), ints({2})});
    NDArray d(meshioplusplus::DType::Float64, {std::size_t{6}, std::size_t{3}});
    for (int i = 0; i < 18; ++i)
        d.As<double>()[i] = Scale * i / 4.0;
    m.AddPointData("DISPLACEMENT", std::move(d));
    m.AddCellData("STRESS:A", {doubles({5.5 * Scale}, {1}), doubles({std::nan("")}, {1}),
                               doubles({std::nan("")}, {1})});
    m.AddCellData("STRESS:X1", {doubles({std::nan("")}, {1}), doubles({7.0}, {1}),
                                doubles({std::nan("")}, {1})});
    m.AddFieldData("nastran:subcase", ints({Subcase}));
    m.AddFieldData("nastran:analysis", ints({Analysis}));
    m.AddFieldData("nastran:mode", ints({static_cast<std::int64_t>(Scale)}));
    return m;
}

}  // namespace

TEST(NastranOp2Writer, FollowsTheHandFramedLayout) {
    const std::string path = temp_op2();
    meshioplusplus::write_nastran_op2(path, oracle_mesh());
    EXPECT_EQ(read_all(path), oracle_bytes());
}

TEST(NastranOp2Writer, ReadsBackWhatTheOracleFileHolds) {
    const std::string path = temp_op2();
    meshioplusplus::write_nastran_op2(path, oracle_mesh());
    const Mesh back = meshioplusplus::read_nastran_op2(path);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "quad");
    EXPECT_EQ(back.CellData("nastran:eid", 0).As<std::int64_t>()[0], 10);
    EXPECT_EQ(back.Region(0).mName, "PSHELL_7");
    const double* d = back.PointData("DISPLACEMENT").As<double>();
    const double* r = back.PointData("DISPLACEMENT_ROT").As<double>();
    for (int n = 0; n < 4; ++n)
        for (int c = 0; c < 3; ++c) {
            EXPECT_NEAR(d[3 * n + c], n + 1 + 0.5 * c, 1e-6);
            EXPECT_NEAR(r[3 * n + c], n + 1 + 0.5 * (c + 3), 1e-6);
        }
    EXPECT_NEAR(back.CellData("STRESS:X1", 0).As<double>()[0], 101.0, 1e-4);
    EXPECT_NEAR(back.CellData("STRESS:VON_MISES2", 0).As<double>()[0], 115.0, 1e-4);
}

TEST(NastranOp2Writer, EveryCardSurvivesARoundTrip) {
    const Mesh mesh = every_card_mesh();
    const std::string path = temp_op2();
    meshioplusplus::write_nastran_op2(path, mesh);
    const Mesh back = meshioplusplus::read_nastran_op2(path);
    ASSERT_EQ(back.NumCellBlocks(), mesh.NumCellBlocks());
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
        const std::string type = mesh.Cells(b).Type();
        std::size_t found = back.NumCellBlocks();
        for (std::size_t k = 0; k < back.NumCellBlocks(); ++k)
            if (back.Cells(k).Type() == type)
                found = k;
        ASSERT_LT(found, back.NumCellBlocks()) << type;
        const auto want = mesh.Cells(b).Conn();
        const auto got = back.Cells(found).Conn();
        ASSERT_EQ(want.Size(), got.Size()) << type;
        for (std::size_t i = 0; i < want.Size(); ++i)
            EXPECT_EQ(want.As<std::int64_t>()[i], got.As<std::int64_t>()[i]) << type << " " << i;
    }
    EXPECT_EQ(back.NumRegions(), 3u);
}

TEST(NastranOp2Writer, ResultsComeBackByElementAndPoint) {
    const std::string path = temp_op2();
    meshioplusplus::write_nastran_op2(path, step_mesh(2.0, 1));
    const Mesh back = meshioplusplus::read_nastran_op2(path);
    ASSERT_EQ(back.NumCellBlocks(), 3u);
    for (std::size_t b = 0; b < 3; ++b) {
        const std::string type = back.Cells(b).Type();
        if (type == "line")
            EXPECT_EQ(back.CellData("STRESS:A", b).As<double>()[0], 11.0);
        if (type == "quad")
            EXPECT_EQ(back.CellData("STRESS:X1", b).As<double>()[0], 7.0);
    }
    const double* d = back.PointData("DISPLACEMENT").As<double>();
    for (int i = 0; i < 18; ++i)
        EXPECT_NEAR(d[i], 2.0 * i / 4.0, 1e-6);
}

TEST(NastranOp2Writer, ASeriesIsTheModelOnceThenABlockPerStep) {
    const std::string path = temp_op2();
    meshioplusplus::NastranOp2SeriesWriter writer(path);
    writer.Write(1.0, step_mesh(1.0, 6));
    // other result arrays, the same step key, and other cells are refused
    Mesh fewer = step_mesh(2.0, 6);
    Mesh stripped = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}}, "line", {{4, 5}});
    stripped.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));
    stripped.AddCellBlock("tetra", mt::conn_from({{0, 1, 3, 4}}));
    stripped.AddFieldData("nastran:analysis", ints({6}));
    EXPECT_THROW(writer.Write(2.0, stripped), WriteError);
    EXPECT_THROW(writer.Write(1.0, step_mesh(1.0, 6)), WriteError);
    EXPECT_THROW(writer.Write(2.0, mt::tri_mesh()), WriteError);
    EXPECT_EQ(writer.NumSteps(), 1u);
    writer.Write(2.0, std::move(fewer));
    writer.Write(4.0, step_mesh(4.0, 6));
    writer.Finalize();
    writer.Finalize();
    EXPECT_EQ(meshioplusplus::nastran_op2_time_values(path), (std::vector<double>{1.0, 2.0, 4.0}));
    ReadOptions last;
    last.mTimeStep = -1;
    const Mesh back = meshioplusplus::read_nastran_op2(path, last);
    EXPECT_EQ(back.FieldData("meshio:time").As<double>()[0], 4.0);
    EXPECT_NEAR(back.PointData("DISPLACEMENT").As<double>()[17], 4.0 * 17 / 4.0, 1e-6);

    meshioplusplus::NastranOp2SeriesWriter empty(temp_op2());
    EXPECT_THROW(empty.Finalize(), WriteError);
}

TEST(NastranOp2Writer, ASeriesIsOneStepsBytesWhenItHoldsOne) {
    Mesh mesh = step_mesh(2.0, 6);
    mesh.AddFieldData("meshio:time", doubles({1.5}, {1}));
    const std::string plain = temp_op2(), series = temp_op2();
    meshioplusplus::write_nastran_op2(plain, mesh);
    meshioplusplus::NastranOp2SeriesWriter writer(series);
    writer.Write(1.5, mesh);
    writer.Finalize();
    EXPECT_EQ(read_all(plain), read_all(series));
}

TEST(NastranOp2Writer, ASeriesWithNoResultsIsTheModel) {
    const std::string path = temp_op2();
    meshioplusplus::NastranOp2SeriesWriter writer(path);
    writer.Write(0.0, every_card_mesh());
    writer.Write(0.0, every_card_mesh());
    writer.Finalize();
    EXPECT_TRUE(meshioplusplus::nastran_op2_time_values(path).empty());
    EXPECT_EQ(meshioplusplus::read_nastran_op2(path).NumCellBlocks(), 15u);
}

TEST(NastranOp2Writer, RefusalsAndDrops) {
    Mesh repeated = step_mesh(1.0, 1);
    repeated.AddCellData("nastran:eid", {ints({1}), ints({1}), ints({2})});
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), repeated), WriteError);
    Mesh nonpositive = step_mesh(1.0, 1);
    nonpositive.AddCellData("nastran:eid", {ints({0}), ints({1}), ints({2})});
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), nonpositive), WriteError);
    Mesh shared = step_mesh(1.0, 1);
    shared.AddCellData("nastran:pid", {ints({4}), ints({1}), ints({1})});
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), shared), WriteError);
    Mesh narrow = step_mesh(1.0, 1);
    narrow.AddPointData("DISPLACEMENT", doubles({1, 2, 3, 4}, {2, 2}));
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), narrow), WriteError);
    Mesh outside = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, "triangle", {{0, 1, 3}});
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), outside), WriteError);
    Mesh code = step_mesh(1.0, 0);
    EXPECT_THROW(meshioplusplus::write_nastran_op2(temp_op2(), code), WriteError);
    EXPECT_THROW(meshioplusplus::write_nastran_op2("/no/such/dir/x.op2", step_mesh(1.0, 1)),
                 WriteError);
    // Arrays with no table are dropped, not written, and the model alone is a file.
    Mesh extra = every_card_mesh();
    NDArray p(meshioplusplus::DType::Float64, {std::size_t{30}});
    extra.AddPointData("pressure", std::move(p));
    const std::string path = temp_op2();
    meshioplusplus::write_nastran_op2(path, extra);
    const Mesh back = meshioplusplus::read_nastran_op2(path);
    EXPECT_FALSE(back.HasPointData("pressure"));
    EXPECT_TRUE(meshioplusplus::nastran_op2_time_values(path).empty());
}

TEST(NastranOp2Writer, TheRegistryAndTheSequenceDriverWriteIt) {
    EXPECT_EQ(meshioplusplus::resolve_write_format("/a/model.op2", ""), "nastran_op2");
    const std::string path = temp_op2("by_name.bin");
    meshioplusplus::registry_writers().at("nastran_op2")(path, oracle_mesh());
    EXPECT_EQ(read_all(path), oracle_bytes());
    EXPECT_EQ(meshioplusplus::sniff_format(path), "nastran_op2");
    std::string why;
    EXPECT_TRUE(meshioplusplus::sequence_write_supports_time("nastran_op2", why));
}
