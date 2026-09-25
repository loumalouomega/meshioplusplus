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
 *        deck route, steps, metadata, sniffing and refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

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
