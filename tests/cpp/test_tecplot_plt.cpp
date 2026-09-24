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
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

namespace {

// A minimal #!TDV112 writer: the fixtures under tests/python/meshes/tecplot/plt
// are TecIO-written; these hand-made files reach the data formats, byte orders
// and refusals TecIO's classic API never writes.
class PltWriter {
public:
    explicit PltWriter(bool BigEndian) : mBig(BigEndian) {}

    void Raw(const std::string& rS) { mOut += rS; }
    void I32(std::int32_t V) { Put(&V, 4); }
    void F32(float V) { Put(&V, 4); }
    void F64(double V) { Put(&V, 8); }
    void I16(std::int16_t V) { Put(&V, 2); }
    void U8(std::uint8_t V) { mOut += static_cast<char>(V); }
    void Str(const std::string& rS) {
        for (char c : rS)
            I32(c);
        I32(0);
    }
    std::string Write(const std::string& rSuffix) const {
        const std::string path = mt::temp_path(rSuffix);
        std::ofstream f(path, std::ios::binary);
        f << mOut;
        return path;
    }
    const std::string& Bytes() const { return mOut; }

private:
    void Put(const void* pV, std::size_t N) {
        char b[8];
        std::memcpy(b, pV, N);
        if (mBig)
            for (std::size_t i = 0; i < N / 2; ++i)
                std::swap(b[i], b[N - 1 - i]);
        mOut.append(b, N);
    }
    bool mBig;
    std::string mOut;
};

// One FE zone (X, Y, V) holding a single triangle, V in data format `Fmt`.
std::string tri_plt(int Fmt, bool BigEndian, int ZoneType = 2,
                    const std::string& rVersion = "112") {
    PltWriter w(BigEndian);
    w.Raw("#!TDV" + rVersion);
    w.I32(1);
    w.I32(0);
    w.Str("t");
    w.I32(3);
    w.Str("X");
    w.Str("Y");
    w.Str("V");
    w.F32(299.0F);
    w.Str("tri");
    w.I32(-1);  // parent
    w.I32(-1);  // strand
    w.F64(0.0);
    w.I32(-1);
    w.I32(ZoneType);
    w.I32(0);  // no var location
    w.I32(0);  // raw face neighbours
    w.I32(0);  // misc face neighbours
    w.I32(3);  // points
    if (ZoneType == 6 || ZoneType == 7) {
        w.I32(3);
        w.I32(6);
        w.I32(0);
        w.I32(0);
    }
    w.I32(1);  // elements
    w.I32(0);
    w.I32(0);
    w.I32(0);
    w.I32(0);  // no aux
    w.F32(357.0F);
    w.F32(299.0F);
    w.I32(2);
    w.I32(1);
    w.I32(Fmt);
    w.I32(0);   // no passive
    w.I32(0);   // no sharing
    w.I32(-1);  // own connectivity
    for (int k = 0; k < 6; ++k)
        w.F64(0.0);  // min/max
    for (double x : {0.0, 1.0, 0.0})
        w.F64(x);
    for (float y : {0.0F, 0.0F, 1.0F})
        w.F32(y);
    for (int v : {3, 250, 7}) {
        if (Fmt == 1)
            w.F32(static_cast<float>(v));
        else if (Fmt == 2)
            w.F64(v);
        else if (Fmt == 3)
            w.I32(v);
        else if (Fmt == 4)
            w.I16(static_cast<std::int16_t>(v));
        else if (Fmt == 5)
            w.U8(static_cast<std::uint8_t>(v));
    }
    w.I32(0);
    w.I32(1);
    w.I32(2);
    return w.Bytes();
}

std::string write_bytes(const std::string& rBytes, const std::string& rSuffix = ".plt") {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream f(path, std::ios::binary);
    f << rBytes;
    return path;
}

void remove_file(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
}

}  // namespace

TEST(TecplotPlt, EveryDataFormatInBothByteOrders) {
    for (bool big : {false, true}) {
        for (int fmt : {1, 2, 3, 4, 5}) {
            const std::string path = write_bytes(tri_plt(fmt, big));
            EXPECT_EQ(meshioplusplus::sniff_format(path), "tecplot");
            const mt::Mesh m = meshioplusplus::read_tecplot(path);
            ASSERT_EQ(m.NumCellBlocks(), 1u);
            EXPECT_EQ(m.Cells(0).Type(), "triangle");
            EXPECT_EQ(m.Cells(0).Conn().As<std::int64_t>()[2], 2);
            EXPECT_EQ(m.Points().Shape()[1], 2u);
            EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.Points(), 2),
                             1.0);  // x of node 1
            EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.PointData("V"), 1), 250.0)
                << "format " << fmt << " big-endian " << big;
            remove_file(path);
        }
    }
}

TEST(TecplotPlt, RefusesWhatItCannotRead) {
    const std::string bit = write_bytes(tri_plt(6, false));
    EXPECT_THROW(meshioplusplus::read_tecplot(bit), meshioplusplus::ReadError);
    const std::string poly = write_bytes(tri_plt(2, false, 6));
    try {
        meshioplusplus::read_tecplot(poly);
        ADD_FAILURE() << "an FEPOLYGON zone was read";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("polygonal"), std::string::npos);
    }
    const std::string old = write_bytes(tri_plt(2, false, 2, "102"));
    EXPECT_THROW(meshioplusplus::read_tecplot(old), meshioplusplus::ReadError);
    const std::string good = tri_plt(2, false);
    const std::string cut = write_bytes(good.substr(0, good.size() - 10));
    EXPECT_THROW(meshioplusplus::read_tecplot(cut), meshioplusplus::ReadError);
    for (const std::string& p : {bit, poly, old, cut})
        remove_file(p);
}

TEST(TecplotPlt, OrderedZonesCellCentredGhostsAndTimeSteps) {
    // One IJK zone (3x2x2) with a cell-centred C stored with ghosts
    // (I*J*(K-1) = 6 values, 2 real), at SOLUTIONTIME 1, and the same zone
    // again at time 2 sharing X/Y/Z.
    PltWriter w(false);
    w.Raw("#!TDV112");
    w.I32(1);
    w.I32(0);
    w.Str("o");
    w.I32(4);
    for (const char* v : {"X", "Y", "Z", "C"})
        w.Str(v);
    for (double t : {1.0, 2.0}) {
        w.F32(299.0F);
        w.Str("block");
        w.I32(-1);
        w.I32(0);  // strand 0
        w.F64(t);
        w.I32(-1);
        w.I32(0);  // ORDERED
        w.I32(1);
        for (int v : {0, 0, 0, 1})
            w.I32(v);
        w.I32(0);
        w.I32(0);
        w.I32(3);
        w.I32(2);
        w.I32(2);
        w.I32(0);
    }
    w.F32(357.0F);
    for (int zone = 0; zone < 2; ++zone) {
        w.F32(299.0F);
        for (int k = 0; k < 4; ++k)
            w.I32(2);
        w.I32(0);
        w.I32(zone);  // zone 2 shares X, Y, Z from zone 1
        if (zone == 1)
            for (int v : {0, 0, 0, -1})
                w.I32(v);
        w.I32(-1);
        for (int k = 0; k < (zone == 0 ? 8 : 2); ++k)
            w.F64(0.0);
        if (zone == 0) {
            for (int axis = 0; axis < 3; ++axis)
                for (int n = 0; n < 12; ++n) {
                    const int i = n % 3, j = (n / 3) % 2, k = n / 6;
                    w.F64(axis == 0 ? i : axis == 1 ? j : k);
                }
        }
        // C over (I=3, J=2, K-1=1): real cells at i < 2, j < 1.
        const double base = zone == 0 ? 10.0 : 20.0;
        for (double v : {base, base + 1, -1.0, -1.0, -1.0, -1.0})
            w.F64(v);
    }
    const std::string path = w.Write(".plt");

    meshioplusplus::ReadOptions opts;
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_tecplot_metadata(path, opts);
    ASSERT_EQ(meta.mTimeValues.size(), 2u);
    EXPECT_EQ(meta.mNumPoints, 12u);
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "hexahedron");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 2u);

    for (int step = 0; step < 2; ++step) {
        opts.mTimeStep = step;
        const mt::Mesh m = meshioplusplus::read_tecplot(path, opts);
        ASSERT_EQ(m.NumCellBlocks(), 1u);
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        const std::vector<std::int64_t> first(c, c + 8);
        EXPECT_EQ(first, (std::vector<std::int64_t>{0, 1, 4, 3, 6, 7, 10, 9}));
        const auto& cd = m.CellData("C", 0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(cd, 0), step == 0 ? 10.0 : 20.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(cd, 1), step == 0 ? 11.0 : 21.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.Points(), 3 * 11 + 0), 2.0);
    }
    remove_file(path);
}

TEST(TecplotPlt, AsciiOrderedFormsMatchTheirConnectivity) {
    const std::string path = mt::temp_path(".dat");
    {
        std::ofstream f(path);
        f << "VARIABLES = \"X\" \"Y\" \"P\"\n"
             "ZONE T=\"j\", I=1, J=3, F=POINT\n"
             "0 0 1\n0 1 2\n0 2 3\n"
             "ZONE T=\"ij\", I=2, J=2, DATAPACKING=POINT\n"
             "0 0 4\n1 0 5\n0 1 6\n1 1 7\n"
             "ZONE T=\"ij2\", I=2, J=2, VARSHARELIST=([1-2])\n"
             "8 9 10 11\n";
    }
    const mt::Mesh m = meshioplusplus::read_tecplot(path);
    ASSERT_EQ(m.NumCellBlocks(), 3u);
    EXPECT_EQ(m.Cells(0).Type(), "line");
    EXPECT_EQ(m.Cells(1).Type(), "quad");
    const std::int64_t* q = m.Cells(1).Conn().As<std::int64_t>();
    EXPECT_EQ(std::vector<std::int64_t>(q, q + 4), (std::vector<std::int64_t>{3, 4, 6, 5}));
    EXPECT_EQ(m.Points().Shape()[0], 11u);  // ij2 owns its P, so its own points
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.PointData("P"), 10), 11.0);
    ASSERT_EQ(m.NumRegions(), 3u);
    EXPECT_TRUE(m.HasRegion("ij2"));
    remove_file(path);
}
