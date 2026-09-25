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
    // An FEPOLYGON zone whose header announces a face map the file lacks.
    const std::string poly = write_bytes(tri_plt(2, false, 6));
    EXPECT_THROW(meshioplusplus::read_tecplot(poly), meshioplusplus::ReadError);
    const std::string old = write_bytes(tri_plt(2, false, 2, "99 "));
    EXPECT_THROW(meshioplusplus::read_tecplot(old), meshioplusplus::ReadError);
    const std::string good = tri_plt(2, false);
    const std::string cut = write_bytes(good.substr(0, good.size() - 10));
    EXPECT_THROW(meshioplusplus::read_tecplot(cut), meshioplusplus::ReadError);
    for (const std::string& p : {bit, poly, old, cut})
        remove_file(p);
}

// A Tecplot 7.5 file (the layout VisIt's test files show): no FileType, a
// 2-D text and a circle in the Tecplot 7 record layouts, an FE-point zone
// (kind 3) whose element type follows its sizes, one extra word before the
// formats and before the 1-based connectivity, no min/max pairs.
TEST(TecplotPlt, Tecplot7PointPackedFeZone) {
    for (bool big : {false, true}) {
        PltWriter w(big);
        w.Raw("#!TDV75 ");
        w.I32(1);
        w.Str("old");
        w.I32(3);
        w.Str("X(M)");  // unit suffixes still name the coordinates
        w.Str("Y(M)");
        w.Str("P");
        w.F32(499.0F);  // text
        w.I32(1);
        w.I32(1);
        w.F64(40.0);
        w.F64(2.5);
        w.I32(0);
        w.I32(1);
        w.F64(3.0);
        w.I32(0);
        w.F64(15.0);
        w.F64(0.1);
        w.I32(0);
        w.I32(7);
        w.F64(90.0);
        w.F64(1.5);
        for (int v : {0, 0, -1, 0})
            w.I32(v);
        w.Str("label");
        w.F32(399.0F);  // a circle
        w.I32(0);
        w.I32(1);
        w.F64(1.0);
        w.F64(2.0);
        for (int v : {0, 0, -1, 0, 7, 0})
            w.I32(v);
        w.I32(3);  // circle
        w.I32(0);  // line pattern
        w.F64(2.0);
        w.F64(0.1);
        w.I32(72);
        w.I32(0);
        w.I32(0);
        w.F64(1.0);
        w.F64(15.0);
        w.Str("");
        w.I32(1);  // float
        w.F32(0.5F);
        w.F32(299.0F);
        w.Str("zone");
        w.I32(3);  // FE point
        w.I32(-1);
        w.I32(4);  // nodes
        w.I32(2);  // elements
        w.I32(1);  // quadrilateral
        w.F32(357.0F);
        w.F32(299.0F);
        w.I32(0);
        for (int k = 0; k < 3; ++k)
            w.I32(1);
        const float xyz[4][3] = {{0, 0, 1}, {1, 0, 2}, {1, 1, 3}, {0, 1, 4}};
        for (const auto& r : xyz)
            for (float v : r)
                w.F32(v);
        w.I32(0);
        for (int v : {1, 2, 3, 4, 1, 3, 4, 2})
            w.I32(v);
        const std::string path = write_bytes(w.Bytes());
        const mt::Mesh m = meshioplusplus::read_tecplot(path);
        ASSERT_EQ(m.NumCellBlocks(), 1u);
        EXPECT_EQ(m.Cells(0).Type(), "quad");
        const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
        EXPECT_EQ(c[0], 0);
        EXPECT_EQ(c[7], 1);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.Points(), 4), 1.0);  // y of node 2
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.PointData("P"), 3), 4.0);
        remove_file(path);
    }
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

namespace {

std::vector<std::int64_t> face_of(const mt::Mesh& rMesh, std::size_t Block, std::size_t Cell,
                                  std::size_t Face) {
    const auto [ptr, size] = rMesh.Cells(Block).Face(Cell, Face);
    return std::vector<std::int64_t>(ptr, ptr + size);
}

}  // namespace

TEST(TecplotPoly, PolyhedralZoneSplitsByNodeCountWithOutwardFaces) {
    // A unit cube and a pyramid on its top face: the shared face is wound
    // outward for the cube (its left element), so the pyramid (the right
    // element) gets it reversed.
    const std::string path = mt::temp_path(".dat");
    {
        std::ofstream f(path);
        f << "VARIABLES = \"X\" \"Y\" \"Z\" \"C\"\n"
             "ZONE T=\"cells\", ZONETYPE=FEPOLYHEDRON, NODES=9, ELEMENTS=2, FACES=10,\n"
             "TOTALNUMFACENODES=36, NUMCONNECTEDBOUNDARYFACES=0, TOTALNUMBOUNDARYCONNECTIONS=0,\n"
             "VARLOCATION=([4]=CELLCENTERED)\n"
             "0 1 1 0 0 1 1 0 0.5\n0 0 1 1 0 0 1 1 0.5\n0 0 0 0 1 1 1 1 2\n10 20\n"
             "4 4 4 4 4 4 3 3 3 3\n"
             "1 4 3 2\n5 6 7 8\n1 2 6 5\n2 3 7 6\n3 4 8 7\n4 1 5 8\n"
             "5 6 9\n6 7 9\n7 8 9\n8 5 9\n"
             "1 1 1 1 1 1 2 2 2 2\n"
             "0 2 0 0 0 0 0 0 0 0\n";
    }
    const mt::Mesh m = meshioplusplus::read_tecplot(path);
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    EXPECT_EQ(m.Cells(0).Type(), "polyhedron8");
    EXPECT_EQ(m.Cells(1).Type(), "polyhedron5");
    ASSERT_TRUE(m.Cells(1).IsPolyhedron());
    EXPECT_EQ(m.Cells(0).NumFaces(0), 6u);
    EXPECT_EQ(m.Cells(1).NumFaces(0), 5u);
    EXPECT_EQ(face_of(m, 0, 0, 1), (std::vector<std::int64_t>{4, 5, 6, 7}));
    EXPECT_EQ(face_of(m, 1, 0, 0), (std::vector<std::int64_t>{7, 6, 5, 4}));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.CellData("C", 1), 0), 20.0);
    ASSERT_EQ(m.NumRegions(), 1u);
    EXPECT_EQ(m.Region(0).NumEntries(), 2u);

    // Written back as an FEPOLYHEDRON zone per block, it reads the same.
    const std::string out = mt::temp_path(".dat");
    meshioplusplus::write_tecplot(out, m);
    const mt::Mesh back = meshioplusplus::read_tecplot(out);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(back.Cells(1).Type(), "polyhedron5");
    for (std::size_t b = 0; b < 2; ++b)
        for (std::size_t f = 0; f < m.Cells(b).NumFaces(0); ++f)
            EXPECT_EQ(face_of(back, b, 0, f), face_of(m, b, 0, f)) << "block " << b;
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.CellData("C", 1), 0), 20.0);
    remove_file(path);
    remove_file(out);
}

TEST(TecplotPoly, BinaryPolygonZoneReadsZeroBasedFaceMaps) {
    // A house: a square (0 1 2 4) and a roof triangle (4 2 3); the shared
    // edge 2-4 has the square on its left and the roof on its right. Binary
    // face maps are 0-based with -1 for "no neighbour".
    PltWriter w(false);
    w.Raw("#!TDV112");
    w.I32(1);
    w.I32(0);
    w.Str("house");
    w.I32(2);
    w.Str("X");
    w.Str("Y");
    w.F32(299.0F);
    w.Str("poly");
    w.I32(-1);
    w.I32(-1);
    w.F64(0.0);
    w.I32(-1);
    w.I32(6);  // FEPOLYGON
    w.I32(0);
    w.I32(0);
    w.I32(0);
    w.I32(5);   // points
    w.I32(6);   // faces
    w.I32(12);  // face nodes
    w.I32(0);   // boundary faces
    w.I32(0);   // boundary connections
    w.I32(2);   // elements
    for (int k = 0; k < 3; ++k)
        w.I32(0);
    w.I32(0);  // no aux
    w.F32(357.0F);
    w.F32(299.0F);
    w.I32(2);
    w.I32(2);
    w.I32(0);
    w.I32(0);
    w.I32(-1);
    for (int k = 0; k < 4; ++k)
        w.F64(0.0);
    for (double x : {0.0, 1.0, 1.0, 0.5, 0.0})
        w.F64(x);
    for (double y : {0.0, 0.0, 1.0, 1.5, 1.0})
        w.F64(y);
    for (int v : {0, 1, 1, 2, 2, 4, 4, 0, 2, 3, 3, 4})
        w.I32(v);
    for (int v : {0, 0, 0, 0, 1, 1})
        w.I32(v);
    for (int v : {-1, -1, 1, -1, -1, -1})
        w.I32(v);
    const std::string path = write_bytes(w.Bytes());
    const mt::Mesh m = meshioplusplus::read_tecplot(path);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    const auto cb = m.Cells(0);
    EXPECT_EQ(cb.Type(), "polygon");
    ASSERT_EQ(cb.NumCells(), 2u);
    EXPECT_EQ(std::vector<std::int64_t>(cb.Row(0), cb.Row(0) + cb.RowSize(0)),
              (std::vector<std::int64_t>{0, 1, 2, 4}));
    EXPECT_EQ(std::vector<std::int64_t>(cb.Row(1), cb.Row(1) + cb.RowSize(1)),
              (std::vector<std::int64_t>{4, 2, 3}));
    remove_file(path);
}
