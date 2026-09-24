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
// Tests for the FEBio `.xplt` reader, on plot files built chunk by chunk here
// (the chunk ids are FEBio's own, FEBioPlot/FEBioPlotFile.h).

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xplt.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;
using namespace mt;

namespace {

// Serialises chunks in either byte order.
struct Plt {
    bool mBig = false;

    std::string U32(std::uint32_t V) const {
        std::string s(4, '\0');
        for (int i = 0; i < 4; ++i)
            s[static_cast<std::size_t>(mBig ? 3 - i : i)] =
                static_cast<char>((V >> (8 * i)) & 0xff);
        return s;
    }
    std::string F32(float V) const {
        std::uint32_t u;
        std::memcpy(&u, &V, 4);
        return U32(u);
    }
    std::string Chunk(std::uint32_t Id, const std::string& rPayload) const {
        return U32(Id) + U32(static_cast<std::uint32_t>(rPayload.size())) + rPayload;
    }
    std::string Floats(const std::vector<float>& rValues) const {
        std::string s;
        for (float v : rValues)
            s += F32(v);
        return s;
    }
    std::string Name64(const std::string& rName) const {
        std::string s = rName;
        s.resize(64, '\0');
        return s;
    }
    std::string Item(std::uint32_t Type, std::uint32_t Fmt, const std::string& rName) const {
        return Chunk(0x01020001, Chunk(0x01020002, U32(Type)) + Chunk(0x01020003, U32(Fmt)) +
                                     Chunk(0x01020004, Name64(rName)));
    }
    std::string Var(std::uint32_t Id, std::uint32_t Region,
                    const std::vector<float>& rValues) const {
        return Chunk(0x02020001, Chunk(0x02020002, U32(Id)) +
                                     Chunk(0x02020003, Chunk(Region, Floats(rValues))));
    }
};

// Two tets sharing a face; one variable of every storage format.
std::vector<std::string> plot_chunks(const Plt& rP, std::uint32_t Version, bool Compressed) {
    const std::string header = rP.Chunk(0x01010000, rP.Chunk(0x01010001, rP.U32(Version)) +
                                                        rP.Chunk(0x01010004, rP.U32(Compressed)));
    const std::string dict =
        rP.Chunk(0x01020000, rP.Chunk(0x01021000, rP.Item(0, 3, "energy")) +
                                 rP.Chunk(0x01023000, rP.Item(1, 0, "displacement")) +
                                 rP.Chunk(0x01024000, rP.Item(0, 1, "p") + rP.Item(0, 2, "m") +
                                                          rP.Item(0, 3, "r") + rP.Item(0, 0, "n") +
                                                          rP.Item(2, 1, "stress")) +
                                 rP.Chunk(0x01025000, rP.Item(0, 1, "traction")));
    std::string coords;
    const float xyz[5][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};
    for (std::uint32_t i = 0; i < 5; ++i)
        coords += rP.U32(10 * (i + 1)) + rP.F32(xyz[i][0]) + rP.F32(xyz[i][1]) + rP.F32(xyz[i][2]);
    const std::string nodes = rP.Chunk(
        0x01041000,
        rP.Chunk(0x01041100, rP.Chunk(0x01041101, rP.U32(5)) + rP.Chunk(0x01041102, rP.U32(3))) +
            rP.Chunk(0x01041200, coords));
    const std::string name = "tets";
    const std::string domain = rP.Chunk(
        0x01042000,
        rP.Chunk(
            0x01042100,
            rP.Chunk(0x01042101, rP.Chunk(0x01042102, rP.U32(2)) + rP.Chunk(0x01042103, rP.U32(4)) +
                                     rP.Chunk(0x01032104, rP.U32(2)) +
                                     rP.Chunk(0x01032105, rP.U32(4) + name)) +
                rP.Chunk(0x01042200, rP.Chunk(0x01042201, rP.U32(7) + rP.U32(0) + rP.U32(1) +
                                                              rP.U32(2) + rP.U32(3)) +
                                         rP.Chunk(0x01042201, rP.U32(8) + rP.U32(1) + rP.U32(3) +
                                                                  rP.U32(2) + rP.U32(4)))));
    const std::string nodeset = rP.Chunk(
        0x01044000,
        rP.Chunk(0x01044100, rP.Chunk(0x01044101, rP.Chunk(0x01044103, rP.U32(4) + "base")) +
                                 rP.Chunk(0x01044200, rP.U32(0) + rP.U32(1) + rP.U32(2))));
    std::vector<std::string> out;
    out.push_back(rP.Chunk(0x01000000, header + dict));
    out.push_back(rP.Chunk(0x01040000, nodes + domain + nodeset));
    for (int s = 0; s < 2; ++s) {
        const float t = 0.5f * static_cast<float>(s + 1);
        std::vector<float> disp(15);
        for (int i = 0; i < 15; ++i)
            disp[static_cast<std::size_t>(i)] = t * static_cast<float>(i);
        const std::string data =
            rP.Chunk(0x02020000,
                     rP.Chunk(0x02020100, rP.Var(1, 0, {10 * t})) +
                         rP.Chunk(0x02020300, rP.Var(1, 0, disp)) +
                         rP.Chunk(0x02020400,
                                  rP.Var(1, 1, {1, 2}) + rP.Var(2, 1, {1, 1, 1, 1, 3, 3, 3, 3}) +
                                      rP.Var(3, 1, {7}) +
                                      // FMT_NODE: domain-local nodes in first-seen order 0,1,2,3,4
                                      rP.Var(4, 1, {0, 1, 2, 3, 4}) +
                                      rP.Var(5, 1, {1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6})) +
                         rP.Chunk(0x02020500, rP.Var(1, 1, {9})));
        out.push_back(
            rP.Chunk(0x02000000, rP.Chunk(0x02010000, rP.Chunk(0x02010002, rP.F32(t)) +
                                                          rP.Chunk(0x02010003, rP.U32(0))) +
                                     data));
    }
    return out;
}

std::string write_plot(const std::vector<std::string>& rChunks, const Plt& rP, bool Compressed,
                       std::size_t Drop = 0) {
    std::string bytes = rP.U32(0x00464542);
    for (std::size_t k = 0; k < rChunks.size(); ++k) {
        if (Compressed && k >= 2) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
            uLongf n = compressBound(static_cast<uLong>(rChunks[k].size()));
            std::string z(n, '\0');
            compress2(reinterpret_cast<Bytef*>(z.data()), &n,
                      reinterpret_cast<const Bytef*>(rChunks[k].data()),
                      static_cast<uLong>(rChunks[k].size()), Z_DEFAULT_COMPRESSION);
            z.resize(n);
            bytes += z;
#endif
        } else {
            bytes += rChunks[k];
        }
    }
    bytes.resize(bytes.size() - Drop);
    const std::string path = temp_path(".xplt");
    std::ofstream(path, std::ios::binary) << bytes;
    return path;
}

void expect_state(const Mesh& rMesh, int Step) {
    const double t = 0.5 * (Step + 1);
    EXPECT_DOUBLE_EQ(detail::read_double(rMesh.FieldData("meshio:time"), 0), t);
    EXPECT_EQ(detail::read_int(rMesh.FieldData("xplt:step"), 0), Step);
    EXPECT_DOUBLE_EQ(detail::read_double(rMesh.FieldData("energy"), 0), 10 * t);
    const NDArray& disp = rMesh.PointData("displacement");
    EXPECT_DOUBLE_EQ(detail::read_double(disp, 14), t * 14);
    const NDArray& p = rMesh.CellData("p", 0);
    EXPECT_EQ(detail::read_double(p, 1), 2.0);
    EXPECT_EQ(detail::read_double(rMesh.CellData("r", 0), 1), 7.0);
    // FMT_MULT averaged over the elements sharing each node: node 0 only in
    // element 7 (1), node 4 only in element 8 (3), node 1 in both (2).
    const NDArray& m = rMesh.PointData("m");
    EXPECT_EQ(detail::read_double(m, 0), 1.0);
    EXPECT_EQ(detail::read_double(m, 1), 2.0);
    EXPECT_EQ(detail::read_double(m, 4), 3.0);
    EXPECT_EQ(detail::read_double(rMesh.PointData("n"), 3), 3.0);
    const NDArray& stress = rMesh.CellData("stress", 0);
    EXPECT_EQ(stress.Shape().size(), 2u);
    EXPECT_EQ(detail::read_double(stress, 11), 6.0);
    EXPECT_FALSE(rMesh.HasCellData("traction"));
}

}  // namespace

TEST(Xplt, ReadsEveryStorageFormat) {
    const Plt p;
    const std::string path = write_plot(plot_chunks(p, 0x35, false), p, false);
    EXPECT_EQ(sniff_format(path), "xplt");
    const Mesh first = read_xplt(path);
    ASSERT_EQ(first.NumCellBlocks(), 1u);
    EXPECT_EQ(first.Cells(0).Type(), "tetra");
    EXPECT_EQ(first.NumPoints(), 5u);
    expect_state(first, 0);
    ReadOptions last;
    last.mTimeStep = -1;
    expect_state(read_xplt(path, last), 1);
    const MeshMetadata meta = read_xplt_metadata(path);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.5, 1.0}));
    bool tets = false, base = false;
    for (std::size_t r = 0; r < first.NumRegions(); ++r) {
        const auto& reg = first.Region(r);
        tets = tets || (reg.mName == "tets" && reg.mTag == 4 && reg.NumEntries() == 2);
        base = base ||
               (reg.mName == "base" && reg.mKind == RegionKind::Point && reg.NumEntries() == 3);
    }
    EXPECT_TRUE(tets);
    EXPECT_TRUE(base);
    ReadOptions narrow;
    narrow.mDataArrays = std::vector<std::string>{"p"};
    const Mesh only_p = read_xplt(path, narrow);
    EXPECT_TRUE(only_p.HasCellData("p"));
    EXPECT_FALSE(only_p.HasPointData("displacement"));
    ReadOptions bad;
    bad.mTimeStep = 2;
    EXPECT_THROW(read_xplt(path, bad), ReadError);
}

TEST(Xplt, ReadsBigEndianFiles) {
    Plt p;
    p.mBig = true;
    const std::string path = write_plot(plot_chunks(p, 0x35, false), p, false);
    EXPECT_EQ(sniff_format(path), "xplt");
    expect_state(read_xplt(path), 0);
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
TEST(Xplt, ReadsCompressedStatesAndDropsATruncatedLastOne) {
    const Plt p;
    const auto chunks = plot_chunks(p, 0x35, true);
    const std::string path = write_plot(chunks, p, true);
    expect_state(read_xplt(path), 0);
    ReadOptions last;
    last.mTimeStep = -1;
    expect_state(read_xplt(path, last), 1);
    // A run killed mid-write: the last stream is cut short.
    const std::string cut = write_plot(chunks, p, true, 10);
    EXPECT_EQ(read_xplt_metadata(cut).mTimeValues, (std::vector<double>{0.5}));
}
#endif

TEST(Xplt, RefusesWhatItCannotRead) {
    const Plt p;
    EXPECT_THROW(read_xplt(write_plot(plot_chunks(p, 0x08, false), p, false)), ReadError);
    const std::string junk = temp_path(".xplt");
    std::ofstream(junk, std::ios::binary) << "not a plot file";
    EXPECT_THROW(read_xplt(junk), ReadError);
}

// The data surface "lid" (the first tet's face 0 1 2) becomes a triangle block
// carrying the surface variable; a second mesh chunk, as FEBio writes after
// remeshing, is the mesh of every state after it.
TEST(Xplt, ReadsSurfaceBlocksAndRemeshedStates) {
    const Plt p;
    auto chunks = plot_chunks(p, 0x35, false);
    const std::string surface =
        p.Chunk(0x01043000,
                p.Chunk(0x01043100,
                        p.Chunk(0x01043101, p.Chunk(0x01043104, p.U32(3) + "lid")) +
                            p.Chunk(0x01043200, p.Chunk(0x01043201, p.U32(1) + p.U32(3) + p.U32(0) +
                                                                        p.U32(1) + p.U32(2)))));
    chunks[1] = p.Chunk(0x01040000, chunks[1].substr(8) + surface);
    // The remeshed copy: same topology, every node moved by +10 in x.
    std::string moved = chunks[1];
    const std::string coords_tag = p.U32(0x01041200);
    const std::size_t at = moved.find(coords_tag) + 8;
    for (std::size_t i = 0; i < 5; ++i) {
        const std::size_t x = at + 16 * i + 4;
        float v;
        std::uint32_t u = 0;
        for (int b = 0; b < 4; ++b)
            u |= static_cast<std::uint32_t>(static_cast<unsigned char>(moved[x + b])) << (8 * b);
        std::memcpy(&v, &u, 4);
        moved.replace(x, 4, p.F32(v + 10));
    }
    chunks.insert(chunks.begin() + 3, moved);
    const std::string path = write_plot(chunks, p, false);

    const Mesh first = read_xplt(path);
    ASSERT_EQ(first.NumCellBlocks(), 2u);
    EXPECT_EQ(first.Cells(1).Type(), "triangle");
    EXPECT_EQ(first.Cells(1).NumCells(), 1u);
    EXPECT_EQ(detail::read_double(first.CellData("traction", 1), 0), 9.0);
    EXPECT_TRUE(std::isnan(detail::read_double(first.CellData("traction", 0), 0)));
    EXPECT_EQ(detail::read_int(first.CellData("xplt:surface", 1), 0), 1);
    bool lid = false;
    for (std::size_t r = 0; r < first.NumRegions(); ++r) {
        const auto& reg = first.Region(r);
        lid = lid || (reg.mName == "lid" && reg.mKind == RegionKind::Cell && reg.NumEntries() == 1);
    }
    EXPECT_TRUE(lid);
    EXPECT_EQ(detail::read_double(first.Points(), 0), 0.0);

    ReadOptions last;
    last.mTimeStep = -1;
    const Mesh second = read_xplt(path, last);
    EXPECT_EQ(detail::read_int(second.FieldData("xplt:step"), 0), 1);
    EXPECT_DOUBLE_EQ(detail::read_double(second.PointData("displacement"), 14), 14.0);
    EXPECT_EQ(detail::read_double(second.CellData("traction", 1), 0), 9.0);
    EXPECT_EQ(detail::read_double(second.Points(), 0), 10.0);
    EXPECT_EQ(read_xplt_metadata(path).mTimeValues, (std::vector<double>{0.5, 1.0}));
}
