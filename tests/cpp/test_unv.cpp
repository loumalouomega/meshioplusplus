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

// Direct tests of the native UNV reader/writer on hand-written datasets, so that a
// broken C++ path cannot hide behind the Python fallback.

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/registry.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadOptions;
using meshioplusplus::RegionKind;
using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

namespace {

std::string write_text(const std::string& rText, const std::string& rSuffix = ".unv") {
    const std::string p = mt::temp_path(rSuffix);
    std::ofstream(p, std::ios::binary) << rText;
    return p;
}

std::string nodes(const std::vector<std::vector<double>>& rPts) {
    std::string out = "    -1\n  2411\n";
    char buf[128];
    for (std::size_t k = 0; k < rPts.size(); ++k) {
        std::snprintf(buf, sizeof(buf), "%10zu%10d%10d%10d\n", k + 1, 1, 1, 11);
        out += buf;
        for (double x : rPts[k]) {
            std::snprintf(buf, sizeof(buf), "%25.16E", x);
            std::string s = buf;
            for (char& c : s)
                if (c == 'E')
                    c = 'D';
            out += s;
        }
        out += "\n";
    }
    return out + "    -1\n";
}

const meshioplusplus::Region* find_region(const Mesh& rMesh, const std::string& rName,
                                          RegionKind Kind) {
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mName == rName && rMesh.Region(i).mKind == Kind)
            return &rMesh.Region(i);
    return nullptr;
}

std::vector<std::int64_t> entries(const meshioplusplus::Region& rR) {
    return std::vector<std::int64_t>(rR.Entries(), rR.Entries() + rR.NumEntries());
}

ReadOptions step(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

// One hexahedron20 in the UNV sandwich order: bottom ring, vertical mid-edges, top ring.
const char* kHex20 =
    "    -1\n  2412\n"
    "         1       116         1         1         7        20\n"
    "         1         9         2        10         3        11         4        12\n"
    "        17        18        19        20         5        13         6        14\n"
    "         7        15         8        16\n"
    "    -1\n";

}  // namespace

TEST(Unv, Hex20SandwichOrderPutsMidNodesOnTheirEdges) {
    // Corners 1-8 of the unit cube, mid-nodes 9-20 at the midpoints the sandwich
    // order names: 9-12 bottom, 13-16 top, 17-20 vertical.
    std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int bottom[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    for (const auto& e : bottom)
        pts.push_back({(pts[e[0]][0] + pts[e[1]][0]) / 2, (pts[e[0]][1] + pts[e[1]][1]) / 2, 0});
    for (const auto& e : bottom)
        pts.push_back({(pts[e[0]][0] + pts[e[1]][0]) / 2, (pts[e[0]][1] + pts[e[1]][1]) / 2, 1});
    for (int k = 0; k < 4; ++k)
        pts.push_back({pts[k][0], pts[k][1], 0.5});
    const std::string p = write_text(nodes(pts) + kHex20);
    const Mesh m = meshioplusplus::read_unv(p);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    ASSERT_EQ(m.Cells(0).Type(), "hexahedron20");
    const NDArray& conn = m.Cells(0).Conn();
    const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                              {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (int k = 0; k < 12; ++k)
        for (int c = 0; c < 3; ++c) {
            const auto node = [&](int i) { return read_int(conn, static_cast<std::size_t>(i)); };
            const double mid = read_double(m.Points(), node(8 + k) * 3 + c);
            const double avg = 0.5 * (read_double(m.Points(), node(edges[k][0]) * 3 + c) +
                                      read_double(m.Points(), node(edges[k][1]) * 3 + c));
            EXPECT_DOUBLE_EQ(mid, avg) << "edge " << k;
        }
    // And the writer puts them back in the same order.
    const std::string q = mt::temp_path(".unv");
    meshioplusplus::write_unv(q, m);
    std::ifstream in(q);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("        17        18        19        20         5        13"),
              std::string::npos);
    std::filesystem::remove(p);
    std::filesystem::remove(q);
}

TEST(Unv, GroupsBecomeRegionsInBothLayouts) {
    const std::string text =
        nodes({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}) +
        "    -1\n  2412\n"
        "         1        91         4         9         7         3\n"
        "         1         2         3\n"
        "    -1\n"
        "    -1\n  2467\n"
        "        10         0         0         0         0         0         0         3\n"
        "mixed\n"
        "         7         1         0         0         7         3         0         0\n"
        "         8         1         0         0\n"
        "    -1\n"
        "    -1\n  2429\n"
        "        11         0         0         0         0         0         0         1\n"
        "pairs\n"
        "         7         2\n"
        "    -1\n";
    const std::string p = write_text(text);
    meshioplusplus::UnvInfo info;
    const Mesh m = meshioplusplus::read_unv(p, info);
    const auto* pts = find_region(m, "mixed", RegionKind::Point);
    const auto* cells = find_region(m, "mixed", RegionKind::Cell);
    const auto* pairs = find_region(m, "pairs", RegionKind::Point);
    ASSERT_TRUE(pts && cells && pairs);
    EXPECT_EQ(entries(*pts), (std::vector<std::int64_t>{0, 2}));
    EXPECT_EQ(entries(*cells), (std::vector<std::int64_t>{0}));
    EXPECT_EQ(cells->mTag, 10);
    EXPECT_EQ(cells->mDim, 2);
    EXPECT_EQ(entries(*pairs), (std::vector<std::int64_t>{1}));
    EXPECT_EQ(info.mPointSets["mixed"], (std::vector<std::int64_t>{0, 2}));
    EXPECT_EQ(read_int(m.CellData("unv:pid", 0), 0), 4);
    EXPECT_EQ(read_int(m.CellData("unv:mid", 0), 0), 9);

    // Written back: one group holding nodes (7) and the element (8), same number.
    const std::string q = mt::temp_path(".unv");
    meshioplusplus::write_unv(q, m);
    const Mesh back = meshioplusplus::read_unv(q);
    const auto* back_cells = find_region(back, "mixed", RegionKind::Cell);
    ASSERT_TRUE(back_cells);
    EXPECT_EQ(back_cells->mTag, 10);
    EXPECT_EQ(entries(*find_region(back, "mixed", RegionKind::Point)),
              (std::vector<std::int64_t>{0, 2}));
    EXPECT_EQ(read_int(back.CellData("unv:mid", 0), 0), 9);
    std::filesystem::remove(p);
    std::filesystem::remove(q);
}

TEST(Unv, ResultsAreStepsWithTensorsReordered) {
    std::string text = nodes({{0, 0, 0}, {1, 0, 0}});
    auto dataset = [](int mode, double freq, double v) {
        char buf[1024];
        std::snprintf(buf, sizeof(buf),
                      "    -1\n  2414\n%10d\nstress\n%10d\nNONE\nNONE\nNONE\nNONE\nNONE\n"
                      "%10d%10d%10d%10d%10d%10d\n"
                      "%10d%10d%10d%10d%10d%10d%10d%10d\n%10d%10d\n"
                      "%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n"
                      "%10d\n%13.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n"
                      "    -1\n",
                      mode, 1, 1, 2, 4, 2, 2, 6, 0, 0, 0, 0, 1, mode, 0, 0, 0, 0, 0.0, freq, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                      // Sxx Sxy Syy Sxz Syz Szz
                      v * 1, v * 4, v * 2, v * 6, v * 5, v * 3);
        return std::string(buf);
    };
    text += dataset(1, 10.0, 1.0) + dataset(2, 20.0, 10.0);
    const std::string p = write_text(text);
    const Mesh second = meshioplusplus::read_unv(p, step(1));
    const NDArray& s = second.PointData("stress");
    for (int k = 0; k < 6; ++k)
        EXPECT_DOUBLE_EQ(read_double(s, k), 10.0 * (k + 1));  // xx yy zz xy yz zx
    for (int k = 6; k < 12; ++k)
        EXPECT_TRUE(std::isnan(read_double(s, k)));  // node 2 has no value
    EXPECT_DOUBLE_EQ(read_double(second.FieldData("meshio:time"), 0), 20.0);
    EXPECT_EQ(read_int(second.FieldData("unv:analysis"), 0), 2);
    EXPECT_EQ(read_int(second.FieldData("unv:step"), 0), 2);
    EXPECT_THROW(meshioplusplus::read_unv(p, step(2)), meshioplusplus::ReadError);

    const auto meta = meshioplusplus::read_unv_metadata(p);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{10.0, 20.0}));
    // The generic registry paths the flat bindings use honour the step too.
    const auto& ex = meshioplusplus::registry_readers_ex();
    ASSERT_TRUE(ex.count("unv"));
    const Mesh via_registry = ex.at("unv")(p, step(-1));
    EXPECT_DOUBLE_EQ(read_double(via_registry.FieldData("meshio:time"), 0), 20.0);
    EXPECT_EQ(meshioplusplus::registry_metadata_readers().at("unv")(p, {}).mTimeValues.size(), 2u);
    std::filesystem::remove(p);
}

TEST(Unv, Legacy55HeaderOfOlderMeshioplusplusStillReads) {
    std::string text = nodes({{0, 0, 0}, {1, 0, 0}});
    text +=
        "    -1\n    55\ntemp\ntemp\ntemp\ntemp\ntemp\n"
        "         1         0         1         0         4         1\n"
        "         0         0         0         0         0         0         0         0\n"
        "         0         0\n"
        "  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00\n"
        "  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00\n"
        "         1\n  1.50000E+00\n         2\n  2.50000E+00\n    -1\n";
    const std::string p = write_text(text);
    const Mesh m = meshioplusplus::read_unv(p);
    EXPECT_DOUBLE_EQ(read_double(m.PointData("temp"), 0), 1.5);
    EXPECT_DOUBLE_EQ(read_double(m.PointData("temp"), 1), 2.5);
    std::filesystem::remove(p);
}

TEST(Unv, Dataset55RecordsSevenAndEightAreCounted) {
    // A pyuff-style normal mode: record 7 is "NINT NRV ints", record 8 the reals.
    std::string text = nodes({{0, 0, 0}, {1, 0, 0}});
    text +=
        "    -1\n    55\nNONE\nNONE\nNONE\nNONE\nNONE\n"
        "         1         2         2         8         2         3\n"
        "         2         4         1         3\n"
        "  1.05000e+01  1.00000e+00  1.00000e-02  0.00000e+00\n"
        "         1\n  1.00000e+00  2.00000e+00  3.00000e+00\n"
        "         2\n  4.00000e+00  5.00000e+00  6.00000e+00\n    -1\n";
    const std::string p = write_text(text);
    const Mesh m = meshioplusplus::read_unv(p);
    const NDArray& d = m.PointData("displacement");  // NONE -> the result type's name
    EXPECT_DOUBLE_EQ(read_double(d, 0), 1.0);
    EXPECT_DOUBLE_EQ(read_double(d, 5), 6.0);
    EXPECT_DOUBLE_EQ(read_double(m.FieldData("meshio:time"), 0), 10.5);
    EXPECT_EQ(read_int(m.FieldData("unv:step"), 0), 3);
    std::filesystem::remove(p);
}

namespace {

std::string function58(bool Binary, int Node, int Dir) {
    // Complex double FRF, three evenly spaced frequencies 0, 5, 10.
    const double re[3] = {1.0 * Node, 2.0 * Node, 3.0 * Node};
    const double im[3] = {-1.0 * Dir, -2.0 * Dir, -3.0 * Dir};
    char rec6[128];
    std::snprintf(rec6, sizeof(rec6), "%5d%10d%5d%10d %-10s%10d%4d %-10s%10d%4d\n", 4, 0, 0, 0,
                  "NONE", Node, Dir, "NONE", 1, -3);
    std::string head = "FRF\nNONE\nNONE\nNONE\nNONE\n" + std::string(rec6) +
                       "         6         3         1  0.00000E+00  5.00000E+00  0.00000E+00\n"
                       "        18    0    0    0 NONE                 NONE                \n"
                       "        12    0    0    0 NONE                 NONE                \n"
                       "        13    0    0    0 NONE                 NONE                \n"
                       "         0    0    0    0 NONE                 NONE                \n";
    if (!Binary) {
        std::string data;
        char buf[64];
        for (int k = 0; k < 3; ++k) {
            std::snprintf(buf, sizeof(buf), "%20.12E%20.12E", re[k], im[k]);
            data += buf;
            if (k % 2 == 1 || k == 2)
                data += "\n";
        }
        return "    -1\n    58\n" + head + data + "    -1\n";
    }
    std::string blob(48, '\0');
    for (int k = 0; k < 3; ++k) {
        std::memcpy(&blob[16 * k], &re[k], 8);
        std::memcpy(&blob[16 * k + 8], &im[k], 8);
    }
    // pyuff declares half of the real size for complex data: the reader must not
    // trust the header's byte count.
    char id[128];
    std::snprintf(id, sizeof(id), "%6d%c%6d%6d%12d%12d%6d%6d%12d%12d\n", 58, 'b', 1, 2, 11, 24, 0,
                  0, 0, 0);
    return "    -1\n" + std::string(id) + head + blob + "    -1\n";
}

}  // namespace

TEST(Unv, Functions58AreAFrequencySequenceAsciiAndBinary) {
    for (bool binary : {false, true}) {
        std::string text = nodes({{0, 0, 0}, {1, 0, 0}});
        for (int dir = 1; dir <= 3; ++dir)
            text += function58(binary, 2, dir);
        const std::string p = write_text(text, ".uff");
        EXPECT_EQ(meshioplusplus::read_unv_metadata(p).mTimeValues,
                  (std::vector<double>{0.0, 5.0, 10.0}));
        const Mesh m = meshioplusplus::read_unv(p, step(2));
        EXPECT_EQ(read_int(m.FieldData("unv:analysis"), 0), 5);
        const NDArray& re = m.PointData("frf_real");
        const NDArray& im = m.PointData("frf_imag");
        EXPECT_TRUE(std::isnan(read_double(re, 0)));  // node 1: no function
        for (int d = 0; d < 3; ++d) {
            // reference direction -3: every value is negated
            EXPECT_DOUBLE_EQ(read_double(re, 3 + d), -6.0) << "binary=" << binary;
            EXPECT_DOUBLE_EQ(read_double(im, 3 + d), 3.0 * (d + 1)) << "binary=" << binary;
        }
        std::filesystem::remove(p);
    }
}

TEST(Unv, Functions58WithoutNodesAreAnError) {
    const std::string p = write_text(function58(false, 2, 1), ".uff");
    EXPECT_THROW(meshioplusplus::read_unv(p), meshioplusplus::ReadError);
    std::filesystem::remove(p);
}

TEST(Unv, CartesianCoordinateSystemsAreApplied) {
    const std::string text =
        "    -1\n  2420\n         1\npart\n"
        "         5         0         8\nlocal\n"
        "   0.0000000000000000E+00  -1.0000000000000000E+00   0.0000000000000000E+00\n"
        "   1.0000000000000000E+00   0.0000000000000000E+00   0.0000000000000000E+00\n"
        "   0.0000000000000000E+00   0.0000000000000000E+00   1.0000000000000000E+00\n"
        "   1.0000000000000000E+01   0.0000000000000000E+00   0.0000000000000000E+00\n"
        "    -1\n"
        "    -1\n  2411\n"
        "         1         5         1        11\n"
        "   1.0000000000000000D+00   0.0000000000000000D+00   0.0000000000000000D+00\n"
        "    -1\n";
    const std::string p = write_text(text);
    const Mesh m = meshioplusplus::read_unv(p);
    EXPECT_DOUBLE_EQ(read_double(m.Points(), 0), 10.0);
    EXPECT_DOUBLE_EQ(read_double(m.Points(), 1), 1.0);
    EXPECT_DOUBLE_EQ(read_double(m.Points(), 2), 0.0);
    std::filesystem::remove(p);
}

TEST(Unv, WrongNodeCountForADescriptorIsSkipped) {
    // gmsh writes its quad9 under the linear quad id 94.
    const std::string text =
        nodes({{0, 0, 0},
               {1, 0, 0},
               {1, 1, 0},
               {0, 1, 0},
               {0.5, 0, 0},
               {1, 0.5, 0},
               {0.5, 1, 0},
               {0, 0.5, 0},
               {0.5, 0.5, 0}}) +
        "    -1\n  2412\n"
        "         1        94         1         1         7         9\n"
        "         1         2         3         4         5         6         7         8\n"
        "         9\n"
        "    -1\n";
    const std::string p = write_text(text);
    EXPECT_EQ(meshioplusplus::read_unv(p).NumCellBlocks(), 0u);
    std::filesystem::remove(p);
}

TEST(Unv, SideRegionsAreDroppedAndOthersBecomeGroups) {
    Mesh m = mt::tet_mesh();
    NDArray pts(DType::Int64, {std::size_t{2}});
    pts.As<std::int64_t>()[0] = 0;
    pts.As<std::int64_t>()[1] = 3;
    m.AddRegion(meshioplusplus::Region("clamped", RegionKind::Point, pts));
    NDArray side(DType::Int64, {std::size_t{1}, std::size_t{2}});
    side.As<std::int64_t>()[0] = 0;
    side.As<std::int64_t>()[1] = 1;
    m.AddRegion(meshioplusplus::Region("wall", RegionKind::Side, 2, -1, side));
    const std::string p = mt::temp_path(".unv");
    // The registry writer has no side channel: regions are its only source of groups.
    meshioplusplus::registry_writers().at("unv")(p, m);
    const Mesh back = meshioplusplus::registry_readers().at("unv")(p);
    ASSERT_TRUE(find_region(back, "clamped", RegionKind::Point));
    EXPECT_EQ(entries(*find_region(back, "clamped", RegionKind::Point)),
              (std::vector<std::int64_t>{0, 3}));
    EXPECT_FALSE(find_region(back, "wall", RegionKind::Side));
    std::filesystem::remove(p);
}
