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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/gmsh.hpp"

namespace {
void rt22(const mt::Mesh& mesh, bool binary) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_gmsh22(p, m, binary); },
                  [](const std::string& p) { return meshioplusplus::read_gmsh(p); }, mesh, ".msh");
}
void rt41(const mt::Mesh& mesh, bool binary) {
    mt::roundtrip([=](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_gmsh41(p, m, binary); },
                  [](const std::string& p) { return meshioplusplus::read_gmsh(p); }, mesh, ".msh");
}

// The 4.1 fixture below is written by hand rather than read from
// tests/python/meshes/: this suite has no test-data path, and more importantly
// a format's own writer is not a sufficient oracle for its reader -- the point
// here is to exercise bytes gmsh produces, including the parts meshio++ never
// emits itself.
//
// Geometry: the unit square as two triangles, with its four boundary curves.
// Physical groups: surface 1 -> tag 7 ("plate"), curve 1 -> tag 8 ("bottom").
// Curve 2 carries elements but NO physical tag, which is what makes the
// untagged-block rule observable.
constexpr const char* kEntitiesAscii = R"($MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
1 8 "bottom"
2 7 "plate"
$EndPhysicalNames
$Entities
4 2 1 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 8 2 1 -2
2 1 0 0 1 1 0 0 2 2 -3
1 0 0 0 1 1 0 1 7 2 1 2
$EndEntities
$Nodes
3 4 1 4
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
2 1 0 2
3
4
1 1 0
0 1 0
$EndNodes
$Elements
3 4 1 5
1 1 1 1
1 1 2
1 2 1 1
2 2 3
2 1 2 2
4 1 2 3
5 1 3 4
$EndElements
)";

/// Write @p rText verbatim to a fresh temp file and return its path.
std::string gmsh_write_fixture(const std::string& rText) {
    const std::string path = mt::temp_path(".msh");
    std::ofstream os(path, std::ios::binary);
    os.write(rText.data(), static_cast<std::streamsize>(rText.size()));
    os.close();
    return path;
}

/// The same fixture in binary, so the `data_size`-width and raw-record paths
/// are covered too. @p rDataSize is 4 or 8 -- real files use both.
std::string gmsh_write_fixture_binary(int rDataSize) {
    std::string b;
    auto put_bytes = [&](const void* p, std::size_t n) {
        b.append(static_cast<const char*>(p), n);
    };
    auto put_sz = [&](std::uint64_t v) { put_bytes(&v, static_cast<std::size_t>(rDataSize)); };
    auto put_i = [&](std::int32_t v) { put_bytes(&v, 4); };
    auto put_d = [&](double v) { put_bytes(&v, 8); };

    b += "$MeshFormat\n4.1 1 " + std::to_string(rDataSize) + "\n";
    put_i(1);
    b += "\n$EndMeshFormat\n";
    b += "$PhysicalNames\n2\n1 8 \"bottom\"\n2 7 \"plate\"\n$EndPhysicalNames\n";

    b += "$Entities\n";
    put_sz(4);
    put_sz(2);
    put_sz(1);
    put_sz(0);
    const double pts[4][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    for (int i = 0; i < 4; ++i) {
        put_i(i + 1);
        for (int c = 0; c < 3; ++c)
            put_d(pts[i][c]);
        put_sz(0);
    }
    // curve 1: physical 8, bounded by points 1 and -2
    put_i(1);
    for (double v : {0.0, 0.0, 0.0, 1.0, 0.0, 0.0})
        put_d(v);
    put_sz(1);
    put_i(8);
    put_sz(2);
    put_i(1);
    put_i(-2);
    // curve 2: no physical tag
    put_i(2);
    for (double v : {1.0, 0.0, 0.0, 1.0, 1.0, 0.0})
        put_d(v);
    put_sz(0);
    put_sz(2);
    put_i(2);
    put_i(-3);
    // surface 1: physical 7
    put_i(1);
    for (double v : {0.0, 0.0, 0.0, 1.0, 1.0, 0.0})
        put_d(v);
    put_sz(1);
    put_i(7);
    put_sz(2);
    put_i(1);
    put_i(2);
    b += "\n$EndEntities\n";

    b += "$Nodes\n";
    put_sz(3);
    put_sz(4);
    put_sz(1);
    put_sz(4);
    const int blk_dim[3] = {0, 0, 2};
    const int blk_tag[3] = {1, 2, 1};
    const std::vector<std::vector<int>> blk_nodes = {{1}, {2}, {3, 4}};
    for (int k = 0; k < 3; ++k) {
        put_i(blk_dim[k]);
        put_i(blk_tag[k]);
        put_i(0);
        put_sz(blk_nodes[static_cast<std::size_t>(k)].size());
        for (int nid : blk_nodes[static_cast<std::size_t>(k)])
            put_sz(static_cast<std::uint64_t>(nid));
        for (int nid : blk_nodes[static_cast<std::size_t>(k)])
            for (int c = 0; c < 3; ++c)
                put_d(pts[nid - 1][c]);
    }
    b += "\n$EndNodes\n";

    b += "$Elements\n";
    put_sz(3);
    put_sz(4);
    put_sz(1);
    put_sz(5);
    // (entity dim, entity tag, gmsh type, count) then (tag, nodes...) rows
    const int hdr[3][4] = {{1, 1, 1, 1}, {1, 2, 1, 1}, {2, 1, 2, 2}};
    const std::vector<std::vector<int>> rows = {{1, 1, 2}, {2, 2, 3}, {4, 1, 2, 3}, {5, 1, 3, 4}};
    std::size_t r = 0;
    for (int k = 0; k < 3; ++k) {
        for (int f = 0; f < 3; ++f)
            put_i(hdr[k][f]);
        put_sz(static_cast<std::uint64_t>(hdr[k][3]));
        for (int e = 0; e < hdr[k][3]; ++e, ++r)
            for (int v : rows[r])
                put_sz(static_cast<std::uint64_t>(v));
    }
    b += "\n$EndElements\n";
    return gmsh_write_fixture(b);
}

/// The block-wise `gmsh:physical` values, as a flat per-block list.
std::vector<std::int64_t> physical_per_block(const meshioplusplus::Mesh& rMesh) {
    std::vector<std::int64_t> out;
    for (std::size_t b = 0; b < rMesh.CellDataNumBlocks("gmsh:physical"); ++b)
        out.push_back(meshioplusplus::detail::read_int(rMesh.CellData("gmsh:physical", b), 0));
    return out;
}

void expect_entities_fixture(const meshioplusplus::Mesh& mesh) {
    EXPECT_EQ(mesh.NumPoints(), 4u);
    ASSERT_EQ(mesh.NumCellBlocks(), 3u);
    EXPECT_EQ(mesh.Cells(0).Type(), "line");
    EXPECT_EQ(mesh.Cells(1).Type(), "line");
    EXPECT_EQ(mesh.Cells(2).Type(), "triangle");

    // The heart of it: physical tags come from $Entities, not from the element
    // rows the way 2.2 does. The untagged curve gets 0, and the array still has
    // one entry per block.
    ASSERT_TRUE(mesh.HasCellData("gmsh:physical"));
    ASSERT_EQ(mesh.CellDataNumBlocks("gmsh:physical"), mesh.NumCellBlocks());
    EXPECT_EQ(physical_per_block(mesh), (std::vector<std::int64_t>{8, 0, 7}));

    // ... and with $PhysicalNames that makes named regions carrying the real
    // dimension and tag.
    ASSERT_EQ(mesh.NumRegions(), 2u);
    const auto& bottom = mesh.HasRegion("bottom") ? mesh.Region(0) : mesh.Region(1);
    const auto& plate = mesh.HasRegion("plate") ? mesh.Region(1) : mesh.Region(0);
    EXPECT_EQ(bottom.mName, "bottom");
    EXPECT_EQ(bottom.mDim, 1);
    EXPECT_EQ(bottom.mTag, 8);
    EXPECT_EQ(bottom.mEntries.Size(), 1u);
    EXPECT_EQ(plate.mName, "plate");
    EXPECT_EQ(plate.mDim, 2);
    EXPECT_EQ(plate.mTag, 7);
    EXPECT_EQ(plate.mEntries.Size(), 2u);
}
}  // namespace

TEST(Gmsh, V22Ascii) {
    rt22(mt::tri_mesh(), false);
    rt22(mt::tet_mesh(), false);
    rt22(mt::hex_mesh(), false);
}
TEST(Gmsh, V22Binary) {
    rt22(mt::tri_mesh(), true);
    rt22(mt::tet_mesh(), true);
}
TEST(Gmsh, V41Ascii) {
    rt41(mt::tri_mesh(), false);
    rt41(mt::tet_mesh(), false);
}
TEST(Gmsh, V41Binary) {
    rt41(mt::hex_mesh(), true);
}
TEST(Gmsh, SecondOrder) {
    rt22(mt::tet10_mesh(), false);
    rt22(mt::hex20_mesh(), false);
}

TEST(Gmsh, EntitiesAscii) {
    const std::string path = gmsh_write_fixture(kEntitiesAscii);
    meshioplusplus::GmshInfo info;
    const meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh(path, info);
    expect_entities_fixture(mesh);

    // Bounding entities are signed, so they ride GmshInfo rather than a Region.
    ASSERT_EQ(info.mBoundingEntities.size(), 3u);
    EXPECT_EQ(info.mBoundingEntities[0], (std::vector<std::int32_t>{1, -2}));
    EXPECT_EQ(info.mBoundingEntities[1], (std::vector<std::int32_t>{2, -3}));
    EXPECT_EQ(info.mBoundingEntities[2], (std::vector<std::int32_t>{1, 2}));
    std::remove(path.c_str());
}

TEST(Gmsh, EntitiesBinary) {
    // Both size_t widths: gmsh writes 4 as readily as 8 (example/example.msh
    // is a `4.1 0 4` file), and assuming 8 silently misreads every count.
    for (int data_size : {4, 8}) {
        const std::string path = gmsh_write_fixture_binary(data_size);
        meshioplusplus::GmshInfo info;
        const meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh(path, info);
        expect_entities_fixture(mesh);
        ASSERT_EQ(info.mBoundingEntities.size(), 3u);
        EXPECT_EQ(info.mBoundingEntities[0], (std::vector<std::int32_t>{1, -2}));
        std::remove(path.c_str());
    }
}

TEST(Gmsh, EntitiesSurviveARoundTrip) {
    for (bool binary : {false, true}) {
        const std::string src = gmsh_write_fixture(kEntitiesAscii);
        meshioplusplus::GmshInfo info;
        const meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh(src, info);

        const std::string out = mt::temp_path(".msh");
        meshioplusplus::write_gmsh41(out, mesh, binary, info);
        meshioplusplus::GmshInfo back_info;
        const meshioplusplus::Mesh back = meshioplusplus::read_gmsh(out, back_info);

        // Physical membership survives a 4.1 round-trip -- which before
        // $Entities was written only 2.2 could do.
        expect_entities_fixture(back);
        EXPECT_EQ(back_info.mBoundingEntities, info.mBoundingEntities);
        mt::expect_mesh_eq(mesh, back);
        std::remove(src.c_str());
        std::remove(out.c_str());
    }
}

TEST(Gmsh, EntitiesMetadataAgreesWithARealRead) {
    const std::string path = gmsh_write_fixture(kEntitiesAscii);
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_gmsh_metadata(path);
    const meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh(path);

    EXPECT_EQ(meta.mNumPoints, mesh.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), mesh.NumCellBlocks());
    // A summary that named different arrays or groups than a real read would be
    // worse than no summary at all.
    EXPECT_NE(std::find(meta.mCellDataNames.begin(), meta.mCellDataNames.end(), "gmsh:physical"),
              meta.mCellDataNames.end());
    ASSERT_EQ(meta.mRegions.size(), mesh.NumRegions());
    for (const auto& rs : meta.mRegions) {
        bool found = false;
        for (std::size_t i = 0; i < mesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = mesh.Region(i);
            if (r.mName != rs.mName)
                continue;
            found = true;
            EXPECT_EQ(rs.mDim, r.mDim);
            EXPECT_EQ(rs.mTag, r.mTag);
            EXPECT_EQ(rs.mNumEntries, r.mEntries.Size());
        }
        EXPECT_TRUE(found) << rs.mName;
    }
    std::remove(path.c_str());
}

TEST(Gmsh, NoPhysicalTagsMeansNoPhysicalCellData) {
    // Dropping $PhysicalNames and every physical tag: synthesizing an all-zero
    // gmsh:physical would invent a group the file does not have.
    std::string text = kEntitiesAscii;
    const std::size_t b = text.find("$PhysicalNames");
    const std::size_t e = text.find("$EndPhysicalNames\n");
    text.erase(b, e + std::strlen("$EndPhysicalNames\n") - b);
    text.replace(text.find("1 0 0 0 1 0 0 1 8 2 1 -2"), std::strlen("1 0 0 0 1 0 0 1 8 2 1 -2"),
                 "1 0 0 0 1 0 0 0 2 1 -2");
    text.replace(text.find("1 0 0 0 1 1 0 1 7 2 1 2"), std::strlen("1 0 0 0 1 1 0 1 7 2 1 2"),
                 "1 0 0 0 1 1 0 0 2 1 2");

    const std::string path = gmsh_write_fixture(text);
    const meshioplusplus::Mesh mesh = meshioplusplus::read_gmsh(path);
    EXPECT_FALSE(mesh.HasCellData("gmsh:physical"));
    EXPECT_TRUE(mesh.HasCellData("gmsh:geometrical"));
    EXPECT_EQ(mesh.NumRegions(), 0u);
    std::remove(path.c_str());
}

TEST(Gmsh, WritingWithoutDimTagsEmitsNoEntities) {
    const std::string path = mt::temp_path(".msh");
    meshioplusplus::write_gmsh41(path, mt::tet_mesh(), false);
    std::ifstream is(path);
    const std::string text((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
    EXPECT_EQ(text.find("$Entities"), std::string::npos);
    std::remove(path.c_str());
}
