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
// Tests for the Kratos MDPA (.mdpa) reader/writer.

// System includes
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::RegionKind;

namespace {

/// Write @p rContents to a fresh temp `.mdpa` file and return its path.
std::string mdpa_temp_file(const std::string& rContents) {
    const std::string path = mt::temp_path(".mdpa");
    std::ofstream out(path);
    out << rContents;
    out.close();
    return path;
}

/// Read the whole file at @p rPath.
std::string mdpa_slurp(const std::string& rPath) {
    std::ifstream in(rPath);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

Mesh mdpa_read_string(const std::string& rContents) {
    const std::string path = mdpa_temp_file(rContents);
    Mesh m = meshioplusplus::read_mdpa(path);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Round trips
// ---------------------------------------------------------------------------

TEST(Mdpa, RoundTripSimpleMeshes) {
    const auto writer = [](const std::string& p, const Mesh& m) {
        meshioplusplus::write_mdpa(p, m);
    };
    const auto reader = [](const std::string& p) { return meshioplusplus::read_mdpa(p); };
    mt::roundtrip(writer, reader, mt::line_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::tri_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::quad_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::tet_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::hex_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::tet10_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::triangle6_mesh(), ".mdpa");
    mt::roundtrip(writer, reader, mt::quad8_mesh(), ".mdpa");
}

// hexahedron20 is one of the two types whose Kratos node order differs from
// meshio's; the permutation must be applied on read and undone on write.
TEST(Mdpa, RoundTripHex20AppliesKratosPermutation) {
    const std::string path = mt::temp_path(".mdpa");
    const Mesh in = mt::hex20_mesh();
    meshioplusplus::write_mdpa(path, in);

    // The file must carry the *Kratos* order, i.e. the mid-edge nodes of the
    // top face (meshio slots 16..19) written in slots 12..15.
    const std::string text = mdpa_slurp(path);
    EXPECT_NE(text.find("Begin Elements Element3D20N"), std::string::npos);

    const Mesh out = meshioplusplus::read_mdpa(path);
    mt::expect_mesh_eq(in, out);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    // A round trip that skipped both permutations would also pass the above,
    // so pin the on-disk order itself: node ids of the single cell, in file
    // order, must be the meshio row permuted by the Kratos table.
    const std::vector<int> kratos = {0,  1, 2,  3,  4,  5,  6,  7,  8,  11,
                                     10, 9, 16, 19, 18, 17, 12, 13, 14, 15};
    const auto cb = in.Cells(0);
    std::string expected;
    for (std::size_t j = 0; j < 20; ++j)
        expected += " " + std::to_string(meshioplusplus::detail::read_int(
                                             cb.Conn(), static_cast<std::size_t>(kratos[j])) +
                                         1);
    EXPECT_NE(text.find(expected), std::string::npos) << text;
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// A mesh mixing volume/surface cells with surface/line cells exercises the
// Elements-vs-Conditions split: a block whose default Kratos element name is
// 2-D is written as a Condition, the rest as Elements.
TEST(Mdpa, MixedElementAndConditionBlocks) {
    const std::string path = mt::temp_path(".mdpa");
    const Mesh in = mt::tri_quad_mesh();  // triangle, quad, triangle
    meshioplusplus::write_mdpa(path, in);
    const std::string text = mdpa_slurp(path);
    EXPECT_NE(text.find("Begin Elements Element3D3N"), std::string::npos) << text;
    EXPECT_NE(text.find("Begin Conditions SurfaceCondition3D4N"), std::string::npos) << text;

    const Mesh out = meshioplusplus::read_mdpa(path);
    ASSERT_EQ(out.NumCellBlocks(), 3u);
    EXPECT_EQ(std::string(out.Cells(0).Type()), "triangle");
    EXPECT_EQ(std::string(out.Cells(1).Type()), "quad");
    EXPECT_EQ(std::string(out.Cells(2).Type()), "triangle");
    mt::expect_mesh_eq(in, out);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, PropertyIdsRoundTripAsGmshPhysical) {
    Mesh in = mt::tri_quad_mesh();
    std::vector<NDArray> tags;
    std::int64_t next = 7;
    for (std::size_t b = 0; b < in.NumCellBlocks(); ++b) {
        NDArray a(meshioplusplus::DType::Int64, {in.Cells(b).NumCells()});
        for (std::size_t r = 0; r < in.Cells(b).NumCells(); ++r)
            a.As<std::int64_t>()[r] = next;
        ++next;
        tags.push_back(std::move(a));
    }
    in.AddCellData("gmsh:physical", std::move(tags));

    const std::string path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(path, in);
    const Mesh out = meshioplusplus::read_mdpa(path);
    ASSERT_TRUE(out.HasCellData("gmsh:physical"));
    ASSERT_EQ(out.CellDataNumBlocks("gmsh:physical"), 3u);
    for (std::size_t b = 0; b < 3; ++b)
        EXPECT_EQ(meshioplusplus::detail::read_int(out.CellData("gmsh:physical", b), 0),
                  7 + static_cast<std::int64_t>(b));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
// SubModelParts <-> regions
// ---------------------------------------------------------------------------

TEST(Mdpa, SubModelPartsBecomeRegions) {
    const Mesh m = mdpa_read_string(R"(Begin Nodes
 1 0.0 0.0 0.0
 2 1.0 0.0 0.0
 3 1.0 1.0 0.0
 4 0.0 1.0 0.0
End Nodes

Begin Elements Element3D3N
 1 0 1 2 3
 2 0 1 3 4
End Elements

Begin Conditions LineCondition3D2N
 1 0 1 2
 2 0 2 3
End Conditions

Begin SubModelPart Domain
    Begin SubModelPartNodes
        1
        2
        3
    End SubModelPartNodes
    Begin SubModelPartElements
        2
    End SubModelPartElements
End SubModelPart

Begin SubModelPart Skin
    Begin SubModelPartConditions
        1
        2
    End SubModelPartConditions
End SubModelPart
)");

    ASSERT_EQ(m.NumCellBlocks(), 2u);
    ASSERT_EQ(m.Cells(0).NumCells(), 2u);  // triangles (elements)
    ASSERT_EQ(m.Cells(1).NumCells(), 2u);  // lines (conditions)

    const std::vector<std::string> names = m.RegionNames();
    ASSERT_EQ(names.size(), 2u);
    EXPECT_EQ(names[0], "Domain");
    EXPECT_EQ(names[1], "Skin");

    const std::size_t ip = m.FindRegion("Domain", RegionKind::Point);
    ASSERT_NE(ip, Mesh::npos);
    ASSERT_EQ(m.Region(ip).NumEntries(), 3u);
    EXPECT_EQ(m.Region(ip).Entries()[0], 0);
    EXPECT_EQ(m.Region(ip).Entries()[2], 2);

    const std::size_t ic = m.FindRegion("Domain", RegionKind::Cell);
    ASSERT_NE(ic, Mesh::npos);
    ASSERT_EQ(m.Region(ic).NumEntries(), 1u);
    EXPECT_EQ(m.Region(ic).Entries()[0], 1);  // element 2 -> global cell 1

    const std::size_t is = m.FindRegion("Skin", RegionKind::Cell);
    ASSERT_NE(is, Mesh::npos);
    ASSERT_EQ(m.Region(is).NumEntries(), 2u);
    // Conditions live in block 1, whose global base is 2.
    EXPECT_EQ(m.Region(is).Entries()[0], 2);
    EXPECT_EQ(m.Region(is).Entries()[1], 3);
}

TEST(Mdpa, RegionsRoundTripThroughSubModelParts) {
    Mesh in = mt::tri_quad_mesh();
    {
        NDArray pts(meshioplusplus::DType::Int64, {3});
        pts.As<std::int64_t>()[0] = 0;
        pts.As<std::int64_t>()[1] = 2;
        pts.As<std::int64_t>()[2] = 5;
        in.AddRegion(meshioplusplus::Region("Support", RegionKind::Point, std::move(pts)));
        NDArray cells(meshioplusplus::DType::Int64, {2});
        cells.As<std::int64_t>()[0] = 0;  // triangle block, row 0 -> element
        cells.As<std::int64_t>()[1] = 2;  // quad block, row 0 -> condition
        in.AddRegion(meshioplusplus::Region("Support", RegionKind::Cell, std::move(cells)));
    }
    const std::string path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(path, in);
    const std::string text = mdpa_slurp(path);
    EXPECT_NE(text.find("Begin SubModelPart Support"), std::string::npos) << text;
    EXPECT_NE(text.find("Begin SubModelPartElements"), std::string::npos) << text;
    EXPECT_NE(text.find("Begin SubModelPartConditions"), std::string::npos) << text;

    const Mesh out = meshioplusplus::read_mdpa(path);
    const std::size_t ip = out.FindRegion("Support", RegionKind::Point);
    const std::size_t ic = out.FindRegion("Support", RegionKind::Cell);
    ASSERT_NE(ip, Mesh::npos);
    ASSERT_NE(ic, Mesh::npos);
    ASSERT_EQ(out.Region(ip).NumEntries(), 3u);
    EXPECT_EQ(out.Region(ip).Entries()[1], 2);
    ASSERT_EQ(out.Region(ic).NumEntries(), 2u);
    EXPECT_EQ(out.Region(ic).Entries()[0], 0);
    EXPECT_EQ(out.Region(ic).Entries()[1], 2);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
// Data blocks
// ---------------------------------------------------------------------------

TEST(Mdpa, NodalDataBecomesPointData) {
    const Mesh m = mdpa_read_string(R"(Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
End Nodes

Begin NodalData TEMPERATURE // scalar
    1 25.5
    2 0 30.5
End NodalData

Begin NodalData DISPLACEMENT[3]
    1 0 0.0 0.0 0.5
    3 1 0.1 0.2 0.3
End NodalData
)");
    ASSERT_TRUE(m.HasPointData("TEMPERATURE"));
    const NDArray& t = m.PointData("TEMPERATURE");
    ASSERT_EQ(t.Size(), 3u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(t, 0), 25.5);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(t, 1), 30.5);
    EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(t, 2)));

    ASSERT_TRUE(m.HasPointData("TEMPERATURE_fixed_status"));
    const NDArray& f = m.PointData("TEMPERATURE_fixed_status");
    EXPECT_EQ(meshioplusplus::detail::read_int(f, 0), -1);
    EXPECT_EQ(meshioplusplus::detail::read_int(f, 1), 0);

    ASSERT_TRUE(m.HasPointData("DISPLACEMENT"));
    const NDArray& d = m.PointData("DISPLACEMENT");
    ASSERT_EQ(d.Size(), 9u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(d, 2), 0.5);
    EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(d, 3)));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(d, 8), 0.3);
    EXPECT_EQ(meshioplusplus::detail::read_int(m.PointData("DISPLACEMENT_fixed_status"), 2), 1);
}

TEST(Mdpa, ElementalAndConditionalDataBecomeCellData) {
    const Mesh m = mdpa_read_string(R"(Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
End Nodes

Begin Elements Element3D3N
1 0 1 2 3
End Elements

Begin Conditions LineCondition3D2N
7 0 1 2
End Conditions

Begin ElementalData STRESS
    1 10.5
End ElementalData

Begin ConditionalData PRESSURE
    7 -5.5
End ConditionalData
)");
    ASSERT_EQ(m.NumCellBlocks(), 2u);
    ASSERT_TRUE(m.HasCellData("STRESS"));
    ASSERT_EQ(m.CellDataNumBlocks("STRESS"), 2u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.CellData("STRESS", 0), 0), 10.5);
    EXPECT_TRUE(std::isnan(meshioplusplus::detail::read_double(m.CellData("STRESS", 1), 0)));
    ASSERT_TRUE(m.HasCellData("PRESSURE"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.CellData("PRESSURE", 1), 0), -5.5);
}

TEST(Mdpa, ModelPartDataBecomesFieldData) {
    const Mesh m = mdpa_read_string(R"(Begin ModelPartData
    // a comment
    DOMAIN_SIZE 3
    GRAVITY_Z -9.81
End ModelPartData

Begin Nodes
1 0.0 0.0 0.0
End Nodes
)");
    ASSERT_TRUE(m.HasFieldData("DOMAIN_SIZE"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.FieldData("DOMAIN_SIZE"), 0), 3.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(m.FieldData("GRAVITY_Z"), 0), -9.81);
}

// Application-specific entity names resolve through their Kratos suffix.
TEST(Mdpa, ResolvesApplicationSpecificEntityNames) {
    const Mesh m = mdpa_read_string(R"(Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 0.0 1.0
End Nodes

Begin Elements SmallDisplacementElement3D4N
1 0 1 2 3 4
End Elements
)");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(m.Cells(0).Type()), "tetra");
}

// ---------------------------------------------------------------------------
// Unsupported constructs must throw by name (so the Python shim can fall back)
// ---------------------------------------------------------------------------

TEST(Mdpa, UnsupportedConstructsThrowByName) {
    struct Case {
        std::string contents;
        std::string needle;
    };
    const std::vector<Case> cases = {
        {"Begin Table 1 TIME VALUE\n 0.0 1.0\nEnd Table\n", "Table"},
        {"Begin Geometries Triangle3D3\n1 1 2 3\nEnd Geometries\n", "Geometries"},
        {"Begin Mesh 1\nEnd Mesh\n", "Mesh"},
        {"Begin ModelPartData\n  NAME \"a string\"\nEnd ModelPartData\n", "ModelPartData"},
        {"Begin Nodes\n1 0 0 0\n", "End Nodes"},
        {"Begin SubModelPart S\n  Begin SubModelPartData\n    K 1\n  End SubModelPartData\n"
         "End SubModelPart\n",
         "SubModelPartData"},
        {"Begin Constraints Foo\nEnd Constraints\n", "Begin Constraints Foo"},
    };
    for (const auto& c : cases) {
        const std::string path = mdpa_temp_file(c.contents);
        try {
            meshioplusplus::read_mdpa(path);
            ADD_FAILURE() << "expected ReadError for: " << c.contents;
        } catch (const meshioplusplus::ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(c.needle), std::string::npos)
                << "message '" << e.what() << "' does not name '" << c.needle << "'";
        }
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

TEST(Mdpa, PropertiesBodyNoLongerThrows) {
    // A non-empty `Begin Properties` used to be on the list above, which made
    // essentially every production deck unreadable. It is now parsed; without
    // an MdpaInfo to hold it the body is dropped with a warning, the way the
    // flat bindings drop MedInfo.
    const Mesh m = mdpa_read_string(
        "Begin Properties 1\n    DENSITY 1.0\nEnd Properties\n"
        "Begin Nodes\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Nodes\n"
        "Begin Elements Element2D3N\n1 1 1 2 3\nEnd Elements\n");
    EXPECT_EQ(m.NumPoints(), 3u);
    EXPECT_EQ(m.NumCellBlocks(), 1u);
}

TEST(Mdpa, UnknownEntityNameThrows) {
    const std::string path = mdpa_temp_file(
        "Begin Nodes\n1 0 0 0\n2 1 0 0\nEnd Nodes\n"
        "Begin Elements TotallyUnknownThing\n1 0 1 2\nEnd Elements\n");
    EXPECT_THROW(meshioplusplus::read_mdpa(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
// Arbitrary node ids (gapped, non-monotonic) -- what a real Kratos deck left by
// a SubModelPart extraction or an entity removal looks like.
// ---------------------------------------------------------------------------

namespace {

/// Gapped AND non-monotonic ids, exercised by every id block that resolves one.
/// Kept textually identical to `GAPPED_NODE_DECK` in tests/python/test_mdpa.py,
/// which is what makes the C++/Python parity test there meaningful.
const char* const kMdpaGappedDeck = R"(Begin Nodes
10 0.0 0.0 0.0
7  1.0 0.0 0.0
42 0.0 1.0 0.0
5  0.0 0.0 1.0
End Nodes

Begin Elements Element3D4N
100 1 10 7 42 5
End Elements

Begin NodalData TEMPERATURE
42 3.5
10 1.5
End NodalData

Begin SubModelPart Gapped
    Begin SubModelPartNodes
        42
        5
    End SubModelPartNodes
End SubModelPart
)";

}  // namespace

TEST(Mdpa, GappedNodeIdsRead) {
    const Mesh m = mdpa_read_string(kMdpaGappedDeck);
    ASSERT_EQ(m.NumPoints(), 4u);

    // Points come back in FILE order, never sorted by id.
    const double* p = m.Points().As<double>();
    const std::vector<double> expected = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (std::size_t i = 0; i < expected.size(); ++i)
        EXPECT_DOUBLE_EQ(p[i], expected[i]) << "coordinate " << i;

    // Connectivity resolved through the map: ids 10/7/42/5 are rows 0/1/2/3.
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    const auto blk = m.Cells(0);
    EXPECT_EQ(std::string(blk.Type()), "tetra");
    ASSERT_EQ(blk.NumCells(), 1u);
    const std::int64_t* c = blk.Conn().As<std::int64_t>();
    for (int j = 0; j < 4; ++j)
        EXPECT_EQ(c[j], j);

    // NodalData keyed by the real file ids.
    ASSERT_TRUE(m.HasPointData("TEMPERATURE"));
    const double* t = m.PointData("TEMPERATURE").As<double>();
    EXPECT_DOUBLE_EQ(t[0], 1.5);  // id 10
    EXPECT_TRUE(std::isnan(t[1]));
    EXPECT_DOUBLE_EQ(t[2], 3.5);  // id 42
    EXPECT_TRUE(std::isnan(t[3]));

    // And so is the SubModelPart node list.
    const std::size_t ip = m.FindRegion("Gapped", RegionKind::Point);
    ASSERT_NE(ip, Mesh::npos);
    ASSERT_EQ(m.Region(ip).NumEntries(), 2u);
    EXPECT_EQ(m.Region(ip).Entries()[0], 2);  // id 42
    EXPECT_EQ(m.Region(ip).Entries()[1], 3);  // id 5
}

TEST(Mdpa, GappedIdsAreNotGatedOnLenient) {
    // Accepting arbitrary ids is strictly more correct, not more lenient: a
    // strict read must succeed and record nothing as skipped.
    const std::string path = mdpa_temp_file(kMdpaGappedDeck);
    meshioplusplus::MdpaInfo info;
    meshioplusplus::ReadOptions opts;  // mLenient stays false
    const Mesh m = meshioplusplus::read_mdpa(path, info, opts);
    EXPECT_EQ(m.NumPoints(), 4u);
    EXPECT_TRUE(info.mSkippedConstructs.empty());
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, NonMonotonicNodeIds) {
    const Mesh m = mdpa_read_string(
        "Begin Nodes\n4 0 0 0\n3 1 0 0\n2 0 1 0\n1 0 0 1\nEnd Nodes\n"
        "Begin Elements Element3D4N\n1 0 1 2 3 4\nEnd Elements\n");
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    // ids 1,2,3,4 are rows 3,2,1,0 -- descending, which is exactly the case a
    // "row == id - 1" reader gets silently backwards.
    for (int j = 0; j < 4; ++j)
        EXPECT_EQ(c[j], 3 - j);
}

TEST(Mdpa, DanglingConnectivityNodeIdThrowsWithTheFileId) {
    const std::string path = mdpa_temp_file(
        "Begin Nodes\n10 0 0 0\n7 1 0 0\n42 0 1 0\n5 0 0 1\nEnd Nodes\n"
        "Begin Elements Element3D4N\n1 0 10 7 42 999\nEnd Elements\n");
    try {
        meshioplusplus::read_mdpa(path);
        ADD_FAILURE() << "expected ReadError for a dangling node id";
    } catch (const meshioplusplus::ReadError& e) {
        const std::string what = e.what();
        // The FILE id, not a row index and not `id + 1`.
        EXPECT_NE(what.find("999"), std::string::npos) << what;
        EXPECT_EQ(what.find("1000"), std::string::npos) << what;
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, DuplicateNodeIdThrows) {
    // Unrepresentable, not merely incomplete: two coordinate rows claim one id.
    const std::string path = mdpa_temp_file("Begin Nodes\n5 0 0 0\n5 1 1 1\nEnd Nodes\n");
    try {
        meshioplusplus::read_mdpa(path);
        ADD_FAILURE() << "expected ReadError for a duplicate node id";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("duplicate node id"), std::string::npos) << e.what();
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, BareCoordinateRowsStillRead) {
    // The id-less `x y z` form, which the reader has always accepted and which
    // nothing covered until the id map made it worth pinning.
    const Mesh m = mdpa_read_string(
        "Begin Nodes\n0 0 0\n1 0 0\n0 1 0\nEnd Nodes\n"
        "Begin Elements Element2D3N\n1 0 1 2 3\nEnd Elements\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    const double* p = m.Points().As<double>();
    EXPECT_DOUBLE_EQ(p[3], 1.0);  // second row's x
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    for (int j = 0; j < 3; ++j)
        EXPECT_EQ(c[j], j);
}

TEST(Mdpa, MixedIdAndBareNodeRows) {
    // An id-less row takes its POSITION as its id -- the documented rule, and
    // the only one that keeps a mixed file whose explicit ids happen to be
    // sequential reading as it always did.
    const Mesh m = mdpa_read_string(
        "Begin Nodes\n5 0 0 0\n1 0 0\n0 1 0\nEnd Nodes\n"
        "Begin Elements Element2D3N\n1 0 5 2 3\nEnd Elements\n");
    ASSERT_EQ(m.NumPoints(), 3u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    const std::int64_t* c = m.Cells(0).Conn().As<std::int64_t>();
    for (int j = 0; j < 3; ++j)
        EXPECT_EQ(c[j], j);
}

TEST(Mdpa, SubModelPartNodesWithUnknownIdIsSkipped) {
    // A set is not connectivity: an unknown member is dropped with a warning,
    // matching what the entity lists already did.
    const Mesh m = mdpa_read_string(
        "Begin Nodes\n10 0 0 0\n7 1 0 0\n42 0 1 0\nEnd Nodes\n"
        "Begin SubModelPart S\n  Begin SubModelPartNodes\n    42\n    999\n"
        "  End SubModelPartNodes\nEnd SubModelPart\n");
    const std::size_t ip = m.FindRegion("S", RegionKind::Point);
    ASSERT_NE(ip, Mesh::npos);
    ASSERT_EQ(m.Region(ip).NumEntries(), 1u);
    EXPECT_EQ(m.Region(ip).Entries()[0], 2);
}

// ---------------------------------------------------------------------------
// Original ids preserved on write (kMdpaIdName / "mdpa:id")
// ---------------------------------------------------------------------------

TEST(Mdpa, GappedIdsRoundTripThroughAWrite) {
    // The write half of roadmap #0: a gapped read attaches point_data/
    // cell_data["mdpa:id"], and the writer must honour it -- so a re-write
    // names the ORIGINAL file ids, not row+1/counter renumbering.
    const std::string in_path = mdpa_temp_file(kMdpaGappedDeck);
    const Mesh m = meshioplusplus::read_mdpa(in_path);
    ASSERT_TRUE(m.HasPointData(meshioplusplus::kMdpaIdName));
    ASSERT_TRUE(m.HasCellData(meshioplusplus::kMdpaIdName));

    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m);
    const std::string text = mdpa_slurp(out_path);
    EXPECT_NE(text.find("Begin Nodes\n 10 "), std::string::npos) << text;
    EXPECT_NE(text.find(" 7 "), std::string::npos) << text;
    EXPECT_NE(text.find(" 42 "), std::string::npos) << text;
    EXPECT_NE(text.find(" 5 "), std::string::npos) << text;
    EXPECT_NE(text.find("Begin Elements Element3D4N\n  100 "), std::string::npos) << text;
    EXPECT_EQ(text.find("Begin Elements Element3D4N\n  1 "), std::string::npos) << text;

    // And the round trip is exact: re-reading the written file reproduces the
    // same mesh, ids included.
    const Mesh out = meshioplusplus::read_mdpa(out_path);
    mt::expect_mesh_eq(m, out);
    ASSERT_TRUE(out.HasPointData(meshioplusplus::kMdpaIdName));
    EXPECT_EQ(meshioplusplus::detail::read_int(out.PointData(meshioplusplus::kMdpaIdName), 0), 10);
    std::error_code ec;
    std::filesystem::remove(in_path, ec);
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, SubModelPartNodeAndElementReferencesUseThePreservedIds) {
    // Point regions AND cell regions (SubModelPart) must reference the SAME
    // preserved ids as the Nodes/Elements blocks they point into -- the file
    // must stay internally consistent whichever numbering is actually used.
    const Mesh m = mdpa_read_string(kMdpaGappedDeck);
    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m);
    const std::string text = mdpa_slurp(out_path);
    const std::size_t smp = text.find("Begin SubModelPart Gapped");
    ASSERT_NE(smp, std::string::npos) << text;
    const std::string tail = text.substr(smp);
    EXPECT_NE(tail.find("        42\n"), std::string::npos) << tail;
    EXPECT_NE(tail.find("        5\n"), std::string::npos) << tail;
    std::error_code ec;
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, SequentialDeckIsNotAffectedByIdPreservation) {
    // The "only when it matters" contract: a plain 1..n deck gets no
    // point_data/cell_data["mdpa:id"] at all, so a re-write stays on the
    // exact old renumbering code path.
    const Mesh m = mdpa_read_string(
        "Begin Nodes\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\nEnd Nodes\n"
        "Begin Elements Element3D4N\n1 0 1 2 3 4\nEnd Elements\n");
    EXPECT_FALSE(m.HasPointData(meshioplusplus::kMdpaIdName));
    EXPECT_FALSE(m.HasCellData(meshioplusplus::kMdpaIdName));

    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m);
    const std::string text = mdpa_slurp(out_path);
    EXPECT_NE(text.find("Begin Nodes\n 1 "), std::string::npos) << text;
    EXPECT_NE(text.find("Begin Elements Element3D4N\n  1 "), std::string::npos) << text;
    std::error_code ec;
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, DuplicateNodeIdInMdpaIdThrowsOnWrite) {
    Mesh m = mt::tri_mesh();
    NDArray ids(meshioplusplus::DType::Int64, {m.NumPoints()});
    std::int64_t* ip = ids.As<std::int64_t>();
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        ip[i] = 5;  // every point claims the same id
    m.AddPointData(meshioplusplus::kMdpaIdName, std::move(ids));
    const std::string out_path = mt::temp_path(".mdpa");
    EXPECT_THROW(meshioplusplus::write_mdpa(out_path, m), meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, DuplicateElementIdInMdpaIdThrowsOnWrite) {
    Mesh m = mt::tri_quad_mesh();  // triangle, quad, triangle -- 3 blocks
    std::vector<NDArray> ids;
    ids.reserve(m.NumCellBlocks());
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        NDArray a(meshioplusplus::DType::Int64, {m.Cells(b).NumCells()});
        std::int64_t* ap = a.As<std::int64_t>();
        for (std::size_t r = 0; r < m.Cells(b).NumCells(); ++r)
            ap[r] = 7;  // every element claims the same id
        ids.push_back(std::move(a));
    }
    m.AddCellData(meshioplusplus::kMdpaIdName, std::move(ids));
    const std::string out_path = mt::temp_path(".mdpa");
    EXPECT_THROW(meshioplusplus::write_mdpa(out_path, m), meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, MismatchedMdpaIdShapeFallsBackRatherThanCrashing) {
    // A wrong-length array is treated as unrelated/stale metadata: the writer
    // falls back to the old renumbering instead of reading out of bounds.
    Mesh m = mt::tri_mesh();
    NDArray ids(meshioplusplus::DType::Int64, {m.NumPoints() - 1});  // too short
    std::int64_t* ip = ids.As<std::int64_t>();
    for (std::size_t i = 0; i < ids.Size(); ++i)
        ip[i] = static_cast<std::int64_t>(i) + 100;
    m.AddPointData(meshioplusplus::kMdpaIdName, std::move(ids));
    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m);
    const std::string text = mdpa_slurp(out_path);
    EXPECT_NE(text.find("Begin Nodes\n 1 "), std::string::npos) << text;
    std::error_code ec;
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, EmptyMeshRoundTrips) {
    const std::string path = mt::temp_path(".mdpa");
    Mesh empty;
    meshioplusplus::write_mdpa(path, empty);
    const Mesh out = meshioplusplus::read_mdpa(path);
    EXPECT_EQ(out.NumPoints(), 0u);
    EXPECT_EQ(out.NumCellBlocks(), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---------------------------------------------------------------------------
// MdpaInfo: production decks (Properties bodies, application entity names)
// ---------------------------------------------------------------------------

namespace {

/// A deck shaped like a real Kratos one: every construct that used to throw.
const char* const kMdpaProductionDeck = R"(Begin ModelPartData
    DOMAIN_SIZE 3
End ModelPartData

Begin Properties 1
    DENSITY 7850.0
    YOUNG_MODULUS 2.1e+11
    POISSON_RATIO 0.3
    CONSTITUTIVE_LAW LinearElastic3DLaw
    LOCAL_AXES [3] (1.0, 0.0, 0.0)
    Begin Table 4 TEMPERATURE YOUNG_MODULUS
        20 2.1e+11
        100 2e+11
    End Table
End Properties

Begin Properties 2
    DENSITY 2700
End Properties

Begin Nodes
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 0.0 0.0 1.0
End Nodes

Begin Elements SmallDisplacementElement3D4N
1 1 1 2 3 4
End Elements

Begin Elements TotalLagrangianElement3D4N
2 2 1 2 3 4
End Elements

Begin Conditions SurfaceCondition3D3N
1 1 1 2 3
End Conditions

Begin SubModelPart Structure
    Begin SubModelPart Structure.Loads
        Begin SubModelPartConditions
            1
        End SubModelPartConditions
    End SubModelPart
    Begin SubModelPartNodes
        1
        2
    End SubModelPartNodes
End SubModelPart
)";

/// Every `Begin Elements|Conditions|Properties` header line of a file, in order.
std::vector<std::string> mdpa_headers(const std::string& rPath) {
    std::vector<std::string> out;
    std::ifstream in(rPath);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("Begin Elements", 0) == 0 || line.rfind("Begin Conditions", 0) == 0 ||
            line.rfind("Begin Properties", 0) == 0)
            out.push_back(line);
    }
    return out;
}

}  // namespace

TEST(Mdpa, ProductionDeckThrowsByDefaultAndNamesTheConstruct) {
    // The Properties bodies no longer throw, but the nested SubModelPart's
    // unsupported siblings still do -- so a strict read of a deck carrying one
    // fails by name, exactly as before.
    const std::string path = mdpa_temp_file(std::string(kMdpaProductionDeck) +
                                            "\nBegin Constraints LinearMasterSlave\n"
                                            "1 1 2\nEnd Constraints\n");
    EXPECT_THROW(meshioplusplus::read_mdpa(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, LenientSkipsUnsupportedBlocksAndRecordsThem) {
    const std::string path = mdpa_temp_file(std::string(kMdpaProductionDeck) +
                                            "\nBegin Constraints LinearMasterSlave\n"
                                            "1 1 2\nEnd Constraints\n"
                                            "\nBegin Geometries Triangle3D3\n"
                                            "1 1 2 3\nEnd Geometries\n");
    meshioplusplus::MdpaInfo info;
    meshioplusplus::ReadOptions opts;
    opts.mLenient = true;
    const Mesh m = meshioplusplus::read_mdpa(path, info, opts);

    EXPECT_EQ(m.NumPoints(), 4u);
    EXPECT_EQ(m.NumCellBlocks(), 3u);
    ASSERT_EQ(info.mSkippedConstructs.size(), 2u);
    EXPECT_NE(info.mSkippedConstructs[0].find("Constraints"), std::string::npos);
    EXPECT_NE(info.mSkippedConstructs[1].find("Geometries"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mdpa, PropertiesBodiesRoundTripThroughMdpaInfo) {
    const std::string path = mdpa_temp_file(kMdpaProductionDeck);
    meshioplusplus::MdpaInfo info;
    const Mesh m = meshioplusplus::read_mdpa(path, info);

    ASSERT_EQ(info.mProperties.size(), 2u);
    EXPECT_EQ(info.mProperties[0].mId, 1);
    EXPECT_EQ(info.mProperties[1].mId, 2);

    const auto& p1 = info.mProperties[0].mValues;
    ASSERT_EQ(p1.size(), 6u);  // five KEY value lines plus the inline table
    EXPECT_EQ(p1[0].mKey, "DENSITY");
    EXPECT_FALSE(p1[0].IsText());
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(p1[0].mValues, 0), 7850.0);
    // A constitutive-law name has no numeric form; it is kept verbatim, which
    // is both lossless and what the pure-Python reference does.
    EXPECT_EQ(p1[3].mKey, "CONSTITUTIVE_LAW");
    EXPECT_TRUE(p1[3].IsText());
    EXPECT_EQ(p1[3].mText, "LinearElastic3DLaw");
    // A bracketed vector likewise stays text -- lossless, not parsed.
    EXPECT_EQ(p1[4].mKey, "LOCAL_AXES");
    EXPECT_TRUE(p1[4].IsText());
    EXPECT_EQ(p1[4].mText, "[3] (1.0, 0.0, 0.0)");
    // The inline table is the one non-scalar that IS parsed.
    const meshioplusplus::PropertyValue* p_table = nullptr;
    for (const auto& v : p1)
        if (v.mIsTable)
            p_table = &v;
    ASSERT_NE(p_table, nullptr);
    EXPECT_EQ(p_table->mKey, "4 TEMPERATURE YOUNG_MODULUS");
    ASSERT_EQ(p_table->mValues.Shape().size(), 2u);
    EXPECT_EQ(p_table->mValues.Shape()[0], 2u);
    EXPECT_EQ(p_table->mValues.Shape()[1], 2u);

    // Write it back and read again: the properties must survive unchanged.
    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m, info);
    meshioplusplus::MdpaInfo info2;
    meshioplusplus::read_mdpa(out_path, info2);
    ASSERT_EQ(info2.mProperties.size(), info.mProperties.size());
    for (std::size_t i = 0; i < info.mProperties.size(); ++i) {
        EXPECT_EQ(info2.mProperties[i].mId, info.mProperties[i].mId);
        ASSERT_EQ(info2.mProperties[i].mValues.size(), info.mProperties[i].mValues.size());
        for (std::size_t j = 0; j < info.mProperties[i].mValues.size(); ++j) {
            const auto& a = info.mProperties[i].mValues[j];
            const auto& b = info2.mProperties[i].mValues[j];
            EXPECT_EQ(b.mKey, a.mKey);
            EXPECT_EQ(b.mText, a.mText);
            EXPECT_EQ(b.mIsTable, a.mIsTable);
            ASSERT_EQ(b.mValues.Size(), a.mValues.Size());
            for (std::size_t k = 0; k < a.mValues.Size(); ++k)
                EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(b.mValues, k),
                                 meshioplusplus::detail::read_double(a.mValues, k));
        }
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, EntityNamesRoundTripAndDegradeWithoutInfo) {
    const std::string path = mdpa_temp_file(kMdpaProductionDeck);
    meshioplusplus::MdpaInfo info;
    const Mesh m = meshioplusplus::read_mdpa(path, info);

    // Two same-type Elements blocks with different Kratos names stay separate,
    // so each keeps its own name rather than collapsing onto the first.
    ASSERT_EQ(m.NumCellBlocks(), 3u);
    ASSERT_EQ(info.mEntityNames.size(), 3u);
    EXPECT_EQ(info.mEntityNames[0].mName, "SmallDisplacementElement3D4N");
    EXPECT_FALSE(info.mEntityNames[0].mIsCondition);
    EXPECT_EQ(info.mEntityNames[1].mName, "TotalLagrangianElement3D4N");
    EXPECT_EQ(info.mEntityNames[2].mName, "SurfaceCondition3D3N");
    EXPECT_TRUE(info.mEntityNames[2].mIsCondition);

    // With the info, the headers come back byte for byte.
    const std::string with_info = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(with_info, m, info);
    const std::vector<std::string> in_headers = mdpa_headers(path);
    const std::vector<std::string> out_headers = mdpa_headers(with_info);
    EXPECT_EQ(out_headers, in_headers);

    // Without it, the names degrade to the canonical ones -- the documented
    // behaviour of the two-argument overload, pinned so it cannot drift.
    const std::string no_info = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(no_info, m);
    const std::vector<std::string> plain = mdpa_headers(no_info);
    ASSERT_EQ(plain.size(), 5u);  // Properties 1, 2 + three entity blocks
    EXPECT_EQ(plain[2], "Begin Elements Element3D4N");
    EXPECT_EQ(plain[3], "Begin Elements Element3D4N");

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(with_info, ec);
    std::filesystem::remove(no_info, ec);
}

TEST(Mdpa, PropertiesDeclaredForEveryReferencedIdWithoutInfo) {
    // The writer used to emit a single hard-coded `Properties 0` while the rows
    // wrote their gmsh:physical value, producing a file that references
    // undeclared properties. Every referenced id must now be declared.
    const std::string path = mdpa_temp_file(kMdpaProductionDeck);
    const Mesh m = meshioplusplus::read_mdpa(path);
    const std::string out_path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out_path, m);
    const std::vector<std::string> headers = mdpa_headers(out_path);
    ASSERT_GE(headers.size(), 2u);
    EXPECT_EQ(headers[0], "Begin Properties 1");
    EXPECT_EQ(headers[1], "Begin Properties 2");
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(out_path, ec);
}

TEST(Mdpa, UntaggedMeshStillEmitsExactlyPropertiesZero) {
    // The byte-identity guard for the change above: a mesh with no
    // gmsh:physical must produce the same two lines it always did.
    const std::string path = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(path, mt::tet_mesh());
    const std::string text = mdpa_slurp(path);
    EXPECT_NE(text.find("Begin Properties 0\nEnd Properties\n"), std::string::npos);
    EXPECT_EQ(text.find("Begin Properties 1"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// --- v9.2.0: Properties bodies ride on the Mesh ----------------------------
//
// Through v9.1.0 they rode the MdpaInfo side channel only, which nothing
// reachable from registry_readers() could ask for -- so a registry-based
// consumer (i.e. every one that does not link formats/mdpa.hpp directly) got
// the ids and no material data, and a write emitted empty blocks.

TEST(Mdpa, PropertiesReachTheMeshWithoutAnMdpaInfo) {
    const std::string path = mt::temp_path(".mdpa");
    {
        std::ofstream f(path);
        f << "Begin Properties 1\n"
             "  DENSITY 7850.0\n"
             "  CONSTITUTIVE_LAW LinearElastic3DLaw\n"
             "End Properties\n\n"
             "Begin Properties 2\n"
             "  DENSITY 2700.0\n"
             "End Properties\n\n"
             "Begin Nodes\n  1 0.0 0.0 0.0\n  2 1.0 0.0 0.0\n  3 0.0 1.0 0.0\n"
             "End Nodes\n\n"
             "Begin Elements Element2D3N\n  1 1 1 2 3\nEnd Elements\n";
    }
    // The plain overload -- exactly what registry_readers() calls.
    const Mesh m = meshioplusplus::read_mdpa(path);
    ASSERT_EQ(m.NumPropertySets(), 2u);
    EXPECT_EQ(m.GetPropertySet(0).mId, 1);
    EXPECT_EQ(m.GetPropertySet(1).mId, 2);
    ASSERT_EQ(m.GetPropertySet(0).mValues.size(), 2u);
    EXPECT_EQ(m.GetPropertySet(0).mValues[0].mKey, "DENSITY");
    EXPECT_DOUBLE_EQ(m.GetPropertySet(0).mValues[0].mValues.As<double>()[0], 7850.0);
    EXPECT_EQ(m.GetPropertySet(0).mValues[1].mKey, "CONSTITUTIVE_LAW");
    EXPECT_EQ(m.GetPropertySet(0).mValues[1].mText, "LinearElastic3DLaw");
    std::remove(path.c_str());
}

TEST(Mdpa, PropertiesRoundTripThroughTheMeshAlone) {
    const std::string path = mt::temp_path(".mdpa");
    {
        std::ofstream f(path);
        f << "Begin Properties 1\n  DENSITY 7850.0\n  CONSTITUTIVE_LAW LinearElastic3DLaw\n"
             "End Properties\n\n"
             "Begin Nodes\n  1 0.0 0.0 0.0\n  2 1.0 0.0 0.0\n  3 0.0 1.0 0.0\n"
             "End Nodes\n\n"
             "Begin Elements Element2D3N\n  1 1 1 2 3\nEnd Elements\n";
    }
    const Mesh m = meshioplusplus::read_mdpa(path);
    const std::string out = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out, m);  // no MdpaInfo anywhere
    const Mesh back = meshioplusplus::read_mdpa(out);

    ASSERT_EQ(back.NumPropertySets(), 1u);
    ASSERT_EQ(back.GetPropertySet(0).mValues.size(), 2u);
    EXPECT_DOUBLE_EQ(back.GetPropertySet(0).mValues[0].mValues.As<double>()[0], 7850.0);
    EXPECT_EQ(back.GetPropertySet(0).mValues[1].mText, "LinearElastic3DLaw");
    std::remove(path.c_str());
    std::remove(out.c_str());
}

TEST(Mdpa, AMeshWithNoPropertySetsWritesExactlyTheOldBytes) {
    // The byte-identity guard for the new writer branch: a mesh carrying
    // gmsh:physical but no property sets must emit what it always did.
    Mesh m = mt::tri_mesh();
    std::vector<NDArray> tags;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        NDArray a(meshioplusplus::DType::Int64, {m.Cells(b).NumCells()});
        for (std::size_t c = 0; c < a.Size(); ++c)
            a.As<std::int64_t>()[c] = 0;
        tags.push_back(std::move(a));
    }
    m.AddCellData("gmsh:physical", std::move(tags));
    ASSERT_EQ(m.NumPropertySets(), 0u);

    const std::string out = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out, m);
    std::ifstream in(out);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("Begin Properties 0\nEnd Properties\n"), std::string::npos);
    EXPECT_EQ(text.find("Begin Properties 1"), std::string::npos);
    std::remove(out.c_str());
}

TEST(Mdpa, AnExplicitMdpaInfoStillWinsOverTheMeshChannel) {
    // MdpaInfo keeps its v9.1.0 meaning: it preserves FILE order, which the
    // mesh channel deliberately does not (it canonicalizes to ascending id).
    const std::string path = mt::temp_path(".mdpa");
    {
        std::ofstream f(path);
        f << "Begin Properties 2\n  DENSITY 2700.0\nEnd Properties\n\n"
             "Begin Properties 1\n  DENSITY 7850.0\nEnd Properties\n\n"
             "Begin Nodes\n  1 0.0 0.0 0.0\n  2 1.0 0.0 0.0\n  3 0.0 1.0 0.0\n"
             "End Nodes\n\n"
             "Begin Elements Element2D3N\n  1 1 1 2 3\nEnd Elements\n";
    }
    meshioplusplus::MdpaInfo info;
    const Mesh m = meshioplusplus::read_mdpa(path, info);
    // The mesh sorted them; the info kept the file's order.
    EXPECT_EQ(m.GetPropertySet(0).mId, 1);
    ASSERT_EQ(info.mProperties.size(), 2u);
    EXPECT_EQ(info.mProperties[0].mId, 2);

    const std::string out = mt::temp_path(".mdpa");
    meshioplusplus::write_mdpa(out, m, info);
    std::ifstream in(out);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_LT(text.find("Begin Properties 2"), text.find("Begin Properties 1"));
    std::remove(path.c_str());
    std::remove(out.c_str());
}
