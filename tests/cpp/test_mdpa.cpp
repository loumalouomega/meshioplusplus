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
        {"Begin Properties 1\n    DENSITY 1.0\nEnd Properties\n", "Properties"},
        {"Begin ModelPartData\n  NAME \"a string\"\nEnd ModelPartData\n", "ModelPartData"},
        {"Begin Nodes\n1 0 0 0\n3 1 1 1\nEnd Nodes\n", "node ids"},
        {"Begin Nodes\n1 0 0 0\n", "End Nodes"},
        {"Begin SubModelPart S\n  Begin SubModelPartData\n    K 1\n  End SubModelPartData\n"
         "End SubModelPart\n",
         "SubModelPartData"},
        {"Begin Constraints Foo\nEnd Constraints\n", "unsupported block"},
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

TEST(Mdpa, UnknownEntityNameThrows) {
    const std::string path = mdpa_temp_file(
        "Begin Nodes\n1 0 0 0\n2 1 0 0\nEnd Nodes\n"
        "Begin Elements TotallyUnknownThing\n1 0 1 2\nEnd Elements\n");
    EXPECT_THROW(meshioplusplus::read_mdpa(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
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
