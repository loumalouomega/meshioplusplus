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
#include <filesystem>
#include <fstream>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xdmf.hpp"

namespace {
void rt(const mt::Mesh& mesh, const std::string& data_format) {
    mt::roundtrip([&](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_xdmf(p, m, data_format, -1); },
                  [](const std::string& p) { return meshioplusplus::read_xdmf(p); }, mesh, ".xdmf");
}
}  // namespace

TEST(Xdmf, Xml) {
    rt(mt::tri_mesh(), "XML");
    rt(mt::tet_mesh(), "XML");
    rt(mt::tri_quad_mesh(), "XML");  // Mixed topology
    rt(mt::tri_mesh_2d(), "XML");    // XY geometry
}
TEST(Xdmf, Binary) {
    rt(mt::tri_mesh(), "Binary");
    rt(mt::hex_mesh(), "Binary");
    rt(mt::tri_quad_mesh(), "Binary");
}
#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(Xdmf, Hdf) {
    rt(mt::tri_mesh(), "HDF");
    rt(mt::tet_mesh(), "HDF");
    rt(mt::tri_quad_mesh(), "HDF");
}
#endif

TEST(Xdmf, ReadRejectsMissingRoot) {
    // Well-formed XML that lacks the <Xdmf> root must raise rather than be
    // treated as an empty mesh.
    std::string path = mt::temp_path(".xdmf");
    {
        std::ofstream f(path);
        f << "<?xml version=\"1.0\"?>\n<NotXdmf/>\n";
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Xdmf, VersionTwoReferencesAndInformationAreNative) {
    const std::string path = mt::temp_path("_v2_reference.xdmf");
    {
        std::ofstream out(path);
        out << R"(<Xdmf Version="2.0"><Domain><Grid>
<Topology TopologyType="Triangle"><DataItem NumberType="Int" Dimensions="1 3" Format="XML">0 1 2</DataItem></Topology>
<Geometry GeometryType="XYZ"><DataItem Reference="XML">/Xdmf/Domain/Grid/Attribute/DataItem</DataItem></Geometry>
<Attribute Name="coords" Center="Node"><DataItem NumberType="Float" Precision="8" Dimensions="3 3" Format="XML">0 0 0 1 0 0 0 1 0</DataItem></Attribute>
<Information><![CDATA[<main><map key="wall" dim="2">7</map></main>]]></Information>
</Grid></Domain></Xdmf>)";
    }
    const auto mesh = meshioplusplus::read_xdmf(path);
    EXPECT_EQ(mesh.NumPoints(), 3);
    EXPECT_EQ(mesh.Cells(0).NumCells(), 1);
    ASSERT_TRUE(mesh.HasFieldData("wall"));
    EXPECT_EQ(mesh.FieldData("wall").As<std::int64_t>()[0], 7);
    const auto meta = meshioplusplus::read_xdmf_metadata(path);
    EXPECT_EQ(meta.mNumPoints, 3);
    ASSERT_EQ(meta.mCellBlocks.size(), 1);
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 1);
    EXPECT_EQ(meta.mFieldDataNames, std::vector<std::string>{"wall"});
    std::filesystem::remove(path);
}
