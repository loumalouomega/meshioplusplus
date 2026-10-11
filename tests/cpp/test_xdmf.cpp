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
#include <iterator>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "../../tests/fuzz/fuzz_bundle.hpp"
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

// A dataset larger than its DataItem's Dimensions is refused before it is read:
// the file's own extent is not bounded by the XML, and a corrupt one declaring
// millions of unwritten chunks made HDF5 1.10's H5Dread crawl (the fuzzer's
// timeout, tests/fuzz/regressions/xdmf).
TEST(Xdmf, HdfDatasetLargerThanItsDeclaredDimensionsIsRefused) {
    const std::string path = mt::temp_path(".xdmf");
    meshioplusplus::write_xdmf(path, mt::tri_mesh(), "HDF", -1);
    std::string text;
    {
        std::ifstream in(path);
        text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    const std::size_t geometry = text.find("<Geometry");
    ASSERT_NE(geometry, std::string::npos);
    const std::string key = "Dimensions=\"";
    const std::size_t begin = text.find(key, geometry);
    ASSERT_NE(begin, std::string::npos);
    const std::size_t value = begin + key.size();
    text.replace(value, text.find('"', value) - value, "1 1");  // stored: at least 3 x 3
    {
        std::ofstream out(path);
        out << text;
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(std::filesystem::path(path).replace_extension(".h5"), ec);
}

// The fuzzer's crash (tests/fuzz/regressions/xdmf/unfiltered-chunk-stored-short):
// an unfiltered chunked dataset whose one stored chunk holds 36 bytes where its 13
// int64 elements take 104. libhdf5 copies the elements out of the stored bytes
// without checking their size (a heap-buffer-overflow inside HDF5 1.10), so the
// reader refuses the dataset before it is read.
TEST(Xdmf, AnUnfilteredChunkStoredInFewerBytesIsAReadError) {
    const std::filesystem::path bundle = std::filesystem::path(__FILE__).parent_path() / ".." /
                                         "fuzz" / "regressions" / "xdmf" /
                                         "unfiltered-chunk-stored-short.bundle";
    ASSERT_TRUE(std::filesystem::exists(bundle)) << bundle;
    std::vector<std::uint8_t> bytes;
    {
        std::ifstream in(bundle, std::ios::binary);
        bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    }
    const std::filesystem::path dir =
        std::filesystem::temp_directory_path() / "meshioplusplus-unfiltered-chunk";
    std::filesystem::create_directories(dir);
    std::filesystem::path xml;
    for (const meshioplusplus::FuzzBundleEntry& entry :
         meshioplusplus::fuzz_bundle_decode(bytes.data(), bytes.size())) {
        const std::filesystem::path path = dir / entry.first;
        std::ofstream out(path, std::ios::binary);
        out.write(reinterpret_cast<const char*>(entry.second.data()),
                  static_cast<std::streamsize>(entry.second.size()));
        if (path.extension() == ".xdmf")
            xml = path;
    }
    ASSERT_FALSE(xml.empty());
    EXPECT_THROW(meshioplusplus::read_xdmf(xml.string()), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}
#endif

// The Dimensions of an XML DataItem name a product the text cannot hold. The
// reader sized the buffer from that product before reading the text, so a
// billion-by-billion item asked the allocator for more than it can hand out
// (the fuzzer's crash, tests/fuzz/regressions/xdmf/xml-dimensions-product-
// overflow). A count above the text's characters is refused before the buffer.
TEST(Xdmf, AnXmlDimensionsProductLargerThanItsTextIsAReadError) {
    const std::filesystem::path input = std::filesystem::path(__FILE__).parent_path() / ".." /
                                        "fuzz" / "regressions" / "xdmf" /
                                        "xml-dimensions-product-overflow";
    ASSERT_TRUE(std::filesystem::exists(input)) << input;
    EXPECT_THROW(meshioplusplus::read_xdmf(input.string()), meshioplusplus::ReadError);
}

// A negative extent cast to size_t wrapped to a huge count, and the product of
// the extents with it, so the reader sized the buffers from a garbage number.
TEST(Xdmf, NegativeDimensionsAreAReadError) {
    const std::string path = mt::temp_path(".xdmf");
    {
        std::ofstream out(path);
        out << "<?xml version=\"1.0\"?>\n<Xdmf Version=\"3.0\"><Domain><Grid>\n"
               "<Topology TopologyType=\"Triangle\" NumberOfElements=\"1\">"
               "<DataItem Dimensions=\"1 3\" Format=\"XML\">0 1 2</DataItem></Topology>\n"
               "<Geometry GeometryType=\"XY\">"
               "<DataItem Dimensions=\"-1 2\" Format=\"XML\">0 0 1 0 0 1</DataItem></Geometry>\n"
               "</Grid></Domain></Xdmf>\n";
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// A Binary payload shorter than its Dimensions declare: the buffer was sized
// from the declared count and the short read left the rest unwritten. The
// count is checked against the file's bytes before the buffer is sized.
TEST(Xdmf, ABinaryPayloadShorterThanItsDimensionsIsAReadError) {
    const std::string bin = mt::temp_path(".bin");
    {
        std::ofstream out(bin, std::ios::binary);
        const float two[2] = {0.0f, 1.0f};
        out.write(reinterpret_cast<const char*>(two), sizeof two);
    }
    const std::string path = mt::temp_path(".xdmf");
    {
        std::ofstream out(path);
        out << "<?xml version=\"1.0\"?>\n<Xdmf Version=\"3.0\"><Domain><Grid>\n"
               "<Topology TopologyType=\"Triangle\" NumberOfElements=\"1\">"
               "<DataItem Dimensions=\"1 3\" Format=\"XML\">0 1 2</DataItem></Topology>\n"
               "<Geometry GeometryType=\"XY\"><DataItem Dimensions=\"3 2\" Precision=\"4\" "
               "Format=\"Binary\">"
            << bin
            << "</DataItem></Geometry>\n"
               "</Grid></Domain></Xdmf>\n";
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(bin, ec);
}

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
