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

// Bundle companion-file fuzzing: the MIOB encoding round-trips, refuses unsafe
// names, and reconstructs XDMF XML/heavy-data bundles through the production
// reader (DataItem/reference paths, not the shared XML parser alone).

// System includes
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../tests/fuzz/fuzz_bundle.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xdmf.hpp"

using meshioplusplus::FuzzBundleEntry;

namespace {

std::vector<std::uint8_t> fuzz_bundle_bytes(const std::string& rText) {
    return std::vector<std::uint8_t>(rText.begin(), rText.end());
}

std::string fuzz_bundle_tmp(const std::string& rName) {
    return (std::filesystem::temp_directory_path() / rName).string();
}

}  // namespace

TEST(FuzzBundle, EncodeDecodeRoundTrip) {
    const std::vector<FuzzBundleEntry> entries = {
        {"input.xdmf", fuzz_bundle_bytes("<Xdmf/>")},
        {"data.h5", fuzz_bundle_bytes("\x89HDF\r\n\x1a\n")},
    };
    const std::vector<std::uint8_t> blob = meshioplusplus::fuzz_bundle_encode(entries);
    EXPECT_TRUE(meshioplusplus::fuzz_bundle_is_bundle(blob.data(), blob.size()));
    const auto back = meshioplusplus::fuzz_bundle_decode(blob.data(), blob.size());
    ASSERT_EQ(back.size(), 2);
    EXPECT_EQ(back[0].first, "input.xdmf");
    EXPECT_EQ(back[1].first, "data.h5");
}

TEST(FuzzBundle, RejectsUnsafeNames) {
    EXPECT_FALSE(meshioplusplus::fuzz_bundle_name_ok("../evil.xdmf"));
    EXPECT_FALSE(meshioplusplus::fuzz_bundle_name_ok("/abs.xdmf"));
    EXPECT_FALSE(meshioplusplus::fuzz_bundle_name_ok("a/../../b"));
    EXPECT_FALSE(meshioplusplus::fuzz_bundle_name_ok(""));
    EXPECT_TRUE(meshioplusplus::fuzz_bundle_name_ok("input.bp/md.idx"));
    const std::vector<FuzzBundleEntry> bad = {{"../evil", fuzz_bundle_bytes("x")}};
    EXPECT_THROW(meshioplusplus::fuzz_bundle_encode(bad), meshioplusplus::ReadError);
}

TEST(FuzzBundle, RejectsTruncatedInputs) {
    const std::vector<std::uint8_t> short_magic = {'M', 'I', 'O'};
    EXPECT_FALSE(meshioplusplus::fuzz_bundle_is_bundle(short_magic.data(), short_magic.size()));
    const std::vector<std::uint8_t> bad_count = {'M', 'I', 'O', 'B', 1, 0, 0};
    EXPECT_THROW(meshioplusplus::fuzz_bundle_decode(bad_count.data(), bad_count.size()),
                 meshioplusplus::ReadError);
}

TEST(FuzzBundle, XdmfXmlBundleReadsThroughProductionReader) {
    const std::string dir = fuzz_bundle_tmp("mio-fuzz-bundle-xdmf-xml");
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    std::filesystem::create_directories(dir, ec);
    const std::string xml = R"(<Xdmf Version="3.0"><Domain><Grid GridType="Uniform">)"
                            R"(<Topology TopologyType="Triangle"><DataItem Dimensions="1 3" )"
                            R"(NumberType="Int" Precision="4" Format="XML">0 1 2</DataItem>)"
                            R"(</Topology><Geometry GeometryType="XYZ"><DataItem Dimensions="3 3" )"
                            R"(NumberType="Float" Precision="8" Format="XML">0 0 0 1 0 0 0 1 0)"
                            R"(</DataItem></Geometry></Grid></Domain></Xdmf>)";
    const std::vector<FuzzBundleEntry> entries = {{"input.xdmf", fuzz_bundle_bytes(xml)}};
    const auto blob = meshioplusplus::fuzz_bundle_encode(entries);
    const auto decoded = meshioplusplus::fuzz_bundle_decode(blob.data(), blob.size());
    const std::string path = dir + "/input.xdmf";
    {
        std::ofstream f(path, std::ios::binary);
        f.write(reinterpret_cast<const char*>(decoded.front().second.data()),
                static_cast<std::streamsize>(decoded.front().second.size()));
    }
    const auto mesh = meshioplusplus::read_xdmf(path);
    EXPECT_EQ(mesh.NumPoints(), 3);
    std::filesystem::remove_all(dir, ec);
}

TEST(FuzzBundle, XdmfReferenceAndCycleRefused) {
    const std::string dir = fuzz_bundle_tmp("mio-fuzz-bundle-xdmf-ref");
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    std::filesystem::create_directories(dir, ec);
    const std::string good = R"(<Xdmf Version="3.0"><Domain><Grid GridType="Uniform">)"
                             R"(<Topology TopologyType="Triangle"><DataItem Dimensions="1 3" )"
                             R"(NumberType="Int" Precision="4" Format="XML">0 1 2</DataItem>)"
                             R"(</Topology><Geometry GeometryType="XYZ">)"
                             R"(<DataItem Reference="/Xdmf/Domain/Grid/Attribute/DataItem"/>)"
                             R"(</Geometry><Attribute Name="pts" Center="Node">)"
                             R"(<DataItem Dimensions="3 3" NumberType="Float" Precision="8" )"
                             R"(Format="XML">0 0 0 1 0 0 0 1 0</DataItem></Attribute>)"
                             R"(</Grid></Domain></Xdmf>)";
    const std::string path = dir + "/input.xdmf";
    {
        std::ofstream f(path);
        f << good;
    }
    EXPECT_EQ(meshioplusplus::read_xdmf(path).NumPoints(), 3);
    const std::string cyclic = R"(<Xdmf Version="3.0"><Domain><Grid GridType="Uniform">)"
                               R"(<Topology TopologyType="Triangle">)"
                               R"(<DataItem Reference="/Xdmf/Domain/Grid/Topology/DataItem"/>)"
                               R"(</Topology></Grid></Domain></Xdmf>)";
    {
        std::ofstream f(path);
        f << cyclic;
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::filesystem::remove_all(dir, ec);
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(FuzzBundle, XdmfHdfBundleReadsHeavyData) {
    mt::roundtrip([&](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_xdmf(p, m, "HDF", -1); },
                  [](const std::string& p) { return meshioplusplus::read_xdmf(p); }, mt::tri_mesh(),
                  ".xdmf");
}
#endif
