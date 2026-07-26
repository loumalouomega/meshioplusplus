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
/**
 * @file test_write_options.cpp
 * @brief `registry_write_ex()` -- the single owner of parameterized writing,
 * shared by both CLIs and the C API's `mio_write_ex`.
 *
 * The rule under test throughout: an option a format cannot honour is an
 * ERROR, never silently ignored.
 */

// System includes
#include <cstdio>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::registry_write_ex;
using meshioplusplus::registry_write_supports;
using meshioplusplus::WriteEncoding;
using meshioplusplus::WriteOptions;

std::string read_all(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

}  // namespace

TEST(WriteOptions, DefaultsMatchTheRegistryWriterExactly) {
    // The all-defaults path must dispatch straight through the registry, so the
    // overwhelmingly common call stays byte-for-byte what it always was.
    const Mesh m = mt::tet_mesh();
    const std::string via_registry = mt::temp_path("_wo_reg.vtu");
    const std::string via_ex = mt::temp_path("_wo_ex.vtu");

    meshioplusplus::registry_writers().at("vtu")(via_registry, m);
    registry_write_ex(via_ex, m, "vtu", WriteOptions{});

    EXPECT_EQ(read_all(via_registry), read_all(via_ex));
    std::remove(via_registry.c_str());
    std::remove(via_ex.c_str());
}

TEST(WriteOptions, EncodingSelectsAsciiOrBinary) {
    const Mesh m = mt::tet_mesh();
    const std::string a = mt::temp_path("_wo_a.vtu");
    const std::string b = mt::temp_path("_wo_b.vtu");

    WriteOptions ascii;
    ascii.mEncoding = WriteEncoding::Ascii;
    WriteOptions binary;
    binary.mEncoding = WriteEncoding::Binary;
    registry_write_ex(a, m, "vtu", ascii);
    registry_write_ex(b, m, "vtu", binary);

    EXPECT_NE(read_all(a), read_all(b));
    // Both must still read back to the same mesh.
    EXPECT_EQ(meshioplusplus::registry_read(a, "vtu", {}).NumPoints(), m.NumPoints());
    EXPECT_EQ(meshioplusplus::registry_read(b, "vtu", {}).NumPoints(), m.NumPoints());
    std::remove(a.c_str());
    std::remove(b.c_str());
}

TEST(WriteOptions, UnsupportedOptionIsAnErrorNotSilentlyIgnored) {
    const Mesh m = mt::tet_mesh();
    std::string why;

    WriteOptions codec_on_gmsh;
    codec_on_gmsh.mCodecSet = true;
    EXPECT_FALSE(registry_write_supports("gmsh", codec_on_gmsh, why));
    EXPECT_NE(why.find("codec"), std::string::npos);
    EXPECT_THROW(registry_write_ex(mt::temp_path("_wo_x.msh"), m, "gmsh", codec_on_gmsh),
                 meshioplusplus::WriteError);

    WriteOptions ff_on_vtu;
    ff_on_vtu.mFloatFormat = ".3f";
    EXPECT_FALSE(registry_write_supports("vtu", ff_on_vtu, why));
    EXPECT_NE(why.find("float-format"), std::string::npos);
}

TEST(WriteOptions, TextOnlyFormatAcceptsAsciiAndRejectsBinary) {
    // "Write this as ASCII" is a sensible thing to ask of a text format -- it is
    // just the normal write -- so it succeeds; BINARY fails by name. This is the
    // set the Python CLI's `ascii` verb accepts, so the two CLIs agree on which
    // files `meshioplusplus ascii` handles.
    const Mesh m = mt::tet_mesh();
    const std::string path = mt::temp_path("_wo_text.mdpa");

    std::string why;
    WriteOptions ascii;
    ascii.mEncoding = WriteEncoding::Ascii;
    EXPECT_TRUE(registry_write_supports("mdpa", ascii, why));
    ASSERT_NO_THROW(registry_write_ex(path, m, "mdpa", ascii));
    EXPECT_FALSE(read_all(path).empty());

    WriteOptions binary;
    binary.mEncoding = WriteEncoding::Binary;
    try {
        registry_write_ex(path, m, "mdpa", binary);
        FAIL() << "expected a WriteError naming the format";
    } catch (const meshioplusplus::WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("mdpa"), std::string::npos);
        EXPECT_NE(std::string(e.what()).find("text-only"), std::string::npos);
    }
    std::remove(path.c_str());
}
