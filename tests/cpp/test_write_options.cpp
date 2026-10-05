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
#include <atomic>
#include <cstdio>
#include <cstring>
#include <filesystem>
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

// The browser viewer's surface path (`convertSurfaceOps(..., {compressVtp: false})`)
// writes its VTP through exactly this call, to skip a deflate the main thread
// would only undo.
TEST(WriteOptions, VtpWithTheNoneCodecIsUncompressedBase64AndReadsBackIdentically) {
    const Mesh m = mt::tri_mesh();
    const std::string raw = mt::temp_path("_wo_vtp_none.vtp");
    const std::string def = mt::temp_path("_wo_vtp_default.vtp");
    WriteOptions opts;
    opts.mEncoding = WriteEncoding::Binary;  // a codec alone would select ASCII
    opts.mCodecSet = true;
    opts.mCodec = meshioplusplus::detail::VtkCodec::None;
    std::string why;
    ASSERT_TRUE(registry_write_supports("vtp", opts, why)) << why;
    registry_write_ex(raw, m, "vtp", opts);
    registry_write_ex(def, m, "vtp", WriteOptions{});

    const std::string raw_text = read_all(raw);
    EXPECT_EQ(raw_text.find("compressor="), std::string::npos);
    EXPECT_NE(raw_text.find("format=\"binary\""), std::string::npos);  // base64, not ASCII
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    EXPECT_NE(read_all(def).find("compressor=\"vtkZLibDataCompressor\""), std::string::npos);
#endif

    const Mesh a = meshioplusplus::registry_read(raw, "vtp", {});
    const Mesh b = meshioplusplus::registry_read(def, "vtp", {});
    ASSERT_EQ(a.NumPoints(), m.NumPoints());
    ASSERT_EQ(a.NumPoints(), b.NumPoints());
    ASSERT_EQ(a.Points().Nbytes(), b.Points().Nbytes());
    EXPECT_EQ(0, std::memcmp(a.Points().Data(), b.Points().Data(), a.Points().Nbytes()));
    ASSERT_EQ(a.NumCellBlocks(), b.NumCellBlocks());
    for (std::size_t i = 0; i < a.NumCellBlocks(); ++i) {
        ASSERT_EQ(a.Cells(i).Conn().Nbytes(), b.Cells(i).Conn().Nbytes());
        EXPECT_EQ(0, std::memcmp(a.Cells(i).Conn().Data(), b.Cells(i).Conn().Data(),
                                 a.Cells(i).Conn().Nbytes()));
    }
    // The trap this guards against: a codec without an encoding writes ASCII.
    WriteOptions codec_only;
    codec_only.mCodecSet = true;
    codec_only.mCodec = meshioplusplus::detail::VtkCodec::None;
    const std::string ascii = mt::temp_path("_wo_vtp_codec_only.vtp");
    registry_write_ex(ascii, m, "vtp", codec_only);
    EXPECT_NE(read_all(ascii).find("format=\"ascii\""), std::string::npos);
    std::remove(ascii.c_str());
    std::remove(raw.c_str());
    std::remove(def.c_str());
}

TEST(WriteOptions, PcdLzfSelectsCompressedBinaryWithoutAnEncodingOverride) {
    const Mesh m = mt::tet_mesh();
    const std::string path = mt::temp_path("_wo_lzf.pcd");
    WriteOptions opts;
    opts.mCodecSet = true;
    opts.mCodec = meshioplusplus::detail::VtkCodec::LZF;
    std::string why;
    EXPECT_TRUE(registry_write_supports("pcd", opts, why)) << why;
    registry_write_ex(path, m, "pcd", opts);
    EXPECT_NE(read_all(path).find("DATA binary_compressed\n"), std::string::npos);
    EXPECT_EQ(meshioplusplus::registry_read(path, "pcd", {}).NumPoints(), m.NumPoints());
    for (const char* fmt : {"vti", "vtu", "vtp", "gmsh"})
        EXPECT_FALSE(registry_write_supports(fmt, opts, why)) << fmt;
    opts.mEncoding = WriteEncoding::Ascii;
    EXPECT_FALSE(registry_write_supports("pcd", opts, why));
    EXPECT_THROW(registry_write_ex(path, m, "pcd", opts), meshioplusplus::WriteError);
    opts.mEncoding = WriteEncoding::Binary;
    EXPECT_TRUE(registry_write_supports("pcd", opts, why));
    opts.mCodec = meshioplusplus::detail::VtkCodec::Zlib;
    EXPECT_FALSE(registry_write_supports("pcd", opts, why));
    std::remove(path.c_str());
}

TEST(WriteOptions, RawAppendedIsVtuOnly) {
    const Mesh m = mt::tet_mesh();
    WriteOptions appended;
    appended.mEncoding = WriteEncoding::RawAppended;
    std::string why;
    EXPECT_TRUE(meshioplusplus::registry_write_supports("vtu", appended, why)) << why;
    for (const char* fmt : {"vtk", "vtp", "gmsh", "stl"}) {
        EXPECT_FALSE(meshioplusplus::registry_write_supports(fmt, appended, why)) << fmt;
        EXPECT_THROW(registry_write_ex(mt::temp_path("_wo_raw_bad"), m, fmt, appended),
                     meshioplusplus::WriteError)
            << fmt;
    }

    const std::string p = mt::temp_path("_wo_raw.vtu");
    registry_write_ex(p, m, "vtu", appended);
    const std::string text = read_all(p);
    EXPECT_NE(text.find("<AppendedData encoding=\"raw\">"), std::string::npos);
    const Mesh back = meshioplusplus::registry_read(p, "vtu", {});
    mt::expect_same_geometry(back, m);
    std::remove(p.c_str());
}

TEST(WriteOptions, OpenfoamEncodingReachesTheBinaryWriter) {
    // roadmap §1.1: openfoam joined the formats with an ASCII/binary variant.
    // Each case gets its own directory: unlike every other format here,
    // openfoam resolves `<parent>/constant/polyMesh` from a `.foam` path, so
    // two cases sharing one parent would collide on the same polyMesh dir.
    static std::atomic<unsigned> counter{0};
    const std::filesystem::path dir_a =
        std::filesystem::temp_directory_path() / ("meshio_wo_a_" + std::to_string(counter++));
    const std::filesystem::path dir_b =
        std::filesystem::temp_directory_path() / ("meshio_wo_b_" + std::to_string(counter++));
    const std::string a = (dir_a / "case.foam").string();
    const std::string b = (dir_b / "case.foam").string();

    const Mesh m = mt::tet_mesh();
    WriteOptions ascii;
    ascii.mEncoding = WriteEncoding::Ascii;
    WriteOptions binary;
    binary.mEncoding = WriteEncoding::Binary;
    registry_write_ex(a, m, "openfoam", ascii);
    registry_write_ex(b, m, "openfoam", binary);

    const std::string a_points = read_all((dir_a / "constant" / "polyMesh" / "points").string());
    const std::string b_points = read_all((dir_b / "constant" / "polyMesh" / "points").string());
    EXPECT_NE(a_points.find("format      ascii;"), std::string::npos);
    EXPECT_NE(b_points.find("format      binary;"), std::string::npos);
    EXPECT_NE(a_points, b_points);

    EXPECT_EQ(meshioplusplus::registry_read(a, "openfoam", {}).NumPoints(), m.NumPoints());
    EXPECT_EQ(meshioplusplus::registry_read(b, "openfoam", {}).NumPoints(), m.NumPoints());

    std::error_code ec;
    std::filesystem::remove_all(dir_a, ec);
    std::filesystem::remove_all(dir_b, ec);
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
