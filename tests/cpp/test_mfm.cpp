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

// System includes
#include <filesystem>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mfm.hpp"

namespace {

void mfm_roundtrip(const mt::Mesh& mesh) {
    mt::roundtrip(
        [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_mfm(p, m, ".16e"); },
        [](const std::string& p) { return meshioplusplus::read_mfm(p); }, mesh, ".mfm");
}

mt::Mesh mfm_read_text(const std::string& rText) {
    const std::string path = mt::temp_path(".mfm");
    {
        auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
        out << rText;
    }
    auto mesh = meshioplusplus::read_mfm(path);
    // Results must own their buffers after a mapped source is released.
    std::filesystem::remove(path);
    {
        auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
        out << std::string(rText.size(), '?');
    }
    std::filesystem::remove(path);
    return mesh;
}

void mfm_expect_triangle(const mt::Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumPoints(), 3u);
    ASSERT_EQ(rMesh.PointDim(), 2u);
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    const auto block = rMesh.Cells(0);
    EXPECT_EQ(block.Type(), "triangle");
    for (std::size_t i = 0; i < 3; ++i)
        EXPECT_EQ(block.Conn().As<std::int64_t>()[i], static_cast<std::int64_t>(i));
    const double expected[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_EQ(rMesh.Points().As<double>()[i], expected[i]);
    EXPECT_EQ(rMesh.CellData("mfm:ref", 0).As<std::int64_t>()[0], 7);
}

}  // namespace

TEST(Mfm, Line) {
    mfm_roundtrip(mt::line_mesh());
}
TEST(Mfm, Triangle) {
    mfm_roundtrip(mt::tri_mesh());
}
TEST(Mfm, Triangle2D) {
    mfm_roundtrip(mt::tri_mesh_2d());
}
TEST(Mfm, Quad) {
    mfm_roundtrip(mt::quad_mesh());
}
TEST(Mfm, Tetra) {
    mfm_roundtrip(mt::tet_mesh());
}
TEST(Mfm, Hexahedron) {
    mfm_roundtrip(mt::hex_mesh());
}
TEST(Mfm, Wedge) {
    mfm_roundtrip(mt::wedge_mesh());
}

TEST(Mfm, RejectsMixedTypes) {
    std::string path = mt::temp_path(".mfm");
    EXPECT_THROW(meshioplusplus::write_mfm(path, mt::tri_quad_mesh(), ".16e"),
                 meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Mfm, WhitespaceAndSourceOwnership) {
    const std::string body = "1 2 3\t0\v0\f0 0 0 0\r\n0 0 1 0 0 1\n7";
    mfm_expect_triangle(mfm_read_text("\n\t\r\n1 3 3 2 3 3 3 1\r\n" + body));
    // Above the default Auto threshold; the result must also own mapped input.
    mfm_expect_triangle(mfm_read_text(std::string(16u << 20, ' ') + "\n1 3 3 2 3 3 3 1\n" + body));
}

TEST(Mfm, LenientTokensAndIgnoredSections) {
    const std::string text =
        "description\n1 3 3 2 3 3 3 1 99 # ignored header tail\n"
        "+1suffix 2.0 3e0 ignored reference tokens are not parsed "
        "nonnumeric 0tail +1suffix nope -0suffix 1.0suffix +7suffix trailing data";
    mfm_expect_triangle(mfm_read_text(text));
}

TEST(Mfm, EmptyMeshAndHeaderWithoutFinalNewline) {
    const auto mesh = mfm_read_text("0 0 0 3 4 4 6 4");
    EXPECT_EQ(mesh.NumPoints(), 0u);
    EXPECT_EQ(mesh.PointDim(), 3u);
    EXPECT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).NumCells(), 0u);
    EXPECT_EQ(mesh.CellData("mfm:ref", 0).Size(), 0u);
}

TEST(Mfm, TruncatedSections) {
    const std::string tokens[] = {"1", "2", "3", "0", "0", "0", "0", "0",
                                  "0", "0", "0", "1", "0", "0", "1", "7"};
    std::string text = "1 3 3 2 3 3 3 1\n";
    const std::string path = mt::temp_path(".mfm");
    for (const auto& token : tokens) {
        {
            auto out = meshioplusplus::detail::make_classic_ofstream(path);
            out << text;
        }
        EXPECT_THROW(meshioplusplus::read_mfm(path), meshioplusplus::ReadError);
        text += token + " ";
    }
    std::filesystem::remove(path);
}
