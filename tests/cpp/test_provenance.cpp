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
#include <fstream>
#include <sstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtu.hpp"

// Asserts the one-line provenance credit every writer emits (see
// detail/provenance.hpp and doc/roadmap.md's "audit and normalize" bullet)
// appears in the *written bytes* -- never inferred by reading the file back
// through our own reader, which cannot tell a right tag from a wrong one
// (the CGNS-ordering-suite rationale: a self-consistent round trip is not an
// oracle for the writer alone).

namespace {

std::string slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(Provenance, ObjCarriesTheTag) {
    std::string path = mt::temp_path(".obj");
    meshioplusplus::write_obj(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    EXPECT_NE(bytes.find(meshioplusplus::detail::kProvenanceTag), std::string::npos);
    // The tag is the very first content line.
    EXPECT_EQ(bytes.rfind(std::string("# ") + meshioplusplus::detail::kProvenanceTag, 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Provenance, OffCarriesTheTag) {
    std::string path = mt::temp_path(".off");
    meshioplusplus::write_off(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    std::string expected = std::string("OFF\n# ") + meshioplusplus::detail::kProvenanceTag;
    EXPECT_EQ(bytes.rfind(expected, 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Provenance, VtuCarriesTheTagAsAnXmlComment) {
    std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, mt::tri_mesh(), /*binary=*/false, /*zlib=*/false);
    std::string bytes = slurp(path);
    std::string expected = std::string("<!--") + meshioplusplus::detail::kProvenanceTag + "-->";
    EXPECT_NE(bytes.find(expected), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The binary STL header is a fixed 80-byte slot: the tag, then NUL padding,
// with no trailing newline -- unlike every other writer here.
TEST(Provenance, StlBinaryHeaderIsTagThenNulPadded) {
    std::string path = mt::temp_path(".stl");
    meshioplusplus::write_stl(path, mt::tri_mesh(), /*binary=*/true, /*skin=*/true);
    std::ifstream f(path, std::ios::binary);
    char header[80];
    f.read(header, 80);
    ASSERT_EQ(f.gcount(), 80);
    std::string tag = meshioplusplus::detail::kProvenanceTag;
    EXPECT_EQ(std::string(header, tag.size()), tag);
    for (std::size_t i = tag.size(); i < 80; ++i)
        EXPECT_EQ(header[i], '\0') << "byte " << i << " of the STL header is not NUL-padded";
    // Must not start with "solid" -- that is the ASCII/binary sniff heuristic
    // both the C++ and Python STL readers use.
    EXPECT_NE(tag.rfind("solid", 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// nastran is the one documented exception: the C++ reader gates on a
// sentinel comment line the Python writer never emits (doc/formats/nastran.md),
// so the sentinel stays first and the tag follows it as a second `$` line --
// both must be present, in that order, for the reader's own gate to keep
// working exactly as before this change.
TEST(Provenance, NastranSentinelPrecedesTheTag) {
    std::string path = mt::temp_path(".bdf");
    meshioplusplus::write_nastran(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    auto sentinel_pos = bytes.find("meshioplusplus-cpp-nastran");
    auto tag_pos = bytes.find(meshioplusplus::detail::kProvenanceTag);
    ASSERT_NE(sentinel_pos, std::string::npos);
    ASSERT_NE(tag_pos, std::string::npos);
    EXPECT_LT(sentinel_pos, tag_pos);
    // The reader must still accept its own output.
    EXPECT_NO_THROW(meshioplusplus::read_nastran(path));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The OpenFOAM banner box is fixed-width; the tag must not push the
// box's closing "|" out of column.
TEST(Provenance, OpenfoamBannerStaysColumnAligned) {
    std::string dir = mt::temp_path("_of");
    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::write_openfoam(dir, mt::hex_mesh(), info);
    std::string points_path = dir + "/constant/polyMesh/points";
    std::string bytes = slurp(points_path);
    auto tag_pos = bytes.find(meshioplusplus::detail::kProvenanceTag);
    ASSERT_NE(tag_pos, std::string::npos);
    auto line_start = bytes.rfind('\n', tag_pos) + 1;
    auto line_end = bytes.find('\n', tag_pos);
    std::string line = bytes.substr(line_start, line_end - line_start);
    // Every banner line (the box border, blank cells, and this one) is the
    // same fixed width, box borders included.
    EXPECT_EQ(line.size(), 79u) << "banner line: " << line;
    EXPECT_EQ(line.front(), '|');
    EXPECT_EQ(line.back(), '|');
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}
