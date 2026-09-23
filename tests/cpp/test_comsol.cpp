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
 * @file test_comsol.cpp
 * @brief COMSOL `.mphtxt`/`.mphbin`: Mesh versions 4 and 2 (parameter records,
 *        up/down pairs, lowest index 1), strings with blanks and '#', Selections
 *        as regions, several Mesh objects, quadratic node order, the binary
 *        twin, the writer's entities and Selections, and the error paths. The
 *        Python tests (tests/python/test_mphtxt.py) cover the generated
 *        fixtures and byte parity between the engines.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;
using meshioplusplus::WriteError;
namespace detail = meshioplusplus::detail;

std::string cs_write_file(const std::string& rBody, const char* pExt = ".mphtxt") {
    const std::string path = mt::temp_path(pExt);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::string cs_read_bytes(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

using CsRegions = std::map<std::tuple<std::string, int>, std::vector<std::int64_t>>;

CsRegions cs_regions(const Mesh& rMesh) {
    CsRegions out;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        out[{reg.mName, reg.mDim}] =
            std::vector<std::int64_t>(reg.Entries(), reg.Entries() + reg.NumEntries());
    }
    return out;
}

std::vector<std::int64_t> cs_geom(const Mesh& rMesh, std::size_t Block) {
    const auto& g = rMesh.CellData("mphtxt:geom", Block);
    std::vector<std::int64_t> out;
    for (std::size_t k = 0; k < g.Size(); ++k)
        out.push_back(detail::read_int(g, k));
    return out;
}

std::vector<std::int64_t> cs_row(const Mesh& rMesh, std::size_t Block, std::size_t Row) {
    const auto cb = rMesh.Cells(Block);
    const std::size_t n = detail::cols(cb.Conn());
    std::vector<std::int64_t> out;
    for (std::size_t j = 0; j < n; ++j)
        out.push_back(detail::read_int(cb.Conn(), Row * n + j));
    return out;
}

// A version 4 file: a tetrahedron (domain 1), a quad boundary (entity 3) and
// two Selections, one label with a blank and a '#'.
const char* kTwoDomains =
    "# comment\n"
    "0 1\n"
    "3 # number of tags\n5 mesh1\n10 mesh1_sel1\n10 mesh1_sel2\n"
    "3 # number of types\n3 obj\n3 obj\n3 obj\n"
    "0 0 1\n4 Mesh # class\n4 # version\n3 # sdim\n"
    "5 # number of mesh vertices\n0 # lowest mesh vertex index\n"
    "0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n"
    "2 # number of element types\n"
    "4 quad # type name\n4\n1\n0 1 2 3\n1\n3\n"
    "3 tet # type name\n4\n1\n0 1 2 4\n1\n1\n"
    "0 0 1\n9 Selection # class\n0\n8 Solid #1 # Label\n5 mesh1\n3\n1\n1\n"
    "0 0 1\n9 Selection\n0\n4 base\n5 mesh1\n2\n1\n3\n";

}  // namespace

TEST(Comsol, ReadsVersion4WithSelections) {
    const Mesh mesh = meshioplusplus::read_mphtxt(cs_write_file(kTwoDomains));
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad");
    // COMSOL's tensor-order quad becomes a ring
    EXPECT_EQ(cs_row(mesh, 0, 0), (std::vector<std::int64_t>{0, 1, 3, 2}));
    EXPECT_EQ(cs_geom(mesh, 0), (std::vector<std::int64_t>{3}));
    EXPECT_EQ(cs_geom(mesh, 1), (std::vector<std::int64_t>{1}));
    const CsRegions expected = {{{"Solid #1", 3}, {1}}, {{"base", 2}, {0}}};
    EXPECT_EQ(cs_regions(mesh), expected);
}

TEST(Comsol, ReadsVersion2ParametersAndUpDown) {
    // lowest index 1; an edg2 whose parameter records hold one value per
    // node, a tri whose records hold three (a 3-D boundary), up/down pairs.
    const std::string body =
        "0 1\n1 5 mesh1\n1 3 obj\n0 0 1\n4 Mesh\n2\n3\n4\n1\n"
        "0 0 0\n1 0 0\n0 1 0\n0.5 0 0\n"
        "2\n"
        "4 edg2\n3\n1\n1 2 4\n3\n1\n0 1 0.5\n1\n0\n0\n"
        "3 tri\n3\n1\n1 2 3\n3\n1\n0 0 7 1 0 7 0 1 7\n1\n5\n1\n0 1\n";
    const Mesh mesh = meshioplusplus::read_mphtxt(cs_write_file(body));
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "line3");
    EXPECT_EQ(cs_row(mesh, 0, 0), (std::vector<std::int64_t>{0, 1, 3}));
    EXPECT_EQ(cs_geom(mesh, 1), (std::vector<std::int64_t>{5}));
}

TEST(Comsol, SeveralMeshObjectsAreMerged) {
    const std::string body =
        "0 1\n3 5 mesh1 5 mesh2 10 mesh2_sel1\n3 3 obj 3 obj 3 obj\n"
        "0 0 1 4 Mesh 4 2 3 0  0 0 1 0 0 1  1 3 tri 3 1 0 1 2 1 1\n"
        "0 0 1 4 Mesh 4 2 3 0  2 0 3 0 2 1  1 3 tri 3 1 0 1 2 1 4\n"
        "0 0 1 9 Selection 0 5 right 5 mesh2 2 1 4\n";
    const Mesh mesh = meshioplusplus::read_mphtxt(cs_write_file(body));
    EXPECT_EQ(mesh.NumPoints(), 6u);
    EXPECT_EQ(cs_row(mesh, 1, 0), (std::vector<std::int64_t>{3, 4, 5}));
    const CsRegions expected = {{{"mesh1", -1}, {0}}, {{"mesh2", -1}, {1}}, {{"right", 2}, {1}}};
    EXPECT_EQ(cs_regions(mesh), expected);
}

TEST(Comsol, QuadraticNodeOrder) {
    // tet2 in COMSOL order: corners, then lattice nodes (z, y, x) lexicographic:
    // m01, m02, m12, m03, m13, m23; meshio++ wants m01, m12, m20, m03, m13, m23.
    const std::string body =
        "0 1\n1 5 mesh1\n1 3 obj\n0 0 1\n4 Mesh\n4\n3\n10\n0\n"
        "0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
        "0.5 0 0\n0 0.5 0\n0.5 0.5 0\n0 0 0.5\n0.5 0 0.5\n0 0.5 0.5\n"
        "1\n4 tet2\n10\n1\n0 1 2 3 4 5 6 7 8 9\n1\n1\n";
    const Mesh mesh = meshioplusplus::read_mphtxt(cs_write_file(body));
    EXPECT_EQ(cs_row(mesh, 0, 0), (std::vector<std::int64_t>{0, 1, 2, 3, 4, 6, 5, 7, 8, 9}));
    // and back: the writer puts COMSOL's order into the file again
    const std::string out = mt::temp_path(".mphtxt");
    meshioplusplus::write_mphtxt(out, mesh);
    EXPECT_NE(cs_read_bytes(out).find("\n0 1 2 3 4 5 6 7 8 9\n"), std::string::npos);
}

TEST(Comsol, BinaryIsTheSameContent) {
    const Mesh text = meshioplusplus::read_mphtxt(cs_write_file(kTwoDomains));
    const std::string path = mt::temp_path(".mphbin");
    meshioplusplus::write_mphbin(path, text);
    const std::string bytes = cs_read_bytes(path);
    // int32 0, 1, then 3 tags; the first "mesh1" as int32 code points
    const std::int32_t head[] = {0, 1, 3, 5, 'm', 'e', 's', 'h', '1'};
    ASSERT_GE(bytes.size(), sizeof(head));
    EXPECT_EQ(std::memcmp(bytes.data(), head, sizeof(head)), 0);  // little-endian host
    const Mesh back = meshioplusplus::read_mphbin(path);
    EXPECT_EQ(cs_regions(back), cs_regions(text));
    EXPECT_EQ(cs_row(back, 0, 0), cs_row(text, 0, 0));
    EXPECT_EQ(cs_geom(back, 0), cs_geom(text, 0));
}

TEST(Comsol, WriterDerivesEntitiesAndSelections) {
    Mesh mesh = mt::tet_mesh();
    auto entries = [](std::vector<std::int64_t> v) {
        meshioplusplus::NDArray a(meshioplusplus::DType::Int64, {v.size()});
        std::copy(v.begin(), v.end(), a.As<std::int64_t>());
        return a;
    };
    mesh.AddRegion(Region("a", RegionKind::Cell, 3, -1, entries({0})));
    mesh.AddRegion(Region("b", RegionKind::Cell, 3, -1, entries({1})));
    mesh.AddRegion(Region("ab", RegionKind::Cell, 3, -1, entries({0, 1})));
    const std::string path = mt::temp_path(".mphtxt");
    meshioplusplus::write_mphtxt(path, mesh);
    const Mesh back = meshioplusplus::read_mphtxt(path);
    // regions in (kind, name, dim, tag) order: "a" is domain 1, "ab" overlaps
    // it and so defines none, "b" is domain 2; all three are Selections
    EXPECT_EQ(cs_geom(back, 0), (std::vector<std::int64_t>{1, 2}));
    const CsRegions expected = {{{"a", 3}, {0}}, {{"ab", 3}, {0, 1}}, {{"b", 3}, {1}}};
    EXPECT_EQ(cs_regions(back), expected);
    EXPECT_NE(cs_read_bytes(path).find("4 # version"), std::string::npos);

    Mesh bad = mt::tet_mesh();
    bad.AddRegion(Region("bad", RegionKind::Cell, entries({2})));
    EXPECT_THROW(meshioplusplus::write_mphtxt(mt::temp_path(".mphtxt"), bad), WriteError);
}

TEST(Comsol, ErrorPaths) {
    EXPECT_THROW(meshioplusplus::read_mphtxt(cs_write_file("0 2\n0\n0\n")), ReadError);
    EXPECT_THROW(meshioplusplus::read_mphtxt(cs_write_file("0 1\n1 5 geom1\n1 3 obj\n0 0 1\n"
                                                           "5 Geom3\n")),
                 ReadError);
    EXPECT_THROW(
        meshioplusplus::read_mphtxt(cs_write_file("0 1\n1 5 mesh1\n1 3 obj\n0 0 1\n4 Mesh\n4\n2\n"
                                                  "3\n0\n0 0\n1 0\n0 1\n1\n3 tri\n3\n1\n0 1 7\n"
                                                  "1\n1\n")),
        ReadError);
    EXPECT_THROW(
        meshioplusplus::read_mphbin(cs_write_file(std::string("\0\0\0\0\1\0\0\0", 8), ".mphbin")),
        ReadError);
}
