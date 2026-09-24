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
 * @file test_libmesh.cpp
 * @brief libMesh `.xda`/`.xdr` reader: an ASCII refined mesh with inline
 *        subdomains and p-levels, the same stream as XDR, the HEX20 node order,
 *        side sets carried to refined children, and the refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/libmesh.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::RegionKind;
namespace detail = meshioplusplus::detail;

std::string write_file(const std::string& rBody, const std::string& rSuffix) {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

// A QUAD4 on [0, 2]^2 refined into four children, with inline subdomain and
// p-level, a sideset on the parent's bottom edge (id 3, named) and a node set.
const char* kAmrXda =
    "libMesh-0.9.6\n"
    "5\t # number of elements\n"
    "9\t # number of nodes\n"
    ".\t # boundary condition specification file\n"
    ".\t # subdomain id specification file\n"
    "n/a\t # processor id specification file\n"
    ".\t # p-level specification file\n"
    "4\t # type size\n0\t # uid size\n0\t # pid size\n4\t # sid size\n4\t # p-level size\n"
    "4\t # eid size\n4\t # side size\n4\t # bid size\n"
    "1\t # subdomain id to name map\n"
    "1\t # vector length\n7\t \n"
    "1\t # vector length\nplate\t \n"
    "1\t # n_elem at level 0\n"
    "5 7 0 0 1 2 3\n"
    "4\t # n_elem at level 1\n"
    "5 0 7 1 0 4 8 7\n"
    "5 0 7 1 4 1 5 8\n"
    "5 0 7 2 8 5 2 6\n"
    "5 0 7 0 7 8 6 3\n"
    "0 0 0\n2 0 0\n2 2 0\n0 2 0\n1 0 0\n2 1 0\n1 2 0\n0 1 0\n1 1 0\n"
    "0\t # presence of unique ids\n"
    "1\t # sideset id to name map\n"
    "1\t # vector length\n3\t \n"
    "1\t # vector length\nbottom\t \n"
    "1\t # number of side boundary conditions\n"
    "0 0 3\n"
    "0\t # nodeset id to name map\n"
    "2\t # number of nodesets\n"
    "0 9\n1 9\n";

// The same value stream, XDR-encoded (every integer 4 bytes: pre-1.3.0).
class Xdr {
public:
    void U32(std::uint32_t v) {
        for (int b = 3; b >= 0; --b)
            mOut += static_cast<char>((v >> (8 * b)) & 0xff);
    }
    void Str(const std::string& s) {
        U32(static_cast<std::uint32_t>(s.size()));
        mOut += s;
        mOut.append((4 - s.size() % 4) % 4, '\0');
    }
    void F64(double d) {
        std::uint64_t v;
        std::memcpy(&v, &d, 8);
        for (int b = 7; b >= 0; --b)
            mOut += static_cast<char>((v >> (8 * b)) & 0xff);
    }
    std::string mOut;
};

std::string amr_xdr() {
    Xdr x;
    x.Str("libMesh-0.9.6");
    x.U32(5);
    x.U32(9);
    for (const char* s : {".", ".", "n/a", "."})
        x.Str(s);
    for (std::uint32_t v : {4u, 0u, 0u, 4u, 4u, 4u, 4u, 4u})
        x.U32(v);
    x.U32(1);  // subdomain names
    x.U32(1);
    x.U32(7);
    x.U32(1);
    x.Str("plate");
    x.U32(1);
    for (std::uint32_t v : {5u, 7u, 0u, 0u, 1u, 2u, 3u})
        x.U32(v);
    x.U32(4);
    const std::uint32_t kids[4][7] = {
        {5, 0, 7, 1, 0, 4, 8}, {5, 0, 7, 1, 4, 1, 5}, {5, 0, 7, 2, 8, 5, 2}, {5, 0, 7, 0, 7, 8, 6}};
    const std::uint32_t last[4] = {7, 8, 6, 3};
    for (int k = 0; k < 4; ++k) {
        for (std::uint32_t v : kids[k])
            x.U32(v);
        x.U32(last[k]);
    }
    const double xy[9][2] = {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 0},
                             {2, 1}, {1, 2}, {0, 1}, {1, 1}};
    for (const auto& p : xy) {
        x.F64(p[0]);
        x.F64(p[1]);
        x.F64(0.0);
    }
    x.U32(0);  // no unique ids
    x.U32(1);
    x.U32(1);
    x.U32(3);
    x.U32(1);
    x.Str("bottom");
    x.U32(1);
    for (std::uint32_t v : {0u, 0u, 3u})
        x.U32(v);
    x.U32(0);  // nodeset names
    x.U32(2);
    for (std::uint32_t v : {0u, 9u, 1u, 9u})
        x.U32(v);
    return x.mOut;
}

void expect_amr(const Mesh& rMesh) {
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "quad");
    EXPECT_EQ(rMesh.Cells(0).NumCells(), 4u);
    EXPECT_EQ(rMesh.NumPoints(), 9u);
    const auto& level = rMesh.CellData("libmesh:level", 0);
    const auto& p = rMesh.CellData("libmesh:p_level", 0);
    for (std::size_t r = 0; r < 4; ++r)
        EXPECT_EQ(detail::read_int(level, r), 1);
    EXPECT_EQ(detail::read_int(p, 2), 2);
    const std::size_t plate = rMesh.FindRegion("plate", RegionKind::Cell);
    ASSERT_NE(plate, Mesh::npos);
    EXPECT_EQ(rMesh.Region(plate).mTag, 7);
    EXPECT_EQ(rMesh.Region(plate).NumEntries(), 4u);
    // The parent's bottom edge reaches the two bottom children.
    const std::size_t bottom = rMesh.FindRegion("bottom", RegionKind::Side);
    ASSERT_NE(bottom, Mesh::npos);
    const auto& side = rMesh.Region(bottom);
    ASSERT_EQ(side.NumEntries(), 2u);
    const auto& pts = rMesh.Points();
    for (std::size_t k = 0; k < 2; ++k) {
        meshioplusplus::CellType type;
        std::vector<std::int64_t> nodes;
        ASSERT_TRUE(detail::facet_nodes(rMesh, side.Entries()[2 * k], side.Entries()[2 * k + 1],
                                        type, nodes));
        for (std::int64_t n : nodes)
            EXPECT_EQ(detail::read_double(pts, static_cast<std::size_t>(n) * 3 + 1), 0.0);
    }
    const std::size_t ns = rMesh.FindRegion("nodeset_9", RegionKind::Point);
    ASSERT_NE(ns, Mesh::npos);
    EXPECT_EQ(rMesh.Region(ns).NumEntries(), 2u);
}

}  // namespace

TEST(LibMesh, ReadsARefinedAsciiMesh) {
    expect_amr(meshioplusplus::read_libmesh(write_file(kAmrXda, ".xda")));
}

TEST(LibMesh, ReadsTheSameStreamAsXdr) {
    expect_amr(meshioplusplus::read_libmesh(write_file(amr_xdr(), ".xdr")));
}

TEST(LibMesh, Hex20IsReorderedIntoMeshioOrder) {
    // libMesh HEX20 on the unit cube: edges 8..11 bottom, 12..15 vertical,
    // 16..19 top (cell_hex20.C).
    const double c[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 5},
                          {2, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {7, 4}};
    std::string body =
        "libMesh-0.7.0+\n1\t #\n20\t #\nn/a\t #\nn/a\t #\nn/a\t #\nn/a\t #\n"
        "1\t # n_elem at level 0\n11";
    for (int k = 0; k < 20; ++k)
        body += " " + std::to_string(k);
    body += "\n";
    auto coord = [&](double x, double y, double z) {
        body += std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
    };
    for (const auto& p : c)
        coord(p[0], p[1], p[2]);
    for (const auto& ed : e)
        coord((c[ed[0]][0] + c[ed[1]][0]) / 2, (c[ed[0]][1] + c[ed[1]][1]) / 2,
              (c[ed[0]][2] + c[ed[1]][2]) / 2);
    const Mesh mesh = meshioplusplus::read_libmesh(write_file(body, ".xda"));
    ASSERT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    const auto& conn = mesh.Cells(0).Conn();
    const auto& pts = mesh.Points();
    // meshio++ slot 12 is the top edge 4-5: (0.5, 0, 1); slot 16 the vertical 0-4.
    const std::size_t n12 = static_cast<std::size_t>(detail::read_int(conn, 12));
    const std::size_t n16 = static_cast<std::size_t>(detail::read_int(conn, 16));
    EXPECT_EQ(detail::read_double(pts, n12 * 3 + 2), 1.0);
    EXPECT_EQ(detail::read_double(pts, n12 * 3 + 0), 0.5);
    EXPECT_EQ(detail::read_double(pts, n16 * 3 + 2), 0.5);
    EXPECT_EQ(detail::read_double(pts, n16 * 3 + 0), 0.0);
}

TEST(LibMesh, RefusesLegacyAndTruncatedFiles) {
    EXPECT_THROW(meshioplusplus::read_libmesh(write_file("DEAL 003:003\n1\n", ".xda")), ReadError);
    const std::string xdr = amr_xdr();
    EXPECT_THROW(meshioplusplus::read_libmesh(write_file(xdr.substr(0, xdr.size() / 2), ".xdr")),
                 ReadError);
    // A node id past the declared count.
    std::string bad = kAmrXda;
    bad.replace(bad.find("5 7 0 0 1 2 3"), 13, "5 7 0 0 1 2 99");
    EXPECT_THROW(meshioplusplus::read_libmesh(write_file(bad, ".xda")), ReadError);
}
