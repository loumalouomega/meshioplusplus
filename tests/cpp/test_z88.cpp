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
 * @file test_z88.cpp
 * @brief Z88 structure files: a hex20 and a tet10 in Z88's node order, the
 *        sibling displacement and stress files, the writer round trip, and the
 *        basename dispatch.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/z88.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
namespace detail = meshioplusplus::detail;
namespace fs = std::filesystem;

fs::path deck_dir() {
    const fs::path dir = fs::path(mt::temp_path("_z88"));
    fs::create_directories(dir);
    return dir;
}

void write(const fs::path& rPath, const std::string& rBody) {
    std::ofstream(rPath, std::ios::binary) << rBody;
}

// A unit hex20 in Z88's order (top face 1-4 first, then 5-8, the top ring, the
// bottom ring, the verticals) and a tet10 (mid-edges 1-2 2-3 3-1 2-4 3-4 1-4).
std::string structure() {
    std::string s = "3 30 2 90 0\n";
    const double top[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    std::vector<std::array<double, 3>> p;
    for (const auto& q : top)
        p.push_back({q[0], q[1], 1.0});
    for (const auto& q : top)
        p.push_back({q[0], q[1], 0.0});
    auto mid = [&](int a, int b) {
        return std::array<double, 3>{(p[a][0] + p[b][0]) / 2, (p[a][1] + p[b][1]) / 2,
                                     (p[a][2] + p[b][2]) / 2};
    };
    for (int k = 0; k < 4; ++k)
        p.push_back(mid(k, (k + 1) % 4));
    for (int k = 0; k < 4; ++k)
        p.push_back(mid(4 + k, 4 + (k + 1) % 4));
    for (int k = 0; k < 4; ++k)
        p.push_back(mid(k, k + 4));
    const std::size_t t0 = p.size();
    p.push_back({2, 0, 0});
    p.push_back({3, 0, 0});
    p.push_back({2, 1, 0});
    p.push_back({2, 0, 1});
    const int te[6][2] = {{0, 1}, {1, 2}, {2, 0}, {1, 3}, {2, 3}, {0, 3}};
    for (const auto& ed : te)
        p.push_back({(p[t0 + ed[0]][0] + p[t0 + ed[1]][0]) / 2,
                     (p[t0 + ed[0]][1] + p[t0 + ed[1]][1]) / 2,
                     (p[t0 + ed[0]][2] + p[t0 + ed[1]][2]) / 2});
    for (std::size_t k = 0; k < p.size(); ++k)
        s += std::to_string(k + 1) + " 3 " + std::to_string(p[k][0]) + " " +
             std::to_string(p[k][1]) + " " + std::to_string(p[k][2]) + "\n";
    s += "1 10\n";
    for (int k = 1; k <= 20; ++k)
        s += std::to_string(k) + (k < 20 ? " " : "\n");
    s += "2 16\n";
    for (int k = 21; k <= 30; ++k)
        s += std::to_string(k) + (k < 30 ? " " : "\n");
    return s;
}

double volume_sign(const Mesh& rMesh, std::size_t Block, const int (&rTet)[4]) {
    const auto& conn = rMesh.Cells(Block).Conn();
    const auto& pts = rMesh.Points();
    double x[4][3];
    for (int c = 0; c < 4; ++c)
        for (int d = 0; d < 3; ++d)
            x[c][d] = detail::read_double(pts, static_cast<std::size_t>(detail::read_int(
                                                   conn, static_cast<std::size_t>(rTet[c]))) *
                                                       3 +
                                                   static_cast<std::size_t>(d));
    double a[3], b[3], c[3];
    for (int d = 0; d < 3; ++d) {
        a[d] = x[1][d] - x[0][d];
        b[d] = x[2][d] - x[0][d];
        c[d] = x[3][d] - x[0][d];
    }
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
           a[2] * (b[0] * c[1] - b[1] * c[0]);
}

}  // namespace

TEST(Z88, ReadsHex20AndTet10InMeshioOrder) {
    const fs::path dir = deck_dir();
    write(dir / "z88i1.txt", structure());
    const Mesh mesh = meshioplusplus::read_z88((dir / "z88i1.txt").string());
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    EXPECT_EQ(mesh.Cells(1).Type(), "tetra10");
    EXPECT_GT(volume_sign(mesh, 0, {0, 1, 3, 4}), 0.0);
    EXPECT_GT(volume_sign(mesh, 1, {0, 1, 2, 3}), 0.0);
    // meshio++ tetra10 slot 7 is the mid-edge 0-3: (2, 0, 0.5).
    const auto& conn = mesh.Cells(1).Conn();
    const std::size_t n7 = static_cast<std::size_t>(detail::read_int(conn, 7));
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), n7 * 3 + 2), 0.5);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.Points(), n7 * 3 + 0), 2.0);
    EXPECT_EQ(detail::read_int(mesh.CellData("z88:type", 0), 0), 10);
}

TEST(Z88, AttachesTheSiblingResults) {
    const fs::path dir = deck_dir();
    write(dir / "z88i1.txt", structure());
    std::string o2 =
        "output file Z88O2.TXT: displacements, computed by Z88R V15OS\n"
        "                       *************\nLoad case: 1\n"
        "node           U(1)              U(2)              U(3)\n\n";
    for (int k = 1; k <= 30; ++k)
        o2 += "    " + std::to_string(k) + "    +1.0000000E-02    +0.0000000E+00    " +
              std::to_string(-0.001 * k) + "\n";
    write(dir / "Z88O2.TXT", o2);
    write(dir / "z88o3.txt",
          "output file Z88O3.TXT : stresses\n\nelement # = 1     type =20-n hexahedron\n"
          " XX YY ZZ SIGXX SIGYY SIGZZ TAUXY TAUYZ TAUZX SIGV\n"
          "  +1.0E+00  +1.0E+00  +1.0E+00  +2.0E+00 +0 +0 +0 +0 +0 +4.0E+00\n"
          "  +1.0E+00  +1.0E+00  +1.0E+00  +4.0E+00 +0 +0 +0 +0 +0 +6.0E+00\n");
    const Mesh mesh = meshioplusplus::read_z88((dir / "z88o2.txt").string());
    ASSERT_TRUE(mesh.HasPointData("U"));
    const auto& u = mesh.PointData("U");
    EXPECT_DOUBLE_EQ(detail::read_double(u, 29 * 3 + 2), -0.03);
    const auto& sig = mesh.CellData("SIG", 0);
    EXPECT_DOUBLE_EQ(detail::read_double(sig, 0), 3.0);
    EXPECT_DOUBLE_EQ(detail::read_double(mesh.CellData("SIGV", 0), 0), 5.0);
    EXPECT_TRUE(std::isnan(detail::read_double(mesh.CellData("SIG", 1), 0)));
    const Mesh bare = meshioplusplus::read_z88((dir / "z88i1.txt").string(), false);
    EXPECT_FALSE(bare.HasPointData("U"));
}

TEST(Z88, WriterRoundTripsAndDispatchIsByName) {
    const fs::path dir = deck_dir();
    write(dir / "z88i1.txt", structure());
    const Mesh mesh = meshioplusplus::read_z88((dir / "z88i1.txt").string(), false);
    const fs::path out = dir / "out";
    fs::create_directories(out);
    meshioplusplus::write_z88((out / "z88i1.txt").string(), mesh, true);
    EXPECT_TRUE(fs::exists(out / "z88i2.txt"));
    const Mesh back = meshioplusplus::read_z88((out / "z88i1.txt").string());
    mt::expect_mesh_eq(mesh, back, 1e-15);
    EXPECT_TRUE(meshioplusplus::is_z88_filename("a/b/Z88I1.TXT"));
    EXPECT_FALSE(meshioplusplus::is_z88_filename("a/b/z88i1.txt.bak"));
    EXPECT_EQ(meshioplusplus::resolve_format((out / "z88i1.txt").string(), ""), "z88");
    const std::string file = (out / "z88i1.txt").string();
    const Mesh via_registry =
        meshioplusplus::registry_read(file, meshioplusplus::resolve_format(file, ""), {});
    EXPECT_EQ(via_registry.Cells(0).Type(), "hexahedron20");
}

TEST(Z88, Errors) {
    const fs::path dir = deck_dir();
    write(dir / "z88i1.txt", "3 1 1 3 0\n1 3 0 0 0\n1 99\n1\n");
    EXPECT_THROW(meshioplusplus::read_z88((dir / "z88i1.txt").string()), ReadError);
    write(dir / "z88i1.txt", "3 2 1 6 0\n1 3 0 0 0\n2 3 1 0 0\n1 17\n1 2\n");
    EXPECT_THROW(meshioplusplus::read_z88((dir / "z88i1.txt").string()), ReadError);
    fs::remove(dir / "z88i1.txt");
    write(dir / "z88o2.txt", "1 0 0 0\n");
    EXPECT_THROW(meshioplusplus::read_z88((dir / "z88o2.txt").string()), ReadError);
}
