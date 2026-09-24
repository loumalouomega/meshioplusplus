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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/z88.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
namespace detail = meshioplusplus::detail;
namespace fs = std::filesystem;

// A fresh directory per call: Z88's file names are fixed, so two tests (or two
// ctest processes, whose `mt::temp_path` counters both start at 0) must never
// share one, or a sibling z88o2.txt left behind is read as this deck's results.
fs::path deck_dir() {
    static std::atomic<unsigned> counter{0};
    const std::string test = ::testing::UnitTest::GetInstance()->current_test_info()->name();
    const fs::path dir = fs::temp_directory_path() /
                         ("meshio_z88_" + test + "_" + std::to_string(std::random_device{}()) +
                          "_" + std::to_string(counter++));
    fs::remove_all(dir);
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

TEST(Z88, ReadsAndWritesTheWholeDeckWithBeamResults) {
    // Two 3-D beams (type 2) and a truss (type 4) along x, clamped at node 1,
    // with every input file and a hand-written z88o3/z88o4 in Z88R's layout.
    const fs::path dir = deck_dir();
    write(dir / "z88i1.txt",
          "3 3 3 18 0\n1 6 0 0 0\n2 6 1 0 0\n3 6 2 0 0\n1 2\n1 2\n2 2\n2 3\n3 4\n1 3\n");
    write(dir / "z88i2.txt", "7\n1 1 2 0\n1 2 2 0\n1 3 2 0\n1 4 2 0\n1 5 2 0\n1 6 2 0\n3 3 1 -5\n");
    write(dir / "z88mat.txt", "2\n1 2 51.txt\n3 3 steel.txt\n");
    write(dir / "51.txt", "210000 0.3\n");
    write(dir / "steel.txt", "200000 0.29\n");
    write(dir / "z88elp.txt", "1\n1 3 100 833 5 833 5 1406 208\n");
    write(dir / "z88int.txt", "1\n1 3 0 0\n");
    write(dir / "z88o3.txt",
          "output file Z88O3.TXT : stresses\n\n"
          "element #=     1 type =3Dbeam   node      1             node      2\n"
          "    SIGXX       TAUXX       SIGZZ1      SIGYY1      SIGZZ2      SIGYY2\n"
          " +1.000E+00  +2.000E+00  +3.000E+00  +4.000E+00  +5.000E+00  +6.000E+00\n\n"
          "element #=     2 type =3Dbeam   node      2             node      3\n"
          "    SIGXX       TAUXX       SIGZZ1      SIGYY1      SIGZZ2      SIGYY2\n"
          " -1.000E+00  +0.000E+00  +0.000E+00  +0.000E+00  +0.000E+00  +0.000E+00\n\n"
          "element # =     3     type = truss in space     SIG =  +7.500E+00\n");
    write(dir / "z88o4.txt",
          "output file Z88O4.TXT : nodal forces\n\nthe nodal sums for each node\n"
          "----------------------------\n  node       F(1)  F(2)  F(3)  F(4)  F(5)  F(6)\n\n"
          "    1 +0.0E+00 +0.0E+00 +5.0E+00 +0.0E+00 +1.0E+01 +0.0E+00\n"
          "    2 +0.0E+00 +0.0E+00 +0.0E+00 +0.0E+00 +0.0E+00 +0.0E+00\n"
          "    3 +0.0E+00 +0.0E+00 -5.0E+00 +0.0E+00 +0.0E+00 +0.0E+00\n");
    const Mesh mesh = meshioplusplus::read_z88((dir / "z88i1.txt").string());
    const auto& sigxx = mesh.CellData("SIGXX", 0);
    EXPECT_EQ(detail::read_double(sigxx, 0), 1.0);
    EXPECT_EQ(detail::read_double(sigxx, 2), 7.5);  // the truss's axial stress
    EXPECT_TRUE(std::isnan(detail::read_double(mesh.CellData("TAUXX", 0), 2)));
    EXPECT_FALSE(mesh.HasCellData("SIG"));
    EXPECT_EQ(detail::read_int(mesh.CellData("z88:material", 0), 2), 2);  // "steel.txt": row 2
    EXPECT_EQ(detail::read_double(mesh.CellData("z88:E", 0), 2), 200000.0);
    const auto& elp = mesh.CellData("z88:elp", 0);
    EXPECT_EQ(detail::read_double(elp, 12 + 5), 1406.0);
    EXPECT_TRUE(std::isnan(detail::read_double(elp, 12 + 7)));
    const auto& u = mesh.PointData("z88:bc:u");
    const auto& f = mesh.PointData("z88:bc:f");
    EXPECT_EQ(detail::read_double(u, 0), 0.0);
    EXPECT_TRUE(std::isnan(detail::read_double(u, 6)));
    EXPECT_EQ(detail::read_double(f, 2 * 6 + 2), -5.0);
    EXPECT_EQ(detail::read_double(mesh.PointData("F"), 4), 10.0);

    // Written back: the same inputs; the material files numbered 51 and 2.
    const fs::path out = deck_dir();
    meshioplusplus::write_z88((out / "z88i1.txt").string(), mesh);
    for (const char* name :
         {"z88i2.txt", "z88mat.txt", "51.txt", "2.txt", "z88elp.txt", "z88int.txt"})
        EXPECT_TRUE(fs::exists(out / name)) << name;
    const Mesh back = meshioplusplus::read_z88((out / "z88i1.txt").string());
    for (const char* name : {"z88:material", "z88:E", "z88:nu", "z88:elp", "z88:int"})
        for (std::size_t k = 0; k < mesh.CellData(name, 0).Size(); ++k) {
            const double a = detail::read_double(mesh.CellData(name, 0), k);
            const double b = detail::read_double(back.CellData(name, 0), k);
            EXPECT_TRUE(a == b || (std::isnan(a) && std::isnan(b))) << name << " " << k;
        }
    for (const char* name : {"z88:bc:u", "z88:bc:f"})
        for (std::size_t k = 0; k < mesh.PointData(name).Size(); ++k) {
            const double a = detail::read_double(mesh.PointData(name), k);
            const double b = detail::read_double(back.PointData(name), k);
            EXPECT_TRUE(a == b || (std::isnan(a) && std::isnan(b))) << name << " " << k;
        }
    fs::remove_all(dir);
    fs::remove_all(out);
}

TEST(Z88, AuroraSetsBecomeRegions) {
    const fs::path dir = deck_dir();
    write(dir / "z88structure.txt",
          "3 4 1 12 0 #AURORA_V2\n1 3 0 0 0\n2 3 1 0 0\n3 3 0 1 0\n"
          "4 3 0 0 1\n1 17\n1 2 3 4\n");
    write(dir / "z88sets.txt",
          "3\n#ELEMENTS MATERIAL 1 1 \"Steel part\"\n 1\n"
          "#SURFACE CONSTRAINT 2 1 \"Set2\"\n 1 0 9 9 9 9\n"
          "#NODES CONSTRAINT 7 2 \"FIX\"\n 1 4\n");
    const Mesh mesh = meshioplusplus::read_z88((dir / "z88structure.txt").string());
    ASSERT_EQ(mesh.NumRegions(), 2u);
    const std::size_t cell = mesh.FindRegion("Steel part", meshioplusplus::RegionKind::Cell);
    ASSERT_NE(cell, Mesh::npos);
    EXPECT_EQ(mesh.Region(cell).mTag, 1);
    const std::size_t fix = mesh.FindRegion("FIX", meshioplusplus::RegionKind::Point);
    ASSERT_NE(fix, Mesh::npos);
    EXPECT_EQ(mesh.Region(fix).mTag, 7);
    EXPECT_EQ(mesh.Region(fix).NumEntries(), 2u);
    EXPECT_EQ(mesh.Region(fix).Entries()[1], 3);
    fs::remove_all(dir);
}

TEST(Z88, Type19LagrangeQuadSurfaceLoadsAndSets) {
    const fs::path dir = deck_dir();
    // One 16-node plate on [0, 3]^2: Z88 node 4i + j + 1 at (i, j).
    std::string s = "2 16 1 48 0\n";
    for (int n = 0; n < 16; ++n)
        s += std::to_string(n + 1) + " 3 " + std::to_string(n / 4) + " " + std::to_string(n % 4) +
             "\n";
    s += "1 19\n";
    for (int n = 0; n < 16; ++n)
        s += std::to_string(n + 1) + (n == 15 ? "\n" : " ");
    write(dir / "z88i1.txt", s);
    write(dir / "z88i5.txt", "2\n1 -2.5\n7 1.0\n");  // element 7 does not exist
    Mesh mesh = meshioplusplus::read_z88((dir / "z88i1.txt").string());
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "VTK_LAGRANGE_QUADRILATERAL");
    ASSERT_EQ(mesh.Cells(0).NodesPerCell(), 16u);
    // VTK's order-3 Lagrange quad: corners, edges along r, s, r, s, interior.
    const int lattice[16][2] = {{0, 0}, {3, 0}, {3, 3}, {0, 3}, {1, 0}, {2, 0}, {3, 1}, {3, 2},
                                {1, 3}, {2, 3}, {0, 1}, {0, 2}, {1, 1}, {2, 1}, {1, 2}, {2, 2}};
    const auto& conn = mesh.Cells(0).Conn();
    const std::size_t dim = mesh.PointDim();
    for (int k = 0; k < 16; ++k) {
        const auto p =
            static_cast<std::size_t>(detail::read_int(conn, static_cast<std::size_t>(k)));
        EXPECT_EQ(detail::read_double(mesh.Points(), dim * p), lattice[k][0]) << k;
        EXPECT_EQ(detail::read_double(mesh.Points(), dim * p + 1), lattice[k][1]) << k;
    }
    const auto& load = mesh.FieldData("z88:surface_load");
    ASSERT_EQ(load.Shape(), (std::vector<std::size_t>{1, 3}));
    EXPECT_EQ(detail::read_double(load, 0), -2.5);
    EXPECT_TRUE(std::isnan(detail::read_double(load, 1)));
    const auto& refs = mesh.FieldData("z88:surface_load:cells");
    EXPECT_EQ(detail::read_int(refs, 0), 0);
    EXPECT_EQ(detail::read_int(refs, 1), -1);

    meshioplusplus::NDArray cells(meshioplusplus::DType::Int64, {1});
    cells.As<std::int64_t>()[0] = 0;
    mesh.AddRegion(
        meshioplusplus::Region("plate", meshioplusplus::RegionKind::Cell, 2, 3, std::move(cells)));
    meshioplusplus::NDArray fix(meshioplusplus::DType::Int64, {2});
    fix.As<std::int64_t>()[0] = 0;
    fix.As<std::int64_t>()[1] = 3;
    mesh.AddRegion(
        meshioplusplus::Region("fix", meshioplusplus::RegionKind::Point, -1, 3, std::move(fix)));
    const fs::path out = dir / "out";
    fs::create_directories(out);
    meshioplusplus::write_z88((out / "z88i1.txt").string(), mesh);
    std::ifstream i5(out / "z88i5.txt");
    const std::string loads((std::istreambuf_iterator<char>(i5)), std::istreambuf_iterator<char>());
    EXPECT_EQ(loads, "1\n1 -2.5000000000000000E+00\n");
    std::ifstream sets_in(out / "z88sets.txt");
    const std::string sets((std::istreambuf_iterator<char>(sets_in)),
                           std::istreambuf_iterator<char>());
    // points first; the cell set, second to claim id 3, gets 1
    EXPECT_EQ(sets,
              "2\n#NODES CONSTRAINT 3 2 \"fix\"\n         1          4 \n"
              "#ELEMENTS MATERIAL 1 1 \"plate\"\n         1 \n");
    const Mesh back = meshioplusplus::read_z88((out / "z88i1.txt").string());
    ASSERT_EQ(back.Cells(0).Type(), "VTK_LAGRANGE_QUADRILATERAL");
    for (std::size_t k = 0; k < 16; ++k)
        EXPECT_EQ(detail::read_int(back.Cells(0).Conn(), k), detail::read_int(conn, k));
    EXPECT_EQ(back.NumRegions(), 2u);
    EXPECT_EQ(detail::read_double(back.FieldData("z88:surface_load"), 0), -2.5);
    fs::remove_all(dir);
}
