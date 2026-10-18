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
 * @file test_lsdyna_d3plot.cpp
 * @brief LS-DYNA d3plot reader: a one-shell family (built here word by word) in
 *        single precision of both byte orders and in double precision, its
 *        second state in `d3plot01`; displacement, stress, plastic strain,
 *        deletion, steps, metadata, file-name resolution and refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/lsdyna_d3plot.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

// Words of one file: ints and floats in a word size and byte order.
struct Words {
    int mWs;
    bool mBig;
    std::string mBytes;

    void Put(std::uint64_t Bits) {
        const int n = mWs;
        for (int k = 0; k < n; ++k) {
            const int shift = mBig ? 8 * (n - 1 - k) : 8 * k;
            mBytes.push_back(static_cast<char>((Bits >> shift) & 0xff));
        }
    }
    void Int(std::int64_t V) {
        Put(mWs == 4 ? static_cast<std::uint32_t>(static_cast<std::int32_t>(V))
                     : static_cast<std::uint64_t>(V));
    }
    void Real(double V) {
        if (mWs == 4)
            Put(std::bit_cast<std::uint32_t>(static_cast<float>(V)));
        else
            Put(std::bit_cast<std::uint64_t>(V));
    }
};

const double kCoords[4][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};

// One shell (a quad on nodes 1..4, part 1), one layer of stress and plastic
// strain, element deletion; state k moves every node by 0.1*k in z, has
// stress component c = 10*k + c and plastic strain 0.5*k, and deletes the
// shell at k = 1.
void state(Words& rW, int K) {
    rW.Real(K * 1.0e-3);             // time
    for (const auto& c : kCoords) {  // current coordinates
        rW.Real(c[0]);
        rW.Real(c[1]);
        rW.Real(c[2] + 0.1 * K);
    }
    for (int c = 0; c < 6; ++c)
        rW.Real(10.0 * K + c);
    rW.Real(0.5 * K);
    rW.Real(K == 0 ? 1.0 : 0.0);  // deletion: the material number while alive
}

// The family as (base file bytes, d3plot01 bytes).
std::pair<std::string, std::string> family(int Ws, bool Big) {
    Words w{Ws, Big, {}};
    for (int i = 0; i < 64; ++i) {
        std::int64_t v = 0;
        switch (i) {
            case 11:
                v = 1;
                break;  // filetype d3plot
            case 15:
                v = 4;
                break;  // NDIM
            case 16:
                v = 4;
                break;  // NUMNP
            case 17:
                v = 6;
                break;  // ICODE
            case 20:
                v = 1;
                break;  // IU
            case 31:
                v = 1;
                break;  // NEL4
            case 32:
                v = 1;
                break;  // NUMMAT4
            case 33:
                v = 7;
                break;  // NV2D: 1 layer x (6 stress + 1 eps)
            case 36:
                v = -10001;
                break;  // MAXINT: 1 layer, element deletion
            case 43:
                v = 1000;
                break;  // IOSHL1
            case 44:
                v = 1000;
                break;  // IOSHL2
            case 51:
                v = 1;
                break;  // NMMAT
            default:
                break;
        }
        if (i < 10)
            w.Put(0x2020202020202020ULL);  // blank title
        else if (i == 14)
            w.Real(971.0);  // version
        else
            w.Int(v);
    }
    for (const auto& c : kCoords)
        for (double x : c)
            w.Real(x);
    for (int n = 1; n <= 4; ++n)
        w.Int(n);
    w.Int(1);  // part
    w.Real(-999999.0);
    state(w, 0);
    Words second{Ws, Big, {}};
    state(second, 1);
    second.Real(-999999.0);  // the end mark: states are counted up to the last non-zero word
    return {w.mBytes, second.mBytes};
}

std::string write_family(const std::pair<std::string, std::string>& rFiles) {
    const std::filesystem::path dir(
        mt::temp_path("_d3plot_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    std::ofstream((dir / "d3plot").string(), std::ios::binary) << rFiles.first;
    std::ofstream((dir / "d3plot01").string(), std::ios::binary) << rFiles.second;
    return (dir / "d3plot").string();
}

ReadOptions step(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

void expect_state(const Mesh& rMesh, int K) {
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "quad");
    ASSERT_EQ(rMesh.NumPoints(), 4u);
    const double* d = rMesh.PointData("displacement").As<double>();
    for (int n = 0; n < 4; ++n)
        EXPECT_NEAR(d[3 * n + 2], 0.1 * K, 1e-6);
    const double* s = rMesh.CellData("stress", 0).As<double>();
    for (int c = 0; c < 6; ++c)
        EXPECT_NEAR(s[c], 10.0 * K + c, 1e-5);
    EXPECT_NEAR(rMesh.CellData("effective_plastic_strain", 0).As<double>()[0], 0.5 * K, 1e-6);
    EXPECT_EQ(rMesh.CellData("lsdyna:alive", 0).As<std::int8_t>()[0], K == 0 ? 1 : 0);
    EXPECT_EQ(rMesh.CellData("lsdyna:part", 0).As<std::int64_t>()[0], 1);
    EXPECT_NEAR(rMesh.FieldData("meshio:time").As<double>()[0], K * 1.0e-3, 1e-9);
}

// One shell on nodes 1-4 (part 1), SPH particles on nodes 5 and 6 (material
// 2), an airbag (id 77) of two particles, a one-segment rigid road (id 9) and a
// reduced rigid body: NDIM 9, two states (the twin of the Python suite's
// _particles_family).
std::string particles_family() {
    Words w{4, false, {}};
    for (int i = 0; i < 64; ++i) {
        std::int64_t v = 0;
        switch (i) {
            case 11:
                v = 1;
                break;
            case 15:
                v = 9;  // NDIM: rigid road and reduced rigid bodies
                break;
            case 16:
                v = 6;
                break;
            case 17:
                v = 6;
                break;
            case 20:
                v = 1;
                break;
            case 31:
                v = 1;
                break;
            case 32:
                v = 1;
                break;
            case 33:
                v = 7;
                break;
            case 36:
                v = -10001;
                break;
            case 37:
                v = 2;  // NMSPH
                break;
            case 43:
            case 44:
                v = 1000;
                break;
            case 51:
                v = 2;
                break;
            case 54:
                v = 1;  // NPEFG: one airbag
                break;
            default:
                break;
        }
        if (i < 10)
            w.Put(0x20202020U);
        else if (i == 14)
            w.Real(971.0);
        else
            w.Int(v);
    }
    for (int v : {11, 1, 1, 6, 1, 1, 1, 1, 0, 1, 0})  // SPH flags
        w.Int(v);
    for (int v : {4, 5, 2, 2})  // airbag: geometry, particle, bag variables, particles
        w.Int(v);
    const char* names[] = {"Start", "Npart", "Bag ID", "NGas",    "GasC ID", "Pos x",
                           "Pos y", "Pos z", "Mass",   "Act Gas", "Bag Vol"};
    for (int t : {1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2})
        w.Int(t);
    for (const char* name : names) {
        std::string n = name;
        n.resize(8, ' ');
        for (char c : n)
            w.Int(c);
    }
    const double coords[6][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {5, 0, 0}, {6, 0, 0}};
    for (const auto& c : coords)
        for (double x : c)
            w.Real(x);
    for (int v : {1, 2, 3, 4, 1})  // the shell
        w.Int(v);
    for (int v : {1, 2, 1, 5, 0})  // a rigid body: part 2, node 5, no active node
        w.Int(v);
    for (int v : {5, 2, 6, 2})  // SPH (node, material)
        w.Int(v);
    for (int v : {1, 2, 77, 1})  // airbag 77: particles 1..2
        w.Int(v);
    for (int v : {4, 1, 1, 1, 101, 102, 103, 104})  // road nodes
        w.Int(v);
    const double road[4][3] = {{0, 0, -1}, {2, 0, -1}, {2, 2, -1}, {0, 2, -1}};
    for (const auto& c : road)
        for (double x : c)
            w.Real(x);
    for (int v : {9, 1, 101, 102, 103, 104})  // road 9, one segment
        w.Int(v);
    w.Real(-999999.0);
    for (int k = 0; k < 2; ++k) {
        w.Real(0.5 * k);
        for (const auto& c : coords) {
            w.Real(c[0]);
            w.Real(c[1]);
            w.Real(c[2] + 0.1 * k);
        }
        for (int c = 0; c < 6; ++c)
            w.Real(10.0 * k + c);
        w.Real(0.25 * k);
        w.Real(1.0);  // the shell lives
        for (int p = 0; p < 2; ++p) {
            w.Real(k == 1 && p == 1 ? -2.0 : 2.0);  // material, negative once deleted
            for (int v = 0; v < 13; ++v)
                w.Real(100.0 * k + 10.0 * p + v);
        }
        w.Int(k == 0 ? 2 : 1);  // airbag: active gas (an integer), volume
        w.Real(3.0 + k);
        for (int p = 0; p < 2; ++p) {  // particles: gas id, x, y, z, mass
            w.Int(1);
            w.Real(10.0 + p);
            w.Real(20.0 + k);
            w.Real(30.0);
            w.Real(0.5 + p);
        }
        for (double v : {0.0, 0.0, -0.1 * k, 0.0, 0.0, -1.0})  // road
            w.Real(v);
        for (int v = 0; v < 12; ++v)  // rigid body: coordinates, rotation
            w.Real(k + v / 10.0);
    }
    return w.mBytes;
}

}  // namespace

TEST(LsdynaD3plot, SphAirbagsRoadsAndRigidBodies) {
    const std::filesystem::path dir(
        mt::temp_path("_d3plot_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "d3plot").string();
    std::ofstream(path, std::ios::binary) << particles_family();
    const Mesh mesh = meshioplusplus::read_lsdyna_d3plot(path, step(1));
    ASSERT_EQ(mesh.NumCellBlocks(), 4u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad");
    EXPECT_EQ(mesh.Cells(1).Type(), "vertex");  // SPH
    EXPECT_EQ(mesh.Cells(2).Type(), "vertex");  // airbag particles
    EXPECT_EQ(mesh.Cells(3).Type(), "quad");    // the road
    ASSERT_EQ(mesh.NumPoints(), 12u);
    const double* p = mesh.Points().As<double>();
    EXPECT_NEAR(p[3 * 6], 10.0, 1e-6);  // particle 1 at x = 10, y = 21
    EXPECT_NEAR(p[3 * 6 + 1], 21.0, 1e-6);
    EXPECT_NEAR(p[3 * 8 + 2], -1.0, 1e-6);  // the road's first node
    EXPECT_TRUE(std::isnan(mesh.PointData("displacement").As<double>()[3 * 6]));
    EXPECT_EQ(mesh.CellData("lsdyna:alive", 1).As<std::int8_t>()[1], 0);  // SPH 2 deleted
    EXPECT_EQ(mesh.CellData("lsdyna:part", 2).As<std::int64_t>()[0], 77);
    EXPECT_EQ(mesh.CellData("lsdyna:part", 3).As<std::int64_t>()[0], 9);
    EXPECT_NEAR(mesh.CellData("sph_radius", 1).As<double>()[1], 110.0, 1e-4);
    EXPECT_EQ(mesh.CellData("airbag_gasc_id", 2).As<double>()[1], 1.0);
    EXPECT_NEAR(mesh.CellData("airbag_mass", 2).As<double>()[1], 1.5, 1e-6);
    EXPECT_EQ(mesh.FieldData("lsdyna:airbag:act_gas").As<double>()[0], 1.0);
    EXPECT_NEAR(mesh.FieldData("lsdyna:road_velocity").As<double>()[2], -1.0, 1e-6);
    EXPECT_NEAR(mesh.FieldData("lsdyna:rigid_body_coordinates").As<double>()[1], 1.1, 1e-6);
    EXPECT_EQ(mesh.FieldData("lsdyna:rigid_body_part").As<std::int64_t>()[0], 2);
    EXPECT_FALSE(mesh.HasFieldData("lsdyna:rigid_body_velocity"));  // reduced
}

TEST(LsdynaD3plot, BothPrecisionsAndByteOrdersReadTheSame) {
    for (const auto& [ws, big] : {std::pair{4, false}, std::pair{4, true}, std::pair{8, false}}) {
        const std::string path = write_family(family(ws, big));
        expect_state(meshioplusplus::read_lsdyna_d3plot(path), 0);
        expect_state(meshioplusplus::read_lsdyna_d3plot(path, step(-1)), 1);
        const auto times = meshioplusplus::lsdyna_d3plot_time_values(path);
        ASSERT_EQ(times.size(), 2u);
        EXPECT_NEAR(times[1], 1.0e-3, 1e-9);
    }
}

TEST(LsdynaD3plot, FoundByNameAndContentWithMetadata) {
    const std::string path = write_family(family(4, false));
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "lsdyna_d3plot");
    EXPECT_EQ(meshioplusplus::sniff_format(path), "lsdyna_d3plot");
    EXPECT_TRUE(meshioplusplus::is_d3plot_filename(path));
    const auto meta = meshioplusplus::registry_read_metadata(path, "lsdyna_d3plot", ReadOptions{});
    EXPECT_EQ(meta.mTimeValues.size(), 2u);
    const auto files = family(8, true);
    EXPECT_TRUE(meshioplusplus::is_d3plot_head(files.first.data(), files.first.size()));
    const std::string text = "*KEYWORD\n*NODE\n       1       0.0       0.0       0.0\n";
    EXPECT_FALSE(meshioplusplus::is_d3plot_head(text.data(), text.size()));
    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"displacement"};
    const Mesh narrowed = meshioplusplus::read_lsdyna_d3plot(path, only);
    EXPECT_TRUE(narrowed.HasPointData("displacement"));
    EXPECT_FALSE(narrowed.HasCellData("stress"));
}

TEST(LsdynaD3plot, MembersTruncationAndOutOfRangeAreRefused) {
    const std::string path = write_family(family(4, false));
    const std::string member = path + "01";
    EXPECT_EQ(meshioplusplus::resolve_format(member, ""), "lsdyna_d3plot");
    EXPECT_THROW(meshioplusplus::read_lsdyna_d3plot(member), ReadError);
    EXPECT_THROW(meshioplusplus::read_lsdyna_d3plot(path, step(2)), ReadError);
    auto files = family(4, false);
    files.first = files.first.substr(0, 300);
    EXPECT_THROW(meshioplusplus::read_lsdyna_d3plot(write_family(files)), ReadError);
}

TEST(LsdynaD3plotWriter, StateSeriesAndSplitRoundTrip) {
    Mesh mesh = mt::hex_mesh();
    meshioplusplus::NDArray u(meshioplusplus::DType::Float64, {mesh.NumPoints(), std::size_t{3}});
    for (std::size_t i = 0; i < u.Size(); ++i)
        u.As<double>()[i] = 0.25 * static_cast<double>(i);
    mesh.AddPointData("displacement", std::move(u));
    meshioplusplus::NDArray s(meshioplusplus::DType::Float64, {std::size_t{1}, std::size_t{6}});
    for (std::size_t i = 0; i < 6; ++i)
        s.As<double>()[i] = 1.5 * static_cast<double>(i);
    mesh.AddCellData("stress", {std::move(s)});
    meshioplusplus::NDArray t(meshioplusplus::DType::Float64, {std::size_t{1}});
    t.As<double>()[0] = 0.5;
    mesh.AddFieldData("meshio:time", std::move(t));

    const std::string dir = mt::temp_path("_d3plot_writer");
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/d3plot";
    meshioplusplus::write_lsdyna_d3plot(path, mesh);
    const Mesh back = meshioplusplus::read_lsdyna_d3plot(path);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    for (std::size_t i = 0; i < 3 * mesh.NumPoints(); ++i)
        EXPECT_NEAR(back.PointData("displacement").As<double>()[i], 0.25 * i, 1e-5);
    for (int c = 0; c < 6; ++c)
        EXPECT_NEAR(back.CellData("stress", 0).As<double>()[c], 1.5 * c, 1e-5);
    EXPECT_NEAR(back.FieldData("meshio:time").As<double>()[0], 0.5, 1e-9);

    // A family split by size: a stale member beside it is replaced.
    std::ofstream(path + "07", std::ios::binary) << "stale";
    meshioplusplus::LsdynaD3plotSeriesWriter writer(path, 600);
    for (int k = 0; k < 4; ++k)
        writer.Write(0.1 * k, mesh);
    writer.Finalize();
    EXPECT_EQ(writer.NumSteps(), 4u);
    EXPECT_TRUE(std::filesystem::exists(path + "01"));
    EXPECT_FALSE(std::filesystem::exists(path + "07"));
    EXPECT_EQ(meshioplusplus::lsdyna_d3plot_time_values(path).size(), 4u);

    meshioplusplus::LsdynaD3plotSeriesWriter mismatch(dir + "/other");
    mismatch.Write(0.0, mesh);
    EXPECT_THROW(mismatch.Write(1.0, mt::tet_mesh()), meshioplusplus::WriteError);
    EXPECT_THROW(meshioplusplus::write_lsdyna_d3plot(dir + "/none", Mesh{}),
                 meshioplusplus::WriteError);
}

// The airbag particles (`Airbag <id>` vertices) and rigid road segments (`Rigid road <id>` quads)
// of a read family are written back into their own sections, not as nodes and SPH particles.
TEST(LsdynaD3plotWriter, AirbagsAndRoadsRoundTrip) {
    const std::filesystem::path dir(
        mt::temp_path("_d3plot_airbag_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    const std::string source = (dir / "d3plot").string();
    std::ofstream(source, std::ios::binary) << particles_family();
    const Mesh mesh = meshioplusplus::read_lsdyna_d3plot(source, step(1));
    ASSERT_EQ(mesh.NumPoints(), 12u);

    std::filesystem::create_directories(dir / "out");
    const std::string out = (dir / "out" / "d3plot").string();
    meshioplusplus::write_lsdyna_d3plot(out, mesh);
    const Mesh back = meshioplusplus::read_lsdyna_d3plot(out);
    ASSERT_EQ(back.NumCellBlocks(), 4u);
    EXPECT_EQ(back.Cells(2).Type(), "vertex");  // the airbag particles
    EXPECT_EQ(back.Cells(3).Type(), "quad");    // the road
    ASSERT_EQ(back.NumPoints(), 12u);
    for (std::size_t i = 0; i < 12 * 3; ++i)
        EXPECT_NEAR(back.Points().As<double>()[i], mesh.Points().As<double>()[i], 1e-6);
    EXPECT_EQ(back.PointData("lsdyna:nid").As<std::int64_t>()[8], -1);
    EXPECT_EQ(back.CellData("lsdyna:part", 2).As<std::int64_t>()[0], 77);
    EXPECT_EQ(back.CellData("lsdyna:part", 3).As<std::int64_t>()[0], 9);
    EXPECT_NEAR(back.CellData("airbag_mass", 2).As<double>()[1], 1.5, 1e-6);
    EXPECT_NEAR(back.FieldData("lsdyna:airbag:bag_vol").As<double>()[0], 4.0, 1e-6);
    EXPECT_NEAR(back.FieldData("lsdyna:road_velocity").As<double>()[2], -1.0, 1e-6);
    EXPECT_FALSE(back.HasFieldData("lsdyna:rigid_body_velocity"));  // NDIM 9: reduced
    EXPECT_TRUE(back.HasFieldData("lsdyna:rigid_body_coordinates"));
}

// The name decides the file type: a `d3part` is file type 5.
TEST(LsdynaD3plotWriter, ADpartCarriesItsFileType) {
    const Mesh mesh = mt::hex_mesh();
    const std::string dir = mt::temp_path("_d3part_writer");
    std::filesystem::create_directories(dir);
    meshioplusplus::write_lsdyna_d3plot(dir + "/d3part", mesh);
    meshioplusplus::write_lsdyna_d3plot(dir + "/d3plot", mesh);
    const auto file_type = [](const std::string& rPath) {
        std::ifstream in(rPath, std::ios::binary);
        char word[4];
        in.seekg(4 * 11);
        in.read(word, 4);
        std::int32_t v;
        std::memcpy(&v, word, 4);
        return v;
    };
    EXPECT_EQ(file_type(dir + "/d3part"), 5);
    EXPECT_EQ(file_type(dir + "/d3plot"), 1);
    EXPECT_EQ(meshioplusplus::read_lsdyna_d3plot(dir + "/d3part").NumCellBlocks(), 1u);
}

// Bricks the reader marks `lsdyna:thick_shell` are written as the thick-shell family.
TEST(LsdynaD3plotWriter, ThickShellsAreAFamilyOfTheirOwn) {
    Mesh mesh = mt::hex_mesh();
    meshioplusplus::NDArray mark(meshioplusplus::DType::Int8, {std::size_t{1}});
    mark.As<std::int8_t>()[0] = 1;
    mesh.AddCellData("lsdyna:thick_shell", {std::move(mark)});
    meshioplusplus::NDArray stress(meshioplusplus::DType::Float64,
                                   {std::size_t{1}, std::size_t{6}});
    for (std::size_t i = 0; i < 6; ++i)
        stress.As<double>()[i] = 0.5 * static_cast<double>(i + 1);
    mesh.AddCellData("stress", {std::move(stress)});
    const std::string dir = mt::temp_path("_d3plot_tshell");
    std::filesystem::create_directories(dir);
    meshioplusplus::write_lsdyna_d3plot(dir + "/d3plot", mesh);
    const Mesh back = meshioplusplus::read_lsdyna_d3plot(dir + "/d3plot");
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
    ASSERT_TRUE(back.HasCellData("lsdyna:thick_shell"));
    EXPECT_EQ(back.CellData("lsdyna:thick_shell", 0).As<std::int8_t>()[0], 1);
    for (int c = 0; c < 6; ++c)
        EXPECT_NEAR(back.CellData("stress", 0).As<double>()[c], 0.5 * (c + 1), 1e-5);
}
