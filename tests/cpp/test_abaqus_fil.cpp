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
 * @file test_abaqus_fil.cpp
 * @brief Abaqus `.fil` reader: one model and two increments written as ASCII and
 *        as little- and big-endian binary (built here record by record), node
 *        and element results, long set labels, steps and metadata. And the writer:
 *        a model read from those files written back in both layouts, results by
 *        node and element, series, refusals, and the registry and sequence
 *        driver.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <variant>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus_fil.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::RegionKind;
using meshioplusplus::WriteError;
namespace detail = meshioplusplus::detail;

// A record item: integer, real or 8-character text.
using Item = std::variant<long long, double, std::string>;
using Record = std::vector<Item>;  // key first

std::vector<Record> records() {
    std::vector<Record> r;
    r.push_back({1921LL, std::string("6.23-1"), std::string("24-Sep-2"), std::string("026"),
                 std::string("12:00:00"), 1LL, 8LL, 1.0});
    // A C3D8R on the unit cube, split over 1900 + 1990.
    r.push_back({1900LL, 1LL, std::string("C3D8R"), 1LL, 2LL, 3LL, 4LL});
    r.push_back({1990LL, 5LL, 6LL, 7LL, 8LL});
    const double c[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int n = 0; n < 8; ++n)
        r.push_back({1901LL, static_cast<long long>(n + 1), c[n][0], c[n][1], c[n][2]});
    // A node set named by a 1940 label longer than 8 characters.
    r.push_back({1931LL, std::string("1"), 1LL, 2LL});
    r.push_back({1932LL, 3LL, 4LL});
    r.push_back(
        {1940LL, 1LL, std::string("ASSEMBLY"), std::string("_BOTTOM_"), std::string("NODES")});
    for (int inc = 1; inc <= 2; ++inc) {
        r.push_back({2000LL, 0.5 * inc, 0.5 * inc, 0.0, 0.0, 1LL, 1LL, static_cast<long long>(inc),
                     0LL, 0.0, 0.0, 0.5});
        r.push_back({1911LL, 0LL, std::string(""), std::string("")});
        r.push_back({1LL, 1LL, 1LL, 0LL, 0LL, std::string(""), 3LL, 3LL, 0LL, 0LL});
        Record s{11LL};
        for (int k = 0; k < 6; ++k)
            s.push_back(100.0 * inc + k);
        r.push_back(s);
        r.push_back({1911LL, 1LL, std::string("")});
        for (int n = 1; n <= 8; ++n)
            r.push_back({101LL, static_cast<long long>(n), 0.001 * n + inc, 0.0, -0.5 * inc});
        r.push_back({2001LL});
    }
    return r;
}

std::string ascii(const std::vector<Record>& rRecords) {
    std::string flat;
    for (const Record& rec : rRecords) {
        flat += '*';
        Record full{static_cast<long long>(rec.size() + 1)};
        full.insert(full.end(), rec.begin(), rec.end());
        for (const Item& it : full) {
            char buf[64];
            if (const long long* i = std::get_if<long long>(&it)) {
                const std::string digits = std::to_string(*i);
                std::snprintf(buf, sizeof(buf), "I%2zu%s", digits.size(), digits.c_str());
            } else if (const double* d = std::get_if<double>(&it)) {
                std::snprintf(buf, sizeof(buf), "D% .15E", *d);
                for (char* q = buf; *q; ++q)
                    if (*q == 'E' && q != buf)
                        *q = 'D';
                buf[0] = 'D';
            } else {
                std::string a = std::get<std::string>(it);
                a.resize(8, ' ');
                std::snprintf(buf, sizeof(buf), "A%s", a.c_str());
            }
            flat += buf;
        }
    }
    std::string out;
    for (std::size_t k = 0; k < flat.size(); k += 80)
        out += flat.substr(k, 80) + "\n";
    return out;
}

std::string ascii() {
    return ascii(records());
}

// The model of records() with what no public file carries: a contact surface
// on the hexahedron's face S1, an eigenfrequency step of two modes, then a
// modal dynamic increment with generalized displacements, total energies, a
// rebar's stress and contact tractions at node 3.
std::vector<Record> extra_records() {
    std::vector<Record> r;
    for (const Record& rec : records()) {
        const long long key = std::get<long long>(rec[0]);
        if (key == 2000)
            break;
        r.push_back(rec);
    }
    r.push_back({1501LL, std::string("CSURF"), 3LL, 1LL, 1LL, 0LL});
    r.push_back({1502LL, 1LL, 1LL, 4LL, 1LL, 2LL, 3LL, 4LL});
    r.push_back({2000LL, 0.0, 0.0, 0.0, 0.0, 41LL, 1LL, 1LL, 0LL, 0.0, 0.0, 0.0});
    for (long long mode = 1; mode <= 2; ++mode) {
        r.push_back({1980LL, mode, 50.0 * static_cast<double>(mode), 2.0, 0.0, 0.5, 0.25});
        r.push_back({1911LL, 1LL, std::string("")});
        for (long long n = 1; n <= 8; ++n)
            r.push_back({101LL, n, static_cast<double>(mode), 0.0, 0.0});
    }
    r.push_back({2001LL});
    r.push_back({2000LL, 1.0, 1.0, 0.0, 0.0, 92LL, 2LL, 1LL, 0LL, 0.0, 0.0, 0.1});
    r.push_back({301LL, 0.25, -0.5});
    r.push_back({1999LL, 1.0, 2.0, 3.0});
    r.push_back({1911LL, 0LL, std::string(""), std::string("")});
    r.push_back({1LL, 1LL, 1LL, 0LL, 3LL, std::string("RB1"), 1LL, 0LL, 0LL, 0LL});
    r.push_back({11LL, 123.0});
    r.push_back({1503LL, 0LL, std::string("CSURF"), std::string("MASTER"), std::string("")});
    r.push_back({1504LL, 3LL, 3LL});
    r.push_back({1511LL, 7.0, 0.5, 0.25});
    r.push_back({2001LL});
    return r;
}

void put_word(std::string& rOut, std::uint64_t Raw, bool Big) {
    for (int b = 0; b < 8; ++b) {
        const int shift = Big ? 8 * (7 - b) : 8 * b;
        rOut += static_cast<char>((Raw >> shift) & 0xff);
    }
}

std::string binary(bool Big) {
    std::string words;
    for (const Record& rec : records()) {
        Record full{static_cast<long long>(rec.size() + 1)};
        full.insert(full.end(), rec.begin(), rec.end());
        for (const Item& it : full) {
            if (const long long* i = std::get_if<long long>(&it)) {
                put_word(words, static_cast<std::uint64_t>(*i), Big);
            } else if (const double* d = std::get_if<double>(&it)) {
                std::uint64_t raw;
                std::memcpy(&raw, d, 8);
                put_word(words, raw, Big);
            } else {
                std::string a = std::get<std::string>(it);
                a.resize(8, ' ');
                words += a;
            }
        }
    }
    words.append((4096 - words.size() % 4096) % 4096, '\0');
    std::string out;
    for (std::size_t b = 0; b < words.size(); b += 4096) {
        std::string marker;
        for (int k = 0; k < 4; ++k)
            marker += static_cast<char>((4096u >> (Big ? 8 * (3 - k) : 8 * k)) & 0xff);
        out += marker + words.substr(b, 4096) + marker;
    }
    return out;
}

std::string write_file(const std::string& rBody) {
    // Tagged: every ctest process starts `mt::temp_path`'s counter at 0.
    const std::string path = mt::temp_path("_" + std::to_string(std::random_device{}()) + ".fil");
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

void expect_model(const Mesh& rMesh, int Increment) {
    ASSERT_EQ(rMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(rMesh.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(rMesh.NumPoints(), 8u);
    EXPECT_EQ(detail::read_double(rMesh.FieldData("meshio:time"), 0), 0.5 * Increment);
    EXPECT_EQ(detail::read_int(rMesh.FieldData("abaqus:increment"), 0), Increment);
    const auto& u = rMesh.PointData("U");
    EXPECT_DOUBLE_EQ(detail::read_double(u, 7 * 3), 0.008 + Increment);
    EXPECT_DOUBLE_EQ(detail::read_double(u, 7 * 3 + 2), -0.5 * Increment);
    // (cells, points * components), in the file's component order.
    const auto& s = rMesh.CellData("S", 0);
    ASSERT_EQ(s.Shape().size(), 2u);
    EXPECT_EQ(s.Shape()[1], 6u);
    EXPECT_DOUBLE_EQ(detail::read_double(s, 5), 100.0 * Increment + 5);
    const auto& layout = rMesh.FieldData("abaqus:layout:S");
    EXPECT_EQ(detail::read_int(layout, 0), 1);
    EXPECT_EQ(detail::read_int(layout, 1), 6);
    const std::size_t set = rMesh.FindRegion("ASSEMBLY_BOTTOM_NODES", RegionKind::Point);
    ASSERT_NE(set, Mesh::npos);
    EXPECT_EQ(rMesh.Region(set).NumEntries(), 4u);
}

ReadOptions step(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

}  // namespace

TEST(AbaqusFil, AsciiAndBothBinaryOrdersReadTheSame) {
    for (const std::string& body : {ascii(), binary(false), binary(true)}) {
        const std::string path = write_file(body);
        expect_model(meshioplusplus::read_abaqus_fil(path), 1);
        expect_model(meshioplusplus::read_abaqus_fil(path, step(-1)), 2);
        EXPECT_EQ(meshioplusplus::abaqus_fil_time_values(path), (std::vector<double>{0.5, 1.0}));
    }
}

TEST(AbaqusFil, StepsMetadataAndSelection) {
    const std::string path = write_file(binary(false));
    EXPECT_THROW(meshioplusplus::read_abaqus_fil(path, step(2)), ReadError);
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "abaqus_fil");
    const auto meta = meshioplusplus::registry_read_metadata(path, "abaqus_fil", ReadOptions{});
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.5, 1.0}));
    ReadOptions only_u;
    only_u.mDataArrays = std::vector<std::string>{"U"};
    const Mesh narrowed = meshioplusplus::read_abaqus_fil(path, only_u);
    EXPECT_TRUE(narrowed.HasPointData("U"));
    EXPECT_FALSE(narrowed.HasCellData("S"));
}

TEST(AbaqusFil, TruncatedFilesAreRefused) {
    const std::string bin = binary(true);
    EXPECT_THROW(meshioplusplus::read_abaqus_fil(write_file(bin.substr(0, 3000))), ReadError);
    const std::string text = ascii();
    EXPECT_THROW(meshioplusplus::read_abaqus_fil(write_file(text.substr(0, 200))), ReadError);
}

TEST(AbaqusFil, ModesContactEnergiesAndRebar) {
    const std::string path = write_file(ascii(extra_records()));
    EXPECT_EQ(meshioplusplus::abaqus_fil_time_values(path).size(), 3u);  // two modes, one increment
    const Mesh mode2 = meshioplusplus::read_abaqus_fil(path, step(1));
    EXPECT_EQ(detail::read_int(mode2.FieldData("abaqus:mode"), 0), 2);
    EXPECT_DOUBLE_EQ(detail::read_double(mode2.FieldData("abaqus:eigenvalue"), 0), 100.0);
    EXPECT_DOUBLE_EQ(detail::read_double(mode2.FieldData("abaqus:effective_mass"), 0), 0.25);
    EXPECT_DOUBLE_EQ(detail::read_double(mode2.PointData("U"), 0), 2.0);
    const std::size_t surf = mode2.FindRegion("CSURF", RegionKind::Side);
    ASSERT_NE(surf, Mesh::npos);
    EXPECT_EQ(mode2.Region(surf).Entries()[1], 4);  // S1 of a hexahedron: meshio++ facet 4
    const Mesh dyn = meshioplusplus::read_abaqus_fil(path, step(2));
    EXPECT_DOUBLE_EQ(detail::read_double(dyn.FieldData("abaqus:GU"), 1), -0.5);
    EXPECT_DOUBLE_EQ(detail::read_double(dyn.FieldData("abaqus:ALLWK"), 0), 3.0);
    EXPECT_DOUBLE_EQ(detail::read_double(dyn.CellData("S@rebar:RB1", 0), 0), 123.0);
    EXPECT_DOUBLE_EQ(detail::read_double(dyn.PointData("CSTRESS"), 2 * 3), 7.0);
}

// ---------------------------------------------------------------------------------------
// The writer
// ---------------------------------------------------------------------------------------

namespace {

NDArray reals(std::initializer_list<double> Values, std::vector<std::size_t> Shape) {
    NDArray a(meshioplusplus::DType::Float64, std::move(Shape));
    std::copy(Values.begin(), Values.end(), a.As<double>());
    return a;
}

NDArray ints(std::initializer_list<std::int64_t> Values) {
    NDArray a(meshioplusplus::DType::Int64, {Values.size()});
    std::copy(Values.begin(), Values.end(), a.As<std::int64_t>());
    return a;
}

std::string read_all(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string temp_fil(const std::string& rName = "written.fil") {
    const std::filesystem::path dir(
        mt::temp_path("_fil_w_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    return (dir / rName).string();
}

// A hexahedron with U per node, S at the centroid and a node set.
Mesh unit_hex(double Scale = 1.0, std::optional<double> Time = std::nullopt) {
    Mesh m = mt::hex_mesh();
    NDArray u(meshioplusplus::DType::Float64, {std::size_t{8}, std::size_t{3}});
    for (int i = 0; i < 24; ++i)
        u.As<double>()[i] = Scale * i / 8.0;
    m.AddPointData("U", std::move(u));
    m.AddCellData("S",
                  {reals({Scale, 2 * Scale, 3 * Scale, 4 * Scale, 5 * Scale, 6 * Scale}, {1, 6})});
    if (Time)
        m.AddFieldData("meshio:time", reals({*Time}, {1}));
    return m;
}

ReadOptions step_of(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

}  // namespace

TEST(AbaqusFilWriter, ResultsComeBackInBothLayouts) {
    for (const bool binary : {false, true}) {
        const std::string path = temp_fil();
        meshioplusplus::write_abaqus_fil(path, unit_hex(2.0, 0.5), binary);
        const std::string bytes = read_all(path);
        if (binary) {
            EXPECT_EQ(bytes.size() % 4104, 0u);  // 512 words of 8 bytes and two 4-byte markers
            EXPECT_NE(bytes[0], '*');
        } else {
            EXPECT_EQ(bytes[0], '*');
        }
        const Mesh back = meshioplusplus::read_abaqus_fil(path);
        ASSERT_EQ(back.NumPoints(), 8u);
        EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
        EXPECT_EQ(back.PointData("abaqus:id").As<std::int64_t>()[7], 8);
        for (int i = 0; i < 24; ++i)
            EXPECT_EQ(back.PointData("U").As<double>()[i], 2.0 * i / 8.0);
        for (int c = 0; c < 6; ++c)
            EXPECT_EQ(back.CellData("S", 0).As<double>()[c], 2.0 * (c + 1));
        EXPECT_EQ(back.FieldData("meshio:time").As<double>()[0], 0.5);
    }
}

TEST(AbaqusFilWriter, TheModelOfTheReaderFixtureSurvivesARewrite) {
    const std::string source = write_file(ascii());
    const Mesh first = meshioplusplus::read_abaqus_fil(source, step_of(1));
    for (const bool binary : {false, true}) {
        const std::string path = temp_fil();
        meshioplusplus::write_abaqus_fil(path, first, binary);
        const Mesh back = meshioplusplus::read_abaqus_fil(path);
        mt::expect_same_geometry(first, back);
        ASSERT_EQ(back.NumRegions(), 1u);
        EXPECT_EQ(back.Region(0).mName, "ASSEMBLY_BOTTOM_NODES");
        EXPECT_EQ(back.Region(0).mEntries.Size(), first.Region(0).mEntries.Size());
        for (int i = 0; i < 24; ++i)
            EXPECT_EQ(back.PointData("U").As<double>()[i], first.PointData("U").As<double>()[i]);
        EXPECT_EQ(back.CellData("S", 0).As<double>()[5], first.CellData("S", 0).As<double>()[5]);
        for (const char* key : {"abaqus:step", "abaqus:increment", "abaqus:procedure"})
            EXPECT_EQ(back.FieldData(key).As<std::int64_t>()[0],
                      first.FieldData(key).As<std::int64_t>()[0]);
    }
}

TEST(AbaqusFilWriter, ASeriesIsTheModelOnceThenAnIncrementPerStep) {
    for (const bool binary : {false, true}) {
        const std::string path = temp_fil();
        meshioplusplus::AbaqusFilSeriesWriter writer(path, binary);
        writer.Write(0.5, unit_hex(1.0));
        EXPECT_THROW(writer.Write(1.0, mt::tri_mesh()), WriteError);
        EXPECT_EQ(writer.NumSteps(), 1u);
        writer.Write(1.0, unit_hex(2.0));
        writer.Write(2.5, unit_hex(3.0));
        writer.Finalize();
        writer.Finalize();
        EXPECT_EQ(meshioplusplus::abaqus_fil_time_values(path),
                  (std::vector<double>{0.5, 1.0, 2.5}));
        const Mesh last = meshioplusplus::read_abaqus_fil(path, step_of(-1));
        EXPECT_EQ(last.PointData("U").As<double>()[23], 3.0 * 23 / 8.0);
        EXPECT_EQ(last.FieldData("abaqus:increment").As<std::int64_t>()[0], 3);
        meshioplusplus::AbaqusFilSeriesWriter empty(temp_fil(), binary);
        EXPECT_THROW(empty.Finalize(), WriteError);
    }
}

TEST(AbaqusFilWriter, ASeriesIsOneStepsBytesWhenItHoldsOne) {
    for (const bool binary : {false, true}) {
        const Mesh mesh = unit_hex(1.5, 0.5);
        const std::string plain = temp_fil(), series = temp_fil();
        meshioplusplus::write_abaqus_fil(plain, mesh, binary);
        meshioplusplus::AbaqusFilSeriesWriter writer(series, binary);
        writer.Write(0.5, mesh);
        writer.Finalize();
        EXPECT_EQ(read_all(plain), read_all(series));
    }
}

TEST(AbaqusFilWriter, IntegrationPointsSectionPointsAndNodesKeepTheirPlaces) {
    Mesh m = mt::tet_mesh();  // two tetrahedra over five points
    m.AddCellData("abaqus:id", {ints({7, 9})});
    m.AddCellData("S", {reals({1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                               13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                              {2, 12})});
    m.AddFieldData("abaqus:layout:S", ints({2, 6}));
    m.AddCellData("MISES@sp3", {reals({10.0, 20.0}, {2})});
    m.AddCellData("NFORC", {reals({1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                                   13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
                                  {2, 12})});
    m.AddFieldData("abaqus:layout:NFORC", ints({4, 3}));
    const std::string path = temp_fil();
    meshioplusplus::write_abaqus_fil(path, m);
    const Mesh back = meshioplusplus::read_abaqus_fil(path);
    EXPECT_EQ(back.CellData("abaqus:id", 0).As<std::int64_t>()[1], 9);
    for (int i = 0; i < 24; ++i) {
        EXPECT_EQ(back.CellData("S", 0).As<double>()[i], i + 1);
        EXPECT_EQ(back.CellData("NFORC", 0).As<double>()[i], i + 1);
    }
    EXPECT_EQ(back.FieldData("abaqus:layout:S").As<std::int64_t>()[0], 2);
    EXPECT_EQ(back.CellData("MISES@sp3", 0).As<double>()[1], 20.0);
}

TEST(AbaqusFilWriter, NaNAndRefusals) {
    Mesh nan = unit_hex();
    NDArray u(meshioplusplus::DType::Float64, {std::size_t{8}, std::size_t{3}});
    for (int i = 0; i < 24; ++i)
        u.As<double>()[i] = i;
    for (int c = 0; c < 3; ++c)
        u.As<double>()[3 + c] = std::nan("");
    u.As<double>()[8] = std::nan("");
    nan.AddPointData("U", std::move(u));
    const std::string path = temp_fil();
    meshioplusplus::write_abaqus_fil(path, nan);
    const Mesh back = meshioplusplus::read_abaqus_fil(path);
    EXPECT_TRUE(std::isnan(back.PointData("U").As<double>()[3]));
    EXPECT_EQ(back.PointData("U").As<double>()[6], 6.0);
    EXPECT_TRUE(std::isnan(back.PointData("U").As<double>()[8]));

    Mesh dup = unit_hex();
    dup.AddPointData("abaqus:id", ints({1, 2, 3, 4, 5, 6, 7, 7}));
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), dup), WriteError);
    Mesh zero = unit_hex();
    zero.AddPointData("abaqus:id", ints({0, 1, 2, 3, 4, 5, 6, 7}));
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), zero), WriteError);
    Mesh layout = unit_hex();
    layout.AddFieldData("abaqus:layout:S", ints({2, 6}));
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), layout), WriteError);
    Mesh huge = unit_hex();
    NDArray h(meshioplusplus::DType::Float64, {std::size_t{8}, std::size_t{3}});
    for (int i = 0; i < 24; ++i)
        h.As<double>()[i] = i == 0 ? 1e150 : 0.0;
    huge.AddPointData("U", std::move(h));
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), huge), WriteError);
    Mesh outside = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, "triangle", {{0, 1, 3}});
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), outside), WriteError);
    EXPECT_THROW(meshioplusplus::write_abaqus_fil("/no/such/dir/x.fil", unit_hex()), WriteError);
}

TEST(AbaqusFilWriter, TheRegistryAndTheSequenceDriverWriteIt) {
    EXPECT_EQ(meshioplusplus::resolve_write_format("/a/job.fil", ""), "abaqus_fil");
    const std::string path = temp_fil("by_name.bin");
    meshioplusplus::registry_writers().at("abaqus_fil")(path, unit_hex());
    EXPECT_EQ(read_all(path)[0], '*');
    EXPECT_EQ(meshioplusplus::sniff_format(path), "abaqus_fil");
    std::string why;
    EXPECT_TRUE(meshioplusplus::sequence_write_supports_time("abaqus_fil", why));
    meshioplusplus::WriteOptions options;
    options.mEncoding = meshioplusplus::WriteEncoding::Binary;
    const std::string bin = temp_fil("binary.fil");
    meshioplusplus::registry_write_ex(bin, unit_hex(), "abaqus_fil", options);
    EXPECT_NE(read_all(bin)[0], '*');
}

namespace {

Mesh contact_hex(bool Pair = true) {
    Mesh m = unit_hex();
    NDArray slave(meshioplusplus::DType::Int64, {std::size_t{1}, std::size_t{2}});
    slave.As<std::int64_t>()[0] = 0;
    slave.As<std::int64_t>()[1] = 0;
    NDArray master(meshioplusplus::DType::Int64, {std::size_t{1}, std::size_t{2}});
    master.As<std::int64_t>()[0] = 0;
    master.As<std::int64_t>()[1] = 1;
    m.AddRegion(meshioplusplus::Region("SLAVE", RegionKind::Side, std::move(slave)));
    if (Pair)
        m.AddRegion(
            meshioplusplus::Region("master_surface_long", RegionKind::Side, std::move(master)));
    NDArray cs(meshioplusplus::DType::Float64, {std::size_t{8}, std::size_t{3}});
    for (int i = 0; i < 24; ++i)
        cs.As<double>()[i] = i < 12 ? std::nan("") : i + 0.5;
    m.AddPointData("CSTRESS", std::move(cs));
    NDArray ar(meshioplusplus::DType::Float64, {std::size_t{8}});
    for (int i = 0; i < 8; ++i)
        ar.As<double>()[i] = i == 5 ? 2.5 : std::nan("");
    m.AddPointData("key_1593", std::move(ar));
    return m;
}

}  // namespace

TEST(AbaqusFilWriter, ContactOutputComesBackInBothLayouts) {
    for (const bool binary : {false, true}) {
        const std::string path = temp_fil();
        meshioplusplus::write_abaqus_fil(path, contact_hex(), binary);
        const Mesh back = meshioplusplus::read_abaqus_fil(path);
        ASSERT_TRUE(back.HasPointData("CSTRESS"));
        ASSERT_TRUE(back.HasPointData("key_1593"));
        for (int i = 0; i < 24; ++i) {
            const double got = back.PointData("CSTRESS").As<double>()[i];
            if (i < 12)
                EXPECT_TRUE(std::isnan(got));
            else
                EXPECT_EQ(got, i + 0.5);
        }
        EXPECT_EQ(back.PointData("key_1593").As<double>()[5], 2.5);
        EXPECT_TRUE(std::isnan(back.PointData("key_1593").As<double>()[4]));
        EXPECT_EQ(back.NumRegions(), 2u);
    }
}

TEST(AbaqusFilWriter, ContactOutputNeedsASurfacePair) {
    const Mesh lone = contact_hex(false);
    const std::string path = temp_fil();
    meshioplusplus::write_abaqus_fil(path, lone);
    EXPECT_FALSE(meshioplusplus::read_abaqus_fil(path).HasPointData("CSTRESS"));
    Mesh bad = contact_hex();
    bad.AddFieldData("abaqus:contact_pair", ints({0, 7}));
    EXPECT_THROW(meshioplusplus::write_abaqus_fil(temp_fil(), bad), WriteError);
}

TEST(AbaqusFilWriter, ModalDynamicsRecordsComeBackInBothLayouts) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const bool binary : {false, true}) {
        Mesh m = unit_hex();
        m.AddFieldData("abaqus:GU", reals({1.5, 2.5, 3.5, 4.5, 5.5, nan}, {2, 3}));
        m.AddFieldData("abaqus:GV", reals({0.25, -1.0}, {2}));
        m.AddFieldData("abaqus:BM", reals({2, 1.0, 0.5, 0.0, 0.0, 0.0, 0.25}, {1, 7}));
        m.AddFieldData("abaqus:KE", reals({7.0}, {1}));
        const std::string path = temp_fil();
        meshioplusplus::write_abaqus_fil(path, m, binary);
        const Mesh back = meshioplusplus::read_abaqus_fil(path);
        const NDArray& gu = back.FieldData("abaqus:GU");
        ASSERT_EQ(gu.Size(), 6u);
        EXPECT_EQ(gu.As<double>()[4], 5.5);
        EXPECT_TRUE(std::isnan(gu.As<double>()[5]));
        EXPECT_EQ(back.FieldData("abaqus:GV").As<double>()[1], -1.0);
        const NDArray& bm = back.FieldData("abaqus:BM");
        ASSERT_EQ(bm.Size(), 7u);
        EXPECT_EQ(bm.As<double>()[0], 2.0);
        EXPECT_EQ(bm.As<double>()[6], 0.25);
        EXPECT_EQ(back.FieldData("abaqus:KE").As<double>()[0], 7.0);
        EXPECT_FALSE(back.HasFieldData("abaqus:GPU"));
    }
}
