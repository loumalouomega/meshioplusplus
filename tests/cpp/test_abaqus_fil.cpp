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
 *        and element results, long set labels, steps and metadata.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus_fil.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::RegionKind;
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

std::string ascii() {
    std::string flat;
    for (const Record& rec : records()) {
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
    const std::string path = mt::temp_path(".fil");
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
    const auto& s = rMesh.CellData("S", 0);
    ASSERT_EQ(s.Shape().size(), 3u);  // (cells, points, components)
    EXPECT_EQ(s.Shape()[2], 6u);
    // Abaqus's 13 and 23 are swapped into meshio++'s yz, zx.
    EXPECT_DOUBLE_EQ(detail::read_double(s, 4), 100.0 * Increment + 5);
    EXPECT_DOUBLE_EQ(detail::read_double(s, 5), 100.0 * Increment + 4);
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
