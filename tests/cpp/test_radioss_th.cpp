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
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/radioss_th.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

// A big-endian Fortran unformatted file, record by record.
class Th {
public:
    // Opens a record; Close() ends it.
    Th& Open() {
        mStart = mBytes.size();
        Word(0);
        return *this;
    }
    Th& Close() {
        const std::uint32_t n = static_cast<std::uint32_t>(mBytes.size() - mStart - 4);
        for (int i = 0; i < 4; ++i)
            mBytes[mStart + i] = static_cast<char>(n >> (24 - 8 * i));
        Word(n);
        return *this;
    }
    Th& Word(std::uint32_t V) {
        for (int i = 0; i < 4; ++i)
            mBytes.push_back(static_cast<char>(V >> (24 - 8 * i)));
        return *this;
    }
    Th& Int(std::int32_t V) { return Word(static_cast<std::uint32_t>(V)); }
    Th& Float(float V) {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &V, 4);
        return Word(bits);
    }
    Th& Text(const std::string& rText, std::size_t Width) {
        std::string t = rText;
        t.resize(Width, ' ');
        mBytes += t;
        return *this;
    }
    Th& Ints(const std::vector<std::int32_t>& rValues) {
        Open();
        for (const std::int32_t v : rValues)
            Int(v);
        return Close();
    }
    std::string mBytes;

private:
    std::size_t mStart = 0;
};

// Globals 1, 2; parts 7 (IE, KE) and 8 (none); subset 0 (KE); TH groups: nodes
// 11 and 12 (3 variables), rigid body 5 (none). Two outputs (t = 0, 0.5) and a
// third cut short. Version 4021 adds the unit factors and 100-character titles.
std::string th_file(std::int32_t Version = 3040) {
    const std::size_t title = Version >= 4021 ? 100 : 40;
    Th t;
    t.Open().Int(Version).Text("run", 80).Close();
    t.Open().Text("Mon Mar 30 16:39:04 2026 RADIOSS", 80).Close();
    if (Version > 3050) {
        t.Ints({1});
        t.Ints({static_cast<std::int32_t>(title)});
        t.Open().Float(1.0f).Float(1000.0f).Float(0.001f).Close();
    }
    t.Ints({2, 1, 0, 1, 2, 2});  // parts, materials, properties, subsets, groups, globals
    t.Ints({1, 2});
    t.Open().Int(7).Text("part7", title).Int(0).Int(0).Int(0).Int(2).Close();
    t.Ints({1, 2});
    t.Open().Int(8).Text("part8", title).Int(0).Int(0).Int(0).Int(0).Close();
    t.Open().Int(1).Text("steel", title).Close();
    t.Open().Int(0).Int(0).Int(0).Int(2).Int(1).Text("GLOBAL MODEL", title).Close();
    t.Ints({7, 8});
    t.Ints({2});
    t.Open().Int(3).Int(0).Int(0).Int(2).Int(3).Text("th_nodes", title).Close();
    t.Open().Int(11).Text("a", title).Close();
    t.Open().Int(12).Text("b", title).Close();
    t.Ints({1, 2, 3});
    t.Open().Int(4).Int(103).Int(0).Int(1).Int(0).Text("TH RBODY", title).Close();
    t.Open().Int(5).Text("body", title).Close();
    for (int k = 0; k < 3; ++k) {
        const float s = static_cast<float>(k);
        t.Open().Float(0.5f * s).Close();
        if (k == 2)
            break;  // cut short
        t.Open().Float(10.0f + s).Float(20.0f + s).Close();
        t.Open().Float(1.0f + s).Float(2.0f + s).Close();
        t.Open().Float(3.0f + s).Close();
        t.Open();
        for (int i = 0; i < 6; ++i)
            t.Float(static_cast<float>(i) + 100.0f * s);
        t.Close();
        t.Open().Close();
    }
    return t.mBytes;
}

std::string write_file(const std::string& rBytes, const std::string& rName = "runT01") {
    const std::filesystem::path dir(
        mt::temp_path("_radioss_th_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    std::ofstream((dir / rName).string(), std::ios::binary) << rBytes;
    return (dir / rName).string();
}

ReadOptions step(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

double scalar(const Mesh& rMesh, const std::string& rName) {
    return rMesh.FieldData(rName).As<double>()[0];
}

}  // namespace

TEST(RadiossTh, OutputsGlobalsPartsSubsetsAndGroups) {
    const std::string path = write_file(th_file());
    EXPECT_EQ(meshioplusplus::radioss_th_time_values(path), (std::vector<double>{0.0, 0.5}));
    const Mesh last = meshioplusplus::read_radioss_th(path, step(-1));
    EXPECT_EQ(last.NumPoints(), 0u);
    EXPECT_EQ(scalar(last, "meshio:time"), 0.5);
    EXPECT_EQ(scalar(last, "radioss_th:global:internal_energy"), 11.0);
    EXPECT_EQ(scalar(last, "radioss_th:global:kinetic_energy"), 21.0);
    EXPECT_EQ(scalar(last, "radioss_th:part:7:IE"), 2.0);
    EXPECT_EQ(scalar(last, "radioss_th:part:7:KE"), 3.0);
    EXPECT_EQ(scalar(last, "radioss_th:subset:0:KE"), 4.0);
    const auto& nodes = last.FieldData("radioss_th:node:3");
    ASSERT_EQ(nodes.Shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_EQ(nodes.As<double>()[5], 105.0);
    EXPECT_EQ(last.FieldData("radioss_th:node:3:ids").As<std::int64_t>()[1], 12);
    EXPECT_EQ(last.FieldData("radioss_th:node:3:variables").As<std::int64_t>()[2], 3);
    EXPECT_EQ(last.FieldData("radioss_th:rbody:4").Shape(), (std::vector<std::size_t>{1, 0}));
    EXPECT_FALSE(last.HasFieldData("radioss_th:unit_factors"));
    EXPECT_THROW(meshioplusplus::read_radioss_th(path, step(2)), ReadError);
}

TEST(RadiossTh, NewerVersionUnitFactorsAndLongTitles) {
    const std::string path = write_file(th_file(4021));
    const Mesh first = meshioplusplus::read_radioss_th(path);
    EXPECT_EQ(scalar(first, "radioss_th:part:7:KE"), 2.0);
    const auto& f = first.FieldData("radioss_th:unit_factors");
    EXPECT_EQ(f.As<double>()[1], 1000.0);
}

TEST(RadiossTh, ResolvedSniffedNarrowedAndRefused) {
    const std::string path = write_file(th_file());
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "radioss_th");
    EXPECT_EQ(meshioplusplus::sniff_format(write_file(th_file(), "history.bin")), "radioss_th");
    EXPECT_TRUE(meshioplusplus::is_radioss_th_filename("/a/crashT02"));
    EXPECT_FALSE(meshioplusplus::is_radioss_th_filename("/a/crashT02.csv"));
    EXPECT_FALSE(meshioplusplus::is_radioss_th_filename("/a/crashA001"));
    const auto meta = meshioplusplus::registry_read_metadata(path, "radioss_th", ReadOptions{});
    EXPECT_EQ(meta.mTimeValues.size(), 2u);
    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"radioss_th:part:7:IE"};
    const Mesh narrowed = meshioplusplus::read_radioss_th(path, only);
    EXPECT_TRUE(narrowed.HasFieldData("radioss_th:part:7:IE"));
    EXPECT_FALSE(narrowed.HasFieldData("radioss_th:node:3"));
    EXPECT_THROW(meshioplusplus::read_radioss_th(write_file(std::string(64, '\0'))), ReadError);
    // the descriptions and one time: no complete output
    const std::string text = th_file();
    EXPECT_THROW(meshioplusplus::read_radioss_th(write_file(text.substr(0, 848))), ReadError);
}
