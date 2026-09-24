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
#include "meshioplusplus/formats/lsdyna_binout.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

// An LSDA file with 8-byte lengths and offsets, little-endian: DATA records,
// then one symbol table of CD and VARIABLE entries.
class Lsda {
public:
    Lsda() {
        const unsigned char header[8] = {8, 8, 8, 1, 1, 1, 0, 0};
        mBytes.append(reinterpret_cast<const char*>(header), 8);
        Command(17, 7);  // the symbol table's offset, patched in Finish
        mOffsetAt = mBytes.size();
        U64(0);
    }

    template <class T>
    void Var(const std::string& rDir, const std::string& rName, int Type,
             const std::vector<T>& rValues) {
        const std::size_t offset = mBytes.size();
        Command(8 + 1 + 1 + 1 + rName.size() + sizeof(T) * rValues.size(), 3);
        mBytes.push_back(static_cast<char>(Type));
        mBytes.push_back(static_cast<char>(rName.size()));
        mBytes += rName;
        for (const T& v : rValues)
            mBytes.append(reinterpret_cast<const char*>(&v), sizeof(T));
        mEntries.push_back({rDir, rName, Type, offset, rValues.size()});
    }

    std::string Finish() {
        const std::size_t table = mBytes.size();
        std::memcpy(mBytes.data() + mOffsetAt, &table, 8);
        Command(9, 5);  // begin
        std::string cwd;
        for (const Entry& e : mEntries) {
            if (e.mDir != cwd) {
                Command(9 + e.mDir.size(), 2);
                mBytes += e.mDir;
                cwd = e.mDir;
            }
            Command(9 + e.mName.size() + 1 + 8 + 8, 4);
            mBytes += e.mName;
            mBytes.push_back(static_cast<char>(e.mType));
            U64(e.mOffset);
            U64(e.mCount);
        }
        Command(17, 6);  // end: then the next table's offset, none
        U64(0);
        return mBytes;
    }

private:
    struct Entry {
        std::string mDir, mName;
        int mType;
        std::size_t mOffset, mCount;
    };
    void U64(std::uint64_t V) { mBytes.append(reinterpret_cast<const char*>(&V), 8); }
    void Command(std::uint64_t Length, int Cmd) {
        U64(Length);
        mBytes.push_back(static_cast<char>(Cmd));
    }
    std::string mBytes;
    std::size_t mOffsetAt = 0;
    std::vector<Entry> mEntries;
};

// nodout: nodes 11 and 12 over two outputs (t = 0, 1), each moved by t in x;
// glstat: one output at t = 0.5.
std::string binout() {
    Lsda l;
    l.Var<std::int32_t>("/nodout/metadata", "ids", 3, {11, 12});
    for (int k = 0; k < 2; ++k) {
        const std::string dir = "/nodout/d00000" + std::to_string(k + 1);
        l.Var<float>(dir, "time", 9, {static_cast<float>(k)});
        l.Var<float>(dir, "x_coordinate", 9, {static_cast<float>(k), 1.0f + k});
        l.Var<float>(dir, "y_coordinate", 9, {0.0f, 0.0f});
        l.Var<float>(dir, "z_coordinate", 9, {0.0f, 0.0f});
        l.Var<float>(dir, "x_displacement", 9, {static_cast<float>(k), static_cast<float>(k)});
        l.Var<float>(dir, "y_displacement", 9, {0.0f, 0.0f});
        l.Var<float>(dir, "z_displacement", 9, {0.0f, 0.0f});
        l.Var<std::int32_t>(dir, "cycle", 3, {100 * k});
    }
    l.Var<float>("/glstat/d000001", "time", 9, {0.5f});
    l.Var<double>("/glstat/d000001", "kinetic_energy", 10, {7.5});
    l.Var<std::int32_t>("/glstat/d000001", "cycle", 3, {50});
    return l.Finish();
}

std::string write_file(const std::string& rBytes, const std::string& rName = "binout") {
    const std::filesystem::path dir(
        mt::temp_path("_binout_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    std::ofstream((dir / rName).string(), std::ios::binary) << rBytes;
    return (dir / rName).string();
}

ReadOptions step(int Step) {
    ReadOptions o;
    o.mTimeStep = Step;
    return o;
}

}  // namespace

TEST(LsdynaBinout, NodoutStepsAndLatestOtherOutputs) {
    const std::string path = write_file(binout());
    EXPECT_EQ(meshioplusplus::lsdyna_binout_time_values(path), (std::vector<double>{0.0, 1.0}));
    const Mesh first = meshioplusplus::read_lsdyna_binout(path);
    ASSERT_EQ(first.NumPoints(), 2u);
    EXPECT_EQ(first.Cells(0).Type(), "vertex");
    EXPECT_EQ(first.PointData("lsdyna:nid").As<std::int64_t>()[1], 12);
    EXPECT_FALSE(first.HasFieldData("binout:glstat:kinetic_energy"));  // none yet at t = 0
    EXPECT_EQ(first.FieldData("binout:nodout:cycle").As<std::int64_t>()[0], 0);

    const Mesh last = meshioplusplus::read_lsdyna_binout(path, step(-1));
    EXPECT_EQ(last.Points().As<double>()[3], 2.0);  // node 12's x at t = 1
    EXPECT_EQ(last.PointData("displacement").As<double>()[0], 1.0);
    EXPECT_EQ(last.FieldData("binout:glstat:kinetic_energy").As<double>()[0], 7.5);
    EXPECT_EQ(last.FieldData("binout:glstat:time").As<double>()[0], 0.5);
    EXPECT_EQ(last.FieldData("binout:glstat:cycle").As<std::int64_t>()[0], 50);
}

TEST(LsdynaBinout, ResolvedSniffedNarrowedAndRefused) {
    const std::string path = write_file(binout());
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "lsdyna_binout");
    EXPECT_EQ(meshioplusplus::sniff_format(write_file(binout(), "results.lsda")), "lsdyna_binout");
    EXPECT_TRUE(meshioplusplus::is_binout_filename("/a/binout0003"));
    EXPECT_FALSE(meshioplusplus::is_binout_filename("/a/binout.txt"));
    const auto meta = meshioplusplus::registry_read_metadata(path, "lsdyna_binout", ReadOptions{});
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.0, 1.0}));
    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"displacement"};
    const Mesh narrowed = meshioplusplus::read_lsdyna_binout(path, only);
    EXPECT_TRUE(narrowed.HasPointData("displacement"));
    EXPECT_FALSE(narrowed.HasPointData("velocity"));
    EXPECT_FALSE(narrowed.HasFieldData("binout:nodout:cycle"));
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(path, step(2)), ReadError);
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(write_file(std::string(64, '\0'))), ReadError);
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(write_file(binout().substr(0, 200))),
                 ReadError);
}
