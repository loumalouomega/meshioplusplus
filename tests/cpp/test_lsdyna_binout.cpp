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
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <optional>
#include <random>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/lsdyna_binout.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::WriteError;

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

namespace {

// Overwrites the eight bytes at At with V (the file is little-endian).
void patch_u64(std::string& rBytes, std::size_t At, std::uint64_t V) {
    std::memcpy(&rBytes[At], &V, 8);
}

}  // namespace

// A table whose next-table offset names itself: a loop, which the reader refuses rather than
// reading the table again forever.
TEST(LsdynaBinout, ATableThatNamesItselfAsTheNextIsRefused) {
    std::string bytes = binout();
    std::uint64_t first = 0;
    std::memcpy(&first, bytes.data() + 17, 8);  // the first table's offset
    patch_u64(bytes, bytes.size() - 8, first);  // the last table's next-table offset
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(write_file(bytes)), ReadError);
}

// The fuzzer's input for the same loop (tests/fuzz/regressions/lsdyna_binout).
TEST(LsdynaBinout, TheFuzzersSymbolTableCycleIsRefused) {
    const std::filesystem::path input = std::filesystem::path(__FILE__).parent_path() / ".." /
                                        "fuzz" / "regressions" / "lsdyna_binout" /
                                        "symbol-table-cycle";
    ASSERT_TRUE(std::filesystem::exists(input));
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(input.string()), ReadError);
}

// Offsets and lengths near 2^64 wrap when added to a position: each one is refused, not read.
TEST(LsdynaBinout, OffsetsAndLengthsNearTwoToTheSixtyFourAreRefused) {
    std::string past_end = binout();
    patch_u64(past_end, 17, 0xFFFFFFFFFFFFFFF8ull);  // the first table's offset
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(write_file(past_end)), ReadError);

    std::string huge = binout();
    std::uint64_t first = 0;
    std::memcpy(&first, huge.data() + 17, 8);
    patch_u64(huge, first + 9, 0xFFFFFFFFFFFFFFFFull);  // the first record's length
    EXPECT_THROW(meshioplusplus::read_lsdyna_binout(write_file(huge)), ReadError);
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

std::string temp_binout(const std::string& rName = "binout") {
    const std::filesystem::path dir(
        mt::temp_path("_binout_w_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    return (dir / rName).string();
}

// Two nodes (ids 11, 12) with a displacement, a cycle count and one glstat output.
Mesh oracle_mesh() {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "vertex", {{0}, {1}});
    m.AddPointData("lsdyna:nid", ints({11, 12}));
    m.AddPointData("displacement", reals({0.5, 0, 0, 0.5, 0, 0}, {2, 3}));
    m.AddFieldData("meshio:time", reals({1.0}, {1}));
    m.AddFieldData("binout:nodout:cycle", ints({100}));
    m.AddFieldData("binout:glstat:time", reals({0.5}, {1}));
    m.AddFieldData("binout:glstat:kinetic_energy", reals({7.5}, {1}));
    return m;
}

std::string oracle_bytes() {
    Lsda l;
    l.Var<std::int32_t>("/nodout/metadata", "ids", 3, {11, 12});
    const std::string dir = "/nodout/d000001";
    l.Var<float>(dir, "time", 9, {1.0f});
    l.Var<float>(dir, "x_coordinate", 9, {0.0f, 1.0f});
    l.Var<float>(dir, "y_coordinate", 9, {0.0f, 0.0f});
    l.Var<float>(dir, "z_coordinate", 9, {0.0f, 0.0f});
    l.Var<float>(dir, "x_displacement", 9, {0.5f, 0.5f});
    l.Var<float>(dir, "y_displacement", 9, {0.0f, 0.0f});
    l.Var<float>(dir, "z_displacement", 9, {0.0f, 0.0f});
    l.Var<std::int32_t>(dir, "cycle", 3, {100});
    l.Var<float>("/glstat/d000001", "time", 9, {0.5f});
    l.Var<float>("/glstat/d000001", "kinetic_energy", 9, {7.5f});
    return l.Finish();
}

// Three nodes over steps, with glstat and elout/beam written less often.
Mesh step_mesh(double T, std::optional<double> Glstat = std::nullopt, bool Temperature = true) {
    Mesh m = mt::make_mesh({{T, 0, 0}, {1 + T, 0, 0}, {1 + T, 1, 0}}, "vertex", {{0}, {1}, {2}});
    m.AddPointData("lsdyna:nid", ints({11, 12, 13}));
    NDArray d(meshioplusplus::DType::Float64, {std::size_t{3}, std::size_t{3}});
    for (int i = 0; i < 9; ++i)
        d.As<double>()[i] = i * T;
    m.AddPointData("displacement", std::move(d));
    if (Temperature)
        m.AddPointData("temperature", reals({T, 2 * T, 3 * T}, {3}));
    if (Glstat) {
        m.AddFieldData("binout:glstat:time", reals({*Glstat}, {1}));
        m.AddFieldData("binout:glstat:kinetic_energy", reals({7.5 * *Glstat}, {1}));
        m.AddFieldData("binout:elout/beam:time", reals({*Glstat}, {1}));
        m.AddFieldData("binout:elout/beam:axial", reals({*Glstat, 2 * *Glstat}, {2}));
        m.AddFieldData("binout:elout/beam:ids", ints({5, 6}));
    }
    return m;
}

std::vector<std::pair<double, Mesh>> series_steps() {
    std::vector<std::pair<double, Mesh>> steps;
    steps.emplace_back(0.0, step_mesh(0.0));
    steps.emplace_back(0.1, step_mesh(0.1, 0.0));
    steps.emplace_back(0.2, step_mesh(0.2, 0.0));
    steps.emplace_back(0.3, step_mesh(0.3, 0.25));
    return steps;
}

}  // namespace

TEST(LsdynaBinoutWriter, FollowsTheHandBuiltLayout) {
    const std::string path = temp_binout();
    meshioplusplus::write_lsdyna_binout(path, oracle_mesh());
    EXPECT_EQ(read_all(path), oracle_bytes());
}

TEST(LsdynaBinoutWriter, ReadsBackWhatTheOracleFileHolds) {
    const std::string path = temp_binout();
    meshioplusplus::write_lsdyna_binout(path, oracle_mesh());
    const Mesh back = meshioplusplus::read_lsdyna_binout(path);
    ASSERT_EQ(back.NumPoints(), 2u);
    EXPECT_EQ(back.PointData("lsdyna:nid").As<std::int64_t>()[1], 12);
    EXPECT_EQ(back.Points().As<double>()[3], 1.0);
    EXPECT_EQ(back.PointData("displacement").As<double>()[0], 0.5);
    EXPECT_EQ(back.FieldData("binout:nodout:cycle").As<std::int64_t>()[0], 100);
    EXPECT_EQ(back.FieldData("binout:glstat:kinetic_energy").As<double>()[0], 7.5);
    EXPECT_EQ(back.FieldData("meshio:time").As<double>()[0], 1.0);
}

TEST(LsdynaBinoutWriter, ASeriesKeepsEachDatabasesOwnCadence) {
    const std::string path = temp_binout();
    meshioplusplus::LsdynaBinoutSeriesWriter writer(path);
    const auto steps = series_steps();
    for (std::size_t i = 0; i < steps.size(); ++i) {
        writer.Write(steps[i].first, steps[i].second);
        EXPECT_EQ(writer.NumSteps(), i + 1);
    }
    writer.Finalize();
    writer.Finalize();
    const auto times = meshioplusplus::lsdyna_binout_time_values(path);
    ASSERT_EQ(times.size(), 4u);
    EXPECT_NEAR(times[3], 0.3, 1e-6);
    const Mesh last = meshioplusplus::read_lsdyna_binout(path, step(-1));
    EXPECT_EQ(last.FieldData("binout:glstat:kinetic_energy").As<double>()[0], 1.875);
    EXPECT_EQ(last.FieldData("binout:glstat:time").As<double>()[0], 0.25);
    EXPECT_EQ(last.FieldData("binout:elout/beam:ids").As<std::int64_t>()[1], 6);
    EXPECT_NEAR(last.PointData("temperature").As<double>()[2], 0.9, 1e-6);
    const Mesh second = meshioplusplus::read_lsdyna_binout(path, step(1));
    EXPECT_EQ(second.FieldData("binout:glstat:kinetic_energy").As<double>()[0], 0.0);
    // glstat has two outputs, nodout four: the file says so
    std::size_t glstat_folders = 0;
    const std::string bytes = read_all(path);
    for (std::size_t at = bytes.find("kinetic_energy"); at != std::string::npos;
         at = bytes.find("kinetic_energy", at + 1))
        ++glstat_folders;
    EXPECT_GE(glstat_folders, 2u);
}

TEST(LsdynaBinoutWriter, ASeriesIsOneStepsBytesWhenItHoldsOne) {
    Mesh mesh = step_mesh(0.5, 0.25);
    mesh.AddFieldData("meshio:time", reals({0.5}, {1}));
    const std::string plain = temp_binout(), series = temp_binout();
    meshioplusplus::write_lsdyna_binout(plain, mesh);
    meshioplusplus::LsdynaBinoutSeriesWriter writer(series);
    writer.Write(0.5, mesh);
    writer.Finalize();
    EXPECT_EQ(read_all(plain), read_all(series));
}

TEST(LsdynaBinoutWriter, ASeriesRefusesOtherNodesAndVariables) {
    meshioplusplus::LsdynaBinoutSeriesWriter writer(temp_binout());
    writer.Write(0.1, step_mesh(0.1, 0.0));
    const Mesh fewer = step_mesh(0.2, 0.0, false);
    EXPECT_THROW(writer.Write(0.2, fewer), WriteError);
    Mesh ids = step_mesh(0.2, 0.0);
    ids.AddPointData("lsdyna:nid", ints({1, 2, 3}));
    EXPECT_THROW(writer.Write(0.2, ids), WriteError);
    Mesh extra = step_mesh(0.2, 0.0);
    extra.AddFieldData("binout:glstat:extra", reals({1.0}, {1}));
    EXPECT_THROW(writer.Write(0.2, extra), WriteError);
    Mesh other = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "vertex", {{0}, {1}});
    EXPECT_THROW(writer.Write(0.2, other), WriteError);
    EXPECT_EQ(writer.NumSteps(), 1u);
    meshioplusplus::LsdynaBinoutSeriesWriter empty(temp_binout());
    EXPECT_THROW(empty.Finalize(), WriteError);
}

TEST(LsdynaBinoutWriter, RefusalsAndDrops) {
    Mesh big = step_mesh(0.5);
    big.AddPointData("lsdyna:nid", ints({1, 2, std::int64_t{1} << 40}));
    EXPECT_THROW(meshioplusplus::write_lsdyna_binout(temp_binout(), big), WriteError);
    for (const char* name : {"binout:d000001:x", "binout:metadata:x", "binout:a//b:x"}) {
        Mesh mesh = step_mesh(0.5);
        mesh.AddFieldData(name, reals({1.0}, {1}));
        EXPECT_THROW(meshioplusplus::write_lsdyna_binout(temp_binout(), mesh), WriteError);
    }
    Mesh nothing;
    nothing.AssignPoints(NDArray(meshioplusplus::DType::Float64, {std::size_t{0}, std::size_t{3}}));
    EXPECT_THROW(meshioplusplus::write_lsdyna_binout(temp_binout(), nothing), WriteError);
    EXPECT_THROW(meshioplusplus::write_lsdyna_binout("/no/such/dir/binout", step_mesh(0.5)),
                 WriteError);
    // A point-free mesh is a database of its own; arrays with no place are dropped.
    Mesh db;
    db.AssignPoints(NDArray(meshioplusplus::DType::Float64, {std::size_t{0}, std::size_t{3}}));
    db.AddFieldData("binout:glstat:kinetic_energy", reals({3.5}, {1}));
    db.AddFieldData("meshio:time", reals({2.0}, {1}));
    db.AddFieldData("note", reals({1.0}, {1}));
    const std::string path = temp_binout();
    meshioplusplus::write_lsdyna_binout(path, db);
    EXPECT_EQ(meshioplusplus::lsdyna_binout_time_values(path), (std::vector<double>{2.0}));
    EXPECT_EQ(meshioplusplus::read_lsdyna_binout(path).NumPoints(), 0u);
}

TEST(LsdynaBinoutWriter, TheRegistryAndTheSequenceDriverWriteIt) {
    EXPECT_EQ(meshioplusplus::resolve_write_format("/a/binout", ""), "lsdyna_binout");
    const std::string path = temp_binout("results.bin");
    meshioplusplus::registry_writers().at("lsdyna_binout")(path, oracle_mesh());
    EXPECT_EQ(read_all(path), oracle_bytes());
    EXPECT_EQ(meshioplusplus::sniff_format(path), "lsdyna_binout");
    std::string why;
    EXPECT_TRUE(meshioplusplus::sequence_write_supports_time("lsdyna_binout", why));
}
