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
#include <iterator>
#include <random>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/version.hpp"
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

// ---------------------------------------------------------------------------
// The writer. `Th` stays the hand-framed oracle: it frames each record itself,
// so the expected bytes below do not go through the framer the writer uses.
// ---------------------------------------------------------------------------

namespace {

using meshioplusplus::NDArray;
using meshioplusplus::WriteError;

void add_field(Mesh& rMesh, const std::string& rName, const std::vector<double>& rValues,
               std::vector<std::size_t> Shape = {}) {
    if (Shape.empty())
        Shape = {rValues.size()};
    NDArray a(meshioplusplus::DType::Float64, std::move(Shape));
    for (std::size_t i = 0; i < rValues.size(); ++i)
        a.As<double>()[i] = rValues[i];
    rMesh.AddFieldData(rName, std::move(a));
}

void add_ints(Mesh& rMesh, const std::string& rName, const std::vector<std::int64_t>& rValues) {
    NDArray a(meshioplusplus::DType::Int64, {rValues.size()});
    for (std::size_t i = 0; i < rValues.size(); ++i)
        a.As<std::int64_t>()[i] = rValues[i];
    rMesh.AddFieldData(rName, std::move(a));
}

// The oracle file's content as field data (see th_file): globals 1, 2; part 7
// (IE, KE); subset 0 (KE); TH groups node 3 (ids 11, 12; variables 1, 2, 3) and
// rigid body 4 (entity 5, no variables).
Mesh sample_mesh(int K, bool Factors = false) {
    Mesh m;
    m.AssignPoints(NDArray(meshioplusplus::DType::Float64, {0, 3}));
    m.AddCellBlock("vertex", NDArray(meshioplusplus::DType::Int64, {0, 1}));
    const double k = K;
    add_field(m, "meshio:time", {0.5 * k});
    add_field(m, "radioss_th:global:internal_energy", {10.0 + k});
    add_field(m, "radioss_th:global:kinetic_energy", {20.0 + k});
    add_field(m, "radioss_th:part:7:IE", {1.0 + k});
    add_field(m, "radioss_th:part:7:KE", {2.0 + k});
    add_field(m, "radioss_th:subset:0:KE", {3.0 + k});
    add_field(
        m, "radioss_th:node:3",
        {100.0 * k, 1 + 100.0 * k, 2 + 100.0 * k, 3 + 100.0 * k, 4 + 100.0 * k, 5 + 100.0 * k},
        {2, 3});
    add_ints(m, "radioss_th:node:3:ids", {11, 12});
    add_ints(m, "radioss_th:node:3:variables", {1, 2, 3});
    add_field(m, "radioss_th:rbody:4", {}, {1, 0});
    add_ints(m, "radioss_th:rbody:4:ids", {5});
    if (Factors)
        add_field(m, "radioss_th:unit_factors", {1.0, 1000.0, 0.001});
    return m;
}

// What the writer must produce for `sample_mesh(k)` for each k in Steps.
std::string expected_bytes(const std::vector<int>& rSteps, bool Factors = false) {
    const std::size_t title = Factors ? 100 : 40;
    const std::string tag = std::string("Written by meshio++ v") + MESHIOPLUSPLUS_VERSION_STRING;
    Th t;
    t.Open().Int(Factors ? 4021 : 3040).Text(tag, 80).Close();
    t.Open().Text("meshio++", 80).Close();
    if (Factors) {
        t.Ints({1});
        t.Ints({static_cast<std::int32_t>(title)});
        t.Open().Float(1.0f).Float(1000.0f).Float(0.001f).Close();
    }
    t.Ints({1, 0, 0, 1, 2, 2});
    t.Ints({1, 2});
    t.Open().Int(7).Text("PART 7", title).Int(0).Int(0).Int(0).Int(2).Close();
    t.Ints({1, 2});
    t.Open().Int(0).Int(0).Int(0).Int(0).Int(1).Text("SUBSET 0", title).Close();
    t.Ints({2});
    t.Open().Int(3).Int(0).Int(0).Int(2).Int(3).Text("TH node 3", title).Close();
    t.Open().Int(11).Text("node 11", title).Close();
    t.Open().Int(12).Text("node 12", title).Close();
    t.Ints({1, 2, 3});
    t.Open().Int(4).Int(103).Int(0).Int(1).Int(0).Text("TH rbody 4", title).Close();
    t.Open().Int(5).Text("rbody 5", title).Close();
    for (const int k : rSteps) {
        const float s = static_cast<float>(k);
        t.Open().Float(0.5f * s).Close();
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

std::string read_bytes(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::string temp_th(const std::string& rName = "runT01") {
    const std::filesystem::path dir(
        mt::temp_path("_radioss_th_w_" + std::to_string(std::random_device{}())));
    std::filesystem::create_directories(dir);
    return (dir / rName).string();
}

}  // namespace

TEST(RadiossThWriter, FollowsTheHandFramedLayout) {
    for (const bool factors : {false, true}) {
        const std::string path = temp_th();
        meshioplusplus::write_radioss_th(path, sample_mesh(1, factors));
        EXPECT_EQ(read_bytes(path), expected_bytes({1}, factors)) << "factors " << factors;
    }
}

TEST(RadiossThWriter, ReadsBackWhatTheOracleFileHolds) {
    // The hand-framed file (version 3040, a material, subset parts) and the
    // written one hold the same variables.
    const Mesh oracle = meshioplusplus::read_radioss_th(write_file(th_file()), step(1));
    const std::string path = temp_th();
    meshioplusplus::write_radioss_th(path, sample_mesh(1));
    const Mesh written = meshioplusplus::read_radioss_th(path);
    EXPECT_EQ(oracle.NumFieldData(), written.NumFieldData());
    for (const std::string& name : oracle.FieldDataNames()) {
        ASSERT_TRUE(written.HasFieldData(name)) << name;
        const auto& a = oracle.FieldData(name);
        const auto& b = written.FieldData(name);
        EXPECT_EQ(a.Dtype(), b.Dtype()) << name;
        EXPECT_EQ(a.Shape(), b.Shape()) << name;
        if (a.Size() > 0)
            EXPECT_EQ(0, std::memcmp(a.Data(), b.Data(),
                                     a.Size() * meshioplusplus::dtype_size(a.Dtype())))
                << name;
    }
}

TEST(RadiossThWriter, DefaultsDropsAndRefusals) {
    // Ids and codes default to 1..n; unnamed codes are var<code>.
    Mesh m;
    add_field(m, "radioss_th:shell:2", {0, 1, 2, 3, 4, 5}, {3, 2});
    add_field(m, "radioss_th:global:var40", {1.5});
    add_field(m, "radioss_th:part:-3:var26", {2.5});
    add_field(m, "unrelated", {1.0});
    const std::string path = temp_th();
    meshioplusplus::write_radioss_th(path, m);
    const Mesh back = meshioplusplus::read_radioss_th(path);
    EXPECT_EQ(back.FieldData("radioss_th:shell:2:ids").As<std::int64_t>()[2], 3);
    EXPECT_EQ(back.FieldData("radioss_th:shell:2:variables").As<std::int64_t>()[1], 2);
    EXPECT_EQ(scalar(back, "radioss_th:global:var40"), 1.5);
    EXPECT_EQ(scalar(back, "radioss_th:part:-3:var26"), 2.5);
    EXPECT_FALSE(back.HasFieldData("unrelated"));

    Mesh flat = sample_mesh(0);
    add_field(flat, "radioss_th:node:3", {1, 2, 3});
    EXPECT_THROW(meshioplusplus::write_radioss_th(temp_th(), flat), WriteError);
    Mesh ids = sample_mesh(0);
    add_ints(ids, "radioss_th:node:3:ids", {1, 2, 3});
    EXPECT_THROW(meshioplusplus::write_radioss_th(temp_th(), ids), WriteError);
    Mesh fractional = sample_mesh(0);
    add_field(fractional, "radioss_th:node:3:ids", {1.5, 2.0});
    EXPECT_THROW(meshioplusplus::write_radioss_th(temp_th(), fractional), WriteError);
    Mesh factors = sample_mesh(0);
    add_field(factors, "radioss_th:unit_factors", {1.0, 2.0});
    EXPECT_THROW(meshioplusplus::write_radioss_th(temp_th(), factors), WriteError);
    EXPECT_THROW(meshioplusplus::write_radioss_th("/no/such/dir/runT01", sample_mesh(0)),
                 WriteError);
}

TEST(RadiossThWriter, ASeriesIsTheDescriptionsOnceThenAnOutputPerStep) {
    const std::string path = temp_th();
    meshioplusplus::RadiossThSeriesWriter writer(path);
    writer.Write(0.0, sample_mesh(0, true));
    // a step with other variables is refused, and does not count
    Mesh other = sample_mesh(1, true);
    Mesh fewer;
    for (const std::string& name : other.FieldDataNames())
        if (name != "radioss_th:part:7:KE")
            fewer.AddFieldData(name, other.FieldData(name));
    EXPECT_THROW(writer.Write(0.5, fewer), WriteError);
    EXPECT_EQ(writer.NumSteps(), 1u);
    // the time is the argument's, not the mesh's
    writer.Write(0.5, sample_mesh(1, true));
    writer.Write(1.0, sample_mesh(2, true));
    writer.Finalize();
    writer.Finalize();
    EXPECT_EQ(read_bytes(path), expected_bytes({0, 1, 2}, true));
    EXPECT_EQ(meshioplusplus::radioss_th_time_values(path), (std::vector<double>{0.0, 0.5, 1.0}));

    meshioplusplus::RadiossThSeriesWriter empty(temp_th());
    EXPECT_THROW(empty.Finalize(), WriteError);

    const std::string timed_path = temp_th();
    meshioplusplus::RadiossThSeriesWriter timed(timed_path);
    timed.Write(7.0, sample_mesh(0));
    timed.Finalize();
    EXPECT_EQ(meshioplusplus::radioss_th_time_values(timed_path), (std::vector<double>{7.0}));
}

TEST(RadiossThWriter, ACutShortWrittenFileKeepsItsCompleteOutputs) {
    const std::string full = expected_bytes({0, 1});
    Th tail;
    tail.Open().Float(1.0f).Close();
    const std::string lone_time = tail.mBytes;
    for (const std::string& extra : {lone_time.substr(0, 7), lone_time, lone_time + lone_time,
                                     lone_time.substr(0, lone_time.size() - 1) + '\x07'}) {
        const std::string path = temp_th();
        std::ofstream(path, std::ios::binary) << full << extra;
        EXPECT_EQ(meshioplusplus::radioss_th_time_values(path), (std::vector<double>{0.0, 0.5}));
    }
    const std::string path = temp_th();
    std::ofstream(path, std::ios::binary) << full.substr(0, full.size() - 20);
    EXPECT_EQ(meshioplusplus::radioss_th_time_values(path), (std::vector<double>{0.0}));
}

TEST(RadiossThWriter, TheRegistryWritesItByNameAndByFormat) {
    EXPECT_EQ(meshioplusplus::resolve_write_format("/a/runT01", ""), "radioss_th");
    const std::string path = temp_th("history.bin");
    meshioplusplus::registry_writers().at("radioss_th")(path, sample_mesh(1));
    EXPECT_EQ(read_bytes(path), expected_bytes({1}));
    EXPECT_EQ(meshioplusplus::sniff_format(path), "radioss_th");
}
