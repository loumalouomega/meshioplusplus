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
 * @file test_radioss_anim.cpp
 * @brief OpenRadioss animation files: a minimal state written byte by byte (a
 *        quad and a degenerate-brick tetra with a scalar, a vector and a 3-D
 *        tensor), its time, parts and alive flags, dispatch and refusals.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <random>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/radioss_anim.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
namespace detail = meshioplusplus::detail;

class Be {
public:
    void Ints(std::initializer_list<std::int32_t> Values) {
        for (std::int32_t v : Values)
            U32(static_cast<std::uint32_t>(v));
    }
    void Floats(std::initializer_list<float> Values) {
        for (float f : Values) {
            std::uint32_t v;
            std::memcpy(&v, &f, 4);
            U32(v);
        }
    }
    void Text(const std::string& rS, std::size_t N) {
        std::string t = rS.substr(0, N - 1);
        t.resize(N, '\0');
        mOut += t;
    }
    void Bytes(std::initializer_list<unsigned char> Values) {
        for (unsigned char c : Values)
            mOut += static_cast<char>(c);
    }
    std::string mOut;

private:
    void U32(std::uint32_t v) {
        for (int b = 3; b >= 0; --b)
            mOut += static_cast<char>((v >> (8 * b)) & 0xff);
    }
};

// Nodes 0..4 of a unit square (z = 0) and an apex; a quad shell and a tetra
// written as the brick 0 1 2 2 4 4 4 4.
std::string state() {
    Be o;
    o.Ints({0x542C});
    o.Floats({0.25f});
    for (const char* t : {"Time=", "ModAnim", "Radioss Run="})
        o.Text(t, 81);
    o.Ints({0, 1, 1, 0, 0, 0, 0, 0, 0, 0});  // node/element ids and 3-D
    o.Ints({5, 1, 1, 0, 1, 1, 0, 0});  // nodes, facets, parts, 0 nodal + 1 element scalar, 1 vector
    o.Floats({0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1});
    o.Ints({0, 1, 2, 3});
    o.Bytes({0xff});
    o.Ints({1});
    o.Text("        7:skin", 50);
    for (int k = 0; k < 15; ++k)
        o.Bytes({0, 0});  // normals
    o.Text("Thickness", 81);
    o.Floats({2.5f});
    o.Text("Velocity", 81);
    o.Floats({1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0});
    o.Ints({10, 11, 12, 13, 14});  // node ids
    o.Ints({500});                 // facet id
    // 3-D: one brick, one part, one scalar, one tensor
    o.Ints({1, 1, 1, 1});
    o.Ints({0, 1, 2, 2, 4, 4, 4, 4});
    o.Bytes({0});  // deleted
    o.Ints({1});
    o.Text("        9:core", 50);
    o.Text("Thickness", 81);
    o.Floats({7.0f});
    o.Text("Stress", 81);
    o.Floats({1, 2, 3, 4, 5, 6});
    o.Ints({900});  // element id
    return o.mOut;
}

std::string write(const std::string& rBody, const std::string& rName) {
    const std::string dir = mt::temp_path("_anim_" + std::to_string(std::random_device{}()));
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/" + rName;
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

}  // namespace

TEST(RadiossAnim, ReadsAState) {
    const std::string path = write(state(), "runA001");
    EXPECT_TRUE(meshioplusplus::is_radioss_anim_filename(path));
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "radioss_anim");
    const Mesh mesh = meshioplusplus::read_radioss_anim(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "quad");
    EXPECT_EQ(mesh.Cells(1).Type(), "tetra");
    EXPECT_FLOAT_EQ(static_cast<float>(mesh.FieldData("meshio:time").As<double>()[0]), 0.25f);
    EXPECT_EQ(detail::read_double(mesh.CellData("Thickness", 0), 0), 2.5);
    EXPECT_EQ(detail::read_double(mesh.CellData("Thickness", 1), 0), 7.0);
    // The quad has no tensor: NaN; the brick's six values as written.
    EXPECT_TRUE(std::isnan(detail::read_double(mesh.CellData("Stress", 0), 0)));
    EXPECT_EQ(detail::read_double(mesh.CellData("Stress", 1), 5), 6.0);
    EXPECT_EQ(detail::read_int(mesh.CellData("radioss:alive", 0), 0), 1);  // 0xff
    EXPECT_EQ(detail::read_int(mesh.CellData("radioss:alive", 1), 0), 0);
    EXPECT_EQ(detail::read_int(mesh.CellData("radioss:part", 1), 0), 9);
    EXPECT_EQ(detail::read_int(mesh.CellData("radioss:element_id", 1), 0), 900);
    EXPECT_EQ(detail::read_int(mesh.PointData("radioss:node_id"), 4), 14);
    EXPECT_EQ(mesh.NumRegions(), 2u);
    const meshioplusplus::MeshMetadata meta = meshioplusplus::read_radioss_anim_metadata(path);
    ASSERT_EQ(meta.mTimeValues.size(), 1u);
    EXPECT_FLOAT_EQ(static_cast<float>(meta.mTimeValues[0]), 0.25f);
}

TEST(RadiossAnim, Refusals) {
    const std::string good = state();
    EXPECT_THROW(
        meshioplusplus::read_radioss_anim(write(good.substr(0, good.size() / 2), "cutA001")),
        ReadError);
    std::string other = good;
    other[3] = 0x2a;  // magic 0x542A
    EXPECT_THROW(meshioplusplus::read_radioss_anim(write(other, "oldA001")), ReadError);
}

TEST(RadiossAnimWriter, StateRoundTripsAndSeriesNumbersFiles) {
    const std::string path = write(state(), "srcA001");
    const Mesh original = meshioplusplus::read_radioss_anim(path);
    const std::string out = (std::filesystem::path(path).parent_path() / "outA001").string();
    meshioplusplus::write_radioss_anim(out, original);
    const Mesh again = meshioplusplus::read_radioss_anim(out);
    ASSERT_EQ(again.NumCellBlocks(), original.NumCellBlocks());
    EXPECT_EQ(again.NumPoints(), original.NumPoints());
    for (std::size_t b = 0; b < original.NumCellBlocks(); ++b)
        EXPECT_EQ(again.Cells(b).Type(), original.Cells(b).Type());
    EXPECT_EQ(detail::read_double(again.CellData("Thickness", 1), 0), 7.0);
    EXPECT_EQ(detail::read_int(again.CellData("radioss:part", 1), 0), 9);
    EXPECT_EQ(detail::read_int(again.PointData("radioss:node_id"), 4), 14);
    EXPECT_EQ(again.NumRegions(), original.NumRegions());

    const std::string first = (std::filesystem::path(path).parent_path() / "serA007").string();
    meshioplusplus::RadiossAnimSeriesWriter writer(first);
    writer.Write(1.5, original);
    writer.Write(2.5, original);
    writer.Finalize();
    EXPECT_EQ(writer.NumSteps(), 2u);
    EXPECT_DOUBLE_EQ(
        meshioplusplus::read_radioss_anim(first).FieldData("meshio:time").As<double>()[0], 1.5);
    EXPECT_TRUE(std::filesystem::exists(
        std::filesystem::path(path).parent_path() / "serA008"));
    // A name that is not an animation file's is a stem.
    meshioplusplus::RadiossAnimSeriesWriter stem(
        (std::filesystem::path(path).parent_path() / "plain").string());
    stem.Write(0.0, original);
    EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(path).parent_path() / "plainA001"));
}

// A mesh of a brick, a quad and a beam with a point region and two time-history cell regions,
// one of which names elements of two families: the lists come back as the same regions.
TEST(RadiossAnimWriter, TimeHistoryRegionsRoundTrip) {
    const std::string path = write(state(), "thA001");
    Mesh mesh = meshioplusplus::read_radioss_anim(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);  // quad (cell 0) and tetra (cell 1)
    const auto indices = [](std::initializer_list<std::int64_t> Values) {
        meshioplusplus::NDArray a(meshioplusplus::DType::Int64, {Values.size()});
        std::copy(Values.begin(), Values.end(), a.As<std::int64_t>());
        return a;
    };
    using meshioplusplus::Region;
    using meshioplusplus::RegionKind;
    mesh.AddRegion(Region("radioss_th:probe", RegionKind::Point, -1, -1, indices({0, 3})));
    mesh.AddRegion(Region("radioss_th:watch", RegionKind::Cell, 3, -1, indices({0, 1})));

    const std::string out = (std::filesystem::path(path).parent_path() / "thoutA001").string();
    meshioplusplus::write_radioss_anim(out, mesh);
    const Mesh again = meshioplusplus::read_radioss_anim(out);
    bool probe = false, watch = false;
    for (std::size_t r = 0; r < again.NumRegions(); ++r) {
        const auto& reg = again.Region(r);
        if (reg.mName == "radioss_th:probe") {
            probe = reg.mKind == RegionKind::Point && reg.NumEntries() == 2 &&
                    detail::read_int(reg.mEntries, 0) == 0 &&
                    detail::read_int(reg.mEntries, 1) == 3;
        } else if (reg.mName == "radioss_th:watch") {
            watch = reg.mKind == RegionKind::Cell && reg.NumEntries() == 2;
        } else {
            EXPECT_NE(reg.mName.rfind("radioss_th:", 0), 0u) << reg.mName;
        }
    }
    EXPECT_TRUE(probe);
    EXPECT_TRUE(watch);
    // The time-history lists are not mistaken for parts: the two parts survive.
    std::size_t parts = 0;
    for (std::size_t r = 0; r < again.NumRegions(); ++r)
        parts += again.Region(r).mName.rfind("radioss_th:", 0) != 0 ? 1 : 0;
    EXPECT_EQ(parts, 2u);
}
