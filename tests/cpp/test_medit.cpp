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
// Binary Medit (.meshb) read and written by the C++ core directly (the Python suite,
// tests/python/test_medit.py, runs both engines on the fixtures).

// System includes
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;

Mesh lists_mesh(std::int64_t Ref) {
    Mesh mesh = mt::tri_mesh();
    meshioplusplus::NDArray pref(meshioplusplus::DType::Int64, {mesh.NumPoints()});
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        pref.As<std::int64_t>()[i] = static_cast<std::int64_t>(i) + 1;
    mesh.AddPointData("medit:ref", std::move(pref));
    meshioplusplus::NDArray cref(meshioplusplus::DType::Int64, {mesh.Cells(0).NumCells()});
    for (std::size_t i = 0; i < cref.Size(); ++i)
        cref.As<std::int64_t>()[i] = i == 0 ? Ref : 3;
    mesh.AppendCellData("medit:ref", std::move(cref));
    meshioplusplus::NDArray corners(meshioplusplus::DType::Int64, {std::size_t{2}});
    corners.As<std::int64_t>()[0] = 0;
    corners.As<std::int64_t>()[1] = 2;
    mesh.AddRegion(meshioplusplus::Region("medit:Corners", meshioplusplus::RegionKind::Point,
                                          std::move(corners)));
    return mesh;
}

}  // namespace

TEST(MeditBinary, RoundTripsRefsAndIdListsInBothWidths) {
    for (const std::int64_t ref : {std::int64_t{7}, std::int64_t{1} << 40}) {
        const Mesh mesh = lists_mesh(ref);
        const std::string path = mt::temp_path(".meshb");
        meshioplusplus::write_medit(path, mesh);
        std::ifstream in(path, std::ios::binary);
        char head[8];
        in.read(head, 8);
        EXPECT_EQ(head[4], ref > 1000 ? 4 : 3);
        const Mesh back = meshioplusplus::read_medit(path);
        ASSERT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).Type(), "triangle");
        ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
        EXPECT_EQ(back.CellData("medit:ref", 0).As<std::int64_t>()[0], ref);
        EXPECT_EQ(back.PointData("medit:ref").As<std::int64_t>()[2], 3);
        ASSERT_EQ(back.NumRegions(), 1u);
        EXPECT_EQ(back.Region(0).mName, "medit:Corners");
        EXPECT_EQ(back.Region(0).NumEntries(), 2u);
    }
}

TEST(MeditBinary, RegistryDispatchesOnTheSuffixAndRefusesDamage) {
    const Mesh mesh = lists_mesh(7);
    const std::string path = mt::temp_path(".meshb");
    meshioplusplus::registry_writers().at("medit")(path, mesh);
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "medit");
    EXPECT_EQ(meshioplusplus::registry_read(path, "medit", {}).NumPoints(), mesh.NumPoints());

    std::ifstream in(path, std::ios::binary);
    const std::string data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const std::size_t cut : {std::size_t{10}, std::size_t{40}, std::size_t{100}}) {
        const std::string cutpath = mt::temp_path(".meshb");
        std::ofstream(cutpath, std::ios::binary) << data.substr(0, cut);
        EXPECT_THROW(meshioplusplus::read_medit_binary(cutpath), meshioplusplus::ReadError);
    }
    std::string wrong = data;
    wrong[0] = 7;
    const std::string wrongpath = mt::temp_path(".meshb");
    std::ofstream(wrongpath, std::ios::binary) << wrong;
    EXPECT_THROW(meshioplusplus::read_medit_binary(wrongpath), meshioplusplus::ReadError);
}
