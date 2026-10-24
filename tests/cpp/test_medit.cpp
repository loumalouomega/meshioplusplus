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
#include "meshioplusplus/detail/value_io.hpp"
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

namespace {

/// The connectivity of a mesh's first block, row-major.
std::vector<std::int64_t> medit_test_conn(const Mesh& rMesh) {
    const auto cb = rMesh.Cells(0);
    std::vector<std::int64_t> out(cb.NumCells() * cb.NodesPerCell());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = meshioplusplus::detail::read_int(cb.Conn(), i);
    return out;
}

}  // namespace

TEST(MeditHighOrder, Triangle6AndTetra10RoundTripInBothEncodings) {
    for (const Mesh& mesh : {mt::triangle6_mesh(), mt::tet10_mesh()}) {
        for (const char* suffix : {".mesh", ".meshb"}) {
            const std::string path = mt::temp_path(suffix);
            meshioplusplus::write_medit(path, mesh);
            const Mesh back = meshioplusplus::read_medit(path);
            ASSERT_EQ(back.NumCellBlocks(), 1u) << suffix;
            EXPECT_EQ(back.Cells(0).Type(), mesh.Cells(0).Type()) << suffix;
            EXPECT_EQ(medit_test_conn(back), medit_test_conn(mesh)) << suffix;
        }
    }
}

TEST(MeditHighOrder, AFileSNodeOrderIsHonouredAndABadOneRefused) {
    // One TrianglesP3 element listed centre first; its ordering (after it in
    // the file) says so. Read back in meshio's order: corners, edge pairs, centre.
    const auto write = [](const std::string& rOrdering) {
        const std::string path = mt::temp_path(".mesh");
        std::ofstream out(path);
        out << "MeshVersionFormatted 2\nDimension 3\nVertices\n10\n";
        for (int i = 0; i < 10; ++i)
            out << i << " 0 0 1\n";
        out << "TrianglesP3\n1\n10 1 2 3 4 5 6 7 8 9 1\n" << rOrdering << "End\n";
        return path;
    };
    const std::string centre_first =
        "TrianglesP3Ordering\n10\n1 1 1\n3 0 0\n0 3 0\n0 0 3\n2 1 0\n1 2 0\n0 2 1\n0 1 2\n1 0 "
        "2\n2 0 1\n";
    EXPECT_EQ(medit_test_conn(meshioplusplus::read_medit(write(centre_first))),
              (std::vector<std::int64_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));
    // Without an ordering, meshio's order is assumed.
    EXPECT_EQ(medit_test_conn(meshioplusplus::read_medit(write(""))),
              (std::vector<std::int64_t>{9, 0, 1, 2, 3, 4, 5, 6, 7, 8}));
    // A tuple listed twice, and one no cubic triangle has, are refused.
    std::string twice = centre_first;
    twice.replace(twice.find("3 0 0"), 5, "0 3 0");
    EXPECT_THROW(meshioplusplus::read_medit(write(twice)), meshioplusplus::ReadError);
    std::string degree = centre_first;
    degree.replace(degree.find("1 1 1"), 5, "2 2 2");
    EXPECT_THROW(meshioplusplus::read_medit(write(degree)), meshioplusplus::ReadError);
}
