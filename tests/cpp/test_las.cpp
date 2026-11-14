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
// LAS point clouds (roadmap 2.4.6.1): written and read back.

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/las.hpp"

TEST(Las, PointsAndFieldsRoundTrip) {
    auto mesh = mt::tri_mesh();
    std::vector<double> cls(mesh.NumPoints());
    for (std::size_t k = 0; k < cls.size(); ++k)
        cls[k] = static_cast<double>(k + 2);
    mesh.AddPointData("classification", mt::data_array(cls));
    const std::string p = mt::temp_path("_cloud.las");
    meshioplusplus::write_las(p, mesh);
    const auto back = meshioplusplus::read_las(p);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    EXPECT_EQ(std::string(back.Cells(0).Type()), "vertex");
    for (std::size_t k = 0; k < mesh.NumPoints() * 3; ++k)
        EXPECT_NEAR(meshioplusplus::detail::read_double(back.Points(), k),
                    meshioplusplus::detail::read_double(mesh.Points(), k), 1e-3);
    EXPECT_EQ(meshioplusplus::detail::read_int(back.PointData("classification"), 1), 3);
    meshioplusplus::LasReadOptions every;
    every.mEvery = 2;
    EXPECT_EQ(meshioplusplus::read_las(p, {}, every).NumPoints(), (mesh.NumPoints() + 1) / 2);
    std::remove(p.c_str());
}
