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
// UGRID NetCDF (roadmap 2.4.3.2): a mixed mesh written and read back.

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ugrid_nc.hpp"

TEST(UgridNc, AMixedMeshRoundTrips) {
    using namespace meshioplusplus;
    Mesh m;
    NDArray pts(DType::Float64, {5, 3});
    const double xyz[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 2, 0.5, 0};
    std::copy(std::begin(xyz), std::end(xyz), pts.As<double>());
    m.AssignPoints(std::move(pts));
    NDArray tri(DType::Int64, {1, 3});
    const std::int64_t t[] = {1, 4, 2};
    std::copy(std::begin(t), std::end(t), tri.As<std::int64_t>());
    m.AddCellBlock("triangle", std::move(tri));
    NDArray quad(DType::Int64, {1, 4});
    for (std::int64_t k = 0; k < 4; ++k)
        quad.As<std::int64_t>()[k] = k;
    m.AddCellBlock("quad", std::move(quad));
    m.AddPointData("h", mt::data_array({0, 1, 2, 3, 4}));
    m.AddCellData("c", {mt::data_array({1}), mt::data_array({2})});

    const std::string p = mt::temp_path("_ugrid.nc");
    write_ugrid_nc(p, m);
    const Mesh back = read_ugrid_nc(p);
    ASSERT_EQ(back.NumPoints(), 5u);
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(std::string(back.Cells(0).Type()), "triangle");
    EXPECT_EQ(std::string(back.Cells(1).Type()), "quad");
    EXPECT_DOUBLE_EQ(detail::read_double(back.PointData("h"), 4), 4.0);
    EXPECT_DOUBLE_EQ(detail::read_double(back.CellData("c", 1), 0), 2.0);
    EXPECT_TRUE(read_ugrid_nc_metadata(p, {}).mTimeValues.empty());
    std::remove(p.c_str());
}

TEST(UgridNc, VolumeCellsAreRefused) {
    const std::string p = mt::temp_path("_ugrid_tet.nc");
    EXPECT_THROW(meshioplusplus::write_ugrid_nc(p, mt::tet_mesh()), meshioplusplus::WriteError);
    std::remove(p.c_str());
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
