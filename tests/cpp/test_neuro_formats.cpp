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
// GIFTI and FreeSurfer surfaces (roadmap 2.4.4): written and read back.

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/freesurfer.hpp"
#include "meshioplusplus/formats/gifti.hpp"

TEST(Gifti, ASurfaceWithDataRoundTrips) {
    auto mesh = mt::tri_mesh();
    std::vector<double> v(mesh.NumPoints());
    for (std::size_t k = 0; k < v.size(); ++k)
        v[k] = 0.25 * static_cast<double>(k);
    mesh.AddPointData("thickness", mt::data_array(v));
    const std::string p = mt::temp_path("_surface.gii");
    meshioplusplus::write_gifti(p, mesh);
    const auto back = meshioplusplus::read_gifti(p);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.PointData("thickness"), 1), 0.25);
    EXPECT_TRUE(back.HasFieldData("gifti:transform"));
    EXPECT_THROW(meshioplusplus::write_gifti(p, mt::tet_mesh()), meshioplusplus::WriteError);
    std::remove(p.c_str());
}

TEST(Freesurfer, ASurfaceRoundTrips) {
    const auto mesh = mt::tri_mesh();
    const std::string p = mt::temp_path("_surface.pial");
    meshioplusplus::write_freesurfer(p, mesh);
    const auto back = meshioplusplus::read_freesurfer(p);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    EXPECT_EQ(back.Cells(0).NumCells(), mesh.Cells(0).NumCells());
    meshioplusplus::FreesurferReadOptions o;
    o.mCurv = {mt::temp_path("_missing.thickness")};
    EXPECT_THROW(meshioplusplus::read_freesurfer(p, o), meshioplusplus::ReadError);
    std::remove(p.c_str());
}
