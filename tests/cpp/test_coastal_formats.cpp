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
// ADCIRC fort.14 and TELEMAC Selafin (roadmap 2.4.3.3, 2.4.3.4): written and read back.

#include <cstdio>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/adcirc.hpp"
#include "meshioplusplus/formats/selafin.hpp"

TEST(Adcirc, AGridWithBoundariesRoundTrips) {
    const std::string p = mt::temp_path("_adcirc.14");
    {
        std::ofstream(p) << "grid\n2 4 = NE NP\n1 0 0 5\n2 1 0 6\n3 1 1 7\n4 0 1 8\n"
                            "1 3 1 2 3\n2 3 1 3 4\n1\n2\n2\n2\n3\n1\n3\n3 20\n3\n4\n1\n";
    }
    const auto mesh = meshioplusplus::read_adcirc(p);
    ASSERT_EQ(mesh.NumPoints(), 4u);
    ASSERT_EQ(mesh.NumRegions(), 2u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.PointData("depth"), 3), 8.0);
    const std::string q = mt::temp_path("_adcirc_out.14");
    meshioplusplus::write_adcirc(q, mesh);
    const auto back = meshioplusplus::read_adcirc(q);
    EXPECT_EQ(back.NumRegions(), 2u);
    EXPECT_EQ(
        meshioplusplus::detail::read_int(back.FieldData("adcirc:boundary:land_boundary_0"), 0), 2);
    EXPECT_THROW(meshioplusplus::write_adcirc(q, mt::tet_mesh()), meshioplusplus::WriteError);
    std::remove(p.c_str());
    std::remove(q.c_str());
}

TEST(Selafin, AWrittenFileReadsBack) {
    auto mesh = mt::tri_mesh();
    std::vector<double> h(mesh.NumPoints());
    for (std::size_t k = 0; k < h.size(); ++k)
        h[k] = 0.5 * static_cast<double>(k);
    mesh.AddPointData("WATER DEPTH", mt::data_array(h));
    const std::string p = mt::temp_path("_selafin.slf");
    meshioplusplus::write_selafin(p, mesh);
    const auto back = meshioplusplus::read_selafin(p);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.PointData("WATER DEPTH"), 1), 0.5);
    EXPECT_EQ(meshioplusplus::read_selafin_metadata(p, {}).mTimeValues.size(), 1u);
    EXPECT_THROW(meshioplusplus::write_selafin(p, mt::tet_mesh()), meshioplusplus::WriteError);
    std::remove(p.c_str());
}
