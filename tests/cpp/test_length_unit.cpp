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
// Length units as metadata (roadmap 7.1.4): the shared key, and the UNV and
// Radioss readers and writers that fill and honour it.

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

std::string lu_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

Mesh lu_tets(double Unit) {
    Mesh m = mt::tet_mesh();
    if (Unit > 0.0)
        detail::set_mesh_length_unit(m, Unit);
    return m;
}

}  // namespace

TEST(LengthUnit, Accessors) {
    Mesh m = mt::tet_mesh();
    EXPECT_LT(detail::mesh_length_unit(m), 0.0);
    detail::set_mesh_length_unit(m, 1e-3);
    EXPECT_EQ(detail::mesh_length_unit(m), 1e-3);
    detail::rescale_mesh_length_unit(m, 1e-3);
    EXPECT_DOUBLE_EQ(detail::mesh_length_unit(m), 1.0);
    Mesh none = mt::tet_mesh();
    detail::rescale_mesh_length_unit(none, 1e-3);
    EXPECT_FALSE(none.HasFieldData(detail::kLengthUnitName));
}

TEST(LengthUnit, UnvWritesAndReadsIt) {
    for (const auto& [unit, code] : {std::pair<double, const char*>{1e-3, "         5"},
                                     {1.0, "         1"},
                                     {0.0254, "         9"}}) {
        const std::string path = mt::temp_path(".unv");
        registry_writers().at("unv")(path, lu_tets(unit));
        EXPECT_NE(lu_slurp(path).find(std::string("\n   164\n") + code), std::string::npos);
        const Mesh back = registry_readers().at("unv")(path);
        EXPECT_NEAR(detail::mesh_length_unit(back), unit, 1e-15 * unit);
        std::filesystem::remove(path);
    }
    const std::string path = mt::temp_path(".unv");
    registry_writers().at("unv")(path, lu_tets(-1.0));
    EXPECT_EQ(lu_slurp(path).find("\n   164\n"), std::string::npos);
    std::filesystem::remove(path);
}

TEST(LengthUnit, RadiossStatesTheWorkUnit) {
    const std::string dir = mt::temp_path("_rad");
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/x_0000.rad";
    registry_writers().at("radioss")(path, lu_tets(1e-3));
    const std::string text = lu_slurp(path);
    EXPECT_NE(text.find(std::string(18, ' ') + "kg" + std::string(18, ' ') + "mm"),
              std::string::npos);
    const Mesh back = registry_readers().at("radioss")(path);
    EXPECT_DOUBLE_EQ(detail::mesh_length_unit(back), 1e-3);
    std::filesystem::remove_all(dir);
}
