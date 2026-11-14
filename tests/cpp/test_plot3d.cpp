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
// PLOT3D (roadmap 2.4.2): the writer's unformatted form read back, and the ASCII
// and raw variants sized from their bytes.

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/plot3d.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

/// A 3x2x2 lattice block with its plot3d:dims.
meshioplusplus::Mesh plot3d_test_block() {
    using namespace meshioplusplus;
    Mesh m;
    NDArray pts(DType::Float64, {12, 3});
    std::size_t p = 0;
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i, ++p) {
                pts.As<double>()[p * 3 + 0] = i;
                pts.As<double>()[p * 3 + 1] = 0.5 * j;
                pts.As<double>()[p * 3 + 2] = 0.25 * k;
            }
    m.AssignPoints(std::move(pts));
    NDArray dims(DType::Int64, {1, 3});
    dims.As<std::int64_t>()[0] = 3;
    dims.As<std::int64_t>()[1] = 2;
    dims.As<std::int64_t>()[2] = 2;
    m.AddFieldData("plot3d:dims", std::move(dims));
    return m;
}

void plot3d_test_write(const std::string& rPath, const std::string& rBytes) {
    std::ofstream(rPath, std::ios::binary) << rBytes;
}

}  // namespace

TEST(Plot3d, TheUnformattedFormRoundTrips) {
    const std::string p = mt::temp_path("_plot3d.x");
    const auto mesh = plot3d_test_block();
    meshioplusplus::write_plot3d(p, mesh);
    const auto back = meshioplusplus::read_plot3d(p);
    ASSERT_EQ(back.NumPoints(), 12u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(back.Cells(0).Type()), "hexahedron");
    EXPECT_EQ(back.Cells(0).NumCells(), 2u);
    for (std::size_t k = 0; k < 36; ++k)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), k),
                         meshioplusplus::detail::read_double(mesh.Points(), k));
    EXPECT_EQ(back.NumRegions(), 1u);
    // The registry reaches it by extension too.
    EXPECT_EQ(
        meshioplusplus::registry_read(p, meshioplusplus::resolve_format(p, ""), {}).NumPoints(),
        12u);
    std::remove(p.c_str());
}

TEST(Plot3d, AsciiAndRawVariantsAreSized) {
    // A 2-D single block, ASCII: dims then x then y.
    const std::string a = mt::temp_path("_plot3d_ascii.p3d");
    plot3d_test_write(a, "2 2\n0 1 0 1\n0 0 1 1\n");
    const auto ascii = meshioplusplus::read_plot3d(a);
    ASSERT_EQ(ascii.NumPoints(), 4u);
    EXPECT_EQ(std::string(ascii.Cells(0).Type()), "quad");
    // The same, raw little-endian float32 with a block count: 4 + 8 + 4*2*4 bytes.
    std::string raw;
    const std::int32_t ints[] = {1, 2, 2};
    raw.append(reinterpret_cast<const char*>(ints), sizeof(ints));
    const float xy[] = {0, 1, 0, 1, 0, 0, 1, 1};
    raw.append(reinterpret_cast<const char*>(xy), sizeof(xy));
    const std::string r = mt::temp_path("_plot3d_raw.x");
    plot3d_test_write(r, raw);
    meshioplusplus::Plot3dReadOptions o;
    o.mPrecision = 4;
    const auto rawmesh = meshioplusplus::read_plot3d(r, o);
    EXPECT_EQ(rawmesh.NumPoints(), 4u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(rawmesh.Points(), 3), 1.0);
    std::remove(a.c_str());
    std::remove(r.c_str());
}

TEST(Plot3d, GarbageAndNonLatticesAreRefused) {
    const std::string g = mt::temp_path("_plot3d_bad.x");
    plot3d_test_write(g, std::string("\x01\x02\x03", 3));
    EXPECT_THROW(meshioplusplus::read_plot3d(g), meshioplusplus::ReadError);
    EXPECT_THROW(meshioplusplus::write_plot3d(g, mt::tri_mesh()), meshioplusplus::WriteError);
    std::remove(g.c_str());
}
