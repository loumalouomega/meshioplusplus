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
// DXF R2000 drawings (roadmap 7.6.4).

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/dxf.hpp"
#include "meshioplusplus/operations/fit_curves.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/registry.hpp"

using namespace meshioplusplus;

namespace {

constexpr double kDxfTestPi = 3.14159265358979323846;

std::string dxf_tmp(const std::string& rName) {
    return (std::filesystem::temp_directory_path() / ("mio_dxf_" + rName)).string();
}

std::string dxf_slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

void dxf_put(const std::string& rPath, const std::string& rText) {
    std::ofstream f(rPath, std::ios::binary);
    f << rText;
}

// A fit_curves sketch of a slot: two lines and two half circles of radius 0.75.
Mesh dxf_slot_sketch() {
    std::vector<std::vector<double>> pts;
    std::vector<std::vector<std::int64_t>> tris;
    std::vector<std::array<double, 2>> outline;
    for (int k = 0; k <= 24; ++k) {
        const double t = -kDxfTestPi / 2 + kDxfTestPi * k / 24;
        outline.push_back({2.0 + 0.75 * std::cos(t), 0.75 * std::sin(t)});
    }
    for (int k = 0; k <= 24; ++k) {
        const double t = kDxfTestPi / 2 + kDxfTestPi * k / 24;
        outline.push_back({-2.0 + 0.75 * std::cos(t), 0.75 * std::sin(t)});
    }
    const auto n = static_cast<std::int64_t>(outline.size());
    for (const auto& r_p : outline)
        pts.push_back({r_p[0], r_p[1], 0.0});
    for (const auto& r_p : outline)
        pts.push_back({r_p[0], r_p[1], 1.0});
    for (std::int64_t i = 0; i < n; ++i) {
        const std::int64_t j = (i + 1) % n;
        tris.push_back({i, j, j + n});
        tris.push_back({i, j + n, i + n});
    }
    SliceOptions so;
    so.mOrigin = {0.0, 0.0, 0.5};
    so.mChain = true;
    FitCurvesOptions fo;
    fo.mTolerance = 1e-6;
    return fit_curves(slice(mt::make_mesh(pts, "triangle", tris), so), fo).mMesh;
}

}  // namespace

TEST(Dxf, ALineMeshRoundTripsExactly) {
    Mesh m = mt::line_mesh();
    NDArray unit = NDArray::Uninit(DType::Float64, {1});
    unit.As<double>()[0] = 0.0254;
    m.AddFieldData("meshio:length_unit", std::move(unit));
    const std::string path = dxf_tmp("lines.dxf");
    write_dxf(path, m);
    const std::string text = dxf_slurp(path);
    EXPECT_NE(text.find("AC1015"), std::string::npos);
    EXPECT_NE(text.find("$INSUNITS\n 70\n1\n"), std::string::npos);  // inches
    EXPECT_EQ(text.rfind("  0\nEOF\n"), text.size() - 8);
    const Mesh back = read_dxf(path);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "line");
    EXPECT_EQ(back.Cells(0).NumCells(), 5u);
    EXPECT_EQ(back.NumPoints(), 4u);
    EXPECT_DOUBLE_EQ(back.FieldData("meshio:length_unit").As<double>()[0], 0.0254);
    EXPECT_TRUE(back.HasCellData("dxf:type"));
    std::filesystem::remove(path);
}

TEST(Dxf, ASketchIsWrittenAsLinesAndArcs) {
    const Mesh sketch = dxf_slot_sketch();
    const std::string path = dxf_tmp("slot.dxf");
    write_dxf(path, sketch);
    const std::string text = dxf_slurp(path);
    std::size_t lines = 0, arcs = 0;
    for (std::size_t p = text.find("\nLINE\n"); p != std::string::npos;
         p = text.find("\nLINE\n", p + 1))
        ++lines;
    for (std::size_t p = text.find("\nARC\n"); p != std::string::npos;
         p = text.find("\nARC\n", p + 1))
        ++arcs;
    EXPECT_EQ(lines, 2u);
    EXPECT_EQ(arcs, 2u);
    // The arcs come back sampled on their radius about (+-2, 0).
    const Mesh back = read_dxf(path);
    const double* p_xyz = back.Points().As<double>();
    for (std::size_t i = 0; i < back.NumPoints(); ++i) {
        const double x = p_xyz[3 * i], y = p_xyz[3 * i + 1];
        if (std::abs(x) <= 2.0)
            continue;  // a straight side
        EXPECT_NEAR(std::hypot(std::abs(x) - 2.0, y), 0.75, 1e-9);
    }
    std::filesystem::remove(path);
}

TEST(Dxf, TheRegistryReadsAndWritesIt) {
    EXPECT_TRUE(registry_readers().count("dxf"));
    EXPECT_TRUE(registry_writers().count("dxf"));
    bool found = false;
    for (const auto& [ext, fmt] : registry_extension_defaults())
        found = found || (ext == ".dxf" && fmt == "dxf");
    EXPECT_TRUE(found);
}

TEST(Dxf, RefusesWhatItDoesNotDraw) {
    const std::string path = dxf_tmp("text.dxf");
    dxf_put(path, "  0\nSECTION\n  2\nENTITIES\n  0\nTEXT\n  1\nhi\n  0\nENDSEC\n  0\nEOF\n");
    EXPECT_THROW(read_dxf(path), ReadError);
    dxf_put(path, "  0\nSECTION\n  2\nENTITIES\n  0\nCIRCLE\n 40\n-1\n  0\nENDSEC\n  0\nEOF\n");
    EXPECT_THROW(read_dxf(path), ReadError);
    dxf_put(path, "  0\nSECTION\n  2\nENTITIES\n  0\nLINE\n 10\nx\n");
    EXPECT_THROW(read_dxf(path), ReadError);
    EXPECT_THROW(write_dxf(path, mt::tri_mesh()), WriteError);
    std::filesystem::remove(path);
}
