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
#include <algorithm>
#include <clocale>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "analytic_surfaces.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

std::string step_text(const Mesh& rMesh, const std::string& rName = "part.step") {
    const std::string path = (std::filesystem::temp_directory_path() / rName).string();
    write_step(path, rMesh);
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    std::filesystem::remove(path);
    return ss.str();
}

std::size_t count(const std::string& rText, const std::string& rEntity) {
    std::size_t n = 0;
    const std::string key = "=" + rEntity + "(";
    for (std::size_t at = rText.find(key); at != std::string::npos; at = rText.find(key, at + 1))
        ++n;
    return n;
}

}  // namespace

TEST(Step, CubeHasSixFacesTwelveEdgesEightVertices) {
    const std::string t = step_text(analytic::box_skin(1));
    EXPECT_EQ(t.rfind("ISO-10303-21;\n", 0), 0u);
    EXPECT_NE(t.find("END-ISO-10303-21;\n"), std::string::npos);
    EXPECT_EQ(count(t, "ADVANCED_FACE"), 6u);
    EXPECT_EQ(count(t, "PLANE"), 6u);
    EXPECT_EQ(count(t, "EDGE_CURVE"), 12u);
    EXPECT_EQ(count(t, "VERTEX_POINT"), 8u);
    EXPECT_EQ(count(t, "ORIENTED_EDGE"), 24u);
    EXPECT_EQ(count(t, "MANIFOLD_SOLID_BREP"), 1u);
    EXPECT_NE(t.find("FILE_NAME('part.step',"), std::string::npos);
}

TEST(Step, VolumeMeshWritesItsSkin) {
    const std::string t = step_text(mt::tet_mesh());
    EXPECT_EQ(count(t, "MANIFOLD_SOLID_BREP"), 1u);
    EXPECT_GE(count(t, "ADVANCED_FACE"), 4u);
}

TEST(Step, OpenSurfaceIsRefusedNamingItsEdges) {
    try {
        step_text(analytic::cylinder());
        FAIL() << "an open surface was written";
    } catch (const WriteError& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("64 open"), std::string::npos) << what;
        EXPECT_NE(what.find("first (0, 1)"), std::string::npos) << what;
    }
}

TEST(Step, InwardSkinIsWrittenReversedAndRaggedPolygonsWork) {
    // The cube's six faces as a ragged polygon block, wound inwards.
    const Mesh box = analytic::box_skin(1);
    const NDArray& c = box.Cells(0).Conn();
    std::vector<std::vector<std::int64_t>> rows;
    for (std::size_t i = 0; i < box.Cells(0).NumCells(); ++i)
        rows.push_back({c.As<std::int64_t>()[i * 4 + 3], c.As<std::int64_t>()[i * 4 + 2],
                        c.As<std::int64_t>()[i * 4 + 1], c.As<std::int64_t>()[i * 4]});
    Mesh m;
    m.AssignPoints(detail::data_owned_copy(box.Points()));
    m.AddPolygonBlock("polygon", rows);
    const std::string inward = step_text(m);
    const std::string outward = step_text(box);
    EXPECT_EQ(count(inward, "ADVANCED_FACE"), 6u);
    // Reversing the inward rings gives the outward cube's face loops back, so
    // the two files are the same bytes.
    EXPECT_EQ(inward, outward);
}

TEST(Step, OutputDoesNotDependOnTheLocale) {
    const std::string before = step_text(analytic::torus(3.0, 1.0, 8, 6));
    const char* old = std::setlocale(LC_NUMERIC, nullptr);
    const std::string saved = old ? old : "C";
    if (!std::setlocale(LC_NUMERIC, "de_DE.UTF-8"))
        GTEST_SKIP() << "no de_DE locale";
    const std::string after = step_text(analytic::torus(3.0, 1.0, 8, 6));
    std::setlocale(LC_NUMERIC, saved.c_str());
    EXPECT_EQ(before, after);
}
