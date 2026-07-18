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

// External includes
#include <gtest/gtest.h>

// System includes
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/formats/tikz.hpp"

namespace {

// Read a whole file into a string.
std::string slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Count non-overlapping occurrences of `rNeedle` in `rHay`.
std::size_t count_occurrences(const std::string& rHay, const std::string& rNeedle) {
    if (rNeedle.empty())
        return 0;
    std::size_t n = 0;
    for (std::size_t pos = rHay.find(rNeedle); pos != std::string::npos;
         pos = rHay.find(rNeedle, pos + rNeedle.size()))
        ++n;
    return n;
}

// A genuinely non-flat 3-D triangle (one vertex off the z=0 plane).
mt::Mesh non_flat_mesh() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 1}}, "triangle", {{0, 1, 2}});
}

}  // namespace

// ---------------------------------------------------------------- SVG --------

TEST(Svg, OnePathPerDrawableCell) {
    struct Case {
        mt::Mesh mesh;
        std::size_t paths;
    };
    // tri_mesh: 2 triangles, quad_mesh: 2 quads, line_mesh: 5 lines.
    for (const auto& c : {Case{mt::tri_mesh(), 2}, Case{mt::tri_mesh_2d(), 2},
                          Case{mt::quad_mesh(), 2}, Case{mt::line_mesh(), 5}}) {
        std::string path = mt::temp_path(".svg");
        meshioplusplus::write_svg(path, c.mesh);
        std::string out = slurp(path);
        EXPECT_NE(out.find("<svg "), std::string::npos);
        EXPECT_NE(out.find("</svg>"), std::string::npos);
        EXPECT_EQ(count_occurrences(out, "<path "), c.paths);
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

TEST(Svg, TriangleAndQuadPathsAreClosed) {
    // triangle/quad paths end with "Z"; line paths do not.
    std::string tp = mt::temp_path(".svg");
    meshioplusplus::write_svg(tp, mt::tri_mesh());
    EXPECT_NE(slurp(tp).find("Z\""), std::string::npos);

    std::string lp = mt::temp_path(".svg");
    meshioplusplus::write_svg(lp, mt::line_mesh());
    EXPECT_EQ(slurp(lp).find("Z\""), std::string::npos);

    std::error_code ec;
    std::filesystem::remove(tp, ec);
    std::filesystem::remove(lp, ec);
}

TEST(Svg, StyleHonoursColours) {
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::tri_mesh(), ".3f", std::nullopt, 100.0, "#ff0000",
                              "#00ff00");
    std::string out = slurp(path);
    EXPECT_NE(out.find("fill: #ff0000"), std::string::npos);
    EXPECT_NE(out.find("stroke: #00ff00"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Svg, ImageWidthScalesCoordinates) {
    // Unit-square mesh, width 1. Default image_width=100 scales coords to 0..100.
    std::string scaled = mt::temp_path(".svg");
    meshioplusplus::write_svg(scaled, mt::tri_mesh());
    EXPECT_NE(slurp(scaled).find("100.000"), std::string::npos);

    // image_width = nullopt -> no scaling; coordinates stay at unit scale.
    std::string raw = mt::temp_path(".svg");
    meshioplusplus::write_svg(raw, mt::tri_mesh(), ".3f", std::nullopt, std::nullopt);
    std::string out = slurp(raw);
    EXPECT_EQ(out.find("100.000"), std::string::npos);
    EXPECT_NE(out.find("1.000"), std::string::npos);

    std::error_code ec;
    std::filesystem::remove(scaled, ec);
    std::filesystem::remove(raw, ec);
}

TEST(Svg, EmptyMeshIsValidWithNoPaths) {
    mt::Mesh empty;
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, empty);
    std::string out = slurp(path);
    EXPECT_NE(out.find("<svg "), std::string::npos);
    EXPECT_EQ(count_occurrences(out, "<path "), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Svg, SkipsUnsupportedCellTypes) {
    // A non-line/triangle/quad block (here a flat, z=0 tetra so the 2-D guard
    // still passes) must be silently dropped rather than drawn.
    mt::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "tetra", {{0, 1, 2, 3}});
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, m);
    EXPECT_EQ(count_occurrences(slurp(path), "<path "), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Svg, NonFlatShellMeshIsProjected) {
    // A non-flat 3-D surface mesh (no volume cells) is projected with the
    // default isometric camera and drawn as-is: 1 path for the 1 triangle.
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, non_flat_mesh());
    EXPECT_EQ(count_occurrences(slurp(path), "<path "), 1u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Svg, VolumeMeshDrawsSkin) {
    // A tetra mesh renders its extracted skin: 6 boundary triangles.
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::tet_mesh());
    EXPECT_EQ(count_occurrences(slurp(path), "<path "), 6u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Svg, CameraAnglesChangeProjection) {
    // A different camera azimuth must actually change the projected output.
    std::string p1 = mt::temp_path(".svg");
    std::string p2 = mt::temp_path(".svg");
    meshioplusplus::write_svg(p1, mt::hex_mesh());
    meshioplusplus::write_svg(p2, mt::hex_mesh(), ".3f", std::nullopt, 100.0, "#c8c5bd", "#000080",
                              /*azimuth=*/10.0, /*elevation=*/60.0);
    const std::string a = slurp(p1);
    const std::string b = slurp(p2);
    EXPECT_EQ(count_occurrences(a, "<path "), 6u);
    EXPECT_EQ(count_occurrences(b, "<path "), 6u);
    EXPECT_NE(a, b);
    std::error_code ec;
    std::filesystem::remove(p1, ec);
    std::filesystem::remove(p2, ec);
}

// --------------------------------------------------------------- TikZ --------

TEST(Tikz, StandaloneDocumentByDefault) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::tri_mesh());
    std::string out = slurp(path);
    EXPECT_NE(out.find("\\documentclass{standalone}"), std::string::npos);
    EXPECT_NE(out.find("\\usepackage{tikz}"), std::string::npos);
    EXPECT_NE(out.find("\\begin{tikzpicture}"), std::string::npos);
    EXPECT_NE(out.find("\\end{tikzpicture}"), std::string::npos);
    EXPECT_EQ(count_occurrences(out, "\\draw"), 2u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tikz, SnippetOmitsDocumentWrapper) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::tri_mesh(), ".6f", /*Standalone=*/false);
    std::string out = slurp(path);
    EXPECT_EQ(out.find("\\documentclass"), std::string::npos);
    EXPECT_NE(out.find("\\begin{tikzpicture}"), std::string::npos);
    EXPECT_EQ(count_occurrences(out, "\\draw"), 2u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tikz, FilledFacesVersusOpenLines) {
    // triangles/quads are filled + closed with cycle; lines are neither.
    std::string tp = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(tp, mt::tri_mesh());
    std::string tri = slurp(tp);
    EXPECT_NE(tri.find("fill=gray!30"), std::string::npos);
    EXPECT_NE(tri.find("-- cycle;"), std::string::npos);

    std::string lp = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(lp, mt::line_mesh());
    std::string line = slurp(lp);
    EXPECT_EQ(line.find("fill="), std::string::npos);
    EXPECT_EQ(line.find("cycle"), std::string::npos);
    EXPECT_EQ(count_occurrences(line, "\\draw"), 5u);

    std::error_code ec;
    std::filesystem::remove(tp, ec);
    std::filesystem::remove(lp, ec);
}

TEST(Tikz, OptionsAppearInOutput) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::tri_mesh(), ".6f", /*Standalone=*/true, "0.4pt", "blue!20",
                               "red", /*rScale=*/2.0);
    std::string out = slurp(path);
    EXPECT_NE(out.find("fill=blue!20"), std::string::npos);
    EXPECT_NE(out.find("draw=red"), std::string::npos);
    EXPECT_NE(out.find("line width=0.4pt"), std::string::npos);
    EXPECT_NE(out.find("[scale=2"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tikz, EmptyMeshHasNoDraws) {
    mt::Mesh empty;
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, empty);
    std::string out = slurp(path);
    EXPECT_NE(out.find("\\begin{tikzpicture}"), std::string::npos);
    EXPECT_EQ(count_occurrences(out, "\\draw"), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tikz, NonFlatShellMeshIsProjected) {
    // A non-flat 3-D surface mesh is projected and drawn as-is.
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, non_flat_mesh());
    std::string out = slurp(path);
    EXPECT_NE(out.find("\\documentclass{standalone}"), std::string::npos);
    EXPECT_EQ(count_occurrences(out, "\\draw"), 1u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tikz, VolumeMeshDrawsSkin) {
    // A hexahedron renders its skin: 6 filled quads, drawn back-to-front.
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::hex_mesh());
    std::string out = slurp(path);
    EXPECT_EQ(count_occurrences(out, "\\draw"), 6u);
    EXPECT_EQ(count_occurrences(out, "-- cycle;"), 6u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
