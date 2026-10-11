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
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/face_color.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"

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
    // still passes) is dropped from the drawing -- but, unlike before, must
    // now warn rather than vanish with no diagnostic at all.
    mt::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "tetra", {{0, 1, 2, 3}});
    std::string path = mt::temp_path(".svg");
    testing::internal::CaptureStderr();
    meshioplusplus::write_svg(path, m);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(count_occurrences(slurp(path), "<path "), 0u);
    EXPECT_NE(err.find("tetra"), std::string::npos) << err;
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

TEST(Tikz, SkipsUnsupportedCellTypesWithAWarning) {
    // Twin of Svg.SkipsUnsupportedCellTypes: a non-line/triangle/quad block
    // is dropped from the drawing but must warn, not vanish silently.
    mt::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "tetra", {{0, 1, 2, 3}});
    std::string path = mt::temp_path(".tikz");
    testing::internal::CaptureStderr();
    meshioplusplus::write_tikz(path, m);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(count_occurrences(slurp(path), "\\draw"), 0u);
    EXPECT_NE(err.find("tetra"), std::string::npos) << err;
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

// ------------------------------------------------------- data-driven colour --

namespace {

// The viridis endpoints, as each format spells them.
constexpr const char* kViridisLoHex = "#440154";
constexpr const char* kViridisHiHex = "#fde725";
constexpr const char* kViridisLoTikz = "{rgb,255:red,68;green,1;blue,84}";
constexpr const char* kViridisHiTikz = "{rgb,255:red,253;green,231;blue,37}";

// A flat two-triangle strip carrying one cell value per triangle.
mt::Mesh two_tri_tagged(double a, double b) {
    mt::Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, "triangle",
                               {{0, 1, 2}, {1, 3, 2}});
    std::vector<meshioplusplus::NDArray> tag;
    tag.push_back(mt::data_array({a, b}));
    m.AddCellData("tag", std::move(tag));
    return m;
}

}  // namespace

TEST(SvgTikzColor, GoldenFillVocabulary) {
    // Pins the exact colour spelling in each format. The TikZ braces matter:
    // without them the inner commas would split the option list.
    std::string svg_path = mt::temp_path(".svg");
    meshioplusplus::write_svg(svg_path, two_tri_tagged(0.0, 1.0), ".3f", std::nullopt, 100.0,
                              "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0, "tag");
    const std::string svg_out = slurp(svg_path);
    EXPECT_NE(svg_out.find(std::string("fill=\"") + kViridisLoHex + "\""), std::string::npos);
    EXPECT_NE(svg_out.find(std::string("fill=\"") + kViridisHiHex + "\""), std::string::npos);

    std::string tikz_path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(tikz_path, two_tri_tagged(0.0, 1.0), ".6f", true, std::nullopt,
                               "gray!30", "black", std::nullopt, 45.0, 35.264389682754654, 0.0,
                               "tag");
    const std::string tikz_out = slurp(tikz_path);
    EXPECT_NE(tikz_out.find(std::string("fill=") + kViridisLoTikz + ", draw=black"),
              std::string::npos);
    EXPECT_NE(tikz_out.find(std::string("fill=") + kViridisHiTikz + ", draw=black"),
              std::string::npos);

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(tikz_path, ec);
}

TEST(SvgTikzColor, ColorByUnsetIsUnchanged) {
    // The regression guard: the colouring parameters at their defaults must
    // reproduce the pre-feature output byte for byte, on both writers.
    for (bool svg : {true, false}) {
        const std::string ext = svg ? ".svg" : ".tikz";
        std::string plain = mt::temp_path(ext);
        std::string explicit_defaults = mt::temp_path(ext);
        if (svg) {
            meshioplusplus::write_svg(plain, mt::data_mesh());
            // colorbar=true is deliberately passed: with no colouring there is
            // nothing to draw a bar for, so it must be ignored.
            meshioplusplus::write_svg(explicit_defaults, mt::data_mesh(), ".3f", std::nullopt,
                                      100.0, "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0,
                                      "", std::nullopt, "viridis", std::nullopt, std::nullopt,
                                      "#808080", true);
        } else {
            meshioplusplus::write_tikz(plain, mt::data_mesh());
            meshioplusplus::write_tikz(explicit_defaults, mt::data_mesh(), ".6f", true,
                                       std::nullopt, "gray!30", "black", std::nullopt, 45.0,
                                       35.264389682754654, 0.0, "", std::nullopt, "viridis",
                                       std::nullopt, std::nullopt, "gray", true);
        }
        EXPECT_EQ(slurp(plain), slurp(explicit_defaults));
        std::error_code ec;
        std::filesystem::remove(plain, ec);
        std::filesystem::remove(explicit_defaults, ec);
    }
}

TEST(SvgTikzColor, NonFiniteValuesUseNanColor) {
    // Every face is NaN, so every face takes nan_color and none is mapped.
    const double nan_v = std::nan("");
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, two_tri_tagged(nan_v, nan_v), ".3f", std::nullopt, 100.0,
                              "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0, "tag",
                              std::nullopt, "viridis", std::nullopt, std::nullopt, "#123456");
    const std::string out = slurp(path);
    EXPECT_EQ(count_occurrences(out, "fill=\"#123456\""), 2u);
    EXPECT_EQ(out.find(kViridisLoHex), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgTikzColor, ColorbarIsAppendOnly) {
    // The colorbar widens only the viewBox; every mesh <path> is untouched.
    std::string without = mt::temp_path(".svg");
    std::string with_bar = mt::temp_path(".svg");
    meshioplusplus::write_svg(without, two_tri_tagged(0.0, 1.0), ".3f", std::nullopt, 100.0,
                              "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0, "tag");
    meshioplusplus::write_svg(with_bar, two_tri_tagged(0.0, 1.0), ".3f", std::nullopt, 100.0,
                              "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0, "tag",
                              std::nullopt, "viridis", std::nullopt, std::nullopt, "#808080", true);
    const std::string a = slurp(without);
    const std::string b = slurp(with_bar);
    EXPECT_EQ(count_occurrences(a, "<rect "), 0u);
    EXPECT_EQ(count_occurrences(b, "<rect "), 32u);
    EXPECT_EQ(count_occurrences(b, "<text "), 2u);
    EXPECT_EQ(count_occurrences(a, "<path "), count_occurrences(b, "<path "));
    // The bar occupies extra viewBox width and nothing else.
    EXPECT_NE(a.find("viewBox=\"0.000 0.000 100.000 100.000\""), std::string::npos);
    EXPECT_NE(b.find("viewBox=\"0.000 0.000 144.000 100.000\""), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(without, ec);
    std::filesystem::remove(with_bar, ec);
}

TEST(SvgTikzColor, InvalidOptionsThrow) {
    const std::string path = mt::temp_path(".svg");
    // Unknown array name.
    EXPECT_THROW(
        meshioplusplus::write_svg(path, mt::data_mesh(), ".3f", std::nullopt, 100.0, "#c8c5bd",
                                  "#000080", 45.0, 35.264389682754654, 0.0, "nope"),
        std::invalid_argument);
    // Unknown colormap.
    EXPECT_THROW(meshioplusplus::write_svg(path, mt::data_mesh(), ".3f", std::nullopt, 100.0,
                                           "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0,
                                           "mat", std::nullopt, "nope"),
                 std::invalid_argument);
    // Component past the end of a scalar array.
    EXPECT_THROW(
        meshioplusplus::write_svg(path, mt::data_mesh(), ".3f", std::nullopt, 100.0, "#c8c5bd",
                                  "#000080", 45.0, 35.264389682754654, 0.0, "mat", 7),
        std::invalid_argument);
    // vmin above vmax.
    EXPECT_THROW(meshioplusplus::write_svg(path, mt::data_mesh(), ".3f", std::nullopt, 100.0,
                                           "#c8c5bd", "#000080", 45.0, 35.264389682754654, 0.0,
                                           "mat", std::nullopt, "viridis", 5.0, 1.0),
                 std::invalid_argument);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgTikzColor, CellDataLandsOnTheOwningSkinFacet) {
    // Two hexes tagged 0 and 1: each contributes 5 of the 10 boundary quads,
    // and every facet must carry ITS OWN hex's colour -- which is exactly what
    // the "surface:parent_cell" provenance is for.
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {2, 0, 0},
                                    {0, 1, 0},
                                    {1, 1, 0},
                                    {2, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {2, 0, 1},
                                    {0, 1, 1},
                                    {1, 1, 1},
                                    {2, 1, 1}}));
    m.AddCellBlock("hexahedron",
                   mt::conn_from({{0, 1, 4, 3, 6, 7, 10, 9}, {1, 2, 5, 4, 7, 8, 11, 10}}));
    std::vector<meshioplusplus::NDArray> tag;
    tag.push_back(mt::data_array({0.0, 1.0}));
    m.AddCellData("tag", std::move(tag));

    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, m, ".3f", std::nullopt, 100.0, "#c8c5bd", "#000080", 45.0,
                              35.264389682754654, 0.0, "tag");
    const std::string out = slurp(path);
    EXPECT_EQ(count_occurrences(out, "<path "), 10u);
    EXPECT_EQ(count_occurrences(out, std::string("fill=\"") + kViridisLoHex + "\""), 5u);
    EXPECT_EQ(count_occurrences(out, std::string("fill=\"") + kViridisHiHex + "\""), 5u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgTikzColor, PointDataAveragesOverCorners) {
    // Corner values 0, 3, 6 -> face value 3, dead centre of an explicit 0..6
    // range, so the face takes the colormap's midpoint colour.
    meshioplusplus::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    m.AddPointData("T", mt::data_array({0.0, 3.0, 6.0}));

    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, m, ".3f", std::nullopt, 100.0, "#c8c5bd", "#000080", 45.0,
                              35.264389682754654, 0.0, "T", std::nullopt, "viridis", 0.0, 6.0);
    const std::string out = slurp(path);
    const meshioplusplus::detail::Rgb mid = meshioplusplus::detail::colormap_lookup(
        meshioplusplus::detail::colormap_table("viridis"), 0.5);
    char buf[16];
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x", static_cast<int>(mid.mR),
                  static_cast<int>(mid.mG), static_cast<int>(mid.mB));
    EXPECT_NE(out.find(std::string("fill=\"") + buf + "\""), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Colormap, LookupClampsAndTableIsCanonical) {
    const std::uint8_t* viridis = meshioplusplus::detail::colormap_table("viridis");
    const meshioplusplus::detail::Rgb lo = meshioplusplus::detail::colormap_lookup(viridis, 0.0);
    const meshioplusplus::detail::Rgb hi = meshioplusplus::detail::colormap_lookup(viridis, 1.0);
    EXPECT_EQ(lo.mR, 68);
    EXPECT_EQ(lo.mG, 1);
    EXPECT_EQ(lo.mB, 84);
    EXPECT_EQ(hi.mR, 253);
    EXPECT_EQ(hi.mG, 231);
    EXPECT_EQ(hi.mB, 37);

    // Out-of-range and NaN clamp rather than reading past the table.
    const meshioplusplus::detail::Rgb under =
        meshioplusplus::detail::colormap_lookup(viridis, -5.0);
    const meshioplusplus::detail::Rgb over = meshioplusplus::detail::colormap_lookup(viridis, 5.0);
    const meshioplusplus::detail::Rgb nan_c =
        meshioplusplus::detail::colormap_lookup(viridis, std::nan(""));
    EXPECT_EQ(under.mR, lo.mR);
    EXPECT_EQ(over.mR, hi.mR);
    EXPECT_EQ(nan_c.mR, lo.mR);

    EXPECT_THROW(meshioplusplus::detail::colormap_table("nope"), std::invalid_argument);
}

TEST(Colormap, ColorParamHandlesADegenerateRange) {
    // An all-constant array gives vmin == vmax; the middle of the colormap is
    // the documented answer, not a division by zero.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::color_param(7.0, 7.0, 7.0), 0.5);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::color_param(-1.0, 0.0, 1.0), 0.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::color_param(2.0, 0.0, 1.0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::color_param(0.25, 0.0, 1.0), 0.25);
}

// ------------------------------------------------------------ readers --------

namespace {

void write_text(const std::string& rPath, const std::string& rText) {
    std::ofstream f(rPath, std::ios::binary);
    f << rText;
}

// A marked SVG document around `rBody` (the marker a planar write leaves).
std::string marked_svg(const std::string& rBody, const std::string& rExtraRoot = "") {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"0 0 10 10\"" +
           rExtraRoot +
           "><!--\nWritten by meshio++ v1.0.0\nDrawing: planar scale=1 offset=0,10\n-->" + rBody +
           "</svg>";
}

// A marked TikZ picture around `rBody`.
std::string marked_tikz(const std::string& rBody) {
    return "% Written by meshio++ v1.0.0\n% Drawing: planar\n\\begin{tikzpicture}\n" + rBody +
           "\n\\end{tikzpicture}\n";
}

// Reads `rText` as `rFormat` and returns the ReadError message ("" if none).
std::string read_error(const std::string& rFormat, const std::string& rText, bool Lenient = false) {
    const std::string path = mt::temp_path(rFormat == "svg" ? ".svg" : ".tikz");
    write_text(path, rText);
    std::string message;
    try {
        meshioplusplus::ReadOptions opts;
        opts.mLenient = Lenient;
        (void)(rFormat == "svg" ? meshioplusplus::read_svg(path, opts)
                                : meshioplusplus::read_tikz(path, opts));
    } catch (const meshioplusplus::ReadError& e) {
        message = e.what();
    }
    std::error_code ec;
    std::filesystem::remove(path, ec);
    return message;
}

double point_at(const mt::Mesh& rMesh, std::size_t i, std::size_t axis) {
    return meshioplusplus::detail::read_double(rMesh.Points(), i * rMesh.PointDim() + axis);
}

std::int64_t cell_node(const mt::Mesh& rMesh, std::size_t block, std::size_t cell,
                       std::size_t node) {
    const auto cb = rMesh.Cells(block);
    return meshioplusplus::detail::read_int(cb.Conn(), cell * cb.NodesPerCell() + node);
}

std::int64_t fill_at(const mt::Mesh& rMesh, const char* pName, std::size_t block, std::size_t cell,
                     std::size_t channel) {
    return meshioplusplus::detail::read_int(rMesh.CellData(pName, block), cell * 3 + channel);
}

}  // namespace

TEST(SvgRead, WritersMarkerNamesTheDrawing) {
    std::string flat = mt::temp_path(".svg");
    meshioplusplus::write_svg(flat, mt::tri_mesh());
    const std::string out = slurp(flat);
    EXPECT_NE(out.find("<!--\nWritten by meshio++ v"), std::string::npos);
    EXPECT_NE(out.find("\nDrawing: planar scale=100 offset=0,1\n-->"), std::string::npos);

    std::string proj = mt::temp_path(".svg");
    meshioplusplus::write_svg(proj, non_flat_mesh());
    EXPECT_NE(slurp(proj).find(
                  "\nDrawing: projected azimuth=45 elevation=35.264389682754654 roll=0 scale="),
              std::string::npos);
    std::error_code ec;
    std::filesystem::remove(flat, ec);
    std::filesystem::remove(proj, ec);
}

TEST(SvgRead, PlanarRoundTripRestoresModelUnits) {
    // 2 x 1 rectangle in two triangles: image_width 100 scales by 50 and the y-flip is
    // about 0 + 1; the marker's scale and offset bring both back (dyadic, so exactly).
    const mt::Mesh in =
        mt::make_mesh({{0, 0}, {2, 0}, {2, 1}, {0, 1}}, "triangle", {{0, 1, 2}, {0, 2, 3}});
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, in);
    const mt::Mesh out = meshioplusplus::read_svg(path);
    ASSERT_EQ(out.NumPoints(), 4u);
    EXPECT_EQ(out.PointDim(), 3u);
    for (std::size_t i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(point_at(out, i, 0), point_at(in, i, 0));
        EXPECT_DOUBLE_EQ(point_at(out, i, 1), point_at(in, i, 1));
        EXPECT_EQ(point_at(out, i, 2), 0.0);
    }
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(out.Cells(0).Type()), "triangle");
    EXPECT_EQ(mt::cell_rows(out), mt::cell_rows(in));
    EXPECT_FALSE(out.HasCellData("svg:fill"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, WeldsOnTheCoordinateText) {
    const mt::Mesh in = mt::tri_mesh();
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, in);
    EXPECT_EQ(meshioplusplus::read_svg(path).NumPoints(), 4u);
    // Two spellings of one number are two points: the weld is on the text.
    write_text(path,
               marked_svg("<path d=\"M 0 0L 1 0L 1 1Z\" /><path d=\"M 0.0 0L 1 0L 1 1Z\" />"));
    EXPECT_EQ(meshioplusplus::read_svg(path).NumPoints(), 4u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, QuadsLinesAndPolygonsKeepTheirType) {
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::quad_mesh());
    const mt::Mesh quads = meshioplusplus::read_svg(path);
    ASSERT_EQ(quads.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(quads.Cells(0).Type()), "quad");
    EXPECT_EQ(quads.Cells(0).NumCells(), 2u);

    meshioplusplus::write_svg(path, mt::line_mesh());
    const mt::Mesh lines = meshioplusplus::read_svg(path);
    ASSERT_EQ(lines.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(lines.Cells(0).Type()), "line");
    EXPECT_EQ(lines.Cells(0).NumCells(), 5u);

    write_text(path, marked_svg("<path d=\"M 0 0L 4 0L 5 2L 2 4L -1 2Z\" />"));
    const mt::Mesh polygon = meshioplusplus::read_svg(path);
    ASSERT_EQ(polygon.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(polygon.Cells(0).Type()), "polygon");
    EXPECT_EQ(polygon.Cells(0).NodesPerCell(), 5u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, ProjectedDrawingComesBackFlat) {
    // A 3-D mesh is drawn projected; what comes back is that picture, with z = 0.
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::tet_mesh());
    const mt::Mesh out = meshioplusplus::read_svg(path);
    EXPECT_GT(out.NumPoints(), 0u);
    EXPECT_GT(out.NumCellBlocks(), 0u);
    for (std::size_t i = 0; i < out.NumPoints(); ++i)
        EXPECT_EQ(point_at(out, i, 2), 0.0);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, FillsComeBackAsRgbAndTheColorbarIsSkipped) {
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, two_tri_tagged(0.0, 1.0), ".3f", std::nullopt, 100.0, "#c8c5bd",
                              "#000080", 45.0, 35.264389682754654, 0.0, "tag", std::nullopt,
                              "viridis", std::nullopt, std::nullopt, "#808080", /*Colorbar=*/true);
    const mt::Mesh out = meshioplusplus::read_svg(path);
    EXPECT_EQ(out.NumPoints(), 4u);  // the 32 swatches added no points
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).NumCells(), 2u);
    ASSERT_TRUE(out.HasCellData("svg:fill"));
    const std::int64_t lo[3] = {68, 1, 84}, hi[3] = {253, 231, 37};
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(fill_at(out, "svg:fill", 0, 0, k), lo[k]);
        EXPECT_EQ(fill_at(out, "svg:fill", 0, 1, k), hi[k]);
    }

    // The narrowing options drop the array.
    meshioplusplus::ReadOptions points_only;
    points_only.mPointsOnly = true;
    EXPECT_FALSE(meshioplusplus::read_svg(path, points_only).HasCellData("svg:fill"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, ACellWithoutAFillIsMinusOne) {
    mt::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}, "triangle", {{0, 1, 2}});
    m.AddCellBlock("line", mt::conn_from({{1, 3}}));
    std::vector<meshioplusplus::NDArray> tag;
    tag.push_back(mt::data_array({0.0}));
    tag.push_back(mt::data_array({1.0}));
    m.AddCellData("tag", std::move(tag));
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, m, ".3f", std::nullopt, 100.0, "#c8c5bd", "#000080", 45.0,
                              35.264389682754654, 0.0, "tag");
    const mt::Mesh out = meshioplusplus::read_svg(path);
    ASSERT_EQ(out.NumCellBlocks(), 2u);
    ASSERT_TRUE(out.HasCellData("svg:fill"));
    EXPECT_NE(fill_at(out, "svg:fill", 0, 0, 0), -1);
    for (std::size_t k = 0; k < 3; ++k)
        EXPECT_EQ(fill_at(out, "svg:fill", 1, 0, k), -1);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, AFileWithoutTheMarkerNeedsStrictFalse) {
    // No marker: refused, unless told to read it as a plain drawing, where the y-flip is
    // undone from the viewBox (here the whole 0..10 box, so y -> 10 - y).
    const std::string text =
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 10 10\">"
        "<path d=\"M 0 10L 10 10L 10 0Z\" /></svg>";
    EXPECT_NE(read_error("svg", text).find("no meshio++ drawing marker"), std::string::npos);
    std::string path = mt::temp_path(".svg");
    write_text(path, text);
    meshioplusplus::ReadOptions lenient;
    lenient.mLenient = true;
    const mt::Mesh out = meshioplusplus::read_svg(path, lenient);
    ASSERT_EQ(out.NumPoints(), 3u);
    EXPECT_EQ(point_at(out, 0, 0), 0.0);
    EXPECT_EQ(point_at(out, 0, 1), 0.0);
    EXPECT_EQ(point_at(out, 1, 0), 10.0);
    EXPECT_EQ(point_at(out, 1, 1), 0.0);
    EXPECT_EQ(point_at(out, 2, 0), 10.0);
    EXPECT_EQ(point_at(out, 2, 1), 10.0);
    // A marker with a Drawing line that is not ours is a ReadError even when lenient.
    const std::string bad =
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1 1\"><!--\nWritten by meshio++ "
        "v1\nDrawing: sideways\n--></svg>";
    EXPECT_NE(read_error("svg", bad, true).find("malformed drawing marker"), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, UnsupportedConstructsAreRefusedByName) {
    struct Case {
        std::string body;
        std::string root;
        std::string expect;
    };
    const std::vector<Case> cases = {
        {"<polygon points=\"0,0 1,0 1,1\" />", "", "<polygon>"},
        {"<path d=\"M 0 0L 1 0L 1 1Z\" transform=\"translate(1)\" />", "", "transform"},
        {"<path d=\"M 0 0C 1 0 1 1 0 1Z\" />", "", "curve command 'C'"},
        {"<path d=\"M 0 0Q 1 0 1 1Z\" />", "", "curve command 'Q'"},
        {"<path d=\"M 0 0A 1 1 0 0 1 1 1Z\" />", "", "curve command 'A'"},
        {"<path d=\"M 0 0l 1 0l 0 1z\" />", "", "relative path command 'l'"},
        {"<path d=\"M 0 0H 1V 1Z\" />", "", "path command 'H'"},
        {"<path d=\"M 0 0L 1 0Z M 2 2L 3 3\" />", "", "subpaths"},
        {"<path d=\"M 0 0 1 0 1 1Z\" />", "", "repeated coordinates"},
        {"<path d=\"M 0 0L 1\" />", "", "missing its coordinates"},
        {"<path />", "", "no d attribute"},
        {"<g><path d=\"M 0 0L 1 0L 1 1Z\" /></g>", "", "<g>"},
        {"<circle cx=\"1\" cy=\"1\" r=\"1\" />", "", "<circle>"},
        {"<path d=\"M 0 0L 1 0L 1 1Z\" />", " transform=\"scale(2)\"", "transform"},
        {"<path d=\"M 0 0L 1 0L 1 1Z\" />", " xmlns:xlink=\"http://www.w3.org/1999/xlink\"",
         "namespace prefixes"},
    };
    for (const Case& c : cases) {
        const std::string message = read_error("svg", marked_svg(c.body, c.root));
        EXPECT_NE(message.find(c.expect), std::string::npos)
            << "body: " << c.body << "\nmessage: " << message;
    }
    // An open path of three vertices is not a line.
    EXPECT_NE(
        read_error("svg", marked_svg("<path d=\"M 0 0L 1 0L 1 1\" />")).find("open path of 3"),
        std::string::npos);
    // A closed path of two is not a face.
    EXPECT_NE(read_error("svg", marked_svg("<path d=\"M 0 0L 1 0Z\" />")).find("closed path of 2"),
              std::string::npos);
}

TEST(SvgRead, MalformedInputIsAlwaysAReadError) {
    for (const std::string& text : std::vector<std::string>{
             "", "not xml", "<svg", "<html></html>", "<svg viewBox=\"0 0 1 1\"></svg>",
             marked_svg("<path d=\"M 1e999 0L 1 0L 1 1Z\" />"),
             "<svg xmlns=\"http://example.org/\" viewBox=\"0 0 1 1\"></svg>",
             "<svg xmlns=\"http://www.w3.org/2000/svg\"><!--\nWritten by meshio++ v1\nDrawing: "
             "planar scale=1 offset=0,0\n--></svg>",
             "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1\"><!--\nWritten by "
             "meshio++ v1\nDrawing: planar scale=1 offset=0,0\n--></svg>",
             "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1 1\"><!--\nWritten by "
             "meshio++ v1\nDrawing: planar scale=0 offset=0,0\n--></svg>",
             "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1 1\"><!--\nWritten by "
             "meshio++ v1\nDrawing: projected azimuth=1\n--></svg>"}) {
        EXPECT_FALSE(read_error("svg", text).empty()) << text;
    }
}

TEST(SvgRead, AnEmptyDrawingReadsAsAnEmptyMesh) {
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::Mesh());
    const mt::Mesh out = meshioplusplus::read_svg(path);
    EXPECT_EQ(out.NumPoints(), 0u);
    EXPECT_EQ(out.NumCellBlocks(), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(SvgRead, TheRegistryReadsItAndPassesTheOptions) {
    std::string path = mt::temp_path(".svg");
    meshioplusplus::write_svg(path, mt::tri_mesh());
    EXPECT_TRUE(meshioplusplus::registry_reader_supports_options("svg"));
    EXPECT_EQ(meshioplusplus::registry_read(path, "svg", meshioplusplus::ReadOptions()).NumPoints(),
              4u);
    EXPECT_EQ(meshioplusplus::registry_readers().at("svg")(path).NumPoints(), 4u);
    write_text(path, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1 1\"></svg>");
    EXPECT_THROW(meshioplusplus::registry_readers().at("svg")(path), meshioplusplus::ReadError);
    meshioplusplus::ReadOptions lenient;
    lenient.mLenient = true;
    EXPECT_NO_THROW(meshioplusplus::registry_read(path, "svg", lenient));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// ---- TikZ ----

TEST(TikzRead, WritersMarkerNamesTheDrawing) {
    std::string flat = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(flat, mt::tri_mesh());
    const std::string out = slurp(flat);
    EXPECT_EQ(out.rfind("% Written by meshio++ v", 0), 0u);
    EXPECT_NE(out.find("\n% Drawing: planar\n\\documentclass{standalone}"), std::string::npos);

    std::string proj = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(proj, non_flat_mesh(), ".6f", false);
    EXPECT_NE(slurp(proj).find("\n% Drawing: projected azimuth=45 elevation=35.264389682754654 "
                               "roll=0\n\\begin{tikzpicture}"),
              std::string::npos);
    std::error_code ec;
    std::filesystem::remove(flat, ec);
    std::filesystem::remove(proj, ec);
}

TEST(TikzRead, PlanarInputComesBackExactly) {
    // Verbatim coordinates: with a lossless float format the points are bit-identical,
    // whether or not the document is wrapped.
    const mt::Mesh in = mt::make_mesh({{0.1, 0.2}, {1.0 / 3.0, 0.0}, {0.7, 0.9}, {-0.25, 1e-7}},
                                      "triangle", {{0, 1, 2}, {0, 2, 3}});
    for (bool standalone : {true, false}) {
        std::string path = mt::temp_path(".tikz");
        meshioplusplus::write_tikz(path, in, ".17g", standalone);
        const mt::Mesh out = meshioplusplus::read_tikz(path);
        ASSERT_EQ(out.NumPoints(), 4u);
        for (std::size_t i = 0; i < 4; ++i) {
            EXPECT_EQ(point_at(out, i, 0), point_at(in, i, 0));
            EXPECT_EQ(point_at(out, i, 1), point_at(in, i, 1));
            EXPECT_EQ(point_at(out, i, 2), 0.0);
        }
        EXPECT_EQ(mt::cell_rows(out), mt::cell_rows(in));
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

TEST(TikzRead, TypesLinesAndTheWeld) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::quad_mesh());
    mt::Mesh out = meshioplusplus::read_tikz(path);
    EXPECT_EQ(out.NumPoints(), 6u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(out.Cells(0).Type()), "quad");

    meshioplusplus::write_tikz(path, mt::line_mesh());
    out = meshioplusplus::read_tikz(path);
    EXPECT_EQ(out.NumPoints(), 4u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(out.Cells(0).Type()), "line");
    EXPECT_EQ(out.Cells(0).NumCells(), 5u);

    write_text(path, marked_tikz("\\draw[fill=gray!30, draw=black] (0,0) -- (4,0) -- (5,2) -- "
                                 "(2,4) -- (-1,2) -- cycle;"));
    out = meshioplusplus::read_tikz(path);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(out.Cells(0).Type()), "polygon");
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(TikzRead, FillsComeBackAsRgbAndTheColorbarIsSkipped) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, two_tri_tagged(0.0, 1.0), ".6f", true, std::nullopt, "gray!30",
                               "black", std::nullopt, 45.0, 35.264389682754654, 0.0, "tag",
                               std::nullopt, "viridis", std::nullopt, std::nullopt, "gray",
                               /*Colorbar=*/true);
    const mt::Mesh out = meshioplusplus::read_tikz(path);
    EXPECT_EQ(out.NumPoints(), 4u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    ASSERT_TRUE(out.HasCellData("tikz:fill"));
    EXPECT_EQ(fill_at(out, "tikz:fill", 0, 0, 0), 68);
    EXPECT_EQ(fill_at(out, "tikz:fill", 0, 1, 2), 37);
    // A flat fill that is not an RGB literal (gray!30) gives no array.
    meshioplusplus::write_tikz(path, mt::tri_mesh());
    EXPECT_FALSE(meshioplusplus::read_tikz(path).HasCellData("tikz:fill"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(TikzRead, ProjectedDrawingComesBackFlat) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::tet_mesh());
    const mt::Mesh out = meshioplusplus::read_tikz(path);
    EXPECT_GT(out.NumPoints(), 0u);
    for (std::size_t i = 0; i < out.NumPoints(); ++i)
        EXPECT_EQ(point_at(out, i, 2), 0.0);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(TikzRead, AFileWithoutTheMarkerNeedsStrictFalse) {
    const std::string text =
        "\\begin{tikzpicture}\n\\draw (0,0) -- (1,0) -- (1,1) -- cycle;\n\\end{tikzpicture}\n";
    EXPECT_NE(read_error("tikz", text).find("no meshio++ drawing marker"), std::string::npos);
    std::string path = mt::temp_path(".tikz");
    write_text(path, text);
    meshioplusplus::ReadOptions lenient;
    lenient.mLenient = true;
    const mt::Mesh out = meshioplusplus::read_tikz(path, lenient);
    ASSERT_EQ(out.NumPoints(), 3u);
    EXPECT_EQ(point_at(out, 2, 1), 1.0);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(TikzRead, UnsupportedConstructsAreRefusedByName) {
    struct Case {
        std::string body;
        std::string expect;
    };
    const std::vector<Case> cases = {
        {"\\path[draw] (0,0) -- (1,1);", "\\path"},
        {"\\foreach \\x in {0,1} \\draw (\\x,0) -- (1,1);", "\\foreach"},
        {"\\draw (0,0) -- (1,1) plot (2,2);", "'plot"},
        {"\\draw (0,0) circle (1);", "'circle"},
        {"\\draw (0,0) .. controls (1,1) and (2,2) .. (3,3);", "'.."},
        {"\\draw (0,0) -- ++(1,0);", "relative"},
        {"\\draw (0,0) -- +(1,0);", "relative"},
        {"\\draw (a) -- (b);", "plain (x,y)"},
        {"\\draw (30:1cm) -- (0,0);", "plain (x,y)"},
        {"\\draw (0,0,0) -- (1,1,1);", "plain (x,y)"},
        {"\\draw (0cm,0cm) -- (1,1);", "plain (x,y)"},
        {"\\draw[->] (0,0) -- (1,1);", "option '->'"},
        {"\\draw[thick] (0,0) -- (1,1);", "option 'thick'"},
        {"\\draw (0,0) -- (1,1) to (2,2);", "'to"},
        {"\\fill (0,0) -- (1,0) -- (1,1) -- cycle;", "colour-bar rectangle"},
        {"\\draw (0,0) -- (1,1)", "unsupported path construct"},
        {"\\draw (0,0) -- (1,1) -- cycle -- (2,2);", "after cycle"},
        {"\\draw (0,0) -- (1,1) -- (2,2);", "open path of 3"},
        {"\\draw (0,0) -- (1,0) -- cycle;", "closed path of 2"},
        {"\\tikzset{x=1cm}", "\\tikzset"},
        {"\\begin{scope}", "'scope'"},
        {"text", "unexpected text"},
    };
    for (const Case& c : cases) {
        const std::string message = read_error("tikz", marked_tikz(c.body));
        EXPECT_NE(message.find(c.expect), std::string::npos)
            << "body: " << c.body << "\nmessage: " << message;
    }
    EXPECT_NE(
        read_error("tikz", "% Written by meshio++ v1\n% Drawing: planar\n").find("no tikzpicture"),
        std::string::npos);
    EXPECT_NE(read_error("tikz",
                         "% Written by meshio++ v1\n% Drawing: "
                         "planar\n\\begin{tikzpicture}[x=1cm]\n\\end{tikzpicture}\n")
                  .find("option 'x'"),
              std::string::npos);
    EXPECT_NE(read_error("tikz",
                         "% Written by meshio++ v1\n% Drawing: "
                         "planar\n\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n")
                  .find("not closed"),
              std::string::npos);
    EXPECT_NE(read_error("tikz",
                         "% Written by meshio++ v1\n% Drawing: projected "
                         "azimuth=1\n\\begin{tikzpicture}\n\\end{tikzpicture}\n")
                  .find("malformed drawing marker"),
              std::string::npos);
}

TEST(TikzRead, CommentsAndEscapesAreHandled) {
    std::string path = mt::temp_path(".tikz");
    write_text(path,
               "% Written by meshio++ v1.0.0\n% Drawing: planar\n\\begin{tikzpicture} % open\n"
               "  \\draw (0,0) -- (1,0) % a comment\n    -- (1,1) -- cycle; % 100\\%\n"
               "\\end{tikzpicture}\n");
    const mt::Mesh out = meshioplusplus::read_tikz(path);
    EXPECT_EQ(out.NumPoints(), 3u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(TikzRead, TheRegistryReadsItAndPassesTheOptions) {
    std::string path = mt::temp_path(".tikz");
    meshioplusplus::write_tikz(path, mt::tri_mesh());
    EXPECT_TRUE(meshioplusplus::registry_reader_supports_options("tikz"));
    EXPECT_EQ(meshioplusplus::registry_readers().at("tikz")(path).NumPoints(), 4u);
    write_text(path, "\\begin{tikzpicture}\\end{tikzpicture}");
    EXPECT_THROW(meshioplusplus::registry_readers().at("tikz")(path), meshioplusplus::ReadError);
    meshioplusplus::ReadOptions lenient;
    lenient.mLenient = true;
    EXPECT_NO_THROW(meshioplusplus::registry_read(path, "tikz", lenient));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
