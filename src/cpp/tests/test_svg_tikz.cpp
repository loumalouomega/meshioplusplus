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

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/face_color.hpp"
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
