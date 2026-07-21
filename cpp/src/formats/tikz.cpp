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

// System includes
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/face_color.hpp"
#include "meshioplusplus/detail/projection.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {

namespace {

// Format a single double with a printf-style spec (spec without leading '%',
// e.g. ".6f").
std::string tikz_fmt_num(double value, const std::string& rSpec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), ("%" + rSpec).c_str(), value);
    return buf;
}

// TikZ's inline colour vocabulary. The braces are required: without them the
// commas inside would split the surrounding key=value option list.
std::string tikz_color_rgb(detail::Rgb c) {
    return "{rgb,255:red," + std::to_string(static_cast<int>(c.mR)) + ";green," +
           std::to_string(static_cast<int>(c.mG)) + ";blue," +
           std::to_string(static_cast<int>(c.mB)) + "}";
}

// The per-path option list for a filled face. Factored out so the constant-fill
// and the per-face-colour paths build it identically -- with coloring off this
// reproduces the previous inline string exactly.
std::string tikz_fill_style(const std::string& rFillColor, const std::string& rDraw,
                            const std::optional<std::string>& rLineWidth) {
    std::string s = "fill=" + rFillColor + ", draw=" + rDraw;
    if (rLineWidth.has_value())
        s += ", line width=" + *rLineWidth;
    return s;
}

// Colorbar geometry, shared by the flat and projected paths: a vertical
// gradient bar to the right of the drawing's bounding box, plus min/max labels.
// Emitted only when the caller asked for it, so the default output is untouched.
constexpr std::size_t kTikzColorbarSegments = 32;

void tikz_append_colorbar(std::vector<std::string>& rLines, const detail::FaceColors& rColors,
                          double MinX, double MaxX, double MinY, double MaxY,
                          const std::string& rFloatFmt) {
    const double w = MaxX - MinX;
    const double h = MaxY - MinY;
    const double gap = w * 0.08;
    const double bar_w = w * 0.06;
    const double x0 = MaxX + gap;
    const double x1 = x0 + bar_w;
    for (std::size_t i = 0; i < kTikzColorbarSegments; ++i) {
        const double y0 =
            MinY + h * (static_cast<double>(i) / static_cast<double>(kTikzColorbarSegments));
        const double y1 =
            MinY + h * (static_cast<double>(i + 1) / static_cast<double>(kTikzColorbarSegments));
        const detail::Rgb c = detail::colormap_lookup(
            rColors.mpTable,
            static_cast<double>(i) / static_cast<double>(kTikzColorbarSegments - 1));
        rLines.push_back("  \\fill[" + tikz_color_rgb(c) + "] (" + tikz_fmt_num(x0, rFloatFmt) +
                         "," + tikz_fmt_num(y0, rFloatFmt) + ") rectangle (" +
                         tikz_fmt_num(x1, rFloatFmt) + "," + tikz_fmt_num(y1, rFloatFmt) + ");");
    }
    rLines.push_back("  \\node[anchor=west, font=\\tiny] at (" + tikz_fmt_num(x1, rFloatFmt) + "," +
                     tikz_fmt_num(MinY, rFloatFmt) + ") {" +
                     tikz_fmt_num(rColors.mVMin, rFloatFmt) + "};");
    rLines.push_back("  \\node[anchor=west, font=\\tiny] at (" + tikz_fmt_num(x1, rFloatFmt) + "," +
                     tikz_fmt_num(MaxY, rFloatFmt) + ") {" +
                     tikz_fmt_num(rColors.mVMax, rFloatFmt) + "};");
}

// 3D rendering path: project the (already skin-extracted) surface mesh with
// the orthographic camera and emit its faces back-to-front as \draw commands.
// Mirrors the flat path's emission; no y-flip (TikZ is y-up).
void tikz_proj_write(const std::string& rPath, const Mesh& rSourceMesh, const Mesh& rDrawMesh,
                     const std::string& rFloatFmt, bool Standalone,
                     const std::optional<std::string>& rLineWidth, const std::string& rFill,
                     const std::string& rDraw, const std::optional<double>& rScale, double azimuth,
                     double elevation, double roll, const detail::ColorSpec& rColorSpec,
                     const std::string& rNanColor) {
    const detail::ProjectedSurface ps =
        detail::project_surface(rDrawMesh, azimuth, elevation, roll);

    const detail::FaceColors colors = detail::resolve_face_colors(
        rColorSpec, rSourceMesh, rDrawMesh, detail::color_faces_from_projection(ps.mFaces));

    const std::string fill_style = tikz_fill_style(rFill, rDraw, rLineWidth);
    std::string line_style = "draw=" + rDraw;
    if (rLineWidth.has_value())
        line_style += ", line width=" + *rLineWidth;

    std::vector<std::string> lines;
    for (std::size_t f = 0; f < ps.mFaces.size(); ++f) {
        const detail::ProjectedFace& face = ps.mFaces[f];
        std::string path;
        for (std::uint8_t k = 0; k < face.mNumNodes; ++k) {
            if (k)
                path += " -- ";
            const std::size_t p = static_cast<std::size_t>(face.mNodes[k]);
            path += "(" + tikz_fmt_num(ps.mX[p], rFloatFmt) + "," +
                    tikz_fmt_num(ps.mY[p], rFloatFmt) + ")";
        }
        if (face.mIsLine) {
            lines.push_back("  \\draw[" + line_style + "] " + path + ";");
        } else if (colors.mActive) {
            const std::optional<detail::Rgb> c = colors.Color(f);
            const std::string fill_color = c.has_value() ? tikz_color_rgb(*c) : rNanColor;
            lines.push_back("  \\draw[" + tikz_fill_style(fill_color, rDraw, rLineWidth) + "] " +
                            path + " -- cycle;");
        } else {
            lines.push_back("  \\draw[" + fill_style + "] " + path + " -- cycle;");
        }
    }

    if (colors.mActive && rColorSpec.mColorbar && !ps.mX.empty()) {
        double min_x = ps.mX[0], max_x = ps.mX[0], min_y = ps.mY[0], max_y = ps.mY[0];
        for (std::size_t i = 1; i < ps.mX.size(); ++i) {
            min_x = std::min(min_x, ps.mX[i]);
            max_x = std::max(max_x, ps.mX[i]);
            min_y = std::min(min_y, ps.mY[i]);
            max_y = std::max(max_y, ps.mY[i]);
        }
        tikz_append_colorbar(lines, colors, min_x, max_x, min_y, max_y, rFloatFmt);
    }

    std::string pic_opts;
    if (rScale.has_value()) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out;
    if (Standalone) {
        out.push_back("\\documentclass{standalone}");
        out.push_back("\\usepackage{tikz}");
        out.push_back("\\begin{document}");
    }
    out.push_back("\\begin{tikzpicture}" + pic_opt_str);
    for (const auto& l : lines)
        out.push_back(l);
    out.push_back("\\end{tikzpicture}");
    if (Standalone)
        out.push_back("\\end{document}");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

}  // namespace

void write_tikz(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
                bool Standalone, const std::optional<std::string>& rLineWidth,
                const std::string& rFill, const std::string& rDraw,
                const std::optional<double>& rScale, double azimuth, double elevation, double roll,
                const std::string& rColorBy, const std::optional<int>& rComponent,
                const std::string& rCmap, const std::optional<double>& rVMin,
                const std::optional<double>& rVMax, const std::string& rNanColor, bool Colorbar) {
    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    detail::ColorSpec spec;
    spec.mColorBy = rColorBy;
    spec.mComponent = rComponent;
    spec.mCmap = rCmap;
    spec.mVMin = rVMin;
    spec.mVMax = rVMax;
    spec.mColorbar = Colorbar;

    // A genuinely non-flat 3D mesh takes the projected-rendering path (skin
    // extraction for volume cells + orthographic camera); a flat one (every
    // z ~ 0) keeps the classic 2D path below, byte-identical to before.
    if (dim == 3) {
        for (std::size_t i = 0; i < num_points; ++i) {
            if (std::fabs(detail::read_double(points, i * dim + 2)) > 1.0e-14) {
                if (has_skinnable_cells(rMesh)) {
                    // Colouring by cell data needs to know which input cell
                    // each skin facet came from; that costs an extra Int64
                    // array, so only ask for it when it will actually be read.
                    const bool need_parents = spec.Active() && rMesh.HasCellData(rColorBy);
                    const Mesh skin =
                        detail::surface_extract(rMesh, /*forceFaceMode=*/true,
                                                /*linearize=*/true, need_parents, "extract_skin");
                    tikz_proj_write(rPath, rMesh, skin, rFloatFmt, Standalone, rLineWidth, rFill,
                                    rDraw, rScale, azimuth, elevation, roll, spec, rNanColor);
                } else {
                    tikz_proj_write(rPath, rMesh, rMesh, rFloatFmt, Standalone, rLineWidth, rFill,
                                    rDraw, rScale, azimuth, elevation, roll, spec, rNanColor);
                }
                return;
            }
        }
    }

    // TikZ/PGF uses the math convention (y-up), so — unlike SVG — no y-flip.
    auto coord = [&](std::int64_t p) {
        const std::size_t idx = static_cast<std::size_t>(p);
        const double px = (0 < dim) ? detail::read_double(points, idx * dim + 0) : 0.0;
        const double py = (1 < dim) ? detail::read_double(points, idx * dim + 1) : 0.0;
        return "(" + tikz_fmt_num(px, rFloatFmt) + "," + tikz_fmt_num(py, rFloatFmt) + ")";
    };

    // Per-path style option lists.
    const std::string fill_style = tikz_fill_style(rFill, rDraw, rLineWidth);
    std::string line_style = "draw=" + rDraw;
    if (rLineWidth.has_value())
        line_style += ", line width=" + *rLineWidth;

    // On the flat path the drawn mesh IS the source mesh, so a face's cell
    // index needs no provenance indirection.
    const detail::FaceColors colors =
        detail::resolve_face_colors(spec, rMesh, rMesh, detail::color_faces_flat(rMesh));

    std::vector<std::string> lines;
    std::size_t face_index = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad")
            continue;

        const NDArray& conn = cb.Conn();
        const std::size_t ncols = detail::cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            std::string path;
            for (std::size_t k = 0; k < ncols; ++k) {
                if (k)
                    path += " -- ";
                path += coord(detail::read_int(conn, r * ncols + k));
            }
            const std::size_t f = face_index++;
            if (type == "line") {
                lines.push_back("  \\draw[" + line_style + "] " + path + ";");
            } else if (colors.mActive) {
                const std::optional<detail::Rgb> c = colors.Color(f);
                const std::string fill_color = c.has_value() ? tikz_color_rgb(*c) : rNanColor;
                lines.push_back("  \\draw[" + tikz_fill_style(fill_color, rDraw, rLineWidth) +
                                "] " + path + " -- cycle;");
            } else {
                lines.push_back("  \\draw[" + fill_style + "] " + path + " -- cycle;");
            }
        }
    }

    if (colors.mActive && spec.mColorbar && num_points > 0) {
        double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
        for (std::size_t i = 0; i < num_points; ++i) {
            const double cx = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
            const double cy = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
            if (i == 0) {
                min_x = max_x = cx;
                min_y = max_y = cy;
            } else {
                min_x = std::min(min_x, cx);
                max_x = std::max(max_x, cx);
                min_y = std::min(min_y, cy);
                max_y = std::max(max_y, cy);
            }
        }
        tikz_append_colorbar(lines, colors, min_x, max_x, min_y, max_y, rFloatFmt);
    }

    // tikzpicture options (scale / line width) — emitted only when set.
    std::string pic_opts;
    if (rScale.has_value()) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out;
    if (Standalone) {
        out.push_back("\\documentclass{standalone}");
        out.push_back("\\usepackage{tikz}");
        out.push_back("\\begin{document}");
    }
    out.push_back("\\begin{tikzpicture}" + pic_opt_str);
    for (const auto& l : lines)
        out.push_back(l);
    out.push_back("\\end{tikzpicture}");
    if (Standalone)
        out.push_back("\\end{document}");

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

}  // namespace meshioplusplus
