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
#include <array>
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
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/typed_view.hpp"
#include "../detail/drawing.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"

namespace meshioplusplus {

namespace {

// Format a single double with a printf-style spec (spec without leading '%',
// e.g. ".6f").
std::string tikz_fmt_num(double value, const std::string& rSpec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), ("%" + rSpec).c_str(), value);
    return buf;
}

// The marker comment that opens the file: the provenance tag, the
// `Drawing: ...` line, then whatever else the provenance block carries, each as a
// `%` line (the twin of `_marker_lines` in `_tikz.py`).
std::vector<std::string> tikz_marker_lines(bool Projected, double Azimuth, double Elevation,
                                           double Roll) {
    std::vector<std::string> lines;
    for (const std::string& line : detail::drawing_comment_lines(
             detail::drawing_marker_line(Projected, Azimuth, Elevation, Roll, nullptr)))
        lines.push_back("% " + line);
    return lines;
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
        detail::snprintf_c(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out = tikz_marker_lines(true, azimuth, elevation, roll);
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

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

// ---------------------------------------------------------------------------
// Reader. Inverts exactly what the writer above emits, and nothing more: see
// doc/formats/tikz.md ("Reading") and detail/drawing.hpp.
// ---------------------------------------------------------------------------

// A name for the construct at the cursor, for a refusal: the next run of
// non-blank characters, cut short.
std::string tikzr_construct_at(std::string_view Text, std::size_t Pos) {
    std::size_t end = Pos;
    while (end < Text.size() && !detail::text_is_blank(Text[end]) && end - Pos < 16)
        ++end;
    return detail::drawing_safe(Text.substr(Pos, end - Pos));
}

// Splits an option list on the commas outside braces.
std::vector<std::string_view> tikzr_split_options(std::string_view Text) {
    std::vector<std::string_view> items;
    int depth = 0;
    std::size_t start = 0;
    for (std::size_t i = 0; i < Text.size(); ++i) {
        if (Text[i] == '{') {
            ++depth;
        } else if (Text[i] == '}') {
            --depth;
        } else if (Text[i] == ',' && depth == 0) {
            items.push_back(Text.substr(start, i - start));
            start = i + 1;
        }
    }
    items.push_back(Text.substr(start));
    return items;
}

std::string_view tikzr_trim(std::string_view Text) {
    while (!Text.empty() && detail::text_is_blank(Text.front()))
        Text.remove_prefix(1);
    while (!Text.empty() && detail::text_is_blank(Text.back()))
        Text.remove_suffix(1);
    return Text;
}

// `{rgb,255:red,R;green,G;blue,B}` -> {R, G, B}; anything else (gray!30, ...) -> nullopt.
std::optional<std::array<int, 3>> tikzr_rgb(std::string_view Value) {
    Value = tikzr_trim(Value);
    if (Value.size() < 2 || Value.front() != '{' || Value.back() != '}')
        return std::nullopt;
    Value = Value.substr(1, Value.size() - 2);
    constexpr std::string_view head = "rgb,255:red,";
    if (Value.substr(0, head.size()) != head)
        return std::nullopt;
    Value.remove_prefix(head.size());
    std::array<int, 3> rgb{};
    const std::array<std::string_view, 3> next = {";green,", ";blue,", ""};
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t end = 0;
        while (end < Value.size() && Value[end] >= '0' && Value[end] <= '9')
            ++end;
        if (end == 0 || end > 3)
            return std::nullopt;
        int v = 0;
        for (std::size_t i = 0; i < end; ++i)
            v = v * 10 + (Value[i] - '0');
        if (v > 255)
            return std::nullopt;
        rgb[k] = v;
        Value.remove_prefix(end);
        if (Value.substr(0, next[k].size()) != next[k])
            return std::nullopt;
        Value.remove_prefix(next[k].size());
    }
    if (!Value.empty())
        return std::nullopt;
    return rgb;
}

// One reading of a figure: the text with its comments cut, and a cursor.
class TikzrReader {
public:
    TikzrReader(std::string_view Body, detail::DrawingBuilder& rBuilder)
        : mText(Body), mBuilder(rBuilder) {}

    void Run() {
        bool in_document = false;
        bool in_picture = false;
        bool seen_picture = false;
        while (true) {
            SkipBlanks();
            if (mPos >= mText.size())
                break;
            if (mText[mPos] != '\\')
                throw ReadError("TikZ: unexpected text '" + tikzr_construct_at(mText, mPos) + "'");
            ++mPos;
            std::size_t end = mPos;
            while (end < mText.size() && ((mText[end] >= 'a' && mText[end] <= 'z') ||
                                          (mText[end] >= 'A' && mText[end] <= 'Z')))
                ++end;
            const std::string_view name = mText.substr(mPos, end - mPos);
            if (name.empty())
                throw ReadError("TikZ: unexpected control sequence '\\" +
                                tikzr_construct_at(mText, mPos) + "'");
            mPos = end;
            if (name == "documentclass" || name == "usepackage") {
                SkipBlanks();
                if (mPos < mText.size() && mText[mPos] == '[')
                    Balanced('[', ']');
                SkipBlanks();
                Balanced('{', '}');
            } else if (name == "begin" || name == "end") {
                SkipBlanks();
                const std::string_view env = Balanced('{', '}');
                if (env == "document") {
                    if (name == "begin" ? in_document : !in_document)
                        throw ReadError("TikZ: unbalanced document environment");
                    in_document = name == "begin";
                } else if (env == "tikzpicture") {
                    if (name == "begin") {
                        if (in_picture || seen_picture)
                            throw ReadError("TikZ: only one tikzpicture is read");
                        SkipBlanks();
                        if (mPos < mText.size() && mText[mPos] == '[')
                            PictureOptions(Balanced('[', ']'));
                        in_picture = seen_picture = true;
                    } else {
                        if (!in_picture)
                            throw ReadError("TikZ: \\end{tikzpicture} without its begin");
                        in_picture = false;
                    }
                } else {
                    throw ReadError("TikZ: the environment '" + detail::drawing_safe(env) +
                                    "' is not supported");
                }
            } else if (name == "draw") {
                RequirePicture(in_picture, name);
                Draw();
            } else if (name == "fill") {
                RequirePicture(in_picture, name);
                ColorBarSwatch();
            } else if (name == "node") {
                RequirePicture(in_picture, name);
                SkipStatement();
            } else {
                throw ReadError("TikZ: the command \\" + std::string(name) + " is not supported");
            }
        }
        if (in_picture)
            throw ReadError("TikZ: the tikzpicture is not closed");
        if (!seen_picture)
            throw ReadError("TikZ: no tikzpicture");
    }

private:
    void SkipBlanks() {
        while (mPos < mText.size() && detail::text_is_blank(mText[mPos]))
            ++mPos;
    }

    static void RequirePicture(bool InPicture, std::string_view Name) {
        if (!InPicture)
            throw ReadError("TikZ: \\" + std::string(Name) + " outside a tikzpicture");
    }

    // The text between the bracket at the cursor and its match (braces nest and
    // hide brackets); the cursor ends just past the closing one.
    std::string_view Balanced(char Open, char Close) {
        if (mPos >= mText.size() || mText[mPos] != Open)
            throw ReadError(std::string("TikZ: expected '") + Open + "'");
        const std::size_t start = mPos + 1;
        int braces = 0;
        for (std::size_t i = start; i < mText.size(); ++i) {
            const char c = mText[i];
            if (c == '{') {
                ++braces;
            } else if (c == '}' && braces > 0) {
                --braces;
            } else if (c == Close && braces == 0) {
                mPos = i + 1;
                return mText.substr(start, i - start);
            }
        }
        throw ReadError(std::string("TikZ: '") + Open + "' is never closed");
    }

    // `scale=` and `line width=` only: anything else changes the drawing.
    static void PictureOptions(std::string_view Options) {
        for (std::string_view item : tikzr_split_options(Options)) {
            item = tikzr_trim(item);
            if (item.empty())
                continue;
            const std::size_t eq = item.find('=');
            const std::string_view key = tikzr_trim(item.substr(0, eq));
            if (key != "scale" && key != "line width")
                throw ReadError("TikZ: the tikzpicture option '" + detail::drawing_safe(key) +
                                "' is not supported");
        }
    }

    // `fill`, `draw` and `line width`; the fill is returned when it is an RGB one.
    static std::optional<std::array<int, 3>> DrawOptions(std::string_view Options) {
        std::optional<std::array<int, 3>> fill;
        for (std::string_view item : tikzr_split_options(Options)) {
            item = tikzr_trim(item);
            if (item.empty())
                continue;
            const std::size_t eq = item.find('=');
            const std::string_view key = tikzr_trim(item.substr(0, eq));
            if (key != "fill" && key != "draw" && key != "line width")
                throw ReadError("TikZ: the \\draw option '" + detail::drawing_safe(key) +
                                "' is not supported");
            if (key == "fill" && eq != std::string_view::npos)
                fill = tikzr_rgb(item.substr(eq + 1));
        }
        return fill;
    }

    // `(x,y)` -> the welded vertex; x and y are plain numbers, nothing else.
    // @p Keep false (a colour-bar corner): checked, not welded into the mesh.
    std::int64_t Coordinate(bool Keep = true) {
        SkipBlanks();
        if (mPos >= mText.size())
            throw ReadError("TikZ: a path ends where a coordinate is expected");
        if (mText[mPos] == '+')
            throw ReadError("TikZ: relative coordinates (+, ++) are not supported");
        if (mText[mPos] != '(')
            throw ReadError("TikZ: unsupported path construct '" + tikzr_construct_at(mText, mPos) +
                            "'");
        const std::size_t close = mText.find(')', mPos);
        if (close == std::string_view::npos)
            throw ReadError("TikZ: a coordinate is never closed");
        const std::string_view inside = mText.substr(mPos + 1, close - mPos - 1);
        mPos = close + 1;
        const std::size_t comma = inside.find(',');
        const std::string_view x = tikzr_trim(inside.substr(0, comma));
        const std::string_view y = comma == std::string_view::npos
                                       ? std::string_view()
                                       : tikzr_trim(inside.substr(comma + 1));
        if (comma == std::string_view::npos || x.empty() || y.empty() ||
            detail::drawing_scan_number(x, 0) != x.size() ||
            detail::drawing_scan_number(y, 0) != y.size())
            throw ReadError("TikZ: the coordinate (" + detail::drawing_safe(inside) +
                            ") is not a plain (x,y) pair of numbers (named, polar, calc, 3-D and "
                            "unit-carrying coordinates are not supported)");
        if (!Keep) {
            detail::drawing_number_value(x, "TikZ");
            detail::drawing_number_value(y, "TikZ");
            return -1;
        }
        return mBuilder.Weld(x, y, [&] {
            return std::make_pair(detail::drawing_number_value(x, "TikZ"),
                                  detail::drawing_number_value(y, "TikZ"));
        });
    }

    // `\draw[opts] (x,y) -- (x,y) ... [-- cycle];`
    void Draw() {
        SkipBlanks();
        std::optional<std::array<int, 3>> fill;
        if (mPos < mText.size() && mText[mPos] == '[')
            fill = DrawOptions(Balanced('[', ']'));
        std::vector<std::int64_t> ids;
        bool closed = false;
        while (true) {
            ids.push_back(Coordinate());
            SkipBlanks();
            if (mPos < mText.size() && mText[mPos] == ';') {
                ++mPos;
                break;
            }
            if (mText.substr(mPos, 2) != "--")
                throw ReadError("TikZ: unsupported path construct '" +
                                tikzr_construct_at(mText, mPos) + "'");
            mPos += 2;
            SkipBlanks();
            if (mText.substr(mPos, 5) == "cycle" &&
                (mPos + 5 >= mText.size() ||
                 !((mText[mPos + 5] >= 'a' && mText[mPos + 5] <= 'z') ||
                   (mText[mPos + 5] >= 'A' && mText[mPos + 5] <= 'Z')))) {
                mPos += 5;
                closed = true;
                SkipBlanks();
                if (mPos >= mText.size() || mText[mPos] != ';')
                    throw ReadError("TikZ: expected ';' after cycle");
                ++mPos;
                break;
            }
        }
        mBuilder.AddCell(ids, closed, fill);
    }

    // `\fill[colour] (x,y) rectangle (x,y);` -- the colour bar's swatches.
    void ColorBarSwatch() {
        SkipBlanks();
        if (mPos < mText.size() && mText[mPos] == '[')
            Balanced('[', ']');
        Coordinate(false);
        SkipBlanks();
        if (mText.substr(mPos, 9) != "rectangle")
            throw ReadError("TikZ: \\fill is only read as a colour-bar rectangle, found '" +
                            tikzr_construct_at(mText, mPos) + "'");
        mPos += 9;
        Coordinate(false);
        SkipBlanks();
        if (mPos >= mText.size() || mText[mPos] != ';')
            throw ReadError("TikZ: expected ';' after the rectangle");
        ++mPos;
    }

    // `\node ... ;` -- the colour bar's labels: up to the next ';' outside braces.
    void SkipStatement() {
        int braces = 0;
        for (; mPos < mText.size(); ++mPos) {
            const char c = mText[mPos];
            if (c == '{') {
                ++braces;
            } else if (c == '}' && braces > 0) {
                --braces;
            } else if (c == ';' && braces == 0) {
                ++mPos;
                return;
            }
        }
        throw ReadError("TikZ: a \\node is never closed");
    }

    std::string_view mText;
    std::size_t mPos = 0;
    detail::DrawingBuilder& mBuilder;
};

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
    const detail::DoubleView point_values(points);
    if (dim == 3) {
        for (std::size_t i = 0; i < num_points; ++i) {
            if (std::fabs(point_values[i * dim + 2]) > 1.0e-14) {
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
        const double px = (0 < dim) ? point_values[idx * dim + 0] : 0.0;
        const double py = (1 < dim) ? point_values[idx * dim + 1] : 0.0;
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
    std::vector<std::string> skipped_types;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad") {
            if (std::find(skipped_types.begin(), skipped_types.end(), type) == skipped_types.end())
                skipped_types.push_back(type);
            continue;
        }

        const NDArray& conn_array = cb.Conn();
        const std::size_t ncols = detail::cols(conn_array);
        const std::size_t n = cb.NumCells();
        const detail::Int64View conn(conn_array);
        for (std::size_t r = 0; r < n; ++r) {
            std::string path;
            for (std::size_t k = 0; k < ncols; ++k) {
                if (k)
                    path += " -- ";
                path += coord(conn[r * ncols + k]);
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

    if (!skipped_types.empty()) {
        std::string joined;
        for (std::size_t i = 0; i < skipped_types.size(); ++i)
            joined += (i ? ", " : "") + skipped_types[i];
        log::warn(
            "TikZ: cell type(s) '{}' are not representable (only line/triangle/quad); "
            "skipping.",
            joined);
    }

    if (colors.mActive && spec.mColorbar && num_points > 0) {
        double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
        for (std::size_t i = 0; i < num_points; ++i) {
            const double cx = (0 < dim) ? point_values[i * dim + 0] : 0.0;
            const double cy = (1 < dim) ? point_values[i * dim + 1] : 0.0;
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
        detail::snprintf_c(buf, sizeof(buf), "scale=%g", *rScale);
        pic_opts = buf;
    }
    if (rLineWidth.has_value()) {
        if (!pic_opts.empty())
            pic_opts += ", ";
        pic_opts += "line width=" + *rLineWidth;
    }
    const std::string pic_opt_str = pic_opts.empty() ? "" : ("[" + pic_opts + "]");

    std::vector<std::string> out = tikz_marker_lines(false, 0.0, 0.0, 0.0);
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

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    for (std::size_t i = 0; i < out.size(); ++i)
        os << out[i] << '\n';
}

Mesh read_tikz(const std::string& rPath, const ReadOptions& rOptions) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);

    // Comments are cut at the first unescaped '%'. The marker is the comment lines
    // that lead the file, before any code.
    std::string body;
    body.reserve(source.Size());
    std::vector<std::string> head;
    bool leading = true;
    for (std::string_view line : detail::split_lines(source.View())) {
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        std::size_t cut = std::string_view::npos;
        for (std::size_t i = 0; i < line.size(); ++i) {
            if (line[i] == '\\') {
                ++i;
            } else if (line[i] == '%') {
                cut = i;
                break;
            }
        }
        const std::string_view code = cut == std::string_view::npos ? line : line.substr(0, cut);
        if (tikzr_trim(code).empty()) {
            if (cut != std::string_view::npos && leading) {
                std::string_view text = line.substr(cut);
                while (!text.empty() && text.front() == '%')
                    text.remove_prefix(1);
                head.emplace_back(tikzr_trim(text));
            }
        } else {
            leading = false;
        }
        body.append(code);
        body.push_back('\n');
    }

    const detail::DrawingMarker marker = detail::parse_drawing_marker(head, "TikZ", false);
    if (!marker.mFound && !rOptions.mLenient)
        throw detail::drawing_unmarked_error("TikZ");

    detail::DrawingBuilder builder("TikZ", "tikz:fill");
    builder.SetWantFill(rOptions.WantsAnyData() && rOptions.WantsArray("tikz:fill"));
    TikzrReader(body, builder).Run();

    // Coordinates are verbatim: no flip, no scale.
    Mesh mesh = builder.Finish(/*Flip=*/false, 1.0, 0.0, 0.0);
    detail::drawing_read_warning("TikZ", marker);
    return mesh;
}

}  // namespace meshioplusplus
