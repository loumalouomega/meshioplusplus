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
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/typed_view.hpp"
#include "../detail/drawing.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"

// External includes
#include "pugixml.hpp"

namespace meshioplusplus {

namespace {

// Format a single double with a printf-style spec (spec without leading '%',
// e.g. ".3f"), mirroring the Python reference's `format(x, float_fmt)`.
std::string svg_fmt_num(double value, const std::string& rSpec) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), ("%" + rSpec).c_str(), value);
    return buf;
}

// SVG's colour vocabulary. Lowercase hex, as Python's "%02x" also produces.
std::string svg_color_hex(detail::Rgb c) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x", static_cast<int>(c.mR), static_cast<int>(c.mG),
                  static_cast<int>(c.mB));
    return buf;
}

// The marker comment right after the opening tag: the provenance tag, the
// `Drawing: ...` line, then whatever else the provenance block carries. Always
// the multi-line form, with "--" broken apart as the XML comment rules need
// (the twin of `render_xml_comment` in `_provenance.py`).
std::string svg_marker_comment(bool Projected, double Azimuth, double Elevation, double Roll,
                               double Scale, double OffsetX, double OffsetY) {
    const double scale_offset[3] = {Scale, OffsetX, OffsetY};
    std::string out = "<!--\n";
    for (std::string line : detail::drawing_comment_lines(
             detail::drawing_marker_line(Projected, Azimuth, Elevation, Roll, scale_offset))) {
        for (std::size_t at = line.find("--"); at != std::string::npos; at = line.find("--", at))
            line.insert(at + 1, 1, ' ');
        out += line;
        out += '\n';
    }
    return out + "-->";
}

// The extra viewBox width a colorbar occupies, as a multiple of the drawing's
// width. Only the viewBox grows: the scaling factor, the stroke width and every
// mesh coordinate are computed from the unmodified width.
constexpr double kSvgColorbarGap = 0.08;
constexpr double kSvgColorbarWidth = 0.06;
constexpr double kSvgColorbarLabels = 0.30;
constexpr std::size_t kSvgColorbarSegments = 32;

// The colorbar: a vertical gradient to the right of the drawing plus min/max
// labels. <rect>, not <path>, because the document-level `path {...}` rule
// would otherwise stroke every swatch. Emitted after the mesh so the painter
// ordering (and a d-attribute-only reading of the output) is undisturbed.
void svg_append_colorbar(std::ostream& rOs, const detail::FaceColors& rColors, double MinX,
                         double MinY, double Width, double Height, const std::string& rFloatFmt) {
    const double x0 = MinX + Width + Width * kSvgColorbarGap;
    const double bar_w = Width * kSvgColorbarWidth;
    for (std::size_t i = 0; i < kSvgColorbarSegments; ++i) {
        const double y0 =
            MinY + Height * (static_cast<double>(i) / static_cast<double>(kSvgColorbarSegments));
        const double seg_h = Height / static_cast<double>(kSvgColorbarSegments);
        // y grows downward in SVG, so the FIRST segment is the top one and
        // must carry the HIGH end of the range.
        const detail::Rgb c = detail::colormap_lookup(
            rColors.mpTable, static_cast<double>(kSvgColorbarSegments - 1 - i) /
                                 static_cast<double>(kSvgColorbarSegments - 1));
        rOs << "<rect x=\"" << svg_fmt_num(x0, rFloatFmt) << "\" y=\"" << svg_fmt_num(y0, rFloatFmt)
            << "\" width=\"" << svg_fmt_num(bar_w, rFloatFmt) << "\" height=\""
            << svg_fmt_num(seg_h, rFloatFmt) << "\" fill=\"" << svg_color_hex(c) << "\" />";
    }
    const double label_x = x0 + bar_w + Width * 0.02;
    const double font_size = Width * 0.04;
    rOs << "<text x=\"" << svg_fmt_num(label_x, rFloatFmt) << "\" y=\""
        << svg_fmt_num(MinY + font_size, rFloatFmt) << "\" font-size=\""
        << svg_fmt_num(font_size, rFloatFmt) << "\">" << svg_fmt_num(rColors.mVMax, rFloatFmt)
        << "</text>";
    rOs << "<text x=\"" << svg_fmt_num(label_x, rFloatFmt) << "\" y=\""
        << svg_fmt_num(MinY + Height, rFloatFmt) << "\" font-size=\""
        << svg_fmt_num(font_size, rFloatFmt) << "\">" << svg_fmt_num(rColors.mVMin, rFloatFmt)
        << "</text>";
}

// 3D rendering path: project the (already skin-extracted) surface mesh with
// the orthographic camera and emit its faces back-to-front. Replicates the
// flat path's bbox -> y-flip -> image_width scaling -> <path> emission, but
// over the sorted projected faces.
void svg_proj_write(const std::string& rPath, const Mesh& rSourceMesh, const Mesh& rDrawMesh,
                    const std::string& rFloatFmt, const std::optional<std::string>& rStrokeWidth,
                    const std::optional<double>& rImageWidth, const std::string& rFill,
                    const std::string& rStroke, double azimuth, double elevation, double roll,
                    const detail::ColorSpec& rColorSpec, const std::string& rNanColor) {
    detail::ProjectedSurface ps = detail::project_surface(rDrawMesh, azimuth, elevation, roll);
    const detail::FaceColors colors = detail::resolve_face_colors(
        rColorSpec, rSourceMesh, rDrawMesh, detail::color_faces_from_projection(ps.mFaces));
    std::vector<double>& x = ps.mX;
    std::vector<double>& y = ps.mY;
    const std::size_t num_points = x.size();

    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    if (num_points > 0) {
        min_x = max_x = x[0];
        min_y = max_y = y[0];
        for (std::size_t i = 1; i < num_points; ++i) {
            min_x = std::min(min_x, x[i]);
            max_x = std::max(max_x, x[i]);
            min_y = std::min(min_y, y[i]);
            max_y = std::max(max_y, y[i]);
        }
    }

    // Flip y (projected math convention y-up -> SVG screen convention y-down).
    const double y_flip = max_y + min_y;
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    double applied_scale = 1.0;
    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
        applied_scale = scaling_factor;
        min_x *= scaling_factor;
        min_y *= scaling_factor;
        width *= scaling_factor;
        height *= scaling_factor;
        for (std::size_t i = 0; i < num_points; ++i) {
            x[i] *= scaling_factor;
            y[i] *= scaling_factor;
        }
    }

    std::string stroke_width;
    if (rStrokeWidth.has_value()) {
        stroke_width = *rStrokeWidth;
    } else {
        char buf[64];
        detail::snprintf_c(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const bool colorbar = colors.mActive && rColorSpec.mColorbar;
    const double view_width =
        colorbar ? width * (1.0 + kSvgColorbarGap + kSvgColorbarWidth + kSvgColorbarLabels) : width;

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\""
       << svg_fmt_num(min_x, rFloatFmt) << ' ' << svg_fmt_num(min_y, rFloatFmt) << ' '
       << svg_fmt_num(view_width, rFloatFmt) << ' ' << svg_fmt_num(height, rFloatFmt) << "\">";

    os << svg_marker_comment(true, azimuth, elevation, roll, applied_scale, 0.0, y_flip);

    os << "<style>path {fill: " << rFill << "; stroke: " << rStroke
       << "; stroke-width: " << stroke_width << "; stroke-linejoin:bevel}</style>";

    for (std::size_t f = 0; f < ps.mFaces.size(); ++f) {
        const detail::ProjectedFace& face = ps.mFaces[f];
        std::string d;
        for (std::uint8_t k = 0; k < face.mNumNodes; ++k) {
            const std::int64_t p = face.mNodes[k];
            d += (k == 0) ? "M " : "L ";
            d += svg_fmt_num(x[static_cast<std::size_t>(p)], rFloatFmt);
            d += ' ';
            d += svg_fmt_num(y[static_cast<std::size_t>(p)], rFloatFmt);
        }
        if (!face.mIsLine)
            d += "Z";
        // A per-path fill attribute ALONE does not override the document-
        // level rule: a presentation attribute has zero CSS specificity and
        // loses to any stylesheet rule, even the plain "path {fill: ...}"
        // type selector above - every cascade-honouring renderer (browsers,
        // cairosvg) paints every face in the flat fallback colour despite
        // this attribute. The redundant inline style wins the cascade and
        // actually colours the face; fill stays for inspection/tests.
        if (colors.mActive && !face.mIsLine) {
            const std::optional<detail::Rgb> c = colors.Color(f);
            const std::string hex_color = c.has_value() ? svg_color_hex(*c) : rNanColor;
            os << "<path d=\"" << d << "\" fill=\"" << hex_color << "\" style=\"fill:" << hex_color
               << "\" />";
        } else {
            os << "<path d=\"" << d << "\" />";
        }
    }

    if (colorbar)
        svg_append_colorbar(os, colors, min_x, min_y, width, height, rFloatFmt);

    os << "</svg>";
}

// ---------------------------------------------------------------------------
// Reader. Inverts exactly what the writer above emits, and nothing more: see
// doc/formats/svg.md ("Reading") and detail/drawing.hpp.
// ---------------------------------------------------------------------------

constexpr const char* kSvgrNamespace = "http://www.w3.org/2000/svg";

// Every element the reader looks at: no transform, and no namespace of its own
// other than SVG's.
void svgr_check_element(const pugi::xml_node& rNode) {
    for (const pugi::xml_attribute attr : rNode.attributes()) {
        const std::string_view name = attr.name();
        if (name == "transform")
            throw ReadError("SVG: a transform on <" + detail::drawing_safe(rNode.name()) +
                            "> is not supported (meshio++ draws absolute coordinates)");
        if (name == "xmlns" && std::string_view(attr.value()) != kSvgrNamespace)
            throw ReadError("SVG: <" + detail::drawing_safe(rNode.name()) +
                            "> is in a namespace other than SVG's");
    }
}

// The four numbers of a viewBox.
std::array<double, 4> svgr_view_box(std::string_view Text) {
    std::array<double, 4> box{};
    std::size_t pos = 0;
    auto skip = [&] {
        while (pos < Text.size() && (detail::text_is_blank(Text[pos]) || Text[pos] == ','))
            ++pos;
    };
    for (double& value : box) {
        skip();
        const std::size_t len = detail::drawing_scan_number(Text, pos);
        if (len == 0)
            throw ReadError("SVG: malformed viewBox");
        value = detail::drawing_number_value(Text.substr(pos, len), "SVG");
        pos += len;
    }
    skip();
    if (pos != Text.size())
        throw ReadError("SVG: malformed viewBox");
    return box;
}

// One path's `d`: `M x y (L x y)* Z?`, absolute, one subpath.
void svgr_read_path(std::string_view D, detail::DrawingBuilder& rBuilder,
                    const std::optional<std::array<int, 3>>& rFill) {
    std::size_t pos = 0;
    auto skip = [&] {
        while (pos < D.size() && (detail::text_is_blank(D[pos]) || D[pos] == ','))
            ++pos;
    };
    auto pair = [&]() -> std::int64_t {
        skip();
        const std::size_t xlen = detail::drawing_scan_number(D, pos);
        if (xlen == 0)
            throw ReadError("SVG: a path command is missing its coordinates");
        const std::string_view x = D.substr(pos, xlen);
        pos += xlen;
        skip();
        const std::size_t ylen = detail::drawing_scan_number(D, pos);
        if (ylen == 0)
            throw ReadError("SVG: a path command is missing its coordinates");
        const std::string_view y = D.substr(pos, ylen);
        pos += ylen;
        return rBuilder.Weld(x, y, [&] {
            return std::make_pair(detail::drawing_number_value(x, "SVG"),
                                  detail::drawing_number_value(y, "SVG"));
        });
    };

    std::vector<std::int64_t> ids;
    bool started = false;
    bool closed = false;
    while (true) {
        skip();
        if (pos >= D.size())
            break;
        const char c = D[pos];
        if (closed)
            throw ReadError("SVG: a path continues after Z (several subpaths are not supported)");
        if (c == 'M') {
            if (started)
                throw ReadError("SVG: a path with several subpaths is not supported");
            ++pos;
            ids.push_back(pair());
            started = true;
        } else if (c == 'L') {
            if (!started)
                throw ReadError("SVG: a path must start with M");
            ++pos;
            ids.push_back(pair());
        } else if (c == 'Z') {
            if (!started)
                throw ReadError("SVG: a path must start with M");
            ++pos;
            closed = true;
        } else if (std::string_view("mlhvcsqtaz").find(c) != std::string_view::npos) {
            throw ReadError(std::string("SVG: the relative path command '") + c +
                            "' is not supported (meshio++ writes absolute M, L and Z)");
        } else if (std::string_view("CSQTA").find(c) != std::string_view::npos) {
            throw ReadError(std::string("SVG: the curve command '") + c + "' is not supported");
        } else if (c == 'H' || c == 'V') {
            throw ReadError(std::string("SVG: the path command '") + c + "' is not supported");
        } else if (detail::drawing_scan_number(D, pos) != 0) {
            throw ReadError("SVG: repeated coordinates after one path command are not supported");
        } else {
            throw ReadError("SVG: unexpected character '" +
                            detail::drawing_safe(std::string_view(&c, 1)) + "' in a path");
        }
    }
    if (!started)
        throw ReadError("SVG: an empty path");
    rBuilder.AddCell(ids, closed, rFill);
}

std::vector<std::string> svgr_split_lines(std::string_view Text) {
    std::vector<std::string> lines;
    for (std::string_view line : detail::split_lines(Text)) {
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        lines.emplace_back(line);
    }
    return lines;
}

}  // namespace

void write_svg(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
               const std::optional<std::string>& rStrokeWidth,
               const std::optional<double>& rImageWidth, const std::string& rFill,
               const std::string& rStroke, double azimuth, double elevation, double roll,
               const std::string& rColorBy, const std::optional<int>& rComponent,
               const std::string& rCmap, const std::optional<double>& rVMin,
               const std::optional<double>& rVMax, const std::string& rNanColor, bool Colorbar) {
    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // A genuinely non-flat 3D mesh takes the projected-rendering path (skin
    // extraction for volume cells + orthographic camera); a flat one (every
    // z ~ 0) keeps the classic 2D path below, byte-identical to before.
    detail::ColorSpec spec;
    spec.mColorBy = rColorBy;
    spec.mComponent = rComponent;
    spec.mCmap = rCmap;
    spec.mVMin = rVMin;
    spec.mVMax = rVMax;
    spec.mColorbar = Colorbar;

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
                    svg_proj_write(rPath, rMesh, skin, rFloatFmt, rStrokeWidth, rImageWidth, rFill,
                                   rStroke, azimuth, elevation, roll, spec, rNanColor);
                } else {
                    svg_proj_write(rPath, rMesh, rMesh, rFloatFmt, rStrokeWidth, rImageWidth, rFill,
                                   rStroke, azimuth, elevation, roll, spec, rNanColor);
                }
                return;
            }
        }
    }

    // Copy the first two coordinate columns.
    std::vector<double> x(num_points), y(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        x[i] = (0 < dim) ? point_values[i * dim + 0] : 0.0;
        y[i] = (1 < dim) ? point_values[i * dim + 1] : 0.0;
    }

    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    if (num_points > 0) {
        min_x = max_x = x[0];
        min_y = max_y = y[0];
        for (std::size_t i = 1; i < num_points; ++i) {
            min_x = std::min(min_x, x[i]);
            max_x = std::max(max_x, x[i]);
            min_y = std::min(min_y, y[i]);
            max_y = std::max(max_y, y[i]);
        }
    }

    // Flip y (mesh math convention y-up -> SVG screen convention y-down).
    const double y_flip = max_y + min_y;
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    double applied_scale = 1.0;
    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
        applied_scale = scaling_factor;
        min_x *= scaling_factor;
        min_y *= scaling_factor;
        width *= scaling_factor;
        height *= scaling_factor;
        for (std::size_t i = 0; i < num_points; ++i) {
            x[i] *= scaling_factor;
            y[i] *= scaling_factor;
        }
    }

    std::string stroke_width;
    if (rStrokeWidth.has_value()) {
        stroke_width = *rStrokeWidth;
    } else {
        char buf[64];
        detail::snprintf_c(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    // On the flat path the drawn mesh IS the source mesh, so a face's cell
    // index needs no provenance indirection.
    const detail::FaceColors colors =
        detail::resolve_face_colors(spec, rMesh, rMesh, detail::color_faces_flat(rMesh));
    const bool colorbar = colors.mActive && spec.mColorbar;
    const double view_width =
        colorbar ? width * (1.0 + kSvgColorbarGap + kSvgColorbarWidth + kSvgColorbarLabels) : width;

    // viewBox: "min_x min_y width height", each float_fmt-formatted.
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\""
       << svg_fmt_num(min_x, rFloatFmt) << ' ' << svg_fmt_num(min_y, rFloatFmt) << ' '
       << svg_fmt_num(view_width, rFloatFmt) << ' ' << svg_fmt_num(height, rFloatFmt) << "\">";

    os << svg_marker_comment(false, 0.0, 0.0, 0.0, applied_scale, 0.0, y_flip);

    // Use path (not polygon): svgo rewrites polygons to paths but drops style.
    os << "<style>path {fill: " << rFill << "; stroke: " << rStroke
       << "; stroke-width: " << stroke_width << "; stroke-linejoin:bevel}</style>";

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
            std::string d;
            for (std::size_t k = 0; k < ncols; ++k) {
                const std::int64_t p = conn[r * ncols + k];
                // "M x y" for the first vertex, "L x y" for the rest — no
                // separating space before the command letter (matches the
                // Python reference's concatenated format strings).
                d += (k == 0) ? "M " : "L ";
                d += svg_fmt_num(x[static_cast<std::size_t>(p)], rFloatFmt);
                d += ' ';
                d += svg_fmt_num(y[static_cast<std::size_t>(p)], rFloatFmt);
            }
            // triangle/quad are closed; line stays open.
            if (type != "line")
                d += "Z";
            const std::size_t f = face_index++;
            if (colors.mActive && type != "line") {
                const std::optional<detail::Rgb> c = colors.Color(f);
                const std::string hex_color = c.has_value() ? svg_color_hex(*c) : rNanColor;
                // See the twin comment in svg_proj_write(): the inline style
                // is what actually wins the CSS cascade; fill stays for
                // inspection/tests.
                os << "<path d=\"" << d << "\" fill=\"" << hex_color
                   << "\" style=\"fill:" << hex_color << "\" />";
            } else {
                os << "<path d=\"" << d << "\" />";
            }
        }
    }

    if (colorbar)
        svg_append_colorbar(os, colors, min_x, min_y, width, height, rFloatFmt);

    os << "</svg>";

    if (!skipped_types.empty()) {
        std::string joined;
        for (std::size_t i = 0; i < skipped_types.size(); ++i)
            joined += (i ? ", " : "") + skipped_types[i];
        log::warn(
            "SVG: cell type(s) '{}' are not representable (only line/triangle/quad); "
            "skipping.",
            joined);
    }
}

Mesh read_svg(const std::string& rPath, const ReadOptions& rOptions) {
    pugi::xml_document doc;
    const pugi::xml_parse_result parsed =
        doc.load_file(rPath.c_str(), pugi::parse_default | pugi::parse_comments);
    if (!parsed)
        throw ReadError("SVG: could not parse " + rPath + " (" + parsed.description() + ")");
    const pugi::xml_node root = doc.document_element();
    if (!root || std::string_view(root.name()) != "svg")
        throw ReadError("SVG: the root element is not <svg>");
    svgr_check_element(root);
    for (const pugi::xml_attribute attr : root.attributes()) {
        if (std::string_view(attr.name()).rfind("xmlns:", 0) == 0)
            throw ReadError("SVG: namespace prefixes are not supported");
    }

    // The marker: the comments that lead the root's children.
    std::vector<std::string> head;
    for (const pugi::xml_node node : root.children()) {
        if (node.type() == pugi::node_element)
            break;
        if (node.type() == pugi::node_comment) {
            std::vector<std::string> lines = svgr_split_lines(node.value());
            head.insert(head.end(), lines.begin(), lines.end());
        }
    }
    const detail::DrawingMarker marker = detail::parse_drawing_marker(head, "SVG", true);
    if (!marker.mFound && !rOptions.mLenient)
        throw detail::drawing_unmarked_error("SVG");

    const pugi::xml_attribute view_box_attr = root.attribute("viewBox");
    if (!view_box_attr)
        throw ReadError("SVG: the <svg> has no viewBox");
    const std::array<double, 4> view_box = svgr_view_box(view_box_attr.value());

    detail::DrawingBuilder builder("SVG", "svg:fill");
    builder.SetWantFill(rOptions.WantsAnyData() && rOptions.WantsArray("svg:fill"));
    for (const pugi::xml_node node : root.children()) {
        if (node.type() != pugi::node_element)
            continue;
        const std::string_view name = node.name();
        svgr_check_element(node);
        if (name == "path") {
            const pugi::xml_attribute d = node.attribute("d");
            if (!d)
                throw ReadError("SVG: a <path> has no d attribute");
            svgr_read_path(d.value(), builder,
                           detail::drawing_hex_color(node.attribute("fill").value()));
        } else if (name == "style" || name == "rect" || name == "text") {
            // The stylesheet is the writer's own; <rect> and <text> are the colour bar.
        } else if (name == "polygon") {
            throw ReadError(
                "SVG: <polygon> is not supported (meshio++ draws faces as <path>; tools such as "
                "svgo rewrite polygons to paths)");
        } else {
            throw ReadError("SVG: the element <" + detail::drawing_safe(name) +
                            "> is not supported");
        }
    }

    // model = (x / S + OX, OY - y / S). Without a marker S = 1 and the y-flip is undone
    // from the viewBox: it was about the middle of the drawing's y extent.
    double scale = 1.0, offset_x = 0.0, offset_y = 2.0 * view_box[1] + view_box[3];
    if (marker.mFound) {
        scale = marker.mScale;
        offset_x = marker.mOffsetX;
        offset_y = marker.mOffsetY;
    }
    Mesh mesh = builder.Finish(/*Flip=*/true, scale, offset_x, offset_y);
    detail::drawing_read_warning("SVG", marker);
    return mesh;
}

}  // namespace meshioplusplus
