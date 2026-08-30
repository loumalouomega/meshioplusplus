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
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/skin.hpp"

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
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
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
        std::snprintf(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const bool colorbar = colors.mActive && rColorSpec.mColorbar;
    const double view_width =
        colorbar ? width * (1.0 + kSvgColorbarGap + kSvgColorbarWidth + kSvgColorbarLabels) : width;

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\""
       << svg_fmt_num(min_x, rFloatFmt) << ' ' << svg_fmt_num(min_y, rFloatFmt) << ' '
       << svg_fmt_num(view_width, rFloatFmt) << ' ' << svg_fmt_num(height, rFloatFmt) << "\">";

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
            os << "<path d=\"" << d << "\" fill=\"" << hex_color << "\" style=\"fill:"
               << hex_color << "\" />";
        } else {
            os << "<path d=\"" << d << "\" />";
        }
    }

    if (colorbar)
        svg_append_colorbar(os, colors, min_x, min_y, width, height, rFloatFmt);

    os << "</svg>";
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
        x[i] = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
        y[i] = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
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
    for (std::size_t i = 0; i < num_points; ++i)
        y[i] = max_y + min_y - y[i];

    double width = max_x - min_x;
    double height = max_y - min_y;

    if (rImageWidth.has_value() && width != 0.0) {
        const double scaling_factor = *rImageWidth / width;
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
        std::snprintf(buf, sizeof(buf), "%g", width / 100.0);
        stroke_width = buf;
    }

    std::ofstream os(rPath, std::ios::binary);
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

    // Use path (not polygon): svgo rewrites polygons to paths but drops style.
    os << "<style>path {fill: " << rFill << "; stroke: " << rStroke
       << "; stroke-width: " << stroke_width << "; stroke-linejoin:bevel}</style>";

    std::size_t face_index = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad")
            continue;

        const NDArray& conn = cb.Conn();
        const std::size_t ncols = detail::cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            std::string d;
            for (std::size_t k = 0; k < ncols; ++k) {
                const std::int64_t p = detail::read_int(conn, r * ncols + k);
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
                os << "<path d=\"" << d << "\" fill=\"" << hex_color << "\" style=\"fill:"
                   << hex_color << "\" />";
            } else {
                os << "<path d=\"" << d << "\" />";
            }
        }
    }

    if (colorbar)
        svg_append_colorbar(os, colors, min_x, min_y, width, height, rFloatFmt);

    os << "</svg>";
}

}  // namespace meshioplusplus
