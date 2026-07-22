// Notebook-only display helpers for the meshio++ C++ (xeus-cpp) examples.
//
// meshio++'s core has no rendering dependency, so these notebooks lean on the
// core's own SVG writer (which already does camera projection and data-driven
// colouring) for mesh renders, plus a couple of hand-rolled SVG primitives
// for the small numeric charts the Python notebooks draw with matplotlib.
// Nothing here is part of the library API -- it exists only to give these
// notebooks something to `xcpp::display()`.
#pragma once

#include "xcpp/xdisplay.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace mio_nb {

// -- rich display plumbing --------------------------------------------------

struct SvgImage {
    std::string svg;
};

inline nl::json mime_bundle_repr(const SvgImage& rImg) {
    auto bundle = nl::json::object();
    bundle["image/svg+xml"] = rImg.svg;
    return bundle;
}

inline std::string read_text_file(const std::string& rPath) {
    std::ifstream f(rPath);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

inline std::string next_temp_path(const std::string& rSuffix) {
    static int counter = 0;
    return "/tmp/meshioplusplus_cpp_nb_" + std::to_string(++counter) + rSuffix;
}

// -- mesh renders, via meshio++'s own SVG writer -----------------------------

// Wraps write_svg: writes to a scratch file, reads it back, and returns it as
// a displayable SvgImage. Mirrors the notebook's `render()` helper in the
// Python examples, minus the PyVista dependency -- the camera + colouring are
// meshio++'s own (see doc/formats/svg.md#data-driven-colouring).
inline SvgImage render(const meshioplusplus::Mesh& rMesh, const std::string& rColorBy = "",
                        const std::string& rCmap = "viridis", bool colorbar = false,
                        std::optional<double> vmin = std::nullopt, std::optional<double> vmax = std::nullopt,
                        double azimuth = 45.0, double elevation = 35.264389682754654,
                        double imageWidth = 900.0) {
    std::string path = next_temp_path(".svg");
    meshioplusplus::write_svg(path, rMesh, ".3f", std::nullopt, imageWidth, "#c8c5bd", "#1f4e79", azimuth,
                               elevation, 0.0, rColorBy, std::nullopt, rCmap, vmin, vmax, "#808080", colorbar);
    return SvgImage{read_text_file(path)};
}

// -- synthetic fixtures ------------------------------------------------------

// An n x n x n hexahedron block on the unit-spaced integer lattice
// [0, n]^3 -- the small, easy-to-read fixture the C++ notebooks build several
// operation demos on (convert_cells, refine, smooth, interpolate), mirroring
// the Python notebooks' `np.meshgrid`-based hex block.
inline meshioplusplus::Mesh hex_block(int n, double spacing = 1.0) {
    using namespace meshioplusplus;
    Mesh m;
    auto pid = [n](int i, int j, int k) { return (i * (n + 1) + j) * (n + 1) + k; };
    std::size_t npts = std::size_t(n + 1) * (n + 1) * (n + 1);
    NDArray p(DType::Float64, {npts, 3});
    double* pd = p.As<double>();
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            for (int k = 0; k <= n; ++k) {
                std::size_t id = std::size_t(pid(i, j, k));
                pd[id * 3 + 0] = i * spacing;
                pd[id * 3 + 1] = j * spacing;
                pd[id * 3 + 2] = k * spacing;
            }
    m.AssignPoints(std::move(p));

    std::size_t ncells = std::size_t(n) * n * n;
    NDArray c(DType::Int64, {ncells, 8});
    std::int64_t* cd = c.As<std::int64_t>();
    std::size_t cell = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                std::int64_t corners[8] = {pid(i, j, k),         pid(i + 1, j, k),         pid(i + 1, j + 1, k),
                                            pid(i, j + 1, k),     pid(i, j, k + 1),         pid(i + 1, j, k + 1),
                                            pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)};
                for (int q = 0; q < 8; ++q) cd[cell * 8 + q] = corners[q];
                ++cell;
            }
    m.AddCellBlock("hexahedron", std::move(c));
    return m;
}

// -- small hand-rolled charts, coloured from meshio++'s own colormap tables -

namespace detail_chart {

inline std::string rgb_hex(meshioplusplus::detail::Rgb c) {
    char buf[8];
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x", c.mR, c.mG, c.mB);
    return std::string(buf);
}

inline std::string svg_escape(const std::string& rIn) {
    std::string out;
    for (char c : rIn) {
        if (c == '&')
            out += "&amp;";
        else if (c == '<')
            out += "&lt;";
        else if (c == '>')
            out += "&gt;";
        else
            out += c;
    }
    return out;
}

}  // namespace detail_chart

// Horizontal bar chart -- e.g. file sizes per format, cell counts per type.
inline SvgImage bar_chart(const std::vector<std::string>& rLabels, const std::vector<double>& rValues,
                           const std::string& rTitle, const std::string& rCmap = "viridis") {
    const std::size_t n = rLabels.size();
    const double barH = 30.0, gap = 12.0, left = 220.0, right = 60.0, top = 50.0, width = 760.0;
    const double maxV = *std::max_element(rValues.begin(), rValues.end());
    const double plotW = width - left - right;
    const double height = top + n * (barH + gap) + 20.0;
    const auto* table = meshioplusplus::detail::colormap_table(rCmap);

    std::ostringstream svg;
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << " " << height
        << "\" font-family=\"sans-serif\" font-size=\"13\">";
    svg << "<text x=\"" << (width / 2.0) << "\" y=\"24\" text-anchor=\"middle\" font-size=\"16\">"
        << detail_chart::svg_escape(rTitle) << "</text>";
    for (std::size_t i = 0; i < n; ++i) {
        double y = top + i * (barH + gap);
        double w = maxV > 0.0 ? (rValues[i] / maxV) * plotW : 0.0;
        auto col = meshioplusplus::detail::colormap_lookup(table, n > 1 ? double(i) / double(n - 1) : 0.0);
        svg << "<text x=\"" << (left - 10.0) << "\" y=\"" << (y + barH * 0.65) << "\" text-anchor=\"end\">"
            << detail_chart::svg_escape(rLabels[i]) << "</text>";
        svg << "<rect x=\"" << left << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << barH
            << "\" fill=\"" << detail_chart::rgb_hex(col) << "\" />";
        std::ostringstream val;
        val.precision(4);
        val << rValues[i];
        svg << "<text x=\"" << (left + w + 8.0) << "\" y=\"" << (y + barH * 0.65) << "\">" << val.str()
            << "</text>";
    }
    svg << "</svg>";
    return SvgImage{svg.str()};
}

// Vertical histogram -- edges.size() == counts.size() + 1.
inline SvgImage histogram_chart(const std::vector<double>& rEdges, const std::vector<std::int64_t>& rCounts,
                                 const std::string& rTitle, const std::string& rXLabel,
                                 const std::string& rColor = "#3b82f6") {
    const std::size_t n = rCounts.size();
    const double width = 760.0, height = 380.0, left = 60.0, right = 20.0, top = 50.0, bottom = 60.0;
    const double plotW = width - left - right, plotH = height - top - bottom;
    const std::int64_t maxC = *std::max_element(rCounts.begin(), rCounts.end());
    const double barW = plotW / double(n);

    std::ostringstream svg;
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " << width << " " << height
        << "\" font-family=\"sans-serif\" font-size=\"12\">";
    svg << "<text x=\"" << (width / 2.0) << "\" y=\"24\" text-anchor=\"middle\" font-size=\"16\">"
        << detail_chart::svg_escape(rTitle) << "</text>";
    svg << "<line x1=\"" << left << "\" y1=\"" << (top + plotH) << "\" x2=\"" << (left + plotW) << "\" y2=\""
        << (top + plotH) << "\" stroke=\"black\" />";
    for (std::size_t i = 0; i < n; ++i) {
        double h = maxC > 0 ? (double(rCounts[i]) / double(maxC)) * plotH : 0.0;
        double x = left + i * barW;
        svg << "<rect x=\"" << (x + barW * 0.05) << "\" y=\"" << (top + plotH - h) << "\" width=\""
            << (barW * 0.9) << "\" height=\"" << h << "\" fill=\"" << rColor << "\" />";
    }
    std::ostringstream lo, hi;
    lo.precision(3);
    hi.precision(3);
    lo << rEdges.front();
    hi << rEdges.back();
    svg << "<text x=\"" << left << "\" y=\"" << (top + plotH + 20.0) << "\">" << lo.str() << "</text>";
    svg << "<text x=\"" << (left + plotW) << "\" y=\"" << (top + plotH + 20.0) << "\" text-anchor=\"end\">"
        << hi.str() << "</text>";
    svg << "<text x=\"" << (width / 2.0) << "\" y=\"" << (height - 10.0) << "\" text-anchor=\"middle\">"
        << detail_chart::svg_escape(rXLabel) << "</text>";
    svg << "</svg>";
    return SvgImage{svg.str()};
}

}  // namespace mio_nb
