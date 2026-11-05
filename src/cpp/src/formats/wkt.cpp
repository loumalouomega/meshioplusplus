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
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/exact_weld.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

std::vector<double> parse_point(const std::string& rS) {
    std::vector<double> p;
    detail::TextStream iss(rS);
    std::string tok;
    while (iss >> tok)
        p.push_back(detail::parse_double(tok));
    return p;
}

}  // namespace

Mesh read_wkt(const std::string& rPath) {
    const detail::FileSource s_source = detail::open_source(rPath, "Could not open file: " + rPath);
    const std::string_view s = s_source.View();

    // Must be a TIN.
    std::size_t tin = s.find("TIN");
    if (tin == std::string::npos)
        throw ReadError("Invalid WKT TIN");

    // Every vertex tuple as written (at most four coordinates, XYZM, padded),
    // its dimensionality, and how many each triangle's linestring holds; the
    // weld that numbers them runs once all are read.
    std::vector<double> coords;
    std::vector<std::size_t> dims;
    std::vector<std::size_t> ring_sizes;
    std::vector<std::array<std::int64_t, 3>> tris;

    // Each triangle's linestring lives at parenthesis depth 3
    // (TIN -> depth 1, triangle -> depth 2, linestring -> depth 3); whitespace
    // may appear between the parens, so track depth rather than match "((".
    std::vector<std::string> triangle_strs;
    {
        int depth = 0;
        std::string cur;
        bool capturing = false;
        for (std::size_t i = tin; i < s.size(); ++i) {
            char c = s[i];
            if (c == '(') {
                ++depth;
                if (depth == 3) {
                    capturing = true;
                    cur.clear();
                }
            } else if (c == ')') {
                if (depth == 3) {
                    capturing = false;
                    triangle_strs.push_back(cur);
                }
                --depth;
            } else if (capturing) {
                cur += c;
            }
        }
    }

    for (const std::string& inner : triangle_strs) {
        // Split the triangle's vertices on commas; each is a coordinate tuple.
        std::size_t count = 0;
        std::size_t start = 0;
        while (start <= inner.size()) {
            std::size_t comma = inner.find(',', start);
            std::string part =
                inner.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
            std::vector<double> pt = parse_point(part);
            if (!pt.empty()) {
                if (pt.size() > 4)
                    throw ReadError("WKT point has more than four coordinates");
                dims.push_back(pt.size());
                pt.resize(4, 0.0);
                coords.insert(coords.end(), pt.begin(), pt.end());
                ++count;
            }
            if (comma == std::string::npos)
                break;
            start = comma + 1;
        }
        ring_sizes.push_back(count);
    }

    // Exact-value weld (`==`: -0.0 meets +0.0, a NaN meets nothing), first
    // appearance first -- the numbering the per-vertex hash map gave.
    const std::size_t nverts = dims.size();
    std::vector<detail::ExactWeldKey> keys(nverts);
    for (std::size_t i = 0; i < nverts; ++i)
        keys[i] = detail::exact_weld_key(coords.data() + i * 4, dims[i],
                                         detail::ExactWeldMode::Values, i);
    const detail::ExactWeld weld = detail::exact_weld_keys(keys);

    std::size_t at = 0;
    for (const std::size_t count : ring_sizes) {
        const std::int64_t* idxs = weld.mIdOfRow.data() + at;
        at += count;
        if (count != 4 || idxs[0] != idxs[3])
            throw ReadError("WKT triangle is not a closed linestring");
        tris.push_back({idxs[0], idxs[1], idxs[2]});
    }

    std::vector<std::vector<double>> points(weld.NumIds());
    for (std::size_t id = 0; id < points.size(); ++id) {
        const std::size_t row = weld.mFirstRow[id];
        points[id].assign(coords.begin() + static_cast<std::ptrdiff_t>(row * 4),
                          coords.begin() + static_cast<std::ptrdiff_t>(row * 4 + dims[row]));
    }

    // Points: all must share a dimensionality.
    std::size_t dim = points.empty() ? 3 : points.front().size();
    for (const auto& p : points)
        if (p.size() != dim)
            throw ReadError("WKT points have mixed dimensionality");

    Mesh mesh;
    NDArray pts(DType::Float64, {points.size(), dim});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t j = 0; j < dim; ++j)
            pp[i * dim + j] = points[i][j];
    mesh.AssignPoints(std::move(pts));

    NDArray data(DType::Int64, {tris.size(), 3});
    std::int64_t* dp = data.As<std::int64_t>();
    for (std::size_t i = 0; i < tris.size(); ++i)
        for (int j = 0; j < 3; ++j)
            dp[i * 3 + j] = tris[i][j];
    mesh.AddCellBlock("triangle", std::move(data));

    return mesh;
}

void write_wkt(const std::string& rPath, const Mesh& rMesh) {
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const detail::DoubleView point_values(points);

    auto point_str = [&](std::int64_t p) {
        std::string out;
        char buf[32];
        for (std::size_t j = 0; j < dim; ++j) {
            detail::snprintf_c(buf, sizeof(buf), "%.17g",
                               point_values[static_cast<std::size_t>(p) * dim + j]);
            if (j)
                out += " ";
            out += buf;
        }
        return out;
    };

    f << "TIN (";
    std::string joiner;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle")
            continue;
        const NDArray& conn_array = cb.Conn();
        const std::size_t ncols = detail::cols(conn_array);
        const std::size_t n = cb.NumCells();
        const detail::Int64View conn(conn_array);
        for (std::size_t r = 0; r < n; ++r) {
            std::int64_t a = conn[r * ncols + 0];
            std::int64_t b = conn[r * ncols + 1];
            std::int64_t c = conn[r * ncols + 2];
            std::string sa = point_str(a);
            f << joiner << "((" << sa << ", " << point_str(b) << ", " << point_str(c) << ", " << sa
              << "))";
            joiner = ", ";
        }
    }
    f << ")";
}

}  // namespace meshioplusplus
