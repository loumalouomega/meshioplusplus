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
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

std::vector<double> parse_point(const std::string& s) {
    std::vector<double> p;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) p.push_back(std::strtod(tok.c_str(), nullptr));
    return p;
}

}  // namespace

Mesh read_wkt(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    // Must be a TIN.
    std::size_t tin = s.find("TIN");
    if (tin == std::string::npos) throw ReadError("Invalid WKT TIN");

    std::map<std::vector<double>, std::int64_t> point_index;  // exact-value dedup
    std::vector<std::vector<double>> points;                  // insertion order
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
                if (depth == 3) { capturing = true; cur.clear(); }
            } else if (c == ')') {
                if (depth == 3) { capturing = false; triangle_strs.push_back(cur); }
                --depth;
            } else if (capturing) {
                cur += c;
            }
        }
    }

    for (const std::string& inner : triangle_strs) {
        // Split the triangle's vertices on commas; each is a coordinate tuple.
        std::vector<std::int64_t> idxs;
        std::size_t start = 0;
        while (start <= inner.size()) {
            std::size_t comma = inner.find(',', start);
            std::string part = inner.substr(
                start, comma == std::string::npos ? std::string::npos : comma - start);
            std::vector<double> pt = parse_point(part);
            if (!pt.empty()) {
                auto it = point_index.find(pt);
                std::int64_t id;
                if (it == point_index.end()) {
                    id = static_cast<std::int64_t>(points.size());
                    point_index.emplace(pt, id);
                    points.push_back(pt);
                } else {
                    id = it->second;
                }
                idxs.push_back(id);
            }
            if (comma == std::string::npos) break;
            start = comma + 1;
        }

        if (idxs.size() != 4 || idxs.front() != idxs.back())
            throw ReadError("WKT triangle is not a closed linestring");
        tris.push_back({idxs[0], idxs[1], idxs[2]});
    }

    // Points: all must share a dimensionality.
    std::size_t dim = points.empty() ? 3 : points.front().size();
    for (const auto& p : points)
        if (p.size() != dim) throw ReadError("WKT points have mixed dimensionality");

    Mesh mesh;
    mesh.points = NDArray(DType::Float64, {points.size(), dim});
    double* pp = mesh.points.as<double>();
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t j = 0; j < dim; ++j) pp[i * dim + j] = points[i][j];

    NDArray data(DType::Int64, {tris.size(), 3});
    std::int64_t* dp = data.as<std::int64_t>();
    for (std::size_t i = 0; i < tris.size(); ++i)
        for (int j = 0; j < 3; ++j) dp[i * 3 + j] = tris[i][j];
    mesh.cells.emplace_back("triangle", std::move(data));

    return mesh;
}

void write_wkt(const std::string& path, const Mesh& mesh) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw WriteError("Could not open file for writing: " + path);

    const std::size_t dim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    auto point_str = [&](std::int64_t p) {
        std::string out;
        char buf[32];
        for (std::size_t j = 0; j < dim; ++j) {
            std::snprintf(buf, sizeof(buf), "%.17g",
                          detail::read_double(mesh.points, static_cast<std::size_t>(p) * dim + j));
            if (j) out += " ";
            out += buf;
        }
        return out;
    };

    f << "TIN (";
    std::string joiner;
    for (const auto& cb : mesh.cells) {
        if (cb.type != "triangle") continue;
        const std::size_t ncols = detail::cols(cb.data);
        const std::size_t n = cb.num_cells();
        for (std::size_t r = 0; r < n; ++r) {
            std::int64_t a = detail::read_int(cb.data, r * ncols + 0);
            std::int64_t b = detail::read_int(cb.data, r * ncols + 1);
            std::int64_t c = detail::read_int(cb.data, r * ncols + 2);
            std::string sa = point_str(a);
            f << joiner << "((" << sa << ", " << point_str(b) << ", " << point_str(c)
              << ", " << sa << "))";
            joiner = ", ";
        }
    }
    f << ")";
}

}  // namespace meshioplusplus
