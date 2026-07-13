#pragma once
//
// C++ analogue of tests/helpers.py: fixture-mesh builders and a round-trip
// helper that writes a mesh to a temp file, reads it back, and asserts the
// points/cells match. Used by the GoogleTest per-format suites.

#include <gtest/gtest.h>

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/mesh.hpp"

namespace mt {

using meshio::CellBlock;
using meshio::DType;
using meshio::Mesh;
using meshio::NDArray;

// ---- builders ----

inline NDArray points_from(const std::vector<std::vector<double>>& pts) {
    std::size_t n = pts.size();
    std::size_t d = n ? pts[0].size() : 3;
    NDArray a(DType::Float64, {n, d});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < d; ++j) a.as<double>()[i * d + j] = pts[i][j];
    return a;
}

inline NDArray conn_from(const std::vector<std::vector<std::int64_t>>& rows) {
    std::size_t n = rows.size();
    std::size_t k = n ? rows[0].size() : 0;
    NDArray a(DType::Int64, {n, k});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < k; ++j) a.as<std::int64_t>()[i * k + j] = rows[i][j];
    return a;
}

inline Mesh make_mesh(std::vector<std::vector<double>> pts, const std::string& type,
                      std::vector<std::vector<std::int64_t>> cells) {
    Mesh m;
    m.points = points_from(pts);
    m.cells.emplace_back(type, conn_from(cells));
    return m;
}

// Fixture meshes mirroring tests/helpers.py (geometry chosen to be valid,
// right-handed volume cells so FLAC3D's determinant reorder round-trips).
inline Mesh line_mesh() {
    return make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "line",
                     {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}});
}
inline Mesh tri_mesh() {
    return make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                     {{0, 1, 2}, {0, 2, 3}});
}
inline Mesh tri_mesh_2d() {
    return make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "triangle",
                     {{0, 1, 2}, {0, 2, 3}});
}
inline Mesh quad_mesh() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {0, 1, 0}}, "quad",
        {{0, 1, 4, 5}, {1, 2, 3, 4}});
}
inline Mesh tet_mesh() {
    return make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 0.5}},
                     "tetra", {{0, 1, 2, 4}, {0, 2, 3, 4}});
}
inline Mesh hex_mesh() {
    return make_mesh({{0, 0, 0},
                      {1, 0, 0},
                      {1, 1, 0},
                      {0, 1, 0},
                      {0, 0, 1},
                      {1, 0, 1},
                      {1, 1, 1},
                      {0, 1, 1}},
                     "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
}
inline Mesh wedge_mesh() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, "wedge",
        {{0, 1, 2, 3, 4, 5}});
}
inline Mesh triangle6_mesh() {
    return make_mesh({{0, 0, 0},
                      {1, 0, 0},
                      {1, 1, 0},
                      {0.5, 0.25, 0},
                      {1.25, 0.5, 0},
                      {0.25, 0.75, 0}},
                     "triangle6", {{0, 1, 2, 3, 4, 5}});
}
inline Mesh quad8_mesh() {
    return make_mesh({{0, 0, 0},
                      {1, 0, 0},
                      {1, 1, 0},
                      {0, 1, 0},
                      {0.5, 0.1, 0},
                      {0.9, 0.5, 0},
                      {0.5, 0.9, 0},
                      {0.1, 0.5, 0}},
                     "quad8", {{0, 1, 2, 3, 4, 5, 6, 7}});
}
inline Mesh tet10_mesh() {
    std::vector<std::vector<double>> p;
    for (int i = 0; i < 10; ++i) p.push_back({0.1 * i, 0.2 * i + 0.05, 0.3 * i + 0.01});
    return make_mesh(p, "tetra10", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}});
}
inline Mesh hex20_mesh() {
    std::vector<std::vector<double>> p;
    for (int i = 0; i < 20; ++i) p.push_back({0.11 * i, 0.07 * i, 0.03 * i});
    std::vector<std::int64_t> row(20);
    for (int i = 0; i < 20; ++i) row[i] = i;
    return make_mesh(p, "hexahedron20", {row});
}
// A hybrid mesh (triangle, quad, triangle blocks).
inline Mesh tri_quad_mesh() {
    Mesh m;
    m.points = points_from(
        {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 1, 0}, {2, 1, 0}, {1, 1, 0}, {0, 1, 0}});
    m.cells.emplace_back("triangle", conn_from({{0, 1, 5}, {0, 5, 6}}));
    m.cells.emplace_back("quad", conn_from({{1, 2, 4, 5}}));
    m.cells.emplace_back("triangle", conn_from({{2, 3, 4}}));
    return m;
}

// ---- comparison / round-trip ----

inline std::string temp_path(const std::string& suffix) {
    static std::atomic<unsigned> counter{0};
    auto dir = std::filesystem::temp_directory_path();
    return (dir / ("meshio_cpp_" + std::to_string(counter++) + suffix)).string();
}

// {type -> multiset of connectivity rows}. Handles formats that merge same-type
// blocks on read; preserves node order within a cell.
inline std::map<std::string, std::multiset<std::vector<std::int64_t>>> cell_rows(
    const Mesh& m) {
    std::map<std::string, std::multiset<std::vector<std::int64_t>>> out;
    for (const auto& cb : m.cells) {
        std::size_t n = cb.num_cells();
        std::size_t k = meshio::detail::cols(cb.data);
        for (std::size_t r = 0; r < n; ++r) {
            std::vector<std::int64_t> row(k);
            for (std::size_t j = 0; j < k; ++j)
                row[j] = meshio::detail::read_int(cb.data, r * k + j);
            out[cb.type].insert(std::move(row));
        }
    }
    return out;
}

inline void expect_points_close(const Mesh& in, const Mesh& out, double atol) {
    ASSERT_EQ(in.num_points(), out.num_points());
    std::size_t din = meshio::detail::cols(in.points);
    std::size_t dout = meshio::detail::cols(out.points);
    ASSERT_GE(dout, din);  // formats may pad 2D -> 3D, never truncate
    for (std::size_t i = 0; i < in.num_points(); ++i)
        for (std::size_t j = 0; j < din; ++j) {
            double a = meshio::detail::read_double(in.points, i * din + j);
            double b = meshio::detail::read_double(out.points, i * dout + j);
            EXPECT_NEAR(a, b, atol) << "point " << i << " comp " << j;
        }
}

inline void expect_mesh_eq(const Mesh& in, const Mesh& out, double atol = 1e-12) {
    expect_points_close(in, out, atol);
    EXPECT_EQ(cell_rows(in), cell_rows(out));
}

using Writer = std::function<void(const std::string&, const Mesh&)>;
using Reader = std::function<Mesh(const std::string&)>;

// Write -> read -> compare, then remove the temp file(s).
inline void roundtrip(const Writer& writer, const Reader& reader, const Mesh& mesh,
                      const std::string& suffix, double atol = 1e-12) {
    std::string path = temp_path(suffix);
    writer(path, mesh);
    Mesh out = reader(path);
    expect_mesh_eq(mesh, out, atol);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

}  // namespace mt
