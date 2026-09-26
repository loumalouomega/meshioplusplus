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
#pragma once

/**
 * @file golden_fixtures.hpp
 * @brief Meshes the golden-digest tests share (`test_facet_tables.cpp`,
 *        `test_op_goldens.cpp`): built from integer arithmetic only, so the
 *        same bits on every platform, and a digest helper.
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

#include "../../src/cpp/benchmark/mesh_digest.hpp"
#include "mesh_fixtures.hpp"

namespace golden {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::bench::MeshDigest;

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#define GOLDEN_PINNED 1
#endif

// Integer LCG offsets scaled by a power of two: bit-identical on every
// platform (no libm).
struct FtJitter {
    std::uint32_t mState = 2463534242u;
    double Next(double amplitude) {
        mState = mState * 1664525u + 1013904223u;
        return amplitude * (static_cast<double>(mState >> 16) / 32768.0 - 1.0);
    }
};

inline NDArray ft_grid_points(std::size_t n, double jitter) {
    const std::size_t np = n + 1;
    NDArray pts = NDArray::Uninit(DType::Float64, {np * np * np, 3});
    double* p = pts.As<double>();
    FtJitter jit;
    for (std::size_t k = 0; k < np; ++k)
        for (std::size_t j = 0; j < np; ++j)
            for (std::size_t i = 0; i < np; ++i) {
                const std::size_t ijk[3] = {i, j, k};
                for (int d = 0; d < 3; ++d) {
                    const bool face = ijk[d] == 0 || ijk[d] == n;
                    const double off = jit.Next(jitter);
                    p[((k * np + j) * np + i) * 3 + d] =
                        (static_cast<double>(ijk[d]) + (face ? 0.0 : off)) / static_cast<double>(n);
                }
            }
    return pts;
}

inline std::int64_t ft_node(std::size_t n, std::size_t i, std::size_t j, std::size_t k) {
    return static_cast<std::int64_t>((k * (n + 1) + j) * (n + 1) + i);
}

// Kuhn tetrahedra, every one positively oriented.
inline Mesh ft_tet_cube(std::size_t n, double jitter) {
    static const int tets[6][4][3] = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 1}, {1, 0, 0}, {1, 1, 1}},
        {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 1}, {0, 0, 1}, {1, 1, 1}},
    };
    NDArray conn = NDArray::Uninit(DType::Int64, {6 * n * n * n, 4});
    std::int64_t* c = conn.As<std::int64_t>();
    for (std::size_t k = 0, t = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                for (int s = 0; s < 6; ++s, ++t)
                    for (int v = 0; v < 4; ++v)
                        c[t * 4 + v] =
                            ft_node(n, i + tets[s][v][0], j + tets[s][v][1], k + tets[s][v][2]);
    Mesh m;
    m.AssignPoints(ft_grid_points(n, jitter));
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

inline std::array<std::int64_t, 8> ft_hex(std::size_t n, std::size_t i, std::size_t j,
                                          std::size_t k) {
    return {ft_node(n, i, j, k),
            ft_node(n, i + 1, j, k),
            ft_node(n, i + 1, j + 1, k),
            ft_node(n, i, j + 1, k),
            ft_node(n, i, j, k + 1),
            ft_node(n, i + 1, j, k + 1),
            ft_node(n, i + 1, j + 1, k + 1),
            ft_node(n, i, j + 1, k + 1)};
}

// Hexahedra; with `polyhedra`, every other one is a polyhedron block of the
// same six quads instead, so hexahedron and polyhedron faces must cancel.
inline Mesh ft_hex_grid(std::size_t n, double jitter, bool polyhedra) {
    std::vector<std::int64_t> hexes;
    std::vector<std::vector<std::vector<std::int64_t>>> polys;
    static const int quads[6][4] = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
                                    {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) {
                const auto h = ft_hex(n, i, j, k);
                if (polyhedra && (i + j + k) % 2 == 1) {
                    std::vector<std::vector<std::int64_t>> faces;
                    for (const auto& q : quads)
                        faces.push_back({h[q[0]], h[q[1]], h[q[2]], h[q[3]]});
                    polys.push_back(std::move(faces));
                } else {
                    hexes.insert(hexes.end(), h.begin(), h.end());
                }
            }
    NDArray conn = NDArray::Uninit(DType::Int64, {hexes.size() / 8, 8});
    std::memcpy(conn.Data(), hexes.data(), hexes.size() * sizeof(std::int64_t));
    Mesh m;
    m.AssignPoints(ft_grid_points(n, jitter));
    m.AddCellBlock("hexahedron", std::move(conn));
    if (!polys.empty())
        m.AddPolyhedronBlock("polyhedron", std::move(polys));
    return m;
}

// Two pentagonal prisms stacked on a shared pentagon: facet keys of five ids.
inline Mesh ft_pentagon_prisms() {
    std::vector<std::vector<double>> pts;
    const double xy[5][2] = {{1.0, 0.0}, {0.25, 0.75}, {-0.75, 0.5}, {-0.75, -0.5}, {0.25, -0.75}};
    for (int layer = 0; layer < 3; ++layer)
        for (const auto& v : xy)
            pts.push_back({v[0], v[1], 0.5 * layer});
    std::vector<std::vector<std::vector<std::int64_t>>> cells;
    for (std::int64_t l = 0; l < 2; ++l) {
        const std::int64_t b = 5 * l, t = 5 * (l + 1);
        std::vector<std::vector<std::int64_t>> faces;
        faces.push_back({b + 4, b + 3, b + 2, b + 1, b + 0});
        faces.push_back({t + 0, t + 1, t + 2, t + 3, t + 4});
        for (std::int64_t e = 0; e < 5; ++e)
            faces.push_back({b + e, b + (e + 1) % 5, t + (e + 1) % 5, t + e});
        cells.push_back(std::move(faces));
    }
    Mesh m;
    m.AssignPoints(mt::points_from(pts));
    m.AddPolyhedronBlock("polyhedron", std::move(cells));
    return m;
}

inline std::uint64_t ft_mesh(const Mesh& rM) {
    MeshDigest d;
    d.Of(rM);
    return d.Value();
}

}  // namespace golden
