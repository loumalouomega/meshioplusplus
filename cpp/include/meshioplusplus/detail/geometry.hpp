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
 * @file geometry.hpp
 * @brief Small, dependency-free 3D vector primitives and mesh-coordinate
 * readers shared by the mesh-operations layer (`operations/quality.cpp`,
 * `operations/surface.cpp`).
 *
 * `Vec3` is a plain `std::array<double, 3>`; every helper is a free inline
 * function in `meshioplusplus::detail` (header-only, so it is exempt from the
 * anonymous-namespace unique-prefix rule the amalgamation imposes on `.cpp`
 * translation units). Coordinates are pulled out of an `NDArray` through
 * `detail::read_double`, so these work regardless of the point/connectivity
 * dtype and under every mesh backend.
 */

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/// A point/vector in 3D. 2D meshes are padded with z = 0 on read.
using Vec3 = std::array<double, 3>;

/** @brief Component-wise difference `a - b`. */
inline Vec3 vec3_sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

/** @brief Component-wise sum `a + b`. */
inline Vec3 vec3_add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

/** @brief Scalar multiple `s * a`. */
inline Vec3 vec3_scale(const Vec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

/** @brief Dot product `a . b`. */
inline double vec3_dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/** @brief Cross product `a x b`. */
inline Vec3 vec3_cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

/** @brief Squared Euclidean length `a . a`. */
inline double vec3_norm_sq(const Vec3& a) {
    return vec3_dot(a, a);
}

/** @brief Euclidean length `|a|`. */
inline double vec3_norm(const Vec3& a) {
    return std::sqrt(vec3_norm_sq(a));
}

/**
 * @brief Unit vector along `a`, or the zero vector when `|a| < eps`.
 * @param a Vector to normalize.
 * @param eps Length below which `a` is treated as degenerate.
 * @return `a / |a|`, or `{0,0,0}` if `|a| < eps`.
 */
inline Vec3 vec3_normalize(const Vec3& a, double eps = 1e-300) {
    const double n = vec3_norm(a);
    if (n < eps)
        return {0.0, 0.0, 0.0};
    return vec3_scale(a, 1.0 / n);
}

/** @brief Scalar triple product `a . (b x c)` (signed volume of the parallelepiped). */
inline double triple_product(const Vec3& a, const Vec3& b, const Vec3& c) {
    return vec3_dot(a, vec3_cross(b, c));
}

/** @brief Determinant of the 3x3 matrix whose columns are `c0`, `c1`, `c2`. */
inline double det3(const Vec3& c0, const Vec3& c1, const Vec3& c2) {
    return triple_product(c0, c1, c2);
}

/**
 * @brief Reads global point @p nodeId as a `Vec3`, padding z = 0 when the mesh
 * is 2D (`pointDim == 2`).
 * @param rPoints The `(num_points, pointDim)` point array.
 * @param pointDim Spatial dimension of the points (2 or 3).
 * @param nodeId Global point index.
 * @return The point's coordinates, with unused components set to 0.
 */
inline Vec3 read_point(const NDArray& rPoints, std::size_t pointDim, std::int64_t nodeId) {
    Vec3 p = {0.0, 0.0, 0.0};
    const std::size_t base = static_cast<std::size_t>(nodeId) * pointDim;
    const std::size_t k = pointDim < 3 ? pointDim : 3;
    for (std::size_t c = 0; c < k; ++c)
        p[c] = read_double(rPoints, base + c);
    return p;
}

/**
 * @brief Reads the first @p n connectivity entries of one cell row into @p rOut
 * as `Vec3` coordinates (used to gather a cell's corner nodes).
 * @param rPoints The point array.
 * @param pointDim Spatial dimension of the points.
 * @param rConn The block connectivity array.
 * @param rowOffset Flat offset of the cell's row (`cell * nodes_per_cell`).
 * @param n Number of leading entries to read (the corner count).
 * @param rOut Cleared and filled with @p n coordinates.
 */
inline void read_corner_coords(const NDArray& rPoints, std::size_t pointDim, const NDArray& rConn,
                               std::size_t rowOffset, std::size_t n, std::vector<Vec3>& rOut) {
    rOut.clear();
    rOut.reserve(n);
    for (std::size_t k = 0; k < n; ++k)
        rOut.push_back(read_point(rPoints, pointDim, read_int(rConn, rowOffset + k)));
}

/**
 * @brief Number of corner (linear-parent) nodes of a cell type. Corners are
 * always the leading connectivity entries in meshio/VTK ordering, so a
 * quadratic cell can be reduced to its linear parent by reading the first
 * `cell_corner_count(type)` nodes.
 * @param type The cell type to query.
 * @return The corner count, or 0 for variable-node-count / unsupported types
 *         (`Polygon`, `Polyhedron`, the VTK Lagrange family, `Custom`) — the
 *         caller must skip those.
 */
inline int cell_corner_count(CellType type) {
    switch (type) {
        case CellType::Vertex:
            return 1;
        case CellType::Line:
        case CellType::Line3:
        case CellType::Line4:
        case CellType::Line5:
        case CellType::Line6:
        case CellType::Line7:
        case CellType::Line8:
        case CellType::Line9:
        case CellType::Line10:
        case CellType::Line11:
            return 2;
        case CellType::Triangle:
        case CellType::Triangle6:
        case CellType::Triangle10:
        case CellType::Triangle15:
        case CellType::Triangle21:
        case CellType::Triangle28:
        case CellType::Triangle36:
        case CellType::Triangle45:
        case CellType::Triangle55:
        case CellType::Triangle66:
            return 3;
        case CellType::Quad:
        case CellType::Quad8:
        case CellType::Quad9:
        case CellType::Quad16:
        case CellType::Quad25:
        case CellType::Quad36:
        case CellType::Quad49:
        case CellType::Quad64:
        case CellType::Quad81:
        case CellType::Quad100:
        case CellType::Quad121:
            return 4;
        case CellType::Tetra:
        case CellType::Tetra10:
        case CellType::Tetra20:
        case CellType::Tetra35:
        case CellType::Tetra56:
        case CellType::Tetra84:
        case CellType::Tetra120:
        case CellType::Tetra165:
        case CellType::Tetra220:
        case CellType::Tetra286:
            return 4;
        case CellType::Hexahedron:
        case CellType::Hexahedron20:
        case CellType::Hexahedron24:
        case CellType::Hexahedron27:
        case CellType::Hexahedron64:
        case CellType::Hexahedron125:
        case CellType::Hexahedron216:
        case CellType::Hexahedron343:
        case CellType::Hexahedron512:
        case CellType::Hexahedron729:
        case CellType::Hexahedron1000:
        case CellType::Hexahedron1331:
            return 8;
        case CellType::Wedge:
        case CellType::Wedge15:
        case CellType::Wedge18:
        case CellType::Wedge40:
        case CellType::Wedge75:
        case CellType::Wedge126:
        case CellType::Wedge196:
        case CellType::Wedge288:
        case CellType::Wedge405:
        case CellType::Wedge550:
            return 6;
        case CellType::Pyramid:
        case CellType::Pyramid13:
        case CellType::Pyramid14:
            return 5;
        default:  // Polygon, Polyhedron, VtkLagrange*, Custom
            return 0;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
