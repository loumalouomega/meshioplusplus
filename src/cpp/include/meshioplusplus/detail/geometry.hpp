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
 * `Vec3` is a plain `std::array<double, 3>`; the `vec3_*`/`triple_product`/
 * `det3` primitives are tiny per-operation arithmetic invoked repeatedly
 * inside per-cell hot loops (surface normals, quality metrics), so they stay
 * `inline` here rather than moving to a `.cpp` — at that call frequency,
 * removing the function-call boundary is what lets the compiler fold/
 * vectorize the surrounding loop. `read_point`/`read_corner_coords`/
 * `cell_corner_count` are each called once per cell (not once per scalar), so
 * their bodies live in `src/cpp/src/detail/geometry.cpp` instead.
 *
 * Coordinates are pulled out of an `NDArray` through `detail::read_double`, so
 * these work regardless of the point/connectivity dtype and under every mesh
 * backend.
 */

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/cell_type.hpp"
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
MESHIOPLUSPLUS_API Vec3 read_point(const NDArray& rPoints, std::size_t pointDim, std::int64_t nodeId);

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
MESHIOPLUSPLUS_API void read_corner_coords(const NDArray& rPoints, std::size_t pointDim, const NDArray& rConn,
                        std::size_t rowOffset, std::size_t n, std::vector<Vec3>& rOut);

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
MESHIOPLUSPLUS_API int cell_corner_count(CellType type);

}  // namespace detail
}  // namespace meshioplusplus
