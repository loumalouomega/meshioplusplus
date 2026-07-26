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
 * @file operations/transform.hpp
 * @brief Dependency-free affine transform of point coordinates: translate,
 * scale (per-axis or uniform), rotate (axis + angle), a general 4x4 affine
 * matrix, and a unit-scale convenience.
 *
 * `transform` returns a new mesh with the transformed points; connectivity,
 * `cell_data`, `field_data`, and (by default) `point_data` are carried through
 * unchanged. Vector/tensor `point_data` can optionally be rotated by the
 * transform's linear part (`rotate_vector_data`, off by default). An
 * orientation-reversing transform (negative determinant of the 3x3 linear
 * block) logs a warning that cell orientation may be flipped.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry.
 */

// System includes
#include <array>
#include <cstddef>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief A 4x4 homogeneous affine transform, stored row-major (element `(r, c)`
 * at `mMatrix[r * 4 + c]`). A point `p` maps to `M * [p.x, p.y, p.z, 1]`.
 */
struct AffineTransform {
    /// Row-major 4x4 matrix; identity by default.
    std::array<double, 16> mMatrix = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
};

/** @brief A pure translation by `(dx, dy, dz)`. */
MESHIOPLUSPLUS_API AffineTransform transform_translation(double dx, double dy, double dz);

/** @brief A per-axis scale `(sx, sy, sz)` about the origin. */
MESHIOPLUSPLUS_API AffineTransform transform_scale(double sx, double sy, double sz);

/** @brief A uniform scale by `factor` about the origin (e.g. unit conversion). */
MESHIOPLUSPLUS_API AffineTransform transform_units(double factor);

/**
 * @brief A rotation of `angle_rad` radians about the axis `(ax, ay, az)`
 * (Rodrigues' formula; the axis is normalized internally).
 * @throws std::invalid_argument if the axis is (near) zero-length.
 */
MESHIOPLUSPLUS_API AffineTransform transform_rotation(double ax, double ay, double az, double angle_rad);

/** @brief Wrap a caller-supplied row-major 4x4 matrix. */
MESHIOPLUSPLUS_API AffineTransform transform_from_matrix(const double* pMatrix16);

/** @brief Compose two transforms: the result applies `rSecond` after `rFirst`. */
MESHIOPLUSPLUS_API AffineTransform transform_compose(const AffineTransform& rSecond, const AffineTransform& rFirst);

/**
 * @brief Apply an affine transform to a mesh's point coordinates.
 * @param rMesh the input mesh (unmodified).
 * @param rXform the affine transform to apply.
 * @param rotate_vector_data when true, rotate vector (trailing dim 3) and
 *        tensor (trailing dim 9) `point_data` by the transform's 3x3 linear
 *        block (`R*v` / `R*A*R^T`); off by default.
 * @return a new mesh with transformed points; connectivity and data preserved.
 */
MESHIOPLUSPLUS_API Mesh transform(const Mesh& rMesh, const AffineTransform& rXform, bool rotate_vector_data = false);

}  // namespace meshioplusplus
