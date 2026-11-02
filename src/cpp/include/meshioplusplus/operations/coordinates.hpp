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
 * @file operations/coordinates.hpp
 * @brief Cylindrical and spherical components about an axis (roadmap 7.6.1):
 * the positions, and vector and tensor fields, of a rotating machine or a
 * pressure vessel read in the components its engineers think in.
 *
 * `convert_coordinates(mesh, {system, axis, origin, reference, arrays})`
 * leaves the points where they are and
 *
 * - records each point's coordinates in the system as point data
 *   (`coords:cylindrical` = `(r, theta, z)`, `coords:spherical` =
 *   `(r, theta, phi)`, see below);
 * - rewrites each named vector (`(n, 3)`) or tensor (`(n, 9)` or `(n, 3, 3)`,
 *   row-major) array, point or cell data, in the local basis of its point (a
 *   cell's at its centroid): `v' = R v` and `T' = R T R^T`, `R`'s rows the
 *   local unit vectors. With `mInverse` the arrays are read as local
 *   components and rewritten in Cartesian ones (`R^T v`, `R^T T R`), so the two
 *   directions round trip.
 *
 * **The systems.** With `a` the unit axis through the origin `o`, `d = p - o`,
 * `z = d . a`, the radial vector `rho = d - z a` and `e1` the reference
 * direction projected off the axis (`(1, 0, 0)` by default, `(0, 1, 0)` when
 * that is along the axis):
 *
 * - cylindrical: `(r, theta, z)` with `r = |rho|`, `theta` the angle from
 *   `e1` to `rho` about `a` in `(-pi, pi]`; basis `e_r = rho / r`,
 *   `e_theta = a x e_r`, `e_z = a`.
 * - spherical (physics convention): `(r, theta, phi)` with `r = |d|`,
 *   `theta` the polar angle from `a` in `[0, pi]`, `phi` the cylindrical
 *   `theta`; basis `e_r = d / r`, `e_theta = cos(theta) e_rho - sin(theta) a`,
 *   `e_phi = a x e_rho`.
 *
 * On the axis (`r = 0` cylindrically, or the spherical origin) the radial
 * direction is undefined: `e1` stands in for it, the angle is 0, and the point
 * is counted in `mNumSingular`.
 *
 * Data-only: points, connectivity, regions and every other array pass through
 * unchanged. A converted array becomes Float64. Both engines (this one and the
 * numpy twin) compute the same values to rounding.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data: Float64 `(n, 3)`, `(r, theta, z)` about the axis.
inline constexpr const char* kCylindricalCoordsName = "coords:cylindrical";
/// Point data: Float64 `(n, 3)`, `(r, theta, phi)` about the axis.
inline constexpr const char* kSphericalCoordsName = "coords:spherical";

/// The curvilinear system.
enum class CoordinateSystem : std::int32_t {
    Cylindrical = 0,
    Spherical = 1,
};

/// What `convert_coordinates` converts, and about what.
struct CoordinateOptions {
    CoordinateSystem mSystem = CoordinateSystem::Cylindrical;
    /// The axis direction (normalised; must not be zero).
    std::array<double, 3> mAxis{0.0, 0.0, 1.0};
    /// A point on the axis (the spherical centre).
    std::array<double, 3> mOrigin{0.0, 0.0, 0.0};
    /// The direction of angle 0, projected off the axis.
    std::array<double, 3> mReference{1.0, 0.0, 0.0};
    /// Vector or tensor arrays (point or cell data) to rewrite in place.
    std::vector<std::string> mArrays;
    /// Read `mArrays` as local components and rewrite them as Cartesian ones.
    bool mInverse = false;
    /// Attach `coords:cylindrical` / `coords:spherical`.
    bool mRecordCoordinates = true;
};

/// What `convert_coordinates` did.
struct CoordinateResult {
    /// The mesh with the coordinates recorded and the arrays rewritten.
    Mesh mMesh;
    /// Arrays rewritten (an array present as point and cell data counts twice).
    std::int64_t mNumArrays = 0;
    /// Points on the axis (or at the spherical centre), where `e1` stood in.
    std::int64_t mNumSingular = 0;
};

/** @brief `"cylindrical"` or `"spherical"`; anything else throws naming them. */
MESHIOPLUSPLUS_API CoordinateSystem coordinate_system_from_name(const std::string& rName);

/**
 * @brief Positions, vector and tensor fields in cylindrical or spherical
 * components about an axis (see the file comment).
 * @throws std::invalid_argument on a zero axis, a named array that is in no
 *         data map or is not a vector or a tensor of 3-D components.
 */
MESHIOPLUSPLUS_API CoordinateResult convert_coordinates(const Mesh& rMesh,
                                                        const CoordinateOptions& rOptions);

}  // namespace meshioplusplus
