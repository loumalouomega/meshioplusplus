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
 * @file detail/surface_fem.hpp
 * @brief The P1 finite elements of a triangle surface: its cotangent stiffness
 * and lumped (one-third area) mass, shared by `geodesic_distance` and
 * `laplacian_eigenmodes` so the two read the same operators.
 *
 * `build_surface_fem` takes the mesh's `triangle` blocks and refuses, naming
 * the operation that fixes it, a volume cell (`extract_surface`), a quad or
 * polygon (`convert_cells(mode="simplexify")`) or a higher-order surface cell
 * (`linearize`); line and vertex blocks are ignored. A triangle whose doubled
 * area is not above `1e-14` times its longest edge squared is degenerate and
 * left out (counted). Everything is serial and in input order, so the
 * operators are the same on every backend.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus {
namespace detail {

/// One surface triangle and its cotangents (at corners 0, 1, 2).
struct SurfaceFemTriangle {
    std::array<std::int64_t, 3> mNodes;
    std::array<double, 3> mCot;
    double mDoubleArea;
};

struct SurfaceFem {
    /// Every mesh point, as three coordinates (2-D points get z = 0).
    std::vector<std::array<double, 3>> mPoints;
    /// The non-degenerate triangles, in input order.
    std::vector<SurfaceFemTriangle> mTriangles;
    /// Triangles left out for (near) zero area.
    std::int64_t mNumDegenerate = 0;
    /// The mean edge length of the kept triangles (each edge once per triangle).
    double mMeanEdge = 0.0;
};

/**
 * @brief The surface of @p rMesh as P1 triangles (see the file comment).
 * @param rOp the operation's name, which every refusal starts with.
 * @throws std::invalid_argument on an out-of-scope cell, a triangle naming a
 *         point outside the mesh, no triangle, or every triangle degenerate.
 */
SurfaceFem build_surface_fem(const Mesh& rMesh, const std::string& rOp);

/**
 * @brief The lumped mass (one third of each triangle's area per corner, into
 * @p rMass, sized to the point count) and the cotangent stiffness triplets
 * (appended to @p rStiffness), in triangle then corner order.
 */
void surface_fem_assemble(const SurfaceFem& rFem, std::vector<double>& rMass,
                          std::vector<SparseTriplet>& rStiffness);

}  // namespace detail
}  // namespace meshioplusplus
