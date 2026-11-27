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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file operations/deform.hpp
 * @brief Parametric Free-Form Deformation (`deform_ffd`, Sederberg & Parry
 * 1986) and Radial Basis Function control-point morphing (`deform_rbf`,
 * de Boer, van der Schoot & Bijl 2007) (roadmap 9.15).
 *
 * Both preserve connectivity, regions, and all point/cell/field data on
 * volume, surface, polyhedral, and high-order meshes alike.
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

/// Radial basis function kernel for `deform_rbf`.
enum class RbfKernel : std::int32_t {
    /// Global 3-D biharmonic / thin-plate spline `phi(r) = r`.
    Tps = 0,
    /// Compactly-supported Wendland C2 `phi(r) = (1 - r/R)_+^4 (4 r/R + 1)`.
    WendlandC2 = 1,
};

/// Parse `"tps"` or `"wendland_c2"`.
MESHIOPLUSPLUS_API RbfKernel rbf_kernel_from_name(const std::string& rName);

/// Options for `deform_ffd`.
struct FfdOptions {
    /// Control-lattice point counts `(nx, ny, nz)` along X, Y, Z (each >= 2).
    std::array<std::int64_t, 3> mGridDims = {2, 2, 2};
    /// B-spline degrees `(px, py, pz)`; -1 selects Bezier (`n - 1`).
    std::array<std::int64_t, 3> mDegrees = {-1, -1, -1};
    /// Control-point displacements, flat `3 * nx * ny * nz` in `(k, j, i)` order (i fastest).
    std::vector<double> mDisplacements;
    /// If true, use `mBoxMin`/`mBoxMax` as the lattice bounds; else use the mesh bounding box.
    bool mHasBounds = false;
    /// Lattice lower corner when `mHasBounds` is true.
    std::array<double, 3> mBoxMin = {0.0, 0.0, 0.0};
    /// Lattice upper corner when `mHasBounds` is true.
    std::array<double, 3> mBoxMax = {1.0, 1.0, 1.0};
    /// Optional Point or Cell region name to restrict deformation (empty = all points).
    std::string mRegion;
};

/// Options for `deform_rbf`.
struct RbfOptions {
    /// Control-point coordinates, flat `3 * M`.
    std::vector<double> mControlPoints;
    /// Prescribed displacements at the control points, flat `3 * M`.
    std::vector<double> mDisplacements;
    /// Radial basis kernel (`Tps` or `WendlandC2`).
    RbfKernel mKernel = RbfKernel::Tps;
    /// Support radius for `WendlandC2`; 0 selects `0.5 * bounding_box_diagonal`.
    double mRadius = 0.0;
    /// Include linear polynomial reproduction (`1, x, y, z`).
    bool mPolynomial = true;
    /// Non-negative Tikhonov regularization added to the RBF diagonal.
    double mRegularization = 0.0;
    /// Optional Point or Cell region name to restrict deformation (empty = all points).
    std::string mRegion;
};

/// Shared result for `deform_ffd` and `deform_rbf`.
struct DeformResult {
    /// Deformed copy of the input mesh.
    Mesh mMesh;
    /// Number of vertices with non-zero displacement.
    std::int64_t mNumNodesMoved = 0;
    /// Number of 3-D cells with negative signed volume after deformation.
    std::int64_t mNumInverted = 0;
    /// Maximum Euclidean displacement across all vertices.
    double mMaxDisplacement = 0.0;
};

/**
 * @brief Deform a mesh with a tensor-product Bezier / B-spline FFD lattice.
 */
MESHIOPLUSPLUS_API DeformResult deform_ffd(const Mesh& rMesh, const FfdOptions& rOptions);

/**
 * @brief Deform a mesh by RBF interpolation of sparse control-point displacements.
 */
MESHIOPLUSPLUS_API DeformResult deform_rbf(const Mesh& rMesh, const RbfOptions& rOptions);

}  // namespace meshioplusplus
