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
 * @file operations/subdivide_surface.hpp
 * @brief Smooth subdivision surfaces: Loop (triangles) and Catmull-Clark
 * (triangles, quads, polygons -> quads) with Hoppe-style crease and boundary
 * preservation (roadmap 9.2).
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Attached Int64 `cell_data` array recording each child's input cell index.
inline constexpr const char* kSubdivideParentCellName = "subdivide:parent_cell";

/// Smooth subdivision scheme for `subdivide_surface`.
enum class SubdivisionScheme : std::int32_t {
    Loop = 0,         ///< Loop triangle subdivision (1987); triangles only.
    CatmullClark = 1  ///< Catmull-Clark subdivision (1978); polygonal -> quads.
};

/// Parse `"loop"` or `"catmull_clark"` / `"catmull-clark"` / `"cc"`.
MESHIOPLUSPLUS_API SubdivisionScheme subdivision_scheme_from_name(const std::string& rName);

/// Options for `subdivide_surface`.
struct SubdivideSurfaceOptions {
    /// Subdivision scheme (`Loop` or `CatmullClark`).
    SubdivisionScheme mScheme = SubdivisionScheme::Loop;
    /// Number of subdivision levels (>= 0; 0 returns a copy).
    int mLevels = 1;
    /// Dihedral angle threshold in degrees for automatic sharp creases (<= 0 or >= 180 disables).
    double mCreaseAngle = 60.0;
    /// Apply crease/boundary rules to open boundary edges (corners pinned).
    bool mPreserveBoundary = true;
    /// Interpolate floating-point `point_data` using the subdivision weights.
    bool mInterpolatePointData = true;
    /// Attach Int64 `subdivide:parent_cell` recording each child's input cell index.
    bool mRecordParentIds = false;
};

/// Result of `subdivide_surface`.
struct SubdivideSurfaceResult {
    Mesh mMesh;
    std::int64_t mNumCreaseEdges = 0;
    std::int64_t mNumBoundaryEdges = 0;
    std::vector<std::string> mWarnings;
};

/// Apply smooth Loop or Catmull-Clark subdivision to a surface mesh.
MESHIOPLUSPLUS_API SubdivideSurfaceResult
subdivide_surface(const Mesh& rMesh, const SubdivideSurfaceOptions& rOptions = {});

}  // namespace meshioplusplus
