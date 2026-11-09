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
 * @file operations/clip.hpp
 * @brief Cut a surface by a plane, keep one side and cap it (roadmap 7.8.3;
 * `cut` is the level-set operation, so this one is `clip`).
 *
 * The surface (triangles, quads and polygons, simplexified as `slice`
 * simplexifies them; a volume block is refused by name) is cut where
 * `d = (x - origin) . normal` changes sign, `d = 0` counting as above, at the
 * crossings `slice` computes (the sorted-edge formula, so the cut boundary is
 * the section `slice` returns, bit for bit). A triangle wholly on the kept
 * side stays, one across the plane is split there (a triangle or two), and
 * the rest go. With `mCap`, each closed loop of the cut is triangulated by
 * the planar CDT (`cdt2d_triangulate_polygon`, no new points; nested loops
 * by the even-odd rule) on the crossings themselves, so a closed part stays
 * closed, its cap facing out of the kept side.
 *
 * Point data is interpolated at the crossings (Float64) and kept elsewhere;
 * cell data is the parent cell's, and NaN (floating) or -1 (integer) on the
 * cap, which the Int64 cell data `clip:cap` marks (1). The result is one
 * `triangle` block.
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

struct ClipOptions {
    std::array<double, 3> mOrigin{0.0, 0.0, 0.0};
    std::array<double, 3> mNormal{0.0, 0.0, 1.0};
    bool mKeepBelow = true;  ///< keep d < 0 (else d >= 0)
    bool mCap = true;
};

struct ClipResult {
    Mesh mMesh;
    std::int64_t mNumKept = 0;   ///< triangles kept whole
    std::int64_t mNumSplit = 0;  ///< triangles split at the plane
    std::int64_t mNumCapTriangles = 0;
    std::int64_t mNumLoops = 0;       ///< closed loops capped
    std::int64_t mNumOpenChains = 0;  ///< open chains of the cut (not capped)
    std::vector<std::string> mWarnings;
};

/**
 * @brief Clip @p rMesh by a plane (see the file comment).
 * @throws std::invalid_argument on a zero normal or a mesh with volume cells.
 */
MESHIOPLUSPLUS_API ClipResult clip(const Mesh& rMesh, const ClipOptions& rOptions = {});

}  // namespace meshioplusplus
