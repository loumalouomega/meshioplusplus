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
 * @file operations/check_tetrahedra.hpp
 * @brief `check_tetrahedra`: are the tetrahedra of a mesh usable?
 *
 * Counts, over the `tetra` and `tetra10` cells (a quadratic tetrahedron by its
 * four corners), the tetrahedra that are inverted (a negative volume beyond
 * the tolerance), degenerate (a volume within it) and the pairs that overlap
 * across a shared face. The tolerance of a tetrahedron is `zero_volume_factor`
 * times the cube of its longest edge plus the effect of six significant
 * digits of written coordinates, so a sliver of a fine mesh is not flagged
 * because of a large coordinate offset. Two tetrahedra overlap when they
 * share a face (and only they use it) and their opposite nodes lie on the same
 * side of it; a conforming mesh of a valid solid has none, whatever its
 * materials. Faces used by three or more tetrahedra are not examined.
 *
 * The mesh is not changed. The same tolerance and overlap count drive the
 * volume fix of `octree_mesh`.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `check_tetrahedra`.
struct TetrahedraCheckOptions {
    /// A volume under this fraction of the longest edge cubed (plus the effect of the written
    /// coordinates' precision) counts as zero; must be finite and non-negative.
    double mZeroVolumeFactor = 1.0e-6;
};

/// The verdict and the counts behind it.
struct TetrahedraCheckResult {
    std::int64_t mNumTetrahedra = 0;  ///< `tetra` and `tetra10` cells examined
    std::int64_t mNumInverted = 0;    ///< volume below minus the tolerance
    std::int64_t mNumDegenerate = 0;  ///< volume within the tolerance
    std::int64_t mNumOverlaps = 0;    ///< pairs overlapping across a shared face
    double mMinVolume = 0.0;          ///< smallest signed volume (0 without tetrahedra)
    double mMaxVolume = 0.0;          ///< largest signed volume (0 without tetrahedra)
    /// No inverted or degenerate tetrahedron and no overlap.
    bool mValid = true;
    /// What the check could not examine, one entry per line.
    std::vector<std::string> mWarnings;
};

MESHIOPLUSPLUS_API TetrahedraCheckResult
check_tetrahedra(const Mesh& rMesh, const TetrahedraCheckOptions& rOptions = {});

}  // namespace meshioplusplus
