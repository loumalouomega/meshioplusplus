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
 * @file operations/flip_normals.hpp
 * @brief Reverse the winding of all cells, or of a cell region (roadmap
 * 7.8.2), with no point moved.
 *
 * A cell is flipped by its type's orientation-reversing node permutation (the
 * table `mirror` uses, `detail/cell_flip.hpp`): a triangle's or polygon's
 * ring is reversed after its first node, a solid's nodes reordered so its
 * volume changes sign, a polyhedron's face rings reversed. `Side` region
 * entries of a flipped cell are moved to the facet the flip gives them; every
 * array is kept as it is. Flipping twice gives the input back bit for bit.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

struct FlipNormalsResult {
    Mesh mMesh;
    std::int64_t mNumFlipped = 0;  ///< cells whose winding was reversed
};

/**
 * @brief Reverse the winding of @p rMesh's cells, or those of the cell region
 * @p rRegion (empty: every cell).
 * @throws std::invalid_argument on an unknown region.
 * @throws Unsupported on a cell type with no orientation flip.
 */
MESHIOPLUSPLUS_API FlipNormalsResult flip_normals(const Mesh& rMesh,
                                                  const std::string& rRegion = "");

}  // namespace meshioplusplus
