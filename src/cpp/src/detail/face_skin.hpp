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
 * @file detail/face_skin.hpp
 * @brief The skin of a cell selection, without the whole face table.
 *
 * A **core-private** header (the `detail/slot_runs.hpp` precedent).
 * `feature_edges` used to build `build_global_faces`' whole face mesh -- every
 * face, its owner and neighbour, the cells' signed face lists -- only to keep
 * the faces with exactly one selected side. This builds the same cell rings
 * and groups them the same way, and keeps the skin alone (roadmap 5.4.4).
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// The skin faces of a selection, as `build_global_faces` numbers and winds them.
struct SelectionSkin {
    std::size_t mNumCells = 0;  ///< volume cells (the face table's compact cells)
    std::vector<std::int64_t> mFaceStart{0};
    std::vector<std::int64_t> mFaceNodes;  ///< each face as its owner winds it
    std::vector<std::uint8_t> mReversed;   ///< wound out of the neighbour instead

    std::size_t NumFaces() const { return mFaceStart.size() - 1; }
};

/**
 * @brief The faces of the volume cells (polyhedra and the skin-supported
 * types) whose owner and neighbour in `build_global_faces(rMesh)` differ in
 * `rSelected` (indexed by global, block-major cell), in that table's face
 * order. A face is wound out of its selected side: the owner's stored ring,
 * or its reverse when `mReversed` is set.
 */
SelectionSkin build_selection_skin(const Mesh& rMesh, const std::vector<std::uint8_t>& rSelected);

}  // namespace detail
}  // namespace meshioplusplus
