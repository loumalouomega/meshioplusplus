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
 * @file cell_faces.hpp
 * @brief Per-cell-type boundary-face topology tables (local node indices of
 * each face of a 3D volume cell), used by the skin extractor (`skin.hpp`).
 *
 * Every face row lists **corner nodes first** (wound so the face normal
 * points *outward*, away from the element interior), then the mid-edge nodes
 * (mid of corner k → corner k+1), then the face-center node — i.e. each row
 * is itself a valid meshio/VTK `triangle6`/`quad8`/`quad9` node ordering.
 *
 * Node numbering is meshio's (= VTK's). Conventions baked in, matching the
 * rest of this repo (see `openfoam.cpp`'s `build_*` orientation checks and
 * the `tests/python/helpers.py` fixtures):
 *  - tetra: base `(0,1,2)` normal points toward apex 3 → outward base is
 *    `(0,2,1)`.
 *  - hexahedron: base `(0,1,2,3)` normal points toward the top `(4,5,6,7)`.
 *  - wedge: base `(0,1,2)` normal points toward the top `(3,4,5)` (the
 *    gmsh-like layout this codebase uses; do NOT trust vtkWedge's own
 *    doc-comment/face array, whose orientation is famously inconsistent).
 *  - pyramid: base `(0,1,2,3)` normal points toward apex 4.
 *  - wedge15 mid-node numbering is pure VTK_QUADRATIC_WEDGE — EnSight's
 *    penta15 involution (see CLAUDE.md) must NOT be applied here.
 *
 * The outward winding of every row is enforced by a gtest invariant
 * (`tests/cpp/test_skin.cpp`): on the reference element, the Newell normal
 * of each face's corner ring must point away from the cell centroid.
 *
 * KEEP IN SYNC: `src/python/meshioplusplus/_skin.py` carries the Python twin of
 * these tables for the pure-Python fallback — any change here must be
 * mirrored there.
 *
 * `cell_faces` is looked up once per cell (not per node/scalar), so its table
 * lookup lives in `src/cpp/src/detail/cell_faces.cpp` rather than inline here.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief One boundary face of a volume cell: its meshio face type and the
 * local node indices into the cell's connectivity row.
 */
struct CellFaceDef {
    /// Face cell type: `Triangle`/`Quad`/`Triangle6`/`Quad8`/`Quad9`.
    CellType mFaceType;
    /// Number of corner nodes (3 or 4) — the leading entries of `mNodes`.
    std::uint8_t mNumCorners;
    /// Total node count of the face (3, 4, 6, 8, or 9).
    std::uint8_t mNumNodes;
    /// Local node indices, corners first (outward winding), then mid-edge
    /// nodes, then the face center; only the first `mNumNodes` are valid.
    std::array<std::uint8_t, 9> mNodes;
};

/**
 * @brief The boundary faces of a volume cell type, or an empty list for
 * types the skin extractor does not support.
 * @param VolumeType The volume cell type to query.
 * @return Reference to the process-wide face table (empty if unsupported).
 */
MESHIOPLUSPLUS_API const std::vector<CellFaceDef>& cell_faces(CellType VolumeType);

/**
 * @brief Whether the skin extractor supports a volume cell type.
 * @param Type The cell type to test.
 * @return `true` when `cell_faces(Type)` is non-empty.
 */
inline bool skin_supported(CellType Type) {
    return !cell_faces(Type).empty();
}

}  // namespace detail
}  // namespace meshioplusplus
