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
 * @file cell_edges.hpp
 * @brief Per-cell-type boundary-edge topology tables (local node indices of
 * each edge of a 2D surface cell), the 2D analogue of `cell_faces.hpp`, used
 * by the surface/boundary extractor (`operations/surface.hpp`) when the input
 * is a surface (max topological dimension 2) mesh.
 *
 * Every edge row lists its two **corner nodes first**, wound consistently with
 * the cell's local corner ring (corner k -> corner k+1), then the mid-edge
 * node for quadratic cells — i.e. each row is itself a valid meshio/VTK
 * `line`/`line3` node ordering.
 *
 * Node numbering is meshio's (= VTK's): for `triangle6` mid node 3 = mid(0,1),
 * 4 = mid(1,2), 5 = mid(2,0); for `quad8`/`quad9` mid node 4 = mid(0,1),
 * 5 = mid(1,2), 6 = mid(2,3), 7 = mid(3,0) (quad9 node 8 is the face center,
 * on no edge).
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief One boundary edge of a surface cell: its meshio edge type and the
 * local node indices into the cell's connectivity row.
 */
struct CellEdgeDef {
    /// Edge cell type: `Line` or `Line3`.
    CellType mEdgeType;
    /// Number of corner nodes (always 2) — the leading entries of `mNodes`.
    std::uint8_t mNumCorners;
    /// Total node count of the edge (2 for `Line`, 3 for `Line3`).
    std::uint8_t mNumNodes;
    /// Local node indices, corners first, then the mid node; only the first
    /// `mNumNodes` are valid.
    std::array<std::uint8_t, 3> mNodes;
};

/**
 * @brief The boundary edges of a surface cell type, or an empty list for types
 * the surface extractor does not support.
 * @param SurfaceType The surface cell type to query.
 * @return Reference to the process-wide edge table (empty if unsupported).
 */
inline const std::vector<CellEdgeDef>& cell_edges(CellType SurfaceType) {
    using CT = CellType;
    static const std::vector<CellEdgeDef> empty = {};
    static const std::vector<CellEdgeDef> triangle = {
        {CT::Line, 2, 2, {0, 1}},
        {CT::Line, 2, 2, {1, 2}},
        {CT::Line, 2, 2, {2, 0}},
    };
    static const std::vector<CellEdgeDef> triangle6 = {
        {CT::Line3, 2, 3, {0, 1, 3}},
        {CT::Line3, 2, 3, {1, 2, 4}},
        {CT::Line3, 2, 3, {2, 0, 5}},
    };
    static const std::vector<CellEdgeDef> quad = {
        {CT::Line, 2, 2, {0, 1}},
        {CT::Line, 2, 2, {1, 2}},
        {CT::Line, 2, 2, {2, 3}},
        {CT::Line, 2, 2, {3, 0}},
    };
    static const std::vector<CellEdgeDef> quad8 = {
        {CT::Line3, 2, 3, {0, 1, 4}},
        {CT::Line3, 2, 3, {1, 2, 5}},
        {CT::Line3, 2, 3, {2, 3, 6}},
        {CT::Line3, 2, 3, {3, 0, 7}},
    };
    switch (SurfaceType) {
        case CT::Triangle:
            return triangle;
        case CT::Triangle6:
            return triangle6;
        case CT::Quad:
            return quad;
        case CT::Quad8:
        case CT::Quad9:  // same edges as quad8; node 8 (center) is on no edge
            return quad8;
        default:
            return empty;
    }
}

/**
 * @brief Whether the surface extractor supports a surface cell type's edges.
 * @param Type The cell type to test.
 * @return `true` when `cell_edges(Type)` is non-empty.
 */
inline bool surface_edge_supported(CellType Type) {
    return !cell_edges(Type).empty();
}

}  // namespace detail
}  // namespace meshioplusplus
