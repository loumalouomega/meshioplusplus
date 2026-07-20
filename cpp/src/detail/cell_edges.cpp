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

// Project includes
#include "meshioplusplus/detail/cell_edges.hpp"

namespace meshioplusplus {
namespace detail {

const std::vector<CellEdgeDef>& cell_edges(CellType SurfaceType) {
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

}  // namespace detail
}  // namespace meshioplusplus
