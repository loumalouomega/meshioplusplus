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
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Project cell_edges()'s rows down to bare corner pairs, so the 2D edge order
// keeps its single owner in cell_edges.cpp.
std::vector<CellEdgePair> subdiv_from_cell_edges(CellType Type) {
    std::vector<CellEdgePair> edges;
    for (const CellEdgeDef& ed : cell_edges(Type))
        edges.push_back({ed.mNodes[0], ed.mNodes[1]});
    return edges;
}

}  // namespace

const std::vector<CellEdgePair>& cell_refine_edges(CellType Type) {
    using CT = CellType;
    static const std::vector<CellEdgePair> empty = {};
    static const std::vector<CellEdgePair> line = {{0, 1}};
    // triangle6 nodes 3-5, quad8/quad9 nodes 4-7.
    static const std::vector<CellEdgePair> triangle = subdiv_from_cell_edges(CT::Triangle);
    static const std::vector<CellEdgePair> quad = subdiv_from_cell_edges(CT::Quad);
    // tetra10 nodes 4-9.
    static const std::vector<CellEdgePair> tetra = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
    // hexahedron20 nodes 8-19: bottom ring, top ring, then the four verticals.
    static const std::vector<CellEdgePair> hexahedron = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
        {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7},
    };
    // wedge15 nodes 6-14: bottom triangle, top triangle, then the three verticals.
    static const std::vector<CellEdgePair> wedge = {
        {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5},
    };
    // pyramid13 nodes 5-12: base ring, then the four apex edges.
    static const std::vector<CellEdgePair> pyramid = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4},
    };
    switch (Type) {
        case CT::Line:
            return line;
        case CT::Triangle:
            return triangle;
        case CT::Quad:
            return quad;
        case CT::Tetra:
            return tetra;
        case CT::Hexahedron:
            return hexahedron;
        case CT::Wedge:
            return wedge;
        case CT::Pyramid:
            return pyramid;
        default:
            return empty;
    }
}

const std::vector<CellQuadFace>& cell_refine_quad_faces(CellType Type) {
    using CT = CellType;
    static const std::vector<CellQuadFace> empty = {};
    // quad9 node 8: the cell's own face.
    static const std::vector<CellQuadFace> quad = {{0, 1, 2, 3}};
    // hexahedron27 nodes 20-25, in cell_faces.hpp's own hexahedron27 order
    // (vtkTriQuadraticHexahedron::Faces): x-min, x-max, y-min, y-max, bottom,
    // top -- row k is node 20+k, which `CellSubdivision.
    // QuadFacesAgreeWithCellFaces` enforces against that table. Reordered in
    // v9.9.0 together with cell_faces.cpp's 20/22/23 correction and
    // refine_templates.cpp's absolute-index references; see cell_faces.cpp.
    static const std::vector<CellQuadFace> hexahedron = {
        {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 1, 2, 3}, {4, 5, 6, 7},
    };
    // wedge18 nodes 15-17: the three quad side faces (the two triangle faces
    // gain no node, which is why a wedge needs no body node either -- see
    // operations/refine.cpp).
    static const std::vector<CellQuadFace> wedge = {
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {2, 0, 3, 5},
    };
    switch (Type) {
        case CT::Quad:
            return quad;
        case CT::Hexahedron:
            return hexahedron;
        case CT::Wedge:
            return wedge;
        default:
            return empty;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
