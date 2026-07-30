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
#include "meshioplusplus/detail/cell_faces.hpp"

namespace meshioplusplus {
namespace detail {

const std::vector<CellFaceDef>& cell_faces(CellType VolumeType) {
    using CT = CellType;
    static const std::vector<CellFaceDef> empty = {};
    static const std::vector<CellFaceDef> tetra = {
        {CT::Triangle, 3, 3, {0, 1, 3}},
        {CT::Triangle, 3, 3, {1, 2, 3}},
        {CT::Triangle, 3, 3, {2, 0, 3}},
        {CT::Triangle, 3, 3, {0, 2, 1}},
    };
    static const std::vector<CellFaceDef> tetra10 = {
        {CT::Triangle6, 3, 6, {0, 1, 3, 4, 8, 7}},
        {CT::Triangle6, 3, 6, {1, 2, 3, 5, 9, 8}},
        {CT::Triangle6, 3, 6, {2, 0, 3, 6, 7, 9}},
        {CT::Triangle6, 3, 6, {0, 2, 1, 6, 5, 4}},
    };
    static const std::vector<CellFaceDef> hexahedron = {
        {CT::Quad, 4, 4, {0, 4, 7, 3}}, {CT::Quad, 4, 4, {1, 2, 6, 5}},
        {CT::Quad, 4, 4, {0, 1, 5, 4}}, {CT::Quad, 4, 4, {3, 7, 6, 2}},
        {CT::Quad, 4, 4, {0, 3, 2, 1}}, {CT::Quad, 4, 4, {4, 5, 6, 7}},
    };
    static const std::vector<CellFaceDef> hexahedron20 = {
        {CT::Quad8, 4, 8, {0, 4, 7, 3, 16, 15, 19, 11}},
        {CT::Quad8, 4, 8, {1, 2, 6, 5, 9, 18, 13, 17}},
        {CT::Quad8, 4, 8, {0, 1, 5, 4, 8, 17, 12, 16}},
        {CT::Quad8, 4, 8, {3, 7, 6, 2, 19, 14, 18, 10}},
        {CT::Quad8, 4, 8, {0, 3, 2, 1, 11, 10, 9, 8}},
        {CT::Quad8, 4, 8, {4, 5, 6, 7, 12, 13, 14, 15}},
    };
    // Mid-face numbering, from vtkTriQuadraticHexahedron::Faces (the authority
    // meshio's hexahedron27 ordering is copied from verbatim -- meshio applies
    // no permutation to this type): 20=(0,4,7,3) x-min, 21=(1,2,6,5) x-max,
    // 22=(0,1,5,4) y-min, 23=(3,7,6,2) y-max, 24=bottom (0,3,2,1),
    // 25=top (4,5,6,7); 26 = body center.
    //
    // Corrected in v9.9.0: 20/22/23 were previously a permuted 3-cycle of the
    // real VTK table, which put extract_surface/extract_skin's quad9 mid-face
    // node at the wrong position for hexahedron27 (never a wrong topology --
    // facet keying uses corners only -- just a wrong quadratic node). Fixed
    // in lockstep with `cell_refine_quad_faces` (cell_subdivision.cpp) and
    // `rtpl_hexahedron_table` (refine_templates.cpp) plus their Python twins,
    // since refine derives node 20+k from the k-th quad-face row; see those
    // files. The format layer (cgns.cpp's HEXA_27 permutation, gmsh.cpp's
    // hexahedron27 perm) was already on this corrected convention, so this
    // makes the codebase self-consistent rather than changing the convention.
    static const std::vector<CellFaceDef> hexahedron27 = {
        {CT::Quad9, 4, 9, {0, 4, 7, 3, 16, 15, 19, 11, 20}},
        {CT::Quad9, 4, 9, {1, 2, 6, 5, 9, 18, 13, 17, 21}},
        {CT::Quad9, 4, 9, {0, 1, 5, 4, 8, 17, 12, 16, 22}},
        {CT::Quad9, 4, 9, {3, 7, 6, 2, 19, 14, 18, 10, 23}},
        {CT::Quad9, 4, 9, {0, 3, 2, 1, 11, 10, 9, 8, 24}},
        {CT::Quad9, 4, 9, {4, 5, 6, 7, 12, 13, 14, 15, 25}},
    };
    static const std::vector<CellFaceDef> wedge = {
        {CT::Triangle, 3, 3, {0, 2, 1}}, {CT::Triangle, 3, 3, {3, 4, 5}},
        {CT::Quad, 4, 4, {0, 1, 4, 3}},  {CT::Quad, 4, 4, {1, 2, 5, 4}},
        {CT::Quad, 4, 4, {2, 0, 3, 5}},
    };
    static const std::vector<CellFaceDef> wedge15 = {
        {CT::Triangle6, 3, 6, {0, 2, 1, 8, 7, 6}},
        {CT::Triangle6, 3, 6, {3, 4, 5, 9, 10, 11}},
        {CT::Quad8, 4, 8, {0, 1, 4, 3, 6, 13, 9, 12}},
        {CT::Quad8, 4, 8, {1, 2, 5, 4, 7, 14, 10, 13}},
        {CT::Quad8, 4, 8, {2, 0, 3, 5, 8, 12, 11, 14}},
    };
    // wedge18 = wedge15 plus a center node on each of the three quad faces:
    // 15=(0,1,4,3), 16=(1,2,5,4), 17=(2,0,3,5) — the quad faces in the order
    // they already appear here. The two triangle faces gain no node.
    static const std::vector<CellFaceDef> wedge18 = {
        {CT::Triangle6, 3, 6, {0, 2, 1, 8, 7, 6}},
        {CT::Triangle6, 3, 6, {3, 4, 5, 9, 10, 11}},
        {CT::Quad9, 4, 9, {0, 1, 4, 3, 6, 13, 9, 12, 15}},
        {CT::Quad9, 4, 9, {1, 2, 5, 4, 7, 14, 10, 13, 16}},
        {CT::Quad9, 4, 9, {2, 0, 3, 5, 8, 12, 11, 14, 17}},
    };
    static const std::vector<CellFaceDef> pyramid = {
        {CT::Quad, 4, 4, {0, 3, 2, 1}},  {CT::Triangle, 3, 3, {0, 1, 4}},
        {CT::Triangle, 3, 3, {1, 2, 4}}, {CT::Triangle, 3, 3, {2, 3, 4}},
        {CT::Triangle, 3, 3, {3, 0, 4}},
    };
    static const std::vector<CellFaceDef> pyramid13 = {
        {CT::Quad8, 4, 8, {0, 3, 2, 1, 8, 7, 6, 5}}, {CT::Triangle6, 3, 6, {0, 1, 4, 5, 10, 9}},
        {CT::Triangle6, 3, 6, {1, 2, 4, 6, 11, 10}}, {CT::Triangle6, 3, 6, {2, 3, 4, 7, 12, 11}},
        {CT::Triangle6, 3, 6, {3, 0, 4, 8, 9, 12}},
    };
    // pyramid14 = pyramid13 plus node 13 at the base-face center.
    static const std::vector<CellFaceDef> pyramid14 = {
        {CT::Quad9, 4, 9, {0, 3, 2, 1, 8, 7, 6, 5, 13}},
        {CT::Triangle6, 3, 6, {0, 1, 4, 5, 10, 9}},
        {CT::Triangle6, 3, 6, {1, 2, 4, 6, 11, 10}},
        {CT::Triangle6, 3, 6, {2, 3, 4, 7, 12, 11}},
        {CT::Triangle6, 3, 6, {3, 0, 4, 8, 9, 12}},
    };
    switch (VolumeType) {
        case CT::Tetra:
            return tetra;
        case CT::Tetra10:
            return tetra10;
        case CT::Hexahedron:
            return hexahedron;
        case CT::Hexahedron20:
            return hexahedron20;
        case CT::Hexahedron27:
            return hexahedron27;
        case CT::Wedge:
            return wedge;
        case CT::Wedge15:
            return wedge15;
        case CT::Wedge18:
            return wedge18;
        case CT::Pyramid:
            return pyramid;
        case CT::Pyramid13:
            return pyramid13;
        case CT::Pyramid14:
            return pyramid14;
        default:
            return empty;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
