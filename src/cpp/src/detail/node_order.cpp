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
// The permutation tables that used to live in formats/med.cpp, formats/frd.cpp
// and formats/unv.cpp, plus Code_Aster's `.mail`, in one place (roadmap §1.21).
// Python twin: src/python/meshioplusplus/_node_order.py.

// System includes
#include <map>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/detail/node_order.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

enum class NodeOrderDirection { ToMeshio, FromMeshio };

struct NodeOrderSource {
    const char* mFormat;
    const char* mCellType;
    NodeOrderDirection mDirection;
    std::vector<int> mTable;
};

// Each table is written in the direction its source gives it; the registry
// derives the other one.
const std::vector<NodeOrderSource>& node_order_sources() {
    using D = NodeOrderDirection;
    static const std::vector<NodeOrderSource> sources = {
        // MED. Corner parts are MED's reversed reference orientation; the
        // mid-edge parts come from MEDCoupling's INTERP_KERNEL/CellModel.cxx
        // edge tables and put every mid-edge slot at its corners' midpoint.
        // hexahedron27/wedge18 extend hexahedron20/wedge15 with the face and
        // body centres Code_Aster's MED reader implies (lrmtyp.F90 composed
        // with the code_aster tables below, each centre placed geometrically);
        // hexahedron27's is not an involution.
        {"med", "tetra", D::ToMeshio, {0, 1, 3, 2}},
        {"med", "pyramid", D::ToMeshio, {0, 3, 2, 1, 4}},
        {"med", "wedge", D::ToMeshio, {3, 4, 5, 0, 1, 2}},
        {"med", "hexahedron", D::ToMeshio, {4, 5, 6, 7, 0, 1, 2, 3}},
        {"med", "tetra10", D::ToMeshio, {0, 1, 3, 2, 4, 8, 7, 6, 5, 9}},
        {"med", "pyramid13", D::ToMeshio, {0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10}},
        {"med", "wedge15", D::ToMeshio, {3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8, 12, 13, 14}},
        {"med",
         "wedge18",
         D::ToMeshio,
         {3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8, 12, 13, 14, 15, 16, 17}},
        {"med", "hexahedron20", D::ToMeshio, {4,  5,  6, 7, 0,  1,  2,  3,  12, 13,
                                              14, 15, 8, 9, 10, 11, 16, 17, 18, 19}},
        {"med", "hexahedron27", D::ToMeshio, {4,  5,  6,  7,  0,  1,  2,  3,  12, 13, 14, 15, 8, 9,
                                              10, 11, 16, 17, 18, 19, 24, 22, 21, 23, 25, 20, 26}},
        // Code_Aster `.mail`. Linear cells, tetra10 and pyramid13 already use
        // meshio++'s order. The solids list the bottom-ring mid-edges, then the
        // vertical ones, then the top ring, then face and body centres. Derived
        // from Code_Aster's gmsh reader (inigms.F90) composed with the gmsh
        // tables, and checked against its MED reader (lrmtyp.F90): both give
        // the same element up to a symmetry of the reference cell.
        {"code_aster", "wedge15", D::ToMeshio, {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}},
        {"code_aster",
         "wedge18",
         D::ToMeshio,
         {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17}},
        {"code_aster", "hexahedron20", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                     10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"code_aster", "hexahedron27", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,
                                                     9,  10, 11, 16, 17, 18, 19, 12, 13,
                                                     14, 15, 24, 22, 21, 23, 20, 25, 26}},
        // FEBio `.feb`/`.xplt`: hex27's mid-height face centres run y-, x+, y+,
        // x- (FECore/FESolidElementShape.cpp FEHex27); every other FEBio type,
        // hex20, penta15, pyra13 and tet10 included, is in meshio++'s order.
        {"febio", "hexahedron27", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,
                                                9,  10, 11, 12, 13, 14, 15, 16, 17,
                                                18, 19, 23, 21, 20, 22, 24, 25, 26}},
        // CalculiX `.frd`: the he20 and pe15 mid-edge groups and the be3 mid
        // node sit elsewhere than in Abaqus order (confirmed against ccx 2.23
        // output for the same `.inp`).
        {"frd", "hexahedron20", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                              10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"frd", "wedge15", D::ToMeshio, {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}},
        {"frd", "line3", D::ToMeshio, {0, 2, 1}},
        // MSC Patran neutral file: hex20 and wedge15 list the bottom ring, the
        // vertical mid-edges, then the top ring (Patran Reference Manual,
        // Element Library); every other Patran shape is in meshio++'s order.
        {"patran", "hexahedron20", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"patran", "wedge15", D::ToMeshio, {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}},
        // Elmer mesh directory: the vertical mid-edge nodes of the 820/827
        // bricks come before the top ring, and the 827 mid-height face centres
        // run y-, x+, y+, x- (ElmerSolver's elements.def reference coordinates;
        // its own VTU writer applies the same permutation). Every other Elmer
        // code is already in meshio++'s order.
        {"elmer", "hexahedron20", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"elmer", "hexahedron27", D::ToMeshio, {0,  1,  2,  3,  4,  5,  6,  7,  8,
                                                9,  10, 11, 16, 17, 18, 19, 12, 13,
                                                14, 15, 23, 21, 20, 22, 24, 25, 26}},
        // COMSOL `.mphtxt`/`.mphbin`: corners in tensor order (x fastest), then
        // every other node of the element's quadratic lattice in lexicographic
        // (z, y, x) order ("Mesh Element Numbering Conventions", COMSOL API
        // guide). Checked against real COMSOL files and against AWS Palace's
        // COMSOL-to-gmsh tables composed with the gmsh ones.
        {"mphtxt", "quad", D::ToMeshio, {0, 1, 3, 2}},
        {"mphtxt", "hexahedron", D::ToMeshio, {0, 1, 3, 2, 4, 5, 7, 6}},
        {"mphtxt", "pyramid", D::ToMeshio, {0, 1, 3, 2, 4}},
        {"mphtxt", "triangle6", D::ToMeshio, {0, 1, 2, 3, 5, 4}},
        {"mphtxt", "quad9", D::ToMeshio, {0, 1, 3, 2, 4, 7, 8, 5, 6}},
        {"mphtxt", "tetra10", D::ToMeshio, {0, 1, 2, 3, 4, 6, 5, 7, 8, 9}},
        {"mphtxt", "hexahedron27", D::ToMeshio, {0,  1,  3,  2,  4,  5,  7,  6,  8,
                                                 11, 12, 9,  22, 25, 26, 23, 13, 15,
                                                 21, 19, 16, 18, 14, 20, 10, 24, 17}},
        {"mphtxt",
         "wedge18",
         D::ToMeshio,
         {0, 1, 2, 3, 4, 5, 6, 8, 7, 15, 17, 16, 9, 11, 14, 10, 13, 12}},
        {"mphtxt", "pyramid14", D::ToMeshio, {0, 1, 3, 2, 4, 5, 8, 9, 6, 10, 11, 13, 12, 7}},
        // I-DEAS UNV: parabolic elements list their mid-side nodes
        // "sandwiched" between the corners of each ring; the solids list the
        // bottom ring, then the vertical mid-edges, then the top ring (pinned
        // against gmsh's .unv/.msh twins and Salome's driver).
        {"unv", "line3", D::FromMeshio, {0, 2, 1}},
        {"unv", "triangle6", D::FromMeshio, {0, 3, 1, 4, 2, 5}},
        {"unv", "quad8", D::FromMeshio, {0, 4, 1, 5, 2, 6, 3, 7}},
        {"unv", "quad9", D::FromMeshio, {0, 4, 1, 5, 2, 6, 3, 7, 8}},
        {"unv", "tetra10", D::FromMeshio, {0, 4, 1, 5, 2, 6, 7, 8, 9, 3}},
        {"unv", "pyramid13", D::FromMeshio, {0, 5, 1, 6, 2, 7, 3, 8, 9, 10, 11, 12, 4}},
        {"unv", "wedge15", D::FromMeshio, {0, 6, 1, 7, 2, 8, 12, 13, 14, 3, 9, 4, 10, 5, 11}},
        {"unv", "hexahedron20", D::FromMeshio, {0,  8,  1, 9,  2, 10, 3, 11, 16, 17,
                                                18, 19, 4, 12, 5, 13, 6, 14, 7,  15}},
    };
    return sources;
}

std::vector<int> node_order_inverse(const std::vector<int>& rTable) {
    std::vector<int> inv(rTable.size());
    for (std::size_t i = 0; i < rTable.size(); ++i)
        inv[static_cast<std::size_t>(rTable[i])] = static_cast<int>(i);
    return inv;
}

using NodeOrderMap = std::map<std::pair<std::string, std::string>, NodeOrder, std::less<>>;

const NodeOrderMap& node_order_map() {
    static const NodeOrderMap m = [] {
        NodeOrderMap out;
        for (const NodeOrderSource& s : node_order_sources()) {
            NodeOrder order;
            if (s.mDirection == NodeOrderDirection::ToMeshio) {
                order.mToMeshio = s.mTable;
                order.mFromMeshio = node_order_inverse(s.mTable);
            } else {
                order.mFromMeshio = s.mTable;
                order.mToMeshio = node_order_inverse(s.mTable);
            }
            out.emplace(std::make_pair(std::string(s.mFormat), std::string(s.mCellType)),
                        std::move(order));
        }
        return out;
    }();
    return m;
}

}  // namespace

const NodeOrder* node_order(std::string_view Format, std::string_view CellType) {
    const NodeOrderMap& m = node_order_map();
    const auto it = m.find(std::make_pair(std::string(Format), std::string(CellType)));
    return it == m.end() ? nullptr : &it->second;
}

std::vector<std::pair<std::string_view, std::string_view>> node_order_keys() {
    std::vector<std::pair<std::string_view, std::string_view>> keys;
    for (const auto& [key, order] : node_order_map())
        keys.emplace_back(key.first, key.second);
    return keys;
}

}  // namespace detail
}  // namespace meshioplusplus
