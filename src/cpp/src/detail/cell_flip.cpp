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
// The orientation flip of each fixed-size cell type and the facet map it
// induces. See cell_flip.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "cell_flip.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"

namespace meshioplusplus::detail {

namespace {

// A higher-order line lists its two ends, then its interior nodes from end 0
// to end 1: reversed, the ends swap and the interior runs backwards.
std::vector<int> cf_line(int N) {
    std::vector<int> p{1, 0};
    for (int i = N - 1; i >= 2; --i)
        p.push_back(i);
    return p;
}

std::vector<int> cf_concat(std::vector<int> Head, const std::vector<int>& rTail) {
    Head.insert(Head.end(), rTail.begin(), rTail.end());
    return Head;
}

// The corner positions of each facet of a cell type, sorted, in facet order.
std::vector<std::vector<int>> cf_facet_corners(CellType Type) {
    std::vector<std::vector<int>> out;
    if (cell_type_dimension(Type) == 3) {
        for (const CellFaceDef& r_f : cell_faces(Type))
            out.emplace_back(r_f.mNodes.begin(), r_f.mNodes.begin() + r_f.mNumCorners);
    } else if (cell_type_dimension(Type) == 2) {
        for (const CellEdgeDef& r_e : cell_edges(Type))
            out.emplace_back(r_e.mNodes.begin(), r_e.mNodes.begin() + r_e.mNumCorners);
    }
    for (auto& r_c : out)
        std::sort(r_c.begin(), r_c.end());
    return out;
}

// Facet j of the flipped cell has the corners P[c] of the original for each of
// its corner slots c; the original facet with that corner set maps to j.
std::vector<int> cf_facet_map(CellType Type, const std::vector<int>& rPerm) {
    const std::vector<std::vector<int>> corners = cf_facet_corners(Type);
    std::vector<int> map(corners.size(), -1);
    for (std::size_t j = 0; j < corners.size(); ++j) {
        std::vector<int> image;
        for (const int c : corners[j])
            image.push_back(rPerm[static_cast<std::size_t>(c)]);
        std::sort(image.begin(), image.end());
        for (std::size_t f = 0; f < corners.size(); ++f)
            if (corners[f] == image)
                map[f] = static_cast<int>(j);
    }
    for (const int m : map)
        if (m < 0)
            throw std::logic_error("meshio++: cell_flip: a facet has no image for " +
                                   cell_type_name(Type));
    return map;
}

struct CfEntry {
    std::vector<int> mPerm;
    std::vector<int> mFacets;
};

const std::map<CellType, CfEntry>& cf_table() {
    static const std::map<CellType, CfEntry> table = [] {
        using CT = CellType;
        const std::vector<int> tri6{0, 2, 1, 5, 4, 3};
        const std::vector<int> quad8{0, 3, 2, 1, 7, 6, 5, 4};
        const std::vector<int> pyr13{0, 3, 2, 1, 4, 8, 7, 6, 5, 9, 12, 11, 10};
        const std::vector<int> wedge15{0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13};
        const std::vector<int> hex20{0,  3,  2,  1,  4,  7,  6,  5,  11, 10,
                                     9,  8,  15, 14, 13, 12, 16, 19, 18, 17};
        std::map<CT, std::vector<int>> perms = {
            {CT::Vertex, {0}},
            {CT::Line, cf_line(2)},
            {CT::Line3, cf_line(3)},
            {CT::Line4, cf_line(4)},
            {CT::Line5, cf_line(5)},
            {CT::Line6, cf_line(6)},
            {CT::Line7, cf_line(7)},
            {CT::Line8, cf_line(8)},
            {CT::Line9, cf_line(9)},
            {CT::Line10, cf_line(10)},
            {CT::Line11, cf_line(11)},
            {CT::Triangle, {0, 2, 1}},
            {CT::Triangle6, tri6},
            {CT::Triangle7, cf_concat(tri6, {6})},
            {CT::Quad, {0, 3, 2, 1}},
            {CT::Quad8, quad8},
            {CT::Quad9, cf_concat(quad8, {8})},
            {CT::Tetra, {0, 2, 1, 3}},
            // Mid-edges 4=(0,1) 5=(1,2) 6=(2,0) 7=(0,3) 8=(1,3) 9=(2,3).
            {CT::Tetra10, {0, 2, 1, 3, 6, 5, 4, 7, 9, 8}},
            {CT::Pyramid, {0, 3, 2, 1, 4}},
            {CT::Pyramid13, pyr13},
            {CT::Pyramid14, cf_concat(pyr13, {13})},
            {CT::Wedge, {0, 2, 1, 3, 5, 4}},
            {CT::Wedge15, wedge15},
            // Face centres 15=(0,1,4,3) 16=(1,2,5,4) 17=(2,0,3,5).
            {CT::Wedge18, cf_concat(wedge15, {17, 16, 15})},
            {CT::Hexahedron, {0, 3, 2, 1, 4, 7, 6, 5}},
            {CT::Hexahedron20, hex20},
            // Face centres 20=x-min 21=x-max 22=y-min 23=y-max 24=bottom
            // 25=top; swapping nodes 1 and 3 swaps the x and y faces.
            {CT::Hexahedron27, cf_concat(hex20, {22, 23, 20, 21, 24, 25, 26})},
        };
        std::map<CT, CfEntry> out;
        for (auto& [type, perm] : perms) {
            CfEntry e;
            e.mFacets = cf_facet_map(type, perm);
            e.mPerm = std::move(perm);
            out.emplace(type, std::move(e));
        }
        return out;
    }();
    return table;
}

}  // namespace

const std::vector<int>* cell_flip_permutation(CellType Type) {
    const auto& r_table = cf_table();
    const auto it = r_table.find(Type);
    return it == r_table.end() ? nullptr : &it->second.mPerm;
}

const std::vector<int>* cell_flip_facet_map(CellType Type) {
    const auto& r_table = cf_table();
    const auto it = r_table.find(Type);
    return it == r_table.end() ? nullptr : &it->second.mFacets;
}

}  // namespace meshioplusplus::detail
