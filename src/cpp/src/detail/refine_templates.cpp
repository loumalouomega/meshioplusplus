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

// The per-(cell type, split-edge mask) subdivision templates behind selective
// refinement. See detail/refine_templates.hpp for the contract; this file is
// the tables plus two pieces of machinery that keep them honest:
//
//  * the tetrahedron's 27 admissible masks are GENERATED from six orbit
//    representatives by the 12 even permutations of its corners, rather than
//    transcribed. An even permutation preserves orientation, so every generated
//    child inherits the representative's positive winding;
//  * the promotion table is built by brute force from the definition -- the
//    intersection of every admissible mask containing m -- rather than from a
//    per-type structural shortcut. That intersection is admissible precisely
//    because each type's admissible set is closed under intersection, which
//    tests/cpp/test_refine_templates.cpp asserts directly.
//
// Anonymous-namespace helpers are `rtpl_`-prefixed: the amalgamation
// concatenates every src/cpp/src/*.cpp into one translation unit.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/refine_templates.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using RtplChildren = std::vector<std::vector<std::uint8_t>>;

/// Everything one cell type contributes: the mask-indexed template table and
/// the promotion table derived from it.
struct RtplTypeTable {
    std::uint8_t mNumCorners = 0;
    std::uint8_t mNumEdges = 0;
    std::uint16_t mFullMask = 0;
    std::vector<RefineMaskTemplate> mMasks;  ///< size 1 << mNumEdges
    std::vector<std::uint16_t> mPromote;     ///< size 1 << mNumEdges
};

const RefineMaskTemplate& rtpl_empty_template() {
    static const RefineMaskTemplate empty;
    return empty;
}

/// The index of edge `(a, b)` in a type's `cell_refine_edges` table.
std::uint8_t rtpl_edge_index(const std::vector<CellEdgePair>& rEdges, std::uint8_t a,
                             std::uint8_t b) {
    for (std::size_t k = 0; k < rEdges.size(); ++k) {
        if ((rEdges[k][0] == a && rEdges[k][1] == b) || (rEdges[k][0] == b && rEdges[k][1] == a))
            return static_cast<std::uint8_t>(k);
    }
    return static_cast<std::uint8_t>(rEdges.size());  // unreachable for a well-formed table
}

/// Fill the promotion table from the definition: the intersection of every
/// admissible mask containing `m`. Closure under intersection (asserted by the
/// tests) is what makes that intersection itself admissible, hence the *least*
/// admissible superset.
void rtpl_build_promote(RtplTypeTable& rTable) {
    const std::size_t n = rTable.mMasks.size();
    std::vector<std::uint16_t> admissible;
    for (std::size_t m = 0; m < n; ++m) {
        if (rTable.mMasks[m].IsAdmissible())
            admissible.push_back(static_cast<std::uint16_t>(m));
    }
    rTable.mPromote.assign(n, rTable.mFullMask);
    for (std::size_t m = 0; m < n; ++m) {
        std::uint16_t acc = rTable.mFullMask;
        for (std::uint16_t a : admissible) {
            if ((a & m) == m)
                acc = static_cast<std::uint16_t>(acc & a);
        }
        rTable.mPromote[m] = acc;
    }
}

// --- the 2D and prismatic tables (hand-written) -------------------------------

RtplTypeTable rtpl_line_table() {
    RtplTypeTable t;
    t.mNumCorners = 2;
    t.mNumEdges = 1;
    t.mFullMask = 0b1;
    t.mMasks.resize(2);
    // line3 layout: 2 = mid(0,1).
    t.mMasks[0b0].mChildren = RtplChildren{{0, 1}};
    t.mMasks[0b1].mChildren = RtplChildren{{0, 2}, {2, 1}};
    rtpl_build_promote(t);
    return t;
}

RtplTypeTable rtpl_triangle_table() {
    RtplTypeTable t;
    t.mNumCorners = 3;
    t.mNumEdges = 3;
    t.mFullMask = 0b111;
    t.mMasks.resize(8);
    // triangle6 layout: 3 = m(0,1), 4 = m(1,2), 5 = m(2,0).
    t.mMasks[0b000].mChildren = RtplChildren{{0, 1, 2}};
    // One split edge: cut from its midpoint to the opposite corner.
    t.mMasks[0b001].mChildren = RtplChildren{{0, 3, 2}, {3, 1, 2}};
    t.mMasks[0b010].mChildren = RtplChildren{{0, 1, 4}, {0, 4, 2}};
    t.mMasks[0b100].mChildren = RtplChildren{{0, 1, 5}, {5, 1, 2}};
    // Two split edges: they share a corner `b`, which is cut off as its own
    // triangle; the remnant quadrilateral (a, m_ab, m_bc, c) takes the diagonal
    // starting at whichever of the two surviving corners a, c has the smaller
    // GLOBAL node id. mChildren is "from a"; mTieA = a, mTieB = c.
    {
        RefineMaskTemplate& m = t.mMasks[0b011];  // corner 1 cut; a = 0, c = 2
        m.mChildren = RtplChildren{{3, 1, 4}, {0, 3, 4}, {0, 4, 2}};
        m.mChildrenAlt = RtplChildren{{3, 1, 4}, {0, 3, 2}, {3, 4, 2}};
        m.mTieA = 0;
        m.mTieB = 2;
    }
    {
        RefineMaskTemplate& m = t.mMasks[0b110];  // corner 2 cut; a = 1, c = 0
        m.mChildren = RtplChildren{{4, 2, 5}, {1, 4, 5}, {1, 5, 0}};
        m.mChildrenAlt = RtplChildren{{4, 2, 5}, {1, 4, 0}, {4, 5, 0}};
        m.mTieA = 1;
        m.mTieB = 0;
    }
    {
        RefineMaskTemplate& m = t.mMasks[0b101];  // corner 0 cut; a = 2, c = 1
        m.mChildren = RtplChildren{{5, 0, 3}, {2, 5, 3}, {2, 3, 1}};
        m.mChildrenAlt = RtplChildren{{5, 0, 3}, {2, 5, 1}, {5, 3, 1}};
        m.mTieA = 2;
        m.mTieB = 1;
    }
    // The standard 1-to-4 split; the central child keeps the parent's winding.
    t.mMasks[0b111].mChildren = RtplChildren{{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};
    rtpl_build_promote(t);
    return t;
}

RtplTypeTable rtpl_quad_table() {
    RtplTypeTable t;
    t.mNumCorners = 4;
    t.mNumEdges = 4;
    t.mFullMask = 0b1111;
    t.mMasks.resize(16);
    // quad9 layout: 4..7 = edge mids (edges 01, 12, 23, 30), 8 = face centre.
    // Only the two OPPOSITE pairs are admissible besides none and all: a
    // quadrangulation of an n-gon satisfies 4Q = B + 2I, so the pentagon left by
    // one split edge and the heptagon left by three cannot be filled with
    // quadrilaterals at any number of interior nodes.
    t.mMasks[0b0000].mChildren = RtplChildren{{0, 1, 2, 3}};
    t.mMasks[0b0101].mChildren = RtplChildren{{0, 4, 6, 3}, {4, 1, 2, 6}};
    t.mMasks[0b1010].mChildren = RtplChildren{{0, 1, 5, 7}, {7, 5, 2, 3}};
    t.mMasks[0b1111].mChildren =
        RtplChildren{{0, 4, 8, 7}, {4, 1, 5, 8}, {8, 5, 2, 6}, {7, 8, 6, 3}};
    rtpl_build_promote(t);
    return t;
}

RtplTypeTable rtpl_wedge_table() {
    RtplTypeTable t;
    t.mNumCorners = 6;
    t.mNumEdges = 9;
    t.mFullMask = 0b111111111;
    t.mMasks.resize(512);
    // wedge18 layout: 6..8 bottom-triangle mids, 9..11 top-triangle mids,
    // 12..14 vertical mids, 15..17 quad-face centres. The six triangle edges
    // form one class and the three verticals another.
    constexpr std::uint16_t tri = 0b000111111;
    constexpr std::uint16_t vert = 0b111000000;
    t.mMasks[0].mChildren = RtplChildren{{0, 1, 2, 3, 4, 5}};
    // Triangles split, full height: the 1-to-4 triangle split extruded.
    t.mMasks[tri].mChildren = RtplChildren{
        {0, 6, 8, 3, 9, 11}, {6, 1, 7, 9, 4, 10}, {8, 7, 2, 11, 10, 5}, {6, 7, 8, 9, 10, 11}};
    // Verticals split only: two stacked wedges.
    t.mMasks[vert].mChildren = RtplChildren{{0, 1, 2, 12, 13, 14}, {12, 13, 14, 3, 4, 5}};
    // Everything split: the mid-level triangle's edge midpoints ARE the three
    // quad-face centres, which is why a wedge needs no body node.
    t.mMasks[t.mFullMask].mChildren =
        RtplChildren{{0, 6, 8, 12, 15, 17},   {6, 1, 7, 15, 13, 16},  {8, 7, 2, 17, 16, 14},
                     {6, 7, 8, 15, 16, 17},   {12, 15, 17, 3, 9, 11}, {15, 13, 16, 9, 4, 10},
                     {17, 16, 14, 11, 10, 5}, {15, 16, 17, 9, 10, 11}};
    rtpl_build_promote(t);
    return t;
}

RtplTypeTable rtpl_hexahedron_table() {
    RtplTypeTable t;
    t.mNumCorners = 8;
    t.mNumEdges = 12;
    t.mFullMask = 0xFFF;
    t.mMasks.resize(4096);
    // hexahedron27 layout: 8..19 edge mids (bottom ring, top ring, verticals),
    // 20..25 face centres, 26 body. The face-centre indices are
    // `cell_refine_quad_faces(Hexahedron)` row order, which since v9.9.0 is
    // cell_faces.hpp's own hexahedron27 order (vtkTriQuadraticHexahedron):
    // 20 = x-min (0,4,7,3), 21 = x-max (1,2,6,5), 22 = y-min (0,1,5,4),
    // 23 = y-max (3,7,6,2), 24 = bottom, 25 = top. The absolute indices below
    // were permuted by the same 3-cycle (20->22, 22->23, 23->20) in that
    // release; they never escape `refine` (output blocks are 8-node), so the
    // renumbering is internal apart from the order new face-centre points are
    // appended in. See cell_faces.cpp for the full derivation.
    //
    // The twelve edges fall into three parallel classes; the admissible masks
    // are their eight unions, so a refinement travels through one dual sheet
    // rather than the whole block.
    constexpr std::uint16_t x = 0b000001010101;  // edges 01, 23, 45, 67
    constexpr std::uint16_t y = 0b000010101010;  // edges 12, 30, 56, 74
    constexpr std::uint16_t z = 0b111100000000;  // the four verticals
    t.mMasks[0].mChildren = RtplChildren{{0, 1, 2, 3, 4, 5, 6, 7}};
    t.mMasks[x].mChildren = RtplChildren{{0, 8, 10, 3, 4, 12, 14, 7}, {8, 1, 2, 10, 12, 5, 6, 14}};
    t.mMasks[y].mChildren = RtplChildren{{0, 1, 9, 11, 4, 5, 13, 15}, {11, 9, 2, 3, 15, 13, 6, 7}};
    t.mMasks[z].mChildren =
        RtplChildren{{0, 1, 2, 3, 16, 17, 18, 19}, {16, 17, 18, 19, 4, 5, 6, 7}};
    // Two classes split: the two faces perpendicular to the untouched direction
    // have all four of their edges split and so carry a centre; no body node.
    // x|y -> bottom/top (24, 25); x|z -> the y-perpendicular pair (22, 23);
    // y|z -> the x-perpendicular pair (21, 20).
    t.mMasks[x | y].mChildren = RtplChildren{{0, 8, 24, 11, 4, 12, 25, 15},
                                             {8, 1, 9, 24, 12, 5, 13, 25},
                                             {24, 9, 2, 10, 25, 13, 6, 14},
                                             {11, 24, 10, 3, 15, 25, 14, 7}};
    t.mMasks[x | z].mChildren = RtplChildren{{0, 8, 10, 3, 16, 22, 23, 19},
                                             {8, 1, 2, 10, 22, 17, 18, 23},
                                             {16, 22, 23, 19, 4, 12, 14, 7},
                                             {22, 17, 18, 23, 12, 5, 6, 14}};
    t.mMasks[y | z].mChildren = RtplChildren{{0, 1, 9, 11, 16, 17, 21, 20},
                                             {11, 9, 2, 3, 20, 21, 18, 19},
                                             {16, 17, 21, 20, 4, 5, 13, 15},
                                             {20, 21, 18, 19, 15, 13, 6, 7}};
    // Everything split: rows follow the parent's own parametric (i,j,k) order
    // over the 3x3x3 lattice, so orientation is preserved by construction.
    t.mMasks[t.mFullMask].mChildren =
        RtplChildren{{0, 8, 24, 11, 16, 22, 26, 20},  {8, 1, 9, 24, 22, 17, 21, 26},
                     {11, 24, 10, 3, 20, 26, 23, 19}, {24, 9, 2, 10, 26, 21, 18, 23},
                     {16, 22, 26, 20, 4, 12, 25, 15}, {22, 17, 21, 26, 12, 5, 13, 25},
                     {20, 26, 23, 19, 15, 25, 14, 7}, {26, 21, 18, 23, 25, 13, 6, 14}};
    rtpl_build_promote(t);
    return t;
}

// --- the tetrahedron table (generated) ---------------------------------------

/// The 12 even permutations of `{0,1,2,3}`, in lexicographic order. Even, so
/// each preserves a tetrahedron's orientation -- which is what lets a
/// representative's positively-wound children be relabelled rather than
/// re-derived.
const std::vector<std::array<std::uint8_t, 4>>& rtpl_even_perms() {
    static const std::vector<std::array<std::uint8_t, 4>> perms = [] {
        // `v` is `int`, not `std::uint8_t`, deliberately: GCC 13's
        // -Wstringop-overflow analysis of std::next_permutation inlined over a
        // byte-sized std::array iterator has a known false-positive blowup
        // (exponential warning/offset enumeration, observed OOM-killing an
        // arm64 CI runner) -- widening the working type sidesteps the buggy
        // codegen path without changing the permutations generated.
        std::vector<std::array<std::uint8_t, 4>> out;
        std::array<int, 4> v{0, 1, 2, 3};
        do {
            int inversions = 0;
            for (std::size_t i = 0; i < 4; ++i)
                for (std::size_t j = i + 1; j < 4; ++j)
                    if (v[i] > v[j])
                        ++inversions;
            if (inversions % 2 == 0)
                out.push_back({static_cast<std::uint8_t>(v[0]), static_cast<std::uint8_t>(v[1]),
                               static_cast<std::uint8_t>(v[2]), static_cast<std::uint8_t>(v[3])});
        } while (std::next_permutation(v.begin(), v.end()));
        return out;
    }();
    return perms;
}

/// Relabel one local node id under a corner permutation: a corner maps through
/// the permutation, an edge slot to the slot of the permuted edge.
std::uint8_t rtpl_map_local(const std::array<std::uint8_t, 4>& rSigma,
                            const std::vector<CellEdgePair>& rEdges, std::uint8_t NumCorners,
                            std::uint8_t Local) {
    if (Local < NumCorners)
        return rSigma[Local];
    const std::uint8_t k = static_cast<std::uint8_t>(Local - NumCorners);
    return static_cast<std::uint8_t>(
        NumCorners + rtpl_edge_index(rEdges, rSigma[rEdges[k][0]], rSigma[rEdges[k][1]]));
}

std::uint16_t rtpl_map_mask(const std::array<std::uint8_t, 4>& rSigma,
                            const std::vector<CellEdgePair>& rEdges, std::uint16_t Mask) {
    std::uint16_t out = 0;
    for (std::size_t k = 0; k < rEdges.size(); ++k) {
        if (Mask & (1u << k))
            out = static_cast<std::uint16_t>(
                out | (1u << rtpl_edge_index(rEdges, rSigma[rEdges[k][0]], rSigma[rEdges[k][1]])));
    }
    return out;
}

RtplChildren rtpl_map_children(const std::array<std::uint8_t, 4>& rSigma,
                               const std::vector<CellEdgePair>& rEdges, std::uint8_t NumCorners,
                               const RtplChildren& rChildren) {
    RtplChildren out;
    out.reserve(rChildren.size());
    for (const std::vector<std::uint8_t>& child : rChildren) {
        std::vector<std::uint8_t> row;
        row.reserve(child.size());
        for (std::uint8_t local : child)
            row.push_back(rtpl_map_local(rSigma, rEdges, NumCorners, local));
        out.push_back(std::move(row));
    }
    return out;
}

RtplTypeTable rtpl_tetra_table() {
    RtplTypeTable t;
    t.mNumCorners = 4;
    t.mNumEdges = 6;
    t.mFullMask = 0b111111;
    t.mMasks.resize(64);

    const std::vector<CellEdgePair>& edges = cell_refine_edges(CellType::Tetra);
    // tetra10 layout: 4 = m(0,1), 5 = m(1,2), 6 = m(0,2), 7 = m(0,3),
    // 8 = m(1,3), 9 = m(2,3).
    const auto bit = [&](std::uint8_t a, std::uint8_t b) {
        return static_cast<std::uint16_t>(1u << rtpl_edge_index(edges, a, b));
    };

    // Six orbit representatives. Everything else in the admissible set is one of
    // these relabelled by an even permutation.
    std::vector<std::pair<std::uint16_t, RefineMaskTemplate>> reps;

    {  // no split edge
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{0, 1, 2, 3}};
        reps.emplace_back(0, std::move(m));
    }
    {  // one split edge (0,1): cut the tetrahedron in half through m(0,1)
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{0, 4, 2, 3}, {4, 1, 2, 3}};
        reps.emplace_back(bit(0, 1), std::move(m));
    }
    {  // two ADJACENT split edges (0,1) and (1,2), meeting at corner 1. Corner
        // 1 becomes its own tetrahedron; the remnant quadrilateral on face
        // (0,1,2) is coned to corner 3 across whichever diagonal starts at the
        // smaller-id surviving corner. mChildren is "from corner 0".
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{4, 1, 5, 3}, {0, 4, 5, 3}, {0, 5, 2, 3}};
        m.mChildrenAlt = RtplChildren{{4, 1, 5, 3}, {0, 4, 2, 3}, {4, 5, 2, 3}};
        m.mTieA = 0;
        m.mTieB = 2;
        reps.emplace_back(static_cast<std::uint16_t>(bit(0, 1) | bit(1, 2)), std::move(m));
    }
    {  // two OPPOSITE split edges (0,1) and (2,3): four tetrahedra around the
        // interior segment m(0,1)-m(2,3), one per edge of the equator.
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{4, 9, 2, 0}, {4, 9, 1, 2}, {4, 9, 3, 1}, {4, 9, 0, 3}};
        reps.emplace_back(static_cast<std::uint16_t>(bit(0, 1) | bit(2, 3)), std::move(m));
    }
    {  // three split edges sharing face (0,1,2): that face takes the 1-to-4
        // triangle split and each sub-triangle is coned to corner 3.
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{0, 4, 6, 3}, {4, 1, 5, 3}, {6, 5, 2, 3}, {4, 5, 6, 3}};
        reps.emplace_back(static_cast<std::uint16_t>(bit(0, 1) | bit(1, 2) | bit(0, 2)),
                          std::move(m));
    }
    {  // every edge split: four corner tetrahedra (each a half-scale homothety
        // about its own vertex), then the residual octahedron split along the
        // fixed interior diagonal 4-9 with the remaining ring 6->7->8->5.
        RefineMaskTemplate m;
        m.mChildren = RtplChildren{{0, 4, 6, 7}, {4, 1, 5, 8}, {6, 5, 2, 9}, {7, 8, 9, 3},
                                   {4, 9, 6, 7}, {4, 9, 7, 8}, {4, 9, 8, 5}, {4, 9, 5, 6}};
        reps.emplace_back(t.mFullMask, std::move(m));
    }

    for (const auto& [rep_mask, rep] : reps) {
        for (const std::array<std::uint8_t, 4>& sigma : rtpl_even_perms()) {
            const std::uint16_t mask = rtpl_map_mask(sigma, edges, rep_mask);
            if (t.mMasks[mask].IsAdmissible())
                continue;  // first permutation reaching a mask wins -- deterministic
            RefineMaskTemplate out;
            out.mChildren = rtpl_map_children(sigma, edges, t.mNumCorners, rep.mChildren);
            if (rep.HasVariant()) {
                out.mChildrenAlt = rtpl_map_children(sigma, edges, t.mNumCorners, rep.mChildrenAlt);
                out.mTieA = sigma[rep.mTieA];
                out.mTieB = sigma[rep.mTieB];
            }
            t.mMasks[mask] = std::move(out);
        }
    }

    rtpl_build_promote(t);
    return t;
}

const std::unordered_map<CellType, RtplTypeTable>& rtpl_tables() {
    static const std::unordered_map<CellType, RtplTypeTable> tables = [] {
        std::unordered_map<CellType, RtplTypeTable> t;
        t.emplace(CellType::Line, rtpl_line_table());
        t.emplace(CellType::Triangle, rtpl_triangle_table());
        t.emplace(CellType::Quad, rtpl_quad_table());
        t.emplace(CellType::Tetra, rtpl_tetra_table());
        t.emplace(CellType::Wedge, rtpl_wedge_table());
        t.emplace(CellType::Hexahedron, rtpl_hexahedron_table());
        return t;
    }();
    return tables;
}

const RtplTypeTable* rtpl_find(CellType Type) {
    const auto& tables = rtpl_tables();
    auto it = tables.find(Type);
    return it == tables.end() ? nullptr : &it->second;
}

}  // namespace

bool refine_type_supported(CellType Type) {
    return rtpl_find(Type) != nullptr;
}

std::uint16_t refine_full_mask(CellType Type) {
    const RtplTypeTable* t = rtpl_find(Type);
    return t == nullptr ? std::uint16_t{0} : t->mFullMask;
}

std::uint16_t refine_promote_mask(CellType Type, std::uint16_t Mask, bool Propagate) {
    const RtplTypeTable* t = rtpl_find(Type);
    if (t == nullptr)
        return 0;
    const std::uint16_t mask = static_cast<std::uint16_t>(Mask & t->mFullMask);
    if (Propagate)
        return mask == 0 ? std::uint16_t{0} : t->mFullMask;
    return t->mPromote[mask];
}

const RefineMaskTemplate& refine_mask_template(CellType Type, std::uint16_t Mask) {
    const RtplTypeTable* t = rtpl_find(Type);
    if (t == nullptr || Mask >= t->mMasks.size())
        return rtpl_empty_template();
    return t->mMasks[Mask];
}

}  // namespace detail
}  // namespace meshioplusplus
