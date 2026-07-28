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

// Properties of the selective-refinement subdivision tables, asserted directly
// rather than sampled through a refined mesh.
//
// Three of these are load-bearing rather than merely nice:
//
//  * closure under intersection is what makes "the smallest admissible superset"
//    well defined, and hence what makes the mesh-wide closure's fixed point
//    unique and independent of traversal order -- the determinism argument;
//  * a template referencing only nodes its own mask creates is what lets
//    refine.cpp derive node existence from the split-edge set alone (an edge
//    node iff split, a face centre iff all four of its edges are, a body iff all
//    twelve), which is in turn what makes two cells sharing an entity agree;
//  * the child facets tiling the parent's boundary is the table-level half of
//    conformity: the mesh-level tests then only have to check that neighbours
//    see the same split edges.

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/refine_templates.hpp"

namespace {

namespace d = meshioplusplus::detail;
using meshioplusplus::CellType;

const std::vector<CellType>& all_types() {
    static const std::vector<CellType> types{CellType::Line,  CellType::Triangle,
                                             CellType::Quad,  CellType::Tetra,
                                             CellType::Wedge, CellType::Hexahedron};
    return types;
}

std::size_t num_edges(CellType Type) {
    return d::cell_refine_edges(Type).size();
}
std::size_t num_faces(CellType Type) {
    return d::cell_refine_quad_faces(Type).size();
}
std::size_t num_corners(CellType Type) {
    return static_cast<std::size_t>(meshioplusplus::cell_type_num_nodes(Type));
}

std::vector<std::uint16_t> admissible_masks(CellType Type) {
    std::vector<std::uint16_t> out;
    for (std::uint32_t m = 0; m < (1u << num_edges(Type)); ++m) {
        if (d::refine_mask_template(Type, static_cast<std::uint16_t>(m)).IsAdmissible())
            out.push_back(static_cast<std::uint16_t>(m));
    }
    return out;
}

// The reference element's corner coordinates, in each type's own node order.
std::vector<std::array<double, 3>> reference_corners(CellType Type) {
    switch (Type) {
        case CellType::Line:
            return {{0, 0, 0}, {1, 0, 0}};
        case CellType::Triangle:
            return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        case CellType::Quad:
            return {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
        case CellType::Tetra:
            return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        case CellType::Wedge:
            return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
        case CellType::Hexahedron:
            return {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        default:
            return {};
    }
}

// The whole flat local node space: corners, then edge midpoints, then quad-face
// centres, then the body centre.
std::vector<std::array<double, 3>> local_coordinates(CellType Type) {
    std::vector<std::array<double, 3>> p = reference_corners(Type);
    const std::size_t nc = p.size();
    const auto mean = [&](const std::vector<std::uint8_t>& rIds) {
        std::array<double, 3> c{0, 0, 0};
        for (std::uint8_t i : rIds)
            for (std::size_t k = 0; k < 3; ++k)
                c[k] += p[i][k];
        for (std::size_t k = 0; k < 3; ++k)
            c[k] /= static_cast<double>(rIds.size());
        return c;
    };
    for (const d::CellEdgePair& e : d::cell_refine_edges(Type))
        p.push_back(mean({e[0], e[1]}));
    for (const d::CellQuadFace& f : d::cell_refine_quad_faces(Type))
        p.push_back(mean({f[0], f[1], f[2], f[3]}));
    std::vector<std::uint8_t> corners(nc);
    for (std::size_t i = 0; i < nc; ++i)
        corners[i] = static_cast<std::uint8_t>(i);
    p.push_back(mean(corners));  // body slot: present in the space, used only by the hexahedron
    return p;
}

double signed_area(const std::vector<std::array<double, 3>>& rP,
                   const std::vector<std::uint8_t>& rPoly) {
    double a = 0.0;
    for (std::size_t i = 0; i < rPoly.size(); ++i) {
        const std::array<double, 3>& u = rP[rPoly[i]];
        const std::array<double, 3>& v = rP[rPoly[(i + 1) % rPoly.size()]];
        a += u[0] * v[1] - v[0] * u[1];
    }
    return 0.5 * a;
}

// The divergence-theorem volume of a closed polyhedron whose faces come from
// cell_faces() (outward-wound), fanned from each face's first corner.
double signed_volume(CellType Type, const std::vector<std::array<double, 3>>& rP,
                     const std::vector<std::uint8_t>& rCell) {
    double v = 0.0;
    for (const d::CellFaceDef& f : d::cell_faces(Type)) {
        for (std::size_t i = 1; i + 1 < f.mNumCorners; ++i) {
            const std::array<double, 3>& a = rP[rCell[f.mNodes[0]]];
            const std::array<double, 3>& b = rP[rCell[f.mNodes[i]]];
            const std::array<double, 3>& c = rP[rCell[f.mNodes[i + 1]]];
            v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                  a[2] * (b[0] * c[1] - b[1] * c[0])) /
                 6.0;
        }
    }
    return v;
}

double measure(CellType Type, const std::vector<std::array<double, 3>>& rP,
               const std::vector<std::uint8_t>& rCell) {
    if (meshioplusplus::cell_type_dimension(Type) == 3)
        return signed_volume(Type, rP, rCell);
    if (meshioplusplus::cell_type_dimension(Type) == 2)
        return signed_area(rP, rCell);
    return rP[rCell[1]][0] - rP[rCell[0]][0];  // line, along x on the reference element
}

// Every local node id the children of `rTpl` reference.
std::set<std::uint8_t> referenced_nodes(const d::RefineMaskTemplate& rTpl) {
    std::set<std::uint8_t> used;
    for (const auto* children : {&rTpl.mChildren, &rTpl.mChildrenAlt})
        for (const std::vector<std::uint8_t>& child : *children)
            for (std::uint8_t n : child)
                used.insert(n);
    return used;
}

// The local node ids a given mask actually creates, by refine.cpp's derived
// rules: corners always, an edge node iff split, a face centre iff all four of
// its edges are split, the body iff every edge is.
std::set<std::uint8_t> nodes_created_by(CellType Type, std::uint16_t Mask) {
    const std::size_t nc = num_corners(Type);
    const std::size_t ne = num_edges(Type);
    const std::vector<d::CellEdgePair>& edges = d::cell_refine_edges(Type);
    std::set<std::uint8_t> present;
    for (std::size_t i = 0; i < nc; ++i)
        present.insert(static_cast<std::uint8_t>(i));
    for (std::size_t k = 0; k < ne; ++k)
        if (Mask & (1u << k))
            present.insert(static_cast<std::uint8_t>(nc + k));
    const std::vector<d::CellQuadFace>& faces = d::cell_refine_quad_faces(Type);
    for (std::size_t f = 0; f < faces.size(); ++f) {
        std::uint16_t need = 0;
        for (std::size_t i = 0; i < 4; ++i) {
            const std::uint8_t a = faces[f][i];
            const std::uint8_t b = faces[f][(i + 1) % 4];
            for (std::size_t k = 0; k < ne; ++k)
                if ((edges[k][0] == a && edges[k][1] == b) ||
                    (edges[k][0] == b && edges[k][1] == a))
                    need = static_cast<std::uint16_t>(need | (1u << k));
        }
        if ((Mask & need) == need)
            present.insert(static_cast<std::uint8_t>(nc + ne + f));
    }
    if (Type == CellType::Hexahedron && Mask == d::refine_full_mask(Type))
        present.insert(static_cast<std::uint8_t>(nc + ne + faces.size()));
    return present;
}

// --- the algebra -------------------------------------------------------------

// The property the whole closure rests on. Without it "the smallest admissible
// superset" is not well defined, the promotion is not a closure operator, and
// the mesh-wide fixed point could depend on the order cells are visited in.
TEST(RefineTemplates, AdmissibleMasksAreClosedUnderIntersection) {
    for (CellType type : all_types()) {
        const std::vector<std::uint16_t> masks = admissible_masks(type);
        ASSERT_FALSE(masks.empty()) << meshioplusplus::cell_type_name(type);
        EXPECT_TRUE(d::refine_mask_template(type, d::refine_full_mask(type)).IsAdmissible())
            << meshioplusplus::cell_type_name(type) << ": the full mask must be admissible";
        for (std::uint16_t a : masks) {
            for (std::uint16_t b : masks) {
                const std::uint16_t both = static_cast<std::uint16_t>(a & b);
                EXPECT_TRUE(d::refine_mask_template(type, both).IsAdmissible())
                    << meshioplusplus::cell_type_name(type) << ": " << a << " & " << b << " = "
                    << both << " is not admissible";
            }
        }
    }
}

TEST(RefineTemplates, PromoteIsAClosureOperatorOntoAdmissibleMasks) {
    for (CellType type : all_types()) {
        const std::uint32_t n = 1u << num_edges(type);
        for (std::uint32_t m = 0; m < n; ++m) {
            const std::uint16_t mask = static_cast<std::uint16_t>(m);
            const std::uint16_t p = d::refine_promote_mask(type, mask, false);
            EXPECT_EQ(p & mask, mask) << "promotion must be a superset";
            EXPECT_TRUE(d::refine_mask_template(type, p).IsAdmissible())
                << "promotion must land on an admissible mask";
            EXPECT_EQ(d::refine_promote_mask(type, p, false), p) << "promotion must be idempotent";
            // The least such superset: nothing admissible sits strictly between.
            for (std::uint16_t a : admissible_masks(type)) {
                if ((a & mask) == mask)
                    EXPECT_EQ(a & p, p)
                        << "found an admissible superset smaller than the promotion";
            }
        }
    }
}

TEST(RefineTemplates, PromoteIsMonotone) {
    for (CellType type : all_types()) {
        const std::uint32_t n = 1u << num_edges(type);
        for (std::uint32_t a = 0; a < n; ++a) {
            for (std::uint32_t b = 0; b < n; ++b) {
                if ((a & b) != a)
                    continue;  // only comparable pairs a <= b
                const std::uint16_t pa =
                    d::refine_promote_mask(type, static_cast<std::uint16_t>(a), false);
                const std::uint16_t pb =
                    d::refine_promote_mask(type, static_cast<std::uint16_t>(b), false);
                EXPECT_EQ(pa & pb, pa) << "promotion must preserve inclusion";
            }
        }
    }
}

TEST(RefineTemplates, PropagateGoesStraightToTheFullSplit) {
    for (CellType type : all_types()) {
        const std::uint16_t full = d::refine_full_mask(type);
        EXPECT_EQ(d::refine_promote_mask(type, 0, true), 0u);
        const std::uint32_t n = 1u << num_edges(type);
        for (std::uint32_t m = 1; m < n; ++m)
            EXPECT_EQ(d::refine_promote_mask(type, static_cast<std::uint16_t>(m), true), full);
    }
}

// --- the tables --------------------------------------------------------------

TEST(RefineTemplates, ChildrenArePositivelyOrientedAndConserveMeasure) {
    for (CellType type : all_types()) {
        const std::vector<std::array<double, 3>> p = local_coordinates(type);
        std::vector<std::uint8_t> parent(num_corners(type));
        for (std::size_t i = 0; i < parent.size(); ++i)
            parent[i] = static_cast<std::uint8_t>(i);
        const double whole = measure(type, p, parent);
        ASSERT_GT(whole, 0.0) << meshioplusplus::cell_type_name(type);

        for (std::uint16_t mask : admissible_masks(type)) {
            const d::RefineMaskTemplate& tpl = d::refine_mask_template(type, mask);
            for (const auto* children : {&tpl.mChildren, &tpl.mChildrenAlt}) {
                if (children->empty())
                    continue;
                double sum = 0.0;
                for (const std::vector<std::uint8_t>& child : *children) {
                    EXPECT_EQ(child.size(), num_corners(type))
                        << meshioplusplus::cell_type_name(type) << " mask " << mask
                        << ": children must keep the parent's cell type";
                    const double m = measure(type, p, child);
                    EXPECT_GT(m, 1e-12) << meshioplusplus::cell_type_name(type) << " mask " << mask
                                        << ": a child is inverted or degenerate";
                    sum += m;
                }
                EXPECT_NEAR(sum, whole, 1e-12) << meshioplusplus::cell_type_name(type) << " mask "
                                               << mask << ": children do not tile the parent";
            }
        }
    }
}

// A template must never reference a node its own mask does not create -- that is
// what lets refine.cpp derive node existence from the split-edge set alone,
// which is in turn why two cells sharing an entity cannot disagree about it.
TEST(RefineTemplates, TemplatesOnlyReferenceNodesTheirMaskCreates) {
    for (CellType type : all_types()) {
        for (std::uint16_t mask : admissible_masks(type)) {
            const std::set<std::uint8_t> created = nodes_created_by(type, mask);
            for (std::uint8_t n : referenced_nodes(d::refine_mask_template(type, mask))) {
                EXPECT_TRUE(created.count(n) != 0)
                    << meshioplusplus::cell_type_name(type) << " mask " << mask
                    << " references local node " << static_cast<int>(n)
                    << ", which that mask does not create";
            }
        }
    }
}

// Every facet of every child is either interior to the parent (shared by exactly
// two children) or lies wholly within one of the parent's own facets. This is
// the table-level half of conformity.
TEST(RefineTemplates, ChildFacetsTileTheParentBoundary) {
    for (CellType type : all_types()) {
        if (meshioplusplus::cell_type_dimension(type) < 2)
            continue;
        const std::size_t nc = num_corners(type);
        const std::size_t ne = num_edges(type);
        const std::vector<d::CellEdgePair>& edges = d::cell_refine_edges(type);

        for (std::uint16_t mask : admissible_masks(type)) {
            const d::RefineMaskTemplate& tpl = d::refine_mask_template(type, mask);
            // The local nodes lying on each of the parent's own facets.
            std::vector<std::set<std::uint8_t>> parent_facets;
            if (meshioplusplus::cell_type_dimension(type) == 2) {
                for (std::size_t k = 0; k < ne; ++k) {
                    std::set<std::uint8_t> s{edges[k][0], edges[k][1]};
                    if (mask & (1u << k))
                        s.insert(static_cast<std::uint8_t>(nc + k));
                    parent_facets.push_back(std::move(s));
                }
            } else {
                for (const d::CellFaceDef& f : d::cell_faces(type)) {
                    std::set<std::uint8_t> s;
                    for (std::size_t i = 0; i < f.mNumCorners; ++i)
                        s.insert(f.mNodes[i]);
                    for (std::size_t k = 0; k < ne; ++k) {
                        if ((mask & (1u << k)) && s.count(edges[k][0]) && s.count(edges[k][1]))
                            s.insert(static_cast<std::uint8_t>(nc + k));
                    }
                    // A quad face's own centre, when it has one.
                    for (std::size_t g = 0; g < num_faces(type); ++g) {
                        const d::CellQuadFace& q = d::cell_refine_quad_faces(type)[g];
                        std::set<std::uint8_t> qs{q[0], q[1], q[2], q[3]};
                        bool same = f.mNumCorners == 4;
                        for (std::size_t i = 0; same && i < f.mNumCorners; ++i)
                            same = qs.count(f.mNodes[i]) != 0;
                        if (same)
                            s.insert(static_cast<std::uint8_t>(nc + ne + g));
                    }
                    parent_facets.push_back(std::move(s));
                }
            }

            for (const auto* children : {&tpl.mChildren, &tpl.mChildrenAlt}) {
                if (children->empty())
                    continue;
                std::map<std::set<std::uint8_t>, int> counts;
                for (const std::vector<std::uint8_t>& child : *children) {
                    if (meshioplusplus::cell_type_dimension(type) == 2) {
                        for (std::size_t i = 0; i < child.size(); ++i)
                            ++counts[{child[i], child[(i + 1) % child.size()]}];
                    } else {
                        for (const d::CellFaceDef& f : d::cell_faces(type)) {
                            std::set<std::uint8_t> key;
                            for (std::size_t i = 0; i < f.mNumCorners; ++i)
                                key.insert(child[f.mNodes[i]]);
                            ++counts[key];
                        }
                    }
                }
                for (const auto& [facet, count] : counts) {
                    EXPECT_LE(count, 2) << meshioplusplus::cell_type_name(type) << " mask " << mask
                                        << ": a facet is shared by more than two children";
                    if (count != 1)
                        continue;
                    bool on_boundary = false;
                    for (const std::set<std::uint8_t>& pf : parent_facets) {
                        bool inside = true;
                        for (std::uint8_t n : facet)
                            inside = inside && pf.count(n) != 0;
                        on_boundary = on_boundary || inside;
                    }
                    EXPECT_TRUE(on_boundary)
                        << meshioplusplus::cell_type_name(type) << " mask " << mask
                        << ": an unshared child facet does not lie on the parent's boundary";
                }
            }
        }
    }
}

// A caller sizes its output from the child count before it can evaluate the
// global-id tie-break, so the two variants must agree about it.
TEST(RefineTemplates, BothVariantsHaveTheSameChildCount) {
    for (CellType type : all_types()) {
        for (std::uint16_t mask : admissible_masks(type)) {
            const d::RefineMaskTemplate& tpl = d::refine_mask_template(type, mask);
            if (!tpl.HasVariant())
                continue;
            EXPECT_EQ(tpl.mChildren.size(), tpl.mChildrenAlt.size());
            EXPECT_NE(tpl.mTieA, tpl.mTieB);
            EXPECT_LT(tpl.mTieA, num_corners(type)) << "the tie-break must be over CORNERS";
            EXPECT_LT(tpl.mTieB, num_corners(type)) << "the tie-break must be over CORNERS";
        }
    }
}

// Only a two-adjacent-edge configuration leaves a quadrilateral remnant needing
// a diagonal; anything else must resolve without consulting global ids, or two
// neighbours numbering a shared face differently could disagree.
TEST(RefineTemplates, OnlyTheAmbiguousMasksCarryAVariant) {
    for (CellType type : all_types()) {
        for (std::uint16_t mask : admissible_masks(type)) {
            const bool has_variant = d::refine_mask_template(type, mask).HasVariant();
            if (!has_variant)
                continue;
            int bits = 0;
            for (std::size_t k = 0; k < num_edges(type); ++k)
                bits += (mask >> k) & 1;
            EXPECT_EQ(bits, 2) << meshioplusplus::cell_type_name(type) << " mask " << mask;
            // The two edges must share a corner -- the adjacent case.
            std::vector<std::uint8_t> ends;
            for (std::size_t k = 0; k < num_edges(type); ++k) {
                if (mask & (1u << k)) {
                    ends.push_back(d::cell_refine_edges(type)[k][0]);
                    ends.push_back(d::cell_refine_edges(type)[k][1]);
                }
            }
            std::set<std::uint8_t> distinct(ends.begin(), ends.end());
            EXPECT_EQ(distinct.size(), 3u) << "the two split edges must share exactly one corner";
        }
    }
}

TEST(RefineTemplates, TetrahedronHasTwentySevenAdmissibleMasks) {
    // 1 empty + 6 singles + 15 pairs + 4 face triples + 1 full. The pairs must
    // ALL be admissible -- both adjacent and opposite -- or the set stops being
    // closed under intersection (two face triples meet in a single edge, and two
    // pairs can meet in anything up to a pair).
    EXPECT_EQ(admissible_masks(CellType::Tetra).size(), 27u);
    EXPECT_EQ(admissible_masks(CellType::Triangle).size(), 8u);
    EXPECT_EQ(admissible_masks(CellType::Quad).size(), 4u);
    EXPECT_EQ(admissible_masks(CellType::Line).size(), 2u);
    EXPECT_EQ(admissible_masks(CellType::Wedge).size(), 4u);
    EXPECT_EQ(admissible_masks(CellType::Hexahedron).size(), 8u);
}

TEST(RefineTemplates, ChildCountsMatchTheDocumentedMatrix) {
    struct Case {
        CellType mType;
        std::uint16_t mMask;
        std::size_t mChildren;
    };
    const std::vector<Case> cases{
        {CellType::Line, 0b1, 2},
        {CellType::Triangle, 0b001, 2},
        {CellType::Triangle, 0b011, 3},
        {CellType::Triangle, 0b111, 4},
        {CellType::Quad, 0b0101, 2},
        {CellType::Quad, 0b1010, 2},
        {CellType::Quad, 0b1111, 4},
        {CellType::Tetra, 0b000001, 2},
        {CellType::Tetra, 0b111111, 8},
        {CellType::Wedge, 0b000111111, 4},
        {CellType::Wedge, 0b111000000, 2},
        {CellType::Wedge, 0b111111111, 8},
        {CellType::Hexahedron, 0b000001010101, 2},
        {CellType::Hexahedron, 0b000010101010, 2},
        {CellType::Hexahedron, 0b111100000000, 2},
        {CellType::Hexahedron, 0b000011111111, 4},
        {CellType::Hexahedron, 0b111111111111, 8},
    };
    for (const Case& c : cases) {
        const d::RefineMaskTemplate& tpl = d::refine_mask_template(c.mType, c.mMask);
        ASSERT_TRUE(tpl.IsAdmissible())
            << meshioplusplus::cell_type_name(c.mType) << " mask " << c.mMask;
        EXPECT_EQ(tpl.NumChildren(), c.mChildren)
            << meshioplusplus::cell_type_name(c.mType) << " mask " << c.mMask;
    }
    // One split edge on a quadrilateral has no all-quad subdivision at all, so
    // it promotes to the opposite pair -- two children, not four.
    EXPECT_EQ(d::refine_promote_mask(CellType::Quad, 0b0001, false), 0b0101u);
    EXPECT_EQ(d::refine_promote_mask(CellType::Quad, 0b0011, false), 0b1111u);
    // A hexahedron promotes per parallel edge class, so one split edge costs a
    // sheet rather than the whole block.
    EXPECT_EQ(d::refine_promote_mask(CellType::Hexahedron, 0b1, false), 0b000001010101u);
    // Three tetrahedron edges meeting at a corner are deliberately inadmissible.
    const std::uint16_t star = 0b0001101;  // edges (0,1), (0,2), (0,3)
    EXPECT_EQ(d::refine_promote_mask(CellType::Tetra, star, false), 0b111111u);
}

}  // namespace
