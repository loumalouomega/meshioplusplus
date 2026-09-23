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
// The mesh of an Ansys MAPDL model, shared by the .cdb reader (formats/ansysinp.cpp)
// and the .rst reader (formats/ansys_rst.cpp). The element categories and the
// degenerate-shape rules follow the open readers pymapdl-reader and mapdl-archive
// (MIT). Python twin: _build in src/python/meshioplusplus/ansysInp/_ansysInp.py.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/ansys_model.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// How an element routine number's nodes become a cell (pymapdl-reader's
// ETYPE_MAP): a point, a line (quadratic when its third node is set), a shell or
// plane (triangle when K == L), a degenerate brick, a native tetrahedron, or a
// line whose extra nodes are orientation nodes.
}  // namespace

AnsysCategory ansys_category(int Routine) {
    static const std::unordered_map<int, AnsysCategory> kTable = [] {
        std::unordered_map<int, AnsysCategory> m;
        for (int n : {7, 21, 71, 175})
            m[n] = AnsysCategory::Point;
        for (int n :
             {1,   3,   4,   8,   10,  11,  12,  14,  16,  17,  18,  20,  23,  24,  31,  32,  33,
              34,  37,  38,  39,  40,  44,  59,  60,  61,  66,  68,  116, 126, 129, 151, 153, 156,
              161, 169, 171, 172, 176, 177, 178, 180, 189, 208, 209, 250, 251, 280, 288, 289, 290})
            m[n] = AnsysCategory::Line;
        for (int n :
             {2,   13,  22,  25,  28,  29,  35,  41,  42,  43,  51,  53,  54,  55,  57,  63,
              67,  75,  77,  78,  79,  81,  82,  83,  88,  91,  93,  99,  106, 115, 118, 121,
              130, 131, 132, 136, 143, 152, 154, 155, 157, 163, 170, 173, 174, 181, 182, 183,
              212, 213, 218, 219, 222, 223, 230, 233, 238, 252, 281, 282, 283, 292, 293})
            m[n] = AnsysCategory::Shell;
        for (int n : {5,   30,  45,  46,  62,  64,  65,  69,  70,  80,  89,  90,  95,  96,
                      97,  100, 101, 102, 103, 104, 105, 107, 108, 117, 120, 122, 164, 185,
                      186, 190, 192, 215, 220, 226, 231, 236, 239, 272, 273, 278, 279})
            m[n] = AnsysCategory::Brick;
        for (int n : {87, 92, 98, 119, 123, 140, 168, 187, 221, 227, 232, 237, 240, 285, 291})
            m[n] = AnsysCategory::Tet;
        for (int n : {188, 214, 216, 217})
            m[n] = AnsysCategory::LinearLine;
        return m;
    }();
    const auto it = kTable.find(Routine);
    return it == kTable.end() ? AnsysCategory::Skip : it->second;
}

// MESH200's shape comes from KEYOPT(1) (pymapdl-reader's MESH200_MAP).
AnsysCategory ansys_mesh200_category(int Keyopt1) {
    if (Keyopt1 >= 0 && Keyopt1 <= 3)
        return AnsysCategory::Line;
    if (Keyopt1 >= 4 && Keyopt1 <= 7)
        return AnsysCategory::Shell;
    if (Keyopt1 == 8 || Keyopt1 == 9)
        return AnsysCategory::Tet;
    if (Keyopt1 == 10 || Keyopt1 == 11)
        return AnsysCategory::Brick;
    return AnsysCategory::Skip;
}

// A cell resolved from one element row: the meshio++ type and, per node of that
// type, the element's slot (a slot past the row's end is a missing midside).
namespace {

struct AmodShape {
    const char* mType = nullptr;
    std::vector<int> mSlots;
};

// Slot maps of the degenerate forms, from the ANSYS brick numbering (corners
// I..P = 0..7; mid-edges Q R S T = 8..11 on IJ JK KL LI, U V W X = 12..15 on MN
// NO OP PM, Y Z A B = 16..19 on IM JN KO LP) into meshio++'s (VTK) orders.
AmodShape amod_resolve(AnsysCategory Category, const std::vector<std::int64_t>& rNodes) {
    const std::size_t n = rNodes.size();
    const auto at = [&](std::size_t k) { return k < n ? rNodes[k] : 0; };
    const auto slots = [](std::initializer_list<int> l) { return std::vector<int>(l); };
    switch (Category) {
        case AnsysCategory::Point:
            if (n >= 1)
                return {"vertex", slots({0})};
            break;
        case AnsysCategory::Line:
            if (n >= 3 && at(2) > 0)
                return {"line3", slots({0, 1, 2})};
            if (n >= 2)
                return {"line", slots({0, 1})};
            break;
        case AnsysCategory::LinearLine:
            if (n >= 2)
                return {"line", slots({0, 1})};
            break;
        case AnsysCategory::Shell:
            if (n == 3)
                return {"triangle", slots({0, 1, 2})};
            if (n == 6)
                return {"triangle6", slots({0, 1, 2, 3, 4, 5})};
            if (n > 5) {  // 8-node (5 is a quad plus an orientation node); absent
                          // trailing midsides are missing ones
                if (at(2) == at(3))
                    return {"triangle6", slots({0, 1, 2, 4, 5, 7})};
                if (at(6) == at(7))  // e.g. a linear contact face with no midsides
                    return {"quad", slots({0, 1, 2, 3})};
                return {"quad8", slots({0, 1, 2, 3, 4, 5, 6, 7})};
            }
            if (n >= 4) {
                if (at(2) == at(3))
                    return {"triangle", slots({0, 1, 2})};
                return {"quad", slots({0, 1, 2, 3})};
            }
            break;
        case AnsysCategory::Tet:
            if (n > 4)
                return {"tetra10", slots({0, 1, 2, 3, 4, 5, 6, 7, 8, 9})};
            if (n >= 4)
                return {"tetra", slots({0, 1, 2, 3})};
            break;
        case AnsysCategory::Brick: {
            if (n < 8)
                break;
            const bool quad = n > 8;
            if (at(6) != at(7))  // hexahedron
                return quad ? AmodShape{"hexahedron20",
                                        slots({0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                               10, 11, 12, 13, 14, 15, 16, 17, 18, 19})}
                            : AmodShape{"hexahedron", slots({0, 1, 2, 3, 4, 5, 6, 7})};
            if (at(5) != at(6))  // wedge: K == L, O == P
                return quad ? AmodShape{"wedge15",
                                        slots({0, 1, 2, 4, 5, 6, 8, 9, 11, 12, 13, 15, 16, 17, 18})}
                            : AmodShape{"wedge", slots({0, 1, 2, 4, 5, 6})};
            if (at(2) != at(3))  // pyramid: M == N == O == P
                return quad ? AmodShape{"pyramid13",
                                        slots({0, 1, 2, 3, 4, 8, 9, 10, 11, 16, 17, 18, 19})}
                            : AmodShape{"pyramid", slots({0, 1, 2, 3, 4})};
            // tetrahedron: K == L, M == N == O == P
            return quad ? AmodShape{"tetra10", slots({0, 1, 2, 4, 8, 9, 11, 16, 17, 18})}
                        : AmodShape{"tetra", slots({0, 1, 2, 4})};
        }
        case AnsysCategory::Skip:
            break;
    }
    return {};
}

// Corner pairs of each mid-edge slot of a meshio++ type (for missing midside
// nodes, written as node 0), indexed by position in the type's node list.
std::vector<std::pair<int, int>> amod_midside_edges(std::string_view Type) {
    if (Type == "line3")
        return {{0, 1}};
    if (Type == "triangle6")
        return {{0, 1}, {1, 2}, {2, 0}};
    if (Type == "quad8")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    if (Type == "tetra10")
        return {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    if (Type == "pyramid13")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    if (Type == "wedge15")
        return {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    if (Type == "hexahedron20")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    return {};
}

std::size_t amod_num_corners(std::string_view Type) {
    static const std::pair<std::string_view, std::size_t> kCorners[] = {
        {"vertex", 1}, {"line", 2},    {"triangle", 3}, {"quad", 4},
        {"tetra", 4},  {"pyramid", 5}, {"wedge", 6},    {"hexahedron", 8},
    };
    for (const auto& [prefix, c] : kCorners)
        if (Type.substr(0, prefix.size()) == prefix)
            return c;
    return 0;
}

}  // namespace

Mesh ansys_build_mesh(const AnsysModel& rModel, bool Lenient, std::string_view Label,
                      AnsysInfo& rInfo,
                      std::unordered_map<std::int64_t, std::int64_t>& rNodeIndex) {
    std::vector<double> coords = rModel.mCoords;
    const std::string label(Label);
    std::unordered_map<std::int64_t, std::int64_t>& node_index = rNodeIndex;
    node_index.clear();
    for (std::size_t k = 0; k < rModel.mNodeIds.size(); ++k)
        node_index.emplace(rModel.mNodeIds[k], static_cast<std::int64_t>(k));
    // Missing midside nodes (node 0) become new points at their edge midpoints.
    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> midsides;
    std::size_t n_missing = 0;

    struct Block {
        std::string mType;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mRoutine, mSlot, mMat, mReal, mSecnum;
        std::size_t Rows() const { return mSlot.size(); }
    };
    std::vector<Block> blocks;
    std::map<std::string, std::size_t> block_of;
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> element_loc;
    std::map<int, std::size_t> skipped;  // routine -> count

    for (const AnsysElement& e : rModel.mElements) {
        const auto rt = rModel.mRoutine.find(e.mSlot);
        if (rt == rModel.mRoutine.end())
            throw ReadError(label + ": element " + std::to_string(e.mId) + " uses element type " +
                            std::to_string(e.mSlot) + ", which no ET or ETBLOCK defines");
        const int routine = rt->second;
        AnsysCategory category = ansys_category(routine);
        if (routine == 200) {
            const auto ko = rModel.mKeyopt.find(e.mSlot);
            int k1 = 0;
            if (ko != rModel.mKeyopt.end())
                if (const auto it = ko->second.find(1); it != ko->second.end())
                    k1 = it->second;
            category = ansys_mesh200_category(k1);
        }
        AmodShape shape = amod_resolve(category, e.mNodes);
        if (shape.mType) {
            // A corner node 0 (TARGE170's pilot and line shapes ...) has no cell.
            const std::size_t n_corners = amod_num_corners(shape.mType);
            for (std::size_t k = 0; k < n_corners && k < shape.mSlots.size(); ++k) {
                const auto slot = static_cast<std::size_t>(shape.mSlots[k]);
                if (slot >= e.mNodes.size() || e.mNodes[slot] == 0)
                    shape.mType = nullptr;
            }
        }
        if (!shape.mType) {
            if (!Lenient)
                throw ReadError(label + ": element " + std::to_string(e.mId) + " (element type " +
                                std::to_string(routine) + ", " + std::to_string(e.mNodes.size()) +
                                " nodes) has no meshio++ cell type; read with lenient to skip it");
            ++skipped[routine];
            continue;
        }
        const auto [it, fresh] = block_of.emplace(shape.mType, blocks.size());
        if (fresh)
            blocks.push_back(Block{shape.mType, {}, {}, {}, {}, {}, {}});
        Block& b = blocks[it->second];
        const std::size_t corners = amod_num_corners(shape.mType);
        const auto edges = amod_midside_edges(shape.mType);
        std::vector<std::int64_t> row;
        for (int slot : shape.mSlots) {
            // A row cut short omits trailing midside nodes: they read as node 0.
            const auto k = static_cast<std::size_t>(slot);
            const std::int64_t id = k < e.mNodes.size() ? e.mNodes[k] : 0;
            if (id == 0 && row.size() >= corners) {
                row.push_back(-1);  // resolved below, once the corners are known
                continue;
            }
            const auto found = node_index.find(id);
            if (found == node_index.end())
                throw ReadError(label + ": element " + std::to_string(e.mId) +
                                " names undefined node " + std::to_string(id));
            row.push_back(found->second);
        }
        for (std::size_t k = corners; k < row.size(); ++k) {
            if (row[k] >= 0)
                continue;
            const auto [a, c] = edges[k - corners];
            std::int64_t p = row[static_cast<std::size_t>(a)], q = row[static_cast<std::size_t>(c)];
            const auto key = std::make_pair(std::min(p, q), std::max(p, q));
            auto [mit, added] = midsides.emplace(key, static_cast<std::int64_t>(coords.size() / 3));
            if (added) {
                for (std::size_t d = 0; d < 3; ++d)
                    coords.push_back(0.5 * (coords[static_cast<std::size_t>(p) * 3 + d] +
                                            coords[static_cast<std::size_t>(q) * 3 + d]));
                ++n_missing;
            }
            row[k] = mit->second;
        }
        element_loc[e.mId] = {it->second, b.Rows()};
        b.mConn.insert(b.mConn.end(), row.begin(), row.end());
        b.mRoutine.push_back(routine);
        b.mSlot.push_back(e.mSlot);
        b.mMat.push_back(e.mMat);
        b.mReal.push_back(e.mReal);
        b.mSecnum.push_back(e.mSecnum);
    }
    for (const auto& [routine, count] : skipped)
        log::warn("{}: {} element(s) of type {} skipped (no meshio++ cell type)", Label, count,
                  routine);
    if (n_missing)
        log::warn("{}: {} missing midside node(s) placed at their edge midpoints", Label,
                  n_missing);

    Mesh mesh;
    NDArray points(DType::Float64, {coords.size() / 3, 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    const auto column = [](const std::vector<std::int64_t>& rValues) {
        NDArray a(DType::Int64, {rValues.size()});
        std::copy(rValues.begin(), rValues.end(), a.As<std::int64_t>());
        return a;
    };
    std::vector<NDArray> routine_data, slot_data, mat_data, real_data, secnum_data;
    std::vector<std::int64_t> bases;
    std::int64_t base = 0;
    for (const Block& b : blocks) {
        const std::size_t k = b.mConn.size() / std::max<std::size_t>(b.Rows(), 1);
        NDArray conn(DType::Int64, {b.Rows(), k});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
        routine_data.push_back(column(b.mRoutine));
        slot_data.push_back(column(b.mSlot));
        mat_data.push_back(column(b.mMat));
        real_data.push_back(column(b.mReal));
        secnum_data.push_back(column(b.mSecnum));
        bases.push_back(base);
        base += static_cast<std::int64_t>(b.Rows());
    }
    if (!blocks.empty()) {
        mesh.AddCellData("ansys:element", std::move(routine_data));
        mesh.AddCellData("ansys:type", std::move(slot_data));
        mesh.AddCellData("ansys:mat", std::move(mat_data));
        mesh.AddCellData("ansys:real", std::move(real_data));
        mesh.AddCellData("ansys:secnum", std::move(secnum_data));
    }

    // Components: point and cell regions, and the legacy side channel.
    for (const AnsysComponent& c : rModel.mComponents) {
        std::vector<std::int64_t> entries;
        if (c.mNodes) {
            for (std::int64_t id : c.mIds)
                if (const auto it = node_index.find(id); it != node_index.end())
                    entries.push_back(it->second);
            rInfo.mPointSets[c.mName] = entries;
        } else {
            std::vector<std::vector<std::int64_t>> per(blocks.size());
            for (std::int64_t id : c.mIds)
                if (const auto it = element_loc.find(id); it != element_loc.end()) {
                    per[it->second.first].push_back(static_cast<std::int64_t>(it->second.second));
                    entries.push_back(bases[it->second.first] +
                                      static_cast<std::int64_t>(it->second.second));
                }
            rInfo.mCellSets[c.mName] = per;
        }
        mesh.AddRegion(Region(c.mName, c.mNodes ? RegionKind::Point : RegionKind::Cell, -1, -1,
                              column(entries)));
    }
    return mesh;
}

}  // namespace detail
}  // namespace meshioplusplus
