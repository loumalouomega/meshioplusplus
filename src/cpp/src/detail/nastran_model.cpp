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
// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Nastran numbers the hex20/wedge15 mid-side nodes bottom, vertical, top;
// meshio++ (VTK) numbers them bottom, top, vertical (the bulk reader's tables).
constexpr int kNmHexa20[20] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                               10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int kNmPenta15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};

constexpr NastranCardSpec kNmCards[] = {
    {"CBAR", "line", 2, nullptr, 0, nullptr},
    {"CBEAM", "line", 2, nullptr, 0, nullptr},
    {"CBUSH", "line", 2, nullptr, 0, nullptr},
    {"CHEXA", "hexahedron", 8, "hexahedron20", 20, kNmHexa20},
    {"CONM2", "vertex", 1, nullptr, 0, nullptr},
    {"CONROD", "line", 2, nullptr, 0, nullptr},
    {"CPENTA", "wedge", 6, "wedge15", 15, kNmPenta15},
    {"CPYRAM", "pyramid", 5, "pyramid13", 13, nullptr},
    {"CQUAD", "quad", 4, "quad9", 9, nullptr},
    {"CQUAD4", "quad", 4, nullptr, 0, nullptr},
    {"CQUAD8", "quad", 4, "quad8", 8, nullptr},
    {"CQUADR", "quad", 4, nullptr, 0, nullptr},
    {"CROD", "line", 2, nullptr, 0, nullptr},
    {"CSHEAR", "quad", 4, nullptr, 0, nullptr},
    {"CTETRA", "tetra", 4, "tetra10", 10, nullptr},
    {"CTRIA3", "triangle", 3, nullptr, 0, nullptr},
    {"CTRIA6", "triangle", 3, "triangle6", 6, nullptr},
    {"CTRIAR", "triangle", 3, nullptr, 0, nullptr},
    {"CTUBE", "line", 2, nullptr, 0, nullptr},
    {"CVISC", "line", 2, nullptr, 0, nullptr},
    {"PLOTEL", "line", 2, nullptr, 0, nullptr},
};

NDArray nm_int_array(const std::vector<std::int64_t>& rV) {
    NDArray a(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

struct NmBlock {
    std::string mCard;
    std::string mType;
    std::size_t mNodes;
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mEid;
    std::vector<std::int64_t> mPid;
};

}  // namespace

const NastranCardSpec* nastran_card_spec(std::string_view Card) {
    for (const NastranCardSpec& spec : kNmCards)
        if (Card == spec.mCard)
            return &spec;
    return nullptr;
}

NastranCells nastran_add_cells(Mesh& rMesh, const std::vector<NastranCardRows>& rCards,
                               const std::unordered_map<std::int64_t, std::size_t>& rGridIndex,
                               const std::unordered_set<std::int64_t>& rScalarPoints,
                               const std::map<std::int64_t, std::string>& rPtype,
                               const std::string& rWho) {
    std::vector<NmBlock> blocks;
    std::size_t dropped = 0;
    for (const NastranCardRows& card : rCards) {
        const NastranCardSpec* spec = nastran_card_spec(card.mCard);
        if (spec == nullptr)
            continue;
        const std::size_t width = card.mWidth;
        const std::size_t n = card.mEid.size();
        NmBlock linear{card.mCard, spec->mLinear, spec->mLinearNodes, {}, {}, {}};
        NmBlock quadratic{
            card.mCard, spec->mQuadratic ? spec->mQuadratic : "", spec->mQuadraticNodes, {}, {},
            {}};
        std::size_t partial = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t* row = card.mNodes.data() + i * width;
            bool quad = false;
            if (spec->mQuadratic != nullptr && width >= spec->mQuadraticNodes) {
                std::size_t given = 0;
                for (std::size_t k = spec->mLinearNodes; k < spec->mQuadraticNodes; ++k)
                    given += row[k] != 0 ? 1 : 0;
                quad = given == spec->mQuadraticNodes - spec->mLinearNodes;
                partial += (given != 0 && !quad) ? 1 : 0;
            }
            NmBlock& b = quad ? quadratic : linear;
            std::vector<std::int64_t> conn(b.mNodes);
            bool ok = true;
            for (std::size_t k = 0; k < b.mNodes && ok; ++k) {
                const std::size_t src = (quad && spec->mPermutation)
                                            ? static_cast<std::size_t>(spec->mPermutation[k])
                                            : k;
                const auto it = rGridIndex.find(row[src]);
                if (it == rGridIndex.end()) {
                    if (rScalarPoints.count(row[src]) == 0)
                        throw ReadError(rWho + ": " + card.mCard + " " +
                                        std::to_string(card.mEid[i]) + " references " +
                                        (row[src] == 0
                                             ? std::string("no node")
                                             : "undefined GRID " + std::to_string(row[src])) +
                                        " as its node " + std::to_string(src + 1));
                    ok = false;  // a scalar point, not a GRID
                } else {
                    conn[k] = static_cast<std::int64_t>(it->second);
                }
            }
            if (!ok) {
                ++dropped;
                continue;
            }
            b.mConn.insert(b.mConn.end(), conn.begin(), conn.end());
            b.mEid.push_back(card.mEid[i]);
            b.mPid.push_back(card.mPid[i]);
        }
        if (partial != 0)
            log::warn("{}: {} {} element(s) have only some mid-side nodes; read as {}", rWho,
                      partial, card.mCard, spec->mLinear);
        if (!linear.mEid.empty())
            blocks.push_back(std::move(linear));
        if (!quadratic.mEid.empty())
            blocks.push_back(std::move(quadratic));
    }
    if (dropped != 0)
        log::warn("{}: skipped {} element(s) that connect scalar points", rWho, dropped);

    NastranCells out;
    {
        std::vector<NDArray> eids;
        std::vector<NDArray> pids;
        std::size_t shared = 0;
        for (NmBlock& b : blocks) {
            out.mOffsets.push_back(out.mNumCells);
            out.mSizes.push_back(b.mEid.size());
            if (b.mCard != "CONM2")
                for (std::size_t i = 0; i < b.mEid.size(); ++i)
                    shared += out.mCellIndex.emplace(b.mEid[i], out.mNumCells + i).second ? 0 : 1;
            out.mNumCells += b.mEid.size();
            NDArray conn(DType::Int64, {b.mEid.size(), b.mNodes});
            std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
            rMesh.AddCellBlock(b.mType, std::move(conn));
            eids.push_back(nm_int_array(b.mEid));
            pids.push_back(nm_int_array(b.mPid));
        }
        if (shared != 0)
            log::warn(
                "{}: {} element id(s) are used by more than one card; their element results "
                "go to the first",
                rWho, shared);
        if (!blocks.empty()) {
            rMesh.AddCellData("nastran:eid", std::move(eids));
            rMesh.AddCellData("nastran:pid", std::move(pids));
        }
    }

    std::map<std::int64_t, std::pair<std::vector<std::int64_t>, int>> by_pid;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const int dim = cell_type_dimension(cell_type_from_name(blocks[b].mType));
        for (std::size_t i = 0; i < blocks[b].mPid.size(); ++i) {
            const std::int64_t p = blocks[b].mPid[i];
            if (p <= 0)
                continue;
            auto& entry = by_pid[p];
            if (entry.first.empty())
                entry.second = dim;
            entry.first.push_back(static_cast<std::int64_t>(out.mOffsets[b] + i));
            entry.second = std::max(entry.second, dim);
        }
    }
    for (auto& [p, entry] : by_pid) {
        const auto it = rPtype.find(p);
        const std::string name =
            (it != rPtype.end() ? it->second : std::string("PID")) + "_" + std::to_string(p);
        rMesh.AddRegion(Region(name, RegionKind::Cell, entry.second, p, nm_int_array(entry.first)));
    }
    return out;
}

void nastran_add_frame(Mesh& rMesh, const std::string& rFrame,
                       const std::vector<std::int64_t>& rValues, const std::string& rWho) {
    const auto nonzero =
        std::count_if(rValues.begin(), rValues.end(), [](std::int64_t c) { return c != 0; });
    if (nonzero == 0)
        return;
    if (rFrame == "CP")
        log::warn(
            "{}: {} GRID(s) have CP != 0; their coordinates are kept in the local system, not "
            "transformed",
            rWho, nonzero);
    else
        log::warn("{}: {} GRID(s) have CD != 0; their results are in the local output system", rWho,
                  nonzero);
    rMesh.AddPointData(std::string("nastran:") + (rFrame == "CP" ? "cp" : "cd"),
                       nm_int_array(rValues));
}

}  // namespace detail
}  // namespace meshioplusplus
