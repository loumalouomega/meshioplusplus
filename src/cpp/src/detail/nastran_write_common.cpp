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
#include <initializer_list>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "nastran_write_common.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

std::string nastran_card_name(std::int64_t Code) {
    static const std::pair<const char*, std::int64_t> kTypes[] = {
        {"CBAR", 24},    {"CBEAM", 54},  {"CBUSH", 26},  {"CDAMP1", 2},  {"CDAMP2", 3},
        {"CELAS1", 6},   {"CELAS2", 7},  {"CGAP", 19},   {"CHEXA", 73},  {"CMASS1", 10},
        {"CMASS2", 11},  {"CONM1", 14},  {"CONM2", 15},  {"CONROD", 16}, {"CPENTA", 41},
        {"CPYRAM", 172}, {"CQUAD", 91},  {"CQUAD4", 51}, {"CQUAD8", 47}, {"CQUADR", 80},
        {"CROD", 30},    {"CSHEAR", 31}, {"CTETRA", 55}, {"CTRIA3", 59}, {"CTRIA6", 48},
        {"CTRIAR", 92},  {"CTUBE", 37},  {"CVISC", 39},  {"PLOTEL", 52},
    };
    for (const auto& [name, number] : kTypes)
        if (number == Code)
            return name;
    return "card " + std::to_string(Code);
}

namespace {

constexpr NastranTypeSpec kNastranWriteTypes[] = {
    {"vertex", "CONM2", 1, -1},    {"line", nullptr, 2, 3},        {"triangle", "CTRIA3", 3, 1},
    {"triangle6", "CTRIA6", 6, 1}, {"quad", "CQUAD4", 4, 1},       {"quad8", "CQUAD8", 8, 1},
    {"quad9", "CQUAD", 9, 1},      {"tetra", "CTETRA", 4, 2},      {"tetra10", "CTETRA", 10, 2},
    {"pyramid", "CPYRAM", 5, 2},   {"pyramid13", "CPYRAM", 13, 2}, {"wedge", "CPENTA", 6, 2},
    {"wedge15", "CPENTA", 15, 2},  {"hexahedron", "CHEXA", 8, 2},  {"hexahedron20", "CHEXA", 20, 2},
};

}  // namespace

const NastranTypeSpec* nastran_type_spec(const std::string& rType) {
    for (const NastranTypeSpec& t : kNastranWriteTypes)
        if (rType == t.mType)
            return &t;
    return nullptr;
}

NastranSignature nastran_signature(const Mesh& rMesh) {
    NastranSignature sig;
    sig.emplace_back("points", rMesh.NumPoints(), 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        sig.emplace_back(cb.Type(), cb.NumCells(), cb.NodesPerCell());
    }
    return sig;
}

std::vector<std::int64_t> nastran_node_ids(const Mesh& rMesh, std::size_t N, const char* pWho) {
    std::vector<std::int64_t> ids(N);
    if (!rMesh.HasPointData("nastran:nid")) {
        for (std::size_t i = 0; i < N; ++i)
            ids[i] = static_cast<std::int64_t>(i + 1);
        return ids;
    }
    const NDArray& a = rMesh.PointData("nastran:nid");
    if (a.Size() != N || is_float_dtype(a.Dtype()))
        nastran_write_fail(pWho, "'nastran:nid' must hold one integer per point");
    const Int64View v(a);
    ids.assign(v.Data(), v.Data() + N);
    std::vector<std::int64_t> sorted = ids;
    std::sort(sorted.begin(), sorted.end());
    if (!sorted.empty() && sorted.front() < 1)
        nastran_write_fail(pWho, "GRID ids must be positive");
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        nastran_write_fail(pWho, "GRID ids are not unique");
    return ids;
}

NastranPlies nastran_plies(const Mesh& rMesh, const char* pWho) {
    NastranPlies out;
    const std::string prefix = "nastran:PCOMP_";
    const std::string suffix = ":PLY";
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.size() <= prefix.size() + suffix.size() || name.rfind(prefix, 0) != 0 ||
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        const std::string digits =
            name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
        if (digits.empty() || digits.size() > 15 ||
            digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        const NDArray& a = rMesh.FieldData(name);
        const std::vector<std::size_t>& shape = a.Shape();
        if (shape.size() != 2 || shape[1] != 4 || shape[0] == 0)
            nastran_write_fail(
                pWho, "'" + name + "' must hold (plies, 4) rows of MID, T, THETA and SOUT");
        const DoubleView view(a);
        auto& plies = out[std::stoll(digits)];
        for (std::size_t r = 0; r < shape[0]; ++r) {
            const double* row = view.Data() + 4 * r;
            plies.push_back({row[0], row[1], row[2], row[3]});
        }
    }
    return out;
}

std::string nastran_join(const std::set<std::string>& rNames) {
    std::string out;
    for (const std::string& n : rNames)
        out += (out.empty() ? "" : ", ") + n;
    return out;
}

void nastran_write_fail(const char* pWho, const std::string& rMessage) {
    throw WriteError(std::string(pWho) + ": " + rMessage);
}

std::int64_t nastran_ptype_pid(const std::string& rName, const std::string& rPtype) {
    if (rName.size() <= rPtype.size() + 1 || rName.compare(0, rPtype.size(), rPtype) != 0 ||
        rName[rPtype.size()] != '_')
        return -1;
    std::int64_t pid = 0;
    for (std::size_t i = rPtype.size() + 1; i < rName.size(); ++i) {
        if (rName[i] < '0' || rName[i] > '9' || pid > 100000000000000000)
            return -1;
        pid = pid * 10 + (rName[i] - '0');
    }
    return pid;
}

std::vector<double> nastran_as_doubles(const NDArray& rA) {
    std::vector<double> out(rA.Size());
    if (rA.Dtype() == DType::Int64) {
        const Int64View view(rA);
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = static_cast<double>(view.Data()[i]);
    } else {
        const DoubleView view(rA);
        std::copy(view.Data(), view.Data() + out.size(), out.begin());
    }
    return out;
}

NastranArrays nastran_named_arrays(const Mesh& rMesh, const std::string& rCard) {
    const std::string kPrefix = "nastran:" + rCard + "_";
    NastranArrays out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind(kPrefix, 0) != 0)
            continue;
        const std::size_t colon = name.find(':', kPrefix.size());
        if (colon == std::string::npos || colon == kPrefix.size())
            continue;
        const std::string digits = name.substr(kPrefix.size(), colon - kPrefix.size());
        const std::string member = name.substr(colon + 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        if (member.empty() ||
            member.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") != std::string::npos)
            continue;
        out[std::stoll(digits)][member] = nastran_as_doubles(rMesh.FieldData(name));
    }
    return out;
}

std::set<std::int64_t> nastran_region_pids(const Mesh& rMesh, const std::string& rPtype) {
    std::set<std::int64_t> out;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const std::int64_t pid = nastran_ptype_pid(rMesh.Region(i).mName, rPtype);
        if (pid >= 0)
            out.insert(pid);
    }
    return out;
}

void nastran_fill_block(const Mesh& rMesh, std::size_t B, const NastranTypeSpec& rSpec,
                        std::size_t NumPoints, std::int64_t& rSequential, const char* pWho,
                        NastranWriteBlock& rBlk) {
    const auto cb = rMesh.Cells(B);
    rBlk.mIndex = B;
    rBlk.mCount = cb.NumCells();
    rBlk.mNodes = cb.NodesPerCell();
    if (rBlk.mCount > 0 && rBlk.mNodes != rSpec.mNodes)
        nastran_write_fail(pWho, "a '" + cb.Type() + "' block needs " +
                                     std::to_string(rSpec.mNodes) + " nodes per cell");
    const Int64View conn(cb.Conn());
    rBlk.mConn.assign(conn.Data(), conn.Data() + rBlk.mCount * rBlk.mNodes);
    for (const std::int64_t node : rBlk.mConn)
        if (node < 0 || static_cast<std::size_t>(node) >= NumPoints)
            nastran_write_fail(
                pWho, "a '" + cb.Type() + "' cell refers to a point the mesh does not have");
    rBlk.mEids.resize(rBlk.mCount);
    rBlk.mPids.resize(rBlk.mCount);
    if (rMesh.HasCellData("nastran:eid")) {
        const NDArray& a = rMesh.CellData("nastran:eid", B);
        if (a.Size() != rBlk.mCount)
            nastran_write_fail(pWho, "'nastran:eid' must have one value per cell");
        const Int64View view(a);
        rBlk.mEids.assign(view.Data(), view.Data() + rBlk.mCount);
    } else {
        for (std::size_t i = 0; i < rBlk.mCount; ++i)
            rBlk.mEids[i] = rSequential + static_cast<std::int64_t>(i);
    }
    rSequential += static_cast<std::int64_t>(rBlk.mCount);
    if (rMesh.HasCellData("nastran:pid")) {
        const NDArray& a = rMesh.CellData("nastran:pid", B);
        if (a.Size() != rBlk.mCount)
            nastran_write_fail(pWho, "'nastran:pid' must have one value per cell");
        const Int64View view(a);
        rBlk.mPids.assign(view.Data(), view.Data() + rBlk.mCount);
    } else {
        rBlk.mPids.assign(rBlk.mCount, rSpec.mPid);
    }
}

void nastran_read_components(const Mesh& rMesh, std::size_t B, const char* pWho,
                             NastranWriteBlock& rBlk) {
    std::vector<std::int64_t> c1(rBlk.mCount, 1);
    std::vector<std::int64_t> c2(rBlk.mCount, 1);
    const auto read_component = [&](const char* pName, std::vector<std::int64_t>& rValues) {
        if (!rMesh.HasCellData(pName))
            return;
        const NDArray& a = rMesh.CellData(pName, B);
        if (a.Size() != rBlk.mCount)
            nastran_write_fail(pWho, std::string("'") + pName + "' must have one value per cell");
        const Int64View view(a);
        rValues.assign(view.Data(), view.Data() + rBlk.mCount);
    };
    read_component("nastran:c1", c1);
    read_component("nastran:c2", c2);
    rBlk.mComps.resize(2 * rBlk.mCount);
    for (std::size_t i = 0; i < rBlk.mCount; ++i) {
        rBlk.mComps[2 * i] = c1[i];
        rBlk.mComps[2 * i + 1] = rBlk.mNodes == 1 ? 0 : c2[i];
    }
}

void nastran_read_field(const Mesh& rMesh, std::size_t B, const char* pName,
                        std::initializer_list<std::size_t> Cards, std::size_t Lo, std::size_t Width,
                        std::size_t BarWidth, const char* pWho, NastranWriteBlock& rBlk) {
    if (!rMesh.HasCellData(pName))
        return;
    const NDArray& a = rMesh.CellData(pName, B);
    if (a.Size() != rBlk.mCount * Width)
        nastran_write_fail(pWho, std::string("'") + pName + "' must have " + std::to_string(Width) +
                                     " value(s) per cell");
    const DoubleView view(a);
    for (std::size_t i = 0; i < rBlk.mCount; ++i)
        if (std::find(Cards.begin(), Cards.end(), rBlk.mCards[i]) != Cards.end())
            for (std::size_t k = 0; k < Width; ++k)
                rBlk.mBar[i * BarWidth + Lo + k] = view.Data()[i * Width + k];
}

void nastran_compact_block(const std::vector<bool>& rKeep, std::size_t Kept, std::size_t BarWidth,
                           NastranWriteBlock& rBlk) {
    const auto compact = [&rKeep](auto& rValues, std::size_t Stride) {
        std::size_t to = 0;
        for (std::size_t i = 0; i < rKeep.size(); ++i)
            if (rKeep[i])
                for (std::size_t s = 0; s < Stride; ++s)
                    rValues[to++] = rValues[i * Stride + s];
        rValues.resize(to);
    };
    compact(rBlk.mConn, rBlk.mNodes);
    compact(rBlk.mEids, 1);
    compact(rBlk.mPids, 1);
    compact(rBlk.mCards, 1);
    compact(rBlk.mComps, 2);
    compact(rBlk.mBar, BarWidth);
    rBlk.mCount = Kept;
}

std::vector<std::vector<NastranWriteRow>> nastran_rows_by_card(
    const std::vector<NastranWriteBlock>& rBlocks, std::size_t NumCards) {
    std::vector<std::vector<NastranWriteRow>> rows(NumCards);
    for (std::size_t bi = 0; bi < rBlocks.size(); ++bi)
        for (std::size_t r = 0; r < rBlocks[bi].mCount; ++r)
            rows[rBlocks[bi].mCards[r]].push_back({rBlocks[bi].mEids[r], bi, r});
    for (auto& card_rows : rows)
        std::stable_sort(
            card_rows.begin(), card_rows.end(),
            [](const NastranWriteRow& rA, const NastranWriteRow& rB) { return rA.mEid < rB.mEid; });
    return rows;
}

}  // namespace detail
}  // namespace meshioplusplus
