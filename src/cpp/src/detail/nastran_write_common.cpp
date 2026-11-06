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
#include <array>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
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

[[noreturn]] void nastran_write_fail(const char* pWho, const std::string& rMessage) {
    throw WriteError(std::string(pWho) + ": " + rMessage);
}

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

}  // namespace detail
}  // namespace meshioplusplus
