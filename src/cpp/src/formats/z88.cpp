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
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/z88.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// Z88 element type -> node count, the meshio++ type it reads as, how many of
// its (leading) nodes that type keeps, and the degrees of freedom per node.
struct Z88Type {
    int mCode;
    int mNodes;
    const char* mCell;
    int mKeep;
    int mDof;
};

const Z88Type* z88_type(std::int64_t Code) {
    static const Z88Type types[] = {
        {1, 8, "hexahedron", 8, 3},   {2, 2, "line", 2, 6},
        {3, 6, "triangle6", 6, 2},    {4, 2, "line", 2, 3},
        {5, 2, "line", 2, 6},         {6, 3, "triangle", 3, 2},
        {7, 8, "quad8", 8, 2},        {8, 8, "quad8", 8, 2},
        {9, 2, "line", 2, 2},         {10, 20, "hexahedron20", 20, 3},
        {11, 12, "quad", 4, 2},       {12, 12, "quad", 4, 2},
        {13, 2, "line", 2, 3},        {14, 6, "triangle6", 6, 2},
        {15, 6, "triangle6", 6, 2},   {16, 10, "tetra10", 10, 3},
        {17, 4, "tetra", 4, 3},       {18, 6, "triangle6", 6, 3},
        {19, 16, "quad", 4, 3},       {20, 8, "quad8", 8, 3},
        {21, 16, "hexahedron", 8, 3}, {22, 12, "wedge", 6, 3},
        {23, 8, "quad8", 8, 6},       {24, 6, "triangle6", 6, 6},
        {25, 2, "line", 2, 6},
    };
    return Code >= 1 && Code <= 25 ? &types[Code - 1] : nullptr;
}

// Types 11/12/19 (cubic and Lagrange quads) and 21/22 (layered shells) keep
// their corners: 11/12 list them first; 19 is a 4x4 lattice row by row (its
// corners 1, 13, 16, 4 run counter-clockwise); 21/22
// list the corners of one layer, then the other.
std::vector<int> z88_corners(int Code) {
    switch (Code) {
        case 19:
            return {0, 12, 15, 3};
        case 21:
            return {0, 1, 2, 3, 8, 9, 10, 11};
        case 22:
            return {0, 1, 2, 6, 7, 8};
        default:
            return {};
    }
}

[[noreturn]] void z88_fail(const std::string& rWhat, std::size_t Line) {
    throw ReadError("Z88: " + rWhat + " (line " + std::to_string(Line) + ")");
}

std::vector<std::string_view> z88_lines(const std::string& rText) {
    std::vector<std::string_view> lines;
    std::size_t pos = 0;
    while (pos < rText.size()) {
        std::size_t eol = rText.find('\n', pos);
        if (eol == std::string::npos)
            eol = rText.size();
        std::string_view line(rText.data() + pos, eol - pos);
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        lines.push_back(line);
        pos = eol + 1;
    }
    return lines;
}

std::vector<std::string_view> z88_tokens(std::string_view Line) {
    std::vector<std::string_view> out;
    std::size_t k = 0;
    while (k < Line.size()) {
        while (k < Line.size() && std::isspace(static_cast<unsigned char>(Line[k])))
            ++k;
        const std::size_t b = k;
        while (k < Line.size() && !std::isspace(static_cast<unsigned char>(Line[k])))
            ++k;
        if (k > b)
            out.push_back(Line.substr(b, k - b));
    }
    return out;
}

bool z88_int(std::string_view Token, std::int64_t& rOut) {
    if (!Token.empty() && Token.front() == '+')
        Token.remove_prefix(1);
    const auto [p, ec] = std::from_chars(Token.data(), Token.data() + Token.size(), rOut);
    return ec == std::errc() && p == Token.data() + Token.size();
}

bool z88_real(std::string_view Token, double& rOut) {
    const std::string t(Token);
    const char* end = nullptr;
    rOut = detail::parse_double(t.c_str(), end);
    return end == t.c_str() + t.size() && !t.empty();
}

// The leading integers of a line (stopping at the first other token).
std::vector<std::int64_t> z88_leading_ints(std::string_view Line) {
    std::vector<std::int64_t> out;
    for (std::string_view t : z88_tokens(Line)) {
        std::int64_t v = 0;
        if (!z88_int(t, v))
            break;
        out.push_back(v);
    }
    return out;
}

std::string z88_read_text(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Z88: cannot open " + rPath);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

std::string z88_lower(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// The file named `rName` (any case) in `rDir`, or an empty path.
std::filesystem::path z88_sibling(const std::filesystem::path& rDir, const std::string& rName) {
    std::error_code ec;
    for (auto it = std::filesystem::directory_iterator(rDir.empty() ? "." : rDir, ec);
         !ec && it != std::filesystem::directory_iterator(); it.increment(ec))
        if (z88_lower(it->path().filename().string()) == rName && it->is_regular_file(ec))
            return it->path();
    return {};
}

struct Z88Element {
    std::int64_t mId;
    const Z88Type* mType;
    std::vector<std::int64_t> mNodes;  // node ids, file order
    std::size_t mLine;
};

// `z88o2.txt`: rows `node u1 u2 [u3 [u4 u5 u6]]`; every other line (the
// language-dependent header) is skipped.
void z88_attach_displacements(Mesh& rMesh, const std::string& rPath,
                              const std::unordered_map<std::int64_t, std::size_t>& rNodeIndex) {
    const std::string text = z88_read_text(rPath);
    std::vector<std::pair<std::size_t, std::vector<double>>> rows;
    std::size_t width = 0;
    for (std::string_view line : z88_lines(text)) {
        const std::vector<std::string_view> t = z88_tokens(line);
        if (t.size() != 3 && t.size() != 4 && t.size() != 7)
            continue;
        std::int64_t id = 0;
        if (!z88_int(t[0], id))
            continue;
        std::vector<double> u(t.size() - 1);
        bool ok = true;
        for (std::size_t k = 1; k < t.size() && ok; ++k)
            ok = z88_real(t[k], u[k - 1]);
        const auto it = rNodeIndex.find(id);
        if (!ok || it == rNodeIndex.end())
            continue;
        width = std::max(width, u.size());
        rows.emplace_back(it->second, std::move(u));
    }
    if (rows.empty()) {
        log::warn("Z88: '{}' holds no displacement rows; ignored", rPath);
        return;
    }
    const std::size_t n = rMesh.NumPoints();
    NDArray u(DType::Float64, {n, width});
    double* d = u.As<double>();
    std::fill(d, d + n * width, std::numeric_limits<double>::quiet_NaN());
    for (const auto& [p, values] : rows)
        std::copy(values.begin(), values.end(), d + p * width);
    rMesh.AddPointData("U", std::move(u));
}

NDArray z88_ids(const std::vector<std::int64_t>& rIds) {
    NDArray a(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

// Accumulates one element's labelled values (sum, count) from z88o3.txt.
using Z88Sums = std::map<std::string, std::pair<double, std::size_t>>;

// Per-block cell data from a per-cell value.
template <class Fill>
std::vector<NDArray> z88_per_block(const std::vector<std::size_t>& rBlockStart, std::size_t Width,
                                   Fill&& rFill) {
    std::vector<NDArray> out;
    for (std::size_t b = 0; b + 1 < rBlockStart.size(); ++b) {
        const std::size_t n = rBlockStart[b + 1] - rBlockStart[b];
        NDArray a(DType::Float64,
                  Width ? std::vector<std::size_t>{n, Width} : std::vector<std::size_t>{n});
        for (std::size_t r = 0; r < n; ++r)
            rFill(rBlockStart[b] + r, a.As<double>() + r * (Width ? Width : 1));
        out.push_back(std::move(a));
    }
    return out;
}

// `z88o3.txt`: per element a header (`element # = N ...`), a line of column
// labels and rows of reals (a truss prints `SIG = v` on its header line). Each
// element keeps the mean of every label over its rows: solids and plane-stress
// elements as the `SIG` tensor, every other label (beams, shafts, trusses,
// tori, plates, shells) as a scalar array of that name; `SIGV` for all.
void z88_attach_stresses(Mesh& rMesh, const std::string& rPath,
                         const std::unordered_map<std::int64_t, std::size_t>& rElementIndex,
                         const std::vector<int>& rCellCode,
                         const std::vector<std::size_t>& rBlockStart) {
    static const std::vector<std::string> solid_labels = {"SIGXX", "SIGYY", "SIGZZ",
                                                          "TAUXY", "TAUYZ", "TAUZX"};
    static const std::vector<std::string> plane_labels = {"SIGXX", "SIGYY", "TAUXY"};
    auto coordinate = [](const std::string& rLabel) {
        return rLabel == "XX" || rLabel == "YY" || rLabel == "ZZ" || rLabel == "RR" ||
               rLabel == "PHI";
    };
    const std::string text = z88_read_text(rPath);
    const std::size_t ncells = rCellCode.size();
    std::vector<Z88Sums> sums(ncells);
    std::int64_t current = -1;
    std::vector<std::string> labels;
    bool have_labels = false;
    for (std::string_view line : z88_lines(text)) {
        const std::size_t hash = line.find('#');
        const std::string lower = z88_lower(std::string(line.substr(0, hash)));
        if (hash != std::string_view::npos && lower.find("element") != std::string::npos) {
            std::string_view rest = line.substr(hash + 1);
            const std::size_t eq = rest.find('=');
            current = -1;
            have_labels = false;
            if (eq != std::string_view::npos) {
                const std::vector<std::int64_t> id = z88_leading_ints(rest.substr(eq + 1));
                if (!id.empty()) {
                    const auto it = rElementIndex.find(id[0]);
                    if (it != rElementIndex.end())
                        current = static_cast<std::int64_t>(it->second);
                }
            }
            const std::size_t at = rest.find("SIG =");
            if (current >= 0 && at != std::string_view::npos) {
                const std::vector<std::string_view> t = z88_tokens(rest.substr(at + 5));
                double v = 0.0;
                if (!t.empty() && z88_real(t[0], v))
                    sums[static_cast<std::size_t>(current)]["SIGXX"] = {v, 1};
            }
            continue;
        }
        if (current < 0)
            continue;
        const std::vector<std::string_view> t = z88_tokens(line);
        if (t.empty())
            continue;
        std::vector<double> v(t.size());
        std::size_t numeric = 0;
        for (std::size_t k = 0; k < t.size(); ++k)
            numeric += z88_real(t[k], v[k]) ? 1 : 0;
        if (numeric == 0) {
            labels.clear();
            for (std::string_view tok : t)
                labels.emplace_back(tok.substr(0, tok.find('(')));
            have_labels = true;
            continue;
        }
        if (numeric != t.size() || !have_labels || t.size() != labels.size())
            continue;
        Z88Sums& s = sums[static_cast<std::size_t>(current)];
        for (std::size_t k = 0; k < labels.size(); ++k) {
            if (coordinate(labels[k]))
                continue;
            auto& acc = s.emplace(labels[k], std::make_pair(0.0, std::size_t(0))).first->second;
            acc.first += v[k];
            ++acc.second;
        }
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto mean = [&](std::size_t c, const std::string& rName) {
        const auto it = sums[c].find(rName);
        return it != sums[c].end() ? it->second.first / static_cast<double>(it->second.second)
                                   : nan;
    };
    std::vector<const std::vector<std::string>*> tensor(ncells, nullptr);
    std::set<std::string> scalars;
    std::set<std::size_t> widths;
    bool any_v = false;
    for (std::size_t c = 0; c < ncells; ++c) {
        const int code = rCellCode[c];
        const bool solid = code == 1 || code == 10 || code == 16 || code == 17;
        const bool plane = code == 3 || code == 7 || code == 11 || code == 14;
        auto all_in = [&](const std::vector<std::string>& rNames) {
            for (const std::string& n : rNames)
                if (!sums[c].count(n))
                    return false;
            return true;
        };
        if (solid && all_in(solid_labels))
            tensor[c] = &solid_labels;
        else if (plane && all_in(plane_labels))
            tensor[c] = &plane_labels;
        if (tensor[c])
            widths.insert(tensor[c]->size());
        for (const auto& [name, acc] : sums[c]) {
            if (name == "SIGV") {
                any_v = true;
                continue;
            }
            if (!tensor[c] ||
                std::find(tensor[c]->begin(), tensor[c]->end(), name) == tensor[c]->end())
                scalars.insert(name);
        }
    }
    if (widths.size() > 1) {
        // Plane and solid elements do not share a file; keep the solids.
        for (auto& t : tensor)
            if (t != &solid_labels)
                t = nullptr;
        widths = {6};
    }
    if (!widths.empty()) {
        const std::size_t w = *widths.begin();
        rMesh.AddCellData("SIG", z88_per_block(rBlockStart, w, [&](std::size_t c, double* pOut) {
                              for (std::size_t k = 0; k < w; ++k)
                                  pOut[k] = tensor[c] ? mean(c, (*tensor[c])[k]) : nan;
                          }));
    }
    if (any_v)
        rMesh.AddCellData("SIGV", z88_per_block(rBlockStart, 0, [&](std::size_t c, double* pOut) {
                              *pOut = mean(c, "SIGV");
                          }));
    for (const std::string& name : scalars) {
        rMesh.AddCellData(name, z88_per_block(rBlockStart, 0, [&](std::size_t c, double* pOut) {
                              const bool in_tensor =
                                  tensor[c] && std::find(tensor[c]->begin(), tensor[c]->end(),
                                                         name) != tensor[c]->end();
                              *pOut = in_tensor ? nan : mean(c, name);
                          }));
    }
}

// `z88o4.txt`: after the per-element blocks, the sums per node (`node F(1) ...
// F(6)`, introduced by "nodal sums" / "aufsummierten").
void z88_attach_forces(Mesh& rMesh, const std::string& rPath,
                       const std::unordered_map<std::int64_t, std::size_t>& rNodeIndex,
                       std::size_t Width) {
    const std::string text = z88_read_text(rPath);
    std::vector<std::pair<std::size_t, std::vector<double>>> rows;
    bool started = false;
    for (std::string_view line : z88_lines(text)) {
        const std::string low = z88_lower(std::string(line));
        if (low.find("nodal sums") != std::string::npos ||
            low.find("aufsummierten") != std::string::npos) {
            started = true;
            continue;
        }
        if (!started)
            continue;
        const std::vector<std::string_view> t = z88_tokens(line);
        if (t.size() != 7)
            continue;
        std::int64_t id = 0;
        std::vector<double> v(6);
        bool ok = z88_int(t[0], id);
        for (std::size_t k = 0; k < 6 && ok; ++k)
            ok = z88_real(t[k + 1], v[k]);
        const auto it = rNodeIndex.find(id);
        if (!ok || it == rNodeIndex.end())
            continue;
        v.resize(Width);
        rows.emplace_back(it->second, std::move(v));
    }
    if (rows.empty())
        return;
    const std::size_t n = rMesh.NumPoints();
    NDArray f(DType::Float64, {n, Width});
    std::fill(f.As<double>(), f.As<double>() + n * Width, std::numeric_limits<double>::quiet_NaN());
    for (const auto& [p, v] : rows)
        std::copy(v.begin(), v.end(), f.As<double>() + p * Width);
    rMesh.AddPointData("F", std::move(f));
}

// A Z88 input file: a count on its first line, then that many rows.
std::vector<std::vector<std::string_view>> z88_count_rows(
    const std::vector<std::string_view>& rLines, const std::string& rPath) {
    std::vector<std::vector<std::string_view>> rows;
    const std::vector<std::int64_t> head =
        rLines.empty() ? std::vector<std::int64_t>{} : z88_leading_ints(rLines[0]);
    if (head.empty() || head[0] < 0) {
        log::warn("Z88: '{}' does not start with a count; ignored", rPath);
        return rows;
    }
    const std::size_t want = static_cast<std::size_t>(head[0]);
    for (std::size_t i = 1; i < rLines.size() && rows.size() < want; ++i) {
        std::vector<std::string_view> t = z88_tokens(rLines[i]);
        if (!t.empty())
            rows.push_back(std::move(t));
    }
    if (rows.size() < want)
        log::warn("Z88: '{}' ends after {} of {} rows", rPath, rows.size(), want);
    return rows;
}

// `z88i2.txt`: `node dof flag value` rows; flag 1 adds a force, flag 2
// prescribes a displacement.
void z88_attach_constraints(Mesh& rMesh, const std::string& rPath,
                            const std::unordered_map<std::int64_t, std::size_t>& rNodeIndex,
                            std::size_t Width) {
    const std::string text = z88_read_text(rPath);
    const auto rows = z88_count_rows(z88_lines(text), rPath);
    if (rows.empty())
        return;
    const std::size_t n = rMesh.NumPoints();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    NDArray u(DType::Float64, {n, Width}), f(DType::Float64, {n, Width});
    std::fill(u.As<double>(), u.As<double>() + n * Width, nan);
    std::fill(f.As<double>(), f.As<double>() + n * Width, nan);
    std::size_t skipped = 0;
    for (const auto& t : rows) {
        std::int64_t node = 0, dof = 0, flag = 0;
        double value = 0.0;
        if (t.size() < 4 || !z88_int(t[0], node) || !z88_int(t[1], dof) || !z88_int(t[2], flag) ||
            !z88_real(t[3], value) || !rNodeIndex.count(node) || dof < 1 ||
            static_cast<std::size_t>(dof) > Width || (flag != 1 && flag != 2)) {
            ++skipped;
            continue;
        }
        const std::size_t at = rNodeIndex.at(node) * Width + static_cast<std::size_t>(dof - 1);
        if (flag == 2) {
            u.As<double>()[at] = value;
        } else {
            double& slot = f.As<double>()[at];
            slot = (std::isnan(slot) ? 0.0 : slot) + value;
        }
    }
    if (skipped)
        log::warn("Z88: {} constraint row(s) of '{}' skipped", skipped, rPath);
    rMesh.AddPointData("z88:bc:u", std::move(u));
    rMesh.AddPointData("z88:bc:f", std::move(f));
}

// Rows `from to ...` over element ids: `rParse(tokens, value)` fills a row's
// values, assigned to every cell in the range; cells no row names stay empty.
template <class Parse>
std::vector<std::vector<double>> z88_ranges(
    const std::vector<std::string_view>& rLines, const std::string& rPath,
    const std::unordered_map<std::int64_t, std::size_t>& rElementIndex, std::size_t NumCells,
    Parse&& rParse) {
    std::vector<std::vector<double>> out(NumCells);
    for (const auto& t : z88_count_rows(rLines, rPath)) {
        std::int64_t from = 0, to = 0;
        std::vector<double> value;
        if (t.size() < 2 || !z88_int(t[0], from) || !z88_int(t[1], to) ||
            !rParse(std::vector<std::string_view>(t.begin() + 2, t.end()), value)) {
            log::warn("Z88: a row of '{}' is malformed; skipped", rPath);
            continue;
        }
        for (std::int64_t id = from; id <= to; ++id) {
            const auto it = rElementIndex.find(id);
            if (it != rElementIndex.end())
                out[it->second] = value;
        }
    }
    return out;
}

// Z88Aurora's `z88sets.txt`: after a count, each set is a header `#KIND
// PURPOSE id count "name"` and then its ids. Element sets become cell regions
// and node sets point regions (tag = set id); surface sets, whose rows do not
// name the structure file's elements, are skipped.
void z88_attach_sets(Mesh& rMesh, const std::string& rPath,
                     const std::unordered_map<std::int64_t, std::size_t>& rNodeIndex,
                     const std::unordered_map<std::int64_t, std::size_t>& rElementIndex,
                     const std::vector<int>& rCellDim) {
    const std::string text = z88_read_text(rPath);
    std::size_t skipped = 0;
    bool open = false;
    std::string kind, name;
    std::int64_t sid = -1;
    std::vector<std::int64_t> ids;
    auto flush = [&]() {
        if (!open)
            return;
        std::set<std::int64_t> members;
        if (kind == "ELEMENTS") {
            int dim = -1;
            for (std::int64_t id : ids) {
                const auto it = rElementIndex.find(id);
                if (it == rElementIndex.end())
                    continue;
                members.insert(static_cast<std::int64_t>(it->second));
                dim = std::max(dim, rCellDim[it->second]);
            }
            rMesh.AddRegion(
                Region(name, RegionKind::Cell, dim, sid,
                       z88_ids(std::vector<std::int64_t>(members.begin(), members.end()))));
        } else if (kind == "NODES") {
            for (std::int64_t id : ids) {
                const auto it = rNodeIndex.find(id);
                if (it != rNodeIndex.end())
                    members.insert(static_cast<std::int64_t>(it->second));
            }
            rMesh.AddRegion(
                Region(name, RegionKind::Point, -1, sid,
                       z88_ids(std::vector<std::int64_t>(members.begin(), members.end()))));
        } else {
            ++skipped;
        }
    };
    for (std::string_view line : z88_lines(text)) {
        const std::vector<std::string_view> t = z88_tokens(line);
        if (t.empty())
            continue;
        if (t[0].front() == '#') {
            flush();
            open = true;
            kind = z88_lower(std::string(t[0].substr(1)));
            for (char& c : kind)
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            sid = -1;
            if (t.size() > 2 && !z88_int(t[2], sid))
                sid = -1;
            const std::size_t q = line.find('"');
            const std::size_t r = line.rfind('"');
            name = q != std::string_view::npos && r > q ? std::string(line.substr(q + 1, r - q - 1))
                                                        : std::string();
            if (name.empty())
                name = "set_" + std::to_string(sid);
            ids.clear();
            continue;
        }
        if (!open)
            continue;
        for (std::string_view tok : t) {
            std::int64_t v = 0;
            if (z88_int(tok, v))
                ids.push_back(v);
        }
    }
    flush();
    if (skipped)
        log::warn("Z88: {} surface set(s) of '{}' skipped", skipped, rPath);
}

// `z88mat.txt` (with the material files it names), `z88elp.txt` and
// `z88int.txt` as `z88:` cell data.
void z88_attach_inputs(Mesh& rMesh, const std::filesystem::path& rDir,
                       const std::unordered_map<std::int64_t, std::size_t>& rElementIndex,
                       const std::vector<std::size_t>& rBlockStart) {
    const std::size_t ncells = rBlockStart.back();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto any = [](const std::vector<std::vector<double>>& rValues) {
        for (const auto& v : rValues)
            if (!v.empty())
                return true;
        return false;
    };
    // One column `Col` of the per-cell values, or the fill where a cell has none.
    auto column = [&](const std::vector<std::vector<double>>& rValues, std::size_t Col,
                      std::size_t Width, double Fill, DType Type) {
        std::vector<NDArray> out;
        for (std::size_t b = 0; b + 1 < rBlockStart.size(); ++b) {
            const std::size_t n = rBlockStart[b + 1] - rBlockStart[b];
            NDArray a(Type,
                      Width ? std::vector<std::size_t>{n, Width} : std::vector<std::size_t>{n});
            const std::size_t w = Width ? Width : 1;
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t k = 0; k < w; ++k) {
                    const auto& v = rValues[rBlockStart[b] + r];
                    const double x = v.empty() ? Fill : v[Col + k];
                    if (Type == DType::Int64)
                        a.As<std::int64_t>()[r * w + k] = static_cast<std::int64_t>(x);
                    else
                        a.As<double>()[r * w + k] = x;
                }
            out.push_back(std::move(a));
        }
        return out;
    };

    const std::filesystem::path mat = z88_sibling(rDir, "z88mat.txt");
    if (!mat.empty()) {
        const std::string text = z88_read_text(mat.string());
        std::size_t k = 0;
        const auto values =
            z88_ranges(z88_lines(text), mat.string(), rElementIndex, ncells,
                       [&](const std::vector<std::string_view>& rT, std::vector<double>& rOut) {
                           ++k;
                           if (rT.empty())
                               return false;
                           const std::string file(rT[0]);
                           const std::string stem = file.substr(0, file.rfind('.'));
                           std::int64_t number = 0;
                           if (!z88_int(stem, number) || number <= 0)
                               number = static_cast<std::int64_t>(k);
                           const std::filesystem::path found = z88_sibling(rDir, z88_lower(file));
                           double e = nan, nu = nan;
                           bool ok = false;
                           if (!found.empty()) {
                               const std::string mtext = z88_read_text(found.string());
                               const std::vector<std::string_view> ml = z88_lines(mtext);
                               const std::vector<std::string_view> t =
                                   ml.empty() ? std::vector<std::string_view>{} : z88_tokens(ml[0]);
                               ok = t.size() >= 2 && z88_real(t[0], e) && z88_real(t[1], nu);
                           }
                           if (!ok) {
                               log::warn("Z88: material file '{}' is missing or malformed", file);
                               e = nu = nan;
                           }
                           rOut = {static_cast<double>(number), e, nu};
                           return true;
                       });
        if (any(values)) {
            rMesh.AddCellData("z88:material", column(values, 0, 0, 0.0, DType::Int64));
            rMesh.AddCellData("z88:E", column(values, 1, 0, nan, DType::Float64));
            rMesh.AddCellData("z88:nu", column(values, 2, 0, nan, DType::Float64));
        }
    }
    const std::filesystem::path elp = z88_sibling(rDir, "z88elp.txt");
    if (!elp.empty()) {
        const std::string text = z88_read_text(elp.string());
        const auto values =
            z88_ranges(z88_lines(text), elp.string(), rElementIndex, ncells,
                       [nan](const std::vector<std::string_view>& rT, std::vector<double>& rOut) {
                           // Fields a row leaves out stay absent (NaN): Z88R leaves them untouched.
                           rOut.assign(12, nan);
                           for (std::size_t k = 0; k < rT.size() && k < 12; ++k)
                               if (!z88_real(rT[k], rOut[k]))
                                   return false;
                           return true;
                       });
        if (any(values))
            rMesh.AddCellData("z88:elp", column(values, 0, 12, nan, DType::Float64));
    }
    const std::filesystem::path integ = z88_sibling(rDir, "z88int.txt");
    if (!integ.empty()) {
        const std::string text = z88_read_text(integ.string());
        const auto values =
            z88_ranges(z88_lines(text), integ.string(), rElementIndex, ncells,
                       [](const std::vector<std::string_view>& rT, std::vector<double>& rOut) {
                           std::int64_t a = 0, b = 0;
                           if (rT.size() < 2 || !z88_int(rT[0], a) || !z88_int(rT[1], b))
                               return false;
                           rOut = {static_cast<double>(a), static_cast<double>(b)};
                           return true;
                       });
        if (any(values))
            rMesh.AddCellData("z88:int", column(values, 0, 2, -1.0, DType::Int64));
    }
}

}  // namespace

bool is_z88_filename(const std::string& rPath) {
    const std::string name = z88_lower(std::filesystem::path(rPath).filename().string());
    return name == "z88i1.txt" || name == "z88structure.txt" || name == "z88o2.txt" ||
           name == "z88o3.txt";
}

Mesh read_z88(const std::string& rPath, bool Results) {
    namespace fs = std::filesystem;
    fs::path structure(rPath);
    const std::string base = z88_lower(structure.filename().string());
    if (base == "z88o2.txt" || base == "z88o3.txt") {
        structure = z88_sibling(structure.parent_path(), "z88i1.txt");
        if (structure.empty())
            structure = z88_sibling(fs::path(rPath).parent_path(), "z88structure.txt");
        if (structure.empty())
            throw ReadError("Z88: no z88i1.txt next to " + rPath);
    }
    const std::string text = z88_read_text(structure.string());
    const std::vector<std::string_view> lines = z88_lines(text);
    std::size_t i = 0;
    while (i < lines.size() && z88_tokens(lines[i]).empty())
        ++i;
    if (i >= lines.size())
        throw ReadError("Z88: " + rPath + " is empty");
    const std::vector<std::int64_t> head = z88_leading_ints(lines[i]);
    if (head.size() < 3 || head[0] < 1 || head[0] > 3 || head[1] < 0 || head[2] < 0)
        z88_fail("the header needs the dimension, node and element counts", i + 1);
    const int ndim = static_cast<int>(head[0]);
    const std::size_t nnodes = static_cast<std::size_t>(head[1]);
    const std::size_t nelem = static_cast<std::size_t>(head[2]);
    // Z88OS v15: ndim nnodes nelem ndof kflag. Z88 <= V13 / Aurora V1: ndim nnodes
    // nelem ndof nmat kflag ibflag ipflag [iqflag] [ihflag].
    const bool legacy = head.size() >= 8;
    const std::int64_t kflag = legacy ? head[5] : (head.size() >= 5 ? head[4] : 0);
    ++i;

    std::vector<double> coords;
    coords.reserve(3 * nnodes);
    std::unordered_map<std::int64_t, std::size_t> node_index;
    node_index.reserve(nnodes);
    std::int64_t node_dof = 0;
    for (std::size_t n = 0; n < nnodes; ++n, ++i) {
        while (i < lines.size() && z88_tokens(lines[i]).empty())
            ++i;
        if (i >= lines.size())
            z88_fail("the file ends before node " + std::to_string(n + 1) + " of " +
                         std::to_string(nnodes),
                     i);
        const std::vector<std::string_view> t = z88_tokens(lines[i]);
        std::int64_t id = 0;
        if (t.size() < static_cast<std::size_t>(2 + ndim) || !z88_int(t[0], id))
            z88_fail("a node line needs its id, degrees of freedom and " + std::to_string(ndim) +
                         " coordinates",
                     i + 1);
        double x[3] = {0.0, 0.0, 0.0};
        for (int d = 0; d < ndim; ++d)
            if (!z88_real(t[2 + static_cast<std::size_t>(d)], x[d]))
                z88_fail("bad coordinate '" + std::string(t[2 + static_cast<std::size_t>(d)]) + "'",
                         i + 1);
        if (kflag == 1) {
            const double r = x[0], phi = x[1] * 3.14159265358979323846 / 180.0;
            x[0] = r * std::cos(phi);
            x[1] = r * std::sin(phi);
        }
        if (!node_index.emplace(id, n).second)
            z88_fail("node " + std::to_string(id) + " is defined twice", i + 1);
        coords.insert(coords.end(), x, x + 3);
        std::int64_t dof = 0;
        if (z88_int(t[1], dof))
            node_dof = std::max(node_dof, dof);
    }
    if (kflag == 1)
        log::warn("Z88: cylindrical input (KFLAG = 1) converted to Cartesian coordinates");

    std::vector<Z88Element> elements;
    elements.reserve(nelem);
    for (std::size_t e = 0; e < nelem; ++e) {
        while (i < lines.size() && z88_tokens(lines[i]).empty())
            ++i;
        if (i >= lines.size())
            z88_fail("the file ends before element " + std::to_string(e + 1) + " of " +
                         std::to_string(nelem),
                     i);
        const std::vector<std::int64_t> h = z88_leading_ints(lines[i]);
        if (h.size() < 2)
            z88_fail("an element starts with its id and type", i + 1);
        const Z88Type* type = z88_type(h[1]);
        if (!type)
            z88_fail(
                "element " + std::to_string(h[0]) + " has unknown type " + std::to_string(h[1]),
                i + 1);
        Z88Element el{h[0], type, {}, i + 1};
        ++i;
        while (el.mNodes.size() < static_cast<std::size_t>(type->mNodes)) {
            if (i >= lines.size())
                z88_fail("element " + std::to_string(el.mId) + " is cut short", el.mLine);
            for (std::int64_t v : z88_leading_ints(lines[i]))
                if (el.mNodes.size() < static_cast<std::size_t>(type->mNodes))
                    el.mNodes.push_back(v);
            ++i;
        }
        elements.push_back(std::move(el));
    }

    Mesh mesh;
    NDArray points(DType::Float64, {nnodes, 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // One block per meshio++ type, in order of first appearance.
    std::vector<std::string> block_types;
    std::map<std::string, std::vector<std::size_t>> by_type;
    std::set<int> reduced;
    for (std::size_t e = 0; e < elements.size(); ++e) {
        const Z88Type* t = elements[e].mType;
        auto [it, fresh] = by_type.emplace(t->mCell, std::vector<std::size_t>{});
        if (fresh)
            block_types.push_back(t->mCell);
        it->second.push_back(e);
        if (t->mKeep < t->mNodes)
            reduced.insert(t->mCode);
    }
    for (int code : reduced) {
        log::warn("Z88: type {} elements keep only their corner nodes", code);
        detail::provenance_note("high-order-dropped", "Z88 type " + std::to_string(code) +
                                                          " elements keep only their corners");
    }
    std::unordered_map<std::int64_t, std::size_t> element_index;
    std::vector<int> cell_code;
    std::vector<std::size_t> block_start{0};
    std::vector<NDArray> type_blocks;
    for (const std::string& ctype : block_types) {
        const std::vector<std::size_t>& members = by_type[ctype];
        const std::size_t k =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(ctype)));
        const detail::NodeOrder* order = detail::node_order("z88", ctype);
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray codes(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const Z88Element& el = elements[members[r]];
            const std::vector<int> corners = z88_corners(el.mType->mCode);
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t src = order ? static_cast<std::size_t>(order->mToMeshio[j]) : j;
                if (!corners.empty())
                    src = static_cast<std::size_t>(
                        corners[order ? static_cast<std::size_t>(order->mToMeshio[j]) : j]);
                const auto it = node_index.find(el.mNodes[src]);
                if (it == node_index.end())
                    z88_fail("element " + std::to_string(el.mId) + " names undefined node " +
                                 std::to_string(el.mNodes[src]),
                             el.mLine);
                c[r * k + j] = static_cast<std::int64_t>(it->second);
            }
            codes.As<std::int64_t>()[r] = el.mType->mCode;
            if (!element_index.emplace(el.mId, cell_code.size()).second)
                z88_fail("element " + std::to_string(el.mId) + " is defined twice", el.mLine);
            cell_code.push_back(el.mType->mCode);
        }
        mesh.AddCellBlock(ctype, std::move(conn));
        type_blocks.push_back(std::move(codes));
        block_start.push_back(cell_code.size());
    }
    if (!type_blocks.empty())
        mesh.AddCellData("z88:type", std::move(type_blocks));

    const std::size_t width = (node_dof == 2 || node_dof == 3 || node_dof == 6)
                                  ? static_cast<std::size_t>(node_dof)
                                  : (ndim == 2 ? 2u : 3u);
    const fs::path dir = structure.parent_path();
    const fs::path i2 = z88_sibling(dir, "z88i2.txt");
    if (!i2.empty())
        z88_attach_constraints(mesh, i2.string(), node_index, width);
    if (!cell_code.empty())
        z88_attach_inputs(mesh, dir, element_index, block_start);
    const fs::path sets = z88_sibling(dir, "z88sets.txt");
    if (!sets.empty()) {
        std::vector<int> cell_dim;
        for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b)
            cell_dim.insert(
                cell_dim.end(), mesh.Cells(b).NumCells(),
                cell_type_dimension(cell_type_from_name(std::string(mesh.Cells(b).Type()))));
        z88_attach_sets(mesh, sets.string(), node_index, element_index, cell_dim);
    }
    if (Results) {
        const fs::path o2 = z88_sibling(dir, "z88o2.txt");
        if (!o2.empty())
            z88_attach_displacements(mesh, o2.string(), node_index);
        const fs::path o4 = z88_sibling(dir, "z88o4.txt");
        if (!o4.empty())
            z88_attach_forces(mesh, o4.string(), node_index, width);
        const fs::path o3 = z88_sibling(dir, "z88o3.txt");
        if (!o3.empty() && !cell_code.empty())
            z88_attach_stresses(mesh, o3.string(), element_index, cell_code, block_start);
    }
    return mesh;
}

// ===========================================================================
// Writer (the Z88OS v15 layout)
// ===========================================================================

namespace {

// Whether a Z88 type only exists in a 3-D structure file.
bool z88_needs_3d(int Code) {
    switch (Code) {
        case 1:
        case 2:
        case 4:
        case 5:
        case 10:
        case 16:
        case 17:
        case 21:
        case 22:
        case 23:
        case 24:
        case 25:
            return true;
        default:
            return false;
    }
}

// The Z88 type a cell of `rCell` is written as, given the file dimension.
int z88_default_code(std::string_view Cell, int Dim) {
    if (Cell == "hexahedron")
        return 1;
    if (Cell == "hexahedron20")
        return 10;
    if (Cell == "tetra")
        return 17;
    if (Cell == "tetra10")
        return 16;
    if (Cell == "triangle6")
        return Dim == 2 ? 14 : 24;
    if (Cell == "quad8")
        return Dim == 2 ? 7 : 23;
    if (Cell == "line")
        return Dim == 2 ? 9 : 4;
    return 0;
}

// The input files the `z88:` arrays describe: `z88i2.txt` from the
// constraints, `z88mat.txt` and one `<n>.txt` per material, `z88elp.txt`,
// `z88int.txt`; element ranges over the written element ids.
std::map<std::string, std::string> z88_deck_files(const Mesh& rMesh,
                                                  const std::vector<std::vector<int>>& rCodes,
                                                  const std::vector<int>& rDof) {
    std::map<std::string, std::string> files;
    char buf[512];
    const std::size_t npts = rMesh.NumPoints();
    const NDArray* u = rMesh.HasPointData("z88:bc:u") ? &rMesh.PointData("z88:bc:u") : nullptr;
    const NDArray* f = rMesh.HasPointData("z88:bc:f") ? &rMesh.PointData("z88:bc:f") : nullptr;
    if (u || f) {
        std::string rows;
        std::size_t count = 0;
        const std::pair<int, const NDArray*> kinds[2] = {{2, u}, {1, f}};
        for (std::size_t p = 0; p < npts; ++p)
            for (int d = 0; d < rDof[p]; ++d)
                for (const auto& [flag, pA] : kinds) {
                    if (!pA || !npts)
                        continue;
                    const std::size_t w = pA->Size() / npts;
                    if (static_cast<std::size_t>(d) >= w)
                        continue;
                    const double v = detail::read_double(*pA, p * w + static_cast<std::size_t>(d));
                    if (std::isnan(v))
                        continue;
                    detail::snprintf_c(buf, sizeof(buf), "%9zu %2d %2d %+.16E\n", p + 1, d + 1,
                                       flag, v);
                    rows += buf;
                    ++count;
                }
        files["z88i2.txt"] = std::to_string(count) + "\n" + rows;
    }

    // An array's rows of the written elements, in element order.
    auto written = [&](const std::string& rName) {
        std::vector<std::vector<double>> out;
        for (std::size_t b = 0; b < rCodes.size(); ++b) {
            const NDArray& a = rMesh.CellData(rName, b);
            const std::size_t n = rCodes[b].size();
            const std::size_t w = n ? a.Size() / n : 0;
            for (std::size_t r = 0; r < n; ++r) {
                if (!rCodes[b][r])
                    continue;
                std::vector<double> row(w);
                for (std::size_t k = 0; k < w; ++k)
                    row[k] = detail::read_double(a, r * w + k);
                out.push_back(std::move(row));
            }
        }
        return out;
    };
    // Runs of equal values as (first id, last id, value); empty values skipped.
    struct Run {
        std::size_t mFirst, mLast;
        std::vector<double> mValue;
    };
    auto ranges = [](const std::vector<std::vector<double>>& rValues) {
        std::vector<Run> runs;
        for (std::size_t k = 0; k < rValues.size(); ++k) {
            if (!runs.empty() && runs.back().mValue == rValues[k] && runs.back().mLast == k)
                runs.back().mLast = k + 1;
            else
                runs.push_back({k + 1, k + 1, rValues[k]});
        }
        std::vector<Run> out;
        for (Run& r : runs)
            if (!r.mValue.empty())
                out.push_back(std::move(r));
        return out;
    };

    if (rMesh.HasCellData("z88:E") && rMesh.HasCellData("z88:nu")) {
        const auto e = written("z88:E");
        const auto nu = written("z88:nu");
        const bool has_want = rMesh.HasCellData("z88:material");
        const auto want = has_want ? written("z88:material") : std::vector<std::vector<double>>{};
        std::vector<std::pair<std::pair<double, double>, std::int64_t>> number;  // insertion order
        auto find = [&](std::pair<double, double> Key) -> std::int64_t {
            for (const auto& [k, n] : number)
                if (k == Key)
                    return n;
            return 0;
        };
        for (std::size_t k = 0; k < e.size(); ++k) {
            const std::pair<double, double> key{e[k][0], nu[k][0]};
            if (std::isnan(key.first) || std::isnan(key.second) || find(key))
                continue;
            std::int64_t w = has_want ? static_cast<std::int64_t>(want[k][0]) : 0;
            std::int64_t used_max = 0;
            bool used = false;
            for (const auto& [kk, n] : number) {
                used_max = std::max(used_max, n);
                used = used || n == w;
            }
            if (w <= 0 || used)
                w = used_max + 1;
            number.emplace_back(key, w);
        }
        std::vector<std::vector<double>> values(e.size());
        for (std::size_t k = 0; k < e.size(); ++k)
            if (!std::isnan(e[k][0]) && !std::isnan(nu[k][0]))
                values[k] = {static_cast<double>(find({e[k][0], nu[k][0]}))};
        const std::vector<Run> runs = ranges(values);
        if (!runs.empty()) {
            std::string body = std::to_string(runs.size()) + "\n";
            for (const Run& r : runs) {
                detail::snprintf_c(buf, sizeof(buf), "%9zu %9zu %lld.txt\n", r.mFirst, r.mLast,
                                   static_cast<long long>(r.mValue[0]));
                body += buf;
            }
            files["z88mat.txt"] = body;
            for (const auto& [key, n] : number) {
                detail::snprintf_c(buf, sizeof(buf), "%+.16E %+.16E\n", key.first, key.second);
                files[std::to_string(n) + ".txt"] = buf;
            }
        }
    }
    if (rMesh.HasCellData("z88:elp")) {
        // The leading fields up to the first NaN, at least the seven section values.
        std::vector<std::vector<double>> values = written("z88:elp");
        for (auto& v : values) {
            std::size_t n = 0;
            while (n < v.size() && !std::isnan(v[n]))
                ++n;
            v.resize(n >= 7 ? n : 0);
        }
        const std::vector<Run> runs = ranges(values);
        if (!runs.empty()) {
            std::string body = std::to_string(runs.size()) + "\n";
            for (const Run& r : runs) {
                detail::snprintf_c(buf, sizeof(buf), "%9zu %9zu", r.mFirst, r.mLast);
                body += buf;
                for (std::size_t k = 0; k < r.mValue.size(); ++k) {
                    if (k == 7)
                        detail::snprintf_c(buf, sizeof(buf), " %lld",
                                           static_cast<long long>(r.mValue[k]));
                    else
                        detail::snprintf_c(buf, sizeof(buf), " %+.16E", r.mValue[k]);
                    body += buf;
                }
                body += "\n";
            }
            files["z88elp.txt"] = body;
        }
    }
    if (rMesh.HasCellData("z88:int")) {
        std::vector<std::vector<double>> values = written("z88:int");
        for (auto& v : values)
            if (v.size() != 2 || v[0] < 0)
                v.clear();
        const std::vector<Run> runs = ranges(values);
        if (!runs.empty()) {
            std::string body = std::to_string(runs.size()) + "\n";
            for (const Run& r : runs) {
                detail::snprintf_c(buf, sizeof(buf), "%9zu %9zu %lld %lld\n", r.mFirst, r.mLast,
                                   static_cast<long long>(r.mValue[0]),
                                   static_cast<long long>(r.mValue[1]));
                body += buf;
            }
            files["z88int.txt"] = body;
        }
    }
    return files;
}

}  // namespace

void write_z88(const std::string& rPath, const Mesh& rMesh, bool Stubs) {
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t npts = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    bool any_3d_cell = false;
    for (const auto cb : rMesh.CellRange())
        if (!cb.IsRagged() && cell_type_dimension(cell_type_from_name(std::string(cb.Type()))) == 3)
            any_3d_cell = true;
    bool flat = pdim < 3;
    if (!flat) {
        flat = true;
        for (std::size_t p = 0; p < npts && flat; ++p)
            flat = detail::read_double(points, p * pdim + 2) == 0.0;
    }
    // Types that only exist in a 3-D file (solids, 3-D beams, trusses and the
    // shaft, shells) make it 3-D even when every z is 0.
    const bool has_type = rMesh.HasCellData("z88:type");
    bool needs_3d = false;
    for (std::size_t b = 0; has_type && b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged())
            continue;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::int64_t want = detail::read_int(rMesh.CellData("z88:type", b), r);
            const Z88Type* t = z88_type(want);
            needs_3d = needs_3d || (t && t->mKeep == t->mNodes && cb.Type() == t->mCell &&
                                    z88_needs_3d(static_cast<int>(want)));
        }
    }
    const int ndim = (flat && !any_3d_cell && !needs_3d) ? 2 : 3;

    // The Z88 type of every block (0: dropped).
    std::vector<std::vector<int>> codes(rMesh.NumCellBlocks());
    std::set<std::string> dropped_types;
    std::size_t nelem = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string cell(cb.Type());
        const int fallback = cb.IsRagged() ? 0 : z88_default_code(cell, ndim);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            int code = fallback;
            if (has_type && !cb.IsRagged()) {
                const std::int64_t want = detail::read_int(rMesh.CellData("z88:type", b), r);
                const Z88Type* t = z88_type(want);
                if (t && t->mKeep == t->mNodes && cell == t->mCell)
                    code = static_cast<int>(want);
            }
            if (!code)
                dropped_types.insert(cell);
            else
                ++nelem;
            codes[b].push_back(code);
        }
    }
    for (const std::string& t : dropped_types) {
        log::warn("Z88 writer: '{}' cells have no Z88 element type here; dropped", t);
        detail::provenance_note("cells-dropped", "Z88 has no element type for '" + t + "' cells");
    }
    if (rMesh.NumRegions()) {
        log::warn("Z88 structure files hold no groups; {} region(s) dropped", rMesh.NumRegions());
        detail::provenance_note("regions-dropped", "a Z88 structure file holds no groups");
    }
    std::size_t deck = 0;
    for (const char* name : {"z88:bc:u", "z88:bc:f"})
        deck += rMesh.HasPointData(name) ? 1 : 0;
    for (const char* name : {"z88:material", "z88:E", "z88:nu", "z88:elp", "z88:int"})
        deck += rMesh.HasCellData(name) ? 1 : 0;
    const std::size_t other_data = rMesh.NumPointData() + rMesh.NumFieldData() +
                                   rMesh.NumCellData() - (has_type ? 1 : 0) - deck;
    if (other_data) {
        log::warn(
            "Z88 decks hold no other data arrays; point, cell and field data other than the z88: "
            "constraints, materials and element parameters dropped");
        detail::provenance_note("data-dropped", "a Z88 deck holds no other data arrays");
    }
    if (!nelem)
        throw WriteError("Z88 writer: no cell has a Z88 element type");

    // Degrees of freedom per node: the most any of its elements needs.
    std::vector<int> dof(npts, ndim == 2 ? 2 : 3);
    std::vector<bool> seen(npts, false);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged())
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            if (!codes[b][r])
                continue;
            const int d = z88_type(codes[b][r])->mDof;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t p = static_cast<std::size_t>(detail::read_int(conn, r * k + j));
                dof[p] = seen[p] ? std::max(dof[p], d) : d;
                seen[p] = true;
            }
        }
    }
    long long total = 0;
    for (int d : dof)
        total += d;

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    std::string out;
    char buf[160];
    detail::snprintf_c(buf, sizeof(buf), "%5d %9zu %9zu %11lld %5d   ", ndim, npts, nelem, total,
                       0);
    out += buf;
    out += detail::provenance_lines(detail::SlotTier::SingleLine)[0];
    out += '\n';
    for (std::size_t p = 0; p < npts; ++p) {
        detail::snprintf_c(buf, sizeof(buf), "%9zu %2d", p + 1, dof[p]);
        out += buf;
        for (int d = 0; d < ndim; ++d) {
            const double v =
                static_cast<std::size_t>(d) < pdim
                    ? detail::read_double(points, p * pdim + static_cast<std::size_t>(d))
                    : 0.0;
            detail::snprintf_c(buf, sizeof(buf), " %+.16E", v);
            out += buf;
        }
        out += '\n';
        if (out.size() > (1u << 20)) {
            f << out;
            out.clear();
        }
    }
    std::size_t id = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged())
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const detail::NodeOrder* order = detail::node_order("z88", std::string(cb.Type()));
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            if (!codes[b][r])
                continue;
            detail::snprintf_c(buf, sizeof(buf), "%9zu %5d\n", ++id, codes[b][r]);
            out += buf;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                detail::snprintf_c(buf, sizeof(buf), "%s%lld", j ? " " : "",
                                   static_cast<long long>(detail::read_int(conn, r * k + src) + 1));
                out += buf;
            }
            out += '\n';
            if (out.size() > (1u << 20)) {
                f << out;
                out.clear();
            }
        }
    }
    f << out;
    if (!f)
        throw WriteError("Z88 writer: failed writing " + rPath);
    std::map<std::string, std::string> files = z88_deck_files(rMesh, codes, dof);
    if (Stubs) {
        files.emplace("z88i2.txt", "0\n");
        files.emplace("z88i5.txt", "0\n");
    }
    const std::filesystem::path dir = std::filesystem::path(rPath).parent_path();
    for (const auto& [name, body] : files) {
        auto s = detail::make_classic_ofstream((dir / name).string(), std::ios::binary);
        if (!s)
            throw WriteError(std::string("Could not open file for writing: ") +
                             (dir / name).string());
        s << body;
        if (!s)
            throw WriteError("Z88 writer: failed writing " + (dir / name).string());
    }
}

}  // namespace meshioplusplus
