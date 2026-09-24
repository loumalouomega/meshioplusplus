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

// `z88o3.txt`: blocks opened by `element # = N` (or `#=`), then rows of reals.
// Solid rows are `x y z` + 6 stresses [+ equivalent], plane-stress rows `x y` +
// 3 stresses [+ equivalent]; the mean over an element's rows is kept.
void z88_attach_stresses(Mesh& rMesh, const std::string& rPath,
                         const std::unordered_map<std::int64_t, std::size_t>& rElementIndex,
                         const std::vector<int>& rCellCode,
                         const std::vector<std::size_t>& rBlockStart) {
    const std::string text = z88_read_text(rPath);
    const std::size_t ncells = rCellCode.size();
    std::vector<std::vector<double>> sum(ncells);
    std::vector<std::vector<double>> sumv(ncells);
    std::vector<std::size_t> count(ncells, 0), countv(ncells, 0);
    std::int64_t current = -1;
    std::size_t skipped = 0;
    int components = 0;
    for (std::string_view line : z88_lines(text)) {
        const std::size_t hash = line.find('#');
        const std::string lower = z88_lower(std::string(line.substr(0, hash)));
        if (hash != std::string_view::npos && lower.find("element") != std::string::npos) {
            std::string_view rest = line.substr(hash + 1);
            const std::size_t eq = rest.find('=');
            current = -1;
            if (eq != std::string_view::npos) {
                const std::vector<std::int64_t> id = z88_leading_ints(rest.substr(eq + 1));
                if (!id.empty()) {
                    const auto it = rElementIndex.find(id[0]);
                    if (it != rElementIndex.end())
                        current = static_cast<std::int64_t>(it->second);
                }
            }
            continue;
        }
        if (current < 0)
            continue;
        const std::vector<std::string_view> t = z88_tokens(line);
        if (t.empty())
            continue;
        std::vector<double> v(t.size());
        bool ok = true;
        for (std::size_t k = 0; k < t.size() && ok; ++k)
            ok = z88_real(t[k], v[k]);
        if (!ok)
            continue;
        const std::size_t c = static_cast<std::size_t>(current);
        const int code = rCellCode[c];
        const bool solid = code == 1 || code == 10 || code == 16 || code == 17;
        const bool plane = code == 3 || code == 7 || code == 11 || code == 14;
        std::size_t first = 0, n = 0;
        if (solid && (v.size() == 9 || v.size() == 10)) {
            first = 3;
            n = 6;
        } else if (plane && (v.size() == 5 || v.size() == 6)) {
            first = 2;
            n = 3;
        } else {
            ++skipped;
            continue;
        }
        if (components && components != static_cast<int>(n)) {
            ++skipped;
            continue;
        }
        components = static_cast<int>(n);
        std::vector<double>& s = sum[c];
        s.resize(n, 0.0);
        for (std::size_t k = 0; k < n; ++k)
            s[k] += v[first + k];
        ++count[c];
        if (v.size() == first + n + 1) {
            sumv[c].resize(1, 0.0);
            sumv[c][0] += v.back();
            ++countv[c];
        }
    }
    if (skipped)
        log::warn(
            "Z88: {} stress row(s) of elements other than solids and plane-stress "
            "elements skipped",
            skipped);
    if (!components)
        return;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t w = static_cast<std::size_t>(components);
    std::vector<NDArray> sig, sigv;
    bool any_v = false;
    for (std::size_t b = 0; b + 1 < rBlockStart.size(); ++b) {
        const std::size_t n = rBlockStart[b + 1] - rBlockStart[b];
        NDArray a(DType::Float64, {n, w});
        NDArray av(DType::Float64, {n});
        for (std::size_t r = 0; r < n; ++r) {
            const std::size_t c = rBlockStart[b] + r;
            for (std::size_t k = 0; k < w; ++k)
                a.As<double>()[r * w + k] =
                    count[c] ? sum[c][k] / static_cast<double>(count[c]) : nan;
            av.As<double>()[r] = countv[c] ? sumv[c][0] / static_cast<double>(countv[c]) : nan;
            any_v = any_v || countv[c] > 0;
        }
        sig.push_back(std::move(a));
        sigv.push_back(std::move(av));
    }
    rMesh.AddCellData("SIG", std::move(sig));
    if (any_v)
        rMesh.AddCellData("SIGV", std::move(sigv));
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

    if (Results) {
        const fs::path dir = structure.parent_path();
        const fs::path o2 = z88_sibling(dir, "z88o2.txt");
        if (!o2.empty())
            z88_attach_displacements(mesh, o2.string(), node_index);
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
    const int ndim = (flat && !any_3d_cell) ? 2 : 3;

    // The Z88 type of every block (0: dropped).
    const bool has_type = rMesh.HasCellData("z88:type");
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
    const std::size_t other_data =
        rMesh.NumPointData() + rMesh.NumFieldData() + rMesh.NumCellData() - (has_type ? 1 : 0);
    if (other_data) {
        log::warn("Z88 structure files hold no data arrays; point, cell and field data dropped");
        detail::provenance_note("data-dropped", "a Z88 structure file holds no data arrays");
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
    if (Stubs) {
        const std::filesystem::path dir = std::filesystem::path(rPath).parent_path();
        for (const char* name : {"z88i2.txt", "z88i5.txt"}) {
            auto s = detail::make_classic_ofstream((dir / name).string(), std::ios::binary);
            if (!s)
                throw WriteError(std::string("Could not open file for writing: ") +
                                 (dir / name).string());
            s << "0\n";
        }
    }
}

}  // namespace meshioplusplus
