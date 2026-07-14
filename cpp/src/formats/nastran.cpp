#include "meshioplusplus/formats/nastran.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kSentinel = "meshioplusplus-cpp-nastran";

const std::unordered_map<std::string, std::string>& nastran_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = {
        {"CTRIA3", "triangle"}, {"CTRIA6", "triangle6"}, {"CQUAD4", "quad"},
        {"CQUAD8", "quad8"},    {"CQUAD9", "quad9"},     {"CTETRA", "tetra"},
        {"CTETRA_", "tetra10"}, {"CPYRA", "pyramid"},    {"CPYRA_", "pyramid13"},
        {"CPENTA", "wedge"},    {"CPENTA_", "wedge15"},  {"CHEXA", "hexahedron"},
        {"CHEXA_", "hexahedron20"}, {"CBAR", "line"}, {"CROD", "line"},
    };
    return m;
}
// meshio -> nastran (matches the Python inverse: last entry per meshio type).
const std::unordered_map<std::string, std::string>& meshio_to_nastran() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "CELAS1"},   {"line", "CBAR"},      {"triangle", "CTRIA3"},
        {"triangle6", "CTRIA6"},{"quad", "CQUAD4"},    {"quad8", "CQUAD8"},
        {"quad9", "CQUAD9"},    {"tetra", "CTETRA"},   {"tetra10", "CTETRA_"},
        {"pyramid", "CPYRA"},   {"pyramid13", "CPYRA_"},{"wedge", "CPENTA"},
        {"wedge15", "CPENTA_"}, {"hexahedron", "CHEXA"},{"hexahedron20", "CHEXA_"},
    };
    return m;
}

// Node reordering between meshio (VTK-like) and Nastran for the few types that
// differ. The given permutation P maps: out[j] = in[P[j]].
const std::vector<int>& reorder_meshio_to_nastran(const std::string& nastran_type) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"CHEXA_", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"CPENTA_", {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}},
    };
    static const std::vector<int> empty;
    auto it = m.find(nastran_type);
    return it == m.end() ? empty : it->second;
}
// Inverse (Nastran -> meshio). CHEXA_/CPENTA_ permutations are involutions.
const std::vector<int>& reorder_nastran_to_meshio(const std::string& nastran_type) {
    return reorder_meshio_to_nastran(nastran_type);
}

std::string nastran_float(double v) {
    if (v == 0.0) return "0.0";
    char buf[40];
    std::string best;
    for (int p = 0; p <= 11; ++p) {
        std::snprintf(buf, sizeof(buf), "%.*E", p, v);
        if (std::strtod(buf, nullptr) == v) { best = buf; break; }
    }
    if (best.empty()) {
        std::snprintf(buf, sizeof(buf), "%.11E", v);
        best = buf;
    }
    std::size_t epos = best.find('E');
    std::string mant = best.substr(0, epos);
    int exp = std::atoi(best.c_str() + epos + 1);
    std::size_t dot = mant.find('.');
    if (dot == std::string::npos) {
        mant += ".";
        dot = mant.size() - 1;
    }
    // trim trailing zeros after the decimal point (keep the dot)
    std::size_t last = mant.size();
    while (last > dot + 1 && mant[last - 1] == '0') --last;
    mant.erase(last);
    std::string es = (exp < 0 ? "-" : "+") + std::to_string(std::abs(exp));
    std::string out = mant + "E" + es;
    // Keep within the 16-char field by shedding mantissa precision if needed.
    while (out.size() > 16 && mant.find('.') != std::string::npos && mant.back() != '.') {
        mant.pop_back();
        out = mant + "E" + es;
    }
    return out;
}

double parse_nastran_float(std::string s) {
    // strip
    std::size_t b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return 0.0;
    std::size_t e = s.find_last_not_of(" \t");
    s = s.substr(b, e - b + 1);
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp != s.c_str() && *endp == '\0') return v;
    // Nastran compressed exponent, e.g. "1.5+1" -> "1.5e+1"
    std::string t;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if ((c == '+' || c == '-') && i > 0 && s[i - 1] != 'e' && s[i - 1] != 'E')
            t += 'e';
        t += c;
    }
    return std::strtod(t.c_str(), nullptr);
}

std::string strip(const std::string& s) {
    std::size_t b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return "";
    std::size_t e = s.find_last_not_of(" \t");
    return s.substr(b, e - b + 1);
}

std::string field(const std::string& line, std::size_t start, std::size_t width) {
    if (start >= line.size()) return "";
    return strip(line.substr(start, width));
}

}  // namespace

void write_nastran(const std::string& path, const Mesh& mesh) {
    std::ofstream os(path);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t n = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    os << "$ " << kSentinel << "\n";
    os << "BEGIN BULK\n";

    // Points: fixed-large GRID*.
    char buf[128];
    for (std::size_t i = 0; i < n; ++i) {
        double xyz[3] = {0, 0, 0};
        for (std::size_t c = 0; c < dim && c < 3; ++c) xyz[c] = detail::read_double(mesh.points, i * dim + c);
        std::string sx = nastran_float(xyz[0]), sy = nastran_float(xyz[1]), sz = nastran_float(xyz[2]);
        std::snprintf(buf, sizeof(buf), "GRID*   %-16d%-16s%16s%16s\n*       %16s\n",
                      static_cast<int>(i + 1), "", sx.c_str(), sy.c_str(), sz.c_str());
        os << buf;
    }

    // Cells: fixed-small element cards (8-char fields), with + continuations.
    const auto& m2n = meshio_to_nastran();
    std::size_t cell_id = 0;
    for (const auto& cb : mesh.cells) {
        auto it = m2n.find(cb.type);
        if (it == m2n.end()) throw WriteError("Nastran writer: unsupported cell type " + cb.type);
        std::string ntype = it->second;
        std::size_t k = cb.data.shape().size() >= 2 ? cb.data.shape()[1] : 1;
        const std::vector<int>& perm = reorder_meshio_to_nastran(ntype);
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            ++cell_id;
            std::vector<long long> nodes(k);
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                nodes[j] = detail::read_int(cb.data, r * k + src) + 1;
            }
            // first line: type, id, ref, up to 6 nodes
            std::snprintf(buf, sizeof(buf), "%-8s%-8d%-8s", ntype.c_str(),
                          static_cast<int>(cell_id), "");
            std::string line = buf;
            std::size_t nipl1 = 6, nipl2 = 14;
            for (std::size_t j = 0; j < k && j < nipl1; ++j) {
                std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                line += buf;
            }
            if (k > nipl1) {
                std::snprintf(buf, sizeof(buf), "+1%-6x", static_cast<unsigned>(cell_id));
                os << line << buf << "\n";
                std::snprintf(buf, sizeof(buf), "+1%-6x", static_cast<unsigned>(cell_id));
                std::string l2 = buf;
                for (std::size_t j = nipl1; j < k && j < nipl2; ++j) {
                    std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                    l2 += buf;
                }
                if (k > nipl2) {
                    std::snprintf(buf, sizeof(buf), "+2%-6x", static_cast<unsigned>(cell_id));
                    os << l2 << buf << "\n";
                    std::snprintf(buf, sizeof(buf), "+2%-6x", static_cast<unsigned>(cell_id));
                    std::string l3 = buf;
                    for (std::size_t j = nipl2; j < k; ++j) {
                        std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                        l3 += buf;
                    }
                    os << l3 << "\n";
                } else {
                    os << l2 << "\n";
                }
            } else {
                os << line << "\n";
            }
        }
    }

    os << "ENDDATA\n";
}

Mesh read_nastran(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw ReadError("Could not open file: " + path);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        lines.push_back(l);
    }

    // Sentinel gate: only parse files this writer produced.
    bool ok = false;
    std::size_t start = 0;
    for (; start < lines.size(); ++start) {
        if (lines[start].find(kSentinel) != std::string::npos) ok = true;
        if (strip(lines[start]).rfind("BEGIN BULK", 0) == 0) { ++start; break; }
    }
    if (!ok) throw ReadError("Not a meshio++-C++ Nastran file");

    const auto& n2m = nastran_to_meshio();
    Mesh mesh;
    std::unordered_map<std::int64_t, std::int64_t> point_ids;
    std::vector<std::array<double, 3>> pts;

    struct Blk { std::string type; int n; std::vector<std::int64_t> conn; std::size_t count = 0; };
    std::vector<Blk> blocks;

    std::size_t i = start;
    while (i < lines.size()) {
        const std::string& line = lines[i];
        std::string s = strip(line);
        if (s.empty() || s[0] == '$' || s.rfind("//", 0) == 0 || s[0] == '#') { ++i; continue; }
        if (s.rfind("ENDDATA", 0) == 0) break;

        std::string kw = field(line, 0, 8);
        if (kw == "GRID*") {
            // line1: 8 + 4x16 (id, ref, x, y); line2: 8 + 16 (z)
            std::int64_t id = std::strtoll(field(line, 8, 16).c_str(), nullptr, 10);
            double x = parse_nastran_float(field(line, 40, 16));
            double y = parse_nastran_float(field(line, 56, 16));
            double z = 0.0;
            if (i + 1 < lines.size()) z = parse_nastran_float(field(lines[i + 1], 8, 16));
            point_ids[id] = static_cast<std::int64_t>(pts.size());
            pts.push_back({x, y, z});
            i += 2;
        } else if (n2m.count(kw)) {
            std::string mtype = n2m.at(kw);
            // gather node fields: first line fields[3..9] (chars 24..72),
            // continuation lines fields[1..9] (chars 8..72).
            std::vector<std::int64_t> nodes;
            // Field 9 (chars 72..80) holds the continuation marker, never a node.
            auto grab = [&](const std::string& ln, std::size_t first_field) {
                for (std::size_t fidx = first_field; fidx < 9; ++fidx) {
                    std::string f = field(ln, fidx * 8, 8);
                    if (!f.empty()) nodes.push_back(std::strtoll(f.c_str(), nullptr, 10));
                }
            };
            grab(line, 3);
            ++i;
            while (i < lines.size() && !lines[i].empty() &&
                   (lines[i][0] == '+' || lines[i][0] == '*')) {
                grab(lines[i], 1);
                ++i;
            }
            int nn = num_nodes_per_cell().count(mtype) ? num_nodes_per_cell().at(mtype) : (int)nodes.size();
            if ((int)nodes.size() != nn)
                throw ReadError("Nastran: node count mismatch for " + kw);
            const std::vector<int>& perm = reorder_nastran_to_meshio(kw);
            if (blocks.empty() || blocks.back().type != mtype) {
                Blk b; b.type = mtype; b.n = nn; blocks.push_back(std::move(b));
            }
            Blk& blk = blocks.back();
            for (int j = 0; j < nn; ++j) {
                int src = perm.empty() ? j : perm[j];
                blk.conn.push_back(nodes[src]);  // 1-based gmsh-ish id
            }
            ++blk.count;
        } else {
            ++i;
        }
    }

    // Points + remap.
    mesh.points = NDArray(DType::Float64, {pts.size(), 3});
    double* pp = mesh.points.as<double>();
    for (std::size_t r = 0; r < pts.size(); ++r)
        for (int c = 0; c < 3; ++c) pp[r * 3 + c] = pts[r][c];

    for (auto& blk : blocks) {
        NDArray data(DType::Int64, {blk.count, static_cast<std::size_t>(blk.n)});
        std::int64_t* dp = data.as<std::int64_t>();
        for (std::size_t idx = 0; idx < blk.conn.size(); ++idx) {
            auto it = point_ids.find(blk.conn[idx]);
            if (it == point_ids.end()) throw ReadError("Nastran: unknown node id");
            dp[idx] = it->second;
        }
        mesh.cells.emplace_back(blk.type, std::move(data));
    }
    return mesh;
}

}  // namespace meshioplusplus
