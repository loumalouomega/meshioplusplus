#include "meshio/formats/su2.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/exceptions.hpp"

namespace meshio {

namespace {

int su2_numnodes(int t) {
    switch (t) {
        case 3: return 2;   // line
        case 5: return 3;   // triangle
        case 9: return 4;   // quad
        case 10: return 4;  // tetra
        case 12: return 8;  // hexahedron
        case 13: return 6;  // wedge
        case 14: return 5;  // pyramid
        default: return 0;
    }
}
std::string su2_to_meshio(int t) {
    switch (t) {
        case 3: return "line";       case 5: return "triangle";
        case 9: return "quad";       case 10: return "tetra";
        case 12: return "hexahedron";case 13: return "wedge";
        case 14: return "pyramid";   default: return "";
    }
}
int meshio_to_su2(const std::string& t) {
    if (t == "line") return 3;        if (t == "triangle") return 5;
    if (t == "quad") return 9;        if (t == "tetra") return 10;
    if (t == "hexahedron") return 12; if (t == "wedge") return 13;
    if (t == "pyramid") return 14;    return -1;
}

std::string strip(const std::string& s) {
    std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    std::size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}
std::vector<std::string> tokens(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

struct Blk {
    std::string type;
    int n = 0;
    std::vector<std::int64_t> conn;
    std::vector<std::int32_t> tag;
    std::size_t count = 0;
};

// Parse `count` element lines (each "vtk_type n0 n1 ... [extra]") into type-
// grouped blocks (sorted by vtk type code, matching numpy.unique), all with
// the given tag.
void read_elem_block(const std::vector<std::string>& lines, std::size_t& li,
                     std::size_t count, std::int32_t tag, std::vector<Blk>& out) {
    std::vector<std::pair<int, std::vector<std::int64_t>>> elems;
    std::set<int> types;
    for (std::size_t e = 0; e < count; ++e) {
        auto t = tokens(lines.at(li++));
        int vt = std::stoi(t[0]);
        int nn = su2_numnodes(vt);
        if (nn == 0) throw ReadError("SU2: unsupported element type " + t[0]);
        std::vector<std::int64_t> nodes(nn);
        for (int j = 0; j < nn; ++j) nodes[j] = std::strtoll(t[1 + j].c_str(), nullptr, 10);
        elems.emplace_back(vt, std::move(nodes));
        types.insert(vt);
    }
    for (int vt : types) {  // std::set is sorted
        Blk b;
        b.type = su2_to_meshio(vt);
        b.n = su2_numnodes(vt);
        for (auto& e : elems) {
            if (e.first != vt) continue;
            b.conn.insert(b.conn.end(), e.second.begin(), e.second.end());
            b.tag.push_back(tag);
            ++b.count;
        }
        out.push_back(std::move(b));
    }
}

}  // namespace

Mesh read_su2(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw ReadError("Could not open file: " + path);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) lines.push_back(l);

    int dim = 0;
    Mesh mesh;
    std::vector<Blk> blocks;
    std::int32_t next_tag_id = 0;

    std::size_t li = 0;
    while (li < lines.size()) {
        std::string line = strip(lines[li]);
        if (line.empty() || line[0] == '%') { ++li; continue; }
        std::size_t eq = line.find('=');
        if (eq == std::string::npos) { ++li; continue; }
        std::string name = strip(line.substr(0, eq));
        std::string rest = strip(line.substr(eq + 1));
        ++li;

        if (name == "NDIME") {
            dim = std::stoi(rest);
            if (dim != 2 && dim != 3) throw ReadError("SU2: invalid NDIME");
        } else if (name == "NPOIN") {
            std::size_t npoin = static_cast<std::size_t>(std::stoll(tokens(rest)[0]));
            mesh.points = NDArray(DType::Float64, {npoin, static_cast<std::size_t>(dim)});
            double* pp = mesh.points.as<double>();
            for (std::size_t i = 0; i < npoin; ++i) {
                auto t = tokens(lines.at(li++));
                for (int c = 0; c < dim; ++c) pp[i * dim + c] = std::strtod(t[c].c_str(), nullptr);
            }
        } else if (name == "NELEM") {
            std::size_t ne = static_cast<std::size_t>(std::stoll(rest));
            read_elem_block(lines, li, ne, 0, blocks);
        } else if (name == "NMARK") {
            // handled implicitly via MARKER_TAG/MARKER_ELEMS
        } else if (name == "MARKER_TAG") {
            try {
                std::size_t pos;
                int v = std::stoi(rest, &pos);
                if (pos == rest.size()) next_tag_id = v;
                else { ++next_tag_id; }
            } catch (...) {
                ++next_tag_id;
            }
        } else if (name == "MARKER_ELEMS") {
            std::size_t ne = static_cast<std::size_t>(std::stoll(rest));
            read_elem_block(lines, li, ne, next_tag_id, blocks);
        }
    }

    // Merge boundary blocks of the same type (lines in 2D; tris/quads in 3D).
    std::vector<std::string> btypes = (dim == 2) ? std::vector<std::string>{"line"}
                                                 : std::vector<std::string>{"triangle", "quad"};
    for (const auto& bt : btypes) {
        int first = -1;
        for (std::size_t i = 0; i < blocks.size(); ++i) {
            if (blocks[i].type != bt) continue;
            if (first < 0) { first = static_cast<int>(i); continue; }
            Blk& dst = blocks[first];
            Blk& src = blocks[i];
            dst.conn.insert(dst.conn.end(), src.conn.begin(), src.conn.end());
            dst.tag.insert(dst.tag.end(), src.tag.begin(), src.tag.end());
            dst.count += src.count;
            src.count = 0;  // mark for removal
            src.conn.clear();
        }
    }

    std::vector<NDArray> tags;
    for (auto& b : blocks) {
        if (b.count == 0) continue;  // merged-away or empty
        NDArray data(DType::Int64, {b.count, static_cast<std::size_t>(b.n)});
        std::memcpy(data.data(), b.conn.data(), b.conn.size() * sizeof(std::int64_t));
        mesh.cells.emplace_back(b.type, std::move(data));
        NDArray tg(DType::Int32, {b.count});
        std::memcpy(tg.data(), b.tag.data(), b.tag.size() * sizeof(std::int32_t));
        tags.push_back(std::move(tg));
    }
    mesh.cell_data.emplace("su2:tag", std::move(tags));
    return mesh;
}

void write_su2(const std::string& path, const Mesh& mesh) {
    std::ofstream os(path);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    const std::size_t npoin = mesh.num_points();

    os << "NDIME= " << dim << "\n";
    os << "NPOIN= " << npoin << "\n";
    char buf[64];
    for (std::size_t i = 0; i < npoin; ++i) {
        for (std::size_t c = 0; c < dim; ++c) {
            std::snprintf(buf, sizeof(buf), "%.16e", detail::read_double(mesh.points, i * dim + c));
            os << buf << (c + 1 == dim ? '\n' : ' ');
        }
    }

    std::vector<std::string> vtypes = (dim == 2)
        ? std::vector<std::string>{"triangle", "quad"}
        : std::vector<std::string>{"tetra", "hexahedron", "wedge", "pyramid"};
    std::vector<std::string> btypes = (dim == 2) ? std::vector<std::string>{"line"}
                                                 : std::vector<std::string>{"triangle", "quad"};
    auto in = [](const std::vector<std::string>& v, const std::string& t) {
        return std::find(v.begin(), v.end(), t) != v.end();
    };

    // Volume cells.
    std::size_t nelem = 0;
    for (const auto& cb : mesh.cells)
        if (in(vtypes, cb.type)) nelem += cb.num_cells();
    os << "NELEM= " << nelem << "\n";
    for (const auto& cb : mesh.cells) {
        if (!in(vtypes, cb.type)) continue;
        int st = meshio_to_su2(cb.type);
        std::size_t k = cb.data.shape().size() >= 2 ? cb.data.shape()[1] : 1;
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            os << st;
            for (std::size_t j = 0; j < k; ++j) os << " " << detail::read_int(cb.data, r * k + j);
            os << "\n";
        }
    }

    // Boundary markers from su2:tag (first int cell_data).
    std::string tag_key;
    for (const auto& kv : mesh.cell_data) {
        if (kv.second.empty()) continue;
        DType t = kv.second.front().dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 ||
            t == DType::Int64 || t == DType::UInt8 || t == DType::UInt16 ||
            t == DType::UInt32 || t == DType::UInt64) { tag_key = kv.first; break; }
    }

    // Collect unique tags (with total counts) over boundary cell blocks.
    std::map<std::int64_t, std::size_t> tag_counts;
    for (std::size_t bi = 0; bi < mesh.cells.size(); ++bi) {
        const CellBlock& cb = mesh.cells[bi];
        if (!in(btypes, cb.type)) continue;
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            std::int64_t tg = 1;
            if (!tag_key.empty()) tg = detail::read_int(mesh.cell_data.at(tag_key)[bi], r);
            ++tag_counts[tg];
        }
    }

    os << "NMARK= " << tag_counts.size() << "\n";
    for (const auto& tc : tag_counts) {
        std::int64_t tag = tc.first;
        os << "MARKER_TAG= " << tag << "\n";
        os << "MARKER_ELEMS= " << tc.second << "\n";
        for (std::size_t bi = 0; bi < mesh.cells.size(); ++bi) {
            const CellBlock& cb = mesh.cells[bi];
            if (!in(btypes, cb.type)) continue;
            int st = meshio_to_su2(cb.type);
            std::size_t k = cb.data.shape().size() >= 2 ? cb.data.shape()[1] : 1;
            for (std::size_t r = 0; r < cb.num_cells(); ++r) {
                std::int64_t tg = 1;
                if (!tag_key.empty()) tg = detail::read_int(mesh.cell_data.at(tag_key)[bi], r);
                if (tg != tag) continue;
                os << st;
                for (std::size_t j = 0; j < k; ++j) os << " " << detail::read_int(cb.data, r * k + j);
                os << "\n";
            }
        }
    }
}

}  // namespace meshio
