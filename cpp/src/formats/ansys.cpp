#include "meshio/formats/ansys.hpp"

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/exceptions.hpp"

namespace meshio {

namespace {

// Cursor over the whole file, mixing line reads with raw binary reads.
struct Buf {
    std::string data;
    std::size_t p = 0;

    bool eof() const { return p >= data.size(); }

    std::string readline() {
        if (p >= data.size()) return "";
        std::size_t nl = data.find('\n', p);
        std::string line;
        if (nl == std::string::npos) {
            line = data.substr(p);
            p = data.size();
        } else {
            line = data.substr(p, nl - p + 1);
            p = nl + 1;
        }
        return line;
    }

    void skip_close(int n) {
        while (n > 0 && p < data.size()) {
            char c = data[p++];
            if (c == '(') ++n;
            else if (c == ')') --n;
        }
    }

    void advance_to(char ch) {
        while (p < data.size() && data[p] != ch) ++p;
        if (p < data.size()) ++p;  // consume it
    }

    const char* raw(std::size_t nbytes) {
        if (p + nbytes > data.size()) throw ReadError("ANSYS: unexpected end of file");
        const char* ptr = data.data() + p;
        p += nbytes;
        return ptr;
    }
};

int count_char(const std::string& s, char c) {
    int n = 0;
    for (char ch : s) if (ch == c) ++n;
    return n;
}

std::string rstrip(const std::string& s) {
    std::size_t b = s.size();
    while (b > 0 && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(0, b);
}

// Parse the bracketed "first second ... " hex group, e.g. "(... (1 1 4 1 3) ...".
std::vector<std::int64_t> parse_header_nums(const std::string& line) {
    std::size_t o1 = line.find('(');
    std::size_t o2 = (o1 == std::string::npos) ? std::string::npos : line.find('(', o1 + 1);
    std::size_t c2 = (o2 == std::string::npos) ? std::string::npos : line.find(')', o2 + 1);
    if (c2 == std::string::npos) throw ReadError("ANSYS: malformed section header");
    std::string nums = line.substr(o2 + 1, c2 - o2 - 1);
    std::vector<std::int64_t> a;
    std::istringstream iss(nums);
    std::string t;
    while (iss >> t) a.push_back(std::strtoll(t.c_str(), nullptr, 16));
    return a;
}

// Leading "(" + ws + digits -> the index string; "" if not a section line.
std::string section_index(const std::string& line) {
    std::size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    if (i >= line.size() || line[i] != '(') return "";
    ++i;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    std::size_t s = i;
    while (i < line.size() && std::isdigit(static_cast<unsigned char>(line[i]))) ++i;
    return line.substr(s, i - s);
}

// "" / "20" / "30" prefix on a 10/12/13 core; returns false if not points/cells/faces.
bool classify(const std::string& idx, int& core, std::string& prefix) {
    static const std::unordered_map<std::string, std::pair<int, std::string>> m = {
        {"10", {10, ""}},   {"2010", {10, "20"}}, {"3010", {10, "30"}},
        {"12", {12, ""}},   {"2012", {12, "20"}}, {"3012", {12, "30"}},
        {"13", {13, ""}},   {"2013", {13, "20"}}, {"3013", {13, "30"}}};
    auto it = m.find(idx);
    if (it == m.end()) return false;
    core = it->second.first;
    prefix = it->second.second;
    return true;
}

const std::unordered_map<int, std::pair<std::string, int>>& cell_type_map() {
    static const std::unordered_map<int, std::pair<std::string, int>> m = {
        {1, {"triangle", 3}}, {2, {"tetra", 4}},   {3, {"quad", 4}},
        {4, {"hexahedron", 8}}, {5, {"pyramid", 5}}, {6, {"wedge", 6}}};
    return m;
}

}  // namespace

Mesh read_ansys(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    Buf buf;
    buf.data.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    std::vector<double> points;  // flat
    std::size_t dim = 3;
    std::int64_t npoints = 0;
    std::int64_t first_point_index_overall = -1;

    struct RawCell { std::string type; std::vector<std::int64_t> data; std::size_t rows; std::size_t cols; };
    std::vector<RawCell> cells;

    while (!buf.eof()) {
        std::string line = buf.readline();
        if (line.empty()) break;
        // blank?
        bool blank = true;
        for (char c : line) if (!std::isspace(static_cast<unsigned char>(c))) { blank = false; break; }
        if (blank) continue;

        std::string idx = section_index(line);
        if (idx.empty()) throw ReadError("ANSYS: expected a section line");

        int core;
        std::string prefix;
        if (idx == "0" || idx == "1" || idx == "2" || idx == "39" || idx == "45") {
            buf.skip_close(count_char(line, '(') - count_char(line, ')'));
            continue;
        }
        if (!classify(idx, core, prefix)) {
            buf.skip_close(count_char(line, '(') - count_char(line, ')'));
            continue;
        }

        // Self-contained declaration line (no data block).
        if (count_char(line, '(') == count_char(line, ')')) continue;

        std::vector<std::int64_t> a = parse_header_nums(line);
        if (a.size() <= 4) throw ReadError("ANSYS: short section header");

        // Position at the data block opener.
        if (rstrip(line).back() != '(') buf.advance_to('(');

        if (core == 10) {
            std::int64_t first = a[1], last = a[2];
            std::int64_t n = last - first + 1;
            int d = static_cast<int>(a[4]);
            if (first_point_index_overall < 0) first_point_index_overall = first;
            if (points.empty()) dim = static_cast<std::size_t>(d);
            if (prefix.empty()) {
                for (std::int64_t k = 0; k < n; ++k) {
                    std::string pl = buf.readline();
                    while (rstrip(pl).empty() && !buf.eof()) pl = buf.readline();
                    std::istringstream iss(pl);
                    for (int c = 0; c < d; ++c) {
                        double v;
                        iss >> v;
                        points.push_back(v);
                    }
                }
            } else {
                std::size_t isz = (prefix == "20") ? 4 : 8;
                const char* ptr = buf.raw(static_cast<std::size_t>(n) * d * isz);
                for (std::int64_t k = 0; k < n * d; ++k) {
                    if (isz == 4) { float f; std::memcpy(&f, ptr + k * 4, 4); points.push_back(f); }
                    else { double db; std::memcpy(&db, ptr + k * 8, 8); points.push_back(db); }
                }
            }
            npoints += n;
            buf.skip_close(2);
        } else if (core == 12) {
            std::int64_t first = a[1], last = a[2];
            std::int64_t zone_type = a[3];
            int element_type = static_cast<int>(a[4]);
            std::int64_t n = last - first + 1;
            if (zone_type == 0) { buf.skip_close(2); continue; }  // dead zone
            auto tit = cell_type_map().find(element_type);
            if (tit == cell_type_map().end())
                throw ReadError("ANSYS: unsupported cell element-type");
            const std::string& key = tit->second.first;
            int npc = tit->second.second;

            std::vector<std::int64_t> cdata(static_cast<std::size_t>(n) * npc);
            if (prefix.empty()) {
                for (std::int64_t k = 0; k < n; ++k) {
                    std::string cl = buf.readline();
                    std::istringstream iss(cl);
                    std::string tok;
                    for (int c = 0; c < npc; ++c) {
                        iss >> tok;
                        cdata[k * npc + c] = std::strtoll(tok.c_str(), nullptr, 16);
                    }
                }
            } else {
                std::size_t isz = (prefix == "20") ? 4 : 8;
                const char* ptr = buf.raw(static_cast<std::size_t>(n) * npc * isz);
                for (std::int64_t k = 0; k < n * npc; ++k) {
                    if (isz == 4) { std::int32_t v; std::memcpy(&v, ptr + k * 4, 4); cdata[k] = v; }
                    else { std::int64_t v; std::memcpy(&v, ptr + k * 8, 8); cdata[k] = v; }
                }
            }
            cells.push_back({key, std::move(cdata), static_cast<std::size_t>(n),
                             static_cast<std::size_t>(npc)});
            buf.skip_close(2);
        } else {  // faces (core == 13) with a data body -> defer to Python
            throw ReadError("ANSYS: face sections handled by Python fallback");
        }
    }

    if (first_point_index_overall < 0) first_point_index_overall = 0;

    Mesh mesh;
    mesh.points = NDArray(DType::Float64,
                          {static_cast<std::size_t>(npoints), dim});
    double* pp = mesh.points.as<double>();
    for (std::size_t i = 0; i < points.size(); ++i) pp[i] = points[i];

    for (auto& rc : cells) {
        NDArray data(DType::Int64, {rc.rows, rc.cols});
        std::int64_t* dp = data.as<std::int64_t>();
        for (std::size_t k = 0; k < rc.data.size(); ++k)
            dp[k] = rc.data[k] - first_point_index_overall;
        mesh.cells.emplace_back(rc.type, std::move(data));
    }

    return mesh;
}

void write_ansys(const std::string& path, const Mesh& mesh, bool binary) {
    std::ofstream fh(path, std::ios::binary);
    if (!fh) throw WriteError("Could not open file for writing: " + path);

    const std::size_t npoints = mesh.num_points();
    const std::size_t dim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    if (dim != 2 && dim != 3)
        throw WriteError("ANSYS: can only write dimension 2 or 3");

    static const std::unordered_map<std::string, int> meshio_to_ansys = {
        {"triangle", 1}, {"tetra", 2}, {"quad", 3},
        {"hexahedron", 4}, {"pyramid", 5}, {"wedge", 6}};

    char hbuf[128];
    fh << "(1 \"meshio C++ core\")\n";
    std::snprintf(hbuf, sizeof(hbuf), "(2 %zu)\n", dim);
    fh << hbuf;

    const std::size_t first_node_index = 1;
    std::snprintf(hbuf, sizeof(hbuf), "(10 (0 %zx %zx 0))\n", first_node_index, npoints);
    fh << hbuf;

    std::size_t total_cells = 0;
    for (const auto& cb : mesh.cells) total_cells += cb.num_cells();
    std::snprintf(hbuf, sizeof(hbuf), "(12 (0 1 %zx 0))\n", total_cells);
    fh << hbuf;

    // Nodes
    const char* nkey = binary ? "3010" : "10";
    std::snprintf(hbuf, sizeof(hbuf), "(%s (1 %zx %zx 1 %zx)(\n", nkey,
                  first_node_index, npoints, dim);
    fh << hbuf;
    if (binary) {
        for (std::size_t i = 0; i < npoints; ++i)
            for (std::size_t c = 0; c < dim; ++c) {
                double v = detail::read_double(mesh.points, i * dim + c);
                fh.write(reinterpret_cast<const char*>(&v), 8);
            }
        fh << "\n)";
        fh << "End of Binary Section 3010)\n";
    } else {
        char cbuf[32];
        for (std::size_t i = 0; i < npoints; ++i) {
            for (std::size_t c = 0; c < dim; ++c) {
                std::snprintf(cbuf, sizeof(cbuf), "%.16e", detail::read_double(mesh.points, i * dim + c));
                fh << cbuf << (c + 1 == dim ? "" : " ");
            }
            fh << "\n";
        }
        fh << "))\n";
    }

    // Cells
    std::size_t first_index = 0;
    for (const auto& cb : mesh.cells) {
        auto it = meshio_to_ansys.find(cb.type);
        if (it == meshio_to_ansys.end())
            throw WriteError("ANSYS: illegal cell type '" + cb.type + "'");
        int ansys_type = it->second;
        std::size_t n = cb.num_cells();
        std::size_t ncols = detail::cols(cb.data);
        std::size_t last_index = first_index + n - 1;
        bool is_i32 = (cb.data.dtype() == DType::Int32);
        const char* ckey = binary ? (is_i32 ? "2012" : "3012") : "12";
        std::snprintf(hbuf, sizeof(hbuf), "(%s (1 %zx %zx 1 %d)(\n", ckey, first_index,
                      last_index, ansys_type);
        fh << hbuf;
        if (binary) {
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t c = 0; c < ncols; ++c) {
                    std::int64_t v = detail::read_int(cb.data, r * ncols + c) + 1;
                    if (is_i32) { std::int32_t v32 = static_cast<std::int32_t>(v); fh.write(reinterpret_cast<const char*>(&v32), 4); }
                    else fh.write(reinterpret_cast<const char*>(&v), 8);
                }
            fh << "\n)";
            std::snprintf(hbuf, sizeof(hbuf), "End of Binary Section %s)\n", ckey);
            fh << hbuf;
        } else {
            char cbuf[24];
            for (std::size_t r = 0; r < n; ++r) {
                for (std::size_t c = 0; c < ncols; ++c) {
                    std::snprintf(cbuf, sizeof(cbuf), "%llx",
                                  static_cast<unsigned long long>(detail::read_int(cb.data, r * ncols + c) + 1));
                    fh << cbuf << (c + 1 == ncols ? "" : " ");
                }
                fh << "\n";
            }
            fh << "))\n";
        }
        first_index = last_index + 1;
    }
}

}  // namespace meshio
