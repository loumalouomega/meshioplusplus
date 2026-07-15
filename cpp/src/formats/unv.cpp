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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// Salome/UNV parabolic node order -> meshio position (0-based).
const std::vector<int>* nd_perm(const std::string& t) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"line3", {0, 2, 1}},
        {"triangle6", {0, 3, 1, 4, 2, 5}},
        {"quad8", {0, 4, 1, 5, 2, 6, 3, 7}},
        {"tetra10", {0, 4, 1, 5, 2, 6, 7, 8, 9, 3}},
        {"hexahedron20",
         {0, 8, 1, 9, 2, 10, 3, 11, 12, 13, 14, 15, 4, 16, 5, 17, 6, 18, 7, 19}}};
    auto it = m.find(t);
    return it == m.end() ? nullptr : &it->second;
}

std::string unv_type(int fedesc) {
    static const std::unordered_map<int, std::string> m = {
        {11, "line"}, {21, "line"}, {22, "line3"}, {24, "line3"},
        {41, "triangle"}, {81, "triangle"}, {91, "triangle"},
        {42, "triangle6"}, {82, "triangle6"}, {92, "triangle6"},
        {44, "quad"}, {84, "quad"}, {94, "quad"}, {122, "quad"},
        {45, "quad8"}, {85, "quad8"}, {95, "quad8"},
        {111, "tetra"}, {118, "tetra10"}, {112, "wedge"},
        {115, "hexahedron"}, {116, "hexahedron20"}};
    auto it = m.find(fedesc);
    return it == m.end() ? std::string() : it->second;
}

// meshio type -> (descriptor, is_beam)
bool meshio_descriptor(const std::string& t, int& desc, bool& beam) {
    static const std::unordered_map<std::string, std::pair<int, bool>> m = {
        {"line", {21, true}}, {"line3", {24, true}},
        {"triangle", {91, false}}, {"triangle6", {92, false}},
        {"quad", {94, false}}, {"quad8", {95, false}},
        {"tetra", {111, false}}, {"tetra10", {118, false}},
        {"wedge", {112, false}}, {"hexahedron", {115, false}},
        {"hexahedron20", {116, false}}};
    auto it = m.find(t);
    if (it == m.end()) return false;
    desc = it->second.first;
    beam = it->second.second;
    return true;
}

bool is_beam(int fedesc) {
    return fedesc == 11 || fedesc == 21 || fedesc == 22 || fedesc == 24;
}

std::vector<std::string> tokens(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

double parse_coord(std::string s) {
    for (char& c : s)
        if (c == 'D' || c == 'd') c = 'E';
    return std::strtod(s.c_str(), nullptr);
}

}  // namespace

Mesh read_unv(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);

    std::vector<std::vector<double>> points;
    std::unordered_map<std::int64_t, std::int64_t> label_to_index;

    struct Group {
        std::string type;
        std::vector<std::vector<std::int64_t>> rows;
        std::vector<std::int64_t> pid;
    };
    std::vector<Group> groups;
    std::unordered_map<std::string, std::size_t> group_index;
    std::size_t dim = 3;

    std::size_t i = 0, n = lines.size();
    auto strip = [](const std::string& s) {
        std::size_t a = s.find_first_not_of(" \t\r");
        std::size_t b = s.find_last_not_of(" \t\r");
        return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
    };

    while (i < n) {
        if (strip(lines[i]) != "-1") { ++i; continue; }
        ++i;
        if (i >= n) break;
        int ds = std::atoi(strip(lines[i]).c_str());
        ++i;
        std::size_t start = i;
        while (i < n && strip(lines[i]) != "-1") ++i;
        std::size_t end = i;  // exclusive
        ++i;                  // skip closing -1

        if (ds == 2411 || ds == 781) {
            std::size_t k = start;
            while (k + 1 < end) {
                auto r1 = tokens(lines[k]);
                if (r1.empty()) { ++k; continue; }
                std::int64_t label = std::strtoll(r1[0].c_str(), nullptr, 10);
                auto co = tokens(lines[k + 1]);
                std::vector<double> p;
                for (const auto& c : co) p.push_back(parse_coord(c));
                if (points.empty() && !p.empty()) dim = p.size();
                label_to_index[label] = static_cast<std::int64_t>(points.size());
                points.push_back(std::move(p));
                k += 2;
            }
        } else if (ds == 2412) {
            std::size_t k = start;
            while (k < end) {
                auto r1 = tokens(lines[k]);
                if (r1.size() < 6) break;
                int fedesc = std::atoi(r1[1].c_str());
                std::int64_t pid = std::strtoll(r1[2].c_str(), nullptr, 10);
                int num_nodes = std::atoi(r1[5].c_str());
                ++k;
                if (is_beam(fedesc)) ++k;  // skip orientation
                std::vector<std::int64_t> nl;
                while (static_cast<int>(nl.size()) < num_nodes && k < end) {
                    for (const auto& v : tokens(lines[k]))
                        nl.push_back(std::strtoll(v.c_str(), nullptr, 10));
                    ++k;
                }
                nl.resize(num_nodes);
                std::string mtype = unv_type(fedesc);
                if (mtype.empty())
                    throw ReadError("UNV: FE descriptor not supported");
                std::vector<std::int64_t> unv_conn(num_nodes);
                for (int j = 0; j < num_nodes; ++j)
                    unv_conn[j] = label_to_index.at(nl[j]);
                const std::vector<int>* nd = nd_perm(mtype);
                std::vector<std::int64_t> conn(num_nodes);
                if (nd)
                    for (int j = 0; j < num_nodes; ++j) conn[(*nd)[j]] = unv_conn[j];
                else
                    conn = unv_conn;

                auto git = group_index.find(mtype);
                if (git == group_index.end()) {
                    group_index[mtype] = groups.size();
                    groups.push_back({mtype, {}, {}});
                    git = group_index.find(mtype);
                }
                groups[git->second].rows.push_back(std::move(conn));
                groups[git->second].pid.push_back(pid);
            }
        } else if (ds == 2467 || ds == 2477 || ds == 2414 || ds == 55 || ds == 56 ||
                   ds == 57) {
            // groups / fields need the Python reference implementation
            throw ReadError("UNV: dataset " + std::to_string(ds) +
                            " handled by Python fallback");
        }
        // other datasets ignored
    }

    Mesh mesh;
    const std::size_t np = points.size();
    mesh.points = NDArray(DType::Float64, {np, dim});
    for (std::size_t r = 0; r < np; ++r)
        for (std::size_t c = 0; c < dim && c < points[r].size(); ++c)
            mesh.points.as<double>()[r * dim + c] = points[r][c];

    std::vector<NDArray> pids;
    for (auto& g : groups) {
        std::size_t ne = g.rows.size();
        std::size_t k = ne ? g.rows[0].size() : 0;
        NDArray data(DType::Int64, {ne, k});
        for (std::size_t r = 0; r < ne; ++r)
            for (std::size_t j = 0; j < k; ++j)
                data.as<std::int64_t>()[r * k + j] = g.rows[r][j];
        mesh.cells.emplace_back(g.type, std::move(data));
        NDArray pd(DType::Int64, {ne});
        for (std::size_t r = 0; r < ne; ++r) pd.as<std::int64_t>()[r] = g.pid[r];
        pids.push_back(std::move(pd));
    }
    if (!pids.empty()) mesh.cell_data.emplace("unv:pid", std::move(pids));
    return mesh;
}

void write_unv(const std::string& path, const Mesh& mesh) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw WriteError("Could not open file for writing: " + path);

    const std::size_t np = mesh.num_points();
    const std::size_t pdim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    // 2411 nodes
    f << "    -1\n  2411\n";
    char buf[128];
    for (std::size_t k = 0; k < np; ++k) {
        std::snprintf(buf, sizeof(buf), "%10zu%10d%10d%10d\n", k + 1, 1, 1, 11);
        f << buf;
        for (int c = 0; c < 3; ++c) {
            double v = c < static_cast<int>(pdim)
                           ? detail::read_double(mesh.points, k * pdim + c)
                           : 0.0;
            std::snprintf(buf, sizeof(buf), "%25.16E", v);
            f << buf;
        }
        f << "\n";
    }
    f << "    -1\n";

    // 2412 elements
    f << "    -1\n  2412\n";
    auto pid_it = mesh.cell_data.find("unv:pid");
    std::int64_t label = 0;
    for (std::size_t bi = 0; bi < mesh.cells.size(); ++bi) {
        const CellBlock& cb = mesh.cells[bi];
        int desc;
        bool beam;
        if (!meshio_descriptor(cb.type, desc, beam))
            throw WriteError("UNV: unsupported cell type " + cb.type);
        const std::vector<int>* nd = nd_perm(cb.type);
        std::size_t ncols = detail::cols(cb.data);
        std::size_t nrows = cb.num_cells();
        const NDArray* pid = (pid_it != mesh.cell_data.end() && bi < pid_it->second.size())
                                 ? &pid_it->second[bi]
                                 : nullptr;
        for (std::size_t r = 0; r < nrows; ++r) {
            ++label;
            std::int64_t pval = pid ? detail::read_int(*pid, r) : 1;
            std::snprintf(buf, sizeof(buf), "%10lld%10d%10lld%10lld%10d%10zu\n",
                          static_cast<long long>(label), desc,
                          static_cast<long long>(pval), static_cast<long long>(pval),
                          11, ncols);
            f << buf;
            if (beam) f << "         0         0         0\n";
            // reorder meshio -> UNV, 1-based, 8 per line
            std::vector<std::int64_t> unv(ncols);
            for (std::size_t j = 0; j < ncols; ++j) {
                std::int64_t node = detail::read_int(cb.data, r * ncols + j) + 1;
                if (nd)
                    unv[j] = 0;  // filled below
                else
                    unv[j] = node;
            }
            if (nd)
                for (std::size_t j = 0; j < ncols; ++j)
                    unv[j] = detail::read_int(cb.data, r * ncols + (*nd)[j]) + 1;
            for (std::size_t j = 0; j < ncols; ++j) {
                std::snprintf(buf, sizeof(buf), "%10lld",
                              static_cast<long long>(unv[j]));
                f << buf;
                if ((j + 1) % 8 == 0 || j + 1 == ncols) f << "\n";
            }
        }
    }
    f << "    -1\n";
}

}  // namespace meshioplusplus
