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
#include <vector>

// Project includes
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

std::string comsol_to_meshio(const std::string& t) {
    static const std::unordered_map<std::string, std::string> m = {
        {"vtx", "vertex"}, {"edg", "line"}, {"tri", "triangle"}, {"quad", "quad"},
        {"tet", "tetra"}, {"prism", "wedge"}, {"pyr", "pyramid"}, {"hex", "hexahedron"},
        {"edg2", "line3"}, {"tri2", "triangle6"}, {"quad2", "quad9"},
        {"tet2", "tetra10"}, {"prism2", "wedge18"}, {"hex2", "hexahedron27"}};
    auto it = m.find(t);
    return it == m.end() ? std::string() : it->second;
}

std::string meshio_to_comsol(const std::string& t) {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "vtx"}, {"line", "edg"}, {"triangle", "tri"}, {"quad", "quad"},
        {"tetra", "tet"}, {"wedge", "prism"}, {"pyramid", "pyr"}, {"hexahedron", "hex"},
        {"line3", "edg2"}, {"triangle6", "tri2"}, {"quad9", "quad2"},
        {"tetra10", "tet2"}, {"wedge18", "prism2"}, {"hexahedron27", "hex2"}};
    auto it = m.find(t);
    return it == m.end() ? std::string() : it->second;
}

const std::vector<int>* perm_of(const std::string& t) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"quad", {0, 1, 3, 2}}, {"hexahedron", {0, 1, 3, 2, 4, 5, 7, 6}}};
    auto it = m.find(t);
    return it == m.end() ? nullptr : &it->second;
}

struct Cursor {
    std::vector<std::string> t;
    std::size_t i = 0;
    const std::string& tok() {
        if (i >= t.size()) throw ReadError("mphtxt: unexpected end of file");
        return t[i++];
    }
    long long integer() { return std::strtoll(tok().c_str(), nullptr, 10); }
    double real() { return std::strtod(tok().c_str(), nullptr); }
    std::string str() {
        integer();  // length prefix
        return tok();
    }
};

}  // namespace

Mesh read_mphtxt(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    Cursor c;
    std::string line;
    while (std::getline(in, line)) {
        std::size_t h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        std::istringstream iss(line);
        std::string w;
        while (iss >> w) c.t.push_back(w);
    }

    c.integer();  // version major
    c.integer();  // version minor
    for (long long k = c.integer(); k > 0; --k) c.str();  // tags
    const long long n_types = c.integer();
    for (long long k = 0; k < n_types; ++k) c.str();  // type names

    Mesh mesh;
    std::vector<NDArray> geom;

    for (long long obj = 0; obj < n_types; ++obj) {
        c.integer(); c.integer(); c.integer();  // object type indices
        c.str();                                 // class name
        c.integer();                             // object version
        const long long sdim = c.integer();
        const long long n_points = c.integer();
        const long long lowest = c.integer();
        mesh.points = NDArray(DType::Float64, {static_cast<std::size_t>(n_points),
                                               static_cast<std::size_t>(sdim)});
        for (long long p = 0; p < n_points * sdim; ++p)
            mesh.points.as<double>()[p] = c.real();

        const long long n_eltypes = c.integer();
        for (long long e = 0; e < n_eltypes; ++e) {
            std::string ctype = c.str();
            std::string mtype = comsol_to_meshio(ctype);
            if (mtype.empty()) throw ReadError("mphtxt: unknown element type " + ctype);
            const long long nn = c.integer();
            const long long ne = c.integer();
            NDArray conn(DType::Int64, {static_cast<std::size_t>(ne),
                                        static_cast<std::size_t>(nn)});
            std::vector<std::int64_t> raw(ne * nn);
            for (long long v = 0; v < ne * nn; ++v) raw[v] = c.integer() - lowest;
            const std::vector<int>* p = perm_of(mtype);
            for (long long r = 0; r < ne; ++r)
                for (long long j = 0; j < nn; ++j)
                    conn.as<std::int64_t>()[r * nn + j] =
                        p ? raw[r * nn + (*p)[j]] : raw[r * nn + j];

            const long long npar_per = c.integer();
            const long long npar = c.integer();
            for (long long v = 0; v < npar * npar_per; ++v) c.tok();
            const long long ngeom = c.integer();
            NDArray g(DType::Int64, {static_cast<std::size_t>(ngeom)});
            for (long long v = 0; v < ngeom; ++v) g.as<std::int64_t>()[v] = c.integer();
            const long long nud = c.integer();
            for (long long v = 0; v < nud * 2; ++v) c.integer();

            mesh.cells.emplace_back(mtype, std::move(conn));
            geom.push_back(std::move(g));
        }
        break;  // first mesh object only
    }

    if (!geom.empty()) mesh.cell_data.emplace("mphtxt:geom", std::move(geom));
    return mesh;
}

void write_mphtxt(const std::string& path, const Mesh& mesh) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw WriteError("Could not open file for writing: " + path);

    const std::size_t sdim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    struct Blk { std::size_t idx; const CellBlock* cb; };
    std::vector<Blk> blocks;
    for (std::size_t k = 0; k < mesh.cells.size(); ++k)
        if (!meshio_to_comsol(mesh.cells[k].type).empty())
            blocks.push_back({k, &mesh.cells[k]});
        else
            throw WriteError("mphtxt: unsupported cell type " + mesh.cells[k].type);

    auto geom_it = mesh.cell_data.find("mphtxt:geom");

    f << "# Created by meshio++ (C++ core)\n\n";
    f << "0 1\n";
    f << "1 # number of tags\n5 mesh1\n";
    f << "1 # number of types\n3 obj\n\n";
    f << "0 0 1\n4 Mesh # class\n2 # version\n";
    f << sdim << " # sdim\n";
    f << mesh.num_points() << " # number of mesh points\n";
    f << "1 # lowest mesh point index\n\n# Mesh point coordinates\n";
    char buf[32];
    for (std::size_t i = 0; i < mesh.num_points(); ++i) {
        for (std::size_t cc = 0; cc < sdim; ++cc) {
            std::snprintf(buf, sizeof(buf), "%.16g", detail::read_double(mesh.points, i * sdim + cc));
            f << buf << (cc + 1 == sdim ? '\n' : ' ');
        }
    }
    f << "\n" << blocks.size() << " # number of element types\n\n";

    int ti = 0;
    for (const auto& b : blocks) {
        const CellBlock& cb = *b.cb;
        std::string ctype = meshio_to_comsol(cb.type);
        const std::vector<int>* p = perm_of(cb.type);
        std::size_t nn = detail::cols(cb.data);
        std::size_t ne = cb.num_cells();
        f << "# Type #" << (++ti) << "\n\n";
        f << ctype.size() << " " << ctype << " # type name\n\n";
        f << nn << " # number of nodes per element\n";
        f << ne << " # number of elements\n# Elements\n";
        for (std::size_t r = 0; r < ne; ++r) {
            for (std::size_t j = 0; j < nn; ++j) {
                std::size_t src = p ? (*p)[j] : j;
                f << (detail::read_int(cb.data, r * nn + src) + 1) << (j + 1 == nn ? '\n' : ' ');
            }
        }
        f << "\n" << nn << " # number of parameter values per element\n";
        f << "0 # number of parameters\n# Parameters\n\n";
        f << ne << " # number of geometric entity indices\n# Geometric entity indices\n";
        const NDArray* g = (geom_it != mesh.cell_data.end() && b.idx < geom_it->second.size())
                               ? &geom_it->second[b.idx]
                               : nullptr;
        for (std::size_t r = 0; r < ne; ++r)
            f << (g ? detail::read_int(*g, r) : 0) << "\n";
        f << "\n0 # number of up/down pairs\n# Up/down\n\n";
    }
}

}  // namespace meshioplusplus
