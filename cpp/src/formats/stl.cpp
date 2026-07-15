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
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// First-occurrence de-duplication of 3-component rows. Returns per-row unique
// index; appends unique rows (raw bytes) to `out_points`.
std::vector<std::int64_t> dedup(const unsigned char* rows, std::size_t nrows,
                                std::size_t isz, std::vector<unsigned char>& out_points) {
    std::unordered_map<std::string, std::int64_t> seen;
    seen.reserve(nrows);
    std::vector<std::int64_t> idx(nrows);
    const std::size_t rowbytes = 3 * isz;
    for (std::size_t i = 0; i < nrows; ++i) {
        std::string key(reinterpret_cast<const char*>(rows) + i * rowbytes, rowbytes);
        auto it = seen.find(key);
        if (it == seen.end()) {
            std::int64_t id = static_cast<std::int64_t>(seen.size());
            seen.emplace(std::move(key), id);
            out_points.insert(out_points.end(), rows + i * rowbytes,
                              rows + (i + 1) * rowbytes);
            idx[i] = id;
        } else {
            idx[i] = it->second;
        }
    }
    return idx;
}

Mesh build_mesh(std::vector<unsigned char>& vert_bytes, DType dt,
                std::vector<unsigned char>* normal_bytes) {
    std::size_t isz = dtype_size(dt);
    std::size_t nverts = vert_bytes.size() / (3 * isz);

    std::vector<unsigned char> point_bytes;
    std::vector<std::int64_t> idx = dedup(vert_bytes.data(), nverts, isz, point_bytes);

    Mesh mesh;
    std::size_t num_unique = point_bytes.size() / (3 * isz);
    mesh.points = NDArray(dt, {num_unique, 3});
    if (!point_bytes.empty())
        std::memcpy(mesh.points.data(), point_bytes.data(), point_bytes.size());

    std::size_t ntri = nverts / 3;
    // An empty STL has no cells (match the Python reader, which returns no
    // cell blocks rather than an empty triangle block).
    if (ntri == 0) return mesh;

    NDArray cells(DType::Int64, {ntri, 3});
    std::int64_t* cp = cells.as<std::int64_t>();
    for (std::size_t i = 0; i < ntri * 3; ++i) cp[i] = idx[i];
    mesh.cells.emplace_back("triangle", std::move(cells));

    if (normal_bytes && !normal_bytes->empty()) {
        NDArray nrm(DType::Float64, {ntri, 3});
        std::memcpy(nrm.data(), normal_bytes->data(), normal_bytes->size());
        mesh.cell_data["facet_normals"].push_back(std::move(nrm));
    }
    return mesh;
}

bool starts_with(const std::string& s, const char* p) {
    return s.rfind(p, 0) == 0;
}

bool is_comment_line(const std::string& s) {
    return starts_with(s, "solid") || starts_with(s, "outer loop") ||
           starts_with(s, "endloop") || starts_with(s, "endfacet") ||
           starts_with(s, "endsolid");
}

std::string lstrip(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    return s.substr(b);
}

Mesh read_ascii(std::ifstream& in) {
    // Collect the last 3 numbers of every non-comment line; rows 0,4,8,... are
    // facet normals, the rest are vertices.
    std::vector<double> data;
    std::string line;
    while (std::getline(in, line)) {
        std::string s = lstrip(line);
        if (s.empty() || is_comment_line(s)) continue;
        std::istringstream iss(s);
        std::vector<std::string> tok;
        std::string t;
        while (iss >> t) tok.push_back(t);
        if (tok.size() < 3) continue;
        for (std::size_t j = tok.size() - 3; j < tok.size(); ++j)
            data.push_back(std::strtod(tok[j].c_str(), nullptr));
    }
    std::size_t nrows = data.size() / 3;
    if (nrows % 4 != 0) throw ReadError("Malformed ascii STL");

    std::vector<unsigned char> verts, normals;
    for (std::size_t r = 0; r < nrows; ++r) {
        const double* row = data.data() + r * 3;
        std::vector<unsigned char>& dst = (r % 4 == 0) ? normals : verts;
        dst.insert(dst.end(), reinterpret_cast<const unsigned char*>(row),
                   reinterpret_cast<const unsigned char*>(row) + 3 * sizeof(double));
    }
    return build_mesh(verts, DType::Float64, &normals);
}

Mesh read_binary(std::ifstream& in, std::uint32_t num_tri) {
    std::vector<unsigned char> verts;
    verts.reserve(num_tri * 9 * sizeof(float));
    unsigned char tri[50];
    for (std::uint32_t i = 0; i < num_tri; ++i) {
        in.read(reinterpret_cast<char*>(tri), 50);
        if (in.gcount() != 50) throw ReadError("Truncated binary STL");
        // bytes [12, 48) are the 9 float32 vertex coords (host is little-endian).
        verts.insert(verts.end(), tri + 12, tri + 48);
    }
    return build_mesh(verts, DType::Float32, nullptr);
}

}  // namespace

Mesh read_stl(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    in.seekg(0, std::ios::end);
    std::streamoff filesize = in.tellg();
    in.seekg(0, std::ios::beg);

    if (filesize < 80) return read_ascii(in);

    char header[80];
    in.read(header, 80);
    std::uint32_t num_tri = 0;
    in.read(reinterpret_cast<char*>(&num_tri), 4);  // little-endian host
    if (static_cast<std::streamoff>(84 + std::uint64_t(num_tri) * 50) == filesize)
        return read_binary(in, num_tri);

    // Fall back to ascii: rewind, skip the first line.
    in.clear();
    in.seekg(0, std::ios::beg);
    std::string first;
    std::getline(in, first);
    return read_ascii(in);
}

namespace {

void gather_triangles(const Mesh& mesh, std::vector<std::array<double, 9>>& tris,
                      std::vector<std::array<double, 3>>& normals) {
    const auto& cd = mesh.cell_data.find("facet_normals");
    bool have_normals = cd != mesh.cell_data.end();
    std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 3;

    std::size_t block = 0;
    for (const auto& cb : mesh.cells) {
        if (cb.type != "triangle") {
            ++block;
            continue;
        }
        std::size_t nc = cb.num_cells();
        const NDArray* nrm = nullptr;
        if (have_normals && block < cd->second.size()) nrm = &cd->second[block];
        for (std::size_t r = 0; r < nc; ++r) {
            std::array<double, 9> tri{};
            double v[3][3];
            for (int k = 0; k < 3; ++k) {
                std::int64_t pi = detail::read_int(cb.data, r * 3 + k);
                for (int c = 0; c < 3; ++c)
                    v[k][c] = (std::size_t(c) < dim)
                                  ? detail::read_double(mesh.points, pi * dim + c)
                                  : 0.0;
                tri[k * 3 + 0] = v[k][0];
                tri[k * 3 + 1] = v[k][1];
                tri[k * 3 + 2] = v[k][2];
            }
            tris.push_back(tri);

            std::array<double, 3> n{};
            if (nrm) {
                for (int c = 0; c < 3; ++c) n[c] = detail::read_double(*nrm, r * 3 + c);
            } else {
                double a[3] = {v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]};
                double b[3] = {v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]};
                n[0] = a[1] * b[2] - a[2] * b[1];
                n[1] = a[2] * b[0] - a[0] * b[2];
                n[2] = a[0] * b[1] - a[1] * b[0];
                double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                if (len > 0) {
                    n[0] /= len;
                    n[1] /= len;
                    n[2] /= len;
                }
            }
            normals.push_back(n);
        }
        ++block;
    }
}

}  // namespace

void write_stl(const std::string& path, const Mesh& mesh, bool binary) {
    std::vector<std::array<double, 9>> tris;
    std::vector<std::array<double, 3>> normals;
    gather_triangles(mesh, tris, normals);

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    if (binary) {
        char header[80];
        std::memset(header, 'X', 80);
        const char* msg = "meshio++ (C++ core) binary STL";
        std::memcpy(header, msg, std::strlen(msg));
        os.write(header, 80);
        std::uint32_t n = static_cast<std::uint32_t>(tris.size());
        os.write(reinterpret_cast<const char*>(&n), 4);
        for (std::size_t i = 0; i < tris.size(); ++i) {
            float buf[12];
            for (int c = 0; c < 3; ++c) buf[c] = static_cast<float>(normals[i][c]);
            for (int c = 0; c < 9; ++c) buf[3 + c] = static_cast<float>(tris[i][c]);
            os.write(reinterpret_cast<const char*>(buf), 48);
            std::uint16_t attr = 0;
            os.write(reinterpret_cast<const char*>(&attr), 2);
        }
    } else {
        auto wr3 = [&](const char* prefix, const double* p) {
            char line[160];
            std::snprintf(line, sizeof(line), "%s %.17g %.17g %.17g\n", prefix, p[0],
                          p[1], p[2]);
            os << line;
        };
        os << "solid\n";
        for (std::size_t i = 0; i < tris.size(); ++i) {
            wr3("facet normal", normals[i].data());
            os << " outer loop\n";
            wr3("  vertex", &tris[i][0]);
            wr3("  vertex", &tris[i][3]);
            wr3("  vertex", &tris[i][6]);
            os << " endloop\nendfacet\n";
        }
        os << "endsolid\n";
    }
}

}  // namespace meshioplusplus
