#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"

namespace meshioplusplus {

namespace {

std::string strip(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

}  // namespace

Mesh read_off(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw ReadError("Could not open file: " + path);

    std::string line;
    if (!std::getline(in, line) || strip(line) != "OFF")
        throw ReadError("Expected the first line to be 'OFF'");

    // Skip comments / blank lines to the counts line.
    std::string counts;
    while (std::getline(in, line)) {
        std::string s = strip(line);
        if (!s.empty() && s[0] != '#') {
            counts = s;
            break;
        }
    }
    std::istringstream cs(counts);
    long long num_verts = 0, num_faces = 0, num_edges = 0;
    cs >> num_verts >> num_faces >> num_edges;

    Mesh mesh;
    mesh.points = NDArray(DType::Float64, {static_cast<std::size_t>(num_verts), 3});
    double* pp = mesh.points.as<double>();
    for (long long i = 0; i < num_verts * 3; ++i) {
        if (!(in >> pp[i])) throw ReadError("OFF: not enough vertex coordinates");
    }

    NDArray cells(DType::Int64, {static_cast<std::size_t>(num_faces), 3});
    std::int64_t* cp = cells.as<std::int64_t>();
    for (long long f = 0; f < num_faces; ++f) {
        long long n;
        if (!(in >> n)) throw ReadError("OFF: not enough faces");
        if (n != 3) throw ReadError("OFF: can only read triangular faces");
        in >> cp[f * 3 + 0] >> cp[f * 3 + 1] >> cp[f * 3 + 2];
    }
    mesh.cells.emplace_back("triangle", std::move(cells));
    return mesh;
}

void write_off(const std::string& path, const Mesh& mesh) {
    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    // Gather triangles (OFF supports triangles only).
    std::vector<std::int64_t> tri;
    std::size_t ntri = 0;
    for (const auto& cb : mesh.cells) {
        if (cb.type != "triangle") continue;
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            for (int k = 0; k < 3; ++k) tri.push_back(detail::read_int(cb.data, r * 3 + k));
            ++ntri;
        }
    }

    os << "OFF\n# Created by meshio++ (C++ core)\n\n";
    os << num_points << ' ' << ntri << " 0\n\n";

    char buf[96];
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? detail::read_double(mesh.points, r * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(mesh.points, r * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(mesh.points, r * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "%.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }
    for (std::size_t t = 0; t < ntri; ++t)
        os << "3 " << tri[t * 3] << ' ' << tri[t * 3 + 1] << ' ' << tri[t * 3 + 2] << '\n';
}

}  // namespace meshioplusplus
