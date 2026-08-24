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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

std::string off_strip(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

std::string off_cell_type_from_count(std::size_t n) {
    switch (n) {
        case 3:
            return "triangle";
        case 4:
            return "quad";
        default:
            return "polygon";
    }
}

}  // namespace

Mesh read_off(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    std::string line;
    if (!std::getline(in, line) || off_strip(line) != "OFF")
        throw ReadError("Expected the first line to be 'OFF'");

    // Skip comments / blank lines to the counts line.
    std::string counts;
    while (std::getline(in, line)) {
        std::string s = off_strip(line);
        if (!s.empty() && s[0] != '#') {
            counts = s;
            break;
        }
    }
    std::istringstream cs(counts);
    long long num_verts = 0, num_faces = 0, num_edges = 0;
    cs >> num_verts >> num_faces >> num_edges;

    Mesh mesh;
    NDArray pts(DType::Float64, {static_cast<std::size_t>(num_verts), 3});
    double* pp = pts.As<double>();
    for (long long i = 0; i < num_verts * 3; ++i) {
        if (!(in >> pp[i]))
            throw ReadError("OFF: not enough vertex coordinates");
    }
    mesh.AssignPoints(std::move(pts));

    if (num_faces == 0) {
        mesh.AddCellBlock("triangle", NDArray(DType::Int64, {0, 3}));
        return mesh;
    }

    // Faces are grouped into blocks by vertex count (a run of same-count faces
    // stays in one block, matching the OBJ reader's approach in this same header):
    // 3 -> triangle, 4 -> quad, anything else -> polygon.
    std::size_t cur_n = 0;
    std::vector<std::int64_t> cur_conn;
    std::size_t cur_count = 0;
    auto flush = [&]() {
        if (cur_count == 0)
            return;
        NDArray data(DType::Int64, {cur_count, cur_n});
        std::memcpy(data.Data(), cur_conn.data(), cur_conn.size() * sizeof(std::int64_t));
        mesh.AddCellBlock(off_cell_type_from_count(cur_n), std::move(data));
        cur_conn.clear();
        cur_count = 0;
    };
    for (long long f = 0; f < num_faces; ++f) {
        long long n;
        if (!(in >> n))
            throw ReadError("OFF: not enough faces");
        if (n < 3)
            throw ReadError("OFF: faces must have at least 3 vertices");
        const std::size_t un = static_cast<std::size_t>(n);
        if (un != cur_n) {
            flush();
            cur_n = un;
        }
        for (std::size_t k = 0; k < un; ++k) {
            std::int64_t idx;
            if (!(in >> idx))
                throw ReadError("OFF: not enough face vertex indices");
            cur_conn.push_back(idx);
        }
        ++cur_count;
    }
    flush();
    return mesh;
}

void write_off(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // OFF represents polygonal faces (triangle/quad/polygon); anything else is
    // skipped with a warning, matching the Python reference writer.
    std::size_t num_faces = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& t = cb.Type();
        if (t != "triangle" && t != "quad" && t != "polygon") {
            log::warn(
                "OFF: '{}' cells are not representable (only triangle/quad/polygon "
                "faces); skipping.",
                t);
            continue;
        }
        num_faces += cb.NumCells();
    }

    os << "OFF\n# " << detail::kProvenanceTag << "\n\n";
    os << num_points << ' ' << num_faces << " 0\n\n";

    char buf[96];
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? detail::read_double(points, r * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(points, r * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(points, r * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "%.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }
    for (const auto cb : rMesh.CellRange()) {
        const std::string& t = cb.Type();
        if (t != "triangle" && t != "quad" && t != "polygon")
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << k;
            for (std::size_t j = 0; j < k; ++j)
                os << ' ' << detail::read_int(conn, r * k + j);
            os << '\n';
        }
    }
}

}  // namespace meshioplusplus
