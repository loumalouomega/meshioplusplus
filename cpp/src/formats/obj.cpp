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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"

namespace meshioplusplus {

namespace {

struct FaceBlock {
    std::size_t size = 0;
    std::vector<std::int64_t> idx;  // flat, 0-based
    std::vector<std::int64_t> gids;
    std::size_t count = 0;
};

std::string cell_type_for(std::size_t n) {
    if (n == 3) return "triangle";
    if (n == 4) return "quad";
    return "polygon";
}

NDArray make_point_data(const std::vector<std::vector<double>>& rows) {
    std::size_t n = rows.size();
    std::size_t nc = n ? rows[0].size() : 0;
    NDArray a(DType::Float64, {n, nc});
    double* p = a.as<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < nc; ++j) p[i * nc + j] = rows[i][j];
    return a;
}

}  // namespace

Mesh read_obj(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw ReadError("Could not open file: " + path);

    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<double>> vn, vt;
    std::vector<FaceBlock> blocks;
    std::int64_t group_id = -1;

    std::string line;
    while (std::getline(in, line)) {
        // strip
        std::size_t b = 0, e = line.size();
        while (b < e && std::isspace(static_cast<unsigned char>(line[b]))) ++b;
        while (e > b && std::isspace(static_cast<unsigned char>(line[e - 1]))) --e;
        if (b == e || line[b] == '#') continue;

        std::istringstream iss(line.substr(b, e - b));
        std::string tag;
        iss >> tag;
        if (tag == "v") {
            std::array<double, 3> p{0, 0, 0};
            iss >> p[0] >> p[1] >> p[2];
            points.push_back(p);
        } else if (tag == "vn") {
            std::vector<double> row;
            double x;
            while (iss >> x) row.push_back(x);
            vn.push_back(row);
        } else if (tag == "vt") {
            std::vector<double> row;
            double x;
            while (iss >> x) row.push_back(x);
            vt.push_back(row);
        } else if (tag == "f") {
            std::vector<std::int64_t> dat;
            std::string item;
            while (iss >> item) {
                std::size_t slash = item.find('/');
                std::string num = (slash == std::string::npos) ? item : item.substr(0, slash);
                dat.push_back(static_cast<std::int64_t>(std::stoll(num)) - 1);
            }
            std::size_t sz = dat.size();
            if (blocks.empty() || (blocks.back().count > 0 && blocks.back().size != sz)) {
                FaceBlock fb;
                fb.size = sz;
                blocks.push_back(std::move(fb));
            }
            FaceBlock& cur = blocks.back();
            if (cur.count == 0) cur.size = sz;
            cur.idx.insert(cur.idx.end(), dat.begin(), dat.end());
            cur.gids.push_back(group_id);
            ++cur.count;
        } else if (tag == "g") {
            FaceBlock fb;
            blocks.push_back(std::move(fb));
            ++group_id;
        }
        // 's' and others: ignored.
    }

    // Drop empty blocks (e.g. from trailing 'g').
    std::vector<FaceBlock> nonempty;
    for (auto& fb : blocks)
        if (fb.count > 0) nonempty.push_back(std::move(fb));

    Mesh mesh;
    std::size_t np = points.size();
    mesh.points = NDArray(DType::Float64, {np, 3});
    double* pp = mesh.points.as<double>();
    for (std::size_t i = 0; i < np; ++i)
        for (int c = 0; c < 3; ++c) pp[i * 3 + c] = points[i][c];

    if (!vt.empty()) mesh.point_data.emplace("obj:vt", make_point_data(vt));
    if (!vn.empty()) mesh.point_data.emplace("obj:vn", make_point_data(vn));

    if (!nonempty.empty()) {
        std::vector<NDArray> gid_blocks;
        for (auto& fb : nonempty) {
            NDArray data(DType::Int64, {fb.count, fb.size});
            std::int64_t* dp = data.as<std::int64_t>();
            for (std::size_t i = 0; i < fb.idx.size(); ++i) dp[i] = fb.idx[i];
            mesh.cells.emplace_back(cell_type_for(fb.size), std::move(data));

            NDArray g(DType::Int64, {fb.count});
            std::int64_t* gp = g.as<std::int64_t>();
            for (std::size_t i = 0; i < fb.count; ++i) gp[i] = fb.gids[i];
            gid_blocks.push_back(std::move(g));
        }
        mesh.cell_data.emplace("obj:group_ids", std::move(gid_blocks));
    }
    return mesh;
}

void write_obj(const std::string& path, const Mesh& mesh) {
    for (const auto& cb : mesh.cells)
        if (cb.type != "triangle" && cb.type != "quad" && cb.type != "polygon")
            throw WriteError("Wavefront .obj files can only contain triangle, quad, "
                             "or polygon cells.");

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    os << "# Created by meshio++ (C++ core)\n";
    char buf[96];
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? detail::read_double(mesh.points, r * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(mesh.points, r * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(mesh.points, r * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "v %.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }

    auto write_pd = [&](const char* key, const char* tag) {
        auto it = mesh.point_data.find(key);
        if (it == mesh.point_data.end()) return;
        const NDArray& d = it->second;
        std::size_t nc = d.shape().size() >= 2 ? d.shape()[1] : 1;
        for (std::size_t r = 0; r < (d.shape().empty() ? 0 : d.shape()[0]); ++r) {
            os << tag;
            for (std::size_t c = 0; c < nc; ++c) {
                std::snprintf(buf, sizeof(buf), " %.17g", detail::read_double(d, r * nc + c));
                os << buf;
            }
            os << '\n';
        }
    };
    write_pd("obj:vn", "vn");
    write_pd("obj:vt", "vt");

    for (const auto& cb : mesh.cells) {
        std::size_t k = cb.data.shape().size() >= 2 ? cb.data.shape()[1] : 1;
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            os << 'f';
            for (std::size_t j = 0; j < k; ++j)
                os << ' ' << (detail::read_int(cb.data, r * k + j) + 1);
            os << '\n';
        }
    }
}

}  // namespace meshioplusplus
