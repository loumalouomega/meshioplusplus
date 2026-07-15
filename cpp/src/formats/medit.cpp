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
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// medit element keyword -> (meshio type, nodes per cell)
const std::unordered_map<std::string, std::pair<std::string, int>>& medit_to_meshio() {
    static const std::unordered_map<std::string, std::pair<std::string, int>> m = {
        {"Edges", {"line", 2}},          {"Triangles", {"triangle", 3}},
        {"Quadrilaterals", {"quad", 4}}, {"Tetrahedra", {"tetra", 4}},
        {"Prisms", {"wedge", 6}},        {"Pyramids", {"pyramid", 5}},
        {"Hexahedra", {"hexahedron", 8}},{"Hexaedra", {"hexahedron", 8}},
    };
    return m;
}

// meshio type -> (medit keyword, nodes per cell), write order.
const std::vector<std::pair<std::string, std::pair<std::string, int>>>&
meshio_to_medit() {
    static const std::vector<std::pair<std::string, std::pair<std::string, int>>> m = {
        {"line", {"Edges", 2}},          {"triangle", {"Triangles", 3}},
        {"quad", {"Quadrilaterals", 4}}, {"tetra", {"Tetrahedra", 4}},
        {"wedge", {"Prisms", 6}},        {"pyramid", {"Pyramids", 5}},
        {"hexahedron", {"Hexahedra", 8}},
    };
    return m;
}

// Whitespace/comment-skipping tokenizer over the whole file.
struct Tokenizer {
    const std::string& buf;
    std::size_t pos = 0;
    explicit Tokenizer(const std::string& b) : buf(b) {}

    bool eof() const { return pos >= buf.size(); }

    void skip_ws() {
        while (pos < buf.size()) {
            char c = buf[pos];
            if (c == '#') {  // comment to end of line
                while (pos < buf.size() && buf[pos] != '\n') ++pos;
            } else if (std::isspace(static_cast<unsigned char>(c))) {
                ++pos;
            } else {
                break;
            }
        }
    }
    std::string next() {
        skip_ws();
        std::size_t start = pos;
        while (pos < buf.size() && !std::isspace(static_cast<unsigned char>(buf[pos])) &&
               buf[pos] != '#')
            ++pos;
        return buf.substr(start, pos - start);
    }
    std::int64_t next_int() { return std::strtoll(next().c_str(), nullptr, 10); }
    double next_double() { return std::strtod(next().c_str(), nullptr); }
    void skip_line() {
        while (pos < buf.size() && buf[pos] != '\n') ++pos;
        if (pos < buf.size()) ++pos;
    }
};

void store_coord(NDArray& a, std::size_t idx, double v) {
    if (a.dtype() == DType::Float32)
        a.as<float>()[idx] = static_cast<float>(v);
    else
        a.as<double>()[idx] = v;
}

const NDArray* pick_first_int(const std::map<std::string, NDArray>& d) {
    for (const auto& kv : d) {
        DType t = kv.second.dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 ||
            t == DType::Int64 || t == DType::UInt8 || t == DType::UInt16 ||
            t == DType::UInt32 || t == DType::UInt64)
            return &kv.second;
    }
    return nullptr;
}

const std::vector<NDArray>* pick_first_int_cell(
    const std::map<std::string, std::vector<NDArray>>& d) {
    for (const auto& kv : d) {
        if (kv.second.empty()) continue;
        DType t = kv.second.front().dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 ||
            t == DType::Int64 || t == DType::UInt8 || t == DType::UInt16 ||
            t == DType::UInt32 || t == DType::UInt64)
            return &kv.second;
    }
    return nullptr;
}

}  // namespace

Mesh read_medit_ascii(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);
    std::string buf((std::istreambuf_iterator<char>(in)),
                    std::istreambuf_iterator<char>());
    Tokenizer tok(buf);

    int dim = 0;
    DType coord_dtype = DType::Float64;
    Mesh mesh;
    std::vector<std::int64_t> point_ref;
    bool have_points = false;

    const auto& e2m = medit_to_meshio();

    while (!tok.eof()) {
        std::string kw = tok.next();
        if (kw.empty()) break;
        if (kw == "MeshVersionFormatted") {
            std::int64_t v = tok.next_int();
            coord_dtype = (v <= 1) ? DType::Float32 : DType::Float64;
        } else if (kw == "Dimension") {
            dim = static_cast<int>(tok.next_int());
        } else if (kw == "Vertices") {
            if (dim <= 0) throw ReadError("Medit: Dimension before Vertices");
            std::int64_t n = tok.next_int();
            mesh.points = NDArray(coord_dtype, {static_cast<std::size_t>(n),
                                                static_cast<std::size_t>(dim)});
            point_ref.resize(n);
            for (std::int64_t i = 0; i < n; ++i) {
                for (int c = 0; c < dim; ++c)
                    store_coord(mesh.points, i * dim + c, tok.next_double());
                point_ref[i] = static_cast<std::int64_t>(tok.next_double());
            }
            have_points = true;
        } else if (e2m.count(kw)) {
            const auto& info = e2m.at(kw);
            const std::string& type = info.first;
            int k = info.second;
            std::int64_t n = tok.next_int();
            NDArray data(DType::Int64, {static_cast<std::size_t>(n),
                                        static_cast<std::size_t>(k)});
            NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
            std::int64_t* dp = data.as<std::int64_t>();
            std::int64_t* rp = ref.as<std::int64_t>();
            for (std::int64_t i = 0; i < n; ++i) {
                for (int j = 0; j < k; ++j) dp[i * k + j] = tok.next_int() - 1;
                rp[i] = tok.next_int();
            }
            mesh.cells.emplace_back(type, std::move(data));
            mesh.cell_data["medit:ref"].push_back(std::move(ref));
        } else if (kw == "Corners") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n; ++i) tok.next();
        } else if (kw == "Normals") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * dim; ++i) tok.next();
        } else if (kw == "NormalAtVertices") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * 2; ++i) tok.next();
        } else if (kw == "SubDomainFromMesh") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * 4; ++i) tok.next();
        } else if (kw == "VertexOnGeometricVertex") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * 2; ++i) tok.next();
        } else if (kw == "VertexOnGeometricEdge") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * 3; ++i) tok.next();
        } else if (kw == "EdgeOnGeometricEdge") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n * 2; ++i) tok.next();
        } else if (kw == "Identifier" || kw == "Geometry") {
            tok.skip_line();
        } else if (kw == "RequiredVertices" || kw == "TangentAtVertices" ||
                   kw == "Tangents" || kw == "Ridges") {
            std::int64_t n = tok.next_int();
            for (std::int64_t i = 0; i < n; ++i) tok.next();
        } else if (kw == "End") {
            break;
        } else {
            throw ReadError("Medit: unknown keyword '" + kw + "'");
        }
    }

    if (!have_points) throw ReadError("Medit: expected Vertices");

    NDArray pr(DType::Int64, {point_ref.size()});
    for (std::size_t i = 0; i < point_ref.size(); ++i) pr.as<std::int64_t>()[i] = point_ref[i];
    mesh.point_data.emplace("medit:ref", std::move(pr));
    return mesh;
}

void write_medit_ascii(const std::string& path, const Mesh& mesh) {
    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t n = mesh.num_points();
    const std::size_t d = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    int version = (mesh.points.dtype() == DType::Float32) ? 1 : 2;

    os << "MeshVersionFormatted " << version << "\n";
    os << "Dimension " << d << "\n";

    // Vertices
    os << "\nVertices\n" << n << "\n";
    const NDArray* vlabels = pick_first_int(mesh.point_data);
    char buf[64];
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t c = 0; c < d; ++c) {
            std::snprintf(buf, sizeof(buf), "%.16e ", detail::read_double(mesh.points, i * d + c));
            os << buf;
        }
        std::int64_t lab = vlabels ? detail::read_int(*vlabels, i) : 1;
        os << lab << "\n";
    }

    // Cells, grouped by medit element keyword.
    const std::vector<NDArray>* clabels = pick_first_int_cell(mesh.cell_data);
    for (const auto& mk : meshio_to_medit()) {
        const std::string& mtype = mk.first;
        const std::string& kw = mk.second.first;
        int k = mk.second.second;
        for (std::size_t ci = 0; ci < mesh.cells.size(); ++ci) {
            const CellBlock& cb = mesh.cells[ci];
            if (cb.type != mtype) continue;
            std::size_t count = cb.num_cells();
            os << "\n" << kw << "\n" << count << "\n";
            const NDArray* lab = (clabels && ci < clabels->size()) ? &(*clabels)[ci] : nullptr;
            for (std::size_t r = 0; r < count; ++r) {
                for (int j = 0; j < k; ++j)
                    os << (detail::read_int(cb.data, r * static_cast<std::size_t>(k) + j) + 1)
                       << " ";
                std::int64_t l = lab ? detail::read_int(*lab, r) : 1;
                os << l << "\n";
            }
        }
    }

    os << "\nEnd\n";
}

}  // namespace meshioplusplus
