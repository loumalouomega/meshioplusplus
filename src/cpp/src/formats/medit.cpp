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
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"

// Project includes (private, not installed)
#include "../detail/row_writer.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// medit element keyword -> (meshio type, nodes per cell)
const std::unordered_map<std::string, std::pair<std::string, int>>& medit_to_meshio() {
    static const std::unordered_map<std::string, std::pair<std::string, int>> m = {
        {"Edges", {"line", 2}},           {"Triangles", {"triangle", 3}},
        {"Quadrilaterals", {"quad", 4}},  {"Tetrahedra", {"tetra", 4}},
        {"Prisms", {"wedge", 6}},         {"Pyramids", {"pyramid", 5}},
        {"Hexahedra", {"hexahedron", 8}}, {"Hexaedra", {"hexahedron", 8}},
    };
    return m;
}

// meshio type -> (medit keyword, nodes per cell), write order.
const std::vector<std::pair<std::string, std::pair<std::string, int>>>& meshio_to_medit() {
    static const std::vector<std::pair<std::string, std::pair<std::string, int>>> m = {
        {"line", {"Edges", 2}},           {"triangle", {"Triangles", 3}},
        {"quad", {"Quadrilaterals", 4}},  {"tetra", {"Tetrahedra", 4}},
        {"wedge", {"Prisms", 6}},         {"pyramid", {"Pyramids", 5}},
        {"hexahedron", {"Hexahedra", 8}},
    };
    return m;
}

// Whitespace/comment-skipping tokenizer over the whole file.
struct Tokenizer {
    std::string_view mBuf;
    std::size_t mPos = 0;
    explicit Tokenizer(std::string_view rB) : mBuf(rB) {}

    bool eof() const { return mPos >= mBuf.size(); }

    void skip_ws() {
        while (mPos < mBuf.size()) {
            char c = mBuf[mPos];
            if (c == '#') {  // comment to end of line
                while (mPos < mBuf.size() && mBuf[mPos] != '\n')
                    ++mPos;
            } else if (std::isspace(static_cast<unsigned char>(c))) {
                ++mPos;
            } else {
                break;
            }
        }
    }
    std::string next() {
        skip_ws();
        std::size_t start = mPos;
        while (mPos < mBuf.size() && !std::isspace(static_cast<unsigned char>(mBuf[mPos])) &&
               mBuf[mPos] != '#')
            ++mPos;
        return std::string(mBuf.substr(start, mPos - start));
    }
    std::int64_t next_int() { return std::strtoll(next().c_str(), nullptr, 10); }
    // A section's entry count: every entry takes at least a byte of the file,
    // so a count beyond its size (or a negative one) is corruption rather than
    // something to allocate or loop over.
    std::int64_t next_count() {
        return static_cast<std::int64_t>(
            detail::checked_count(next_int(), mBuf.size(), "Medit", "entry"));
    }
    double next_double() { return detail::parse_double(next()); }
    // Tokens on the line of the next token, without consuming anything.
    std::size_t tokens_on_next_line() {
        const std::size_t saved = mPos;
        skip_ws();
        std::size_t count = 0;
        while (mPos < mBuf.size() && mBuf[mPos] != '\n' && mBuf[mPos] != '#') {
            if (std::isspace(static_cast<unsigned char>(mBuf[mPos]))) {
                ++mPos;
                continue;
            }
            ++count;
            while (mPos < mBuf.size() && !std::isspace(static_cast<unsigned char>(mBuf[mPos])) &&
                   mBuf[mPos] != '#')
                ++mPos;
        }
        mPos = saved;
        return count;
    }
    void skip_line() {
        while (mPos < mBuf.size() && mBuf[mPos] != '\n')
            ++mPos;
        if (mPos < mBuf.size())
            ++mPos;
    }
};

void store_coord(NDArray& rA, std::size_t idx, double v) {
    if (rA.Dtype() == DType::Float32)
        rA.As<float>()[idx] = static_cast<float>(v);
    else
        rA.As<double>()[idx] = v;
}

// Both pickers iterate in sorted key order so the "first int field" chosen is
// stable regardless of the backend's storage order.
const NDArray* pick_first_int(const Mesh& rMesh) {
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& v = rMesh.PointData(name);
        DType t = v.Dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
            t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64)
            return &v;
    }
    return nullptr;
}

// Name of the first (sorted) integer cell-data field, or "" if none.
std::string pick_first_int_cell(const Mesh& rMesh) {
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        DType t = rMesh.CellData(name, 0).Dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
            t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64)
            return name;
    }
    return "";
}

}  // namespace

Mesh read_medit_ascii(const std::string& rPath) {
    const detail::FileSource buf_source =
        detail::open_source(rPath, "Could not open file: " + rPath);
    const std::string_view buf = buf_source.View();
    Tokenizer tok(buf);

    int dim = 0;
    DType coord_dtype = DType::Float64;
    Mesh mesh;
    std::vector<std::int64_t> point_ref;
    bool have_points = false;

    const auto& e2m = medit_to_meshio();

    while (!tok.eof()) {
        std::string kw = tok.next();
        if (kw.empty())
            break;
        if (kw == "MeshVersionFormatted") {
            std::int64_t v = tok.next_int();
            coord_dtype = (v <= 1) ? DType::Float32 : DType::Float64;
        } else if (kw == "Dimension") {
            dim = static_cast<int>(tok.next_int());
        } else if (kw == "Vertices") {
            std::int64_t n = tok.next_count();
            // No `Dimension` keyword (FEconv writes none): a vertex row holds
            // the coordinates and a reference.
            if (dim <= 0)
                dim = n > 0 ? static_cast<int>(tok.tokens_on_next_line()) - 1 : 3;
            if (dim < 1 || dim > 3)
                throw ReadError("Medit: cannot tell the dimension from the Vertices rows");
            NDArray pts(coord_dtype, {static_cast<std::size_t>(n), static_cast<std::size_t>(dim)});
            point_ref.resize(n);
            for (std::int64_t i = 0; i < n; ++i) {
                for (int c = 0; c < dim; ++c)
                    store_coord(pts, i * dim + c, tok.next_double());
                point_ref[i] = detail::checked_integer<std::int64_t>(tok.next_double(), "Medit");
            }
            mesh.AssignPoints(std::move(pts));
            have_points = true;
        } else if (e2m.count(kw)) {
            const auto& info = e2m.at(kw);
            const std::string& type = info.first;
            int k = info.second;
            std::int64_t n = tok.next_count();
            NDArray data(DType::Int64, {static_cast<std::size_t>(n), static_cast<std::size_t>(k)});
            NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
            std::int64_t* dp = data.As<std::int64_t>();
            std::int64_t* rp = ref.As<std::int64_t>();
            for (std::int64_t i = 0; i < n; ++i) {
                for (int j = 0; j < k; ++j)
                    dp[i * k + j] = tok.next_int() - 1;
                rp[i] = tok.next_int();
            }
            mesh.AddCellBlock(type, std::move(data));
            mesh.AppendCellData("medit:ref", std::move(ref));
        } else if (kw == "Corners") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n; ++i)
                tok.next();
        } else if (kw == "Normals") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * dim; ++i)
                tok.next();
        } else if (kw == "NormalAtVertices") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * 2; ++i)
                tok.next();
        } else if (kw == "SubDomainFromMesh") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * 4; ++i)
                tok.next();
        } else if (kw == "VertexOnGeometricVertex") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * 2; ++i)
                tok.next();
        } else if (kw == "VertexOnGeometricEdge") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * 3; ++i)
                tok.next();
        } else if (kw == "EdgeOnGeometricEdge") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n * 2; ++i)
                tok.next();
        } else if (kw == "Identifier" || kw == "Geometry") {
            tok.skip_line();
        } else if (kw == "RequiredVertices" || kw == "TangentAtVertices" || kw == "Tangents" ||
                   kw == "Ridges") {
            std::int64_t n = tok.next_count();
            for (std::int64_t i = 0; i < n; ++i)
                tok.next();
        } else if (kw == "End") {
            break;
        } else {
            throw ReadError("Medit: unknown keyword '" + kw + "'");
        }
    }

    if (!have_points)
        throw ReadError("Medit: expected Vertices");

    NDArray pr(DType::Int64, {point_ref.size()});
    for (std::size_t i = 0; i < point_ref.size(); ++i)
        pr.As<std::int64_t>()[i] = point_ref[i];
    mesh.AddPointData("medit:ref", std::move(pr));
    return mesh;
}

void write_medit_ascii(const std::string& rPath, const Mesh& rMesh) {
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t d = rMesh.PointDim();
    int version = (points.Dtype() == DType::Float32) ? 1 : 2;

    os << "MeshVersionFormatted " << version << "\n";
    os << "Dimension " << d << "\n";

    // Vertices
    os << "\nVertices\n" << n << "\n";
    const NDArray* vlabels = pick_first_int(rMesh);
    // Rows formatted in parallel chunks (row_writer.hpp), byte for byte.
    {
        const detail::DoubleView pv(points);
        const std::optional<detail::Int64View> labels =
            vlabels ? std::optional<detail::Int64View>(std::in_place, *vlabels) : std::nullopt;
        const detail::CNumber num;
        detail::write_row_chunks(os, n,
                                 [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                                     for (std::size_t i = First; i < Last; ++i) {
                                         for (std::size_t c = 0; c < d; ++c)
                                             num.Append(rBuf, "%.16e ", pv[i * d + c]);
                                         detail::append_int(rBuf, labels ? (*labels)[i] : 1);
                                         rBuf += '\n';
                                     }
                                 });
    }

    // Cells, grouped by medit element keyword.
    const std::string clabel_key = pick_first_int_cell(rMesh);
    for (const auto& mk : meshio_to_medit()) {
        const std::string& mtype = mk.first;
        const std::string& kw = mk.second.first;
        int k = mk.second.second;
        for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
            const auto cb = rMesh.Cells(ci);
            if (cb.Type() != mtype)
                continue;
            std::size_t count = cb.NumCells();
            os << "\n" << kw << "\n" << count << "\n";
            const NDArray* lab = (!clabel_key.empty() && ci < rMesh.CellDataNumBlocks(clabel_key))
                                     ? &rMesh.CellData(clabel_key, ci)
                                     : nullptr;
            const detail::Int64View conn(cb.Conn());
            const std::optional<detail::Int64View> labels =
                lab ? std::optional<detail::Int64View>(std::in_place, *lab) : std::nullopt;
            const std::size_t kk = static_cast<std::size_t>(k);
            detail::write_row_chunks(os, count,
                                     [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                                         for (std::size_t r = First; r < Last; ++r) {
                                             for (std::size_t j = 0; j < kk; ++j) {
                                                 detail::append_int(rBuf, conn[r * kk + j] + 1);
                                                 rBuf += ' ';
                                             }
                                             detail::append_int(rBuf, labels ? (*labels)[r] : 1);
                                             rBuf += '\n';
                                         }
                                     });
        }
    }

    os << "\nEnd\n";
}

}  // namespace meshioplusplus
