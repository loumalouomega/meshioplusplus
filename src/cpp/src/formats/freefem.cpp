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
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// Next non-blank line's whitespace tokens, as views into the text @p rCursor
// walks (which the caller keeps alive).
bool next_tokens(detail::TextCursor& rCursor, std::vector<std::string_view>& rOut) {
    while (!rCursor.AtEnd()) {
        detail::split_blanks(rCursor.Line(), rOut);
        if (!rOut.empty())
            return true;
    }
    return false;
}

}  // namespace

Mesh read_freefem(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    detail::TextCursor cursor(source.View());

    std::vector<std::string_view> tok;
    if (!next_tokens(cursor, tok) || tok.size() != 3)
        throw ReadError("FreeFem: expected a 3-integer header");
    // Every vertex and element is a line of at least two bytes.
    const std::size_t max_rows = source.Size() / 2;
    const auto nver = static_cast<std::int64_t>(detail::checked_count(
        detail::strtoll_token(tok[0]), max_rows, "FreeFem", "vertex"));
    const auto n1 = static_cast<std::int64_t>(
        detail::checked_count(std::max<long long>(0, detail::strtoll_token(tok[1])),
                              max_rows, "FreeFem", "element"));
    const auto n2 = static_cast<std::int64_t>(
        detail::checked_count(std::max<long long>(0, detail::strtoll_token(tok[2])),
                              max_rows, "FreeFem", "element"));

    if (!next_tokens(cursor, tok))
        throw ReadError("FreeFem: missing vertices");
    const int dim = static_cast<int>(tok.size()) - 1;
    if (dim != 2 && dim != 3)
        throw ReadError("FreeFem: bad vertex dimension");

    Mesh mesh;
    NDArray pts(DType::Float64, {static_cast<std::size_t>(nver), static_cast<std::size_t>(dim)});
    NDArray pref(DType::Int64, {static_cast<std::size_t>(nver)});
    for (std::int64_t i = 0; i < nver; ++i) {
        if (i > 0 && !next_tokens(cursor, tok))
            throw ReadError("FreeFem: truncated vertices");
        detail::need_tokens(tok, static_cast<std::size_t>(dim) + 1, "FreeFem");
        for (int c = 0; c < dim; ++c)
            pts.As<double>()[i * dim + c] = detail::parse_double_prefix(tok[c]);
        pref.As<std::int64_t>()[i] = detail::strtoll_token(tok[dim]);
    }
    mesh.AssignPoints(std::move(pts));
    mesh.AddPointData("freefem:ref", std::move(pref));

    const char* t1 = dim == 2 ? "triangle" : "tetra";
    const int lnv1 = dim == 2 ? 3 : 4;
    const char* t2 = dim == 2 ? "line" : "triangle";
    const int lnv2 = dim == 2 ? 2 : 3;

    std::vector<NDArray> cell_refs;
    auto read_block = [&](std::int64_t n, const char* type, int lnv) {
        if (n <= 0)
            return;
        NDArray data(DType::Int64, {static_cast<std::size_t>(n), static_cast<std::size_t>(lnv)});
        NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
        for (std::int64_t k = 0; k < n; ++k) {
            if (!next_tokens(cursor, tok))
                throw ReadError("FreeFem: truncated elements");
            detail::need_tokens(tok, static_cast<std::size_t>(lnv) + 1, "FreeFem");
            for (int j = 0; j < lnv; ++j)
                data.As<std::int64_t>()[k * lnv + j] =
                    detail::zero_based(detail::strtoll_token(tok[j]));
            ref.As<std::int64_t>()[k] = detail::strtoll_token(tok[lnv]);
        }
        mesh.AddCellBlock(type, std::move(data));
        cell_refs.push_back(std::move(ref));
    };
    read_block(n1, t1, lnv1);
    read_block(n2, t2, lnv2);
    if (!cell_refs.empty())
        mesh.AddCellData("freefem:ref", std::move(cell_refs));

    return mesh;
}

void write_freefem(const std::string& rPath, const Mesh& rMesh) {
    const int dim = static_cast<int>(rMesh.PointDim());
    if (dim != 2 && dim != 3)
        throw WriteError("FreeFem: can only write 2D/3D meshes");

    const std::string t1 = dim == 2 ? "triangle" : "tetra";
    const std::string t2 = dim == 2 ? "line" : "triangle";

    // Reject unsupported cell types so the shim falls back to Python (which
    // warns and skips). This keeps behaviour identical to the reference impl.
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() != t1 && cb.Type() != t2)
            throw WriteError("FreeFem: unsupported cell type " + cb.Type());

    const bool has_ref = rMesh.HasCellData("freefem:ref");

    struct Row {
        Mesh::CellView mCb;
        const NDArray* mRef;
    };
    std::vector<Row> b1, b2;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const NDArray* ref = (has_ref && i < rMesh.CellDataNumBlocks("freefem:ref"))
                                 ? &rMesh.CellData("freefem:ref", i)
                                 : nullptr;
        if (rMesh.Cells(i).Type() == t1)
            b1.push_back({rMesh.Cells(i), ref});
        else
            b2.push_back({rMesh.Cells(i), ref});
    }
    auto count = [](const std::vector<Row>& b) {
        std::size_t n = 0;
        for (const auto& r : b)
            n += r.mCb.NumCells();
        return n;
    };

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t nver = rMesh.NumPoints();
    f << nver << " " << count(b1) << " " << count(b2) << "\n";

    const NDArray* pref =
        rMesh.HasPointData("freefem:ref") ? &rMesh.PointData("freefem:ref") : nullptr;

    const NDArray& points = rMesh.Points();
    char buf[32];
    const detail::DoubleView point_values(points);
    std::optional<detail::Int64View> point_refs;
    if (pref)
        point_refs.emplace(*pref);
    for (std::size_t i = 0; i < nver; ++i) {
        for (int c = 0; c < dim; ++c) {
            detail::snprintf_c(buf, sizeof(buf), "%.16e", point_values[i * dim + c]);
            f << buf << " ";
        }
        f << (point_refs ? (*point_refs)[i] : 0) << "\n";
    }
    auto write_block = [&](const std::vector<Row>& b, int lnv) {
        for (const auto& r : b) {
            std::size_t n = r.mCb.NumCells();
            const NDArray& conn_array = r.mCb.Conn();
            std::size_t k = detail::cols(conn_array);
            const detail::Int64View conn(conn_array);
            std::optional<detail::Int64View> refs;
            if (r.mRef)
                refs.emplace(*r.mRef);
            for (std::size_t rr = 0; rr < n; ++rr) {
                for (int j = 0; j < lnv && static_cast<std::size_t>(j) < k; ++j)
                    f << (conn[rr * k + j] + 1) << " ";
                f << (refs ? (*refs)[rr] : 0) << "\n";
            }
        }
    };
    write_block(b1, dim == 2 ? 3 : 4);
    write_block(b2, dim == 2 ? 2 : 3);
}

}  // namespace meshioplusplus
