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
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

struct FaceBlock {
    std::size_t mSize = 0;
    std::vector<std::int64_t> mIdx;  // flat, 0-based
    std::vector<std::int64_t> mGids;
    std::size_t mCount = 0;
};

// The `vn` or `vt` rows of a file, flat. Rows may differ in width while
// reading; one that does is refused when the array is made.
struct AttributeRows {
    std::vector<double> mValues;
    std::size_t mRows = 0;
    std::size_t mWidth = 0;
    bool mUniform = true;

    void AddRow(const double* pValues, std::size_t Count) {
        if (mRows == 0)
            mWidth = Count;
        else if (Count != mWidth)
            mUniform = false;
        mValues.insert(mValues.end(), pValues, pValues + Count);
        ++mRows;
    }
};

std::string cell_type_for(std::size_t n) {
    if (n == 3)
        return "triangle";
    if (n == 4)
        return "quad";
    return "polygon";
}

NDArray make_point_data(const AttributeRows& rRows) {
    if (!rRows.mUniform)
        throw ReadError("OBJ: rows of one attribute with different lengths");
    NDArray a(DType::Float64, {rRows.mRows, rRows.mWidth});
    if (!rRows.mValues.empty())
        std::memcpy(a.As<double>(), rRows.mValues.data(), rRows.mValues.size() * sizeof(double));
    return a;
}

// `std::stoll(std::string(Token), nullptr, 10)` over a token that holds no
// blank: the leading integer, an optional sign included, with whatever follows
// ignored. Like `stoll`, a token with no digit throws `std::invalid_argument`
// and a value outside `int64` throws `std::out_of_range`.
std::int64_t obj_stoll(std::string_view Token) {
    std::size_t i = 0;
    bool negative = false;
    if (i < Token.size() && (Token[i] == '+' || Token[i] == '-')) {
        negative = Token[i] == '-';
        ++i;
    }
    if (i >= Token.size() || Token[i] < '0' || Token[i] > '9')
        throw std::invalid_argument("stoll");
    std::uint64_t magnitude = 0;
    const auto parsed =
        std::from_chars(Token.data() + i, Token.data() + Token.size(), magnitude);
    constexpr auto maximum = static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if (parsed.ec == std::errc::result_out_of_range || magnitude > maximum + (negative ? 1u : 0u))
        throw std::out_of_range("stoll");
    return negative ? static_cast<std::int64_t>(0 - magnitude)
                    : static_cast<std::int64_t>(magnitude);
}

// `v`, `vn` and `vt` carry numbers until the first token that is not one.
void obj_read_numbers(std::string_view Rest, std::vector<double>& rOut) {
    detail::TextStream iss(Rest);
    double x;
    while (iss >> x)
        rOut.push_back(x);
}

}  // namespace

Mesh read_obj(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    detail::TextCursor cursor(source.View());

    std::vector<double> points;  // flat x, y, z
    AttributeRows vn, vt;
    std::vector<FaceBlock> blocks;
    std::int64_t group_id = -1;
    std::vector<double> row;
    std::vector<std::int64_t> dat;

    while (!cursor.AtEnd()) {
        std::string_view line = cursor.Line();
        std::size_t b = 0, e = line.size();
        while (b < e && detail::text_is_blank(line[b]))
            ++b;
        while (e > b && detail::text_is_blank(line[e - 1]))
            --e;
        if (b == e || line[b] == '#')
            continue;
        line = line.substr(b, e - b);

        // The tag is the first blank-separated token; the rest is its payload.
        std::size_t t = 0;
        while (t < line.size() && !detail::text_is_blank(line[t]))
            ++t;
        const std::string_view tag = line.substr(0, t);
        const std::string_view rest = line.substr(t);
        if (tag == "v") {
            // A vertex short of three numbers keeps zeros for the rest, and a
            // value that does not parse ends the line, as a stream would.
            std::array<double, 3> p{0, 0, 0};
            detail::TextStream iss(rest);
            iss >> p[0] >> p[1] >> p[2];
            points.insert(points.end(), p.begin(), p.end());
        } else if (tag == "vn" || tag == "vt") {
            row.clear();
            obj_read_numbers(rest, row);
            (tag == "vn" ? vn : vt).AddRow(row.data(), row.size());
        } else if (tag == "f") {
            dat.clear();
            std::size_t i = 0;
            while (true) {
                while (i < rest.size() && detail::text_is_blank(rest[i]))
                    ++i;
                if (i >= rest.size())
                    break;
                const std::size_t first = i;
                std::size_t slash = std::string_view::npos;
                while (i < rest.size() && !detail::text_is_blank(rest[i])) {
                    if (rest[i] == '/' && slash == std::string_view::npos)
                        slash = i;
                    ++i;
                }
                const std::size_t last = slash == std::string_view::npos ? i : slash;
                dat.push_back(obj_stoll(rest.substr(first, last - first)) - 1);
            }
            std::size_t sz = dat.size();
            if (blocks.empty() || (blocks.back().mCount > 0 && blocks.back().mSize != sz)) {
                FaceBlock fb;
                fb.mSize = sz;
                blocks.push_back(std::move(fb));
            }
            FaceBlock& cur = blocks.back();
            if (cur.mCount == 0)
                cur.mSize = sz;
            cur.mIdx.insert(cur.mIdx.end(), dat.begin(), dat.end());
            cur.mGids.push_back(group_id);
            ++cur.mCount;
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
        if (fb.mCount > 0)
            nonempty.push_back(std::move(fb));

    Mesh mesh;
    std::size_t np = points.size() / 3;
    NDArray pts(DType::Float64, {np, 3});
    if (!points.empty())
        std::memcpy(pts.As<double>(), points.data(), points.size() * sizeof(double));
    mesh.AssignPoints(std::move(pts));

    if (vt.mRows > 0)
        mesh.AddPointData("obj:vt", make_point_data(vt));
    if (vn.mRows > 0)
        mesh.AddPointData("obj:vn", make_point_data(vn));

    if (!nonempty.empty()) {
        std::vector<NDArray> gid_blocks;
        for (auto& fb : nonempty) {
            NDArray data(DType::Int64, {fb.mCount, fb.mSize});
            if (!fb.mIdx.empty())
                std::memcpy(data.As<std::int64_t>(), fb.mIdx.data(),
                            fb.mIdx.size() * sizeof(std::int64_t));
            mesh.AddCellBlock(cell_type_for(fb.mSize), std::move(data));

            NDArray g(DType::Int64, {fb.mCount});
            std::memcpy(g.As<std::int64_t>(), fb.mGids.data(), fb.mCount * sizeof(std::int64_t));
            gid_blocks.push_back(std::move(g));
        }
        mesh.AddCellData("obj:group_ids", std::move(gid_blocks));
    }
    return mesh;
}

void write_obj(const std::string& rPath, const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() != "triangle" && cb.Type() != "quad" && cb.Type() != "polygon")
            throw WriteError(
                "Wavefront .obj files can only contain triangle, quad, "
                "or polygon cells.");

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    os << detail::provenance_render_lines(detail::SlotTier::Block, "# ");
    char buf[96];
    const detail::DoubleView point_values(points);
    for (std::size_t r = 0; r < num_points; ++r) {
        double x = (0 < dim) ? point_values[r * dim + 0] : 0.0;
        double y = (1 < dim) ? point_values[r * dim + 1] : 0.0;
        double z = (2 < dim) ? point_values[r * dim + 2] : 0.0;
        detail::snprintf_c(buf, sizeof(buf), "v %.17g %.17g %.17g\n", x, y, z);
        os << buf;
    }

    auto write_pd = [&](const char* key, const char* tag) {
        if (!rMesh.HasPointData(key))
            return;
        const NDArray& d = rMesh.PointData(key);
        std::size_t nc = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        const detail::DoubleView values(d);
        for (std::size_t r = 0; r < (d.Shape().empty() ? 0 : d.Shape()[0]); ++r) {
            os << tag;
            for (std::size_t c = 0; c < nc; ++c) {
                detail::snprintf_c(buf, sizeof(buf), " %.17g", values[r * nc + c]);
                os << buf;
            }
            os << '\n';
        }
    };
    write_pd("obj:vn", "vn");
    write_pd("obj:vt", "vt");

    for (const auto cb : rMesh.CellRange()) {
        const NDArray& conn_array = cb.Conn();
        std::size_t k = conn_array.Shape().size() >= 2 ? conn_array.Shape()[1] : 1;
        const detail::Int64View conn(conn_array);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << 'f';
            for (std::size_t j = 0; j < k; ++j)
                os << ' ' << (conn[r * k + j] + 1);
            os << '\n';
        }
    }
}

}  // namespace meshioplusplus
