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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/open_source.hpp"

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

namespace meshioplusplus {

namespace {

// netgen cell node count -> meshio type, per topological dimension.
const std::unordered_map<int, std::string>& netgen_type(int dim) {
    static const std::unordered_map<int, std::string> d0 = {{1, "vertex"}};
    static const std::unordered_map<int, std::string> d1 = {{2, "line"}};
    static const std::unordered_map<int, std::string> d2 = {
        {3, "triangle"}, {6, "triangle6"}, {4, "quad"}, {8, "quad8"}};
    static const std::unordered_map<int, std::string> d3 = {
        {4, "tetra"},    {5, "pyramid"},    {6, "wedge"},    {8, "hexahedron"},
        {10, "tetra10"}, {13, "pyramid13"}, {15, "wedge15"}, {20, "hexahedron20"}};
    switch (dim) {
        case 0:
            return d0;
        case 1:
            return d1;
        case 2:
            return d2;
        default:
            return d3;
    }
}

// netgen -> meshio node permutation: meshio[i] = netgen[pmap[i]].
const std::unordered_map<std::string, std::vector<int>>& n2m_pmap() {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"vertex", {0}},
        {"line", {0, 1}},
        {"triangle", {0, 1, 2}},
        {"triangle6", {0, 1, 2, 5, 3, 4}},
        {"quad", {0, 1, 2, 3}},
        {"quad8", {0, 1, 2, 3, 4, 7, 5, 6}},
        {"tetra", {0, 2, 1, 3}},
        {"tetra10", {0, 2, 1, 3, 5, 7, 4, 6, 9, 8}},
        {"pyramid", {0, 3, 2, 1, 4}},
        {"pyramid13", {0, 3, 2, 1, 4, 7, 6, 8, 5, 9, 12, 11, 10}},
        {"wedge", {0, 2, 1, 3, 5, 4}},
        {"wedge15", {0, 2, 1, 3, 5, 4, 7, 8, 6, 13, 14, 12, 9, 11, 10}},
        {"hexahedron", {0, 3, 2, 1, 4, 7, 6, 5}},
        {"hexahedron20", {0, 3, 2, 1, 4, 7, 6, 5, 10, 9, 11, 8, 16, 19, 18, 17, 14, 13, 15, 12}},
    };
    return m;
}

// meshio -> netgen node permutation (inverse of n2m_pmap).
const std::unordered_map<std::string, std::vector<int>>& m2n_pmap() {
    static const std::unordered_map<std::string, std::vector<int>> m = [] {
        std::unordered_map<std::string, std::vector<int>> out;
        for (const auto& kv : n2m_pmap()) {
            const auto& p = kv.second;
            std::vector<int> inv(p.size());
            for (std::size_t i = 0; i < p.size(); ++i)
                inv[p[i]] = static_cast<int>(i);
            out.emplace(kv.first, std::move(inv));
        }
        return out;
    }();
    return m;
}

int topo_dim(const std::string& rType) {
    auto it = topological_dimension().find(rType);
    return it == topological_dimension().end() ? -1 : it->second;
}

std::vector<std::string_view> netgen_split_ws(std::string_view rS) {
    return detail::split_blanks(rS);
}

std::string_view netgen_strip(std::string_view rS) {
    std::size_t a = 0, b = rS.size();
    while (a < b && std::isspace(static_cast<unsigned char>(rS[a])))
        ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(rS[b - 1])))
        --b;
    return rS.substr(a, b - a);
}

// Cursor over the file's lines, with comment/blank handling like the Python
// reader's _fast_forward_over_blank_lines.
struct LineCursor {
    std::vector<std::string_view> mLines;
    std::size_t mPos = 0;
    std::size_t mBytes = 0;

    explicit LineCursor(std::string_view Text) : mLines(detail::split_lines(Text)) {
        for (const auto line : mLines) {
            mBytes += line.size() + 1;
        }
    }

    bool Eof() const { return mPos >= mLines.size(); }

    // Next non-blank, non-comment line (stripped). Sets is_eof when exhausted.
    std::string_view NextReal(bool& rIsEof) {
        while (mPos < mLines.size()) {
            const auto s = netgen_strip(mLines[mPos++]);
            if (!s.empty() && s[0] != '#') {
                rIsEof = false;
                return s;
            }
        }
        rIsEof = true;
        return "";
    }

    // Next line raw (stripped), used for count lines that directly follow a
    // keyword; skips any stray blank/comment lines defensively.
    std::string_view NextCount() {
        bool eof = false;
        return NextReal(eof);
    }
};

struct NetgenRawBlock {
    std::string mType;
    std::vector<std::vector<std::int64_t>> mRows;  // meshio node order, 0-based
    std::vector<std::int64_t> mIndex;
};

std::int64_t netgen_integer(std::string_view rToken) {
    std::int64_t value = 0;
    if (!detail::parse_int_token(rToken, value))
        throw ReadError("Netgen: invalid integer '" + std::string(rToken) + "'");
    return value;
}

void read_cells(LineCursor& rC, const std::string& rSection, std::vector<NetgenRawBlock>& rBlocks,
                bool TwoLines) {
    int dim, pi0, i_index, fixed_nump = -1;
    if (rSection == "pointelements") {
        dim = 0;
        pi0 = 0;
        i_index = 1;
        fixed_nump = 1;
    } else if (rSection.rfind("edgesegments", 0) == 0) {
        dim = 1;
        pi0 = 2;
        i_index = 0;
        fixed_nump = 2;
    } else if (rSection.rfind("surfaceelements", 0) == 0) {
        dim = 2;
        pi0 = 5;
        i_index = 1;
    } else if (rSection == "volumeelements") {
        dim = 3;
        pi0 = 2;
        i_index = 0;
    } else {
        throw ReadError("Netgen: unknown cell section '" + rSection + "'");
    }

    const std::size_t num_cells =
        detail::checked_count(netgen_integer(rC.NextCount()), rC.mBytes, "Netgen", "cell");
    const auto& tmap = netgen_type(dim);

    for (std::size_t k = 0; k < num_cells; ++k) {
        bool eof = false;
        const auto line = rC.NextReal(eof);
        if (eof)
            throw ReadError("Netgen: unexpected end of file in " + rSection);
        const auto data = netgen_split_ws(line);
        // The node count sits at a fixed column; check the row reaches it.
        detail::need_tokens(data, dim == 2 ? 5 : (dim == 3 ? 2 : 0), "Netgen");

        std::int64_t nump = fixed_nump;
        if (dim == 2)
            nump = netgen_integer(data[4]);
        else if (dim == 3)
            nump = netgen_integer(data[1]);

        auto tit = nump >= 1 && nump <= 20 ? tmap.find(static_cast<int>(nump)) : tmap.end();
        if (tit != tmap.end())
            detail::need_tokens(
                data, static_cast<std::size_t>(std::max<std::int64_t>(i_index + 1, pi0 + nump)),
                "Netgen");
        std::int64_t index = tit == tmap.end() ? 0 : netgen_integer(data[i_index]);
        if (tit == tmap.end())
            throw ReadError("Netgen: unsupported element with " + std::to_string(nump) + " nodes");
        const std::string& type = tit->second;

        std::vector<std::int64_t> pi(nump);
        for (int j = 0; j < nump; ++j)
            pi[j] = netgen_integer(data[pi0 + j]);

        if (rBlocks.empty() || rBlocks.back().mType != type) {
            rBlocks.push_back(NetgenRawBlock{type, {}, {}});
        }
        rBlocks.back().mRows.push_back(std::move(pi));
        rBlocks.back().mIndex.push_back(index);
        if (TwoLines && rSection == "edgesegmentsgi2") {
            rC.NextReal(eof);
            if (eof)
                throw ReadError("Netgen: unexpected end of file in two-line edge data");
        }
    }
}

}  // namespace

Mesh read_netgen(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    std::string bytes;
    const bool gzip = rPath.size() >= 7 && rPath.compare(rPath.size() - 7, 7, ".vol.gz") == 0;
    if (gzip) {
        const std::string_view compressed = source.View();
        std::size_t pos = 0;
        do {
            std::size_t consumed = 0;
            bytes += detail::zlib_inflate(std::string_view(compressed).substr(pos), 31, &consumed,
                                          "Netgen");
            pos += consumed;
        } while (pos < compressed.size());
    }
    LineCursor c(gzip ? std::string_view(bytes) : source.View());

    bool eof = false;
    std::string line(c.NextReal(eof));
    if (line != "mesh3d")
        throw ReadError("Not a valid Netgen mesh");

    int dimension = 3;
    std::vector<double> raw_points;  // flat, 3 per point
    std::int64_t num_points = 0;
    std::vector<NetgenRawBlock> blocks;
    std::map<std::string, NDArray> fields;
    bool two_lines = false;
    const std::map<std::string, int> codims = {
        {"materials", 0}, {"bcnames", 1}, {"cd2names", 2}, {"cd3names", 3}};
    auto count = [&]() {
        return detail::checked_count(netgen_integer(c.NextCount()), c.mBytes, "Netgen", "section");
    };

    while (true) {
        line = c.NextReal(eof);
        if (eof)
            break;
        if (line == "dimension") {
            dimension = static_cast<int>(detail::strtoll_token(c.NextCount()));
            if (dimension < 1 || dimension > 3)
                throw ReadError("Netgen: dimension must be 1, 2 or 3");
        } else if (line == "geomtype") {
            c.NextCount();  // value; ignored
        } else if (line == "points") {
            // A point row is at least a few bytes: bound the count by the file.
            num_points = static_cast<std::int64_t>(detail::checked_count(
                detail::strtoll_token(c.NextCount()), c.mBytes, "Netgen", "point"));
            raw_points.resize(static_cast<std::size_t>(num_points) * 3, 0.0);
            for (std::int64_t i = 0; i < num_points; ++i) {
                const auto pl = c.NextReal(eof);
                if (eof)
                    throw ReadError("Netgen: unexpected EOF in points");
                const auto toks = netgen_split_ws(pl);
                for (int j = 0; j < 3 && j < static_cast<int>(toks.size()); ++j)
                    raw_points[i * 3 + j] = detail::parse_double_prefix(toks[j]);
            }
        } else if (line == "pointelements" || line == "edgesegments" || line == "edgesegmentsgi" ||
                   line == "surfaceelements" || line == "surfaceelementsgi" ||
                   line == "surfaceelementsuv" || line == "volumeelements") {
            read_cells(c, line, blocks, two_lines);
        } else if (line == "edgesegmentsgi2") {
            read_cells(c, line, blocks, two_lines);
        } else if (netgen_split_ws(line) ==
                   std::vector<std::string_view>{"surf1", "surf2", "p1", "p2"}) {
            two_lines = true;
        } else if (codims.count(line)) {
            const int edim = dimension - codims.at(line);
            const std::size_t n = count();
            for (std::size_t i = 0; i < n; ++i) {
                if (c.Eof())
                    throw ReadError("Netgen: unexpected EOF in name table");
                const auto tokens = netgen_split_ws(c.mLines[c.mPos++]);
                if (tokens.size() != 2)
                    continue;  // an unnamed slot, as in the Python reference
                NDArray data(DType::Int64, {2});
                if (!detail::parse_int_token(tokens[0], data.As<std::int64_t>()[0]))
                    throw ReadError("Netgen: invalid name-table index");
                data.As<std::int64_t>()[1] = edim;
                fields.insert_or_assign(std::string(tokens[1]), std::move(data));
            }
        } else if (line == "identifications" || line == "identificationtypes") {
            const std::string key = "netgen:" + line;
            const std::size_t n = count();
            const bool pairs = line == "identifications";
            NDArray data(DType::Int64,
                         pairs ? std::vector<std::size_t>{n, 3} : std::vector<std::size_t>{1, n});
            if (n) {
                for (std::size_t i = 0; i < (pairs ? n : 1); ++i) {
                    const auto tokens = netgen_split_ws(c.NextReal(eof));
                    const std::size_t width = pairs ? 3 : n;
                    if (eof || tokens.size() != width)
                        throw ReadError("Netgen: malformed periodic table");
                    for (std::size_t j = 0; j < width; ++j)
                        if (!detail::parse_int_token(tokens[j],
                                                     data.As<std::int64_t>()[i * width + j]))
                            throw ReadError("Netgen: invalid periodic-table integer");
                }
            }
            fields.insert_or_assign(key, std::move(data));
        } else if (line == "face_colours" || line == "singular_edge_left" ||
                   line == "singular_edge_right" || line == "singular_face_inside" ||
                   line == "singular_face_outside" || line == "singular_points") {
            const std::size_t n = count();
            for (std::size_t i = 0; i < n; ++i) {
                if (c.Eof())
                    throw ReadError("Netgen: unexpected EOF in auxiliary section");
                ++c.mPos;
            }
        } else if (line == "endmesh") {
            break;
        } else {
            throw ReadError("Netgen: unknown token '" + line + "'");
        }
    }

    Mesh mesh;
    NDArray pts(DType::Float64,
                {static_cast<std::size_t>(num_points), static_cast<std::size_t>(dimension)});
    double* pp = pts.As<double>();
    for (std::int64_t i = 0; i < num_points; ++i)
        for (int j = 0; j < dimension; ++j)
            pp[i * dimension + j] = raw_points[i * 3 + j];
    mesh.AssignPoints(std::move(pts));

    std::vector<NDArray> index_blocks;
    for (auto& b : blocks) {
        const std::vector<int>& pmap = n2m_pmap().at(b.mType);
        std::size_t n = b.mRows.size();
        std::size_t k = pmap.size();
        NDArray data(DType::Int64, {n, k});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t r = 0; r < n; ++r)
            for (std::size_t j = 0; j < k; ++j)
                dp[r * k + j] = b.mRows[r][pmap[j]] - 1;
        mesh.AddCellBlock(b.mType, std::move(data));

        NDArray idx(DType::Int64, {n});
        for (std::size_t r = 0; r < n; ++r)
            idx.As<std::int64_t>()[r] = b.mIndex[r];
        index_blocks.push_back(std::move(idx));
    }
    mesh.AddCellData("netgen:index", std::move(index_blocks));
    for (auto& [name, data] : fields)
        mesh.AddFieldData(name, std::move(data));

    return mesh;
}

namespace {

void write_block(std::ostream& rOs, Mesh::CellView cb, const NDArray* pIndex) {
    if (cb.NumCells() == 0)
        return;
    int dim = topo_dim(cb.Type());
    const std::vector<int>& pmap = m2n_pmap().at(cb.Type());
    const int np = static_cast<int>(pmap.size());

    std::vector<std::int64_t> pre, post;
    int i_index = 0;
    if (dim == 0) {
        post = {1};
        i_index = 1;
    } else if (dim == 1) {
        pre = {1, 0};
        post = {-1, -1, 0, 0, 1, 0, 1, 0};
    } else if (dim == 2) {
        pre = {1, 1, 0, 0, np};
        i_index = 1;
    } else {  // dim == 3
        pre = {1, np};
    }

    const NDArray& conn = cb.Conn();
    const std::size_t n = cb.NumCells();
    for (std::size_t r = 0; r < n; ++r) {
        std::vector<std::int64_t> cols;
        cols.reserve(pre.size() + np + post.size());
        for (auto v : pre)
            cols.push_back(v);
        for (int j = 0; j < np; ++j)
            cols.push_back(detail::read_int(conn, r * np + pmap[j]) + 1);
        for (auto v : post)
            cols.push_back(v);
        if (pIndex)
            cols[i_index] = detail::read_int(*pIndex, r);

        for (std::size_t j = 0; j < cols.size(); ++j)
            rOs << cols[j] << (j + 1 == cols.size() ? '\n' : ' ');
    }
}

}  // namespace

void write_netgen(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt) {
    const bool gzip = rPath.size() >= 7 && rPath.compare(rPath.size() - 7, 7, ".vol.gz") == 0;
    if (gzip && !detail::zlib_available())
        throw WriteError("Netgen: gzip needs -DMESHIOPLUSPLUS_WITH_ZLIB=ON");
    auto file = detail::make_classic_ofstream();
    auto buffer = detail::make_classic_ostringstream();
    if (!gzip)
        file.open(rPath, std::ios::binary);
    if (!gzip && !file)
        throw WriteError("Could not open file for writing: " + rPath);
    std::ostream& f = gzip ? static_cast<std::ostream&>(buffer) : static_cast<std::ostream&>(file);

    const NDArray& points = rMesh.Points();
    const int dimension = points.Shape().size() >= 2 ? static_cast<int>(points.Shape()[1]) : 3;

    // Pick the single integer cell index, preferring "netgen:index".
    bool have_index = false;
    std::string index_key;
    if (rMesh.HasCellData("netgen:index")) {
        have_index = true;
        index_key = "netgen:index";
    } else {
        for (const auto& name : rMesh.CellDataNames()) {
            if (rMesh.CellDataNumBlocks(name) == 0)
                continue;
            DType t = rMesh.CellData(name, 0).Dtype();
            if (t != DType::Float32 && t != DType::Float64) {
                have_index = true;
                index_key = name;
                break;
            }
        }
    }
    auto index_for = [&](std::size_t ci) -> const NDArray* {
        if (!have_index || ci >= rMesh.CellDataNumBlocks(index_key))
            return nullptr;
        return &rMesh.CellData(index_key, ci);
    };

    std::int64_t per_dim[4] = {0, 0, 0, 0};
    for (const auto cb : rMesh.CellRange()) {
        int d = topo_dim(cb.Type());
        if (d >= 0 && d <= 3)
            per_dim[d] += static_cast<std::int64_t>(cb.NumCells());
    }

    f << detail::provenance_render_lines(detail::SlotTier::Block, "# ");
    f << "mesh3d\n\n";
    f << "dimension\n" << dimension << "\n\n";
    f << "geomtype\n0\n";

    f << "\n# surfnr    bcnr   domin  domout      np      p1      p2      p3\n";
    f << "surfaceelements\n" << per_dim[2] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 2)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n#  matnr      np      p1      p2      p3      p4\n";
    f << "volumeelements\n" << per_dim[3] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 3)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n# surfid  0   p1   p2   trignum1    trignum2   domin/surfnr1    "
         "domout/surfnr2   ednr1   dist1   ednr2   dist2\n";
    f << "edgesegmentsgi2\n" << per_dim[1] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 1)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    f << "\n#          X             Y             Z\n";
    f << "points\n" << rMesh.NumPoints() << "\n";
    std::string fmt = "%" + rFloatFmt;
    char buf[64];
    const std::size_t npts = rMesh.NumPoints();
    for (std::size_t i = 0; i < npts; ++i) {
        for (int j = 0; j < 3; ++j) {
            double v = (j < dimension) ? detail::read_double(points, i * dimension + j) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), fmt.c_str(), v);
            f << buf << (j == 2 ? '\n' : ' ');
        }
    }

    f << "\n#          pnum             index\n";
    f << "pointelements\n" << per_dim[0] << "\n";
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci)
        if (topo_dim(rMesh.Cells(ci).Type()) == 0)
            write_block(f, rMesh.Cells(ci), index_for(ci));

    for (const char* key : {"netgen:identifications", "netgen:identificationtypes"}) {
        if (!rMesh.HasFieldData(key))
            continue;
        const NDArray& data = rMesh.FieldData(key);
        const bool pairs = std::string(key) == "netgen:identifications";
        if (pairs && (data.Ndim() != 2 || data.Shape()[1] != 3))
            throw WriteError("Netgen: identifications must have shape (n,3)");
        const std::size_t n = pairs ? data.Shape()[0] : data.Size();
        f << '\n' << (pairs ? "identifications" : "identificationtypes") << '\n' << n << '\n';
        for (std::size_t i = 0; i < data.Size(); ++i)
            f << detail::read_int(data, i) << ((pairs ? (i % 3 == 2) : (i + 1 == n)) ? '\n' : ' ');
    }
    const char* codim_names[] = {"materials", "bcnames", "cd2names", "cd3names"};
    for (int codim = 0; codim <= dimension; ++codim) {
        std::map<std::int64_t, std::string> names;
        for (const std::string& name : rMesh.FieldDataNames()) {
            if (name == "netgen:identifications" || name == "netgen:identificationtypes")
                continue;
            const NDArray& data = rMesh.FieldData(name);
            if (data.Size() != 2)
                throw WriteError("Netgen: name-table field '" + name +
                                 "' must hold [id, dimension]");
            if (detail::read_int(data, 1) == dimension - codim)
                names[detail::read_int(data, 0)] = name;
        }
        if (names.empty()) {
            for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
                if (topo_dim(rMesh.Cells(ci).Type()) != dimension - codim)
                    continue;
                const NDArray* idx = index_for(ci);
                if (idx)
                    for (std::size_t i = 0; i < idx->Size(); ++i) {
                        const auto id = detail::read_int(*idx, i);
                        names[id] = "cd" + std::to_string(codim) + "_" + std::to_string(id);
                    }
            }
        }
        if (names.empty())
            continue;
        const std::int64_t max_id = names.rbegin()->first;
        if (max_id > 10000000)
            throw WriteError("Netgen: name-table id exceeds the 10000000-entry budget");
        f << '\n' << codim_names[codim] << '\n' << max_id << '\n';
        for (std::int64_t id = 1; id <= max_id; ++id) {
            const auto it = names.find(id);
            f << id << ' ' << (it == names.end() ? "" : it->second) << '\n';
        }
    }
    f << "\nendmesh\n";
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    if (gzip) {
        const std::string text = buffer.str();
        gzFile out = gzopen(rPath.c_str(), "wb");
        if (!out)
            throw WriteError("Netgen: could not open gzip output " + rPath);
        bool ok = true;
        for (std::size_t pos = 0; pos < text.size();) {
            const unsigned n =
                static_cast<unsigned>(std::min<std::size_t>(text.size() - pos, 1 << 20));
            if (gzwrite(out, text.data() + pos, n) != static_cast<int>(n)) {
                ok = false;
                break;
            }
            pos += n;
        }
        if (gzclose(out) != Z_OK || !ok)
            throw WriteError("Netgen: failed writing gzip output " + rPath);
    }
#endif
    if (!f)
        throw WriteError("Netgen: failed writing " + rPath);
}

}  // namespace meshioplusplus
