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
#include <algorithm>
#include <bit>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"

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
    std::string_view next() {
        skip_ws();
        std::size_t start = mPos;
        while (mPos < mBuf.size() && !std::isspace(static_cast<unsigned char>(mBuf[mPos])) &&
               mBuf[mPos] != '#')
            ++mPos;
        return mBuf.substr(start, mPos - start);
    }
    std::int64_t next_int() { return detail::strtoll_token(next()); }
    // A section's entry count: every entry takes at least a byte of the file,
    // so a count beyond its size (or a negative one) is corruption rather than
    // something to allocate or loop over.
    std::int64_t next_count() {
        return static_cast<std::int64_t>(
            detail::checked_count(next_int(), mBuf.size(), "Medit", "entry"));
    }
    double next_double() { return detail::parse_double_prefix(next()); }
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

// An id-list keyword and the region it becomes: Corners and RequiredVertices list
// vertices; Ridges and RequiredEdges list entries of the Edges section, and
// RequiredTriangles of the Triangles section.
struct MeditIdList {
    const char* mKeyword;
    const char* mRegion;
    RegionKind mKind;
    const char* mCellType;  // the type whose section the ids index; null for points
};

const MeditIdList kMeditIdLists[] = {
    {"Corners", "medit:Corners", RegionKind::Point, nullptr},
    {"RequiredVertices", "medit:RequiredVertices", RegionKind::Point, nullptr},
    {"Ridges", "medit:Ridges", RegionKind::Cell, "line"},
    {"RequiredEdges", "medit:RequiredEdges", RegionKind::Cell, "line"},
    {"RequiredTriangles", "medit:RequiredTriangles", RegionKind::Cell, "triangle"},
};

const MeditIdList* medit_id_list(const std::string& rKeyword) {
    for (const MeditIdList& l : kMeditIdLists)
        if (rKeyword == l.mKeyword)
            return &l;
    return nullptr;
}

/// The id-list sections to write: each region named like an id-list keyword, as 0-based ids
/// within the section (cell regions: within the first block of the section's type).
std::vector<std::pair<const MeditIdList*, std::vector<std::int64_t>>> medit_id_list_sections(
    const Mesh& rMesh) {
    std::vector<std::pair<const MeditIdList*, std::vector<std::int64_t>>> out;
    std::vector<std::int64_t> block_start(rMesh.NumCellBlocks() + 1, 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        block_start[b + 1] = block_start[b] + static_cast<std::int64_t>(rMesh.Cells(b).NumCells());
    for (const MeditIdList& list : kMeditIdLists) {
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
            const Region& reg = rMesh.Region(r);
            if (reg.mName != list.mRegion || reg.mKind != list.mKind)
                continue;
            std::vector<std::int64_t> ids(reg.Entries(), reg.Entries() + reg.NumEntries());
            if (list.mKind == RegionKind::Cell) {
                std::size_t block = rMesh.NumCellBlocks();
                for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
                    if (rMesh.Cells(b).Type() == list.mCellType) {
                        block = b;
                        break;
                    }
                if (block == rMesh.NumCellBlocks()) {
                    log::warn("Medit: region '{}' needs {} cells; not written", list.mRegion,
                              list.mCellType);
                    break;
                }
                const std::size_t before = ids.size();
                ids.erase(std::remove_if(ids.begin(), ids.end(),
                                         [&](std::int64_t Id) {
                                             return Id < block_start[block] ||
                                                    Id >= block_start[block + 1];
                                         }),
                          ids.end());
                if (ids.size() != before)
                    log::warn(
                        "Medit: region '{}' holds cells outside the first {} block; they are "
                        "not written",
                        list.mRegion, list.mCellType);
                for (std::int64_t& id : ids)
                    id -= block_start[block];
            }
            out.emplace_back(&list, std::move(ids));
            break;
        }
    }
    return out;
}

/// The region an id-list section becomes: ids checked against the section they index.
void medit_add_id_list_region(Mesh& rMesh, const MeditIdList& rList,
                              std::vector<std::int64_t> Ids) {
    std::int64_t base = 0, limit = static_cast<std::int64_t>(rMesh.NumPoints());
    if (rList.mKind == RegionKind::Cell) {
        bool found = false;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks() && !found; ++b) {
            if (rMesh.Cells(b).Type() == rList.mCellType) {
                limit = static_cast<std::int64_t>(rMesh.Cells(b).NumCells());
                found = true;
            } else {
                base += static_cast<std::int64_t>(rMesh.Cells(b).NumCells());
            }
        }
        if (!found) {
            log::warn("Medit: {} lists {} cells the file does not have", rList.mKeyword,
                      rList.mCellType);
            return;
        }
    }
    for (std::int64_t& id : Ids) {
        if (id < 0 || id >= limit)
            throw ReadError(std::string("Medit: ") + rList.mKeyword +
                            " names an entry the file does not have");
        id += base;
    }
    NDArray entries = NDArray::Uninit(DType::Int64, {Ids.size()});
    std::copy(Ids.begin(), Ids.end(), entries.As<std::int64_t>());
    rMesh.AddRegion(Region(rList.mRegion, rList.mKind, std::move(entries)));
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
    std::vector<std::pair<const MeditIdList*, std::vector<std::int64_t>>> id_lists;

    while (!tok.eof()) {
        std::string kw(tok.next());
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
        } else if (const MeditIdList* id_list = medit_id_list(kw)) {
            // Vertex or element ids (1-based), kept as a region.
            const std::int64_t n = tok.next_count();
            std::vector<std::int64_t> ids(static_cast<std::size_t>(n));
            for (std::int64_t& id : ids)
                id = tok.next_int() - 1;
            id_lists.emplace_back(id_list, std::move(ids));
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
        } else if (kw == "TangentAtVertices" || kw == "Tangents") {
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

    for (auto& [list, ids] : id_lists)
        medit_add_id_list_region(mesh, *list, std::move(ids));
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
        // Int64View is non-copyable: emplace rather than a ternary, which MSVC copies.
        std::optional<detail::Int64View> labels;
        if (vlabels)
            labels.emplace(*vlabels);
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
            std::optional<detail::Int64View> labels;
            if (lab)
                labels.emplace(*lab);
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

    // Regions named like the id-list keywords, written back as those keywords.
    for (const auto& [list, ids] : medit_id_list_sections(rMesh)) {
        os << "\n" << list->mKeyword << "\n" << ids.size() << "\n";
        for (const std::int64_t id : ids)
            os << (id + 1) << "\n";
    }

    os << "\nEnd\n";
}

// --- binary (.meshb, libMeshb GMF) ----------------------------------------------------------

namespace {

constexpr int kMbEnd = 54;
constexpr int kMbMaxCode = 198;
constexpr int kMbVertices = 4;

bool medit_is_reserved(int Code) {
    return Code == 0 || Code == 2 || Code == 11 || Code == 12 || Code == 51;
}

// The element sections: libMeshb code -> (meshio type, nodes per cell).
struct MeditBinCell {
    int mCode;
    const char* mType;
    int mNodes;
};
const MeditBinCell kMeditBinCells[] = {
    {5, "line", 2},  {6, "triangle", 3},    {7, "quad", 4},     {8, "tetra", 4},
    {9, "wedge", 6}, {10, "hexahedron", 8}, {49, "pyramid", 5},
};

// The id-list sections: libMeshb code -> keyword.
struct MeditBinIdList {
    int mCode;
    const char* mKeyword;
};
const MeditBinIdList kMeditBinIdLists[] = {
    {13, "Corners"},           {14, "Ridges"}, {15, "RequiredVertices"}, {16, "RequiredEdges"},
    {17, "RequiredTriangles"},
};

/// Little- or big-endian reads at a byte offset, bounds-checked.
struct MeditBinCursor {
    const char* mData;
    std::size_t mSize;
    std::size_t mAt = 0;
    bool mSwap = false;

    [[noreturn]] static void Fail(const std::string& rWhy) { throw ReadError("Medit: " + rWhy); }

    std::size_t Left() const { return mSize - mAt; }

    std::int64_t Int(std::size_t Width) {
        if (Left() < Width)
            Fail("the file is truncated");
        std::int64_t v;
        if (Width == 4) {
            std::uint32_t u;
            std::memcpy(&u, mData + mAt, 4);
            if (mSwap)
                u = detail::bswap32(u);
            v = static_cast<std::int32_t>(u);
        } else {
            std::uint64_t u;
            std::memcpy(&u, mData + mAt, 8);
            if (mSwap)
                u = detail::bswap64(u);
            v = static_cast<std::int64_t>(u);
        }
        mAt += Width;
        return v;
    }

    double Real(std::size_t Width) {
        if (Left() < Width)
            Fail("the file is truncated");
        double v;
        if (Width == 4) {
            std::uint32_t u;
            std::memcpy(&u, mData + mAt, 4);
            if (mSwap)
                u = detail::bswap32(u);
            v = std::bit_cast<float>(u);
        } else {
            std::uint64_t u;
            std::memcpy(&u, mData + mAt, 8);
            if (mSwap)
                u = detail::bswap64(u);
            v = std::bit_cast<double>(u);
        }
        mAt += Width;
        return v;
    }
};

void medit_put_int(std::string& rOut, std::int64_t V, std::size_t Width) {
    const auto u = static_cast<std::uint64_t>(V);
    for (std::size_t k = 0; k < Width; ++k)
        rOut.push_back(static_cast<char>((u >> (8 * k)) & 0xff));
}

void medit_put_double(std::string& rOut, double V) {
    medit_put_int(rOut, static_cast<std::int64_t>(std::bit_cast<std::uint64_t>(V)), 8);
}

}  // namespace

Mesh read_medit_binary(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    const std::string_view buf = source.View();
    MeditBinCursor c{buf.data(), buf.size()};

    const std::int64_t code = c.Int(4);
    if (code == 16777216)
        c.mSwap = true;
    else if (code != 1)
        throw ReadError("Medit: invalid code");
    const std::int64_t version = c.Int(4);
    if (version < 1 || version > 4)
        throw ReadError("Medit: invalid version");
    const std::size_t isize = version == 4 ? 8 : 4;
    const std::size_t fsize = version == 1 ? 4 : 8;
    const std::size_t psize = version >= 3 ? 8 : 4;
    const std::int64_t field = c.Int(4);
    if (field != 3)
        throw ReadError("Medit: invalid dimension code : " + std::to_string(field) +
                        " it should be 3");
    c.Int(psize);
    const std::int64_t dim = c.Int(4);
    if (dim != 2 && dim != 3)
        throw ReadError("Medit: invalid mesh dimension : " + std::to_string(dim));

    Mesh mesh;
    bool have_points = false;
    std::vector<std::int64_t> point_ref;
    std::vector<std::pair<const MeditIdList*, std::vector<std::int64_t>>> id_lists;
    const DType fdtype = version == 1 ? DType::Float32 : DType::Float64;

    while (true) {
        if (c.Left() < 4) {
            log::warn("Medit: end-of-file reached before GmfEnd keyword");
            break;
        }
        const std::int64_t kw = c.Int(4);
        if (kw < 0 || kw > kMbMaxCode)
            throw ReadError("Medit: unsupported field " + std::to_string(kw));
        if (kw == kMbEnd)
            break;
        if (medit_is_reserved(static_cast<int>(kw)))
            continue;
        const std::int64_t next = c.Int(psize);
        const std::size_t section_start = c.mAt;

        const MeditBinCell* cell = nullptr;
        for (const MeditBinCell& k : kMeditBinCells)
            if (k.mCode == kw)
                cell = &k;
        const MeditBinIdList* idl = nullptr;
        for (const MeditBinIdList& k : kMeditBinIdLists)
            if (k.mCode == kw)
                idl = &k;

        if (kw == kMbVertices) {
            const std::size_t row = static_cast<std::size_t>(dim) * fsize + isize;
            const std::size_t n =
                detail::checked_count(c.Int(isize), c.Left() / row, "Medit", "vertex");
            NDArray pts(fdtype, {n, static_cast<std::size_t>(dim)});
            point_ref.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::int64_t d = 0; d < dim; ++d) {
                    const double v = c.Real(fsize);
                    if (fdtype == DType::Float32)
                        pts.As<float>()[i * dim + d] = static_cast<float>(v);
                    else
                        pts.As<double>()[i * dim + d] = v;
                }
                point_ref[i] = c.Int(isize);
            }
            mesh.AssignPoints(std::move(pts));
            have_points = true;
        } else if (cell) {
            const auto k = static_cast<std::size_t>(cell->mNodes);
            const std::size_t row = (k + 1) * isize;
            const std::size_t n =
                detail::checked_count(c.Int(isize), c.Left() / row, "Medit", "cell");
            NDArray data(DType::Int64, {n, k});
            NDArray ref(DType::Int64, {n});
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t j = 0; j < k; ++j)
                    data.As<std::int64_t>()[i * k + j] = c.Int(isize) - 1;
                ref.As<std::int64_t>()[i] = c.Int(isize);
            }
            mesh.AddCellBlock(cell->mType, std::move(data));
            mesh.AppendCellData("medit:ref", std::move(ref));
        } else if (idl) {
            const std::size_t n =
                detail::checked_count(c.Int(isize), c.Left() / isize, "Medit", "entry");
            std::vector<std::int64_t> ids(n);
            for (std::int64_t& id : ids)
                id = c.Int(isize) - 1;
            id_lists.emplace_back(medit_id_list(idl->mKeyword), std::move(ids));
        } else {
            // A section meshio++ has no home for: skipped through the position of the next one.
            if (next <= static_cast<std::int64_t>(section_start) ||
                static_cast<std::uint64_t>(next) > buf.size())
                throw ReadError("Medit: cannot skip keyword " + std::to_string(kw) +
                                ": its next-section position is outside the file");
            log::warn("meshio++ doesn't know Medit keyword {}. Skipping.", kw);
            c.mAt = static_cast<std::size_t>(next);
        }
    }

    if (!have_points)
        throw ReadError("Medit: expected Vertices");
    NDArray pr(DType::Int64, {point_ref.size()});
    for (std::size_t i = 0; i < point_ref.size(); ++i)
        pr.As<std::int64_t>()[i] = point_ref[i];
    mesh.AddPointData("medit:ref", std::move(pr));
    for (auto& [list, ids] : id_lists)
        medit_add_id_list_region(mesh, *list, std::move(ids));
    return mesh;
}

void write_medit_binary(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t d = rMesh.PointDim();
    if (d != 2 && d != 3)
        throw WriteError("Medit: a binary file holds 2-D or 3-D points, not " + std::to_string(d) +
                         "-D");
    const NDArray* vlabels = pick_first_int(rMesh);
    if (vlabels && rMesh.PointDataNames().size() > 1)
        log::warn("Medit can only write one point data array. Picking the first integer one.");
    const std::string clabel_key = pick_first_int_cell(rMesh);

    // Which blocks are written, and does everything fit 32-bit integers?
    struct Block {
        std::size_t mIndex;
        int mCode;
        std::size_t mNodes;
    };
    std::vector<Block> blocks;
    bool big = n > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    const auto fits = [](std::int64_t V) {
        return V >= std::numeric_limits<std::int32_t>::min() &&
               V <= std::numeric_limits<std::int32_t>::max();
    };
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const MeditBinCell* found = nullptr;
        for (const MeditBinCell& k : kMeditBinCells)
            if (cb.Type() == k.mType)
                found = &k;
        if (!found) {
            log::warn("MEDIT's mesh format doesn't know {} cells. Skipping.",
                      std::string(cb.Type()));
            continue;
        }
        blocks.push_back({b, found->mCode, static_cast<std::size_t>(found->mNodes)});
        const detail::Int64View conn(cb.Conn());
        for (std::size_t i = 0; i < cb.NumCells() * cb.NodesPerCell() && !big; ++i)
            big = !fits(conn[i] + 1);
        if (!clabel_key.empty() && b < rMesh.CellDataNumBlocks(clabel_key)) {
            const NDArray& lab = rMesh.CellData(clabel_key, b);
            const detail::Int64View lv(lab);
            for (std::size_t i = 0; i < lab.Size() && !big; ++i)
                big = !fits(lv[i]);
        }
    }
    std::optional<detail::Int64View> plabels;
    if (vlabels) {
        plabels.emplace(*vlabels);
        for (std::size_t i = 0; i < n && !big; ++i)
            big = !fits((*plabels)[i]);
    }
    const auto sections = medit_id_list_sections(rMesh);
    for (const auto& [list, ids] : sections)
        for (std::size_t i = 0; i < ids.size() && !big; ++i)
            big = !fits(ids[i] + 1);

    const std::size_t isize = big ? 8 : 4;
    const std::size_t psize = 8;
    std::string out;
    medit_put_int(out, 1, 4);
    medit_put_int(out, big ? 4 : 3, 4);
    medit_put_int(out, 3, 4);  // GmfDimension
    std::int64_t pos = 4 * 4 + static_cast<std::int64_t>(psize);
    medit_put_int(out, pos, psize);
    medit_put_int(out, static_cast<std::int64_t>(d), 4);

    const auto header = [&](int Code, std::int64_t Count, std::int64_t Bytes) {
        pos += Bytes + 4 + static_cast<std::int64_t>(psize + isize);
        medit_put_int(out, Code, 4);
        medit_put_int(out, pos, psize);
        medit_put_int(out, Count, isize);
    };

    header(kMbVertices, static_cast<std::int64_t>(n),
           static_cast<std::int64_t>(n * (d * 8 + isize)));
    {
        const detail::DoubleView pv(rMesh.Points());
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t c = 0; c < d; ++c)
                medit_put_double(out, pv[i * d + c]);
            medit_put_int(out, plabels ? (*plabels)[i] : 1, isize);
        }
    }
    for (const Block& blk : blocks) {
        const auto cb = rMesh.Cells(blk.mIndex);
        const std::size_t count = cb.NumCells();
        header(blk.mCode, static_cast<std::int64_t>(count),
               static_cast<std::int64_t>(count * (blk.mNodes + 1) * isize));
        const detail::Int64View conn(cb.Conn());
        std::optional<detail::Int64View> labels;
        if (!clabel_key.empty() && blk.mIndex < rMesh.CellDataNumBlocks(clabel_key))
            labels.emplace(rMesh.CellData(clabel_key, blk.mIndex));
        for (std::size_t r = 0; r < count; ++r) {
            for (std::size_t j = 0; j < blk.mNodes; ++j)
                medit_put_int(out, conn[r * blk.mNodes + j] + 1, isize);
            medit_put_int(out, labels ? (*labels)[r] : 1, isize);
        }
    }
    for (const auto& [list, ids] : sections) {
        int code = 0;
        for (const MeditBinIdList& k : kMeditBinIdLists)
            if (std::string(k.mKeyword) == list->mKeyword)
                code = k.mCode;
        header(code, static_cast<std::int64_t>(ids.size()),
               static_cast<std::int64_t>(ids.size() * isize));
        for (const std::int64_t id : ids)
            medit_put_int(out, id + 1, isize);
    }
    medit_put_int(out, kMbEnd, 4);
    medit_put_int(out, 0, psize);

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    os.close();
    if (!os)
        throw WriteError("Medit: failed writing " + rPath);
}

Mesh read_medit(const std::string& rPath) {
    return !rPath.empty() && rPath.back() == 'b' ? read_medit_binary(rPath)
                                                 : read_medit_ascii(rPath);
}

void write_medit(const std::string& rPath, const Mesh& rMesh) {
    if (!rPath.empty() && rPath.back() == 'b')
        write_medit_binary(rPath, rMesh);
    else
        write_medit_ascii(rPath, rMesh);
}

}  // namespace meshioplusplus
