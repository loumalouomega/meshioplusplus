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
#include <array>
#include <bit>
#include <cctype>
#include <cmath>
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
#include "meshioplusplus/ndarray.hpp"
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
        {"Edges", {"line", 2}},
        {"Triangles", {"triangle", 3}},
        {"Quadrilaterals", {"quad", 4}},
        {"Tetrahedra", {"tetra", 4}},
        {"Prisms", {"wedge", 6}},
        {"Pyramids", {"pyramid", 5}},
        {"Hexahedra", {"hexahedron", 8}},
        {"Hexaedra", {"hexahedron", 8}},
        {"EdgesP2", {"line3", 3}},
        {"EdgesP3", {"line4", 4}},
        {"TrianglesP2", {"triangle6", 6}},
        {"TrianglesP3", {"triangle10", 10}},
        {"TetrahedraP2", {"tetra10", 10}},
    };
    return m;
}

// meshio type -> (medit keyword, nodes per cell), write order.
const std::vector<std::pair<std::string, std::pair<std::string, int>>>& meshio_to_medit() {
    static const std::vector<std::pair<std::string, std::pair<std::string, int>>> m = {
        {"line", {"Edges", 2}},
        {"triangle", {"Triangles", 3}},
        {"quad", {"Quadrilaterals", 4}},
        {"tetra", {"Tetrahedra", 4}},
        {"wedge", {"Prisms", 6}},
        {"pyramid", {"Pyramids", 5}},
        {"hexahedron", {"Hexahedra", 8}},
        // High-order sections, each preceded by its ordering (see MeditHighOrder).
        {"line3", {"EdgesP2", 3}},
        {"line4", {"EdgesP3", 4}},
        {"triangle6", {"TrianglesP2", 6}},
        {"triangle10", {"TrianglesP3", 10}},
        {"tetra10", {"TetrahedraP2", 10}},
    };
    return m;
}

/**
 * libMeshb's high-order elements. The format defines no node order for them: a
 * file states its own in a `<Keyword>Ordering` section, one integer tuple per
 * node, and a program maps it onto its own order by matching tuples
 * (`GmfSetHONodesOrdering`). meshio++'s tuples, in meshio's node order: for a
 * simplex the barycentric multiplicities, component k counting vertex k; for
 * an edge the node's position from the first vertex. The writers emit them;
 * the readers match a file's against them, and read a section with no
 * ordering in meshio's order (libMeshb's own example writes TrianglesP2 so).
 */
struct MeditHighOrder {
    const char* mType;
    const char* mKeyword;
    int mCode;
    const char* mOrderingKeyword;
    int mOrderingCode;
    int mNodes;
    int mCoords;
    std::array<int, 40> mTuples;  ///< mNodes * mCoords, meshio order
};

const MeditHighOrder kMeditHighOrder[] = {
    {"line3", "EdgesP2", 25, "EdgesP2Ordering", 127, 3, 1, {0, 2, 1}},
    {"line4", "EdgesP3", 92, "EdgesP3Ordering", 128, 4, 1, {0, 3, 1, 2}},
    {"triangle6",
     "TrianglesP2",
     24,
     "TrianglesP2Ordering",
     129,
     6,
     3,
     {2, 0, 0, 0, 2, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 1, 0, 1}},
    {"triangle10", "TrianglesP3", 90, "TrianglesP3Ordering", 130, 10, 3, {3, 0, 0, 0, 3, 0, 0, 0,
                                                                          3, 2, 1, 0, 1, 2, 0, 0,
                                                                          2, 1, 0, 1, 2, 1, 0, 2,
                                                                          2, 0, 1, 1, 1, 1}},
    {"tetra10", "TetrahedraP2", 30, "TetrahedraP2Ordering", 133, 10, 4, {2, 0, 0, 0, 0, 2, 0, 0,
                                                                         0, 0, 2, 0, 0, 0, 0, 2,
                                                                         1, 1, 0, 0, 0, 1, 1, 0,
                                                                         1, 0, 1, 0, 1, 0, 0, 1,
                                                                         0, 1, 0, 1, 0, 0, 1, 1}},
};

const MeditHighOrder* medit_high_order_of_type(std::string_view Type) {
    for (const MeditHighOrder& h : kMeditHighOrder)
        if (Type == h.mType)
            return &h;
    return nullptr;
}

const MeditHighOrder* medit_high_order_of_ordering(std::string_view Keyword) {
    for (const MeditHighOrder& h : kMeditHighOrder)
        if (Keyword == h.mOrderingKeyword)
            return &h;
    return nullptr;
}

/// The node permutation (file position -> meshio position) an ordering
/// section's tuples give; a table that is not a permutation of meshio++'s
/// tuples (a tuple it does not have, one listed twice) is refused rather than
/// guessed at.
std::vector<int> medit_ordering_permutation(const MeditHighOrder& rH,
                                            const std::vector<std::int64_t>& rTuples) {
    const auto nodes = static_cast<std::size_t>(rH.mNodes);
    const auto coords = static_cast<std::size_t>(rH.mCoords);
    std::vector<int> perm(nodes, -1);
    std::vector<bool> used(nodes, false);
    for (std::size_t f = 0; f < nodes; ++f) {
        for (std::size_t m = 0; m < nodes; ++m) {
            bool same = true;
            for (std::size_t k = 0; k < coords && same; ++k)
                same = rTuples[f * coords + k] == rH.mTuples[m * coords + k];
            if (same) {
                if (used[m])
                    throw ReadError(std::string("Medit: ") + rH.mOrderingKeyword +
                                    " lists a node twice");
                used[m] = true;
                perm[f] = static_cast<int>(m);
                break;
            }
        }
        if (perm[f] < 0)
            throw ReadError(std::string("Medit: ") + rH.mOrderingKeyword + " lists a node no " +
                            rH.mKeyword + " element has");
    }
    return perm;
}

/// A cell section as read, held back until the file's ordering sections
/// (which may come after it) are known.
struct MeditStagedBlock {
    std::string mType;
    NDArray mConn;
    NDArray mRef;
};

/// Puts the high-order blocks' nodes in meshio's order by the file's ordering
/// sections (keyed by meshio type), then adds every block to the mesh in file order.
void medit_add_blocks(Mesh& rMesh, std::vector<MeditStagedBlock>& rBlocks,
                      const std::unordered_map<std::string, std::vector<int>>& rPerms) {
    for (MeditStagedBlock& blk : rBlocks) {
        const auto it = rPerms.find(blk.mType);
        if (it != rPerms.end()) {
            const std::vector<int>& perm = it->second;
            const std::size_t k = perm.size();
            std::int64_t* conn = blk.mConn.As<std::int64_t>();
            const std::size_t rows = blk.mConn.Shape()[0];
            std::vector<std::int64_t> row(k);
            for (std::size_t r = 0; r < rows; ++r) {
                for (std::size_t j = 0; j < k; ++j)
                    row[static_cast<std::size_t>(perm[j])] = conn[r * k + j];
                std::copy(row.begin(), row.end(), conn + r * k);
            }
        }
        rMesh.AddCellBlock(blk.mType, std::move(blk.mConn));
        rMesh.AppendCellData("medit:ref", std::move(blk.mRef));
    }
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

// A solution keyword (`SolAt*`): the libMeshb code and the cell type its rows
// belong to (null: one row per vertex).
struct MeditSolKeyword {
    const char* mKeyword;
    int mCode;
    const char* mCellType;
};

const MeditSolKeyword kMeditSolKeywords[] = {
    {"SolAtVertices", 62, nullptr},       {"SolAtEdges", 63, "line"},
    {"SolAtTriangles", 64, "triangle"},   {"SolAtQuadrilaterals", 65, "quad"},
    {"SolAtTetrahedra", 66, "tetra"},     {"SolAtPrisms", 67, "wedge"},
    {"SolAtHexahedra", 68, "hexahedron"}, {"SolAtPyramids", 26, "pyramid"},
};

const MeditSolKeyword* medit_sol_keyword(std::string_view Keyword) {
    for (const MeditSolKeyword& k : kMeditSolKeywords)
        if (Keyword == k.mKeyword)
            return &k;
    return nullptr;
}

/// Values per row of one field of libMeshb solution type `Type`: a scalar, a
/// vector, a symmetric matrix (lower triangle) or a full matrix.
std::size_t medit_sol_components(std::int64_t Type, std::int64_t Dim) {
    switch (Type) {
        case 1:
            return 1;
        case 2:
            return static_cast<std::size_t>(Dim);
        case 3:
            return static_cast<std::size_t>(Dim * (Dim + 1) / 2);
        case 4:
            return static_cast<std::size_t>(Dim * Dim);
        default:
            throw ReadError("Medit: unknown solution type " + std::to_string(Type));
    }
}

/// The solution type holding `Components` values in `Dim`-D, or 0.
int medit_sol_type_of(std::size_t Components, std::size_t Dim) {
    for (int t = 1; t <= 4; ++t)
        if (medit_sol_components(t, static_cast<std::int64_t>(Dim)) == Components)
            return t;
    return 0;
}

/// One `SolAt*` section: its fields' types and `mRows x mTotal` values.
struct MeditSolSection {
    const MeditSolKeyword* mKeyword = nullptr;
    std::vector<int> mTypes;
    std::size_t mRows = 0;
    std::size_t mTotal = 0;
    DType mDtype = DType::Float64;
    std::vector<double> mValues;
};

std::size_t medit_sol_total(const std::vector<int>& rTypes, std::int64_t Dim) {
    if (rTypes.empty())
        throw ReadError("Medit: a solution section declares no field");
    std::size_t total = 0;
    for (const int t : rTypes)
        total += medit_sol_components(t, Dim);
    return total;
}

/// An ascii `SolAt*` section, its keyword already consumed.
MeditSolSection medit_sol_read_ascii_section(Tokenizer& rTok, const MeditSolKeyword& rKeyword,
                                             int Dim, DType Dtype) {
    if (Dim != 2 && Dim != 3)
        throw ReadError(std::string("Medit: ") + rKeyword.mKeyword + " before Dimension");
    MeditSolSection section;
    section.mKeyword = &rKeyword;
    section.mDtype = Dtype;
    section.mRows = static_cast<std::size_t>(rTok.next_count());
    const std::int64_t num_types = rTok.next_count();
    for (std::int64_t i = 0; i < num_types; ++i)
        section.mTypes.push_back(static_cast<int>(rTok.next_int()));
    section.mTotal = medit_sol_total(section.mTypes, Dim);
    // Every value takes at least two bytes of the file.
    if (section.mRows > rTok.mBuf.size() / section.mTotal)
        throw ReadError(std::string("Medit: ") + rKeyword.mKeyword + " is truncated");
    section.mValues.resize(section.mRows * section.mTotal);
    for (double& v : section.mValues) {
        const std::string_view token = rTok.next();
        if (token.empty())
            throw ReadError(std::string("Medit: ") + rKeyword.mKeyword + " is truncated");
        v = detail::parse_double_prefix(token);
    }
    return section;
}

/// The section's fields: `(rows,)` for a scalar, else `(rows, components)`.
std::vector<NDArray> medit_sol_fields(const MeditSolSection& rSection, int Dim) {
    std::vector<NDArray> out;
    std::size_t offset = 0;
    for (const int t : rSection.mTypes) {
        const std::size_t c = medit_sol_components(t, Dim);
        NDArray field = c == 1 ? NDArray(rSection.mDtype, {rSection.mRows})
                               : NDArray(rSection.mDtype, {rSection.mRows, c});
        for (std::size_t r = 0; r < rSection.mRows; ++r)
            for (std::size_t j = 0; j < c; ++j)
                store_coord(field, r * c + j, rSection.mValues[r * rSection.mTotal + offset + j]);
        out.push_back(std::move(field));
        offset += c;
    }
    return out;
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
    std::vector<MeditStagedBlock> blocks;
    std::unordered_map<std::string, std::vector<int>> perms;

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
            blocks.push_back({type, std::move(data), std::move(ref)});
        } else if (const MeditHighOrder* ho = medit_high_order_of_ordering(kw)) {
            // One tuple per node; anything else is not this element's ordering.
            if (tok.next_count() != ho->mNodes)
                throw ReadError(std::string("Medit: ") + ho->mOrderingKeyword + " lists " +
                                "a number of nodes other than " + std::to_string(ho->mNodes));
            std::vector<std::int64_t> tuples(static_cast<std::size_t>(ho->mNodes * ho->mCoords));
            for (std::int64_t& v : tuples)
                v = tok.next_int();
            perms[ho->mType] = medit_ordering_permutation(*ho, tuples);
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
        } else if (const MeditSolKeyword* sol = medit_sol_keyword(kw)) {
            // A solution inside the mesh file: its home is a `.sol` beside it.
            medit_sol_read_ascii_section(tok, *sol, dim, coord_dtype);
            log::warn(
                "Medit: {} in a mesh file is skipped; keep the solution in a .sol beside the "
                "mesh and read it with medit.read(path, sol=...)",
                kw);
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
    medit_add_blocks(mesh, blocks, perms);

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

    // Cells, grouped by medit element keyword; a type with none is said so.
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        const std::string_view type = rMesh.Cells(ci).Type();
        if (std::none_of(meshio_to_medit().begin(), meshio_to_medit().end(),
                         [&](const auto& rMk) { return rMk.first == type; }))
            log::warn("MEDIT's mesh format doesn't know {} cells. Skipping.", std::string(type));
    }
    const std::string clabel_key = pick_first_int_cell(rMesh);
    for (const auto& mk : meshio_to_medit()) {
        const std::string& mtype = mk.first;
        const std::string& kw = mk.second.first;
        int k = mk.second.second;
        bool ordered = false;
        for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
            const auto cb = rMesh.Cells(ci);
            if (cb.Type() != mtype)
                continue;
            if (const MeditHighOrder* ho = medit_high_order_of_type(mtype); ho && !ordered) {
                ordered = true;
                os << "\n" << ho->mOrderingKeyword << "\n" << ho->mNodes << "\n";
                for (int i = 0; i < ho->mNodes; ++i) {
                    for (int j = 0; j < ho->mCoords; ++j)
                        os << (j ? " " : "")
                           << ho->mTuples[static_cast<std::size_t>(i * ho->mCoords + j)];
                    os << "\n";
                }
            }
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
    {5, "line", 2},   {6, "triangle", 3},    {7, "quad", 4},         {8, "tetra", 4},
    {9, "wedge", 6},  {10, "hexahedron", 8}, {49, "pyramid", 5},     {25, "line3", 3},
    {92, "line4", 4}, {24, "triangle6", 6},  {90, "triangle10", 10}, {30, "tetra10", 10},
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
    std::vector<MeditStagedBlock> blocks;
    std::unordered_map<std::string, std::vector<int>> perms;
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
        const MeditHighOrder* ordering = nullptr;
        for (const MeditHighOrder& h : kMeditHighOrder)
            if (h.mOrderingCode == kw)
                ordering = &h;

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
            blocks.push_back({cell->mType, std::move(data), std::move(ref)});
        } else if (ordering) {
            // One tuple per node; anything else is not this element's ordering.
            if (c.Int(isize) != ordering->mNodes)
                throw ReadError(std::string("Medit: ") + ordering->mOrderingKeyword +
                                " lists a number of nodes other than " +
                                std::to_string(ordering->mNodes));
            std::vector<std::int64_t> tuples(
                static_cast<std::size_t>(ordering->mNodes * ordering->mCoords));
            for (std::int64_t& v : tuples)
                v = c.Int(isize);
            perms[ordering->mType] = medit_ordering_permutation(*ordering, tuples);
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
    medit_add_blocks(mesh, blocks, perms);
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
    std::vector<const MeditHighOrder*> ordered;
    for (const Block& blk : blocks) {
        const auto cb = rMesh.Cells(blk.mIndex);
        // A high-order type's ordering, once, before its first section.
        if (const MeditHighOrder* ho = medit_high_order_of_type(cb.Type());
            ho && std::find(ordered.begin(), ordered.end(), ho) == ordered.end()) {
            ordered.push_back(ho);
            const auto values = static_cast<std::size_t>(ho->mNodes * ho->mCoords);
            header(ho->mOrderingCode, ho->mNodes, static_cast<std::int64_t>(values * isize));
            for (std::size_t i = 0; i < values; ++i)
                medit_put_int(out, ho->mTuples[i], isize);
        }
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

// --- solutions (.sol ascii, .solb binary) ----------------------------------------------------

namespace {

using MeditSolFile = std::pair<int, std::vector<MeditSolSection>>;

MeditSolFile medit_sol_read_ascii(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    Tokenizer tok(source.View());
    int dim = 0;
    DType dtype = DType::Float64;
    std::vector<MeditSolSection> sections;
    while (!tok.eof()) {
        const std::string kw(tok.next());
        if (kw.empty())
            break;
        if (kw == "MeshVersionFormatted") {
            dtype = tok.next_int() <= 1 ? DType::Float32 : DType::Float64;
        } else if (kw == "Dimension") {
            dim = static_cast<int>(tok.next_int());
            if (dim != 2 && dim != 3)
                throw ReadError("Medit: invalid solution dimension " + std::to_string(dim));
        } else if (const MeditSolKeyword* sol = medit_sol_keyword(kw)) {
            sections.push_back(medit_sol_read_ascii_section(tok, *sol, dim, dtype));
        } else if (kw == "Time" || kw == "Iterations") {
            tok.next();
        } else if (kw == "End") {
            break;
        } else {
            throw ReadError("Medit: unknown solution keyword '" + kw + "'");
        }
    }
    if (dim == 0)
        throw ReadError("Medit: expected Dimension");
    return {dim, std::move(sections)};
}

MeditSolFile medit_sol_read_binary(const std::string& rPath) {
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
    if (c.Int(4) != 3)
        throw ReadError("Medit: expected Dimension");
    c.Int(psize);
    const std::int64_t dim = c.Int(4);
    if (dim != 2 && dim != 3)
        throw ReadError("Medit: invalid solution dimension " + std::to_string(dim));

    std::vector<MeditSolSection> sections;
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
        const MeditSolKeyword* sol = nullptr;
        for (const MeditSolKeyword& k : kMeditSolKeywords)
            if (k.mCode == kw)
                sol = &k;
        if (sol) {
            MeditSolSection section;
            section.mKeyword = sol;
            section.mDtype = version == 1 ? DType::Float32 : DType::Float64;
            const std::int64_t rows = c.Int(isize);
            const std::int64_t num_types = c.Int(4);
            if (rows < 0 || num_types < 0 || static_cast<std::uint64_t>(num_types) > c.Left() / 4)
                throw ReadError(std::string("Medit: ") + sol->mKeyword + " has an impossible size");
            for (std::int64_t i = 0; i < num_types; ++i)
                section.mTypes.push_back(static_cast<int>(c.Int(4)));
            section.mTotal = medit_sol_total(section.mTypes, dim);
            section.mRows = detail::checked_count(rows, c.Left() / (section.mTotal * fsize),
                                                  "Medit", "solution row");
            section.mValues.resize(section.mRows * section.mTotal);
            for (double& v : section.mValues)
                v = c.Real(fsize);
            sections.push_back(std::move(section));
        } else {
            if (next <= static_cast<std::int64_t>(section_start) ||
                static_cast<std::uint64_t>(next) > buf.size())
                throw ReadError("Medit: cannot skip keyword " + std::to_string(kw) +
                                ": its next-section position is outside the file");
            log::warn("meshio++ doesn't know Medit keyword {}. Skipping.", kw);
            c.mAt = static_cast<std::size_t>(next);
        }
    }
    return {static_cast<int>(dim), std::move(sections)};
}

MeditSolFile medit_sol_read_file(const std::string& rPath) {
    return !rPath.empty() && rPath.back() == 'b' ? medit_sol_read_binary(rPath)
                                                 : medit_sol_read_ascii(rPath);
}

/// The vertex fields, named `medit:sol:<i>` in file order, and their row count.
std::vector<NDArray> medit_sol_vertex_fields(const MeditSolFile& rFile,
                                             std::optional<std::size_t>& rCount) {
    std::vector<NDArray> out;
    for (const MeditSolSection& s : rFile.second) {
        if (s.mKeyword->mCellType)
            continue;
        if (rCount && *rCount != s.mRows)
            throw ReadError("Medit: two SolAtVertices sections of different sizes");
        rCount = s.mRows;
        for (NDArray& f : medit_sol_fields(s, rFile.first))
            out.push_back(std::move(f));
    }
    return out;
}

/// The trailing shape of a field: what follows its row count.
std::vector<std::size_t> medit_sol_trailing(const NDArray& rA) {
    return std::vector<std::size_t>(rA.Shape().begin() + 1, rA.Shape().end());
}

NDArray medit_sol_rows(const NDArray& rField, std::size_t First, std::size_t Count) {
    std::vector<std::size_t> shape = rField.Shape();
    shape[0] = Count;
    NDArray out(rField.Dtype(), shape);
    const std::size_t width = rField.Size() / std::max<std::size_t>(rField.Shape()[0], 1);
    const std::size_t bytes = rField.Dtype() == DType::Float32 ? 4 : 8;
    std::memcpy(out.Data(), rField.Data() + First * width * bytes, Count * width * bytes);
    return out;
}

void medit_sol_attach(Mesh& rMesh, const MeditSolFile& rFile) {
    const int dim = rFile.first;
    if (rMesh.PointDim() != static_cast<std::size_t>(dim))
        throw ReadError("Medit: the solution is " + std::to_string(dim) + "-D and the mesh " +
                        std::to_string(rMesh.PointDim()) + "-D");
    std::optional<std::size_t> count;
    std::vector<NDArray> point_fields = medit_sol_vertex_fields(rFile, count);
    if (count && *count != rMesh.NumPoints())
        throw ReadError("Medit: SolAtVertices has " + std::to_string(*count) +
                        " rows and the mesh " + std::to_string(rMesh.NumPoints()) + " points");
    for (std::size_t i = 0; i < point_fields.size(); ++i)
        rMesh.AddPointData("medit:sol:" + std::to_string(i), std::move(point_fields[i]));

    const std::size_t num_blocks = rMesh.NumCellBlocks();
    // name -> one array per block (absent: the fill), in the order first met
    std::vector<std::pair<std::string, std::vector<std::optional<NDArray>>>> cell_fields;
    for (const MeditSolSection& s : rFile.second) {
        if (!s.mKeyword->mCellType)
            continue;
        std::vector<std::size_t> blocks;
        std::size_t total = 0;
        for (std::size_t b = 0; b < num_blocks; ++b) {
            const auto cb = rMesh.Cells(b);
            if (cb.Type() == s.mKeyword->mCellType) {
                blocks.push_back(b);
                total += cb.NumCells();
            }
        }
        if (total != s.mRows)
            throw ReadError(std::string("Medit: ") + s.mKeyword->mKeyword + " has " +
                            std::to_string(s.mRows) + " rows and the mesh " +
                            std::to_string(total) + " " + s.mKeyword->mCellType + " cells");
        std::vector<NDArray> fields = medit_sol_fields(s, dim);
        for (std::size_t j = 0; j < fields.size(); ++j) {
            std::string name = "medit:sol:" + std::to_string(j);
            auto find = [&](const std::string& rName) {
                return std::find_if(cell_fields.begin(), cell_fields.end(),
                                    [&](const auto& rEntry) { return rEntry.first == rName; });
            };
            auto known = find(name);
            if (known != cell_fields.end()) {
                const auto other = std::find_if(known->second.begin(), known->second.end(),
                                                [](const auto& rA) { return rA.has_value(); });
                if (medit_sol_trailing(**other) != medit_sol_trailing(fields[j]))
                    name += std::string(":") + s.mKeyword->mCellType;
            }
            known = find(name);
            if (known == cell_fields.end()) {
                cell_fields.emplace_back(name, std::vector<std::optional<NDArray>>(num_blocks));
                known = cell_fields.end() - 1;
            }
            std::size_t start = 0;
            for (const std::size_t b : blocks) {
                const std::size_t n = rMesh.Cells(b).NumCells();
                known->second[b] = medit_sol_rows(fields[j], start, n);
                start += n;
            }
        }
    }
    for (auto& [name, per_block] : cell_fields) {
        // The fill's dtype and trailing shape, taken before the arrays move out.
        const NDArray& sample = **std::find_if(per_block.begin(), per_block.end(),
                                               [](const auto& rA) { return rA.has_value(); });
        const DType dtype = sample.Dtype();
        const std::vector<std::size_t> sample_shape = sample.Shape();
        std::vector<NDArray> arrays;
        for (std::size_t b = 0; b < num_blocks; ++b) {
            if (per_block[b]) {
                arrays.push_back(std::move(*per_block[b]));
                continue;
            }
            std::vector<std::size_t> shape = sample_shape;
            shape[0] = rMesh.Cells(b).NumCells();
            NDArray fill(dtype, shape);
            for (std::size_t i = 0; i < fill.Size(); ++i)
                store_coord(fill, i, std::numeric_limits<double>::quiet_NaN());
            arrays.push_back(std::move(fill));
        }
        rMesh.AddCellData(name, std::move(arrays));
    }
}

std::size_t medit_sol_dimension(const Mesh& rMesh) {
    const std::size_t d = rMesh.PointDim();
    if (d == 2 || d == 3)
        return d;
    if (rMesh.HasFieldData("medit:dimension")) {
        const NDArray& f = rMesh.FieldData("medit:dimension");
        if (f.Size() > 0) {
            const double v = detail::DoubleView(f)[0];
            if (v == 2 || v == 3)
                return static_cast<std::size_t>(v);
        }
    }
    throw WriteError(
        "Medit: a solution is 2-D or 3-D; give the mesh 2-D or 3-D points or "
        "field_data['medit:dimension']");
}

/// The solution type of one array and its row width, or type 0 (with a warning).
std::pair<int, std::size_t> medit_sol_classify(const std::string& rName, const NDArray& rA,
                                               std::size_t Dim) {
    if (!detail::is_float_dtype(rA.Dtype()) || rA.Ndim() == 0) {
        if (rName != "medit:ref")
            log::warn("Medit: a solution holds reals; '{}' is not written", rName);
        return {0, 0};
    }
    const std::size_t width = rA.Shape()[0] ? rA.Size() / rA.Shape()[0] : 1;
    const int type = medit_sol_type_of(width, Dim);
    if (!type)
        log::warn("Medit: '{}' has {} components, which no {}-D solution type holds; not written",
                  rName, width, Dim);
    return {type, width};
}

/// Append the columns of `rA` (rows x width) to a section being built.
void medit_sol_append(MeditSolSection& rSection, int Type, std::size_t Width,
                      const std::vector<const NDArray*>& rParts) {
    std::size_t rows = 0;
    for (const NDArray* p : rParts)
        rows += Width ? p->Size() / Width : 0;
    if (rSection.mTypes.empty())
        rSection.mRows = rows;
    std::vector<double> merged(rows * (rSection.mTotal + Width));
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t j = 0; j < rSection.mTotal; ++j)
            merged[r * (rSection.mTotal + Width) + j] = rSection.mValues[r * rSection.mTotal + j];
    std::size_t r = 0;
    for (const NDArray* p : rParts) {
        const detail::DoubleView v(*p);
        const std::size_t n = Width ? p->Size() / Width : 0;
        for (std::size_t i = 0; i < n; ++i, ++r)
            for (std::size_t j = 0; j < Width; ++j)
                merged[r * (rSection.mTotal + Width) + rSection.mTotal + j] = v[i * Width + j];
    }
    rSection.mValues = std::move(merged);
    rSection.mTotal += Width;
    rSection.mTypes.push_back(Type);
}

std::vector<MeditSolSection> medit_sol_sections_of(const Mesh& rMesh, std::size_t Dim) {
    std::vector<MeditSolSection> sections;
    MeditSolSection vertices;
    vertices.mKeyword = &kMeditSolKeywords[0];
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const auto [type, width] = medit_sol_classify(name, a, Dim);
        if (type)
            medit_sol_append(vertices, type, width, {&a});
    }
    if (!vertices.mTypes.empty())
        sections.push_back(std::move(vertices));

    const std::size_t num_blocks = rMesh.NumCellBlocks();
    std::vector<std::string> names;
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name == "medit:ref" || rMesh.CellDataNumBlocks(name) == 0)
            continue;
        bool real = true;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            real = real && detail::is_float_dtype(rMesh.CellData(name, b).Dtype());
        if (!real) {
            log::warn("Medit: a solution holds reals; '{}' is not written", name);
            continue;
        }
        names.push_back(name);
        for (std::size_t b = 0; b < num_blocks && b < rMesh.CellDataNumBlocks(name); ++b) {
            const std::string type(rMesh.Cells(b).Type());
            if (std::any_of(std::begin(kMeditSolKeywords), std::end(kMeditSolKeywords),
                            [&](const MeditSolKeyword& rK) {
                                return rK.mCellType && type == rK.mCellType;
                            }))
                continue;
            const NDArray& a = rMesh.CellData(name, b);
            const detail::DoubleView v(a);
            bool all_nan = true;
            for (std::size_t i = 0; i < a.Size() && all_nan; ++i)
                all_nan = std::isnan(v[i]);
            if (!all_nan)
                log::warn("Medit: a solution has no {} section; '{}' skips it", type, name);
        }
    }
    for (const MeditSolKeyword& k : kMeditSolKeywords) {
        if (!k.mCellType)
            continue;
        std::vector<std::size_t> blocks;
        for (std::size_t b = 0; b < num_blocks; ++b)
            if (rMesh.Cells(b).Type() == k.mCellType)
                blocks.push_back(b);
        if (blocks.empty())
            continue;
        MeditSolSection section;
        section.mKeyword = &k;
        for (const std::string& name : names) {
            std::vector<const NDArray*> parts;
            std::size_t size = 0;
            bool all_nan = true;
            for (const std::size_t b : blocks) {
                if (b >= rMesh.CellDataNumBlocks(name))
                    throw WriteError("Medit: cell data '" + name + "' misses a block");
                const NDArray& a = rMesh.CellData(name, b);
                parts.push_back(&a);
                size += a.Size();
                const detail::DoubleView v(a);
                for (std::size_t i = 0; i < a.Size() && all_nan; ++i)
                    all_nan = std::isnan(v[i]);
            }
            if (size && all_nan)
                continue;  // the fill of a field other cell types carry
            const auto [type, width] = medit_sol_classify(name, *parts[0], Dim);
            if (type)
                medit_sol_append(section, type, width, parts);
        }
        if (!section.mTypes.empty())
            sections.push_back(std::move(section));
    }
    if (sections.empty())
        throw WriteError("Medit: no floating-point point or cell data to write as a solution");
    return sections;
}

void medit_sol_write_ascii(const std::string& rPath, std::size_t Dim,
                           const std::vector<MeditSolSection>& rSections) {
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os << "MeshVersionFormatted 2\nDimension " << Dim << "\n";
    const detail::CNumber num;
    for (const MeditSolSection& s : rSections) {
        os << "\n" << s.mKeyword->mKeyword << "\n" << s.mRows << "\n" << s.mTypes.size();
        for (const int t : s.mTypes)
            os << " " << t;
        os << "\n";
        detail::write_row_chunks(
            os, s.mRows, [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t r = First; r < Last; ++r) {
                    for (std::size_t j = 0; j < s.mTotal; ++j)
                        num.Append(rBuf, j ? " %.16e" : "%.16e", s.mValues[r * s.mTotal + j]);
                    rBuf += '\n';
                }
            });
    }
    os << "\nEnd\n";
    os.close();
    if (!os)
        throw WriteError("Medit: failed writing " + rPath);
}

void medit_sol_write_binary(const std::string& rPath, std::size_t Dim,
                            const std::vector<MeditSolSection>& rSections) {
    bool big = false;
    for (const MeditSolSection& s : rSections)
        big = big || s.mRows > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    const std::size_t isize = big ? 8 : 4;
    const std::size_t psize = 8;
    std::string out;
    medit_put_int(out, 1, 4);
    medit_put_int(out, big ? 4 : 3, 4);
    medit_put_int(out, 3, 4);  // GmfDimension
    std::int64_t pos = 4 * 4 + static_cast<std::int64_t>(psize);
    medit_put_int(out, pos, psize);
    medit_put_int(out, static_cast<std::int64_t>(Dim), 4);
    for (const MeditSolSection& s : rSections) {
        const std::size_t bytes = 4 + 4 * s.mTypes.size() + s.mRows * s.mTotal * 8;
        pos += static_cast<std::int64_t>(bytes + 4 + psize + isize);
        medit_put_int(out, s.mKeyword->mCode, 4);
        medit_put_int(out, pos, psize);
        medit_put_int(out, static_cast<std::int64_t>(s.mRows), isize);
        medit_put_int(out, static_cast<std::int64_t>(s.mTypes.size()), 4);
        for (const int t : s.mTypes)
            medit_put_int(out, t, 4);
        for (const double v : s.mValues)
            medit_put_double(out, v);
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

}  // namespace

Mesh read_medit_sol(const std::string& rPath) {
    const MeditSolFile file = medit_sol_read_file(rPath);
    std::optional<std::size_t> count;
    std::vector<NDArray> fields = medit_sol_vertex_fields(file, count);
    for (const MeditSolSection& s : file.second)
        if (s.mKeyword->mCellType)
            log::warn(
                "Medit: {} needs the mesh it belongs to; read it with medit.read(mesh_path, "
                "sol=...). Skipping.",
                s.mKeyword->mKeyword);
    Mesh mesh;
    mesh.AssignPoints(NDArray(DType::Float64, {count.value_or(0), 0}));
    for (std::size_t i = 0; i < fields.size(); ++i)
        mesh.AddPointData("medit:sol:" + std::to_string(i), std::move(fields[i]));
    NDArray dimension(DType::Int64, {1});
    dimension.As<std::int64_t>()[0] = file.first;
    mesh.AddFieldData("medit:dimension", std::move(dimension));
    return mesh;
}

void read_medit_sol_into(const std::string& rPath, Mesh& rMesh) {
    medit_sol_attach(rMesh, medit_sol_read_file(rPath));
}

void write_medit_sol(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t dim = medit_sol_dimension(rMesh);
    const std::vector<MeditSolSection> sections = medit_sol_sections_of(rMesh, dim);
    if (!rPath.empty() && rPath.back() == 'b')
        medit_sol_write_binary(rPath, dim, sections);
    else
        medit_sol_write_ascii(rPath, dim, sections);
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
