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
// Elmer mesh directory reader and writer. The layouts follow ElmerGrid's own
// writer and ElmerSolver's reader; see doc/formats/elmer.md.
// Python twin: src/python/meshioplusplus/elmer/_elmer.py.

// System includes
#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/elmer.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

// Elmer type code <-> meshio++ cell type. Codes not listed (102 periodic
// pairs, 412/416 quads, 614 pyramid, the higher-order families) have no
// meshio++ counterpart with the same nodes.
struct ElmType {
    int mCode;
    const char* mType;
};

constexpr ElmType kElmTypes[] = {
    {101, "vertex"},   {202, "line"},       {203, "line3"},        {204, "line4"},
    {303, "triangle"}, {306, "triangle6"},  {310, "triangle10"},   {404, "quad"},
    {408, "quad8"},    {409, "quad9"},      {504, "tetra"},        {510, "tetra10"},
    {605, "pyramid"},  {613, "pyramid13"},  {706, "wedge"},        {715, "wedge15"},
    {718, "wedge18"},  {808, "hexahedron"}, {820, "hexahedron20"}, {827, "hexahedron27"},
};

const char* elm_type_of(int Code) {
    for (const ElmType& t : kElmTypes)
        if (t.mCode == Code)
            return t.mType;
    return nullptr;
}

int elm_code_of(std::string_view Type) {
    for (const ElmType& t : kElmTypes)
        if (Type == t.mType)
            return t.mCode;
    return -1;
}

// Corner count of a cell type's family: the nodes a facet lookup keys on.
std::size_t elm_num_corners(std::string_view Type) {
    static constexpr std::pair<std::string_view, std::size_t> kFamilies[] = {
        {"vertex", 1}, {"line", 2},    {"triangle", 3}, {"quad", 4},
        {"tetra", 4},  {"pyramid", 5}, {"wedge", 6},    {"hexahedron", 8},
    };
    for (const auto& [prefix, n] : kFamilies)
        if (Type.substr(0, prefix.size()) == prefix)
            return n;
    return 0;
}

[[noreturn]] void elm_fail(const fs::path& rFile, std::size_t Line, const std::string& rWhat) {
    throw ReadError("Elmer mesh: " + rFile.string() + ":" + std::to_string(Line) + ": " + rWhat);
}

// Whitespace-separated tokens of one line.
void elm_split(std::string_view Line, std::vector<std::string_view>& rTokens) {
    rTokens.clear();
    std::size_t i = 0;
    while (i < Line.size()) {
        while (i < Line.size() && (Line[i] == ' ' || Line[i] == '\t' || Line[i] == '\r'))
            ++i;
        const std::size_t start = i;
        while (i < Line.size() && Line[i] != ' ' && Line[i] != '\t' && Line[i] != '\r')
            ++i;
        if (i > start)
            rTokens.push_back(Line.substr(start, i - start));
    }
}

std::optional<std::int64_t> elm_int(std::string_view Token) {
    std::int64_t v = 0;
    const char* end = Token.data() + Token.size();
    const auto [ptr, ec] = std::from_chars(Token.data(), end, v);
    if (ec != std::errc() || ptr != end)
        return std::nullopt;
    return v;
}

// An element id, possibly `id/part` (ElmerGrid's halo elements name the part
// that owns them).
struct ElmId {
    std::int64_t mId = 0;
    int mOwner = -1;  // 1-based owning part, or -1 when not given
};

std::optional<ElmId> elm_id(std::string_view Token) {
    const std::size_t slash = Token.find('/');
    ElmId out;
    const auto id = elm_int(Token.substr(0, slash));
    if (!id)
        return std::nullopt;
    out.mId = *id;
    if (slash != std::string_view::npos) {
        const auto owner = elm_int(Token.substr(slash + 1));
        if (!owner)
            return std::nullopt;
        out.mOwner = static_cast<int>(*owner);
    }
    return out;
}

// Calls `rOnLine(tokens, line_number)` for every non-blank line of a file.
template <class F>
void elm_for_each_line(const fs::path& rFile, F&& rOnLine) {
    auto in = detail::make_classic_ifstream(rFile.string());
    if (!in)
        throw ReadError("Elmer mesh: cannot open " + rFile.string());
    std::string line;
    std::vector<std::string_view> tokens;
    std::size_t number = 0;
    while (std::getline(in, line)) {
        ++number;
        elm_split(line, tokens);
        if (!tokens.empty())
            rOnLine(tokens, number);
    }
}

// Every element and boundary element read, in file order, before blocks are
// built. Node ids are the file's.
struct ElmRecords {
    std::vector<std::int64_t> mIds;
    std::vector<std::int64_t> mTags;  // body or boundary id
    std::vector<int> mCodes;
    std::vector<int> mParts;             // 0-based part, or -1 for a serial mesh
    std::vector<std::int64_t> mParents;  // a boundary element's first parent, 0 if none
    std::vector<std::size_t> mOffsets{0};
    std::vector<std::int64_t> mNodes;

    std::size_t Size() const { return mIds.size(); }
};

struct ElmMesh {
    std::vector<std::int64_t> mNodeIds;
    std::vector<double> mCoords;
    std::unordered_map<std::int64_t, std::int64_t> mNodeIndex;
    ElmRecords mBulk;
    ElmRecords mBoundary;
    std::map<std::int64_t, std::string> mBodyNames;
    std::map<std::int64_t, std::string> mBoundaryNames;
    bool mHasParts = false;
};

void elm_read_nodes(const fs::path& rFile, ElmMesh& rMesh) {
    elm_for_each_line(rFile, [&](const std::vector<std::string_view>& rTok, std::size_t Line) {
        if (rTok.size() < 5)
            elm_fail(rFile, Line, "a node line needs `id part x y z`");
        const auto id = elm_int(rTok[0]);
        if (!id)
            elm_fail(rFile, Line, "bad node id '" + std::string(rTok[0]) + "'");
        double xyz[3];
        for (std::size_t d = 0; d < 3; ++d) {
            const std::string text(rTok[2 + d]);
            const char* end = nullptr;
            xyz[d] = detail::parse_double(text.c_str(), end);
            if (end != text.c_str() + text.size())
                elm_fail(rFile, Line, "bad coordinate '" + text + "'");
        }
        // A shared node is listed by every part that uses it: keep the first.
        if (!rMesh.mNodeIndex.emplace(*id, static_cast<std::int64_t>(rMesh.mNodeIds.size())).second)
            return;
        rMesh.mNodeIds.push_back(*id);
        rMesh.mCoords.insert(rMesh.mCoords.end(), xyz, xyz + 3);
    });
}

// `mesh.elements` (`id body type nodes`) or `mesh.boundary` (`id boundary
// parent1 parent2 type nodes`). `Part` is the 0-based part the file belongs
// to, -1 for a serial mesh. Bulk elements are keyed by id, boundary elements by
// (id, nodes), so a merged partitioned mesh holds each once.
void elm_read_elements(const fs::path& rFile, bool Boundary, int Part, bool Lenient,
                       ElmRecords& rOut, std::unordered_map<std::int64_t, std::size_t>& rSeen,
                       std::set<std::pair<std::int64_t, std::vector<std::int64_t>>>& rSeenSides) {
    const std::size_t lead = Boundary ? 5 : 3;
    std::size_t skipped = 0;
    std::vector<std::int64_t> nodes;
    elm_for_each_line(rFile, [&](const std::vector<std::string_view>& rTok, std::size_t Line) {
        if (rTok.size() < lead)
            elm_fail(rFile, Line,
                     Boundary ? "a boundary line needs `id boundary parent1 parent2 type nodes`"
                              : "an element line needs `id body type nodes`");
        const auto id = elm_id(rTok[0]);
        const auto tag = elm_int(rTok[1]);
        const auto code = elm_int(rTok[lead - 1]);
        if (!id || !tag || !code)
            elm_fail(rFile, Line, "malformed line");
        const std::size_t count = static_cast<std::size_t>(*code % 100);
        if (rTok.size() != lead + count)
            elm_fail(
                rFile, Line,
                "type " + std::to_string(*code) + " needs " + std::to_string(count) + " nodes");
        if (!elm_type_of(static_cast<int>(*code))) {
            if (!Lenient)
                elm_fail(rFile, Line,
                         "element type " + std::to_string(*code) + " has no meshio++ cell type");
            ++skipped;
            return;
        }
        nodes.clear();
        for (std::size_t k = 0; k < count; ++k) {
            const auto n = elm_int(rTok[lead + k]);
            if (!n)
                elm_fail(rFile, Line, "bad node id '" + std::string(rTok[lead + k]) + "'");
            nodes.push_back(*n);
        }
        const int part = id->mOwner > 0 ? id->mOwner - 1 : Part;
        if (Part >= 0) {
            if (Boundary) {
                std::vector<std::int64_t> key = nodes;
                std::sort(key.begin(), key.end());
                if (!rSeenSides.emplace(*tag, std::move(key)).second)
                    return;
            } else {
                const auto [it, fresh] = rSeen.emplace(id->mId, rOut.Size());
                if (!fresh) {
                    // A halo copy names its owner; the owner's own copy has no
                    // `/part`. Either way the owner wins.
                    if (id->mOwner < 0)
                        rOut.mParts[it->second] = Part;
                    return;
                }
            }
        }
        rOut.mIds.push_back(id->mId);
        rOut.mTags.push_back(*tag);
        rOut.mCodes.push_back(static_cast<int>(*code));
        rOut.mParts.push_back(part);
        rOut.mParents.push_back(Boundary ? elm_int(rTok[2]).value_or(0) : 0);
        rOut.mNodes.insert(rOut.mNodes.end(), nodes.begin(), nodes.end());
        rOut.mOffsets.push_back(rOut.mNodes.size());
    });
    if (skipped)
        log::warn("Elmer mesh: {} element(s) of types with no meshio++ cell type skipped in {}",
                  skipped, rFile.string());
}

void elm_read_names(const fs::path& rFile, ElmMesh& rMesh) {
    std::error_code ec;
    if (!fs::is_regular_file(rFile, ec))
        return;
    auto in = detail::make_classic_ifstream(rFile.string());
    std::string line;
    bool bodies = true;
    while (std::getline(in, line)) {
        std::string lower = line;
        for (char& c : lower)
            if (c >= 'A' && c <= 'Z')
                c = static_cast<char>(c - 'A' + 'a');
        const std::size_t dollar = line.find('$');
        const std::size_t equals = line.find('=', dollar == std::string::npos ? 0 : dollar);
        if (dollar == std::string::npos || equals == std::string::npos) {
            if (lower.find("names for bound") != std::string::npos)
                bodies = false;
            else if (lower.find("names for bod") != std::string::npos)
                bodies = true;
            continue;
        }
        std::string name = line.substr(dollar + 1, equals - dollar - 1);
        const std::size_t first = name.find_first_not_of(" \t");
        const std::size_t last = name.find_last_not_of(" \t");
        name = first == std::string::npos ? std::string() : name.substr(first, last - first + 1);
        std::vector<std::string_view> tokens;
        elm_split(std::string_view(line).substr(equals + 1), tokens);
        const auto id = tokens.empty() ? std::nullopt : elm_int(tokens[0]);
        if (!id || name.empty())
            continue;
        (bodies ? rMesh.mBodyNames : rMesh.mBoundaryNames)[*id] = name;
    }
}

// The 1-based part directories' files, `part.1.*` ... while they exist.
std::size_t elm_count_parts(const fs::path& rDir) {
    std::error_code ec;
    std::size_t n = 0;
    while (fs::is_regular_file(rDir / ("part." + std::to_string(n + 1) + ".header"), ec))
        ++n;
    return n;
}

std::vector<fs::path> elm_partition_dirs(const fs::path& rDir) {
    std::vector<fs::path> out;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rDir, ec)) {
        const std::string name = entry.path().filename().string();
        if (name.rfind("partitioning.", 0) == 0 && entry.is_directory(ec) &&
            elm_count_parts(entry.path()) > 0)
            out.push_back(entry.path());
    }
    std::sort(out.begin(), out.end());
    return out;
}

void elm_check_text(const fs::path& rDir, const std::string& rStem) {
    std::error_code ec;
    if (fs::is_regular_file(rDir / rStem, ec))
        return;
    if (fs::is_regular_file(rDir / (rStem + ".bin"), ec))
        throw ReadError("Elmer mesh: " + (rDir / rStem).string() +
                        " is only present in ElmerGrid's binary form, which is not supported; "
                        "write the mesh without -bin");
    throw ReadError("Elmer mesh: missing " + (rDir / rStem).string());
}

void elm_read_serial(const fs::path& rDir, bool Lenient, ElmMesh& rMesh) {
    for (const char* stem : {"mesh.nodes", "mesh.elements", "mesh.boundary"})
        elm_check_text(rDir, stem);
    std::unordered_map<std::int64_t, std::size_t> seen;
    std::set<std::pair<std::int64_t, std::vector<std::int64_t>>> seen_sides;
    elm_read_nodes(rDir / "mesh.nodes", rMesh);
    elm_read_elements(rDir / "mesh.elements", false, -1, Lenient, rMesh.mBulk, seen, seen_sides);
    elm_read_elements(rDir / "mesh.boundary", true, -1, Lenient, rMesh.mBoundary, seen, seen_sides);
    elm_read_names(rDir / "mesh.names", rMesh);
}

// Parts `First`..`Last` (0-based, inclusive) of a partitioning directory, merged.
void elm_read_parts(const fs::path& rDir, std::size_t First, std::size_t Last, bool Lenient,
                    ElmMesh& rMesh) {
    std::unordered_map<std::int64_t, std::size_t> seen;
    std::set<std::pair<std::int64_t, std::vector<std::int64_t>>> seen_sides;
    for (std::size_t p = First; p <= Last; ++p) {
        const std::string stem = "part." + std::to_string(p + 1);
        for (const char* ext : {".nodes", ".elements", ".boundary"})
            elm_check_text(rDir, stem + ext);
        elm_read_nodes(rDir / (stem + ".nodes"), rMesh);
        elm_read_elements(rDir / (stem + ".elements"), false, static_cast<int>(p), Lenient,
                          rMesh.mBulk, seen, seen_sides);
        elm_read_elements(rDir / (stem + ".boundary"), true, static_cast<int>(p), Lenient,
                          rMesh.mBoundary, seen, seen_sides);
    }
    rMesh.mHasParts = true;
    // ElmerGrid keeps mesh.names next to the partitioning directory.
    elm_read_names(rDir.parent_path() / "mesh.names", rMesh);
}

// Labels a serial mesh with the parts of its one partitioning directory: bulk
// elements by id, boundary elements by their first parent's part (-1 with none).
void elm_label_serial(const fs::path& rPartDir, ElmMesh& rMesh) {
    const std::size_t parts = elm_count_parts(rPartDir);
    std::unordered_map<std::int64_t, int> part_of;
    for (std::size_t p = 0; p < parts; ++p) {
        const fs::path file = rPartDir / ("part." + std::to_string(p + 1) + ".elements");
        std::error_code ec;
        if (!fs::is_regular_file(file, ec))
            return;
        elm_for_each_line(file, [&](const std::vector<std::string_view>& rTok, std::size_t) {
            const auto id = elm_id(rTok[0]);
            if (!id)
                return;
            const int owner = id->mOwner > 0 ? id->mOwner - 1 : static_cast<int>(p);
            if (id->mOwner < 0 || part_of.find(id->mId) == part_of.end())
                part_of[id->mId] = owner;
        });
    }
    std::vector<int> labels(rMesh.mBulk.Size(), -1);
    for (std::size_t e = 0; e < rMesh.mBulk.Size(); ++e) {
        const auto it = part_of.find(rMesh.mBulk.mIds[e]);
        if (it == part_of.end()) {
            log::warn("Elmer mesh: {} does not match the serial mesh; no partition labels",
                      rPartDir.string());
            return;
        }
        labels[e] = it->second;
    }
    rMesh.mBulk.mParts = labels;
    for (std::size_t e = 0; e < rMesh.mBoundary.Size(); ++e) {
        const auto it = part_of.find(rMesh.mBoundary.mParents[e]);
        rMesh.mBoundary.mParts[e] = it == part_of.end() ? -1 : it->second;
    }
    rMesh.mHasParts = true;
}

Mesh elm_build(const ElmMesh& rIn, const fs::path& rWhere) {
    Mesh mesh;
    NDArray points(DType::Float64, {rIn.mNodeIds.size(), 3});
    std::copy(rIn.mCoords.begin(), rIn.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    std::vector<NDArray> part_blocks;
    std::map<std::int64_t, std::vector<std::int64_t>> body_cells, boundary_cells;
    std::map<std::int64_t, int> body_dim, boundary_dim;
    std::int64_t global = 0;
    for (const bool boundary : {false, true}) {
        const ElmRecords& rec = boundary ? rIn.mBoundary : rIn.mBulk;
        std::vector<int> codes;
        for (int c : rec.mCodes)
            if (std::find(codes.begin(), codes.end(), c) == codes.end())
                codes.push_back(c);
        for (const int code : codes) {
            const std::string type = elm_type_of(code);
            const std::size_t k = static_cast<std::size_t>(code % 100);
            const detail::NodeOrder* order = detail::node_order("elmer", type);
            const int dim = cell_type_dimension(cell_type_from_name(type));
            std::vector<std::size_t> rows;
            for (std::size_t e = 0; e < rec.Size(); ++e)
                if (rec.mCodes[e] == code)
                    rows.push_back(e);
            NDArray conn(DType::Int64, {rows.size(), k});
            NDArray part(DType::Int64, {rows.size()});
            std::int64_t* c = conn.As<std::int64_t>();
            std::int64_t* pp = part.As<std::int64_t>();
            for (std::size_t r = 0; r < rows.size(); ++r) {
                const std::size_t e = rows[r];
                const std::int64_t* src = rec.mNodes.data() + rec.mOffsets[e];
                for (std::size_t j = 0; j < k; ++j) {
                    const std::size_t from =
                        order ? static_cast<std::size_t>(order->mToMeshio[j]) : j;
                    const auto it = rIn.mNodeIndex.find(src[from]);
                    if (it == rIn.mNodeIndex.end())
                        throw ReadError("Elmer mesh: " + rWhere.string() + ": " +
                                        std::string(boundary ? "boundary element " : "element ") +
                                        std::to_string(rec.mIds[e]) + " names undefined node " +
                                        std::to_string(src[from]));
                    c[r * k + j] = it->second;
                }
                pp[r] = rec.mParts[e];
                auto& cells = boundary ? boundary_cells : body_cells;
                auto& dims = boundary ? boundary_dim : body_dim;
                cells[rec.mTags[e]].push_back(global);
                const auto d = dims.emplace(rec.mTags[e], dim);
                d.first->second = std::max(d.first->second, dim);
                ++global;
            }
            mesh.AddCellBlock(type, std::move(conn));
            part_blocks.push_back(std::move(part));
        }
    }
    if (rIn.mHasParts)
        mesh.AddCellData("partition:part", std::move(part_blocks));

    // A name used for both a body and a boundary would make two regions that
    // differ only in their tags; prefix both.
    std::set<std::string> body_names, clash;
    for (const auto& [id, name] : rIn.mBodyNames)
        if (body_cells.count(id))
            body_names.insert(name);
    for (const auto& [id, name] : rIn.mBoundaryNames)
        if (boundary_cells.count(id) && body_names.count(name))
            clash.insert(name);
    for (const bool boundary : {false, true}) {
        const auto& cells = boundary ? boundary_cells : body_cells;
        const auto& names = boundary ? rIn.mBoundaryNames : rIn.mBodyNames;
        const auto& dims = boundary ? boundary_dim : body_dim;
        for (const auto& [id, list] : cells) {
            const auto named = names.find(id);
            std::string name;
            if (named == names.end())
                name = (boundary ? "boundary_" : "body_") + std::to_string(id);
            else if (clash.count(named->second))
                name = (boundary ? "boundary:" : "body:") + named->second;
            else
                name = named->second;
            NDArray entries(DType::Int64, {list.size()});
            std::copy(list.begin(), list.end(), entries.As<std::int64_t>());
            mesh.AddRegion(Region(name, RegionKind::Cell, dims.at(id), id, std::move(entries)));
        }
    }
    return mesh;
}

}  // namespace

Mesh read_elmer(const std::string& rPath, const ReadOptions& rOptions) {
    std::error_code ec;
    fs::path dir(rPath);
    if (dir.filename() == "mesh.header" && fs::is_regular_file(dir, ec))
        dir = dir.has_parent_path() ? dir.parent_path() : fs::path(".");
    if (!fs::is_directory(dir, ec))
        throw ReadError("Elmer mesh: '" + rPath + "' is not a mesh directory");

    ElmMesh mesh;
    const bool serial = fs::is_regular_file(dir / "mesh.header", ec);
    fs::path part_dir;
    if (dir.filename().string().rfind("partitioning.", 0) == 0 && elm_count_parts(dir) > 0) {
        part_dir = dir;
    } else {
        const std::vector<fs::path> parts = elm_partition_dirs(dir);
        if (parts.size() == 1)
            part_dir = parts.front();
        else if (parts.size() > 1 && (!serial || rOptions.mPieceSet))
            throw ReadError("Elmer mesh: " + dir.string() +
                            " holds several partitioning directories; name one of them");
        else if (parts.size() > 1)
            log::warn("Elmer mesh: {} holds several partitioning directories; no partition labels",
                      dir.string());
        if (!serial && part_dir.empty())
            throw ReadError("Elmer mesh: " + dir.string() +
                            " holds neither mesh.header nor a partitioning directory");
    }

    if (rOptions.mPieceSet) {
        if (part_dir.empty())
            throw ReadError("Elmer mesh: " + dir.string() +
                            " is not partitioned; it has no pieces");
        const std::size_t piece = rOptions.ResolvePiece(elm_count_parts(part_dir));
        elm_read_parts(part_dir, piece, piece, rOptions.mLenient, mesh);
        return elm_build(mesh, part_dir);
    }
    if (serial) {
        elm_read_serial(dir, rOptions.mLenient, mesh);
        if (!part_dir.empty())
            elm_label_serial(part_dir, mesh);
        return elm_build(mesh, dir);
    }
    elm_read_parts(part_dir, 0, elm_count_parts(part_dir) - 1, rOptions.mLenient, mesh);
    return elm_build(mesh, part_dir);
}

namespace {

// One group of cells written as a body or a boundary: its id and name (empty
// for an unnamed default).
struct ElmGroup {
    std::int64_t mId;
    std::string mName;
};

// `$` and `=` delimit a mesh.names line, and a line break ends it.
std::string elm_clean_name(const std::string& rName) {
    std::string out = rName;
    for (char& c : out)
        if (c == '$' || c == '=' || c == '\n' || c == '\r' || c == '!')
            c = '_';
    return out;
}

void elm_append_int(std::string& rOut, std::int64_t Value) {
    char buf[24];
    const auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), Value);
    rOut.append(buf, ptr);
}

// "%-6d": left-aligned, padded to six columns (ElmerGrid's header layout).
void elm_append_padded(std::string& rOut, std::int64_t Value) {
    const std::size_t start = rOut.size();
    elm_append_int(rOut, Value);
    while (rOut.size() - start < 6)
        rOut += ' ';
}

void elm_write_file(const fs::path& rFile, const std::string& rText) {
    auto out = detail::make_classic_ofstream(rFile.string(), std::ios::binary);
    if (!out)
        throw WriteError("Elmer mesh: cannot write " + rFile.string());
    out << rText;
}

// Hands out ids: a region's tag when positive and not taken, otherwise the next
// free id above everything handed out so far.
struct ElmIds {
    std::set<std::int64_t> mTaken;

    std::int64_t Take(std::int64_t Tag) {
        if (Tag > 0 && mTaken.insert(Tag).second)
            return Tag;
        std::int64_t id = mTaken.empty() ? 1 : *mTaken.rbegin() + 1;
        mTaken.insert(id);
        return id;
    }
};

}  // namespace

void write_elmer(const std::string& rPath, const Mesh& rMesh) {
    const fs::path dir(rPath);
    const std::size_t n_blocks = rMesh.NumCellBlocks();
    int bulk_dim = -1;
    std::vector<int> codes(n_blocks), dims(n_blocks);
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        codes[b] = cb.IsRagged() ? -1 : elm_code_of(type);
        if (codes[b] < 0)
            throw WriteError("Elmer mesh writer: cell type '" + type +
                             "' has no Elmer element code");
        dims[b] = cell_type_dimension(cell_type_from_name(type));
        if (cb.NumCells() > 0)
            bulk_dim = std::max(bulk_dim, dims[b]);
    }
    if (bulk_dim < 1)
        throw WriteError("Elmer mesh writer: the mesh has no cells of dimension 1 or more");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::int64_t n_cells = detail::total_cells(bases);
    std::vector<bool> is_bulk(static_cast<std::size_t>(n_cells), false);
    for (std::size_t b = 0; b < n_blocks; ++b)
        if (dims[b] == bulk_dim)
            std::fill(is_bulk.begin() + bases[b], is_bulk.begin() + bases[b + 1], true);

    // --- ids: regions first, in the mesh's (kind, name, dim, tag) order ----------
    std::vector<std::int64_t> cell_id(static_cast<std::size_t>(n_cells), 0);
    ElmIds body_ids, boundary_ids;
    std::vector<ElmGroup> bodies, boundaries;
    std::size_t overlaps = 0, point_regions = 0;
    std::vector<std::pair<std::int64_t, std::int64_t>> side_facets;  // (cell, facet)
    std::vector<std::int64_t> side_ids;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const std::int64_t* e = reg.Entries();
        const std::size_t n = reg.NumEntries();
        if (reg.mKind == RegionKind::Point) {
            ++point_regions;
            continue;
        }
        if (reg.mKind == RegionKind::Side) {
            if (n == 0)
                continue;
            const std::int64_t id = boundary_ids.Take(reg.mTag);
            boundaries.push_back({id, reg.mName});
            for (std::size_t j = 0; j < n; ++j) {
                side_facets.emplace_back(e[2 * j], e[2 * j + 1]);
                side_ids.push_back(id);
            }
            continue;
        }
        // A cell region may hold bulk cells, lower-dimensional cells, or both.
        std::vector<std::int64_t> bulk, lower;
        for (std::size_t j = 0; j < n; ++j)
            if (e[j] >= 0 && e[j] < n_cells)
                (is_bulk[static_cast<std::size_t>(e[j])] ? bulk : lower).push_back(e[j]);
        for (const bool boundary : {false, true}) {
            const auto& list = boundary ? lower : bulk;
            if (list.empty())
                continue;
            std::vector<std::int64_t> free_cells;
            for (std::int64_t g : list) {
                if (cell_id[static_cast<std::size_t>(g)] != 0)
                    ++overlaps;
                else
                    free_cells.push_back(g);
            }
            if (free_cells.empty())
                continue;
            std::string name = reg.mName;
            if (!bulk.empty() && !lower.empty() && boundary)
                name += "_boundary";
            const std::int64_t id = (boundary ? boundary_ids : body_ids).Take(reg.mTag);
            (boundary ? boundaries : bodies).push_back({id, name});
            for (std::int64_t g : free_cells)
                cell_id[static_cast<std::size_t>(g)] = id;
        }
    }
    // Cells in no region: one fresh id per block.
    for (std::size_t b = 0; b < n_blocks; ++b) {
        std::int64_t id = 0;
        for (std::int64_t g = bases[b]; g < bases[b + 1]; ++g) {
            if (cell_id[static_cast<std::size_t>(g)] != 0)
                continue;
            if (id == 0) {
                id = (dims[b] == bulk_dim ? body_ids : boundary_ids).Take(-1);
                (dims[b] == bulk_dim ? bodies : boundaries).push_back({id, ""});
            }
            cell_id[static_cast<std::size_t>(g)] = id;
        }
    }

    // --- notes, rendered into mesh.names' provenance block ----------------------
    if (overlaps)
        log::warn(
            "Elmer mesh writer: {} cell(s) belong to more than one region; the first "
            "region (in name order) wins",
            overlaps);
    if (point_regions) {
        log::warn("Elmer mesh writer: {} point region(s) have no Elmer equivalent and were dropped",
                  point_regions);
        detail::provenance_note("regions-dropped", std::to_string(point_regions) +
                                                       " point region(s) have no Elmer equivalent");
    }
    if (rMesh.NumPointData() + rMesh.NumCellData() + rMesh.NumFieldData() > 0) {
        log::warn(
            "Elmer mesh writer: an Elmer mesh holds no data arrays; point, cell and field "
            "data dropped");
        detail::provenance_note("data-dropped", "an Elmer mesh directory holds no data arrays");
    }

    // --- parents -----------------------------------------------------------------
    // Bulk element numbers, 1-based, in block order.
    std::vector<std::int64_t> element_no(static_cast<std::size_t>(n_cells), 0);
    std::int64_t n_bulk = 0;
    for (std::int64_t g = 0; g < n_cells; ++g)
        if (is_bulk[static_cast<std::size_t>(g)])
            element_no[static_cast<std::size_t>(g)] = ++n_bulk;

    detail::FacetIndexOptions facet_options;
    facet_options.mSolidFaces = bulk_dim == 3;
    facet_options.mSurfaceEdges = bulk_dim == 2;
    std::optional<detail::FacetIndex> facets;
    std::vector<std::vector<std::int64_t>> node_cells;  // node -> bulk cells, built on demand
    std::size_t orphans = 0;
    const auto bulk_parent = [&](std::int64_t g) -> std::int64_t {
        return g >= 0 && is_bulk[static_cast<std::size_t>(g)]
                   ? element_no[static_cast<std::size_t>(g)]
                   : 0;
    };
    // Parents of the boundary element with these corner nodes and dimension.
    const auto parents_of = [&](const std::vector<std::int64_t>& rCorners,
                                int Dim) -> std::pair<std::int64_t, std::int64_t> {
        if (Dim == bulk_dim - 1 && Dim >= 1) {
            if (!facets)
                facets.emplace(rMesh, facet_options);
            const detail::FacetHit* hit = facets->Find(rCorners.data(), rCorners.size());
            if (!hit)
                return {0, 0};
            return {bulk_parent(hit->mFirst.mCell),
                    hit->mCount > 1 ? bulk_parent(hit->mSecond.mCell) : 0};
        }
        // An edge of a solid, or a point: the first bulk cell holding every corner.
        if (node_cells.empty()) {
            node_cells.resize(rMesh.NumPoints());
            for (std::size_t b = 0; b < n_blocks; ++b) {
                if (dims[b] != bulk_dim)
                    continue;
                const auto cb = rMesh.Cells(b);
                const NDArray& conn = cb.Conn();
                const std::size_t k = cb.NodesPerCell();
                for (std::size_t r = 0; r < cb.NumCells(); ++r)
                    for (std::size_t j = 0; j < k; ++j) {
                        const auto p = static_cast<std::size_t>(detail::read_int(conn, r * k + j));
                        auto& list = node_cells[p];
                        const std::int64_t g = bases[b] + static_cast<std::int64_t>(r);
                        if (list.empty() || list.back() != g)
                            list.push_back(g);
                    }
            }
        }
        for (std::int64_t g : node_cells[static_cast<std::size_t>(rCorners[0])]) {
            bool all = true;
            for (std::size_t c = 1; c < rCorners.size() && all; ++c) {
                const auto& list = node_cells[static_cast<std::size_t>(rCorners[c])];
                all = std::find(list.begin(), list.end(), g) != list.end();
            }
            if (all)
                return {element_no[static_cast<std::size_t>(g)], 0};
        }
        return {0, 0};
    };

    // --- files ---------------------------------------------------------------------
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (!fs::is_directory(dir, ec))
        throw WriteError("Elmer mesh writer: cannot create directory " + dir.string());
    for (const auto& entry : fs::directory_iterator(dir, ec))
        if (entry.path().filename().string().rfind("partitioning.", 0) == 0)
            log::warn(
                "Elmer mesh writer: {} is left over from an earlier mesh and no longer "
                "matches it",
                entry.path().string());

    std::map<int, std::int64_t> type_counts;
    std::string out;
    char buf[32];
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t npts = rMesh.NumPoints();
    for (std::size_t p = 0; p < npts; ++p) {
        elm_append_int(out, static_cast<std::int64_t>(p + 1));
        out += " -1";
        for (std::size_t d = 0; d < 3; ++d) {
            const double v = d < pdim ? detail::read_double(points, p * pdim + d) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), " %.17g", v);
            out += buf;
        }
        out += '\n';
    }
    elm_write_file(dir / "mesh.nodes", out);

    // Bulk elements, then boundary elements.
    std::string elements, boundary;
    std::int64_t n_boundary = 0;
    std::vector<std::int64_t> corners, row;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const std::string type(cb.Type());
        const detail::NodeOrder* order = detail::node_order("elmer", type);
        const std::size_t n_corners = elm_num_corners(type);
        const bool bulk = dims[b] == bulk_dim;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::int64_t g = bases[b] + static_cast<std::int64_t>(r);
            std::string& text = bulk ? elements : boundary;
            if (bulk) {
                elm_append_int(text, element_no[static_cast<std::size_t>(g)]);
                text += ' ';
                elm_append_int(text, cell_id[static_cast<std::size_t>(g)]);
            } else {
                corners.clear();
                for (std::size_t c = 0; c < n_corners; ++c)
                    corners.push_back(detail::read_int(conn, r * k + c));
                const auto [p1, p2] = parents_of(corners, dims[b]);
                orphans += p1 == 0 ? 1 : 0;
                elm_append_int(text, ++n_boundary);
                text += ' ';
                elm_append_int(text, cell_id[static_cast<std::size_t>(g)]);
                text += ' ';
                elm_append_int(text, p1);
                text += ' ';
                elm_append_int(text, p2);
            }
            text += ' ';
            elm_append_int(text, codes[b]);
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                text += ' ';
                elm_append_int(text, detail::read_int(conn, r * k + src) + 1);
            }
            text += '\n';
            ++type_counts[codes[b]];
        }
    }
    // Side regions: each facet one more boundary element.
    std::size_t bad_facets = 0;
    CellType facet_type = CellType::Custom;
    for (std::size_t s = 0; s < side_facets.size(); ++s) {
        if (!detail::facet_nodes(rMesh, side_facets[s].first, side_facets[s].second, facet_type,
                                 row)) {
            ++bad_facets;
            continue;
        }
        const std::string type = cell_type_name(facet_type);
        const int code = elm_code_of(type);
        const detail::NodeOrder* order = detail::node_order("elmer", type);
        const std::size_t n_corners = elm_num_corners(type);
        corners.assign(row.begin(), row.begin() + static_cast<std::ptrdiff_t>(n_corners));
        const auto [p1, p2] = parents_of(corners, cell_type_dimension(facet_type));
        orphans += p1 == 0 ? 1 : 0;
        elm_append_int(boundary, ++n_boundary);
        boundary += ' ';
        elm_append_int(boundary, side_ids[s]);
        boundary += ' ';
        elm_append_int(boundary, p1);
        boundary += ' ';
        elm_append_int(boundary, p2);
        boundary += ' ';
        elm_append_int(boundary, code);
        for (std::size_t j = 0; j < row.size(); ++j) {
            const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
            boundary += ' ';
            elm_append_int(boundary, row[src] + 1);
        }
        boundary += '\n';
        ++type_counts[code];
    }
    if (bad_facets)
        log::warn("Elmer mesh writer: {} side region entr(ies) name no facet and were dropped",
                  bad_facets);
    if (orphans)
        log::warn(
            "Elmer mesh writer: {} boundary element(s) lie on no bulk element; written "
            "with parent 0",
            orphans);
    elm_write_file(dir / "mesh.elements", elements);
    elm_write_file(dir / "mesh.boundary", boundary);

    out.clear();
    elm_append_padded(out, static_cast<std::int64_t>(npts));
    out += ' ';
    elm_append_padded(out, n_bulk);
    out += ' ';
    elm_append_padded(out, n_boundary);
    out += '\n';
    elm_append_padded(out, static_cast<std::int64_t>(type_counts.size()));
    out += '\n';
    for (const auto& [code, count] : type_counts) {
        elm_append_padded(out, code);
        out += ' ';
        elm_append_padded(out, count);
        out += '\n';
    }
    elm_write_file(dir / "mesh.header", out);

    // mesh.names: the provenance block, then the named bodies and boundaries.
    out = detail::provenance_render_lines(detail::SlotTier::Block, "! ");
    for (char& c : out)
        if (c == '$')
            c = '_';
    std::set<std::string> body_names;
    for (const ElmGroup& g : bodies)
        if (!g.mName.empty())
            body_names.insert(elm_clean_name(g.mName));
    for (const bool is_boundary : {false, true}) {
        std::vector<ElmGroup> list = is_boundary ? boundaries : bodies;
        std::sort(list.begin(), list.end(),
                  [](const ElmGroup& a, const ElmGroup& b) { return a.mId < b.mId; });
        out += is_boundary ? "! ----- names for boundaries -----\n"
                           : "! ----- names for bodies -----\n";
        for (const ElmGroup& g : list) {
            if (g.mName.empty())
                continue;
            std::string name = elm_clean_name(g.mName);
            // ElmerSolver matches a name against bodies and boundaries alike.
            if (is_boundary && body_names.count(name))
                name += "_boundary";
            out += "$ " + name + " = ";
            elm_append_int(out, g.mId);
            out += '\n';
        }
    }
    elm_write_file(dir / "mesh.names", out);
}

}  // namespace meshioplusplus
