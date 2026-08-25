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
// The GiD postprocess READER. See gid.hpp for the format's documentation.
//
// Deliberately a separate translation unit from gid.cpp (the vtu.cpp /
// vtu_read.cpp split precedent), and deliberately OUTSIDE that file's
// MESHIOPLUSPLUS_HAS_GIDPOST guard: gidpost is a write-only library, so
// reading needs none of it. The reader's real dependencies are per flavour --
// ASCII needs nothing, the gzip-wrapped flavours need zlib, HDF5 needs HDF5 --
// which is what makes `gid` readable in strictly more build configurations
// than it is writable.

// System includes
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"

#include "gid_common.hpp"

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/detail/hdf5_util.hpp"
#endif

namespace meshioplusplus {

namespace {

using gid_detail::gid_ascii_paths;
using gid_detail::gid_has_suffix;
using gid_detail::gid_resolve_mode;

// ---------------------------------------------------------------------------
// Cell-type mapping: (GiD ElemType name, Nnode) -> meshio++ type.
//
// The inverse of gid.cpp's own write table, and deliberately no larger. GiD's
// legal Nnode values per type (from gidpost's own ValidateConnectivity) are
// wider than what meshio++ can map: Hexahedra also admits 27, Prism 15 and
// Pyramid 13. Those three are exactly the orderings the writer refuses because
// they were never independently verified against GiD's own geometry -- so the
// reader refuses them too, by name. Reading a permutation we could not write
// would be the same unverified guess in the other direction.
//
// Sphere/Circle are refused for a different reason: they have no meshio++
// counterpart, AND their element rows are not all-integer (they embed a radius
// and, for Circle, a normal vector), so they need a different row parser
// entirely rather than merely a type name.
struct GidReadType {
    const char* mGidName;
    int mNumNodes;
    const char* mMeshioName;
};

const std::vector<GidReadType>& gid_read_type_table() {
    static const std::vector<GidReadType> table = {
        {"Point", 1, "vertex"},          {"Linear", 2, "line"},
        {"Linear", 3, "line3"},          {"Triangle", 3, "triangle"},
        {"Triangle", 6, "triangle6"},    {"Quadrilateral", 4, "quad"},
        {"Quadrilateral", 8, "quad8"},   {"Quadrilateral", 9, "quad9"},
        {"Tetrahedra", 4, "tetra"},      {"Tetrahedra", 10, "tetra10"},
        {"Hexahedra", 8, "hexahedron"},  {"Hexahedra", 20, "hexahedron20"},
        {"Prism", 6, "wedge"},           {"Pyramid", 5, "pyramid"},
    };
    return table;
}

std::string gid_meshio_type(const std::string& rGidName, int nnode) {
    for (const GidReadType& e : gid_read_type_table())
        if (rGidName == e.mGidName && nnode == e.mNumNodes)
            return e.mMeshioName;
    throw ReadError("GiD: element type '" + rGidName + "' with Nnode " + std::to_string(nnode) +
                    " has no verified meshio++ mapping");
}

// ---------------------------------------------------------------------------
// Tokenization. std::strtod / std::strtoll throughout, never std::from_chars --
// its floating-point overload is a real Emscripten/libc++ hazard (the same rule
// dex.cpp / wkt.cpp / nastran.cpp follow).

double gid_to_double(const std::string& rTok, const char* pWhat) {
    const char* start = rTok.c_str();
    char* end = nullptr;
    const double v = std::strtod(start, &end);
    if (end == start)
        throw ReadError(std::string("GiD: expected a number for ") + pWhat + ", got '" + rTok +
                        "'");
    return v;
}

std::int64_t gid_to_int(const std::string& rTok, const char* pWhat) {
    const char* start = rTok.c_str();
    char* end = nullptr;
    const long long v = std::strtoll(start, &end, 10);
    if (end == start)
        throw ReadError(std::string("GiD: expected an integer for ") + pWhat + ", got '" + rTok +
                        "'");
    return static_cast<std::int64_t>(v);
}

/// Splits on whitespace, keeping `"quoted strings"` (which gidpost guarantees
/// contain no embedded quote: change_quotes() rewrites any `"` inside a user
/// string to `'` BEFORE embedding it, so "up to the next quote" is exact, not
/// a heuristic). The quotes themselves are stripped from the returned token.
std::vector<std::string> gid_split(const std::string& rLine) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < rLine.size()) {
        while (i < rLine.size() && std::isspace(static_cast<unsigned char>(rLine[i])))
            ++i;
        if (i >= rLine.size())
            break;
        if (rLine[i] == '"') {
            const std::size_t close = rLine.find('"', i + 1);
            if (close == std::string::npos)
                throw ReadError("GiD: unterminated quoted string: " + rLine);
            out.push_back(rLine.substr(i + 1, close - i - 1));
            i = close + 1;
        } else {
            const std::size_t start = i;
            while (i < rLine.size() && !std::isspace(static_cast<unsigned char>(rLine[i])))
                ++i;
            out.push_back(rLine.substr(start, i - start));
        }
    }
    return out;
}

/// True when the line's first non-space character is `#`.
///
/// gidpost emits `# <Name>: <Value>` user attributes (where meshio++'s own
/// provenance block goes) and `# color <r> <g> <b>`, at file scope, after a
/// MESH header and after a Result header. There is no lexer-level restriction
/// on where they appear, so the reader simply drops them everywhere. The
/// provenance block is recovered separately, from the file's bytes, by
/// detail::read_provenance_lines -- not from here.
bool gid_is_comment(const std::string& rLine) {
    for (char ch : rLine) {
        if (std::isspace(static_cast<unsigned char>(ch)))
            continue;
        return ch == '#';
    }
    return false;  // blank
}

bool gid_is_blank(const std::string& rLine) {
    for (char ch : rLine)
        if (!std::isspace(static_cast<unsigned char>(ch)))
            return false;
    return true;
}

// ---------------------------------------------------------------------------
// The line cursor. ASCII and gzip-inflated-ASCII share it verbatim -- the
// gzipped flavour (gidpost's GiD_PostAsciiZipped) is the SAME text, so
// inflating up front is the whole of its support.
class GidLineCursor {
public:
    explicit GidLineCursor(std::string_view text) : mText(text) {}

    /// Next line that is neither blank nor a comment; empty when exhausted.
    bool Next(std::string& rOut) {
        while (mPos < mText.size()) {
            std::size_t nl = mText.find('\n', mPos);
            if (nl == std::string_view::npos)
                nl = mText.size();
            std::string_view raw = mText.substr(mPos, nl - mPos);
            mPos = nl + 1;
            if (!raw.empty() && raw.back() == '\r')
                raw.remove_suffix(1);
            // A Values row may legitimately begin with whitespace -- see the
            // id-suppression rule in gid_read_values -- so leading space must
            // NOT disqualify a line here; only blank and comment lines are
            // dropped.
            std::string line(raw);
            if (gid_is_blank(line) || gid_is_comment(line))
                continue;
            mRawStartedWithSpace =
                !raw.empty() && std::isspace(static_cast<unsigned char>(raw.front()));
            rOut = std::move(line);
            return true;
        }
        return false;
    }

    /// Whether the line just returned by Next() began with whitespace.
    bool LastStartedWithSpace() const { return mRawStartedWithSpace; }

private:
    std::string_view mText;
    std::size_t mPos = 0;
    bool mRawStartedWithSpace = false;
};

// ---------------------------------------------------------------------------
// Staged mesh, assembled across every MESH block before anything is committed.

struct GidBlock {
    std::string mMeshName;
    std::string mMeshioType;
    int mNumNodes = 0;
    std::vector<std::int64_t> mConn;      // global node ids, row-major
    std::vector<std::int64_t> mElemIds;   // as written; may restart per block
    std::vector<std::int64_t> mMaterial;  // empty when the file carried none
};

struct GidStaged {
    // One GLOBAL node table keyed by id, de-duplicated. Real files repeat the
    // full table in every MESH block (Kratos does); meshio++'s own writer
    // writes it once and emits empty Coordinates pairs thereafter. Both are
    // legal and both land here identically. std::map, not unordered_map:
    // iteration order defines the emitted point order, so it must be stable.
    std::map<std::int64_t, std::array<double, 3>> mNodes;
    std::vector<GidBlock> mBlocks;
    bool mAnyThirdCoord = false;  // false only when every row carried 2 coords
};

/// One result, still in file terms (ids not yet resolved to rows).
struct GidResult {
    std::string mName;
    std::string mAnalysis;
    double mStep = 0.0;
    std::string mType;      // Scalar / Vector / Matrix / ...
    std::string mLocation;  // OnNodes / OnGaussPoints / OnNurbs*
    std::string mGaussName;
    std::size_t mNumComponents = 0;
    std::vector<std::int64_t> mIds;
    std::vector<double> mValues;  // mIds.size() * mNumComponents
};

/// A `GaussPoints` declaration: what it is called, and which mesh it is on.
struct GidGaussSet {
    std::string mMeshName;
    int mNumPoints = 1;
};

// ---------------------------------------------------------------------------
// Geometry parsing (.post.msh).

void gid_read_coordinates(GidLineCursor& rCur, GidStaged& rStaged) {
    std::string line;
    while (rCur.Next(line)) {
        if (line.rfind("End Coordinates", 0) == 0)
            return;
        const std::vector<std::string> tok = gid_split(line);
        // gidpost's block writer always emits 3 coordinates, but the per-node
        // GiD_WriteCoordinates2D ASCII path emits only 2 -- so count, never
        // assume 4 tokens.
        if (tok.size() < 3 || tok.size() > 4)
            throw ReadError("GiD: malformed coordinate row: " + line);
        const std::int64_t id = gid_to_int(tok[0], "a node id");
        std::array<double, 3> xyz{0.0, 0.0, 0.0};
        for (std::size_t c = 0; c + 1 < tok.size(); ++c)
            xyz[c] = gid_to_double(tok[c + 1], "a node coordinate");
        if (tok.size() == 4)
            rStaged.mAnyThirdCoord = true;
        // Repeats are expected (the Kratos case) -- keep the first, which
        // makes a repeated-but-identical table a no-op.
        rStaged.mNodes.emplace(id, xyz);
    }
    throw ReadError("GiD: unterminated Coordinates block");
}

void gid_read_elements(GidLineCursor& rCur, GidBlock& rBlock) {
    const std::size_t nn = static_cast<std::size_t>(rBlock.mNumNodes);
    std::string line;
    bool material_seen = false;
    while (rCur.Next(line)) {
        if (line.rfind("End Elements", 0) == 0) {
            if (!material_seen)
                rBlock.mMaterial.clear();
            return;
        }
        const std::vector<std::string> tok = gid_split(line);
        // THE material-column ambiguity. There is no separator between the
        // connectivity and an optional trailing material id, so Nnode -- from
        // this block's own MESH header -- is the only disambiguator:
        //   1 + Nnode      tokens -> id + connectivity
        //   1 + Nnode + 1  tokens -> id + connectivity + material
        // Anything else is malformed, and saying so beats guessing.
        const bool has_mat = tok.size() == nn + 2;
        if (tok.size() != nn + 1 && !has_mat)
            throw ReadError("GiD: element row has " + std::to_string(tok.size()) +
                            " values; expected " + std::to_string(nn + 1) + " or " +
                            std::to_string(nn + 2) + " for Nnode " + std::to_string(nn) + ": " +
                            line);
        rBlock.mElemIds.push_back(gid_to_int(tok[0], "an element id"));
        for (std::size_t k = 0; k < nn; ++k)
            rBlock.mConn.push_back(gid_to_int(tok[k + 1], "an element node"));
        if (has_mat) {
            material_seen = true;
            rBlock.mMaterial.push_back(gid_to_int(tok[nn + 1], "a material id"));
        } else {
            rBlock.mMaterial.push_back(0);
        }
    }
    throw ReadError("GiD: unterminated Elements block");
}

/// Skips to a matching `End <what>` line, for constructs meshio++ understands
/// structurally but does not map (`Group`, `ResultRangesTable`, `OnGroup`).
void gid_skip_to_end(GidLineCursor& rCur, const char* pWhat) {
    const std::string terminator = std::string("End ") + pWhat;
    std::string line;
    while (rCur.Next(line))
        if (line.rfind(terminator, 0) == 0)
            return;
    throw ReadError(std::string("GiD: unterminated ") + pWhat + " block");
}

void gid_parse_mesh_text(std::string_view text, GidStaged& rStaged) {
    GidLineCursor cur(text);
    std::string line;
    while (cur.Next(line)) {
        const std::vector<std::string> tok = gid_split(line);
        if (tok.empty())
            continue;

        if (tok[0] == "MESH") {
            // MESH "<name>" dimension <2|3> ElemType <T> Nnode <K>
            GidBlock block;
            std::string gid_type;
            int nnode = 0;
            block.mMeshName = tok.size() > 1 ? tok[1] : std::string();
            for (std::size_t i = 2; i + 1 < tok.size(); ++i) {
                if (tok[i] == "ElemType")
                    gid_type = tok[i + 1];
                else if (tok[i] == "Nnode")
                    nnode = static_cast<int>(gid_to_int(tok[i + 1], "Nnode"));
            }
            if (gid_type.empty() || nnode <= 0)
                throw ReadError("GiD: malformed MESH header: " + line);
            block.mNumNodes = nnode;
            block.mMeshioType = gid_meshio_type(gid_type, nnode);
            rStaged.mBlocks.push_back(std::move(block));
            continue;
        }

        if (tok[0] == "Coordinates") {
            gid_read_coordinates(cur, rStaged);
            continue;
        }
        if (tok[0] == "Elements") {
            if (rStaged.mBlocks.empty())
                throw ReadError("GiD: Elements block before any MESH header");
            gid_read_elements(cur, rStaged.mBlocks.back());
            continue;
        }
        if (tok[0] == "Group") {
            // Mesh groups have no meshio++ counterpart yet (a documented
            // roadmap item); the MESH blocks inside are read normally, so only
            // the wrapper is ignored.
            continue;
        }
        if (tok[0] == "End" || tok[0] == "Unit")
            continue;  // "End Group"/"End Mesh"-style closers, and mesh units
        // Anything else at mesh scope is a construct we do not map; ignoring it
        // is safe because every block we DO map is self-delimiting.
    }
}

// ---------------------------------------------------------------------------
// Results parsing (.post.res).

/// Legal component counts per GiD result type, from gidpost's own table.
/// Used only to sanity-check the count inferred from the row width, because
/// the `Result` header carries no `:N` dimension suffix -- `Vector` alone does
/// not say whether rows hold 2, 3 or 4 components.
bool gid_dim_is_legal(const std::string& rType, std::size_t k) {
    if (rType == "Scalar")
        return k == 1;
    if (rType == "Vector")
        return k == 2 || k == 3 || k == 4;
    if (rType == "Matrix")
        return k == 3 || k == 6;
    if (rType == "PlainDeformationMatrix")
        return k == 4;
    if (rType == "MainMatrix")
        return k == 12;
    if (rType == "LocalAxes")
        return k == 3;
    if (rType == "ComplexScalar")
        return k == 2;
    if (rType == "ComplexVector")
        return k == 4 || k == 6;
    if (rType == "ComplexMatrix")
        return k == 6 || k == 12;
    return true;  // unknown type: accept whatever the rows say
}

/**
 * @brief Reads a `Values` ... `End Values` block.
 *
 * Handles gidpost's **id-suppression rule**: inside a Values block the row id
 * is omitted whenever it equals the previous row's, so a result with G>1 Gauss
 * points per element writes the element id on the first row only and the next
 * G-1 rows begin with whitespace. A row that carries fewer tokens than
 * `1 + ncomp` -- or that begins with whitespace once the width is known -- is
 * therefore a continuation of the preceding id, not a new one.
 *
 * The component count itself is inferred from the FIRST row's width, since the
 * header cannot supply it.
 */
void gid_read_values(GidLineCursor& rCur, GidResult& rResult) {
    std::string line;
    bool first = true;
    std::int64_t last_id = 0;
    while (rCur.Next(line)) {
        if (line.rfind("End Values", 0) == 0)
            return;
        const std::vector<std::string> tok = gid_split(line);
        if (tok.empty())
            continue;

        if (first) {
            if (tok.size() < 2)
                throw ReadError("GiD: malformed first Values row: " + line);
            rResult.mNumComponents = tok.size() - 1;
            if (!gid_dim_is_legal(rResult.mType, rResult.mNumComponents))
                throw ReadError("GiD: result '" + rResult.mName + "' declares type " +
                                rResult.mType + " but its rows carry " +
                                std::to_string(rResult.mNumComponents) +
                                " components, which that type does not admit");
            first = false;
        }

        const std::size_t k = rResult.mNumComponents;
        const bool continuation = tok.size() == k && rCur.LastStartedWithSpace();
        if (!continuation && tok.size() != k + 1)
            throw ReadError("GiD: Values row has " + std::to_string(tok.size()) +
                            " tokens; expected " + std::to_string(k + 1) + " (or " +
                            std::to_string(k) + " for a suppressed repeated id): " + line);

        const std::size_t base = continuation ? 0 : 1;
        if (!continuation)
            last_id = gid_to_int(tok[0], "a result id");
        rResult.mIds.push_back(last_id);
        for (std::size_t c = 0; c < k; ++c)
            rResult.mValues.push_back(gid_to_double(tok[base + c], "a result value"));
    }
    throw ReadError("GiD: unterminated Values block");
}

void gid_parse_res_text(std::string_view text, std::vector<GidResult>& rResults,
                        std::unordered_map<std::string, GidGaussSet>& rGauss) {
    GidLineCursor cur(text);
    std::string line;
    while (cur.Next(line)) {
        const std::vector<std::string> tok = gid_split(line);
        if (tok.empty())
            continue;

        if (tok[0] == "GaussPoints") {
            // GaussPoints "<gp>" ElemType <T> ["<meshname>"]
            // There is NO OnMeshName keyword -- the mesh name is a bare
            // trailing quoted string, which is why it is found positionally.
            GidGaussSet set;
            const std::string gp_name = tok.size() > 1 ? tok[1] : std::string();
            if (tok.size() >= 5)
                set.mMeshName = tok[4];
            std::string inner;
            while (cur.Next(inner)) {
                if (inner.rfind("End GaussPoints", 0) == 0)
                    break;
                const std::vector<std::string> it = gid_split(inner);
                if (it.size() >= 5 && it[0] == "Number" && it[1] == "Of")
                    set.mNumPoints = static_cast<int>(gid_to_int(it[4], "a Gauss point count"));
            }
            rGauss[gp_name] = set;
            continue;
        }

        if (tok[0] == "Result") {
            // Result "<name>" "<analysis>" <step> <Type> <Location> ["<gp>"]
            if (tok.size() < 6)
                throw ReadError("GiD: malformed Result header: " + line);
            GidResult res;
            res.mName = tok[1];
            res.mAnalysis = tok[2];
            res.mStep = gid_to_double(tok[3], "a result step");
            res.mType = tok[4];
            res.mLocation = tok[5];
            if (tok.size() >= 7)
                res.mGaussName = tok[6];
            // Optional ResultRangesTable / ComponentNames / Unit lines may sit
            // between the header and Values, in that fixed order.
            std::string inner;
            while (cur.Next(inner)) {
                const std::vector<std::string> it = gid_split(inner);
                if (!it.empty() && it[0] == "Values") {
                    gid_read_values(cur, res);
                    break;
                }
                if (!it.empty() &&
                    (it[0] == "ResultRangesTable" || it[0] == "ComponentNames" || it[0] == "Unit"))
                    continue;
                throw ReadError("GiD: unexpected line in result '" + res.mName + "': " + inner);
            }
            rResults.push_back(std::move(res));
            continue;
        }

        if (tok[0] == "ResultRangesTable") {
            gid_skip_to_end(cur, "ResultRangesTable");
            continue;
        }
        if (tok[0] == "ResultGroup") {
            // A ResultGroup packs several results into one wide Values row.
            // meshio++ has never written one and unpacking it needs the
            // per-member ResultDescription dims; refuse by name rather than
            // mis-associate columns.
            throw ReadError(
                "GiD: ResultGroup blocks are not supported by this reader (a documented gap)");
        }
        if (tok[0] == "OnGroup") {
            gid_skip_to_end(cur, "OnGroup");
            continue;
        }
        // Everything else at result scope (the file magic, End markers, ...)
        // is skipped.
    }
}

// ---------------------------------------------------------------------------
// Assembly: staged file structures -> Mesh.

void gid_apply_results(Mesh& rMesh, const GidStaged& rStaged,
                       const std::map<std::int64_t, std::size_t>& rNodeRow,
                       const std::vector<GidResult>& rResults,
                       const std::unordered_map<std::string, GidGaussSet>& rGauss,
                       int time_step) {
    // Distinct step values, in first-seen order, so mTimeStep can select one.
    std::vector<double> steps;
    for (const GidResult& r : rResults) {
        bool seen = false;
        for (double s : steps)
            seen = seen || s == r.mStep;
        if (!seen)
            steps.push_back(r.mStep);
    }
    double wanted = steps.empty() ? 0.0 : steps.front();
    if (!steps.empty()) {
        std::int64_t idx = time_step;
        if (idx < 0)
            idx += static_cast<std::int64_t>(steps.size());
        if (idx < 0 || idx >= static_cast<std::int64_t>(steps.size()))
            throw ReadError("GiD: time step " + std::to_string(time_step) + " out of range (" +
                            std::to_string(steps.size()) + " available)");
        wanted = steps[static_cast<std::size_t>(idx)];
    }

    // Per-block element-id -> row. Element ids may RESTART per block in real
    // files (Kratos does exactly that), so this must never be one global map.
    // Same-named results are accumulated here and committed once, in
    // first-seen order, so a cell_data array split across blocks survives.
    std::map<std::string, std::vector<NDArray>> pending;
    std::vector<std::string> order;

    std::vector<std::map<std::int64_t, std::size_t>> elem_rows(rStaged.mBlocks.size());
    for (std::size_t b = 0; b < rStaged.mBlocks.size(); ++b)
        for (std::size_t r = 0; r < rStaged.mBlocks[b].mElemIds.size(); ++r)
            elem_rows[b].emplace(rStaged.mBlocks[b].mElemIds[r], r);

    for (const GidResult& res : rResults) {
        if (res.mStep != wanted)
            continue;
        if (res.mNumComponents == 0)
            continue;

        if (res.mLocation == "OnNodes") {
            NDArray arr(DType::Float64, res.mNumComponents == 1
                                            ? std::vector<std::size_t>{rNodeRow.size()}
                                            : std::vector<std::size_t>{rNodeRow.size(),
                                                                       res.mNumComponents});
            double* dst = arr.As<double>();
            for (std::size_t i = 0; i < res.mIds.size(); ++i) {
                auto it = rNodeRow.find(res.mIds[i]);
                if (it == rNodeRow.end())
                    continue;  // a value for a node the geometry never defined
                for (std::size_t c = 0; c < res.mNumComponents; ++c)
                    dst[it->second * res.mNumComponents + c] =
                        res.mValues[i * res.mNumComponents + c];
            }
            rMesh.AddPointData(res.mName, std::move(arr));
            continue;
        }

        if (res.mLocation != "OnGaussPoints") {
            log::warn("gid: result '{}' has location '{}', which has no meshio++ counterpart "
                      "-- skipped",
                      res.mName, res.mLocation);
            continue;
        }

        auto gp = rGauss.find(res.mGaussName);
        if (gp == rGauss.end()) {
            log::warn("gid: result '{}' names Gauss-point set '{}', which the file never "
                      "declares -- skipped",
                      res.mName, res.mGaussName);
            continue;
        }
        if (gp->second.mNumPoints != 1) {
            // meshio++'s cell_data is (n,)/(n,k), never per-node-within-cell --
            // the same structural limit MED's ELNO/ELGA already documents.
            // Averaging or taking the first point would invent data.
            log::warn("gid: result '{}' has {} Gauss points per element; meshio++'s cell_data "
                      "cannot represent per-point values (the MED ELNO/ELGA limit) -- skipped",
                      res.mName, gp->second.mNumPoints);
            continue;
        }

        std::size_t block = rStaged.mBlocks.size();
        for (std::size_t b = 0; b < rStaged.mBlocks.size(); ++b)
            if (rStaged.mBlocks[b].mMeshName == gp->second.mMeshName)
                block = b;
        if (block == rStaged.mBlocks.size()) {
            log::warn("gid: Gauss-point set '{}' is on mesh '{}', which matches no MESH block "
                      "-- result '{}' skipped",
                      res.mGaussName, gp->second.mMeshName, res.mName);
            continue;
        }

        // A cell_data array spanning several cell blocks is written as several
        // Result blocks SHARING ONE NAME, each against its own block's
        // Gauss-point set. So same-named results must be MERGED into one array
        // set here: AddCellData is insert-or-assign, and calling it once per
        // Result block would leave every block but the last zeroed.
        auto slot = pending.find(res.mName);
        if (slot == pending.end()) {
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b < rStaged.mBlocks.size(); ++b) {
                const std::size_t n = rStaged.mBlocks[b].mElemIds.size();
                blocks.emplace_back(DType::Float64,
                                    res.mNumComponents == 1
                                        ? std::vector<std::size_t>{n}
                                        : std::vector<std::size_t>{n, res.mNumComponents});
            }
            slot = pending.emplace(res.mName, std::move(blocks)).first;
            order.push_back(res.mName);
        }
        if (slot->second[block].Shape().size() >= 2 &&
            slot->second[block].Shape()[1] != res.mNumComponents) {
            log::warn("gid: result '{}' has {} components on mesh '{}' but a different width "
                      "elsewhere -- this block skipped",
                      res.mName, res.mNumComponents, gp->second.mMeshName);
            continue;
        }
        double* dst = slot->second[block].As<double>();
        for (std::size_t i = 0; i < res.mIds.size(); ++i) {
            auto it = elem_rows[block].find(res.mIds[i]);
            if (it == elem_rows[block].end())
                continue;  // a value for an element this block never defined
            for (std::size_t c = 0; c < res.mNumComponents; ++c)
                dst[it->second * res.mNumComponents + c] = res.mValues[i * res.mNumComponents + c];
        }
    }

    for (const std::string& name : order)
        rMesh.AddCellData(name, std::move(pending[name]));
}

Mesh gid_assemble(const GidStaged& rStaged, const std::vector<GidResult>& rResults,
                  const std::unordered_map<std::string, GidGaussSet>& rGauss, int time_step) {
    Mesh mesh;

    // Points, in ascending node id (std::map order). Ids may be gapped and
    // non-contiguous, so connectivity is remapped through this table -- the
    // abaqus/mdpa precedent.
    const std::size_t dim = rStaged.mAnyThirdCoord ? 3 : 2;
    NDArray points(DType::Float64, {rStaged.mNodes.size(), dim});
    double* pp = points.As<double>();
    std::map<std::int64_t, std::size_t> node_row;
    {
        std::size_t row = 0;
        for (const auto& kv : rStaged.mNodes) {
            node_row.emplace(kv.first, row);
            for (std::size_t c = 0; c < dim; ++c)
                pp[row * dim + c] = kv.second[c];
            ++row;
        }
    }
    mesh.AssignPoints(std::move(points));

    for (const GidBlock& block : rStaged.mBlocks) {
        const std::size_t nn = static_cast<std::size_t>(block.mNumNodes);
        const std::size_t ncells = block.mElemIds.size();
        NDArray conn(DType::Int64, {ncells, nn});
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < ncells * nn; ++i) {
            auto it = node_row.find(block.mConn[i]);
            if (it == node_row.end())
                throw ReadError("GiD: element in mesh '" + block.mMeshName +
                                "' references node id " + std::to_string(block.mConn[i]) +
                                ", which no Coordinates block defines");
            cp[i] = static_cast<std::int64_t>(it->second);
        }
        mesh.AddCellBlock(block.mMeshioType, std::move(conn));
    }

    // The material column round-trips through "gmsh:physical" -- the exact
    // inverse of the key write_gid consumes.
    // Only when the file actually carries material information. The binary and
    // HDF5 flavours ALWAYS emit a material column (unlike ASCII, which omits it
    // when there is nothing to say), so an all-zero column means "no materials"
    // -- surfacing it as a gmsh:physical array would invent data the source
    // mesh never had, and would make the three flavours disagree on round trip.
    bool any_material = false;
    for (const GidBlock& b : rStaged.mBlocks)
        for (std::int64_t v : b.mMaterial)
            any_material = any_material || v != 0;
    if (any_material) {
        std::vector<NDArray> mats;
        for (const GidBlock& b : rStaged.mBlocks) {
            NDArray a(DType::Int64, {b.mElemIds.size()});
            std::int64_t* ap = a.As<std::int64_t>();
            for (std::size_t i = 0; i < b.mElemIds.size(); ++i)
                ap[i] = i < b.mMaterial.size() ? b.mMaterial[i] : 0;
            mats.push_back(std::move(a));
        }
        mesh.AddCellData("gmsh:physical", std::move(mats));
    }

    gid_apply_results(mesh, rStaged, node_row, rResults, rGauss, time_step);
    return mesh;
}

// ---------------------------------------------------------------------------
// Flavour detection and decompression.

bool gid_looks_gzip(std::string_view bytes) {
    return bytes.size() >= 2 && static_cast<unsigned char>(bytes[0]) == 0x1f &&
           static_cast<unsigned char>(bytes[1]) == 0x8b;
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
/// Inflates a whole gzip stream. `uncompress()` cannot serve here: it speaks
/// the zlib wrapper, not gzip, and the decompressed size is not known up front.
std::string gid_gunzip(std::string_view bytes) {
    z_stream strm{};
    if (inflateInit2(&strm, 15 + 16) != Z_OK)
        throw ReadError("GiD: could not initialize gzip decompression");
    strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(bytes.data()));
    strm.avail_in = static_cast<uInt>(bytes.size());

    std::string out;
    std::vector<char> chunk(1 << 16);
    int rc = Z_OK;
    do {
        strm.next_out = reinterpret_cast<Bytef*>(chunk.data());
        strm.avail_out = static_cast<uInt>(chunk.size());
        rc = inflate(&strm, Z_NO_FLUSH);
        if (rc != Z_OK && rc != Z_STREAM_END) {
            inflateEnd(&strm);
            throw ReadError("GiD: gzip decompression failed");
        }
        out.append(chunk.data(), chunk.size() - strm.avail_out);
    } while (rc != Z_STREAM_END);
    inflateEnd(&strm);
    return out;
}
#endif  // MESHIOPLUSPLUS_HAS_ZLIB

#ifdef MESHIOPLUSPLUS_HAS_HDF5
// ---------------------------------------------------------------------------
// HDF5 flavour. A genuinely different traversal from the two text-shaped ones,
// so it gets its own function rather than a third cursor: gidpost stores the
// same logical model as HDF5 groups whose ATTRIBUTES carry what the ASCII
// header line carries. Layout (verified against a real written file):
//
//   /                       attrs: "GiD Post Results File", WriteStatus
//   /Meshes/<n>             attrs: Name, Dimension, ElemType, Nnode
//   /Meshes/<n>/Coordinates dataset
//   /Meshes/<n>/Elements    dataset
//   /GaussPoints/<n>        attrs: Name, ElemType, MeshName, GP_number, ...
//   /Results/<n>            attrs: Name, Analysis, Step, ResultType, ResultLocation
//
// Group names are integers as STRINGS ("1", "2", ...), so they are visited in
// HDF5 link-creation order, matching the order the writer emitted them.
/// Reads one column-group ("Coordinates", "Elements", a Results group) as a
/// list of 1-D column datasets named "1", "2", ... in creation order.
std::vector<NDArray> gid_h5_columns(hid_t loc) {
    std::vector<NDArray> cols;
    for (const std::string& key : h5::group_links_crt(loc))
        cols.push_back(h5::read_dataset(loc, key));
    return cols;
}

Mesh gid_read_hdf5(const std::string& rPath, const ReadOptions& rOptions) {
    h5::Hid file = h5::open_file_read(rPath);

    GidStaged staged;
    staged.mAnyThirdCoord = true;  // the HDF5 writer always stores three coordinate columns

    if (h5::exists(file, "Meshes")) {
        h5::Hid meshes = h5::open_group(file, "Meshes");
        for (const std::string& key : h5::group_links_crt(meshes)) {
            h5::Hid g = h5::open_group(meshes, key);
            GidBlock block;
            block.mMeshName = h5::has_attr(g, "Name") ? h5::read_attr_string(g, "Name") : key;
            // gidpost stores every attribute as a STRING, even the numeric ones.
            const std::string etype =
                h5::has_attr(g, "ElemType") ? h5::read_attr_string(g, "ElemType") : std::string();
            const std::string nnode_s =
                h5::has_attr(g, "Nnode") ? h5::read_attr_string(g, "Nnode") : std::string("0");
            const int nnode = static_cast<int>(gid_to_int(nnode_s, "Nnode"));
            block.mNumNodes = nnode;
            block.mMeshioType = gid_meshio_type(etype, nnode);

            if (h5::exists(g, "Coordinates")) {
                h5::Hid cg = h5::open_group(g, "Coordinates");
                const std::vector<NDArray> cols = gid_h5_columns(cg);
                // Column-oriented: 1 = ids, 2..4 = x, y, z.
                if (cols.size() >= 2) {
                    const std::size_t n = cols[0].Size();
                    for (std::size_t r = 0; r < n; ++r) {
                        std::array<double, 3> xyz{0.0, 0.0, 0.0};
                        for (std::size_t k = 0; k + 1 < cols.size() && k < 3; ++k)
                            xyz[k] = detail::read_double(cols[k + 1], r);
                        staged.mNodes.emplace(detail::read_int(cols[0], r), xyz);
                    }
                }
            }

            if (h5::exists(g, "Elements")) {
                h5::Hid eg = h5::open_group(g, "Elements");
                const std::vector<NDArray> cols = gid_h5_columns(eg);
                // 1 = element ids, then Nnode connectivity columns, then --
                // always, as in the binary flavour -- one material column.
                if (static_cast<int>(cols.size()) >= nnode + 1) {
                    const std::size_t n = cols[0].Size();
                    const bool has_mat = static_cast<int>(cols.size()) >= nnode + 2;
                    for (std::size_t r = 0; r < n; ++r) {
                        block.mElemIds.push_back(detail::read_int(cols[0], r));
                        for (int k = 0; k < nnode; ++k)
                            block.mConn.push_back(
                                detail::read_int(cols[static_cast<std::size_t>(k) + 1], r));
                        block.mMaterial.push_back(
                            has_mat ? detail::read_int(cols[static_cast<std::size_t>(nnode) + 1], r)
                                    : 0);
                    }
                    if (!has_mat)
                        block.mMaterial.clear();
                }
            }
            staged.mBlocks.push_back(std::move(block));
        }
    }

    std::unordered_map<std::string, GidGaussSet> gauss;
    if (h5::exists(file, "GaussPoints")) {
        h5::Hid gps = h5::open_group(file, "GaussPoints");
        for (const std::string& key : h5::group_links_crt(gps)) {
            h5::Hid g = h5::open_group(gps, key);
            GidGaussSet set;
            if (h5::has_attr(g, "MeshName"))
                set.mMeshName = h5::read_attr_string(g, "MeshName");
            if (h5::has_attr(g, "GP_number"))
                set.mNumPoints = static_cast<int>(
                    gid_to_int(h5::read_attr_string(g, "GP_number"), "GP_number"));
            gauss[h5::has_attr(g, "Name") ? h5::read_attr_string(g, "Name") : key] = set;
        }
    }

    std::vector<GidResult> results;
    if (h5::exists(file, "Results")) {
        h5::Hid res_root = h5::open_group(file, "Results");
        for (const std::string& key : h5::group_links_crt(res_root)) {
            h5::Hid g = h5::open_group(res_root, key);
            GidResult res;
            res.mName = h5::has_attr(g, "Name") ? h5::read_attr_string(g, "Name") : key;
            if (h5::has_attr(g, "Analysis"))
                res.mAnalysis = h5::read_attr_string(g, "Analysis");
            if (h5::has_attr(g, "Step"))
                res.mStep = gid_to_double(h5::read_attr_string(g, "Step"), "a step");
            if (h5::has_attr(g, "ResultType"))
                res.mType = h5::read_attr_string(g, "ResultType");
            if (h5::has_attr(g, "ResultLocation"))
                res.mLocation = h5::read_attr_string(g, "ResultLocation");
            if (h5::has_attr(g, "GaussPointsName"))
                res.mGaussName = h5::read_attr_string(g, "GaussPointsName");

            // 1 = ids, then one dataset per component.
            const std::vector<NDArray> cols = gid_h5_columns(g);
            if (cols.size() >= 2) {
                res.mNumComponents = cols.size() - 1;
                const std::size_t n = cols[0].Size();
                for (std::size_t r = 0; r < n; ++r) {
                    res.mIds.push_back(detail::read_int(cols[0], r));
                    for (std::size_t k = 0; k < res.mNumComponents; ++k)
                        res.mValues.push_back(detail::read_double(cols[k + 1], r));
                }
            }
            results.push_back(std::move(res));
        }
    }

    return gid_assemble(staged, results, gauss, rOptions.mTimeStep);
}
#endif  // MESHIOPLUSPLUS_HAS_HDF5

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
// ---------------------------------------------------------------------------
// Compressed-binary flavour (.post.bin).
//
// A gzip stream wrapping gidpost's own record layout (gidpostInt.c's
// CPostBinary):
//
//   int  0x91d              magic, and the endianness detector
//   then a stream of:
//     string  -> int length (INCLUDING the NUL) + that many bytes
//     int     -> raw 4-byte int
//     real    -> raw 4-byte FLOAT (gidpost narrows every double on the way out,
//                which is why binary round trips are float-precision, not a bug)
//
// The keyword vocabulary is the ASCII one with a distinct spelling for the
// three block openers: "Coordinates -1 Indexed", "Elements -1 Indexed",
// "Values -1 Indexed". `End Values` is preceded by a bare int -1.
class GidBinaryCursor {
public:
    explicit GidBinaryCursor(const std::string& rBytes) : mBytes(rBytes) {
        std::int32_t magic = 0;
        if (!RawInt(magic))
            throw ReadError("GiD: truncated binary file (no magic)");
        if (magic == 0x91d)
            mSwap = false;
        else if (ByteSwap(magic) == 0x91d)
            mSwap = true;  // foreign byte order, the ensight precedent
        else
            throw ReadError("GiD: not a GiD binary post file (bad magic)");
    }

    bool AtEnd() const { return mPos >= mBytes.size(); }

    std::int32_t NextInt() {
        std::int32_t v = 0;
        if (!RawInt(v))
            throw ReadError("GiD: truncated binary file (expected an int)");
        return mSwap ? ByteSwap(v) : v;
    }

    double NextReal() {
        float f = 0.0F;
        if (mPos + sizeof(float) > mBytes.size())
            throw ReadError("GiD: truncated binary file (expected a real)");
        std::memcpy(&f, mBytes.data() + mPos, sizeof(float));
        mPos += sizeof(float);
        if (mSwap) {
            std::int32_t bits = 0;
            std::memcpy(&bits, &f, sizeof(bits));
            bits = ByteSwap(bits);
            std::memcpy(&f, &bits, sizeof(f));
        }
        return static_cast<double>(f);
    }

    /// Reads a length-prefixed string record.
    std::string NextString() {
        const std::int32_t len = NextInt();
        if (len < 0 || mPos + static_cast<std::size_t>(len) > mBytes.size())
            throw ReadError("GiD: truncated binary file (bad string length)");
        std::string s(mBytes.data() + mPos, static_cast<std::size_t>(len));
        mPos += static_cast<std::size_t>(len);
        if (!s.empty() && s.back() == '\0')
            s.pop_back();  // the length includes the terminator
        return s;
    }

    std::size_t Tell() const { return mPos; }
    void Seek(std::size_t p) { mPos = p; }

    /**
     * @brief Consumes the next record iff it is the string @p pExpected.
     *
     * Binary blocks are terminated by a **string record** ("End Coordinates",
     * "End Elements", "End Values"), not by a sentinel value -- verified
     * against real inflated bytes. A row id and a string's length prefix are
     * both plain ints, so the only sound way to tell a further row from the
     * terminator is to try to decode the terminator and rewind on failure.
     * Comparing the FULL content (not just a plausible length) is what makes
     * this exact rather than heuristic: a false positive would need float
     * payload bytes to literally spell "End Coordinates".
     */
    bool TryTerminator(const char* pExpected) {
        const std::size_t save = mPos;
        std::int32_t len = 0;
        if (!RawInt(len)) {
            mPos = save;
            return false;
        }
        if (mSwap)
            len = ByteSwap(len);
        if (len <= 0 || mPos + static_cast<std::size_t>(len) > mBytes.size()) {
            mPos = save;
            return false;
        }
        std::string s(mBytes.data() + mPos, static_cast<std::size_t>(len));
        if (!s.empty() && s.back() == '\0')
            s.pop_back();
        if (s != pExpected) {
            mPos = save;
            return false;
        }
        mPos += static_cast<std::size_t>(len);
        return true;
    }

    /// Consumes the next int iff it equals @p want (the `-1` that precedes
    /// binary's "End Values" record).
    bool TryInt(std::int32_t want) {
        const std::size_t save = mPos;
        std::int32_t v = 0;
        if (!RawInt(v)) {
            mPos = save;
            return false;
        }
        if (mSwap)
            v = ByteSwap(v);
        if (v != want) {
            mPos = save;
            return false;
        }
        return true;
    }

private:
    bool RawInt(std::int32_t& rOut) {
        if (mPos + sizeof(std::int32_t) > mBytes.size())
            return false;
        std::memcpy(&rOut, mBytes.data() + mPos, sizeof(std::int32_t));
        mPos += sizeof(std::int32_t);
        return true;
    }
    static std::int32_t ByteSwap(std::int32_t v) {
        std::uint32_t u = static_cast<std::uint32_t>(v);
        u = ((u & 0x000000FFU) << 24) | ((u & 0x0000FF00U) << 8) | ((u & 0x00FF0000U) >> 8) |
            ((u & 0xFF000000U) >> 24);
        return static_cast<std::int32_t>(u);
    }

    const std::string& mBytes;
    std::size_t mPos = 0;
    bool mSwap = false;
};

/// The binary block openers carry a `-1 Indexed` suffix the ASCII ones do not.
bool gid_binary_keyword_is(const std::string& rRecord, const char* pKeyword) {
    const std::string kw(pKeyword);
    return rRecord == kw || rRecord.rfind(kw + " ", 0) == 0;
}

Mesh gid_read_binary(const std::string& rBytes, const ReadOptions& rOptions) {
    GidBinaryCursor cur(rBytes);
    GidStaged staged;
    staged.mAnyThirdCoord = true;  // the binary writer always emits three columns
    std::vector<GidResult> results;
    std::unordered_map<std::string, GidGaussSet> gauss;

    while (!cur.AtEnd()) {
        std::string rec;
        try {
            rec = cur.NextString();
        } catch (const ReadError&) {
            break;  // trailing padding: a clean end of stream
        }
        const std::vector<std::string> tok = gid_split(rec);
        if (tok.empty())
            continue;

        if (tok[0] == "MESH") {
            GidBlock block;
            std::string gid_type;
            int nnode = 0;
            block.mMeshName = tok.size() > 1 ? tok[1] : std::string();
            for (std::size_t i = 2; i + 1 < tok.size(); ++i) {
                if (tok[i] == "ElemType")
                    gid_type = tok[i + 1];
                else if (tok[i] == "Nnode")
                    nnode = static_cast<int>(gid_to_int(tok[i + 1], "Nnode"));
            }
            if (gid_type.empty() || nnode <= 0)
                throw ReadError("GiD: malformed binary MESH header: " + rec);
            block.mNumNodes = nnode;
            block.mMeshioType = gid_meshio_type(gid_type, nnode);
            staged.mBlocks.push_back(std::move(block));
            continue;
        }

        if (gid_binary_keyword_is(rec, "Coordinates")) {
            while (!cur.TryTerminator("End Coordinates")) {
                const std::int32_t id = cur.NextInt();
                std::array<double, 3> xyz{0.0, 0.0, 0.0};
                for (std::size_t k = 0; k < 3; ++k)
                    xyz[k] = cur.NextReal();
                staged.mNodes.emplace(static_cast<std::int64_t>(id), xyz);
            }
            continue;
        }

        if (gid_binary_keyword_is(rec, "Elements")) {
            if (staged.mBlocks.empty())
                throw ReadError("GiD: binary Elements block before any MESH header");
            GidBlock& block = staged.mBlocks.back();
            const int nn = block.mNumNodes;
            while (!cur.TryTerminator("End Elements")) {
                block.mElemIds.push_back(static_cast<std::int64_t>(cur.NextInt()));
                for (int k = 0; k < nn; ++k)
                    block.mConn.push_back(static_cast<std::int64_t>(cur.NextInt()));
                // gidpost's binary element writer ALWAYS emits the material
                // column (verified in the inflated bytes), so unlike ASCII
                // there is no row-width ambiguity to resolve here.
                block.mMaterial.push_back(static_cast<std::int64_t>(cur.NextInt()));
            }
            continue;
        }

        if (tok[0] == "GaussPoints") {
            GidGaussSet set;
            const std::string gp_name = tok.size() > 1 ? tok[1] : std::string();
            if (tok.size() >= 5)
                set.mMeshName = tok[4];
            while (!cur.AtEnd()) {
                const std::size_t save = cur.Tell();
                std::string inner;
                try {
                    inner = cur.NextString();
                } catch (const ReadError&) {
                    break;
                }
                if (inner.rfind("End GaussPoints", 0) == 0)
                    break;
                const std::vector<std::string> it = gid_split(inner);
                if (it.size() >= 5 && it[0] == "Number" && it[1] == "Of")
                    set.mNumPoints = static_cast<int>(gid_to_int(it[4], "a Gauss point count"));
                else if (it.empty())
                    cur.Seek(save);
            }
            gauss[gp_name] = set;
            continue;
        }

        if (tok[0] == "Result") {
            if (tok.size() < 6)
                throw ReadError("GiD: malformed binary Result header: " + rec);
            GidResult res;
            res.mName = tok[1];
            res.mAnalysis = tok[2];
            res.mStep = gid_to_double(tok[3], "a result step");
            res.mType = tok[4];
            res.mLocation = tok[5];
            if (tok.size() >= 7)
                res.mGaussName = tok[6];
            // The binary stream carries no component count either, and unlike
            // ASCII there is no row width to infer it from -- so the declared
            // type's own canonical dimension is the only source.
            res.mNumComponents = res.mType == "Scalar"   ? 1
                                 : res.mType == "Vector" ? 3
                                 : res.mType == "Matrix" ? 6
                                                         : 0;
            if (res.mNumComponents == 0)
                throw ReadError("GiD: binary result '" + res.mName + "' has type '" + res.mType +
                                "', whose component count this reader cannot infer");
            std::string opener;
            while (!cur.AtEnd()) {
                opener = cur.NextString();
                if (gid_binary_keyword_is(opener, "Values"))
                    break;
            }
            while (!cur.TryInt(-1) && !cur.TryTerminator("End Values")) {
                res.mIds.push_back(static_cast<std::int64_t>(cur.NextInt()));
                for (std::size_t k = 0; k < res.mNumComponents; ++k)
                    res.mValues.push_back(cur.NextReal());
            }
            cur.TryTerminator("End Values");  // consume it if the -1 ended the loop
            results.push_back(std::move(res));
            continue;
        }
        // Anything else is a keyword record we do not map.
    }

    return gid_assemble(staged, results, gauss, rOptions.mTimeStep);
}
#endif  // MESHIOPLUSPLUS_HAS_ZLIB

std::string gid_missing_flavour_message(GidMode mode) {
    if (mode == GidMode::Hdf5)
        return "meshio++: reading the 'gid' HDF5 flavour needs a build with "
               "-DMESHIOPLUSPLUS_WITH_HDF5=ON";
    return "meshio++: reading the 'gid' compressed flavours needs a build with "
           "-DMESHIOPLUSPLUS_WITH_ZLIB=ON";
}

}  // namespace

// ---------------------------------------------------------------------------
// Public entry points.

bool gid_readable(GidMode mode) {
    switch (mode) {
        case GidMode::Hdf5:
#ifdef MESHIOPLUSPLUS_HAS_HDF5
            return true;
#else
            return false;
#endif
        case GidMode::Binary:
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
            return true;
#else
            return false;
#endif
        default:
            return true;  // Ascii (and Auto, which resolves to it)
    }
}

Mesh read_gid(const std::string& rPath, const ReadOptions& rOptions) {
    const GidMode resolved = gid_resolve_mode(rPath, GidMode::Auto);

    if (resolved == GidMode::Hdf5) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        return gid_read_hdf5(rPath, rOptions);
#else
        throw ReadError(gid_missing_flavour_message(GidMode::Hdf5));
#endif
    }

    if (resolved == GidMode::Binary) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        const detail::FileSource src(rPath, rOptions.mMmap);
        return gid_read_binary(gid_gunzip(src.View()), rOptions);
#else
        throw ReadError(gid_missing_flavour_message(GidMode::Binary));
#endif
    }

    // ASCII (possibly gzipped: GiD_PostAsciiZipped writes the same text
    // through gzprintf, and the extension cannot say so -- only the bytes can).
    const auto paths = gid_ascii_paths(rPath);

    const detail::FileSource mesh_src = [&] {
        try {
            return detail::FileSource(paths.first, rOptions.mMmap);
        } catch (const ReadError&) {
            throw ReadError("GiD: could not open geometry file: " + paths.first);
        }
    }();

    std::string mesh_owned;
    std::string_view mesh_text = mesh_src.View();
    if (gid_looks_gzip(mesh_text)) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
        mesh_owned = gid_gunzip(mesh_text);
        mesh_text = mesh_owned;
#else
        throw ReadError(gid_missing_flavour_message(GidMode::Binary));
#endif
    }

    GidStaged staged;
    gid_parse_mesh_text(mesh_text, staged);

    // The results sibling is OPTIONAL (the triangle .node/.ele precedent): a
    // geometry file with no results reads back as geometry only.
    std::vector<GidResult> results;
    std::unordered_map<std::string, GidGaussSet> gauss;
    try {
        const detail::FileSource res_src(paths.second, rOptions.mMmap);
        std::string res_owned;
        std::string_view res_text = res_src.View();
        if (gid_looks_gzip(res_text)) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
            res_owned = gid_gunzip(res_text);
            res_text = res_owned;
#else
            throw ReadError(gid_missing_flavour_message(GidMode::Binary));
#endif
        }
        gid_parse_res_text(res_text, results, gauss);
    } catch (const ReadError&) {
        // Absent or unreadable results file: geometry only. A malformed one is
        // deliberately not fatal either -- the geometry is still good, and the
        // alternative is refusing a file GiD itself opens.
        results.clear();
        gauss.clear();
    }

    return gid_assemble(staged, results, gauss, rOptions.mTimeStep);
}

MeshMetadata read_gid_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    // Deliberately minimal: this declines (by throwing) for everything except
    // the plain ASCII flavour, and registry_read_metadata then falls back to a
    // full read. Declining costs a slower answer, never a failed one.
    if (gid_resolve_mode(rPath, GidMode::Auto) != GidMode::Ascii)
        throw ReadError("GiD: no cheap metadata path for this flavour");

    const auto paths = gid_ascii_paths(rPath);
    const detail::FileSource src(paths.first, rOptions.mMmap);
    if (gid_looks_gzip(src.View()))
        throw ReadError("GiD: no cheap metadata path for a gzipped file");

    // Walk MESH headers and count rows, without building any array.
    GidStaged staged;
    gid_parse_mesh_text(src.View(), staged);

    MeshMetadata meta;
    meta.mNumPoints = staged.mNodes.size();
    meta.mPointDim = staged.mAnyThirdCoord ? 3 : 2;
    for (const GidBlock& b : staged.mBlocks) {
        CellBlockInfo info;
        info.mType = b.mMeshioType;
        info.mNumCells = b.mElemIds.size();
        info.mNodesPerCell = static_cast<std::size_t>(b.mNumNodes);
        meta.mCellBlocks.push_back(info);
    }
    return meta;
}

}  // namespace meshioplusplus
