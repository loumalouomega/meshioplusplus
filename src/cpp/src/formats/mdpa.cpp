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
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Lexing helpers. Every anonymous-namespace symbol here is `mdpa_`-prefixed
// because the single-header amalgamation concatenates all of src/cpp/src.
// ---------------------------------------------------------------------------

std::string_view mdpa_strip(std::string_view S) {
    const std::size_t b = S.find_first_not_of(" \t\r\n");
    if (b == std::string_view::npos)
        return {};
    const std::size_t e = S.find_last_not_of(" \t\r\n");
    return S.substr(b, e - b + 1);
}

/// The line with any `//` comment removed, then stripped.
std::string_view mdpa_clean(std::string_view S) {
    const std::size_t c = S.find("//");
    return mdpa_strip(c == std::string_view::npos ? S : S.substr(0, c));
}

/// The line's blank-separated tokens, as views into it (detail/text_cursor.hpp).
std::vector<std::string_view> mdpa_tokens(std::string_view S) {
    return detail::split_blanks(S);
}

bool mdpa_starts_with(std::string_view S, std::string_view Prefix) {
    return S.size() >= Prefix.size() && S.compare(0, Prefix.size(), Prefix) == 0;
}

/// strtoll(…, 10) over the whole token.
bool mdpa_parse_int(std::string_view S, std::int64_t& rOut) {
    return detail::parse_int_token(S, rOut);
}

/// Membership ids must not saturate: that would name a different entity.
bool mdpa_parse_member_id(std::string_view S, std::int64_t& rOut) {
    if (!S.empty() && S.front() == '+') {
        S.remove_prefix(1);
        if (S.empty() || S.front() < '0' || S.front() > '9')
            return false;
    }
    if (S.empty())
        return false;
    const auto result = std::from_chars(S.data(), S.data() + S.size(), rOut);
    return result.ec == std::errc{} && result.ptr == S.data() + S.size();
}

/// parse_double over the whole token.
bool mdpa_parse_double(std::string_view S, double& rOut) {
    return detail::parse_double_token(S, rOut);
}

// ---------------------------------------------------------------------------
// Kratos <-> VTK node orderings.
//
// `mdpa_kratos_node_order(type)[i]` is the meshio/VTK slot the node in Kratos
// slot `i` belongs in. Reading therefore scatters (`out[table[i]] = in[i]`) and
// writing gathers (`out[j] = in[table[j]]`). The tables are the "mdpa" entries
// of the node-ordering registry (detail/node_order.cpp), read from Kratos's own
// geometry classes.
// ---------------------------------------------------------------------------

const std::vector<int>& mdpa_kratos_node_order(CellType type) {
    static const std::vector<int> none;
    const detail::NodeOrder* order = detail::node_order("mdpa", cell_type_name(type));
    return order ? order->mFromMeshio : none;
}

/**
 * @brief Resolve a Kratos entity name to a meshio cell type.
 *
 * Exact lookup first (`kratos_names.hpp` owns the tables), then the longest
 * suffix that resolves — which is what makes application-specific names such
 * as `SmallDisplacementElement3D4N` work, mirroring the Python reference's
 * substring fallback for the cases that occur in practice.
 */
CellType mdpa_entity_cell_type(const std::string& rName) {
    // kratos_names.hpp owns both the tables and the exact-then-longest-suffix
    // rule, so this and ModelPart::CreateNewElement cannot disagree about which
    // names a deck may use.
    const CellType t = cell_type_from_kratos_name_or_suffix(rName);
    if (t != CellType::Custom || rName.size() < 3)
        return t;
    // Then the Python reference's last resort: the longest Kratos name found
    // anywhere inside it (`Triangle2D3N` -> `Triangle2D3`). Plain meshio names
    // (`quad`, `line`) are not Kratos names and never match here; a wrong
    // guess still fails the per-row node-count check.
    for (std::size_t len = rName.size() - 1; len >= 2; --len)
        for (std::size_t i = 0; i + len <= rName.size(); ++i) {
            const std::string sub = rName.substr(i, len);
            if (cell_type_from_name(sub) != CellType::Custom)
                continue;
            const CellType s = cell_type_from_kratos_name(sub);
            if (s != CellType::Custom)
                return s;
        }
    return CellType::Custom;
}

// ---------------------------------------------------------------------------
// Reader staging
// ---------------------------------------------------------------------------

struct MdpaBlock {
    std::string mType;
    /// The Kratos entity name the file spelled, kept so it can be written back.
    std::string mEntityName;
    std::size_t mNodes = 0;
    bool mIsCondition = false;
    /// Raw *file* node ids until `mdpa_read_impl`'s materialize pass resolves
    /// them to point rows -- the Nodes block is not required to precede this
    /// one, and deferring keeps the file id available for the error message.
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mProps;
    /// Each row's raw file id, in append order -- captured unconditionally
    /// (parallel to `mProps`) so the materialize pass can decide whether it is
    /// worth keeping without having re-derived it.
    std::vector<std::int64_t> mFileIds;
    std::size_t mCount = 0;
};

/// One parsed row of a `NodalData` / `ElementalData` / `ConditionalData` block.
struct MdpaDataRow {
    std::size_t mBlock = 0;  ///< cell block index (0 for nodal data)
    std::size_t mRow = 0;    ///< point index, or row within the cell block
    int mFixed = -1;         ///< the optional leading 0/1 column, -1 if absent
    std::vector<double> mValues;
};

/// A cursor over the file's lines, so every block parser advances one index.
using MdpaCursor = detail::RecordCursor<std::string_view>;

/// The refusal's tail when the caller could have kept the construct in an MdpaInfo.
constexpr const char* kMdpaNeedsInfo =
    " needs an MdpaInfo to be kept (read_mdpa(path, info), or mio_read_with_info on the flat "
    "ABI); set ReadOptions::mLenient to skip it instead";
/// The refusal's tail for what not even an MdpaInfo holds.
constexpr const char* kMdpaUnsupported =
    " is not supported by the C++ reader (set ReadOptions::mLenient to skip it instead)";

/// Consume the rest of a block, ignoring blank/comment-only lines.
void mdpa_consume_block(MdpaCursor& rCur, const std::string& rEnd) {
    while (!rCur.Done())
        if (mdpa_clean(rCur.Next()) == rEnd)
            return;
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief Reject a construct this reader cannot represent, or skip it.
 *
 * Strict (the default) throws naming it, so a caller with a Python fallback can
 * take it and a caller without one at least learns what was in the file.
 * `Lenient` warns, records it in @p pInfo and consumes the block instead --
 * which is the only way a production deck reads at all where there is no Python
 * (the C API, Fortran, Julia, R, WASM, the native CLI).
 */
void mdpa_reject_or_skip(MdpaCursor& rCur, const std::string& rEnd, const std::string& rWhat,
                         bool Lenient, MdpaInfo* pInfo) {
    if (!Lenient)
        throw ReadError("MDPA: " + rWhat + (pInfo ? kMdpaUnsupported : kMdpaNeedsInfo));
    log::warn("mdpa: skipping {} (ReadOptions::mLenient)", rWhat);
    if (pInfo)
        pInfo->mSkippedConstructs.push_back(rWhat);
    mdpa_consume_block(rCur, rEnd);
}

/// Consume a block that must be empty; a non-empty one goes through the above.
void mdpa_expect_empty_block(MdpaCursor& rCur, const std::string& rEnd, const std::string& rWhat,
                             bool Lenient, MdpaInfo* pInfo) {
    while (!rCur.Done()) {
        const std::string_view line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == rEnd)
            return;
        if (!Lenient)
            throw ReadError("MDPA: " + rWhat + " (offending line: '" + std::string(line) + "')" +
                            (pInfo ? kMdpaUnsupported : kMdpaNeedsInfo));
        log::warn("mdpa: skipping {} (ReadOptions::mLenient)", rWhat);
        if (pInfo)
            pInfo->mSkippedConstructs.push_back(rWhat);
        mdpa_consume_block(rCur, rEnd);
        return;
    }
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief One `KEY value` line: a number becomes a Float64 `{1}`, anything else text.
 *
 * Key plus the rest of the line, the Python reference's `split(None, 1)`, so
 * both readers agree on where the value starts. False for a valueless line.
 */
bool mdpa_parse_kv_line(std::string_view Line, PropertyValue& rOut) {
    const std::size_t sep = Line.find_first_of(" \t");
    if (sep == std::string_view::npos)
        return false;
    rOut = PropertyValue{};
    rOut.mKey = Line.substr(0, sep);
    const std::string_view rest = mdpa_strip(Line.substr(sep + 1));
    double scalar = 0.0;
    if (mdpa_parse_double(rest, scalar)) {
        NDArray a(DType::Float64, {1});
        a.As<double>()[0] = scalar;
        rOut.mValues = std::move(a);
    } else {
        rOut.mText = rest;
    }
    return true;
}

/**
 * @brief Parse an inline `Begin Table <args>` inside a `Properties` body.
 *
 * @p rHeader is the whole header line; everything after `Begin Table` becomes
 * the entry's `mKey` verbatim, so the block re-emits with its id and variable
 * names unchanged. Rows are whitespace-separated numbers; the first usable row
 * fixes the column count and a row that disagrees is warned about and skipped.
 */
PropertyValue mdpa_parse_property_table(MdpaCursor& rCur, std::string_view rHeader) {
    PropertyValue out;
    out.mIsTable = true;
    out.mKey = mdpa_strip(rHeader.substr(std::string("Begin Table").size()));

    std::vector<double> values;
    std::size_t ncols = 0;
    bool terminated = false;
    while (!rCur.Done()) {
        const std::string_view line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == "End Table") {
            terminated = true;
            break;
        }
        const std::vector<std::string_view> toks = mdpa_tokens(line);
        std::vector<double> row;
        row.reserve(toks.size());
        bool ok = true;
        for (const std::string_view tok : toks) {
            double v = 0.0;
            if (!mdpa_parse_double(tok, v)) {
                ok = false;
                break;
            }
            row.push_back(v);
        }
        if (!ok) {
            log::warn("mdpa: skipping non-numeric Table row: {}", line);
            continue;
        }
        if (ncols == 0)
            ncols = row.size();
        if (row.size() != ncols) {
            log::warn("mdpa: skipping Table row with {} values (expected {}): {}", row.size(),
                      ncols, line);
            continue;
        }
        values.insert(values.end(), row.begin(), row.end());
    }
    if (!terminated)
        throw ReadError("MDPA: EOF before 'End Table'");

    const std::size_t nrows = ncols ? values.size() / ncols : 0;
    NDArray a(DType::Float64, {nrows, ncols});
    double* p = a.As<double>();
    for (std::size_t i = 0; i < values.size(); ++i)
        p[i] = values[i];
    out.mValues = std::move(a);
    return out;
}

/**
 * @brief Parse a `Begin Properties <id>` body.
 *
 * Never throws on content: a plain number becomes a Float64 scalar, an inline
 * table an `(n, k)` array, and everything else -- a constitutive-law name, a
 * bracketed vector or matrix -- is kept verbatim as text, which is both
 * lossless and what the pure-Python reference does.
 */
PropertySet mdpa_parse_properties(MdpaCursor& rCur, std::string_view rHeader) {
    PropertySet out;
    const std::vector<std::string_view> head = mdpa_tokens(rHeader);
    if (head.size() < 3 || !mdpa_parse_int(head[2], out.mId)) {
        log::warn("mdpa: Properties block with no readable id, using 0: {}", rHeader);
        out.mId = 0;
    }

    while (true) {
        if (rCur.Done())
            throw ReadError("MDPA: EOF before 'End Properties'");
        const std::string_view line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == "End Properties")
            break;
        if (mdpa_starts_with(line, "Begin Table")) {
            out.mValues.push_back(mdpa_parse_property_table(rCur, line));
            continue;
        }
        PropertyValue v;
        if (!mdpa_parse_kv_line(line, v)) {
            log::warn("mdpa: skipping valueless Properties line: {}", line);
            continue;
        }
        out.mValues.push_back(std::move(v));
    }
    return out;
}

/// A `KEY value` block body (`SubModelPartData`, `MeshData`) up to @p rEnd.
std::vector<PropertyValue> mdpa_parse_kv_block(MdpaCursor& rCur, const std::string& rEnd) {
    std::vector<PropertyValue> out;
    while (!rCur.Done()) {
        const std::string_view line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == rEnd)
            return out;
        PropertyValue v;
        if (!mdpa_parse_kv_line(line, v)) {
            log::warn("mdpa: skipping valueless line in {} block: {}", rEnd.substr(4), line);
            continue;
        }
        out.push_back(std::move(v));
    }
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief Keep a block verbatim: its header, its body lines and its terminator.
 *
 * The body is the file's own text (a trailing `\r` aside), so the block writes
 * back exactly as it was read; only the terminator is matched comment-free.
 */
MdpaRawBlock mdpa_capture_raw_block(MdpaCursor& rCur, std::string_view Header,
                                    const std::string& rEnd) {
    MdpaRawBlock out;
    out.mHeader = Header;
    out.mEnd = rEnd;
    while (!rCur.Done()) {
        std::string_view raw = rCur.Next();
        if (mdpa_clean(raw) == rEnd)
            return out;
        if (!raw.empty() && raw.back() == '\r')
            raw.remove_suffix(1);
        out.mBody.append(raw.data(), raw.size());
        out.mBody += '\n';
    }
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief Parse a `*Data` block body, mirroring `_parse_generic_data_block`.
 *
 * @param rCur cursor positioned just after the `Begin ...` header
 * @param rEnd the terminating token, e.g. `"End NodalData"`
 * @param nodal whether the optional leading `fixed` column may appear
 * @param rResolve maps a 1-based entity id to `(block, row)`; returns false to
 *        skip the line (unknown id)
 * @param rRows out: the parsed rows
 * @param rHasFixed out: whether any row carried a `fixed` column
 * @return the number of value components (0 for a membership flag), or -1 when
 *         the block held no usable line
 */
int mdpa_parse_data_block(
    MdpaCursor& rCur, const std::string& rEnd, bool nodal,
    const std::function<bool(std::int64_t, std::size_t&, std::size_t&)>& rResolve,
    std::vector<MdpaDataRow>& rRows, bool& rHasFixed) {
    int nc = -1;
    bool terminated = false;
    std::vector<std::string_view> toks;  // reused: one allocation per block, not per row
    while (!rCur.Done()) {
        const std::string_view raw = mdpa_strip(rCur.Next());
        if (raw == rEnd) {
            terminated = true;
            break;
        }
        const std::string_view line = mdpa_clean(raw);
        if (line.empty())
            continue;
        if (line == rEnd) {
            terminated = true;
            break;
        }
        detail::split_blanks(line, toks);
        std::int64_t id = 0;
        if (toks.empty())
            continue;
        if (!mdpa_parse_int(toks[0], id)) {
            log::warn("mdpa: skipping data line with non-integer id: {}", line);
            continue;
        }
        std::size_t block = 0, row = 0;
        if (!rResolve(id, block, row)) {
            log::warn("mdpa: skipping data for unknown entity id {}", id);
            continue;
        }

        // Split the remaining tokens into an optional `fixed` flag and values.
        const std::size_t rest = toks.size() - 1;
        int fixed = -1;
        std::size_t first_value = 1;
        std::int64_t flag = 0;
        const bool flag_like =
            nodal && rest > 0 && mdpa_parse_int(toks[1], flag) && (flag == 0 || flag == 1);
        if (nc < 0) {  // first usable line defines the layout
            if (flag_like && rest > 1) {
                fixed = static_cast<int>(flag);
                first_value = 2;
            }
            nc = static_cast<int>(toks.size() - first_value);
        } else {
            const std::size_t want = static_cast<std::size_t>(nc);
            if (flag_like && rest == want + 1) {
                fixed = static_cast<int>(flag);
                first_value = 2;
            } else if (rest != want) {
                log::warn("mdpa: skipping data line with {} values (expected {}): {}", rest, want,
                          line);
                continue;
            }
        }
        if (fixed >= 0)
            rHasFixed = true;

        MdpaDataRow r;
        r.mBlock = block;
        r.mRow = row;
        r.mFixed = fixed;
        bool ok = true;
        for (std::size_t j = first_value; j < toks.size(); ++j) {
            double v = 0.0;
            if (!mdpa_parse_double(toks[j], v)) {
                ok = false;
                break;
            }
            r.mValues.push_back(v);
        }
        if (!ok) {
            log::warn("mdpa: skipping data line with non-numeric value: {}", line);
            continue;
        }
        rRows.push_back(std::move(r));
    }
    if (!terminated)
        throw ReadError("MDPA: EOF before '" + rEnd + "'");
    return nc;
}

/// Build one data array of @p count rows and @p nc components from @p rRows.
NDArray mdpa_data_array(const std::vector<MdpaDataRow>& rRows, std::size_t block, std::size_t count,
                        int nc) {
    if (nc == 0) {  // membership flag
        NDArray a(DType::Int64, {count});
        std::int64_t* p = a.As<std::int64_t>();
        for (std::size_t i = 0; i < count; ++i)
            p[i] = 0;
        for (const auto& r : rRows)
            if (r.mBlock == block && r.mRow < count)
                p[r.mRow] = 1;
        return a;
    }
    const std::size_t k = static_cast<std::size_t>(nc);
    std::vector<std::size_t> shape;
    if (k == 1)
        shape = {count};
    else
        shape = {count, k};
    NDArray a(DType::Float64, shape);
    double* p = a.As<double>();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < count * k; ++i)
        p[i] = nan;
    for (const auto& r : rRows) {
        if (r.mBlock != block || r.mRow >= count || r.mValues.size() != k)
            continue;
        for (std::size_t j = 0; j < k; ++j)
            p[r.mRow * k + j] = r.mValues[j];
    }
    return a;
}

}  // namespace

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

namespace {

/**
 * @brief The one parse body behind all three `read_mdpa` overloads.
 *
 * @param rPath   file to read
 * @param Lenient `ReadOptions::mLenient`
 * @param pInfo   where to put what the `Mesh` cannot hold, or null to drop it
 */
Mesh mdpa_read_impl(const std::string& rPath, bool Lenient, MdpaInfo* pInfo) {
    // The file read once; its lines are views into it (detail/text_cursor.hpp),
    // not a string each.
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    const std::vector<std::string_view> lines = detail::split_lines(source.View());

    MdpaCursor cur{lines};

    std::vector<double> coords;  // flat (n, 3)
    std::size_t num_points = 0;
    std::vector<MdpaBlock> blocks;
    // The Properties bodies, staged for the Mesh. Kept in file order here; the
    // mesh canonicalizes them to ascending id, which is why MdpaInfo remains
    // the way to preserve an unusual declaration order verbatim.
    std::vector<PropertySet> property_sets;
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> element_ids,
        condition_ids;
    // The node half of the id maps above -- but materialized LAZILY. Kratos node
    // ids are 1..n in file order in the overwhelming majority of decks, and those
    // are the million-node ones; until an id turns up that is not `row + 1`,
    // "row == id - 1" IS the map and building one would be pure overhead. The
    // moment one does, the identity entries read so far are back-filled and the
    // map takes over for the rest of the file. This mirrors abaqus.cpp's
    // `mPointIds` and unv.cpp's `label_to_index`, minus their unconditional cost.
    std::unordered_map<std::int64_t, std::size_t> node_ids;  // file id -> point row
    bool node_ids_dense = true;  // ids so far are exactly 1..num_points, in order
    // Every node's raw file id, in file (row) order -- captured unconditionally
    // (cheap: one push_back per row already being appended) so it is available
    // at materialize time regardless of whether `node_ids_dense` ever flips.
    std::vector<std::int64_t> raw_node_ids;
    // Same "dense" idea as node ids, but tracked directly rather than lazily:
    // elements and conditions each have their own independent 1-based counter,
    // spanning every block of that kind in file order (not per block), which is
    // exactly how the writer numbers them. The moment either counter's next
    // expected value disagrees with a row's actual id, ids are no longer the
    // trivial "renumber from 1" case and the original ones are worth keeping.
    std::int64_t next_element_id = 1, next_condition_id = 1;
    bool entities_dense = true;
    std::map<std::string, NDArray> field_data;

    // Staged data arrays, materialized after every block is known.
    struct StagedData {
        std::string mName;
        int mComponents = 0;
        bool mHasFixed = false;
        std::vector<MdpaDataRow> mRows;
    };
    std::vector<StagedData> point_data, cell_data;

    // SubModelParts, accumulated by (hierarchical) name.
    struct StagedSmp {
        std::vector<std::int64_t> mNodes;
        std::vector<std::pair<std::size_t, std::size_t>> mCells;
    };
    std::map<std::string, StagedSmp> smps;
    std::vector<std::string> smp_stack;

    // Side-channel content, staged only when there is an MdpaInfo to receive
    // it. Geometry and Mesh-block node references stay raw file ids until the
    // materialize pass, for the same reason MdpaBlock::mConn does.
    struct StagedGeometry {
        std::string mName;
        CellType mType = CellType::Custom;
        std::size_t mNodes = 0;
        std::vector<std::int64_t> mConn;  // raw file node ids, meshio order
        std::vector<std::int64_t> mIds;
    };
    std::vector<StagedGeometry> geometries;
    std::vector<std::vector<std::int64_t>> mesh_block_nodes;      // raw ids, per MdpaMeshBlock
    std::unordered_map<std::string, std::size_t> smp_info_index;  // name -> mSubModelParts slot
    auto smp_info = [&](const std::string& rName) -> MdpaSubModelPart& {
        const auto it = smp_info_index.find(rName);
        if (it != smp_info_index.end())
            return pInfo->mSubModelParts[it->second];
        smp_info_index.emplace(rName, pInfo->mSubModelParts.size());
        pInfo->mSubModelParts.push_back(MdpaSubModelPart{rName, {}, {}, {}, {}});
        return pInfo->mSubModelParts.back();
    };

    auto smp_name = [&]() -> std::string {
        std::string out;
        for (std::size_t i = 0; i < smp_stack.size(); ++i) {
            if (i)
                out += "/";
            out += smp_stack[i];
        }
        return out;
    };

    auto read_id_list = [&](const std::string& rEnd, std::vector<std::int64_t>& rOut,
                            bool Strict = false) {
        while (!cur.Done()) {
            const std::string_view line = mdpa_clean(cur.Next());
            if (line.empty())
                continue;
            if (line == rEnd)
                return;
            std::int64_t v = 0;
            if (!(Strict ? mdpa_parse_member_id(line, v) : mdpa_parse_int(line, v))) {
                if (Strict)
                    throw ReadError("MDPA: non-integer id in " + rEnd + ": " + std::string(line));
                log::warn("mdpa: skipping non-integer id in {}: {}", rEnd, line);
                continue;
            }
            rOut.push_back(v);
        }
        throw ReadError("MDPA: EOF before '" + rEnd + "'");
    };

    // File node id -> point row. False when the file defines no such node. The
    // dense branch is the pre-v9.13.0 arithmetic, bit for bit.
    auto node_row = [&](std::int64_t Id, std::size_t& rRow) -> bool {
        if (node_ids_dense) {
            if (Id < 1 || static_cast<std::size_t>(Id) > num_points)
                return false;
            rRow = static_cast<std::size_t>(Id) - 1;
            return true;
        }
        const auto it = node_ids.find(Id);
        if (it == node_ids.end())
            return false;
        rRow = it->second;
        return true;
    };

    while (!cur.Done()) {
        const std::string_view raw = mdpa_strip(cur.Next());
        const std::string_view line = mdpa_clean(raw);
        if (line.empty())
            continue;

        if (line == "Begin ModelPartData") {
            while (true) {
                if (cur.Done())
                    throw ReadError("MDPA: EOF before 'End ModelPartData'");
                const std::string_view e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End ModelPartData")
                    break;
                PropertyValue v;
                if (!mdpa_parse_kv_line(e, v)) {
                    log::warn("mdpa: skipping malformed ModelPartData line: {}", e);
                    continue;
                }
                if (v.IsText()) {
                    // NDArray has no string dtype: the value goes to the side
                    // channel, or -- without one -- is refused by name.
                    if (pInfo) {
                        pInfo->mModelPartData.push_back(std::move(v));
                        continue;
                    }
                    const std::string what =
                        "a non-numeric ModelPartData value for '" + v.mKey + "'";
                    if (!Lenient)
                        throw ReadError("MDPA: " + what + kMdpaNeedsInfo);
                    log::warn("mdpa: skipping {} (ReadOptions::mLenient)", what);
                    continue;
                }
                field_data[v.mKey] = std::move(v.mValues);
            }
        } else if (line == "Begin Nodes") {
            if (num_points)
                throw ReadError("MDPA: more than one Nodes block");
            bool terminated = false;
            // The format carries no node count: the block's lines (views, so
            // cheap to look ahead over) bound it, and size the coordinates.
            {
                std::size_t rows = 0;
                for (std::size_t k = cur.Pos(); k < lines.size(); ++k) {
                    const std::string_view ahead = mdpa_clean(lines[k]);
                    if (ahead == "End Nodes")
                        break;
                    rows += ahead.empty() ? 0 : 1;
                }
                coords.reserve(coords.size() + 3 * rows);
                raw_node_ids.reserve(raw_node_ids.size() + rows);
            }
            std::vector<std::string_view> t;  // reused: one allocation per block, not per row
            while (!cur.Done()) {
                const std::string_view e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End Nodes") {
                    terminated = true;
                    break;
                }
                detail::split_blanks(e, t);
                if (t.size() < 3)
                    throw ReadError("MDPA: node line with fewer than 3 coordinates: " +
                                    std::string(e));
                // An id-less row (`x y z`) takes its position as its id, which is
                // exactly what "connectivity is 1-based into row order" already
                // meant for a fully id-less file -- so such a file never leaves
                // the dense path and reads byte-identically to before.
                std::int64_t id = static_cast<std::int64_t>(num_points) + 1;
                if (t.size() >= 4 && !mdpa_parse_int(t[0], id))
                    throw ReadError("MDPA: non-integer node id: " + std::string(e));
                if (node_ids_dense && id != static_cast<std::int64_t>(num_points) + 1) {
                    node_ids.reserve(num_points * 2 + 16);
                    for (std::size_t r = 0; r < num_points; ++r)
                        node_ids.emplace(static_cast<std::int64_t>(r) + 1, r);
                    node_ids_dense = false;
                }
                // A duplicate is unrepresentable -- two coordinate rows would
                // claim one id -- so it throws, always. The dense path cannot
                // produce one by construction (`id == row + 1` strictly
                // increases), so the check only runs where it can fire.
                if (!node_ids_dense && !node_ids.emplace(id, num_points).second)
                    throw ReadError("MDPA: duplicate node id " + std::to_string(id));
                raw_node_ids.push_back(id);
                for (std::size_t c = t.size() - 3; c < t.size(); ++c) {
                    double v = 0.0;
                    if (!mdpa_parse_double(t[c], v))
                        throw ReadError("MDPA: non-numeric node coordinate: " + std::string(e));
                    coords.push_back(v);
                }
                ++num_points;
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before 'End Nodes'");
        } else if (mdpa_starts_with(line, "Begin Elements") ||
                   mdpa_starts_with(line, "Begin Conditions")) {
            const bool is_condition = mdpa_starts_with(line, "Begin Conditions");
            const std::string end_token = is_condition ? "End Conditions" : "End Elements";
            const std::vector<std::string_view> head = mdpa_tokens(line);
            const std::string entity_name = head.size() >= 3 ? std::string(head[2]) : std::string();
            const CellType type = mdpa_entity_cell_type(entity_name);
            if (type == CellType::Custom)
                throw ReadError("MDPA: unknown Kratos entity name '" + entity_name + "'");
            const int nn = cell_type_num_nodes(type);
            if (nn <= 0)
                throw ReadError("MDPA: entity '" + entity_name +
                                "' maps to a variable-node-count cell type");
            const std::string& type_name = cell_type_name(type);
            const std::vector<int>& order = mdpa_kratos_node_order(type);

            bool terminated = false;
            std::vector<std::string_view> t;  // reused: one allocation per block, not per row
            while (!cur.Done()) {
                const std::string_view e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == end_token) {
                    terminated = true;
                    break;
                }
                if (mdpa_starts_with(e, "End "))
                    throw ReadError("MDPA: expected '" + end_token + "', got '" + std::string(e) +
                                    "'");
                detail::split_blanks(e, t);
                if (static_cast<int>(t.size()) != nn + 2)
                    throw ReadError("MDPA: " + entity_name + " row with " +
                                    std::to_string(t.size() >= 2 ? t.size() - 2 : 0) +
                                    " nodes (expected " + std::to_string(nn) +
                                    "): " + std::string(e));
                std::int64_t id = 0, prop = 0;
                if (!mdpa_parse_int(t[0], id) || !mdpa_parse_int(t[1], prop))
                    throw ReadError("MDPA: non-integer id/property in: " + std::string(e));

                // The Kratos *name* is part of the split key, not just the cell
                // type: two adjacent SmallDisplacementElement3D4N and
                // TotalLagrangianElement3D4N blocks are both `tetra`, and
                // merging them would leave one name to write both back under --
                // exactly the silent degradation MdpaInfo exists to stop.
                if (blocks.empty() || blocks.back().mType != type_name ||
                    blocks.back().mIsCondition != is_condition ||
                    blocks.back().mEntityName != entity_name) {
                    MdpaBlock b;
                    b.mType = type_name;
                    b.mEntityName = entity_name;
                    b.mNodes = static_cast<std::size_t>(nn);
                    b.mIsCondition = is_condition;
                    blocks.push_back(std::move(b));
                }
                MdpaBlock& blk = blocks.back();
                const std::size_t base = blk.mConn.size();
                blk.mConn.resize(base + static_cast<std::size_t>(nn));
                for (int j = 0; j < nn; ++j) {
                    std::int64_t node = 0;
                    if (!mdpa_parse_int(t[static_cast<std::size_t>(j) + 2], node))
                        throw ReadError("MDPA: non-integer node id in: " + std::string(e));
                    const std::size_t slot =
                        order.empty()
                            ? static_cast<std::size_t>(j)
                            : static_cast<std::size_t>(order[static_cast<std::size_t>(j)]);
                    // The raw file id; resolved to a point row in the
                    // materialize pass below (see MdpaBlock::mConn).
                    blk.mConn[base + slot] = node;
                }
                blk.mProps.push_back(prop);
                blk.mFileIds.push_back(id);
                std::int64_t& next_id = is_condition ? next_condition_id : next_element_id;
                if (id != next_id)
                    entities_dense = false;
                ++next_id;
                const std::size_t row = blk.mCount++;
                auto& id_map = is_condition ? condition_ids : element_ids;
                id_map[id] = {blocks.size() - 1, row};
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before '" + end_token + "'");
        } else if (mdpa_starts_with(line, "Begin Properties")) {
            // Parsed unconditionally, not gated on mLenient: this is a pure
            // de-throwing, so no read that used to succeed changes.
            PropertySet ps = mdpa_parse_properties(cur, line);
            // Onto the mesh, so a registry-based consumer gets it. Through
            // v9.1.0 this rode the MdpaInfo side channel only, which nothing
            // reachable from registry_readers() could ask for -- so the values
            // were unreachable from every consumer that did not link
            // formats/mdpa.hpp and call read_mdpa directly.
            property_sets.push_back(ps);
            if (pInfo)
                pInfo->mProperties.push_back(std::move(ps));
        } else if (mdpa_starts_with(line, "Begin NodalData")) {
            const std::vector<std::string_view> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed NodalData header: " + std::string(line));
            if (num_points == 0)
                throw ReadError("MDPA: NodalData before Nodes");
            std::string name(head[2]);
            const std::size_t br = name.find('[');
            if (br != std::string::npos)
                name = name.substr(0, br);
            StagedData sd;
            sd.mName = name;
            const int nc = mdpa_parse_data_block(
                cur, "End NodalData", /*nodal=*/true,
                [&](std::int64_t id, std::size_t& block, std::size_t& row) {
                    block = 0;
                    return node_row(id, row);
                },
                sd.mRows, sd.mHasFixed);
            if (nc < 0) {
                log::warn("mdpa: NodalData block '{}' held no usable line", name);
                continue;
            }
            sd.mComponents = nc;
            point_data.push_back(std::move(sd));
        } else if (mdpa_starts_with(line, "Begin ElementalData") ||
                   mdpa_starts_with(line, "Begin ConditionalData")) {
            const bool elemental = mdpa_starts_with(line, "Begin ElementalData");
            const std::string end_token = elemental ? "End ElementalData" : "End ConditionalData";
            const std::vector<std::string_view> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed " + end_token.substr(4) +
                                " header: " + std::string(line));
            std::string name(head[2]);
            const std::size_t br = name.find('[');
            if (br != std::string::npos)
                name = name.substr(0, br);
            if (name == "gmsh:physical")
                name += "_data";  // never shadow the property-id array
            const auto& id_map = elemental ? element_ids : condition_ids;
            StagedData sd;
            sd.mName = name;
            const int nc = mdpa_parse_data_block(
                cur, end_token, /*nodal=*/false,
                [&](std::int64_t id, std::size_t& block, std::size_t& row) {
                    auto it = id_map.find(id);
                    if (it == id_map.end())
                        return false;
                    block = it->second.first;
                    row = it->second.second;
                    return true;
                },
                sd.mRows, sd.mHasFixed);
            if (nc < 0) {
                log::warn("mdpa: {} block '{}' held no usable line", end_token.substr(4), name);
                continue;
            }
            sd.mComponents = nc;
            cell_data.push_back(std::move(sd));
        } else if (mdpa_starts_with(line, "Begin SubModelPartData")) {
            if (pInfo && !smp_stack.empty()) {
                std::vector<PropertyValue> data = mdpa_parse_kv_block(cur, "End SubModelPartData");
                if (!data.empty()) {
                    MdpaSubModelPart& r_smp = smp_info(smp_name());
                    for (PropertyValue& r_v : data)
                        r_smp.mData.push_back(std::move(r_v));
                }
            } else {
                mdpa_expect_empty_block(cur, "End SubModelPartData",
                                        "a non-empty SubModelPartData block", Lenient, pInfo);
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPartTables")) {
            if (pInfo && !smp_stack.empty()) {
                std::vector<std::int64_t> ids;
                read_id_list("End SubModelPartTables", ids);
                if (!ids.empty()) {
                    MdpaSubModelPart& r_smp = smp_info(smp_name());
                    r_smp.mTables.insert(r_smp.mTables.end(), ids.begin(), ids.end());
                }
            } else {
                mdpa_expect_empty_block(cur, "End SubModelPartTables",
                                        "a non-empty SubModelPartTables block", Lenient, pInfo);
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPartGeometries") ||
                   mdpa_starts_with(line, "Begin SubModelPartConstraints")) {
            const bool geometry = mdpa_starts_with(line, "Begin SubModelPartGeometries");
            const std::string tag = geometry ? "SubModelPartGeometries" : "SubModelPartConstraints";
            if (smp_stack.empty())
                throw ReadError("MDPA: " + tag + " outside a SubModelPart");
            std::vector<std::int64_t> ids;
            read_id_list("End " + tag, ids, /*Strict=*/true);
            if (!ids.empty()) {
                if (pInfo) {
                    MdpaSubModelPart& r_smp = smp_info(smp_name());
                    auto& out = geometry ? r_smp.mGeometryIds : r_smp.mConstraintIds;
                    out.insert(out.end(), ids.begin(), ids.end());
                } else if (!Lenient) {
                    throw ReadError("MDPA: a non-empty " + tag + " block" + kMdpaNeedsInfo);
                } else {
                    log::warn("mdpa: skipping a non-empty {} block (ReadOptions::mLenient)", tag);
                }
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPartNodes")) {
            if (smp_stack.empty())
                throw ReadError("MDPA: SubModelPartNodes outside a SubModelPart");
            std::vector<std::int64_t> ids;
            read_id_list("End SubModelPartNodes", ids);
            auto& smp = smps[smp_name()];
            for (std::int64_t id : ids) {
                std::size_t row = 0;
                if (!node_row(id, row)) {
                    log::warn("mdpa: SubModelPart references unknown node id {}", id);
                    continue;
                }
                smp.mNodes.push_back(static_cast<std::int64_t>(row));
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPartElements") ||
                   mdpa_starts_with(line, "Begin SubModelPartConditions")) {
            const bool elemental = mdpa_starts_with(line, "Begin SubModelPartElements");
            if (smp_stack.empty())
                throw ReadError("MDPA: SubModelPart entity list outside a SubModelPart");
            std::vector<std::int64_t> ids;
            read_id_list(elemental ? "End SubModelPartElements" : "End SubModelPartConditions",
                         ids);
            const auto& id_map = elemental ? element_ids : condition_ids;
            auto& smp = smps[smp_name()];
            for (std::int64_t id : ids) {
                auto it = id_map.find(id);
                if (it == id_map.end()) {
                    log::warn("mdpa: SubModelPart references unknown entity id {}", id);
                    continue;
                }
                smp.mCells.push_back(it->second);
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPart")) {
            const std::vector<std::string_view> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed SubModelPart header: " + std::string(line));
            smp_stack.emplace_back(head[2]);
            smps[smp_name()];  // an entity-less SubModelPart is still a group
        } else if (line == "End SubModelPart") {
            if (smp_stack.empty())
                throw ReadError("MDPA: 'End SubModelPart' without a matching 'Begin'");
            smp_stack.pop_back();
        } else if (mdpa_starts_with(line, "Begin Table")) {
            if (pInfo)
                pInfo->mTables.push_back(mdpa_parse_property_table(cur, line));
            else
                mdpa_reject_or_skip(cur, "End Table", "a top-level Table block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin Geometries")) {
            if (!pInfo) {
                mdpa_reject_or_skip(cur, "End Geometries", "a Geometries block", Lenient, pInfo);
                continue;
            }
            const std::vector<std::string_view> head = mdpa_tokens(line);
            const std::string name = head.size() >= 3 ? std::string(head[2]) : std::string();
            const CellType type = mdpa_entity_cell_type(name);
            if (type == CellType::Custom)
                throw ReadError("MDPA: unknown Kratos geometry name '" + name + "'");
            const int nn = cell_type_num_nodes(type);
            if (nn <= 0)
                throw ReadError("MDPA: geometry '" + name +
                                "' maps to a variable-node-count cell type");
            const std::vector<int>& order = mdpa_kratos_node_order(type);
            if (geometries.empty() || geometries.back().mName != name) {
                StagedGeometry g;
                g.mName = name;
                g.mType = type;
                g.mNodes = static_cast<std::size_t>(nn);
                geometries.push_back(std::move(g));
            }
            StagedGeometry& r_geo = geometries.back();
            bool terminated = false;
            std::vector<std::string_view> t;
            while (!cur.Done()) {
                const std::string_view e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End Geometries") {
                    terminated = true;
                    break;
                }
                detail::split_blanks(e, t);
                if (static_cast<int>(t.size()) != nn + 1)
                    throw ReadError("MDPA: " + name + " row with " +
                                    std::to_string(t.empty() ? 0 : t.size() - 1) +
                                    " nodes (expected " + std::to_string(nn) +
                                    "): " + std::string(e));
                std::int64_t id = 0;
                if (!mdpa_parse_int(t[0], id))
                    throw ReadError("MDPA: non-integer geometry id in: " + std::string(e));
                const std::size_t base = r_geo.mConn.size();
                r_geo.mConn.resize(base + static_cast<std::size_t>(nn));
                for (int j = 0; j < nn; ++j) {
                    std::int64_t node = 0;
                    if (!mdpa_parse_int(t[static_cast<std::size_t>(j) + 1], node))
                        throw ReadError("MDPA: non-integer node id in: " + std::string(e));
                    const std::size_t slot =
                        order.empty()
                            ? static_cast<std::size_t>(j)
                            : static_cast<std::size_t>(order[static_cast<std::size_t>(j)]);
                    r_geo.mConn[base + slot] = node;
                }
                r_geo.mIds.push_back(id);
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before 'End Geometries'");
        } else if (mdpa_starts_with(line, "Begin Mesh")) {
            if (!pInfo) {
                mdpa_reject_or_skip(cur, "End Mesh", "a Mesh block", Lenient, pInfo);
                continue;
            }
            // Kratos reserves mesh 0 for the model part itself; the Python
            // reference warns and skips such a header, and so does this.
            const std::vector<std::string_view> head = mdpa_tokens(line);
            std::int64_t mesh_id = 0;
            if (head.size() < 3 || !mdpa_parse_int(head[2], mesh_id) || mesh_id == 0) {
                log::warn("mdpa: skipping Mesh block with a missing, non-integer or 0 id: {}",
                          line);
                mdpa_consume_block(cur, "End Mesh");
                continue;
            }
            MdpaMeshBlock mb;
            mb.mId = mesh_id;
            std::vector<std::int64_t> nodes;
            bool terminated = false;
            while (!cur.Done()) {
                const std::string_view e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End Mesh") {
                    terminated = true;
                    break;
                }
                if (mdpa_starts_with(e, "Begin MeshData")) {
                    std::vector<PropertyValue> data = mdpa_parse_kv_block(cur, "End MeshData");
                    for (PropertyValue& r_v : data)
                        mb.mData.push_back(std::move(r_v));
                } else if (mdpa_starts_with(e, "Begin MeshNodes")) {
                    read_id_list("End MeshNodes", nodes);
                } else if (mdpa_starts_with(e, "Begin MeshElements")) {
                    read_id_list("End MeshElements", mb.mElementIds);
                } else if (mdpa_starts_with(e, "Begin MeshConditions")) {
                    read_id_list("End MeshConditions", mb.mConditionIds);
                } else {
                    log::warn("mdpa: skipping unknown line in Mesh {}: {}", mesh_id, e);
                }
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before 'End Mesh'");
            pInfo->mMeshBlocks.push_back(std::move(mb));
            mesh_block_nodes.push_back(std::move(nodes));
        } else if (mdpa_starts_with(line, "Begin ")) {
            // Everything unrecognized, which is how `Begin Constraints` and any
            // block a future Kratos adds are covered without a case each. The
            // terminator is the header's first word after `Begin`, so a nested
            // `End <other>` cannot end the scan early. A top-level one is kept
            // verbatim when there is an MdpaInfo to hold it.
            const std::vector<std::string_view> head = mdpa_tokens(line);
            const std::string end_token =
                "End " + (head.size() >= 2 ? std::string(head[1]) : std::string());
            if (pInfo && smp_stack.empty())
                pInfo->mRawBlocks.push_back(mdpa_capture_raw_block(cur, line, end_token));
            else
                mdpa_reject_or_skip(cur, end_token, "the block '" + std::string(line) + "'",
                                    Lenient, pInfo);
        } else {
            throw ReadError("MDPA: unexpected line outside a block: '" + std::string(line) + "'");
        }
    }
    if (!smp_stack.empty())
        throw ReadError("MDPA: EOF before 'End SubModelPart'");

    // ---- materialize ------------------------------------------------------
    Mesh mesh;
    {
        NDArray pts(DType::Float64, {num_points, 3});
        double* pp = pts.As<double>();
        for (std::size_t i = 0; i < coords.size(); ++i)
            pp[i] = coords[i];
        mesh.AssignPoints(std::move(pts));
    }
    // Attach original node ids ONLY when they weren't already the trivial
    // `1..n` renumbering the writer would produce anyway -- so a sequential (or
    // id-less) deck's `Mesh` is untouched by this feature and a re-write is
    // byte-identical to before. `write_mdpa` looks for this exact name.
    if (!node_ids_dense) {
        NDArray ids(DType::Int64, {num_points});
        std::int64_t* ip = ids.As<std::int64_t>();
        for (std::size_t i = 0; i < num_points; ++i)
            ip[i] = raw_node_ids[i];
        mesh.AddPointData(kMdpaIdName, std::move(ids));
    }

    std::vector<NDArray> props;
    std::vector<std::size_t> block_base(blocks.size(), 0);
    std::size_t running = 0;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        MdpaBlock& blk = blocks[b];
        block_base[b] = running;
        running += blk.mCount;
        NDArray conn(DType::Int64, {blk.mCount, blk.mNodes});
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < blk.mConn.size(); ++i) {
            // Resolve the raw file ids staged above. The message names the FILE
            // id, never a row: with arbitrary ids "the file has N nodes" is no
            // longer the criterion (id 500 can be valid in a 4-node deck), so it
            // is reported as context instead.
            std::size_t row = 0;
            if (!node_row(blk.mConn[i], row))
                throw ReadError("MDPA: connectivity refers to node id " +
                                std::to_string(blk.mConn[i]) +
                                ", which the file's Nodes block does not define (" +
                                std::to_string(num_points) + " nodes read)");
            cp[i] = static_cast<std::int64_t>(row);
        }
        mesh.AddCellBlock(blk.mType, std::move(conn));
        if (pInfo)
            pInfo->mEntityNames.push_back(MdpaEntityName{blk.mEntityName, blk.mIsCondition});
        NDArray tag(DType::Int64, {blk.mCount});
        std::int64_t* tp = tag.As<std::int64_t>();
        for (std::size_t i = 0; i < blk.mProps.size(); ++i)
            tp[i] = blk.mProps[i];
        props.push_back(std::move(tag));
    }
    if (!blocks.empty())
        mesh.AddCellData("gmsh:physical", std::move(props));
    // Same "only when it matters" rule as the node ids above: elements and
    // conditions each have their own independent 1-based file-order counter,
    // and only when EITHER disagreed with a trivial renumbering is the
    // original id worth carrying -- so a fresh write of an untouched deck
    // stays byte-identical, and a gapped/reclassified one round-trips.
    if (!blocks.empty() && !entities_dense) {
        std::vector<NDArray> ids;
        ids.reserve(blocks.size());
        for (const MdpaBlock& blk : blocks) {
            NDArray a(DType::Int64, {blk.mCount});
            std::int64_t* ap = a.As<std::int64_t>();
            for (std::size_t i = 0; i < blk.mFileIds.size(); ++i)
                ap[i] = blk.mFileIds[i];
            ids.push_back(std::move(a));
        }
        mesh.AddCellData(kMdpaIdName, std::move(ids));
    }

    for (auto& fd : field_data)
        mesh.AddFieldData(fd.first, std::move(fd.second));

    for (const auto& sd : point_data) {
        mesh.AddPointData(sd.mName, mdpa_data_array(sd.mRows, 0, num_points, sd.mComponents));
        if (sd.mHasFixed) {
            NDArray fx(DType::Int64, {num_points});
            std::int64_t* fp = fx.As<std::int64_t>();
            for (std::size_t i = 0; i < num_points; ++i)
                fp[i] = -1;
            for (const auto& r : sd.mRows)
                if (r.mFixed >= 0 && r.mRow < num_points)
                    fp[r.mRow] = r.mFixed;
            mesh.AddPointData(sd.mName + "_fixed_status", std::move(fx));
        }
    }
    for (const auto& sd : cell_data) {
        std::vector<NDArray> arrays;
        arrays.reserve(blocks.size());
        for (std::size_t b = 0; b < blocks.size(); ++b)
            arrays.push_back(mdpa_data_array(sd.mRows, b, blocks[b].mCount, sd.mComponents));
        mesh.AddCellData(sd.mName, std::move(arrays));
    }

    for (const auto& smp : smps) {
        if (smp.second.mNodes.empty() && smp.second.mCells.empty()) {
            // An entity-less SubModelPart is still a named group: carry it as
            // an empty Point region rather than losing the name.
            mesh.AddRegion(Region(smp.first, RegionKind::Point, NDArray(DType::Int64, {0})));
            continue;
        }
        if (!smp.second.mNodes.empty()) {
            NDArray e(DType::Int64, {smp.second.mNodes.size()});
            std::int64_t* p = e.As<std::int64_t>();
            for (std::size_t i = 0; i < smp.second.mNodes.size(); ++i)
                p[i] = smp.second.mNodes[i];
            mesh.AddRegion(Region(smp.first, RegionKind::Point, std::move(e)));
        }
        if (!smp.second.mCells.empty()) {
            NDArray e(DType::Int64, {smp.second.mCells.size()});
            std::int64_t* p = e.As<std::int64_t>();
            for (std::size_t i = 0; i < smp.second.mCells.size(); ++i)
                p[i] = static_cast<std::int64_t>(block_base[smp.second.mCells[i].first] +
                                                 smp.second.mCells[i].second);
            mesh.AddRegion(Region(smp.first, RegionKind::Cell, std::move(e)));
        }
    }
    if (pInfo) {
        for (const StagedGeometry& r_geo : geometries) {
            MdpaGeometryBlock gb;
            gb.mName = r_geo.mName;
            gb.mType = cell_type_name(r_geo.mType);
            gb.mConn = NDArray(DType::Int64, {r_geo.mIds.size(), r_geo.mNodes});
            std::int64_t* gp = gb.mConn.As<std::int64_t>();
            for (std::size_t i = 0; i < r_geo.mConn.size(); ++i) {
                std::size_t row = 0;
                if (!node_row(r_geo.mConn[i], row))
                    throw ReadError("MDPA: geometry refers to node id " +
                                    std::to_string(r_geo.mConn[i]) +
                                    ", which the file's Nodes block does not define");
                gp[i] = static_cast<std::int64_t>(row);
            }
            gb.mIds = r_geo.mIds;
            pInfo->mGeometries.push_back(std::move(gb));
        }
        // Mesh-block nodes resolve like SubModelPart nodes: an id the file
        // never defined is warned about and dropped (the Python reference's
        // `_node_rows`).
        for (std::size_t k = 0; k < mesh_block_nodes.size(); ++k) {
            std::vector<std::int64_t>& r_rows = pInfo->mMeshBlocks[k].mNodes;
            for (std::int64_t id : mesh_block_nodes[k]) {
                std::size_t row = 0;
                if (!node_row(id, row)) {
                    log::warn("mdpa: Mesh {} references unknown node id {}",
                              pInfo->mMeshBlocks[k].mId, id);
                    continue;
                }
                r_rows.push_back(static_cast<std::int64_t>(row));
            }
        }
    }

    // Properties last: they are keyed by id, so order relative to the cell
    // blocks and regions above does not matter.
    for (PropertySet& r_ps : property_sets)
        mesh.AddPropertySet(std::move(r_ps));

    return mesh;
}

}  // namespace

Mesh read_mdpa(const std::string& rPath) {
    return mdpa_read_impl(rPath, /*Lenient=*/false, /*pInfo=*/nullptr);
}

Mesh read_mdpa(const std::string& rPath, const ReadOptions& rOptions) {
    return mdpa_read_impl(rPath, rOptions.mLenient, /*pInfo=*/nullptr);
}

Mesh read_mdpa(const std::string& rPath, MdpaInfo& rInfo, const ReadOptions& rOptions) {
    rInfo = MdpaInfo{};
    return mdpa_read_impl(rPath, rOptions.mLenient, &rInfo);
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

namespace {

bool mdpa_is_int_dtype(DType dt) {
    return dt == DType::Int8 || dt == DType::Int16 || dt == DType::Int32 || dt == DType::Int64 ||
           dt == DType::UInt8 || dt == DType::UInt16 || dt == DType::UInt32 || dt == DType::UInt64;
}

/// One value of @p rArray, formatted the way MDPA spells numbers.
std::string mdpa_format_value(const NDArray& rArray, std::size_t index) {
    char buf[64];
    if (mdpa_is_int_dtype(rArray.Dtype())) {
        std::snprintf(buf, sizeof(buf), "%lld",
                      static_cast<long long>(detail::read_int(rArray, index)));
    } else {
        detail::snprintf_c(buf, sizeof(buf), "%.16g", detail::read_double(rArray, index));
    }
    return buf;
}

/// A data array formatted as `mdpa_format_value` formats it, with the dtype
/// switch taken once for the array instead of once per value.
class MdpaValues {
public:
    explicit MdpaValues(const NDArray& rArray) {
        if (mdpa_is_int_dtype(rArray.Dtype()))
            mInts.emplace(rArray);
        else
            mDoubles.emplace(rArray);
    }

    /// `std::isnan(detail::read_double(rArray, index))`: an integer is never NaN.
    bool IsNan(std::size_t index) const { return mDoubles && std::isnan((*mDoubles)[index]); }

    std::string Format(std::size_t index) const {
        char buf[64];
        if (mInts)
            std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>((*mInts)[index]));
        else
            detail::snprintf_c(buf, sizeof(buf), "%.16g", (*mDoubles)[index]);
        return buf;
    }

private:
    std::optional<detail::Int64View> mInts;
    std::optional<detail::DoubleView> mDoubles;
};

/// Number of trailing components of a data array (1 for a 1-D array).
std::size_t mdpa_components(const NDArray& rArray, std::size_t rows) {
    if (rows == 0)
        return 1;
    return rArray.Size() / rows;
}

bool mdpa_skip_point_data(const std::string& rName) {
    if (rName == kMdpaIdName)
        return true;
    const std::string suffix = "_fixed_status";
    if (rName.size() >= suffix.size() &&
        rName.compare(rName.size() - suffix.size(), suffix.size(), suffix) == 0)
        return true;
    return mdpa_starts_with(rName, "gmsh:");
}

bool mdpa_skip_cell_data(const std::string& rName) {
    if (rName == kMdpaIdName)
        return true;
    const std::string suffix = "_tag";
    if (rName.size() >= suffix.size() &&
        rName.compare(rName.size() - suffix.size(), suffix.size(), suffix) == 0)
        return true;
    return mdpa_starts_with(rName, "gmsh:");
}

/**
 * @brief Emit one `Begin Table <mKey>` block at @p rIndent.
 *
 * `mKey` holds the header's arguments verbatim (id + variable names), so the
 * block comes back out exactly as it went in.
 */
void mdpa_write_table(std::ostream& rOs, const PropertyValue& rTable, const char* pIndent) {
    rOs << pIndent << "Begin Table " << rTable.mKey << "\n";
    const std::size_t ncols = rTable.mValues.Shape().size() >= 2 ? rTable.mValues.Shape()[1] : 1;
    const std::size_t nrows = ncols ? rTable.mValues.Size() / ncols : 0;
    const MdpaValues values(rTable.mValues);
    for (std::size_t r = 0; r < nrows; ++r) {
        rOs << pIndent << "  ";
        for (std::size_t c = 0; c < ncols; ++c)
            rOs << " " << values.Format(r * ncols + c);
        rOs << "\n";
    }
    rOs << pIndent << "End Table\n";
}

/// Emit one `KEY value` line of a Properties/ModelPartData/*Data body.
void mdpa_write_kv(std::ostream& rOs, const PropertyValue& rValue, const char* pIndent) {
    rOs << pIndent << rValue.mKey << " ";
    if (rValue.IsText()) {
        rOs << rValue.mText;
    } else {
        const MdpaValues values(rValue.mValues);
        for (std::size_t i = 0; i < rValue.mValues.Size(); ++i) {
            if (i)
                rOs << " ";
            rOs << values.Format(i);
        }
    }
    rOs << "\n";
}

/// Emit one `Begin Properties <id>` block, bodies included.
void mdpa_write_properties(std::ostream& rOs, const PropertySet& rSet) {
    rOs << "Begin Properties " << rSet.mId << "\n";
    for (const PropertyValue& v : rSet.mValues) {
        if (v.mIsTable)
            mdpa_write_table(rOs, v, "  ");
        else
            mdpa_write_kv(rOs, v, "  ");
    }
    rOs << "End Properties\n\n";
}

/// Emit an id list sub-block (`SubModelPartNodes`, `MeshElements`, ...), if non-empty.
void mdpa_write_id_list(std::ostream& rOs, const char* pTag, const std::vector<std::int64_t>& rIds,
                        const std::string& rIndent = "    ") {
    if (rIds.empty())
        return;
    rOs << rIndent << "Begin " << pTag << "\n";
    for (std::int64_t id : rIds)
        rOs << rIndent << "    " << id << "\n";
    rOs << rIndent << "End " << pTag << "\n";
}

}  // namespace

void write_mdpa(const std::string& rPath, const Mesh& rMesh) {
    write_mdpa(rPath, rMesh, MdpaInfo{});
}

void write_mdpa(const std::string& rPath, const Mesh& rMesh, const MdpaInfo& rInfo) {
    // No provenance slot in this format: drop the notes this write raises on
    // the way out rather than let them reach the next file written.
    const detail::ProvenanceSlotlessWrite slotless;
    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    // ---- per-block decisions (entity kind, Kratos name, written ids) -------
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<bool> is_condition(nblocks, false);
    std::vector<std::string> entity_name(nblocks);
    std::vector<std::vector<std::int64_t>> written_ids(nblocks);
    std::vector<std::size_t> block_base(nblocks, 0);
    {
        // Honour cell_data["mdpa:id"] (kMdpaIdName) when the mesh carries one
        // array per block AND every array's row count matches its block's cell
        // count -- anything short of that is treated as unrelated/stale
        // metadata and falls back to the old renumbering, exactly like the
        // node-id check below. An id that survives is still validated for
        // uniqueness (elements and conditions each have their own Kratos
        // namespace, so a collision is only checked within its own kind) --
        // writing a duplicate would silently produce an invalid Kratos deck.
        bool preserve_entity_ids =
            rMesh.HasCellData(kMdpaIdName) && rMesh.CellDataNumBlocks(kMdpaIdName) == nblocks;
        for (std::size_t b = 0; preserve_entity_ids && b < nblocks; ++b)
            if (rMesh.CellData(kMdpaIdName, b).Size() != rMesh.Cells(b).NumCells())
                preserve_entity_ids = false;

        std::int64_t next_element = 1, next_condition = 1;
        std::unordered_set<std::int64_t> seen_element_ids, seen_condition_ids;
        std::size_t running = 0;
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            if (cb.IsRagged())
                throw WriteError("MDPA: ragged/polyhedron cell blocks are not supported (block " +
                                 std::to_string(b) + ", type '" + std::string(cb.Type()) + "')");
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            if (type == CellType::Custom)
                throw WriteError("MDPA: no Kratos entity name for cell type '" +
                                 std::string(cb.Type()) + "'");
            // The Python reference's rule for a mesh with no physical tags: a
            // block whose default Kratos *element* name is a 2-D one is written
            // as a Condition, everything else as an Element.
            is_condition[b] = kratos_element_name(type).find("2D") != std::string::npos;
            entity_name[b] =
                is_condition[b] ? kratos_condition_name(type) : kratos_element_name(type);
            // A name the reader kept wins over the derived one, which is what
            // makes an application-specific SmallDisplacementElement3D4N
            // survive a round trip instead of collapsing to Element3D4N. The
            // recorded kind comes with it: inferring "Condition" from the name
            // would be a second guess on top of the cell-type heuristic.
            if (b < rInfo.mEntityNames.size() && !rInfo.mEntityNames[b].mName.empty()) {
                entity_name[b] = rInfo.mEntityNames[b].mName;
                is_condition[b] = rInfo.mEntityNames[b].mIsCondition;
            }
            block_base[b] = running;
            running += cb.NumCells();
            written_ids[b].resize(cb.NumCells());
            std::unordered_set<std::int64_t>& seen =
                is_condition[b] ? seen_condition_ids : seen_element_ids;
            std::optional<detail::Int64View> ids;
            if (preserve_entity_ids)
                ids.emplace(rMesh.CellData(kMdpaIdName, b));
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t wid =
                    ids ? (*ids)[r] : (is_condition[b] ? next_condition++ : next_element++);
                if (!seen.insert(wid).second)
                    throw WriteError("MDPA: duplicate " +
                                     std::string(is_condition[b] ? "condition" : "element") +
                                     " id " + std::to_string(wid) + " in cell_data['" +
                                     std::string(kMdpaIdName) + "']");
                written_ids[b][r] = wid;
            }
        }
    }

    // ---- ModelPartData ----------------------------------------------------
    os << "Begin ModelPartData\n";
    for (const auto& name : rMesh.FieldDataNames()) {
        const NDArray& a = rMesh.FieldData(name);
        if (a.Size() != 1) {
            log::warn("mdpa: field_data '{}' has {} values; only scalars are written", name,
                      a.Size());
            detail::provenance_note("data-dropped", "field_data '" + name +
                                                        "' not written -- MDPA's ModelPartData "
                                                        "holds scalars only");
            continue;
        }
        os << "    " << name << " " << mdpa_format_value(a, 0) << "\n";
    }
    for (const PropertyValue& r_v : rInfo.mModelPartData) {
        if (rMesh.HasFieldData(r_v.mKey)) {
            log::warn(
                "mdpa: ModelPartData '{}' is both field_data and MdpaInfo text; the "
                "field_data value is written",
                r_v.mKey);
            continue;
        }
        mdpa_write_kv(os, r_v, "    ");
    }
    os << "End ModelPartData\n\n";

    // ---- Properties -------------------------------------------------------
    // With an MdpaInfo the blocks come back with their bodies. Without one,
    // every id the entity rows below will actually reference gets an empty
    // block, ascending: the rows have always written their `gmsh:physical`
    // value as the property id, so hard-coding a single `Properties 0` left a
    // tagged mesh referencing undeclared properties, which Kratos's own
    // ModelPartIO rejects. A mesh whose ids are all 0 -- every mesh with no
    // `gmsh:physical` -- still emits exactly the old two lines.
    const bool has_props =
        rMesh.HasCellData("gmsh:physical") && rMesh.CellDataNumBlocks("gmsh:physical") == nblocks;
    std::set<std::int64_t> referenced_ids;
    if (has_props) {
        for (std::size_t b = 0; b < nblocks; ++b) {
            const NDArray& tags = rMesh.CellData("gmsh:physical", b);
            const detail::Int64View tag_values(tags);
            for (std::size_t r = 0; r < tags.Size(); ++r)
                referenced_ids.insert(tag_values[r]);
        }
    }
    if (!rInfo.mProperties.empty()) {
        // An explicit MdpaInfo wins, and keeps the caller's order verbatim --
        // which is the one thing the mesh channel cannot do, since it
        // canonicalizes to ascending id.
        for (const PropertySet& ps : rInfo.mProperties)
            mdpa_write_properties(os, ps);
    } else if (rMesh.NumPropertySets() > 0) {
        // Bodies carried on the mesh (v9.2.0): what read_mdpa now stores, so a
        // registry-driven mdpa -> mdpa round trip keeps its material data
        // instead of emitting empty blocks.
        for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i) {
            mdpa_write_properties(os, rMesh.GetPropertySet(i));
            referenced_ids.erase(rMesh.GetPropertySet(i).mId);
        }
        // Any id the rows reference but no set covers still has to be declared:
        // Kratos's own ModelPartIO rejects a row naming an undeclared id.
        for (std::int64_t id : referenced_ids)
            os << "Begin Properties " << id << "\nEnd Properties\n\n";
    } else {
        if (referenced_ids.empty())
            referenced_ids.insert(0);
        for (std::int64_t id : referenced_ids)
            os << "Begin Properties " << id << "\nEnd Properties\n\n";
    }

    // ---- Tables -----------------------------------------------------------
    for (const PropertyValue& r_table : rInfo.mTables) {
        mdpa_write_table(os, r_table, "");
        os << "\n";
    }

    // ---- Nodes ------------------------------------------------------------
    os << "Begin Nodes\n";
    std::vector<std::int64_t> written_node_ids;
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = rMesh.PointDim();
        const std::size_t np = rMesh.NumPoints();
        // Honour point_data["mdpa:id"] (kMdpaIdName) when present and the right
        // length; anything short of that (absent, wrong size, wrong dtype) is
        // treated as unrelated metadata and falls back to the old row+1
        // numbering. Values are validated for uniqueness -- a duplicate would
        // silently produce an ambiguous file.
        const bool preserve_node_ids =
            rMesh.HasPointData(kMdpaIdName) && rMesh.PointData(kMdpaIdName).Size() == np;
        written_node_ids.resize(np);
        std::unordered_set<std::int64_t> seen_node_ids;
        const detail::DoubleView point_values(points);
        std::optional<detail::Int64View> node_ids;
        if (preserve_node_ids)
            node_ids.emplace(rMesh.PointData(kMdpaIdName));
        char buf[64];
        for (std::size_t i = 0; i < np; ++i) {
            const std::int64_t id = node_ids ? (*node_ids)[i] : static_cast<std::int64_t>(i) + 1;
            if (!seen_node_ids.insert(id).second)
                throw WriteError("MDPA: duplicate node id " + std::to_string(id) +
                                 " in point_data['" + std::string(kMdpaIdName) + "']");
            written_node_ids[i] = id;
            os << " " << id;
            for (std::size_t c = 0; c < 3; ++c) {
                const double v = c < dim ? point_values[i * dim + c] : 0.0;
                detail::snprintf_c(buf, sizeof(buf), "%.16e", v);
                os << " " << buf;
            }
            os << "\n";
        }
    }
    os << "End Nodes\n\n";

    // ---- Elements / Conditions -------------------------------------------
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const std::vector<int>& order = mdpa_kratos_node_order(type);
        const std::string kind = is_condition[b] ? "Conditions" : "Elements";
        os << "Begin " << kind << " " << entity_name[b] << "\n";
        const detail::Int64View conn(cb.Conn());
        const std::size_t k = cb.NodesPerCell();
        std::optional<detail::Int64View> props;
        if (has_props)
            props.emplace(rMesh.CellData("gmsh:physical", b));
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::int64_t prop = 0;
            if (props)
                prop = (*props)[r];
            os << "  " << written_ids[b][r] << " " << prop;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t slot = order.empty() ? j : static_cast<std::size_t>(order[j]);
                // Through `written_node_ids`, not a bare `+ 1`: connectivity
                // must name whichever node numbering was actually written
                // (preserved or row+1), the same rule the Nodes block itself
                // and the SubModelPart node lists follow.
                const std::size_t row = static_cast<std::size_t>(conn[r * k + slot]);
                os << " " << written_node_ids[row];
            }
            os << "\n";
        }
        os << "End " << kind << "\n\n";
    }

    // ---- Geometries -------------------------------------------------------
    const std::size_t np = rMesh.NumPoints();
    std::int64_t next_geometry = 1;
    for (const MdpaGeometryBlock& r_geo : rInfo.mGeometries) {
        const CellType type = cell_type_from_name(r_geo.mType);
        if (type == CellType::Custom)
            throw WriteError("MDPA: geometry block of unknown cell type '" + r_geo.mType + "'");
        const std::size_t k = r_geo.mConn.Shape().size() == 2 ? r_geo.mConn.Shape()[1] : 0;
        if (static_cast<int>(k) != cell_type_num_nodes(type))
            throw WriteError("MDPA: geometry block '" + r_geo.mType + "' has " + std::to_string(k) +
                             " nodes per row");
        const std::size_t rows = r_geo.mConn.Size() / k;
        const bool keep_ids = r_geo.mIds.size() == rows;
        if (!keep_ids && !r_geo.mIds.empty())
            log::warn("mdpa: geometry block '{}' has {} ids for {} rows; renumbering", r_geo.mType,
                      r_geo.mIds.size(), rows);
        const std::vector<int>& order = mdpa_kratos_node_order(type);
        const detail::Int64View geo_conn(r_geo.mConn);
        os << "Begin Geometries "
           << (r_geo.mName.empty() ? kratos_geometry_name(type) : r_geo.mName) << "\n";
        for (std::size_t r = 0; r < rows; ++r) {
            os << "  " << (keep_ids ? r_geo.mIds[r] : next_geometry++);
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t slot = order.empty() ? j : static_cast<std::size_t>(order[j]);
                const std::int64_t row = geo_conn[r * k + slot];
                if (row < 0 || static_cast<std::size_t>(row) >= np)
                    throw WriteError("MDPA: geometry row names point " + std::to_string(row) +
                                     " of a mesh with " + std::to_string(np) + " points");
                os << " " << written_node_ids[static_cast<std::size_t>(row)];
            }
            os << "\n";
        }
        os << "End Geometries\n\n";
    }

    // ---- NodalData --------------------------------------------------------
    for (const auto& name : rMesh.PointDataNames()) {
        if (mdpa_skip_point_data(name))
            continue;
        const NDArray& a = rMesh.PointData(name);
        const std::size_t nc = mdpa_components(a, np);
        const std::string fixed_name = name + "_fixed_status";
        const bool has_fixed = rMesh.HasPointData(fixed_name);
        const MdpaValues values(a);
        std::optional<detail::Int64View> fixed;
        if (has_fixed)
            fixed.emplace(rMesh.PointData(fixed_name));
        os << "Begin NodalData " << name << "\n";
        for (std::size_t i = 0; i < np; ++i) {
            bool all_nan = true;
            for (std::size_t j = 0; j < nc; ++j)
                if (!values.IsNan(i * nc + j))
                    all_nan = false;
            if (all_nan)
                continue;
            os << "  " << written_node_ids[i];
            if (fixed) {
                const std::int64_t f = (*fixed)[i];
                if (f >= 0)
                    os << " " << f;
            }
            for (std::size_t j = 0; j < nc; ++j)
                os << " " << values.Format(i * nc + j);
            os << "\n";
        }
        os << "End NodalData\n\n";
    }

    // ---- ElementalData / ConditionalData ----------------------------------
    for (const auto& name : rMesh.CellDataNames()) {
        if (mdpa_skip_cell_data(name))
            continue;
        if (rMesh.CellDataNumBlocks(name) != nblocks)
            continue;
        for (int pass = 0; pass < 2; ++pass) {
            const bool conditions = pass == 1;
            auto body = detail::make_classic_ostringstream();
            for (std::size_t b = 0; b < nblocks; ++b) {
                if (is_condition[b] != conditions)
                    continue;
                const auto cb = rMesh.Cells(b);
                const NDArray& a = rMesh.CellData(name, b);
                const std::size_t nc = mdpa_components(a, cb.NumCells());
                const MdpaValues values(a);
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    bool all_nan = true;
                    for (std::size_t j = 0; j < nc; ++j)
                        if (!values.IsNan(r * nc + j))
                            all_nan = false;
                    if (all_nan)
                        continue;
                    body << "  " << written_ids[b][r];
                    for (std::size_t j = 0; j < nc; ++j)
                        body << " " << values.Format(r * nc + j);
                    body << "\n";
                }
            }
            if (body.str().empty())
                continue;
            const std::string kind = conditions ? "ConditionalData" : "ElementalData";
            os << "Begin " << kind << " " << name << "\n" << body.str() << "End " << kind << "\n\n";
        }
    }

    // ---- SubModelParts from named regions ---------------------------------
    // Side-channel-only parts are retained. Hierarchical region keys are
    // emitted as nested blocks, not a single Kratos name containing '/'.
    std::unordered_map<std::string, const MdpaSubModelPart*> smp_extras;
    for (const MdpaSubModelPart& r_smp : rInfo.mSubModelParts)
        smp_extras.emplace(r_smp.mName, &r_smp);
    auto write_smp_extras = [&](const std::string& rName, const std::string& rIndent) {
        const auto it = smp_extras.find(rName);
        if (it == smp_extras.end())
            return;
        if (!it->second->mData.empty()) {
            os << rIndent << "Begin SubModelPartData\n";
            for (const PropertyValue& r_v : it->second->mData)
                mdpa_write_kv(os, r_v, (rIndent + "    ").c_str());
            os << rIndent << "End SubModelPartData\n";
        }
        mdpa_write_id_list(os, "SubModelPartTables", it->second->mTables, rIndent);
        mdpa_write_id_list(os, "SubModelPartGeometries", it->second->mGeometryIds, rIndent);
        mdpa_write_id_list(os, "SubModelPartConstraints", it->second->mConstraintIds, rIndent);
    };
    std::unordered_map<std::string, std::vector<std::string>> smp_children;
    std::unordered_set<std::string> smp_seen;
    auto add_smp = [&](const std::string& rName) {
        if (rName.empty() || rName.front() == '/' || rName.back() == '/' ||
            rName.find("//") != std::string::npos)
            throw WriteError("MDPA: invalid SubModelPart hierarchy name '" + rName + "'");
        std::string parent;
        std::size_t start = 0;
        while (start < rName.size()) {
            const std::size_t slash = rName.find('/', start);
            const std::string path = rName.substr(0, slash);
            if (smp_seen.insert(path).second)
                smp_children[parent].push_back(path);
            parent = path;
            if (slash == std::string::npos)
                break;
            start = slash + 1;
        }
    };
    for (const auto& name : rMesh.RegionNames()) {
        bool has_membership = false;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const auto& region = rMesh.Region(i);
            if (region.mName != name)
                continue;
            if (region.mKind == RegionKind::Side) {
                log::warn("mdpa: dropping side region '{}' (MDPA has no facet sets)", name);
                detail::provenance_note(
                    "regions-dropped",
                    "side region '" + name + "' dropped -- MDPA has no facet sets");
            } else {
                has_membership = true;
            }
        }
        if (has_membership)
            add_smp(name);
    }
    for (const MdpaSubModelPart& r_smp : rInfo.mSubModelParts)
        add_smp(r_smp.mName);
    auto write_smp = [&](const std::string& name, std::size_t depth) {
        std::vector<std::int64_t> nodes;
        std::vector<std::int64_t> elements, conditions;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = rMesh.Region(i);
            if (r.mName != name)
                continue;
            if (r.mKind == RegionKind::Side) {
                continue;
            }
            const std::int64_t* e = r.Entries();
            for (std::size_t j = 0; j < r.NumEntries(); ++j) {
                if (r.mKind == RegionKind::Point) {
                    // `written_node_ids` already reflects whichever numbering
                    // was actually written (preserved or row+1), so this stays
                    // consistent with the Nodes block above with no extra work.
                    nodes.push_back(written_node_ids[static_cast<std::size_t>(e[j])]);
                    continue;
                }
                // Cell region: global block-major index -> (block, row) -> id.
                const std::size_t g = static_cast<std::size_t>(e[j]);
                for (std::size_t b = 0; b < nblocks; ++b) {
                    const std::size_t count = rMesh.Cells(b).NumCells();
                    if (g < block_base[b] || g >= block_base[b] + count)
                        continue;
                    const std::int64_t id = written_ids[b][g - block_base[b]];
                    (is_condition[b] ? conditions : elements).push_back(id);
                    break;
                }
            }
        }
        const std::string indent(depth * 4, ' ');
        const std::string body_indent = indent + "    ";
        const std::size_t slash = name.rfind('/');
        const std::string leaf = slash == std::string::npos ? name : name.substr(slash + 1);
        os << indent << "Begin SubModelPart " << leaf << "\n";
        write_smp_extras(name, body_indent);
        mdpa_write_id_list(os, "SubModelPartNodes", nodes, body_indent);
        mdpa_write_id_list(os, "SubModelPartElements", elements, body_indent);
        mdpa_write_id_list(os, "SubModelPartConditions", conditions, body_indent);
    };
    // Explicit traversal stack: a valid deeply nested deck must not consume
    // one C++ call frame per part.
    struct MdpaSmpWriteFrame {
        std::string mName;
        std::size_t mDepth;
        bool mClose;
    };
    std::vector<MdpaSmpWriteFrame> smp_frames;
    const auto& roots = smp_children[""];
    for (auto it = roots.rbegin(); it != roots.rend(); ++it)
        smp_frames.push_back({*it, 0, false});
    while (!smp_frames.empty()) {
        const MdpaSmpWriteFrame frame = std::move(smp_frames.back());
        smp_frames.pop_back();
        if (frame.mClose) {
            os << std::string(frame.mDepth * 4, ' ') << "End SubModelPart\n"
               << (frame.mDepth == 0 ? "\n" : "");
            continue;
        }
        write_smp(frame.mName, frame.mDepth);
        smp_frames.push_back({frame.mName, frame.mDepth, true});
        const auto children = smp_children.find(frame.mName);
        if (children != smp_children.end())
            for (auto it = children->second.rbegin(); it != children->second.rend(); ++it)
                smp_frames.push_back({*it, frame.mDepth + 1, false});
    }

    // ---- Mesh blocks ------------------------------------------------------
    // Nodes go through the written numbering; element and condition members
    // are the file's own ids, written verbatim (the Python reference's rule).
    for (const MdpaMeshBlock& r_mb : rInfo.mMeshBlocks) {
        os << "Begin Mesh " << r_mb.mId << "\n";
        if (!r_mb.mData.empty()) {
            os << "    Begin MeshData\n";
            for (const PropertyValue& r_v : r_mb.mData)
                mdpa_write_kv(os, r_v, "        ");
            os << "    End MeshData\n";
        }
        std::vector<std::int64_t> node_ids;
        node_ids.reserve(r_mb.mNodes.size());
        for (std::int64_t row : r_mb.mNodes) {
            if (row < 0 || static_cast<std::size_t>(row) >= np)
                throw WriteError("MDPA: Mesh " + std::to_string(r_mb.mId) + " names point " +
                                 std::to_string(row) + " of a mesh with " + std::to_string(np) +
                                 " points");
            node_ids.push_back(written_node_ids[static_cast<std::size_t>(row)]);
        }
        mdpa_write_id_list(os, "MeshNodes", node_ids);
        mdpa_write_id_list(os, "MeshElements", r_mb.mElementIds);
        mdpa_write_id_list(os, "MeshConditions", r_mb.mConditionIds);
        os << "End Mesh\n\n";
    }

    // ---- Raw blocks, verbatim ---------------------------------------------
    for (const MdpaRawBlock& r_raw : rInfo.mRawBlocks)
        os << r_raw.mHeader << "\n" << r_raw.mBody << r_raw.mEnd << "\n\n";
}

}  // namespace meshioplusplus
