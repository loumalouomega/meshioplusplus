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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Lexing helpers. Every anonymous-namespace symbol here is `mdpa_`-prefixed
// because the single-header amalgamation concatenates all of src/cpp/src.
// ---------------------------------------------------------------------------

std::string mdpa_strip(const std::string& rS) {
    const std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    const std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}

/// The line with any `//` comment removed, then stripped.
std::string mdpa_clean(const std::string& rS) {
    const std::size_t c = rS.find("//");
    return mdpa_strip(c == std::string::npos ? rS : rS.substr(0, c));
}

std::vector<std::string> mdpa_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

bool mdpa_starts_with(const std::string& rS, const std::string& rPrefix) {
    return rS.size() >= rPrefix.size() && rS.compare(0, rPrefix.size(), rPrefix) == 0;
}

bool mdpa_parse_int(const std::string& rS, std::int64_t& rOut) {
    if (rS.empty())
        return false;
    char* end = nullptr;
    const long long v = std::strtoll(rS.c_str(), &end, 10);
    if (end != rS.c_str() + rS.size())
        return false;
    rOut = static_cast<std::int64_t>(v);
    return true;
}

bool mdpa_parse_double(const std::string& rS, double& rOut) {
    if (rS.empty())
        return false;
    char* end = nullptr;
    const double v = std::strtod(rS.c_str(), &end);
    if (end != rS.c_str() + rS.size())
        return false;
    rOut = v;
    return true;
}

// ---------------------------------------------------------------------------
// Kratos <-> VTK node orderings.
//
// `mdpa_kratos_node_order(type)[i]` is the meshio/VTK slot the node in Kratos
// slot `i` belongs in. Reading therefore scatters (`out[table[i]] = in[i]`) and
// writing gathers (`out[j] = in[table[j]]`) — the exact inverse pair the Python
// reference builds with `np.argsort`.
// ---------------------------------------------------------------------------

const std::vector<int>& mdpa_kratos_node_order(CellType type) {
    static const std::vector<int> h20 = {0,  1, 2,  3,  4,  5,  6,  7,  8,  11,
                                         10, 9, 16, 19, 18, 17, 12, 13, 14, 15};
    static const std::vector<int> h27 = {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 10, 9,  16, 19,
                                         18, 17, 12, 15, 14, 13, 20, 23, 21, 24, 22, 25, 26};
    static const std::vector<int> none;
    if (type == CellType::Hexahedron20)
        return h20;
    if (type == CellType::Hexahedron27)
        return h27;
    return none;
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
    return cell_type_from_kratos_name_or_suffix(rName);
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
struct MdpaCursor {
    const std::vector<std::string>* mpLines = nullptr;
    std::size_t mIndex = 0;

    bool Done() const { return mIndex >= mpLines->size(); }
    const std::string& Next() { return (*mpLines)[mIndex++]; }
};

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
        throw ReadError("MDPA: " + rWhat +
                        " is not supported by the C++ reader (set ReadOptions::mLenient to skip "
                        "it instead)");
    log::warn("mdpa: skipping {} (ReadOptions::mLenient)", rWhat);
    if (pInfo)
        pInfo->mSkippedConstructs.push_back(rWhat);
    mdpa_consume_block(rCur, rEnd);
}

/// Consume a block that must be empty; a non-empty one goes through the above.
void mdpa_expect_empty_block(MdpaCursor& rCur, const std::string& rEnd, const std::string& rWhat,
                             bool Lenient, MdpaInfo* pInfo) {
    while (!rCur.Done()) {
        const std::string line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == rEnd)
            return;
        if (!Lenient)
            throw ReadError("MDPA: " + rWhat +
                            " is not supported by the C++ reader (offending line: '" + line +
                            "'; set ReadOptions::mLenient to skip it instead)");
        log::warn("mdpa: skipping {} (ReadOptions::mLenient)", rWhat);
        if (pInfo)
            pInfo->mSkippedConstructs.push_back(rWhat);
        mdpa_consume_block(rCur, rEnd);
        return;
    }
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief Parse an inline `Begin Table <args>` inside a `Properties` body.
 *
 * @p rHeader is the whole header line; everything after `Begin Table` becomes
 * the entry's `mKey` verbatim, so the block re-emits with its id and variable
 * names unchanged. Rows are whitespace-separated numbers; the first usable row
 * fixes the column count and a row that disagrees is warned about and skipped.
 */
PropertyValue mdpa_parse_property_table(MdpaCursor& rCur, const std::string& rHeader) {
    PropertyValue out;
    out.mIsTable = true;
    out.mKey = mdpa_strip(rHeader.substr(std::string("Begin Table").size()));

    std::vector<double> values;
    std::size_t ncols = 0;
    bool terminated = false;
    while (!rCur.Done()) {
        const std::string line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == "End Table") {
            terminated = true;
            break;
        }
        const std::vector<std::string> toks = mdpa_tokens(line);
        std::vector<double> row;
        row.reserve(toks.size());
        bool ok = true;
        for (const std::string& tok : toks) {
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
PropertySet mdpa_parse_properties(MdpaCursor& rCur, const std::string& rHeader) {
    PropertySet out;
    const std::vector<std::string> head = mdpa_tokens(rHeader);
    if (head.size() < 3 || !mdpa_parse_int(head[2], out.mId)) {
        log::warn("mdpa: Properties block with no readable id, using 0: {}", rHeader);
        out.mId = 0;
    }

    while (true) {
        if (rCur.Done())
            throw ReadError("MDPA: EOF before 'End Properties'");
        const std::string line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == "End Properties")
            break;
        if (mdpa_starts_with(line, "Begin Table")) {
            out.mValues.push_back(mdpa_parse_property_table(rCur, line));
            continue;
        }
        // key + the rest of the line, exactly the Python reference's
        // `split(None, 1)`, so both readers agree on where the value starts.
        const std::size_t sep = line.find_first_of(" \t");
        if (sep == std::string::npos) {
            log::warn("mdpa: skipping valueless Properties line: {}", line);
            continue;
        }
        PropertyValue v;
        v.mKey = line.substr(0, sep);
        const std::string rest = mdpa_strip(line.substr(sep + 1));
        double scalar = 0.0;
        if (mdpa_parse_double(rest, scalar)) {
            NDArray a(DType::Float64, {1});
            a.As<double>()[0] = scalar;
            v.mValues = std::move(a);
        } else {
            v.mText = rest;
        }
        out.mValues.push_back(std::move(v));
    }
    return out;
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
    while (!rCur.Done()) {
        const std::string raw = mdpa_strip(rCur.Next());
        if (raw == rEnd) {
            terminated = true;
            break;
        }
        const std::string line = mdpa_clean(raw);
        if (line.empty())
            continue;
        if (line == rEnd) {
            terminated = true;
            break;
        }
        const std::vector<std::string> toks = mdpa_tokens(line);
        std::int64_t id = 0;
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
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l))
        lines.push_back(l);

    MdpaCursor cur{&lines, 0};

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

    auto smp_name = [&]() -> std::string {
        std::string out;
        for (std::size_t i = 0; i < smp_stack.size(); ++i) {
            if (i)
                out += "/";
            out += smp_stack[i];
        }
        return out;
    };

    auto read_id_list = [&](const std::string& rEnd, std::vector<std::int64_t>& rOut) {
        while (!cur.Done()) {
            const std::string line = mdpa_clean(cur.Next());
            if (line.empty())
                continue;
            if (line == rEnd)
                return;
            std::int64_t v = 0;
            if (!mdpa_parse_int(line, v)) {
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
        const std::string raw = mdpa_strip(cur.Next());
        const std::string line = mdpa_clean(raw);
        if (line.empty())
            continue;

        if (line == "Begin ModelPartData") {
            while (true) {
                if (cur.Done())
                    throw ReadError("MDPA: EOF before 'End ModelPartData'");
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End ModelPartData")
                    break;
                const std::vector<std::string> t = mdpa_tokens(e);
                if (t.size() < 2) {
                    log::warn("mdpa: skipping malformed ModelPartData line: {}", e);
                    continue;
                }
                double v = 0.0;
                if (t.size() != 2 || !mdpa_parse_double(t[1], v)) {
                    const std::string what = "a non-numeric ModelPartData value for '" + t[0] + "'";
                    if (!Lenient)
                        throw ReadError("MDPA: " + what +
                                        " is not supported by the C++ reader (set "
                                        "ReadOptions::mLenient to skip it instead)");
                    log::warn("mdpa: skipping {} (ReadOptions::mLenient)", what);
                    if (pInfo)
                        pInfo->mSkippedConstructs.push_back(what);
                    continue;
                }
                NDArray a(DType::Float64, {1});
                a.As<double>()[0] = v;
                field_data[t[0]] = std::move(a);
            }
        } else if (line == "Begin Nodes") {
            if (num_points)
                throw ReadError("MDPA: more than one Nodes block");
            bool terminated = false;
            while (!cur.Done()) {
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End Nodes") {
                    terminated = true;
                    break;
                }
                const std::vector<std::string> t = mdpa_tokens(e);
                if (t.size() < 3)
                    throw ReadError("MDPA: node line with fewer than 3 coordinates: " + e);
                // An id-less row (`x y z`) takes its position as its id, which is
                // exactly what "connectivity is 1-based into row order" already
                // meant for a fully id-less file -- so such a file never leaves
                // the dense path and reads byte-identically to before.
                std::int64_t id = static_cast<std::int64_t>(num_points) + 1;
                if (t.size() >= 4 && !mdpa_parse_int(t[0], id))
                    throw ReadError("MDPA: non-integer node id: " + e);
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
                        throw ReadError("MDPA: non-numeric node coordinate: " + e);
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
            const std::vector<std::string> head = mdpa_tokens(line);
            const std::string entity_name = head.size() >= 3 ? head[2] : std::string();
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
            while (!cur.Done()) {
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == end_token) {
                    terminated = true;
                    break;
                }
                if (mdpa_starts_with(e, "End "))
                    throw ReadError("MDPA: expected '" + end_token + "', got '" + e + "'");
                const std::vector<std::string> t = mdpa_tokens(e);
                if (static_cast<int>(t.size()) != nn + 2)
                    throw ReadError("MDPA: " + entity_name + " row with " +
                                    std::to_string(t.size() >= 2 ? t.size() - 2 : 0) +
                                    " nodes (expected " + std::to_string(nn) + "): " + e);
                std::int64_t id = 0, prop = 0;
                if (!mdpa_parse_int(t[0], id) || !mdpa_parse_int(t[1], prop))
                    throw ReadError("MDPA: non-integer id/property in: " + e);

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
                        throw ReadError("MDPA: non-integer node id in: " + e);
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
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed NodalData header: " + line);
            if (num_points == 0)
                throw ReadError("MDPA: NodalData before Nodes");
            std::string name = head[2];
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
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed " + end_token.substr(4) + " header: " + line);
            std::string name = head[2];
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
            mdpa_expect_empty_block(cur, "End SubModelPartData",
                                    "a non-empty SubModelPartData block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin SubModelPartTables")) {
            mdpa_expect_empty_block(cur, "End SubModelPartTables",
                                    "a non-empty SubModelPartTables block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin SubModelPartGeometries")) {
            mdpa_expect_empty_block(cur, "End SubModelPartGeometries",
                                    "a non-empty SubModelPartGeometries block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin SubModelPartConstraints")) {
            mdpa_expect_empty_block(cur, "End SubModelPartConstraints",
                                    "a non-empty SubModelPartConstraints block", Lenient, pInfo);
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
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed SubModelPart header: " + line);
            smp_stack.push_back(head[2]);
            smps[smp_name()];  // an entity-less SubModelPart is still a group
        } else if (line == "End SubModelPart") {
            if (smp_stack.empty())
                throw ReadError("MDPA: 'End SubModelPart' without a matching 'Begin'");
            smp_stack.pop_back();
        } else if (mdpa_starts_with(line, "Begin Table")) {
            mdpa_reject_or_skip(cur, "End Table", "a top-level Table block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin Geometries")) {
            mdpa_reject_or_skip(cur, "End Geometries", "a Geometries block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin Mesh")) {
            mdpa_reject_or_skip(cur, "End Mesh", "a Mesh block", Lenient, pInfo);
        } else if (mdpa_starts_with(line, "Begin ")) {
            // Everything unrecognized, which is how `Begin Constraints` and any
            // block a future Kratos adds are covered without a case each. The
            // terminator is the header's first word after `Begin`, so a nested
            // `End <other>` cannot end the scan early.
            const std::vector<std::string> head = mdpa_tokens(line);
            const std::string end_token = "End " + (head.size() >= 2 ? head[1] : std::string());
            mdpa_reject_or_skip(cur, end_token, "the block '" + line + "'", Lenient, pInfo);
        } else {
            throw ReadError("MDPA: unexpected line outside a block: '" + line + "'");
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
        std::snprintf(buf, sizeof(buf), "%.16g", detail::read_double(rArray, index));
    }
    return buf;
}

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

/// Emit one `Begin Properties <id>` block, bodies included.
void mdpa_write_properties(std::ostream& rOs, const PropertySet& rSet) {
    rOs << "Begin Properties " << rSet.mId << "\n";
    for (const PropertyValue& v : rSet.mValues) {
        if (v.mIsTable) {
            // mKey holds the header's arguments verbatim (id + variable names),
            // so the block comes back out exactly as it went in.
            rOs << "  Begin Table " << v.mKey << "\n";
            const std::size_t ncols = v.mValues.Shape().size() >= 2 ? v.mValues.Shape()[1] : 1;
            const std::size_t nrows = ncols ? v.mValues.Size() / ncols : 0;
            for (std::size_t r = 0; r < nrows; ++r) {
                rOs << "   ";
                for (std::size_t c = 0; c < ncols; ++c)
                    rOs << " " << mdpa_format_value(v.mValues, r * ncols + c);
                rOs << "\n";
            }
            rOs << "  End Table\n";
            continue;
        }
        rOs << "  " << v.mKey << " ";
        if (v.IsText()) {
            rOs << v.mText;
        } else {
            for (std::size_t i = 0; i < v.mValues.Size(); ++i) {
                if (i)
                    rOs << " ";
                rOs << mdpa_format_value(v.mValues, i);
            }
        }
        rOs << "\n";
    }
    rOs << "End Properties\n\n";
}

}  // namespace

void write_mdpa(const std::string& rPath, const Mesh& rMesh) {
    write_mdpa(rPath, rMesh, MdpaInfo{});
}

void write_mdpa(const std::string& rPath, const Mesh& rMesh, const MdpaInfo& rInfo) {
    std::ofstream os(rPath);
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
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t wid =
                    preserve_entity_ids ? detail::read_int(rMesh.CellData(kMdpaIdName, b), r)
                                        : (is_condition[b] ? next_condition++ : next_element++);
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
            for (std::size_t r = 0; r < tags.Size(); ++r)
                referenced_ids.insert(detail::read_int(tags, r));
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
        char buf[64];
        for (std::size_t i = 0; i < np; ++i) {
            const std::int64_t id = preserve_node_ids
                                        ? detail::read_int(rMesh.PointData(kMdpaIdName), i)
                                        : static_cast<std::int64_t>(i) + 1;
            if (!seen_node_ids.insert(id).second)
                throw WriteError("MDPA: duplicate node id " + std::to_string(id) +
                                 " in point_data['" + std::string(kMdpaIdName) + "']");
            written_node_ids[i] = id;
            os << " " << id;
            for (std::size_t c = 0; c < 3; ++c) {
                const double v = c < dim ? detail::read_double(points, i * dim + c) : 0.0;
                std::snprintf(buf, sizeof(buf), "%.16e", v);
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
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::int64_t prop = 0;
            if (has_props)
                prop = detail::read_int(rMesh.CellData("gmsh:physical", b), r);
            os << "  " << written_ids[b][r] << " " << prop;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t slot = order.empty() ? j : static_cast<std::size_t>(order[j]);
                // Through `written_node_ids`, not a bare `+ 1`: connectivity
                // must name whichever node numbering was actually written
                // (preserved or row+1), the same rule the Nodes block itself
                // and the SubModelPart node lists follow.
                const std::size_t row =
                    static_cast<std::size_t>(detail::read_int(conn, r * k + slot));
                os << " " << written_node_ids[row];
            }
            os << "\n";
        }
        os << "End " << kind << "\n\n";
    }

    // ---- NodalData --------------------------------------------------------
    const std::size_t np = rMesh.NumPoints();
    for (const auto& name : rMesh.PointDataNames()) {
        if (mdpa_skip_point_data(name))
            continue;
        const NDArray& a = rMesh.PointData(name);
        const std::size_t nc = mdpa_components(a, np);
        const std::string fixed_name = name + "_fixed_status";
        const bool has_fixed = rMesh.HasPointData(fixed_name);
        os << "Begin NodalData " << name << "\n";
        for (std::size_t i = 0; i < np; ++i) {
            bool all_nan = true;
            for (std::size_t j = 0; j < nc; ++j)
                if (!std::isnan(detail::read_double(a, i * nc + j)))
                    all_nan = false;
            if (all_nan)
                continue;
            os << "  " << written_node_ids[i];
            if (has_fixed) {
                const std::int64_t f = detail::read_int(rMesh.PointData(fixed_name), i);
                if (f >= 0)
                    os << " " << f;
            }
            for (std::size_t j = 0; j < nc; ++j)
                os << " " << mdpa_format_value(a, i * nc + j);
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
            std::ostringstream body;
            for (std::size_t b = 0; b < nblocks; ++b) {
                if (is_condition[b] != conditions)
                    continue;
                const auto cb = rMesh.Cells(b);
                const NDArray& a = rMesh.CellData(name, b);
                const std::size_t nc = mdpa_components(a, cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    bool all_nan = true;
                    for (std::size_t j = 0; j < nc; ++j)
                        if (!std::isnan(detail::read_double(a, r * nc + j)))
                            all_nan = false;
                    if (all_nan)
                        continue;
                    body << "  " << written_ids[b][r];
                    for (std::size_t j = 0; j < nc; ++j)
                        body << " " << mdpa_format_value(a, r * nc + j);
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
    for (const auto& name : rMesh.RegionNames()) {
        std::vector<std::int64_t> nodes;
        std::vector<std::int64_t> elements, conditions;
        bool any = false;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = rMesh.Region(i);
            if (r.mName != name)
                continue;
            if (r.mKind == RegionKind::Side) {
                log::warn("mdpa: dropping side region '{}' (MDPA has no facet sets)", name);
                detail::provenance_note("regions-dropped", "side region '" + name +
                                                               "' dropped -- MDPA has no facet "
                                                               "sets");
                continue;
            }
            any = true;
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
        if (!any)
            continue;
        os << "Begin SubModelPart " << name << "\n";
        auto emit = [&](const char* pTag, const std::vector<std::int64_t>& rIds) {
            if (rIds.empty())
                return;
            os << "    Begin SubModelPart" << pTag << "\n";
            for (std::int64_t id : rIds)
                os << "        " << id << "\n";
            os << "    End SubModelPart" << pTag << "\n";
        };
        emit("Nodes", nodes);
        emit("Elements", elements);
        emit("Conditions", conditions);
        os << "End SubModelPart\n\n";
    }
}

}  // namespace meshioplusplus
