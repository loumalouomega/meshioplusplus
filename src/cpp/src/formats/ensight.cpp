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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Shared tables
// ---------------------------------------------------------------------------

// EnSight Gold element keyword <-> meshio cell type (fixed-node types only;
// nsided/nfaced are handled separately as ragged blocks).
struct EnsightTypeEntry {
    const char* mKeyword;
    const char* mMeshioType;
    int mNumNodes;
};

const std::vector<EnsightTypeEntry>& ensight_type_table() {
    static const std::vector<EnsightTypeEntry> table = {
        {"point", "vertex", 1},         {"bar2", "line", 2},
        {"bar3", "line3", 3},           {"tria3", "triangle", 3},
        {"tria6", "triangle6", 6},      {"quad4", "quad", 4},
        {"quad8", "quad8", 8},          {"tetra4", "tetra", 4},
        {"tetra10", "tetra10", 10},     {"pyramid5", "pyramid", 5},
        {"pyramid13", "pyramid13", 13}, {"penta6", "wedge", 6},
        {"penta15", "wedge15", 15},     {"hexa8", "hexahedron", 8},
        {"hexa20", "hexahedron20", 20}};
    return table;
}

const EnsightTypeEntry* ensight_entry_from_keyword(const std::string& rKeyword) {
    static const std::unordered_map<std::string, const EnsightTypeEntry*> map = [] {
        std::unordered_map<std::string, const EnsightTypeEntry*> m;
        for (const auto& e : ensight_type_table())
            m.emplace(e.mKeyword, &e);
        return m;
    }();
    auto it = map.find(rKeyword);
    return it == map.end() ? nullptr : it->second;
}

const EnsightTypeEntry* ensight_entry_from_meshio(const std::string& rType) {
    static const std::unordered_map<std::string, const EnsightTypeEntry*> map = [] {
        std::unordered_map<std::string, const EnsightTypeEntry*> m;
        for (const auto& e : ensight_type_table())
            m.emplace(e.mMeshioType, &e);
        return m;
    }();
    auto it = map.find(rType);
    return it == map.end() ? nullptr : it->second;
}

// meshio <-> EnSight node-order permutation. meshio ordering is VTK ordering
// for every supported type; the only difference is the prism triangle
// winding of penta15 vs wedge15 — the same involution VTK's EnSight readers
// apply, so one table serves both directions (result[j] = source index for
// position j). Kratos's penta15/hexa20 group swaps fix Kratos-specific
// ordering and must NOT be ported here.
const std::vector<int>* ensight_permutation(const std::string& rMeshioType) {
    static const std::vector<int> wedge15 = {0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13};
    if (rMeshioType == "wedge15")
        return &wedge15;
    return nullptr;
}

// ---------------------------------------------------------------------------
// Small string / path helpers
// ---------------------------------------------------------------------------

bool ensight_has_suffix(const std::string& rPath, const std::string& rSuffix) {
    if (rPath.size() < rSuffix.size())
        return false;
    return rPath.compare(rPath.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
}

std::string ensight_trim(const std::string& rLine) {
    std::size_t b = 0, e = rLine.size();
    while (b < e && (std::isspace(static_cast<unsigned char>(rLine[b])) || rLine[b] == '\0'))
        ++b;
    while (e > b &&
           (std::isspace(static_cast<unsigned char>(rLine[e - 1])) || rLine[e - 1] == '\0'))
        --e;
    return rLine.substr(b, e - b);
}

bool ensight_starts_with(std::string_view str, const char* pPrefix) {
    return str.rfind(pPrefix, 0) == 0;
}

std::string ensight_dirname(const std::string& rPath) {
    std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? std::string() : rPath.substr(0, slash + 1);
}

std::string ensight_basename(const std::string& rPath) {
    std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? rPath : rPath.substr(slash + 1);
}

// "<stem>.case" or "<stem>.geo" -> {case path, geo path}.
std::pair<std::string, std::string> ensight_case_geo_paths(const std::string& rPath, bool& rOk) {
    rOk = true;
    if (ensight_has_suffix(rPath, ".case"))
        return {rPath, rPath.substr(0, rPath.size() - 5) + ".geo"};
    if (ensight_has_suffix(rPath, ".geo"))
        return {rPath.substr(0, rPath.size() - 4) + ".case", rPath};
    rOk = false;
    return {"", ""};
}

/**
 * @brief Whole-file access, mapped where that pays (detail/file_source.hpp).
 *
 * Returns the source itself rather than a string so the caller controls its
 * lifetime: the cursors below hold a view into it, so it must outlive them.
 */
detail::FileSource ensight_read_whole_file(const std::string& rPath, const char* pWhat) {
    try {
        return detail::FileSource(rPath);
    } catch (const ReadError&) {
        throw ReadError(std::string("EnSight: could not open ") + pWhat + ": " + rPath);
    }
}

// ---------------------------------------------------------------------------
// Geometry cursors: one record/number stream over ASCII or C-binary data
// ---------------------------------------------------------------------------

class EnsightCursor {
public:
    virtual ~EnsightCursor() = default;
    /// True once only trailing whitespace/padding remains.
    virtual bool AtEnd() = 0;
    /// Consume and return the next string record (trimmed line / 80-char record).
    virtual std::string NextRecord() = 0;
    /// Return the next string record without consuming it ("" at end).
    virtual std::string PeekRecord() = 0;
    virtual std::int64_t NextInt() = 0;
    virtual void ReadInts(std::size_t n, std::int64_t* pDst) = 0;
    virtual void ReadFloats(std::size_t n, double* pDst) = 0;
    virtual void SkipInts(std::size_t n) = 0;
    /// Binary-only hook: peek the next int32 and enable byte-swapping when it
    /// is implausible as-is but plausible swapped. With PreferSmaller (used at
    /// the tiny part-number record, where e.g. bswap(1) = 16777216 is still
    /// "plausible"), the smaller of two plausible interpretations wins.
    /// No-op for ASCII.
    virtual void CheckSwap(std::int64_t Max, bool PreferSmaller) {
        (void)Max;
        (void)PreferSmaller;
    }
};

class EnsightAsciiCursor final : public EnsightCursor {
public:
    explicit EnsightAsciiCursor(std::string_view text) : mText(text) {}

    bool AtEnd() override {
        std::size_t p = mPos;
        while (p < mText.size() &&
               (std::isspace(static_cast<unsigned char>(mText[p])) || mText[p] == '\0'))
            ++p;
        return p >= mText.size();
    }

    std::string NextRecord() override {
        while (mPos < mText.size()) {
            std::size_t eol = mText.find('\n', mPos);
            if (eol == std::string::npos)
                eol = mText.size();
            std::string line = ensight_trim(std::string(mText.substr(mPos, eol - mPos)));
            mPos = eol < mText.size() ? eol + 1 : eol;
            if (!line.empty())
                return line;
        }
        throw ReadError("EnSight: unexpected end of geometry file");
    }

    std::string PeekRecord() override {
        if (AtEnd())
            return "";
        const std::size_t saved = mPos;
        std::string rec = NextRecord();
        mPos = saved;
        return rec;
    }

    std::int64_t NextInt() override {
        const char* start = mText.data() + mPos;
        char* end = nullptr;
        const std::int64_t v = std::strtoll(start, &end, 10);
        if (end == start)
            throw ReadError("EnSight: expected an integer in geometry file");
        mPos = static_cast<std::size_t>(end - mText.data());
        return v;
    }

    void ReadInts(std::size_t n, std::int64_t* pDst) override {
        for (std::size_t i = 0; i < n; ++i)
            pDst[i] = NextInt();
    }

    void ReadFloats(std::size_t n, double* pDst) override {
        for (std::size_t i = 0; i < n; ++i) {
            const char* start = mText.data() + mPos;
            char* end = nullptr;
            pDst[i] = std::strtod(start, &end);
            if (end == start)
                throw ReadError("EnSight: expected a number in geometry file");
            mPos = static_cast<std::size_t>(end - mText.data());
        }
    }

    void SkipInts(std::size_t n) override {
        for (std::size_t i = 0; i < n; ++i)
            NextInt();
    }

private:
    std::string_view mText;
    std::size_t mPos = 0;
};

class EnsightBinaryCursor final : public EnsightCursor {
public:
    // Text is the whole file; the cursor starts after the leading
    // "C Binary" 80-char record.
    explicit EnsightBinaryCursor(std::string_view data) : mData(data), mPos(80) {}

    bool AtEnd() override { return mPos >= mData.size(); }

    std::string NextRecord() override {
        if (mPos + 80 > mData.size())
            throw ReadError("EnSight: truncated binary geometry file");
        std::string rec(mData.data() + mPos, 80);
        mPos += 80;
        const std::size_t nul = rec.find('\0');
        if (nul != std::string::npos)
            rec.resize(nul);
        return ensight_trim(rec);
    }

    std::string PeekRecord() override {
        if (mPos + 80 > mData.size())
            return "";
        const std::size_t saved = mPos;
        std::string rec = NextRecord();
        mPos = saved;
        return rec;
    }

    std::int64_t NextInt() override {
        std::int32_t v;
        Require(4);
        std::memcpy(&v, mData.data() + mPos, 4);
        mPos += 4;
        if (mSwap)
            detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
        return v;
    }

    void ReadInts(std::size_t n, std::int64_t* pDst) override {
        Require(4 * n);
        const char* src = mData.data() + mPos;
        for (std::size_t i = 0; i < n; ++i) {
            std::int32_t v;
            std::memcpy(&v, src + 4 * i, 4);
            if (mSwap)
                detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
            pDst[i] = v;
        }
        mPos += 4 * n;
    }

    void ReadFloats(std::size_t n, double* pDst) override {
        Require(4 * n);
        const char* src = mData.data() + mPos;
        for (std::size_t i = 0; i < n; ++i) {
            float v;
            std::memcpy(&v, src + 4 * i, 4);
            if (mSwap)
                detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
            pDst[i] = v;
        }
        mPos += 4 * n;
    }

    void SkipInts(std::size_t n) override {
        Require(4 * n);
        mPos += 4 * n;
    }

    void CheckSwap(std::int64_t Max, bool PreferSmaller) override {
        if (mPos + 4 > mData.size())
            return;
        std::int32_t v;
        std::memcpy(&v, mData.data() + mPos, 4);
        if (mSwap)
            detail::bswap_inplace(reinterpret_cast<char*>(&v), 4);
        std::int32_t s = v;
        detail::bswap_inplace(reinterpret_cast<char*>(&s), 4);
        const bool v_ok = v >= 0 && v <= Max;
        const bool s_ok = s >= 0 && s <= Max;
        if ((!v_ok && s_ok) || (PreferSmaller && v_ok && s_ok && s < v))
            mSwap = !mSwap;
    }

private:
    void Require(std::size_t n) {
        if (mPos + n > mData.size())
            throw ReadError("EnSight: truncated binary geometry file");
    }

    std::string_view mData;
    std::size_t mPos;
    bool mSwap = false;
};

// ---------------------------------------------------------------------------
// .case parsing
// ---------------------------------------------------------------------------

/// One `VARIABLE` section entry -- `scalar per node:`/`vector per node:`/
/// `scalar per element:`/`vector per element:` are the four kinds this
/// reader understands; anything else (complex/tensor variables, `per
/// measured node`, constants) is recorded but never read, matching the
/// roadmap's variable-reading scope.
struct EnsightVariableEntry {
    std::string mKind;
    std::string mName;
    std::string mFilePattern;  // relative to the case file's directory; may contain '*'
};

/// A `.case` file's GEOMETRY/TIME/VARIABLE sections, resolved against
/// `rCasePath`'s directory. `mTimeValues` is the *first* `time set:` found
/// (real-world Gold case files overwhelmingly have exactly one); a file with
/// several is read against that one, which is a documented narrowing, not a
/// silent wrong answer -- every variable's `[ts]` is otherwise ignored.
struct EnsightCaseInfo {
    std::string mGeoPath;
    bool mGeoIsWildcard = false;
    std::vector<double> mTimeValues;
    long mFileNameStart = 0;
    long mFileNameIncrement = 1;
    std::vector<EnsightVariableEntry> mVariables;
};

/// Splits a whitespace-separated record and drops its leading run of pure
/// integer tokens (a `[ts] [fs]` prefix) -- the same rule `model:`/variable
/// lines both use to make the leading timeset/fileset optional.
std::vector<std::string> ensight_tokens_after_leading_ints(const std::string& rValue) {
    std::istringstream toks(rValue);
    std::vector<std::string> tokens;
    std::string tok;
    while (toks >> tok)
        tokens.push_back(tok);
    std::size_t first = 0;
    while (first < tokens.size()) {
        char* end = nullptr;
        (void)std::strtoll(tokens[first].c_str(), &end, 10);
        if (end == tokens[first].c_str() || *end != '\0')
            break;
        ++first;
    }
    tokens.erase(tokens.begin(), tokens.begin() + static_cast<std::ptrdiff_t>(first));
    return tokens;
}

/// Parse the .case file: FORMAT/GEOMETRY (as before), plus TIME and
/// VARIABLE, needed for `ReadOptions::mTimeStep` and variable-file reading.
EnsightCaseInfo ensight_parse_case(const std::string& rCasePath) {
    const detail::FileSource source = ensight_read_whole_file(rCasePath, "case file");
    const std::string data(source.View());  // small text file; parsed via istringstream

    std::string section;
    std::string format_type;
    std::string model_value;
    EnsightCaseInfo info;
    bool in_time_values = false;
    bool have_time_set = false;
    long num_steps = -1;
    std::istringstream stream(data);
    std::string raw;
    while (std::getline(stream, raw)) {
        std::string line = ensight_trim(raw);
        if (line.empty() || line[0] == '#') {
            in_time_values = false;
            continue;
        }
        if (line == "FORMAT" || line == "GEOMETRY" || line == "VARIABLE" || line == "TIME" ||
            line == "FILE" || line == "MATERIAL" || line == "SCRIPTS") {
            section = line;
            in_time_values = false;
            continue;
        }
        if (section == "FORMAT" && ensight_starts_with(line, "type:")) {
            format_type = ensight_trim(line.substr(5));
        } else if (section == "GEOMETRY" && ensight_starts_with(line, "model:")) {
            model_value = ensight_trim(line.substr(6));
        } else if (section == "VARIABLE") {
            static const char* kKinds[] = {"scalar per node:", "vector per node:",
                                           "scalar per element:", "vector per element:"};
            for (const char* kind : kKinds) {
                if (!ensight_starts_with(line, kind))
                    continue;
                const std::string kind_str(kind, std::strlen(kind) - 1);  // drop trailing ':'
                const std::vector<std::string> toks =
                    ensight_tokens_after_leading_ints(line.substr(std::strlen(kind)));
                if (toks.size() < 2)
                    throw ReadError("EnSight: malformed '" + kind_str + "' line: " + line);
                EnsightVariableEntry entry;
                entry.mFilePattern = toks.back();
                std::string joined;
                for (std::size_t i = 0; i + 1 < toks.size(); ++i)
                    joined += (i ? " " : "") + toks[i];
                entry.mName = joined;
                entry.mKind = kind_str;
                info.mVariables.push_back(std::move(entry));
                break;
            }
        } else if (section == "TIME") {
            if (ensight_starts_with(line, "time set:")) {
                // A second time set: only the first is honoured (see
                // EnsightCaseInfo's own doc comment) -- clearing `section`
                // stops every TIME branch below from matching until the next
                // recognized section keyword resets it.
                if (have_time_set) {
                    section.clear();
                    continue;
                }
                have_time_set = true;
            } else if (ensight_starts_with(line, "number of steps:")) {
                num_steps =
                    std::strtol(line.c_str() + std::strlen("number of steps:"), nullptr, 10);
            } else if (ensight_starts_with(line, "filename start number:")) {
                info.mFileNameStart = std::strtol(
                    line.c_str() + std::strlen("filename start number:"), nullptr, 10);
            } else if (ensight_starts_with(line, "filename increment:")) {
                info.mFileNameIncrement =
                    std::strtol(line.c_str() + std::strlen("filename increment:"), nullptr, 10);
            } else if (ensight_starts_with(line, "time values:")) {
                in_time_values = true;
                const std::string rest = ensight_trim(line.substr(std::strlen("time values:")));
                std::istringstream iss(rest);
                double v;
                while (iss >> v)
                    info.mTimeValues.push_back(v);
            } else if (in_time_values) {
                std::istringstream iss(line);
                double v;
                while (iss >> v)
                    info.mTimeValues.push_back(v);
                if (num_steps >= 0 &&
                    info.mTimeValues.size() >= static_cast<std::size_t>(num_steps))
                    in_time_values = false;
            }
        }
    }

    if (format_type.find("ensight gold") == std::string::npos)
        throw ReadError("EnSight: case file is not 'type: ensight gold' (got '" + format_type +
                        "')");
    if (model_value.empty())
        throw ReadError("EnSight: case file has no GEOMETRY 'model:' entry");

    const std::vector<std::string> tokens = ensight_tokens_after_leading_ints(model_value);
    if (tokens.empty())
        throw ReadError("EnSight: malformed 'model:' line in case file");
    const std::string& filename = tokens.front();
    info.mGeoIsWildcard = filename.find('*') != std::string::npos;
    if (info.mGeoIsWildcard)
        throw ReadError("EnSight: transient (wildcard) geometry is not supported");
    info.mGeoPath = ensight_dirname(rCasePath) + filename;

    return info;
}

/// Replaces `filename`'s run of `*` characters with `Number`, zero-padded to
/// the run's own width -- the EnSight Gold filename-templating convention
/// TIME's `filename start number:`/`filename increment:` resolve into.
std::string ensight_resolve_wildcard(const std::string& rPattern, long Number) {
    const std::size_t star = rPattern.find('*');
    if (star == std::string::npos)
        return rPattern;
    std::size_t width = 0;
    while (star + width < rPattern.size() && rPattern[star + width] == '*')
        ++width;
    std::ostringstream num;
    num << std::setfill('0') << std::setw(static_cast<int>(width)) << Number;
    std::string digits = num.str();
    if (digits.size() > width)
        digits = digits.substr(digits.size() - width);  // Number overflowed the field width
    return rPattern.substr(0, star) + digits + rPattern.substr(star + width);
}

// ---------------------------------------------------------------------------
// Geometry reading
// ---------------------------------------------------------------------------

// A staged cell block (added to the mesh only after AssignPoints).
struct EnsightBlock {
    std::string mType;
    NDArray mConn{DType::Int64, {}};                                       // rectangular
    std::vector<std::vector<std::int64_t>> mPolygonRows;                   // nsided
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;  // nfaced
    int mKind = 0;  // 0 rectangular, 1 polygon, 2 polyhedron
    std::int64_t mPartId = 0;
    std::size_t mNumCells = 0;
};

// One part's node range (for per-node variables) and its blocks' cell ranges
// (for per-element variables, which are per-BLOCK in meshio++ but per-PART
// per-element-type in EnSight). A variable file's own "part"/id sections and
// "coordinates"/<element type> sections are structured identically to the
// geometry file's, in the same part/type order, which is what lets a
// variable file be read against this layout with no name matching at all.
struct EnsightPartLayout {
    std::int64_t mPartId = 0;
    std::size_t mPointOffset = 0, mNumPoints = 0;
    // One entry per element-type section this part had, in file order; each
    // names the Mesh cell-block index it landed in and how many cells.
    std::vector<std::pair<std::size_t, std::size_t>> mBlocks;  // (block index, num cells)
};

// "given" and "ignore" both put id arrays in the file; only the presence
// matters — Gold connectivity is positional, so ids are always skipped.
bool ensight_ids_in_file(const std::string& rRecord, const char* pWhat) {
    // rRecord is e.g. "node id assign"; the mode is the last token.
    std::istringstream iss(rRecord);
    std::string tok, mode;
    while (iss >> tok)
        mode = tok;
    if (mode == "given" || mode == "ignore")
        return true;
    if (mode == "off" || mode == "assign")
        return false;
    throw ReadError(std::string("EnSight: malformed '") + pWhat + " id' record: " + rRecord);
}

Mesh ensight_parse_geo(EnsightCursor& rCur, std::vector<EnsightPartLayout>* pLayout = nullptr) {
    constexpr std::int64_t plausible_max = 100000000;  // generous id/count bound

    rCur.NextRecord();  // description line 1
    rCur.NextRecord();  // description line 2
    std::string node_id_rec = rCur.NextRecord();
    if (!ensight_starts_with(node_id_rec, "node id"))
        throw ReadError("EnSight: expected 'node id' record, got: " + node_id_rec);
    std::string elem_id_rec = rCur.NextRecord();
    if (!ensight_starts_with(elem_id_rec, "element id"))
        throw ReadError("EnSight: expected 'element id' record, got: " + elem_id_rec);
    const bool node_ids_in_file = ensight_ids_in_file(node_id_rec, "node");
    const bool elem_ids_in_file = ensight_ids_in_file(elem_id_rec, "element");

    if (ensight_starts_with(rCur.PeekRecord(), "extents")) {
        rCur.NextRecord();
        double extents[6];
        rCur.ReadFloats(6, extents);
    }

    std::vector<double> coords;  // xyz-interleaved, all parts concatenated
    std::vector<EnsightBlock> blocks;
    std::int64_t num_parts = 0;
    // Point range per part, for a variable file's per-node sections; filled
    // regardless of pLayout (cheap) and only exposed through it.
    std::vector<std::pair<std::int64_t, std::pair<std::size_t, std::size_t>>> part_point_ranges;

    while (!rCur.AtEnd()) {
        std::string rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "part"))
            throw ReadError("EnSight: expected 'part' record, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/true);
        const std::int64_t part_id = rCur.NextInt();
        ++num_parts;
        rCur.NextRecord();  // part description

        rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "coordinates"))
            throw ReadError("EnSight: expected 'coordinates' record, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/false);
        const std::int64_t nn = rCur.NextInt();
        if (nn < 0)
            throw ReadError("EnSight: negative node count");
        const std::int64_t point_offset = static_cast<std::int64_t>(coords.size() / 3);
        part_point_ranges.emplace_back(
            part_id,
            std::make_pair(static_cast<std::size_t>(point_offset), static_cast<std::size_t>(nn)));

        if (node_ids_in_file)
            rCur.SkipInts(static_cast<std::size_t>(nn));
        std::vector<double> x(static_cast<std::size_t>(nn));
        std::vector<double> y(static_cast<std::size_t>(nn));
        std::vector<double> z(static_cast<std::size_t>(nn));
        rCur.ReadFloats(static_cast<std::size_t>(nn), x.data());
        rCur.ReadFloats(static_cast<std::size_t>(nn), y.data());
        rCur.ReadFloats(static_cast<std::size_t>(nn), z.data());
        coords.reserve(coords.size() + static_cast<std::size_t>(nn) * 3);
        for (std::int64_t i = 0; i < nn; ++i) {
            coords.push_back(x[static_cast<std::size_t>(i)]);
            coords.push_back(y[static_cast<std::size_t>(i)]);
            coords.push_back(z[static_cast<std::size_t>(i)]);
        }

        // 1-based positional index within this part -> global 0-based index.
        auto resolve = [&](std::int64_t v) -> std::int64_t {
            const std::int64_t local = v - 1;
            if (local < 0 || local >= nn)
                throw ReadError("EnSight: connectivity index out of range");
            return local + point_offset;
        };

        // Element sections until the next part / EOF.
        while (!rCur.AtEnd()) {
            std::string kw = rCur.PeekRecord();
            if (kw.empty() || ensight_starts_with(kw, "part"))
                break;
            rCur.NextRecord();
            // Ghost-cell sections carry the same data as their base type.
            if (ensight_starts_with(kw, "g_"))
                kw = kw.substr(2);

            rCur.CheckSwap(plausible_max, /*PreferSmaller=*/false);
            const std::int64_t ne = rCur.NextInt();
            if (ne < 0)
                throw ReadError("EnSight: negative element count");
            if (elem_ids_in_file)
                rCur.SkipInts(static_cast<std::size_t>(ne));

            if (kw == "nsided") {
                std::vector<std::int64_t> sizes(static_cast<std::size_t>(ne));
                rCur.ReadInts(sizes.size(), sizes.data());
                std::vector<std::vector<std::int64_t>> rows(static_cast<std::size_t>(ne));
                std::vector<std::int64_t> flat;
                std::size_t total = 0;
                for (auto s : sizes)
                    total += static_cast<std::size_t>(s);
                flat.resize(total);
                rCur.ReadInts(total, flat.data());
                std::size_t at = 0;
                for (std::size_t c = 0; c < rows.size(); ++c) {
                    rows[c].resize(static_cast<std::size_t>(sizes[c]));
                    for (std::size_t j = 0; j < rows[c].size(); ++j)
                        rows[c][j] = resolve(flat[at++]);
                }
                EnsightBlock b;
                b.mType = "polygon";
                b.mKind = 1;
                b.mPartId = part_id;
                b.mNumCells = rows.size();
                b.mPolygonRows = std::move(rows);
                blocks.push_back(std::move(b));
            } else if (kw == "nfaced") {
                std::vector<std::int64_t> nfaces(static_cast<std::size_t>(ne));
                rCur.ReadInts(nfaces.size(), nfaces.data());
                std::size_t total_faces = 0;
                for (auto f : nfaces)
                    total_faces += static_cast<std::size_t>(f);
                std::vector<std::int64_t> fsizes(total_faces);
                rCur.ReadInts(total_faces, fsizes.data());
                std::size_t total_nodes = 0;
                for (auto s : fsizes)
                    total_nodes += static_cast<std::size_t>(s);
                std::vector<std::int64_t> flat(total_nodes);
                rCur.ReadInts(total_nodes, flat.data());

                std::vector<std::vector<std::vector<std::int64_t>>> cells(
                    static_cast<std::size_t>(ne));
                std::size_t face_at = 0, node_at = 0;
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    cells[c].resize(static_cast<std::size_t>(nfaces[c]));
                    for (auto& face : cells[c]) {
                        face.resize(static_cast<std::size_t>(fsizes[face_at++]));
                        for (auto& v : face)
                            v = resolve(flat[node_at++]);
                    }
                }
                // Group by unique node count into "polyhedron<N>" blocks (the
                // openfoam convention), preserving first-seen order.
                std::vector<std::size_t> node_counts(cells.size());
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    std::vector<std::int64_t> uniq;
                    for (const auto& face : cells[c])
                        uniq.insert(uniq.end(), face.begin(), face.end());
                    std::sort(uniq.begin(), uniq.end());
                    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                    node_counts[c] = uniq.size();
                }
                std::vector<std::size_t> group_order;
                std::map<std::size_t, std::vector<std::size_t>> groups;
                for (std::size_t c = 0; c < cells.size(); ++c) {
                    if (groups.find(node_counts[c]) == groups.end())
                        group_order.push_back(node_counts[c]);
                    groups[node_counts[c]].push_back(c);
                }
                for (std::size_t n : group_order) {
                    std::vector<std::vector<std::vector<std::int64_t>>> group_cells;
                    for (std::size_t c : groups[n])
                        group_cells.push_back(std::move(cells[c]));
                    EnsightBlock b;
                    b.mType = "polyhedron" + std::to_string(n);
                    b.mKind = 2;
                    b.mPartId = part_id;
                    b.mNumCells = group_cells.size();
                    b.mPolyhedronCells = std::move(group_cells);
                    blocks.push_back(std::move(b));
                }
            } else {
                const EnsightTypeEntry* entry = ensight_entry_from_keyword(kw);
                if (entry == nullptr)
                    throw ReadError("EnSight: unsupported element keyword: " + kw);
                const std::size_t npc = static_cast<std::size_t>(entry->mNumNodes);
                NDArray conn = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ne), npc});
                std::int64_t* cp = conn.As<std::int64_t>();
                rCur.ReadInts(static_cast<std::size_t>(ne) * npc, cp);
                const std::vector<int>* perm = ensight_permutation(entry->mMeshioType);
                if (perm != nullptr) {
                    std::vector<std::int64_t> tmp(npc);
                    for (std::int64_t r = 0; r < ne; ++r) {
                        std::int64_t* row = cp + r * static_cast<std::int64_t>(npc);
                        for (std::size_t j = 0; j < npc; ++j)
                            tmp[j] = row[(*perm)[j]];
                        for (std::size_t j = 0; j < npc; ++j)
                            row[j] = tmp[j];
                    }
                }
                const std::size_t nvals = static_cast<std::size_t>(ne) * npc;
                for (std::size_t k = 0; k < nvals; ++k)
                    cp[k] = resolve(cp[k]);
                EnsightBlock b;
                b.mType = entry->mMeshioType;
                b.mKind = 0;
                b.mPartId = part_id;
                b.mNumCells = static_cast<std::size_t>(ne);
                b.mConn = std::move(conn);
                blocks.push_back(std::move(b));
            }
        }
    }

    Mesh mesh;
    const std::size_t npoints = coords.size() / 3;
    NDArray pts = NDArray::Uninit(DType::Float64, {npoints, 3});
    if (npoints > 0)
        std::memcpy(pts.Data(), coords.data(), npoints * 3 * sizeof(double));
    mesh.AssignPoints(std::move(pts));

    for (auto& b : blocks) {
        if (b.mKind == 0)
            mesh.AddCellBlock(b.mType, std::move(b.mConn));
        else if (b.mKind == 1)
            mesh.AddPolygonBlock(b.mType, std::move(b.mPolygonRows));
        else
            mesh.AddPolyhedronBlock(b.mType, std::move(b.mPolyhedronCells));
    }

    if (pLayout != nullptr) {
        // Part order: first-seen, matching part_point_ranges/the file itself.
        for (const auto& [pid, range] : part_point_ranges) {
            EnsightPartLayout pl;
            pl.mPartId = pid;
            pl.mPointOffset = range.first;
            pl.mNumPoints = range.second;
            for (std::size_t i = 0; i < blocks.size(); ++i)
                if (blocks[i].mPartId == pid)
                    pl.mBlocks.emplace_back(i, blocks[i].mNumCells);
            pLayout->push_back(std::move(pl));
        }
    }

    if (num_parts >= 2) {
        std::vector<NDArray> tags;
        tags.reserve(blocks.size());
        for (const auto& b : blocks) {
            NDArray a = NDArray::Uninit(DType::Int64, {b.mNumCells});
            std::int64_t* ap = a.As<std::int64_t>();
            for (std::size_t i = 0; i < b.mNumCells; ++i)
                ap[i] = b.mPartId;
            tags.push_back(std::move(a));
        }
        mesh.AddCellData("ensight:part", std::move(tags));
    }

    return mesh;
}

// ---------------------------------------------------------------------------
// Variable (point_data/cell_data) reading
// ---------------------------------------------------------------------------

/**
 * @brief Reads one EnSight Gold variable file against the geometry's own
 *        part/block layout.
 *
 * A variable file mirrors the geometry file's structure exactly -- the same
 * `part`/id sequence, the same per-part `coordinates` (per-node) or
 * element-type-keyword (per-element) sections in the same order -- but
 * carries no counts of its own (a part's point/cell counts are only ever
 * given once, in the geometry file) and no `node id`/`element id` header.
 * Multi-component (vector) data is component-major (every X, then every Y,
 * then every Z), the same convention geometry coordinates use.
 *
 * @param rCur cursor over the variable file (ascii or binary; the caller
 *        selects the concrete cursor type, as `read_ensight` does for the
 *        geometry file)
 * @param PerNode point (`coordinates`) sections when true, element-type
 *        sections when false
 * @param NumComponents 1 for a scalar variable, 3 for a vector
 * @param rLayout the geometry's own per-part layout, in file order
 * @param TotalPoints total point count (only used when `PerNode`)
 * @param pPointOut filled when `PerNode`; otherwise untouched
 * @param pCellOut filled (one entry per Mesh cell block that has one) when
 *        `!PerNode`; otherwise untouched. Must already have as many entries
 *        as the mesh has cell blocks.
 * @throws ReadError if the file's part/type sequence does not match the
 *         geometry's own.
 */
void ensight_read_variable_file(EnsightCursor& rCur, bool PerNode, std::size_t NumComponents,
                                const std::vector<EnsightPartLayout>& rLayout,
                                std::size_t TotalPoints, NDArray* pPointOut,
                                std::vector<NDArray>* pCellOut) {
    constexpr std::int64_t plausible_max = 100000000;
    rCur.NextRecord();  // description line

    NDArray point_out;
    if (PerNode)
        point_out = NDArray::Uninit(DType::Float64, NumComponents == 1
                                                         ? std::vector<std::size_t>{TotalPoints}
                                                         : std::vector<std::size_t>{TotalPoints,
                                                                                    NumComponents});
    double* pp = PerNode ? point_out.As<double>() : nullptr;

    for (const EnsightPartLayout& part : rLayout) {
        std::string rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "part"))
            throw ReadError("EnSight: expected 'part' record in variable file, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/true);
        const std::int64_t pid = rCur.NextInt();
        if (pid != part.mPartId)
            throw ReadError(
                "EnSight: variable file's part sequence does not match the geometry's");

        if (PerNode) {
            rec = rCur.NextRecord();
            if (!ensight_starts_with(rec, "coordinates"))
                throw ReadError("EnSight: expected 'coordinates' record in variable file, got: " +
                                rec);
            std::vector<double> comp(part.mNumPoints);
            for (std::size_t c = 0; c < NumComponents; ++c) {
                rCur.ReadFloats(part.mNumPoints, comp.data());
                for (std::size_t i = 0; i < part.mNumPoints; ++i)
                    pp[(part.mPointOffset + i) * NumComponents + c] = comp[i];
            }
            continue;
        }

        for (const auto& [block_index, num_cells] : part.mBlocks) {
            std::string kw = rCur.NextRecord();  // element type; not re-validated by name
            (void)kw;
            NDArray block(DType::Float64, NumComponents == 1
                                              ? std::vector<std::size_t>{num_cells}
                                              : std::vector<std::size_t>{num_cells,
                                                                         NumComponents});
            double* bp = block.As<double>();
            std::vector<double> comp(num_cells);
            for (std::size_t c = 0; c < NumComponents; ++c) {
                rCur.ReadFloats(num_cells, comp.data());
                for (std::size_t i = 0; i < num_cells; ++i)
                    bp[i * NumComponents + c] = comp[i];
            }
            (*pCellOut)[block_index] = std::move(block);
        }
    }

    if (PerNode)
        *pPointOut = std::move(point_out);
}

/// Dispatches to the ascii/binary cursor, mirroring read_ensight's own
/// geometry-file dispatch.
void ensight_read_variable_file_auto(const std::string& rPath, bool PerNode,
                                     std::size_t NumComponents,
                                     const std::vector<EnsightPartLayout>& rLayout,
                                     std::size_t TotalPoints, NDArray* pPointOut,
                                     std::vector<NDArray>* pCellOut) {
    const detail::FileSource source = ensight_read_whole_file(rPath, "variable file");
    const std::string_view data = source.View();
    if (ensight_starts_with(data, "Fortran Binary"))
        throw ReadError("EnSight: Fortran-binary variable files are not supported");
    if (data.size() >= 80 && ensight_starts_with(data, "C Binary")) {
        EnsightBinaryCursor cur(data);
        ensight_read_variable_file(cur, PerNode, NumComponents, rLayout, TotalPoints, pPointOut,
                                   pCellOut);
        return;
    }
    EnsightAsciiCursor cur(data);
    ensight_read_variable_file(cur, PerNode, NumComponents, rLayout, TotalPoints, pPointOut,
                               pCellOut);
}

}  // namespace

Mesh read_ensight(const std::string& rPath) { return read_ensight(rPath, ReadOptions{}); }

Mesh read_ensight(const std::string& rPath, const ReadOptions& rOptions) {
    const bool have_case = ensight_has_suffix(rPath, ".case");
    EnsightCaseInfo case_info;
    std::string geo_path = rPath;
    if (have_case) {
        case_info = ensight_parse_case(rPath);
        geo_path = case_info.mGeoPath;
    }

    // The source outlives both cursors below, which only hold views into it.
    const detail::FileSource source = ensight_read_whole_file(geo_path, "geometry file");
    const std::string_view data = source.View();
    if (ensight_starts_with(data, "Fortran Binary"))
        throw ReadError("EnSight: Fortran-binary geometry files are not supported");

    std::vector<EnsightPartLayout> layout;
    Mesh mesh;
    if (data.size() >= 80 && ensight_starts_with(data, "C Binary")) {
        EnsightBinaryCursor cur(data);
        mesh = ensight_parse_geo(cur, &layout);
    } else {
        EnsightAsciiCursor cur(data);
        mesh = ensight_parse_geo(cur, &layout);
    }

    if (!have_case || case_info.mVariables.empty() || !rOptions.WantsAnyData())
        return mesh;

    // Which step's variable files to read. A file with no TIME section (the
    // static-geometry-plus-static-variables case) has exactly one step; a
    // non-default mTimeStep against it is a real request this reader cannot
    // honour, so it is refused rather than silently answering step 0.
    std::size_t step = 0;
    if (!case_info.mTimeValues.empty()) {
        step = rOptions.ResolveTimeStep(case_info.mTimeValues.size());
    } else if (rOptions.mTimeStep != 0) {
        throw ReadError(
            "EnSight: mTimeStep requested but the case file has no TIME section to resolve it "
            "against");
    }
    const long file_number =
        case_info.mFileNameStart + static_cast<long>(step) * case_info.mFileNameIncrement;

    const std::string dir = ensight_dirname(rPath);
    for (const EnsightVariableEntry& var : case_info.mVariables) {
        if (!rOptions.WantsArray(var.mName))
            continue;
        const bool per_node = var.mKind.find("per node") != std::string::npos;
        const bool per_element = var.mKind.find("per element") != std::string::npos;
        if (!per_node && !per_element)
            continue;  // a kind this reader does not (yet) understand
        const std::size_t ncomp = ensight_starts_with(var.mKind, "vector") ? 3 : 1;
        const std::string resolved = var.mFilePattern.find('*') != std::string::npos
                                         ? ensight_resolve_wildcard(var.mFilePattern, file_number)
                                         : var.mFilePattern;
        const std::string var_path = dir + resolved;

        if (per_node) {
            NDArray arr;
            ensight_read_variable_file_auto(var_path, true, ncomp, layout, mesh.NumPoints(), &arr,
                                            nullptr);
            mesh.AddPointData(var.mName, std::move(arr));
        } else {
            std::vector<NDArray> blocks(mesh.NumCellBlocks());
            ensight_read_variable_file_auto(var_path, false, ncomp, layout, 0, nullptr, &blocks);
            mesh.AddCellData(var.mName, std::move(blocks));
        }
    }

    return mesh;
}

MeshMetadata read_ensight_metadata(const std::string& rPath, const ReadOptions& /*rOptions*/) {
    if (!ensight_has_suffix(rPath, ".case"))
        throw ReadError("EnSight: metadata needs a .case file (a bare geometry file has no TIME)");
    const EnsightCaseInfo info = ensight_parse_case(rPath);

    // No native header-only shape scan (unlike CGNS/Gmsh 4.1): the same
    // full-read-plus-override shape Exodus's own metadata has.
    MeshMetadata meta = metadata_from_mesh(read_ensight(rPath, ReadOptions{}));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "ensight";
    meta.mTimeValues = info.mTimeValues;
    return meta;
}

namespace {

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void ensight_append_str80(std::vector<char>& rOut, const std::string& rStr) {
    char buf[80] = {};
    rStr.copy(buf, std::min<std::size_t>(rStr.size(), 79));
    rOut.insert(rOut.end(), buf, buf + 80);
}

void ensight_append_i32(std::vector<char>& rOut, std::int64_t v) {
    const std::int32_t i = static_cast<std::int32_t>(v);
    const char* p = reinterpret_cast<const char*>(&i);
    rOut.insert(rOut.end(), p, p + 4);
}

// Validate the mesh and return one keyword entry per cell block.
// The two ragged keywords have no fixed node count, so they cannot live in
// ensight_type_table() (which is keyed on a meshio type name and carries one).
// `mNumNodes < 0` is the marker the geo writers branch on.
const EnsightTypeEntry kEnsightNsided = {"nsided", "", -1};
const EnsightTypeEntry kEnsightNfaced = {"nfaced", "", -2};

std::vector<const EnsightTypeEntry*> ensight_writable_blocks(const Mesh& rMesh) {
    std::vector<const EnsightTypeEntry*> entries;
    for (const auto cb : rMesh.CellRange()) {
        // Ragged blocks are EnSight's own `nsided` / `nfaced`, whose wire
        // format is the exact inverse of what the reader already parses: a
        // per-cell count run, then (for nfaced) a per-face size run, then the
        // node ids. No orientation contract, no global face table.
        if (cb.IsPolyhedron()) {
            entries.push_back(&kEnsightNfaced);
            continue;
        }
        if (cb.IsRagged()) {
            entries.push_back(&kEnsightNsided);
            continue;
        }
        const EnsightTypeEntry* entry = ensight_entry_from_meshio(cb.Type());
        if (entry == nullptr)
            throw WriteError("EnSight: cell type '" + cb.Type() + "' has no EnSight keyword");
        entries.push_back(entry);
    }
    return entries;
}

void ensight_write_geo_ascii(std::ostream& rOs, const Mesh& rMesh,
                             const std::vector<const EnsightTypeEntry*>& rEntries) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    std::string out;
    out.reserve(200 + np * 42);
    out += "EnSight Gold Geometry File\n";
    out += detail::provenance_lines(detail::SlotTier::Bounded)[0] + "\n";
    out += "node id assign\n";
    out += "element id assign\n";
    out += "part\n";

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%10d\n", 1);
    out += buf;
    out += "Mesh\n";
    out += "coordinates\n";
    std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(np));
    out += buf;
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t i = 0; i < np; ++i) {
            const double v = c < dim ? detail::read_double(points, i * dim + c) : 0.0;
            std::snprintf(buf, sizeof(buf), "%12.5e\n", v);
            out += buf;
        }
    }

    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const EnsightTypeEntry* entry = rEntries[bi];
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t ne = cb.NumCells();
        const NDArray& conn = cb.Conn();
        const std::vector<int>* perm = ensight_permutation(cb.Type());

        out += entry->mKeyword;
        out += "\n";
        std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(ne));
        out += buf;
        if (entry->mNumNodes == -2) {
            // nfaced: faces-per-cell, then nodes-per-face, then the node ids --
            // the three runs read_faces consumes, in the same order.
            for (std::size_t r = 0; r < ne; ++r) {
                std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(cb.NumFaces(r)));
                out += buf;
            }
            for (std::size_t r = 0; r < ne; ++r)
                for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                    std::snprintf(buf, sizeof(buf), "%10lld\n",
                                  static_cast<long long>(cb.Face(r, f).second));
                    out += buf;
                }
            for (std::size_t r = 0; r < ne; ++r)
                for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                    const auto face = cb.Face(r, f);
                    for (std::size_t k = 0; k < face.second; ++k) {
                        std::snprintf(buf, sizeof(buf), "%10lld",
                                      static_cast<long long>(face.first[k]) + 1);
                        out += buf;
                    }
                    out += "\n";
                }
            continue;
        }
        if (entry->mNumNodes == -1) {
            // nsided: nodes-per-cell, then the node ids.
            for (std::size_t r = 0; r < ne; ++r) {
                std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(cb.RowSize(r)));
                out += buf;
            }
            for (std::size_t r = 0; r < ne; ++r) {
                for (std::size_t k = 0; k < cb.RowSize(r); ++k) {
                    std::snprintf(buf, sizeof(buf), "%10lld",
                                  static_cast<long long>(cb.Row(r)[k]) + 1);
                    out += buf;
                }
                out += "\n";
            }
            continue;
        }
        for (std::size_t r = 0; r < ne; ++r) {
            for (std::size_t j = 0; j < npc; ++j) {
                const std::size_t src = perm != nullptr ? static_cast<std::size_t>((*perm)[j]) : j;
                const long long v =
                    static_cast<long long>(detail::read_int(conn, r * npc + src)) + 1;
                std::snprintf(buf, sizeof(buf), "%10lld", v);
                out += buf;
            }
            out += "\n";
        }
    }

    rOs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

void ensight_write_geo_binary(std::ostream& rOs, const Mesh& rMesh,
                              const std::vector<const EnsightTypeEntry*>& rEntries) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    constexpr std::size_t i32_max =
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    if (np > i32_max)
        throw WriteError("EnSight: mesh too large for 32-bit binary EnSight output");

    std::vector<char> out;
    out.reserve(80 * 8 + np * 12 + 64);
    ensight_append_str80(out, "C Binary");
    ensight_append_str80(out, "EnSight Gold Geometry File");
    ensight_append_str80(out, detail::provenance_lines(detail::SlotTier::Bounded)[0]);
    ensight_append_str80(out, "node id assign");
    ensight_append_str80(out, "element id assign");
    ensight_append_str80(out, "part");
    ensight_append_i32(out, 1);
    ensight_append_str80(out, "Mesh");
    ensight_append_str80(out, "coordinates");
    ensight_append_i32(out, static_cast<std::int64_t>(np));
    {
        std::vector<float> col(np * 3);
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t i = 0; i < np; ++i)
                col[c * np + i] =
                    c < dim ? static_cast<float>(detail::read_double(points, i * dim + c)) : 0.0f;
        const char* p = reinterpret_cast<const char*>(col.data());
        out.insert(out.end(), p, p + col.size() * sizeof(float));
    }

    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const EnsightTypeEntry* entry = rEntries[bi];
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t ne = cb.NumCells();
        if (ne > i32_max)
            throw WriteError("EnSight: mesh too large for 32-bit binary EnSight output");
        const NDArray& conn = cb.Conn();
        const std::vector<int>* perm = ensight_permutation(cb.Type());

        ensight_append_str80(out, entry->mKeyword);
        ensight_append_i32(out, static_cast<std::int64_t>(ne));
        if (entry->mNumNodes < 0) {
            // nsided / nfaced: the same three (or two) runs as the ASCII path,
            // as int32 -- which is exactly what read_binary_faces consumes.
            std::vector<std::int32_t> counts;
            std::vector<std::int32_t> sizes;
            std::vector<std::int32_t> nodes;
            const bool faced = entry->mNumNodes == -2;
            for (std::size_t r = 0; r < ne; ++r) {
                if (faced) {
                    counts.push_back(static_cast<std::int32_t>(cb.NumFaces(r)));
                    for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                        const auto face = cb.Face(r, f);
                        sizes.push_back(static_cast<std::int32_t>(face.second));
                        for (std::size_t k = 0; k < face.second; ++k)
                            nodes.push_back(static_cast<std::int32_t>(face.first[k]) + 1);
                    }
                } else {
                    counts.push_back(static_cast<std::int32_t>(cb.RowSize(r)));
                    for (std::size_t k = 0; k < cb.RowSize(r); ++k)
                        nodes.push_back(static_cast<std::int32_t>(cb.Row(r)[k]) + 1);
                }
            }
            auto append = [&out](const std::vector<std::int32_t>& v) {
                if (v.empty())
                    return;
                const char* q = reinterpret_cast<const char*>(v.data());
                out.insert(out.end(), q, q + v.size() * sizeof(std::int32_t));
            };
            append(counts);
            if (faced)
                append(sizes);
            append(nodes);
            continue;
        }
        std::vector<std::int32_t> flat(ne * npc);
        for (std::size_t r = 0; r < ne; ++r)
            for (std::size_t j = 0; j < npc; ++j) {
                const std::size_t src = perm != nullptr ? static_cast<std::size_t>((*perm)[j]) : j;
                flat[r * npc + j] =
                    static_cast<std::int32_t>(detail::read_int(conn, r * npc + src)) + 1;
            }
        const char* p = reinterpret_cast<const char*>(flat.data());
        out.insert(out.end(), p, p + flat.size() * sizeof(std::int32_t));
    }

    rOs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

}  // namespace

void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary) {
    bool ok = false;
    auto paths = ensight_case_geo_paths(rPath, ok);
    if (!ok)
        throw WriteError("EnSight: must specify a .case or .geo file");
    const std::string& case_path = paths.first;
    const std::string& geo_path = paths.second;

    if (rMesh.PointDim() > 3)
        throw WriteError("EnSight: points must have at most three components");
    const std::vector<const EnsightTypeEntry*> entries = ensight_writable_blocks(rMesh);

    {
        std::ofstream cf(case_path, std::ios::binary);
        if (!cf)
            throw WriteError("Could not open file for writing: " + case_path);
        std::string out;
        out += "FORMAT\n";
        out += "type: ensight gold\n";
        out += "\n";
        out += "GEOMETRY\n";
        out += "model: " + ensight_basename(geo_path) + "\n";
        cf.write(out.data(), static_cast<std::streamsize>(out.size()));
    }

    std::ofstream gf(geo_path, std::ios::binary);
    if (!gf)
        throw WriteError("Could not open file for writing: " + geo_path);
    if (binary)
        ensight_write_geo_binary(gf, rMesh, entries);
    else
        ensight_write_geo_ascii(gf, rMesh, entries);
}

}  // namespace meshioplusplus
