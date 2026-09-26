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
#include <bit>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
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
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

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
            const char* end = nullptr;
            pDst[i] = detail::parse_double(start, end);
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
    // Data is the whole file; a geometry cursor starts after the leading
    // "C Binary" 80-char record, a variable cursor at its description.
    // `Swap` seeds the byte order when the file says it (Fortran markers).
    explicit EnsightBinaryCursor(std::string_view data, std::size_t Start = 80, bool Swap = false)
        : mData(data), mPos(Start), mSwap(Swap) {}

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
    bool mSwap;
};

// ---------------------------------------------------------------------------
// Fortran binary: the C-binary stream with every WRITE framed as a Fortran
// sequential unformatted record (its byte length before and after it, 4 or 8
// bytes, in the writer's byte order). Joining the payloads gives the C-binary
// stream back, whatever the producer's record boundaries were.
// ---------------------------------------------------------------------------

enum class EnsightEncoding { Ascii, CBinary, Fortran };

struct EnsightUnframed {
    std::string mBytes;
    bool mSwap = false;  ///< the markers (so the values) are in the other byte order
};

/// The joined payloads of a Fortran-framed `Data`, or nothing when it is not
/// framed. `rWhat` names the file in errors.
std::optional<EnsightUnframed> ensight_unframe_fortran(std::string_view Data,
                                                       const std::string& rWhat) {
    const auto layout = detail::sniff_fortran_records(Data.data(), Data.size());
    if (!layout)
        return std::nullopt;
    const std::vector<detail::FortranRecord> records =
        detail::fortran_records(Data.data(), Data.size(), *layout, "EnSight " + rWhat);
    EnsightUnframed out;
    std::size_t total = 0;
    for (const detail::FortranRecord& r : records)
        total += r.mSize;
    out.mBytes.reserve(total);
    for (const detail::FortranRecord& r : records)
        out.mBytes.append(Data.data() + r.mOffset, r.mSize);
    out.mSwap = layout->mBigEndian != (std::endian::native == std::endian::big);
    return out;
}

/// Whether `Data` is a Fortran-binary geometry file: its first record is the
/// 80-byte "Fortran Binary" string.
bool ensight_is_fortran_geometry(std::string_view Data) {
    const auto layout = detail::sniff_fortran_records(Data.data(), Data.size());
    if (!layout)
        return false;
    const std::size_t m = static_cast<std::size_t>(layout->mMarkerBytes);
    return Data.size() >= m + 80 && ensight_starts_with(Data.substr(m), "Fortran Binary");
}

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
    std::vector<std::pair<std::string, double>> mConstants;  // "constant per case:" entries
};

/// Splits a whitespace-separated record and drops its leading run of pure
/// integer tokens (a `[ts] [fs]` prefix) -- the same rule `model:`/variable
/// lines both use to make the leading timeset/fileset optional.
std::vector<std::string> ensight_tokens_after_leading_ints(const std::string& rValue) {
    auto toks = detail::make_classic_istringstream(rValue);
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
    auto stream = detail::make_classic_istringstream(data);
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
        } else if (section == "VARIABLE" && ensight_starts_with(line, "constant per case:")) {
            // Inline, not a file reference: `[ts] [fs] <name> <value>`.
            const std::vector<std::string> toks =
                ensight_tokens_after_leading_ints(line.substr(std::strlen("constant per case:")));
            if (toks.size() < 2)
                throw ReadError("EnSight: malformed 'constant per case' line: " + line);
            std::string joined;
            for (std::size_t i = 0; i + 1 < toks.size(); ++i)
                joined += (i ? " " : "") + toks[i];
            info.mConstants.emplace_back(joined, detail::parse_double(toks.back()));
        } else if (section == "VARIABLE") {
            static const char* kKinds[] = {
                "scalar per node:",       "vector per node:",       "tensor symm per node:",
                "tensor asym per node:",  "scalar per element:",    "vector per element:",
                "tensor symm per element:", "tensor asym per element:"};
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
                info.mFileNameStart =
                    std::strtol(line.c_str() + std::strlen("filename start number:"), nullptr, 10);
            } else if (ensight_starts_with(line, "filename increment:")) {
                info.mFileNameIncrement =
                    std::strtol(line.c_str() + std::strlen("filename increment:"), nullptr, 10);
            } else if (ensight_starts_with(line, "time values:")) {
                in_time_values = true;
                const std::string rest = ensight_trim(line.substr(std::strlen("time values:")));
                auto iss = detail::make_classic_istringstream(rest);
                double v;
                while (iss >> v)
                    info.mTimeValues.push_back(v);
            } else if (in_time_values) {
                auto iss = detail::make_classic_istringstream(line);
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
    auto num = detail::make_classic_ostringstream();
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
    auto iss = detail::make_classic_istringstream(rRecord);
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
        point_out = NDArray::Uninit(DType::Float64,
                                    NumComponents == 1
                                        ? std::vector<std::size_t>{TotalPoints}
                                        : std::vector<std::size_t>{TotalPoints, NumComponents});
    double* pp = PerNode ? point_out.As<double>() : nullptr;

    for (const EnsightPartLayout& part : rLayout) {
        std::string rec = rCur.NextRecord();
        if (!ensight_starts_with(rec, "part"))
            throw ReadError("EnSight: expected 'part' record in variable file, got: " + rec);
        rCur.CheckSwap(plausible_max, /*PreferSmaller=*/true);
        const std::int64_t pid = rCur.NextInt();
        if (pid != part.mPartId)
            throw ReadError("EnSight: variable file's part sequence does not match the geometry's");

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
                                              : std::vector<std::size_t>{num_cells, NumComponents});
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

/// A binary variable file starts with its description record; files meshio++
/// wrote before v16.17.0 put a "C Binary" record before it, which no other
/// reader expects. Where the binary cursor starts in `Data`.
std::size_t ensight_variable_start(std::string_view Data) {
    const bool standard = Data.size() >= 160 && ensight_starts_with(Data.substr(80), "part");
    return !standard && ensight_starts_with(Data, "C Binary") ? 80 : 0;
}

/// Whether a variable file is plain text: no NUL in its first KiB and its
/// description line ends within 80 characters. A binary variable file's
/// description is an 80-byte record, so its second record starts at byte 80.
bool ensight_variable_looks_ascii(std::string_view Data) {
    const std::string_view head = Data.substr(0, std::min<std::size_t>(Data.size(), 1024));
    if (head.find('\0') != std::string_view::npos)
        return false;
    const std::size_t eol = head.find('\n');
    return eol != std::string_view::npos && eol < 81;
}

/// Reads a variable file in the encoding of its geometry file, as EnSight
/// defines it (a variable file carries no format record of its own). A text
/// variable file next to a binary geometry file, which hand-made cases hold
/// and meshio++ read before v16.17.0, is still read as text.
void ensight_read_variable_file_auto(const std::string& rPath, EnsightEncoding Encoding,
                                     bool PerNode, std::size_t NumComponents,
                                     const std::vector<EnsightPartLayout>& rLayout,
                                     std::size_t TotalPoints, NDArray* pPointOut,
                                     std::vector<NDArray>* pCellOut) {
    const detail::FileSource source = ensight_read_whole_file(rPath, "variable file");
    const std::string_view data = source.View();
    const bool binary_layout =
        (data.size() >= 160 && ensight_starts_with(data.substr(80), "part") &&
         data.substr(0, 80).find('\n') == std::string_view::npos) ||
        ensight_starts_with(data, "C Binary");
    if (Encoding != EnsightEncoding::Ascii && !binary_layout && ensight_variable_looks_ascii(data))
        Encoding = EnsightEncoding::Ascii;
    if (Encoding == EnsightEncoding::Fortran) {
        const auto unframed = ensight_unframe_fortran(data, "variable file");
        if (!unframed)
            throw ReadError("EnSight: variable file '" + rPath +
                            "' is not Fortran binary like its geometry file");
        EnsightBinaryCursor cur(unframed->mBytes, 0, unframed->mSwap);
        ensight_read_variable_file(cur, PerNode, NumComponents, rLayout, TotalPoints, pPointOut,
                                   pCellOut);
        return;
    }
    if (Encoding == EnsightEncoding::CBinary) {
        EnsightBinaryCursor cur(data, ensight_variable_start(data));
        ensight_read_variable_file(cur, PerNode, NumComponents, rLayout, TotalPoints, pPointOut,
                                   pCellOut);
        return;
    }
    EnsightAsciiCursor cur(data);
    ensight_read_variable_file(cur, PerNode, NumComponents, rLayout, TotalPoints, pPointOut,
                               pCellOut);
}

/// Component count for a `.case` `VARIABLE` kind: 1 scalar, 3 vector, 6
/// tensor symm, 9 tensor asym.
std::size_t ensight_variable_ncomp(const std::string& rKind) {
    if (ensight_starts_with(rKind, "vector"))
        return 3;
    if (ensight_starts_with(rKind, "tensor symm"))
        return 6;
    if (ensight_starts_with(rKind, "tensor asym"))
        return 9;
    return 1;
}

/**
 * @brief Reorders a `tensor symm` array's last two components in place.
 *
 * EnSight Gold's own file order for `tensor symm` is `11 22 33 12 13 23`
 * (xx, yy, zz, xy, xz, yz); meshio++'s six-component symmetric-tensor
 * convention (see `doc/mesh_data_model.md`) is `xx yy zz xy yz zx` -- the
 * same six values, with the last two swapped (`zx` and `xz` are the same
 * component of a symmetric tensor). Verified empirically against
 * ParaView's own EnSight Gold reader (see `doc/formats/ensight.md`); VTK's
 * internal tensor6 order is `xx yy zz xy yz xz`, so `vtkEnSightGoldReader`
 * performs the identical swap on its own read. The swap is its own
 * inverse, so this one function serves both read and write.
 * @param pData Row-major `(Rows, 6)` buffer, reordered in place.
 * @param Rows Number of tensor entries (points or cells).
 */
void ensight_swap_tensor_symm_last_two(double* pData, std::size_t Rows) {
    for (std::size_t r = 0; r < Rows; ++r)
        std::swap(pData[r * 6 + 4], pData[r * 6 + 5]);
}

}  // namespace

Mesh read_ensight(const std::string& rPath) {
    return read_ensight(rPath, ReadOptions{});
}

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

    std::vector<EnsightPartLayout> layout;
    Mesh mesh;
    EnsightEncoding encoding = EnsightEncoding::Ascii;
    if (ensight_is_fortran_geometry(data)) {
        encoding = EnsightEncoding::Fortran;
        const auto unframed = ensight_unframe_fortran(data, "geometry file");
        EnsightBinaryCursor cur(unframed->mBytes, 80, unframed->mSwap);
        mesh = ensight_parse_geo(cur, &layout);
    } else if (data.size() >= 80 && ensight_starts_with(data, "C Binary")) {
        encoding = EnsightEncoding::CBinary;
        EnsightBinaryCursor cur(data);
        mesh = ensight_parse_geo(cur, &layout);
    } else {
        EnsightAsciiCursor cur(data);
        mesh = ensight_parse_geo(cur, &layout);
    }

    if (!have_case || !rOptions.WantsAnyData() ||
        (case_info.mVariables.empty() && case_info.mConstants.empty()))
        return mesh;

    for (const auto& [name, value] : case_info.mConstants) {
        if (!rOptions.WantsArray(name))
            continue;
        NDArray arr(DType::Float64, {std::size_t{1}});
        *arr.As<double>() = value;
        mesh.AddFieldData(name, std::move(arr));
    }
    if (case_info.mVariables.empty())
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
        const std::size_t ncomp = ensight_variable_ncomp(var.mKind);
        const bool tensor_symm = ensight_starts_with(var.mKind, "tensor symm");
        const std::string resolved = var.mFilePattern.find('*') != std::string::npos
                                         ? ensight_resolve_wildcard(var.mFilePattern, file_number)
                                         : var.mFilePattern;
        const std::string var_path = dir + resolved;

        if (per_node) {
            NDArray arr;
            ensight_read_variable_file_auto(var_path, encoding, true, ncomp, layout,
                                            mesh.NumPoints(), &arr, nullptr);
            if (tensor_symm)
                ensight_swap_tensor_symm_last_two(arr.As<double>(), mesh.NumPoints());
            mesh.AddPointData(var.mName, std::move(arr));
        } else {
            std::vector<NDArray> blocks(mesh.NumCellBlocks());
            ensight_read_variable_file_auto(var_path, encoding, false, ncomp, layout, 0, nullptr,
                                            &blocks);
            if (tensor_symm)
                for (NDArray& blk : blocks)
                    if (blk.Size() > 0)
                        ensight_swap_tensor_symm_last_two(blk.As<double>(), blk.Shape()[0]);
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

/// A binary EnSight file as its records: C binary writes them back to back,
/// Fortran binary frames each one as the WRITE of a Fortran producer (its byte
/// length before and after it, 4 bytes in the native byte order). A record is
/// one string, one count, or one array: the layout vtkEnSightGoldBinaryReader
/// reads Fortran files in (each coordinate component and each variable
/// component is its own record).
class EnsightRecordWriter {
public:
    explicit EnsightRecordWriter(bool Fortran) : mFortran(Fortran) {}

    void Str80(const std::string& rStr) {
        char buf[80] = {};
        rStr.copy(buf, std::min<std::size_t>(rStr.size(), 79));
        Record(buf, 80);
    }

    void Int(std::int64_t Value) {
        const std::int32_t i = static_cast<std::int32_t>(Value);
        Record(reinterpret_cast<const char*>(&i), 4);
    }

    void Ints(const std::vector<std::int32_t>& rValues) {
        Record(reinterpret_cast<const char*>(rValues.data()), rValues.size() * 4);
    }

    void Floats(const float* pValues, std::size_t Count) {
        Record(reinterpret_cast<const char*>(pValues), Count * 4);
    }

    const std::vector<char>& Bytes() const { return mOut; }
    std::vector<char>& Bytes() { return mOut; }

private:
    void Record(const char* pData, std::size_t Size) {
        if (mFortran) {
            if (Size > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
                throw WriteError("EnSight: a Fortran-binary record over 2 GiB cannot be written");
            Marker(Size);
        }
        mOut.insert(mOut.end(), pData, pData + Size);
        if (mFortran)
            Marker(Size);
    }

    void Marker(std::size_t Size) {
        const std::int32_t n = static_cast<std::int32_t>(Size);
        const char* p = reinterpret_cast<const char*>(&n);
        mOut.insert(mOut.end(), p, p + 4);
    }

    bool mFortran;
    std::vector<char> mOut;
};

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
            detail::snprintf_c(buf, sizeof(buf), "%12.5e\n", v);
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
                              const std::vector<const EnsightTypeEntry*>& rEntries, bool Fortran) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    constexpr std::size_t i32_max =
        static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());
    if (np > i32_max)
        throw WriteError("EnSight: mesh too large for 32-bit binary EnSight output");

    EnsightRecordWriter out(Fortran);
    out.Bytes().reserve(80 * 8 + np * 12 + 64);
    out.Str80(Fortran ? "Fortran Binary" : "C Binary");
    out.Str80("EnSight Gold Geometry File");
    out.Str80(detail::provenance_lines(detail::SlotTier::Bounded)[0]);
    out.Str80("node id assign");
    out.Str80("element id assign");
    out.Str80("part");
    out.Int(1);
    out.Str80("Mesh");
    out.Str80("coordinates");
    out.Int(static_cast<std::int64_t>(np));
    {
        std::vector<float> col(np * 3);
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t i = 0; i < np; ++i)
                col[c * np + i] =
                    c < dim ? static_cast<float>(detail::read_double(points, i * dim + c)) : 0.0f;
        for (std::size_t c = 0; c < 3; ++c)
            out.Floats(col.data() + c * np, np);
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

        out.Str80(entry->mKeyword);
        out.Int(static_cast<std::int64_t>(ne));
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
                if (!v.empty())
                    out.Ints(v);
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
        out.Ints(flat);
    }

    rOs.write(out.Bytes().data(), static_cast<std::streamsize>(out.Bytes().size()));
}

// ---------------------------------------------------------------------------
// VARIABLE section writing
// ---------------------------------------------------------------------------

// One point_data/cell_data array this write collected. `mNumComponents` is
// the *written* component count (2 pads to 3, matching the geometry writer's
// own 2D-coordinate padding); `mKind` is the exact `VARIABLE` line keyword.
struct EnsightVariableToWrite {
    std::string mName;
    std::string mKind;  // "scalar", "vector", "tensor symm" or "tensor asym"
    std::size_t mNumComponents;
    bool mPerNode;
};

// EnSight's kind word and written component count for a data array's actual
// component count; `false` when the count has no EnSight representation.
bool ensight_kind_for_ncomp(std::size_t NumComponents, std::string& rKind,
                            std::size_t& rWritten) {
    switch (NumComponents) {
        case 1:
            rKind = "scalar";
            rWritten = 1;
            return true;
        case 2:
        case 3:
            rKind = "vector";
            rWritten = 3;
            return true;
        case 6:
            rKind = "tensor symm";
            rWritten = 6;
            return true;
        case 9:
            rKind = "tensor asym";
            rWritten = 9;
            return true;
        default:
            return false;
    }
}

// The variable-file extension this write uses for a kind/location pair --
// arbitrary (the `.case` file's own declared kind is what the reader goes
// by, not the extension), but kept distinct per kind for readability on
// disk.
std::string ensight_variable_extension(const std::string& rKind, bool PerNode) {
    std::string ext = "scl";
    if (rKind == "vector")
        ext = "vec";
    else if (rKind == "tensor symm")
        ext = "tsym";
    else if (rKind == "tensor asym")
        ext = "tasym";
    return PerNode ? ext : ("e" + ext);
}

/**
 * @brief One written component's values, in EnSight file order.
 *
 * `tensor symm`'s last two components are transposed relative to meshio++'s
 * own `xx yy zz xy yz zx` convention -- `ensight_swap_tensor_symm_last_two`
 * documents why, and undoes the same transposition on read; a missing
 * trailing component (the vector-padding case) reads as `0.0`.
 * @param rMesh The mesh being written.
 * @param rVar Which array, kind and location.
 * @param Comp 0-based component index, in EnSight file order.
 * @param BlockIndex Cell block index; ignored when `rVar.mPerNode`.
 * @return One value per point (or per cell of that block).
 */
std::vector<double> ensight_variable_column(const Mesh& rMesh, const EnsightVariableToWrite& rVar,
                                            std::size_t Comp, std::size_t BlockIndex) {
    const std::size_t mio_comp = (rVar.mKind == "tensor symm" && (Comp == 4 || Comp == 5))
                                     ? (Comp == 4 ? 5 : 4)
                                     : Comp;
    const NDArray& arr =
        rVar.mPerNode ? rMesh.PointData(rVar.mName) : rMesh.CellData(rVar.mName, BlockIndex);
    const std::size_t stored_ncomp = arr.Shape().size() >= 2 ? arr.Shape()[1] : 1;
    const std::size_t n =
        rVar.mPerNode ? rMesh.NumPoints() : rMesh.Cells(BlockIndex).NumCells();
    std::vector<double> col(n);
    for (std::size_t i = 0; i < n; ++i)
        col[i] =
            mio_comp < stored_ncomp ? detail::read_double(arr, i * stored_ncomp + mio_comp) : 0.0;
    return col;
}

void ensight_write_variable_ascii(std::ostream& rOs, const Mesh& rMesh,
                                  const std::vector<const EnsightTypeEntry*>& rEntries,
                                  const EnsightVariableToWrite& rVar) {
    std::string out;
    out += "variable\n";  // description line; not re-read
    out += "part\n";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%10d\n", 1);
    out += buf;

    auto write_col = [&](const std::vector<double>& col) {
        for (double v : col) {
            detail::snprintf_c(buf, sizeof(buf), "%12.5e\n", v);
            out += buf;
        }
    };

    if (rVar.mPerNode) {
        out += "coordinates\n";
        for (std::size_t c = 0; c < rVar.mNumComponents; ++c)
            write_col(ensight_variable_column(rMesh, rVar, c, 0));
    } else {
        for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
            out += rEntries[bi]->mKeyword;
            out += "\n";
            for (std::size_t c = 0; c < rVar.mNumComponents; ++c)
                write_col(ensight_variable_column(rMesh, rVar, c, bi));
        }
    }
    rOs.write(out.data(), static_cast<std::streamsize>(out.size()));
}

// A variable file has no format record: it starts with its description, and
// its encoding is its geometry file's (until v16.17.0 a "C Binary" record came
// first, which kept VTK and ParaView from reading the variables).
void ensight_write_variable_binary(std::ostream& rOs, const Mesh& rMesh,
                                   const std::vector<const EnsightTypeEntry*>& rEntries,
                                   const EnsightVariableToWrite& rVar, bool Fortran) {
    EnsightRecordWriter out(Fortran);
    out.Str80("variable");
    out.Str80("part");
    out.Int(1);

    auto append_col = [&](const std::vector<double>& col) {
        std::vector<float> f(col.size());
        for (std::size_t i = 0; i < col.size(); ++i)
            f[i] = static_cast<float>(col[i]);
        out.Floats(f.data(), f.size());
    };

    if (rVar.mPerNode) {
        out.Str80("coordinates");
        for (std::size_t c = 0; c < rVar.mNumComponents; ++c)
            append_col(ensight_variable_column(rMesh, rVar, c, 0));
    } else {
        for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
            out.Str80(rEntries[bi]->mKeyword);
            for (std::size_t c = 0; c < rVar.mNumComponents; ++c)
                append_col(ensight_variable_column(rMesh, rVar, c, bi));
        }
    }
    rOs.write(out.Bytes().data(), static_cast<std::streamsize>(out.Bytes().size()));
}

/// Scans `rMesh`'s data maps for what this writer can express: `point_data`
/// and `cell_data` arrays of 1, 2, 3, 6 or 9 components (anything else is
/// skipped with a warning; a `cell_data` array missing from any cell block
/// is skipped too, since a variable file must cover every block the
/// geometry file does), and single-scalar `field_data` as `constant per
/// case` entries (anything else skipped with a warning).
void ensight_collect_variables(const Mesh& rMesh, std::vector<EnsightVariableToWrite>& rVars,
                               std::vector<std::pair<std::string, double>>& rConstants) {
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& arr = rMesh.PointData(name);
        const std::size_t ncomp = arr.Shape().size() >= 2 ? arr.Shape()[1] : 1;
        std::string kind;
        std::size_t written = 0;
        if (!ensight_kind_for_ncomp(ncomp, kind, written)) {
            log::warn(
                "EnSight: skipping point data '{}' with {} components (1, 2, 3, 6 or 9 are "
                "supported)",
                name, ncomp);
            continue;
        }
        rVars.push_back({name, kind, written, true});
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks()) {
            log::warn("EnSight: skipping cell data '{}': not present on every cell block", name);
            continue;
        }
        const NDArray& first = rMesh.CellData(name, 0);
        const std::size_t ncomp = first.Shape().size() >= 2 ? first.Shape()[1] : 1;
        std::string kind;
        std::size_t written = 0;
        if (!ensight_kind_for_ncomp(ncomp, kind, written)) {
            log::warn(
                "EnSight: skipping cell data '{}' with {} components (1, 2, 3, 6 or 9 are "
                "supported)",
                name, ncomp);
            continue;
        }
        rVars.push_back({name, kind, written, false});
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        const NDArray& arr = rMesh.FieldData(name);
        if (arr.Size() != 1) {
            log::warn(
                "EnSight: skipping field data '{}': only a single scalar can become a "
                "'constant per case'",
                name);
            continue;
        }
        rConstants.emplace_back(name, detail::read_double(arr, 0));
    }
}

}  // namespace

void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary) {
    write_ensight(rPath, rMesh, binary, /*fortran=*/false);
}

void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary, bool fortran) {
    if (fortran && !binary)
        throw WriteError("EnSight: Fortran binary is a binary encoding (binary=true)");
    bool ok = false;
    auto paths = ensight_case_geo_paths(rPath, ok);
    if (!ok)
        throw WriteError("EnSight: must specify a .case or .geo file");
    const std::string& case_path = paths.first;
    const std::string& geo_path = paths.second;
    const std::string base = ensight_basename(geo_path);
    const std::string dir = ensight_dirname(geo_path);

    if (rMesh.PointDim() > 3)
        throw WriteError("EnSight: points must have at most three components");
    const std::vector<const EnsightTypeEntry*> entries = ensight_writable_blocks(rMesh);

    std::vector<EnsightVariableToWrite> vars;
    std::vector<std::pair<std::string, double>> constants;
    ensight_collect_variables(rMesh, vars, constants);

    {
        auto cf = detail::make_classic_ofstream(case_path, std::ios::binary);
        if (!cf)
            throw WriteError("Could not open file for writing: " + case_path);
        std::string out;
        out += "FORMAT\n";
        out += "type: ensight gold\n";
        out += "\n";
        out += "GEOMETRY\n";
        out += "model: " + base + "\n";
        if (!vars.empty() || !constants.empty()) {
            out += "\n";
            out += "VARIABLE\n";
            for (const EnsightVariableToWrite& v : vars) {
                const std::string ext = ensight_variable_extension(v.mKind, v.mPerNode);
                out += v.mKind + " per " + (v.mPerNode ? "node" : "element") + ": " + v.mName +
                       " " + v.mName + "." + ext + "\n";
            }
            for (const auto& [name, value] : constants) {
                char buf[40];
                detail::snprintf_c(buf, sizeof(buf), "%.17g", value);
                out += "constant per case: " + name + " " + buf + "\n";
            }
        }
        cf.write(out.data(), static_cast<std::streamsize>(out.size()));
    }

    auto gf = detail::make_classic_ofstream(geo_path, std::ios::binary);
    if (!gf)
        throw WriteError("Could not open file for writing: " + geo_path);
    if (binary)
        ensight_write_geo_binary(gf, rMesh, entries, fortran);
    else
        ensight_write_geo_ascii(gf, rMesh, entries);

    for (const EnsightVariableToWrite& v : vars) {
        const std::string ext = ensight_variable_extension(v.mKind, v.mPerNode);
        const std::string var_path = dir + v.mName + "." + ext;
        auto vf = detail::make_classic_ofstream(var_path, std::ios::binary);
        if (!vf)
            throw WriteError("Could not open file for writing: " + var_path);
        if (binary)
            ensight_write_variable_binary(vf, rMesh, entries, v, fortran);
        else
            ensight_write_variable_ascii(vf, rMesh, entries, v);
    }
}

}  // namespace meshioplusplus
