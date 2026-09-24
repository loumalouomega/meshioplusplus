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
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// Written as the first comment line of every file the C++ writer produces. The
// reader no longer needs it (it reads any bulk-data deck); it stays so older
// meshio++ releases, whose reader was gated on it, can still read new output.
constexpr const char* kSentinel = "meshioplusplus-cpp-nastran";

// An element card: its meshio type and how many node fields it holds. A fixed
// count (> 0) reads exactly that many fields and ignores the rest of the card
// (THETA, ZOFFS, thicknesses, an orientation vector ...); 0 marks a solid whose
// linear or quadratic form is told apart by the number of node fields given.
struct NasElement {
    const char* mType;
    int mNodes;
};

const std::unordered_map<std::string, NasElement>& nas_elements() {
    static const std::unordered_map<std::string, NasElement> m = {
        {"CELAS1", {"vertex", 1}},    {"CBEAM", {"line", 2}},        {"CBUSH", {"line", 2}},
        {"CBUSH1D", {"line", 2}},     {"CROD", {"line", 2}},         {"CGAP", {"line", 2}},
        {"CBAR", {"line", 2}},        {"CTRIAR", {"triangle", 3}},   {"CTRIA3", {"triangle", 3}},
        {"CTRAX6", {"triangle6", 6}}, {"CTRIAX6", {"triangle6", 6}}, {"CTRIA6", {"triangle6", 6}},
        {"CQUADR", {"quad", 4}},      {"CSHEAR", {"quad", 4}},       {"CQUAD4", {"quad", 4}},
        {"CQUAD8", {"quad8", 8}},     {"CQUAD9", {"quad9", 9}},      {"CTETRA", {"tetra", 0}},
        {"CPYRAM", {"pyramid", 0}},   {"CPYRA", {"pyramid", 0}},     {"CPENTA", {"wedge", 0}},
        {"CHEXA", {"hexahedron", 0}},
    };
    return m;
}

// A solid card's linear and quadratic meshio types and node counts.
struct NasSolid {
    const char* mLinear;
    int mLinearNodes;
    const char* mQuadratic;
    int mQuadraticNodes;
};

const NasSolid* nas_solid(const std::string& rType) {
    static const NasSolid solids[] = {
        {"tetra", 4, "tetra10", 10},
        {"pyramid", 5, "pyramid13", 13},
        {"wedge", 6, "wedge15", 15},
        {"hexahedron", 8, "hexahedron20", 20},
    };
    for (const NasSolid& s : solids)
        if (rType == s.mLinear)
            return &s;
    return nullptr;
}

// meshio type -> the card the writer emits.
const std::unordered_map<std::string, std::string>& meshio_to_nastran() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "CELAS1"},    {"line", "CBAR"},        {"triangle", "CTRIA3"},
        {"triangle6", "CTRIA6"}, {"quad", "CQUAD4"},      {"quad8", "CQUAD8"},
        {"quad9", "CQUAD9"},     {"tetra", "CTETRA"},     {"tetra10", "CTETRA"},
        {"pyramid", "CPYRA"},    {"pyramid13", "CPYRA"},  {"wedge", "CPENTA"},
        {"wedge15", "CPENTA"},   {"hexahedron", "CHEXA"}, {"hexahedron20", "CHEXA"},
    };
    return m;
}

// meshio slot j holds Nastran slot P[j]. hexahedron20 and wedge15 put the
// vertical mid-edges before the top ring in Nastran; both are involutions.
const std::vector<int>& nas_perm_hex20() {
    static const std::vector<int> p = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                       10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
    return p;
}
const std::vector<int>& nas_perm_wedge15() {
    static const std::vector<int> p = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};
    return p;
}
// CTRIAX6/CTRAX6 list corner, mid, corner, mid, corner, mid.
const std::vector<int>& nas_perm_triax6_read() {
    static const std::vector<int> p = {0, 2, 4, 1, 3, 5};
    return p;
}

const std::vector<int>* nas_read_perm(const std::string& rCard, const std::string& rType) {
    if (rCard == "CTRIAX6" || rCard == "CTRAX6")
        return &nas_perm_triax6_read();
    if (rType == "hexahedron20")
        return &nas_perm_hex20();
    if (rType == "wedge15")
        return &nas_perm_wedge15();
    return nullptr;
}

const std::vector<int>* nas_write_perm(const std::string& rType) {
    if (rType == "hexahedron20")
        return &nas_perm_hex20();
    if (rType == "wedge15")
        return &nas_perm_wedge15();
    return nullptr;
}

// Cards read silently although they carry nothing meshio++ keeps: properties,
// materials, loads, constraints, coordinate systems, tables and solution
// parameters. Anything else not read is counted and named in one warning.
bool nas_is_quietly_skipped(const std::string& rKeyword) {
    static const char* const prefixes[] = {
        "P",      "MAT",   "SPC",   "MPC",   "FORCE", "MOMENT", "LOAD",   "TEMP",  "GRAV",
        "RFORCE", "ACCEL", "CORD",  "TABLE", "EIGR",  "EIGC",   "NLPARM", "TSTEP", "FREQ",
        "SUPORT", "DAREA", "DLOAD", "RLOAD", "TLOAD", "SPOINT", "ASET",   "OMIT",  "INCLUDE",
    };
    for (const char* p : prefixes)
        if (rKeyword.rfind(p, 0) == 0)
            return true;
    return false;
}

std::string nas_strip(const std::string& rS) {
    const std::size_t b = rS.find_first_not_of(" \t");
    if (b == std::string::npos)
        return "";
    const std::size_t e = rS.find_last_not_of(" \t");
    return rS.substr(b, e - b + 1);
}

bool nas_is_comment(const std::string& rLine) {
    return rLine.size() < 3 || rLine[0] == '$' || rLine[0] == '#' || rLine.rfind("//", 0) == 0;
}

// One raw field of a card line; `mNone` marks a continuation marker, which is
// dropped when the card's fields are flattened.
struct NasChunk {
    std::string mText;
    bool mNone = false;
};
using NasChunks = std::vector<NasChunk>;

bool nas_is_free(const std::string& rLine) {
    return rLine.find(',') != std::string::npos;
}

// A free-field line splits on commas; a fixed-field line into (at most ten)
// 8-column fields, the tenth being the continuation marker.
NasChunks nas_chunk_line(const std::string& rLine) {
    NasChunks out;
    if (nas_is_free(rLine)) {
        std::size_t start = 0;
        while (true) {
            const std::size_t comma = rLine.find(',', start);
            if (comma == std::string::npos) {
                out.push_back({rLine.substr(start)});
                break;
            }
            out.push_back({rLine.substr(start, comma - start)});
            start = comma + 1;
        }
        return out;
    }
    for (std::size_t i = 0; i < rLine.size() && out.size() < 10; i += 8)
        out.push_back({rLine.substr(i, 8)});
    return out;
}

// Large-field lines hold 8 + 4x16 + 8 columns: re-merge each pair of 8-column
// chunks into one 16-column field.
NasChunks nas_merge_large(const NasChunks& rC) {
    NasChunks d;
    d.push_back(rC[0]);
    for (std::size_t k = 1; k <= 7 && k < rC.size(); k += 2) {
        NasChunk f = rC[k];
        if (k + 1 < rC.size() && !rC[k + 1].mNone)
            f.mText += rC[k + 1].mText;
        d.push_back(f);
    }
    if (rC.size() > 9)
        d.push_back(rC[9]);
    return d;
}

// The logical cards of the bulk section: each is the flattened, stripped list
// of its fields, continuation lines merged in.
std::vector<std::vector<std::string>> nas_cards(const std::vector<std::string>& rLines) {
    std::vector<std::vector<std::string>> cards;
    const std::string blank8(8, ' ');
    std::size_t i = 0;
    while (i < rLines.size()) {
        std::vector<NasChunks> chunks;
        chunks.push_back(nas_chunk_line(rLines[i]));
        const bool free = nas_is_free(rLines[i]);
        ++i;
        while (i < rLines.size()) {
            const std::string& next = rLines[i];
            if (next[0] == '+' || next[0] == '*') {
                if (chunks.back().size() == 10)
                    chunks.back().back().mNone = true;
                NasChunks c = nas_chunk_line(next);
                c[0].mNone = true;
                chunks.push_back(std::move(c));
                ++i;
            } else if (chunks.back().size() == 10 && !chunks.back().back().mNone &&
                       chunks.back().back().mText == blank8) {
                // Implicit continuation: a blank tenth field followed by a line
                // whose first field is blank too.
                NasChunks c = nas_chunk_line(next);
                if (!c.empty() && c[0].mText == blank8) {
                    chunks.back()[9].mNone = true;
                    c[0].mNone = true;
                    chunks.push_back(std::move(c));
                    ++i;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        const std::string head = nas_strip(chunks[0][0].mText);
        if (!free && !head.empty() && head.back() == '*')
            for (NasChunks& c : chunks)
                c = nas_merge_large(c);
        std::vector<std::string> fields;
        for (const NasChunks& c : chunks)
            for (const NasChunk& f : c)
                if (!f.mNone)
                    fields.push_back(nas_strip(f.mText));
        cards.push_back(std::move(fields));
    }
    return cards;
}

const std::string& nas_field(const std::vector<std::string>& rFields, std::size_t k) {
    static const std::string empty;
    return k < rFields.size() ? rFields[k] : empty;
}

std::int64_t nas_int(const std::string& rText, const std::string& rCard) {
    return detail::card_to_int(rText, " in a " + rCard + " card", "Nastran");
}

double nas_real(const std::string& rText, const std::string& rCard) {
    return detail::card_to_real(rText, " in a " + rCard + " card", "Nastran");
}

// A list of ids with `a THRU b` ranges, as $HMMOVE lines and SET cards hold them.
// Explicit ids land in `rIds`; each range in `rRanges`.
void nas_parse_id_list(const std::vector<std::string>& rTokens, std::size_t First,
                       const std::string& rCard, std::vector<std::int64_t>& rIds,
                       std::vector<std::pair<std::int64_t, std::int64_t>>& rRanges) {
    std::vector<std::string> t;
    for (std::size_t k = First; k < rTokens.size(); ++k)
        if (!rTokens[k].empty())
            t.push_back(rTokens[k]);
    for (std::size_t k = 0; k < t.size(); ++k) {
        if (k + 2 < t.size() && t[k + 1] == "THRU") {
            rRanges.emplace_back(nas_int(t[k], rCard), nas_int(t[k + 2], rCard));
            k += 2;
        } else {
            rIds.push_back(nas_int(t[k], rCard));
        }
    }
}

// Members of a HyperMesh component or a SET, as parsed.
struct NasGroup {
    std::vector<std::int64_t> mIds;
    std::vector<std::pair<std::int64_t, std::int64_t>> mRanges;
};

// The HyperMesh comment cards: component membership ($HMMOVE plus the `$`
// id lines after it), component names ($HMNAME COMP) and SET names ($HMSET).
struct NasHyperMesh {
    std::map<std::int64_t, NasGroup> mComponents;
    std::map<std::int64_t, std::string> mComponentNames;
    // Component id -> the property id `$HMNAME COMP <id>"name" <pid> "type"`
    // records after the name, when it does.
    std::map<std::int64_t, std::int64_t> mComponentProperties;
    std::map<std::int64_t, std::string> mSetNames;
    std::int64_t mActive = -1;
    bool mHasActive = false;

    // The quoted name after position `Pos`, or nothing when there is none.
    static bool QuotedAfter(const std::string& rLine, std::size_t Pos, std::int64_t& rId,
                            std::string& rName, std::size_t* pEnd = nullptr) {
        std::size_t k = Pos;
        while (k < rLine.size() && rLine[k] == ' ')
            ++k;
        const std::size_t digits = k;
        while (k < rLine.size() && rLine[k] >= '0' && rLine[k] <= '9')
            ++k;
        if (k == digits)
            return false;
        rId = std::strtoll(rLine.substr(digits, k - digits).c_str(), nullptr, 10);
        const std::size_t open = rLine.find('"', k);
        if (open == std::string::npos)
            return false;
        const std::size_t close = rLine.find('"', open + 1);
        if (close == std::string::npos)
            return false;
        rName = rLine.substr(open + 1, close - open - 1);
        if (pEnd)
            *pEnd = close + 1;
        return true;
    }

    // `$` followed by blanks, then 8-column fields of ids and THRU.
    static bool IdLine(const std::string& rLine, std::vector<std::string>& rTokens) {
        if (rLine.empty() || rLine[0] != '$')
            return false;
        for (std::size_t k = 1; k < 8 && k < rLine.size(); ++k)
            if (rLine[k] != ' ')
                return false;
        bool any = false;
        for (std::size_t k = 8; k < rLine.size(); k += 8) {
            const std::string f = nas_strip(rLine.substr(k, 8));
            if (f.empty())
                continue;
            if (f != "THRU" && f.find_first_not_of("0123456789") != std::string::npos)
                return false;
            any = any || f != "THRU";
            rTokens.push_back(f);
        }
        return any;
    }

    void Feed(const std::string& rLine) {
        if (rLine.rfind("$HMMOVE", 0) == 0) {
            const std::string id = nas_strip(rLine.substr(7, 9));
            mHasActive = !id.empty() && id.find_first_not_of("0123456789") == std::string::npos;
            if (mHasActive) {
                mActive = std::strtoll(id.c_str(), nullptr, 10);
                mComponents[mActive];
            }
            return;
        }
        std::vector<std::string> tokens;
        if (mHasActive && IdLine(rLine, tokens)) {
            NasGroup& g = mComponents[mActive];
            nas_parse_id_list(tokens, 0, "$HMMOVE", g.mIds, g.mRanges);
            return;
        }
        mHasActive = false;
        std::int64_t id = 0;
        std::string name;
        std::size_t end = 0;
        if (rLine.rfind("$HMNAME COMP ", 0) == 0 && QuotedAfter(rLine, 13, id, name, &end)) {
            mComponentNames[id] = name;
            std::size_t k = end;
            while (k < rLine.size() && rLine[k] == ' ')
                ++k;
            std::size_t e = k;
            while (e < rLine.size() && rLine[e] >= '0' && rLine[e] <= '9')
                ++e;
            if (e > k)
                mComponentProperties[id] =
                    std::strtoll(rLine.substr(k, e - k).c_str(), nullptr, 10);
        } else if (rLine.rfind("$HMSET ", 0) == 0) {
            // $HMSET <id> <type> "name"
            std::size_t k = 7;
            while (k < rLine.size() && rLine[k] == ' ')
                ++k;
            std::size_t e = k;
            while (e < rLine.size() && rLine[e] >= '0' && rLine[e] <= '9')
                ++e;
            if (e > k && QuotedAfter(rLine, e, id, name))
                mSetNames[std::strtoll(rLine.substr(k, e - k).c_str(), nullptr, 10)] = name;
        }
    }
};

// The global indices a group names. Ids no entity has are counted in
// `rMissing` when listed explicitly; a range only picks the ids that exist.
std::vector<std::int64_t> nas_resolve(const NasGroup& rGroup,
                                      const std::unordered_map<std::int64_t, std::int64_t>& rIndex,
                                      const std::map<std::int64_t, std::int64_t>& rSorted,
                                      std::size_t& rMissing) {
    std::vector<std::int64_t> out;
    for (std::int64_t id : rGroup.mIds) {
        auto it = rIndex.find(id);
        if (it == rIndex.end())
            ++rMissing;
        else
            out.push_back(it->second);
    }
    for (const auto& r : rGroup.mRanges)
        for (auto it = rSorted.lower_bound(r.first); it != rSorted.end() && it->first <= r.second;
             ++it)
            out.push_back(it->second);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::string nastran_float(double v) {
    if (v == 0.0)
        return "0.0";
    char buf[40];
    std::string best;
    for (int p = 0; p <= 11; ++p) {
        detail::snprintf_c(buf, sizeof(buf), "%.*E", p, v);
        const char* end = nullptr;
        if (detail::parse_double(buf, end) == v) {
            best = buf;
            break;
        }
    }
    if (best.empty()) {
        detail::snprintf_c(buf, sizeof(buf), "%.11E", v);
        best = buf;
    }
    std::size_t epos = best.find('E');
    std::string mant = best.substr(0, epos);
    int exp = std::atoi(best.c_str() + epos + 1);
    std::size_t dot = mant.find('.');
    if (dot == std::string::npos) {
        mant += ".";
        dot = mant.size() - 1;
    }
    // trim trailing zeros after the decimal point (keep the dot)
    std::size_t last = mant.size();
    while (last > dot + 1 && mant[last - 1] == '0')
        --last;
    mant.erase(last);
    std::string es = (exp < 0 ? "-" : "+") + std::to_string(std::abs(exp));
    std::string out = mant + "E" + es;
    // Keep within the 16-char field by shedding mantissa precision if needed.
    while (out.size() > 16 && mant.find('.') != std::string::npos && mant.back() != '.') {
        mant.pop_back();
        out = mant + "E" + es;
    }
    return out;
}

// The HyperMesh comment block that records disjoint cell regions as
// components: `$HMMOVE` with the element ids (runs as `a THRU b`), then
// `$HMNAME COMP`. Comments only, so no solver sees them. The Python writer
// emits the same bytes (nastran/_nastran.py, `_hypermesh_block`).
std::string nas_hypermesh_block(const Mesh& rMesh) {
    std::vector<const Region*> kept;
    std::vector<char> owned(
        static_cast<std::size_t>(detail::total_cells(detail::block_bases(rMesh))), 0);
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind != RegionKind::Cell) {
            log::warn("Nastran: {} region '{}' dropped; HyperMesh components hold cells only",
                      reg.mKind == RegionKind::Point ? "point" : "side", reg.mName);
            continue;
        }
        const std::int64_t* e = reg.Entries();
        for (std::size_t j = 0; j < reg.NumEntries(); ++j)
            if (e[j] < 0 || static_cast<std::size_t>(e[j]) >= owned.size())
                throw WriteError("Nastran writer: cell region '" + reg.mName + "' names cell " +
                                 std::to_string(e[j]) + " of " + std::to_string(owned.size()));
        bool overlaps = false;
        for (std::size_t j = 0; j < reg.NumEntries() && !overlaps; ++j)
            overlaps = owned[static_cast<std::size_t>(e[j])] != 0;
        if (overlaps) {
            log::warn(
                "Nastran: cell region '{}' overlaps an earlier one and is dropped; a "
                "HyperMesh component owns each element once",
                reg.mName);
            continue;
        }
        for (std::size_t j = 0; j < reg.NumEntries(); ++j)
            owned[static_cast<std::size_t>(e[j])] = 1;
        kept.push_back(&reg);
    }
    if (kept.empty())
        return "";
    // Keep the regions' tags as component ids when they are usable ones.
    bool tags_ok = true;
    std::vector<std::int64_t> seen;
    for (const Region* reg : kept) {
        if (reg->mTag <= 0 || std::find(seen.begin(), seen.end(), reg->mTag) != seen.end())
            tags_ok = false;
        seen.push_back(reg->mTag);
    }
    std::string out;
    char buf[64];
    for (std::size_t k = 0; k < kept.size(); ++k) {
        const Region& reg = *kept[k];
        const std::int64_t comp = tags_ok ? reg.mTag : static_cast<std::int64_t>(k + 1);
        std::string name = reg.mName;
        std::replace(name.begin(), name.end(), '"', '\'');
        out += "$\n";
        if (reg.NumEntries() > 0) {
            std::snprintf(buf, sizeof(buf), "$HMMOVE %8lld\n", static_cast<long long>(comp));
            out += buf;
            const std::int64_t* e = reg.Entries();
            std::string singles;
            int nsingles = 0;
            auto flush = [&]() {
                if (nsingles > 0)
                    out += "$       " + singles + "\n";
                singles.clear();
                nsingles = 0;
            };
            std::size_t j = 0;
            while (j < reg.NumEntries()) {
                std::size_t run = j;
                while (run + 1 < reg.NumEntries() && e[run + 1] == e[run] + 1)
                    ++run;
                if (run > j) {
                    flush();
                    std::snprintf(buf, sizeof(buf), "$       %8lldTHRU    %8lld\n",
                                  static_cast<long long>(e[j] + 1),
                                  static_cast<long long>(e[run] + 1));
                    out += buf;
                } else {
                    std::snprintf(buf, sizeof(buf), "%8lld", static_cast<long long>(e[j] + 1));
                    singles += buf;
                    if (++nsingles == 8)
                        flush();
                }
                j = run + 1;
            }
            flush();
        }
        std::snprintf(buf, sizeof(buf), "$HMNAME COMP%20lld", static_cast<long long>(comp));
        out += buf;
        out += "\"" + name + "\"\n";
    }
    return out;
}

}  // namespace

void write_nastran(const std::string& rPath, const Mesh& rMesh) {
    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    const NDArray* point_refs =
        rMesh.HasPointData("nastran:ref") ? &rMesh.PointData("nastran:ref") : nullptr;
    const bool cell_refs = rMesh.HasCellData("nastran:ref");

    // Validate before writing anything.
    const auto& m2n = meshio_to_nastran();
    for (const auto cb : rMesh.CellRange())
        if (!m2n.count(cb.Type()))
            throw WriteError("Nastran writer: unsupported cell type " + cb.Type());

    os << "$ " << kSentinel << "\n";
    os << detail::provenance_render_lines(detail::SlotTier::Block, "$ ");
    os << "BEGIN BULK\n";

    // Points: fixed-large GRID*. A zero reference field is written blank.
    char buf[128];
    for (std::size_t i = 0; i < n; ++i) {
        double xyz[3] = {0, 0, 0};
        for (std::size_t c = 0; c < dim && c < 3; ++c)
            xyz[c] = detail::read_double(points, i * dim + c);
        std::string sx = nastran_float(xyz[0]), sy = nastran_float(xyz[1]),
                    sz = nastran_float(xyz[2]);
        std::string ref;
        if (point_refs) {
            const long long v = detail::read_int(*point_refs, i);
            if (v != 0)
                ref = std::to_string(v);
        }
        std::snprintf(buf, sizeof(buf), "GRID*   %-16d%-16s%16s%16s\n*       %16s\n",
                      static_cast<int>(i + 1), ref.c_str(), sx.c_str(), sy.c_str(), sz.c_str());
        os << buf;
    }

    // Cells: fixed-small element cards (8-char fields), with + continuations.
    std::size_t cell_id = 0;
    std::size_t block = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& ntype = m2n.at(cb.Type());
        const NDArray* refs = cell_refs ? &rMesh.CellData("nastran:ref", block) : nullptr;
        ++block;
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const std::vector<int>* perm = nas_write_perm(cb.Type());
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            ++cell_id;
            std::vector<long long> nodes(k);
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t src = perm ? static_cast<std::size_t>((*perm)[j]) : j;
                nodes[j] = detail::read_int(conn, r * k + src) + 1;
            }
            std::string ref;
            if (refs) {
                const long long v = detail::read_int(*refs, r);
                if (v != 0)
                    ref = std::to_string(v);
            }
            // first line: type, id, ref, up to 6 nodes
            std::snprintf(buf, sizeof(buf), "%-8s%-8d%-8s", ntype.c_str(),
                          static_cast<int>(cell_id), ref.c_str());
            std::string line = buf;
            std::size_t nipl1 = 6, nipl2 = 14;
            for (std::size_t j = 0; j < k && j < nipl1; ++j) {
                std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                line += buf;
            }
            if (k > nipl1) {
                std::snprintf(buf, sizeof(buf), "+1%-6x", static_cast<unsigned>(cell_id));
                os << line << buf << "\n";
                std::snprintf(buf, sizeof(buf), "+1%-6x", static_cast<unsigned>(cell_id));
                std::string l2 = buf;
                for (std::size_t j = nipl1; j < k && j < nipl2; ++j) {
                    std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                    l2 += buf;
                }
                if (k > nipl2) {
                    std::snprintf(buf, sizeof(buf), "+2%-6x", static_cast<unsigned>(cell_id));
                    os << l2 << buf << "\n";
                    std::snprintf(buf, sizeof(buf), "+2%-6x", static_cast<unsigned>(cell_id));
                    std::string l3 = buf;
                    for (std::size_t j = nipl2; j < k; ++j) {
                        std::snprintf(buf, sizeof(buf), "%-8lld", nodes[j]);
                        l3 += buf;
                    }
                    os << l3 << "\n";
                } else {
                    os << l2 << "\n";
                }
            } else {
                os << line << "\n";
            }
        }
    }

    os << nas_hypermesh_block(rMesh);
    os << "ENDDATA\n";
    if (!os)
        throw WriteError("Nastran writer: failed writing " + rPath);
}

namespace {

// The bulk reader proper; also hands back the GRID ids (point order) and the
// element ids (global cell order) when asked.
Mesh nas_read(const std::string& rPath, std::vector<std::int64_t>* pGridIds,
              std::vector<std::int64_t>* pCellIds) {
    auto in = detail::make_classic_ifstream(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);

    // Everything before BEGIN BULK (executive and case control, I/O options) is
    // skipped; comment lines feed the HyperMesh parser; ENDDATA ends the deck.
    std::vector<std::string> lines;
    NasHyperMesh hm;
    std::string l;
    bool bulk = false;
    while (std::getline(in, l)) {
        if (!l.empty() && l.back() == '\r')
            l.pop_back();
        if (!bulk) {
            bulk = nas_strip(l).rfind("BEGIN BULK", 0) == 0;
            continue;
        }
        if (l.rfind("ENDDATA", 0) == 0)
            break;
        if (!l.empty() && l[0] == '$')
            hm.Feed(l);
        else
            hm.mHasActive = false;
        if (!nas_is_comment(l))
            lines.push_back(l);
    }
    if (!bulk)
        throw ReadError("Nastran: \"BEGIN BULK\" statement not found in " + rPath);

    struct Blk {
        std::string mType;
        int mN = 0;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mRefs;
        std::vector<std::int64_t> mIds;
    };
    std::vector<Blk> blocks;
    std::vector<double> pts;
    std::vector<std::int64_t> point_refs;
    bool any_point_ref = false, any_cell_ref = false;
    std::unordered_map<std::int64_t, std::int64_t> point_index, cell_index;
    std::vector<std::int64_t> cell_dims;
    std::vector<std::int64_t> cell_pids;
    std::map<std::int64_t, std::pair<std::string, NasGroup>> sets;
    std::map<std::string, std::size_t> skipped;
    const auto& elements = nas_elements();

    for (const std::vector<std::string>& f : nas_cards(lines)) {
        std::string kw = f.empty() ? std::string() : f[0];
        if (!kw.empty() && kw.back() == '*')
            kw.pop_back();
        if (kw == "GRID") {
            const std::int64_t id = nas_int(nas_field(f, 1), kw);
            const std::string& ref = nas_field(f, 2);
            any_point_ref = any_point_ref || !ref.empty();
            point_refs.push_back(nas_int(ref, kw));
            point_index[id] = static_cast<std::int64_t>(point_refs.size() - 1);
            if (pGridIds)
                pGridIds->push_back(id);
            for (std::size_t c = 3; c < 6; ++c)
                pts.push_back(nas_real(nas_field(f, c), kw));
            continue;
        }
        auto el = elements.find(kw);
        if (el != elements.end()) {
            const std::int64_t id = nas_int(nas_field(f, 1), kw);
            const std::string& ref = nas_field(f, 2);
            std::string type = el->second.mType;
            std::vector<std::int64_t> nodes;
            if (el->second.mNodes > 0) {
                for (int j = 0; j < el->second.mNodes; ++j) {
                    const std::string& t = nas_field(f, 3 + static_cast<std::size_t>(j));
                    if (t.empty())
                        throw ReadError("Nastran: " + kw + " " + std::to_string(id) +
                                        " is missing node " + std::to_string(j + 1));
                    nodes.push_back(nas_int(t, kw));
                }
            } else {
                for (std::size_t j = 3; j < f.size(); ++j)
                    if (!f[j].empty())
                        nodes.push_back(nas_int(f[j], kw));
                const NasSolid* s = nas_solid(type);
                const int nn = static_cast<int>(nodes.size());
                if (nn == s->mQuadraticNodes)
                    type = s->mQuadratic;
                else if (nn != s->mLinearNodes)
                    throw ReadError("Nastran: " + kw + " " + std::to_string(id) + " has " +
                                    std::to_string(nn) + " nodes; expected " +
                                    std::to_string(s->mLinearNodes) + " or " +
                                    std::to_string(s->mQuadraticNodes));
            }
            if (const std::vector<int>* perm = nas_read_perm(kw, type)) {
                std::vector<std::int64_t> p(nodes.size());
                for (std::size_t j = 0; j < p.size(); ++j)
                    p[j] = nodes[static_cast<std::size_t>((*perm)[j])];
                nodes = std::move(p);
            }
            if (blocks.empty() || blocks.back().mType != type) {
                Blk b;
                b.mType = type;
                b.mN = static_cast<int>(nodes.size());
                blocks.push_back(std::move(b));
            }
            Blk& blk = blocks.back();
            blk.mConn.insert(blk.mConn.end(), nodes.begin(), nodes.end());
            any_cell_ref = any_cell_ref || !ref.empty();
            blk.mRefs.push_back(nas_int(ref, kw));
            blk.mIds.push_back(id);
            cell_pids.push_back(blk.mRefs.back());
            cell_index[id] = static_cast<std::int64_t>(cell_dims.size());
            cell_dims.push_back(cell_type_dimension(cell_type_from_name(type)));
            continue;
        }
        if (kw == "SET") {
            // OptiStruct: SET, id, GRID|ELEM, LIST, ids (with THRU ranges).
            const std::int64_t id = nas_int(nas_field(f, 1), kw);
            const std::string& kind = nas_field(f, 2);
            const std::string& sub = nas_field(f, 3);
            if ((kind != "GRID" && kind != "ELEM") || sub != "LIST") {
                log::warn("Nastran: SET {} of type '{} {}' is not read; skipped", id, kind, sub);
                continue;
            }
            auto& entry = sets[id];
            entry.first = kind;
            nas_parse_id_list(f, 4, kw, entry.second.mIds, entry.second.mRanges);
            continue;
        }
        if (!nas_is_quietly_skipped(kw))
            ++skipped[kw];
    }

    if (!skipped.empty()) {
        std::string list;
        std::size_t total = 0;
        for (const auto& [k, c] : skipped) {
            list += (list.empty() ? "" : ", ") + k + " (" + std::to_string(c) + ")";
            total += c;
        }
        log::warn("Nastran: skipped {} card(s) meshio++ does not read: {}", total, list);
    }

    Mesh mesh;
    const std::size_t np = point_refs.size();
    NDArray points(DType::Float64, {np, 3});
    std::copy(pts.begin(), pts.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    std::vector<NDArray> cell_ref_arrays;
    for (Blk& blk : blocks) {
        const std::size_t count = blk.mRefs.size();
        NDArray data(DType::Int64, {count, static_cast<std::size_t>(blk.mN)});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t idx = 0; idx < blk.mConn.size(); ++idx) {
            auto it = point_index.find(blk.mConn[idx]);
            if (it == point_index.end())
                throw ReadError("Nastran: an element references grid " +
                                std::to_string(blk.mConn[idx]) + ", which is not defined");
            dp[idx] = it->second;
        }
        mesh.AddCellBlock(blk.mType, std::move(data));
        NDArray refs(DType::Int64, {count});
        std::copy(blk.mRefs.begin(), blk.mRefs.end(), refs.As<std::int64_t>());
        cell_ref_arrays.push_back(std::move(refs));
    }
    if (any_point_ref) {
        NDArray refs(DType::Int64, {np});
        std::copy(point_refs.begin(), point_refs.end(), refs.As<std::int64_t>());
        mesh.AddPointData("nastran:ref", std::move(refs));
    }
    if (any_cell_ref)
        mesh.AddCellData("nastran:ref", std::move(cell_ref_arrays));

    // Regions: HyperMesh components (tag = component id), then SET cards
    // (tag = set id). Names resolve last: $HMNAME may follow the elements.
    const std::map<std::int64_t, std::int64_t> sorted_cells(cell_index.begin(), cell_index.end());
    const std::map<std::int64_t, std::int64_t> sorted_points(point_index.begin(),
                                                             point_index.end());
    std::size_t missing = 0;
    auto cell_dim = [&](const std::vector<std::int64_t>& rIds) -> int {
        int dim = -1;
        for (std::size_t k = 0; k < rIds.size(); ++k) {
            const int d = static_cast<int>(cell_dims[static_cast<std::size_t>(rIds[k])]);
            if (k == 0)
                dim = d;
            else if (d != dim)
                return -1;
        }
        return dim;
    };
    auto add = [&](const std::string& rName, RegionKind Kind, int Dim, std::int64_t Tag,
                   const std::vector<std::int64_t>& rIds) {
        NDArray entries(DType::Int64, {rIds.size()});
        std::copy(rIds.begin(), rIds.end(), entries.As<std::int64_t>());
        mesh.AddRegion(Region(rName, Kind, Dim, Tag, std::move(entries)));
    };
    for (const auto& [comp, _] : hm.mComponentNames)
        hm.mComponents[comp];
    std::map<std::int64_t, std::vector<std::int64_t>> members;
    std::vector<char> moved(cell_dims.size(), 0);
    for (const auto& [comp, group] : hm.mComponents) {
        members[comp] = nas_resolve(group, cell_index, sorted_cells, missing);
        for (std::int64_t c : members[comp])
            moved[static_cast<std::size_t>(c)] = 1;
    }
    // An element no $HMMOVE lists belongs to the component whose recorded
    // property is its PID (HyperMesh writes $HMMOVE only for the others).
    if (!hm.mComponentProperties.empty()) {
        std::unordered_map<std::int64_t, std::int64_t> by_property;
        for (const auto& [comp, pid] : hm.mComponentProperties)
            by_property.emplace(pid, comp);
        bool added = false;
        for (std::size_t c = 0; c < cell_pids.size(); ++c) {
            if (moved[c])
                continue;
            auto it = by_property.find(cell_pids[c]);
            if (it != by_property.end()) {
                members[it->second].push_back(static_cast<std::int64_t>(c));
                added = true;
            }
        }
        if (added)
            for (auto& [comp, ids] : members)
                std::sort(ids.begin(), ids.end());
    }
    for (const auto& [comp, ids] : members) {
        const auto name_it = hm.mComponentNames.find(comp);
        const std::string name = name_it != hm.mComponentNames.end()
                                     ? name_it->second
                                     : "component_" + std::to_string(comp);
        add(name, RegionKind::Cell, cell_dim(ids), comp, ids);
    }
    for (const auto& [id, entry] : sets) {
        const auto name_it = hm.mSetNames.find(id);
        const std::string name =
            name_it != hm.mSetNames.end() ? name_it->second : "set_" + std::to_string(id);
        if (entry.first == "GRID") {
            const std::vector<std::int64_t> ids =
                nas_resolve(entry.second, point_index, sorted_points, missing);
            add(name, RegionKind::Point, -1, id, ids);
        } else {
            const std::vector<std::int64_t> ids =
                nas_resolve(entry.second, cell_index, sorted_cells, missing);
            add(name, RegionKind::Cell, cell_dim(ids), id, ids);
        }
    }
    if (missing > 0)
        log::warn("Nastran: {} region member id(s) name no grid or element; dropped", missing);
    if (pCellIds)
        for (const Blk& blk : blocks)
            pCellIds->insert(pCellIds->end(), blk.mIds.begin(), blk.mIds.end());
    return mesh;
}

}  // namespace

Mesh read_nastran(const std::string& rPath) {
    return nas_read(rPath, nullptr, nullptr);
}

namespace detail {

Mesh nastran_read_deck(const std::string& rPath, std::vector<std::int64_t>& rGridIds,
                       std::vector<std::int64_t>& rCellIds) {
    rGridIds.clear();
    rCellIds.clear();
    return nas_read(rPath, &rGridIds, &rCellIds);
}

}  // namespace detail

}  // namespace meshioplusplus
