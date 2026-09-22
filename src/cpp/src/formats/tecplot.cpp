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
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

namespace meshioplusplus {

namespace {

std::string tecplot_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string tecplot_strip(const std::string& rS) {
    std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}
std::vector<std::string> tecplot_tokens(const std::string& rS) {
    std::vector<std::string> out;
    auto iss = detail::make_classic_istringstream(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}
bool is_float_token(const std::string& rS) {
    if (rS.empty())
        return false;
    const char* endp = nullptr;
    detail::parse_double(rS.c_str(), endp);
    return endp == rS.c_str() + rS.size();
}

/// Splits a joined header (e.g. `ZONE T = "VARLOCATION" N=4 ...`) into tokens,
/// keeping a quoted `"..."` span or a parenthesized `(...)` span as one
/// token (quotes/parens included) rather than letting their commas, `=` or
/// spaces fall through to the general splitter -- a title can legitimately
/// contain any of those (a real adversarial fixture: `T = "VARLOCATION"`,
/// where a naive substring search for the VARLOCATION *field* would instead
/// find this title). `=` becomes its own one-character token; everything
/// else splits on whitespace/commas exactly like tecplot_tokens.
std::vector<std::string> tecplot_header_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::size_t i = 0;
    const std::size_t n = rS.size();
    while (i < n) {
        const char c = rS[i];
        if (c == ' ' || c == '\t' || c == ',') {
            ++i;
        } else if (c == '=') {
            out.emplace_back("=");
            ++i;
        } else if (c == '"') {
            std::size_t j = rS.find('"', i + 1);
            if (j == std::string::npos)
                j = n - 1;
            out.push_back(rS.substr(i, j - i + 1));
            i = j + 1;
        } else if (c == '(') {
            std::size_t j = rS.find(')', i + 1);
            if (j == std::string::npos)
                j = n - 1;
            out.push_back(rS.substr(i, j - i + 1));
            i = j + 1;
        } else {
            std::size_t j = i;
            while (j < n && rS[j] != ' ' && rS[j] != '\t' && rS[j] != ',' && rS[j] != '=' &&
                  rS[j] != '"' && rS[j] != '(')
                ++j;
            out.push_back(rS.substr(i, j - i));
            i = j;
        }
    }
    return out;
}

/// Strips the surrounding quotes of a token tecplot_header_tokens quoted
/// (a no-op on a token that was never quoted).
std::string tecplot_unquote(const std::string& rTok) {
    if (rTok.size() >= 2 && rTok.front() == '"' && rTok.back() == '"')
        return rTok.substr(1, rTok.size() - 2);
    return rTok;
}

/// Splits `rInner` on commas that are not inside a `[...]` range, so
/// `"[1,2]=1,[4-6]=2"` splits into `{"[1,2]=1", "[4-6]=2"}` rather than
/// breaking the first range apart. Used for VARLOCATION/VARSHARELIST's
/// `([range]=target, ...)` value and for PASSIVEVARLIST's `([range])`.
std::vector<std::string> tecplot_split_entries(const std::string& rInner) {
    std::vector<std::string> out;
    std::string cur;
    int depth = 0;
    for (char c : rInner) {
        if (c == '[')
            ++depth;
        else if (c == ']')
            --depth;
        if (c == ',' && depth == 0) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty())
        out.push_back(cur);
    return out;
}

/// Expands a `[a,b,c-d,...]` (or bare, bracket-less) range spec into 0-based
/// indices. Tecplot variable/zone numbers are 1-based in the file.
std::vector<std::size_t> tecplot_parse_ranges(std::string rSpec) {
    if (!rSpec.empty() && rSpec.front() == '[')
        rSpec = rSpec.substr(1);
    if (!rSpec.empty() && rSpec.back() == ']')
        rSpec.pop_back();
    std::vector<std::size_t> out;
    for (const std::string& part : tecplot_split_entries(rSpec)) {
        const std::size_t dash = part.find('-');
        if (dash == std::string::npos) {
            if (!part.empty())
                out.push_back(static_cast<std::size_t>(std::stoul(part) - 1));
        } else {
            const int a = std::stoi(part.substr(0, dash));
            const int b = std::stoi(part.substr(dash + 1));
            for (int k = a; k <= b; ++k)
                out.push_back(static_cast<std::size_t>(k - 1));
        }
    }
    return out;
}

std::string tecplot_to_meshio(const std::string& rZ) {
    std::string u = tecplot_upper(rZ);
    if (u == "LINESEG" || u == "FELINESEG")
        return "line";
    if (u == "TRIANGLE" || u == "FETRIANGLE")
        return "triangle";
    if (u == "QUADRILATERAL" || u == "FEQUADRILATERAL")
        return "quad";
    if (u == "TETRAHEDRON" || u == "FETETRAHEDRON")
        return "tetra";
    if (u == "BRICK" || u == "FEBRICK")
        return "hexahedron";
    return "";
}
std::string meshio_to_tecplot(const std::string& rM) {
    if (rM == "line")
        return "FELINESEG";
    if (rM == "triangle")
        return "FETRIANGLE";
    if (rM == "quad")
        return "FEQUADRILATERAL";
    if (rM == "tetra")
        return "FETETRAHEDRON";
    if (rM == "pyramid" || rM == "wedge" || rM == "hexahedron")
        return "FEBRICK";
    return "";
}
const std::vector<int>& tecplot_order(const std::string& rM) {
    static const std::map<std::string, std::vector<int>> o = {
        {"line", {0, 1}},
        {"triangle", {0, 1, 2}},
        {"quad", {0, 1, 2, 3}},
        {"tetra", {0, 1, 2, 3}},
        {"pyramid", {0, 1, 2, 3, 4, 4, 4, 4}},
        {"wedge", {0, 1, 4, 3, 2, 2, 5, 5}},
        {"hexahedron", {0, 1, 2, 3, 4, 5, 6, 7}},
    };
    static const std::vector<int> empty;
    auto it = o.find(rM);
    return it == o.end() ? empty : it->second;
}

}  // namespace

namespace {

// One ZONE header's parsed fields, its data section's line range, and its
// transient identity (SOLUTIONTIME/STRANDID). Shared by read_tecplot (which
// decodes exactly one zone's data) and read_tecplot_metadata (which decodes
// none): both start from the same tecplot_scan_zones pass, so which zone the
// data-decoding step picks and which zones the metadata's timeline lists can
// never drift against each other.
struct TecplotZoneHeader {
    std::map<std::string, std::string> mFields;  // NODES/N/ELEMENTS/E/DATAPACKING/ZONETYPE/F/ET/NV
    std::string mVarloc;
    std::string mTitle;  // T="...", unquoted; empty if not given
    bool mHasSolutionTime = false;
    double mSolutionTime = 0.0;
    bool mHasStrandId = false;
    int mStrandId = 0;
    std::size_t mDataStart = 0;
    std::size_t mNumNodes = 0;
    std::size_t mNumCells = 0;
    // Variable index (0-based) -> the 1-based zone number it is shared from
    // (VARSHARELIST=([m,n]=z)): that variable has no data of its own here.
    std::map<std::size_t, std::size_t> mVarShareZone;
    // Variable indices with no data anywhere for this zone (PASSIVEVARLIST).
    std::set<std::size_t> mPassiveVars;
    // 1-based zone number this zone's connectivity is shared from
    // (CONNECTIVITYSHAREZONE=z), or 0 for its own.
    std::size_t mConnShareZone = 0;
};

std::string tecplot_zone_field(const TecplotZoneHeader& rZ, const char* pA, const char* pB) {
    auto it = rZ.mFields.find(pA);
    if (it != rZ.mFields.end())
        return it->second;
    it = rZ.mFields.find(pB);
    return it != rZ.mFields.end() ? it->second : std::string();
}

// The zone's element type and, for FEBLOCK data, which variables are
// cell-centered -- the same derivation whether or not this zone's data is
// ever decoded.
void tecplot_zone_format(const TecplotZoneHeader& rZ, std::size_t NumVariables, bool& rFeblock,
                         std::string& rZtype, std::vector<int>& rCellCentered) {
    std::string fmt;
    if (rZ.mFields.count("F")) {
        fmt = tecplot_upper(rZ.mFields.at("F"));
        rZtype = rZ.mFields.count("ET") ? rZ.mFields.at("ET") : "";
    } else {
        fmt = "FE" + tecplot_upper(tecplot_zone_field(rZ, "DATAPACKING", ""));
        rZtype = tecplot_zone_field(rZ, "ZONETYPE", "");
    }
    rFeblock = (fmt == "FEBLOCK");

    rCellCentered.assign(NumVariables, 0);
    if (!rFeblock)
        return;
    if (rZ.mFields.count("NV")) {
        int nv = std::stoi(rZ.mFields.at("NV"));
        for (std::size_t k = static_cast<std::size_t>(nv); k < NumVariables; ++k)
            rCellCentered[k] = 1;
    } else if (!rZ.mVarloc.empty()) {
        std::string vc = rZ.mVarloc.substr(1, rZ.mVarloc.size() - 2);  // strip outer ()
        for (const std::string& entry : tecplot_split_entries(vc)) {
            const std::size_t eq = entry.find('=');
            if (eq == std::string::npos)
                continue;
            const std::string rng = entry.substr(0, eq);
            const std::string loc = tecplot_upper(entry.substr(eq + 1));
            if (loc != "CELLCENTERED")
                continue;
            for (std::size_t idx : tecplot_parse_ranges(rng))
                if (idx < NumVariables)
                    rCellCentered[idx] = 1;
        }
    }
}

// How many numeric tokens this zone's data block holds, in file order --
// FEBLOCK is one run per *owned, active* variable (cell-centered ones
// NumCells long, the rest NumNodes long); a variable this zone shares from an
// earlier one (VARSHARELIST) or carries nowhere (PASSIVEVARLIST) has no data
// here at all, so it contributes nothing to the budget. POINT/FEPOINT has no
// such sharing (it is NumNodes rows of NumVariables each, always).
std::size_t tecplot_zone_data_token_count(const TecplotZoneHeader& rZ, std::size_t NumVariables,
                                          bool Feblock, const std::vector<int>& rCellCentered) {
    if (!Feblock)
        return rZ.mNumNodes * NumVariables;
    std::size_t total = 0;
    for (std::size_t k = 0; k < NumVariables; ++k) {
        if (rZ.mVarShareZone.count(k) || rZ.mPassiveVars.count(k))
            continue;
        total += rCellCentered[k] ? rZ.mNumCells : rZ.mNumNodes;
    }
    return total;
}

// One pass over every line, splitting VARIABLES from every ZONE header (not
// just the first) and locating each zone's data section by actually counting
// off its own token budget plus its connectivity lines -- Tecplot ASCII has
// no fixed tokens-per-line convention, so this is the only reliable way to
// find where one zone's data ends and the next one's header begins.
std::vector<TecplotZoneHeader> tecplot_scan_zones(const std::vector<std::string>& rLines,
                                                  std::vector<std::string>& rVariables) {
    std::vector<TecplotZoneHeader> zones;
    std::size_t i = 0;
    for (; i < rLines.size(); ++i) {
        std::string u = tecplot_upper(rLines[i]);
        if (u.rfind("VARIABLES", 0) == 0) {
            std::string joined = rLines[i];
            while (i + 1 < rLines.size() && tecplot_strip(rLines[i + 1])[0] == '"')
                joined += " " + rLines[++i];
            std::string rhs = joined.substr(joined.find('=') + 1);
            std::size_t p = 0;
            while (p < rhs.size()) {
                if (rhs[p] == '"') {
                    std::size_t q = rhs.find('"', p + 1);
                    rVariables.push_back(rhs.substr(p + 1, q - p - 1));
                    p = q + 1;
                } else if (std::isspace((unsigned char)rhs[p]) || rhs[p] == ',') {
                    ++p;
                } else {
                    std::size_t q = p;
                    while (q < rhs.size() && !std::isspace((unsigned char)rhs[q]) && rhs[q] != ',')
                        ++q;
                    rVariables.push_back(rhs.substr(p, q - p));
                    p = q;
                }
            }
            continue;
        }
        if (u.rfind("ZONE", 0) != 0)
            continue;

        TecplotZoneHeader z;
        std::string joined = rLines[i];
        while (i + 1 < rLines.size() && !is_float_token(tecplot_tokens(rLines[i + 1])[0]))
            joined += " " + rLines[++i];
        z.mDataStart = i + 1;

        // tk[0] is the "ZONE" keyword itself; the rest is KEY = VALUE triples,
        // where a quoted or parenthesized VALUE is already one token (so a
        // title like T = "VARLOCATION" cannot be mistaken for the field of
        // the same name).
        const std::vector<std::string> tk = tecplot_header_tokens(joined);
        for (std::size_t k = 1; k + 2 < tk.size(); k += 3) {
            if (tk[k + 1] != "=") {
                ++k;  // resync: not a KEY = VALUE triple after all
                continue;
            }
            const std::string key = tecplot_upper(tk[k]);
            const std::string& val = tk[k + 2];
            if (key == "NODES" || key == "N" || key == "ELEMENTS" || key == "E" ||
                key == "DATAPACKING" || key == "ZONETYPE" || key == "F" || key == "ET" ||
                key == "NV") {
                z.mFields[key] = val;
            } else if (key == "T") {
                z.mTitle = tecplot_unquote(val);
            } else if (key == "SOLUTIONTIME") {
                z.mSolutionTime = detail::parse_double(val);
                z.mHasSolutionTime = true;
            } else if (key == "STRANDID") {
                z.mStrandId = std::stoi(val);
                z.mHasStrandId = true;
            } else if (key == "VARLOCATION") {
                z.mVarloc = val;
                z.mVarloc.erase(std::remove(z.mVarloc.begin(), z.mVarloc.end(), ' '),
                                z.mVarloc.end());
            } else if (key == "VARSHARELIST") {
                std::string inner = val.substr(1, val.size() - 2);  // strip outer ()
                for (const std::string& entry : tecplot_split_entries(inner)) {
                    const std::size_t eq = entry.find('=');
                    if (eq == std::string::npos)
                        continue;
                    const std::size_t zoneno =
                        static_cast<std::size_t>(std::stoul(entry.substr(eq + 1)));
                    for (std::size_t idx : tecplot_parse_ranges(entry.substr(0, eq)))
                        z.mVarShareZone[idx] = zoneno;
                }
            } else if (key == "PASSIVEVARLIST") {
                std::string inner = val.substr(1, val.size() - 2);  // strip outer ()
                for (std::size_t idx : tecplot_parse_ranges(inner))
                    z.mPassiveVars.insert(idx);
            } else if (key == "CONNECTIVITYSHAREZONE") {
                z.mConnShareZone = static_cast<std::size_t>(std::stoul(val));
            }
        }
        z.mNumNodes = std::stoull(tecplot_zone_field(z, "NODES", "N"));
        z.mNumCells = std::stoull(tecplot_zone_field(z, "ELEMENTS", "E"));

        bool feblock = false;
        std::string ztype;
        std::vector<int> cell_centered;
        tecplot_zone_format(z, rVariables.size(), feblock, ztype, cell_centered);
        const std::size_t want =
            tecplot_zone_data_token_count(z, rVariables.size(), feblock, cell_centered);

        std::size_t li = z.mDataStart, got = 0;
        while (got < want && li < rLines.size()) {
            got += tecplot_tokens(rLines[li]).size();
            ++li;
        }
        if (z.mConnShareZone == 0)
            li += z.mNumCells;  // one connectivity line per cell -- none if shared
        zones.push_back(z);
        i = li - 1;  // the for-loop's ++i resumes scanning right after
    }
    if (rVariables.empty())
        throw ReadError("Tecplot: no VARIABLES");
    if (zones.empty())
        throw ReadError("Tecplot: no ZONE");
    return zones;
}

// One entry per step, each the zone indices (file order preserved) that make
// up that step's mesh: several parts at once, exactly like read_tecplot's
// eventual Mesh. A non-transient file (no zone carries SOLUTIONTIME) is a
// single "step" listing every zone -- the several-static-parts case. A
// transient file's steps are its distinct SOLUTIONTIME values, ascending,
// each grouping every zone at that time sharing zones[0]'s STRANDID when any
// zone carries one (SOLUTIONTIME with no STRANDID at all groups every zone
// together); a zone with no SOLUTIONTIME in an otherwise transient file
// belongs to no step and is dropped.
std::vector<std::vector<std::size_t>> tecplot_timeline(
    const std::vector<TecplotZoneHeader>& rZones) {
    if (!rZones[0].mHasSolutionTime) {
        std::vector<std::size_t> all(rZones.size());
        for (std::size_t k = 0; k < rZones.size(); ++k)
            all[k] = k;
        return {std::move(all)};
    }
    std::map<double, std::vector<std::size_t>> by_time;
    for (std::size_t k = 0; k < rZones.size(); ++k) {
        if (!rZones[k].mHasSolutionTime)
            continue;
        if (rZones[0].mHasStrandId && rZones[k].mHasStrandId &&
            rZones[k].mStrandId != rZones[0].mStrandId)
            continue;
        by_time[rZones[k].mSolutionTime].push_back(k);
    }
    std::vector<std::vector<std::size_t>> steps;
    steps.reserve(by_time.size());
    for (auto& [time, idxs] : by_time)  // std::map is sorted by key (time)
        steps.push_back(std::move(idxs));
    return steps;
}

/// One zone's decoded data columns, connectivity and topology -- shared by
/// every step that references it (a zone can be VARSHARELIST'd or
/// CONNECTIVITYSHAREZONE'd from more than one later zone), so it is decoded
/// once and cached.
struct TecplotDecodedZone {
    std::vector<std::vector<double>> mCols;  // per-variable; length NumNodes or NumCells
    std::vector<int> mCellCentered;
    std::string mMeshioType;
    NDArray mConn;  // (NumCells, NodesPerCell) int64, 0-based
    std::size_t mNumNodes = 0;
    std::size_t mNumCells = 0;
};

/// Decodes zone `ZoneIdx`, resolving VARSHARELIST/CONNECTIVITYSHAREZONE by
/// recursively decoding (and caching) whichever earlier zone they name.
const TecplotDecodedZone& tecplot_decode_zone(std::size_t ZoneIdx,
                                              const std::vector<TecplotZoneHeader>& rZones,
                                              const std::vector<std::string>& rVariables,
                                              const std::vector<std::string>& rLines,
                                              std::map<std::size_t, TecplotDecodedZone>& rCache) {
    const auto found = rCache.find(ZoneIdx);
    if (found != rCache.end())
        return found->second;

    const TecplotZoneHeader& z = rZones[ZoneIdx];
    bool feblock = false;
    std::string ztype;
    std::vector<int> cell_centered;
    tecplot_zone_format(z, rVariables.size(), feblock, ztype, cell_centered);
    const std::string mtype = tecplot_to_meshio(ztype);
    if (mtype.empty())
        throw ReadError("Tecplot: unsupported zone type " + ztype);

    std::vector<std::size_t> ndata(rVariables.size());
    for (std::size_t k = 0; k < rVariables.size(); ++k)
        ndata[k] = cell_centered[k] ? z.mNumCells : z.mNumNodes;
    const std::size_t want =
        tecplot_zone_data_token_count(z, rVariables.size(), feblock, cell_centered);
    const std::size_t want_full = feblock ? want : z.mNumNodes * rVariables.size();

    std::vector<double> flat;
    flat.reserve(want_full);
    std::size_t li = z.mDataStart;
    while (flat.size() < want_full && li < rLines.size()) {
        for (const auto& t : tecplot_tokens(rLines[li]))
            flat.push_back(detail::parse_double(t));
        ++li;
    }

    std::vector<std::vector<double>> cols(rVariables.size());
    if (feblock) {
        std::size_t off = 0;
        for (std::size_t k = 0; k < rVariables.size(); ++k) {
            const auto share = z.mVarShareZone.find(k);
            if (share != z.mVarShareZone.end()) {
                const TecplotDecodedZone& src =
                    tecplot_decode_zone(share->second - 1, rZones, rVariables, rLines, rCache);
                cols[k] = src.mCols[k];
            } else if (z.mPassiveVars.count(k)) {
                cols[k].assign(ndata[k], std::numeric_limits<double>::quiet_NaN());
            } else {
                cols[k].assign(flat.begin() + static_cast<std::ptrdiff_t>(off),
                               flat.begin() + static_cast<std::ptrdiff_t>(off + ndata[k]));
                off += ndata[k];
            }
        }
    } else {
        const std::size_t nv = rVariables.size();
        for (std::size_t k = 0; k < nv; ++k)
            cols[k].resize(z.mNumNodes);
        for (std::size_t r = 0; r < z.mNumNodes; ++r)
            for (std::size_t k = 0; k < nv; ++k)
                cols[k][r] = flat[r * nv + k];
    }

    std::size_t nn;
    if (mtype == "line")
        nn = 2;
    else if (mtype == "triangle")
        nn = 3;
    else if (mtype == "quad" || mtype == "tetra")
        nn = 4;
    else
        nn = 8;

    NDArray conn;
    if (z.mConnShareZone != 0) {
        const TecplotDecodedZone& src =
            tecplot_decode_zone(z.mConnShareZone - 1, rZones, rVariables, rLines, rCache);
        conn = src.mConn;
    } else {
        conn = NDArray(DType::Int64, {z.mNumCells, nn});
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t c = 0; c < z.mNumCells; ++c) {
            auto t = tecplot_tokens(rLines.at(li++));
            for (std::size_t j = 0; j < nn; ++j)
                cp[c * nn + j] = std::strtoll(t[j].c_str(), nullptr, 10) - 1;
        }
    }

    TecplotDecodedZone result;
    result.mCols = std::move(cols);
    result.mCellCentered = std::move(cell_centered);
    result.mMeshioType = mtype;
    result.mConn = std::move(conn);
    result.mNumNodes = z.mNumNodes;
    result.mNumCells = z.mNumCells;
    const auto [inserted, _] = rCache.emplace(ZoneIdx, std::move(result));
    return inserted->second;
}

/// Resolves the zone whose data literally *is* `ZoneIdx`'s point set: itself,
/// unless X (and Y, and Z when 3-D) are all VARSHARELIST'd from the very same
/// earlier zone, in which case the chain is followed to that zone's own
/// origin. A zone that shares only some of its coordinate variables, or
/// shares them from different zones, is not literally the same point set --
/// it owns its own (possibly duplicated) points, same as no sharing at all.
/// Tecplot itself requires VARSHARELIST partners to agree on NODES, so this
/// only matters for zones that concatenate several element-type zones over
/// one shared point cloud (a common hybrid-mesh export shape), not for the
/// zones-as-a-timeline case, where sharing spans different steps.
std::size_t tecplot_points_origin_zone(std::size_t ZoneIdx,
                                       const std::vector<TecplotZoneHeader>& rZones, int Xi,
                                       int Yi, int Zi) {
    std::set<std::size_t> seen;
    std::size_t cur = ZoneIdx;
    while (seen.insert(cur).second) {
        const TecplotZoneHeader& z = rZones[cur];
        const auto it_x = z.mVarShareZone.find(static_cast<std::size_t>(Xi));
        if (it_x == z.mVarShareZone.end())
            return cur;
        const auto it_y = z.mVarShareZone.find(static_cast<std::size_t>(Yi));
        if (it_y == z.mVarShareZone.end() || it_y->second != it_x->second)
            return cur;
        if (Zi >= 0) {
            const auto it_z = z.mVarShareZone.find(static_cast<std::size_t>(Zi));
            if (it_z == z.mVarShareZone.end() || it_z->second != it_x->second)
                return cur;
        }
        cur = it_x->second - 1;
    }
    return ZoneIdx;  // cycle guard: fall back to owning its own points
}

/// The 0-based (X, Y, Z) variable indices, Z absent (-1) for a 2D file. The
/// same three variables for every zone in the file: VARIABLES is per-file.
void tecplot_xyz_indices(const std::vector<std::string>& rVariables, int& rXi, int& rYi, int& rZi) {
    rXi = rYi = rZi = -1;
    for (std::size_t k = 0; k < rVariables.size(); ++k) {
        const std::string v = tecplot_upper(rVariables[k]);
        if (v == "X")
            rXi = static_cast<int>(k);
        else if (v == "Y")
            rYi = static_cast<int>(k);
        else if (v == "Z")
            rZi = static_cast<int>(k);
    }
}

/// Builds one step's Mesh by concatenating every zone in `rZoneIdxs`: one
/// cell block per zone (never merged by type -- each zone is its own named
/// part), points offset per zone (zones are not welded), point/cell data
/// concatenated per variable (a PASSIVEVARLIST zone already decoded to NaN
/// for that variable, so every zone contributes a value), `tecplot:zone`
/// naming each cell's zone, and one Cell region per zone (named by its own
/// title when it has one, de-duplicated, else `zone_<i>`).
Mesh tecplot_build_step_mesh(const std::vector<std::size_t>& rZoneIdxs,
                             const std::vector<TecplotZoneHeader>& rZones,
                             const std::vector<std::string>& rVariables,
                             const std::vector<std::string>& rLines,
                             std::map<std::size_t, TecplotDecodedZone>& rCache) {
    int xi = -1, yi = -1, zi = -1;
    tecplot_xyz_indices(rVariables, xi, yi, zi);
    const std::size_t ndim = (zi >= 0) ? 3 : 2;

    std::vector<const TecplotDecodedZone*> decoded;
    decoded.reserve(rZoneIdxs.size());
    for (std::size_t idx : rZoneIdxs) {
        const TecplotDecodedZone& d =
            tecplot_decode_zone(idx, rZones, rVariables, rLines, rCache);
        decoded.push_back(&d);
    }

    // A zone whose X (and Y, Z) are all VARSHARELIST'd from an earlier zone
    // in this same step literally reuses that zone's points -- it gets no
    // new point block and its connectivity is offset the same as its
    // origin's, instead of being appended and offset as an independent part.
    std::map<std::size_t, std::size_t> idx_to_pos;
    for (std::size_t k = 0; k < rZoneIdxs.size(); ++k)
        idx_to_pos.emplace(rZoneIdxs[k], k);
    std::vector<std::size_t> point_offset(decoded.size());
    std::vector<bool> owns_points(decoded.size());
    std::size_t total_points = 0;
    for (std::size_t k = 0; k < decoded.size(); ++k) {
        const std::size_t origin = tecplot_points_origin_zone(rZoneIdxs[k], rZones, xi, yi, zi);
        const auto it = idx_to_pos.find(origin);
        if (origin != rZoneIdxs[k] && it != idx_to_pos.end() && it->second < k) {
            point_offset[k] = point_offset[it->second];
            owns_points[k] = false;
        } else {
            point_offset[k] = total_points;
            owns_points[k] = true;
            total_points += decoded[k]->mNumNodes;
        }
    }

    Mesh mesh;
    NDArray pts(DType::Float64, {total_points, ndim});
    double* pp = pts.As<double>();
    for (std::size_t k = 0; k < decoded.size(); ++k) {
        if (!owns_points[k])
            continue;
        const TecplotDecodedZone& d = *decoded[k];
        const std::size_t poff = point_offset[k];
        for (std::size_t r = 0; r < d.mNumNodes; ++r) {
            pp[(poff + r) * ndim + 0] = d.mCols[static_cast<std::size_t>(xi)][r];
            pp[(poff + r) * ndim + 1] = d.mCols[static_cast<std::size_t>(yi)][r];
            if (zi >= 0)
                pp[(poff + r) * ndim + 2] = d.mCols[static_cast<std::size_t>(zi)][r];
        }
    }
    mesh.AssignPoints(std::move(pts));

    for (std::size_t k = 0; k < decoded.size(); ++k) {
        const TecplotDecodedZone& d = *decoded[k];
        NDArray conn = d.mConn;  // deep copy: offset in place below
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t j = 0; j < conn.Size(); ++j)
            cp[j] += static_cast<std::int64_t>(point_offset[k]);
        mesh.AddCellBlock(d.mMeshioType, std::move(conn));
    }

    for (std::size_t k = 0; k < rVariables.size(); ++k) {
        if (static_cast<int>(k) == xi || static_cast<int>(k) == yi || static_cast<int>(k) == zi)
            continue;
        if (decoded[0]->mCellCentered[k]) {
            std::vector<NDArray> blk;
            blk.reserve(decoded.size());
            for (const TecplotDecodedZone* d : decoded) {
                NDArray arr(DType::Float64, {d->mCols[k].size()});
                std::memcpy(arr.Data(), d->mCols[k].data(), d->mCols[k].size() * sizeof(double));
                blk.push_back(std::move(arr));
            }
            mesh.AddCellData(rVariables[k], std::move(blk));
        } else {
            std::vector<double> col;
            for (const TecplotDecodedZone* d : decoded)
                col.insert(col.end(), d->mCols[k].begin(), d->mCols[k].end());
            NDArray arr(DType::Float64, {col.size()});
            std::memcpy(arr.Data(), col.data(), col.size() * sizeof(double));
            mesh.AddPointData(rVariables[k], std::move(arr));
        }
    }

    std::vector<NDArray> zone_ids;
    zone_ids.reserve(decoded.size());
    for (std::size_t k = 0; k < decoded.size(); ++k) {
        NDArray zn(DType::Int64, {decoded[k]->mNumCells});
        std::int64_t* zp = zn.As<std::int64_t>();
        std::fill(zp, zp + decoded[k]->mNumCells, static_cast<std::int64_t>(rZoneIdxs[k]));
        zone_ids.push_back(std::move(zn));
    }
    mesh.AddCellData("tecplot:zone", std::move(zone_ids));

    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    std::vector<std::string> used_names;
    for (std::size_t k = 0; k < rZoneIdxs.size(); ++k) {
        std::string name = rZones[rZoneIdxs[k]].mTitle;
        if (name.empty())
            name = "zone_" + std::to_string(k);
        std::string unique = name;
        for (int suffix = 2; std::find(used_names.begin(), used_names.end(), unique) !=
                             used_names.end();
            ++suffix)
            unique = name + "_" + std::to_string(suffix);
        used_names.push_back(unique);

        NDArray entries(DType::Int64, {decoded[k]->mNumCells});
        std::int64_t* ep = entries.As<std::int64_t>();
        for (std::size_t r = 0; r < decoded[k]->mNumCells; ++r)
            ep[r] = detail::block_row_to_global(bases, k, static_cast<std::int64_t>(r));
        mesh.AddRegion(Region(unique, RegionKind::Cell, -1, static_cast<std::int64_t>(rZoneIdxs[k]),
                              std::move(entries)));
    }
    return mesh;
}

}  // namespace

MeshMetadata read_tecplot_metadata(const std::string& rPath, const ReadOptions& /*rOptions*/) {
    auto in = detail::make_classic_ifstream(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        std::string s = tecplot_strip(l);
        if (s.empty() || s[0] == '#')
            continue;
        lines.push_back(s);
    }

    std::vector<std::string> variables;
    const std::vector<TecplotZoneHeader> zones = tecplot_scan_zones(lines, variables);
    const std::vector<std::vector<std::size_t>> timeline = tecplot_timeline(zones);

    MeshMetadata meta;
    meta.mFormat = "tecplot";
    const std::vector<std::size_t>& first_step = timeline[0];
    int xi = -1, yi = -1, zi = -1;
    tecplot_xyz_indices(variables, xi, yi, zi);
    std::map<std::size_t, std::size_t> idx_to_pos;
    for (std::size_t k = 0; k < first_step.size(); ++k)
        idx_to_pos.emplace(first_step[k], k);
    std::size_t total_points = 0;
    for (std::size_t k = 0; k < first_step.size(); ++k) {
        const std::size_t idx = first_step[k];
        const std::size_t origin = tecplot_points_origin_zone(idx, zones, xi, yi, zi);
        const auto it = idx_to_pos.find(origin);
        const bool shares_earlier =
            origin != idx && it != idx_to_pos.end() && it->second < k;
        if (!shares_earlier)
            total_points += zones[idx].mNumNodes;
        bool feblock = false;
        std::string ztype;
        std::vector<int> cell_centered;
        tecplot_zone_format(zones[idx], variables.size(), feblock, ztype, cell_centered);
        CellBlockInfo block;
        block.mType = tecplot_to_meshio(ztype);
        block.mNumCells = zones[idx].mNumCells;
        meta.mCellBlocks.push_back(std::move(block));
    }
    meta.mNumPoints = total_points;
    meta.mPointDim = 0;  // not knowable without decoding X/Y/Z columns
    if (zones[0].mHasSolutionTime)
        for (const std::vector<std::size_t>& step : timeline)
            meta.mTimeValues.push_back(zones[step[0]].mSolutionTime);
    return meta;
}

Mesh read_tecplot(const std::string& rPath, const ReadOptions& rOptions) {
    auto in = detail::make_classic_ifstream(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        std::string s = tecplot_strip(l);
        if (s.empty() || s[0] == '#')
            continue;
        lines.push_back(s);
    }

    std::vector<std::string> variables;
    const std::vector<TecplotZoneHeader> zones = tecplot_scan_zones(lines, variables);
    const std::vector<std::vector<std::size_t>> timeline = tecplot_timeline(zones);
    const std::size_t step = rOptions.ResolveTimeStep(timeline.size());
    std::map<std::size_t, TecplotDecodedZone> cache;
    return tecplot_build_step_mesh(timeline[step], zones, variables, lines, cache);
}

Mesh read_tecplot(const std::string& rPath) {
    return read_tecplot(rPath, ReadOptions{});
}

// Writes one Tecplot ZONE per cell block (no single-type restriction). Zone
// 1 carries the coordinates and every nodal field; later zones reuse them
// through VARSHARELIST rather than duplicating the (unwelded, shared) point
// array. Cell-centred variables are per-block: a zone that owns no values
// for a given variable declares it PASSIVEVARLIST instead of writing zeros.
void write_tecplot(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::size_t> blocks;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const auto cb = rMesh.Cells(i);
        if (!meshio_to_tecplot(cb.Type()).empty())
            blocks.push_back(i);
        else
            log::warn("tecplot: skipping unsupported cell type '{}'", cb.Type());
    }
    if (blocks.empty())
        throw WriteError("Tecplot writer found no supported cell blocks");

    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t dim = rMesh.PointDim();
    const std::size_t num_nodes = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();

    // Shared (nodal) variables: X, Y, (Z), then every point-data variable.
    std::vector<std::string> variables = {"X", "Y"};
    if (dim == 3)
        variables.push_back("Z");
    std::vector<std::vector<double>> shared_data;
    auto push_point_col = [&](std::size_t comp) {
        std::vector<double> col(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            col[r] = detail::read_double(points, r * dim + comp);
        shared_data.push_back(std::move(col));
    };
    push_point_col(0);
    push_point_col(1);
    if (dim == 3)
        push_point_col(2);
    for (const auto& k : rMesh.PointDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        const NDArray& v = rMesh.PointData(k);
        std::size_t ncomp = v.Shape().size() >= 2 ? v.Shape()[1] : 1;
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            std::vector<double> col(num_nodes);
            for (std::size_t r = 0; r < num_nodes; ++r)
                col[r] = detail::read_double(v, r * ncomp + c);
            shared_data.push_back(std::move(col));
        }
    }
    const std::size_t num_shared = variables.size();

    // Cell-centred variables (component-expanded); `mKey`/`mComp`/`mNumComp`
    // let each zone look its own values up (or find them absent) later.
    struct CellVarComp {
        std::string mKey;
        std::size_t mComp;
        std::size_t mNumComp;
    };
    std::vector<CellVarComp> cell_vars;
    for (const auto& k : rMesh.CellDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        if (rMesh.CellDataNumBlocks(k) == 0)
            continue;
        std::size_t ncomp = 1;
        for (std::size_t b : blocks) {
            if (b < rMesh.CellDataNumBlocks(k)) {
                const NDArray& vv = rMesh.CellData(k, b);
                ncomp = vv.Shape().size() >= 2 ? vv.Shape()[1] : 1;
                break;
            }
        }
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            cell_vars.push_back({k, c, ncomp});
        }
    }

    os << "TITLE = \"" << detail::provenance_lines(detail::SlotTier::SingleLine)[0] << "\"\n";
    os << "VARIABLES = ";
    for (std::size_t k = 0; k < variables.size(); ++k)
        os << (k ? ", " : "") << "\"" << variables[k] << "\"";
    os << "\n";

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    char buf[40];
    auto write_column = [&](const std::vector<double>& col) {
        for (std::size_t i = 0; i < col.size(); ++i) {
            detail::snprintf_c(buf, sizeof(buf), "%.17g", col[i]);
            os << buf << ((i + 1) % 20 == 0 || i + 1 == col.size() ? '\n' : ' ');
        }
        if (col.empty())
            os << "\n";
    };

    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        const std::size_t block = blocks[bi];
        const auto cb = rMesh.Cells(block);
        const std::string ztype = meshio_to_tecplot(cb.Type());
        const std::vector<int>& order = tecplot_order(cb.Type());
        const std::size_t num_cells = cb.NumCells();

        // Zone title: an exactly-matching Cell region's name, else block_<i>.
        std::string title = "block_" + std::to_string(block);
        for (std::size_t ri = 0; ri < rMesh.NumRegions(); ++ri) {
            const meshioplusplus::Region& reg = rMesh.Region(ri);
            if (reg.mKind != RegionKind::Cell || reg.NumEntries() == 0)
                continue;
            if (reg.NumEntries() !=
                static_cast<std::size_t>(bases[block + 1] - bases[block]))
                continue;
            const std::int64_t* e = reg.Entries();
            if (e[0] == bases[block] && e[reg.NumEntries() - 1] == bases[block + 1] - 1) {
                title = reg.mName;
                break;
            }
        }

        // Which cell variables this zone owns vs. leaves passive.
        std::vector<std::size_t> present, passive;
        for (std::size_t j = 0; j < cell_vars.size(); ++j) {
            if (block < rMesh.CellDataNumBlocks(cell_vars[j].mKey))
                present.push_back(j);
            else
                passive.push_back(j);
        }

        os << "ZONE T = \"" << title << "\", NODES = " << num_nodes
           << ", ELEMENTS = " << num_cells << ",\n";
        os << "DATAPACKING = BLOCK, ZONETYPE = " << ztype;
        if (bi > 0)
            os << ",\nVARSHARELIST = ([1-" << num_shared << "] = 1)";
        if (!present.empty()) {
            os << ",\nVARLOCATION = ([";
            for (std::size_t i = 0; i < present.size(); ++i)
                os << (i ? "," : "") << (num_shared + present[i] + 1);
            os << "] = CELLCENTERED)";
        }
        if (!passive.empty()) {
            os << ",\nPASSIVEVARLIST = ([";
            for (std::size_t i = 0; i < passive.size(); ++i)
                os << (i ? "," : "") << (num_shared + passive[i] + 1);
            os << "])";
        }
        os << "\n";

        if (bi == 0)
            for (const auto& col : shared_data)
                write_column(col);
        for (std::size_t j : present) {
            const CellVarComp& cv = cell_vars[j];
            const NDArray& vv = rMesh.CellData(cv.mKey, block);
            std::vector<double> col(vv.Shape()[0]);
            for (std::size_t r = 0; r < vv.Shape()[0]; ++r)
                col[r] = detail::read_double(vv, r * cv.mNumComp + cv.mComp);
            write_column(col);
        }

        const NDArray& conn = cb.Conn();
        const std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < num_cells; ++r) {
            for (std::size_t j = 0; j < order.size(); ++j) {
                std::size_t src = static_cast<std::size_t>(order[j]);
                if (src >= k)
                    src = k - 1;
                os << (detail::read_int(conn, r * k + src) + 1)
                   << (j + 1 == order.size() ? '\n' : ' ');
            }
        }
    }
}

}  // namespace meshioplusplus
