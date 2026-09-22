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
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus {

namespace {

constexpr double kUnvPi = 3.14159265358979323846;
constexpr double kUnvNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------
// Element tables
// ---------------------------------------------------------------------------

// UNV node order -> meshio position (0-based): meshio_conn[perm[i]] = unv_conn[i].
// Parabolic elements list their mid-side nodes "sandwiched" between the corners
// of each ring; the solids list the bottom ring, then the vertical mid-edges,
// then the top ring (pinned against gmsh's .unv/.msh twins and Salome's driver).
const std::vector<int>* unv_perm(const std::string& rType) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"line3", {0, 2, 1}},
        {"triangle6", {0, 3, 1, 4, 2, 5}},
        {"quad8", {0, 4, 1, 5, 2, 6, 3, 7}},
        {"quad9", {0, 4, 1, 5, 2, 6, 3, 7, 8}},
        {"tetra10", {0, 4, 1, 5, 2, 6, 7, 8, 9, 3}},
        {"pyramid13", {0, 5, 1, 6, 2, 7, 3, 8, 9, 10, 11, 12, 4}},
        {"wedge15", {0, 6, 1, 7, 2, 8, 12, 13, 14, 3, 9, 4, 10, 5, 11}},
        {"hexahedron20", {0, 8, 1, 9, 2, 10, 3, 11, 16, 17, 18, 19, 4, 12, 5, 13, 6, 14, 7, 15}}};
    auto it = m.find(rType);
    return it == m.end() ? nullptr : &it->second;
}

bool unv_is_beam(int FeId) {
    return FeId == 11 || (FeId >= 21 && FeId <= 25);
}

// FE descriptor id + node count -> meshio cell type; empty when unsupported.
// The node count must match the type: gmsh, for one, writes its quad9/hex27 under
// the linear ids 94/115 in its own order, which is refused rather than misread.
std::string unv_cell_type(int FeId, std::size_t NumNodes) {
    auto pick = [&](const char* pType, std::size_t N) {
        return NumNodes == N ? std::string(pType) : std::string();
    };
    if (unv_is_beam(FeId))
        return NumNodes == 2 ? "line" : (NumNodes == 3 ? "line3" : std::string());
    switch (FeId) {
        case 41:
        case 51:
        case 61:
        case 74:
        case 81:
        case 91:
            return pick("triangle", 3);
        case 42:
        case 52:
        case 62:
        case 72:
        case 82:
        case 92:
            return pick("triangle6", 6);
        case 44:
        case 54:
        case 64:
        case 71:
        case 84:
        case 94:
        case 122:
            return pick("quad", 4);
        case 45:
        case 55:
        case 65:
        case 75:
        case 85:
        case 95:
            return NumNodes == 9 ? "quad9" : pick("quad8", 8);
        case 111:
            return pick("tetra", 4);
        case 118:
            return pick("tetra10", 10);
        case 112:
            return pick("wedge", 6);
        case 113:
            return pick("wedge15", 15);
        case 115:
            return pick("hexahedron", 8);
        case 116:
            return pick("hexahedron20", 20);
        case 119:
        case 312:
            return pick("pyramid", 5);
        case 114:
            return pick("pyramid13", 13);
        default:
            return std::string();
    }
}

// meshio type -> (descriptor, is_beam) used on write.
bool unv_descriptor(const std::string& rType, int& rDesc, bool& rBeam) {
    static const std::unordered_map<std::string, std::pair<int, bool>> m = {
        {"line", {21, true}},           {"line3", {24, true}},     {"triangle", {91, false}},
        {"triangle6", {92, false}},     {"quad", {94, false}},     {"quad8", {95, false}},
        {"quad9", {95, false}},         {"tetra", {111, false}},   {"tetra10", {118, false}},
        {"wedge", {112, false}},        {"wedge15", {113, false}}, {"hexahedron", {115, false}},
        {"hexahedron20", {116, false}}, {"pyramid", {312, false}}, {"pyramid13", {114, false}}};
    auto it = m.find(rType);
    if (it == m.end())
        return false;
    rDesc = it->second.first;
    rBeam = it->second.second;
    return true;
}

// Permanent-group datasets and whether each entity is a (type, tag, leaf, component)
// quadruple (true) or a bare (type, tag) pair (false).
bool unv_group_layout(int Id, bool& rQuad) {
    switch (Id) {
        case 2417:
        case 2429:
        case 2430:
        case 2432:
            rQuad = false;
            return true;
        case 2435:
        case 2452:
        case 2467:
        case 2477:
            rQuad = true;
            return true;
        default:
            return false;
    }
}

// Result-type code (2414 record 9 field 4 / 55 record 6 field 4) -> a data name used
// when the dataset's own name is blank or "NONE".
std::string unv_result_type_name(int Code) {
    static const std::unordered_map<int, const char*> m = {{1, "general"},
                                                           {2, "stress"},
                                                           {3, "strain"},
                                                           {4, "element_force"},
                                                           {5, "temperature"},
                                                           {6, "heat_flux"},
                                                           {7, "strain_energy"},
                                                           {8, "displacement"},
                                                           {9, "reaction_force"},
                                                           {10, "kinetic_energy"},
                                                           {11, "velocity"},
                                                           {12, "acceleration"},
                                                           {13, "strain_energy_density"},
                                                           {14, "kinetic_energy_density"},
                                                           {15, "hydrostatic_pressure"},
                                                           {16, "heat_gradient"},
                                                           {17, "code_checking_value"},
                                                           {18, "pressure_coefficient"}};
    auto it = m.find(Code);
    return it == m.end() ? std::string("unv:field") : std::string(it->second);
}

// Dataset-58 function type -> data name.
std::string unv_function_name(int Type) {
    static const char* const kNames[] = {"function",
                                         "time_response",
                                         "auto_spectrum",
                                         "cross_spectrum",
                                         "frf",
                                         "transmissibility",
                                         "coherence",
                                         "auto_correlation",
                                         "cross_correlation",
                                         "psd",
                                         "esd",
                                         "pdf",
                                         "spectrum",
                                         "cumulative_frequency_distribution",
                                         "peaks_valley",
                                         "stress_cycles",
                                         "strain_cycles",
                                         "orbit",
                                         "mode_indicator_function",
                                         "force_pattern",
                                         "partial_power",
                                         "partial_coherence",
                                         "eigenvalue",
                                         "eigenvector",
                                         "shock_response_spectrum",
                                         "fir_filter",
                                         "multiple_coherence",
                                         "order_function"};
    if (Type < 0 || Type >= static_cast<int>(sizeof(kNames) / sizeof(kNames[0])))
        return "function_" + std::to_string(Type);
    return kNames[Type];
}

// Component count -> data characteristic (1 scalar, 2 3-DOF vector, 4 symmetric
// tensor, 5 general tensor); 0 when unrecognized.
int unv_data_char(std::size_t NumComps) {
    switch (NumComps) {
        case 1:
            return 1;
        case 3:
            return 2;
        case 6:
            return 4;
        case 9:
            return 5;
        default:
            return 0;
    }
}

// File tensor order -> meshio order. A symmetric tensor is stored
// Sxx Sxy Syy Sxz Syz Szz (meshio: xx yy zz xy yz zx); a general one column by
// column (meshio: row-major).
void unv_tensor_to_meshio(int DataChar, std::size_t NumComps, double* pValues) {
    if (DataChar == 4 && NumComps == 6) {
        const double f[6] = {pValues[0], pValues[1], pValues[2],
                             pValues[3], pValues[4], pValues[5]};
        const int from[6] = {0, 2, 5, 1, 4, 3};
        for (int k = 0; k < 6; ++k)
            pValues[k] = f[from[k]];
    } else if (DataChar == 5 && NumComps == 9) {
        double f[9];
        std::copy(pValues, pValues + 9, f);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                pValues[r * 3 + c] = f[c * 3 + r];
    }
}

void unv_tensor_from_meshio(std::size_t NumComps, double* pValues) {
    if (NumComps == 6) {
        const double m[6] = {pValues[0], pValues[1], pValues[2],
                             pValues[3], pValues[4], pValues[5]};
        const int from[6] = {0, 3, 1, 5, 4, 2};
        for (int k = 0; k < 6; ++k)
            pValues[k] = m[from[k]];
    } else if (NumComps == 9) {
        double m[9];
        std::copy(pValues, pValues + 9, m);
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                pValues[c * 3 + r] = m[r * 3 + c];
    }
}

// ---------------------------------------------------------------------------
// Text helpers
// ---------------------------------------------------------------------------

std::string_view unv_strip(std::string_view s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string_view::npos)
        return {};
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::vector<std::string_view> unv_split(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r'))
            ++i;
        std::size_t j = i;
        while (j < s.size() && s[j] != ' ' && s[j] != '\t' && s[j] != '\r')
            ++j;
        if (j > i)
            out.push_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}

std::int64_t unv_int(std::string_view t) {
    if (!t.empty() && t.front() == '+')
        t.remove_prefix(1);
    std::int64_t v = 0;
    auto res = std::from_chars(t.data(), t.data() + t.size(), v);
    if (res.ec != std::errc())
        throw ReadError("UNV: expected an integer, got '" + std::string(t) + "'");
    return v;
}

std::vector<std::int64_t> unv_ints(std::string_view line) {
    std::vector<std::int64_t> out;
    for (auto t : unv_split(line))
        out.push_back(unv_int(t));
    return out;
}

double unv_real(std::string_view t) {
    std::string s(t);
    for (char& c : s)
        if (c == 'D' || c == 'd')
            c = 'E';
    const char* end = nullptr;
    const double v = detail::parse_double(s.c_str(), end);
    if (end == s.c_str())
        throw ReadError("UNV: expected a real number, got '" + s + "'");
    return v;
}

// Real tokens of a line. Fixed-width fields (E13.5 and friends) run together when a
// value is negative and fills its field, so a sign that does not follow an exponent
// letter also starts a new number.
std::vector<double> unv_reals(std::string_view line) {
    std::vector<double> out;
    for (auto t : unv_split(line)) {
        std::size_t start = 0;
        for (std::size_t k = 1; k < t.size(); ++k) {
            const char c = t[k];
            const char p = t[k - 1];
            if ((c == '-' || c == '+') && p != 'E' && p != 'e' && p != 'D' && p != 'd') {
                out.push_back(unv_real(t.substr(start, k - start)));
                start = k;
            }
        }
        out.push_back(unv_real(t.substr(start)));
    }
    return out;
}

std::string unv_name(std::string_view s) {
    return std::string(unv_strip(s));
}

bool unv_is_none(const std::string& rName) {
    if (rName.size() != 4)
        return rName.empty();
    std::string up = rName;
    for (char& c : up)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return up == "NONE";
}

// ---------------------------------------------------------------------------
// Dataset splitting (byte-safe: 58b carries raw binary data)
// ---------------------------------------------------------------------------

struct UnvDataset {
    int mId = 0;
    bool mBinary = false;
    int mByteOrder = 1;  // 1 little endian, 2 big endian
    int mFpFormat = 2;   // 2 IEEE 754
    std::vector<std::string_view> mLines;
    std::string_view mBlob;
};

class UnvLineReader {
public:
    explicit UnvLineReader(std::string_view data) : mData(data) {}
    bool AtEnd() const { return mPos >= mData.size(); }
    std::size_t Pos() const { return mPos; }
    void Seek(std::size_t pos) { mPos = pos; }
    std::string_view Next() {
        const std::size_t eol = mData.find('\n', mPos);
        const std::size_t end = eol == std::string_view::npos ? mData.size() : eol;
        std::string_view line = mData.substr(mPos, end - mPos);
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        mPos = eol == std::string_view::npos ? mData.size() : eol + 1;
        return line;
    }
    std::string_view Data() const { return mData; }

private:
    std::string_view mData;
    std::size_t mPos = 0;
};

// Bytes of a 58b data block, from record 7 (ordinate type, count, spacing); -1 when
// the record cannot be read.
std::int64_t unv_58b_bytes(const std::vector<std::string_view>& rLines) {
    if (rLines.size() < 7)
        return -1;
    auto r7 = unv_split(rLines[6]);
    if (r7.size() < 3)
        return -1;
    try {
        const std::int64_t type = unv_int(r7[0]);
        const std::int64_t npts = unv_int(r7[1]);
        const bool even = unv_int(r7[2]) == 1;
        if (npts < 0 || (type != 2 && type != 4 && type != 5 && type != 6))
            return -1;
        const std::int64_t per = ((type == 5 || type == 6) ? 2 : 1) + (even ? 0 : 1);
        return npts * per * ((type == 4 || type == 6) ? 8 : 4);
    } catch (const ReadError&) {
        return -1;
    }
}

std::vector<UnvDataset> unv_split_datasets(std::string_view data) {
    std::vector<UnvDataset> out;
    UnvLineReader in(data);
    bool warned_58b_size = false;
    while (!in.AtEnd()) {
        if (unv_strip(in.Next()) != "-1")
            continue;
        // Consecutive "-1" lines: the second one opens the dataset.
        std::string_view header;
        do {
            if (in.AtEnd())
                return out;
            header = in.Next();
        } while (unv_strip(header) == "-1" || unv_strip(header).empty());
        auto tokens = unv_split(header);
        UnvDataset ds;
        std::string_view id = tokens[0];
        if (!id.empty() && (id.back() == 'b' || id.back() == 'B')) {
            ds.mBinary = true;
            id.remove_suffix(1);
        }
        ds.mId = static_cast<int>(unv_int(id));
        if (ds.mBinary) {
            if (tokens.size() < 5)
                throw ReadError("UNV: malformed binary dataset header '" + std::string(header) +
                                "'");
            ds.mByteOrder = static_cast<int>(unv_int(tokens[1]));
            ds.mFpFormat = static_cast<int>(unv_int(tokens[2]));
            const std::int64_t n_ascii = unv_int(tokens[3]);
            std::int64_t n_bytes = unv_int(tokens[4]);
            for (std::int64_t k = 0; k < n_ascii && !in.AtEnd(); ++k)
                ds.mLines.push_back(in.Next());
            // A 58b's size follows from its record 7; some writers (pyuff among them)
            // declare half of it for complex data, so the record wins.
            const std::int64_t expected = ds.mId == 58 ? unv_58b_bytes(ds.mLines) : -1;
            if (expected >= 0 && expected != n_bytes) {
                // pyuff declares exactly half the size of complex data; any other
                // mismatch deserves a warning.
                if (2 * n_bytes == expected) {
                    log::debug("UNV: dataset 58b declares {} bytes, record 7 describes {}",
                               n_bytes, expected);
                } else if (!warned_58b_size) {
                    log::warn(
                        "UNV: dataset 58b declares {} bytes but its record 7 describes {}; "
                        "using the record",
                        n_bytes, expected);
                    warned_58b_size = true;
                }
                n_bytes = expected;
            }
            if (n_bytes < 0 || in.Pos() + static_cast<std::size_t>(n_bytes) > data.size())
                throw ReadError("UNV: binary dataset " + std::to_string(ds.mId) + " is truncated");
            ds.mBlob = data.substr(in.Pos(), static_cast<std::size_t>(n_bytes));
            in.Seek(in.Pos() + static_cast<std::size_t>(n_bytes));
            // The closing "-1" follows the binary data directly (no newline first).
            const std::string_view rest = in.Next();
            if (unv_strip(rest) != "-1")
                log::warn("UNV: binary dataset {} is not closed by '-1'", ds.mId);
        } else {
            while (!in.AtEnd()) {
                const std::string_view line = in.Next();
                if (unv_strip(line) == "-1")
                    break;
                ds.mLines.push_back(line);
            }
        }
        out.push_back(std::move(ds));
    }
    return out;
}

// ---------------------------------------------------------------------------
// Parsed file
// ---------------------------------------------------------------------------

struct UnvNode {
    std::int64_t mLabel = 0;
    std::int64_t mCs = 0;
    std::size_t mNumCoords = 3;
    double mX[3] = {0.0, 0.0, 0.0};
};

struct UnvElement {
    std::int64_t mLabel = 0;
    int mFeId = 0;
    std::int64_t mPid = 0;
    std::int64_t mMid = 0;
    std::vector<std::int64_t> mNodes;
};

struct UnvGroup {
    std::int64_t mNumber = 0;
    std::string mName;
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mElements;
};

struct UnvCs {
    int mType = 0;  // 0 Cartesian, 1 cylindrical, 2 spherical
    double mM[4][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
};

struct UnvStepKey {
    int mAnalysis = 0;
    std::int64_t mId = 0;
    double mValue = 0.0;
    bool operator==(const UnvStepKey& rOther) const {
        return mAnalysis == rOther.mAnalysis && mId == rOther.mId && mValue == rOther.mValue;
    }
};

// A results dataset (2414, 55, 56): values per entity label.
struct UnvResult {
    int mDataset = 2414;
    int mLocation = 1;  // 1 nodes, 2 elements
    std::string mName;
    int mDataChar = 0;
    std::size_t mNumComps = 1;
    bool mComplex = false;
    UnvStepKey mKey;
    std::vector<std::pair<std::int64_t, std::vector<double>>> mValues;
};

// One dataset-58 function.
struct UnvFunction {
    int mType = 0;
    std::int64_t mLoadCase = 0;
    std::int64_t mRspNode = 0;
    int mRspDir = 0;
    std::int64_t mRefNode = 0;
    int mRefDir = 0;
    int mAbscissaType = 0;
    bool mComplex = false;
    std::vector<double> mX;
    std::vector<double> mRe;
    std::vector<double> mIm;
};

struct UnvFile {
    std::vector<UnvNode> mNodes;
    std::vector<UnvElement> mElements;
    std::vector<UnvGroup> mGroups;
    std::unordered_map<std::int64_t, UnvCs> mCs;
    bool mHasUnits = false;
    std::int64_t mUnitsCode = 0;
    std::vector<double> mUnitFactors;
    std::vector<UnvResult> mResults;
    std::vector<UnvFunction> mFunctions;
};

// Gather `Count` reals starting at line `rK` (advancing it).
std::vector<double> unv_take_reals(const std::vector<std::string_view>& rLines, std::size_t& rK,
                                   std::size_t Count) {
    std::vector<double> vals;
    while (vals.size() < Count && rK < rLines.size()) {
        for (double v : unv_reals(rLines[rK]))
            vals.push_back(v);
        ++rK;
    }
    if (vals.size() < Count)
        throw ReadError("UNV: dataset ends inside a data record");
    vals.resize(Count);
    return vals;
}

std::vector<std::int64_t> unv_take_ints(const std::vector<std::string_view>& rLines,
                                        std::size_t& rK, std::size_t Count) {
    std::vector<std::int64_t> vals;
    while (vals.size() < Count && rK < rLines.size()) {
        for (std::int64_t v : unv_ints(rLines[rK]))
            vals.push_back(v);
        ++rK;
    }
    if (vals.size() < Count)
        throw ReadError("UNV: dataset ends inside an integer record");
    vals.resize(Count);
    return vals;
}

void unv_parse_nodes(const UnvDataset& rDs, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    std::size_t k = 0;
    while (k < lines.size()) {
        if (unv_strip(lines[k]).empty()) {
            ++k;
            continue;
        }
        UnvNode node;
        if (rDs.mId == 15) {
            // 4I10,1P3E13.5 on one line.
            const std::string_view line = lines[k++];
            auto ints = unv_split(line.substr(0, std::min<std::size_t>(40, line.size())));
            if (ints.size() < 2)
                throw ReadError("UNV: malformed dataset-15 node record");
            node.mLabel = unv_int(ints[0]);
            node.mCs = unv_int(ints[1]);
            auto xs = line.size() > 40 ? unv_reals(line.substr(40)) : std::vector<double>{};
            node.mNumCoords = std::min<std::size_t>(3, xs.size());
            for (std::size_t c = 0; c < node.mNumCoords; ++c)
                node.mX[c] = xs[c];
        } else {
            // 2411 / 781: `label def_cs disp_cs colour`, then the coordinates.
            auto r1 = unv_ints(lines[k]);
            if (r1.empty() || k + 1 >= lines.size())
                throw ReadError("UNV: malformed node record in dataset " + std::to_string(rDs.mId));
            node.mLabel = r1[0];
            node.mCs = r1.size() > 1 ? r1[1] : 0;
            auto xs = unv_reals(lines[k + 1]);
            node.mNumCoords = std::min<std::size_t>(3, xs.size());
            for (std::size_t c = 0; c < node.mNumCoords; ++c)
                node.mX[c] = xs[c];
            k += 2;
        }
        rFile.mNodes.push_back(node);
    }
}

void unv_parse_elements(const UnvDataset& rDs, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    std::size_t k = 0;
    while (k < lines.size()) {
        if (unv_strip(lines[k]).empty()) {
            ++k;
            continue;
        }
        auto r1 = unv_ints(lines[k++]);
        UnvElement el;
        std::int64_t num_nodes = 0;
        if (rDs.mId == 780) {
            // label, FE id, phys bin, phys prop, mat bin, mat prop, colour, node count
            if (r1.size() < 8)
                throw ReadError("UNV: malformed dataset-780 element record");
            el.mLabel = r1[0];
            el.mFeId = static_cast<int>(r1[1]);
            el.mPid = r1[3];
            el.mMid = r1[5];
            num_nodes = r1[7];
        } else {
            if (r1.size() < 6)
                throw ReadError("UNV: malformed dataset-2412 element record");
            el.mLabel = r1[0];
            el.mFeId = static_cast<int>(r1[1]);
            el.mPid = r1[2];
            el.mMid = r1[3];
            num_nodes = r1[5];
        }
        if (num_nodes < 0)
            throw ReadError("UNV: negative node count on element " + std::to_string(el.mLabel));
        if (unv_is_beam(el.mFeId))
            ++k;  // orientation node and cross sections
        el.mNodes = unv_take_ints(lines, k, static_cast<std::size_t>(num_nodes));
        rFile.mElements.push_back(std::move(el));
    }
}

void unv_parse_groups(const UnvDataset& rDs, bool Quad, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    const std::size_t stride = Quad ? 4 : 2;
    std::size_t k = 0;
    while (k < lines.size()) {
        if (unv_strip(lines[k]).empty()) {
            ++k;
            continue;
        }
        auto r1 = unv_ints(lines[k++]);
        if (r1.size() < 8)
            throw ReadError("UNV: malformed group record in dataset " + std::to_string(rDs.mId));
        if (r1[7] < 0)
            throw ReadError("UNV: negative entity count in dataset " + std::to_string(rDs.mId));
        UnvGroup g;
        g.mNumber = r1[0];
        g.mName = k < lines.size() ? unv_name(lines[k]) : std::string();
        ++k;
        auto vals = unv_take_ints(lines, k, stride * static_cast<std::size_t>(r1[7]));
        for (std::size_t e = 0; e < static_cast<std::size_t>(r1[7]); ++e) {
            const std::int64_t type = vals[stride * e];
            const std::int64_t tag = vals[stride * e + 1];
            if (type == 7)
                g.mNodes.push_back(tag);
            else if (type == 8)
                g.mElements.push_back(tag);
        }
        rFile.mGroups.push_back(std::move(g));
    }
}

void unv_parse_units(const UnvDataset& rDs, UnvFile& rFile) {
    if (rDs.mLines.empty())
        return;
    auto r1 = unv_split(rDs.mLines[0]);
    if (r1.empty())
        return;
    rFile.mHasUnits = true;
    rFile.mUnitsCode = unv_int(r1[0]);
    std::size_t k = 1;
    rFile.mUnitFactors.clear();
    while (k < rDs.mLines.size() && rFile.mUnitFactors.size() < 4) {
        for (double v : unv_reals(rDs.mLines[k]))
            rFile.mUnitFactors.push_back(v);
        ++k;
    }
    rFile.mUnitFactors.resize(4, 0.0);
}

void unv_parse_cs(const UnvDataset& rDs, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    std::size_t k = 2;  // part UID, part name
    while (k < lines.size()) {
        if (unv_strip(lines[k]).empty()) {
            ++k;
            continue;
        }
        auto r3 = unv_ints(lines[k]);
        if (r3.size() < 2 || k + 2 >= lines.size())
            break;
        UnvCs cs;
        cs.mType = static_cast<int>(r3[1]);
        std::size_t j = k + 2;
        auto m = unv_take_reals(lines, j, 12);
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 3; ++c)
                cs.mM[r][c] = m[r * 3 + c];
        rFile.mCs[r3[0]] = cs;
        k = j;
    }
}

// Step key of a 2414 dataset from records 9-13.
UnvStepKey unv_key_2414(int Analysis, const std::vector<std::int64_t>& rInts,
                        const std::vector<double>& rReals) {
    auto i = [&](std::size_t n) { return n < rInts.size() ? rInts[n] : 0; };
    auto r = [&](std::size_t n) { return n < rReals.size() ? rReals[n] : 0.0; };
    UnvStepKey key;
    key.mAnalysis = Analysis;
    switch (Analysis) {
        case 2:  // normal mode: mode number, frequency
            key.mId = i(5);
            key.mValue = r(1);
            break;
        case 3:
        case 7:  // complex eigenvalue: mode number, |Im(lambda)| / 2 pi
            key.mId = i(5);
            key.mValue = std::abs(r(7)) / (2.0 * kUnvPi);
            break;
        case 4:  // transient: time step number, time
        case 9:  // static non-linear
            key.mId = i(6);
            key.mValue = r(0);
            break;
        case 5:  // frequency response: frequency number, frequency
            key.mId = i(7);
            key.mValue = r(1);
            break;
        case 6:  // buckling: mode number, eigenvalue
            key.mId = i(5);
            key.mValue = r(2);
            break;
        default:  // static / unknown: load set
            key.mId = i(4);
            key.mValue = r(0) != 0.0 ? r(0) : static_cast<double>(i(4));
            break;
    }
    return key;
}

// Step key of a 55/56 dataset from records 7-8.
UnvStepKey unv_key_55(int Analysis, const std::vector<std::int64_t>& rInts,
                      const std::vector<double>& rReals) {
    auto i = [&](std::size_t n) { return n < rInts.size() ? rInts[n] : 0; };
    auto r = [&](std::size_t n) { return n < rReals.size() ? rReals[n] : 0.0; };
    UnvStepKey key;
    key.mAnalysis = Analysis;
    switch (Analysis) {
        case 2:
        case 4:
        case 5:
            key.mId = i(1);
            key.mValue = r(0);
            break;
        case 3:
        case 7:
            key.mId = i(1);
            key.mValue = std::abs(r(1)) / (2.0 * kUnvPi);
            break;
        case 6:
            key.mId = i(0);
            key.mValue = r(0);
            break;
        default:
            key.mId = i(0);
            key.mValue = r(0) != 0.0 ? r(0) : static_cast<double>(i(0));
            break;
    }
    return key;
}

void unv_parse_result(const UnvDataset& rDs, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    UnvResult res;
    res.mDataset = rDs.mId;
    std::size_t k = 0;
    int data_type = 2;
    std::size_t ndv = 1;
    int result_type = 0;
    bool legacy = false;
    if (rDs.mId == 2414) {
        if (lines.size() < 13)
            throw ReadError("UNV: dataset 2414 has a truncated header");
        res.mName = unv_name(lines[1]);
        auto loc = unv_ints(lines[2]);
        res.mLocation = loc.empty() ? 1 : static_cast<int>(loc[0]);
        auto r9 = unv_ints(lines[8]);
        if (r9.size() < 6)
            throw ReadError("UNV: dataset 2414 record 9 needs six integers");
        const int analysis = static_cast<int>(r9[1]);
        res.mDataChar = static_cast<int>(r9[2]);
        result_type = static_cast<int>(r9[3]);
        data_type = static_cast<int>(r9[4]);
        ndv = static_cast<std::size_t>(std::max<std::int64_t>(r9[5], 0));
        std::vector<std::int64_t> ints = unv_ints(lines[9]);
        for (std::int64_t v : unv_ints(lines[10]))
            ints.push_back(v);
        std::vector<double> reals = unv_reals(lines[11]);
        for (double v : unv_reals(lines[12]))
            reals.push_back(v);
        res.mKey = unv_key_2414(analysis, ints, reals);
        k = 13;
    } else {
        // 55 (data at nodes), 56 (data at elements): ID lines, record 6, then records 7-8
        // sized by their own NINT/NRV counts.
        if (lines.size() < 7)
            throw ReadError("UNV: dataset " + std::to_string(rDs.mId) + " has a truncated header");
        res.mName = unv_name(lines[0]);
        res.mLocation = rDs.mId == 55 ? 1 : (rDs.mId == 56 ? 2 : 3);
        auto r6 = unv_ints(lines[5]);
        if (r6.size() < 6)
            throw ReadError("UNV: dataset " + std::to_string(rDs.mId) +
                            " record 6 needs six integers");
        const int analysis = static_cast<int>(r6[1]);
        res.mDataChar = static_cast<int>(r6[2]);
        result_type = static_cast<int>(r6[3]);
        data_type = static_cast<int>(r6[4]);
        ndv = static_cast<std::size_t>(std::max<std::int64_t>(r6[5], 0));
        k = 6;
        auto head = unv_ints(lines[k]);
        if (head.size() >= 2 && head[0] == 0 && head[1] == 0) {
            // meshio++ <= 15.5 wrote four zero lines here (8 ints, 2 ints, 2 x 6 reals), and
            // its "57" held element data; a valid record 7 has NINT >= 1.
            log::warn("UNV: dataset {} uses the meshio++ <= 15.5 header layout", rDs.mId);
            legacy = true;
            k += 4;
            res.mKey = UnvStepKey{analysis, 0, 0.0};
        } else {
            if (head.size() < 2)
                throw ReadError("UNV: dataset " + std::to_string(rDs.mId) +
                                " record 7 needs NINT and NRV");
            const std::size_t nint = static_cast<std::size_t>(std::max<std::int64_t>(head[0], 0));
            const std::size_t nrv = static_cast<std::size_t>(std::max<std::int64_t>(head[1], 0));
            std::vector<std::int64_t> ints(head.begin() + 2, head.end());
            ++k;
            while (ints.size() < nint && k < lines.size())
                for (std::int64_t v : unv_ints(lines[k++]))
                    ints.push_back(v);
            std::vector<double> reals = unv_take_reals(lines, k, nrv);
            res.mKey = unv_key_55(analysis, ints, reals);
        }
    }
    if (rDs.mId == 57 && legacy)
        res.mLocation = 2;
    if (unv_is_none(res.mName))
        res.mName = unv_result_type_name(result_type);

    if (res.mLocation != 1 && res.mLocation != 2) {
        log::warn(
            "UNV: skipping '{}' (dataset {}, location {}): only data at nodes and on "
            "elements is read",
            res.mName, rDs.mId, res.mLocation);
        return;
    }
    if (data_type != 1 && data_type != 2 && data_type != 4 && data_type != 5 && data_type != 6) {
        log::warn("UNV: skipping '{}' (dataset {}): unknown data type {}", res.mName, rDs.mId,
                  data_type);
        return;
    }
    if (ndv == 0)
        return;
    res.mComplex = data_type == 5 || data_type == 6;
    res.mNumComps = ndv;
    const std::size_t width = res.mComplex ? 2 : 1;
    bool warned_layers = false;
    while (k < lines.size()) {
        if (unv_strip(lines[k]).empty()) {
            ++k;
            continue;
        }
        auto rec = unv_ints(lines[k++]);
        std::size_t count = ndv;
        if (res.mLocation == 2 && rec.size() >= 2 && rec[1] > 0)
            count = static_cast<std::size_t>(rec[1]);  // NDVAL of this element
        auto vals = unv_take_reals(lines, k, count * width);
        if (count != ndv) {
            if (!warned_layers)
                log::warn(
                    "UNV: '{}' has {} values on an element but {} per component set; "
                    "keeping the first",
                    res.mName, count, ndv);
            warned_layers = true;
            vals.resize(ndv * width, kUnvNaN);
        }
        res.mValues.emplace_back(rec[0], std::move(vals));
    }
    rFile.mResults.push_back(std::move(res));
}

void unv_parse_function(const UnvDataset& rDs, UnvFile& rFile) {
    const auto& lines = rDs.mLines;
    if (lines.size() < 11)
        throw ReadError("UNV: dataset 58 has a truncated header");
    UnvFunction fn;
    // Record 6: Format(2(I5,I10),2(1X,10A1,I10,I4)); names may hold spaces, so slice.
    const std::string rec6(lines[5]);
    auto field = [&](std::size_t a, std::size_t n) {
        return a < rec6.size() ? unv_strip(std::string_view(rec6).substr(a, n))
                               : std::string_view();
    };
    auto int_field = [&](std::size_t a, std::size_t n) {
        auto f = field(a, n);
        return f.empty() ? std::int64_t{0} : unv_int(f);
    };
    if (rec6.size() >= 80) {
        fn.mType = static_cast<int>(int_field(0, 5));
        fn.mLoadCase = int_field(20, 10);
        fn.mRspNode = int_field(41, 10);
        fn.mRspDir = static_cast<int>(int_field(51, 4));
        fn.mRefNode = int_field(66, 10);
        fn.mRefDir = static_cast<int>(int_field(76, 4));
    } else {
        auto t = unv_split(rec6);
        if (t.size() < 10)
            throw ReadError("UNV: dataset 58 record 6 is malformed");
        fn.mType = static_cast<int>(unv_int(t[0]));
        fn.mLoadCase = unv_int(t[3]);
        fn.mRspNode = unv_int(t[5]);
        fn.mRspDir = static_cast<int>(unv_int(t[6]));
        fn.mRefNode = unv_int(t[8]);
        fn.mRefDir = static_cast<int>(unv_int(t[9]));
    }
    auto r7 = unv_split(lines[6]);
    if (r7.size() < 3)
        throw ReadError("UNV: dataset 58 record 7 is malformed");
    const int ord_type = static_cast<int>(unv_int(r7[0]));
    const std::size_t npts = static_cast<std::size_t>(std::max<std::int64_t>(unv_int(r7[1]), 0));
    const bool even = unv_int(r7[2]) == 1;
    std::vector<double> r7_reals;
    for (std::size_t n = 3; n < r7.size(); ++n)
        for (double v : unv_reals(r7[n]))
            r7_reals.push_back(v);
    r7_reals.resize(3, 0.0);
    auto r8 = unv_split(lines[7]);
    fn.mAbscissaType = r8.empty() ? 0 : static_cast<int>(unv_int(r8[0]));
    if (ord_type != 2 && ord_type != 4 && ord_type != 5 && ord_type != 6)
        throw ReadError("UNV: dataset 58 ordinate data type " + std::to_string(ord_type) +
                        " is not supported");
    fn.mComplex = ord_type == 5 || ord_type == 6;
    const bool dbl = ord_type == 4 || ord_type == 6;
    const std::size_t per = (fn.mComplex ? 2 : 1) + (even ? 0 : 1);
    const std::size_t count = npts * per;

    std::vector<double> values;
    if (rDs.mBinary) {
        if (rDs.mFpFormat != 2)
            throw ReadError("UNV: dataset 58b floating-point format " +
                            std::to_string(rDs.mFpFormat) + " is not IEEE 754");
        const bool swap = (rDs.mByteOrder == 2) != (std::endian::native == std::endian::big);
        const char* p = rDs.mBlob.data();
        const std::size_t size = rDs.mBlob.size();
        values.reserve(count);
        // Uneven double-precision data keeps a single-precision abscissa in ASCII (E13.5),
        // but in binary every value has the ordinate's width.
        const std::size_t w = dbl ? 8 : 4;
        if (size < count * w)
            throw ReadError("UNV: dataset 58b holds fewer bytes than its header declares");
        for (std::size_t n = 0; n < count; ++n) {
            char buf[8];
            std::memcpy(buf, p + n * w, w);
            if (swap)
                detail::bswap_inplace(buf, static_cast<int>(w));
            if (dbl) {
                double v;
                std::memcpy(&v, buf, 8);
                values.push_back(v);
            } else {
                float v;
                std::memcpy(&v, buf, 4);
                values.push_back(static_cast<double>(v));
            }
        }
    } else {
        std::size_t k = 11;
        values = unv_take_reals(lines, k, count);
    }
    fn.mX.resize(npts);
    fn.mRe.resize(npts);
    if (fn.mComplex)
        fn.mIm.resize(npts);
    for (std::size_t n = 0; n < npts; ++n) {
        const double* v = values.data() + n * per;
        std::size_t o = 0;
        if (even)
            fn.mX[n] = r7_reals[0] + static_cast<double>(n) * r7_reals[1];
        else
            fn.mX[n] = v[o++];
        fn.mRe[n] = v[o++];
        if (fn.mComplex)
            fn.mIm[n] = v[o++];
    }
    rFile.mFunctions.push_back(std::move(fn));
}

UnvFile unv_parse(const std::string& rPath) {
    detail::FileSource src(rPath);
    const std::string_view data = src.View();
    UnvFile file;
    for (const UnvDataset& ds : unv_split_datasets(data)) {
        bool quad = false;
        if (ds.mId == 2411 || ds.mId == 781 || ds.mId == 15)
            unv_parse_nodes(ds, file);
        else if (ds.mId == 2412 || ds.mId == 780)
            unv_parse_elements(ds, file);
        else if (unv_group_layout(ds.mId, quad))
            unv_parse_groups(ds, quad, file);
        else if (ds.mId == 164)
            unv_parse_units(ds, file);
        else if (ds.mId == 2420)
            unv_parse_cs(ds, file);
        else if (ds.mId == 2414 || ds.mId == 55 || ds.mId == 56 || ds.mId == 57)
            unv_parse_result(ds, file);
        else if (ds.mId == 58)
            unv_parse_function(ds, file);
        else
            log::debug("UNV: dataset {} ignored", ds.mId);
    }
    return file;
}

// ---------------------------------------------------------------------------
// Steps
// ---------------------------------------------------------------------------

// Dataset-58 functions of one kind: (type, load case, reference node, reference
// direction) and one abscissa grid.
struct UnvFunctionGroup {
    std::string mName;
    int mAnalysis = 0;
    std::vector<double> mX;
    std::vector<std::size_t> mFunctions;
};

struct UnvStep {
    UnvStepKey mKey;
    std::vector<std::size_t> mResults;
    std::vector<std::pair<std::size_t, std::size_t>> mSamples;  // (group, sample index)
};

const char* unv_dir_name(int Dir) {
    static const char* const kNames[] = {"", "x", "y", "z", "rx", "ry", "rz"};
    const int a = std::abs(Dir);
    return a <= 6 ? kNames[a] : "";
}

std::vector<UnvFunctionGroup> unv_group_functions(const UnvFile& rFile) {
    std::vector<UnvFunctionGroup> groups;
    struct Kind {
        int mType;
        std::int64_t mLoadCase, mRefNode;
        int mRefAxis;
        std::vector<std::size_t> mGroups;
    };
    std::vector<Kind> kinds;
    for (std::size_t f = 0; f < rFile.mFunctions.size(); ++f) {
        const UnvFunction& fn = rFile.mFunctions[f];
        const int ref_axis = std::abs(fn.mRefDir);
        auto kit = std::find_if(kinds.begin(), kinds.end(), [&](const Kind& k) {
            return k.mType == fn.mType && k.mLoadCase == fn.mLoadCase &&
                   k.mRefNode == fn.mRefNode && k.mRefAxis == ref_axis;
        });
        if (kit == kinds.end()) {
            kinds.push_back({fn.mType, fn.mLoadCase, fn.mRefNode, ref_axis, {}});
            kit = kinds.end() - 1;
        }
        std::size_t g = groups.size();
        for (std::size_t cand : kit->mGroups)
            if (groups[cand].mX == fn.mX) {
                g = cand;
                break;
            }
        if (g == groups.size()) {
            if (!kit->mGroups.empty())
                log::warn(
                    "UNV: {} functions of one kind do not share an abscissa grid; each grid "
                    "becomes its own steps",
                    unv_function_name(fn.mType));
            UnvFunctionGroup grp;
            grp.mX = fn.mX;
            grp.mAnalysis = fn.mAbscissaType == 18 ? 5 : (fn.mAbscissaType == 17 ? 4 : 0);
            groups.push_back(std::move(grp));
            kit->mGroups.push_back(g);
        }
        groups[g].mFunctions.push_back(f);
    }
    // Names: the function type, qualified by the reference only when two kinds share it.
    std::map<std::string, int> type_count;
    for (const Kind& k : kinds)
        ++type_count[unv_function_name(k.mType)];
    std::unordered_set<std::string> used;
    for (const Kind& k : kinds) {
        std::string base = unv_function_name(k.mType);
        if (type_count[base] > 1)
            base += "_ref" + std::to_string(k.mRefNode) + unv_dir_name(k.mRefAxis);
        for (std::size_t n = 0; n < k.mGroups.size(); ++n) {
            std::string name = base;
            for (int s = 2; used.count(name); ++s)
                name = base + "_" + std::to_string(s);
            used.insert(name);
            groups[k.mGroups[n]].mName = name;
        }
    }
    return groups;
}

std::vector<UnvStep> unv_steps(const UnvFile& rFile, const std::vector<UnvFunctionGroup>& rGroups) {
    std::vector<UnvStep> steps;
    auto step_of = [&](const UnvStepKey& rKey) -> UnvStep& {
        for (UnvStep& s : steps)
            if (s.mKey == rKey)
                return s;
        steps.push_back(UnvStep{rKey, {}, {}});
        return steps.back();
    };
    for (std::size_t r = 0; r < rFile.mResults.size(); ++r)
        step_of(rFile.mResults[r].mKey).mResults.push_back(r);
    for (std::size_t g = 0; g < rGroups.size(); ++g)
        for (std::size_t n = 0; n < rGroups[g].mX.size(); ++n) {
            UnvStepKey key{rGroups[g].mAnalysis, static_cast<std::int64_t>(n + 1),
                           rGroups[g].mX[n]};
            step_of(key).mSamples.emplace_back(g, n);
        }
    return steps;
}

NDArray unv_scalar(DType Type, double Value) {
    NDArray out(Type, {std::size_t{1}});
    if (Type == DType::Float64)
        *out.As<double>() = Value;
    else
        *out.As<std::int64_t>() = static_cast<std::int64_t>(Value);
    return out;
}

NDArray unv_nan_array(std::size_t N, std::size_t NumComps) {
    NDArray out =
        NumComps == 1 ? NDArray(DType::Float64, {N}) : NDArray(DType::Float64, {N, NumComps});
    std::fill(out.As<double>(), out.As<double>() + N * NumComps, kUnvNaN);
    return out;
}

// ---------------------------------------------------------------------------
// Mesh assembly
// ---------------------------------------------------------------------------

Mesh unv_build(const UnvFile& rFile, UnvInfo& rInfo, const ReadOptions& rOpts,
               std::vector<double>* pTimes) {
    Mesh mesh;

    // Points, moved into the global system where a node was defined in a local one.
    std::unordered_map<std::int64_t, std::int64_t> node_index;
    std::size_t dim =
        rFile.mNodes.empty() ? 3 : std::max<std::size_t>(rFile.mNodes[0].mNumCoords, 1);
    std::vector<std::array<double, 3>> coords;
    coords.reserve(rFile.mNodes.size());
    std::set<std::int64_t> warned_cs;
    for (const UnvNode& node : rFile.mNodes) {
        std::array<double, 3> x = {node.mX[0], node.mX[1], node.mX[2]};
        auto cit = rFile.mCs.find(node.mCs);
        if (cit != rFile.mCs.end()) {
            const UnvCs& cs = cit->second;
            if (cs.mType == 0) {
                std::array<double, 3> g;
                for (int r = 0; r < 3; ++r)
                    g[r] =
                        cs.mM[r][0] * x[0] + cs.mM[r][1] * x[1] + cs.mM[r][2] * x[2] + cs.mM[3][r];
                x = g;
                if (dim < 3 && (x[2] != 0.0 || (dim < 2 && x[1] != 0.0)))
                    dim = 3;
            } else if (warned_cs.insert(node.mCs).second) {
                log::warn(
                    "UNV: coordinate system {} is {}; its nodes are left in local "
                    "coordinates",
                    node.mCs, cs.mType == 1 ? "cylindrical" : "spherical");
            }
        } else if (node.mCs > 1 && warned_cs.insert(node.mCs).second) {
            log::warn("UNV: coordinate system {} is not defined; its nodes are taken as global",
                      node.mCs);
        }
        auto [it, fresh] =
            node_index.emplace(node.mLabel, static_cast<std::int64_t>(coords.size()));
        if (fresh) {
            coords.push_back(x);
        } else {
            log::warn("UNV: node {} is defined twice; the last definition wins", node.mLabel);
            coords[static_cast<std::size_t>(it->second)] = x;
        }
    }
    const std::size_t np = coords.size();
    NDArray pts(DType::Float64, {np, dim});
    for (std::size_t r = 0; r < np; ++r)
        for (std::size_t c = 0; c < dim; ++c)
            pts.As<double>()[r * dim + c] = coords[r][c];
    mesh.AssignPoints(std::move(pts));

    // Cells, one block per type in first-appearance order.
    struct Block {
        std::string mType;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mPid, mMid;
        std::size_t mNumNodes = 0;
    };
    std::vector<Block> blocks;
    std::unordered_map<std::string, std::size_t> block_of;
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> elem_ref;
    std::set<std::pair<int, std::size_t>> warned_fe;
    for (const UnvElement& el : rFile.mElements) {
        const std::string type = unv_cell_type(el.mFeId, el.mNodes.size());
        if (type.empty()) {
            if (warned_fe.emplace(el.mFeId, el.mNodes.size()).second)
                log::warn(
                    "UNV: FE descriptor {} with {} nodes is not supported; skipping those "
                    "elements",
                    el.mFeId, el.mNodes.size());
            continue;
        }
        auto [bit, fresh] = block_of.emplace(type, blocks.size());
        if (fresh)
            blocks.push_back(Block{type, {}, {}, {}, el.mNodes.size()});
        Block& b = blocks[bit->second];
        const std::vector<int>* perm = unv_perm(type);
        const std::size_t base = b.mConn.size();
        b.mConn.resize(base + b.mNumNodes);
        for (std::size_t j = 0; j < b.mNumNodes; ++j) {
            auto nit = node_index.find(el.mNodes[j]);
            if (nit == node_index.end())
                throw ReadError("UNV: element " + std::to_string(el.mLabel) +
                                " references undefined node " + std::to_string(el.mNodes[j]));
            b.mConn[base + (perm ? static_cast<std::size_t>((*perm)[j]) : j)] = nit->second;
        }
        elem_ref[el.mLabel] = {bit->second, b.mPid.size()};
        b.mPid.push_back(el.mPid);
        b.mMid.push_back(el.mMid);
    }
    std::vector<std::size_t> block_sizes;
    std::vector<NDArray> pids, mids;
    for (Block& b : blocks) {
        const std::size_t ne = b.mPid.size();
        NDArray conn(DType::Int64, {ne, b.mNumNodes});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
        NDArray pid(DType::Int64, {ne});
        std::copy(b.mPid.begin(), b.mPid.end(), pid.As<std::int64_t>());
        pids.push_back(std::move(pid));
        NDArray mid(DType::Int64, {ne});
        std::copy(b.mMid.begin(), b.mMid.end(), mid.As<std::int64_t>());
        mids.push_back(std::move(mid));
        block_sizes.push_back(ne);
    }
    if (!blocks.empty()) {
        mesh.AddCellData("unv:pid", std::move(pids));
        mesh.AddCellData("unv:mid", std::move(mids));
    }

    // Groups: node members -> a Point region, element members -> a Cell region.
    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    for (const UnvGroup& g : rFile.mGroups) {
        std::size_t missing = 0;
        std::vector<std::int64_t> pts_idx;
        for (std::int64_t label : g.mNodes) {
            auto it = node_index.find(label);
            if (it == node_index.end())
                ++missing;
            else
                pts_idx.push_back(it->second);
        }
        std::vector<std::int64_t> cells;
        std::vector<std::vector<std::int64_t>> per_block(blocks.size());
        int cdim = -1;
        for (std::int64_t label : g.mElements) {
            auto it = elem_ref.find(label);
            if (it == elem_ref.end()) {
                ++missing;
                continue;
            }
            const auto [b, row] = it->second;
            cells.push_back(detail::block_row_to_global(bases, b, static_cast<std::int64_t>(row)));
            per_block[b].push_back(static_cast<std::int64_t>(row));
            cdim = std::max(cdim, cell_type_dimension(cell_type_from_name(blocks[b].mType)));
        }
        if (missing)
            log::warn("UNV: group '{}' names {} entities that were not read; they are dropped",
                      g.mName, missing);
        auto entries = [](const std::vector<std::int64_t>& rV) {
            NDArray out(DType::Int64, {rV.size()});
            std::copy(rV.begin(), rV.end(), out.As<std::int64_t>());
            return out;
        };
        if (!g.mNodes.empty()) {
            mesh.AddRegion(Region(g.mName, RegionKind::Point, -1, g.mNumber, entries(pts_idx)));
            rInfo.mPointSets[g.mName] = pts_idx;
        }
        if (!g.mElements.empty() || g.mNodes.empty()) {
            mesh.AddRegion(Region(g.mName, RegionKind::Cell, cdim, g.mNumber, entries(cells)));
            rInfo.mCellSets[g.mName] = std::move(per_block);
        }
    }

    if (rFile.mHasUnits) {
        mesh.AddFieldData("unv:units",
                          unv_scalar(DType::Int64, static_cast<double>(rFile.mUnitsCode)));
        NDArray factors(DType::Float64, {std::size_t{4}});
        std::copy(rFile.mUnitFactors.begin(), rFile.mUnitFactors.end(), factors.As<double>());
        mesh.AddFieldData("unv:unit_factors", std::move(factors));
    }

    // Results: one step of the sequence.
    if (!rFile.mFunctions.empty() && np == 0)
        throw ReadError(
            "UNV: the file holds dataset-58 functions but no nodes (datasets "
            "15/781/2411) to attach them to");
    const std::vector<UnvFunctionGroup> groups = unv_group_functions(rFile);
    const std::vector<UnvStep> steps = unv_steps(rFile, groups);
    if (pTimes) {
        pTimes->clear();
        for (const UnvStep& s : steps)
            pTimes->push_back(s.mKey.mValue);
    }
    if (steps.empty()) {
        rOpts.ResolveTimeStep(0);
        return mesh;
    }
    const UnvStep& step = steps[rOpts.ResolveTimeStep(steps.size())];
    if (steps.size() > 1 || step.mKey.mAnalysis != 0) {
        mesh.AddFieldData(kSequenceTimeKey, unv_scalar(DType::Float64, step.mKey.mValue));
        mesh.AddFieldData("unv:analysis",
                          unv_scalar(DType::Int64, static_cast<double>(step.mKey.mAnalysis)));
        mesh.AddFieldData("unv:step", unv_scalar(DType::Int64, static_cast<double>(step.mKey.mId)));
    }
    if (!rOpts.WantsAnyData())
        return mesh;

    std::unordered_set<std::string> used = {"unv:pid", "unv:mid"};
    auto unique = [&](const std::string& rBase) {
        std::string name = rBase;
        for (int k = 2; used.count(name); ++k)
            name = rBase + "_" + std::to_string(k);
        used.insert(name);
        return name;
    };

    for (std::size_t r : step.mResults) {
        const UnvResult& res = rFile.mResults[r];
        const std::size_t nc = res.mNumComps;
        const std::string base = unique(res.mName);
        const int nparts = res.mComplex ? 2 : 1;
        for (int part = 0; part < nparts; ++part) {
            const std::string name =
                res.mComplex ? (part == 0 ? base + "_real" : base + "_imag") : base;
            if (!rOpts.WantsArray(name))
                continue;
            std::vector<double> buf(nc);
            auto fill = [&](const std::vector<double>& rVals) {
                for (std::size_t c = 0; c < nc; ++c)
                    buf[c] = res.mComplex ? rVals[2 * c + part] : rVals[c];
                unv_tensor_to_meshio(res.mDataChar, nc, buf.data());
            };
            if (res.mLocation == 1) {
                NDArray arr = unv_nan_array(np, nc);
                for (const auto& [label, vals] : res.mValues) {
                    auto it = node_index.find(label);
                    if (it == node_index.end())
                        continue;
                    fill(vals);
                    std::copy(buf.begin(), buf.end(),
                              arr.As<double>() + static_cast<std::size_t>(it->second) * nc);
                }
                mesh.AddPointData(name, std::move(arr));
            } else if (!block_sizes.empty()) {
                std::vector<NDArray> arrs;
                for (std::size_t ne : block_sizes)
                    arrs.push_back(unv_nan_array(ne, nc));
                for (const auto& [label, vals] : res.mValues) {
                    auto it = elem_ref.find(label);
                    if (it == elem_ref.end())
                        continue;
                    fill(vals);
                    std::copy(buf.begin(), buf.end(),
                              arrs[it->second.first].As<double>() + it->second.second * nc);
                }
                mesh.AddCellData(name, std::move(arrs));
            }
        }
    }

    // Dataset-58 samples of this step, one set of arrays per function group.
    std::set<std::int64_t> warned_nodes;
    for (const auto& [g, n] : step.mSamples) {
        const UnvFunctionGroup& grp = groups[g];
        bool has_trans = false, has_rot = false, has_scalar = false, is_complex = false;
        for (std::size_t f : grp.mFunctions) {
            const int a = std::abs(rFile.mFunctions[f].mRspDir);
            has_trans = has_trans || (a >= 1 && a <= 3);
            has_rot = has_rot || (a >= 4 && a <= 6);
            has_scalar = has_scalar || a == 0 || a > 6;
            is_complex = is_complex || rFile.mFunctions[f].mComplex;
        }
        struct Target {
            std::string mName;
            std::size_t mNumComps;
            int mKind;  // 0 scalar, 1 translation, 2 rotation
        };
        std::vector<Target> targets;
        const std::string base = unique(grp.mName);
        if (has_trans)
            targets.push_back({base, 3, 1});
        if (has_rot)
            targets.push_back({unique(grp.mName + "_rot"), 3, 2});
        if (has_scalar)
            targets.push_back({has_trans ? unique(grp.mName + "_scalar") : base, 1, 0});
        for (const Target& t : targets) {
            const int nparts = is_complex ? 2 : 1;
            for (int part = 0; part < nparts; ++part) {
                const std::string name =
                    is_complex ? (part == 0 ? t.mName + "_real" : t.mName + "_imag") : t.mName;
                if (!rOpts.WantsArray(name))
                    continue;
                NDArray arr = unv_nan_array(np, t.mNumComps);
                for (std::size_t f : grp.mFunctions) {
                    const UnvFunction& fn = rFile.mFunctions[f];
                    const int a = std::abs(fn.mRspDir);
                    const int kind = (a >= 1 && a <= 3) ? 1 : ((a >= 4 && a <= 6) ? 2 : 0);
                    if (kind != t.mKind)
                        continue;
                    auto it = node_index.find(fn.mRspNode);
                    if (it == node_index.end()) {
                        if (warned_nodes.insert(fn.mRspNode).second)
                            log::warn(
                                "UNV: dataset-58 response node {} is not defined; its "
                                "functions are dropped",
                                fn.mRspNode);
                        continue;
                    }
                    const double sign =
                        (fn.mRspDir < 0 ? -1.0 : 1.0) * (fn.mRefDir < 0 ? -1.0 : 1.0);
                    const double v = part == 0 ? fn.mRe[n] : (fn.mComplex ? fn.mIm[n] : 0.0);
                    const std::size_t comp = kind == 0 ? 0 : static_cast<std::size_t>((a - 1) % 3);
                    arr.As<double>()[static_cast<std::size_t>(it->second) * t.mNumComps + comp] =
                        sign * v;
                }
                mesh.AddPointData(name, std::move(arr));
            }
        }
    }
    return mesh;
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

std::int64_t unv_field_int(const Mesh& rMesh, const char* pName, std::int64_t Default) {
    if (!rMesh.HasFieldData(pName))
        return Default;
    const NDArray& a = rMesh.FieldData(pName);
    return a.Size() ? detail::read_int(a, 0) : Default;
}

double unv_field_real(const Mesh& rMesh, const char* pName, double Default) {
    if (!rMesh.HasFieldData(pName))
        return Default;
    const NDArray& a = rMesh.FieldData(pName);
    return a.Size() ? detail::read_double(a, 0) : Default;
}

void unv_write_reals(std::ofstream& rF, const double* pValues, std::size_t N, const char* pFmt,
                     std::size_t PerLine) {
    char buf[64];
    for (std::size_t i = 0; i < N; ++i) {
        detail::snprintf_c(buf, sizeof(buf), pFmt, pValues[i]);
        rF << buf;
        if ((i + 1) % PerLine == 0 || i + 1 == N)
            rF << "\n";
    }
}

}  // namespace

Mesh read_unv(const std::string& rPath, UnvInfo& rInfo, const ReadOptions& rOpts) {
    const UnvFile file = unv_parse(rPath);
    return unv_build(file, rInfo, rOpts, nullptr);
}

Mesh read_unv(const std::string& rPath, const ReadOptions& rOpts) {
    UnvInfo info;
    return read_unv(rPath, info, rOpts);
}

Mesh read_unv(const std::string& rPath, UnvInfo& rInfo) {
    return read_unv(rPath, rInfo, ReadOptions{});
}

Mesh read_unv(const std::string& rPath) {
    UnvInfo info;
    return read_unv(rPath, info, ReadOptions{});
}

MeshMetadata read_unv_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    const UnvFile file = unv_parse(rPath);
    UnvInfo info;
    std::vector<double> times;
    MeshMetadata meta = metadata_from_mesh(unv_build(file, info, ReadOptions{}, &times));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "unv";
    meta.mTimeValues = std::move(times);
    return meta;
}

void write_unv(const std::string& rPath, const Mesh& rMesh, bool code_aster, int node_dataset) {
    UnvInfo info;
    write_unv(rPath, rMesh, info, code_aster, node_dataset);
}

void write_unv(const std::string& rPath, const Mesh& rMesh, const UnvInfo& rInfo, bool code_aster,
               int node_dataset) {
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t np = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();

    if (node_dataset != 2411 && node_dataset != 781)
        node_dataset = 2411;

    char buf[160];

    // 164 units, when the mesh carries them.
    if (rMesh.HasFieldData("unv:units") && rMesh.HasFieldData("unv:unit_factors")) {
        const NDArray& fac = rMesh.FieldData("unv:unit_factors");
        double v[4] = {1.0, 1.0, 1.0, 0.0};
        for (std::size_t i = 0; i < 4 && i < fac.Size(); ++i)
            v[i] = detail::read_double(fac, i);
        std::snprintf(buf, sizeof(buf), "    -1\n   164\n%10lld%20s%10d\n",
                      static_cast<long long>(unv_field_int(rMesh, "unv:units", 1)), "", 2);
        f << buf;
        for (int i = 0; i < 4; ++i) {
            detail::snprintf_c(buf, sizeof(buf), "%25.16E", v[i]);
            f << buf;
            if (i == 2 || i == 3)
                f << "\n";
        }
        f << "    -1\n";
    }

    // 2411 / 781 nodes
    std::snprintf(buf, sizeof(buf), "    -1\n%6d\n", node_dataset);
    f << buf;
    for (std::size_t k = 0; k < np; ++k) {
        std::snprintf(buf, sizeof(buf), "%10zu%10d%10d%10d\n", k + 1, 1, 1, 11);
        f << buf;
        for (int c = 0; c < 3; ++c) {
            double v = c < static_cast<int>(pdim) ? detail::read_double(points, k * pdim + c) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), "%25.16E", v);
            f << buf;
        }
        f << "\n";
    }
    f << "    -1\n";

    // 2412 elements
    f << "    -1\n  2412\n";
    const bool has_pid = rMesh.HasCellData("unv:pid");
    const bool has_mid = rMesh.HasCellData("unv:mid");
    std::int64_t label = 0;
    const std::size_t nblocks = rMesh.NumCellBlocks();
    // per-block 1-based element labels (empty for skipped blocks) so groups and element
    // data can resolve (block, local) -> label.
    std::vector<std::vector<std::int64_t>> block_labels(nblocks);
    for (std::size_t bi = 0; bi < nblocks; ++bi) {
        const auto cb = rMesh.Cells(bi);
        int desc;
        bool beam;
        if (!unv_descriptor(cb.Type(), desc, beam)) {
            log::warn("UNV does not support '{}' cells. Skipping.", cb.Type());
            detail::provenance_note(
                "cells-dropped",
                "cell block(s) of type " + std::string(cb.Type()) + " have no UNV equivalent");
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::vector<int>* perm = unv_perm(cb.Type());
        const std::size_t ncols = detail::cols(conn);
        const std::size_t nrows = cb.NumCells();
        block_labels[bi].reserve(nrows);
        const NDArray* pid = (has_pid && bi < rMesh.CellDataNumBlocks("unv:pid"))
                                 ? &rMesh.CellData("unv:pid", bi)
                                 : nullptr;
        const NDArray* mid = (has_mid && bi < rMesh.CellDataNumBlocks("unv:mid"))
                                 ? &rMesh.CellData("unv:mid", bi)
                                 : nullptr;
        std::vector<std::int64_t> unv(ncols);
        for (std::size_t r = 0; r < nrows; ++r) {
            ++label;
            block_labels[bi].push_back(label);
            const std::int64_t pval = pid ? detail::read_int(*pid, r) : 1;
            const std::int64_t mval = mid ? detail::read_int(*mid, r) : pval;
            std::snprintf(buf, sizeof(buf), "%10lld%10d%10lld%10lld%10d%10zu\n",
                          static_cast<long long>(label), desc, static_cast<long long>(pval),
                          static_cast<long long>(mval), 11, ncols);
            f << buf;
            if (beam)
                f << "         0         1         1\n";
            // meshio -> UNV order, 1-based, 8 per line
            for (std::size_t j = 0; j < ncols; ++j)
                unv[j] = detail::read_int(conn, r * ncols + (perm ? (*perm)[j] : j)) + 1;
            for (std::size_t j = 0; j < ncols; ++j) {
                std::snprintf(buf, sizeof(buf), "%10lld", static_cast<long long>(unv[j]));
                f << buf;
                if ((j + 1) % 8 == 0 || j + 1 == ncols)
                    f << "\n";
            }
        }
    }
    f << "    -1\n";

    // 2467 permanent groups from the mesh's regions (and any UnvInfo sets, which replace
    // a region of the same name and kind). A point and a cell region sharing a name are
    // one group.
    struct OutGroup {
        std::int64_t mNumber = -1;
        std::vector<std::int64_t> mNodes, mElements;
        bool mHasNodes = false, mHasElements = false;
    };
    std::map<std::string, OutGroup> out_groups;
    std::vector<std::string> order;
    auto group = [&](const std::string& rName) -> OutGroup& {
        auto [it, fresh] = out_groups.emplace(rName, OutGroup{});
        if (fresh)
            order.push_back(rName);
        return it->second;
    };
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::size_t side_dropped = 0, cells_dropped = 0;
    for (std::size_t ri = 0; ri < rMesh.NumRegions(); ++ri) {
        const meshioplusplus::Region& reg = rMesh.Region(ri);
        if (reg.mKind == RegionKind::Side) {
            ++side_dropped;
            continue;
        }
        if (reg.mKind == RegionKind::Point && rInfo.mPointSets.count(reg.mName))
            continue;
        if (reg.mKind == RegionKind::Cell && rInfo.mCellSets.count(reg.mName))
            continue;
        OutGroup& g = group(reg.mName);
        if (g.mNumber < 0 && reg.mTag > 0)
            g.mNumber = reg.mTag;
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i) {
            if (reg.mKind == RegionKind::Point) {
                g.mHasNodes = true;
                g.mNodes.push_back(e[i] + 1);
            } else {
                g.mHasElements = true;
                auto [b, row] = detail::global_to_block_row(bases, e[i]);
                if (b >= nblocks || static_cast<std::size_t>(row) >= block_labels[b].size()) {
                    ++cells_dropped;
                    continue;
                }
                g.mElements.push_back(block_labels[b][static_cast<std::size_t>(row)]);
            }
        }
        if (reg.mKind == RegionKind::Cell)
            g.mHasElements = true;
    }
    for (const auto& [name, idx] : rInfo.mPointSets) {
        OutGroup& g = group(name);
        g.mHasNodes = true;
        for (std::int64_t i : idx)
            g.mNodes.push_back(i + 1);
    }
    for (const auto& [name, sets] : rInfo.mCellSets) {
        OutGroup& g = group(name);
        g.mHasElements = true;
        for (std::size_t bi = 0; bi < sets.size() && bi < block_labels.size(); ++bi)
            for (std::int64_t local : sets[bi]) {
                if (local >= 0 && local < static_cast<std::int64_t>(block_labels[bi].size()))
                    g.mElements.push_back(block_labels[bi][static_cast<std::size_t>(local)]);
                else
                    ++cells_dropped;
            }
    }
    if (side_dropped) {
        log::warn("UNV has no facet groups; {} side region(s) dropped", side_dropped);
        detail::provenance_note("regions-dropped",
                                std::to_string(side_dropped) + " side region(s) have no UNV group");
    }
    if (cells_dropped)
        log::warn("UNV: {} group member(s) sit in cell blocks UNV cannot hold; dropped",
                  cells_dropped);
    if (!order.empty()) {
        std::set<std::int64_t> taken;
        for (const std::string& name : order) {
            OutGroup& g = out_groups[name];
            if (g.mNumber > 0 && !taken.insert(g.mNumber).second)
                g.mNumber = -1;
        }
        std::int64_t next = 1;
        for (const std::string& name : order) {
            OutGroup& g = out_groups[name];
            if (g.mNumber > 0)
                continue;
            while (taken.count(next))
                ++next;
            g.mNumber = next;
            taken.insert(next);
        }
        std::vector<std::string> by_number = order;
        std::stable_sort(by_number.begin(), by_number.end(),
                         [&](const std::string& a, const std::string& b) {
                             return out_groups[a].mNumber < out_groups[b].mNumber;
                         });
        f << "    -1\n  2467\n";
        for (const std::string& name : by_number) {
            const OutGroup& g = out_groups[name];
            std::snprintf(buf, sizeof(buf), "%10lld%10d%10d%10d%10d%10d%10d%10zu\n",
                          static_cast<long long>(g.mNumber), 0, 0, 0, 0, 0, 0,
                          g.mNodes.size() + g.mElements.size());
            f << buf << name << "\n";
            std::size_t col = 0;
            auto entity = [&](int type, std::int64_t tag) {
                std::snprintf(buf, sizeof(buf), "%10d%10lld%10d%10d", type,
                              static_cast<long long>(tag), 0, 0);
                f << buf;
                if (++col == 2) {
                    f << "\n";
                    col = 0;
                }
            };
            for (std::int64_t t : g.mNodes)
                entity(7, t);
            for (std::int64_t t : g.mElements)
                entity(8, t);
            if (col != 0)
                f << "\n";
        }
        f << "    -1\n";
    }

    // Results: point_data -> data at nodes, cell_data -> data on elements, as one step
    // described by the mesh's `unv:analysis` / `unv:step` / `meshio:time` field data.
    const int analysis = static_cast<int>(unv_field_int(rMesh, "unv:analysis", 0));
    const std::int64_t step_id = unv_field_int(rMesh, "unv:step", 0);
    const double time = unv_field_real(rMesh, kSequenceTimeKey, 0.0);
    auto write_field = [&](int field_id, const std::string& rName, int location, std::size_t nc,
                           std::vector<std::int64_t> rLabels, std::vector<double> rFlat) {
        // An entity with no value (NaN everywhere, as the reader pads) gets no record.
        {
            std::vector<std::int64_t> keep_labels;
            std::vector<double> keep_flat;
            for (std::size_t r = 0; r < rLabels.size(); ++r) {
                const double* row = rFlat.data() + r * nc;
                if (std::all_of(row, row + nc, [](double v) { return std::isnan(v); }))
                    continue;
                keep_labels.push_back(rLabels[r]);
                keep_flat.insert(keep_flat.end(), row, row + nc);
            }
            if (keep_labels.size() != rLabels.size()) {
                rLabels.swap(keep_labels);
                rFlat.swap(keep_flat);
            }
        }
        const int ch = unv_data_char(nc);
        if (ch == 4 || ch == 5)
            for (std::size_t r = 0; r < rLabels.size(); ++r)
                unv_tensor_from_meshio(nc, rFlat.data() + r * nc);
        if (code_aster) {
            // Legacy 55 (nodes) / 56 (elements), single precision E13.5 as the spec reads.
            std::snprintf(buf, sizeof(buf), "    -1\n%6d\n", location == 1 ? 55 : 56);
            f << buf;
            f << rName << "\n";
            for (int i = 0; i < 4; ++i)
                f << "NONE\n";
            std::snprintf(buf, sizeof(buf), "%10d%10d%10d%10d%10d%10zu\n", 1, analysis, ch, 0, 2,
                          nc);
            f << buf;
            std::vector<std::int64_t> ints;
            std::vector<double> reals;
            switch (analysis) {
                case 2:
                    ints = {1, step_id};
                    reals = {time, 0.0, 0.0, 0.0};
                    break;
                case 3:
                case 7:
                    ints = {1, step_id};
                    reals = {0.0, 2.0 * kUnvPi * time, 0.0, 0.0, 0.0, 0.0};
                    break;
                case 4:
                case 5:
                    ints = {1, step_id};
                    reals = {time};
                    break;
                default:
                    ints = {step_id};
                    reals = {time};
                    break;
            }
            std::snprintf(buf, sizeof(buf), "%10zu%10zu", ints.size(), reals.size());
            f << buf;
            for (std::int64_t v : ints) {
                std::snprintf(buf, sizeof(buf), "%10lld", static_cast<long long>(v));
                f << buf;
            }
            f << "\n";
            unv_write_reals(f, reals.data(), reals.size(), "%13.5E", 6);
            for (std::size_t r = 0; r < rLabels.size(); ++r) {
                if (location == 1)
                    std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(rLabels[r]));
                else
                    std::snprintf(buf, sizeof(buf), "%10lld%10zu\n",
                                  static_cast<long long>(rLabels[r]), nc);
                f << buf;
                unv_write_reals(f, rFlat.data() + r * nc, nc, "%13.5E", 6);
            }
        } else {
            f << "    -1\n  2414\n";
            std::snprintf(buf, sizeof(buf), "%10d\n", field_id);
            f << buf << rName << "\n";
            std::snprintf(buf, sizeof(buf), "%10d\n", location);
            f << buf;
            for (int i = 0; i < 5; ++i)
                f << "NONE\n";
            std::snprintf(buf, sizeof(buf), "%10d%10d%10d%10d%10d%10zu\n", 1, analysis, ch, 0, 4,
                          nc);
            f << buf;
            std::int64_t r10[8] = {0, 0, 0, 0, 0, 0, 0, 0};
            double r12[12] = {0.0};
            switch (analysis) {
                case 2:
                    r10[5] = step_id;
                    r12[1] = time;
                    break;
                case 3:
                case 7:
                    r10[5] = step_id;
                    r12[7] = 2.0 * kUnvPi * time;
                    break;
                case 4:
                case 9:
                    r10[6] = step_id;
                    r12[0] = time;
                    break;
                case 5:
                    r10[7] = step_id;
                    r12[1] = time;
                    break;
                case 6:
                    r10[5] = step_id;
                    r12[2] = time;
                    break;
                default:
                    r10[4] = step_id;
                    r12[0] = time;
                    break;
            }
            for (std::int64_t v : r10) {
                std::snprintf(buf, sizeof(buf), "%10lld", static_cast<long long>(v));
                f << buf;
            }
            f << "\n         0         0\n";
            unv_write_reals(f, r12, 6, "%13.5E", 6);
            unv_write_reals(f, r12 + 6, 6, "%13.5E", 6);
            for (std::size_t r = 0; r < rLabels.size(); ++r) {
                if (location == 1)
                    std::snprintf(buf, sizeof(buf), "%10lld\n", static_cast<long long>(rLabels[r]));
                else
                    std::snprintf(buf, sizeof(buf), "%10lld%10zu\n",
                                  static_cast<long long>(rLabels[r]), nc);
                f << buf;
                // One line per entity, whatever its length: pyuff pairs each element
                // record with exactly one value line.
                unv_write_reals(f, rFlat.data() + r * nc, nc, "%20.12E", nc);
            }
        }
        f << "    -1\n";
    };

    int field_id = 0;
    std::vector<std::int64_t> node_labels(np);
    for (std::size_t k = 0; k < np; ++k)
        node_labels[k] = static_cast<std::int64_t>(k + 1);
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& arr = rMesh.PointData(name);
        const std::size_t nc = np ? arr.Size() / np : 0;
        if (nc == 0)
            continue;
        std::vector<double> flat(np * nc);
        for (std::size_t i = 0; i < np * nc; ++i)
            flat[i] = detail::read_double(arr, i);
        write_field(++field_id, name, 1, nc, node_labels, flat);
    }
    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "unv:pid" || name == "unv:mid")
            continue;
        std::vector<std::int64_t> labels;
        std::vector<double> flat;
        std::size_t nc = 0;
        for (std::size_t bi = 0; bi < nblocks && bi < rMesh.CellDataNumBlocks(name); ++bi) {
            if (block_labels[bi].empty())
                continue;
            const NDArray& blk = rMesh.CellData(name, bi);
            const std::size_t ne = block_labels[bi].size();
            const std::size_t bnc = ne ? blk.Size() / ne : 0;
            if (bnc == 0)
                continue;
            if (nc == 0)
                nc = bnc;
            if (bnc != nc)
                continue;
            for (std::size_t r = 0; r < ne; ++r) {
                labels.push_back(block_labels[bi][r]);
                for (std::size_t c = 0; c < nc; ++c)
                    flat.push_back(detail::read_double(blk, r * nc + c));
            }
        }
        if (nc == 0 || labels.empty())
            continue;
        write_field(++field_id, name, 2, nc, labels, flat);
    }
}

}  // namespace meshioplusplus
