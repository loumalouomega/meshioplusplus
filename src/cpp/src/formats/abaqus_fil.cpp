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
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/abaqus_fil.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/abaqus_types.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"
#include "abaqus_face.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

// --- words and records ---------------------------------------------------------------

// One 8-byte word. ASCII items carry their tag ('I', 'D', 'A'); binary words
// are raw ('B') and are read as whatever the record layout says they are.
struct FilWord {
    char mTag;
    std::uint64_t mRaw;  // 'I': the int64; 'D': the double's bits; 'A'/'B': the bytes
};

struct FilRecord {
    std::int64_t mKey = 0;
    std::vector<FilWord> mWords;  // the attributes (after length and key)
};

struct FilData {
    std::vector<FilRecord> mRecords;
    bool mSwap = false;  // binary words in the other byte order
};

std::uint64_t fil_bytes_to_raw(const char* p) {
    std::uint64_t v;
    std::memcpy(&v, p, 8);
    return v;
}

std::int64_t fil_word_int(const FilWord& rWord, bool Swap) {
    switch (rWord.mTag) {
        case 'I':
            return static_cast<std::int64_t>(rWord.mRaw);
        case 'D':
            return static_cast<std::int64_t>(std::bit_cast<double>(rWord.mRaw));
        case 'A': {
            char s[9] = {};
            std::memcpy(s, &rWord.mRaw, 8);
            return std::atoll(s);
        }
        default: {
            // An 8-byte integer, or a 4-byte one in the word's first bytes
            // (Fortran EQUIVALENCE of an INTEGER array onto the REAL*8 one).
            const std::uint64_t v = Swap ? detail::bswap64(rWord.mRaw) : rWord.mRaw;
            const auto wide = static_cast<std::int64_t>(v);
            if (wide >= -2147483648LL && wide <= 2147483647LL)
                return wide;
            std::uint32_t lo;
            std::memcpy(&lo, &rWord.mRaw, 4);
            if (Swap)
                lo = detail::bswap32(lo);
            return static_cast<std::int32_t>(lo);
        }
    }
}

double fil_word_real(const FilWord& rWord, bool Swap) {
    switch (rWord.mTag) {
        case 'I':
            return static_cast<double>(static_cast<std::int64_t>(rWord.mRaw));
        case 'D':
            return std::bit_cast<double>(rWord.mRaw);
        case 'A':
            return std::numeric_limits<double>::quiet_NaN();
        default:
            return std::bit_cast<double>(Swap ? detail::bswap64(rWord.mRaw) : rWord.mRaw);
    }
}

std::string fil_word_text(const FilWord& rWord, bool Swap) {
    if (rWord.mTag == 'I' || rWord.mTag == 'D')
        return std::to_string(fil_word_int(rWord, Swap));
    char s[8];
    std::memcpy(s, &rWord.mRaw, 8);
    return std::string(s, 8);
}

std::string fil_trim(std::string s) {
    const std::size_t b = s.find_first_not_of(" \t\0", 0, 3);
    if (b == std::string::npos)
        return {};
    const std::size_t e = s.find_last_not_of(" \t\0", std::string::npos, 3);
    return s.substr(b, e - b + 1);
}

// Binary: the 512-word blocks' payloads, joined, cut into records.
void fil_parse_binary(std::string_view rText, FilData& rOut) {
    std::string words;
    if (const auto layout = detail::sniff_fortran_records(rText.data(), rText.size())) {
        for (const auto& r :
             detail::fortran_records(rText.data(), rText.size(), *layout, "Abaqus .fil"))
            words.append(rText, r.mOffset, r.mSize);
        rOut.mSwap = layout->mBigEndian != (std::endian::native == std::endian::big);
    } else {
        // Blocks written without record markers.
        words = rText;
    }
    if (words.size() % 8 != 0)
        throw ReadError("Abaqus .fil: binary payload is not a whole number of 8-byte words");
    const std::size_t n = words.size() / 8;
    std::size_t w = 0;
    while (w < n) {
        const FilWord len_word{'B', fil_bytes_to_raw(words.data() + 8 * w)};
        const std::int64_t len = fil_word_int(len_word, rOut.mSwap);
        if (len == 0) {
            // Zero padding after the last record of a block.
            ++w;
            continue;
        }
        if (len < 2 || static_cast<std::uint64_t>(len) > n - w)
            throw ReadError("Abaqus .fil: record at word " + std::to_string(w) +
                            " has invalid length " + std::to_string(len));
        FilRecord rec;
        rec.mKey =
            fil_word_int(FilWord{'B', fil_bytes_to_raw(words.data() + 8 * (w + 1))}, rOut.mSwap);
        for (std::int64_t k = 2; k < len; ++k)
            rec.mWords.push_back(
                {'B', fil_bytes_to_raw(words.data() + 8 * (w + static_cast<std::size_t>(k)))});
        rOut.mRecords.push_back(std::move(rec));
        w += static_cast<std::size_t>(len);
    }
}

// ASCII: line breaks dropped, then item by item from each `*`.
void fil_parse_ascii(std::string_view rText, FilData& rOut) {
    std::string s;
    s.reserve(rText.size());
    for (char c : rText)
        if (c != '\n' && c != '\r')
            s += c;
    std::size_t pos = 0;
    auto fail = [&](const std::string& rWhat) {
        throw ReadError("Abaqus .fil: " + rWhat + " (character " + std::to_string(pos) + ")");
    };
    auto item = [&]() -> FilWord {
        if (pos >= s.size())
            fail("the file ends inside a record");
        const char tag = s[pos];
        if (tag == 'I') {
            if (pos + 3 > s.size())
                fail("truncated integer item");
            const std::string width = s.substr(pos + 1, 2);
            const std::size_t first = width.find_first_not_of(' ');
            if (first == std::string::npos)
                fail("bad integer width");
            const long n = std::atol(width.c_str() + first);
            if (n < 1 || pos + 3 + static_cast<std::size_t>(n) > s.size())
                fail("bad integer width");
            const std::string digits = s.substr(pos + 3, static_cast<std::size_t>(n));
            for (std::size_t k = 0; k < digits.size(); ++k)
                if (!(std::isdigit(static_cast<unsigned char>(digits[k])) ||
                      (k == 0 && (digits[k] == '-' || digits[k] == ' '))))
                    fail("bad integer '" + digits + "'");
            pos += 3 + static_cast<std::size_t>(n);
            return {'I', static_cast<std::uint64_t>(std::atoll(digits.c_str()))};
        }
        if (tag == 'D' || tag == 'E') {
            if (pos + 23 > s.size())
                fail("truncated real item");
            std::string num = s.substr(pos + 1, 22);
            for (char& c : num)
                if (c == 'D' || c == 'd')
                    c = 'E';
            // D22.15 with a three-digit exponent drops the letter: 1.0+100.
            if (num.find('E') == std::string::npos) {
                const std::size_t sign = num.find_last_of("+-");
                if (sign != std::string::npos && sign > 2)
                    num.insert(sign, 1, 'E');
            }
            const char* end = nullptr;
            const double v = detail::parse_double(num.c_str(), end);
            if (end == num.c_str())
                fail("bad real '" + num + "'");
            pos += 23;
            return {'D', std::bit_cast<std::uint64_t>(v)};
        }
        if (tag == 'A') {
            if (pos + 9 > s.size())
                fail("truncated text item");
            char buf[8];
            std::memcpy(buf, s.data() + pos + 1, 8);
            pos += 9;
            return {'A', fil_bytes_to_raw(buf)};
        }
        fail(std::string("unknown item tag '") + tag + "'");
        return FilWord{'I', 0};
    };
    while (true) {
        pos = s.find('*', pos);
        if (pos == std::string::npos)
            break;
        ++pos;
        const FilWord len_word = item();
        if (len_word.mTag != 'I')
            fail("a record starts with its length");
        const auto len = static_cast<std::int64_t>(len_word.mRaw);
        if (len < 2)
            fail("record length " + std::to_string(len));
        FilRecord rec;
        rec.mKey = fil_word_int(item(), false);
        for (std::int64_t k = 2; k < len; ++k)
            rec.mWords.push_back(item());
        rOut.mRecords.push_back(std::move(rec));
    }
}

FilData fil_parse(const std::string& rPath) {
    const detail::FileSource text_source =
        detail::open_source(rPath, "Abaqus .fil: cannot open " + rPath);
    const std::string_view text = text_source.View();
    FilData out;
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first != std::string::npos && text[first] == '*')
        fil_parse_ascii(text, out);
    else
        fil_parse_binary(text, out);
    return out;
}

// --- names ------------------------------------------------------------------------

const char* fil_nodal_name(std::int64_t Key) {
    static const std::unordered_map<std::int64_t, const char*> m = {
        {101, "U"},     {102, "V"},    {103, "A"},     {104, "RF"},   {105, "EPOT"}, {106, "CF"},
        {107, "COORD"}, {108, "POR"},  {109, "RVF"},   {110, "RVT"},  {113, "TU"},   {114, "TV"},
        {115, "TA"},    {119, "RCHG"}, {120, "CECHG"}, {136, "PCAV"}, {137, "CVOL"}, {145, "VF"},
        {146, "TF"},    {151, "PABS"}, {201, "NT"},    {204, "RFL"},  {206, "CFL"},  {214, "RFLE"},
        {221, "NNC"},   {320, "CFF"},
    };
    const auto it = m.find(Key);
    return it == m.end() ? nullptr : it->second;
}

const char* fil_element_name(std::int64_t Key) {
    static const std::unordered_map<std::int64_t, const char*> m = {
        {2, "TEMP"},
        {3, "LOADS"},
        {4, "FLUXS"},
        {5, "SDV"},
        {6, "VOIDR"},
        {7, "FOUND"},
        {8, "COORD"},
        {9, "FV"},
        {10, "NFLUX"},
        {11, "S"},
        {12, "SINV"},
        {13, "SF"},
        {14, "ENER"},
        {15, "NFORC"},
        {17, "JK"},
        {18, "POR"},
        {19, "ELEN"},
        {21, "E"},
        {22, "PE"},
        {23, "CE"},
        {24, "IE"},
        {25, "EE"},
        {26, "CRACK"},
        {27, "STH"},
        {28, "HFL"},
        {29, "SE"},
        {30, "DG"},
        {31, "CONF"},
        {32, "SJP"},
        {35, "SAT"},
        {36, "SS"},
        {38, "CONC"},
        {39, "MFL"},
        {42, "SPE"},
        {45, "PEQC"},
        {47, "SEPE"},
        {48, "TSHR"},
        {50, "EPG"},
        {51, "EFLX"},
        {61, "STATUS"},
        {73, "PEEQ"},
        {74, "PRESS"},
        {75, "MISES"},
        {76, "IVOL"},
        {77, "SVOL"},
        {78, "EVOL"},
        {83, "SSAVG"},
        {86, "ALPHA"},
        {87, "UVARM"},
        {88, "THE"},
        {89, "LE"},
        {90, "NE"},
        {91, "ER"},
        {401, "SP"},
        {402, "ALPHAP"},
        {403, "EP"},
        {404, "NEP"},
        {405, "LEP"},
        {406, "ERP"},
        {407, "DGP"},
        {408, "EEP"},
        {409, "IEP"},
        {410, "THEP"},
        {411, "PEP"},
        {412, "CEP"},
        // Abaqus/Explicit only
        {421, "CKE"},
        {422, "CKLE"},
        {423, "CKLS"},
        {424, "CKSTAT"},
        {441, "CKEMAG"},
        {476, "EMSF"},
        {477, "EDT"},
        {507, "CFAILST"},
        {559, "CDMG"},
        {560, "CDIF"},
        {561, "CDIM"},
        {562, "CDIP"},
    };
    const auto it = m.find(Key);
    return it == m.end() ? nullptr : it->second;
}

// Keys read elsewhere (modal, energies, contact, element matrices) or not at
// all (contour integrals, cavity radiation, section output...): none is a
// nodal or element result.
bool fil_is_skipped_key(std::int64_t Key) {
    return (Key >= 301 && Key <= 310) || Key >= 1000;
}

// Records 301-310, per increment of a mode-based dynamic step.
const char* fil_modal_name(std::int64_t Key) {
    static const char* const names[] = {"GU",  "GV",  "GA",  "BM", "GPU",
                                        "GPV", "GPA", "SNE", "KE", "T"};
    return Key >= 301 && Key <= 310 ? names[Key - 301] : nullptr;
}

// Record 1999's attributes: Standard's and, where they differ, Explicit's
// (7 unused, 9 ALLDC, 15 DMASS, 17-18 heat energies).
const char* fil_energy_name(std::size_t Index, bool Explicit) {
    static const char* const standard[] = {"ALLKE", "ALLSE", "ALLWK", "ALLPD", "ALLCD", "ALLVD",
                                           "ALLKL", "ALLAE", "ALLQB", "ALLEE", "ALLIE", "ETOTAL",
                                           "ALLFD", "ALLJD", "ALLSD", "ALLDMD"};
    static const char* const explicit_[] = {"ALLKE", "ALLSE", "ALLWK", "ALLPD",  "ALLCD",  "ALLVD",
                                            nullptr, "ALLAE", "ALLDC", nullptr,  "ALLIE",  "ETOTAL",
                                            "ALLFD", nullptr, "DMASS", "ALLDMD", "ALLIHE", "ALLHF"};
    if (Explicit)
        return Index < 18 ? explicit_[Index] : nullptr;
    return Index < 16 ? standard[Index] : nullptr;
}

// Contact output (records 1511-1578, after a 1504 node header).
const char* fil_contact_name(std::int64_t Key) {
    static const std::unordered_map<std::int64_t, const char*> m = {
        {1511, "CSTRESS"}, {1512, "CDSTRESS"}, {1521, "CDISP"}, {1522, "CFN"},   {1523, "CFS"},
        {1524, "CAREA"},   {1526, "CMN"},      {1527, "CMS"},   {1528, "HFL"},   {1529, "HFLA"},
        {1530, "HTL"},     {1531, "HTLA"},     {1532, "SFDR"},  {1533, "SFDRA"}, {1534, "SFDRT"},
        {1535, "SFDRTA"},  {1536, "WEIGHT"},   {1537, "SJD"},   {1538, "SJDA"},  {1539, "SJDT"},
        {1540, "SJDTA"},   {1541, "ECD"},      {1542, "ECDA"},  {1543, "ECDT"},  {1544, "ECDTA"},
        {1545, "PFL"},     {1546, "PFLA"},     {1547, "PTL"},   {1548, "PTLA"},  {1549, "TPFL"},
        {1550, "TPTL"},    {1570, "DBT"},      {1571, "DBSF"},  {1572, "DBS"},   {1573, "XN"},
        {1574, "XS"},      {1575, "CFT"},      {1576, "CMT"},   {1577, "XT"},    {1578, "CTRQ"},
        {1592, "PPRESS"},
    };
    const auto it = m.find(Key);
    return it == m.end() ? nullptr : it->second;
}

// Element matrix records (1011-1031): the field data stem of each.
const char* fil_matrix_name(std::int64_t Key) {
    switch (Key) {
        case 1011:
        case 1012:
            return "stiffness";
        case 1021:
        case 1022:
            return "mass";
        case 1031:
            return "load";
        default:
            return nullptr;
    }
}

// --- the model --------------------------------------------------------------------

struct FilIncrement {
    std::size_t mBegin = 0, mEnd = 0;  // record range after the 2000 record
    double mTotalTime = 0.0, mStepTime = 0.0;
    std::int64_t mProcedure = 0, mStep = 0, mIncrement = 0;
};

// An eigenvalue step writes one increment whose modes each start with a 1980
// record: every mode becomes an increment of its own, from its 1980 record to
// the next one.
std::vector<FilIncrement> fil_split_modes(const FilData& rData,
                                          const std::vector<FilIncrement>& rIncrements) {
    std::vector<FilIncrement> out;
    for (const FilIncrement& inc : rIncrements) {
        std::vector<std::size_t> modes;
        for (std::size_t r = inc.mBegin; r < inc.mEnd; ++r)
            if (rData.mRecords[r].mKey == 1980)
                modes.push_back(r);
        if (modes.empty()) {
            out.push_back(inc);
            continue;
        }
        for (std::size_t k = 0; k < modes.size(); ++k) {
            FilIncrement m = inc;
            m.mBegin = modes[k];
            m.mEnd = k + 1 < modes.size() ? modes[k + 1] : inc.mEnd;
            out.push_back(m);
        }
    }
    return out;
}

std::vector<FilIncrement> fil_increments(const FilData& rData) {
    std::vector<FilIncrement> out;
    for (std::size_t r = 0; r < rData.mRecords.size(); ++r) {
        const FilRecord& rec = rData.mRecords[r];
        if (rec.mKey == 2000) {
            FilIncrement inc;
            inc.mBegin = r + 1;
            const auto& w = rec.mWords;
            const bool s = rData.mSwap;
            if (w.size() > 0)
                inc.mTotalTime = fil_word_real(w[0], s);
            if (w.size() > 1)
                inc.mStepTime = fil_word_real(w[1], s);
            if (w.size() > 4)
                inc.mProcedure = fil_word_int(w[4], s);
            if (w.size() > 5)
                inc.mStep = fil_word_int(w[5], s);
            if (w.size() > 6)
                inc.mIncrement = fil_word_int(w[6], s);
            if (!out.empty() && out.back().mEnd == 0)
                out.back().mEnd = r;
            out.push_back(inc);
        } else if (rec.mKey == 2001 && !out.empty() && out.back().mEnd == 0) {
            out.back().mEnd = r;
        }
    }
    if (!out.empty() && out.back().mEnd == 0)
        out.back().mEnd = rData.mRecords.size();
    return fil_split_modes(rData, out);
}

NDArray fil_scalar(DType Type, double Value) {
    NDArray a(Type, {});
    if (Type == DType::Int64)
        a.As<std::int64_t>()[0] = static_cast<std::int64_t>(Value);
    else
        a.As<double>()[0] = Value;
    return a;
}

NDArray fil_ids(const std::vector<std::int64_t>& rIds) {
    NDArray a(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

struct FilSet {
    std::string mName;
    bool mNodes;
    std::vector<std::int64_t> mLabels;
};

bool fil_is_explicit(std::int64_t Procedure) {
    return Procedure == 17 || Procedure == 21 || Procedure == 74;
}

// A contact surface (records 1501, 1502): its facets as (element, face key).
struct FilSurface {
    std::string mName;
    std::int64_t mType = 1;  // 1 deformable, 2 rigid
    std::vector<std::pair<std::int64_t, std::int64_t>> mFacets;
};

// Element matrices of one kind (records 1011-1031): the values of every
// element one after the other, and per element (label, offset, count, flag).
struct FilMatrices {
    std::vector<double> mValues;
    std::vector<std::int64_t> mIndex;
};

// One element result: the values of every (element, point) the step lists.
struct FilElementField {
    std::string mName;
    int mLocation;
    std::unordered_map<std::int64_t, std::map<std::int64_t, std::vector<double>>> mValues;
};

}  // namespace

Mesh read_abaqus_fil(const std::string& rPath, const ReadOptions& rOpts) {
    const FilData data = fil_parse(rPath);
    const bool sw = data.mSwap;

    // --- model records --------------------------------------------------------------
    std::vector<std::int64_t> node_labels;
    std::vector<double> coords;
    struct FilElement {
        std::int64_t mLabel;
        std::string mType;
        std::vector<std::int64_t> mNodes;
    };
    std::vector<FilElement> elements;
    std::vector<FilSet> sets;
    std::unordered_map<std::int64_t, std::string> labels;  // 1940
    std::vector<FilSurface> surfaces;
    std::map<std::string, FilMatrices> matrices;  // stiffness, mass, load, dofs
    std::int64_t matrix_element = 0, last_matrix_key = 0;
    for (const FilRecord& rec : data.mRecords) {
        const auto& w = rec.mWords;
        if (rec.mKey != last_matrix_key && (rec.mKey < 1001 || rec.mKey > 1043))
            last_matrix_key = 0;
        switch (rec.mKey) {
            case 1501: {  // contact surface: name, dimension, type, facet count, node
                if (w.empty())
                    break;
                FilSurface surf{fil_trim(fil_word_text(w[0], sw)), 1, {}};
                if (w.size() > 2)
                    surf.mType = fil_word_int(w[2], sw);
                surfaces.push_back(std::move(surf));
                break;
            }
            case 1502:  // a facet: element, face key, node count, nodes
                if (!surfaces.empty() && w.size() >= 2)
                    surfaces.back().mFacets.emplace_back(fil_word_int(w[0], sw),
                                                         fil_word_int(w[1], sw));
                break;
            case 1001:  // element matrix header: element, type, nodes
                matrix_element = w.empty() ? 0 : fil_word_int(w[0], sw);
                last_matrix_key = 0;
                break;
            case 1002:
            case 1011:
            case 1012:
            case 1021:
            case 1022:
            case 1031: {
                const std::string kind =
                    rec.mKey == 1002 ? "matrix_dofs" : fil_matrix_name(rec.mKey);
                FilMatrices& m = matrices[kind];
                std::size_t first = 0;
                std::int64_t flag = rec.mKey == 1011 || rec.mKey == 1021 ? 1 : 0;
                if (rec.mKey == 1031) {  // a load case, then the loads
                    flag = w.empty() ? 0 : fil_word_int(w[0], sw);
                    first = 1;
                }
                const bool continued = rec.mKey == last_matrix_key && !m.mIndex.empty();
                if (!continued) {
                    m.mIndex.insert(
                        m.mIndex.end(),
                        {matrix_element, static_cast<std::int64_t>(m.mValues.size()), 0, flag});
                }
                for (std::size_t k = continued ? 0 : first; k < w.size(); ++k)
                    m.mValues.push_back(rec.mKey == 1002
                                            ? static_cast<double>(fil_word_int(w[k], sw))
                                            : fil_word_real(w[k], sw));
                m.mIndex[m.mIndex.size() - 2] =
                    static_cast<std::int64_t>(m.mValues.size()) - m.mIndex[m.mIndex.size() - 3];
                last_matrix_key = rec.mKey;
                break;
            }
            case 1901: {
                if (w.empty())
                    break;
                node_labels.push_back(fil_word_int(w[0], sw));
                for (std::size_t d = 0; d < 3; ++d)
                    coords.push_back(d + 1 < w.size() ? fil_word_real(w[d + 1], sw) : 0.0);
                break;
            }
            case 1900: {
                if (w.size() < 2)
                    break;
                FilElement el{fil_word_int(w[0], sw), fil_trim(fil_word_text(w[1], sw)), {}};
                for (char& c : el.mType)
                    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                for (std::size_t k = 2; k < w.size(); ++k)
                    el.mNodes.push_back(fil_word_int(w[k], sw));
                elements.push_back(std::move(el));
                break;
            }
            case 1990:
                if (!elements.empty())
                    for (const FilWord& x : w)
                        elements.back().mNodes.push_back(fil_word_int(x, sw));
                break;
            case 1931:
            case 1933: {
                if (w.empty())
                    break;
                FilSet set{fil_trim(fil_word_text(w[0], sw)), rec.mKey == 1931, {}};
                for (std::size_t k = 1; k < w.size(); ++k)
                    set.mLabels.push_back(fil_word_int(w[k], sw));
                sets.push_back(std::move(set));
                break;
            }
            case 1932:
            case 1934:
                if (!sets.empty())
                    for (const FilWord& x : w)
                        sets.back().mLabels.push_back(fil_word_int(x, sw));
                break;
            case 1940: {
                if (w.empty())
                    break;
                std::string label;
                for (std::size_t k = 1; k < w.size(); ++k)
                    label += fil_word_text(w[k], sw);
                labels[fil_word_int(w[0], sw)] = fil_trim(label);
                break;
            }
            default:
                break;
        }
    }

    Mesh mesh;
    std::unordered_map<std::int64_t, std::size_t> node_index;
    for (std::size_t p = 0; p < node_labels.size(); ++p)
        if (!node_index.emplace(node_labels[p], p).second)
            throw ReadError("Abaqus .fil: node " + std::to_string(node_labels[p]) +
                            " is defined twice");
    NDArray points(DType::Float64, {node_labels.size(), 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    mesh.AddPointData("abaqus:id", fil_ids(node_labels));

    // Cells: one block per meshio++ type, in order of first appearance.
    std::vector<std::string> block_types;
    std::map<std::string, std::vector<std::size_t>> by_type;
    std::map<std::string, std::size_t> skipped_types;
    std::vector<std::string> cell_type_of(elements.size());
    for (std::size_t e = 0; e < elements.size(); ++e) {
        const std::string t =
            detail::abaqus_cell_type(elements[e].mType, elements[e].mNodes.size());
        if (t.empty()) {
            ++skipped_types[elements[e].mType];
            continue;
        }
        auto [it, fresh] = by_type.emplace(t, std::vector<std::size_t>{});
        if (fresh)
            block_types.push_back(t);
        it->second.push_back(e);
        cell_type_of[e] = t;
    }
    for (const auto& [t, n] : skipped_types)
        log::warn("Abaqus .fil: skipping {} element(s) of type {} (no meshio++ equivalent)", n, t);
    std::unordered_map<std::int64_t, std::size_t> element_index;  // label -> global cell
    std::vector<std::size_t> block_start{0};
    std::vector<std::vector<std::int64_t>> cell_nodes;  // global cell -> node labels
    std::vector<int> cell_dim;
    std::vector<NDArray> id_blocks;
    for (const std::string& t : block_types) {
        const std::vector<std::size_t>& members = by_type[t];
        const std::size_t k = static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(t)));
        NDArray conn(DType::Int64, {members.size(), k});
        std::vector<std::int64_t> ids;
        for (std::size_t r = 0; r < members.size(); ++r) {
            const FilElement& el = elements[members[r]];
            for (std::size_t j = 0; j < k; ++j) {
                const auto it = node_index.find(el.mNodes[j]);
                if (it == node_index.end())
                    throw ReadError("Abaqus .fil: element " + std::to_string(el.mLabel) +
                                    " names undefined node " + std::to_string(el.mNodes[j]));
                conn.As<std::int64_t>()[r * k + j] = static_cast<std::int64_t>(it->second);
            }
            if (!element_index.emplace(el.mLabel, cell_nodes.size()).second)
                throw ReadError("Abaqus .fil: element " + std::to_string(el.mLabel) +
                                " is defined twice");
            cell_nodes.push_back(el.mNodes);
            cell_dim.push_back(cell_type_dimension(cell_type_from_name(t)));
            ids.push_back(el.mLabel);
        }
        mesh.AddCellBlock(t, std::move(conn));
        id_blocks.push_back(fil_ids(ids));
        block_start.push_back(cell_nodes.size());
    }
    if (!id_blocks.empty())
        mesh.AddCellData("abaqus:id", std::move(id_blocks));

    // Sets -> regions; a digits-only name is a 1940 cross-reference.
    std::size_t missing = 0;
    for (const FilSet& set : sets) {
        std::string name = set.mName;
        if (!name.empty() && name.find_first_not_of("0123456789") == std::string::npos) {
            const auto it = labels.find(std::stoll(name));
            if (it != labels.end())
                name = it->second;
        }
        std::vector<std::int64_t> entries;
        int dim = -1;
        for (std::int64_t label : set.mLabels) {
            if (set.mNodes) {
                const auto it = node_index.find(label);
                if (it == node_index.end())
                    ++missing;
                else
                    entries.push_back(static_cast<std::int64_t>(it->second));
            } else {
                const auto it = element_index.find(label);
                if (it == element_index.end()) {
                    ++missing;
                } else {
                    entries.push_back(static_cast<std::int64_t>(it->second));
                    dim = std::max(dim, cell_dim[it->second]);
                }
            }
        }
        mesh.AddRegion(Region(name, set.mNodes ? RegionKind::Point : RegionKind::Cell,
                              set.mNodes ? -1 : dim, -1, fil_ids(entries)));
    }
    if (missing)
        log::warn("Abaqus .fil: sets name {} node(s) or element(s) that are not in the mesh",
                  missing);

    // Contact surfaces -> side regions (cell, local facet); rigid surfaces and
    // facets of unknown elements are left out.
    std::size_t unplaced = 0;
    for (const FilSurface& surf : surfaces) {
        std::string name = surf.mName;
        if (!name.empty() && name.find_first_not_of("0123456789") == std::string::npos) {
            const auto it = labels.find(std::stoll(name));
            if (it != labels.end())
                name = it->second;
        }
        std::vector<std::int64_t> entries;
        int dim = -1;
        for (const auto& [label, face] : surf.mFacets) {
            const auto it = element_index.find(label);
            if (it == element_index.end() || face < 1 || face > 8) {
                ++unplaced;
                continue;
            }
            const std::size_t c = it->second;
            const std::size_t b = static_cast<std::size_t>(
                std::upper_bound(block_start.begin(), block_start.end(), c) - block_start.begin() -
                1);
            const std::string key = face == 7   ? "SPOS"
                                    : face == 8 ? "SNEG"
                                                : "S" + std::to_string(face);
            const int facet = detail::abaqus_face_index(block_types[b], key);
            if (facet < 0) {
                ++unplaced;
                continue;
            }
            entries.push_back(static_cast<std::int64_t>(c));
            entries.push_back(facet);
            dim = std::max(dim, cell_dim[c] - 1);
        }
        NDArray a(DType::Int64, {entries.size() / 2, 2});
        std::copy(entries.begin(), entries.end(), a.As<std::int64_t>());
        mesh.AddRegion(Region(name, RegionKind::Side, dim, -1, std::move(a)));
    }
    if (unplaced)
        log::warn("Abaqus .fil: {} contact surface facet(s) name no element face of the mesh",
                  unplaced);
    for (const auto& [kind, m] : matrices) {
        const std::size_t rows = m.mIndex.size() / 4;
        const bool dofs = kind == "matrix_dofs";
        NDArray values(dofs ? DType::Int64 : DType::Float64, {m.mValues.size()});
        for (std::size_t k = 0; k < m.mValues.size(); ++k) {
            if (dofs)
                values.As<std::int64_t>()[k] = static_cast<std::int64_t>(m.mValues[k]);
            else
                values.As<double>()[k] = m.mValues[k];
        }
        NDArray index(DType::Int64, {rows, 4});
        std::copy(m.mIndex.begin(), m.mIndex.end(), index.As<std::int64_t>());
        mesh.AddFieldData("abaqus:" + kind, std::move(values));
        mesh.AddFieldData("abaqus:" + kind + ":index", std::move(index));
    }

    // --- the selected increment ----------------------------------------------------
    const std::vector<FilIncrement> increments = fil_increments(data);
    if (increments.empty()) {
        rOpts.ResolveTimeStep(0);
        return mesh;
    }
    const FilIncrement& inc = increments[rOpts.ResolveTimeStep(increments.size())];
    mesh.AddFieldData(kSequenceTimeKey, fil_scalar(DType::Float64, inc.mTotalTime));
    mesh.AddFieldData("abaqus:step", fil_scalar(DType::Int64, static_cast<double>(inc.mStep)));
    mesh.AddFieldData("abaqus:increment",
                      fil_scalar(DType::Int64, static_cast<double>(inc.mIncrement)));
    mesh.AddFieldData("abaqus:step_time", fil_scalar(DType::Float64, inc.mStepTime));
    mesh.AddFieldData("abaqus:procedure",
                      fil_scalar(DType::Int64, static_cast<double>(inc.mProcedure)));
    if (!rOpts.WantsAnyData())
        return mesh;

    // Nodal fields: name -> node label -> values.
    std::vector<std::string> nodal_order;
    std::map<std::string, std::unordered_map<std::int64_t, std::vector<double>>> nodal;
    // Element fields keyed by (name, location, section point).
    std::vector<FilElementField> fields;
    std::map<std::pair<std::string, int>, std::size_t> field_of;
    std::set<std::int64_t> skipped_keys;
    bool nodal_mode = false;
    struct Header {
        std::int64_t mElement, mPoint, mSection;
        int mLocation;
    };
    std::optional<Header> header;
    std::string rebar;  // the rebar the element header names
    const bool explicit_run = fil_is_explicit(inc.mProcedure);
    std::map<std::string, std::vector<std::vector<double>>> modal_rows;  // 301-310
    std::vector<std::string> modal_order;
    std::int64_t contact_node = 0;
    bool contact = false;
    for (std::size_t r = inc.mBegin; r < inc.mEnd; ++r) {
        const FilRecord& rec = data.mRecords[r];
        const auto& w = rec.mWords;
        if (rec.mKey == 1980) {  // a mode: number, eigenvalue, mass, damping, factors
            if (w.size() > 0)
                mesh.AddFieldData(
                    "abaqus:mode",
                    fil_scalar(DType::Int64, static_cast<double>(fil_word_int(w[0], sw))));
            if (w.size() > 1)
                mesh.AddFieldData("abaqus:eigenvalue",
                                  fil_scalar(DType::Float64, fil_word_real(w[1], sw)));
            if (w.size() > 2)
                mesh.AddFieldData("abaqus:generalized_mass",
                                  fil_scalar(DType::Float64, fil_word_real(w[2], sw)));
            if (w.size() > 3)
                mesh.AddFieldData("abaqus:composite_damping",
                                  fil_scalar(DType::Float64, fil_word_real(w[3], sw)));
            std::vector<double> factor, mass;
            for (std::size_t k = 4; k + 1 < w.size(); k += 2) {
                factor.push_back(fil_word_real(w[k], sw));
                mass.push_back(fil_word_real(w[k + 1], sw));
            }
            if (!factor.empty()) {
                NDArray f(DType::Float64, {factor.size()}), m(DType::Float64, {mass.size()});
                std::copy(factor.begin(), factor.end(), f.As<double>());
                std::copy(mass.begin(), mass.end(), m.As<double>());
                mesh.AddFieldData("abaqus:participation_factor", std::move(f));
                mesh.AddFieldData("abaqus:effective_mass", std::move(m));
            }
            continue;
        }
        if (rec.mKey == 1999) {  // total energies
            for (std::size_t k = 0; k < w.size(); ++k)
                if (const char* e = fil_energy_name(k, explicit_run))
                    mesh.AddFieldData(std::string("abaqus:") + e,
                                      fil_scalar(DType::Float64, fil_word_real(w[k], sw)));
            continue;
        }
        if (const char* g = fil_modal_name(rec.mKey)) {  // generalized quantities per mode
            std::vector<double> row;
            for (std::size_t k = 0; k < w.size(); ++k)
                if (!(rec.mKey == 304 && k == 7))  // BM's base name
                    row.push_back(rec.mKey == 304 && k == 0
                                      ? static_cast<double>(fil_word_int(w[k], sw))
                                      : fil_word_real(w[k], sw));
            auto [it, fresh] = modal_rows.emplace(g, std::vector<std::vector<double>>{});
            if (fresh)
                modal_order.push_back(g);
            it->second.push_back(std::move(row));
            continue;
        }
        if (rec.mKey == 1503) {  // contact output request: the node records follow
            contact = true;
            continue;
        }
        if (rec.mKey == 1504) {
            contact_node = w.empty() ? 0 : fil_word_int(w[0], sw);
            continue;
        }
        if (contact && rec.mKey >= 1505 && rec.mKey <= 1599) {
            const char* known = fil_contact_name(rec.mKey);
            const std::string name = known ? known : "key_" + std::to_string(rec.mKey);
            if (!rOpts.WantsArray(name))
                continue;
            std::vector<double> v;
            for (const FilWord& x : w)
                v.push_back(fil_word_real(x, sw));
            auto [it, fresh] = nodal.emplace(name, decltype(nodal)::mapped_type{});
            if (fresh)
                nodal_order.push_back(name);
            it->second[contact_node] = std::move(v);
            continue;
        }
        if (rec.mKey == 1911) {
            nodal_mode = !w.empty() && fil_word_int(w[0], sw) == 1;
            header.reset();
            contact = false;
            continue;
        }
        if (rec.mKey == 1) {
            if (w.size() < 4) {
                header.reset();
                continue;
            }
            header = Header{fil_word_int(w[0], sw), fil_word_int(w[1], sw), fil_word_int(w[2], sw),
                            static_cast<int>(fil_word_int(w[3], sw))};
            rebar = header->mLocation == 3 && w.size() > 4 ? fil_trim(fil_word_text(w[4], sw)) : "";
            if (!rebar.empty() && rebar.find_first_not_of("0123456789") == std::string::npos) {
                const auto it = labels.find(std::stoll(rebar));
                if (it != labels.end())
                    rebar = it->second;
            }
            continue;
        }
        if ((rec.mKey >= 1900 && rec.mKey <= 2001) || fil_is_skipped_key(rec.mKey)) {
            // model records, and those read from the whole file (surfaces,
            // element matrices) are not skipped results
            const bool whole_file =
                (rec.mKey >= 1001 && rec.mKey <= 1043) || rec.mKey == 1501 || rec.mKey == 1502;
            if ((rec.mKey < 1900 || rec.mKey > 2001) && !whole_file)
                skipped_keys.insert(rec.mKey);
            continue;
        }
        if (nodal_mode) {
            if (w.empty())
                continue;
            const char* known = fil_nodal_name(rec.mKey);
            const std::string name = known ? known : "key_" + std::to_string(rec.mKey);
            if (!rOpts.WantsArray(name))
                continue;
            std::vector<double> v;
            for (std::size_t k = 1; k < w.size(); ++k)
                v.push_back(fil_word_real(w[k], sw));
            auto [it, fresh] = nodal.emplace(name, decltype(nodal)::mapped_type{});
            if (fresh)
                nodal_order.push_back(name);
            it->second[fil_word_int(w[0], sw)] = std::move(v);
            continue;
        }
        if (!header) {
            skipped_keys.insert(rec.mKey);
            continue;
        }
        const char* known =
            rec.mKey == 79 ? (explicit_run ? "ERV" : "RATIO") : fil_element_name(rec.mKey);
        std::string name = known ? known : "key_" + std::to_string(rec.mKey);
        // Rebar values are per integration point of the element, named by the rebar.
        if (header->mLocation == 3)
            name += "@rebar:" + rebar;
        // Continuum elements write section point 0, shells and beams 1..n: the
        // first one shares the plain name, the others are `@sp<k>`.
        if (header->mSection > 1)
            name += "@sp" + std::to_string(header->mSection);
        if (!rOpts.WantsArray(name))
            continue;
        const int location =
            header->mLocation == 3 ? 0 : (header->mLocation == 5 ? 1 : header->mLocation);
        const auto key = std::make_pair(name, location);
        auto [fit, ffresh] = field_of.emplace(key, fields.size());
        if (ffresh)
            fields.push_back({name, location, {}});
        std::vector<double> v;
        for (const FilWord& x : w)
            v.push_back(fil_word_real(x, sw));
        fields[fit->second].mValues[header->mElement][header->mPoint] = std::move(v);
    }
    if (!skipped_keys.empty())
        log::warn(
            "Abaqus .fil: {} record key(s) outside the nodal and element results skipped "
            "(first: {})",
            skipped_keys.size(), *skipped_keys.begin());
    for (const std::string& g : modal_order) {  // one row per record, NaN-padded
        const auto& rows = modal_rows[g];
        std::size_t width = 0;
        for (const auto& row : rows)
            width = std::max(width, row.size());
        NDArray a(DType::Float64, rows.size() == 1 ? std::vector<std::size_t>{width}
                                                   : std::vector<std::size_t>{rows.size(), width});
        double* d = a.As<double>();
        std::fill(d, d + rows.size() * width, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t k = 0; k < rows.size(); ++k)
            std::copy(rows[k].begin(), rows[k].end(), d + k * width);
        mesh.AddFieldData("abaqus:" + g, std::move(a));
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t npts = node_labels.size();
    auto add_point_field =
        [&](const std::string& rName,
            const std::unordered_map<std::int64_t, std::vector<double>>& rValues) {
            std::size_t width = 0;
            for (const auto& [label, v] : rValues)
                width = std::max(width, v.size());
            if (!width)
                return;
            NDArray a(DType::Float64, width == 1 ? std::vector<std::size_t>{npts}
                                                 : std::vector<std::size_t>{npts, width});
            double* d = a.As<double>();
            std::fill(d, d + npts * width, nan);
            for (const auto& [label, v] : rValues) {
                const auto it = node_index.find(label);
                if (it != node_index.end())
                    std::copy(v.begin(), v.end(), d + it->second * width);
            }
            std::string name = rName;
            while (mesh.HasPointData(name))
                name += "@avg";
            mesh.AddPointData(name, std::move(a));
        };
    for (const std::string& name : nodal_order)
        add_point_field(name, nodal[name]);

    const std::size_t nblocks = block_types.size();
    for (const FilElementField& f : fields) {
        if (f.mLocation == 4) {  // averaged at the nodes: the header names the node
            std::unordered_map<std::int64_t, std::vector<double>> per_node;
            for (const auto& [label, points] : f.mValues)
                if (!points.empty())
                    per_node[label] = points.begin()->second;
            add_point_field(f.mName, per_node);
            continue;
        }
        if (!nblocks)
            continue;
        // One rectangular array with the same width in every block (what every
        // writer can hold): per-point data is flattened point-major, column
        // `point * W + component`, NaN where a block has fewer points or
        // components.
        const bool per_point = f.mLocation == 0 || f.mLocation == 2;
        std::size_t width = 0, count = 1;
        if (f.mLocation == 2)
            for (std::size_t b = 0; b < nblocks; ++b)
                count = std::max(count, mesh.Cells(b).NodesPerCell());
        std::vector<std::vector<const std::map<std::int64_t, std::vector<double>>*>> rows(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b)
            rows[b].assign(block_start[b + 1] - block_start[b], nullptr);
        for (const auto& [label, points] : f.mValues) {
            const auto it = element_index.find(label);
            if (it == element_index.end())
                continue;
            const std::size_t c = it->second;
            const std::size_t b = static_cast<std::size_t>(
                std::upper_bound(block_start.begin(), block_start.end(), c) - block_start.begin() -
                1);
            rows[b][c - block_start[b]] = &points;
            for (const auto& [pt, v] : points) {
                width = std::max(width, v.size());
                if (f.mLocation == 0)
                    count =
                        std::max(count, static_cast<std::size_t>(std::max<std::int64_t>(pt, 1)));
            }
        }
        if (!width)
            continue;
        // Point numbers and value counts come from the records: a corrupt one
        // must not size an array of billions of columns.
        if (count > (std::size_t{1} << 16) || width > (std::size_t{1} << 16))
            throw ReadError("Abaqus .fil: implausible integration point or component count");
        const std::size_t pts = per_point ? count : 1;
        const std::size_t cols = pts * width;
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < nblocks; ++b) {
            const std::size_t n = rows[b].size();
            NDArray a(DType::Float64,
                      cols == 1 ? std::vector<std::size_t>{n} : std::vector<std::size_t>{n, cols});
            double* d = a.As<double>();
            std::fill(d, d + n * cols, nan);
            for (std::size_t r = 0; r < n; ++r) {
                if (!rows[b][r])
                    continue;
                const std::vector<std::int64_t>& cn = cell_nodes[block_start[b] + r];
                for (const auto& [pt, v] : *rows[b][r]) {
                    std::size_t slot = 0;
                    if (f.mLocation == 0) {
                        slot = static_cast<std::size_t>(std::max<std::int64_t>(pt, 1) - 1);
                    } else if (f.mLocation == 2) {
                        const auto at = std::find(cn.begin(), cn.end(), pt);
                        if (at == cn.end())
                            continue;
                        slot = static_cast<std::size_t>(at - cn.begin());
                    }
                    if (slot >= pts)
                        continue;
                    std::copy(v.begin(),
                              v.begin() + static_cast<std::ptrdiff_t>(std::min(v.size(), width)),
                              d + r * cols + slot * width);
                }
            }
            blocks.push_back(std::move(a));
        }
        std::string name = f.mName;
        while (mesh.HasCellData(name))
            name += f.mLocation == 0 ? "@ip" : "@el";
        mesh.AddCellData(name, std::move(blocks));
        if (per_point) {
            // How to unflatten it: (points, components).
            NDArray layout(DType::Int64, {2});
            layout.As<std::int64_t>()[0] = static_cast<std::int64_t>(pts);
            layout.As<std::int64_t>()[1] = static_cast<std::int64_t>(width);
            mesh.AddFieldData("abaqus:layout:" + name, std::move(layout));
        }
    }
    return mesh;
}

std::vector<double> abaqus_fil_time_values(const std::string& rPath) {
    const FilData data = fil_parse(rPath);
    std::vector<double> out;
    for (const FilIncrement& inc : fil_increments(data))
        out.push_back(inc.mTotalTime);
    return out;
}

MeshMetadata read_abaqus_fil_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    MeshMetadata meta = metadata_from_mesh(read_abaqus_fil(rPath, ReadOptions{}));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "abaqus_fil";
    meta.mTimeValues = abaqus_fil_time_values(rPath);
    return meta;
}

}  // namespace meshioplusplus
