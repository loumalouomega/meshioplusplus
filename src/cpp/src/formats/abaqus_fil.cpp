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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
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
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"
#include "abaqus_face.hpp"
#include "../detail/open_source.hpp"
#include "../detail/typed_view.hpp"

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
            v.reserve(w.size());
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
        v.reserve(w.size());
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

// ---------------------------------------------------------------------------------------
// The writer: the reader's inverse for the model, the nodal and element results, in the
// ASCII layout or the binary one (512-word blocks, little-endian). The twin of
// abaqus_fil/_write.py, which writes the same bytes.
// ---------------------------------------------------------------------------------------

namespace {

[[noreturn]] void filw_fail(const std::string& rMessage) {
    throw WriteError("Abaqus .fil: " + rMessage);
}

/// One item of a record: an integer word, a real word or an 8-character text.
struct FilwItem {
    char mTag = 'I';
    std::int64_t mInt = 0;
    double mReal = 0.0;
    std::string mText;
};
using FilwRecord = std::vector<FilwItem>;

FilwItem filw_i(std::int64_t V) {
    FilwItem it;
    it.mTag = 'I';
    it.mInt = V;
    return it;
}

FilwItem filw_d(double V) {
    if (V != V || std::isinf(V))
        filw_fail("a value is not finite");
    if (V != 0.0 && !(std::fabs(V) >= 1e-99 && std::fabs(V) < 1e100))
        filw_fail("a value is outside the range a D22.15 item holds");
    FilwItem it;
    it.mTag = 'D';
    it.mReal = V;
    return it;
}

/// The text as the file holds it: ASCII, a `?` for every other byte.
std::string filw_ascii(const std::string& rText) {
    std::string out = rText;
    for (char& c : out)
        if (static_cast<unsigned char>(c) >= 0x80)
            c = '?';
    return out;
}

FilwItem filw_a(const std::string& rText) {
    std::string t = filw_ascii(rText);
    t.resize(8, ' ');
    FilwItem it;
    it.mTag = 'A';
    it.mText = t;
    return it;
}

void filw_text_items(FilwRecord& rOut, const std::string& rText, std::size_t Words = 10) {
    std::string t = filw_ascii(rText);
    t.resize(8 * Words, ' ');
    for (std::size_t k = 0; k < Words; ++k)
        rOut.push_back(filw_a(t.substr(8 * k, 8)));
}

FilwRecord filw_rec(std::int64_t Key) {
    return FilwRecord{filw_i(Key)};
}

std::string filw_ascii_bytes(const std::vector<FilwRecord>& rRecords) {
    std::string text, line;
    const auto put = [&](const std::string& rPiece) {
        line += rPiece;
        while (line.size() >= 80) {
            text += line.substr(0, 80);
            text += '\n';
            line.erase(0, 80);
        }
    };
    for (const FilwRecord& r : rRecords) {
        put("*");
        std::vector<FilwItem> items;
        items.push_back(filw_i(static_cast<std::int64_t>(r.size()) + 1));
        items.insert(items.end(), r.begin(), r.end());
        for (const FilwItem& it : items) {
            if (it.mTag == 'I') {
                const std::string digits = std::to_string(it.mInt);
                std::string count = std::to_string(digits.size());
                if (count.size() < 2)
                    count.insert(0, 2 - count.size(), ' ');
                put("I" + count + digits);
            } else if (it.mTag == 'D') {
                char buf[64];
                detail::snprintf_c(buf, sizeof(buf), "% .15E", it.mReal);
                std::string v = buf;
                std::replace(v.begin(), v.end(), 'E', 'D');
                put("D" + v);
            } else {
                put("A" + it.mText);
            }
        }
        if (r[0].mInt == 2001) {
            line.resize(80, ' ');
            text += line + "\n" + std::string(80, ' ') + "\n";
            line.clear();
        }
    }
    if (!line.empty())
        text += line + "\n";
    return text;
}

void filw_put_le(std::string& rOut, std::uint64_t V, std::size_t N) {
    for (std::size_t i = 0; i < N; ++i)
        rOut.push_back(static_cast<char>((V >> (8 * i)) & 0xffu));
}

void filw_put_be(std::string& rOut, std::uint64_t V, std::size_t N) {
    for (std::size_t i = N; i-- > 0;)
        rOut.push_back(static_cast<char>((V >> (8 * i)) & 0xffu));
}

std::string filw_binary_bytes(const std::vector<FilwRecord>& rRecords, bool BigEndian = false,
                              int MarkerSize = 4) {
    if (MarkerSize != 4 && MarkerSize != 8)
        filw_fail("the record marker size must be 4 or 8");
    const auto put = [BigEndian](std::string& rOut, std::uint64_t V, std::size_t N) {
        if (BigEndian)
            filw_put_be(rOut, V, N);
        else
            filw_put_le(rOut, V, N);
    };
    constexpr std::size_t kBlock = 512;
    std::vector<std::string> words;
    for (const FilwRecord& r : rRecords) {
        std::vector<std::string> rec_words;
        std::string w;
        put(w, static_cast<std::uint64_t>(r.size() + 1), 8);
        rec_words.push_back(w);
        for (const FilwItem& it : r) {
            std::string x;
            if (it.mTag == 'I') {
                put(x, static_cast<std::uint64_t>(it.mInt), 8);
            } else if (it.mTag == 'D') {
                put(x, std::bit_cast<std::uint64_t>(it.mReal), 8);
            } else {
                x = it.mText;
            }
            rec_words.push_back(x);
        }
        if (r[0].mInt == 2001) {
            // The last record of an increment pads to the end of its block.
            const std::size_t pad = (kBlock - (words.size() + rec_words.size()) % kBlock) % kBlock;
            std::string len;
            put(len, rec_words.size() + pad, 8);
            rec_words[0] = len;
            for (std::size_t k = 0; k < pad; ++k)
                rec_words.push_back(std::string(8, '\0'));
        }
        words.insert(words.end(), rec_words.begin(), rec_words.end());
    }
    const std::size_t tail = (kBlock - words.size() % kBlock) % kBlock;
    for (std::size_t k = 0; k < tail; ++k)
        words.push_back(std::string(8, '\0'));
    std::string marker;
    put(marker, 8 * kBlock, static_cast<std::size_t>(MarkerSize));
    std::string out;
    for (std::size_t b = 0; b < words.size(); b += kBlock) {
        out += marker;
        for (std::size_t k = b; k < b + kBlock; ++k)
            out += words[k];
        out += marker;
    }
    return out;
}

const std::map<std::string, std::string>& filw_types() {
    static const std::map<std::string, std::string> types = {
        {"line", "T3D2"},     {"line3", "T3D3"},       {"triangle", "S3"},
        {"quad", "S4R"},      {"triangle6", "STRI65"}, {"quad8", "S8R"},
        {"quad9", "S9R5"},    {"tetra", "C3D4"},       {"tetra10", "C3D10"},
        {"pyramid", "C3D5"},  {"pyramid13", "C3D13"},  {"wedge", "C3D6"},
        {"wedge15", "C3D15"}, {"hexahedron", "C3D8R"}, {"hexahedron20", "C3D20R"},
    };
    return types;
}

bool filw_is_shell(const std::string& rType) {
    return rType == "triangle" || rType == "quad" || rType == "triangle6" || rType == "quad8" ||
           rType == "quad9";
}

/// The Abaqus key of a result name: the reader's names inverted (`key_<n>` is `<n>`).
std::optional<std::int64_t> filw_key_of(const std::string& rName, bool Nodal) {
    static const std::unordered_map<std::string, std::int64_t> nodal = [] {
        std::unordered_map<std::string, std::int64_t> m;
        for (std::int64_t k = 1; k < 1000; ++k)
            if (const char* n = fil_nodal_name(k))
                m.emplace(n, k);
        return m;
    }();
    static const std::unordered_map<std::string, std::int64_t> element = [] {
        std::unordered_map<std::string, std::int64_t> m;
        for (std::int64_t k = 1; k < 1000; ++k)
            if (const char* n = fil_element_name(k))
                m.emplace(n, k);
        m["RATIO"] = 79;
        m["ERV"] = 79;
        return m;
    }();
    const auto& table = Nodal ? nodal : element;
    const auto it = table.find(rName);
    if (it != table.end())
        return it->second;
    if (rName.size() > 4 && rName.compare(0, 4, "key_") == 0 && rName.size() < 12 &&
        std::all_of(rName.begin() + 4, rName.end(), [](char C) { return C >= '0' && C <= '9'; }))
        return std::stoll(rName.substr(4));
    return std::nullopt;
}

/// The contact output record (1505-1599) a point array is: the reader's contact names, and
/// `key_<n>` for a key it has no name for.
std::optional<std::int64_t> filw_contact_key(const std::string& rName) {
    static const std::unordered_map<std::string, std::int64_t> named = [] {
        std::unordered_map<std::string, std::int64_t> m;
        for (std::int64_t k = 1505; k <= 1599; ++k)
            if (const char* n = fil_contact_name(k))
                m.emplace(n, k);
        return m;
    }();
    const auto it = named.find(rName);
    if (it != named.end())
        return it->second;
    if (rName.size() > 4 && rName.compare(0, 4, "key_") == 0 && rName.size() < 12 &&
        std::all_of(rName.begin() + 4, rName.end(), [](char C) { return C >= '0' && C <= '9'; })) {
        const std::int64_t key = std::stoll(rName.substr(4));
        if (key >= 1505 && key <= 1599 && fil_contact_name(key) == nullptr)
            return key;
    }
    return std::nullopt;
}

/// The (slave, master) side regions contact output is written for, as region indices: the two
/// indices of `abaqus:contact_pair`, else the first two side regions; none without a pair.
std::optional<std::pair<std::size_t, std::size_t>> filw_contact_pair(const Mesh& rMesh) {
    if (rMesh.HasFieldData("abaqus:contact_pair")) {
        const NDArray& a = rMesh.FieldData("abaqus:contact_pair");
        const detail::Int64View v(a);
        bool ok = a.Size() == 2;
        for (std::size_t k = 0; ok && k < 2; ++k)
            ok = v.Data()[k] >= 0 && static_cast<std::size_t>(v.Data()[k]) < rMesh.NumRegions() &&
                 rMesh.Region(static_cast<std::size_t>(v.Data()[k])).mKind == RegionKind::Side;
        if (!ok)
            filw_fail("'abaqus:contact_pair' must hold the indices of two side regions");
        return std::make_pair(static_cast<std::size_t>(v.Data()[0]),
                              static_cast<std::size_t>(v.Data()[1]));
    }
    std::vector<std::size_t> sides;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mKind == RegionKind::Side)
            sides.push_back(i);
    if (sides.size() < 2)
        return std::nullopt;
    return std::make_pair(sides[0], sides[1]);
}

std::vector<std::int64_t> filw_labels(const std::vector<std::int64_t>* pGiven, std::size_t Count,
                                      const char* pWhat) {
    std::vector<std::int64_t> labels;
    if (pGiven == nullptr) {
        for (std::size_t i = 0; i < Count; ++i)
            labels.push_back(static_cast<std::int64_t>(i + 1));
        return labels;
    }
    labels = *pGiven;
    if (labels.size() != Count)
        filw_fail(std::string("'abaqus:id' must hold one label per ") + pWhat);
    for (const std::int64_t v : labels)
        if (v <= 0)
            filw_fail(std::string(pWhat) + " labels must be positive");
    std::vector<std::int64_t> sorted = labels;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        filw_fail(std::string(pWhat) + " labels are not unique");
    return labels;
}

using FilwSignature = std::vector<std::tuple<std::string, std::size_t, std::size_t>>;

FilwSignature filw_signature(const Mesh& rMesh) {
    FilwSignature sig;
    sig.emplace_back("points", rMesh.NumPoints(), 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        sig.emplace_back(cb.Type(), cb.NumCells(), cb.NodesPerCell());
    }
    return sig;
}

/// The mesh as the file holds it: records, labels and the written cells.
class FilwModel {
public:
    explicit FilwModel(const Mesh& rMesh) {
        mSignature = filw_signature(rMesh);
        // contact output refers to its surface pair by a 1940 label
        for (const std::string& name : rMesh.PointDataNames())
            if (filw_contact_key(name)) {
                if (const auto pair = filw_contact_pair(rMesh)) {
                    mContactNames.insert(rMesh.Region(pair->first).mName);
                    mContactNames.insert(rMesh.Region(pair->second).mName);
                }
                break;
            }
        const std::size_t n = rMesh.NumPoints();
        const std::size_t pdim = rMesh.PointDim();
        mDim = pdim == 2 ? 2 : 3;
        std::vector<std::int64_t> node_ids;
        const bool have_ids = rMesh.HasPointData("abaqus:id");
        if (have_ids) {
            const detail::Int64View v(rMesh.PointData("abaqus:id"));
            node_ids.assign(v.Data(), v.Data() + rMesh.PointData("abaqus:id").Size());
        }
        mNodeLabels = filw_labels(have_ids ? &node_ids : nullptr, n, "node");
        std::vector<double> points(n * mDim, 0.0);
        {
            const detail::DoubleView pts(rMesh.Points());
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t k = 0; k < mDim && k < pdim; ++k)
                    points[i * mDim + k] = pts.Data()[i * pdim + k];
        }
        std::size_t total = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            total += rMesh.Cells(b).NumCells();
        std::vector<std::int64_t> cell_ids;
        const bool have_cell_ids = rMesh.HasCellData("abaqus:id");
        if (have_cell_ids)
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                const NDArray& a = rMesh.CellData("abaqus:id", b);
                const detail::Int64View v(a);
                cell_ids.insert(cell_ids.end(), v.Data(), v.Data() + a.Size());
            }
        mCellLabels = filw_labels(have_cell_ids ? &cell_ids : nullptr, total, "element");
        std::size_t g = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const auto it = filw_types().find(cb.Type());
            for (std::size_t r = 0; r < cb.NumCells(); ++r, ++g) {
                if (it == filw_types().end())
                    ++mSkipped[cb.Type()];
                else
                    mCells.push_back({g, it->second, b, r});
            }
        }
        for (const Cell& c : mCells)
            mWritten.insert(c.mGlobal);
        BuildRecords(rMesh, points);
    }

    struct Cell {
        std::size_t mGlobal;
        std::string mName;
        std::size_t mBlock, mRow;
    };
    const FilwSignature& Signature() const { return mSignature; }
    const std::vector<FilwRecord>& Records() const { return mRecords; }
    const std::vector<std::int64_t>& NodeLabels() const { return mNodeLabels; }
    const std::vector<std::int64_t>& CellLabels() const { return mCellLabels; }
    const std::vector<Cell>& Cells() const { return mCells; }
    const std::map<std::string, std::size_t>& Skipped() const { return mSkipped; }
    std::set<std::string> mDropped;
    /// The 1940 label of each surface (side region) the contact pair names.
    std::map<std::string, std::int64_t> mSurfaceLabels;

private:
    std::set<std::string> mContactNames;

    void BuildRecords(const Mesh& rMesh, const std::vector<double>& rPoints) {
        FilwRecord head = filw_rec(1921);
        for (int k = 0; k < 4; ++k)
            head.push_back(k == 0 ? filw_a("meshio++") : filw_a(""));
        head.push_back(filw_i(static_cast<std::int64_t>(mCells.size())));
        head.push_back(filw_i(static_cast<std::int64_t>(mNodeLabels.size())));
        head.push_back(filw_d(1.0));
        mRecords.push_back(std::move(head));
        FilwRecord title = filw_rec(1922);
        filw_text_items(title, "meshio++ results file");
        mRecords.push_back(std::move(title));
        // the active degrees of freedom (1, 2 and, in 3-D, 3); 35 words
        FilwRecord dof = filw_rec(1902);
        const std::int64_t dofs[3] = {1, 2, mDim == 3 ? 3 : 0};
        for (const std::int64_t d : dofs)
            dof.push_back(filw_i(d));
        for (int k = 0; k < 32; ++k)
            dof.push_back(filw_i(0));
        mRecords.push_back(std::move(dof));
        constexpr std::size_t kChunk = 10;
        for (const Cell& c : mCells) {
            const auto cb = rMesh.Cells(c.mBlock);
            const detail::Int64View conn(cb.Conn());
            const std::size_t k = cb.NodesPerCell();
            std::vector<std::int64_t> labels;
            for (std::size_t j = 0; j < k; ++j) {
                const std::int64_t p = conn.Data()[c.mRow * k + j];
                if (p < 0 || static_cast<std::size_t>(p) >= mNodeLabels.size())
                    filw_fail("a cell refers to a point the mesh does not have");
                labels.push_back(mNodeLabels[static_cast<std::size_t>(p)]);
            }
            FilwRecord r = filw_rec(1900);
            r.push_back(filw_i(mCellLabels[c.mGlobal]));
            r.push_back(filw_a(c.mName));
            for (std::size_t j = 0; j < std::min(kChunk, labels.size()); ++j)
                r.push_back(filw_i(labels[j]));
            mRecords.push_back(std::move(r));
            for (std::size_t j = kChunk; j < labels.size(); j += kChunk) {
                FilwRecord more = filw_rec(1990);
                for (std::size_t m = j; m < std::min(j + kChunk, labels.size()); ++m)
                    more.push_back(filw_i(labels[m]));
                mRecords.push_back(std::move(more));
            }
        }
        for (std::size_t i = 0; i < mNodeLabels.size(); ++i) {
            FilwRecord r = filw_rec(1901);
            r.push_back(filw_i(mNodeLabels[i]));
            for (std::size_t k = 0; k < mDim; ++k)
                r.push_back(filw_d(rPoints[i * mDim + k]));
            mRecords.push_back(std::move(r));
        }
        Sets(rMesh);
        Matrices(rMesh);
    }

    void Sets(const Mesh& rMesh) {
        struct Ref {
            std::string mKind, mName;
            std::size_t mIndex;
        };
        std::vector<Ref> order;
        std::set<std::int64_t> used_digits;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const std::string& name = rMesh.Region(i).mName;
            order.push_back({region_kind_name(rMesh.Region(i).mKind), name, i});
            if (!name.empty() && name.size() < 18 &&
                std::all_of(name.begin(), name.end(), [](char C) { return C >= '0' && C <= '9'; }))
                used_digits.insert(std::stoll(name));
        }
        std::stable_sort(order.begin(), order.end(), [](const Ref& A, const Ref& B) {
            return std::tie(A.mKind, A.mName) < std::tie(B.mKind, B.mName);
        });
        std::vector<FilwRecord> refs;
        std::int64_t label = 0;
        constexpr std::size_t kChunk = 10;
        for (const Ref& ref : order) {
            const auto& region = rMesh.Region(ref.mIndex);
            std::string name = filw_ascii(region.mName);
            if (name.size() > 80)
                filw_fail("region name '" + name.substr(0, 20) +
                          "...' is longer than 80 characters");
            const bool side = region.mKind == RegionKind::Side;
            if (name.size() > 8 || (side && mContactNames.count(region.mName) != 0)) {
                ++label;
                while (used_digits.count(label) != 0)
                    ++label;
                FilwRecord r = filw_rec(1940);
                r.push_back(filw_i(label));
                filw_text_items(r, name);
                refs.push_back(std::move(r));
                if (side)
                    mSurfaceLabels[region.mName] = label;
                name = std::to_string(label);
            }
            if (region.mKind == RegionKind::Side) {
                Surface(rMesh, region, name);
                continue;
            }
            const detail::Int64View entries(region.mEntries);
            std::vector<std::int64_t> values;
            const bool point = region.mKind == RegionKind::Point;
            for (std::size_t k = 0; k < region.mEntries.Size(); ++k) {
                const std::int64_t e = entries.Data()[k];
                if (point) {
                    if (e < 0 || static_cast<std::size_t>(e) >= mNodeLabels.size())
                        filw_fail("region '" + region.mName +
                                  "' names a point the mesh does not have");
                    values.push_back(mNodeLabels[static_cast<std::size_t>(e)]);
                } else if (e >= 0 && mWritten.count(static_cast<std::size_t>(e)) != 0) {
                    values.push_back(mCellLabels[static_cast<std::size_t>(e)]);
                }
            }
            FilwRecord r = filw_rec(point ? 1931 : 1933);
            r.push_back(filw_a(name));
            for (std::size_t k = 0; k < std::min<std::size_t>(5, values.size()); ++k)
                r.push_back(filw_i(values[k]));
            mRecords.push_back(std::move(r));
            for (std::size_t k = 5; k < values.size(); k += kChunk) {
                FilwRecord more = filw_rec(point ? 1932 : 1934);
                for (std::size_t m = k; m < std::min(k + kChunk, values.size()); ++m)
                    more.push_back(filw_i(values[m]));
                mRecords.push_back(std::move(more));
            }
        }
        mRecords.insert(mRecords.end(), refs.begin(), refs.end());
    }

    /// A side region as a contact surface: a 1501 record and a 1502 record per (element, face).
    void Surface(const Mesh& rMesh, const Region& rRegion, const std::string& rName) {
        std::unordered_map<std::size_t, std::size_t> block_of;
        for (const Cell& c : mCells)
            block_of[c.mGlobal] = c.mBlock;
        FilwRecord head = filw_rec(1501);
        head.push_back(filw_a(rName));
        head.push_back(filw_i(static_cast<std::int64_t>(mDim)));
        head.push_back(filw_i(1));
        std::vector<FilwRecord> faces;  // the facet count and slave count follow the type
        std::size_t dropped = 0;
        const detail::Int64View entries(rRegion.mEntries);
        for (std::size_t k = 0; k + 1 < rRegion.mEntries.Size(); k += 2) {
            const auto cell = static_cast<std::size_t>(entries.Data()[k]);
            const std::int64_t facet = entries.Data()[k + 1];
            const auto it = block_of.find(cell);
            if (it == block_of.end()) {
                ++dropped;
                continue;
            }
            const std::string label =
                detail::abaqus_face_name(std::string(rMesh.Cells(it->second).Type()), facet);
            if (label.empty()) {
                ++dropped;
                continue;
            }
            FilwRecord r = filw_rec(1502);
            r.push_back(filw_i(mCellLabels[cell]));
            r.push_back(filw_i(std::stoll(label.substr(1))));
            faces.push_back(std::move(r));
        }
        if (dropped)
            mDropped.insert(std::to_string(dropped) + " facet(s) of surface '" + rRegion.mName +
                            "'");
        head.push_back(filw_i(static_cast<std::int64_t>(faces.size())));
        head.push_back(filw_i(0));
        mRecords.push_back(std::move(head));
        mRecords.insert(mRecords.end(), faces.begin(), faces.end());
    }

    /// The element matrices: per row of an `abaqus:<kind>:index` table (element, offset, length,
    /// flag) a 1001 record and its values.
    void Matrices(const Mesh& rMesh) {
        static const struct {
            const char* mKind;
            std::int64_t mKeyFlag, mKeyPlain;
        } kinds[] = {{"matrix_dofs", 1002, 1002},
                     {"stiffness", 1011, 1012},
                     {"mass", 1021, 1022},
                     {"load", 1031, 1031}};
        for (const auto& kind : kinds) {
            const std::string values_name = std::string("abaqus:") + kind.mKind;
            const std::string index_name = values_name + ":index";
            if (!rMesh.HasFieldData(values_name) || !rMesh.HasFieldData(index_name))
                continue;
            const NDArray& values = rMesh.FieldData(values_name);
            const NDArray& index = rMesh.FieldData(index_name);
            for (std::size_t row = 0; row + 3 < index.Size(); row += 4) {
                const std::int64_t element = detail::read_int(index, row);
                const std::int64_t offset = detail::read_int(index, row + 1);
                const std::int64_t length = detail::read_int(index, row + 2);
                const std::int64_t flag = detail::read_int(index, row + 3);
                if (offset < 0 || length < 0 ||
                    static_cast<std::size_t>(offset + length) > values.Size())
                    filw_fail("'" + index_name + "' reaches past '" + values_name + "'");
                FilwRecord head = filw_rec(1001);
                head.push_back(filw_i(element));
                mRecords.push_back(std::move(head));
                const std::string kind_name = kind.mKind;
                FilwRecord r = filw_rec(kind_name == "matrix_dofs" ? 1002
                                        : kind_name == "load"      ? 1031
                                        : flag == 1                ? kind.mKeyFlag
                                                                   : kind.mKeyPlain);
                if (kind_name == "load")
                    r.push_back(filw_i(flag));
                for (std::int64_t k = 0; k < length; ++k) {
                    const auto at = static_cast<std::size_t>(offset + k);
                    r.push_back(kind_name == "matrix_dofs"
                                    ? filw_i(detail::read_int(values, at))
                                    : filw_d(detail::read_double(values, at)));
                }
                mRecords.push_back(std::move(r));
            }
        }
    }

    FilwSignature mSignature;
    std::size_t mDim = 3;
    std::vector<std::int64_t> mNodeLabels, mCellLabels;
    std::vector<Cell> mCells;
    std::set<std::size_t> mWritten;
    std::map<std::string, std::size_t> mSkipped;
    std::vector<FilwRecord> mRecords;
};

/// A record's values: trailing NaN dropped, an inner NaN written as 0.
std::vector<double> filw_trim(const double* pValues, std::size_t N) {
    std::vector<double> v(pValues, pValues + N);
    while (!v.empty() && v.back() != v.back())
        v.pop_back();
    for (double& x : v)
        if (x != x)
            x = 0.0;
    return v;
}

struct FilwName {
    std::string mBase;
    std::optional<std::int64_t> mSection;
    std::optional<std::string> mRebar;
    std::string mCollision;
    bool mValid = true;
};

/// `(base, section point, rebar, collision suffix)` of a cell array name.
FilwName filw_split_name(const std::string& rName) {
    FilwName out;
    const std::size_t at = rName.find('@');
    out.mBase = rName.substr(0, at);
    if (at == std::string::npos)
        return out;
    std::size_t start = at;
    while (start < rName.size()) {
        const std::size_t next = rName.find('@', start + 1);
        const std::string part =
            rName.substr(start + 1, (next == std::string::npos ? rName.size() : next) - start - 1);
        start = next == std::string::npos ? rName.size() : next;
        if (part.empty())
            continue;
        if (part.compare(0, 2, "sp") == 0 && part.size() > 2 &&
            std::all_of(part.begin() + 2, part.end(),
                        [](char C) { return C >= '0' && C <= '9'; }) &&
            part.size() < 12)
            out.mSection = std::stoll(part.substr(2));
        else if (part.compare(0, 6, "rebar:") == 0)
            out.mRebar = part.substr(6);
        else if (part == "ip" || part == "el")
            out.mCollision = part;
        else {
            out.mBase = rName;
            out.mSection.reset();
            out.mRebar.reset();
            out.mCollision.clear();
            return out;
        }
    }
    return out;
}

void filw_ndi_nshr(std::size_t Count, std::int64_t& rNdi, std::int64_t& rNshr) {
    rNdi = Count == 6 ? 3 : Count == 4 ? 3 : Count == 3 ? 2 : 0;
    rNshr = Count == 6 ? 3 : Count == 4 ? 1 : Count == 3 ? 1 : 0;
}

bool filw_ends_with(const std::string& rS, const std::string& rSuffix) {
    return rS.size() >= rSuffix.size() &&
           rS.compare(rS.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
}

std::optional<double> filw_scalar(const Mesh& rMesh, const std::string& rName) {
    if (!rMesh.HasFieldData(rName))
        return std::nullopt;
    const NDArray& a = rMesh.FieldData(rName);
    if (a.Size() == 0)
        return std::nullopt;
    return detail::DoubleView(a).Data()[0];
}

std::optional<std::int64_t> filw_int_scalar(const Mesh& rMesh, const std::string& rName) {
    if (!rMesh.HasFieldData(rName))
        return std::nullopt;
    const NDArray& a = rMesh.FieldData(rName);
    if (a.Size() == 0)
        return std::nullopt;
    return detail::Int64View(a).Data()[0];
}

/// One increment: a 2000 record, the step's results and a 2001 record.
struct FilwIncrement {
    std::vector<FilwRecord> mRecords;
    std::set<std::string> mDropped;

    FilwIncrement(const FilwModel& rModel, const Mesh& rMesh, double Time, std::size_t Index) {
        const double step_time = filw_scalar(rMesh, "abaqus:step_time").value_or(Time);
        const std::int64_t procedure = filw_int_scalar(rMesh, "abaqus:procedure").value_or(1);
        const std::int64_t step = filw_int_scalar(rMesh, "abaqus:step").value_or(1);
        const std::int64_t increment = filw_int_scalar(rMesh, "abaqus:increment")
                                           .value_or(static_cast<std::int64_t>(Index + 1));
        FilwRecord head = filw_rec(2000);
        head.push_back(filw_d(Time));
        head.push_back(filw_d(step_time));
        head.push_back(filw_d(0.0));
        head.push_back(filw_d(0.0));
        head.push_back(filw_i(procedure));
        head.push_back(filw_i(step));
        head.push_back(filw_i(increment));
        head.push_back(filw_i(0));
        head.push_back(filw_d(0.0));
        head.push_back(filw_d(0.0));
        head.push_back(filw_d(1.0));
        filw_text_items(head, "");
        mRecords.push_back(std::move(head));
        Mode(rMesh);
        Energies(rMesh, procedure);
        Modal(rMesh);
        const std::vector<FilwRecord> element = Element(rModel, rMesh);
        const std::vector<FilwRecord> nodal = Nodal(rModel, rMesh);
        const std::vector<FilwRecord> contact = Contact(rModel, rMesh);
        if (!element.empty()) {
            FilwRecord request = filw_rec(1911);
            request.push_back(filw_i(0));
            request.push_back(filw_a(""));
            request.push_back(filw_a(""));
            mRecords.push_back(std::move(request));
            mRecords.insert(mRecords.end(), element.begin(), element.end());
        }
        if (!nodal.empty()) {
            FilwRecord request = filw_rec(1911);
            request.push_back(filw_i(1));
            request.push_back(filw_a(""));
            mRecords.push_back(std::move(request));
            mRecords.insert(mRecords.end(), nodal.begin(), nodal.end());
        }
        mRecords.insert(mRecords.end(), contact.begin(), contact.end());
        mRecords.push_back(filw_rec(2001));
        static const std::set<std::string> known = {
            kSequenceTimeKey,           "abaqus:step",
            "abaqus:increment",         "abaqus:step_time",
            "abaqus:procedure",         "abaqus:mode",
            "abaqus:eigenvalue",        "abaqus:generalized_mass",
            "abaqus:composite_damping", "abaqus:participation_factor",
            "abaqus:effective_mass",    "abaqus:contact_pair"};
        std::set<std::string> handled;
        for (const bool explicit_names : {false, true})
            for (std::size_t k = 0; k < 18; ++k)
                if (const char* n = fil_energy_name(k, explicit_names))
                    handled.insert(std::string("abaqus:") + n);
        for (std::int64_t key = 301; key <= 310; ++key)
            handled.insert(std::string("abaqus:") + fil_modal_name(key));
        for (const char* kind : {"matrix_dofs", "stiffness", "mass", "load"}) {
            handled.insert(std::string("abaqus:") + kind);
            handled.insert(std::string("abaqus:") + kind + ":index");
        }
        for (const std::string& name : rMesh.FieldDataNames())
            if (known.count(name) == 0 && handled.count(name) == 0 &&
                name.compare(0, 14, "abaqus:layout:") != 0)
                mDropped.insert(name);
    }

private:
    /// Record 1999, the total energies, positional: those the step's procedure names, up to the
    /// last one present (an absent one in between is 0).
    void Energies(const Mesh& rMesh, std::int64_t Procedure) {
        const bool explicit_names = Procedure == 17 || Procedure == 21 || Procedure == 74;
        std::vector<std::pair<std::size_t, double>> present;
        for (std::size_t k = 0; k < 18; ++k) {
            const char* name = fil_energy_name(k, explicit_names);
            if (!name)
                continue;
            const std::string key = std::string("abaqus:") + name;
            if (rMesh.HasFieldData(key) && rMesh.FieldData(key).Size() > 0)
                present.emplace_back(k, detail::read_double(rMesh.FieldData(key), 0));
        }
        if (present.empty())
            return;
        std::vector<double> values(present.back().first + 1, 0.0);
        for (const auto& [k, v] : present)
            values[k] = v;
        FilwRecord r = filw_rec(1999);
        for (const double v : values)
            r.push_back(filw_d(v));
        mRecords.push_back(std::move(r));
    }

    /// Records 301-310, the generalized quantities of a mode-based dynamic step: one record per
    /// row of `abaqus:GU` and the like (trailing NaN trimmed, an inner NaN written as 0). BM's
    /// row is the motion type and six components, followed by the base name (blank).
    void Modal(const Mesh& rMesh) {
        for (std::int64_t key = 301; key <= 310; ++key) {
            const std::string name = std::string("abaqus:") + fil_modal_name(key);
            if (!rMesh.HasFieldData(name))
                continue;
            const NDArray& a = rMesh.FieldData(name);
            if (a.Shape().size() > 2)
                filw_fail("'" + name + "' must be one row or a matrix of rows");
            const std::size_t width = a.Shape().empty() ? 1 : a.Shape().back();
            const std::size_t rows = width == 0 ? 0 : a.Size() / width;
            const detail::DoubleView v(a);
            for (std::size_t row = 0; row < rows; ++row) {
                std::vector<double> values = filw_trim(v.Data() + row * width, width);
                if (values.empty())
                    continue;
                FilwRecord r = filw_rec(key);
                if (key == 304) {
                    values.resize(std::min<std::size_t>(values.size(), 7));
                    r.push_back(filw_i(static_cast<std::int64_t>(values[0])));
                    for (std::size_t k = 1; k < values.size(); ++k)
                        r.push_back(filw_d(values[k]));
                    if (values.size() == 7)
                        r.push_back(filw_a(""));
                } else {
                    for (const double x : values)
                        r.push_back(filw_d(x));
                }
                mRecords.push_back(std::move(r));
            }
        }
    }

    void Mode(const Mesh& rMesh) {
        const auto mode = filw_int_scalar(rMesh, "abaqus:mode");
        if (!mode)
            return;
        std::vector<double> factors, masses;
        if (rMesh.HasFieldData("abaqus:participation_factor")) {
            const NDArray& a = rMesh.FieldData("abaqus:participation_factor");
            const detail::DoubleView v(a);
            factors.assign(v.Data(), v.Data() + a.Size());
        }
        if (rMesh.HasFieldData("abaqus:effective_mass")) {
            const NDArray& a = rMesh.FieldData("abaqus:effective_mass");
            const detail::DoubleView v(a);
            masses.assign(v.Data(), v.Data() + a.Size());
        }
        if (factors.size() != masses.size())
            filw_fail("'abaqus:participation_factor' and 'abaqus:effective_mass' differ");
        FilwRecord r = filw_rec(1980);
        r.push_back(filw_i(*mode));
        for (const char* name :
             {"abaqus:eigenvalue", "abaqus:generalized_mass", "abaqus:composite_damping"})
            r.push_back(filw_d(filw_scalar(rMesh, name).value_or(0.0)));
        for (std::size_t k = 0; k < factors.size(); ++k) {
            r.push_back(filw_d(factors[k]));
            r.push_back(filw_d(masses[k]));
        }
        mRecords.push_back(std::move(r));
    }

    std::vector<FilwRecord> Nodal(const FilwModel& rModel, const Mesh& rMesh) {
        std::vector<FilwRecord> out;
        const std::size_t n = rModel.NodeLabels().size();
        for (const std::string& name : rMesh.PointDataNames()) {
            if (name == "abaqus:id" || filw_ends_with(name, "@avg"))
                continue;
            if (filw_contact_key(name))
                continue;  // contact output, written by Contact
            const auto key = filw_key_of(name, true);
            if (!key) {
                if (!filw_key_of(name, false))
                    mDropped.insert(name);
                continue;
            }
            const NDArray& a = rMesh.PointData(name);
            if (a.Size() == 0 || a.Shape().empty() || a.Shape()[0] != n) {
                mDropped.insert(name);
                continue;
            }
            const detail::DoubleView v(a);
            const std::size_t cols = a.Size() / n;
            for (std::size_t i = 0; i < n; ++i) {
                const std::vector<double> values = filw_trim(v.Data() + i * cols, cols);
                if (values.empty())
                    continue;
                FilwRecord r = filw_rec(*key);
                r.push_back(filw_i(rModel.NodeLabels()[i]));
                for (const double x : values)
                    r.push_back(filw_d(x));
                out.push_back(std::move(r));
            }
        }
        return out;
    }

    /// Contact output: a 1503 request naming the slave and master surface by their 1940 labels,
    /// then per node holding data a 1504 header and a record per array.
    std::vector<FilwRecord> Contact(const FilwModel& rModel, const Mesh& rMesh) {
        std::vector<std::pair<std::int64_t, std::string>> arrays;
        for (const std::string& name : rMesh.PointDataNames())
            if (const auto key = filw_contact_key(name))
                arrays.emplace_back(*key, name);
        std::vector<FilwRecord> out;
        if (arrays.empty())
            return out;
        std::sort(arrays.begin(), arrays.end());
        std::int64_t labels[2] = {0, 0};
        if (const auto pair = filw_contact_pair(rMesh)) {
            const auto a = rModel.mSurfaceLabels.find(rMesh.Region(pair->first).mName);
            const auto b = rModel.mSurfaceLabels.find(rMesh.Region(pair->second).mName);
            if (a != rModel.mSurfaceLabels.end() && b != rModel.mSurfaceLabels.end()) {
                labels[0] = a->second;
                labels[1] = b->second;
            }
        }
        if (labels[0] == 0) {
            for (const auto& entry : arrays)
                mDropped.insert(entry.second);
            return out;
        }
        const std::size_t n = rModel.NodeLabels().size();
        std::vector<const NDArray*> data;
        std::vector<std::int64_t> keys;
        for (const auto& [key, name] : arrays) {
            const NDArray& a = rMesh.PointData(name);
            if (a.Size() == 0 || a.Shape().empty() || a.Shape()[0] != n) {
                mDropped.insert(name);
                continue;
            }
            data.push_back(&a);
            keys.push_back(key);
        }
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<FilwRecord> items;
            for (std::size_t k = 0; k < data.size(); ++k) {
                const std::size_t cols = data[k]->Size() / n;
                const detail::DoubleView v(*data[k]);
                const std::vector<double> values = filw_trim(v.Data() + i * cols, cols);
                if (values.empty())
                    continue;
                FilwRecord r = filw_rec(keys[k]);
                for (const double x : values)
                    r.push_back(filw_d(x));
                items.push_back(std::move(r));
            }
            if (items.empty())
                continue;
            FilwRecord h = filw_rec(1504);
            h.push_back(filw_i(rModel.NodeLabels()[i]));
            h.push_back(filw_i(static_cast<std::int64_t>(items.front().size()) - 1));
            out.push_back(std::move(h));
            for (FilwRecord& r : items)
                out.push_back(std::move(r));
        }
        if (!out.empty()) {
            FilwRecord request = filw_rec(1503);
            request.push_back(filw_i(0));
            request.push_back(filw_a(std::to_string(labels[0])));
            request.push_back(filw_a(std::to_string(labels[1])));
            request.push_back(filw_a(""));
            out.insert(out.begin(), std::move(request));
        }
        return out;
    }

    std::vector<FilwRecord> Element(const FilwModel& rModel, const Mesh& rMesh) {
        std::map<std::size_t, std::size_t> index;  // global cell -> position
        for (std::size_t k = 0; k < rModel.Cells().size(); ++k)
            index[rModel.Cells()[k].mGlobal] = k;
        using SlotKey = std::tuple<std::size_t, int, std::int64_t, std::int64_t, std::string>;
        std::map<SlotKey, std::vector<std::pair<std::int64_t, std::vector<double>>>> slots;
        std::vector<std::size_t> offsets{0};
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            offsets.push_back(offsets.back() + rMesh.Cells(b).NumCells());
        for (const std::string& name : rMesh.CellDataNames()) {
            if (name == "abaqus:id")
                continue;
            const FilwName parsed = filw_split_name(name);
            const auto key = filw_key_of(parsed.mBase, false);
            if (!key) {
                mDropped.insert(name);
                continue;
            }
            std::vector<std::int64_t> layout;
            if (rMesh.HasFieldData("abaqus:layout:" + name)) {
                const NDArray& l = rMesh.FieldData("abaqus:layout:" + name);
                const detail::Int64View v(l);
                layout.assign(v.Data(), v.Data() + l.Size());
                if (layout.size() != 2 || layout[0] < 1 || layout[1] < 1)
                    filw_fail("'abaqus:layout:" + name + "' must be (points, components)");
            }
            int location = 1;
            if (parsed.mRebar)
                location = 3;
            else if (!layout.empty())
                location = (parsed.mBase == "NFORC" || parsed.mCollision == "el") ? 2 : 0;
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                const auto cb = rMesh.Cells(b);
                const std::size_t n = cb.NumCells();
                const NDArray& a = rMesh.CellData(name, b);
                if (n == 0 || a.Shape().empty() || a.Shape()[0] != n)
                    continue;
                const std::size_t cols = a.Size() / n;
                const detail::DoubleView v(a);
                const std::int64_t sec =
                    parsed.mSection ? *parsed.mSection : (filw_is_shell(cb.Type()) ? 1 : 0);
                if (!layout.empty() && cols != static_cast<std::size_t>(layout[0]) *
                                                   static_cast<std::size_t>(layout[1]))
                    filw_fail("'" + name + "' does not match its 'abaqus:layout:" + name + "'");
                const detail::Int64View conn(cb.Conn());
                for (std::size_t r = 0; r < n; ++r) {
                    const auto at = index.find(offsets[b] + r);
                    if (at == index.end())
                        continue;
                    const double* row = v.Data() + r * cols;
                    if (layout.empty()) {
                        const std::vector<double> values = filw_trim(row, cols);
                        if (!values.empty())
                            slots[{at->second, location, sec, 0, ""}].emplace_back(*key, values);
                        continue;
                    }
                    for (std::int64_t p = 0; p < layout[0]; ++p) {
                        const std::vector<double> values =
                            filw_trim(row + p * layout[1], static_cast<std::size_t>(layout[1]));
                        if (values.empty())
                            continue;
                        std::int64_t point = p + 1;
                        if (location == 2) {
                            if (static_cast<std::size_t>(p) >= cb.NodesPerCell())
                                continue;
                            point = rModel.NodeLabels()[static_cast<std::size_t>(
                                conn.Data()[r * cb.NodesPerCell() + static_cast<std::size_t>(p)])];
                        }
                        slots[{at->second, location, sec, point, parsed.mRebar.value_or("")}]
                            .emplace_back(*key, values);
                    }
                }
            }
        }
        std::vector<FilwRecord> out;
        for (const auto& [slot, items] : slots) {
            const auto& [position, location, sec, point, rebar] = slot;
            const std::size_t g = rModel.Cells()[position].mGlobal;
            std::int64_t ndi = 0, nshr = 0;
            for (const auto& item : items) {
                const std::size_t c = item.second.size();
                if (c == 3 || c == 4 || c == 6) {
                    filw_ndi_nshr(c, ndi, nshr);
                    break;
                }
            }
            FilwRecord h = filw_rec(1);
            h.push_back(filw_i(rModel.CellLabels()[g]));
            h.push_back(filw_i(point));
            h.push_back(filw_i(sec));
            h.push_back(filw_i(location));
            h.push_back(filw_a(rebar));
            h.push_back(filw_i(ndi));
            h.push_back(filw_i(nshr));
            h.push_back(filw_i(0));
            h.push_back(filw_i(0));
            out.push_back(std::move(h));
            for (const auto& [key, values] : items) {
                FilwRecord r = filw_rec(key);
                for (const double x : values)
                    r.push_back(filw_d(x));
                out.push_back(std::move(r));
            }
        }
        Averaged(rModel, rMesh, out);
        return out;
    }

    /// Point data named like an element result: averaged at the nodes (location 4).
    void Averaged(const FilwModel& rModel, const Mesh& rMesh, std::vector<FilwRecord>& rOut) {
        const std::size_t n = rModel.NodeLabels().size();
        for (const std::string& name : rMesh.PointDataNames()) {
            std::string base = name;
            while (filw_ends_with(base, "@avg"))
                base.erase(base.size() - 4);
            if (base == "abaqus:id")
                continue;
            if (base == name && filw_contact_key(base))
                continue;  // contact output
            if (base == name && filw_key_of(base, true))
                continue;  // a nodal result
            const auto key = filw_key_of(base, false);
            if (!key)
                continue;
            const NDArray& a = rMesh.PointData(name);
            if (a.Size() == 0 || a.Shape().empty() || a.Shape()[0] != n)
                continue;
            const detail::DoubleView v(a);
            const std::size_t cols = a.Size() / n;
            for (std::size_t i = 0; i < n; ++i) {
                const std::vector<double> values = filw_trim(v.Data() + i * cols, cols);
                if (values.empty())
                    continue;
                FilwRecord h = filw_rec(1);
                h.push_back(filw_i(rModel.NodeLabels()[i]));
                for (int k = 0; k < 2; ++k)
                    h.push_back(filw_i(0));
                h.push_back(filw_i(4));
                h.push_back(filw_a(""));
                for (int k = 0; k < 4; ++k)
                    h.push_back(filw_i(0));
                rOut.push_back(std::move(h));
                FilwRecord r = filw_rec(*key);
                for (const double x : values)
                    r.push_back(filw_d(x));
                rOut.push_back(std::move(r));
            }
        }
    }
};

}  // namespace

namespace {

std::string filw_encode(const AbaqusFilWriteOptions& rOptions,
                        const std::vector<FilwRecord>& rRecords) {
    return rOptions.mBinary
               ? filw_binary_bytes(rRecords, rOptions.mBigEndian, rOptions.mMarkerSize)
               : filw_ascii_bytes(rRecords);
}

bool filw_has_results(const Mesh& rMesh) {
    if (rMesh.HasFieldData(kSequenceTimeKey))
        return true;
    for (const std::string& n : rMesh.PointDataNames())
        if (n != "abaqus:id")
            return true;
    for (const std::string& n : rMesh.CellDataNames())
        if (n != "abaqus:id")
            return true;
    return false;
}

double filw_time_of(const Mesh& rMesh) {
    return filw_scalar(rMesh, kSequenceTimeKey).value_or(0.0);
}

void filw_warn_dropped(const FilwModel& rModel, const std::vector<FilwIncrement>& rIncrements) {
    std::set<std::string> names = rModel.mDropped;
    for (const FilwIncrement& inc : rIncrements)
        names.insert(inc.mDropped.begin(), inc.mDropped.end());
    if (!names.empty()) {
        std::string list;
        for (const std::string& n : names)
            list += (list.empty() ? "" : ", ") + n;
        log::warn("Abaqus .fil: arrays and sets without a place in the file are not written: {}",
                  list);
    }
    if (!rModel.Skipped().empty()) {
        std::string list;
        for (const auto& [type, count] : rModel.Skipped())
            list += (list.empty() ? "" : ", ") + type;
        log::warn("Abaqus .fil: cell blocks with no Abaqus element type are not written: {}", list);
    }
}

void filw_save(const std::string& rPath, const std::string& rBytes) {
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("Could not open file for writing: " + rPath);
    out.write(rBytes.data(), static_cast<std::streamsize>(rBytes.size()));
    out.close();
    if (!out)
        filw_fail("failed writing " + rPath);
}

}  // namespace

void write_abaqus_fil(const std::string& rPath, const Mesh& rMesh, bool Binary) {
    AbaqusFilWriteOptions options;
    options.mBinary = Binary;
    write_abaqus_fil(rPath, rMesh, options);
}

void write_abaqus_fil(const std::string& rPath, const Mesh& rMesh,
                      const AbaqusFilWriteOptions& rOptions) {
    const detail::ProvenanceSlotlessWrite slotless;
    const FilwModel model(rMesh);
    std::vector<FilwIncrement> increments;
    if (filw_has_results(rMesh))
        increments.emplace_back(model, rMesh, filw_time_of(rMesh), 0);
    filw_warn_dropped(model, increments);
    std::vector<FilwRecord> records = model.Records();
    for (const FilwIncrement& inc : increments)
        records.insert(records.end(), inc.mRecords.begin(), inc.mRecords.end());
    filw_save(rPath, filw_encode(rOptions, records));
}

struct AbaqusFilSeriesWriter::Impl {
    std::string mPath;
    AbaqusFilWriteOptions mOptions;
    std::optional<FilwModel> mModel;
    std::vector<FilwIncrement> mIncrements;
    bool mDone = false;
};

AbaqusFilSeriesWriter::AbaqusFilSeriesWriter(const std::string& rPath, bool Binary)
    : AbaqusFilSeriesWriter(rPath, [&] {
          AbaqusFilWriteOptions options;
          options.mBinary = Binary;
          return options;
      }()) {}

AbaqusFilSeriesWriter::AbaqusFilSeriesWriter(const std::string& rPath,
                                             const AbaqusFilWriteOptions& rOptions)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
    mpImpl->mOptions = rOptions;
    // Fail now, not when the last step is written, if the path cannot be written.
    auto probe = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!probe)
        throw WriteError("Could not open file for writing: " + rPath);
}

AbaqusFilSeriesWriter::~AbaqusFilSeriesWriter() = default;
AbaqusFilSeriesWriter::AbaqusFilSeriesWriter(AbaqusFilSeriesWriter&&) noexcept = default;
AbaqusFilSeriesWriter& AbaqusFilSeriesWriter::operator=(AbaqusFilSeriesWriter&&) noexcept = default;

void AbaqusFilSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        filw_fail("series: the writer was moved from");
    Impl& s = *mpImpl;
    const std::size_t index = s.mIncrements.size();
    if (!s.mModel) {
        s.mModel.emplace(rMesh);
    } else if (s.mModel->Signature() != filw_signature(rMesh)) {
        filw_fail("series: step " + std::to_string(index) +
                  "'s points or cells differ from the first step's; a results file holds one "
                  "model -- write one file per step with '{step}'");
    }
    s.mIncrements.emplace_back(*s.mModel, rMesh, Time, index);
}

std::size_t AbaqusFilSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mIncrements.size() : 0;
}

void AbaqusFilSeriesWriter::Finalize() {
    if (!mpImpl || mpImpl->mDone)
        return;
    if (mpImpl->mIncrements.empty())
        filw_fail("series: no step was written to " + mpImpl->mPath);
    mpImpl->mDone = true;
    filw_warn_dropped(*mpImpl->mModel, mpImpl->mIncrements);
    std::vector<FilwRecord> records = mpImpl->mModel->Records();
    for (const FilwIncrement& inc : mpImpl->mIncrements)
        records.insert(records.end(), inc.mRecords.begin(), inc.mRecords.end());
    filw_save(mpImpl->mPath, filw_encode(mpImpl->mOptions, records));
}

}  // namespace meshioplusplus
