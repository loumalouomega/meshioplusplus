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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

const std::string kOp2Who = "Nastran OP2";
const double kOp2Nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t kOp2Npos = std::numeric_limits<std::size_t>::max();

[[noreturn]] void op2_fail(const std::string& rMessage) {
    throw ReadError(kOp2Who + ": " + rMessage);
}

// --- words ---------------------------------------------------------------------------

struct Op2Words {
    int mWs = 4;
    bool mSwap = false;

    std::int64_t Int(const char* p) const {
        if (mWs == 4) {
            std::uint32_t u;
            std::memcpy(&u, p, 4);
            if (mSwap)
                u = detail::bswap32(u);
            return static_cast<std::int32_t>(u);
        }
        std::uint64_t u;
        std::memcpy(&u, p, 8);
        if (mSwap)
            u = detail::bswap64(u);
        return static_cast<std::int64_t>(u);
    }

    double Float(const char* p) const {
        if (mWs == 4) {
            std::uint32_t u;
            std::memcpy(&u, p, 4);
            if (mSwap)
                u = detail::bswap32(u);
            return static_cast<double>(std::bit_cast<float>(u));
        }
        return Double(p);
    }

    double Double(const char* p) const {
        std::uint64_t u;
        std::memcpy(&u, p, 8);
        if (mSwap)
            u = detail::bswap64(u);
        return std::bit_cast<double>(u);
    }

    std::vector<std::int64_t> Ints(const std::string& rRaw) const {
        const std::size_t n = rRaw.size() / static_cast<std::size_t>(mWs);
        std::vector<std::int64_t> out(n);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = Int(rRaw.data() + i * static_cast<std::size_t>(mWs));
        return out;
    }
};

// --- framing -------------------------------------------------------------------------

// A data record: one or more Fortran blocks, joined on demand.
struct Op2Record {
    std::vector<std::pair<std::size_t, std::size_t>> mSpans;  // (offset, size)

    std::size_t Size() const {
        std::size_t n = 0;
        for (const auto& s : mSpans)
            n += s.second;
        return n;
    }

    std::string Join(const char* pBase) const {
        std::string out;
        out.reserve(Size());
        for (const auto& [off, size] : mSpans)
            out.append(pBase + off, size);
        return out;
    }
};

class Op2Stream {
public:
    Op2Stream(const char* pData, std::size_t Size) : mpData(pData) {
        const auto layout = detail::sniff_fortran_records(pData, Size);
        if (!layout)
            op2_fail("the file is not Fortran unformatted (no record framing)");
        mBlocks = detail::fortran_records(pData, Size, *layout, kOp2Who);
        const std::size_t first = mBlocks.front().mSize;
        if (first != 4 && first != 8)
            op2_fail("the first record is not a one-word marker");
        mWords.mWs = static_cast<int>(first);
        mWords.mSwap = layout->mBigEndian != (std::endian::native == std::endian::big);
    }

    const Op2Words& Words() const { return mWords; }
    const char* Base() const { return mpData; }
    bool Done() const { return mI >= mBlocks.size(); }

    std::int64_t Peek() const { return MarkerValue(mI); }

    std::int64_t Marker() {
        const std::int64_t v = MarkerValue(mI);
        ++mI;
        return v;
    }

    std::pair<std::size_t, std::size_t> Block() {
        if (Done())
            op2_fail("the file is truncated");
        const auto& b = mBlocks[mI++];
        return {b.mOffset, b.mSize};
    }

    Op2Record Record() {
        const std::int64_t n = Marker();
        if (n <= 0)
            op2_fail("expected a record, found marker " + std::to_string(n));
        Op2Record r;
        r.mSpans.push_back(Block());
        while (!Done()) {
            const auto& b = mBlocks[mI];
            if (b.mSize != static_cast<std::size_t>(mWords.mWs) || MarkerValue(mI) <= 0)
                break;
            ++mI;
            r.mSpans.push_back(Block());
        }
        return r;
    }

private:
    std::int64_t MarkerValue(std::size_t Index) const {
        if (Index >= mBlocks.size())
            op2_fail("the file is truncated");
        const auto& b = mBlocks[Index];
        if (b.mSize != static_cast<std::size_t>(mWords.mWs))
            op2_fail("expected a marker at byte " + std::to_string(b.mOffset - 4) + ", found a " +
                     std::to_string(b.mSize) + "-byte record");
        return mWords.Int(mpData + b.mOffset);
    }

    const char* mpData;
    std::vector<detail::FortranRecord> mBlocks;
    Op2Words mWords;
    std::size_t mI = 0;
};

// A Nastran string: 8-byte words hold 4 characters then 4 blanks in NX 64-bit
// files (interlaced) and 8 characters in MSC ones.
std::string op2_text(const std::string& rRaw, int Ws) {
    if (Ws == 8 && rRaw.size() % 8 == 0 && rRaw.size() >= 16) {
        bool interlaced = rRaw.compare(8, 4, "    ") != 0;
        for (std::size_t k = 0; interlaced && k < rRaw.size(); k += 8)
            interlaced = rRaw.compare(k + 4, 4, "    ") == 0;
        if (interlaced) {
            std::string out;
            for (std::size_t k = 0; k < rRaw.size(); k += 8)
                out.append(rRaw, k, 4);
            return out;
        }
    }
    return rRaw;
}

std::string op2_strip(const std::string& rS) {
    const char* ws = " \t\r\n\v\f";
    const auto first = rS.find_first_not_of(ws);
    if (first == std::string::npos)
        return {};
    return rS.substr(first, rS.find_last_not_of(ws) - first + 1);
}

struct Op2Table {
    std::string mName;
    std::vector<std::pair<std::int64_t, Op2Record>> mRecords;  // (mark, record)
};

std::vector<Op2Table> op2_parse(Op2Stream& rS) {
    if (rS.Peek() == 3) {  // PARAM,POST,-1 header: date, tape code, version
        rS.Marker();
        rS.Block();
        rS.Marker();
        rS.Block();
        rS.Record();
        if (rS.Marker() != -1 || rS.Marker() != 0)
            op2_fail("the file header is not closed by the markers -1, 0");
    }
    std::vector<Op2Table> tables;
    while (!rS.Done()) {
        std::int64_t m = rS.Peek();
        if (m == 0) {
            rS.Marker();
            continue;
        }
        if (m < 0)
            op2_fail("expected a table name, found marker " + std::to_string(m));
        Op2Table t;
        t.mName = op2_strip(op2_text(rS.Record().Join(rS.Base()), rS.Words().mWs).substr(0, 8));
        std::int64_t mark = 0;
        while (true) {
            if (rS.Done())
                op2_fail("table " + t.mName + " is truncated");
            m = rS.Peek();
            if (m == 0) {
                rS.Marker();
                break;
            }
            if (m == -1) {
                rS.Marker();
                mark = -1;
                continue;
            }
            if (m < 0) {
                const std::int64_t a = rS.Marker(), b = rS.Marker(), c = rS.Marker();
                if (b != 1 || (c != 0 && c != 1))
                    op2_fail("table " + t.mName + ": malformed marker triple (" +
                             std::to_string(a) + ", " + std::to_string(b) + ", " +
                             std::to_string(c) + ")");
                mark = m;
                continue;
            }
            t.mRecords.emplace_back(mark, rS.Record());
        }
        tables.push_back(std::move(t));
    }
    return tables;
}

// --- geometry ------------------------------------------------------------------------

struct Op2ElementSpec {
    int mK1, mK2;
    const char* mCard;
    int mSizes[3];
    int mNumSizes;
    int mFirst;    // word of the first node
    int mCount;    // nodes
    int mPidWord;  // -1: none
};

// Where versions wrote different sizes (NX and MSC CQUAD4: 14 and 15 words) the
// first size whose entries validate (increasing ids, corners that are GRIDs)
// wins, in pyNastran's order.
constexpr Op2ElementSpec kOp2Elements[] = {
    {2408, 24, "CBAR", {16}, 1, 2, 2, 1},
    {5408, 54, "CBEAM", {18}, 1, 2, 2, 1},
    {2608, 26, "CBUSH", {14}, 1, 2, 2, 1},
    {1501, 15, "CONM2", {13}, 1, 1, 1, -1},
    {1601, 16, "CONROD", {8}, 1, 1, 2, -1},
    {4108, 41, "CPENTA", {17}, 1, 2, 15, 1},
    {17200, 172, "CPYRAM", {16}, 1, 2, 13, 1},
    {9108, 91, "CQUAD", {11}, 1, 2, 9, 1},
    {2958, 51, "CQUAD4", {14, 15}, 2, 2, 4, 1},
    {4701, 47, "CQUAD8", {17, 16, 18}, 3, 2, 8, 1},
    {8009, 80, "CQUADR", {14, 15}, 2, 2, 4, 1},
    {3001, 30, "CROD", {4}, 1, 2, 2, 1},
    {3101, 31, "CSHEAR", {6}, 1, 2, 4, 1},
    {7308, 73, "CHEXA", {22}, 1, 2, 20, 1},
    {5508, 55, "CTETRA", {12}, 1, 2, 10, 1},
    {5959, 59, "CTRIA3", {13, 14}, 2, 2, 3, 1},
    {4801, 48, "CTRIA6", {13, 14, 15}, 3, 2, 6, 1},
    {9200, 92, "CTRIAR", {13, 14}, 2, 2, 3, 1},
    {3701, 37, "CTUBE", {4}, 1, 2, 2, 1},
    {3901, 39, "CVISC", {4}, 1, 2, 2, 1},
    {5201, 52, "PLOTEL", {3}, 1, 1, 2, -1},
    // NX 2019 and later write CQUAD4/CTRIA3 as these (pyNastran's CQUADRN/CTRIARN)
    {15401, 154, "CQUAD4", {14, 15}, 2, 2, 4, 1},
    {15301, 153, "CTRIA3", {13, 14}, 2, 2, 3, 1},
};

// Element records without a cell type, named in the warning.
constexpr std::tuple<int, int, const char*> kOp2NamedElements[] = {
    {601, 6, "CELAS1"},   {701, 7, "CELAS2"},     {801, 8, "CELAS3"},    {901, 9, "CELAS4"},
    {201, 2, "CDAMP1"},   {301, 3, "CDAMP2"},     {401, 4, "CDAMP3"},    {501, 5, "CDAMP4"},
    {1001, 10, "CMASS1"}, {1101, 11, "CMASS2"},   {1201, 12, "CMASS3"},  {1301, 13, "CMASS4"},
    {1401, 14, "CONM1"},  {1908, 19, "CGAP"},     {5608, 56, "CBUSH1D"}, {6108, 61, "CTRIAX6"},
    {9008, 90, "CQUADX"}, {10108, 101, "CTRIAX"},
};

enum class Op2PropLayout { Fixed, Sizes, Pcomp, Pcompg, Terminated };

struct Op2PropertySpec {
    int mK1, mK2;
    const char* mCard;
    Op2PropLayout mLayout;
    int mSizes[4];
    int mNumSizes;
};

constexpr Op2PropertySpec kOp2Properties[] = {
    {2302, 23, "PSHELL", Op2PropLayout::Fixed, {11}, 1},
    {2402, 24, "PSOLID", Op2PropLayout::Fixed, {7}, 1},
    {902, 9, "PROD", Op2PropLayout::Fixed, {6}, 1},
    {1002, 10, "PSHEAR", Op2PropLayout::Fixed, {6}, 1},
    {1602, 16, "PTUBE", Op2PropLayout::Fixed, {5}, 1},
    {52, 20, "PBAR", Op2PropLayout::Fixed, {19}, 1},
    {5402, 54, "PBEAM", Op2PropLayout::Fixed, {197}, 1},
    {302, 3, "PELAS", Op2PropLayout::Fixed, {4}, 1},
    {1802, 18, "PVISC", Op2PropLayout::Fixed, {3}, 1},
    {202, 2, "PDAMP", Op2PropLayout::Fixed, {2}, 1},
    {402, 4, "PMASS", Op2PropLayout::Fixed, {2}, 1},
    {2102, 21, "PGAP", Op2PropLayout::Fixed, {11}, 1},
    {4706, 47, "PLSOLID", Op2PropLayout::Fixed, {7}, 1},
    {4606, 46, "PLPLANE", Op2PropLayout::Fixed, {11}, 1},
    {1402, 14, "PBUSH", Op2PropLayout::Sizes, {18, 23, 24, 27}, 4},
    {2706, 27, "PCOMP", Op2PropLayout::Pcomp, {}, 0},
    {15006, 150, "PCOMPG", Op2PropLayout::Pcompg, {}, 0},
    {9102, 91, "PBARL", Op2PropLayout::Terminated, {}, 0},
    {9202, 92, "PBEAML", Op2PropLayout::Terminated, {}, 0},
};

// A card record: its key and its raw bytes (key included).
struct Op2CardRecord {
    std::int64_t mK1, mK2, mK3;
    std::string mRaw;
};

std::vector<Op2CardRecord> op2_geometry_records(const Op2Stream& rS,
                                                const std::vector<Op2Table>& rTables,
                                                const std::vector<std::string>& rPrefixes) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    std::vector<Op2CardRecord> out;
    for (const Op2Table& t : rTables) {
        bool match = false;
        for (const std::string& p : rPrefixes)
            match = match || t.mName.compare(0, p.size(), p) == 0;
        if (!match)
            continue;
        for (const auto& [mark, rec] : t.mRecords) {
            const std::size_t size = rec.Size();
            if (mark > -3 || size < 3 * ws || size % ws)
                continue;
            std::string raw = rec.Join(rS.Base());
            out.push_back({w.Int(raw.data()), w.Int(raw.data() + ws), w.Int(raw.data() + 2 * ws),
                           std::move(raw)});
        }
    }
    return out;
}

struct Op2Grids {
    std::vector<std::int64_t> mIds, mCp, mCd;
    std::vector<double> mXyz;
    std::unordered_set<std::int64_t> mScalarPoints;
    std::vector<detail::NastranCoordCard> mCords;
};

// CORD1R/C/S entries are `cid, type, 1|2, g1, g2, g3` (6 words); CORD2R/C/S
// entries `cid, type, 2, rid, a1..c3` (13 words, the coordinates in the
// file's precision) or, in 32-bit files, 22 words with the coordinates as
// doubles (NX's GEOM1N). Returns false when no layout fits.
bool op2_cord_rows(const Op2Words& rW, const std::string& rRaw, int Type, bool ByGrids,
                   std::vector<detail::NastranCoordCard>& rOut) {
    const auto ws = static_cast<std::size_t>(rW.mWs);
    const char* body = rRaw.data() + 3 * ws;
    const std::size_t nwords = (rRaw.size() - 3 * ws) / ws;
    const std::vector<std::size_t> layouts =
        ByGrids ? std::vector<std::size_t>{6} : std::vector<std::size_t>{13, 22};
    for (std::size_t size : layouts) {
        if ((size == 22 && ws != 4) || nwords == 0 || nwords % size)
            continue;
        std::vector<detail::NastranCoordCard> cards(nwords / size);
        bool ok = true;
        for (std::size_t i = 0; i < cards.size() && ok; ++i) {
            const char* e = body + i * size * ws;
            detail::NastranCoordCard& c = cards[i];
            c.mCid = rW.Int(e);
            c.mType = Type;
            c.mByGrids = ByGrids;
            ok = c.mCid > 0;
            if (ByGrids) {
                // The two flag words vary (pyNastran has seen 1 or 2 for a CORD1C).
                for (int k = 0; k < 3; ++k) {
                    c.mGrids[k] = rW.Int(e + (3 + k) * ws);
                    ok = ok && c.mGrids[k] > 0;
                }
            } else {
                ok = ok && rW.Int(e + ws) == Type;
                c.mRid = rW.Int(e + 3 * ws);
                for (std::size_t k = 0; k < 9; ++k)
                    c.mAbc[k] = size == 13 ? rW.Float(e + (4 + k) * ws) : rW.Double(e + 16 + 8 * k);
                ok = ok && c.mRid >= 0;
            }
        }
        if (!ok)
            continue;
        rOut.insert(rOut.end(), cards.begin(), cards.end());
        return true;
    }
    return false;
}

// Entries are `id, cp, x, y, z, cd, ps, seid`: 8 words, the coordinates in the
// file's precision; or, in 32-bit files, 11 words with the coordinates as
// doubles (NX writes them under the key (4501, 45, 1120001)).
bool op2_grid_rows(const Op2Words& rW, const std::string& rRaw, std::int64_t K3, Op2Grids& rOut) {
    const auto ws = static_cast<std::size_t>(rW.mWs);
    const char* body = rRaw.data() + 3 * ws;
    const std::size_t nwords = (rRaw.size() - 3 * ws) / ws;
    const std::vector<std::size_t> layouts =
        K3 != 1120001 ? std::vector<std::size_t>{8, 11} : std::vector<std::size_t>{11, 8};
    for (std::size_t size : layouts) {
        if (size == 11 && ws != 4)
            continue;
        if (nwords == 0 || nwords % size)
            continue;
        const std::size_t n = nwords / size;
        const std::size_t stride = size * ws;
        std::vector<std::int64_t> id(n), cp(n), cd(n);
        std::vector<double> xyz(3 * n);
        bool ok = true;
        for (std::size_t i = 0; i < n && ok; ++i) {
            const char* e = body + i * stride;
            id[i] = rW.Int(e);
            cp[i] = rW.Int(e + ws);
            for (std::size_t k = 0; k < 3; ++k)
                xyz[3 * i + k] = size == 8 ? rW.Float(e + (2 + k) * ws) : rW.Double(e + 8 + 8 * k);
            cd[i] = rW.Int(e + (size == 8 ? 5 * ws : 32));
            ok = id[i] > 0 && (i == 0 || id[i] > id[i - 1]) && cd[i] >= 0;
        }
        if (!ok)
            continue;
        rOut.mIds.insert(rOut.mIds.end(), id.begin(), id.end());
        rOut.mCp.insert(rOut.mCp.end(), cp.begin(), cp.end());
        rOut.mCd.insert(rOut.mCd.end(), cd.begin(), cd.end());
        rOut.mXyz.insert(rOut.mXyz.end(), xyz.begin(), xyz.end());
        return true;
    }
    return false;
}

Op2Grids op2_read_grids(const Op2Stream& rS, const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    Op2Grids g;
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"GEOM1"})) {
        if (r.mK1 == 4501 && r.mK2 == 45) {
            if (!op2_grid_rows(w, r.mRaw, r.mK3, g))
                log::warn("{}: a GRID record of {} words matches no GRID layout; skipped", kOp2Who,
                          r.mRaw.size() / static_cast<std::size_t>(w.mWs) - 3);
            continue;
        }
        // CORD1C/R/S and CORD2C/R/S, keyed (1701, 17), (1801, 18), (1901, 19),
        // (2001, 20), (2101, 21), (2201, 22).
        static const std::map<std::pair<std::int64_t, std::int64_t>, std::pair<int, bool>> kCords =
            {{{1701, 17}, {2, true}},  {{1801, 18}, {1, true}},  {{1901, 19}, {3, true}},
             {{2001, 20}, {2, false}}, {{2101, 21}, {1, false}}, {{2201, 22}, {3, false}}};
        const auto cord = kCords.find({r.mK1, r.mK2});
        if (cord != kCords.end() &&
            !op2_cord_rows(w, r.mRaw, cord->second.first, cord->second.second, g.mCords))
            log::warn(
                "{}: a coordinate system record ({}, {}) of {} words matches no layout; "
                "skipped",
                kOp2Who, r.mK1, r.mK2, r.mRaw.size() / static_cast<std::size_t>(w.mWs) - 3);
    }
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"GEOM2", "GEOM1"}))
        if ((r.mK1 == 5551 && r.mK2 == 49) || (r.mK1 == 707 && r.mK2 == 7)) {
            const auto words = w.Ints(r.mRaw);
            for (std::size_t i = 3; i < words.size(); ++i)
                g.mScalarPoints.insert(words[i]);
        }
    return g;
}

std::vector<detail::NastranCardRows> op2_read_elements(
    const Op2Stream& rS, const std::vector<Op2Table>& rTables,
    const std::unordered_set<std::int64_t>& rGrids) {
    const Op2Words& w = rS.Words();
    std::map<std::string, detail::NastranCardRows> cards;
    std::vector<std::string> skipped;
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"GEOM2"})) {
        if ((r.mK1 == 65535 && r.mK2 == 65535) || (r.mK1 == 5551 && r.mK2 == 49) ||
            (r.mK1 == 707 && r.mK2 == 7))
            continue;
        const Op2ElementSpec* spec = nullptr;
        for (const Op2ElementSpec& e : kOp2Elements)
            if (e.mK1 == r.mK1 && e.mK2 == r.mK2) {
                spec = &e;
                break;
            }
        if (spec == nullptr) {
            std::string name = "(" + std::to_string(r.mK1) + "," + std::to_string(r.mK2) + "," +
                               std::to_string(r.mK3) + ")";
            for (const auto& [k1, k2, n] : kOp2NamedElements)
                if (k1 == r.mK1 && k2 == r.mK2)
                    name = n;
            if (std::find(skipped.begin(), skipped.end(), name) == skipped.end())
                skipped.push_back(name);
            continue;
        }
        const std::vector<std::int64_t> all = w.Ints(r.mRaw);
        const std::vector<std::int64_t> words(all.begin() + 3, all.end());
        const std::size_t corners = detail::nastran_card_spec(spec->mCard)->mLinearNodes;
        const auto first = static_cast<std::size_t>(spec->mFirst);
        std::size_t chosen = 0;
        for (int k = 0; k < spec->mNumSizes && chosen == 0; ++k) {
            const auto size = static_cast<std::size_t>(spec->mSizes[k]);
            if (words.size() % size)
                continue;
            const std::size_t n = words.size() / size;
            bool ok = n > 0;
            for (std::size_t i = 0; i < n && ok; ++i) {
                const std::int64_t* row = words.data() + i * size;
                ok = row[0] > 0 && (i == 0 || row[0] > words[(i - 1) * size]);
                for (std::size_t c = 0; c < corners && ok; ++c)
                    ok = row[first + c] > 0 && rGrids.count(row[first + c]) != 0;
            }
            if (ok)
                chosen = size;
        }
        if (chosen == 0) {
            log::warn("{}: a {} record of {} words matches no known layout; skipped", kOp2Who,
                      spec->mCard, words.size());
            continue;
        }
        detail::NastranCardRows& rows = cards[spec->mCard];
        rows.mCard = spec->mCard;
        rows.mWidth = static_cast<std::size_t>(spec->mCount);
        for (std::size_t i = 0; i < words.size() / chosen; ++i) {
            const std::int64_t* row = words.data() + i * chosen;
            rows.mEid.push_back(row[0]);
            rows.mPid.push_back(spec->mPidWord >= 0 ? row[spec->mPidWord] : -1);
            rows.mNodes.insert(rows.mNodes.end(), row + first, row + first + spec->mCount);
        }
    }
    if (!skipped.empty()) {
        std::string list;
        for (const std::string& s : skipped)
            list += (list.empty() ? "" : ", ") + s;
        log::warn("{}: skipped element records with no cell type: {}", kOp2Who, list);
    }
    std::vector<detail::NastranCardRows> out;
    for (auto& [card, rows] : cards)
        out.push_back(std::move(rows));
    return out;
}

// NX also writes each PCOMP as a PSHELL whose material ids are 100000000 and
// up; those PSHELL entries are ignored, and a PCOMPG wins over a PCOMP of the
// same id, as pyNastran reads them.
std::map<std::int64_t, std::string> op2_read_properties(const Op2Stream& rS,
                                                        const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    std::map<std::int64_t, std::vector<std::string>> found;
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        const Op2PropertySpec* spec = nullptr;
        for (const Op2PropertySpec& p : kOp2Properties)
            if (p.mK1 == r.mK1 && p.mK2 == r.mK2) {
                spec = &p;
                break;
            }
        if (spec == nullptr)
            continue;
        const std::vector<std::int64_t> all = w.Ints(r.mRaw);
        const std::vector<std::int64_t> words(all.begin() + 3, all.end());
        const std::size_t n = words.size();
        std::vector<std::int64_t> pids;
        switch (spec->mLayout) {
            case Op2PropLayout::Fixed: {
                const auto size = static_cast<std::size_t>(spec->mSizes[0]);
                if (n % size == 0)
                    for (std::size_t k = 0; k < n; k += size)
                        if (std::string(spec->mCard) != "PSHELL" || words[k + 1] < 100000000)
                            pids.push_back(words[k]);
                break;
            }
            case Op2PropLayout::Sizes:
                for (int s = 0; s < spec->mNumSizes; ++s) {
                    const auto size = static_cast<std::size_t>(spec->mSizes[s]);
                    if (n % size)
                        continue;
                    bool ok = true;
                    for (std::size_t k = 0; k < n && ok; k += size)
                        ok = words[k] > 0 && (k == 0 || words[k] > words[k - size]);
                    if (ok) {
                        for (std::size_t k = 0; k < n; k += size)
                            pids.push_back(words[k]);
                        break;
                    }
                }
                break;
            case Op2PropLayout::Pcomp: {
                std::size_t k = 0;
                while (k + 8 <= n) {
                    const std::size_t nlayers = static_cast<std::size_t>(std::abs(words[k + 1]));
                    if (words[k] <= 0 || nlayers == 0 || k + 8 + 4 * nlayers > n)
                        break;
                    pids.push_back(words[k]);
                    k += 8 + 4 * nlayers;
                }
                break;
            }
            case Op2PropLayout::Pcompg: {
                std::size_t k = 0;
                while (k + 8 <= n) {
                    if (words[k] <= 0)
                        break;
                    const std::int64_t pid = words[k];
                    k += 8;
                    const auto terminator = [&](std::size_t At) {
                        for (std::size_t j = 0; j < 5; ++j)
                            if (words[At + j] != -1)
                                return false;
                        return true;
                    };
                    while (k + 5 <= n && !terminator(k))
                        k += 5;
                    if (k + 5 > n)
                        break;
                    pids.push_back(pid);
                    k += 5;
                }
                break;
            }
            case Op2PropLayout::Terminated: {
                std::size_t k = 0;
                while (k + 6 <= n) {
                    if (words[k] <= 0 || (!pids.empty() && words[k] <= pids.back()))
                        break;
                    std::size_t end = k + 6;
                    while (end < n && words[end] != -1)
                        ++end;
                    if (end >= n)
                        break;
                    pids.push_back(words[k]);
                    k = end + 1;
                }
                break;
            }
        }
        for (std::int64_t p : pids)
            found[p].push_back(spec->mCard);
    }
    std::map<std::int64_t, std::string> ptype;
    for (const auto& [pid, list] : found) {
        std::string chosen = list.front();
        for (const char* priority : {"PCOMPG", "PCOMP"})
            if (std::find(list.begin(), list.end(), priority) != list.end()) {
                chosen = priority;
                break;
            }
        ptype[pid] = chosen;
    }
    return ptype;
}

// --- results -------------------------------------------------------------------------

bool op2_starts(const std::string& rName, std::initializer_list<const char*> Prefixes) {
    for (const char* p : Prefixes)
        if (rName.compare(0, std::strlen(p), p) == 0)
            return true;
    return false;
}

const char* op2_nodal_name(std::int64_t TableCode) {
    switch (TableCode) {
        case 1:
            return "DISPLACEMENT";
        case 7:
            return "EIGENVECTOR";
        case 10:
            return "VELOCITY";
        case 11:
            return "ACCELERATION";
        case 3:
            return "SPC_FORCE";
        case 39:
            return "MPC_FORCE";
        case 2:
            return "APPLIED_LOAD";
        default:
            return nullptr;
    }
}

std::string op2_element_type_name(std::int64_t Type) {
    static const std::map<std::int64_t, const char*> names = {{1, "CROD"},
                                                              {2, "CBEAM"},
                                                              {3, "CTUBE"},
                                                              {4, "CSHEAR"},
                                                              {10, "CONROD"},
                                                              {11, "CELAS1"},
                                                              {12, "CELAS2"},
                                                              {13, "CELAS3"},
                                                              {14, "CELAS4"},
                                                              {33, "CQUAD4"},
                                                              {34, "CBAR"},
                                                              {39, "CTETRA"},
                                                              {64, "CQUAD8"},
                                                              {67, "CHEXA"},
                                                              {68, "CPENTA"},
                                                              {70, "CTRIAR"},
                                                              {74, "CTRIA3"},
                                                              {75, "CTRIA6"},
                                                              {82, "CQUADR"},
                                                              {95, "CQUAD4 composite"},
                                                              {96, "CQUAD8 composite"},
                                                              {97, "CTRIA3 composite"},
                                                              {98, "CTRIA6 composite"},
                                                              {100, "CBAR stations"},
                                                              {102, "CBUSH"},
                                                              {144, "CQUAD4 corner"},
                                                              {255, "CPYRAM"}};
    const auto it = names.find(Type);
    return it != names.end() ? it->second : "type " + std::to_string(Type);
}

using Op2Layout = std::vector<std::pair<std::size_t, std::string>>;  // (word, member)

// How an element table's entries hold their values.
enum class Op2Values {
    Row,      // one entry per element, members at fixed words
    Ply,      // one entry per ply of a composite element: word 1 is the ply
    Station,  // one entry per station along a bar (CBAR type 100), in order
    Blocks,   // one entry per element made of node blocks, each led by its GRID:
              // the centre then the corners, or a beam's stations
};

struct Op2ElementLayout {
    Op2Values mKind = Op2Values::Row;
    Op2Layout mMembers;       // Blocks: the word within a block (its GRID is word 0)
    std::size_t mFirst = 0;   // Blocks: word of the first block
    std::size_t mBlock = 0;   // Blocks: words per block
    std::size_t mBlocks = 1;  // Blocks: blocks per element
    bool mStations = false;   // Blocks: beam stations rather than centre + corners
};

std::optional<Op2ElementLayout> op2_element_layout(std::int64_t Type, std::int64_t NumWide,
                                                   std::int64_t SCode) {
    static const char* plate[16] = {"FD1",    "X1",     "Y1",     "TXY1", "ANGLE1", "MAJOR1",
                                    "MINOR1", nullptr,  "FD2",    "X2",   "Y2",     "TXY2",
                                    "ANGLE2", "MAJOR2", "MINOR2", nullptr};
    static const char* bar[15] = {"X1A", "X2A", "X3A", "X4A", "AX",   "MAXA", "MINA", "MST",
                                  "X1B", "X2B", "X3B", "X4B", "MAXB", "MINB", "MSC"};
    static const std::pair<std::size_t, const char*> solid[] = {
        {1, "X"},     {2, "TXY"},  {3, "PRINCIPAL_A"}, {7, "PRESSURE"},
        {8, nullptr}, {9, "Y"},    {10, "TYZ"},        {11, "PRINCIPAL_B"},
        {15, "Z"},    {16, "TZX"}, {17, "PRINCIPAL_C"}};
    const std::string vm = (SCode & 1) ? "VON_MISES" : "MAX_SHEAR";
    const auto corner_nodes = [](std::int64_t T) -> std::int64_t {
        switch (T) {
            case 64:
            case 82:
            case 144:
                return 5;
            case 70:
            case 75:
                return 4;
            default:
                return 0;
        }
    };
    const auto solid_nodes = [](std::int64_t T) -> std::int64_t {
        switch (T) {
            case 39:
                return 5;
            case 67:
                return 9;
            case 68:
                return 7;
            case 255:
                return 6;
            default:
                return 0;
        }
    };
    Op2ElementLayout out;
    auto rows = [&](Op2Layout rMembers) {
        out.mMembers = std::move(rMembers);
        return out;
    };
    if ((Type == 1 || Type == 10) && NumWide == 5)
        return rows({{1, "A"}, {2, "MSA"}, {3, "T"}, {4, "MST"}});
    if (Type == 3 && NumWide == 5)
        return rows({{1, "AS"}, {2, "MSA"}, {3, "TS"}, {4, "MST"}});
    if (Type == 4 && NumWide == 4)
        return rows({{1, "TMAX"}, {2, "TAVG"}, {3, "MS"}});
    if (Type == 34 && NumWide == 16) {
        for (std::size_t k = 0; k < 15; ++k)
            out.mMembers.emplace_back(1 + k, bar[k]);
        return out;
    }
    // CBAR stations (100): sd, the four fibres, axial, max, min, margin.
    if (Type == 100 && NumWide == 10) {
        out.mKind = Op2Values::Station;
        return rows({{1, "SD"},
                     {2, "XC"},
                     {3, "XD"},
                     {4, "XE"},
                     {5, "XF"},
                     {6, "AX"},
                     {7, "MAX"},
                     {8, "MIN"},
                     {9, "MS"}});
    }
    // Composite shells (QUAD4, QUAD8, TRIA3, TRIA6; NX QUADR, TRIAR): one entry
    // per ply, the MSC HDF5 names.
    if ((Type == 95 || Type == 96 || Type == 97 || Type == 98 || Type == 232 || Type == 233) &&
        NumWide == 11) {
        out.mKind = Op2Values::Ply;
        return rows({{2, "X1"},
                     {3, "Y1"},
                     {4, "T1"},
                     {5, "L1"},
                     {6, "L2"},
                     {7, "ANGLE"},
                     {8, "MAJOR"},
                     {9, "MINOR"},
                     {10, vm}});
    }
    // CBEAM (2): 11 stations of grid, sd, the four fibres, max, min and margins.
    if (Type == 2 && NumWide == 111) {
        out.mKind = Op2Values::Blocks;
        out.mFirst = 1;
        out.mBlock = 10;
        out.mBlocks = 11;
        out.mStations = true;
        return rows({{1, "SD"},
                     {2, "XC"},
                     {3, "XD"},
                     {4, "XE"},
                     {5, "XF"},
                     {6, "MAX"},
                     {7, "MIN"},
                     {8, "MST"},
                     {9, "MSC"}});
    }
    if ((Type == 33 || Type == 74) && NumWide == 17) {
        for (std::size_t k = 0; k < 16; ++k)
            out.mMembers.emplace_back(
                1 + k, plate[k] ? std::string(plate[k]) : vm + std::to_string(1 + k / 8));
        return out;
    }
    const std::int64_t cn = corner_nodes(Type);
    if (cn && NumWide == 2 + 17 * cn) {
        out.mKind = Op2Values::Blocks;
        out.mFirst = 2;
        out.mBlock = 17;
        out.mBlocks = static_cast<std::size_t>(cn);
        for (std::size_t k = 0; k < 16; ++k)
            out.mMembers.emplace_back(
                1 + k, plate[k] ? std::string(plate[k]) : vm + std::to_string(1 + k / 8));
        return out;
    }
    const std::int64_t sn = solid_nodes(Type);
    if (sn && NumWide == 4 + 21 * sn) {
        out.mKind = Op2Values::Blocks;
        out.mFirst = 4;
        out.mBlock = 21;
        out.mBlocks = static_cast<std::size_t>(sn);
        const std::string octa = (SCode & 1) ? "VON_MISES" : "OCT_SHEAR";
        for (const auto& [word, name] : solid)
            out.mMembers.emplace_back(word, name ? std::string(name) : octa);
        return out;
    }
    return std::nullopt;
}

double op2_time_of(std::int64_t Analysis, double W5, double W6) {
    switch (Analysis) {
        case 5:
        case 6:
        case 10:
        case 12:
            return W5;
        case 2:
        case 8:
        case 9:
            return W6;
        default:
            return 0.0;
    }
}

struct Op2Step {
    std::int64_t mSubcase, mAnalysis, mMode;
    double mTime;
};

struct Op2Block {
    std::size_t mStep;
    bool mNodal;
    bool mBasic = false;  // BOUG*: already in the basic system
    std::string mBase;  // nodal: point data name; element: STRESS/STRAIN
    Op2ElementLayout mLayout;
    std::size_t mNumWide = 8;
    bool mGridForce = false;  // OGPFB: grid point forces
    const Op2Record* mRecord;
};

}  // namespace

namespace {

const detail::FileSource& op2_nonempty(const detail::FileSource& rSource,
                                       const std::string& rPath) {
    if (rSource.Size() == 0)
        op2_fail("'" + rPath + "' is empty");
    return rSource;
}

// The whole file: its tables, and the result blocks of the supported tables.
class Op2Reader {
public:
    explicit Op2Reader(const std::string& rPath)
        : mPath(rPath),
          mSource(rPath),
          mStream(op2_nonempty(mSource, rPath).Data(), mSource.Size()) {
        mTables = op2_parse(mStream);
        ScanResults();
    }

    const std::string mPath;
    detail::FileSource mSource;
    Op2Stream mStream;
    std::vector<Op2Table> mTables;
    std::vector<Op2Step> mSteps;
    std::vector<Op2Block> mBlocks;
    std::vector<std::string> mSkipped;

private:
    void Skip(const std::string& rReason) {
        if (std::find(mSkipped.begin(), mSkipped.end(), rReason) == mSkipped.end())
            mSkipped.push_back(rReason);
    }

    void ScanResults() {
        const Op2Words& w = mStream.Words();
        const auto ws = static_cast<std::size_t>(w.mWs);
        std::map<std::tuple<std::int64_t, std::int64_t, std::int64_t>, std::size_t> step_index;
        struct Deferred {
            std::tuple<std::int64_t, std::int64_t, std::int64_t> mKey;
            Op2Block mBlock;
            std::string mTable;
        };
        std::vector<Deferred> deferred;
        for (const Op2Table& t : mTables) {
            const bool nodal = op2_starts(t.mName, {"OUG", "BOUG", "OQG", "OQMG", "OPG"});
            const bool elemental = op2_starts(t.mName, {"OES", "OSTR"});
            const bool grid_force = op2_starts(t.mName, {"OGPF"});
            if (!nodal && !elemental && !grid_force) {
                if (!t.mName.empty() && t.mName[0] == 'O')
                    Skip(t.mName);
                continue;
            }
            std::string header;
            for (const auto& [mark, rec] : t.mRecords) {
                if (mark > -3)
                    continue;
                if (rec.Size() == 146 * ws) {
                    header = rec.Join(mStream.Base());
                    continue;
                }
                if (header.empty())
                    continue;
                const char* h = header.data();
                const auto word = [&](std::size_t I) { return w.Int(h + I * ws); };
                const std::int64_t approach = word(0), tcode = word(1), etype = word(2),
                                   subcase = word(3);
                const std::int64_t device = approach % 10;
                const std::int64_t analysis = (approach - device) / 10;
                const std::int64_t table_code = tcode % 1000;
                const std::int64_t sort_code = tcode / 1000;
                const std::int64_t num_wide = word(9);
                const std::int64_t s_code = word(10);
                const std::int64_t thermal = word(22);
                if (sort_code & 1) {
                    Skip(t.mName + " (complex)");
                    continue;
                }
                if (sort_code & 4) {
                    Skip(t.mName + " (random)");
                    continue;
                }
                if (sort_code & 2) {
                    Skip(t.mName + " (SORT2)");
                    continue;
                }
                Op2Block block{0, nodal, op2_starts(t.mName, {"BOUG"}), {}, {}, 8, false, &rec};
                if (grid_force) {
                    if (table_code != 19 || num_wide != 10) {
                        Skip(t.mName + " (table code " + std::to_string(table_code) + ", " +
                             std::to_string(num_wide) + " words)");
                        continue;
                    }
                    block.mGridForce = true;
                    block.mBase = "GRID_FORCE";
                    block.mNumWide = 10;
                } else if (nodal) {
                    const char* base = op2_nodal_name(table_code);
                    // MPC forces share the SPC forces' table code; the name tells.
                    if (op2_starts(t.mName, {"OQMG"}) && (table_code == 3 || table_code == 39))
                        base = "MPC_FORCE";
                    if (base == nullptr) {
                        Skip(t.mName + " (table code " + std::to_string(table_code) + ")");
                        continue;
                    }
                    if (num_wide != 8) {
                        Skip(t.mName + " (" + std::to_string(num_wide) + " words per node)");
                        continue;
                    }
                    if (thermal == 1) {
                        if (table_code != 1) {
                            Skip(t.mName + " (thermal table code " + std::to_string(table_code) +
                                 ")");
                            continue;
                        }
                        base = "TEMPERATURE";
                    }
                    block.mBase = base;
                } else {
                    if (table_code != 5) {
                        Skip(t.mName + " (table code " + std::to_string(table_code) + ")");
                        continue;
                    }
                    auto layout = op2_element_layout(etype, num_wide, s_code);
                    if (!layout) {
                        Skip(t.mName + " " + op2_element_type_name(etype));
                        continue;
                    }
                    block.mBase = (s_code & 8) ? "STRAIN" : "STRESS";
                    block.mLayout = std::move(*layout);
                    block.mNumWide = static_cast<std::size_t>(num_wide);
                }
                const std::int64_t w5 = word(4);
                const auto key = std::make_tuple(subcase, analysis, w5);
                auto it = step_index.find(key);
                if (grid_force) {
                    deferred.push_back({key, std::move(block), t.mName});
                    continue;
                }
                if (it == step_index.end()) {
                    it = step_index.emplace(key, mSteps.size()).first;
                    const bool moded = analysis == 2 || analysis == 8 || analysis == 9;
                    mSteps.push_back(
                        {subcase, analysis, moded ? w5 : 0,
                         op2_time_of(analysis, w.Float(h + 4 * ws), w.Float(h + 5 * ws))});
                }
                block.mStep = it->second;
                mBlocks.push_back(std::move(block));
            }
        }
        // Grid point forces join a step the other tables made, never a new one:
        // MSC writes 0 in their word 5 where the other tables of a static step
        // write the load set, and a buckling run's forces carry analysis 2.
        for (Deferred& d : deferred) {
            auto it = step_index.find(d.mKey);
            if (it == step_index.end()) {
                std::size_t same = 0;
                for (auto s = step_index.begin(); s != step_index.end(); ++s)
                    if (std::get<0>(s->first) == std::get<0>(d.mKey) &&
                        std::get<1>(s->first) == std::get<1>(d.mKey)) {
                        it = s;
                        ++same;
                    }
                if (same != 1) {
                    Skip(d.mTable + " (no matching step)");
                    continue;
                }
            }
            d.mBlock.mStep = it->second;
            mBlocks.push_back(std::move(d.mBlock));
        }
    }
};

NDArray op2_int_array(const std::vector<std::int64_t>& rV) {
    NDArray a(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

NDArray op2_scalar(double V, DType T) {
    NDArray a(T, {1});
    if (T == DType::Float64)
        a.As<double>()[0] = V;
    else
        a.As<std::int64_t>()[0] = static_cast<std::int64_t>(V);
    return a;
}

// The input deck beside an OP2, and the paths tried.
std::pair<std::string, std::vector<std::string>> op2_sibling_deck(const std::string& rPath) {
    // os.path.splitext: the last dot of the file name, leading dots excluded.
    const std::size_t slash = rPath.find_last_of("/\\");
    const std::size_t name = slash == std::string::npos ? 0 : slash + 1;
    std::size_t lead = name;
    while (lead < rPath.size() && rPath[lead] == '.')
        ++lead;
    const std::size_t dot = rPath.find_last_of('.');
    const std::string stem =
        (dot != std::string::npos && dot > lead) ? rPath.substr(0, dot) : rPath;
    std::vector<std::string> tried;
    for (const char* ext : {".bdf", ".dat", ".nas", ".blk", ".BDF", ".DAT", ".NAS", ".BLK"}) {
        tried.push_back(stem + ext);
        std::error_code ec;
        if (fs::is_regular_file(tried.back(), ec))
            return {tried.back(), tried};
    }
    return {std::string(), tried};
}

struct Op2Model {
    Mesh mMesh;
    detail::NastranCoordSystems mSystems;
    std::vector<std::int64_t> mCd;  // per point; empty for a deck-built model
    std::unordered_map<std::int64_t, std::size_t> mGridIndex;
    std::unordered_map<std::int64_t, std::size_t> mCellIndex;
    std::vector<std::size_t> mOffsets, mSizes;
};

Op2Model op2_build_mesh(const Op2Reader& rR) {
    Op2Model m;
    Op2Grids g = op2_read_grids(rR.mStream, rR.mTables);
    if (g.mIds.empty()) {
        auto [deck, tried] = op2_sibling_deck(rR.mPath);
        if (deck.empty()) {
            std::string list;
            for (const std::string& t : tried)
                list += (list.empty() ? "" : ", ") + t;
            op2_fail(
                "the file has no GEOM1 GRID records (rerun with PARAM,POST,-1 or provide the "
                "input deck beside it); looked for " +
                list);
        }
        std::vector<std::int64_t> grid_ids, cell_ids;
        m.mMesh = detail::nastran_read_deck(deck, grid_ids, cell_ids);
        for (std::size_t i = 0; i < grid_ids.size(); ++i)
            m.mGridIndex.emplace(grid_ids[i], i);
        std::vector<NDArray> eids;
        std::size_t base = 0;
        for (const auto cb : m.mMesh.CellRange()) {
            const std::size_t n = cb.NumCells();
            m.mOffsets.push_back(base);
            m.mSizes.push_back(n);
            std::vector<std::int64_t> ids(cell_ids.begin() + static_cast<std::ptrdiff_t>(base),
                                          cell_ids.begin() + static_cast<std::ptrdiff_t>(base + n));
            for (std::size_t i = 0; i < n; ++i)
                m.mCellIndex.emplace(ids[i], base + i);
            eids.push_back(op2_int_array(ids));
            base += n;
        }
        if (!eids.empty())
            m.mMesh.AddCellData("nastran:eid", std::move(eids));
        return m;
    }
    for (std::size_t i = 0; i < g.mIds.size(); ++i)
        if (!m.mGridIndex.emplace(g.mIds[i], i).second)
            op2_fail("GRID " + std::to_string(g.mIds[i]) + " is defined twice");
    NDArray points(DType::Float64, {g.mIds.size(), 3});
    std::copy(g.mXyz.begin(), g.mXyz.end(), points.As<double>());
    m.mMesh.AssignPoints(std::move(points));
    m.mSystems = detail::nastran_apply_frames(m.mMesh, g.mCords, g.mIds, g.mCp, g.mCd, kOp2Who);
    m.mCd = g.mCd;
    const std::unordered_set<std::int64_t> grids(g.mIds.begin(), g.mIds.end());
    const auto cards = op2_read_elements(rR.mStream, rR.mTables, grids);
    const detail::NastranCells cells =
        detail::nastran_add_cells(m.mMesh, cards, m.mGridIndex, g.mScalarPoints,
                                  op2_read_properties(rR.mStream, rR.mTables), kOp2Who);
    m.mCellIndex = cells.mCellIndex;
    m.mOffsets = cells.mOffsets;
    m.mSizes = cells.mSizes;
    return m;
}

void op2_warn_skipped(const Op2Reader& rR) {
    if (rR.mSkipped.empty())
        return;
    std::string list;
    for (const std::string& s : rR.mSkipped)
        list += (list.empty() ? "" : ", ") + s;
    log::warn("{}: result tables not read: {}", kOp2Who, list);
}

}  // namespace

Mesh read_nastran_op2(const std::string& rPath, const ReadOptions& rOpts) {
    const Op2Reader r(rPath);
    Op2Model model = op2_build_mesh(r);
    Mesh& mesh = model.mMesh;
    const std::size_t n = r.mSteps.size();
    if (n == 0) {
        if (rOpts.mTimeStep != 0 && rOpts.mTimeStep != -1)
            throw ReadError("time step " + std::to_string(rOpts.mTimeStep) +
                            " is out of range: the file has no steps");
        op2_warn_skipped(r);
        return mesh;
    }
    const std::size_t index = rOpts.ResolveTimeStep(n);
    const Op2Step& step = r.mSteps[index];
    mesh.AddFieldData(kSequenceTimeKey, op2_scalar(step.mTime, DType::Float64));
    mesh.AddFieldData("nastran:subcase",
                      op2_scalar(static_cast<double>(step.mSubcase), DType::Int64));
    mesh.AddFieldData("nastran:analysis",
                      op2_scalar(static_cast<double>(step.mAnalysis), DType::Int64));
    mesh.AddFieldData("nastran:mode", op2_scalar(static_cast<double>(step.mMode), DType::Int64));
    op2_warn_skipped(r);
    if (rOpts.mPointsOnly)
        return mesh;

    const Op2Words& w = r.mStream.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    const std::size_t npts = mesh.NumPoints();
    std::size_t ncells = 0;
    for (std::size_t s : model.mSizes)
        ncells += s;
    // name -> (components, values)
    std::map<std::string, std::pair<std::size_t, std::vector<double>>> point_arrays;
    std::map<std::string, std::vector<double>> cell_arrays;
    const NDArray basic_points = mesh.Points();
    const double* basic = basic_points.As<double>();
    // Multi-valued cell arrays (per ply, station, corner or element node):
    // (cell, column, value), laid out once every block is read; the first value
    // of a cell and column wins, as for the centre values.
    struct Wide {
        std::vector<std::size_t> mCell, mCol;
        std::vector<double> mValue;
    };
    std::map<std::string, Wide> wide;
    auto push = [&](const std::string& rName, std::size_t Cell, std::size_t Col, double V) {
        Wide& e = wide[rName];
        e.mCell.push_back(Cell);
        e.mCol.push_back(Col);
        e.mValue.push_back(V);
    };
    // The position of point `Point` in the connectivity of global cell `Cell`.
    auto node_position = [&](std::size_t Cell, std::size_t Point) -> std::size_t {
        const std::size_t b = static_cast<std::size_t>(
            std::upper_bound(model.mOffsets.begin(), model.mOffsets.end(), Cell) -
            model.mOffsets.begin() - 1);
        const NDArray& conn = mesh.Cells(b).Conn();
        const std::size_t width = conn.Shape()[1];
        const std::size_t row = Cell - model.mOffsets[b];
        for (std::size_t k = 0; k < width; ++k)
            if (static_cast<std::size_t>(detail::read_int(conn, row * width + k)) == Point)
                return k;
        return kOp2Npos;
    };
    for (const Op2Block& b : r.mBlocks) {
        if (b.mStep != index)
            continue;
        const std::string raw = b.mRecord->Join(r.mStream.Base());
        const std::size_t nwords = raw.size() / ws;
        if (b.mNodal) {
            if (nwords % 8)
                op2_fail("a " + b.mBase + " record is not a whole number of rows");
            std::vector<std::tuple<std::string, std::size_t, std::size_t>> outputs;
            if (b.mBase == "TEMPERATURE")
                outputs.emplace_back(b.mBase, 2, 1);
            else {
                outputs.emplace_back(b.mBase, 2, 3);
                outputs.emplace_back(b.mBase + "_ROT", 5, 3);
            }
            for (const auto& [name, c0, nc] : outputs) {
                if (!rOpts.WantsArray(name))
                    continue;
                auto it = point_arrays.find(name);
                if (it == point_arrays.end())
                    it = point_arrays
                             .emplace(name,
                                      std::make_pair(nc, std::vector<double>(npts * nc, kOp2Nan)))
                             .first;
                std::vector<double>& values = it->second.second;
                for (std::size_t r0 = 0; r0 < nwords / 8; ++r0) {
                    const char* row = raw.data() + r0 * 8 * ws;
                    const auto p = model.mGridIndex.find(w.Int(row) / 10);
                    if (p == model.mGridIndex.end() || !std::isnan(values[p->second * nc]))
                        continue;
                    double* dst = values.data() + p->second * nc;
                    for (std::size_t c = 0; c < nc; ++c)
                        dst[c] = w.Float(row + (c0 + c) * ws);
                    // Results are in the GRID's output system (CD) unless the table is BOUG*.
                    if (nc == 3 && !b.mBasic && !model.mCd.empty())
                        detail::nastran_rotate_to_basic(model.mSystems, model.mCd[p->second],
                                                        basic + 3 * p->second, dst, 1);
                }
            }
        } else if (b.mGridForce) {
            // Entries: grid*10+device, element (0 for the totals, applied loads,
            // SPC forces, ...), its 8-character name, F1 F2 F3 M1 M2 M3 in the
            // GRID's output system (CD). Element rows become (cells, nodes) in
            // the cell's node order; the others point data named after the label.
            if (nwords % 10)
                op2_fail("a " + b.mBase + " record is not a whole number of entries");
            static const char* members[6] = {"F1", "F2", "F3", "M1", "M2", "M3"};
            for (std::size_t r0 = 0; r0 < nwords / 10; ++r0) {
                const char* row = raw.data() + r0 * 10 * ws;
                const auto p = model.mGridIndex.find(w.Int(row) / 10);
                if (p == model.mGridIndex.end())
                    continue;
                double v[6];
                for (std::size_t c = 0; c < 6; ++c)
                    v[c] = w.Float(row + (4 + c) * ws);
                if (!model.mCd.empty())
                    detail::nastran_rotate_to_basic(model.mSystems, model.mCd[p->second],
                                                    basic + 3 * p->second, v, 2);
                const std::int64_t eid = w.Int(row + ws);
                if (eid > 0) {
                    const auto c = model.mCellIndex.find(eid);
                    if (c == model.mCellIndex.end())
                        continue;
                    const std::size_t pos = node_position(c->second, p->second);
                    if (pos == kOp2Npos)
                        continue;
                    for (std::size_t k = 0; k < 6; ++k) {
                        const std::string name = b.mBase + ":" + members[k];
                        if (rOpts.WantsArray(name))
                            push(name, c->second, pos, v[k]);
                    }
                    continue;
                }
                std::string label;
                for (std::size_t k = 0; k < 2 * ws; ++k) {
                    const char ch = row[2 * ws + k];
                    if (ch != ' ' && ch != '*' && ch != '\0')
                        label += ch;
                }
                for (std::size_t k = 0; k < 6; ++k) {
                    const std::string name = b.mBase + ":" + label + ":" + members[k];
                    if (!rOpts.WantsArray(name))
                        continue;
                    auto it = point_arrays.find(name);
                    if (it == point_arrays.end())
                        it = point_arrays
                                 .emplace(name, std::make_pair(std::size_t{1},
                                                               std::vector<double>(npts, kOp2Nan)))
                                 .first;
                    if (std::isnan(it->second.second[p->second]))
                        it->second.second[p->second] = v[k];
                }
            }
        } else {
            const std::size_t nw = b.mNumWide;
            if (nwords % nw)
                op2_fail("a " + b.mBase + " record is not a whole number of elements");
            const Op2ElementLayout& lay = b.mLayout;
            std::unordered_map<std::int64_t, std::size_t> stations;  // CBAR: stations seen
            for (std::size_t r0 = 0; r0 < nwords / nw; ++r0) {
                const char* row = raw.data() + r0 * nw * ws;
                const std::int64_t eid = w.Int(row) / 10;
                const auto c = model.mCellIndex.find(eid);
                const std::size_t station = lay.mKind == Op2Values::Station ? stations[eid]++ : 0;
                if (c == model.mCellIndex.end())
                    continue;
                const std::size_t cell = c->second;
                if (lay.mKind == Op2Values::Ply || lay.mKind == Op2Values::Station) {
                    const bool ply = lay.mKind == Op2Values::Ply;
                    const std::int64_t ply_id = ply ? w.Int(row + ws) : 1;
                    if (ply_id < 1)
                        continue;
                    const std::size_t col = ply ? static_cast<std::size_t>(ply_id - 1) : station;
                    const char* suffix = ply ? "@ply" : "@station";
                    for (const auto& [word, member] : lay.mMembers) {
                        const std::string name = b.mBase + ":" + member + suffix;
                        if (rOpts.WantsArray(name))
                            push(name, cell, col, w.Float(row + word * ws));
                    }
                    continue;
                }
                // Row: the values; Blocks: the first block's (the centre, or end A).
                const std::size_t base = lay.mKind == Op2Values::Blocks ? lay.mFirst : 0;
                for (const auto& [word, member] : lay.mMembers) {
                    const std::string name = b.mBase + ":" + member;
                    if (!rOpts.WantsArray(name))
                        continue;
                    auto it = cell_arrays.find(name);
                    if (it == cell_arrays.end())
                        it = cell_arrays.emplace(name, std::vector<double>(ncells, kOp2Nan)).first;
                    if (std::isnan(it->second[cell]))
                        it->second[cell] = w.Float(row + (base + word) * ws);
                }
                if (lay.mKind != Op2Values::Blocks)
                    continue;
                for (std::size_t k = lay.mStations ? 0 : 1; k < lay.mBlocks; ++k) {
                    const char* block = row + (lay.mFirst + k * lay.mBlock) * ws;
                    std::size_t col = k;
                    // A beam station with no GRID and no distance was not output.
                    if (lay.mStations && k > 0 && w.Int(block) == 0 && w.Float(block + ws) == 0.0)
                        continue;
                    if (!lay.mStations) {
                        const auto g = model.mGridIndex.find(w.Int(block));
                        if (g == model.mGridIndex.end())
                            continue;
                        col = node_position(cell, g->second);
                        if (col == kOp2Npos)
                            continue;
                    }
                    const char* suffix = lay.mStations ? "@station" : "@corner";
                    for (const auto& [word, member] : lay.mMembers) {
                        const std::string name = b.mBase + ":" + member + suffix;
                        if (rOpts.WantsArray(name))
                            push(name, cell, col, w.Float(block + word * ws));
                    }
                }
            }
        }
    }
    for (auto& [name, entry] : point_arrays) {
        const auto& [nc, values] = entry;
        NDArray a = nc == 1 ? NDArray(DType::Float64, {npts}) : NDArray(DType::Float64, {npts, nc});
        std::copy(values.begin(), values.end(), a.As<double>());
        mesh.AddPointData(name, std::move(a));
    }
    for (auto& [name, values] : cell_arrays) {
        std::vector<NDArray> per_block;
        for (std::size_t b = 0; b < model.mSizes.size(); ++b) {
            NDArray a(DType::Float64, {model.mSizes[b]});
            std::copy(
                values.begin() + static_cast<std::ptrdiff_t>(model.mOffsets[b]),
                values.begin() + static_cast<std::ptrdiff_t>(model.mOffsets[b] + model.mSizes[b]),
                a.As<double>());
            per_block.push_back(std::move(a));
        }
        mesh.AddCellData(name, std::move(per_block));
    }
    for (const auto& [name, e] : wide) {
        std::size_t width = 0;
        for (std::size_t col : e.mCol)
            width = std::max(width, col + 1);
        std::vector<NDArray> per_block;
        for (std::size_t b = 0; b < model.mSizes.size(); ++b) {
            NDArray a(DType::Float64, {model.mSizes[b], width});
            std::fill(a.As<double>(), a.As<double>() + a.Size(), kOp2Nan);
            per_block.push_back(std::move(a));
        }
        for (std::size_t i = 0; i < e.mCell.size(); ++i) {
            const std::size_t b = static_cast<std::size_t>(
                std::upper_bound(model.mOffsets.begin(), model.mOffsets.end(), e.mCell[i]) -
                model.mOffsets.begin() - 1);
            double& slot =
                per_block[b].As<double>()[(e.mCell[i] - model.mOffsets[b]) * width + e.mCol[i]];
            if (std::isnan(slot))
                slot = e.mValue[i];
        }
        mesh.AddCellData(name, std::move(per_block));
        NDArray layout(DType::Int64, {std::size_t{2}});
        layout.As<std::int64_t>()[0] = static_cast<std::int64_t>(width);
        layout.As<std::int64_t>()[1] = 1;
        mesh.AddFieldData("nastran:layout:" + name, std::move(layout));
    }
    return mesh;
}

std::vector<double> nastran_op2_time_values(const std::string& rPath) {
    const Op2Reader r(rPath);
    std::vector<double> out;
    for (const Op2Step& s : r.mSteps)
        out.push_back(s.mTime);
    return out;
}

MeshMetadata read_nastran_op2_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    ReadOptions options = rOpts;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(read_nastran_op2(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "nastran_op2";
    meta.mTimeValues = nastran_op2_time_values(rPath);
    return meta;
}

}  // namespace meshioplusplus
