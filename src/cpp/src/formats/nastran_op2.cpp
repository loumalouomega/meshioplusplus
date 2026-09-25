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
    {201, 2, "CDAMP1", {6}, 1, 2, 2, 1},
    {301, 3, "CDAMP2", {6}, 1, 2, 2, -1},
    {601, 6, "CELAS1", {6}, 1, 2, 2, 1},
    {701, 7, "CELAS2", {8}, 1, 2, 2, -1},
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
    {801, 8, "CELAS3"},   {901, 9, "CELAS4"},     {401, 4, "CDAMP3"},    {501, 5, "CDAMP4"},
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
            // A fluid (acoustic) GRID's CD is -1: it has no output system.
            ok = id[i] > 0 && (i == 0 || id[i] > id[i - 1]) && cd[i] >= -1;
            cd[i] = std::max<std::int64_t>(cd[i], 0);
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

// The first record of mark `Mark` of the first table whose name starts with
// `Prefix`, or nullptr.
const Op2Record* op2_table_record(const std::vector<Op2Table>& rTables, const char* pPrefix,
                                  std::int64_t Mark) {
    for (const Op2Table& t : rTables) {
        if (t.mName.compare(0, std::strlen(pPrefix), pPrefix) != 0)
            continue;
        for (const auto& [mark, rec] : t.mRecords)
            if (mark == Mark)
                return &rec;
        return nullptr;
    }
    return nullptr;
}

// The grid points of a file without GEOM1, from its basic grid point table
// (BGPDT/BGPDTS): basic coordinates and output system per point. BGPDTS rows
// are `cd, x, y, z` in internal order, named by EQEXIN (pairs `id, internal
// sequence` and `id, 10 * sil + type`: type 2 is a scalar point); NX's BGPDT
// rows carry the id: `cd, sil, id, 61, ps, 0, x, y, z` with the coordinates
// as doubles. Returns false when there is none or it matches no layout.
bool op2_bgpdt_grids(const Op2Stream& rS, const std::vector<Op2Table>& rTables, Op2Grids& rOut) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    const Op2Record* header = op2_table_record(rTables, "BGPDT", -1);
    const Op2Record* data = op2_table_record(rTables, "BGPDT", -3);
    if (header == nullptr || data == nullptr || header->Size() < 2 * ws)
        return false;
    const std::string head = header->Join(rS.Base());
    const std::string raw = data->Join(rS.Base());
    const std::int64_t n = w.Int(head.data() + ws);
    if (n <= 0 || raw.size() % static_cast<std::size_t>(n))
        return false;
    const auto count = static_cast<std::size_t>(n);
    const std::size_t row = raw.size() / count;  // bytes per point
    std::vector<std::int64_t> ids(count, 0);
    std::vector<char> scalar(count, 0);
    if (row == 4 * ws) {
        const Op2Record* order = op2_table_record(rTables, "EQEXIN", -3);
        const Op2Record* kinds = op2_table_record(rTables, "EQEXIN", -4);
        if (order == nullptr || kinds == nullptr)
            return false;
        const auto seq = w.Ints(order->Join(rS.Base()));
        const auto type = w.Ints(kinds->Join(rS.Base()));
        std::unordered_map<std::int64_t, std::int64_t> type_of;
        for (std::size_t i = 0; i + 1 < type.size(); i += 2)
            type_of[type[i]] = type[i + 1] % 10;
        for (std::size_t i = 0; i + 1 < seq.size(); i += 2) {
            const std::int64_t k = seq[i + 1] - 1;
            if (k < 0 || k >= n)
                return false;
            ids[static_cast<std::size_t>(k)] = seq[i];
            scalar[static_cast<std::size_t>(k)] = type_of[seq[i]] == 2;
        }
    } else if (row != 12 * 4 && !(ws == 8 && row == 9 * 8)) {
        return false;
    }
    for (std::size_t k = 0; k < count; ++k) {
        const char* e = raw.data() + k * row;
        const std::int64_t id = row == 4 * ws ? ids[k] : w.Int(e + 2 * ws);
        if (id <= 0)
            return false;
        if (scalar[k]) {
            rOut.mScalarPoints.insert(id);
            continue;
        }
        rOut.mIds.push_back(id);
        rOut.mCp.push_back(0);
        rOut.mCd.push_back(std::max<std::int64_t>(w.Int(e), 0));
        for (std::size_t d = 0; d < 3; ++d)
            rOut.mXyz.push_back(row == 4 * ws ? w.Float(e + (1 + d) * ws)
                                              : w.Double(e + 6 * ws + 8 * d));
    }
    return !rOut.mIds.empty();
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
        // A spring or damper's ends may be scalar points, or 0 (grounded).
        const bool scalar_ends = std::strncmp(spec->mCard, "CELAS", 5) == 0 ||
                                 std::strncmp(spec->mCard, "CDAMP", 5) == 0;
        const std::size_t corners =
            scalar_ends ? 0 : detail::nastran_card_spec(spec->mCard)->mLinearNodes;
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
                                                              {20, "CDAMP1"},
                                                              {21, "CDAMP2"},
                                                              {22, "CDAMP3"},
                                                              {23, "CDAMP4"},
                                                              {24, "CVISC"},
                                                              {38, "CGAP"},
                                                              {53, "CTRIAX6"},
                                                              {69, "CBEND"},
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
                                                              {107, "CHBDYE"},
                                                              {108, "CHBDYG"},
                                                              {109, "CHBDYP"},
                                                              {110, "CONV"},
                                                              {126, "CFAST"},
                                                              {227, "CTRIAR"},
                                                              {228, "CQUADR"},
                                                              {232, "CQUADR composite"},
                                                              {233, "CTRIAR composite"},
                                                              {280, "CBEAR"},
                                                              {300, "CHEXA"},
                                                              {301, "CPENTA"},
                                                              {302, "CTETRA"},
                                                              {303, "CPYRAM"},
                                                              {144, "CQUAD4 corner"},
                                                              {255, "CPYRAM"}};
    const auto it = names.find(Type);
    return it != names.end() ? it->second : "type " + std::to_string(Type);
}

// One value of an element entry: its word (within the entry, or within a node
// block), the word of its imaginary part in a complex table (kOp2Npos: real),
// and its member name.
struct Op2Member {
    std::size_t mWord;
    std::size_t mImag;
    std::string mName;
};
using Op2Layout = std::vector<Op2Member>;

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
    bool mCentre = true;      // Blocks: the first block is the element centre
};

// What an element table holds: stresses or strains (OES, OSTR), forces (OEF),
// heat fluxes (a heat-transfer OEF) or energies (ONRGY, OEKE).
enum class Op2Family { Stress, Force, Flux, Energy };

// Real members `rNames` at words First, First + 1 ...
Op2Layout op2_real(const std::vector<std::string>& rNames, std::size_t First) {
    Op2Layout out;
    for (std::size_t k = 0; k < rNames.size(); ++k)
        out.push_back({First + k, kOp2Npos, rNames[k]});
    return out;
}

// Complex members, all real parts from First, then all imaginary parts.
Op2Layout op2_block(const std::vector<std::string>& rNames, std::size_t First) {
    Op2Layout out;
    for (std::size_t k = 0; k < rNames.size(); ++k)
        out.push_back({First + k, First + rNames.size() + k, rNames[k]});
    return out;
}

// Complex members as (real, imaginary) pairs from First.
Op2Layout op2_pairs(const std::vector<std::string>& rNames, std::size_t First) {
    Op2Layout out;
    for (std::size_t k = 0; k < rNames.size(); ++k)
        out.push_back({First + 2 * k, First + 2 * k + 1, rNames[k]});
    return out;
}

Op2Layout op2_join(Op2Layout A, const Op2Layout& rB) {
    A.insert(A.end(), rB.begin(), rB.end());
    return A;
}

std::optional<Op2ElementLayout> op2_element_layout(Op2Family Family, std::int64_t Type,
                                                   std::int64_t NumWide, std::int64_t SCode,
                                                   bool Complex, bool Random) {
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
    const auto is = [&](std::initializer_list<std::int64_t> Types) {
        return std::find(Types.begin(), Types.end(), Type) != Types.end();
    };
    Op2ElementLayout out;
    const auto rows = [&](Op2Layout Members) {
        out.mMembers = std::move(Members);
        return out;
    };
    const auto blocks = [&](Op2Layout Members, std::size_t First, std::size_t Block,
                            std::size_t Count) {
        out.mKind = Op2Values::Blocks;
        out.mFirst = First;
        out.mBlock = Block;
        out.mBlocks = Count;
        out.mMembers = std::move(Members);
        return out;
    };
    const std::int64_t cn = corner_nodes(Type), sn = solid_nodes(Type);
    const std::vector<std::string> shell_force = {"MX",  "MY",   "MXY", "BMX",
                                                  "BMY", "BMXY", "TX",  "TY"};
    const std::vector<std::string> bush = {"FX", "FY", "FZ", "MX", "MY", "MZ"};

    if (Family == Op2Family::Energy) {
        if (!Complex && NumWide == 4)
            return rows(op2_real({"ENERGY", "PCT", "DEN"}, 1));
        if (Complex && NumWide == 5)
            return rows(op2_join({{1, 2, "ENERGY"}}, op2_real({"PCT", "DEN"}, 3)));
        return std::nullopt;
    }

    if (Family == Op2Family::Flux) {
        const std::vector<std::string> flux = {"XGRAD", "YGRAD", "ZGRAD",
                                               "XFLUX", "YFLUX", "ZFLUX"};
        if (is({1, 2, 3, 10, 34, 69, 33, 53, 64, 74, 75}) && NumWide == 9)
            return rows(op2_real(flux, 3));
        if (is({39, 67, 68}) && NumWide == 10)
            return rows(op2_real(flux, 3));
        if (is({107, 108, 109}) && NumWide == 8)
            return rows(op2_real({"FAPPLIED", "FREECONV", "FORCECONV", "FRAD", "FTOTAL"}, 3));
        if (Type == 110 && NumWide == 4)
            return rows({{1, kOp2Npos, "FREECONV"}, {3, kOp2Npos, "FREECONVK"}});
        return std::nullopt;
    }

    if (Family == Op2Family::Force) {
        // Random force tables use the real layouts.
        const auto real_or_complex = [&](const std::vector<std::string>& rNames,
                                         std::int64_t RealWide) -> std::optional<Op2ElementLayout> {
            if (!Complex && NumWide == RealWide)
                return rows(op2_real(rNames, 1));
            if (Complex && NumWide == 1 + 2 * static_cast<std::int64_t>(rNames.size()))
                return rows(op2_block(rNames, 1));
            return std::nullopt;
        };
        if (is({1, 3, 10, 24}))
            return real_or_complex({"AF", "TRQ"}, 3);
        if (is({11, 12, 13, 14, 20, 21, 22, 23}))
            return real_or_complex({"F"}, 2);
        if (Type == 34)
            return real_or_complex({"BM1A", "BM2A", "BM1B", "BM2B", "TS1", "TS2", "AF", "TRQ"}, 9);
        if (Type == 100 && !Complex && NumWide == 8) {
            out.mKind = Op2Values::Station;
            return rows(op2_real({"SD", "BM1", "BM2", "TS1", "TS2", "AF", "TRQ"}, 1));
        }
        if (Type == 4) {
            const std::vector<std::string> forces = {"F41", "F21", "F12", "F32",
                                                     "F23", "F43", "F34", "F14"};
            const std::vector<std::string> kicks = {"KF1", "S12", "KF2", "S23",
                                                    "KF3", "S34", "KF4", "S41"};
            if (!Complex && NumWide == 17)
                return rows(op2_join(op2_real(forces, 1), op2_real(kicks, 9)));
            // The forces' real then imaginary parts, then the kick forces' and
            // shear flows' (pyNastran pairs words 1-16 with 17-32 instead).
            if (Complex && NumWide == 33)
                return rows(op2_join(op2_block(forces, 1), op2_block(kicks, 17)));
            return std::nullopt;
        }
        if (is({33, 74, 227, 228}))
            return real_or_complex(shell_force, 9);
        if (is({102, 126, 280}))
            return real_or_complex(bush, 7);
        if (Type == 38 && !Complex && NumWide == 9)
            return rows(op2_real({"FX", "SFY", "SFZ", "U", "V", "W", "SV", "SW"}, 1));
        if (Type == 2) {
            const std::vector<std::string> beam = {"BM1", "BM2",  "TS1", "TS2",
                                                   "AF",  "TTRQ", "WTRQ"};
            out.mStations = true;
            if (!Complex && NumWide == 100)
                return blocks(op2_join({{1, kOp2Npos, "SD"}}, op2_real(beam, 2)), 1, 9, 11);
            if (Complex && NumWide == 177)
                return blocks(op2_join({{1, kOp2Npos, "SD"}}, op2_block(beam, 2)), 1, 16, 11);
            return std::nullopt;
        }
        if (cn && !Complex && NumWide == 2 + 9 * cn)
            return blocks(op2_real(shell_force, 1), 2, 9, static_cast<std::size_t>(cn));
        if (cn && Complex && NumWide == 2 + 17 * cn)
            return blocks(op2_block(shell_force, 1), 2, 17, static_cast<std::size_t>(cn));
        return std::nullopt;
    }

    // Stresses and strains.
    const Op2Layout plate_c =
        op2_join(op2_join({{1, kOp2Npos, "FD1"}}, op2_pairs({"X1", "Y1", "TXY1"}, 2)),
                 op2_join({{8, kOp2Npos, "FD2"}}, op2_pairs({"X2", "Y2", "TXY2"}, 9)));
    const Op2Layout plate_cvm =
        op2_join(op2_join(op2_join({{1, kOp2Npos, "FD1"}}, op2_pairs({"X1", "Y1", "TXY1"}, 2)),
                          {{8, kOp2Npos, "VON_MISES1"}, {9, kOp2Npos, "FD2"}}),
                 op2_join(op2_pairs({"X2", "Y2", "TXY2"}, 10), {{16, kOp2Npos, "VON_MISES2"}}));
    const std::vector<std::string> tensor = {"X", "Y", "Z", "TXY", "TYZ", "TZX"};
    const std::vector<std::string> ply = {"X1", "Y1", "T1", "L1", "L2"};
    const std::vector<std::string> bush_stress = {"TX", "TY", "TZ", "RX", "RY", "RZ"};
    if (Random) {
        // Magnitudes only (PSD, RMS ...): pyNastran's random layouts.
        const std::vector<std::string> rp = {"FD1", "X1", "Y1", "TXY1", "FD2", "X2", "Y2", "TXY2"};
        const std::vector<std::string> rpvm = {"FD1", "X1", "Y1", "TXY1", "VON_MISES1",
                                               "FD2", "X2", "Y2", "TXY2", "VON_MISES2"};
        if (is({1, 10}) && NumWide == 3)
            return rows(op2_real({"A", "T"}, 1));
        if (Type == 3 && NumWide == 3)
            return rows(op2_real({"AS", "TS"}, 1));
        if (Type == 4 && NumWide == 3)
            return rows(op2_real({"TMAX", "TAVG"}, 1));
        if (Type == 34 && NumWide == 10)
            return rows(
                op2_real({"X1A", "X2A", "X3A", "X4A", "AX", "X1B", "X2B", "X3B", "X4B"}, 1));
        if (Type == 2 && NumWide == 67) {
            out.mStations = true;
            return blocks(op2_real({"SD", "XC", "XD", "XE", "XF"}, 1), 1, 6, 11);
        }
        if (is({33, 74, 227, 228}) && NumWide == 9)
            return rows(op2_real(rp, 1));
        if (is({33, 74, 227, 228}) && NumWide == 11)
            return rows(op2_real(rpvm, 1));
        if (cn && NumWide == 2 + 9 * cn)
            return blocks(op2_real(rp, 1), 2, 9, static_cast<std::size_t>(cn));
        if (cn && NumWide == 2 + 11 * cn)
            return blocks(op2_real(rpvm, 1), 2, 11, static_cast<std::size_t>(cn));
        if (sn && NumWide == 4 + 7 * sn)
            return blocks(op2_real(tensor, 1), 4, 7, static_cast<std::size_t>(sn));
        if (sn && NumWide == 4 + 8 * sn)
            return blocks(op2_join(op2_real(tensor, 1), {{7, kOp2Npos, "VON_MISES"}}), 4, 8,
                          static_cast<std::size_t>(sn));
        if (is({95, 96, 97, 98, 232, 233}) && (NumWide == 7 || NumWide == 8)) {
            out.mKind = Op2Values::Ply;
            Op2Layout m = op2_real(ply, 2);
            if (NumWide == 8)
                m.push_back({7, kOp2Npos, "VON_MISES"});
            return rows(std::move(m));
        }
        if (Type == 102 && NumWide == 7)
            return rows(op2_real(bush_stress, 1));
        return std::nullopt;
    }
    if (Complex) {
        if (is({1, 10}) && NumWide == 5)
            return rows(op2_pairs({"A", "T"}, 1));
        if (Type == 3 && NumWide == 5)
            return rows(op2_pairs({"AS", "TS"}, 1));
        if (is({11, 12, 13, 14}) && NumWide == 3)
            return rows({{1, 2, "S"}});
        if (Type == 4 && NumWide == 5)
            return rows(op2_pairs({"TMAX", "TAVG"}, 1));
        if (Type == 34 && NumWide == 19)
            return rows(op2_join(op2_block({"X1A", "X2A", "X3A", "X4A", "AX"}, 1),
                                 op2_block({"X1B", "X2B", "X3B", "X4B"}, 11)));
        if (Type == 2 && NumWide == 111) {
            out.mStations = true;
            return blocks(op2_join({{1, kOp2Npos, "SD"}}, op2_block({"XC", "XD", "XE", "XF"}, 2)),
                          1, 10, 11);
        }
        if (is({33, 74, 227, 228}) && NumWide == 15)
            return rows(plate_c);
        if (is({33, 74, 227, 228}) && NumWide == 17)
            return rows(plate_cvm);
        if (cn && NumWide == 2 + 15 * cn)
            return blocks(plate_c, 2, 15, static_cast<std::size_t>(cn));
        if (cn && NumWide == 2 + 17 * cn)
            return blocks(plate_cvm, 2, 17, static_cast<std::size_t>(cn));
        if (sn && NumWide == 4 + 13 * sn)
            return blocks(op2_block(tensor, 1), 4, 13, static_cast<std::size_t>(sn));
        if (sn && NumWide == 4 + 14 * sn)
            return blocks(op2_join(op2_block(tensor, 1), {{13, kOp2Npos, "VON_MISES"}}), 4, 14,
                          static_cast<std::size_t>(sn));
        if (is({95, 96, 97, 98, 232, 233}) && (NumWide == 12 || NumWide == 13)) {
            out.mKind = Op2Values::Ply;
            Op2Layout m = op2_block(ply, 2);
            if (NumWide == 13)
                m.push_back({12, kOp2Npos, "VON_MISES"});
            return rows(std::move(m));
        }
        if (Type == 102 && NumWide == 13)
            return rows(op2_block(bush_stress, 1));
        return std::nullopt;
    }
    if ((Type == 1 || Type == 10) && NumWide == 5)
        return rows(op2_real({"A", "MSA", "T", "MST"}, 1));
    if (Type == 3 && NumWide == 5)
        return rows(op2_real({"AS", "MSA", "TS", "MST"}, 1));
    if (Type == 4 && NumWide == 4)
        return rows(op2_real({"TMAX", "TAVG", "MS"}, 1));
    if (is({11, 12, 13, 14}) && NumWide == 2)
        return rows(op2_real({"S"}, 1));
    if (Type == 102 && NumWide == 7)
        return rows(op2_real(bush_stress, 1));
    if (Type == 34 && NumWide == 16) {
        for (std::size_t k = 0; k < 15; ++k)
            out.mMembers.push_back({1 + k, kOp2Npos, bar[k]});
        return out;
    }
    // CBAR stations (100): sd, the four fibres, axial, max, min, margin.
    if (Type == 100 && NumWide == 10) {
        out.mKind = Op2Values::Station;
        return rows(op2_real({"SD", "XC", "XD", "XE", "XF", "AX", "MAX", "MIN", "MS"}, 1));
    }
    // Composite shells (QUAD4, QUAD8, TRIA3, TRIA6; NX QUADR, TRIAR): one entry
    // per ply, the MSC HDF5 names.
    if (is({95, 96, 97, 98, 232, 233}) && NumWide == 11) {
        out.mKind = Op2Values::Ply;
        return rows(op2_real({"X1", "Y1", "T1", "L1", "L2", "ANGLE", "MAJOR", "MINOR", vm}, 2));
    }
    // CBEAM (2): 11 stations of grid, sd, the four fibres, max, min and margins.
    if (Type == 2 && NumWide == 111) {
        out.mStations = true;
        return blocks(op2_real({"SD", "XC", "XD", "XE", "XF", "MAX", "MIN", "MST", "MSC"}, 1), 1,
                      10, 11);
    }
    Op2Layout shell;
    for (std::size_t k = 0; k < 16; ++k)
        shell.push_back(
            {1 + k, kOp2Npos, plate[k] ? std::string(plate[k]) : vm + std::to_string(1 + k / 8)});
    if (is({33, 74}) && NumWide == 17)
        return rows(shell);
    if (cn && NumWide == 2 + 17 * cn)
        return blocks(shell, 2, 17, static_cast<std::size_t>(cn));
    if (sn && NumWide == 4 + 21 * sn) {
        const std::string octa = (SCode & 1) ? "VON_MISES" : "OCT_SHEAR";
        Op2Layout m;
        for (const auto& [word, name] : solid)
            m.push_back({word, kOp2Npos, name ? std::string(name) : octa});
        return blocks(std::move(m), 4, 21, static_cast<std::size_t>(sn));
    }
    // NX's newer solids (CHEXA 300, CPENTA 301, CTETRA 302, CPYRAM 303): the
    // corners only, no centre.
    const std::int64_t nx = Type == 300   ? 8
                            : Type == 301 ? 6
                            : Type == 302 ? 4
                            : Type == 303 ? 5
                                          : 0;
    if (nx && NumWide == 3 + 8 * nx) {
        out.mCentre = false;
        return blocks(op2_join(op2_real(tensor, 1), {{7, kOp2Npos, "VON_MISES"}}), 3, 8,
                      static_cast<std::size_t>(nx));
    }
    return std::nullopt;
}

// A complex value's real and imaginary parts from its two words: themselves,
// or a magnitude and a phase in degrees.
std::pair<double, double> op2_complex(double A, double B, bool MagPhase) {
    if (!MagPhase)
        return {A, B};
    constexpr double kDegree = 3.14159265358979323846 / 180.0;
    return {A * std::cos(B * kDegree), A * std::sin(B * kDegree)};
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
    std::int64_t mW5;    // its word 5 (a SORT2 row's first word): the key of its rows
    double mEigi = 0.0;  // a complex mode's imaginary eigenvalue
};

struct Op2Block {
    std::size_t mStep = 0;
    bool mNodal = false;
    bool mBasic = false;  // BOUG*: already in the basic system
    std::string mBase;    // nodal: point data name; element: STRESS/STRAIN
    Op2ElementLayout mLayout;
    std::size_t mNumWide = 8;
    bool mGridForce = false;  // OGPFB: grid point forces
    const Op2Record* mRecord = nullptr;
    bool mComplex = false;   // real and imaginary parts (format 2) ...
    bool mMagPhase = false;  // ... or magnitude and phase in degrees (format 3)
    std::string mSuffix;     // a random table's: _PSD, _ATO, _RMS, _NO, _CRM
    // SORT2: one entity's rows over every step, each led by the step's word 5;
    // the entity's id * 10 + device, as the file writes it.
    bool mSort2 = false;
    std::int64_t mSubcase = 0, mAnalysis = 0;
    std::string mEntity;
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
            std::optional<Op2Step> mOwn;  // the step it makes when none matches
        };
        std::vector<Deferred> deferred;
        // Complex eigenvalues (CLAMA): mode -> (real, imaginary part); rows are
        // `mode, order, eigr, eigi, frequency, damping`.
        std::map<std::int64_t, std::pair<double, double>> clama;
        for (const Op2Table& t : mTables) {
            if (!op2_starts(t.mName, {"CLAMA"}))
                continue;
            for (const auto& [mark, rec] : t.mRecords) {
                if (mark > -3 || rec.Size() == 146 * ws || rec.Size() % (6 * ws))
                    continue;
                const std::string raw = rec.Join(mStream.Base());
                for (std::size_t r0 = 0; r0 < raw.size(); r0 += 6 * ws)
                    clama.emplace(w.Int(raw.data() + r0),
                                  std::make_pair(w.Float(raw.data() + r0 + 2 * ws),
                                                 w.Float(raw.data() + r0 + 3 * ws)));
            }
        }
        for (const Op2Table& t : mTables) {
            const bool nodal = op2_starts(t.mName, {"OUG", "BOUG", "OQG", "OQMG", "OPG"});
            const bool elemental =
                op2_starts(t.mName, {"OES", "OSTR", "OEF", "HOEF", "ONRGY", "OEKE"});
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
                const std::int64_t format_code = word(8);
                const std::int64_t num_wide = word(9);
                const std::int64_t s_code = word(10);
                const std::int64_t thermal = word(22);
                const bool complex = (sort_code & 1) != 0;
                const bool sort2 = (sort_code & 2) != 0;
                // Random tables: the hundreds of the table code name the
                // quantity (5 CRM, 6 PSD, 7 ATO, 8 RMS, 9 NO), the rest what.
                const bool random = (sort_code & 4) != 0 || table_code >= 500;
                std::int64_t what = table_code;
                std::string suffix;
                if (random) {
                    static const char* const kRandom[] = {"_CRM", "_PSD", "_ATO", "_RMS", "_NO"};
                    const std::int64_t kind = table_code / 100;
                    if (kind < 5 || kind > 9 || complex || grid_force) {
                        Skip(t.mName + " (random)");
                        continue;
                    }
                    suffix = kRandom[kind - 5];
                    what = table_code % 100;
                }
                if (complex && grid_force) {
                    Skip(t.mName + " (complex)");
                    continue;
                }
                Op2Block block;
                block.mNodal = nodal;
                block.mBasic = op2_starts(t.mName, {"BOUG"});
                block.mRecord = &rec;
                block.mComplex = complex;
                block.mMagPhase = complex && format_code == 3;
                block.mSuffix = suffix;
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
                    const char* base = op2_nodal_name(what);
                    // MPC forces share the SPC forces' table code; the name tells.
                    if (op2_starts(t.mName, {"OQMG"}) && (what == 3 || what == 39))
                        base = "MPC_FORCE";
                    if (base == nullptr) {
                        Skip(t.mName + " (table code " + std::to_string(table_code) + ")");
                        continue;
                    }
                    if (num_wide != (complex ? 14 : 8)) {
                        Skip(t.mName + " (" + std::to_string(num_wide) + " words per node)");
                        continue;
                    }
                    if (thermal == 1) {
                        if (what != 1 || complex) {
                            Skip(t.mName + " (thermal table code " + std::to_string(table_code) +
                                 ")");
                            continue;
                        }
                        base = "TEMPERATURE";
                    }
                    block.mBase = base;
                    block.mNumWide = static_cast<std::size_t>(num_wide);
                } else {
                    Op2Family family = Op2Family::Stress;
                    if (what == 5) {
                        block.mBase = (s_code & 8) ? "STRAIN" : "STRESS";
                    } else if (what == 4) {
                        family = thermal == 1 ? Op2Family::Flux : Op2Family::Force;
                        block.mBase = thermal == 1 ? "HEAT_FLUX" : "ELEMENT_FORCE";
                    } else if (what == 18 || what == 36) {
                        family = Op2Family::Energy;
                        block.mBase = what == 18 ? "ENERGY" : "KINETIC_ENERGY";
                    } else {
                        Skip(t.mName + " (table code " + std::to_string(table_code) + ")");
                        continue;
                    }
                    auto layout = op2_element_layout(family, etype, num_wide, s_code, complex,
                                                     random && family == Op2Family::Stress);
                    if (!layout) {
                        // An energy table's word 3 is the total energy, not a type.
                        Skip(t.mName + " " +
                             (family == Op2Family::Energy ? std::string("energy")
                                                          : op2_element_type_name(etype)) +
                             (complex ? " complex" : "") + " (" + std::to_string(num_wide) +
                             " words)");
                        continue;
                    }
                    block.mLayout = std::move(*layout);
                    block.mNumWide = static_cast<std::size_t>(num_wide);
                }
                const std::int64_t w5 = word(4);
                const bool moded = analysis == 2 || analysis == 8 || analysis == 9;
                if (sort2) {
                    if (grid_force) {
                        Skip(t.mName + " (SORT2)");
                        continue;
                    }
                    // Every row makes (or joins) the step its first word names.
                    const std::string raw = rec.Join(mStream.Base());
                    const std::size_t width = block.mNumWide * ws;
                    if (width == 0 || raw.size() % width)
                        op2_fail("a SORT2 " + t.mName + " record is not a whole number of rows");
                    for (std::size_t r0 = 0; r0 < raw.size(); r0 += width) {
                        const std::int64_t key5 = w.Int(raw.data() + r0);
                        const auto key = std::make_tuple(subcase, analysis, key5);
                        if (step_index.count(key))
                            continue;
                        step_index.emplace(key, mSteps.size());
                        const double time = w.Float(raw.data() + r0);
                        mSteps.push_back({subcase, analysis, moded ? key5 : 0,
                                          op2_time_of(analysis, time, time), key5});
                    }
                    block.mSort2 = true;
                    block.mSubcase = subcase;
                    block.mAnalysis = analysis;
                    block.mEntity.assign(h + 4 * ws, ws);
                    mBlocks.push_back(std::move(block));
                    continue;
                }
                const auto key = std::make_tuple(subcase, analysis, w5);
                auto it = step_index.find(key);
                const Op2Step own{subcase, analysis, moded ? w5 : 0,
                                  op2_time_of(analysis, w.Float(h + 4 * ws), w.Float(h + 5 * ws)),
                                  w5};
                if (grid_force) {
                    deferred.push_back({key, std::move(block), t.mName, std::nullopt});
                    continue;
                }
                // Energies too: ONRGY writes 0 in word 5 where the other tables
                // of a static step write the load set.
                if (block.mBase == "ENERGY" || block.mBase == "KINETIC_ENERGY") {
                    deferred.push_back({key, std::move(block), t.mName, own});
                    continue;
                }
                if (it == step_index.end()) {
                    it = step_index.emplace(key, mSteps.size()).first;
                    mSteps.push_back(own);
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
                if (same == 0 && d.mOwn) {
                    it = step_index.emplace(d.mKey, mSteps.size()).first;
                    mSteps.push_back(*d.mOwn);
                } else if (same != 1) {
                    Skip(d.mTable + " (no matching step)");
                    continue;
                }
            }
            d.mBlock.mStep = it->second;
            mBlocks.push_back(std::move(d.mBlock));
        }
        for (Op2Step& st : mSteps)
            if (st.mAnalysis == 9)
                if (const auto it = clama.find(st.mMode); it != clama.end()) {
                    st.mTime = it->second.first;
                    st.mEigi = it->second.second;
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
    bool from_bgpdt = false;
    if (g.mIds.empty()) {
        auto [deck, tried] = op2_sibling_deck(rR.mPath);
        if (deck.empty() && op2_bgpdt_grids(rR.mStream, rR.mTables, g)) {
            from_bgpdt = true;
        } else if (deck.empty()) {
            std::string list;
            for (const std::string& t : tried)
                list += (list.empty() ? "" : ", ") + t;
            op2_fail(
                "the file has no GEOM1 GRID records and no basic grid point table (rerun with "
                "PARAM,POST,-1 or provide the input deck beside it); looked for " +
                list);
        }
    }
    if (g.mIds.empty()) {
        auto [deck, tried] = op2_sibling_deck(rR.mPath);
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
    // A file can repeat a GRID in a second GEOM1 table (a restart's): kept
    // once when both definitions agree.
    {
        Op2Grids unique;
        unique.mScalarPoints = g.mScalarPoints;
        unique.mCords = g.mCords;
        std::unordered_map<std::int64_t, std::size_t> seen;
        for (std::size_t i = 0; i < g.mIds.size(); ++i) {
            const auto [it, fresh] = seen.emplace(g.mIds[i], i);
            if (!fresh) {
                const std::size_t j = it->second;
                const bool same = g.mCp[i] == g.mCp[j] && g.mCd[i] == g.mCd[j] &&
                                  std::equal(&g.mXyz[3 * i], &g.mXyz[3 * i + 3], &g.mXyz[3 * j]);
                if (!same)
                    op2_fail("GRID " + std::to_string(g.mIds[i]) + " is defined twice");
                continue;
            }
            unique.mIds.push_back(g.mIds[i]);
            unique.mCp.push_back(g.mCp[i]);
            unique.mCd.push_back(g.mCd[i]);
            unique.mXyz.insert(unique.mXyz.end(), &g.mXyz[3 * i], &g.mXyz[3 * i + 3]);
        }
        g = std::move(unique);
    }
    for (std::size_t i = 0; i < g.mIds.size(); ++i)
        m.mGridIndex.emplace(g.mIds[i], i);
    NDArray points(DType::Float64, {g.mIds.size(), 3});
    std::copy(g.mXyz.begin(), g.mXyz.end(), points.As<double>());
    m.mMesh.AssignPoints(std::move(points));
    if (from_bgpdt) {
        log::warn(
            "{}: no GEOM1 table and no input deck beside the file; the {} points come from "
            "the basic grid point table (BGPDT)",
            kOp2Who, g.mIds.size());
        if (std::any_of(g.mCd.begin(), g.mCd.end(), [](std::int64_t C) { return C > 0; }))
            log::warn(
                "{}: without GEOM1 the output coordinate systems are unknown; results "
                "stay in them",
                kOp2Who);
        std::fill(g.mCd.begin(), g.mCd.end(), 0);
    }
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
        return std::move(mesh);  // a reference into the local model
    }
    const std::size_t index = rOpts.ResolveTimeStep(n);
    const Op2Step& step = r.mSteps[index];
    mesh.AddFieldData(kSequenceTimeKey, op2_scalar(step.mTime, DType::Float64));
    mesh.AddFieldData("nastran:subcase",
                      op2_scalar(static_cast<double>(step.mSubcase), DType::Int64));
    mesh.AddFieldData("nastran:analysis",
                      op2_scalar(static_cast<double>(step.mAnalysis), DType::Int64));
    mesh.AddFieldData("nastran:mode", op2_scalar(static_cast<double>(step.mMode), DType::Int64));
    if (step.mAnalysis == 9)
        mesh.AddFieldData("nastran:eigi", op2_scalar(step.mEigi, DType::Float64));
    op2_warn_skipped(r);
    if (rOpts.mPointsOnly)
        return std::move(mesh);  // a reference into the local model

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
        std::string raw;
        if (b.mSort2) {
            // This step's rows of the entity, rewritten as SORT1 rows.
            if (b.mSubcase != step.mSubcase || b.mAnalysis != step.mAnalysis)
                continue;
            const std::string all = b.mRecord->Join(r.mStream.Base());
            const std::size_t width = b.mNumWide * ws;
            for (std::size_t r0 = 0; r0 + width <= all.size(); r0 += width)
                if (w.Int(all.data() + r0) == step.mW5) {
                    raw += b.mEntity;
                    raw.append(all, r0 + ws, width - ws);
                }
            if (raw.empty())
                continue;
        } else {
            if (b.mStep != index)
                continue;
            raw = b.mRecord->Join(r.mStream.Base());
        }
        const std::size_t nwords = raw.size() / ws;
        if (b.mNodal) {
            const std::size_t nw = b.mNumWide;
            if (nwords % nw)
                op2_fail("a " + b.mBase + " record is not a whole number of rows");
            // (name, first word of the value or real part, of the imaginary
            // part, components, part: 0 the value, 1 the real, 2 the imaginary)
            struct Output {
                std::string mName;
                std::size_t mRe, mIm, mCount;
                int mPart;
            };
            std::vector<Output> outputs;
            if (b.mBase == "TEMPERATURE") {
                outputs.push_back({b.mBase + b.mSuffix, 2, 0, 1, 0});
            } else if (b.mComplex) {
                outputs.push_back({b.mBase + "_real", 2, 8, 3, 1});
                outputs.push_back({b.mBase + "_imag", 2, 8, 3, 2});
                outputs.push_back({b.mBase + "_ROT_real", 5, 11, 3, 1});
                outputs.push_back({b.mBase + "_ROT_imag", 5, 11, 3, 2});
            } else {
                outputs.push_back({b.mBase + b.mSuffix, 2, 0, 3, 0});
                outputs.push_back({b.mBase + "_ROT" + b.mSuffix, 5, 0, 3, 0});
            }
            for (const Output& o : outputs) {
                if (!rOpts.WantsArray(o.mName))
                    continue;
                const std::size_t nc = o.mCount;
                auto it = point_arrays.find(o.mName);
                if (it == point_arrays.end())
                    it = point_arrays
                             .emplace(o.mName,
                                      std::make_pair(nc, std::vector<double>(npts * nc, kOp2Nan)))
                             .first;
                std::vector<double>& values = it->second.second;
                for (std::size_t r0 = 0; r0 < nwords / nw; ++r0) {
                    const char* row = raw.data() + r0 * nw * ws;
                    const auto p = model.mGridIndex.find(w.Int(row) / 10);
                    if (p == model.mGridIndex.end() || !std::isnan(values[p->second * nc]))
                        continue;
                    double* dst = values.data() + p->second * nc;
                    for (std::size_t c = 0; c < nc; ++c) {
                        const double a = w.Float(row + (o.mRe + c) * ws);
                        if (o.mPart == 0) {
                            dst[c] = a;
                            continue;
                        }
                        const double bb = w.Float(row + (o.mIm + c) * ws);
                        const auto [re, im] = op2_complex(a, bb, b.mMagPhase);
                        dst[c] = o.mPart == 1 ? re : im;
                    }
                    // Results are in the GRID's output system (CD) unless the
                    // table is BOUG*; random ones (spectral densities, RMS ...)
                    // are not vectors and stay there.
                    if (nc == 3 && !b.mBasic && !model.mCd.empty() && b.mSuffix.empty())
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
            // A member's values at pBase: itself (with a random table's suffix),
            // or its real and imaginary parts.
            const auto values = [&](const char* pBase, const Op2Member& rM, const auto& rEmit) {
                const double a = w.Float(pBase + rM.mWord * ws);
                if (rM.mImag == kOp2Npos) {
                    rEmit(rM.mName + b.mSuffix, a);
                    return;
                }
                const auto [re, im] = op2_complex(a, w.Float(pBase + rM.mImag * ws), b.mMagPhase);
                rEmit(rM.mName + "_real", re);
                rEmit(rM.mName + "_imag", im);
            };
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
                    const std::string suffix = ply ? "@ply" : "@station";
                    for (const Op2Member& m : lay.mMembers)
                        values(row, m, [&](const std::string& rMember, double V) {
                            const std::string name = b.mBase + ":" + rMember + suffix;
                            if (rOpts.WantsArray(name))
                                push(name, cell, col, V);
                        });
                    continue;
                }
                // Row: the values; Blocks: the first block's (the centre, or end A).
                if (lay.mKind != Op2Values::Blocks || lay.mCentre) {
                    const char* at = row + (lay.mKind == Op2Values::Blocks ? lay.mFirst * ws : 0);
                    for (const Op2Member& m : lay.mMembers)
                        values(at, m, [&](const std::string& rMember, double V) {
                            const std::string name = b.mBase + ":" + rMember;
                            if (!rOpts.WantsArray(name))
                                return;
                            auto it = cell_arrays.find(name);
                            if (it == cell_arrays.end())
                                it = cell_arrays.emplace(name, std::vector<double>(ncells, kOp2Nan))
                                         .first;
                            if (std::isnan(it->second[cell]))
                                it->second[cell] = V;
                        });
                }
                if (lay.mKind != Op2Values::Blocks)
                    continue;
                for (std::size_t k = (lay.mStations || !lay.mCentre) ? 0 : 1; k < lay.mBlocks;
                     ++k) {
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
                    const std::string suffix = lay.mStations ? "@station" : "@corner";
                    for (const Op2Member& m : lay.mMembers)
                        values(block, m, [&](const std::string& rMember, double V) {
                            const std::string name = b.mBase + ":" + rMember + suffix;
                            if (rOpts.WantsArray(name))
                                push(name, cell, col, V);
                        });
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
    return std::move(mesh);  // a reference into the local model
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
