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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "../detail/typed_view.hpp"

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
    {1001, 10, "CMASS1", {6}, 1, 2, 2, 1},
    {1101, 11, "CMASS2", {6}, 1, 2, 2, -1},
    {1401, 14, "CONM1", {24}, 1, 1, 1, -1},
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
    {1908, 19, "CGAP", {9}, 1, 2, 2, 1},
    {5201, 52, "PLOTEL", {3}, 1, 1, 2, -1},
    // NX 2019 and later write CQUAD4/CTRIA3 as these (pyNastran's CQUADRN/CTRIARN)
    {15401, 154, "CQUAD4", {14, 15}, 2, 2, 4, 1},
    {15301, 153, "CTRIA3", {13, 14}, 2, 2, 3, 1},
};

// Element records without a cell type, named in the warning.
constexpr std::tuple<int, int, const char*> kOp2NamedElements[] = {
    {801, 8, "CELAS3"},   {901, 9, "CELAS4"},     {401, 4, "CDAMP3"},    {501, 5, "CDAMP4"},
    {1201, 12, "CMASS3"}, {1301, 13, "CMASS4"},   {5608, 56, "CBUSH1D"}, {6108, 61, "CTRIAX6"},
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
    const std::unordered_set<std::int64_t>& rGrids, std::vector<std::vector<std::int64_t>>& rComps,
    std::vector<std::vector<double>>& rBar) {
    const Op2Words& w = rS.Words();
    std::map<std::string, detail::NastranCardRows> cards;
    // Per card, two components per element: a spring's or damper's at each end.
    std::map<std::string, std::vector<std::int64_t>> comps;
    // Per CBAR, twelve doubles per element: its orientation, G0, pins and offsets.
    std::map<std::string, std::vector<double>> bars;
    std::vector<std::string> skipped;
    std::size_t grid_oriented = 0;     // CBEAMs whose X vector is not read
    std::size_t element_oriented = 0;  // CBUSHes oriented by their element system
    std::size_t conm1_oriented = 0;    // CONM1s whose mass matrix is in another system
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
                                 std::strncmp(spec->mCard, "CDAMP", 5) == 0 ||
                                 std::strncmp(spec->mCard, "CMASS", 5) == 0;
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
        std::vector<std::int64_t>& card_comps = comps[spec->mCard];
        std::vector<double>& card_bars = bars[spec->mCard];
        for (std::size_t i = 0; i < words.size() / chosen; ++i) {
            const std::int64_t* row = words.data() + i * chosen;
            rows.mEid.push_back(row[0]);
            rows.mPid.push_back(spec->mPidWord >= 0 ? row[spec->mPidWord] : -1);
            rows.mNodes.insert(rows.mNodes.end(), row + first, row + first + spec->mCount);
            // The component numbers follow the two ends of a spring or damper.
            const std::size_t after = first + static_cast<std::size_t>(spec->mCount);
            card_comps.push_back(scalar_ends ? row[after] : 0);
            card_comps.push_back(scalar_ends ? row[after + 1] : 0);
            // A float word of this row; the record's three key words come before the row.
            const std::size_t base = 3 + i * chosen;
            auto fword = [&](std::size_t Word) {
                return w.Float(r.mRaw.data() + (base + Word) * static_cast<std::size_t>(w.mWs));
            };
            const std::size_t row_start = card_bars.size();
            if (std::strcmp(spec->mCard, "CBAR") == 0) {
                // X1..X3 at words 4 to 6, or a G0 at 4 with the flag 2 at 7; PA and PB at
                // 8 and 9; the offsets W1A..W3A and W1B..W3B at 10 to 15.
                const bool g0 = row[7] == 2;
                for (std::size_t k = 0; k < 3; ++k)
                    card_bars.push_back(g0 ? 0.0 : fword(4 + k));
                card_bars.push_back(g0 ? static_cast<double>(row[4]) : 0.0);
                card_bars.push_back(static_cast<double>(row[8]));
                card_bars.push_back(static_cast<double>(row[9]));
                for (std::size_t k = 0; k < 6; ++k)
                    card_bars.push_back(fword(10 + k));
            } else if (std::strcmp(spec->mCard, "CBEAM") == 0) {
                // X1..X3 at words 6 to 8 when the flag at 9 names a vector (its low two bits:
                // 0 basic, 1 global), the G0 at word 6 when it is 2; PA and PB at 10 and 11;
                // the offsets at 12 to 17.
                const std::int64_t low = row[9] & 3;
                if (low == 3)
                    ++grid_oriented;
                for (std::size_t k = 0; k < 3; ++k)
                    card_bars.push_back(low <= 1 ? fword(6 + k) : 0.0);
                card_bars.push_back(low == 2 ? static_cast<double>(row[6]) : 0.0);
                card_bars.push_back(static_cast<double>(row[10]));
                card_bars.push_back(static_cast<double>(row[11]));
                for (std::size_t k = 0; k < 6; ++k)
                    card_bars.push_back(fword(12 + k));
            } else if (std::strcmp(spec->mCard, "CGAP") == 0) {
                // X1..X3 at words 4 to 6 when the flag at 7 is not 2, the G0 at word 4 when it is.
                const bool g0 = row[7] == 2;
                for (std::size_t k = 0; k < 3; ++k)
                    card_bars.push_back(g0 ? 0.0 : fword(4 + k));
                card_bars.push_back(g0 ? static_cast<double>(row[4]) : 0.0);
                card_bars.insert(card_bars.end(), 8, 0.0);
            } else if (std::strcmp(spec->mCard, "CBUSH") == 0) {
                // X1..X3 at words 4 to 6 when the flag at 7 is 0 or 1, the G0 at word 4 when it
                // is 2; a flag of -1 orients the bush by its element system, which is not read.
                const std::int64_t flag = row[7];
                if (flag == -1)
                    ++element_oriented;
                for (std::size_t k = 0; k < 3; ++k)
                    card_bars.push_back((flag == 0 || flag == 1) ? fword(4 + k) : 0.0);
                card_bars.push_back(flag == 2 ? static_cast<double>(row[4]) : 0.0);
                card_bars.insert(card_bars.end(), 8, 0.0);
            } else if (std::strcmp(spec->mCard, "CELAS2") == 0) {
                // K at word 1, GE at 6 and S at 7.
                card_bars.push_back(fword(1));
                card_bars.push_back(fword(6));
                card_bars.push_back(fword(7));
                card_bars.insert(card_bars.end(), 9, 0.0);
            } else if (std::strcmp(spec->mCard, "CDAMP2") == 0 ||
                       std::strcmp(spec->mCard, "CMASS2") == 0) {
                // B (a damper's) or M (a mass's) at word 1.
                card_bars.push_back(fword(1));
                card_bars.insert(card_bars.end(), 11, 0.0);
            } else if (std::strcmp(spec->mCard, "CONM1") == 0) {
                // The mass matrix M11, M21, M22, M31, ... M66 at words 3 to 23; a coordinate
                // system CID (word 2) other than the basic one is not read.
                if (row[2] > 0)
                    ++conm1_oriented;
                for (std::size_t k = 0; k < 21; ++k)
                    card_bars.push_back(fword(3 + k));
            }
            // Every row of a card with fields takes 21 columns, a CONM1's mass terms.
            if (card_bars.size() > row_start)
                card_bars.resize(row_start + 21, 0.0);
        }
    }
    if (!skipped.empty()) {
        std::string list;
        for (const std::string& s : skipped)
            list += (list.empty() ? "" : ", ") + s;
        log::warn("{}: skipped element records with no cell type: {}", kOp2Who, list);
    }
    if (grid_oriented > 0)
        log::warn(
            "{}: {} CBEAM element(s) are oriented by a local system; their X "
            "vector is read as zeros",
            kOp2Who, grid_oriented);
    if (conm1_oriented > 0)
        log::warn(
            "{}: {} CONM1 element(s) have their mass matrix in a coordinate system other than "
            "the basic one; it is read as given",
            kOp2Who, conm1_oriented);
    if (element_oriented > 0)
        log::warn(
            "{}: {} CBUSH element(s) are oriented by their element coordinate system; "
            "their X vector is read as zeros",
            kOp2Who, element_oriented);
    std::vector<detail::NastranCardRows> out;
    out.reserve(cards.size());
    rComps.clear();
    rBar.clear();
    for (auto& [card, rows] : cards) {
        out.push_back(std::move(rows));
        rComps.push_back(std::move(comps[card]));
        rBar.push_back(std::move(bars[card]));
    }
    return out;
}

// The property and material values the field data can set, as
// `nastran:<CARD>_<id>:<MEMBER>`: each member's word in the EPT/MPT record, default and
// whether it is an integer (the MID members).
struct Op2wValueMember {
    const char* mName;
    std::size_t mWord;
    double mDefault;
    bool mInt;
};
struct Op2wValueCard {
    const char* mCard;
    std::vector<Op2wValueMember> mMembers;
};
const std::vector<Op2wValueCard>& op2w_value_cards() {
    static const std::vector<Op2wValueCard> cards = {
        {"PSHELL", {{"MID1", 1, 1, true}, {"T", 2, 1.0, false}, {"NSM", 7, 0.0, false}}},
        {"PROD",
         {{"MID", 1, 1, true},
          {"A", 2, 1.0, false},
          {"J", 3, 0.0, false},
          {"C", 4, 0.0, false},
          {"NSM", 5, 0.0, false}}},
        {"PTUBE",
         {{"MID", 1, 1, true},
          {"OD", 2, 1.0, false},
          {"T", 3, 0.0, false},
          {"NSM", 4, 0.0, false}}},
        {"PSOLID", {{"MID", 1, 1, true}}},
        {"PMASS", {{"M", 1, 1.0, false}}},
        {"PCOMP",
         {{"Z0", 2, 0.0, false},
          {"NSM", 3, 0.0, false},
          {"SB", 4, 0.0, false},
          {"FT", 5, 0.0, true},
          {"TREF", 6, 0.0, false},
          {"GE", 7, 0.0, false},
          {"SYM", 1, 0.0, true}}},
        {"PBAR",
         {{"MID", 1, 1, true},
          {"A", 2, 1.0, false},
          {"I1", 3, 0.0, false},
          {"I2", 4, 0.0, false},
          {"J", 5, 0.0, false},
          {"NSM", 6, 0.0, false},
          {"C1", 8, 0.0, false},
          {"C2", 9, 0.0, false},
          {"D1", 10, 0.0, false},
          {"D2", 11, 0.0, false},
          {"E1", 12, 0.0, false},
          {"E2", 13, 0.0, false},
          {"F1", 14, 0.0, false},
          {"F2", 15, 0.0, false},
          {"K1", 16, 0.0, false},
          {"K2", 17, 0.0, false},
          {"I12", 18, 0.0, false}}},
        {"PELAS", {{"K", 1, 1.0, false}, {"GE", 2, 0.0, false}, {"S", 3, 0.0, false}}},
        {"PDAMP", {{"B", 1, 1.0, false}}},
        {"PVISC", {{"CE", 1, 1.0, false}, {"CR", 2, 0.0, false}}},
        {"PGAP",
         {{"U0", 1, 0.0, false},
          {"F0", 2, 0.0, false},
          {"KA", 3, 1.0e8, false},
          {"KB", 4, 1.0e-6, false},
          {"KT", 5, 0.0, false},
          {"MU1", 6, 0.0, false},
          {"MU2", 7, 0.0, false},
          {"TMAX", 8, 0.0, false},
          {"MAR", 9, 100.0, false},
          {"TRMIN", 10, 0.001, false}}},
        {"PSHEAR",
         {{"MID", 1, 1, true},
          {"T", 2, 1.0, false},
          {"NSM", 3, 0.0, false},
          {"F1", 4, 0.0, false},
          {"F2", 5, 0.0, false}}},
        {"MAT1",
         {{"E", 1, 1.0, false},
          {"G", 2, 0.0, false},
          {"NU", 3, 0.3, false},
          {"RHO", 4, 0.0, false},
          {"A", 5, 0.0, false},
          {"TREF", 6, 0.0, false},
          {"GE", 7, 0.0, false}}},
        // The anisotropic shell material: the stiffness matrix, the thermal expansion and the
        // stress limits (the MAT2 record is 17 words).
        {"MAT2",
         {{"G11", 1, 1.0, false},
          {"G12", 2, 0.0, false},
          {"G13", 3, 0.0, false},
          {"G22", 4, 1.0, false},
          {"G23", 5, 0.0, false},
          {"G33", 6, 1.0, false},
          {"RHO", 7, 0.0, false},
          {"A1", 8, 0.0, false},
          {"A2", 9, 0.0, false},
          {"A3", 10, 0.0, false},
          {"TREF", 11, 0.0, false},
          {"GE", 12, 0.0, false},
          {"ST", 13, 0.0, false},
          {"SC", 14, 0.0, false},
          {"SS", 15, 0.0, false},
          {"MCSID", 16, 0, true}}},
        // The thermal conductivity material (the MAT4 record is 11 words).
        {"MAT4",
         {{"K", 1, 1.0, false},
          {"CP", 2, 0.0, false},
          {"RHO", 3, 1.0, false},
          {"H", 4, 0.0, false},
          {"MU", 5, 0.0, false},
          {"HGEN", 6, 1.0, false},
          {"REFENTH", 7, 0.0, false},
          {"TCH", 8, 0.0, false},
          {"TDELTA", 9, 0.0, false},
          {"QLAT", 10, 0.0, false}}},
        // The anisotropic thermal conductivity material (the MAT5 record is 10 words).
        {"MAT5",
         {{"K1", 1, 1.0, false},
          {"K2", 2, 0.0, false},
          {"K3", 3, 0.0, false},
          {"K4", 4, 1.0, false},
          {"K5", 5, 0.0, false},
          {"K6", 6, 1.0, false},
          {"CP", 7, 0.0, false},
          {"RHO", 8, 1.0, false},
          {"HGEN", 9, 1.0, false}}},
        // The anisotropic solid material: the upper triangle of the stiffness matrix by rows,
        // the density, the thermal expansion, TREF and GE (the MAT9 record is 35 words).
        {"MAT9", {{"G11", 1, 1.0, false},  {"G12", 2, 0.0, false},   {"G13", 3, 0.0, false},
                  {"G14", 4, 0.0, false},  {"G15", 5, 0.0, false},   {"G16", 6, 0.0, false},
                  {"G22", 7, 1.0, false},  {"G23", 8, 0.0, false},   {"G24", 9, 0.0, false},
                  {"G25", 10, 0.0, false}, {"G26", 11, 0.0, false},  {"G33", 12, 1.0, false},
                  {"G34", 13, 0.0, false}, {"G35", 14, 0.0, false},  {"G36", 15, 0.0, false},
                  {"G44", 16, 1.0, false}, {"G45", 17, 0.0, false},  {"G46", 18, 0.0, false},
                  {"G55", 19, 1.0, false}, {"G56", 20, 0.0, false},  {"G66", 21, 1.0, false},
                  {"RHO", 22, 0.0, false}, {"A1", 23, 0.0, false},   {"A2", 24, 0.0, false},
                  {"A3", 25, 0.0, false},  {"A4", 26, 0.0, false},   {"A5", 27, 0.0, false},
                  {"A6", 28, 0.0, false},  {"TREF", 29, 0.0, false}, {"GE", 30, 0.0, false}}},
        // The orthotropic ply material (the MAT8 record is 19 words).
        {"MAT8",
         {{"E1", 1, 1.0, false},
          {"E2", 2, 1.0, false},
          {"NU12", 3, 0.3, false},
          {"G12", 4, 0.0, false},
          {"G1Z", 5, 0.0, false},
          {"G2Z", 6, 0.0, false},
          {"RHO", 7, 0.0, false},
          {"A1", 8, 0.0, false},
          {"A2", 9, 0.0, false},
          {"TREF", 10, 0.0, false},
          {"XT", 11, 0.0, false},
          {"XC", 12, 0.0, false},
          {"YT", 13, 0.0, false},
          {"YC", 14, 0.0, false},
          {"S", 15, 0.0, false},
          {"GE", 16, 0.0, false},
          {"F12", 17, 0.0, false},
          {"STRN", 18, 0.0, false}}},
    };
    return cards;
}

/// The members of the value card named @p pCard.
const Op2wValueCard& op2w_value_card(const char* pCard) {
    for (const Op2wValueCard& c : op2w_value_cards())
        if (std::string_view(c.mCard) == pCard)
            return c;
    throw std::logic_error(std::string("no value card ") + pCard);
}

/// The cross-section types of PBARL and PBEAML and their dimensions per section.
const std::map<std::string, std::size_t>& op2_section_dims() {
    static const std::map<std::string, std::size_t> dims = {
        {"ROD", 1},  {"TUBE", 2},  {"TUBE2", 2}, {"L", 4},     {"I", 6},  {"CHAN", 4},
        {"T", 4},    {"BOX", 4},   {"BAR", 2},   {"CROSS", 4}, {"H", 4},  {"T1", 4},
        {"I1", 4},   {"CHAN1", 4}, {"Z", 4},     {"CHAN2", 4}, {"T2", 4}, {"BOX1", 6},
        {"HEXA", 3}, {"HAT", 4},   {"HAT1", 5},  {"DBOX", 10}};
    return dims;
}

/// A text of up to eight characters as a `uint8` array, without its padding.
NDArray op2_text_array(const std::string& rText) {
    NDArray a(DType::UInt8, {rText.size()});
    std::copy(rText.begin(), rText.end(), a.As<std::uint8_t>());
    return a;
}

/// `nastran:<CARD>_<id>:<MEMBER>` field data of the PBARL, PBEAML and PCOMPG records. A PBARL or
/// PBEAML is the id, MID, GROUP and TYPE (two words each, text) and the values up to a -1: a
/// PBARL's dimensions then NSM, a PBEAML's stations as SO, X/XB, the dimensions and NSM. They
/// give `MID` (int64), `GROUP` and `TYPE` (`uint8` text) and `DIM` (`(n)`, or `(stations, n)`),
/// `NSM`, and a PBEAML's `SO` and `XXB`. A PCOMPG is pid, LAMOPT, Z0, NSM, SB, FT, TREF and GE,
/// then the plies as GPLYID, MID, T, THETA and SOUT up to five -1: `LAMOPT` and `FT` (int64),
/// the scalars, and `PLY`, `(n, 5)` rows. A file of 8-byte words gives none (its texts are
/// four characters a word). The first record of an id wins.
void op2_add_line_values(Mesh& rMesh, const Op2Stream& rS, const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    if (w.mWs != 4)
        return;
    const auto put_floats = [&](const std::string& rName, const std::vector<double>& rValues,
                                std::vector<std::size_t> Shape) {
        NDArray a(DType::Float64, std::move(Shape));
        std::copy(rValues.begin(), rValues.end(), a.As<double>());
        rMesh.AddFieldData(rName, std::move(a));
    };
    const auto put_int = [&](const std::string& rName, std::int64_t Value) {
        NDArray a(DType::Int64, {1});
        a.As<std::int64_t>()[0] = Value;
        rMesh.AddFieldData(rName, std::move(a));
    };
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        const bool pcompg = r.mK1 == 15006 && r.mK2 == 150;
        const bool pbarl = r.mK1 == 9102 && r.mK2 == 91;
        const bool pbeaml = r.mK1 == 9202 && r.mK2 == 92;
        if (!pcompg && !pbarl && !pbeaml)
            continue;
        const std::size_t body = 12;
        if (r.mRaw.size() <= body || (r.mRaw.size() - body) % 4)
            continue;
        const std::size_t total = (r.mRaw.size() - body) / 4;
        const char* base = r.mRaw.data() + body;
        const auto int_at = [&](std::size_t k) { return w.Int(base + k * 4); };
        const auto float_at = [&](std::size_t k) { return w.Float(base + k * 4); };
        std::size_t k = 0;
        if (pcompg) {
            while (k + 8 <= total && int_at(k) > 0) {
                const std::int64_t pid = int_at(k);
                const std::size_t start = k + 8;
                std::size_t at = start;
                const auto terminator = [&](std::size_t At) {
                    for (std::size_t j = 0; j < 5; ++j)
                        if (int_at(At + j) != -1)
                            return false;
                    return true;
                };
                while (at + 5 <= total && !terminator(at))
                    at += 5;
                if (at + 5 > total)
                    break;
                const std::string prefix = "nastran:PCOMPG_" + std::to_string(pid) + ":";
                if (!rMesh.HasFieldData(prefix + "PLY")) {
                    const std::size_t count = (at - start) / 5;
                    std::vector<double> plies(count * 5);
                    for (std::size_t j = 0; j < count; ++j) {
                        const std::size_t p = start + 5 * j;
                        plies[5 * j] = static_cast<double>(int_at(p));
                        plies[5 * j + 1] = static_cast<double>(int_at(p + 1));
                        plies[5 * j + 2] = float_at(p + 2);
                        plies[5 * j + 3] = float_at(p + 3);
                        plies[5 * j + 4] = static_cast<double>(int_at(p + 4));
                    }
                    put_floats(prefix + "PLY", plies, {count, 5});
                    put_int(prefix + "LAMOPT", int_at(k + 1));
                    put_int(prefix + "FT", int_at(k + 5));
                    for (const auto& [member, word] : {std::pair<const char*, std::size_t>{"Z0", 2},
                                                       {"NSM", 3},
                                                       {"SB", 4},
                                                       {"TREF", 6},
                                                       {"GE", 7}})
                        put_floats(prefix + member, {float_at(k + word)}, {1});
                }
                k = at + 5;
            }
            continue;
        }
        const std::string card = pbarl ? "PBARL" : "PBEAML";
        while (k + 6 <= total && int_at(k) > 0) {
            std::size_t end = k + 6;
            while (end < total && int_at(end) != -1)
                ++end;
            if (end >= total)
                break;
            const std::int64_t pid = int_at(k);
            const std::int64_t mid = int_at(k + 1);
            const std::string prefix = "nastran:" + card + "_" + std::to_string(pid) + ":";
            const auto text_at = [&](std::size_t At) {
                std::string t(base + At * 4, 8);
                const std::size_t last = t.find_last_not_of(std::string(" \0", 2));
                return last == std::string::npos ? std::string() : t.substr(0, last + 1);
            };
            const std::string group = text_at(k + 2);
            std::string section = text_at(k + 4);
            for (char& c : section)
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            std::vector<double> values;
            for (std::size_t j = k + 6; j < end; ++j)
                values.push_back(float_at(j));
            k = end + 1;
            if (rMesh.HasFieldData(prefix + "TYPE"))
                continue;
            if (pbarl) {
                if (values.empty())
                    continue;
                const double nsm = values.back();
                values.pop_back();
                put_floats(prefix + "DIM", values, {values.size()});
                put_floats(prefix + "NSM", {nsm}, {1});
            } else {
                const auto it = op2_section_dims().find(section);
                if (it == op2_section_dims().end() || values.empty() ||
                    values.size() % (it->second + 3) != 0)
                    continue;
                const std::size_t ndim = it->second;
                const std::size_t n = values.size() / (ndim + 3);
                std::vector<double> so(n), xxb(n), nsm(n), dim(n * ndim);
                for (std::size_t j = 0; j < n; ++j) {
                    const double* row = values.data() + j * (ndim + 3);
                    so[j] = row[0];
                    xxb[j] = row[1];
                    std::copy(row + 2, row + 2 + ndim, dim.begin() + static_cast<long>(j * ndim));
                    nsm[j] = row[ndim + 2];
                }
                put_floats(prefix + "SO", so, {n});
                put_floats(prefix + "XXB", xxb, {n});
                put_floats(prefix + "DIM", dim, {n, ndim});
                put_floats(prefix + "NSM", nsm, {n});
            }
            put_int(prefix + "MID", mid);
            rMesh.AddFieldData(prefix + "GROUP", op2_text_array(group));
            rMesh.AddFieldData(prefix + "TYPE", op2_text_array(section));
        }
    }
}

/// `nastran:<CARD>_<id>:<MEMBER>` field data of the PSHELL, PROD, PTUBE, PSOLID and MAT1
/// records (the members `op2w_value_cards` lists).
/// `nastran:PCOMP_<id>:<MEMBER>` field data of the PCOMP records: the scalars Z0, NSM, SB,
/// TREF and GE, the failure theory FT and SYM (1 for a symmetric layup), and PLY, the
/// plies as (n, 4) rows of MID, T, THETA and SOUT. The first record of an id wins.
void op2_add_pcomp_values(Mesh& rMesh, const Op2Stream& rS, const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    const auto put_float = [&](const std::string& rName, double Value) {
        if (rMesh.HasFieldData(rName))
            return;
        NDArray a(DType::Float64, {1});
        a.As<double>()[0] = Value;
        rMesh.AddFieldData(rName, std::move(a));
    };
    const auto put_int = [&](const std::string& rName, std::int64_t Value) {
        if (rMesh.HasFieldData(rName))
            return;
        NDArray a(DType::Int64, {1});
        a.As<std::int64_t>()[0] = Value;
        rMesh.AddFieldData(rName, std::move(a));
    };
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        if (r.mK1 != 2706 || r.mK2 != 27)
            continue;
        const std::size_t body = 3 * ws;
        if (r.mRaw.size() <= body || (r.mRaw.size() - body) % ws)
            continue;
        const std::size_t total = (r.mRaw.size() - body) / ws;
        const char* base = r.mRaw.data() + body;
        const auto int_at = [&](std::size_t k) { return w.Int(base + k * ws); };
        const auto float_at = [&](std::size_t k) { return w.Float(base + k * ws); };
        std::size_t k = 0;
        while (k + 8 <= total) {
            const std::int64_t pid = int_at(k);
            const std::int64_t nlayers = int_at(k + 1);
            const auto count = static_cast<std::size_t>(nlayers < 0 ? -nlayers : nlayers);
            if (pid <= 0 || count == 0 || k + 8 + 4 * count > total)
                break;
            const std::string prefix = "nastran:PCOMP_" + std::to_string(pid);
            put_float(prefix + ":Z0", float_at(k + 2));
            put_float(prefix + ":NSM", float_at(k + 3));
            put_float(prefix + ":SB", float_at(k + 4));
            put_int(prefix + ":FT", int_at(k + 5));
            put_float(prefix + ":TREF", float_at(k + 6));
            put_float(prefix + ":GE", float_at(k + 7));
            put_int(prefix + ":SYM", nlayers < 0 ? 1 : 0);
            if (!rMesh.HasFieldData(prefix + ":PLY")) {
                NDArray plies(DType::Float64, {count, 4});
                double* out = plies.As<double>();
                for (std::size_t j = 0; j < count; ++j) {
                    const std::size_t at = k + 8 + 4 * j;
                    out[4 * j] = static_cast<double>(int_at(at));
                    out[4 * j + 1] = float_at(at + 1);
                    out[4 * j + 2] = float_at(at + 2);
                    out[4 * j + 3] = static_cast<double>(int_at(at + 3));
                }
                rMesh.AddFieldData(prefix + ":PLY", std::move(plies));
            }
            k += 8 + 4 * count;
        }
    }
}

/// `nastran:PBUSH_<id>:<MEMBER>` field data of the PBUSH records: K and B as (6,) arrays, GE
/// as one value (the 18-word NX layout) or six (the 23-word MSC layout), and SA, ST, EA and ET
/// as scalars. A record of another size is not read. The first record of an id wins.
void op2_add_bush_values(Mesh& rMesh, const Op2Stream& rS, const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    static const char* const kScalar[] = {"SA", "ST", "EA", "ET"};
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        if (r.mK1 != 1402 || r.mK2 != 14)
            continue;
        const std::size_t body = 3 * ws;
        if (r.mRaw.size() != body + 18 * ws && r.mRaw.size() != body + 23 * ws)
            continue;
        const std::size_t words = (r.mRaw.size() - body) / ws;
        const char* base = r.mRaw.data() + body;
        const std::string prefix = "nastran:PBUSH_" + std::to_string(w.Int(base)) + ":";
        if (rMesh.HasFieldData(prefix + "K"))
            continue;
        const auto put_array = [&](const char* pMember, std::size_t First, std::size_t Count) {
            NDArray a(DType::Float64, {Count});
            double* out = a.As<double>();
            for (std::size_t k = 0; k < Count; ++k)
                out[k] = w.Float(base + (First + k) * ws);
            rMesh.AddFieldData(prefix + pMember, std::move(a));
        };
        put_array("K", 1, 6);
        put_array("B", 7, 6);
        put_array("GE", 13, words - 17);
        for (std::size_t k = 0; k < 4; ++k)
            put_array(kScalar[k], words - 4 + k, 1);
    }
}

/// The word of station @p I (of @p N) in a PBEAM record's eleven station slots: the first
/// station in slot 0, the last in slot 10, and any between them in the slots after the first.
std::size_t op2_pbeam_slot(std::size_t I, std::size_t N) {
    return I == 0 ? 0 : (I + 1 == N ? 10 : I);
}

/// `nastran:PBEAM_<id>:<MEMBER>` field data of the PBEAM records: MID, NSEG and CCF as
/// scalars, each station's SO, XXB, A, I1, I2, I12, J and NSM as (n,) arrays and the points
/// C1, C2, D1 and D2 as (n, 2) rows, for the n = NSEG + 1 stations. The first record of an
/// id wins.
void op2_add_pbeam_values(Mesh& rMesh, const Op2Stream& rS, const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    static const char* const kStation[] = {"SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM"};
    static const char* const kPoint[] = {"C1", "C2", "D1", "D2"};
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        if (r.mK1 != 5402 || r.mK2 != 54)
            continue;
        const std::size_t body = 3 * ws;
        if (r.mRaw.size() != body + 197 * ws)
            continue;
        const char* base = r.mRaw.data() + body;
        const std::int64_t pid = w.Int(base);
        const std::int64_t mid = w.Int(base + ws);
        const std::int64_t nseg = w.Int(base + 2 * ws);
        const std::int64_t ccf = w.Int(base + 3 * ws);
        if (nseg < 1 || nseg > 10)
            continue;
        const auto n = static_cast<std::size_t>(nseg) + 1;
        const std::string prefix = "nastran:PBEAM_" + std::to_string(pid) + ":";
        if (rMesh.HasFieldData(prefix + "MID"))
            continue;
        const auto put_int = [&](const std::string& rName, std::int64_t Value) {
            NDArray a(DType::Int64, {1});
            a.As<std::int64_t>()[0] = Value;
            rMesh.AddFieldData(rName, std::move(a));
        };
        put_int(prefix + "MID", mid);
        put_int(prefix + "NSEG", nseg);
        put_int(prefix + "CCF", ccf);
        const auto station = [&](std::size_t I, std::size_t Word) {
            return w.Float(base + (5 + 16 * op2_pbeam_slot(I, n) + Word) * ws);
        };
        for (std::size_t m = 0; m < 8; ++m) {
            NDArray a(DType::Float64, {n});
            double* out = a.As<double>();
            for (std::size_t i = 0; i < n; ++i)
                out[i] = station(i, m);
            rMesh.AddFieldData(prefix + kStation[m], std::move(a));
        }
        for (std::size_t p = 0; p < 4; ++p) {
            NDArray a(DType::Float64, {n, 2});
            double* out = a.As<double>();
            for (std::size_t i = 0; i < n; ++i) {
                out[2 * i] = station(i, 8 + 2 * p);
                out[2 * i + 1] = station(i, 9 + 2 * p);
            }
            rMesh.AddFieldData(prefix + kPoint[p], std::move(a));
        }
    }
}

void op2_add_property_values(Mesh& rMesh, const Op2Stream& rS,
                             const std::vector<Op2Table>& rTables) {
    const Op2Words& w = rS.Words();
    const auto ws = static_cast<std::size_t>(w.mWs);
    auto take = [&](const std::string& rCard, const std::string& rRaw, std::size_t Size) {
        if (rRaw.size() <= 3 * ws || (rRaw.size() - 3 * ws) % (Size * ws))
            return;
        const std::size_t rows = (rRaw.size() - 3 * ws) / (Size * ws);
        for (const Op2wValueCard& c : op2w_value_cards()) {
            if (rCard != c.mCard)
                continue;
            for (std::size_t r = 0; r < rows; ++r) {
                const char* row = rRaw.data() + 3 * ws + r * Size * ws;
                const std::int64_t ident = w.Int(row);
                if (rCard == "PSHELL" && w.Int(row + ws) >= 100000000)
                    continue;  // NX's PSHELL twin of a PCOMP
                for (const Op2wValueMember& m : c.mMembers) {
                    const std::string name =
                        "nastran:" + rCard + "_" + std::to_string(ident) + ":" + m.mName;
                    if (rMesh.HasFieldData(name))
                        continue;  // the first record wins
                    const char* at = row + m.mWord * ws;
                    NDArray a(m.mInt ? DType::Int64 : DType::Float64, {1});
                    if (m.mInt)
                        a.As<std::int64_t>()[0] = w.Int(at);
                    else
                        a.As<double>()[0] = w.Float(at);
                    rMesh.AddFieldData(name, std::move(a));
                }
            }
        }
    };
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"EPT"})) {
        for (const Op2PropertySpec& p : kOp2Properties)
            if (p.mK1 == r.mK1 && p.mK2 == r.mK2) {
                const std::string card = p.mCard;
                if (card == "PSHELL" || card == "PROD" || card == "PTUBE" || card == "PSOLID" ||
                    card == "PBAR" || card == "PELAS" || card == "PDAMP" || card == "PVISC" ||
                    card == "PGAP" || card == "PSHEAR" || card == "PMASS")
                    take(card, r.mRaw, static_cast<std::size_t>(p.mSizes[0]));
                break;
            }
    }
    for (const Op2CardRecord& r : op2_geometry_records(rS, rTables, {"MPT", "MAT"}))
        if (r.mK1 == 103 && r.mK2 == 1)
            take("MAT1", r.mRaw, 12);
        else if (r.mK1 == 203 && r.mK2 == 2)
            take("MAT2", r.mRaw, 17);
        else if (r.mK1 == 2103 && r.mK2 == 21)
            take("MAT4", r.mRaw, 11);
        else if (r.mK1 == 2203 && r.mK2 == 22)
            take("MAT5", r.mRaw, 10);
        else if (r.mK1 == 2503 && r.mK2 == 25)
            take("MAT8", r.mRaw, 19);
        else if (r.mK1 == 2603 && r.mK2 == 26)
            take("MAT9", r.mRaw, 35);
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

// The node blocks (the centre, then the corners) of a plate type that reports its
// corners and of a solid type; 0 for any other type.
std::int64_t op2_corner_nodes(std::int64_t Type) {
    switch (Type) {
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
}

std::int64_t op2_solid_nodes(std::int64_t Type) {
    switch (Type) {
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
    const std::int64_t cn = op2_corner_nodes(Type), sn = op2_solid_nodes(Type);
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

// GRID ids that are not 1 to n, as `nastran:nid` point data, from the model's GRID index
// (GRID id -> point index).
void op2_add_node_ids(Mesh& rMesh,
                      const std::unordered_map<std::int64_t, std::size_t>& rGridIndex) {
    const std::size_t n = rMesh.NumPoints();
    std::vector<std::int64_t> ids(n, 0);
    for (const auto& [id, index] : rGridIndex)
        if (index < n)
            ids[index] = id;
    for (std::size_t i = 0; i < n; ++i)
        if (ids[i] != static_cast<std::int64_t>(i + 1)) {
            rMesh.AddPointData("nastran:nid", op2_int_array(ids));
            return;
        }
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
    std::vector<std::int64_t> mCd;       // per point; empty for a deck-built model
    std::vector<std::int64_t> mCordIds;  // the id of every coordinate system card the file holds
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
    for (const detail::NastranCoordCard& c : g.mCords)
        m.mCordIds.push_back(c.mCid);
    m.mCd = g.mCd;
    const std::unordered_set<std::int64_t> grids(g.mIds.begin(), g.mIds.end());
    std::vector<std::vector<std::int64_t>> comps;
    std::vector<std::vector<double>> bars;
    const auto cards = op2_read_elements(rR.mStream, rR.mTables, grids, comps, bars);
    detail::NastranCellIds ids;
    const detail::NastranCells cells = detail::nastran_add_cells(
        m.mMesh, cards, m.mGridIndex, g.mScalarPoints, op2_read_properties(rR.mStream, rR.mTables),
        kOp2Who, comps, bars, &ids);
    if (!ids.mCard.empty()) {
        m.mMesh.AddCellData("nastran:card", std::move(ids.mCard));
        m.mMesh.AddCellData("nastran:c1", std::move(ids.mC1));
        m.mMesh.AddCellData("nastran:c2", std::move(ids.mC2));
    }
    if (ids.mHasBar) {
        m.mMesh.AddCellData("nastran:x", std::move(ids.mX));
        m.mMesh.AddCellData("nastran:g0", std::move(ids.mG0));
        m.mMesh.AddCellData("nastran:pin", std::move(ids.mPin));
        m.mMesh.AddCellData("nastran:w", std::move(ids.mW));
    }
    if (ids.mHasSpring) {
        m.mMesh.AddCellData("nastran:k", std::move(ids.mK));
        m.mMesh.AddCellData("nastran:ge", std::move(ids.mGe));
        m.mMesh.AddCellData("nastran:s", std::move(ids.mS));
    }
    if (ids.mHasDamper)
        m.mMesh.AddCellData("nastran:b", std::move(ids.mB));
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

// A 3-vector as a float64 array.
NDArray op2_vec3(const double* pV) {
    NDArray a(DType::Float64, {3});
    std::copy(pV, pV + 3, a.As<double>());
    return a;
}

// The coordinate systems the file defines, as `nastran:CORD2_<id>:*` field data: the type,
// the origin A and the axes X, Y and Z, all in basic coordinates. A write puts GRIDs back in
// them from these.
void op2_add_cord_systems(Mesh& rMesh, const detail::NastranCoordSystems& rSystems,
                          const std::vector<std::int64_t>& rIds) {
    for (const std::int64_t cid : std::set<std::int64_t>(rIds.begin(), rIds.end())) {
        if (cid == 0 || !rSystems.Has(cid))
            continue;
        const detail::NastranCoordSystems::System& s = rSystems.Get(cid);
        const std::string prefix = "nastran:CORD2_" + std::to_string(cid) + ":";
        rMesh.AddFieldData(prefix + "TYPE", op2_scalar(static_cast<double>(s.mType), DType::Int64));
        rMesh.AddFieldData(prefix + "A", op2_vec3(s.mOrigin));
        rMesh.AddFieldData(prefix + "X", op2_vec3(s.mAxes));
        rMesh.AddFieldData(prefix + "Y", op2_vec3(s.mAxes + 3));
        rMesh.AddFieldData(prefix + "Z", op2_vec3(s.mAxes + 6));
    }
}

}  // namespace

Mesh read_nastran_op2(const std::string& rPath, const ReadOptions& rOpts) {
    const Op2Reader r(rPath);
    Op2Model model = op2_build_mesh(r);
    Mesh& mesh = model.mMesh;
    // GRID ids that are not 1 to n, as point data: only when point data is wanted
    // (`points_only` drops it, and `arrays` names what to keep).
    if (rOpts.WantsAnyData() && rOpts.WantsArray("nastran:nid"))
        op2_add_node_ids(mesh, model.mGridIndex);
    op2_add_property_values(mesh, r.mStream, r.mTables);
    op2_add_pcomp_values(mesh, r.mStream, r.mTables);
    op2_add_pbeam_values(mesh, r.mStream, r.mTables);
    op2_add_bush_values(mesh, r.mStream, r.mTables);
    op2_add_line_values(mesh, r.mStream, r.mTables);
    op2_add_cord_systems(mesh, model.mSystems, model.mCordIds);
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
        // Columns are plies, stations or element nodes: a column index far
        // beyond the values read is a corrupt record, not a width to allocate.
        if (width > e.mValue.size() + 4096)
            op2_fail("a multi-valued element table names column " + std::to_string(width - 1));
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
    out.reserve(r.mSteps.size());
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

// ---------------------------------------------------------------------------------------
// The writer: 32-bit little-endian MSC OP2, SORT1 real. The twin of
// nastran_op2/_write.py, which writes the same bytes.
// ---------------------------------------------------------------------------------------

namespace {

constexpr detail::FortranRecordLayout kOp2wLayout{4, false};
constexpr std::int32_t kOp2wDevice = 1;
constexpr const char* kOp2wVersion = "XXXXXXXX";
constexpr const char* kOp2wTapeId = "NASTRAN FORT TAPE ID CODE - ";

[[noreturn]] void op2w_fail(const std::string& rMessage) {
    throw WriteError(kOp2Who + ": " + rMessage);
}

/// A double as the file's single precision, rounded like numpy's cast (an
/// overflow is an infinity).
float op2w_float(double V) {
    if (V != V)
        return std::numeric_limits<float>::quiet_NaN();
    if (std::fabs(V) >= 3.4028235677973366e38)
        return V < 0 ? -std::numeric_limits<float>::infinity()
                     : std::numeric_limits<float>::infinity();
    return static_cast<float>(V);
}

void op2w_i32(std::string& rBuf, std::int64_t V) {
    detail::fortran_put_int32(rBuf, static_cast<std::int32_t>(V), false);
}

void op2w_f32(std::string& rBuf, double V) {
    detail::fortran_put_float32(rBuf, op2w_float(V), false);
}

/// A row of @p Size words, filled by index.
struct Op2wRowBuf {
    explicit Op2wRowBuf(std::size_t Size) : mWords(Size, 0u) {}
    void I(std::size_t At, std::int64_t V) { mWords[At] = static_cast<std::uint32_t>(V); }
    void F(std::size_t At, double V) { mWords[At] = std::bit_cast<std::uint32_t>(op2w_float(V)); }
    void Append(std::string& rOut) const {
        for (const std::uint32_t w : mWords)
            op2w_i32(rOut, static_cast<std::int32_t>(w));
    }
    std::vector<std::uint32_t> mWords;
};

void op2w_marker(std::string& rOut, std::int32_t V) {
    std::string payload;
    detail::fortran_put_int32(payload, V, false);
    detail::fortran_append_record(rOut, payload.data(), payload.size(), kOp2wLayout);
}

void op2w_record(std::string& rOut, const std::string& rPayload) {
    op2w_marker(rOut, static_cast<std::int32_t>(rPayload.size() / 4));
    detail::fortran_append_record(rOut, rPayload.data(), rPayload.size(), kOp2wLayout);
}

std::string op2w_text(const std::string& rText) {
    std::string t = rText.substr(0, 8);
    t.resize(8, ' ');
    return t;
}

std::string op2w_key(std::int32_t K1, std::int32_t K2, std::int32_t K3) {
    std::string key;
    op2w_i32(key, K1);
    op2w_i32(key, K2);
    op2w_i32(key, K3);
    return key;
}

/// A named table: its name, the 7-word header, then every record behind a marker
/// triple, and the closing triple and 0.
std::string op2w_table(const std::string& rName, const std::vector<std::string>& rRecords) {
    std::string out;
    op2w_record(out, op2w_text(rName));
    op2w_marker(out, -1);
    std::string head;
    op2w_i32(head, 101);
    for (int i = 0; i < 6; ++i)
        op2w_i32(head, 0);
    op2w_record(out, head);
    std::int32_t k = -2;
    const auto triple = [&] {
        op2w_marker(out, k);
        op2w_marker(out, 1);
        op2w_marker(out, 0);
    };
    triple();
    op2w_record(out, op2w_text(rName));
    for (const std::string& r : rRecords) {
        --k;
        triple();
        op2w_record(out, r);
    }
    --k;
    triple();
    op2w_marker(out, 0);
    return out;
}

std::string op2w_file_header() {
    std::string out;
    op2w_marker(out, 3);
    std::string date;
    op2w_i32(date, 1);
    op2w_i32(date, 1);
    op2w_i32(date, 0);
    detail::fortran_append_record(out, date.data(), date.size(), kOp2wLayout);
    op2w_marker(out, 7);
    const std::string tape = kOp2wTapeId;
    detail::fortran_append_record(out, tape.data(), tape.size(), kOp2wLayout);
    op2w_record(out, op2w_text(kOp2wVersion));
    op2w_marker(out, -1);
    op2w_marker(out, 0);
    return out;
}

// -- the element cards ---------------------------------------------------------------

enum Op2wProp {
    kOp2wNoProp = -1,
    kOp2wProd = 0,
    kOp2wPshell,
    kOp2wPsolid,
    kOp2wPtube,
    kOp2wPelas,
    kOp2wPdamp,
    kOp2wPbar,
    kOp2wPcomp,
    kOp2wPbeam,
    kOp2wPbush,
    kOp2wPvisc,
    kOp2wPshear,
    kOp2wPgap,
    kOp2wPmass,
    kOp2wPbarl,
    kOp2wPbeaml,
    kOp2wPcompg
};
// The property cards' names, in the order of `Op2wProp`.
constexpr const char* kOp2wPropNames[] = {"PROD", "PSHELL", "PSOLID", "PTUBE",  "PELAS", "PDAMP",
                                          "PBAR", "PCOMP",  "PBEAM",  "PBUSH",  "PVISC", "PSHEAR",
                                          "PGAP", "PMASS",  "PBARL",  "PBEAML", "PCOMPG"};

constexpr int kOp2wHexa20[20] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int kOp2wPenta15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};

struct Op2wCard {
    const char* mName;
    std::int32_t mK1, mK2, mK3;
    std::size_t mSize;    // record words
    std::size_t mLinear;  // nodes of the linear card
    std::size_t mFirst;   // first node word
    const int* mPerm;     // to Nastran's order for the quadratic cell
    Op2wProp mProp;
    std::size_t mCorners;  // corner nodes with results (solids)
};

// In card-name order: the order of the GEOM2 records.
constexpr Op2wCard kOp2wCards[] = {
    {"CBAR", 2408, 24, 180, 16, 2, 2, nullptr, kOp2wPbar, 0},
    {"CBEAM", 5408, 54, 261, 18, 2, 2, nullptr, kOp2wPbeam, 0},
    {"CBUSH", 2608, 26, 60, 14, 2, 2, nullptr, kOp2wPbush, 0},
    {"CDAMP1", 201, 2, 69, 6, 2, 2, nullptr, kOp2wPdamp, 0},
    {"CDAMP2", 301, 3, 70, 6, 2, 2, nullptr, kOp2wNoProp, 0},
    {"CELAS1", 601, 6, 73, 6, 2, 2, nullptr, kOp2wPelas, 0},
    {"CELAS2", 701, 7, 74, 8, 2, 2, nullptr, kOp2wNoProp, 0},
    {"CGAP", 1908, 19, 104, 9, 2, 2, nullptr, kOp2wPgap, 0},
    {"CHEXA", 7308, 73, 253, 22, 8, 2, kOp2wHexa20, kOp2wPsolid, 8},
    {"CMASS1", 1001, 10, 65, 6, 2, 2, nullptr, kOp2wPmass, 0},
    {"CMASS2", 1101, 11, 66, 6, 2, 2, nullptr, kOp2wNoProp, 0},
    {"CONM1", 1401, 14, 63, 24, 1, 1, nullptr, kOp2wNoProp, 0},
    {"CONM2", 1501, 15, 64, 13, 1, 1, nullptr, kOp2wNoProp, 0},
    {"CONROD", 1601, 16, 47, 8, 2, 1, nullptr, kOp2wNoProp, 0},
    {"CPENTA", 4108, 41, 280, 17, 6, 2, kOp2wPenta15, kOp2wPsolid, 6},
    {"CPYRAM", 17200, 172, 1000, 16, 5, 2, nullptr, kOp2wPsolid, 5},
    {"CQUAD", 9108, 91, 507, 11, 4, 2, nullptr, kOp2wPshell, 0},
    {"CQUAD4", 2958, 51, 177, 14, 4, 2, nullptr, kOp2wPshell, 0},
    {"CQUAD8", 4701, 47, 326, 17, 4, 2, nullptr, kOp2wPshell, 0},
    {"CROD", 3001, 30, 48, 4, 2, 2, nullptr, kOp2wProd, 0},
    {"CSHEAR", 3101, 31, 61, 6, 4, 2, nullptr, kOp2wPshear, 0},
    {"CTETRA", 5508, 55, 217, 12, 4, 2, nullptr, kOp2wPsolid, 4},
    {"CTRIA3", 5959, 59, 282, 13, 3, 2, nullptr, kOp2wPshell, 0},
    {"CTRIA6", 4801, 48, 327, 14, 3, 2, nullptr, kOp2wPshell, 0},
    {"CTUBE", 3701, 37, 49, 4, 2, 2, nullptr, kOp2wPtube, 0},
    {"CVISC", 3901, 39, 50, 4, 2, 2, nullptr, kOp2wPvisc, 0},
    {"PLOTEL", 5201, 52, 11, 3, 2, 1, nullptr, kOp2wNoProp, 0},
};
constexpr std::size_t kOp2wNumCards = sizeof(kOp2wCards) / sizeof(kOp2wCards[0]);

std::size_t op2w_card_index(std::string_view Name) {
    for (std::size_t i = 0; i < kOp2wNumCards; ++i)
        if (Name == kOp2wCards[i].mName)
            return i;
    return kOp2wNumCards;
}

/// The card a line or vertex takes from its `nastran:card` number, or nullptr when the
/// writer does not write that card yet (CBEAM, CBUSH, ...).
const char* op2w_line_card(std::int64_t Code) {
    switch (Code) {
        case 30:
            return "CROD";
        case 16:
            return "CONROD";
        case 37:
            return "CTUBE";
        case 6:
            return "CELAS1";
        case 2:
            return "CDAMP1";
        case 15:
            return "CONM2";
        case 10:
            return "CMASS1";
        case 11:
            return "CMASS2";
        case 14:
            return "CONM1";
        case 24:
            return "CBAR";
        case 54:
            return "CBEAM";
        case 26:
            return "CBUSH";
        case 3:
            return "CDAMP2";
        case 7:
            return "CELAS2";
        case 39:
            return "CVISC";
        case 52:
            return "PLOTEL";
        case 19:
            return "CGAP";
        default:
            return nullptr;
    }
}

/// Whether a cell block of @p rType may hold an element of card @p Card.
bool op2w_block_holds(const std::string& rType, std::string_view Card) {
    if (rType == "vertex")
        return Card == "CONM2" || Card == "CELAS1" || Card == "CDAMP1" || Card == "CELAS2" ||
               Card == "CDAMP2" || Card == "CMASS1" || Card == "CMASS2" || Card == "CONM1";
    if (rType == "line")
        return Card == "CROD" || Card == "CONROD" || Card == "CTUBE" || Card == "CELAS1" ||
               Card == "CDAMP1" || Card == "CBAR" || Card == "CBEAM" || Card == "CBUSH" ||
               Card == "CELAS2" || Card == "CDAMP2" || Card == "CVISC" || Card == "PLOTEL" ||
               Card == "CGAP" || Card == "CMASS1" || Card == "CMASS2";
    return false;
}

/// The card a `nastran:card` number names (its GEOM2 type), for the warnings.
std::string op2w_card_name(std::int64_t Code) {
    static const std::pair<const char*, std::int64_t> kTypes[] = {
        {"CBAR", 24},    {"CBEAM", 54},  {"CBUSH", 26},  {"CDAMP1", 2},  {"CDAMP2", 3},
        {"CELAS1", 6},   {"CELAS2", 7},  {"CGAP", 19},   {"CHEXA", 73},  {"CMASS1", 10},
        {"CMASS2", 11},  {"CONM1", 14},  {"CONM2", 15},  {"CONROD", 16}, {"CPENTA", 41},
        {"CPYRAM", 172}, {"CQUAD", 91},  {"CQUAD4", 51}, {"CQUAD8", 47}, {"CQUADR", 80},
        {"CROD", 30},    {"CSHEAR", 31}, {"CTETRA", 55}, {"CTRIA3", 59}, {"CTRIA6", 48},
        {"CTRIAR", 92},  {"CTUBE", 37},  {"CVISC", 39},  {"PLOTEL", 52},
    };
    for (const auto& [name, number] : kTypes)
        if (number == Code)
            return name;
    return "card " + std::to_string(Code);
}

struct Op2wTypeSpec {
    const char* mType;
    const char* mCard;  // null: a line takes CROD, CTUBE or CONROD by its property
    std::size_t mNodes;
    std::int64_t mPid;  // the property id of a block with no 'nastran:pid'
};
constexpr Op2wTypeSpec kOp2wTypes[] = {
    {"vertex", "CONM2", 1, -1},    {"line", nullptr, 2, 3},        {"triangle", "CTRIA3", 3, 1},
    {"triangle6", "CTRIA6", 6, 1}, {"quad", "CQUAD4", 4, 1},       {"quad8", "CQUAD8", 8, 1},
    {"quad9", "CQUAD", 9, 1},      {"tetra", "CTETRA", 4, 2},      {"tetra10", "CTETRA", 10, 2},
    {"pyramid", "CPYRAM", 5, 2},   {"pyramid13", "CPYRAM", 13, 2}, {"wedge", "CPENTA", 6, 2},
    {"wedge15", "CPENTA", 15, 2},  {"hexahedron", "CHEXA", 8, 2},  {"hexahedron20", "CHEXA", 20, 2},
};

const Op2wTypeSpec* op2w_type_spec(const std::string& rType) {
    for (const Op2wTypeSpec& t : kOp2wTypes)
        if (rType == t.mType)
            return &t;
    return nullptr;
}

/// The columns of a row's card fields: a CONM1's 21 mass terms are the widest.
constexpr std::size_t kOp2wBarWidth = 21;

/// The mesh block (and cells) behind one element row.
struct Op2wBlock {
    std::size_t mIndex = 0;  // the mesh's block
    std::size_t mCount = 0;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;  // 0-based, mCount x mNodes
    std::vector<std::int64_t> mEids;
    std::vector<std::int64_t> mPids;
    std::vector<std::size_t> mCards;   // card index per cell
    std::vector<std::int64_t> mComps;  // two per cell: a spring's or damper's components
    std::vector<double> mBar;          // twelve per cell: a CBAR's fields, a CELAS2's K, GE
                                       // and S, or a CDAMP2's B (the card's first columns)
};

struct Op2wRow {
    std::int64_t mEid;
    std::size_t mBlock;  // into Op2wModel::mBlocks
    std::size_t mRow;
};

/// What every step of a series shares: the points' count and the cell blocks'
/// types and shapes.
using Op2wSignature = std::vector<std::tuple<std::string, std::size_t, std::size_t>>;

Op2wSignature op2w_signature(const Mesh& rMesh) {
    Op2wSignature sig;
    sig.emplace_back("points", rMesh.NumPoints(), 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        sig.emplace_back(cb.Type(), cb.NumCells(), cb.NodesPerCell());
    }
    return sig;
}

/// The digits of "<PTYPE>_<digits>" when the prefix is @p rPtype, else -1.
std::int64_t op2w_ptype_pid(const std::string& rName, const std::string& rPtype) {
    if (rName.size() <= rPtype.size() + 1 || rName.compare(0, rPtype.size(), rPtype) != 0 ||
        rName[rPtype.size()] != '_')
        return -1;
    std::int64_t pid = 0;
    for (std::size_t i = rPtype.size() + 1; i < rName.size(); ++i) {
        if (rName[i] < '0' || rName[i] > '9' || pid > 100000000000000000)
            return -1;
        pid = pid * 10 + (rName[i] - '0');
    }
    return pid;
}

/// The GRID id of each point: `point_data["nastran:nid"]`, else 1 to n.
std::vector<std::int64_t> op2w_node_ids(const Mesh& rMesh, std::size_t N) {
    std::vector<std::int64_t> ids(N);
    if (!rMesh.HasPointData("nastran:nid")) {
        for (std::size_t i = 0; i < N; ++i)
            ids[i] = static_cast<std::int64_t>(i + 1);
        return ids;
    }
    const NDArray& a = rMesh.PointData("nastran:nid");
    if (a.Size() != N || detail::is_float_dtype(a.Dtype()))
        op2w_fail("'nastran:nid' must hold one integer per point");
    const detail::Int64View v(a);
    ids.assign(v.Data(), v.Data() + N);
    std::vector<std::int64_t> sorted = ids;
    std::sort(sorted.begin(), sorted.end());
    if (!sorted.empty() && sorted.front() < 1)
        op2w_fail("GRID ids must be positive");
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        op2w_fail("GRID ids are not unique");
    return ids;
}

using Op2wValues = std::map<std::pair<std::string, std::int64_t>, std::map<std::string, double>>;

/// `{(card, id): {member: value}}` from the `nastran:<CARD>_<id>:<MEMBER>` field data; keys
/// that name no card or member here are left alone.
/// The plies of each PCOMP id: `nastran:PCOMP_<id>:PLY`, rows of MID, T, THETA and SOUT.
using Op2wPlies = std::map<std::int64_t, std::vector<std::array<double, 4>>>;

Op2wPlies op2w_plies(const Mesh& rMesh) {
    Op2wPlies out;
    const std::string prefix = "nastran:PCOMP_";
    const std::string suffix = ":PLY";
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.size() <= prefix.size() + suffix.size() || name.rfind(prefix, 0) != 0 ||
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        const std::string digits =
            name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
        if (digits.empty() || digits.size() > 15 ||
            digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        const NDArray& a = rMesh.FieldData(name);
        const std::vector<std::size_t>& shape = a.Shape();
        if (shape.size() != 2 || shape[1] != 4 || shape[0] == 0)
            op2w_fail("'" + name + "' must hold (plies, 4) rows of MID, T, THETA and SOUT");
        const detail::DoubleView view(a);
        auto& plies = out[std::stoll(digits)];
        for (std::size_t r = 0; r < shape[0]; ++r) {
            const double* row = view.Data() + 4 * r;
            plies.push_back({row[0], row[1], row[2], row[3]});
        }
    }
    return out;
}

Op2wValues op2w_property_values(const Mesh& rMesh) {
    Op2wValues out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind("nastran:", 0) != 0)
            continue;
        const std::size_t colon = name.rfind(':');
        const std::size_t under = name.rfind('_', colon);
        if (colon <= 8 || under == std::string::npos || under < 8 || under + 1 >= colon)
            continue;
        const std::string card = name.substr(8, under - 8);
        const std::string member = name.substr(colon + 1);
        const std::string digits = name.substr(under + 1, colon - under - 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        for (const Op2wValueCard& c : op2w_value_cards()) {
            if (card != c.mCard)
                continue;
            for (const Op2wValueMember& m : c.mMembers)
                if (member == m.mName) {
                    const NDArray& a = rMesh.FieldData(name);
                    if (a.Size() != 1)
                        op2w_fail("'" + name + "' must hold one value");
                    out[{card, std::stoll(digits)}][member] =
                        m.mInt ? static_cast<double>(detail::Int64View(a).Data()[0])
                               : detail::DoubleView(a).Data()[0];
                }
        }
    }
    return out;
}

/// `nastran:<CARD>_<id>:<MEMBER>` field data of the PBEAM or PBUSH records, as values per member.
using Op2wArrays = std::map<std::int64_t, std::map<std::string, std::vector<double>>>;

std::vector<double> op2w_as_doubles(const NDArray& rA) {
    std::vector<double> out(rA.Size());
    if (rA.Dtype() == DType::Int64) {
        const detail::Int64View view(rA);
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = static_cast<double>(view.Data()[i]);
    } else {
        const detail::DoubleView view(rA);
        std::copy(view.Data(), view.Data() + out.size(), out.begin());
    }
    return out;
}

Op2wArrays op2w_named_arrays(const Mesh& rMesh, const std::string& rCard) {
    const std::string kPrefix = "nastran:" + rCard + "_";
    Op2wArrays out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind(kPrefix, 0) != 0)
            continue;
        const std::size_t colon = name.find(':', kPrefix.size());
        if (colon == std::string::npos || colon == kPrefix.size())
            continue;
        const std::string digits = name.substr(kPrefix.size(), colon - kPrefix.size());
        const std::string member = name.substr(colon + 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        if (member.empty() ||
            member.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") != std::string::npos)
            continue;
        out[std::stoll(digits)][member] = op2w_as_doubles(rMesh.FieldData(name));
    }
    return out;
}

/// A PBEAM's scalar (MID, NSEG or CCF) from the beams' field data, or @p Default.
double op2w_array_scalar(const Op2wArrays& rArrays, std::int64_t Pid, const char* pMember,
                         double Default) {
    const auto beam = rArrays.find(Pid);
    if (beam == rArrays.end())
        return Default;
    const auto it = beam->second.find(pMember);
    if (it == beam->second.end() || it->second.empty())
        return Default;
    return it->second[0];
}

/// A station member the field data leaves out: the defaults of a unit section, and XXB
/// spaced evenly from 0 to 1 (as numpy's linspace spaces it).
std::vector<double> op2w_beam_default(const std::string& rMember, std::size_t N) {
    if (rMember == "XXB") {
        std::vector<double> out(N);
        const double step = 1.0 / static_cast<double>(N - 1);
        for (std::size_t i = 0; i < N; ++i)
            out[i] = static_cast<double>(i) * step;
        out.back() = 1.0;
        return out;
    }
    const bool unit =
        rMember == "SO" || rMember == "A" || rMember == "I1" || rMember == "I2" || rMember == "J";
    return std::vector<double>(N, unit ? 1.0 : 0.0);
}

/// The 197-word PBEAM record of one id: the header (pid, MID, segments, CCF and a zero), then
/// the stations in the slots `op2_pbeam_slot` names, from the id's field data.
std::string op2w_pbeam_row(std::int64_t Pid, const Op2wArrays& rArrays) {
    const auto mid = static_cast<std::int64_t>(op2w_array_scalar(rArrays, Pid, "MID", 1));
    const auto nseg = static_cast<std::int64_t>(op2w_array_scalar(rArrays, Pid, "NSEG", 1));
    const auto ccf = static_cast<std::int64_t>(op2w_array_scalar(rArrays, Pid, "CCF", 1));
    if (nseg < 1 || nseg > 10)
        op2w_fail("PBEAM " + std::to_string(Pid) + " has " + std::to_string(nseg) +
                  " segments; a PBEAM holds 1 to 10");
    const auto n = static_cast<std::size_t>(nseg) + 1;
    static const char* const kStation[] = {"SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM"};
    static const char* const kPoint[] = {"C1", "C2", "D1", "D2"};
    static const std::map<std::string, std::vector<double>> kNone;
    const auto beam = rArrays.find(Pid);
    const auto& given = beam == rArrays.end() ? kNone : beam->second;
    const auto given_values = [&](const char* pName, std::size_t Count,
                                  std::vector<double> Default) {
        const auto it = given.find(pName);
        if (it == given.end())
            return Default;
        if (it->second.size() != Count)
            op2w_fail("'nastran:PBEAM_" + std::to_string(Pid) + ":" + pName + "' must have " +
                      std::to_string(Count) + " value(s)");
        return it->second;
    };
    std::vector<std::vector<double>> stations;
    for (const char* pName : kStation)
        stations.push_back(given_values(pName, n, op2w_beam_default(pName, n)));
    std::vector<std::vector<double>> points;
    for (const char* pName : kPoint)
        points.push_back(given_values(pName, 2 * n, std::vector<double>(2 * n, 0.0)));
    Op2wRowBuf w(197);
    w.I(0, Pid);
    w.I(1, mid);
    w.I(2, nseg);
    w.I(3, ccf);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t at = 5 + 16 * op2_pbeam_slot(i, n);
        for (std::size_t m = 0; m < 8; ++m)
            w.F(at + m, stations[m][i]);
        for (std::size_t p = 0; p < 4; ++p) {
            w.F(at + 8 + 2 * p, points[p][2 * i]);
            w.F(at + 9 + 2 * p, points[p][2 * i + 1]);
        }
    }
    std::string out;
    w.Append(out);
    return out;
}

/// The 23-word PBUSH record of one id: the id, K, B and GE (six each) and SA, ST, EA and ET,
/// from the id's field data. An id with none gets K of 1, B of 0 and GE of 0, and stress and
/// strain coefficients of 1. A K, B or GE given as one value is written to all six.
std::string op2w_pbush_row(std::int64_t Pid, const Op2wArrays& rArrays) {
    static const std::map<std::string, std::vector<double>> kNone;
    const auto bush = rArrays.find(Pid);
    const auto& given = bush == rArrays.end() ? kNone : bush->second;
    Op2wRowBuf w(23);
    w.I(0, Pid);
    const auto six = [&](const char* pMember, double Default, std::size_t First) {
        const auto it = given.find(pMember);
        std::vector<double> values = it == given.end() ? std::vector<double>{Default} : it->second;
        if (values.size() == 1) {
            // Copy the value first: assign() would read it from a buffer it may free.
            const double value = values[0];
            values.assign(6, value);
        }
        if (values.size() != 6)
            op2w_fail("'nastran:PBUSH_" + std::to_string(Pid) + ":" + pMember +
                      "' must have 1 or 6 value(s)");
        for (std::size_t k = 0; k < 6; ++k)
            w.F(First + k, values[k]);
    };
    six("K", 1.0, 1);
    six("B", 0.0, 7);
    six("GE", 0.0, 13);
    const char* const kScalar[] = {"SA", "ST", "EA", "ET"};
    for (std::size_t k = 0; k < 4; ++k) {
        const auto it = given.find(kScalar[k]);
        const std::vector<double> value = it == given.end() ? std::vector<double>{1.0} : it->second;
        if (value.size() != 1)
            op2w_fail(std::string("'nastran:PBUSH_") + std::to_string(Pid) + ":" + kScalar[k] +
                      "' must have 1 value");
        w.F(19 + k, value[0]);
    }
    std::string out;
    w.Append(out);
    return out;
}

/// A coordinate system a GRID's CP names, as the writer needs it: its type (1 rectangular,
/// 2 cylindrical, 3 spherical), its origin A and its axes X, Y and Z, all in basic
/// coordinates, from the `nastran:CORD2_<id>:*` field data.
struct Op2wCord {
    std::int64_t mCid = 0;
    int mType = 1;
    double mA[3] = {0, 0, 0};
    double mX[3] = {1, 0, 0};
    double mY[3] = {0, 1, 0};
    double mZ[3] = {0, 0, 1};
};

// Degrees to radians: the reader's constant, so both engines round the same way.
constexpr double kOp2wDegree = 3.14159265358979323846 / 180.0;

/// The basic point @p pBasic in the local coordinates of system @p rC: the inverse of the
/// reader's transform, in the operation order of the Python twin.
void op2w_to_local(const Op2wCord& rC, const double* pBasic, double* pOut) {
    const double d0 = pBasic[0] - rC.mA[0];
    const double d1 = pBasic[1] - rC.mA[1];
    const double d2 = pBasic[2] - rC.mA[2];
    const double xl = d0 * rC.mX[0] + d1 * rC.mX[1] + d2 * rC.mX[2];
    const double yl = d0 * rC.mY[0] + d1 * rC.mY[1] + d2 * rC.mY[2];
    const double zl = d0 * rC.mZ[0] + d1 * rC.mZ[1] + d2 * rC.mZ[2];
    if (rC.mType == 1) {
        pOut[0] = xl;
        pOut[1] = yl;
        pOut[2] = zl;
    } else if (rC.mType == 2) {
        pOut[0] = std::sqrt(xl * xl + yl * yl);
        pOut[1] = std::atan2(yl, xl) / kOp2wDegree;
        pOut[2] = zl;
    } else {
        const double r = std::sqrt(xl * xl + yl * yl + zl * zl);
        pOut[0] = r;
        pOut[1] = r > 0.0 ? std::acos(std::min(1.0, std::max(-1.0, zl / r))) / kOp2wDegree : 0.0;
        pOut[2] = std::atan2(yl, xl) / kOp2wDegree;
    }
}

/// The basic vector @p pV at the basic point @p pPoint in the components of system @p rC: the
/// inverse of the reader's rotation, in the operation order of the Python twin.
void op2w_vector_to_local(const Op2wCord& rC, const double* pPoint, const double* pV,
                          double* pOut) {
    const double w0 = rC.mX[0] * pV[0] + rC.mX[1] * pV[1] + rC.mX[2] * pV[2];
    const double w1 = rC.mY[0] * pV[0] + rC.mY[1] * pV[1] + rC.mY[2] * pV[2];
    const double w2 = rC.mZ[0] * pV[0] + rC.mZ[1] * pV[1] + rC.mZ[2] * pV[2];
    if (rC.mType == 1) {
        pOut[0] = w0;
        pOut[1] = w1;
        pOut[2] = w2;
        return;
    }
    const double d0 = pPoint[0] - rC.mA[0];
    const double d1 = pPoint[1] - rC.mA[1];
    const double d2 = pPoint[2] - rC.mA[2];
    const double l0 = rC.mX[0] * d0 + rC.mX[1] * d1 + rC.mX[2] * d2;
    const double l1 = rC.mY[0] * d0 + rC.mY[1] * d1 + rC.mY[2] * d2;
    const double l2 = rC.mZ[0] * d0 + rC.mZ[1] * d1 + rC.mZ[2] * d2;
    const double ph = std::atan2(l1, l0);
    const double cp = std::cos(ph);
    const double sp = std::sin(ph);
    if (rC.mType == 2) {
        pOut[0] = w0 * cp + w1 * sp;
        pOut[1] = -w0 * sp + w1 * cp;
        pOut[2] = w2;
        return;
    }
    const double th = std::atan2(std::sqrt(l0 * l0 + l1 * l1), l2);
    const double ct = std::cos(th);
    const double st = std::sin(th);
    const double er[3] = {st * cp, st * sp, ct};
    const double et[3] = {ct * cp, ct * sp, -st};
    const double ep[3] = {-sp, cp, 0.0};
    pOut[0] = w0 * er[0] + w1 * er[1] + w2 * er[2];
    pOut[1] = w0 * et[0] + w1 * et[1] + w2 * et[2];
    pOut[2] = w0 * ep[0] + w1 * ep[1] + w2 * ep[2];
}

/// System @p Cid's type, origin and axes, from its `nastran:CORD2_<id>:*` field data.
Op2wCord op2w_cord(const Mesh& rMesh, std::int64_t Cid) {
    const std::string p = "nastran:CORD2_" + std::to_string(Cid) + ":";
    for (const char* member : {"TYPE", "A", "X", "Y", "Z"})
        if (!rMesh.HasFieldData(p + member))
            op2w_fail("GRIDs are in coordinate system " + std::to_string(Cid) + ", which has no '" +
                      p + "*' field data");
    Op2wCord c;
    c.mCid = Cid;
    const NDArray& t = rMesh.FieldData(p + "TYPE");
    if (t.Size() != 1)
        op2w_fail("'" + p + "TYPE' must have one value");
    const detail::DoubleView kind(t);
    if (kind.Data()[0] != 1.0 && kind.Data()[0] != 2.0 && kind.Data()[0] != 3.0)
        op2w_fail("'" + p + "TYPE' must be 1, 2 or 3");
    c.mType = static_cast<int>(kind.Data()[0]);
    const auto vector = [&](const char* member, double* pOut) {
        const NDArray& a = rMesh.FieldData(p + member);
        if (a.Size() != 3)
            op2w_fail("'" + p + member + "' must have 3 values");
        const detail::DoubleView view(a);
        std::copy(view.Data(), view.Data() + 3, pOut);
    };
    vector("A", c.mA);
    vector("X", c.mX);
    vector("Y", c.mY);
    vector("Z", c.mZ);
    return c;
}

/// Each point's coordinates in its own system (basic for CP 0), each point's CP and CD
/// (the system its results are in) and the systems they name, by id: what GEOM1 holds.
void op2w_local_points(const Mesh& rMesh, std::size_t N, const std::vector<double>& rBasic,
                       std::vector<std::int64_t>& rCp, std::vector<std::int64_t>& rCd,
                       std::vector<Op2wCord>& rCords, std::vector<double>& rLocal) {
    rCp.assign(N, 0);
    if (rMesh.HasPointData("nastran:cp")) {
        const NDArray& a = rMesh.PointData("nastran:cp");
        if (a.Size() != N)
            op2w_fail("'nastran:cp' must hold one integer per point");
        const detail::Int64View view(a);
        rCp.assign(view.Data(), view.Data() + N);
    }
    rCd.assign(N, 0);
    if (rMesh.HasPointData("nastran:cd")) {
        const NDArray& a = rMesh.PointData("nastran:cd");
        if (a.Size() != N)
            op2w_fail("'nastran:cd' must hold one integer per point");
        const detail::Int64View view(a);
        rCd.assign(view.Data(), view.Data() + N);
    }
    std::set<std::int64_t> ids;
    for (const std::int64_t cid : rCp) {
        if (cid < 0)
            op2w_fail("coordinate system ids must be positive");
        if (cid != 0)
            ids.insert(cid);
    }
    for (const std::int64_t cid : rCd) {
        if (cid < 0)
            op2w_fail("coordinate system ids must be positive");
        if (cid != 0)
            ids.insert(cid);
    }
    rCords.clear();
    std::map<std::int64_t, std::size_t> at;
    for (const std::int64_t cid : ids) {
        at[cid] = rCords.size();
        rCords.push_back(op2w_cord(rMesh, cid));
    }
    rLocal = rBasic;
    for (std::size_t i = 0; i < N; ++i)
        if (rCp[i] != 0)
            op2w_to_local(rCords[at[rCp[i]]], rBasic.data() + 3 * i, rLocal.data() + 3 * i);
}

/// The characters of a text member of the field data, without its padding.
std::string op2w_text_of(const std::vector<double>& rChars) {
    std::string t;
    for (const double c : rChars)
        t.push_back(static_cast<char>(static_cast<int>(c)));
    const std::size_t first = t.find_first_not_of(std::string(" \0", 2));
    if (first == std::string::npos)
        return {};
    return t.substr(first, t.find_last_not_of(std::string(" \0", 2)) - first + 1);
}

/// The id, MID, GROUP and TYPE of a PBARL or PBEAML (six words), and the section's number of
/// dimensions.
std::size_t op2w_section_header(const char* pCard, std::int64_t Pid,
                                const std::map<std::string, std::vector<double>>& rGiven,
                                Op2wRowBuf& rW) {
    std::string section = op2w_text_of(rGiven.at("TYPE"));
    for (char& c : section)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    const auto it = op2_section_dims().find(section);
    if (it == op2_section_dims().end())
        op2w_fail(std::string(pCard) + " " + std::to_string(Pid) + ": the section type '" +
                  section + "' is not one of the cards'");
    std::string group = "MSCBML0";
    if (const auto g = rGiven.find("GROUP"); g != rGiven.end())
        group = op2w_text_of(g->second);
    rW.I(0, Pid);
    const auto mid = rGiven.find("MID");
    rW.I(1, mid != rGiven.end() && !mid->second.empty() ? static_cast<std::int64_t>(mid->second[0])
                                                        : 1);
    const auto text_word = [](const std::string& rText, std::size_t K) {
        std::string t = rText.substr(0, 8);
        t.resize(8, ' ');
        std::uint32_t w = 0;
        std::memcpy(&w, t.data() + 4 * K, 4);
        return w;
    };
    rW.mWords[2] = text_word(group, 0);
    rW.mWords[3] = text_word(group, 1);
    rW.mWords[4] = text_word(section, 0);
    rW.mWords[5] = text_word(section, 1);
    return it->second;
}

/// A PBARL record: the header, the section's dimensions, NSM and a -1.
std::string op2w_pbarl_row(std::int64_t Pid,
                           const std::map<std::string, std::vector<double>>& rGiven) {
    Op2wRowBuf head(6);
    const std::size_t ndim = op2w_section_header("PBARL", Pid, rGiven, head);
    const auto dim = rGiven.find("DIM");
    const std::size_t have = dim == rGiven.end() ? 0 : dim->second.size();
    if (have != ndim)
        op2w_fail("'nastran:PBARL_" + std::to_string(Pid) + ":DIM' must have " +
                  std::to_string(ndim) + " value(s) for a '" + op2w_text_of(rGiven.at("TYPE")) +
                  "' section");
    Op2wRowBuf body(ndim + 2);
    for (std::size_t k = 0; k < ndim; ++k)
        body.F(k, dim->second[k]);
    const auto nsm = rGiven.find("NSM");
    body.F(ndim, nsm != rGiven.end() && !nsm->second.empty() ? nsm->second[0] : 0.0);
    body.I(ndim + 1, -1);
    std::string out;
    head.Append(out);
    body.Append(out);
    return out;
}

/// A PBEAML record: the header, per station SO, X/XB, the section's dimensions and NSM, and a -1.
std::string op2w_pbeaml_row(std::int64_t Pid,
                            const std::map<std::string, std::vector<double>>& rGiven) {
    Op2wRowBuf head(6);
    const std::size_t ndim = op2w_section_header("PBEAML", Pid, rGiven, head);
    const auto dim = rGiven.find("DIM");
    const std::size_t have = dim == rGiven.end() ? 0 : dim->second.size();
    if (have == 0 || have % ndim != 0)
        op2w_fail("'nastran:PBEAML_" + std::to_string(Pid) + ":DIM' must hold (stations, " +
                  std::to_string(ndim) + ") rows for a '" + op2w_text_of(rGiven.at("TYPE")) +
                  "' section");
    const std::size_t n = have / ndim;
    const auto column = [&](const char* pMember, std::vector<double> Default) {
        const auto it = rGiven.find(pMember);
        if (it == rGiven.end())
            return Default;
        if (it->second.size() != n)
            op2w_fail("'nastran:PBEAML_" + std::to_string(Pid) + ":" + pMember + "' must have " +
                      std::to_string(n) + " value(s)");
        return it->second;
    };
    std::vector<double> xxb(n, 0.0);
    if (n > 1)
        xxb = op2w_beam_default("XXB", n);
    const std::vector<double> so = column("SO", std::vector<double>(n, 0.0));
    const std::vector<double> x = column("XXB", xxb);
    const std::vector<double> nsm = column("NSM", std::vector<double>(n, 0.0));
    Op2wRowBuf body(n * (ndim + 3) + 1);
    for (std::size_t j = 0; j < n; ++j) {
        const std::size_t at = j * (ndim + 3);
        body.F(at, so[j]);
        body.F(at + 1, x[j]);
        for (std::size_t k = 0; k < ndim; ++k)
            body.F(at + 2 + k, dim->second[j * ndim + k]);
        body.F(at + 2 + ndim, nsm[j]);
    }
    body.I(n * (ndim + 3), -1);
    std::string out;
    head.Append(out);
    body.Append(out);
    return out;
}

/// A PCOMPG record: the header (pid, LAMOPT, Z0, NSM, SB, FT, TREF, GE), the plies (GPLYID, MID,
/// T, THETA, SOUT) and five -1.
std::string op2w_pcompg_row(std::int64_t Pid,
                            const std::map<std::string, std::vector<double>>& rGiven) {
    const std::vector<double>& flat = rGiven.at("PLY");
    if (flat.empty() || flat.size() % 5 != 0)
        op2w_fail("'nastran:PCOMPG_" + std::to_string(Pid) +
                  ":PLY' must hold (plies, 5) rows of GPLYID, MID, T, THETA and SOUT");
    const std::size_t count = flat.size() / 5;
    const auto value = [&](const char* pMember) {
        const auto it = rGiven.find(pMember);
        return it != rGiven.end() && !it->second.empty() ? it->second[0] : 0.0;
    };
    Op2wRowBuf w(8 + 5 * count + 5);
    w.I(0, Pid);
    w.I(1, static_cast<std::int64_t>(value("LAMOPT")));
    w.F(2, value("Z0"));
    w.F(3, value("NSM"));
    w.F(4, value("SB"));
    w.I(5, static_cast<std::int64_t>(value("FT")));
    w.F(6, value("TREF"));
    w.F(7, value("GE"));
    for (std::size_t j = 0; j < count; ++j) {
        const std::size_t at = 8 + 5 * j;
        w.I(at, static_cast<std::int64_t>(flat[5 * j]));
        w.I(at + 1, static_cast<std::int64_t>(flat[5 * j + 1]));
        w.F(at + 2, flat[5 * j + 2]);
        w.F(at + 3, flat[5 * j + 3]);
        w.I(at + 4, static_cast<std::int64_t>(flat[5 * j + 4]));
    }
    for (std::size_t j = 0; j < 5; ++j)
        w.I(8 + 5 * count + j, -1);
    std::string out;
    w.Append(out);
    return out;
}

class Op2wModel {
public:
    explicit Op2wModel(const Mesh& rMesh) {
        mSignature = op2w_signature(rMesh);
        mValues = op2w_property_values(rMesh);
        mPlies = op2w_plies(rMesh);
        mArrays = op2w_named_arrays(rMesh, "PBEAM");
        mBushes = op2w_named_arrays(rMesh, "PBUSH");
        mPbarl = op2w_named_arrays(rMesh, "PBARL");
        mPbeaml = op2w_named_arrays(rMesh, "PBEAML");
        mPcompg = op2w_named_arrays(rMesh, "PCOMPG");
        mNpoints = rMesh.NumPoints();
        mNodeIds = op2w_node_ids(rMesh, mNpoints);
        {
            const NDArray& pts = rMesh.Points();
            const std::size_t dim = rMesh.PointDim();
            const detail::DoubleView view(pts);
            mPoints.assign(mNpoints * 3, 0.0);
            for (std::size_t i = 0; i < mNpoints; ++i)
                for (std::size_t k = 0; k < 3 && k < dim; ++k)
                    mPoints[3 * i + k] = view.Data()[i * dim + k];
        }
        op2w_local_points(rMesh, mNpoints, mPoints, mCp, mCd, mCords, mLocal);
        for (std::size_t k = 0; k < mCords.size(); ++k)
            mCordAt[mCords[k].mCid] = k;
        std::set<std::int64_t> ptube;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const std::int64_t pid = op2w_ptype_pid(rMesh.Region(i).mName, "PTUBE");
            if (pid >= 0)
                ptube.insert(pid);
        }
        const bool have_eid = rMesh.HasCellData("nastran:eid");
        const bool have_pid = rMesh.HasCellData("nastran:pid");
        const bool have_card = rMesh.HasCellData("nastran:card");
        std::int64_t sequential = 1;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const Op2wTypeSpec* spec = op2w_type_spec(cb.Type());
            if (spec == nullptr) {
                if (cb.NumCells() > 0)
                    mSkipped.insert(cb.Type());
                continue;
            }
            Op2wBlock blk;
            blk.mIndex = b;
            blk.mCount = cb.NumCells();
            blk.mNodes = cb.NodesPerCell();
            if (blk.mCount > 0 && blk.mNodes != spec->mNodes)
                op2w_fail("a '" + cb.Type() + "' block needs " + std::to_string(spec->mNodes) +
                          " nodes per cell");
            const detail::Int64View conn(cb.Conn());
            blk.mConn.assign(conn.Data(), conn.Data() + blk.mCount * blk.mNodes);
            for (const std::int64_t node : blk.mConn)
                if (node < 0 || static_cast<std::size_t>(node) >= mNpoints)
                    op2w_fail("a '" + cb.Type() +
                              "' cell refers to a point the mesh does not have");
            blk.mEids.resize(blk.mCount);
            blk.mPids.resize(blk.mCount);
            if (have_eid) {
                const NDArray& a = rMesh.CellData("nastran:eid", b);
                if (a.Size() != blk.mCount)
                    op2w_fail("'nastran:eid' must have one value per cell");
                const detail::Int64View view(a);
                blk.mEids.assign(view.Data(), view.Data() + blk.mCount);
            } else {
                for (std::size_t i = 0; i < blk.mCount; ++i)
                    blk.mEids[i] = sequential + static_cast<std::int64_t>(i);
            }
            sequential += static_cast<std::int64_t>(blk.mCount);
            if (have_pid) {
                const NDArray& a = rMesh.CellData("nastran:pid", b);
                if (a.Size() != blk.mCount)
                    op2w_fail("'nastran:pid' must have one value per cell");
                const detail::Int64View view(a);
                blk.mPids.assign(view.Data(), view.Data() + blk.mCount);
            } else {
                blk.mPids.assign(blk.mCount, spec->mPid);
            }
            blk.mCards.resize(blk.mCount);
            // A line or vertex takes the card its `nastran:card` number names; a card the
            // writer does not write yet leaves its elements out.
            const bool line_or_vertex = cb.Type() == "line" || cb.Type() == "vertex";
            std::vector<std::int64_t> codes(blk.mCount, 0);
            if (have_card) {
                const NDArray& a = rMesh.CellData("nastran:card", b);
                if (a.Size() != blk.mCount)
                    op2w_fail("'nastran:card' must have one value per cell");
                const detail::Int64View view(a);
                codes.assign(view.Data(), view.Data() + blk.mCount);
            }
            std::vector<bool> keep(blk.mCount, true);
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                if (have_card && line_or_vertex) {
                    const char* name = op2w_line_card(codes[i]);
                    if (name == nullptr) {
                        keep[i] = false;
                        mUnwritten.insert(op2w_card_name(codes[i]));
                    } else if (!op2w_block_holds(cb.Type(), name)) {
                        op2w_fail("a '" + cb.Type() + "' block cannot hold " + name + " cards");
                    } else {
                        blk.mCards[i] = op2w_card_index(name);
                    }
                } else if (have_card && codes[i] == 31) {
                    // A quad is a CSHEAR; no other cell type holds one.
                    if (cb.Type() == "quad") {
                        blk.mCards[i] = op2w_card_index("CSHEAR");
                    } else {
                        keep[i] = false;
                        mUnwritten.insert("CSHEAR");
                    }
                } else if (spec->mCard != nullptr) {
                    blk.mCards[i] = op2w_card_index(spec->mCard);
                } else {
                    const std::int64_t pid = blk.mPids[i];
                    blk.mCards[i] = op2w_card_index(pid <= 0                ? "CONROD"
                                                    : ptube.count(pid) != 0 ? "CTUBE"
                                                                            : "CROD");
                }
            }
            const std::size_t cbar = op2w_card_index("CBAR");
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cbar && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    op2w_fail("a CBAR needs its orientation: nastran:x, or nastran:g0 for a G0");
            const std::size_t cbeam = op2w_card_index("CBEAM");
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cbeam && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    op2w_fail("a CBEAM needs its orientation: nastran:x, or nastran:g0 for a G0");
            const std::size_t cbush = op2w_card_index("CBUSH");
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cbush && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    op2w_fail("a CBUSH needs its orientation: nastran:x, or nastran:g0 for a G0");
            const std::size_t cgap = op2w_card_index("CGAP");
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cgap && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    op2w_fail("a CGAP needs its orientation: nastran:x, or nastran:g0 for a G0");
            const std::size_t celas2 = op2w_card_index("CELAS2");
            const std::size_t cdamp2 = op2w_card_index("CDAMP2");
            const std::size_t cmass2 = op2w_card_index("CMASS2");
            const std::size_t conm1 = op2w_card_index("CONM1");
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                if (keep[i] && blk.mCards[i] == celas2 && !rMesh.HasCellData("nastran:k"))
                    op2w_fail("a CELAS2 needs its value: nastran:k");
                if (keep[i] && blk.mCards[i] == cdamp2 && !rMesh.HasCellData("nastran:b"))
                    op2w_fail("a CDAMP2 needs its value: nastran:b");
                if (keep[i] && blk.mCards[i] == cmass2 && !rMesh.HasCellData("nastran:m"))
                    op2w_fail("a CMASS2 needs its value: nastran:m");
                if (keep[i] && blk.mCards[i] == conm1 && !rMesh.HasCellData("nastran:conm1"))
                    op2w_fail("a CONM1 needs its value: nastran:conm1");
            }
            std::size_t kept = 0;
            for (const bool k : keep)
                kept += k ? 1 : 0;
            if (kept == 0)
                continue;
            // A spring's or damper's components: `nastran:c1` and `nastran:c2`, else 1 at
            // each end. A vertex's second end is ground, so 0.
            std::vector<std::int64_t> c1(blk.mCount, 1);
            std::vector<std::int64_t> c2(blk.mCount, 1);
            auto read_component = [&](const char* pName, std::vector<std::int64_t>& rValues) {
                if (!rMesh.HasCellData(pName))
                    return;
                const NDArray& a = rMesh.CellData(pName, b);
                if (a.Size() != blk.mCount)
                    op2w_fail(std::string("'") + pName + "' must have one value per cell");
                const detail::Int64View view(a);
                rValues.assign(view.Data(), view.Data() + blk.mCount);
            };
            read_component("nastran:c1", c1);
            read_component("nastran:c2", c2);
            // The card fields of each row, from the cell data of its own card only: a CBAR's
            // X1..X3, G0, PA and PB, W1A..W3A and W1B..W3B; a CELAS2's K, GE and S; a CDAMP2's B.
            // A CELAS2's K and a CDAMP2's B share column 0, so a row takes only its own card's.
            blk.mBar.assign(blk.mCount * kOp2wBarWidth, 0.0);
            const auto read_field = [&](const char* pName, std::initializer_list<std::size_t> Cards,
                                        std::size_t Lo, std::size_t Width) {
                if (!rMesh.HasCellData(pName))
                    return;
                const NDArray& a = rMesh.CellData(pName, b);
                if (a.Size() != blk.mCount * Width)
                    op2w_fail(std::string("'") + pName + "' must have " + std::to_string(Width) +
                              " value(s) per cell");
                const detail::DoubleView view(a);
                for (std::size_t i = 0; i < blk.mCount; ++i)
                    if (std::find(Cards.begin(), Cards.end(), blk.mCards[i]) != Cards.end())
                        for (std::size_t k = 0; k < Width; ++k)
                            blk.mBar[i * kOp2wBarWidth + Lo + k] = view.Data()[i * Width + k];
            };
            // A CBAR and a CBEAM share the orientation, G0, pin and offset fields.
            read_field("nastran:x", {cbar, cbeam, cbush, cgap}, 0, 3);
            read_field("nastran:g0", {cbar, cbeam, cbush, cgap}, 3, 1);
            read_field("nastran:pin", {cbar, cbeam}, 4, 2);
            read_field("nastran:w", {cbar, cbeam}, 6, 6);
            read_field("nastran:k", {celas2}, 0, 1);
            read_field("nastran:ge", {celas2}, 1, 1);
            read_field("nastran:s", {celas2}, 2, 1);
            read_field("nastran:b", {cdamp2}, 0, 1);
            read_field("nastran:m", {cmass2}, 0, 1);
            read_field("nastran:conm1", {conm1}, 0, 21);
            blk.mComps.resize(2 * blk.mCount);
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                blk.mComps[2 * i] = c1[i];
                blk.mComps[2 * i + 1] = blk.mNodes == 1 ? 0 : c2[i];
            }
            if (kept < blk.mCount) {
                auto compact = [&keep](auto& rValues, std::size_t Stride) {
                    std::size_t to = 0;
                    for (std::size_t i = 0; i < keep.size(); ++i)
                        if (keep[i])
                            for (std::size_t s = 0; s < Stride; ++s)
                                rValues[to++] = rValues[i * Stride + s];
                    rValues.resize(to);
                };
                compact(blk.mConn, blk.mNodes);
                compact(blk.mEids, 1);
                compact(blk.mPids, 1);
                compact(blk.mCards, 1);
                compact(blk.mComps, 2);
                compact(blk.mBar, kOp2wBarWidth);
                blk.mCount = kept;
            }
            mBlocks.push_back(std::move(blk));
        }
        std::vector<std::int64_t> every;
        for (const Op2wBlock& blk : mBlocks) {
            for (const std::int64_t e : blk.mEids)
                if (e <= 0)
                    op2w_fail("element ids must be positive");
            every.insert(every.end(), blk.mEids.begin(), blk.mEids.end());
        }
        std::sort(every.begin(), every.end());
        if (std::adjacent_find(every.begin(), every.end()) != every.end())
            op2w_fail("element ids are not unique");
        // Cards in name order, each one's rows by element id.
        mRows.resize(kOp2wNumCards);
        for (std::size_t bi = 0; bi < mBlocks.size(); ++bi)
            for (std::size_t r = 0; r < mBlocks[bi].mCount; ++r)
                mRows[mBlocks[bi].mCards[r]].push_back({mBlocks[bi].mEids[r], bi, r});
        for (auto& rows : mRows)
            std::stable_sort(rows.begin(), rows.end(), [](const Op2wRow& rA, const Op2wRow& rB) {
                return rA.mEid < rB.mEid;
            });
        // Property ids by family.
        for (const Op2wBlock& blk : mBlocks)
            for (std::size_t r = 0; r < blk.mCount; ++r) {
                const Op2wProp prop = kOp2wCards[blk.mCards[r]].mProp;
                if (prop != kOp2wNoProp)
                    mProps[prop].insert(blk.mPids[r]);
            }
        const char* const* const names = kOp2wPropNames;
        // A shell whose property has plies is a PCOMP, not a PSHELL.
        for (auto it = mProps[kOp2wPshell].begin(); it != mProps[kOp2wPshell].end();) {
            if (mPlies.count(*it) != 0) {
                mProps[kOp2wPcomp].insert(*it);
                it = mProps[kOp2wPshell].erase(it);
            } else {
                ++it;
            }
        }
        // A bar, beam or shell whose id has the line-property field data is a PBARL, PBEAML or
        // PCOMPG, not a PBAR, PBEAM, PCOMP or PSHELL.
        const auto replace = [&](const Op2wArrays& rArrays, const char* pKey, Op2wProp To,
                                 std::initializer_list<Op2wProp> From) {
            for (const auto& [pid, members] : rArrays) {
                if (members.count(pKey) == 0)
                    continue;
                bool any = false;
                for (const Op2wProp f : From)
                    any = any || mProps[f].count(pid) != 0;
                if (!any)
                    continue;
                for (const Op2wProp f : From)
                    mProps[f].erase(pid);
                mProps[To].insert(pid);
            }
        };
        replace(mPbarl, "TYPE", kOp2wPbarl, {kOp2wPbar});
        replace(mPbeaml, "TYPE", kOp2wPbeaml, {kOp2wPbeam});
        replace(mPcompg, "PLY", kOp2wPcompg, {kOp2wPshell, kOp2wPcomp});
        std::map<std::int64_t, int> used;
        for (const int prop :
             {kOp2wPshell, kOp2wPsolid, kOp2wProd, kOp2wPtube, kOp2wPbar, kOp2wPelas, kOp2wPdamp,
              kOp2wPcomp, kOp2wPbeam, kOp2wPbush, kOp2wPvisc, kOp2wPshear, kOp2wPgap, kOp2wPmass,
              kOp2wPbarl, kOp2wPbeaml, kOp2wPcompg})
            for (const std::int64_t p : mProps[prop]) {
                const auto it = used.find(p);
                if (it != used.end())
                    op2w_fail("property id " + std::to_string(p) + " is used by " +
                              names[it->second] + " and " + names[prop] +
                              " elements: shells, solids and rods need ids of their own");
                used[p] = prop;
            }
        for (const auto& [pid, prop] : used)
            if (pid <= 0)
                op2w_fail("property ids must be positive");
        MaterialCards();  // refuses a material id given as two cards
    }

    const Op2wSignature& Signature() const { return mSignature; }

    /// The `(points, 3)` basic vectors in @p rValues, rotated in place into each point's
    /// output system CD.
    void ToOutput(std::vector<double>& rValues) const {
        for (std::size_t i = 0; i < mNpoints; ++i)
            if (mCd[i] != 0) {
                double out[3];
                op2w_vector_to_local(mCords[mCordAt.at(mCd[i])], mPoints.data() + 3 * i,
                                     rValues.data() + 3 * i, out);
                std::copy(out, out + 3, rValues.data() + 3 * i);
            }
    }
    bool HasOutputSystems() const {
        return std::any_of(mCd.begin(), mCd.end(), [](std::int64_t C) { return C != 0; });
    }
    const std::vector<std::int64_t>& Cd() const { return mCd; }
    const Op2wCord& OutputSystem(std::size_t Point) const { return mCords[mCordAt.at(mCd[Point])]; }
    const std::vector<double>& BasicPoints() const { return mPoints; }
    std::size_t NumPoints() const { return mNpoints; }
    const std::vector<std::int64_t>& NodeIds() const { return mNodeIds; }
    const std::vector<Op2wBlock>& Blocks() const { return mBlocks; }
    const std::vector<Op2wRow>& Rows(std::size_t Card) const { return mRows[Card]; }
    const std::set<std::string>& Skipped() const { return mSkipped; }
    const std::set<std::string>& Unwritten() const { return mUnwritten; }
    bool HasProperties() const {
        for (const auto& s : mProps)
            if (!s.empty())
                return true;
        return false;
    }

    std::string Geom1() const {
        // The coordinate systems the GRIDs are in: CORD2C, then CORD2R, then CORD2S, by id.
        static const struct {
            int mType;
            std::int32_t mK1, mK2, mK3;
        } kCords[] = {{2, 2001, 20, 9}, {1, 2101, 21, 8}, {3, 2201, 22, 10}};
        std::vector<std::string> records;
        for (const auto& k : kCords) {
            std::string rows;
            for (const Op2wCord& c : mCords) {
                if (c.mType != k.mType)
                    continue;
                Op2wRowBuf w(13);
                w.I(0, c.mCid);
                w.I(1, c.mType);
                w.I(2, 2);
                w.I(3, 0);  // RID 0: the points are in basic coordinates
                for (std::size_t j = 0; j < 3; ++j) {
                    w.F(4 + j, c.mA[j]);
                    w.F(7 + j, c.mA[j] + c.mZ[j]);   // B: a point on +z
                    w.F(10 + j, c.mA[j] + c.mX[j]);  // C: a point in the xz plane
                }
                w.Append(rows);
            }
            if (!rows.empty())
                records.push_back(op2w_key(k.mK1, k.mK2, k.mK3) + rows);
        }
        std::string rows;
        std::vector<std::size_t> order(mNpoints);  // GEOM1 holds the GRIDs by id
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(), [this](std::size_t A, std::size_t B) {
            return mNodeIds[A] < mNodeIds[B];
        });
        for (const std::size_t i : order) {
            Op2wRowBuf w(8);
            w.I(0, mNodeIds[i]);
            w.I(1, mCp[i]);  // CP: the system the coordinates are in
            w.I(5, mCd[i]);  // CD: the system the results are in
            for (std::size_t k = 0; k < 3; ++k)
                w.F(2 + k, mLocal[3 * i + k]);
            w.Append(rows);
        }
        records.push_back(op2w_key(4501, 45, 1) + rows);
        records.push_back(op2w_key(65535, 65535, 65535));
        return op2w_table("GEOM1", records);
    }

    std::string Geom2() const {
        std::vector<std::string> records;
        for (std::size_t c = 0; c < kOp2wNumCards; ++c) {
            if (mRows[c].empty())
                continue;
            const Op2wCard& card = kOp2wCards[c];
            std::string rows;
            for (const Op2wRow& row : mRows[c]) {
                const Op2wBlock& blk = mBlocks[row.mBlock];
                Op2wRowBuf w(card.mSize);
                w.I(0, row.mEid);
                std::vector<std::int64_t> nodes(blk.mNodes);
                const std::int64_t* conn = blk.mConn.data() + row.mRow * blk.mNodes;
                for (std::size_t k = 0; k < blk.mNodes; ++k) {
                    const std::size_t at =
                        (blk.mNodes > card.mLinear && card.mPerm != nullptr) ? card.mPerm[k] : k;
                    nodes[at] = mNodeIds[static_cast<std::size_t>(conn[k])];
                }
                if (std::string_view(card.mName) == "CONM1") {
                    // GRID, the coordinate system CID (0) and the 21 mass terms
                    w.I(1, nodes[0]);
                    for (std::size_t k = 0; k < 21; ++k)
                        w.F(3 + k, blk.mBar[row.mRow * kOp2wBarWidth + k]);
                } else if (std::string_view(card.mName) == "CONM2") {
                    w.I(1, nodes[0]);
                } else if (std::string_view(card.mName) == "PLOTEL") {
                    w.I(1, nodes[0]);  // no property: the two nodes at words 1 and 2
                    w.I(2, nodes[1]);
                } else if (std::string_view(card.mName) == "CONROD") {
                    w.I(1, nodes[0]);
                    w.I(2, nodes[1]);
                    w.I(3, 1);    // MID
                    w.F(4, 1.0);  // A
                } else {
                    w.I(1, blk.mPids[row.mRow]);
                    for (std::size_t k = 0; k < nodes.size(); ++k)
                        w.I(2 + k, nodes[k]);
                    if (std::string_view(card.mName) == "CBAR") {
                        // X1..X3 at 4 to 6, or the G0 at 4 with the flag 2; the flag at 7;
                        // PA and PB at 8 and 9; the offsets at 10 to 15.
                        const double* f = blk.mBar.data() + row.mRow * kOp2wBarWidth;
                        const auto g0 = static_cast<std::int64_t>(f[3]);
                        if (g0 != 0) {
                            w.I(4, g0);
                        } else {
                            w.F(4, f[0]);
                            w.F(5, f[1]);
                            w.F(6, f[2]);
                        }
                        w.I(7, g0 != 0 ? 2 : 0);
                        w.I(8, static_cast<std::int64_t>(f[4]));
                        w.I(9, static_cast<std::int64_t>(f[5]));
                        for (std::size_t k = 0; k < 6; ++k)
                            w.F(10 + k, f[6 + k]);
                    }
                    if (std::string_view(card.mName) == "CBEAM") {
                        // X1..X3 at 6 to 8 (flag 1: global) or the G0 at 6 with zero 7 and 8
                        // (flag 2), the flag at 9, PA and PB at 10 and 11, the offsets at 12 to
                        // 17; the scalar ends SA and SB stay zero.
                        const double* f = blk.mBar.data() + row.mRow * kOp2wBarWidth;
                        const auto g0 = static_cast<std::int64_t>(f[3]);
                        if (g0 != 0) {
                            w.I(6, g0);
                        } else {
                            for (std::size_t k = 0; k < 3; ++k)
                                w.F(6 + k, f[k]);
                        }
                        w.I(9, g0 != 0 ? 2 : 1);
                        w.I(10, static_cast<std::int64_t>(f[4]));
                        w.I(11, static_cast<std::int64_t>(f[5]));
                        for (std::size_t k = 0; k < 6; ++k)
                            w.F(12 + k, f[6 + k]);
                    }
                    if (std::string_view(card.mName) == "CBUSH") {
                        // X1..X3 at 4 to 6, or the G0 at 4 (flag 2); the flag at 7; the element
                        // system 0 and the stress coefficient 0.5; output system -1.
                        const double* f = blk.mBar.data() + row.mRow * kOp2wBarWidth;
                        const auto g0 = static_cast<std::int64_t>(f[3]);
                        if (g0 != 0) {
                            w.I(4, g0);
                        } else {
                            w.F(4, f[0]);
                            w.F(5, f[1]);
                            w.F(6, f[2]);
                        }
                        w.I(7, g0 != 0 ? 2 : 0);
                        w.F(9, 0.5);
                        w.I(10, -1);
                    }
                    if (std::string_view(card.mName) == "CGAP") {
                        // X1..X3 at 4 to 6, or the G0 at 4 (flag 2); the flag at 7 (1 for a
                        // vector); the coordinate system CID (8) stays 0.
                        const double* f = blk.mBar.data() + row.mRow * kOp2wBarWidth;
                        const auto g0 = static_cast<std::int64_t>(f[3]);
                        if (g0 != 0) {
                            w.I(4, g0);
                        } else {
                            w.F(4, f[0]);
                            w.F(5, f[1]);
                            w.F(6, f[2]);
                        }
                        w.I(7, g0 != 0 ? 2 : 1);
                    }
                    if (std::string_view(card.mName) == "CELAS1" ||
                        std::string_view(card.mName) == "CDAMP1" ||
                        std::string_view(card.mName) == "CELAS2" ||
                        std::string_view(card.mName) == "CDAMP2" ||
                        std::string_view(card.mName) == "CMASS1" ||
                        std::string_view(card.mName) == "CMASS2") {
                        w.I(4, blk.mComps[2 * row.mRow]);
                        w.I(5, blk.mComps[2 * row.mRow + 1]);
                    }
                    if (std::string_view(card.mName) == "CELAS2") {
                        // K at 1, in place of the property id; GE at 6 and S at 7.
                        const double* f = blk.mBar.data() + row.mRow * kOp2wBarWidth;
                        w.F(1, f[0]);
                        w.F(6, f[1]);
                        w.F(7, f[2]);
                    } else if (std::string_view(card.mName) == "CDAMP2" ||
                               std::string_view(card.mName) == "CMASS2") {
                        // B (a damper's) or M (a mass's), in place of the property id
                        w.F(1, blk.mBar[row.mRow * kOp2wBarWidth]);
                    }
                }
                w.Append(rows);
            }
            records.push_back(op2w_key(card.mK1, card.mK2, card.mK3) + rows);
        }
        records.push_back(op2w_key(65535, 65535, 65535));
        return op2w_table("GEOM2", records);
    }

    std::string Ept() const {
        struct Prop {
            Op2wProp mProp;
            const char* mCard;
            std::int32_t mK1, mK2, mK3;
            std::size_t mSize;
        };
        static const Prop props[] = {
            {kOp2wProd, "PROD", 902, 9, 29, 6},        {kOp2wPshell, "PSHELL", 2302, 23, 283, 11},
            {kOp2wPsolid, "PSOLID", 2402, 24, 281, 7}, {kOp2wPtube, "PTUBE", 1602, 16, 30, 5},
            {kOp2wPbar, "PBAR", 52, 20, 181, 19},      {kOp2wPelas, "PELAS", 302, 3, 46, 4},
            {kOp2wPdamp, "PDAMP", 202, 2, 45, 2},      {kOp2wPvisc, "PVISC", 1802, 18, 31, 3},
            {kOp2wPgap, "PGAP", 2102, 21, 121, 11},    {kOp2wPshear, "PSHEAR", 1002, 10, 42, 6},
            {kOp2wPmass, "PMASS", 402, 4, 44, 2}};
        std::vector<std::string> records;
        for (const Prop& p : props) {
            if (mProps[p.mProp].empty())
                continue;
            std::string rows;
            for (const std::int64_t pid : mProps[p.mProp]) {
                Op2wRowBuf w(p.mSize);
                w.I(0, pid);
                const auto given = mValues.find({p.mCard, pid});
                for (const Op2wValueCard& c : op2w_value_cards()) {
                    if (std::string_view(c.mCard) != p.mCard)
                        continue;
                    for (const Op2wValueMember& m : c.mMembers) {
                        double v = m.mDefault;
                        if (given != mValues.end()) {
                            const auto it = given->second.find(m.mName);
                            if (it != given->second.end())
                                v = it->second;
                        }
                        if (m.mInt)
                            w.I(m.mWord, static_cast<std::int64_t>(v));
                        else
                            w.F(m.mWord, v);
                    }
                }
                if (p.mProp == kOp2wPsolid)
                    w.I(6, 0x43454D53);  // FCTN 'SMEC'
                w.Append(rows);
            }
            records.push_back(op2w_key(p.mK1, p.mK2, p.mK3) + rows);
        }
        if (!mProps[kOp2wPcomp].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPcomp]) {
                const auto given = mValues.find({"PCOMP", pid});
                const auto value = [&](const char* pMember, double Default) {
                    if (given != mValues.end()) {
                        const auto it = given->second.find(pMember);
                        if (it != given->second.end())
                            return it->second;
                    }
                    return Default;
                };
                const auto& plies = mPlies.at(pid);
                const auto count = static_cast<std::int64_t>(plies.size());
                op2w_i32(rows, pid);
                op2w_i32(rows, value("SYM", 0) != 0 ? -count : count);
                op2w_f32(rows, value("Z0", 0.0));
                op2w_f32(rows, value("NSM", 0.0));
                op2w_f32(rows, value("SB", 0.0));
                op2w_i32(rows, static_cast<std::int64_t>(value("FT", 0.0)));
                op2w_f32(rows, value("TREF", 0.0));
                op2w_f32(rows, value("GE", 0.0));
                for (const std::array<double, 4>& ply : plies) {
                    op2w_i32(rows, static_cast<std::int64_t>(ply[0]));
                    op2w_f32(rows, ply[1]);
                    op2w_f32(rows, ply[2]);
                    op2w_i32(rows, static_cast<std::int64_t>(ply[3]));
                }
            }
            records.push_back(op2w_key(2706, 27, 287) + rows);
        }
        if (!mProps[kOp2wPbeam].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPbeam])
                rows += op2w_pbeam_row(pid, mArrays);
            records.push_back(op2w_key(5402, 54, 262) + rows);
        }
        if (!mProps[kOp2wPbush].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPbush])
                rows += op2w_pbush_row(pid, mBushes);
            records.push_back(op2w_key(1402, 14, 37) + rows);
        }
        if (!mProps[kOp2wPbarl].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPbarl])
                rows += op2w_pbarl_row(pid, mPbarl.at(pid));
            records.push_back(op2w_key(9102, 91, 52) + rows);
        }
        if (!mProps[kOp2wPbeaml].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPbeaml])
                rows += op2w_pbeaml_row(pid, mPbeaml.at(pid));
            records.push_back(op2w_key(9202, 92, 53) + rows);
        }
        if (!mProps[kOp2wPcompg].empty()) {
            std::string rows;
            for (const std::int64_t pid : mProps[kOp2wPcompg])
                rows += op2w_pcompg_row(pid, mPcompg.at(pid));
            records.push_back(op2w_key(15006, 150, 604) + rows);
        }
        records.push_back(op2w_key(65535, 65535, 65535));
        return op2w_table("EPT", records);
    }

    /// The card of each material id the field data gives (`nastran:MAT1_<id>:*`, `MAT2` or
    /// `MAT8`); an id given as two cards is an error. An id with none is a MAT1.
    std::map<std::int64_t, std::string> MaterialCards() const {
        std::map<std::int64_t, std::string> cards;
        for (const auto& [key, members] : mValues) {
            if (key.first != "MAT1" && key.first != "MAT2" && key.first != "MAT4" &&
                key.first != "MAT5" && key.first != "MAT8" && key.first != "MAT9")
                continue;
            const auto [it, added] = cards.emplace(key.second, key.first);
            if (!added)
                op2w_fail("material id " + std::to_string(key.second) + " is given as both " +
                          it->second + " and " + key.first + " field data");
        }
        return cards;
    }

    /// A MAT1 per material the properties name (unit modulus unless the field data gives
    /// its values), and per `nastran:MAT1_<id>:*` given; a material given as `nastran:MAT2_<id>:*`
    /// or `nastran:MAT8_<id>:*` is written as that card instead.
    std::string Mpt() const {
        const std::map<std::int64_t, std::string> given_cards = MaterialCards();
        std::set<std::int64_t> ids;
        for (const auto& [pid, plies] : mPlies)
            if (mProps[kOp2wPcomp].count(pid) != 0)
                for (const std::array<double, 4>& ply : plies)
                    ids.insert(static_cast<std::int64_t>(ply[0]));
        for (const int prop :
             {kOp2wPshell, kOp2wPsolid, kOp2wProd, kOp2wPtube, kOp2wPbar, kOp2wPshear}) {
            for (const std::int64_t pid : mProps[prop]) {
                std::int64_t mid = 1;
                const auto given = mValues.find({kOp2wPropNames[prop], pid});
                if (given != mValues.end()) {
                    const auto it = given->second.find(prop == kOp2wPshell ? "MID1" : "MID");
                    if (it != given->second.end())
                        mid = static_cast<std::int64_t>(it->second);
                }
                ids.insert(mid);
            }
        }
        for (const std::int64_t pid : mProps[kOp2wPbeam])
            ids.insert(static_cast<std::int64_t>(op2w_array_scalar(mArrays, pid, "MID", 1)));
        for (const std::int64_t pid : mProps[kOp2wPbarl])
            ids.insert(static_cast<std::int64_t>(op2w_array_scalar(mPbarl, pid, "MID", 1)));
        for (const std::int64_t pid : mProps[kOp2wPbeaml])
            ids.insert(static_cast<std::int64_t>(op2w_array_scalar(mPbeaml, pid, "MID", 1)));
        for (const std::int64_t pid : mProps[kOp2wPcompg]) {
            const std::vector<double>& plies = mPcompg.at(pid).at("PLY");
            for (std::size_t j = 1; j < plies.size(); j += 5)
                ids.insert(static_cast<std::int64_t>(plies[j]));
        }
        for (const auto& [mid, card] : given_cards)
            ids.insert(mid);
        if (ids.empty())
            ids.insert(1);
        static const struct {
            const char* mCard;
            std::int32_t mK1, mK2, mK3;
            std::size_t mSize;
        } kMaterials[] = {{"MAT1", 103, 1, 77, 12},    {"MAT2", 203, 2, 78, 17},
                          {"MAT4", 2103, 21, 234, 11}, {"MAT5", 2203, 22, 235, 10},
                          {"MAT8", 2503, 25, 288, 19}, {"MAT9", 2603, 26, 300, 35}};
        std::vector<std::string> records;
        for (const auto& k : kMaterials) {
            std::string rows;
            for (const std::int64_t mid : ids) {
                const auto at = given_cards.find(mid);
                if ((at == given_cards.end() ? "MAT1" : at->second) != k.mCard)
                    continue;
                Op2wRowBuf w(k.mSize);
                w.I(0, mid);
                const auto given = mValues.find({k.mCard, mid});
                for (const Op2wValueMember& m : op2w_value_card(k.mCard).mMembers) {
                    double v = m.mDefault;
                    if (given != mValues.end()) {
                        const auto it = given->second.find(m.mName);
                        if (it != given->second.end())
                            v = it->second;
                    }
                    if (m.mInt)
                        w.I(m.mWord, static_cast<std::int64_t>(v));
                    else
                        w.F(m.mWord, v);
                }
                w.Append(rows);
            }
            if (!rows.empty())
                records.push_back(op2w_key(k.mK1, k.mK2, k.mK3) + rows);
        }
        records.push_back(op2w_key(65535, 65535, 65535));
        return op2w_table("MPT", records);
    }

private:
    Op2wSignature mSignature;
    std::size_t mNpoints = 0;
    std::vector<double> mPoints;
    std::vector<std::int64_t> mCp;       // each point's coordinate system (0: basic)
    std::vector<std::int64_t> mCd;       // the system each point's results are in (0: basic)
    std::map<std::int64_t, std::size_t> mCordAt;  // a system's id -> its place in mCords
    std::vector<double> mLocal;          // each point's coordinates in that system
    std::vector<Op2wCord> mCords;        // the systems the CPs name, by id
    std::vector<std::int64_t> mNodeIds;  // the GRID id of each point
    std::vector<Op2wBlock> mBlocks;
    std::vector<std::vector<Op2wRow>> mRows;
    std::set<std::int64_t> mProps[17];
    std::set<std::string> mSkipped;
    std::set<std::string> mUnwritten;  // element cards the file does not hold yet
    Op2wValues mValues;
    Op2wPlies mPlies;
    Op2wArrays mArrays;
    Op2wArrays mBushes;
    Op2wArrays mPbarl;  // the line-property field data, by id: texts as their characters
    Op2wArrays mPbeaml;
    Op2wArrays mPcompg;
};

// -- steps -----------------------------------------------------------------------------

std::int64_t op2w_scalar_int(const Mesh& rMesh, const std::string& rName, std::int64_t Default) {
    if (!rMesh.HasFieldData(rName))
        return Default;
    const NDArray& a = rMesh.FieldData(rName);
    if (a.Size() == 0)
        return Default;
    return detail::Int64View(a).Data()[0];
}

double op2w_time_of(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kSequenceTimeKey))
        return 0.0;
    const NDArray& a = rMesh.FieldData(kSequenceTimeKey);
    return a.Size() == 0 ? 0.0 : detail::DoubleView(a).Data()[0];
}

/// The step words every result header carries.
struct Op2wHeader {
    Op2wHeader(const Mesh& rMesh, double Time) {
        mTime = static_cast<double>(op2w_float(Time));
        mSubcase = op2w_scalar_int(rMesh, "nastran:subcase", 1);
        mAnalysis = op2w_scalar_int(rMesh, "nastran:analysis", mTime == 0.0 ? 1 : 6);
        mMode = op2w_scalar_int(rMesh, "nastran:mode", 1);
        if (mSubcase <= 0 || mAnalysis <= 0 || mAnalysis > 99)
            op2w_fail("'nastran:subcase' and 'nastran:analysis' must be positive codes");
        mW5 = 1;
        mW6 = 0.0;
        if (mAnalysis == 5 || mAnalysis == 6 || mAnalysis == 10 || mAnalysis == 12) {
            mW5 = std::bit_cast<std::int32_t>(op2w_float(mTime));
        } else if (mAnalysis == 2 || mAnalysis == 8 || mAnalysis == 9) {
            mW5 = mMode;
            mW6 = mTime;
        }
    }

    std::tuple<std::int64_t, std::int64_t, std::int64_t> Key() const {
        return {mSubcase, mAnalysis, mW5};
    }

    std::string Words(std::int32_t TCode, std::int32_t EType, std::size_t NumWide,
                      std::int32_t SCode = 0, std::int32_t Thermal = 0,
                      std::int32_t FormatCode = 1) const {
        Op2wRowBuf w(146);
        for (std::size_t i = 50; i < 146; ++i)
            w.mWords[i] = 0x20202020u;  // the title, subtitle and label: blanks
        w.I(0, mAnalysis * 10 + kOp2wDevice);
        w.I(1, TCode);
        w.I(2, EType);
        w.I(3, mSubcase);
        w.mWords[4] = static_cast<std::uint32_t>(mW5);
        w.F(5, mW6);
        w.I(8, FormatCode);
        w.I(9, static_cast<std::int64_t>(NumWide));
        w.I(10, SCode);
        w.I(22, Thermal);
        std::string out;
        w.Append(out);
        return out;
    }

    double mTime = 0.0;
    std::int64_t mSubcase = 1, mAnalysis = 1, mMode = 1, mW5 = 1;
    double mW6 = 0.0;
};

/// One step's records per result table, and the (table, table code, element
/// type) of every block.
struct Op2wStep {
    std::map<std::string, std::vector<std::string>> mTables;
    std::vector<std::tuple<std::string, std::int32_t, std::int32_t>> mLayout;
    std::tuple<std::int64_t, std::int64_t, std::int64_t> mKey;
};

struct Op2wNodal {
    const char* mName;
    const char* mTable;
    std::int32_t mTCode;
    bool mThermal;
};
constexpr Op2wNodal kOp2wNodal[] = {
    {"DISPLACEMENT", "OUGV1", 1, false}, {"EIGENVECTOR", "OUGV1", 7, false},
    {"VELOCITY", "OUGV1", 10, false},    {"ACCELERATION", "OUGV1", 11, false},
    {"TEMPERATURE", "OUGV1", 1, true},   {"SPC_FORCE", "OQG1", 3, false},
    {"MPC_FORCE", "OQMG1", 3, false},    {"APPLIED_LOAD", "OPG1", 2, false},
};

/// A cell-data array of the cells of a card, in the card's row order: one value per
/// row (the first column of a wider array), or @p Width per row for a per-location
/// array. Empty when no block has it.
std::vector<double> op2w_gather(const Mesh& rMesh, const Op2wModel& rModel,
                                const std::string& rName, const std::vector<Op2wRow>& rRows,
                                std::size_t Width) {
    if (!rMesh.HasCellData(rName))
        return {};
    const std::size_t stride = Width == 0 ? 1 : Width;
    std::vector<double> out(rRows.size() * stride, std::numeric_limits<double>::quiet_NaN());
    std::vector<char> used(rModel.Blocks().size(), 0);
    for (const Op2wRow& row : rRows)
        used[row.mBlock] = 1;
    std::vector<std::unique_ptr<detail::DoubleView>> views(rModel.Blocks().size());
    std::vector<std::size_t> columns(rModel.Blocks().size(), 0);
    for (std::size_t bi = 0; bi < rModel.Blocks().size(); ++bi) {
        const Op2wBlock& blk = rModel.Blocks()[bi];
        if (blk.mCount == 0 || !used[bi])
            continue;
        const NDArray& a = rMesh.CellData(rName, blk.mIndex);
        const auto& shape = a.Shape();
        if (Width == 0) {
            if (a.Size() < blk.mCount || a.Size() % blk.mCount != 0)
                op2w_fail("'" + rName + "' must have one value per cell");
            columns[bi] = a.Size() / blk.mCount;
        } else {
            if (shape.size() != 2 || shape[0] != blk.mCount)
                continue;
            columns[bi] = shape[1];
        }
        views[bi] = std::make_unique<detail::DoubleView>(a);
    }
    for (std::size_t r = 0; r < rRows.size(); ++r) {
        const std::size_t bi = rRows[r].mBlock;
        if (!views[bi])
            continue;
        const double* row = views[bi]->Data() + rRows[r].mRow * columns[bi];
        for (std::size_t k = 0; k < (Width == 0 ? 1 : std::min(Width, columns[bi])); ++k)
            out[r * stride + k] = row[k];
    }
    return out;
}

// -- element result tables, driven by the reader's layouts --------------------------
//
// The reader decodes an element table through `op2_element_layout` (which words hold
// which member, per element type, word count and complex/random flag); the writer uses
// the same table the other way round, so what one reads the other writes. An element
// type is chosen from the arrays the mesh has: for a card, the first of its alternatives
// (element type, whether its location arrays are required) the mesh has data for.

constexpr std::uint32_t kOp2wCen = 0x2F4E4543u;   // 'CEN/'
constexpr std::uint32_t kOp2wGrid = 0x44495247u;  // 'GRID'

struct Op2wResultAlt {
    std::int32_t mEtype;
    bool mStrict;  // needs its location arrays (`@corner`, `@station`, `@ply`)
};

struct Op2wResultCard {
    const char* mCard;
    std::vector<Op2wResultAlt> mAlts;
};

/// The cards that have element results, in file order, with the element types of their
/// tables.
const std::vector<Op2wResultCard>& op2w_result_cards() {
    static const std::vector<Op2wResultCard> cards = {
        {"CROD", {{1, false}}},
        {"CTUBE", {{3, false}}},
        {"CONROD", {{10, false}}},
        {"CQUAD4", {{95, true}, {144, true}, {33, false}}},
        {"CTETRA", {{39, false}}},
        {"CHEXA", {{67, false}}},
        {"CPENTA", {{68, false}}},
        {"CTRIA3", {{97, true}, {74, false}}},
        {"CPYRAM", {{255, false}}},
        {"CQUAD8", {{96, true}, {64, false}}},
        {"CTRIA6", {{98, true}, {75, false}}},
        {"CBAR", {{100, true}, {34, false}}},
        {"CBEAM", {{2, false}}},
        {"CBUSH", {{102, false}}},
        {"CSHEAR", {{4, false}}},
        {"CELAS1", {{11, false}}},
        {"CELAS2", {{12, false}}},
        {"CDAMP1", {{20, false}}},
        {"CDAMP2", {{21, false}}},
        {"CVISC", {{24, false}}},
    };
    return cards;
}

struct Op2wResultGroup {
    const char* mName;  // the cell-data prefix
    Op2Family mFamily;
    std::int32_t mBase;  // the table code
    const char* mTable;
    const char* mRandomPrefix;  // a random table's name starts with it
    std::int32_t mThermal;
};
constexpr Op2wResultGroup kOp2wResultGroups[] = {
    {"ELEMENT_FORCE", Op2Family::Force, 4, "OEF1X", "OEF", 0},
    {"STRESS", Op2Family::Stress, 5, "OES1X1", "OES", 0},
    {"STRAIN", Op2Family::Stress, 5, "OSTR1X", "OSTR", 0},
    {"HEAT_FLUX", Op2Family::Flux, 4, "HOEF1", "", 1},
    {"ENERGY", Op2Family::Energy, 18, "ONRGY1", "", 0},
    {"KINETIC_ENERGY", Op2Family::Energy, 36, "OEKE1", "", 0},
};

// The random tables' suffixes and the hundreds of their table codes, in file order.
constexpr std::pair<const char*, std::int32_t> kOp2wRandomKinds[] = {
    {"_CRM", 5}, {"_PSD", 6}, {"_ATO", 7}, {"_RMS", 8}, {"_NO", 9}};

/// The result tables in file order: the model's, then the random ones.
const std::vector<std::string>& op2w_table_order() {
    static const std::vector<std::string> tables = [] {
        std::vector<std::string> out = {"OUGV1",  "OQG1",   "OQMG1",  "OPG1",  "OEF1X", "HOEF1",
                                        "OES1X1", "OSTR1X", "ONRGY1", "OEKE1", "OGPFB1"};
        for (const char* prefix : {"OES", "OSTR", "OEF"})
            for (const auto& kind : kOp2wRandomKinds)
                out.push_back(std::string(prefix) + (kind.first + 1) + "1");
        return out;
    }();
    return tables;
}

struct Op2wCandidate {
    Op2ElementLayout mLayout;
    std::size_t mNumWide;
    std::int32_t mSCode;
};

/// The layouts a (family, element type, complex, random) can take, scanned once from
/// the reader's table: the stress code with von Mises first, so that it wins a tie.
const std::vector<Op2wCandidate>& op2w_candidates(Op2Family Family, std::int32_t Etype,
                                                  bool Complex, bool Random) {
    static std::mutex lock;
    static std::map<std::tuple<int, std::int32_t, bool, bool>, std::vector<Op2wCandidate>> cache;
    const std::lock_guard<std::mutex> guard(lock);
    const auto key = std::make_tuple(static_cast<int>(Family), Etype, Complex, Random);
    auto it = cache.find(key);
    if (it != cache.end())
        return it->second;
    std::vector<Op2wCandidate> found;
    for (const std::int32_t s_code : {1, 0})
        for (std::size_t wide = 2; wide < 400; ++wide) {
            auto layout = op2_element_layout(Family, Etype, static_cast<std::int64_t>(wide), s_code,
                                             Complex, Random);
            if (layout)
                found.push_back({std::move(*layout), wide, s_code});
        }
    return cache.emplace(key, std::move(found)).first->second;
}

/// The suffix of a layout's location arrays, or nullptr for plain rows.
const char* op2w_location(const Op2ElementLayout& rLayout) {
    if (rLayout.mKind == Op2Values::Ply)
        return "@ply";
    if (rLayout.mKind == Op2Values::Station || rLayout.mStations)
        return "@station";
    if (rLayout.mKind == Op2Values::Blocks)
        return "@corner";
    return nullptr;
}

/// The array names of a layout member and the word of each: itself (with a random
/// table's suffix) or its real and imaginary parts.
std::vector<std::pair<std::string, std::size_t>> op2w_names_of(const Op2Member& rMember,
                                                               const std::string& rSuffix) {
    if (rMember.mImag == kOp2Npos)
        return {{rMember.mName + rSuffix, rMember.mWord}};
    return {{rMember.mName + "_real", rMember.mWord}, {rMember.mName + "_imag", rMember.mImag}};
}

/// A per-location array, `mWidth` columns per row; zero columns: the mesh has none.
struct Op2wWide {
    std::vector<double> mData;
    std::size_t mWidth = 0;
    double At(std::size_t Row, std::size_t Col) const { return mData[Row * mWidth + Col]; }
};

Op2wWide op2w_wide_gather(const Mesh& rMesh, const Op2wModel& rModel, const std::string& rName,
                          const std::vector<Op2wRow>& rRows) {
    Op2wWide out;
    if (!rMesh.HasCellData(rName))
        return out;
    std::vector<char> used(rModel.Blocks().size(), 0);
    for (const Op2wRow& row : rRows)
        used[row.mBlock] = 1;
    for (std::size_t bi = 0; bi < rModel.Blocks().size(); ++bi) {
        if (!used[bi])
            continue;
        const NDArray& a = rMesh.CellData(rName, rModel.Blocks()[bi].mIndex);
        const auto& shape = a.Shape();
        const std::size_t width = shape.size() == 2 ? shape[1] : (a.Size() ? 1 : 0);
        out.mWidth = std::max(out.mWidth, width);
    }
    if (out.mWidth == 0)
        return out;
    out.mData = op2w_gather(rMesh, rModel, rName, rRows, out.mWidth);
    return out;
}

double op2w_zero_nan(double V) {
    return std::isnan(V) ? 0.0 : V;
}

/// One table block of a card in one element type; false when the mesh holds none of
/// its arrays.
bool op2w_card_table(const Mesh& rMesh, const Op2wModel& rModel, const Op2wHeader& rHeader,
                     Op2wStep& rStep, std::set<std::string>& rTaken, const Op2wResultGroup& rGroup,
                     const char* pCard, std::int32_t Etype, bool Strict, const std::string& rSuffix,
                     bool Complex, bool Random) {
    const std::size_t card = op2w_card_index(pCard);
    const std::vector<Op2wRow>& rows = rModel.Rows(card);
    if (rows.empty())
        return false;
    const bool stress_like = rGroup.mFamily == Op2Family::Stress;
    const std::vector<Op2wCandidate>& candidates =
        op2w_candidates(rGroup.mFamily, Etype, Complex, Random && stress_like);
    if (candidates.empty())
        return false;
    const std::string prefix = std::string(rGroup.mName) + ":";
    std::map<std::string, std::vector<double>> plain_cache;
    std::map<std::string, Op2wWide> wide_cache;
    const auto plain = [&](const std::string& rName) -> const std::vector<double>& {
        auto it = plain_cache.find(rName);
        if (it == plain_cache.end())
            it = plain_cache.emplace(rName, op2w_gather(rMesh, rModel, prefix + rName, rows, 0))
                     .first;
        return it->second;
    };
    const auto wide = [&](const std::string& rName) -> const Op2wWide& {
        auto it = wide_cache.find(rName);
        if (it == wide_cache.end())
            it = wide_cache.emplace(rName, op2w_wide_gather(rMesh, rModel, prefix + rName, rows))
                     .first;
        return it->second;
    };

    const Op2wCandidate* best = nullptr;
    int best_score = 0;
    for (const Op2wCandidate& cand : candidates) {
        const Op2ElementLayout& lay = cand.mLayout;
        const char* loc = op2w_location(lay);
        const bool is_plain = !loc || (lay.mKind == Op2Values::Blocks && !lay.mStations);
        int score = 0, own = 0;
        bool have_loc = false;
        for (const Op2Member& member : lay.mMembers)
            for (const auto& [name, word] : op2w_names_of(member, rSuffix)) {
                // a complex or random table needs an array of its own: the real-only
                // members of a complex layout (the fibre distances) are shared
                const bool is_own =
                    !(Complex || Random) || name != member.mName || !rSuffix.empty();
                if (is_plain && !plain(name).empty()) {
                    ++score;
                    own += is_own;
                }
                if (loc && wide(name + loc).mWidth != 0) {
                    ++score;
                    own += is_own;
                    have_loc = true;
                }
            }
        if ((Strict && !have_loc) || !score || !own)
            continue;
        if (!best || score > best_score) {
            best = &cand;
            best_score = score;
        }
    }
    if (!best)
        return false;
    const Op2ElementLayout& lay = best->mLayout;
    const std::size_t num_wide = best->mNumWide;
    const char* loc = op2w_location(lay);
    const bool by_entry = lay.mKind == Op2Values::Station || lay.mKind == Op2Values::Ply;
    std::vector<std::pair<std::string, std::size_t>> members;
    for (const Op2Member& member : lay.mMembers)
        for (auto& pair : op2w_names_of(member, rSuffix))
            members.push_back(std::move(pair));
    std::map<std::string, const std::vector<double>*> plain_values;
    std::map<std::string, const Op2wWide*> loc_values;
    for (const auto& [name, word] : members) {
        const std::vector<double>& a = plain(name);
        if (!a.empty() &&
            (lay.mKind == Op2Values::Row || (lay.mKind == Op2Values::Blocks && !lay.mStations))) {
            plain_values[name] = &a;
            rTaken.insert(prefix + name);
        }
        if (loc) {
            const Op2wWide& b = wide(name + loc);
            if (b.mWidth != 0) {
                loc_values[name] = &b;
                rTaken.insert(prefix + name + loc);
            }
        }
    }

    const std::size_t count = rows.size();
    std::string data;
    double energy_total = 0.0;  // an energy table's header holds the sum of its energies
    if (by_entry) {
        // one entry per station or ply with any value, in element order
        std::size_t most = 0;
        for (const auto& [name, b] : loc_values)
            most = std::max(most, b->mWidth);
        std::size_t kept = 0;
        for (std::size_t r = 0; r < count; ++r)
            for (std::size_t c = 0; c < most; ++c) {
                bool any = false;
                for (const auto& [name, b] : loc_values)
                    any = any || (c < b->mWidth && !std::isnan(b->At(r, c)));
                if (!any)
                    continue;
                ++kept;
                Op2wRowBuf w(num_wide);
                w.I(0, rows[r].mEid * 10 + kOp2wDevice);
                if (lay.mKind == Op2Values::Ply)
                    w.I(1, static_cast<std::int64_t>(c) + 1);
                for (const auto& [name, word] : members) {
                    const auto it = loc_values.find(name);
                    if (it == loc_values.end() || c >= it->second->mWidth)
                        continue;
                    w.F(word, op2w_zero_nan(it->second->At(r, c)));
                }
                w.Append(data);
            }
        if (kept == 0)
            return false;
    } else {
        std::vector<std::size_t> keep;
        for (std::size_t r = 0; r < count; ++r) {
            bool shown = false;
            for (const auto& [name, a] : plain_values)
                shown = shown || !std::isnan((*a)[r]);
            for (const auto& [name, b] : loc_values)
                for (std::size_t c = 0; c < b->mWidth; ++c)
                    shown = shown || !std::isnan(b->At(r, c));
            if (shown)
                keep.push_back(r);
        }
        if (keep.empty())
            return false;
        const std::int64_t solid_nodes = op2_solid_nodes(Etype);
        const std::int64_t plate_nodes = op2_corner_nodes(Etype);
        const std::size_t ncorner =
            lay.mStations ? 0
                          : static_cast<std::size_t>((solid_nodes ? solid_nodes : plate_nodes) - 1);
        for (const std::size_t r : keep) {
            Op2wRowBuf w(num_wide);
            w.I(0, rows[r].mEid * 10 + kOp2wDevice);
            if (lay.mKind == Op2Values::Row) {
                for (const auto& [name, word] : members) {
                    const auto it = plain_values.find(name);
                    if (it != plain_values.end())
                        w.F(word, op2w_zero_nan((*it->second)[r]));
                }
                if (rGroup.mFamily == Op2Family::Energy) {
                    const auto it = plain_values.find(members[0].first);
                    if (it != plain_values.end())
                        energy_total += op2w_zero_nan((*it->second)[r]);
                }
                w.Append(data);
                continue;
            }
            const Op2wBlock& blk = rModel.Blocks()[rows[r].mBlock];
            const std::int64_t* conn = blk.mConn.data() + rows[r].mRow * blk.mNodes;
            const auto node_id = [&](std::size_t K) {
                return rModel.NodeIds()[static_cast<std::size_t>(conn[K])];
            };
            if (plate_nodes)
                w.mWords[1] = kOp2wCen;
            if (solid_nodes) {
                w.mWords[2] = kOp2wGrid;
                w.I(3, static_cast<std::int64_t>(ncorner));
            }
            for (std::size_t k = 0; k < lay.mBlocks; ++k) {
                const std::size_t base = lay.mFirst + k * lay.mBlock;
                if (lay.mStations) {
                    // the ends carry their GRIDs, the stations between them none
                    if (k == 0)
                        w.I(base, node_id(0));
                    else if (k == lay.mBlocks - 1)
                        w.I(base, node_id(1));
                } else if (k == 0) {  // the centre
                    w.mWords[base] = solid_nodes ? kOp2wCen : 0u;
                } else {
                    w.I(base, node_id(k - 1));
                }
                for (const auto& [name, word] : members) {
                    if (lay.mStations) {
                        const auto it = loc_values.find(name);
                        if (it != loc_values.end() && k < it->second->mWidth)
                            w.F(base + word, op2w_zero_nan(it->second->At(r, k)));
                        continue;
                    }
                    double v = 0.0;
                    const auto pit = plain_values.find(name);
                    if (pit != plain_values.end())
                        v = op2w_zero_nan((*pit->second)[r]);
                    if (k > 0) {
                        const auto lit = loc_values.find(name);
                        if (lit != loc_values.end() && k - 1 < lit->second->mWidth) {
                            const double c = lit->second->At(r, k - 1);
                            if (!std::isnan(c))
                                v = c;
                        }
                    }
                    w.F(base + word, v);
                }
            }
            w.Append(data);
        }
    }
    std::int32_t tcode = rGroup.mBase;
    std::string name = rGroup.mTable;
    if (Random) {
        for (const auto& kind : kOp2wRandomKinds)
            if (rSuffix == kind.first) {
                tcode += 4000 + 100 * kind.second;
                name = std::string(rGroup.mRandomPrefix) + (kind.first + 1) + "1";
            }
    } else if (Complex) {
        tcode += 1000;
    }
    const std::int32_t s_out =
        stress_like ? ((std::string_view(rGroup.mName) == "STRAIN" ? 10 : 0) | (best->mSCode & 1))
                    : 0;
    rStep.mLayout.emplace_back(name, tcode, Etype);
    auto& records = rStep.mTables[name];
    const std::int32_t header_type = rGroup.mFamily == Op2Family::Energy
                                         ? std::bit_cast<std::int32_t>(op2w_float(energy_total))
                                         : Etype;
    records.push_back(
        rHeader.Words(tcode, header_type, num_wide, s_out, rGroup.mThermal, Complex ? 2 : 1));
    records.push_back(std::move(data));
    return true;
}

/// Every element result table of the step: the forces, stresses, strains and heat
/// fluxes of the cards that have results, in the reader's layouts.
void op2w_result_tables(const Mesh& rMesh, const Op2wModel& rModel, const Op2wHeader& rHeader,
                        Op2wStep& rStep, std::set<std::string>& rTaken) {
    for (const Op2wResultGroup& group : kOp2wResultGroups) {
        struct Variant {
            std::string mSuffix;
            bool mComplex, mRandom;
        };
        std::vector<Variant> variants = {{"", false, false}, {"", true, false}};
        if (group.mFamily == Op2Family::Stress || group.mFamily == Op2Family::Force)
            for (const auto& kind : kOp2wRandomKinds)
                variants.push_back({kind.first, false, true});
        for (const Variant& variant : variants)
            for (const Op2wResultCard& card : op2w_result_cards()) {
                // an energy table holds any element
                static const std::vector<Op2wResultAlt> any = {{0, false}};
                for (const Op2wResultAlt& alt :
                     group.mFamily == Op2Family::Energy ? any : card.mAlts)
                    if (op2w_card_table(rMesh, rModel, rHeader, rStep, rTaken, group, card.mCard,
                                        alt.mEtype, alt.mStrict, variant.mSuffix, variant.mComplex,
                                        variant.mRandom))
                        break;
            }
    }
}

constexpr const char* kOp2wGpfMembers[6] = {"F1", "F2", "F3", "M1", "M2", "M3"};

/// Word @p K (0 or 1) of an 8-character text, blank padded.
std::uint32_t op2w_text_word(const std::string& rText, std::size_t K) {
    std::string t = rText.substr(0, 8);
    t.resize(8, ' ');
    std::uint32_t w = 0;
    std::memcpy(&w, t.data() + 4 * K, 4);
    return w;
}

/// Grid point forces, one `OGPFB1` record: a row per node of each element with
/// `GRID_FORCE:F1`...`M3` cell data (`(cells, nodes)` per block, a column per node of
/// the connectivity), then a row per point of each `GRID_FORCE:<label>:<M>` point data
/// (label `TOTALS`, `APP-LOAD`, `F-OF-SPC`...) with element 0.
void op2w_grid_force_tables(const Mesh& rMesh, const Op2wModel& rModel, const Op2wHeader& rHeader,
                            Op2wStep& rStep, std::set<std::string>& rTaken) {
    std::string data;
    std::size_t kept = 0;
    // The forces and moments (basic) in the output system of the row's point.
    const auto output = [&](std::size_t Point, double* pV) {
        if (rModel.Cd()[Point] == 0)
            return;
        const Op2wCord& system = rModel.OutputSystem(Point);
        const double* at = rModel.BasicPoints().data() + 3 * Point;
        double out[3];
        op2w_vector_to_local(system, at, pV, out);
        std::copy(out, out + 3, pV);
        op2w_vector_to_local(system, at, pV + 3, out);
        std::copy(out, out + 3, pV + 3);
    };
    for (std::size_t card = 0; card < kOp2wNumCards; ++card) {
        const std::vector<Op2wRow>& rows = rModel.Rows(card);
        if (rows.empty())
            continue;
        std::array<Op2wWide, 6> arrays;
        std::size_t width = 0;
        for (std::size_t k = 0; k < 6; ++k) {
            arrays[k] = op2w_wide_gather(rMesh, rModel,
                                         std::string("GRID_FORCE:") + kOp2wGpfMembers[k], rows);
            if (arrays[k].mWidth != 0) {
                rTaken.insert(std::string("GRID_FORCE:") + kOp2wGpfMembers[k]);
                width = std::max(width, arrays[k].mWidth);
            }
        }
        if (width == 0)
            continue;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const Op2wBlock& blk = rModel.Blocks()[rows[r].mBlock];
            const std::int64_t* conn = blk.mConn.data() + rows[r].mRow * blk.mNodes;
            for (std::size_t c = 0; c < std::min(width, blk.mNodes); ++c) {
                bool shown = false;
                for (const Op2wWide& a : arrays)
                    shown = shown || (c < a.mWidth && !std::isnan(a.At(r, c)));
                if (!shown)
                    continue;
                ++kept;
                Op2wRowBuf w(10);
                const auto point = static_cast<std::size_t>(conn[c]);
                w.I(0, rModel.NodeIds()[point] * 10 + kOp2wDevice);
                w.I(1, rows[r].mEid);
                w.mWords[2] = op2w_text_word(kOp2wCards[card].mName, 0);
                w.mWords[3] = op2w_text_word(kOp2wCards[card].mName, 1);
                double v[6];
                for (std::size_t k = 0; k < 6; ++k)
                    v[k] = c < arrays[k].mWidth ? op2w_zero_nan(arrays[k].At(r, c)) : 0.0;
                output(point, v);
                for (std::size_t k = 0; k < 6; ++k)
                    w.F(4 + k, v[k]);
                w.Append(data);
            }
        }
    }
    std::set<std::string> labels;
    for (const std::string& n : rMesh.PointDataNames()) {
        if (n.rfind("GRID_FORCE:", 0) != 0)
            continue;
        const std::size_t last = n.rfind(':');
        const std::size_t first = n.find(':');
        if (last == first || n.find(':', first + 1) != last)
            continue;
        const std::string member = n.substr(last + 1);
        if (std::find(std::begin(kOp2wGpfMembers), std::end(kOp2wGpfMembers), member) !=
            std::end(kOp2wGpfMembers))
            labels.insert(n.substr(first + 1, last - first - 1));
    }
    const std::size_t npts = rModel.NumPoints();
    std::vector<std::size_t> by_id(npts);  // the rows are by GRID id
    std::iota(by_id.begin(), by_id.end(), std::size_t{0});
    std::stable_sort(by_id.begin(), by_id.end(), [&rModel](std::size_t A, std::size_t B) {
        return rModel.NodeIds()[A] < rModel.NodeIds()[B];
    });
    for (const std::string& label : labels) {
        std::array<std::optional<detail::DoubleView>, 6> columns;
        for (std::size_t k = 0; k < 6; ++k) {
            const std::string name = "GRID_FORCE:" + label + ":" + kOp2wGpfMembers[k];
            if (!rMesh.HasPointData(name))
                continue;
            rTaken.insert(name);
            const NDArray& a = rMesh.PointData(name);
            if (a.Size() != npts)
                op2w_fail("'" + name + "' must have one value per point");
            columns[k].emplace(a);
        }
        const std::string text = label == "TOTALS" ? "*TOTALS*" : label;
        for (const std::size_t i : by_id) {
            bool shown = false;
            for (const auto& c : columns)
                shown = shown || (c && !std::isnan(c->Data()[i]));
            if (!shown)
                continue;
            ++kept;
            Op2wRowBuf w(10);
            w.I(0, rModel.NodeIds()[i] * 10 + kOp2wDevice);
            w.mWords[2] = op2w_text_word(text, 0);
            w.mWords[3] = op2w_text_word(text, 1);
            double v[6];
            for (std::size_t k = 0; k < 6; ++k)
                v[k] = columns[k] ? op2w_zero_nan(columns[k]->Data()[i]) : 0.0;
            output(i, v);
            for (std::size_t k = 0; k < 6; ++k)
                w.F(4 + k, v[k]);
            w.Append(data);
        }
    }
    if (kept == 0)
        return;
    rStep.mLayout.emplace_back("OGPFB1", 19, 0);
    auto& records = rStep.mTables["OGPFB1"];
    records.push_back(rHeader.Words(19, 0, 10));
    records.push_back(std::move(data));
}

Op2wStep op2w_step_tables(const Mesh& rMesh, const Op2wModel& rModel, double Time,
                          std::set<std::string>& rTaken) {
    const Op2wHeader header(rMesh, Time);
    Op2wStep step;
    step.mKey = header.Key();
    const std::size_t npts = rModel.NumPoints();
    for (const Op2wNodal& nodal : kOp2wNodal) {
        if (!rMesh.HasPointData(nodal.mName))
            continue;
        rTaken.insert(nodal.mName);
        const NDArray& base_array = rMesh.PointData(nodal.mName);
        const std::size_t ncomp = nodal.mThermal ? 1 : 3;
        if (base_array.Size() != npts * ncomp)
            op2w_fail(std::string("'") + nodal.mName + "' must have " + std::to_string(ncomp) +
                      " component(s) per point");
        const detail::DoubleView base_view(base_array);
        std::vector<double> base(base_view.Data(), base_view.Data() + npts * ncomp);
        if (!nodal.mThermal)
            rModel.ToOutput(base);  // in each point's output system
        std::vector<double> rot;
        bool has_rot = false;
        if (!nodal.mThermal) {
            const std::string rot_name = std::string(nodal.mName) + "_ROT";
            if (rMesh.HasPointData(rot_name)) {
                rTaken.insert(rot_name);
                const NDArray& rot_array = rMesh.PointData(rot_name);
                if (rot_array.Size() != npts * 3)
                    op2w_fail("'" + rot_name + "' must have 3 components per point");
                const detail::DoubleView rot_view(rot_array);
                rot.assign(rot_view.Data(), rot_view.Data() + npts * 3);
                rModel.ToOutput(rot);
                has_rot = true;
            }
        }
        std::string data;
        std::size_t kept = 0;
        std::vector<std::size_t> by_id(npts);  // the rows are by GRID id
        std::iota(by_id.begin(), by_id.end(), std::size_t{0});
        std::stable_sort(by_id.begin(), by_id.end(), [&rModel](std::size_t A, std::size_t B) {
            return rModel.NodeIds()[A] < rModel.NodeIds()[B];
        });
        for (const std::size_t i : by_id) {
            bool all_nan = true;
            for (std::size_t k = 0; k < ncomp; ++k)
                all_nan = all_nan && std::isnan(base[i * ncomp + k]);
            if (all_nan)
                continue;
            ++kept;
            Op2wRowBuf w(8);
            w.I(0, rModel.NodeIds()[i] * 10 + kOp2wDevice);
            w.I(1, 1);
            for (std::size_t k = 0; k < ncomp; ++k)
                w.F(2 + k, base[i * ncomp + k]);
            if (has_rot)
                for (std::size_t k = 0; k < 3; ++k)
                    w.F(5 + k, rot[i * 3 + k]);
            w.Append(data);
        }
        if (kept == 0)
            continue;
        step.mLayout.emplace_back(nodal.mTable, nodal.mTCode, 0);
        auto& records = step.mTables[nodal.mTable];
        records.push_back(header.Words(nodal.mTCode, 0, 8, 0, nodal.mThermal ? 1 : 0));
        records.push_back(std::move(data));
    }
    // Complex nodal tables: `<name>_real`/`_imag` and their `_ROT_` rotations, table code
    // 1000 + the real one, format code 2 (real/imaginary).
    for (const Op2wNodal& nodal : kOp2wNodal) {
        const std::string name = nodal.mName;
        if (nodal.mThermal || !rMesh.HasPointData(name + "_real") ||
            !rMesh.HasPointData(name + "_imag"))
            continue;
        std::array<std::vector<double>, 4> parts;  // empty: the mesh has no such part
        static const char* const suffixes[] = {"_real", "_imag", "_ROT_real", "_ROT_imag"};
        for (std::size_t k = 0; k < 4; ++k) {
            if (!rMesh.HasPointData(name + suffixes[k]))
                continue;
            rTaken.insert(name + suffixes[k]);
            const NDArray& a = rMesh.PointData(name + suffixes[k]);
            if (a.Size() != npts * 3)
                op2w_fail("'" + name + suffixes[k] + "' must have 3 components per point");
            const detail::DoubleView view(a);
            parts[k].assign(view.Data(), view.Data() + npts * 3);
            rModel.ToOutput(parts[k]);
        }
        std::string data;
        std::size_t kept = 0;
        std::vector<std::size_t> by_id(npts);
        std::iota(by_id.begin(), by_id.end(), std::size_t{0});
        std::stable_sort(by_id.begin(), by_id.end(), [&rModel](std::size_t A, std::size_t B) {
            return rModel.NodeIds()[A] < rModel.NodeIds()[B];
        });
        for (const std::size_t i : by_id) {
            bool all_nan = true;
            for (std::size_t k = 0; k < 3; ++k)
                all_nan =
                    all_nan && std::isnan(parts[0][i * 3 + k]) && std::isnan(parts[1][i * 3 + k]);
            if (all_nan)
                continue;
            ++kept;
            Op2wRowBuf w(14);
            w.I(0, rModel.NodeIds()[i] * 10 + kOp2wDevice);
            w.I(1, 1);
            for (std::size_t k = 0; k < 3; ++k) {
                w.F(2 + k, parts[0][i * 3 + k]);
                w.F(8 + k, parts[1][i * 3 + k]);
                if (!parts[2].empty())
                    w.F(5 + k, parts[2][i * 3 + k]);
                if (!parts[3].empty())
                    w.F(11 + k, parts[3][i * 3 + k]);
            }
            w.Append(data);
        }
        if (kept == 0)
            continue;
        step.mLayout.emplace_back(nodal.mTable, 1000 + nodal.mTCode, 0);
        auto& records = step.mTables[nodal.mTable];
        records.push_back(header.Words(1000 + nodal.mTCode, 0, 14, 0, 0, 2));
        records.push_back(std::move(data));
    }
    op2w_result_tables(rMesh, rModel, header, step, rTaken);
    op2w_grid_force_tables(rMesh, rModel, header, step, rTaken);
    return step;
}

void op2w_warn_dropped(const Mesh& rMesh, const Op2wModel& rModel,
                       const std::set<std::string>& rTaken) {
    std::set<std::string> names;
    for (const std::string& n : rMesh.PointDataNames())
        if (rTaken.count(n) == 0 && n != "nastran:cp" && n != "nastran:cd" && n != "nastran:nid")
            names.insert(n);
    for (const std::string& n : rMesh.CellDataNames())
        if (rTaken.count(n) == 0 && n != "nastran:eid" && n != "nastran:pid" &&
            n != "nastran:card" && n != "nastran:c1" && n != "nastran:c2" && n != "nastran:x" &&
            n != "nastran:g0" && n != "nastran:pin" && n != "nastran:w" && n != "nastran:k" &&
            n != "nastran:ge" && n != "nastran:s" && n != "nastran:b" && n != "nastran:m" &&
            n != "nastran:conm1")
            names.insert(n);
    for (const std::string& n : rMesh.FieldDataNames())
        if (n != kSequenceTimeKey && n.rfind("nastran:", 0) != 0)
            names.insert(n);
    if (!names.empty()) {
        std::string list;
        for (const std::string& n : names)
            list += (list.empty() ? "" : ", ") + n;
        log::warn("Nastran OP2: arrays without a place in the file are not written: {}", list);
    }
    if (!rModel.Skipped().empty()) {
        std::string list;
        for (const std::string& n : rModel.Skipped())
            list += (list.empty() ? "" : ", ") + n;
        log::warn("Nastran OP2: cell blocks with no element card are not written: {}", list);
    }
    if (!rModel.Unwritten().empty()) {
        std::string list;
        for (const std::string& n : rModel.Unwritten())
            list += (list.empty() ? "" : ", ") + n;
        log::warn("Nastran OP2: element cards not written yet, their elements are left out: {}",
                  list);
    }
}

/// The whole file: header, model tables, then every result table with all steps'
/// blocks, step by step.
std::string op2w_assemble(const Op2wModel& rModel, const std::vector<Op2wStep>& rSteps) {
    std::string out = op2w_file_header();
    out += rModel.Geom1();
    out += rModel.Geom2();
    out += rModel.Ept();
    if (rModel.HasProperties())
        out += rModel.Mpt();
    for (const std::string& table : op2w_table_order()) {
        std::vector<std::string> records;
        for (const Op2wStep& step : rSteps) {
            const auto it = step.mTables.find(table);
            if (it != step.mTables.end())
                records.insert(records.end(), it->second.begin(), it->second.end());
        }
        if (!records.empty())
            out += op2w_table(table, records);
    }
    op2w_marker(out, 0);  // the end of the file
    return out;
}

void op2w_save(const std::string& rPath, const std::string& rBytes) {
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("Could not open file for writing: " + rPath);
    out.write(rBytes.data(), static_cast<std::streamsize>(rBytes.size()));
    out.close();
    if (!out)
        op2w_fail("failed writing " + rPath);
}

}  // namespace

void write_nastran_op2(const std::string& rPath, const Mesh& rMesh) {
    const detail::ProvenanceSlotlessWrite slotless;
    const Op2wModel model(rMesh);
    std::set<std::string> taken;
    const Op2wStep step = op2w_step_tables(rMesh, model, op2w_time_of(rMesh), taken);
    op2w_warn_dropped(rMesh, model, taken);
    op2w_save(rPath, op2w_assemble(model, {step}));
}

struct NastranOp2SeriesWriter::Impl {
    std::string mPath;
    std::optional<Op2wModel> mModel;
    std::vector<Op2wStep> mSteps;
    bool mDone = false;
};

NastranOp2SeriesWriter::NastranOp2SeriesWriter(const std::string& rPath)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
    // Fail now, not when the last step is written, if the path cannot be written.
    auto probe = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!probe)
        throw WriteError("Could not open file for writing: " + rPath);
}

NastranOp2SeriesWriter::~NastranOp2SeriesWriter() = default;
NastranOp2SeriesWriter::NastranOp2SeriesWriter(NastranOp2SeriesWriter&&) noexcept = default;
NastranOp2SeriesWriter& NastranOp2SeriesWriter::operator=(NastranOp2SeriesWriter&&) noexcept =
    default;

void NastranOp2SeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        op2w_fail("series: the writer was moved from");
    Impl& s = *mpImpl;
    const std::size_t index = s.mSteps.size();
    if (!s.mModel) {
        s.mModel.emplace(rMesh);
    } else if (s.mModel->Signature() != op2w_signature(rMesh)) {
        op2w_fail("series: step " + std::to_string(index) +
                  "'s points or cells differ from the first step's; an OP2 file holds one "
                  "model -- write one file per step with '{step}'");
    }
    std::set<std::string> taken;
    Op2wStep step = op2w_step_tables(rMesh, *s.mModel, Time, taken);
    if (s.mSteps.empty()) {
        op2w_warn_dropped(rMesh, *s.mModel, taken);
    } else if (step.mLayout != s.mSteps.front().mLayout) {
        op2w_fail("series: step " + std::to_string(index) +
                  "'s result arrays differ from the first step's; a result table holds one "
                  "block per step");
    }
    // A step with no result table is no step of the file: its key is free.
    if (!step.mLayout.empty())
        for (const Op2wStep& earlier : s.mSteps)
            if (earlier.mKey == step.mKey)
                op2w_fail("series: step " + std::to_string(index) +
                          " has the subcase, analysis and load set, mode or time of an earlier "
                          "step; give it a 'nastran:subcase' or a different time");
    s.mSteps.push_back(std::move(step));
}

std::size_t NastranOp2SeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps.size() : 0;
}

void NastranOp2SeriesWriter::Finalize() {
    if (!mpImpl || mpImpl->mDone)
        return;
    if (mpImpl->mSteps.empty())
        op2w_fail("series: no step was written to " + mpImpl->mPath);
    mpImpl->mDone = true;
    op2w_save(mpImpl->mPath, op2w_assemble(*mpImpl->mModel, mpImpl->mSteps));
}

}  // namespace meshioplusplus
