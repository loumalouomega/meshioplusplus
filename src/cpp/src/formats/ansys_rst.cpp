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
// Ansys MAPDL binary results (.rst, .rth) reader. The record layout follows
// Ansys's fdresu.inc as the open reader pymapdl-reader (MIT) documents it; see
// doc/formats/ansys_rst.md. Python twin:
// src/python/meshioplusplus/ansys_rst/_ansys_rst.py.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ansys_rst.hpp"
#include "meshioplusplus/detail/ansys_model.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

namespace {

// Record flags (the flags word's top byte).
constexpr unsigned kRstBsparse = 0x08;
constexpr unsigned kRstWsparse = 0x10;
constexpr unsigned kRstZlib = 0x20;
constexpr unsigned kRstPrec = 0x40;
constexpr unsigned kRstInt = 0x80;

// MAPDL's "no value".
const double kRstUndefined = std::ldexp(1.0, 100);

[[noreturn]] void rst_fail(const std::string& rWhat) {
    throw ReadError("Ansys .rst: " + rWhat);
}

// One decoded record: its values (integers are exact in a double) and the
// pointer of the record after it.
struct RstRecord {
    std::vector<double> mValues;
    std::uint64_t mNext = 0;

    std::int64_t Int(std::size_t K) const {
        return K < mValues.size() ? static_cast<std::int64_t>(mValues[K]) : 0;
    }
};

class RstFile {
public:
    RstFile(const std::string& rPath, const ReadOptions& rOptions)
        : mSource(rPath, rOptions.mMmap), mWords(mSource.Size() / 4) {}

    std::int32_t Word(std::uint64_t K) const {
        std::int32_t v;
        std::memcpy(&v, mSource.Data() + K * 4, 4);
        return v;
    }

    RstRecord Record(std::uint64_t Ptr) const {
        if (Ptr + 2 > mWords)
            rst_fail("record pointer " + std::to_string(Ptr) + " is outside the file");
        const std::int32_t n = Word(Ptr);
        const unsigned flags = (static_cast<std::uint32_t>(Word(Ptr + 1)) >> 24) & 0xFFu;
        const std::uint64_t end = Ptr + 2 + static_cast<std::uint64_t>(n < 0 ? 0 : n);
        if (n < 0 || end > mWords)
            rst_fail("the record at word " + std::to_string(Ptr) +
                     " runs past the end of the file");
        if (flags & kRstZlib)
            rst_fail("zlib-compressed records are not supported; write the file with /FCOMP,RST,0");
        const bool is_int = (flags & kRstInt) != 0, prec = (flags & kRstPrec) != 0;
        const std::size_t item = is_int ? (prec ? 2 : 4) : (prec ? 4 : 8);
        const char* raw = mSource.Data() + (Ptr + 2) * 4;
        const std::size_t bytes = static_cast<std::size_t>(n) * 4;
        const auto value = [&](std::size_t Offset) -> double {
            if (is_int && prec) {
                std::int16_t v;
                std::memcpy(&v, raw + Offset, 2);
                return v;
            }
            if (is_int) {
                std::int32_t v;
                std::memcpy(&v, raw + Offset, 4);
                return v;
            }
            if (prec) {
                float v;
                std::memcpy(&v, raw + Offset, 4);
                return v;
            }
            double v;
            std::memcpy(&v, raw + Offset, 8);
            return v;
        };
        const auto word = [&](std::size_t K) {
            std::int32_t v;
            std::memcpy(&v, raw + K * 4, 4);
            return v;
        };
        const std::string where = " at word " + std::to_string(Ptr);

        RstRecord out;
        out.mNext = Ptr + static_cast<std::uint64_t>(n) + 3;
        if (flags & kRstBsparse) {
            // A size, a bit mask, then the values of the set bits.
            if (bytes < 8)
                rst_fail("the bit-sparse record" + where + " is truncated");
            const std::int32_t size = word(0);
            const auto bits = static_cast<std::uint32_t>(word(1));
            if (size < 0 || size > 32)
                rst_fail("the bit-sparse record" + where + " has size " + std::to_string(size));
            out.mValues.assign(static_cast<std::size_t>(size), 0.0);
            std::size_t k = 0;
            for (int i = 0; i < size; ++i)
                if ((bits >> i) & 1u) {
                    if (8 + (k + 1) * item > bytes)
                        rst_fail("the bit-sparse record" + where + " is truncated");
                    out.mValues[static_cast<std::size_t>(i)] = value(8 + k * item);
                    ++k;
                }
        } else if (flags & kRstWsparse) {
            // A size, a window count, then windows: an isolated value (loc > 0),
            // a run of values, or a repeated constant.
            if (item == 2)
                rst_fail("the record" + where + " is a windowed-sparse int16 record");
            const std::size_t n_words = bytes / 4;
            if (n_words < 2)
                rst_fail("the windowed-sparse record" + where + " is truncated");
            const std::int32_t size = word(0), n_windows = word(1);
            const std::size_t shift = item / 4;
            out.mValues.assign(static_cast<std::size_t>(std::max(size, 0)), 0.0);
            std::size_t pos = 2;
            const auto take = [&](std::size_t Count) {
                if (pos + shift * Count > n_words)
                    rst_fail("the windowed-sparse record" + where + " is truncated");
                const std::size_t at = pos;
                pos += shift * Count;
                return at;
            };
            for (std::int32_t w = 0; w < n_windows; ++w) {
                if (pos >= n_words)
                    rst_fail("the windowed-sparse record" + where + " is truncated");
                const std::int32_t loc = word(pos++);
                if (loc > 0) {
                    if (loc >= size)
                        rst_fail("the windowed-sparse record" + where + " is corrupt");
                    out.mValues[static_cast<std::size_t>(loc)] = value(take(1) * 4);
                    continue;
                }
                if (pos >= n_words)
                    rst_fail("the windowed-sparse record" + where + " is truncated");
                const std::int64_t start = -static_cast<std::int64_t>(loc);
                const std::int32_t length = word(pos++);
                const std::int64_t count = length < 0 ? -static_cast<std::int64_t>(length) : length;
                if (start + count > size)
                    rst_fail("the windowed-sparse record" + where + " is corrupt");
                if (length > 0) {
                    const std::size_t at = take(static_cast<std::size_t>(count));
                    for (std::int64_t k = 0; k < count; ++k)
                        out.mValues[static_cast<std::size_t>(start + k)] =
                            value((at + static_cast<std::size_t>(k) * shift) * 4);
                } else {
                    const double v = value(take(1) * 4);
                    for (std::int64_t k = 0; k < count; ++k)
                        out.mValues[static_cast<std::size_t>(start + k)] = v;
                }
            }
        } else {
            out.mValues.resize(bytes / item);
            for (std::size_t k = 0; k < out.mValues.size(); ++k)
                out.mValues[k] = value(k * item);
        }
        return out;
    }

private:
    detail::FileSource mSource;
    std::uint64_t mWords;
};

// A 64-bit word pointer split over rHeader[Lo] (low, unsigned) and rHeader[Hi];
// an older file's shorter header has only the 32-bit rHeader[Fallback] (or
// rHeader[Lo] itself).
std::uint64_t rst_pointer(const RstRecord& rHeader, std::size_t Lo, std::size_t Hi,
                          std::size_t Fallback = std::numeric_limits<std::size_t>::max()) {
    const auto low = [&](std::size_t K) {
        return static_cast<std::uint64_t>(static_cast<std::uint32_t>(rHeader.Int(K)));
    };
    if (Hi < rHeader.mValues.size())
        return low(Lo) | (static_cast<std::uint64_t>(rHeader.Int(Hi)) << 32);
    if (Fallback != std::numeric_limits<std::size_t>::max())
        Lo = Fallback;
    return Lo < rHeader.mValues.size() ? low(Lo) : 0;
}

// A name stored as 4-character words, each byte-reversed; cut at the first NUL
// and stripped of spaces and control characters.
std::string rst_name(const RstRecord& rTable, std::size_t Begin, std::size_t End) {
    std::string out;
    for (std::size_t k = Begin; k < End; ++k) {
        const auto v = static_cast<std::uint32_t>(rTable.Int(k));
        for (int shift : {24, 16, 8, 0})
            out += static_cast<char>((v >> shift) & 0xFFu);
    }
    out = out.substr(0, out.find('\0'));
    const auto blank = [](char c) { return static_cast<unsigned char>(c) <= 32; };
    while (!out.empty() && blank(out.back()))
        out.pop_back();
    std::size_t a = 0;
    while (a < out.size() && blank(out[a]))
        ++a;
    return out.substr(a);
}

struct RstSet {
    std::uint64_t mPointer = 0;
    double mTime = 0.0;
    std::int64_t mLoadStep = 0, mSubstep = 0, mCumulative = 0;
};

// The headers and the per-set pointers of a results file.
struct RstResults {
    std::int64_t mNumNodes = 0, mNumSectors = 0;
    std::vector<std::int64_t> mNeqv;
    std::vector<RstSet> mSets;
    std::uint64_t mGeometry = 0;

    explicit RstResults(const RstFile& rFile) {
        const RstRecord standard = rFile.Record(0);
        if (standard.mValues.size() < 2 || standard.Int(0) != 12)
            rst_fail("not a MAPDL results file (its standard header does not name file 12)");
        const RstRecord h = rFile.Record(standard.mNext);
        mNumNodes = h.Int(2);
        const std::int64_t resmax = h.Int(3), nsets = h.Int(8);
        mNumSectors = h.Int(20);
        const std::int64_t global_nnod = h.Int(48);
        if (global_nnod != 0 && global_nnod != mNumNodes)
            rst_fail("a partial result file of a distributed solve (" + std::to_string(mNumNodes) +
                     " of " + std::to_string(global_nnod) +
                     " nodes); read the combined file (RESCOMBINE)");
        const RstRecord neqv = rFile.Record(rst_pointer(h, 14, 45));
        for (std::size_t k = 0; k < neqv.mValues.size() && k < static_cast<std::size_t>(mNumNodes);
             ++k)
            mNeqv.push_back(neqv.Int(k));
        if (nsets > 0) {
            const RstRecord dsi = rFile.Record(rst_pointer(h, 10, 40));
            const RstRecord tim = rFile.Record(rst_pointer(h, 11, 41));
            const RstRecord lsp = rFile.Record(rst_pointer(h, 12, 42));
            for (std::int64_t i = 0; i < nsets; ++i) {
                const auto k = static_cast<std::size_t>(i);
                RstSet s;
                s.mPointer =
                    static_cast<std::uint64_t>(static_cast<std::uint32_t>(dsi.Int(k))) |
                    (static_cast<std::uint64_t>(dsi.Int(static_cast<std::size_t>(resmax) + k))
                     << 32);
                if (k >= tim.mValues.size())
                    rst_fail("the time table is shorter than the " + std::to_string(nsets) +
                             " result sets");
                s.mTime = tim.mValues[k];
                s.mLoadStep = lsp.Int(3 * k);
                s.mSubstep = lsp.Int(3 * k + 1);
                s.mCumulative = lsp.Int(3 * k + 2);
                mSets.push_back(s);
            }
        }
        mGeometry = rst_pointer(h, 15, 46);
    }
};

struct RstModel {
    Mesh mMesh;
    std::unordered_map<std::int64_t, std::int64_t> mNodeIndex;
    std::vector<double> mAngles;  // THXY, THYZ, THZX per node, in degrees
};

RstModel rst_model(const RstFile& rFile, const RstResults& rResults, bool Lenient) {
    const RstRecord g = rFile.Record(rResults.mGeometry);
    const std::int64_t maxety = g.Int(1), nnod = g.Int(3), nelm = g.Int(4);
    const std::uint64_t ptr_ety = rst_pointer(g, 20, 21, 6);
    const std::uint64_t ptr_loc = rst_pointer(g, 26, 27, 8);
    const std::uint64_t ptr_eid = rst_pointer(g, 28, 29, 10);

    detail::AnsysModel model;
    std::unordered_map<int, std::int64_t> nodelm;
    // Element types: an index record of offsets from ptrETY (in the record after
    // it since 2021R1's mapFlag), each to one type's description.
    const RstRecord table = rFile.Record(ptr_ety);
    std::vector<std::int64_t> offsets;
    if (g.Int(64)) {
        const RstRecord map = rFile.Record(table.mNext);
        for (std::size_t k = 0; k < map.mValues.size(); ++k)
            offsets.push_back(map.Int(k));
    } else {
        for (std::size_t k = 0; k < table.mValues.size() && k < static_cast<std::size_t>(maxety);
             ++k)
            if (table.Int(k))
                offsets.push_back(table.Int(k));
    }
    for (std::int64_t off : offsets) {
        const RstRecord info = rFile.Record(ptr_ety + static_cast<std::uint64_t>(off));
        if (info.mValues.size() < 2)
            continue;
        const int slot = static_cast<int>(info.Int(0));
        model.mRoutine[slot] = static_cast<int>(info.Int(1));
        model.mKeyopt[slot][1] = static_cast<int>(info.Int(2));
        nodelm[slot] = info.Int(60);
    }

    std::vector<double> angles;
    std::uint64_t ptr = ptr_loc;
    for (std::int64_t k = 0; k < nnod; ++k) {
        const RstRecord node = rFile.Record(ptr);
        ptr = node.mNext;
        const auto at = [&](std::size_t K) {
            return K < node.mValues.size() ? node.mValues[K] : 0.0;
        };
        model.mNodeIds.push_back(static_cast<std::int64_t>(at(0)));
        for (std::size_t d = 1; d <= 3; ++d)
            model.mCoords.push_back(at(d));
        for (std::size_t d = 4; d <= 6; ++d)
            angles.push_back(at(d));
    }

    if (nelm > 0) {
        // The element index holds 64-bit offsets from ptrEID in 32-bit words.
        const RstRecord index = rFile.Record(ptr_eid);
        for (std::int64_t e = 0; e < nelm; ++e) {
            const auto k = static_cast<std::size_t>(e);
            const std::int64_t off = static_cast<std::int64_t>(
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(index.Int(2 * k))) |
                (static_cast<std::uint64_t>(index.Int(2 * k + 1)) << 32));
            const std::uint64_t at = ptr_eid + static_cast<std::uint64_t>(off);
            const RstRecord rec = rFile.Record(at);
            if (rec.mValues.size() < 10)
                rst_fail("the element record at word " + std::to_string(at) + " is truncated");
            detail::AnsysElement el;
            el.mMat = rec.Int(0);
            el.mSlot = static_cast<int>(rec.Int(1));
            el.mReal = rec.Int(2);
            el.mSecnum = rec.Int(3);
            el.mId = rec.Int(8);
            std::size_t n = rec.mValues.size() - 10;
            if (const auto it = nodelm.find(el.mSlot); it != nodelm.end() && it->second > 0)
                n = std::min(n, static_cast<std::size_t>(it->second));
            for (std::size_t j = 0; j < n; ++j)
                el.mNodes.push_back(rec.Int(10 + j));
            model.mElements.push_back(std::move(el));
        }
    }

    ptr = rst_pointer(g, 50, 51);
    const std::int64_t n_components = ptr ? g.Int(48) : 0;
    for (std::int64_t c = 0; c < n_components; ++c) {
        const RstRecord comp = rFile.Record(ptr);
        ptr = comp.mNext;
        const std::int64_t kind = comp.Int(0);
        if (comp.mValues.size() < 9 || (kind != 1 && kind != 2))
            continue;
        detail::AnsysComponent out;
        out.mName = rst_name(comp, 1, 9);
        out.mNodes = kind == 1;
        // Ids with -last closing a run opened by the value before it.
        for (std::size_t k = 9; k < comp.mValues.size(); ++k) {
            const std::int64_t v = comp.Int(k);
            if (v > 0)
                out.mIds.push_back(v);
            else if (v < 0 && !out.mIds.empty())
                for (std::int64_t id = out.mIds.back() + 1; id <= -v; ++id)
                    out.mIds.push_back(id);
        }
        model.mComponents.push_back(std::move(out));
    }

    RstModel out;
    AnsysInfo info;
    out.mMesh = detail::ansys_build_mesh(model, Lenient, "Ansys .rst", info, out.mNodeIndex);
    out.mAngles = std::move(angles);
    return out;
}

// Nodal to global: a node's axes are the global ones turned about Z by THXY,
// then about the new X by THYZ, then about the newest Y by THZX, so
// v_global = Rz Rx Ry v_nodal.
void rst_rotate(double* pV, const double* pAngles) {
    const double xy = pAngles[0], yz = pAngles[1], zx = pAngles[2];
    if (xy == 0.0 && yz == 0.0 && zx == 0.0)
        return;
    const double deg = std::acos(-1.0) / 180.0;
    double x = pV[0], y = pV[1], z = pV[2];
    double c = std::cos(zx * deg), s = std::sin(zx * deg);
    double nx = c * x + s * z, nz = -s * x + c * z;
    x = nx;
    z = nz;
    c = std::cos(yz * deg);
    s = std::sin(yz * deg);
    double ny = c * y - s * z;
    nz = s * y + c * z;
    y = ny;
    z = nz;
    c = std::cos(xy * deg);
    s = std::sin(xy * deg);
    nx = c * x - s * y;
    ny = s * x + c * y;
    pV[0] = nx;
    pV[1] = ny;
    pV[2] = z;
}

const char* rst_scalar_name(std::int64_t Code) {
    switch (Code) {
        case 16:
            return "WARP";
        case 17:
            return "CONC";
        case 18:
            return "HDSP";
        case 19:
            return "PRES";
        case 20:
            return "TEMP";
        case 21:
            return "VOLT";
        case 22:
            return "MAG";
        case 23:
            return "ENKE";
        case 24:
            return "ENDS";
        case 25:
            return "EMF";
        case 26:
            return "CURR";
        default:
            return nullptr;
    }
}

void rst_solution(const RstFile& rFile, const RstResults& rResults, std::size_t Index,
                  RstModel& rModel, const ReadOptions& rOptions) {
    const std::uint64_t base = rResults.mSets[Index].mPointer;
    const RstRecord s = rFile.Record(base);
    const std::int64_t nnod = s.Int(2), numdof = s.Int(19);
    std::vector<std::int64_t> dofs;
    for (std::int64_t k = 0; k < numdof; ++k)
        dofs.push_back(s.Int(20 + static_cast<std::size_t>(k)));
    const std::int64_t sumdof = numdof + s.Int(97);
    const std::uint64_t ptr_nsl = rst_pointer(s, 104, 105, 10);
    if (!ptr_nsl || numdof <= 0)
        return;
    const RstRecord values = rFile.Record(base + ptr_nsl);
    const auto width = static_cast<std::size_t>(sumdof);
    const std::size_t rows = std::min(static_cast<std::size_t>(std::max<std::int64_t>(nnod, 0)),
                                      values.mValues.size() / width);
    std::vector<std::int64_t> numbers;
    if (rows < static_cast<std::size_t>(nnod)) {
        // Not every node has a solution: the next record lists which.
        const RstRecord which = rFile.Record(values.mNext);
        if (which.mValues.size() < rows)
            rst_fail("result set " + std::to_string(Index + 1) + " has a corrupt node index");
        for (std::size_t r = 0; r < rows; ++r) {
            const std::int64_t k = which.Int(r) - 1;
            if (k < 0 || k >= static_cast<std::int64_t>(rResults.mNeqv.size()))
                rst_fail("result set " + std::to_string(Index + 1) + " has a corrupt node index");
            numbers.push_back(rResults.mNeqv[static_cast<std::size_t>(k)]);
        }
    } else {
        for (std::size_t r = 0; r < rows; ++r)
            numbers.push_back(r < rResults.mNeqv.size() ? rResults.mNeqv[r] : -1);
    }
    std::vector<std::int64_t> points(rows, -1);
    for (std::size_t r = 0; r < rows; ++r)
        if (const auto it = rModel.mNodeIndex.find(numbers[r]); it != rModel.mNodeIndex.end())
            points[r] = it->second;
    const auto cell = [&](std::size_t R, std::size_t K) {
        const double v = values.mValues[R * width + K];
        return std::fabs(v) == kRstUndefined ? std::numeric_limits<double>::quiet_NaN() : v;
    };

    std::unordered_map<std::int64_t, std::size_t> column;
    for (std::size_t k = 0; k < dofs.size(); ++k)
        column.emplace(dofs[k], k);
    const std::size_t n_points = rModel.mMesh.NumPoints();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    static const std::pair<const char*, std::int64_t> kVectors[] = {
        {"U", 1}, {"ROT", 4}, {"A", 7}, {"V", 10}};
    for (const auto& [name, first] : kVectors) {
        const bool any = column.count(first) || column.count(first + 1) || column.count(first + 2);
        if (!any || !rOptions.WantsArray(name))
            continue;
        NDArray out(DType::Float64, {n_points, 3});
        double* p = out.As<double>();
        std::fill(p, p + n_points * 3, nan);
        for (std::size_t r = 0; r < rows; ++r) {
            if (points[r] < 0)
                continue;
            double v[3] = {0.0, 0.0, 0.0};
            for (std::int64_t d = 0; d < 3; ++d)
                if (const auto it = column.find(first + d); it != column.end())
                    v[d] = cell(r, it->second);
            const auto point = static_cast<std::size_t>(points[r]);
            if (point * 3 + 2 < rModel.mAngles.size())
                rst_rotate(v, &rModel.mAngles[point * 3]);
            std::copy(v, v + 3, p + point * 3);
        }
        rModel.mMesh.AddPointData(name, std::move(out));
    }
    for (std::size_t k = 0; k < dofs.size(); ++k) {
        const std::int64_t code = dofs[k];
        if (code >= 1 && code <= 12)
            continue;
        const char* known = rst_scalar_name(code);
        const std::string name = known ? known : "DOF" + std::to_string(code);
        if (!rOptions.WantsArray(name))
            continue;
        NDArray out(DType::Float64, {n_points});
        double* p = out.As<double>();
        std::fill(p, p + n_points, nan);
        for (std::size_t r = 0; r < rows; ++r)
            if (points[r] >= 0)
                p[points[r]] = cell(r, k);
        rModel.mMesh.AddPointData(name, std::move(out));
    }
}

}  // namespace

Mesh read_ansys_rst(const std::string& rPath, const ReadOptions& rOptions) {
    const RstFile file(rPath, rOptions);
    const RstResults results(file);
    if (results.mNumSectors > 1)
        log::warn("Ansys .rst: a cyclic-symmetry model ({} sectors); only the base sector is read",
                  results.mNumSectors);
    RstModel model = rst_model(file, results, rOptions.mLenient);
    const std::size_t n = results.mSets.size();
    if (n == 0) {
        if (rOptions.mTimeStep != 0 && rOptions.mTimeStep != -1)
            rst_fail("time step " + std::to_string(rOptions.mTimeStep) +
                     " is out of range: the file has no result sets");
        return std::move(model.mMesh);
    }
    const std::size_t index = rOptions.ResolveTimeStep(n);
    const RstSet& set = results.mSets[index];
    const auto scalar = [](DType T, double V) {
        NDArray a(T, {1});
        if (T == DType::Float64)
            a.As<double>()[0] = V;
        else
            a.As<std::int64_t>()[0] = static_cast<std::int64_t>(V);
        return a;
    };
    Mesh& mesh = model.mMesh;
    mesh.AddFieldData("meshio:time", scalar(DType::Float64, set.mTime));
    mesh.AddFieldData("ansys:load_step", scalar(DType::Int64, static_cast<double>(set.mLoadStep)));
    mesh.AddFieldData("ansys:substep", scalar(DType::Int64, static_cast<double>(set.mSubstep)));
    mesh.AddFieldData("ansys:cumulative",
                      scalar(DType::Int64, static_cast<double>(set.mCumulative)));
    if (!rOptions.mPointsOnly)
        rst_solution(file, results, index, model, rOptions);
    return std::move(model.mMesh);
}

MeshMetadata read_ansys_rst_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    ReadOptions options = rOptions;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    options.mLenient = true;  // a summary skips (and warns about) cell-less elements
    MeshMetadata meta = metadata_from_mesh(read_ansys_rst(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "ansys_rst";
    const RstFile file(rPath, rOptions);
    for (const RstSet& s : RstResults(file).mSets)
        meta.mTimeValues.push_back(s.mTime);
    return meta;
}

}  // namespace meshioplusplus
