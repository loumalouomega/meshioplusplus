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
// Ansys MAPDL binary results (.rst, .rth) reader: the mesh, each result set's
// nodal solution, reaction forces and element nodal stresses, strains and
// forces; a distributed solve's partial files merged; a static cyclic model's
// full rotor. The record layout follows Ansys's fdresu.inc as the open reader
// pymapdl-reader (MIT) documents it; see doc/formats/ansys_rst.md. Python twin:
// src/python/meshioplusplus/ansys_rst/_ansys_rst.py.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
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
    std::int64_t mNumNodes = 0, mNumSectors = 0, mKan = 0;
    std::int64_t mCsEls = 0, mCsCord = 0, mGlobalNodes = 0;
    bool mSparseEns = true;  // else ENS holds 11 items per node
    bool mDistributed = false;
    std::uint64_t mPtrGnod = 0;
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
        mKan = h.Int(7);
        mCsEls = h.Int(18);
        mNumSectors = h.Int(20);
        mCsCord = h.Int(21);
        mSparseEns = h.Int(39) != 0;
        mGlobalNodes = h.Int(48);
        mPtrGnod = rst_pointer(h, 49, 50);
        mDistributed = mGlobalNodes != 0 && mGlobalNodes != mNumNodes;
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

// Where an element type's results sit: nodes with element nodal forces (item 63),
// nodes with stresses and strains (item 94: the corners), and the surfaces a
// layered shell stores them for (SHELL181/281 with KEYOPT(8) = 0: bottom, top).
struct RstLayout {
    std::int64_t mNodfor = 0, mNodstr = 0;
    int mLayers = 1;
};

// One results file of a read: its records, headers and model.
struct RstPart {
    RstFile mFile;
    RstResults mResults;
    detail::AnsysModel mModel;
    std::vector<double> mAngles;  // THXY, THYZ, THZX per node, in degrees
    std::map<int, RstLayout> mLayout;
    std::uint64_t mPtrCsy = 0;
    std::int64_t mMaxCsy = 0;
    bool mMapFlag = false;

    RstPart(const std::string& rPath, const ReadOptions& rOptions)
        : mFile(rPath, rOptions), mResults(mFile) {}

    void ReadModel();
    std::pair<std::array<double, 9>, std::array<double, 3>> CoordinateSystem(
        std::int64_t Number) const;
};

void RstPart::ReadModel() {
    const RstFile& file = mFile;
    const RstRecord g = file.Record(mResults.mGeometry);
    const std::int64_t maxety = g.Int(1), nnod = g.Int(3), nelm = g.Int(4);
    const std::uint64_t ptr_ety = rst_pointer(g, 20, 21, 6);
    const std::uint64_t ptr_loc = rst_pointer(g, 26, 27, 8);
    const std::uint64_t ptr_eid = rst_pointer(g, 28, 29, 10);
    mMapFlag = g.Int(64) != 0;
    mPtrCsy = rst_pointer(g, 24, 25, 9);
    mMaxCsy = g.Int(5);

    std::unordered_map<int, std::int64_t> nodelm;
    // Element types: an index record of offsets from ptrETY (in the record after
    // it since 2021R1's mapFlag), each to one type's description.
    const RstRecord table = file.Record(ptr_ety);
    std::vector<std::int64_t> offsets;
    if (mMapFlag) {
        const RstRecord map = file.Record(table.mNext);
        for (std::size_t k = 0; k < map.mValues.size(); ++k)
            offsets.push_back(map.Int(k));
    } else {
        for (std::size_t k = 0; k < table.mValues.size() && k < static_cast<std::size_t>(maxety);
             ++k)
            if (table.Int(k))
                offsets.push_back(table.Int(k));
    }
    for (std::int64_t off : offsets) {
        const RstRecord info = file.Record(ptr_ety + static_cast<std::uint64_t>(off));
        if (info.mValues.size() < 2)
            continue;
        const int slot = static_cast<int>(info.Int(0));
        const int routine = static_cast<int>(info.Int(1));
        mModel.mRoutine[slot] = routine;
        mModel.mKeyopt[slot][1] = static_cast<int>(info.Int(2));
        nodelm[slot] = info.Int(60);
        RstLayout layout;
        layout.mNodfor = info.Int(62);
        layout.mNodstr = info.Int(93);
        layout.mLayers = (routine == 181 || routine == 281) && info.Int(9) == 0 ? 2 : 1;
        mLayout[slot] = layout;
    }

    std::uint64_t ptr = ptr_loc;
    for (std::int64_t k = 0; k < nnod; ++k) {
        const RstRecord node = file.Record(ptr);
        ptr = node.mNext;
        const auto at = [&](std::size_t K) {
            return K < node.mValues.size() ? node.mValues[K] : 0.0;
        };
        mModel.mNodeIds.push_back(static_cast<std::int64_t>(at(0)));
        for (std::size_t d = 1; d <= 3; ++d)
            mModel.mCoords.push_back(at(d));
        for (std::size_t d = 4; d <= 6; ++d)
            mAngles.push_back(at(d));
    }

    if (nelm > 0) {
        // The element index holds 64-bit offsets from ptrEID in 32-bit words.
        const RstRecord index = file.Record(ptr_eid);
        for (std::int64_t e = 0; e < nelm; ++e) {
            const auto k = static_cast<std::size_t>(e);
            const std::int64_t off = static_cast<std::int64_t>(
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(index.Int(2 * k))) |
                (static_cast<std::uint64_t>(index.Int(2 * k + 1)) << 32));
            const std::uint64_t at = ptr_eid + static_cast<std::uint64_t>(off);
            const RstRecord rec = file.Record(at);
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
            mModel.mElements.push_back(std::move(el));
        }
    }

    ptr = rst_pointer(g, 50, 51);
    const std::int64_t n_components = ptr ? g.Int(48) : 0;
    for (std::int64_t c = 0; c < n_components; ++c) {
        const RstRecord comp = file.Record(ptr);
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
        mModel.mComponents.push_back(std::move(out));
    }
}

// (axes, origin) of local coordinate system Number: the rows of axes are its X,
// Y and Z axes in global coordinates.
std::pair<std::array<double, 9>, std::array<double, 3>> RstPart::CoordinateSystem(
    std::int64_t Number) const {
    if (mPtrCsy) {
        const RstRecord table = mFile.Record(mPtrCsy);
        std::vector<std::int64_t> offsets;
        if (mMapFlag) {
            const RstRecord map = mFile.Record(table.mNext);
            for (std::size_t k = 0; k < map.mValues.size(); ++k)
                offsets.push_back(map.Int(k));
        } else {
            for (std::size_t k = 0;
                 k < table.mValues.size() &&
                 k < static_cast<std::size_t>(std::max<std::int64_t>(mMaxCsy, 0));
                 ++k)
                offsets.push_back(table.Int(k));
        }
        for (std::int64_t off : offsets) {
            if (!off)
                continue;
            const RstRecord data = mFile.Record(mPtrCsy + static_cast<std::uint64_t>(off));
            if (data.mValues.size() >= 22 && data.Int(21) == Number) {
                std::array<double, 9> axes{};
                std::array<double, 3> origin{};
                std::copy(data.mValues.begin(), data.mValues.begin() + 9, axes.begin());
                std::copy(data.mValues.begin() + 9, data.mValues.begin() + 12, origin.begin());
                return {axes, origin};
            }
        }
    }
    rst_fail("coordinate system " + std::to_string(Number) +
             " (the cyclic axis) is not in the file");
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

// A DOF's label, for RF_<label> reaction names.
std::string rst_dof_label(std::int64_t Code) {
    static const char* const kLabels[] = {"UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ",
                                          "AX", "AY", "AZ", "VX",   "VY",   "VZ"};
    if (Code >= 1 && Code <= 12)
        return kLabels[Code - 1];
    if (const char* known = rst_scalar_name(Code))
        return known;
    return "DOF" + std::to_string(Code);
}

// The element solution: each element's pointer table (ptrESL) has one entry per
// kind of record, in this order (fdresu.inc).
constexpr std::size_t kRstEnf = 1, kRstEns = 2, kRstEel = 5, kRstEpl = 6, kRstEcr = 7, kRstEth = 8,
                      kRstEul = 9;

// Element nodal tensors: name, table entry, items per node (the first six are
// the components xx yy zz xy yz xz), stress (else engineering strain).
struct RstTensor {
    const char* mName;
    std::size_t mEntry;
    std::size_t mItems;
    bool mStress;
};
constexpr RstTensor kRstTensors[] = {
    {"S", kRstEns, 6, true},     {"EPEL", kRstEel, 7, false}, {"EPPL", kRstEpl, 7, false},
    {"EPCR", kRstEcr, 7, false}, {"EPTH", kRstEth, 8, false},
};
constexpr const char* kRstTop = "@top";  // a layered shell's top surface

bool rst_is_tensor(const std::string& rName, bool* pStress) {
    for (const RstTensor& t : kRstTensors)
        if (rName == t.mName || rName == std::string(t.mName) + kRstTop) {
            if (pStress)
                *pStress = t.mStress;
            return true;
        }
    return false;
}

bool rst_is_vector(const std::string& rName) {
    return rName == "U" || rName == "ROT" || rName == "A" || rName == "V" || rName == "RF" ||
           rName == "RMOM";
}

// Q T Q^T for one tensor xx yy zz xy yz xz (in place); a strain's shear
// components are engineering strains.
void rst_rotate_tensor(double* pT, const double* pQ, bool Stress) {
    const double shear = Stress ? 1.0 : 0.5;
    const double t[9] = {pT[0],         shear * pT[3], shear * pT[5], shear * pT[3], pT[1],
                         shear * pT[4], shear * pT[5], shear * pT[4], pT[2]};
    double qt[9], out[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double v = 0.0;
            for (int k = 0; k < 3; ++k)
                v += pQ[i * 3 + k] * t[k * 3 + j];
            qt[i * 3 + j] = v;
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double v = 0.0;
            for (int k = 0; k < 3; ++k)
                v += qt[i * 3 + k] * pQ[j * 3 + k];
            out[i * 3 + j] = v;
        }
    pT[0] = out[0];
    pT[1] = out[4];
    pT[2] = out[8];
    pT[3] = out[1] / shear;
    pT[4] = out[5] / shear;
    pT[5] = out[2] / shear;
}

// MAPDL's element rotation (THXY, THYZ, THZX in degrees, 3-1-2) as the matrix
// Q = R^T taking a tensor from the element system to the global one.
std::array<double, 9> rst_euler_to_global(const double* pAngles) {
    const double deg = std::acos(-1.0) / 180.0;
    const double c1 = std::cos(pAngles[0] * deg), c2 = std::cos(pAngles[1] * deg),
                 c3 = std::cos(pAngles[2] * deg);
    const double s1 = std::sin(pAngles[0] * deg), s2 = std::sin(pAngles[1] * deg),
                 s3 = std::sin(pAngles[2] * deg);
    const double r[9] = {
        c1 * c3 - s1 * s2 * s3, s1 * c3 + c1 * s2 * s3, -s3 * c2, -s1 * c2, c1 * c2, s2,
        c1 * s3 + s1 * s2 * c3, s1 * s3 - c1 * c3 * s2, c2 * c3};
    return {r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]};
}

// Rodrigues: the rotation by Theta about unit vector pAxis.
std::array<double, 9> rst_axis_rotation(const double* pAxis, double Theta) {
    const double x = pAxis[0], y = pAxis[1], z = pAxis[2];
    const double c = std::cos(Theta), s = std::sin(Theta), t = 1.0 - c;
    return {c + t * x * x,     t * x * y - s * z, t * x * z + s * y,
            t * x * y + s * z, c + t * y * y,     t * y * z - s * x,
            t * x * z - s * y, t * y * z + s * x, c + t * z * z};
}

// The model of one results file, or of a distributed solve's files merged by
// node number (each element lives in one file).
struct RstModel {
    std::vector<std::unique_ptr<RstPart>> mParts;
    std::vector<std::size_t> mOffsets;  // each part's first element in the merged list
    detail::AnsysModel mModel;
    std::map<int, RstLayout> mLayout;
    std::vector<double> mAngles;
    Mesh mMesh;
    std::unordered_map<std::int64_t, std::int64_t> mNodeIndex;
    std::vector<detail::AnsysCellLocation> mCells;
};

// The results file, with its partial siblings when it is the main file of a
// distributed solve (file0.rst beside file1.rst ...).
std::vector<std::unique_ptr<RstPart>> rst_open(const std::string& rPath,
                                               const ReadOptions& rOptions) {
    std::vector<std::unique_ptr<RstPart>> parts;
    parts.push_back(std::make_unique<RstPart>(rPath, rOptions));
    const RstResults& main = parts[0]->mResults;
    if (!main.mDistributed)
        return parts;
    if (!main.mPtrGnod)
        rst_fail("a partial result file of a distributed solve (" + std::to_string(main.mNumNodes) +
                 " of " + std::to_string(main.mGlobalNodes) +
                 " nodes) that is not the main one; read the file whose name ends in 0 (it "
                 "finds the others), or the combined file (RESCOMBINE)");
    const std::filesystem::path path(rPath);
    const std::string stem = path.stem().string();
    const std::string ext = path.extension().string();
    if (stem.empty() || stem.back() != '0')
        rst_fail("the main file of a distributed solve must be named <job>0" + ext +
                 " to find its partial files");
    const std::string job = stem.substr(0, stem.size() - 1);
    for (;;) {
        const std::filesystem::path sibling =
            path.parent_path() / (job + std::to_string(parts.size()) + ext);
        std::error_code ec;
        if (!std::filesystem::exists(sibling, ec))
            break;
        auto part = std::make_unique<RstPart>(sibling.string(), rOptions);
        if (!part->mResults.mDistributed || part->mResults.mSets.size() != main.mSets.size())
            rst_fail(sibling.filename().string() + " is not a partial file of the same solve");
        parts.push_back(std::move(part));
    }
    return parts;
}

RstModel rst_model(std::vector<std::unique_ptr<RstPart>> Parts, bool Lenient) {
    RstModel out;
    out.mParts = std::move(Parts);
    std::unordered_map<std::int64_t, std::size_t> seen;
    std::vector<std::pair<std::string, bool>> component_order;
    std::map<std::pair<std::string, bool>, std::vector<std::int64_t>> components;
    for (auto& part : out.mParts) {
        part->ReadModel();
        const detail::AnsysModel& m = part->mModel;
        out.mOffsets.push_back(out.mModel.mElements.size());
        out.mModel.mElements.insert(out.mModel.mElements.end(), m.mElements.begin(),
                                    m.mElements.end());
        for (const auto& [slot, routine] : m.mRoutine)
            out.mModel.mRoutine.emplace(slot, routine);
        for (const auto& [slot, keys] : m.mKeyopt)
            out.mModel.mKeyopt.emplace(slot, keys);
        for (const auto& [slot, layout] : part->mLayout)
            out.mLayout.emplace(slot, layout);
        for (std::size_t k = 0; k < m.mNodeIds.size(); ++k) {
            if (!seen.emplace(m.mNodeIds[k], out.mModel.mNodeIds.size()).second)
                continue;
            out.mModel.mNodeIds.push_back(m.mNodeIds[k]);
            for (std::size_t d = 0; d < 3; ++d) {
                out.mModel.mCoords.push_back(m.mCoords[k * 3 + d]);
                out.mAngles.push_back(part->mAngles[k * 3 + d]);
            }
        }
        for (const detail::AnsysComponent& c : m.mComponents) {
            const auto key = std::make_pair(c.mName, c.mNodes);
            auto [it, fresh] = components.emplace(key, std::vector<std::int64_t>{});
            if (fresh) {
                component_order.push_back(key);
                it->second = c.mIds;
                continue;
            }
            std::unordered_set<std::int64_t> known(it->second.begin(), it->second.end());
            for (std::int64_t id : c.mIds)
                if (known.insert(id).second)
                    it->second.push_back(id);
        }
    }
    if (out.mParts.size() > 1) {
        const RstPart& main = *out.mParts[0];
        const RstRecord gnod = main.mFile.Record(main.mResults.mPtrGnod);
        std::size_t missing = 0;
        for (std::size_t k = 0; k < gnod.mValues.size(); ++k)
            if (!seen.count(gnod.Int(k)))
                ++missing;
        if (missing)
            rst_fail("the partial files hold " + std::to_string(seen.size()) + " of the " +
                     std::to_string(seen.size() + missing) +
                     " nodes of the distributed solve; one is missing");
    }
    for (const auto& key : component_order) {
        detail::AnsysComponent c;
        c.mName = key.first;
        c.mNodes = key.second;
        c.mIds = components[key];
        out.mModel.mComponents.push_back(std::move(c));
    }
    AnsysInfo info;
    out.mMesh = detail::ansys_build_mesh(out.mModel, Lenient, "Ansys .rst", info, out.mNodeIndex,
                                         out.mCells);
    return out;
}

std::int64_t rst_point_of(const RstModel& rModel, std::int64_t Number) {
    const auto it = rModel.mNodeIndex.find(Number);
    return it == rModel.mNodeIndex.end() ? -1 : it->second;
}

NDArray rst_filled(std::size_t Rows, std::size_t Cols) {
    NDArray out = Cols ? NDArray(DType::Float64, {Rows, Cols}) : NDArray(DType::Float64, {Rows});
    double* p = out.As<double>();
    std::fill(p, p + out.Size(), std::numeric_limits<double>::quiet_NaN());
    return out;
}

// The nodal DOF solution of result set Index, from every part.
void rst_solution(RstModel& rModel, std::size_t Index, const ReadOptions& rOptions) {
    const std::size_t n_points = rModel.mMesh.NumPoints();
    std::map<std::string, NDArray> vectors, scalars;
    for (const auto& part : rModel.mParts) {
        const RstFile& file = part->mFile;
        const RstResults& results = part->mResults;
        const std::uint64_t base = results.mSets[Index].mPointer;
        const RstRecord s = file.Record(base);
        const std::int64_t nnod = s.Int(2), numdof = s.Int(19);
        std::vector<std::int64_t> dofs;
        for (std::int64_t k = 0; k < numdof; ++k)
            dofs.push_back(s.Int(20 + static_cast<std::size_t>(k)));
        const std::int64_t sumdof = numdof + s.Int(97);
        const std::uint64_t ptr_nsl = rst_pointer(s, 104, 105, 10);
        if (!ptr_nsl || numdof <= 0)
            continue;
        const RstRecord values = file.Record(base + ptr_nsl);
        const auto width = static_cast<std::size_t>(sumdof);
        const std::size_t rows = std::min(static_cast<std::size_t>(std::max<std::int64_t>(nnod, 0)),
                                          values.mValues.size() / width);
        std::vector<std::int64_t> numbers;
        if (rows < static_cast<std::size_t>(nnod)) {
            // Not every node has a solution: the next record lists which.
            const RstRecord which = file.Record(values.mNext);
            if (which.mValues.size() < rows)
                rst_fail("result set " + std::to_string(Index + 1) + " has a corrupt node index");
            for (std::size_t r = 0; r < rows; ++r) {
                const std::int64_t k = which.Int(r) - 1;
                if (k < 0 || k >= static_cast<std::int64_t>(results.mNeqv.size()))
                    rst_fail("result set " + std::to_string(Index + 1) +
                             " has a corrupt node index");
                numbers.push_back(results.mNeqv[static_cast<std::size_t>(k)]);
            }
        } else {
            for (std::size_t r = 0; r < rows; ++r)
                numbers.push_back(r < results.mNeqv.size() ? results.mNeqv[r] : -1);
        }
        std::vector<std::int64_t> points(rows, -1);
        for (std::size_t r = 0; r < rows; ++r)
            points[r] = rst_point_of(rModel, numbers[r]);
        const auto cell = [&](std::size_t R, std::size_t K) {
            const double v = values.mValues[R * width + K];
            return std::fabs(v) == kRstUndefined ? std::numeric_limits<double>::quiet_NaN() : v;
        };

        std::unordered_map<std::int64_t, std::size_t> column;
        for (std::size_t k = 0; k < dofs.size(); ++k)
            column.emplace(dofs[k], k);
        static const std::pair<const char*, std::int64_t> kVectors[] = {
            {"U", 1}, {"ROT", 4}, {"A", 7}, {"V", 10}};
        for (const auto& [name, first] : kVectors) {
            const bool any =
                column.count(first) || column.count(first + 1) || column.count(first + 2);
            if (!any || !rOptions.WantsArray(name))
                continue;
            auto it = vectors.find(name);
            if (it == vectors.end())
                it = vectors.emplace(name, rst_filled(n_points, 3)).first;
            double* p = it->second.As<double>();
            for (std::size_t r = 0; r < rows; ++r) {
                if (points[r] < 0)
                    continue;
                double v[3] = {0.0, 0.0, 0.0};
                for (std::int64_t d = 0; d < 3; ++d)
                    if (const auto c = column.find(first + d); c != column.end())
                        v[d] = cell(r, c->second);
                const auto point = static_cast<std::size_t>(points[r]);
                rst_rotate(v, &rModel.mAngles[point * 3]);
                std::copy(v, v + 3, p + point * 3);
            }
        }
        for (std::size_t k = 0; k < dofs.size(); ++k) {
            const std::int64_t code = dofs[k];
            if (code >= 1 && code <= 12)
                continue;
            const char* known = rst_scalar_name(code);
            const std::string name = known ? known : "DOF" + std::to_string(code);
            if (!rOptions.WantsArray(name))
                continue;
            auto it = scalars.find(name);
            if (it == scalars.end())
                it = scalars.emplace(name, rst_filled(n_points, 0)).first;
            double* p = it->second.As<double>();
            for (std::size_t r = 0; r < rows; ++r)
                if (points[r] >= 0)
                    p[points[r]] = cell(r, k);
        }
    }
    for (auto& [name, data] : vectors)
        rModel.mMesh.AddPointData(name, std::move(data));
    for (auto& [name, data] : scalars)
        rModel.mMesh.AddPointData(name, std::move(data));
}

// The reaction forces of result set Index: RF and RMOM (global axes) and
// RF_<DOF label>.
void rst_reactions(RstModel& rModel, std::size_t Index, const ReadOptions& rOptions) {
    const std::size_t n_points = rModel.mMesh.NumPoints();
    std::map<std::string, NDArray> out;
    for (const auto& part : rModel.mParts) {
        const RstFile& file = part->mFile;
        const RstResults& results = part->mResults;
        const std::uint64_t base = results.mSets[Index].mPointer;
        const RstRecord s = file.Record(base);
        const std::int64_t nrf = s.Int(7), numdof = s.Int(19);
        const std::uint64_t ptr_rf = rst_pointer(s, 106, 107, 12);
        if (nrf <= 0 || !ptr_rf || numdof <= 0)
            continue;
        std::vector<std::int64_t> dofs;
        for (std::int64_t k = 0; k < numdof; ++k)
            dofs.push_back(s.Int(20 + static_cast<std::size_t>(k)));
        // (N - 1) * numdof + k: N the node's position in the nodal equivalence
        // table, k the DOF's position in the set's DOF list (1-based).
        const RstRecord table = file.Record(base + ptr_rf);
        const RstRecord values = file.Record(table.mNext);
        const auto n = static_cast<std::size_t>(nrf);
        if (table.mValues.size() < 2 * n || values.mValues.size() < n)
            rst_fail("result set " + std::to_string(Index + 1) +
                     " has a truncated reaction record");
        for (std::size_t r = 0; r < n; ++r) {
            const std::int64_t index = static_cast<std::int64_t>(
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(table.Int(2 * r))) |
                (static_cast<std::uint64_t>(table.Int(2 * r + 1)) << 32));
            const std::int64_t position = (index - 1) / numdof;
            const auto k = static_cast<std::size_t>((index - 1) % numdof);
            if (index < 1 || position >= static_cast<std::int64_t>(results.mNeqv.size()))
                rst_fail("result set " + std::to_string(Index + 1) +
                         " has a corrupt reaction record");
            const std::int64_t point =
                rst_point_of(rModel, results.mNeqv[static_cast<std::size_t>(position)]);
            if (point < 0)
                continue;
            const std::int64_t code = dofs[k];
            std::string name;
            std::size_t component = 0, width = 0;
            if (code >= 1 && code <= 3) {
                name = "RF";
                component = static_cast<std::size_t>(code - 1);
                width = 3;
            } else if (code >= 4 && code <= 6) {
                name = "RMOM";
                component = static_cast<std::size_t>(code - 4);
                width = 3;
            } else {
                name = "RF_" + rst_dof_label(code);
            }
            if (!rOptions.WantsArray(name))
                continue;
            auto it = out.find(name);
            if (it == out.end())
                it = out.emplace(name, rst_filled(n_points, width)).first;
            // A distributed solve's files each hold their share of a reaction at
            // the nodes they have in common: summed.
            double* p = it->second.As<double>();
            if (width) {
                double* row = p + static_cast<std::size_t>(point) * 3;
                for (std::size_t d = 0; d < 3; ++d)
                    if (std::isnan(row[d]))
                        row[d] = 0.0;
                row[component] += values.mValues[r];
            } else {
                p[point] = std::isnan(p[point]) ? values.mValues[r] : p[point] + values.mValues[r];
            }
        }
    }
    for (const char* name : {"RF", "RMOM"}) {
        const auto it = out.find(name);
        if (it == out.end())
            continue;
        double* p = it->second.As<double>();
        for (std::size_t point = 0; point < n_points; ++point)
            if (!std::isnan(p[point * 3]))
                rst_rotate(p + point * 3, &rModel.mAngles[point * 3]);
    }
    for (auto& [name, data] : out)
        rModel.mMesh.AddPointData(name, std::move(data));
}

// Entry Entry of the pointer table at word Table: false when the element has
// no such record. A negative entry -n stands for a record of n zeros that is
// not written.
bool rst_element_record(const RstFile& rFile, std::uint64_t Table, const RstRecord& rPointers,
                        std::size_t Entry, std::vector<double>& rOut) {
    if (Entry >= rPointers.mValues.size())
        return false;
    const std::int64_t ptr = rPointers.Int(Entry);
    if (ptr == 0)
        return false;
    if (ptr < 0) {
        rOut.assign(static_cast<std::size_t>(-ptr), 0.0);
        return true;
    }
    rOut = rFile.Record(Table + static_cast<std::uint64_t>(ptr)).mValues;
    return true;
}

bool rst_no_result_cell(const std::string& rType) {
    return rType == "vertex" || rType == "line" || rType == "line3";
}

// Element nodal stresses and strains (cell data per element node, and averaged
// at the corner nodes as point data) and element nodal forces of result set
// Index.
void rst_elements(RstModel& rModel, std::size_t Index, const ReadOptions& rOptions,
                  std::vector<std::int64_t>& rEnfDofs) {
    Mesh& mesh = rModel.mMesh;
    const std::size_t n_points = mesh.NumPoints();
    const std::size_t n_blocks = mesh.NumCellBlocks();
    std::vector<std::string> types(n_blocks);
    std::vector<std::size_t> rows_of(n_blocks), width_of(n_blocks);
    std::vector<const std::int64_t*> conn(n_blocks);
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto view = mesh.Cells(b);
        types[b] = view.Type();
        rows_of[b] = view.NumCells();
        width_of[b] = view.NodesPerCell();
        conn[b] = view.Conn().template As<std::int64_t>();
    }
    struct Accumulator {
        std::vector<NDArray> mCells;
        std::vector<double> mSum;
        std::vector<std::int64_t> mCount;
    };
    std::map<std::string, Accumulator> tensors;
    const auto accumulator = [&](const std::string& rName) -> Accumulator& {
        auto it = tensors.find(rName);
        if (it != tensors.end())
            return it->second;
        Accumulator a;
        for (std::size_t b = 0; b < n_blocks; ++b) {
            NDArray arr(DType::Float64, {rows_of[b], width_of[b], 6});
            std::fill(arr.As<double>(), arr.As<double>() + arr.Size(),
                      std::numeric_limits<double>::quiet_NaN());
            a.mCells.push_back(std::move(arr));
        }
        a.mSum.assign(n_points * 6, 0.0);
        a.mCount.assign(n_points, 0);
        return tensors.emplace(rName, std::move(a)).first->second;
    };
    bool want_any_tensor = false;
    for (const RstTensor& t : kRstTensors)
        want_any_tensor = want_any_tensor || rOptions.WantsArray(t.mName) ||
                          rOptions.WantsArray(std::string(t.mName) + kRstTop);
    const bool want_enf = rOptions.WantsArray("ENF");
    if (!want_any_tensor && !want_enf)
        return;
    std::vector<NDArray> enf;
    std::size_t enf_width = 0;
    std::vector<double> values, rotation;
    for (std::size_t p = 0; p < rModel.mParts.size(); ++p) {
        const RstPart& part = *rModel.mParts[p];
        const RstFile& file = part.mFile;
        const RstResults& results = part.mResults;
        const std::uint64_t base = results.mSets[Index].mPointer;
        const RstRecord s = file.Record(base);
        const std::uint64_t ptr_esl = rst_pointer(s, 118, 119, 11);
        if (!ptr_esl)
            continue;
        std::vector<std::int64_t> dofs;
        for (std::int64_t k = 0; k < s.Int(19); ++k)
            dofs.push_back(s.Int(20 + static_cast<std::size_t>(k)));
        const RstRecord index = file.Record(base + ptr_esl);
        const std::size_t n_elements = index.mValues.size() / 2;
        for (std::size_t k = 0; k < n_elements; ++k) {
            const std::int64_t off = static_cast<std::int64_t>(
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(index.Int(2 * k))) |
                (static_cast<std::uint64_t>(index.Int(2 * k + 1)) << 32));
            const std::size_t pos = rModel.mOffsets[p] + k;
            if (off == 0 || pos >= rModel.mCells.size() || rModel.mCells[pos].mBlock < 0)
                continue;
            const detail::AnsysCellLocation& loc = rModel.mCells[pos];
            const auto b = static_cast<std::size_t>(loc.mBlock);
            const auto row = static_cast<std::size_t>(loc.mRow);
            if (rst_no_result_cell(types[b]))
                continue;
            RstLayout layout;
            if (const auto it = rModel.mLayout.find(rModel.mModel.mElements[pos].mSlot);
                it != rModel.mLayout.end())
                layout = it->second;
            const std::uint64_t table = base + ptr_esl + static_cast<std::uint64_t>(off);
            const RstRecord pointers = file.Record(table);
            const std::int64_t* cell_nodes = conn[b] + row * width_of[b];
            const auto nodstr = static_cast<std::size_t>(std::max<std::int64_t>(layout.mNodstr, 0));
            bool have_rotation = false;
            for (const RstTensor& t : kRstTensors) {
                if (!rOptions.WantsArray(t.mName) &&
                    !rOptions.WantsArray(std::string(t.mName) + kRstTop))
                    continue;
                const std::size_t items =
                    t.mEntry == kRstEns && !results.mSparseEns ? 11 : t.mItems;
                if (!nodstr || !rst_element_record(file, table, pointers, t.mEntry, values))
                    continue;
                std::size_t n_layers = static_cast<std::size_t>(layout.mLayers);
                if (values.size() < nodstr * n_layers * items) {
                    n_layers = 1;
                    if (values.size() < nodstr * items)
                        continue;
                }
                // The first six items of each node: xx yy zz xy yz xz.
                std::vector<double> nodes(nodstr * n_layers * 6);
                for (std::size_t i = 0; i < nodstr * n_layers; ++i)
                    for (std::size_t c = 0; c < 6; ++c) {
                        const double v = values[i * items + c];
                        nodes[i * 6 + c] = std::fabs(v) == kRstUndefined
                                               ? std::numeric_limits<double>::quiet_NaN()
                                               : v;
                    }
                if (!have_rotation) {
                    if (!rst_element_record(file, table, pointers, kRstEul, rotation))
                        rotation.clear();
                    have_rotation = true;
                }
                if (rotation.size() >= 3) {
                    const std::size_t per_node = rotation.size() / 3;
                    for (std::size_t i = 0; i < nodstr * n_layers; ++i) {
                        const std::size_t j = per_node == 1 ? 0 : i % per_node;
                        const double* angles = &rotation[3 * j];
                        if (angles[0] == 0.0 && angles[1] == 0.0 && angles[2] == 0.0)
                            continue;
                        const auto q = rst_euler_to_global(angles);
                        rst_rotate_tensor(&nodes[i * 6], q.data(), t.mStress);
                    }
                }
                for (std::size_t layer = 0; layer < n_layers; ++layer) {
                    const std::string key = std::string(t.mName) + (layer ? kRstTop : "");
                    if (!rOptions.WantsArray(key))
                        continue;
                    Accumulator& acc = accumulator(key);
                    double* target = acc.mCells[b].As<double>() + row * width_of[b] * 6;
                    for (std::size_t j = 0; j < loc.mSlots.size() && j < width_of[b]; ++j) {
                        const auto slot = static_cast<std::size_t>(loc.mSlots[j]);
                        if (slot >= nodstr)
                            continue;
                        const double* v = &nodes[(layer * nodstr + slot) * 6];
                        std::copy(v, v + 6, target + j * 6);
                        bool finite = true;
                        for (std::size_t c = 0; c < 6; ++c)
                            finite = finite && std::isfinite(v[c]);
                        if (finite) {
                            const auto point = static_cast<std::size_t>(cell_nodes[j]);
                            for (std::size_t c = 0; c < 6; ++c)
                                acc.mSum[point * 6 + c] += v[c];
                            ++acc.mCount[point];
                        }
                    }
                }
            }
            const auto nodfor = static_cast<std::size_t>(std::max<std::int64_t>(layout.mNodfor, 0));
            if (want_enf && nodfor && !dofs.empty() &&
                rst_element_record(file, table, pointers, kRstEnf, values) &&
                values.size() >= nodfor * dofs.size()) {
                if (enf.empty()) {
                    enf_width = dofs.size();
                    rEnfDofs = dofs;
                    for (std::size_t bb = 0; bb < n_blocks; ++bb)
                        enf.push_back([&] {
                            NDArray arr(DType::Float64, {rows_of[bb], width_of[bb], enf_width});
                            std::fill(arr.As<double>(), arr.As<double>() + arr.Size(),
                                      std::numeric_limits<double>::quiet_NaN());
                            return arr;
                        }());
                }
                if (dofs.size() != enf_width)
                    continue;
                double* target = enf[b].As<double>() + row * width_of[b] * enf_width;
                for (std::size_t j = 0; j < loc.mSlots.size() && j < width_of[b]; ++j) {
                    const auto slot = static_cast<std::size_t>(loc.mSlots[j]);
                    if (slot < nodfor)
                        std::copy(&values[slot * enf_width], &values[slot * enf_width] + enf_width,
                                  target + j * enf_width);
                }
            }
        }
    }
    for (auto& [name, acc] : tensors) {
        NDArray mean(DType::Float64, {n_points, 6});
        double* m = mean.As<double>();
        for (std::size_t point = 0; point < n_points; ++point)
            for (std::size_t c = 0; c < 6; ++c)
                m[point * 6 + c] = acc.mCount[point] ? acc.mSum[point * 6 + c] /
                                                           static_cast<double>(acc.mCount[point])
                                                     : std::numeric_limits<double>::quiet_NaN();
        mesh.AddCellData(name, std::move(acc.mCells));
        mesh.AddPointData(name, std::move(mean));
    }
    if (!enf.empty())
        mesh.AddCellData("ENF", std::move(enf));
}

// The full rotor: the base sector's cells (element numbers up to csEls) and
// their points repeated round the cyclic axis, results rotated with them.
// Coincident nodes on the sector boundaries are not merged.
Mesh rst_expand_cyclic(const RstModel& rModel, const std::vector<std::int64_t>& rEnfDofs) {
    const RstPart& part = *rModel.mParts[0];
    const RstResults& results = part.mResults;
    const Mesh& mesh = rModel.mMesh;
    const auto n = static_cast<std::size_t>(results.mNumSectors);
    double axis[3] = {0.0, 0.0, 1.0}, origin[3] = {0.0, 0.0, 0.0};
    if (results.mCsCord > 1) {
        const auto [axes, o] = part.CoordinateSystem(results.mCsCord);
        const double norm = std::sqrt(axes[6] * axes[6] + axes[7] * axes[7] + axes[8] * axes[8]);
        for (std::size_t d = 0; d < 3; ++d) {
            axis[d] = axes[6 + d] / norm;
            origin[d] = o[d];
        }
    }
    const std::size_t n_blocks = mesh.NumCellBlocks();
    std::vector<std::vector<char>> keep(n_blocks);
    for (std::size_t b = 0; b < n_blocks; ++b)
        keep[b].assign(mesh.Cells(b).NumCells(), 0);
    for (std::size_t pos = 0; pos < rModel.mCells.size(); ++pos) {
        const detail::AnsysCellLocation& loc = rModel.mCells[pos];
        if (loc.mBlock >= 0 && rModel.mModel.mElements[pos].mId <= results.mCsEls)
            keep[static_cast<std::size_t>(loc.mBlock)][static_cast<std::size_t>(loc.mRow)] = 1;
    }
    const std::size_t n_points = mesh.NumPoints();
    std::vector<char> used(n_points, 0);
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto view = mesh.Cells(b);
        const std::int64_t* c = view.Conn().template As<std::int64_t>();
        const std::size_t w = view.NodesPerCell();
        for (std::size_t r = 0; r < view.NumCells(); ++r)
            if (keep[b][r])
                for (std::size_t j = 0; j < w; ++j)
                    if (c[r * w + j] >= 0)
                        used[static_cast<std::size_t>(c[r * w + j])] = 1;
    }
    std::vector<std::int64_t> new_point(n_points, -1), old_point;
    for (std::size_t p = 0; p < n_points; ++p)
        if (used[p]) {
            new_point[p] = static_cast<std::int64_t>(old_point.size());
            old_point.push_back(static_cast<std::int64_t>(p));
        }
    const std::size_t n_pts = old_point.size();
    std::vector<std::array<double, 9>> rotations;
    const double pi = std::acos(-1.0);
    for (std::size_t i = 0; i < n; ++i)
        rotations.push_back(
            rst_axis_rotation(axis, 2.0 * pi * static_cast<double>(i) / static_cast<double>(n)));
    const auto rotate = [](const std::array<double, 9>& rQ, const double* pIn, double* pOut) {
        for (std::size_t d = 0; d < 3; ++d)
            pOut[d] = rQ[d * 3] * pIn[0] + rQ[d * 3 + 1] * pIn[1] + rQ[d * 3 + 2] * pIn[2];
    };

    Mesh out;
    {
        const double* src = mesh.Points().template As<double>();
        NDArray points(DType::Float64, {n * n_pts, 3});
        double* dst = points.As<double>();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < n_pts; ++k) {
                const double* p = src + static_cast<std::size_t>(old_point[k]) * 3;
                const double rel[3] = {p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]};
                double* q = dst + (i * n_pts + k) * 3;
                rotate(rotations[i], rel, q);
                for (std::size_t d = 0; d < 3; ++d)
                    q[d] += origin[d];
            }
        out.AssignPoints(std::move(points));
    }
    // Kept blocks and, per kept block, the base mesh's global cell indices.
    std::vector<std::size_t> kept_blocks;
    std::vector<std::vector<std::size_t>> kept_rows;
    std::vector<std::int64_t> new_cell;
    std::int64_t base_cells = 0;
    std::vector<std::int64_t> block_start(n_blocks, 0);
    for (std::size_t b = 0; b < n_blocks; ++b) {
        block_start[b] = base_cells;
        base_cells += static_cast<std::int64_t>(keep[b].size());
    }
    new_cell.assign(static_cast<std::size_t>(base_cells), -1);
    std::vector<std::int64_t> expanded_start;
    std::int64_t expanded = 0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        std::vector<std::size_t> rows;
        for (std::size_t r = 0; r < keep[b].size(); ++r)
            if (keep[b][r])
                rows.push_back(r);
        if (rows.empty())
            continue;
        for (std::size_t k = 0; k < rows.size(); ++k)
            new_cell[static_cast<std::size_t>(block_start[b]) + rows[k]] =
                static_cast<std::int64_t>(k);
        const auto view = mesh.Cells(b);
        const std::size_t w = view.NodesPerCell();
        const std::int64_t* c = view.Conn().template As<std::int64_t>();
        NDArray conn(DType::Int64, {n * rows.size(), w});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < rows.size(); ++k)
                for (std::size_t j = 0; j < w; ++j)
                    dst[(i * rows.size() + k) * w + j] =
                        new_point[static_cast<std::size_t>(c[rows[k] * w + j])] +
                        static_cast<std::int64_t>(i * n_pts);
        out.AddCellBlock(view.Type(), std::move(conn));
        kept_blocks.push_back(b);
        kept_rows.push_back(std::move(rows));
        expanded_start.push_back(expanded);
        expanded += static_cast<std::int64_t>(n * kept_rows.back().size());
    }

    // Values rotated with their sector: vectors, tensors, and ENF's UX..UZ and
    // ROTX..ROTZ columns.
    const auto spin = [&](const std::string& rName, std::size_t I, double* pData, std::size_t Count,
                          std::size_t Width, const std::vector<std::int64_t>* pDofs) {
        const auto& q = rotations[I];
        bool stress = false;
        if (Width == 6 && rst_is_tensor(rName, &stress)) {
            for (std::size_t k = 0; k < Count; ++k)
                rst_rotate_tensor(pData + k * 6, q.data(), stress);
        } else if (Width == 3 && rst_is_vector(rName)) {
            for (std::size_t k = 0; k < Count; ++k) {
                double v[3];
                rotate(q, pData + k * 3, v);
                std::copy(v, v + 3, pData + k * 3);
            }
        } else if (pDofs) {
            for (std::int64_t first : {1, 4}) {
                std::size_t idx[3];
                bool all = true;
                for (std::int64_t d = 0; d < 3; ++d) {
                    const auto it = std::find(pDofs->begin(), pDofs->end(), first + d);
                    all = all && it != pDofs->end();
                    if (all)
                        idx[d] = static_cast<std::size_t>(it - pDofs->begin());
                }
                if (!all)
                    continue;
                for (std::size_t k = 0; k < Count; ++k) {
                    double* row = pData + k * Width;
                    const double in[3] = {row[idx[0]], row[idx[1]], row[idx[2]]};
                    double v[3];
                    rotate(q, in, v);
                    for (std::size_t d = 0; d < 3; ++d)
                        row[idx[d]] = v[d];
                }
            }
        }
    };
    for (const std::string& name : mesh.PointDataNames()) {
        const NDArray& src = mesh.PointData(name);
        const std::size_t width = src.Size() / std::max<std::size_t>(n_points, 1);
        std::vector<std::size_t> shape = src.Shape();
        shape[0] = n * n_pts;
        NDArray dst(src.Dtype(), shape);
        const std::size_t item = dtype_size(src.Dtype()) * width;
        const auto* s = reinterpret_cast<const unsigned char*>(src.Data());
        auto* d = reinterpret_cast<unsigned char*>(dst.Data());
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < n_pts; ++k)
                std::memcpy(d + (i * n_pts + k) * item,
                            s + static_cast<std::size_t>(old_point[k]) * item, item);
            if (src.Dtype() == DType::Float64)
                spin(name, i, dst.As<double>() + i * n_pts * width, n_pts, width, nullptr);
        }
        out.AddPointData(name, std::move(dst));
    }
    for (const std::string& name : mesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        for (std::size_t kb = 0; kb < kept_blocks.size(); ++kb) {
            const std::size_t b = kept_blocks[kb];
            const NDArray& src = mesh.CellData(name, b);
            const std::size_t rows = keep[b].size();
            const std::size_t width = src.Size() / std::max<std::size_t>(rows, 1);
            std::vector<std::size_t> shape = src.Shape();
            shape[0] = n * kept_rows[kb].size();
            NDArray dst(src.Dtype(), shape);
            const std::size_t item = dtype_size(src.Dtype()) * width;
            const auto* s = reinterpret_cast<const unsigned char*>(src.Data());
            auto* d = reinterpret_cast<unsigned char*>(dst.Data());
            const std::size_t count = kept_rows[kb].size();
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t k = 0; k < count; ++k)
                    std::memcpy(d + (i * count + k) * item, s + kept_rows[kb][k] * item, item);
                if (src.Dtype() == DType::Float64) {
                    // Per element node: width = nodes * components.
                    const std::size_t comps = shape.size() >= 3 ? shape.back() : width;
                    spin(name, i, dst.As<double>() + i * count * width, count * (width / comps),
                         comps, name == "ENF" ? &rEnfDofs : nullptr);
                }
            }
            blocks.push_back(std::move(dst));
        }
        out.AddCellData(name, std::move(blocks));
    }
    {
        std::vector<NDArray> sector;
        for (std::size_t kb = 0; kb < kept_blocks.size(); ++kb) {
            const std::size_t count = kept_rows[kb].size();
            NDArray a(DType::Int64, {n * count});
            for (std::size_t i = 0; i < n; ++i)
                std::fill(a.As<std::int64_t>() + i * count, a.As<std::int64_t>() + (i + 1) * count,
                          static_cast<std::int64_t>(i));
            sector.push_back(std::move(a));
        }
        out.AddCellData("ansys:sector", std::move(sector));
    }
    for (const std::string& name : mesh.FieldDataNames())
        out.AddFieldData(name, mesh.FieldData(name));
    {
        NDArray sectors(DType::Int64, {1});
        sectors.As<std::int64_t>()[0] = static_cast<std::int64_t>(n);
        out.AddFieldData("ansys:sectors", std::move(sectors));
    }
    // Block of each base cell among the kept blocks.
    std::vector<std::int64_t> kept_index(n_blocks, -1);
    for (std::size_t kb = 0; kb < kept_blocks.size(); ++kb)
        kept_index[kept_blocks[kb]] = static_cast<std::int64_t>(kb);
    for (std::size_t r = 0; r < mesh.NumRegions(); ++r) {
        const meshioplusplus::Region& region = mesh.Region(r);
        const std::int64_t* e = region.mEntries.template As<std::int64_t>();
        const std::size_t count = region.mEntries.Size();
        std::vector<std::int64_t> entries;
        if (region.mKind == RegionKind::Point) {
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t k = 0; k < count; ++k)
                    if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < n_points &&
                        new_point[static_cast<std::size_t>(e[k])] >= 0)
                        entries.push_back(new_point[static_cast<std::size_t>(e[k])] +
                                          static_cast<std::int64_t>(i * n_pts));
        } else if (region.mKind == RegionKind::Cell) {
            for (std::size_t k = 0; k < count; ++k) {
                if (e[k] < 0 || e[k] >= base_cells)
                    continue;
                const std::int64_t local = new_cell[static_cast<std::size_t>(e[k])];
                if (local < 0)
                    continue;
                const auto b = static_cast<std::size_t>(
                    std::upper_bound(block_start.begin(), block_start.end(), e[k]) -
                    block_start.begin() - 1);
                const auto kb = static_cast<std::size_t>(kept_index[b]);
                const auto size = static_cast<std::int64_t>(kept_rows[kb].size());
                for (std::size_t i = 0; i < n; ++i)
                    entries.push_back(expanded_start[kb] + static_cast<std::int64_t>(i) * size +
                                      local);
            }
        } else {
            continue;
        }
        NDArray arr(DType::Int64, {entries.size()});
        std::copy(entries.begin(), entries.end(), arr.As<std::int64_t>());
        out.AddRegion(meshioplusplus::Region(region.mName, region.mKind, region.mDim, region.mTag,
                                             std::move(arr)));
    }
    return out;
}

Mesh rst_read(const std::string& rPath, const ReadOptions& rOptions, bool Cyclic) {
    std::vector<std::unique_ptr<RstPart>> parts = rst_open(rPath, rOptions);
    const RstResults& main = parts[0]->mResults;
    if (Cyclic) {
        if (main.mNumSectors <= 1)
            rst_fail("not a cyclic-symmetry model; read it as ansys_rst");
        if (main.mKan != 0)
            rst_fail("only static cyclic-symmetry results can be expanded (this is analysis type " +
                     std::to_string(main.mKan) + "); read the base sector as ansys_rst");
    } else if (main.mNumSectors > 1) {
        log::warn(
            "Ansys .rst: a cyclic-symmetry model ({} sectors); only the base sector is read "
            "(ansys_rst_cyclic expands a static one)",
            main.mNumSectors);
    }
    RstModel model = rst_model(std::move(parts), rOptions.mLenient);
    const RstResults& results = model.mParts[0]->mResults;
    std::vector<std::int64_t> enf_dofs;
    const std::size_t n = results.mSets.size();
    if (n == 0) {
        if (rOptions.mTimeStep != 0 && rOptions.mTimeStep != -1)
            rst_fail("time step " + std::to_string(rOptions.mTimeStep) +
                     " is out of range: the file has no result sets");
        return Cyclic ? rst_expand_cyclic(model, enf_dofs) : std::move(model.mMesh);
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
    if (!rOptions.mPointsOnly) {
        rst_solution(model, index, rOptions);
        rst_reactions(model, index, rOptions);
        rst_elements(model, index, rOptions, enf_dofs);
    }
    return Cyclic ? rst_expand_cyclic(model, enf_dofs) : std::move(model.mMesh);
}

MeshMetadata rst_metadata(const std::string& rPath, const ReadOptions& rOptions, bool Cyclic,
                          const char* pFormat) {
    ReadOptions options = rOptions;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    options.mLenient = true;  // a summary skips (and warns about) cell-less elements
    MeshMetadata meta = metadata_from_mesh(rst_read(rPath, options, Cyclic));
    meta.mFellBackToFullRead = true;
    meta.mFormat = pFormat;
    const RstFile file(rPath, rOptions);
    for (const RstSet& s : RstResults(file).mSets)
        meta.mTimeValues.push_back(s.mTime);
    return meta;
}

}  // namespace

Mesh read_ansys_rst(const std::string& rPath, const ReadOptions& rOptions) {
    return rst_read(rPath, rOptions, false);
}

Mesh read_ansys_rst_cyclic(const std::string& rPath, const ReadOptions& rOptions) {
    return rst_read(rPath, rOptions, true);
}

MeshMetadata read_ansys_rst_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    return rst_metadata(rPath, rOptions, false, "ansys_rst");
}

MeshMetadata read_ansys_rst_cyclic_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    return rst_metadata(rPath, rOptions, true, "ansys_rst_cyclic");
}

}  // namespace meshioplusplus
