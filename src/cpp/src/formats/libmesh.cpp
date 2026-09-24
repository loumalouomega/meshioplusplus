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
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/libmesh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/binary_stream.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kLmWhat = "libMesh";

// --- element types ----------------------------------------------------------------

// libMesh's ElemType (include/enums/enum_elem_type.h): node count, the meshio++
// type it reads as, how many of its nodes that type keeps, and its base shape
// for the side tables.
enum class LmShape { None, Point, Edge, Tri, Quad, Tet, Hex, Prism, Pyramid };

struct LmType {
    int mNodes;         // nodes in the file (0: not a readable type)
    const char* mCell;  // meshio++ type, nullptr: skipped
    int mKeep;          // leading file nodes the meshio++ type keeps
    LmShape mShape;
    const char* mName;  // libMesh's name, for messages
};

const LmType* lm_type(std::uint64_t Code) {
    static const LmType types[] = {
        {2, "line", 2, LmShape::Edge, "EDGE2"},                // 0
        {3, "line3", 3, LmShape::Edge, "EDGE3"},               // 1
        {4, "line4", 4, LmShape::Edge, "EDGE4"},               // 2
        {3, "triangle", 3, LmShape::Tri, "TRI3"},              // 3
        {6, "triangle6", 6, LmShape::Tri, "TRI6"},             // 4
        {4, "quad", 4, LmShape::Quad, "QUAD4"},                // 5
        {8, "quad8", 8, LmShape::Quad, "QUAD8"},               // 6
        {9, "quad9", 9, LmShape::Quad, "QUAD9"},               // 7
        {4, "tetra", 4, LmShape::Tet, "TET4"},                 // 8
        {10, "tetra10", 10, LmShape::Tet, "TET10"},            // 9
        {8, "hexahedron", 8, LmShape::Hex, "HEX8"},            // 10
        {20, "hexahedron20", 20, LmShape::Hex, "HEX20"},       // 11
        {27, "hexahedron27", 27, LmShape::Hex, "HEX27"},       // 12
        {6, "wedge", 6, LmShape::Prism, "PRISM6"},             // 13
        {15, "wedge15", 15, LmShape::Prism, "PRISM15"},        // 14
        {18, "wedge18", 18, LmShape::Prism, "PRISM18"},        // 15
        {5, "pyramid", 5, LmShape::Pyramid, "PYRAMID5"},       // 16
        {13, "pyramid13", 13, LmShape::Pyramid, "PYRAMID13"},  // 17
        {14, "pyramid14", 14, LmShape::Pyramid, "PYRAMID14"},  // 18
        {2, nullptr, 0, LmShape::None, "INFEDGE2"},            // 19
        {4, nullptr, 0, LmShape::None, "INFQUAD4"},            // 20
        {6, nullptr, 0, LmShape::None, "INFQUAD6"},            // 21
        {8, nullptr, 0, LmShape::None, "INFHEX8"},             // 22
        {16, nullptr, 0, LmShape::None, "INFHEX16"},           // 23
        {18, nullptr, 0, LmShape::None, "INFHEX18"},           // 24
        {6, nullptr, 0, LmShape::None, "INFPRISM6"},           // 25
        {12, nullptr, 0, LmShape::None, "INFPRISM12"},         // 26
        {1, "vertex", 1, LmShape::Point, "NODEELEM"},          // 27
        {0, nullptr, 0, LmShape::None, "REMOTEELEM"},          // 28
        {3, "triangle", 3, LmShape::Tri, "TRI3SUBDIVISION"},   // 29
        {3, "triangle", 3, LmShape::Tri, "TRISHELL3"},         // 30
        {4, "quad", 4, LmShape::Quad, "QUADSHELL4"},           // 31
        {8, "quad8", 8, LmShape::Quad, "QUADSHELL8"},          // 32
        {7, "triangle7", 7, LmShape::Tri, "TRI7"},             // 33
        {14, "tetra10", 10, LmShape::Tet, "TET14"},            // 34
        {20, "wedge18", 18, LmShape::Prism, "PRISM20"},        // 35
        {21, "wedge18", 18, LmShape::Prism, "PRISM21"},        // 36
        {18, "pyramid14", 14, LmShape::Pyramid, "PYRAMID18"},  // 37
        {9, "quad9", 9, LmShape::Quad, "QUADSHELL9"},          // 38
        {0, nullptr, 0, LmShape::None, "C0POLYGON"},           // 39
        {0, nullptr, 0, LmShape::None, "C0POLYHEDRON"},        // 40
    };
    return Code < std::size(types) ? &types[Code] : nullptr;
}

// Corner nodes of each side, libMesh's side numbering (Hex8::side_nodes_map ...).
// Corners are in the same slots in libMesh and meshio++ for every shape.
const std::vector<std::vector<int>>& lm_sides(LmShape Shape) {
    static const std::vector<std::vector<int>> none;
    static const std::vector<std::vector<int>> tri = {{0, 1}, {1, 2}, {2, 0}};
    static const std::vector<std::vector<int>> quad = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const std::vector<std::vector<int>> tet = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    static const std::vector<std::vector<int>> hex = {{0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5},
                                                      {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
    static const std::vector<std::vector<int>> prism = {
        {0, 2, 1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {3, 4, 5}};
    static const std::vector<std::vector<int>> pyramid = {
        {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}, {0, 3, 2, 1}};
    switch (Shape) {
        case LmShape::Tri:
            return tri;
        case LmShape::Quad:
            return quad;
        case LmShape::Tet:
            return tet;
        case LmShape::Hex:
            return hex;
        case LmShape::Prism:
            return prism;
        case LmShape::Pyramid:
            return pyramid;
        default:
            return none;
    }
}

// --- the value stream ---------------------------------------------------------------

// libMesh's `Xdr` in READ (ASCII) or DECODE (XDR) mode. ASCII follows
// xdr_cxx.C: a scalar is `>> value` and the rest of its line is dropped (the
// `\t # comment`), a string is the rest of the current line up to the first tab,
// a vector is a scalar length, its items, and the rest of that line, and a
// stream is bare whitespace-separated values.
class LmStream {
public:
    LmStream(const std::string& rText, bool Xdr)
        : mText(rText), mXdr(Xdr), mBin(rText.data(), rText.size(), true, kLmWhat) {}

    bool Xdr() const { return mXdr; }

    std::string String() {
        if (mXdr) {
            const std::uint32_t n = mBin.U32();
            std::string s = mBin.Bytes(n);
            mBin.Skip((4 - n % 4) % 4);
            const std::size_t nul = s.find('\0');
            if (nul != std::string::npos)
                s.resize(nul);
            return s;
        }
        const std::size_t eol = LineEnd();
        std::string s = mText.substr(mPos, eol - mPos);
        mPos = eol < mText.size() ? eol + 1 : eol;
        const std::size_t tab = s.find('\t');
        if (tab != std::string::npos)
            s.resize(tab);
        if (!s.empty() && s.back() == '\r')
            s.pop_back();
        return s;
    }

    // A `data()` integer of `Width` bytes (4 or 8).
    std::int64_t Scalar(int Width) {
        if (mXdr)
            return BinInt(Width);
        const std::int64_t v = Int();
        SkipLine();
        return v;
    }

    // A `data_stream()` integer of `Width` bytes.
    std::int64_t StreamInt(int Width) { return mXdr ? BinInt(Width) : Int(); }

    double StreamReal() {
        if (mXdr)
            return mBin.F64();
        const std::string_view t = Token();
        const std::string tok(t);
        const char* end = nullptr;
        const double v = detail::parse_double(tok.c_str(), end);
        if (end != tok.c_str() + tok.size()) {
            std::string low;
            for (char c : tok)
                low += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (low.find("nan") != std::string::npos)
                return std::numeric_limits<double>::quiet_NaN();
            Fail("bad real '" + tok + "'");
        }
        return v;
    }

    // `data(std::vector<T>)`: a 4-byte length, then its items.
    std::vector<std::int64_t> IntVector(int Width) {
        const std::int64_t n = Scalar(4);
        std::vector<std::int64_t> out;
        for (std::int64_t k = 0; k < n; ++k)
            out.push_back(StreamInt(Width));
        if (!mXdr)
            SkipLine();
        return out;
    }

    std::vector<std::string> StringVector() {
        const std::int64_t n = Scalar(4);
        std::vector<std::string> out;
        for (std::int64_t k = 0; k < n; ++k) {
            if (mXdr) {
                out.push_back(String());
            } else {
                out.emplace_back(Token());
            }
        }
        if (!mXdr)
            SkipLine();
        return out;
    }

    [[noreturn]] void Fail(const std::string& rWhat) const {
        throw ReadError(std::string("libMesh: ") + rWhat +
                        (mXdr ? " (byte " + std::to_string(mBin.Offset()) + ")"
                              : " (line " + std::to_string(LineNumber()) + ")"));
    }

private:
    std::int64_t BinInt(int Width) {
        return Width == 8 ? static_cast<std::int64_t>(mBin.U64())
                          : static_cast<std::int64_t>(mBin.U32());
    }

    std::size_t LineEnd() const {
        const std::size_t eol = mText.find('\n', mPos);
        return eol == std::string::npos ? mText.size() : eol;
    }

    void SkipLine() {
        const std::size_t eol = LineEnd();
        mPos = eol < mText.size() ? eol + 1 : eol;
    }

    std::string_view Token() {
        while (mPos < mText.size() && std::isspace(static_cast<unsigned char>(mText[mPos])))
            ++mPos;
        if (mPos >= mText.size())
            Fail("file ends early");
        const std::size_t b = mPos;
        while (mPos < mText.size() && !std::isspace(static_cast<unsigned char>(mText[mPos])))
            ++mPos;
        return std::string_view(mText).substr(b, mPos - b);
    }

    std::int64_t Int() {
        const std::string_view t = Token();
        std::int64_t v = 0;
        const auto [p, ec] = std::from_chars(t.data(), t.data() + t.size(), v);
        if (ec != std::errc() || p != t.data() + t.size())
            Fail("bad integer '" + std::string(t) + "'");
        return v;
    }

    std::size_t LineNumber() const {
        return 1 + static_cast<std::size_t>(std::count(
                       mText.begin(),
                       mText.begin() + static_cast<std::ptrdiff_t>(std::min(mPos, mText.size())),
                       '\n'));
    }

    const std::string& mText;
    bool mXdr;
    std::size_t mPos = 0;
    detail::ByteCursor mBin;
};

// libMesh's version_at_least_* tests: substring checks for the known versions.
bool lm_has_any(const std::string& rVersion, std::initializer_list<const char*> Versions) {
    for (const char* v : Versions)
        if (rVersion.find(v) != std::string::npos)
            return true;
    return false;
}

struct LmElement {
    const LmType* mType;
    std::uint64_t mCode;
    std::int64_t mParent;  // file index, -1 at level 0
    std::int64_t mSubdomain;
    std::int64_t mPLevel;
    int mLevel;
    std::vector<std::int64_t> mNodes;
};

struct LmBoundary {
    std::int64_t mElem, mSide, mId;
};

struct LmFile {
    std::vector<LmElement> mElements;
    std::vector<double> mCoords;  // 3 per node id
    std::map<std::int64_t, std::string> mSubdomainNames, mSidesetNames, mNodesetNames;
    std::vector<LmBoundary> mSides;
    std::vector<std::pair<std::int64_t, std::int64_t>> mNodesets;  // (node, id)
    std::size_t mSkippedEdgeBcs = 0;
    bool mInlinePLevel = false;
};

std::map<std::int64_t, std::string> lm_name_map(LmStream& rIo, int HeaderWidth) {
    std::map<std::int64_t, std::string> out;
    const std::int64_t n = rIo.Scalar(HeaderWidth);
    if (n == 0)
        return out;
    const std::vector<std::int64_t> ids = rIo.IntVector(HeaderWidth);
    const std::vector<std::string> names = rIo.StringVector();
    for (std::size_t k = 0; k < ids.size() && k < names.size(); ++k)
        out.emplace(ids[k], names[k]);
    return out;
}

LmFile lm_parse(const std::string& rText, bool Xdr) {
    LmStream io(rText, Xdr);
    LmFile f;
    const std::string version = io.String();
    if (version.find("libMesh") == std::string::npos)
        throw ReadError("libMesh: '" + version +
                        "' is a legacy (pre-libMesh) mesh file, which libMesh itself no longer "
                        "reads");
    const bool v092 = lm_has_any(version, {"0.9.2", "0.9.6", "1.1.0", "1.3.0", "1.8.0"});
    const bool v096 = lm_has_any(version, {"0.9.6", "1.1.0", "1.3.0", "1.8.0"});
    const bool v110 = lm_has_any(version, {"1.1.0", "1.3.0", "1.8.0"});
    const bool v130 = lm_has_any(version, {"1.3.0", "1.8.0"});
    const bool v180 = lm_has_any(version, {"1.8.0"});
    const int hw = v130 ? 8 : 4;  // header integers: uint64 from 1.3.0 on

    const std::int64_t n_elem = io.Scalar(hw);
    const std::int64_t n_nodes = io.Scalar(hw);
    const std::string bc_file = io.String();
    const std::string sid_file = io.String();
    const std::string pid_file = io.String();
    const std::string pl_file = io.String();
    if (n_elem < 0 || n_nodes < 0)
        io.Fail("negative element or node count");

    std::int64_t sizes[8] = {8, 0, 0, 0, 0, 0, 0, 0};  // type uid pid sid p eid side bid
    if (v092)
        for (std::int64_t& s : sizes)
            s = io.Scalar(hw);
    const std::int64_t field_width = sizes[0];
    // Pre-1.3.0 files were written with 32-bit connectivity whatever "type size"
    // says (XdrIO::read); from 1.3.0 on "type size" is the width.
    const int tw = (!v130 || field_width == 4) ? 4 : 8;
    const bool read_uid = v092 && sizes[1] != 0;

    std::size_t n_elem_ints = 0, n_node_ints = 0;
    if (v180) {
        io.Scalar(hw);  // extra integer size
        n_node_ints = io.StringVector().size();
        n_elem_ints = io.StringVector().size();
        const std::vector<std::int64_t> codes = io.IntVector(tw);
        for (std::size_t k = 0; k < codes.size(); ++k)
            io.IntVector(tw);
    }
    if (v092)
        f.mSubdomainNames = lm_name_map(io, hw);
    if (n_elem == 0)
        return f;

    const bool read_pid = pid_file == ".";
    const bool read_sid = sid_file == ".";
    const bool read_p = pl_file == ".";
    f.mInlinePLevel = read_p;
    f.mElements.reserve(static_cast<std::size_t>(n_elem));
    std::int64_t at_level = 0, done_at_level = 0;
    int level = -1;
    for (std::int64_t e = 0; e < n_elem; ++e, ++done_at_level) {
        if (done_at_level == at_level) {
            at_level = io.Scalar(tw);
            done_at_level = 0;
            ++level;
        }
        LmElement el{};
        el.mCode = static_cast<std::uint64_t>(io.StreamInt(tw));
        el.mType = lm_type(el.mCode);
        if (!el.mType || el.mType->mNodes == 0)
            io.Fail("element " + std::to_string(e) + " has unsupported type " +
                    std::to_string(el.mCode) +
                    (el.mType ? std::string(" (") + el.mType->mName + ")" : std::string()));
        if (read_uid)
            io.StreamInt(tw);
        el.mParent = -1;
        if (level > 0) {
            el.mParent = io.StreamInt(tw);
            if (el.mParent < 0 || el.mParent >= e)
                io.Fail("element " + std::to_string(e) + " names parent " +
                        std::to_string(el.mParent) + ", which is not an earlier element");
        }
        if (read_pid)
            io.StreamInt(tw);
        el.mSubdomain = read_sid ? io.StreamInt(tw) : 0;
        el.mPLevel = read_p ? io.StreamInt(tw) : 0;
        el.mLevel = level;
        el.mNodes.resize(static_cast<std::size_t>(el.mType->mNodes));
        for (std::int64_t& n : el.mNodes) {
            n = io.StreamInt(tw);
            if (n < 0 || n >= n_nodes)
                io.Fail("element " + std::to_string(e) + " names node " + std::to_string(n) +
                        " of " + std::to_string(n_nodes));
        }
        for (std::size_t k = 0; k < n_elem_ints; ++k)
            io.StreamInt(tw);
        f.mElements.push_back(std::move(el));
    }

    f.mCoords.resize(3 * static_cast<std::size_t>(n_nodes));
    for (double& c : f.mCoords)
        c = io.StreamReal();
    if (v096 && io.Scalar(4) != 0)  // "presence of unique ids" (an unsigned short)
        for (std::int64_t k = 0; k < n_nodes; ++k)
            io.StreamInt(field_width == 8 ? 8 : 4);
    for (std::size_t k = 0; k < n_node_ints * static_cast<std::size_t>(n_nodes); ++k)
        io.StreamInt(tw);

    if (bc_file == "n/a")
        return f;
    auto read_triples = [&](std::vector<LmBoundary>* pOut) {
        std::map<std::int64_t, std::string> names;
        if (v092)
            names = lm_name_map(io, hw);
        const std::int64_t n = io.Scalar(hw);
        for (std::int64_t k = 0; k < n; ++k) {
            LmBoundary b{};
            b.mElem = io.StreamInt(tw);
            b.mSide = io.StreamInt(tw);
            b.mId = io.StreamInt(tw);
            if (pOut)
                pOut->push_back(b);
        }
        return std::make_pair(names, static_cast<std::size_t>(n));
    };
    f.mSidesetNames = read_triples(&f.mSides).first;
    if (v092) {
        f.mNodesetNames = lm_name_map(io, hw);
        const std::int64_t n = io.Scalar(hw);
        for (std::int64_t k = 0; k < n; ++k) {
            const std::int64_t node = io.StreamInt(tw);
            const std::int64_t id = io.StreamInt(tw);
            f.mNodesets.emplace_back(node, id);
        }
    }
    if (v110) {
        f.mSkippedEdgeBcs += read_triples(nullptr).second;  // edge
        f.mSkippedEdgeBcs += read_triples(nullptr).second;  // shell face
    }
    return f;
}

// --- geometry for carrying a side down to refined children ------------------------

using LmPoint = std::array<double, 3>;

LmPoint lm_sub(const LmPoint& a, const LmPoint& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
double lm_dot(const LmPoint& a, const LmPoint& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
LmPoint lm_cross(const LmPoint& a, const LmPoint& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Whether `rQ` lies on the (flat) side whose corners are `rP`: on the segment
// for two corners, else on the plane and inside the triangle or quad.
bool lm_on_side(const std::vector<LmPoint>& rP, const LmPoint& rQ) {
    double scale = 0.0;
    for (std::size_t k = 0; k < rP.size(); ++k) {
        const LmPoint d = lm_sub(rP[(k + 1) % rP.size()], rP[k]);
        scale = std::max(scale, std::sqrt(lm_dot(d, d)));
    }
    const double tol = 1e-6 * (scale > 0 ? scale : 1.0);
    if (rP.size() == 2) {
        const LmPoint d = lm_sub(rP[1], rP[0]);
        const LmPoint q = lm_sub(rQ, rP[0]);
        const double len2 = lm_dot(d, d);
        if (len2 == 0.0)
            return false;
        const double t = lm_dot(q, d) / len2;
        const LmPoint off = lm_cross(q, d);
        return t >= -1e-6 && t <= 1 + 1e-6 && std::sqrt(lm_dot(off, off) / len2) <= tol;
    }
    auto in_triangle = [&](const LmPoint& a, const LmPoint& b, const LmPoint& c) {
        const LmPoint n = lm_cross(lm_sub(b, a), lm_sub(c, a));
        const double n2 = lm_dot(n, n);
        if (n2 == 0.0)
            return false;
        if (std::abs(lm_dot(lm_sub(rQ, a), n)) / std::sqrt(n2) > tol)
            return false;
        const double eps = -1e-6 * n2;
        return lm_dot(lm_cross(lm_sub(b, a), lm_sub(rQ, a)), n) >= eps &&
               lm_dot(lm_cross(lm_sub(c, b), lm_sub(rQ, b)), n) >= eps &&
               lm_dot(lm_cross(lm_sub(a, c), lm_sub(rQ, c)), n) >= eps;
    };
    if (rP.size() == 3)
        return in_triangle(rP[0], rP[1], rP[2]);
    return in_triangle(rP[0], rP[1], rP[2]) || in_triangle(rP[0], rP[2], rP[3]);
}

NDArray lm_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride = 1) {
    NDArray a(DType::Int64, Stride == 1 ? std::vector<std::size_t>{rIds.size()}
                                        : std::vector<std::size_t>{rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

Mesh read_libmesh(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("libMesh: cannot open " + rPath);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    // XDR starts with the version string's 4-byte big-endian length; ASCII with
    // the version text itself.
    bool xdr = false;
    if (text.size() >= 8 && text.compare(0, 7, "libMesh") != 0) {
        const auto b = [&](std::size_t k) { return static_cast<unsigned char>(text[k]); };
        const std::uint32_t len = (std::uint32_t(b(0)) << 24) | (std::uint32_t(b(1)) << 16) |
                                  (std::uint32_t(b(2)) << 8) | std::uint32_t(b(3));
        xdr = len > 0 && len < 256 && text.size() >= 4 + len && text.compare(4, 7, "libMesh") == 0;
    }
    const LmFile f = lm_parse(text, xdr);

    // Active elements: the ones no element names as its parent.
    const std::size_t ne = f.mElements.size();
    std::vector<bool> active(ne, true);
    std::vector<std::vector<std::size_t>> children(ne);
    int max_level = 0;
    for (std::size_t e = 0; e < ne; ++e) {
        const LmElement& el = f.mElements[e];
        max_level = std::max(max_level, el.mLevel);
        if (el.mParent >= 0) {
            active[static_cast<std::size_t>(el.mParent)] = false;
            children[static_cast<std::size_t>(el.mParent)].push_back(e);
        }
    }

    // Points: every node id an element uses; unused ids (NaN in the file) go.
    const std::size_t nn = f.mCoords.size() / 3;
    std::vector<std::int64_t> node_index(nn, -1);
    for (const LmElement& el : f.mElements)
        for (std::int64_t n : el.mNodes)
            node_index[static_cast<std::size_t>(n)] = 0;
    std::vector<std::int64_t> kept;
    for (std::size_t n = 0; n < nn; ++n)
        if (node_index[n] == 0) {
            node_index[n] = static_cast<std::int64_t>(kept.size());
            kept.push_back(static_cast<std::int64_t>(n));
        }
    Mesh mesh;
    NDArray points(DType::Float64, {kept.size(), 3});
    double* pts = points.As<double>();
    for (std::size_t p = 0; p < kept.size(); ++p)
        for (std::size_t d = 0; d < 3; ++d)
            pts[3 * p + d] = f.mCoords[3 * static_cast<std::size_t>(kept[p]) + d];
    mesh.AssignPoints(std::move(points));
    if (kept.size() != nn)
        mesh.AddPointData("libmesh:id", lm_ids(kept));

    // Cells: active elements, one block per meshio++ type in order of first use.
    std::vector<std::string> block_types;
    std::map<std::string, std::vector<std::size_t>> by_type;
    std::map<std::string, std::size_t> dropped_nodes, skipped;
    for (std::size_t e = 0; e < ne; ++e) {
        if (!active[e])
            continue;
        const LmType* t = f.mElements[e].mType;
        if (!t->mCell) {
            ++skipped[t->mName];
            continue;
        }
        auto [it, fresh] = by_type.emplace(t->mCell, std::vector<std::size_t>{});
        if (fresh)
            block_types.push_back(t->mCell);
        it->second.push_back(e);
        if (t->mKeep < t->mNodes)
            ++dropped_nodes[t->mName];
    }
    for (const auto& [name, count] : skipped)
        log::warn("libMesh: skipping {} {} element(s) (no meshio++ equivalent)", count, name);
    for (const auto& [name, count] : dropped_nodes) {
        log::warn("libMesh: {} {} element(s) keep only the nodes of the nearest meshio++ type",
                  count, name);
        detail::provenance_note("high-order-dropped", std::to_string(count) + " " + name +
                                                          " element(s) lost their extra nodes");
    }

    std::vector<std::int64_t> cell_of(ne, -1);  // file element -> global cell
    std::vector<int> cell_dim;
    std::vector<NDArray> sid_blocks, level_blocks, p_blocks;
    std::int64_t next_cell = 0;
    for (const std::string& type : block_types) {
        const std::vector<std::size_t>& members = by_type[type];
        const int k = cell_type_num_nodes(cell_type_from_name(type));
        const int dim = cell_type_dimension(cell_type_from_name(type));
        const detail::NodeOrder* order = detail::node_order("libmesh", type);
        NDArray conn(DType::Int64, {members.size(), static_cast<std::size_t>(k)});
        NDArray sid(DType::Int64, {members.size()});
        NDArray lvl(DType::Int64, {members.size()});
        NDArray pl(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const LmElement& el = f.mElements[members[r]];
            for (int j = 0; j < k; ++j) {
                const int src = order ? order->mToMeshio[static_cast<std::size_t>(j)] : j;
                c[r * static_cast<std::size_t>(k) + static_cast<std::size_t>(j)] =
                    node_index[static_cast<std::size_t>(el.mNodes[static_cast<std::size_t>(src)])];
            }
            sid.As<std::int64_t>()[r] = el.mSubdomain;
            lvl.As<std::int64_t>()[r] = el.mLevel;
            pl.As<std::int64_t>()[r] = el.mPLevel;
            cell_of[members[r]] = next_cell++;
            cell_dim.push_back(dim);
        }
        mesh.AddCellBlock(type, std::move(conn));
        sid_blocks.push_back(std::move(sid));
        level_blocks.push_back(std::move(lvl));
        p_blocks.push_back(std::move(pl));
    }
    if (block_types.empty())
        return mesh;
    mesh.AddCellData("libmesh:subdomain", std::move(sid_blocks));
    if (max_level > 0)
        mesh.AddCellData("libmesh:level", std::move(level_blocks));
    if (f.mInlinePLevel)
        mesh.AddCellData("libmesh:p_level", std::move(p_blocks));

    // Subdomains -> cell regions.
    std::map<std::int64_t, std::vector<std::int64_t>> by_sid;
    std::map<std::int64_t, int> sid_dim;
    for (std::size_t e = 0; e < ne; ++e) {
        if (cell_of[e] < 0)
            continue;
        const std::int64_t sid = f.mElements[e].mSubdomain;
        by_sid[sid].push_back(cell_of[e]);
        int& d = sid_dim.emplace(sid, -1).first->second;
        d = std::max(d, cell_dim[static_cast<std::size_t>(cell_of[e])]);
    }
    for (const auto& [sid, cells] : by_sid) {
        const auto name = f.mSubdomainNames.find(sid);
        mesh.AddRegion(Region(
            name != f.mSubdomainNames.end() ? name->second : "subdomain_" + std::to_string(sid),
            RegionKind::Cell, sid_dim[sid], sid, lm_ids(cells)));
    }

    // Side sets -> side regions, through the facet each side's corners name.
    if (!f.mSides.empty()) {
        const detail::FacetIndex facets(mesh);
        auto corners_of = [&](std::size_t e, std::size_t s) {
            const LmElement& el = f.mElements[e];
            std::vector<std::int64_t> out;
            for (int k : lm_sides(el.mType->mShape)[s])
                out.push_back(
                    node_index[static_cast<std::size_t>(el.mNodes[static_cast<std::size_t>(k)])]);
            return out;
        };
        auto point_of = [&](std::int64_t file_node) {
            const std::size_t n = static_cast<std::size_t>(file_node);
            return LmPoint{f.mCoords[3 * n], f.mCoords[3 * n + 1], f.mCoords[3 * n + 2]};
        };
        std::map<std::int64_t, std::vector<std::int64_t>> by_id;
        std::map<std::int64_t, int> id_dim;
        std::size_t unmatched = 0, point_sides = 0;
        auto add = [&](std::int64_t id, std::size_t e, std::size_t s) {
            const std::vector<std::int64_t> corners = corners_of(e, s);
            const detail::FacetHit* hit = facets.Find(corners.data(), corners.size());
            const std::int64_t cell = cell_of[e];
            const detail::FacetOwner* owner = nullptr;
            if (hit && hit->mFirst.mCell == cell)
                owner = &hit->mFirst;
            else if (hit && hit->mSecond.mCell == cell)
                owner = &hit->mSecond;
            if (!owner) {
                ++unmatched;
                return;
            }
            std::vector<std::int64_t>& v = by_id[id];
            v.push_back(owner->mCell);
            v.push_back(owner->mFacet);
            int& d = id_dim.emplace(id, -1).first->second;
            d = std::max(d, cell_dim[static_cast<std::size_t>(cell)] - 1);
        };
        for (const LmBoundary& b : f.mSides) {
            if (b.mElem < 0 || static_cast<std::size_t>(b.mElem) >= ne) {
                ++unmatched;
                continue;
            }
            const std::size_t e = static_cast<std::size_t>(b.mElem);
            const LmElement& el = f.mElements[e];
            const auto& sides = lm_sides(el.mType->mShape);
            if (el.mType->mShape == LmShape::Edge) {
                ++point_sides;
                continue;
            }
            if (b.mSide < 0 || static_cast<std::size_t>(b.mSide) >= sides.size() ||
                !el.mType->mCell) {
                ++unmatched;
                continue;
            }
            const std::size_t s = static_cast<std::size_t>(b.mSide);
            if (active[e]) {
                add(b.mId, e, s);
                continue;
            }
            // A refined element: its active descendants' sides that lie on it.
            std::vector<LmPoint> side_pts;
            for (int k : sides[s])
                side_pts.push_back(point_of(el.mNodes[static_cast<std::size_t>(k)]));
            std::vector<std::size_t> stack = children[e];
            while (!stack.empty()) {
                const std::size_t c = stack.back();
                stack.pop_back();
                if (!active[c]) {
                    stack.insert(stack.end(), children[c].begin(), children[c].end());
                    continue;
                }
                const LmElement& ch = f.mElements[c];
                if (!ch.mType->mCell)
                    continue;
                const auto& csides = lm_sides(ch.mType->mShape);
                for (std::size_t cs = 0; cs < csides.size(); ++cs) {
                    bool on = true;
                    for (int k : csides[cs])
                        if (!lm_on_side(side_pts,
                                        point_of(ch.mNodes[static_cast<std::size_t>(k)]))) {
                            on = false;
                            break;
                        }
                    if (on)
                        add(b.mId, c, cs);
                }
            }
        }
        if (unmatched)
            log::warn("libMesh: {} boundary side(s) match no cell facet and are skipped",
                      unmatched);
        if (point_sides)
            log::warn(
                "libMesh: {} boundary side(s) of line elements (end points) have no side "
                "region form and are skipped",
                point_sides);
        for (auto& [id, entries] : by_id) {
            // Refined parents can list one child side twice (two set entries).
            std::set<std::pair<std::int64_t, std::int64_t>> unique;
            for (std::size_t k = 0; k + 1 < entries.size(); k += 2)
                unique.emplace(entries[k], entries[k + 1]);
            std::vector<std::int64_t> flat;
            for (const auto& [c, s] : unique) {
                flat.push_back(c);
                flat.push_back(s);
            }
            const auto name = f.mSidesetNames.find(id);
            mesh.AddRegion(Region(
                name != f.mSidesetNames.end() ? name->second : "boundary_" + std::to_string(id),
                RegionKind::Side, id_dim[id], id, lm_ids(flat, 2)));
        }
    }

    // Node sets -> point regions.
    std::map<std::int64_t, std::vector<std::int64_t>> by_nodeset;
    for (const auto& [node, id] : f.mNodesets)
        if (node >= 0 && static_cast<std::size_t>(node) < nn &&
            node_index[static_cast<std::size_t>(node)] >= 0)
            by_nodeset[id].push_back(node_index[static_cast<std::size_t>(node)]);
    for (const auto& [id, pts_in] : by_nodeset) {
        const auto name = f.mNodesetNames.find(id);
        mesh.AddRegion(
            Region(name != f.mNodesetNames.end() ? name->second : "nodeset_" + std::to_string(id),
                   RegionKind::Point, -1, id, lm_ids(pts_in)));
    }
    if (f.mSkippedEdgeBcs)
        log::warn("libMesh: {} edge/shell-face boundary condition(s) skipped", f.mSkippedEdgeBcs);
    return mesh;
}

}  // namespace meshioplusplus
