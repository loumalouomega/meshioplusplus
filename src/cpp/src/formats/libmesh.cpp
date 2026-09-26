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
#include <cstring>
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
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

// External includes
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
#include <bzlib.h>
#endif

namespace meshioplusplus {

namespace {

constexpr const char* kLmWhat = "libMesh";

// --- compression: libMesh's `.gz` (gzip) and `.bz2` (bzip2) files ---------------

// gzip, laid out as Python's gzip.compress(data, mtime=0) lays it out (raw
// deflate at level 9, mtime 0, XFL 2, OS 255), so both engines write the same
// bytes when they share a zlib.
std::string lm_gzip(const std::string& rIn) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    z_stream z{};
    if (deflateInit2(&z, Z_BEST_COMPRESSION, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        throw WriteError("libMesh: zlib deflateInit2 failed");
    std::string out("\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff", 10);
    out.resize(10 + deflateBound(&z, static_cast<uLong>(rIn.size())));
    z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(rIn.data()));
    z.avail_in = static_cast<uInt>(rIn.size());
    z.next_out = reinterpret_cast<Bytef*>(&out[10]);
    z.avail_out = static_cast<uInt>(out.size() - 10);
    const int rc = deflate(&z, Z_FINISH);
    const std::size_t produced = 10 + z.total_out;
    deflateEnd(&z);
    if (rc != Z_STREAM_END)
        throw WriteError("libMesh: zlib deflate failed");
    out.resize(produced);
    const uLong crc =
        crc32(0L, reinterpret_cast<const Bytef*>(rIn.data()), static_cast<uInt>(rIn.size()));
    for (const std::uint32_t v :
         {static_cast<std::uint32_t>(crc), static_cast<std::uint32_t>(rIn.size())})
        for (int b = 0; b < 4; ++b)
            out += static_cast<char>((v >> (8 * b)) & 0xff);
    return out;
#else
    (void)rIn;
    throw WriteError("libMesh: this build has no zlib (MESHIOPLUSPLUS_WITH_ZLIB) to write .gz");
#endif
}

std::string lm_bzip2(const std::string& rIn) {
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
    // The bound bzip2 documents: 1% larger plus 600 bytes.
    unsigned int size = static_cast<unsigned int>(rIn.size() + rIn.size() / 100 + 601);
    std::string out(size, '\0');
    const int rc = BZ2_bzBuffToBuffCompress(&out[0], &size, const_cast<char*>(rIn.data()),
                                            static_cast<unsigned int>(rIn.size()), 9, 0, 0);
    if (rc != BZ_OK)
        throw WriteError("libMesh: bzip2 compression failed (" + std::to_string(rc) + ")");
    out.resize(size);
    return out;
#else
    (void)rIn;
    throw WriteError("libMesh: this build has no bzip2 (MESHIOPLUSPLUS_WITH_BZIP2) to write .bz2");
#endif
}

// Every bzip2 stream of `rIn`, one after the other (as `bzip2 -d` reads them).
std::string lm_bunzip2(const std::string& rIn, const std::string& rPath) {
#ifdef MESHIOPLUSPLUS_HAS_BZIP2
    std::string out;
    std::size_t pos = 0;
    char buf[1 << 16];
    while (pos < rIn.size()) {
        bz_stream z{};
        if (BZ2_bzDecompressInit(&z, 0, 0) != BZ_OK)
            throw ReadError("libMesh: bzip2 initialisation failed");
        z.next_in = const_cast<char*>(rIn.data() + pos);
        z.avail_in = static_cast<unsigned int>(rIn.size() - pos);
        int rc = BZ_OK;
        while (rc == BZ_OK) {
            z.next_out = buf;
            z.avail_out = sizeof(buf);
            rc = BZ2_bzDecompress(&z);
            out.append(buf, sizeof(buf) - z.avail_out);
            if (rc == BZ_OK && z.avail_in == 0 && z.avail_out != 0)
                rc = BZ_UNEXPECTED_EOF;
        }
        pos = rIn.size() - z.avail_in;
        BZ2_bzDecompressEnd(&z);
        if (rc != BZ_STREAM_END)
            throw ReadError("libMesh: " + rPath + " is a corrupt or truncated bzip2 stream");
    }
    return out;
#else
    (void)rIn;
    throw ReadError("libMesh: " + rPath +
                    " is bzip2-compressed and this build has no bzip2 "
                    "(MESHIOPLUSPLUS_WITH_BZIP2; the Python reader inflates it)");
#endif
}

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

// Corner nodes of each edge, libMesh's edge numbering (Hex8::edge_nodes_map ...).
// A quadratic element's mid-edge node for edge k is node `vertices + k`.
const std::vector<std::array<int, 2>>& lm_edges(LmShape Shape) {
    static const std::vector<std::array<int, 2>> none;
    static const std::vector<std::array<int, 2>> tri = {{0, 1}, {1, 2}, {2, 0}};
    static const std::vector<std::array<int, 2>> quad = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const std::vector<std::array<int, 2>> tet = {{0, 1}, {1, 2}, {0, 2},
                                                        {0, 3}, {1, 3}, {2, 3}};
    static const std::vector<std::array<int, 2>> hex = {{0, 1}, {1, 2}, {2, 3}, {0, 3},
                                                        {0, 4}, {1, 5}, {2, 6}, {3, 7},
                                                        {4, 5}, {5, 6}, {6, 7}, {4, 7}};
    static const std::vector<std::array<int, 2>> prism = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 4},
                                                          {2, 5}, {3, 4}, {4, 5}, {3, 5}};
    static const std::vector<std::array<int, 2>> pyramid = {{0, 1}, {1, 2}, {2, 3}, {0, 3},
                                                            {0, 4}, {1, 4}, {2, 4}, {3, 4}};
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

int lm_vertices(LmShape Shape) {
    switch (Shape) {
        case LmShape::Point:
            return 1;
        case LmShape::Edge:
            return 2;
        case LmShape::Tri:
            return 3;
        case LmShape::Quad:
        case LmShape::Tet:
            return 4;
        case LmShape::Pyramid:
            return 5;
        case LmShape::Prism:
            return 6;
        case LmShape::Hex:
            return 8;
        default:
            return 0;
    }
}

// Region-name suffixes for the sets that are not sides: an edge set's line
// cells, and a shell face (0 or 1) of 2-D elements.
constexpr const char* kLmEdgeSuffix = ":edge";
constexpr const char* kLmShellfaceSuffix = ":shellface";

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

    /** @brief The whole input's size: the bound on any count it declares. */
    std::size_t Size() const { return mText.size(); }

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
    std::vector<LmBoundary> mSides, mEdges, mShellfaces;
    std::vector<std::pair<std::int64_t, std::int64_t>> mNodesets;  // (node, id)
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
    // Each element takes at least a byte of the file and each node three
    // coordinates of at least a byte: counts beyond that are corruption, not
    // something to reserve memory for.
    if (n_elem < 0 || static_cast<std::uint64_t>(n_elem) > io.Size() || n_nodes < 0 ||
        static_cast<std::uint64_t>(n_nodes) > io.Size() / 3)
        io.Fail("element or node count larger than the file");
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
    // Side, edge and shell-face sets share one id space, and each repeats the
    // sideset name map.
    auto read_triples = [&](std::vector<LmBoundary>& rOut) {
        if (v092)
            f.mSidesetNames.merge(lm_name_map(io, hw));
        const std::int64_t n = io.Scalar(hw);
        for (std::int64_t k = 0; k < n; ++k) {
            LmBoundary b{};
            b.mElem = io.StreamInt(tw);
            b.mSide = io.StreamInt(tw);
            b.mId = io.StreamInt(tw);
            rOut.push_back(b);
        }
    };
    read_triples(f.mSides);
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
        read_triples(f.mEdges);
        read_triples(f.mShellfaces);
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

// --- writing --------------------------------------------------------------------------

// libMesh's `Xdr` in WRITE (ASCII) or ENCODE (XDR) mode, the reverse of
// `LmStream`: a scalar is `value\t # comment`, a vector its length line then
// `item\t ` per item, bulk data bare values; XDR is big-endian with 4-byte
// lengths and strings padded to 4 bytes.
class LmOut {
public:
    explicit LmOut(bool Xdr) : mXdr(Xdr) {}

    void String(const std::string& rValue, const char* pComment = "") {
        if (mXdr) {
            BinString(rValue);
            return;
        }
        mOut += rValue;
        Comment(pComment);
    }

    void Scalar(std::int64_t Value, const char* pComment, int Width = 8) {
        if (mXdr) {
            BinInt(Value, Width);
            return;
        }
        mOut += std::to_string(Value);
        Comment(pComment);
    }

    void IntVector(const std::vector<std::int64_t>& rValues, const char* pComment = "") {
        if (mXdr) {
            BinInt(static_cast<std::int64_t>(rValues.size()), 4);
            for (std::int64_t v : rValues)
                BinInt(v, 8);
            return;
        }
        Scalar(static_cast<std::int64_t>(rValues.size()), "# vector length");
        for (std::int64_t v : rValues)
            mOut += std::to_string(v) + "\t ";
        Comment(pComment);
    }

    void StringVector(const std::vector<std::string>& rValues, const char* pComment = "") {
        if (mXdr) {
            BinInt(static_cast<std::int64_t>(rValues.size()), 4);
            for (const std::string& v : rValues)
                BinString(v);
            return;
        }
        Scalar(static_cast<std::int64_t>(rValues.size()), "# vector length");
        for (const std::string& v : rValues)
            mOut += v + "\t ";
        Comment(pComment);
    }

    // `data_stream` of 8-byte integers, one line in ASCII.
    void Ints(const std::vector<std::int64_t>& rValues) {
        if (mXdr) {
            for (std::int64_t v : rValues)
                BinInt(v, 8);
            return;
        }
        for (std::size_t k = 0; k < rValues.size(); ++k) {
            if (k)
                mOut += ' ';
            mOut += std::to_string(rValues[k]);
        }
        mOut += '\n';
    }

    // `data_stream` of reals, three to a line in ASCII (`%.17e`, as libMesh's
    // `std::scientific` with `max_digits10`).
    void Reals(const double* pValues, std::size_t N) {
        char buf[40];
        for (std::size_t k = 0; k < N; ++k) {
            if (mXdr) {
                std::uint64_t bits = 0;
                std::memcpy(&bits, &pValues[k], 8);
                BinInt(static_cast<std::int64_t>(bits), 8);
                continue;
            }
            detail::snprintf_c(buf, sizeof(buf), "%.17e", pValues[k]);
            mOut += buf;
            mOut += (k % 3 == 2 || k + 1 == N) ? '\n' : ' ';
        }
    }

    const std::string& Data() const { return mOut; }

private:
    void Comment(const char* pComment) {
        if (*pComment) {
            mOut += "\t ";
            mOut += pComment;
        }
        mOut += '\n';
    }

    void BinInt(std::int64_t Value, int Width) {
        const std::uint64_t v = static_cast<std::uint64_t>(Value);
        for (int b = Width - 1; b >= 0; --b)
            mOut += static_cast<char>((v >> (8 * b)) & 0xff);
    }

    void BinString(const std::string& rValue) {
        BinInt(static_cast<std::int64_t>(rValue.size()), 4);
        mOut += rValue;
        mOut.append((4 - rValue.size() % 4) % 4, '\0');
    }

    bool mXdr;
    std::string mOut;
};

// meshio++ type -> libMesh ElemType, for the types libMesh stores as they are.
int lm_code(std::string_view Type) {
    static const std::pair<const char*, int> codes[] = {
        {"line", 0},          {"line3", 1},      {"line4", 2},       {"triangle", 3},
        {"triangle6", 4},     {"quad", 5},       {"quad8", 6},       {"quad9", 7},
        {"tetra", 8},         {"tetra10", 9},    {"hexahedron", 10}, {"hexahedron20", 11},
        {"hexahedron27", 12}, {"wedge", 13},     {"wedge15", 14},    {"wedge18", 15},
        {"pyramid", 16},      {"pyramid13", 17}, {"pyramid14", 18},  {"vertex", 27},
        {"triangle7", 33},
    };
    for (const auto& [name, code] : codes)
        if (Type == name)
            return code;
    return -1;
}

bool lm_ends_with(const std::string& rName, std::string_view Suffix) {
    return rName.size() > Suffix.size() &&
           rName.compare(rName.size() - Suffix.size(), Suffix.size(), Suffix) == 0;
}

// A shell-face region's face (0 or 1) and base name, else -1.
int lm_shellface(const std::string& rName, std::string* pBase) {
    for (int k = 0; k < 2; ++k)
        if (lm_ends_with(rName, std::string(kLmShellfaceSuffix) + std::to_string(k))) {
            if (pBase)
                *pBase =
                    rName.substr(0, rName.size() - std::string_view(kLmShellfaceSuffix).size() - 1);
            return k;
        }
    return -1;
}

}  // namespace

Mesh read_libmesh(const std::string& rPath) {
    // One bulk read; the text is owned, since a compressed file is replaced
    // by its inflated bytes.
    std::string text(detail::open_source(rPath, "libMesh: cannot open " + rPath).View());
    // libMesh writes `.xda.gz`/`.xdr.gz` through gzip and `.bz2` through bzip2.
    if (text.size() >= 2 && static_cast<unsigned char>(text[0]) == 0x1f &&
        static_cast<unsigned char>(text[1]) == 0x8b)
        text = detail::zlib_inflate(text, 15 + 16, nullptr, kLmWhat);
    else if (text.compare(0, 3, "BZh") == 0)
        text = lm_bunzip2(text, rPath);

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

    // Edge sets -> line cells (the element's edge corners, plus its mid-edge
    // node when it is quadratic), shared by every set naming that edge. They
    // are not libMesh elements: their subdomain (and level, p-level) is -1.
    std::map<std::int64_t, std::set<std::int64_t>> by_edge_set;
    if (!f.mEdges.empty()) {
        std::map<std::vector<std::int64_t>, std::pair<int, std::size_t>> seen;
        std::vector<std::int64_t> rows[2];  // line, line3
        std::size_t unmatched = 0;
        std::vector<std::pair<std::int64_t, std::pair<int, std::size_t>>> members;
        for (const LmBoundary& b : f.mEdges) {
            if (b.mElem < 0 || static_cast<std::size_t>(b.mElem) >= ne) {
                ++unmatched;
                continue;
            }
            const LmElement& el = f.mElements[static_cast<std::size_t>(b.mElem)];
            const auto& edges = lm_edges(el.mType->mShape);
            if (b.mSide < 0 || static_cast<std::size_t>(b.mSide) >= edges.size()) {
                ++unmatched;
                continue;
            }
            const std::size_t k = static_cast<std::size_t>(b.mSide);
            const std::size_t nv = static_cast<std::size_t>(lm_vertices(el.mType->mShape));
            std::vector<std::int64_t> row;
            for (int c : edges[k])
                row.push_back(
                    node_index[static_cast<std::size_t>(el.mNodes[static_cast<std::size_t>(c)])]);
            if (el.mNodes.size() >= nv + edges.size())
                row.push_back(node_index[static_cast<std::size_t>(el.mNodes[nv + k])]);
            std::vector<std::int64_t> key = row;
            std::sort(key.begin(), key.end());
            const int kind = row.size() == 3 ? 1 : 0;
            auto [it, fresh] =
                seen.emplace(key, std::make_pair(kind, rows[kind].size() / (2 + kind)));
            if (fresh)
                rows[kind].insert(rows[kind].end(), row.begin(), row.end());
            members.emplace_back(b.mId, it->second);
        }
        std::size_t first[2] = {0, 0};
        for (int kind = 0; kind < 2; ++kind) {
            const std::size_t width = static_cast<std::size_t>(2 + kind);
            const std::size_t n = rows[kind].size() / width;
            if (n == 0)
                continue;
            first[kind] = static_cast<std::size_t>(next_cell);
            NDArray conn(DType::Int64, {n, width});
            std::copy(rows[kind].begin(), rows[kind].end(), conn.As<std::int64_t>());
            mesh.AddCellBlock(kind ? "line3" : "line", std::move(conn));
            NDArray none(DType::Int64, {n});
            std::fill(none.As<std::int64_t>(), none.As<std::int64_t>() + n, -1);
            sid_blocks.push_back(none);
            level_blocks.push_back(none);
            p_blocks.push_back(std::move(none));
            next_cell += static_cast<std::int64_t>(n);
            cell_dim.insert(cell_dim.end(), n, 1);
        }
        for (const auto& [id, at] : members)
            by_edge_set[id].insert(static_cast<std::int64_t>(first[at.first] + at.second));
        if (unmatched)
            log::warn("libMesh: {} edge boundary condition(s) name no element edge and are skipped",
                      unmatched);
    }
    mesh.AddCellData("libmesh:subdomain", std::move(sid_blocks));
    if (max_level > 0) {
        mesh.AddCellData("libmesh:level", std::move(level_blocks));
        // The refinement tree, every element in file order, so the writer can
        // write it back: `libmesh:tree` (cell or -1, parent row or -1, libMesh
        // type, subdomain, p-level) and `libmesh:tree:nodes` (its nodes as
        // point indices, in libMesh's order, -1 past them).
        std::size_t width = 0;
        for (const LmElement& el : f.mElements)
            width = std::max(width, el.mNodes.size());
        NDArray tree(DType::Int64, {ne, 5});
        NDArray tree_nodes(DType::Int64, {ne, width});
        std::int64_t* t = tree.As<std::int64_t>();
        std::int64_t* tn = tree_nodes.As<std::int64_t>();
        std::fill(tn, tn + ne * width, -1);
        for (std::size_t e = 0; e < ne; ++e) {
            const LmElement& el = f.mElements[e];
            t[5 * e] = cell_of[e];
            t[5 * e + 1] = el.mParent;
            t[5 * e + 2] = static_cast<std::int64_t>(el.mCode);
            t[5 * e + 3] = el.mSubdomain;
            t[5 * e + 4] = el.mPLevel;
            for (std::size_t k = 0; k < el.mNodes.size(); ++k)
                tn[e * width + k] = node_index[static_cast<std::size_t>(el.mNodes[k])];
        }
        mesh.AddFieldData("libmesh:tree", std::move(tree));
        mesh.AddFieldData("libmesh:tree:nodes", std::move(tree_nodes));
    }
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

    // Edge sets -> cell regions on their line cells.
    auto set_name = [&](std::int64_t id) {
        const auto name = f.mSidesetNames.find(id);
        return name != f.mSidesetNames.end() ? name->second : "boundary_" + std::to_string(id);
    };
    for (const auto& [id, cells] : by_edge_set)
        mesh.AddRegion(Region(set_name(id) + kLmEdgeSuffix, RegionKind::Cell, 1, id,
                              lm_ids(std::vector<std::int64_t>(cells.begin(), cells.end()))));

    // Shell-face sets -> cell regions `<name>:shellface<k>` on the 2-D cells
    // (a refined element's face is carried to all its active descendants).
    if (!f.mShellfaces.empty()) {
        std::map<std::pair<std::int64_t, std::int64_t>, std::set<std::int64_t>> by_face;
        std::size_t unmatched = 0;
        for (const LmBoundary& b : f.mShellfaces) {
            if (b.mElem < 0 || static_cast<std::size_t>(b.mElem) >= ne || b.mSide < 0 ||
                b.mSide > 1) {
                ++unmatched;
                continue;
            }
            std::set<std::int64_t>& cells = by_face[{b.mId, b.mSide}];
            std::vector<std::size_t> stack{static_cast<std::size_t>(b.mElem)};
            while (!stack.empty()) {
                const std::size_t e = stack.back();
                stack.pop_back();
                if (!active[e]) {
                    stack.insert(stack.end(), children[e].begin(), children[e].end());
                    continue;
                }
                if (cell_of[e] >= 0 && cell_dim[static_cast<std::size_t>(cell_of[e])] == 2)
                    cells.insert(cell_of[e]);
                else
                    ++unmatched;
            }
        }
        if (unmatched)
            log::warn(
                "libMesh: {} shell-face boundary condition(s) name no 2-D element and are "
                "skipped",
                unmatched);
        for (const auto& [key, cells] : by_face)
            if (!cells.empty())
                mesh.AddRegion(
                    Region(set_name(key.first) + kLmShellfaceSuffix + std::to_string(key.second),
                           RegionKind::Cell, 2, key.first,
                           lm_ids(std::vector<std::int64_t>(cells.begin(), cells.end()))));
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
    return mesh;
}

void write_libmesh(const std::string& rPath, const Mesh& rMesh) {
    // No provenance slot in this format: drop the notes this write raises on
    // the way out rather than let them reach the next file written.
    const detail::ProvenanceSlotlessWrite slotless;
    std::string lower;
    for (char c : rPath)
        lower += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    // A trailing `.gz`/`.bz2` compresses an ASCII stream, as libMesh does: its
    // XDR files ignore the suffix (`Xdr` opens them with plain stdio), so a
    // `.xdr.gz`/`.xdr.bz2` is plain XDR, which libMesh then reads.
    bool gz = lm_ends_with(lower, ".gz");
    bool bz2 = lm_ends_with(lower, ".bz2");
    if (gz)
        lower.resize(lower.size() - 3);
    if (bz2)
        lower.resize(lower.size() - 4);
    const bool xdr = lm_ends_with(lower, ".xdr");
    gz = gz && !xdr;
    bz2 = bz2 && !xdr;

    // Cells: global (block-major) offsets and each cell's dimension.
    const std::size_t nb = rMesh.NumCellBlocks();
    std::vector<std::size_t> start{0};
    std::vector<int> block_dim(nb, -1);
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cb = rMesh.Cells(b);
        start.push_back(start.back() + cb.NumCells());
        block_dim[b] = cell_type_dimension(cell_type_from_name(std::string(cb.Type())));
    }
    const std::size_t ncells = start.back();
    auto block_of = [&](std::int64_t Cell) {
        return static_cast<std::size_t>(
            std::upper_bound(start.begin(), start.end(), static_cast<std::size_t>(Cell)) -
            start.begin() - 1);
    };
    auto dim_of = [&](std::int64_t Cell) { return block_dim[block_of(Cell)]; };

    // Region roles: edge sets and shell faces (by name suffix), side sets,
    // node sets, and the remaining cell regions as subdomains.
    std::vector<const Region*> edge_sets, shell_sets, side_sets, node_sets, subdomains;
    std::vector<bool> placeholder(ncells, false);  // an edge set's line cells
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const std::int64_t* e = reg.mEntries.As<std::int64_t>();
        const std::size_t n = reg.mEntries.Size();
        if (reg.mKind == RegionKind::Side) {
            side_sets.push_back(&reg);
            continue;
        }
        if (reg.mKind == RegionKind::Point) {
            node_sets.push_back(&reg);
            continue;
        }
        bool valid = true;
        for (std::size_t k = 0; k < n; ++k)
            valid = valid && e[k] >= 0 && static_cast<std::size_t>(e[k]) < ncells;
        auto all_dim = [&](int Dim) {
            for (std::size_t k = 0; k < n; ++k)
                if (dim_of(e[k]) != Dim)
                    return false;
            return true;
        };
        if (valid && n && lm_ends_with(reg.mName, kLmEdgeSuffix) && all_dim(1)) {
            edge_sets.push_back(&reg);
            for (std::size_t k = 0; k < n; ++k)
                placeholder[static_cast<std::size_t>(e[k])] = true;
        } else if (valid && n && lm_shellface(reg.mName, nullptr) >= 0 && all_dim(2)) {
            shell_sets.push_back(&reg);
        } else if (valid) {
            subdomains.push_back(&reg);
        }
    }

    // Elements: every cell libMesh has a type for, in cell order.
    struct LmOutElem {
        int mCode;
        std::int64_t mCell;
        std::vector<std::int64_t> mNodes;  // point indices, libMesh order
    };
    std::vector<LmOutElem> elems;
    std::vector<std::int64_t> elem_of(ncells, -1);
    std::map<std::string, std::size_t> dropped;
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const int code = cb.IsRagged() ? -1 : lm_code(type);
        if (code < 0) {
            if (cb.NumCells())
                dropped[type] += cb.NumCells();
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const detail::NodeOrder* order = detail::node_order("libmesh", type);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::size_t g = start[b] + r;
            if (placeholder[g])
                continue;
            LmOutElem el{code, static_cast<std::int64_t>(g), {}};
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                el.mNodes.push_back(detail::read_int(conn, r * k + src));
            }
            elem_of[g] = static_cast<std::int64_t>(elems.size());
            elems.push_back(std::move(el));
        }
    }
    for (const auto& [type, count] : dropped) {
        log::warn("libMesh: {} '{}' cell(s) have no libMesh element type and are dropped", count,
                  type);
        detail::provenance_note("cells-dropped", std::to_string(count) + " '" + type +
                                                     "' cell(s) have no libMesh element type");
    }

    // Node ids: `libmesh:id` when it is a valid numbering, else the point index.
    const std::size_t np = rMesh.NumPoints();
    std::vector<std::int64_t> node_id(np);
    for (std::size_t p = 0; p < np; ++p)
        node_id[p] = static_cast<std::int64_t>(p);
    std::int64_t max_node_id = static_cast<std::int64_t>(np);
    if (rMesh.HasPointData("libmesh:id")) {
        const NDArray& ids = rMesh.PointData("libmesh:id");
        std::set<std::int64_t> unique;
        std::vector<std::int64_t> v(np);
        bool ok = ids.Size() == np;
        for (std::size_t p = 0; ok && p < np; ++p) {
            v[p] = detail::read_int(ids, p);
            ok = v[p] >= 0 && unique.insert(v[p]).second;
        }
        if (ok) {
            node_id = std::move(v);
            max_node_id = np ? *unique.rbegin() + 1 : 0;
        }
    }

    // Subdomain ids: `libmesh:subdomain`, else the first cell region holding
    // the cell (its tag, or a fresh id), else 0.
    std::vector<std::int64_t> sid(ncells, 0);
    std::map<std::int64_t, std::string> subdomain_names;
    if (rMesh.HasCellData("libmesh:subdomain")) {
        for (std::size_t b = 0; b < nb; ++b)
            for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r)
                sid[start[b] + r] = detail::read_int(rMesh.CellData("libmesh:subdomain", b), r);
        for (const Region* pReg : subdomains)
            if (pReg->mTag >= 0 && pReg->mName != "subdomain_" + std::to_string(pReg->mTag))
                subdomain_names.emplace(pReg->mTag, pReg->mName);
    } else {
        std::int64_t next = 0;
        for (const Region* pReg : subdomains)
            next = std::max(next, pReg->mTag + 1);
        std::vector<bool> assigned(ncells, false);
        for (const Region* pReg : subdomains) {
            const std::int64_t id = pReg->mTag >= 0 ? pReg->mTag : next++;
            const std::int64_t* e = pReg->mEntries.As<std::int64_t>();
            for (std::size_t k = 0; k < pReg->mEntries.Size(); ++k) {
                const std::size_t c = static_cast<std::size_t>(e[k]);
                if (!assigned[c]) {
                    assigned[c] = true;
                    sid[c] = id;
                }
            }
            if (pReg->mName != "subdomain_" + std::to_string(id))
                subdomain_names.emplace(id, pReg->mName);
        }
    }
    for (const LmOutElem& el : elems)
        if (sid[static_cast<std::size_t>(el.mCell)] < 0 ||
            sid[static_cast<std::size_t>(el.mCell)] > 65534)
            throw WriteError("libMesh: subdomain id " +
                             std::to_string(sid[static_cast<std::size_t>(el.mCell)]) +
                             " is outside libMesh's 0..65534");
    const bool write_p = rMesh.HasCellData("libmesh:p_level");

    // Boundary ids: one id space for side, edge and shell-face sets.
    std::int64_t next_bid = 0;
    for (const auto* pList : {&side_sets, &edge_sets, &shell_sets})
        for (const Region* pReg : *pList)
            next_bid = std::max(next_bid, pReg->mTag + 1);
    std::map<const Region*, std::int64_t> bid;
    std::map<std::int64_t, std::string> sideset_names;
    auto boundary_id = [&](const Region* pReg, const std::string& rBase) {
        auto it = bid.find(pReg);
        if (it != bid.end())
            return it->second;
        const std::int64_t id = pReg->mTag >= 0 ? pReg->mTag : next_bid++;
        bid.emplace(pReg, id);
        if (rBase != "boundary_" + std::to_string(id))
            sideset_names.emplace(id, rBase);
        return id;
    };

    // Side sets: (element, libMesh side, id), matched by the facet's corners.
    std::set<std::array<std::int64_t, 3>> sides;
    std::size_t sides_lost = 0;
    for (const Region* pReg : side_sets) {
        const std::int64_t id = boundary_id(pReg, pReg->mName);
        const std::int64_t* e = pReg->mEntries.As<std::int64_t>();
        for (std::size_t k = 0; k + 1 < pReg->mEntries.Size(); k += 2) {
            CellType ftype{};
            std::vector<std::int64_t> fnodes;
            if (e[k] < 0 || static_cast<std::size_t>(e[k]) >= ncells ||
                elem_of[static_cast<std::size_t>(e[k])] < 0 ||
                !detail::facet_nodes(rMesh, e[k], e[k + 1], ftype, fnodes)) {
                ++sides_lost;
                continue;
            }
            const LmOutElem& el =
                elems[static_cast<std::size_t>(elem_of[static_cast<std::size_t>(e[k])])];
            const std::string& fname = cell_type_name(ftype);
            const std::size_t corners = cell_type_dimension(ftype) == 1   ? 2
                                        : fname.rfind("triangle", 0) == 0 ? 3
                                                                          : 4;
            fnodes.resize(std::min(corners, fnodes.size()));
            std::sort(fnodes.begin(), fnodes.end());
            const auto& table = lm_sides(lm_type(static_cast<std::uint64_t>(el.mCode))->mShape);
            bool found = false;
            for (std::size_t s = 0; s < table.size() && !found; ++s) {
                std::vector<std::int64_t> key;
                for (int c : table[s])
                    key.push_back(el.mNodes[static_cast<std::size_t>(c)]);
                std::sort(key.begin(), key.end());
                if (key == fnodes) {
                    sides.insert({elem_of[static_cast<std::size_t>(e[k])],
                                  static_cast<std::int64_t>(s), id});
                    found = true;
                }
            }
            if (!found)
                ++sides_lost;
        }
    }

    // Edge sets: each line cell on the first element holding that edge.
    std::set<std::array<std::int64_t, 3>> edges;
    std::size_t edges_lost = 0;
    // Edge cells no active element holds (a refined element's whole edge):
    // (corner, corner, id), matched against the tree's level-0 elements.
    std::vector<std::array<std::int64_t, 3>> edges_unmatched;
    if (!edge_sets.empty()) {
        std::map<std::pair<std::int64_t, std::int64_t>, std::pair<std::int64_t, std::int64_t>>
            edge_owner;
        for (std::size_t i = 0; i < elems.size(); ++i) {
            const auto& table =
                lm_edges(lm_type(static_cast<std::uint64_t>(elems[i].mCode))->mShape);
            for (std::size_t k = 0; k < table.size(); ++k) {
                std::int64_t a = elems[i].mNodes[static_cast<std::size_t>(table[k][0])];
                std::int64_t b = elems[i].mNodes[static_cast<std::size_t>(table[k][1])];
                edge_owner.emplace(
                    std::make_pair(std::min(a, b), std::max(a, b)),
                    std::make_pair(static_cast<std::int64_t>(i), static_cast<std::int64_t>(k)));
            }
        }
        for (const Region* pReg : edge_sets) {
            const std::string base =
                pReg->mName.substr(0, pReg->mName.size() - std::string_view(kLmEdgeSuffix).size());
            const std::int64_t id = boundary_id(pReg, base);
            const std::int64_t* e = pReg->mEntries.As<std::int64_t>();
            for (std::size_t k = 0; k < pReg->mEntries.Size(); ++k) {
                const std::size_t b = block_of(e[k]);
                const NDArray& conn = rMesh.Cells(b).Conn();
                const std::size_t w = rMesh.Cells(b).NodesPerCell();
                const std::size_t r = static_cast<std::size_t>(e[k]) - start[b];
                const std::int64_t n0 = detail::read_int(conn, r * w);
                const std::int64_t n1 = detail::read_int(conn, r * w + 1);
                const auto it = edge_owner.find({std::min(n0, n1), std::max(n0, n1)});
                if (it == edge_owner.end()) {
                    edges_unmatched.push_back({std::min(n0, n1), std::max(n0, n1), id});
                    continue;
                }
                edges.insert({it->second.first, it->second.second, id});
            }
        }
    }

    // Shell faces.
    std::set<std::array<std::int64_t, 3>> shellfaces;
    for (const Region* pReg : shell_sets) {
        std::string base;
        const int face = lm_shellface(pReg->mName, &base);
        const std::int64_t id = boundary_id(pReg, base);
        const std::int64_t* e = pReg->mEntries.As<std::int64_t>();
        for (std::size_t k = 0; k < pReg->mEntries.Size(); ++k)
            if (elem_of[static_cast<std::size_t>(e[k])] >= 0)
                shellfaces.insert({elem_of[static_cast<std::size_t>(e[k])], face, id});
    }

    // Node sets.
    std::int64_t next_nid = 0;
    for (const Region* pReg : node_sets)
        next_nid = std::max(next_nid, pReg->mTag + 1);
    std::set<std::array<std::int64_t, 2>> nodesets;
    std::map<std::int64_t, std::string> nodeset_names;
    for (const Region* pReg : node_sets) {
        const std::int64_t id = pReg->mTag >= 0 ? pReg->mTag : next_nid++;
        if (pReg->mName != "nodeset_" + std::to_string(id))
            nodeset_names.emplace(id, pReg->mName);
        const std::int64_t* e = pReg->mEntries.As<std::int64_t>();
        for (std::size_t k = 0; k < pReg->mEntries.Size(); ++k)
            if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < np)
                nodesets.insert({node_id[static_cast<std::size_t>(e[k])], id});
    }
    // The refinement tree (`libmesh:tree`, as the reader keeps it), when it
    // still describes these cells: each written cell is exactly one leaf of
    // the tree, with the same type and nodes; else the mesh is written flat.
    struct LmTreeRow {
        std::int64_t mCell, mParent, mCode, mSid, mP;
        int mLevel;
        std::vector<std::int64_t> mNodes;
    };
    std::vector<LmTreeRow> tree;
    std::vector<std::int64_t> row_of_elem(elems.size(), -1);
    bool use_tree = false;
    if (rMesh.HasFieldData("libmesh:tree") && rMesh.HasFieldData("libmesh:tree:nodes")) {
        const NDArray& t = rMesh.FieldData("libmesh:tree");
        const NDArray& tn = rMesh.FieldData("libmesh:tree:nodes");
        const std::size_t nt = t.Shape().size() == 2 && t.Shape()[1] == 5 ? t.Shape()[0] : 0;
        const std::size_t w = tn.Shape().size() == 2 && tn.Shape()[0] == nt ? tn.Shape()[1] : 0;
        bool ok = nt > 0 && w > 0;
        std::vector<bool> has_child(nt, false);
        int last_level = 0;
        for (std::size_t r = 0; ok && r < nt; ++r) {
            LmTreeRow row{detail::read_int(t, 5 * r),
                          detail::read_int(t, 5 * r + 1),
                          detail::read_int(t, 5 * r + 2),
                          detail::read_int(t, 5 * r + 3),
                          detail::read_int(t, 5 * r + 4),
                          0,
                          {}};
            const LmType* type =
                row.mCode >= 0 ? lm_type(static_cast<std::uint64_t>(row.mCode)) : nullptr;
            ok = type && type->mNodes > 0 && row.mParent >= -1 &&
                 row.mParent < static_cast<std::int64_t>(r) && row.mSid >= 0 && row.mSid <= 65534;
            for (std::size_t k = 0; ok && k < w; ++k) {
                const std::int64_t v = detail::read_int(tn, r * w + k);
                if (v < 0)
                    break;
                ok = static_cast<std::size_t>(v) < np;
                row.mNodes.push_back(v);
            }
            if (!ok || row.mNodes.size() != static_cast<std::size_t>(type->mNodes)) {
                ok = false;
                break;
            }
            if (row.mParent >= 0) {
                row.mLevel = tree[static_cast<std::size_t>(row.mParent)].mLevel + 1;
                has_child[static_cast<std::size_t>(row.mParent)] = true;
            }
            ok = row.mLevel >= last_level;
            last_level = row.mLevel;
            tree.push_back(std::move(row));
        }
        std::vector<bool> named(ncells, false);
        for (std::size_t r = 0; ok && r < nt; ++r) {
            const LmTreeRow& row = tree[r];
            if (has_child[r]) {
                ok = row.mCell == -1;
                continue;
            }
            ok = row.mCell >= 0 && static_cast<std::size_t>(row.mCell) < ncells &&
                 elem_of[static_cast<std::size_t>(row.mCell)] >= 0 &&
                 !named[static_cast<std::size_t>(row.mCell)];
            if (!ok)
                break;
            named[static_cast<std::size_t>(row.mCell)] = true;
            const std::size_t i =
                static_cast<std::size_t>(elem_of[static_cast<std::size_t>(row.mCell)]);
            ok = elems[i].mCode == row.mCode && elems[i].mNodes == row.mNodes;
            row_of_elem[i] = static_cast<std::int64_t>(r);
        }
        for (std::size_t i = 0; ok && i < elems.size(); ++i)
            ok = row_of_elem[i] >= 0;
        if (!ok) {
            log::warn(
                "libMesh: `libmesh:tree` no longer matches the cells; the mesh is written "
                "flat, as its active cells");
            detail::provenance_note("refinement-tree-dropped",
                                    "libmesh:tree does not match the cells");
            tree.clear();
        }
        use_tree = ok;
    }
    if (use_tree) {
        // libMesh keeps boundary ids on level-0 elements: a set's entries are
        // lifted to each level-0 side (edge, shell face) whose active pieces
        // are all in the set.
        const std::size_t nt = tree.size();
        std::vector<std::vector<std::size_t>> kids(nt);
        std::vector<bool> leaf(nt, true);
        for (std::size_t r = 0; r < nt; ++r)
            if (tree[r].mParent >= 0) {
                kids[static_cast<std::size_t>(tree[r].mParent)].push_back(r);
                leaf[static_cast<std::size_t>(tree[r].mParent)] = false;
            }
        const NDArray& pts_in = rMesh.Points();
        const std::size_t dim = rMesh.PointDim();
        auto point = [&](std::int64_t P) {
            LmPoint q{0.0, 0.0, 0.0};
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                q[d] = detail::read_double(pts_in, static_cast<std::size_t>(P) * dim + d);
            return q;
        };
        auto leaves = [&](std::size_t Root) {
            std::vector<std::size_t> out, stack{Root};
            while (!stack.empty()) {
                const std::size_t r = stack.back();
                stack.pop_back();
                if (leaf[r])
                    out.push_back(r);
                else
                    stack.insert(stack.end(), kids[r].rbegin(), kids[r].rend());
            }
            return out;
        };
        auto shape_of = [&](std::size_t R) {
            return lm_type(static_cast<std::uint64_t>(tree[R].mCode))->mShape;
        };
        // A shape's sides, or its edges (corner pairs), as corner lists.
        auto table_of = [&](std::size_t R, bool Edges) {
            if (!Edges)
                return lm_sides(shape_of(R));
            std::vector<std::vector<int>> out;
            for (const auto& e : lm_edges(shape_of(R)))
                out.push_back({e[0], e[1]});
            return out;
        };
        // (row, local index) pieces of the leaves under Root lying on the
        // polygon (or segment) `rOn`, from the `rTable` of each leaf's shape.
        using Pieces = std::vector<std::pair<std::int64_t, std::int64_t>>;
        auto pieces_on = [&](std::size_t Root, const std::vector<LmPoint>& rOn, bool Edges) {
            Pieces out;
            for (const std::size_t c : leaves(Root)) {
                const auto table = table_of(c, Edges);
                for (std::size_t k = 0; k < table.size(); ++k) {
                    bool on = true;
                    for (int v : table[k])
                        on = on &&
                             lm_on_side(rOn, point(tree[c].mNodes[static_cast<std::size_t>(v)]));
                    if (on)
                        out.emplace_back(static_cast<std::int64_t>(c),
                                         static_cast<std::int64_t>(k));
                }
            }
            return out;
        };
        std::size_t lifted_lost = 0;
        auto lift = [&](const std::set<std::array<std::int64_t, 3>>& rIn, bool Edges) {
            std::map<std::int64_t, std::set<std::pair<std::int64_t, std::int64_t>>> want;
            for (const auto& t : rIn)
                want[t[2]].emplace(row_of_elem[static_cast<std::size_t>(t[0])], t[1]);
            std::set<std::array<std::int64_t, 3>> out;
            for (const auto& [id, pieces] : want) {
                std::set<std::pair<std::int64_t, std::int64_t>> covered;
                for (std::size_t r0 = 0; r0 < nt && tree[r0].mParent < 0; ++r0) {
                    const auto table = table_of(r0, Edges);
                    for (std::size_t k = 0; k < table.size(); ++k) {
                        std::vector<LmPoint> on;
                        for (int v : table[k])
                            on.push_back(point(tree[r0].mNodes[static_cast<std::size_t>(v)]));
                        const Pieces d = pieces_on(r0, on, Edges);
                        bool all = !d.empty();
                        for (const auto& piece : d)
                            all = all && pieces.count(piece);
                        if (!all)
                            continue;
                        out.insert(
                            {static_cast<std::int64_t>(r0), static_cast<std::int64_t>(k), id});
                        covered.insert(d.begin(), d.end());
                    }
                }
                for (const auto& piece : pieces)
                    lifted_lost += covered.count(piece) ? 0 : 1;
            }
            return out;
        };
        sides = lift(sides, false);
        edges = lift(edges, true);
        // Whole edges of refined level-0 elements.
        std::map<std::pair<std::int64_t, std::int64_t>, std::pair<std::int64_t, std::int64_t>>
            root_edge;
        for (std::size_t r0 = 0; r0 < nt && tree[r0].mParent < 0; ++r0) {
            const auto table = table_of(r0, true);
            for (std::size_t k = 0; k < table.size(); ++k) {
                const std::int64_t a = tree[r0].mNodes[static_cast<std::size_t>(table[k][0])];
                const std::int64_t b = tree[r0].mNodes[static_cast<std::size_t>(table[k][1])];
                root_edge.emplace(
                    std::make_pair(std::min(a, b), std::max(a, b)),
                    std::make_pair(static_cast<std::int64_t>(r0), static_cast<std::int64_t>(k)));
            }
        }
        std::vector<std::array<std::int64_t, 3>> still;
        for (const auto& u : edges_unmatched) {
            const auto it = root_edge.find({u[0], u[1]});
            if (it == root_edge.end())
                still.push_back(u);
            else
                edges.insert({it->second.first, it->second.second, u[2]});
        }
        edges_unmatched = std::move(still);
        {
            std::map<std::pair<std::int64_t, std::int64_t>, std::set<std::int64_t>> want;
            for (const auto& t : shellfaces)
                want[{t[2], t[1]}].insert(row_of_elem[static_cast<std::size_t>(t[0])]);
            std::set<std::array<std::int64_t, 3>> out;
            for (const auto& [key, rows] : want) {
                std::set<std::int64_t> covered;
                for (std::size_t r0 = 0; r0 < nt && tree[r0].mParent < 0; ++r0) {
                    bool all = true;
                    const std::vector<std::size_t> under = leaves(r0);
                    for (const std::size_t c : under)
                        all = all && rows.count(static_cast<std::int64_t>(c));
                    if (!all)
                        continue;
                    out.insert({static_cast<std::int64_t>(r0), key.second, key.first});
                    covered.insert(under.begin(), under.end());
                }
                for (const std::int64_t r : rows)
                    lifted_lost += covered.count(r) ? 0 : 1;
            }
            shellfaces = std::move(out);
        }
        if (lifted_lost) {
            log::warn(
                "libMesh: {} boundary set entries cover only part of a level-0 element's "
                "side, edge or face and are dropped (libMesh keeps boundary ids on level-0 "
                "elements)",
                lifted_lost);
            detail::provenance_note(
                "regions-dropped",
                std::to_string(lifted_lost) + " boundary entries are not whole level-0 sides");
        }
    }
    edges_lost += edges_unmatched.size();
    if (sides_lost || edges_lost) {
        log::warn(
            "libMesh: {} side and {} edge set entries match no written element and are "
            "dropped",
            sides_lost, edges_lost);
        detail::provenance_note("regions-dropped",
                                std::to_string(sides_lost + edges_lost) +
                                    " side/edge set entries match no libMesh element");
    }
    const bool bcs = !sides.empty() || !edges.empty() || !shellfaces.empty() || !nodesets.empty();

    // The stream, as XdrIO::write lays it out (libMesh-1.8.0, 8-byte ids).
    LmOut io(xdr);
    io.String("libMesh-1.8.0");
    io.Scalar(static_cast<std::int64_t>(use_tree ? tree.size() : elems.size()),
              "# number of elements");
    io.Scalar(max_node_id, "# number of nodes");
    io.String(bcs ? "." : "n/a", "# boundary condition specification file");
    io.String(".", "# subdomain id specification file");
    io.String("n/a", "# processor id specification file");
    io.String(write_p ? "." : "n/a", "# p-level specification file");
    io.Scalar(8, "# type size");
    io.Scalar(0, "# uid size");
    io.Scalar(0, "# pid size");
    io.Scalar(8, "# sid size");
    io.Scalar(write_p ? 8 : 0, "# p-level size");
    io.Scalar(bcs ? 8 : 0, "# eid size");
    io.Scalar(bcs ? 8 : 0, "# side size");
    io.Scalar(bcs ? 8 : 0, "# bid size");
    io.Scalar(0, "# extra integer size");
    io.StringVector({}, "# node integer names");
    io.StringVector({}, "# elem integer names");
    io.IntVector({}, "# elemset codes");

    auto name_map = [&](const std::map<std::int64_t, std::string>& rNames, const char* pComment) {
        io.Scalar(static_cast<std::int64_t>(rNames.size()), pComment);
        if (rNames.empty())
            return;
        std::vector<std::int64_t> ids;
        std::vector<std::string> names;
        for (const auto& [id, name] : rNames) {
            ids.push_back(id);
            names.push_back(name);
        }
        io.IntVector(ids);
        io.StringVector(names);
    };
    name_map(subdomain_names, "# subdomain id to name map");

    if (use_tree) {
        // Level by level, each element after its parent; the leaves' subdomain
        // and p-level from the cells, their ancestors' from the tree.
        for (std::size_t r = 0; r < tree.size(); ++r) {
            const LmTreeRow& row = tree[r];
            if (r == 0 || row.mLevel != tree[r - 1].mLevel) {
                std::size_t n = 0;
                for (std::size_t q = r; q < tree.size() && tree[q].mLevel == row.mLevel; ++q)
                    ++n;
                const std::string legend = "# n_elem at level " + std::to_string(row.mLevel) +
                                           ", [ type " + (row.mLevel ? "parent " : "") + "sid " +
                                           (write_p ? "p_level " : "") + "(n0 ... nN-1) ]";
                io.Scalar(static_cast<std::int64_t>(n), legend.c_str());
            }
            std::vector<std::int64_t> rec{row.mCode};
            if (row.mLevel)
                rec.push_back(row.mParent);
            if (row.mCell >= 0) {
                const std::size_t c = static_cast<std::size_t>(row.mCell);
                rec.push_back(sid[c]);
                if (write_p) {
                    const std::size_t b = block_of(row.mCell);
                    rec.push_back(
                        detail::read_int(rMesh.CellData("libmesh:p_level", b), c - start[b]));
                }
            } else {
                rec.push_back(row.mSid);
                if (write_p)
                    rec.push_back(row.mP);
            }
            for (std::int64_t n : row.mNodes)
                rec.push_back(node_id[static_cast<std::size_t>(n)]);
            io.Ints(rec);
        }
    }
    if (!use_tree && !elems.empty()) {
        const std::string legend = std::string("# n_elem at level 0, [ type sid ") +
                                   (write_p ? "p_level " : "") + "(n0 ... nN-1) ]";
        io.Scalar(static_cast<std::int64_t>(elems.size()), legend.c_str());
    }
    for (const LmOutElem& el : use_tree ? std::vector<LmOutElem>{} : elems) {
        const std::size_t c = static_cast<std::size_t>(el.mCell);
        std::vector<std::int64_t> rec{el.mCode, sid[c]};
        if (write_p) {
            const std::size_t b = block_of(el.mCell);
            rec.push_back(detail::read_int(rMesh.CellData("libmesh:p_level", b), c - start[b]));
        }
        for (std::int64_t n : el.mNodes)
            rec.push_back(node_id[static_cast<std::size_t>(n)]);
        io.Ints(rec);
    }

    // Unused node ids: NaN in XDR, as libMesh writes them; 0 in ASCII, where
    // libMesh's `>>` cannot read back the `nan` it writes (it only loads nodes
    // elements use, so the value is never looked at).
    std::vector<double> coords(3 * static_cast<std::size_t>(max_node_id),
                               xdr ? std::numeric_limits<double>::quiet_NaN() : 0.0);
    const NDArray& pts = rMesh.Points();
    const std::size_t pd = rMesh.PointDim();
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t d = 0; d < 3; ++d)
            coords[3 * static_cast<std::size_t>(node_id[p]) + d] =
                d < pd ? detail::read_double(pts, p * pd + d) : 0.0;
    io.Reals(coords.data(), coords.size());
    io.Scalar(0, "# presence of unique ids", 4);

    auto triples = [&](const std::set<std::array<std::int64_t, 3>>& rSet, const char* pComment) {
        name_map(sideset_names, "# sideset id to name map");
        io.Scalar(static_cast<std::int64_t>(rSet.size()), pComment);
        for (const auto& t : rSet)
            io.Ints({t[0], t[1], t[2]});
    };
    triples(sides, "# number of side boundary conditions");
    name_map(nodeset_names, "# nodeset id to name map");
    io.Scalar(static_cast<std::int64_t>(nodesets.size()), "# number of nodesets");
    for (const auto& t : nodesets)
        io.Ints({t[0], t[1]});
    triples(edges, "# number of edge boundary conditions");
    triples(shellfaces, "# number of shellface boundary conditions");

    const std::string data = gz ? lm_gzip(io.Data()) : bz2 ? lm_bzip2(io.Data()) : io.Data();
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("libMesh: cannot open " + rPath + " for writing");
    out.write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!out)
        throw WriteError("libMesh: failed writing " + rPath);
}

}  // namespace meshioplusplus
