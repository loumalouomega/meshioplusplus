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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Geometry tables. MFEM's local edge and face lists (fem/geom.cpp and
// Mesh::GetElementToFaceTable) in MFEM's vertex order; they decide the order
// in which edges and faces are first met, hence the degree-of-freedom numbering.
// ---------------------------------------------------------------------------

struct MfGeom {
    const char* mLinear;
    const char* mQuadratic;  // nullptr: no order-2 support
    int mDim;
    std::size_t mNumVertices;
    std::vector<std::array<int, 2>> mEdges;
    std::vector<std::vector<int>> mFaces;
    bool mInterior;  // one order-2 interior dof when this is an element
};

const std::vector<MfGeom>& mf_geoms() {
    static const std::vector<MfGeom> g = {
        {"vertex", "vertex", 0, 1, {}, {}, false},
        {"line", "line3", 1, 2, {}, {}, true},
        {"triangle", "triangle6", 2, 3, {{0, 1}, {1, 2}, {2, 0}}, {}, false},
        {"quad", "quad9", 2, 4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, {}, true},
        {"tetra",
         "tetra10",
         3,
         4,
         {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}},
         {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}},
         false},
        {"hexahedron",
         "hexahedron27",
         3,
         8,
         {{0, 1},
          {1, 2},
          {3, 2},
          {0, 3},
          {4, 5},
          {5, 6},
          {7, 6},
          {4, 7},
          {0, 4},
          {1, 5},
          {2, 6},
          {3, 7}},
         {{3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}},
         true},
        {"wedge",
         "wedge18",
         3,
         6,
         {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}},
         {{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}},
         false},
        {"pyramid",
         nullptr,
         3,
         5,
         {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 4}, {2, 4}, {3, 4}},
         {{3, 2, 1, 0}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}},
         false},
    };
    return g;
}

// The MFEM vertex that meshio++ corner `K` of a `Geom` cell comes from: always
// `K`. MFEM's reference cells are meshio++'s, the prism included (its base
// triangle turns towards the top face, as gmsh's does); MFEM's own VTK export
// reverses its prisms (mesh/vtk.cpp's PrismMap) only because classic VTK winds
// the wedge the other way round.
std::size_t mf_vtk_corner(int /*Geom*/, std::size_t K) {
    return K;
}

// The non-corner nodes of each order-2 meshio++ type, by the corners (meshio++
// local indices) they sit between: 2 = an edge, 4 = a quad face, all corners =
// the cell centre.
const std::vector<std::vector<int>>& mf_slots(std::string_view Type) {
    static const std::map<std::string_view, std::vector<std::vector<int>>> slots = {
        {"vertex", {}},
        {"line3", {{0, 1}}},
        {"triangle6", {{0, 1}, {1, 2}, {2, 0}}},
        {"quad9", {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 1, 2, 3}}},
        {"tetra10", {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}}},
        {"wedge18",
         {{0, 1},
          {1, 2},
          {2, 0},
          {3, 4},
          {4, 5},
          {5, 3},
          {0, 3},
          {1, 4},
          {2, 5},
          {0, 1, 4, 3},
          {1, 2, 5, 4},
          {2, 0, 3, 5}}},
        {"hexahedron27",
         {{0, 1},
          {1, 2},
          {2, 3},
          {3, 0},
          {4, 5},
          {5, 6},
          {6, 7},
          {7, 4},
          {0, 4},
          {1, 5},
          {2, 6},
          {3, 7},
          {0, 4, 7, 3},
          {1, 2, 6, 5},
          {0, 1, 5, 4},
          {3, 2, 6, 7},
          {0, 1, 2, 3},
          {4, 5, 6, 7},
          {0, 1, 2, 3, 4, 5, 6, 7}}},
    };
    static const std::vector<std::vector<int>> none;
    const auto it = slots.find(Type);
    return it == slots.end() ? none : it->second;
}

int mf_geom_of_type(std::string_view Type) {
    // Serendipity and quadratic pyramids: written by completing them (or as
    // their corners); never read.
    if (Type == "quad8")
        return 3;
    if (Type == "hexahedron20")
        return 5;
    if (Type == "wedge15")
        return 6;
    if (Type == "pyramid13" || Type == "pyramid14")
        return 7;
    const auto& g = mf_geoms();
    for (std::size_t k = 0; k < g.size(); ++k)
        if (Type == g[k].mLinear || (g[k].mQuadratic && Type == g[k].mQuadratic))
            return static_cast<int>(k);
    return -1;
}

using MfKey = std::vector<std::int64_t>;  // sorted vertex ids

struct MfKeyHash {
    std::size_t operator()(const MfKey& rKey) const {
        std::size_t h = rKey.size();
        for (std::int64_t v : rKey)
            h ^= std::hash<std::int64_t>()(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

MfKey mf_key(std::vector<std::int64_t> v) {
    std::sort(v.begin(), v.end());
    return v;
}

struct MfElement {
    std::int64_t mAttribute;
    int mGeom;
    std::vector<std::int64_t> mVertices;  // MFEM order
    std::size_t mLine;
};

// MFEM's order-2 H1 numbering for a set of elements: vertices, edges and (3-D)
// faces in order of first appearance, then element interiors.
struct MfNumbering {
    std::size_t mNumVertices = 0;
    std::unordered_map<MfKey, std::size_t, MfKeyHash> mEdges;  // key -> dof
    std::unordered_map<MfKey, std::size_t, MfKeyHash> mFaces;  // quad faces only
    std::vector<std::int64_t> mInterior;                       // element -> dof, or -1
    std::vector<MfKey> mNodeKeys;                              // dof -> vertices it averages
    std::size_t Size() const { return mNodeKeys.size(); }
};

MfNumbering mf_number(const std::vector<MfElement>& rElements, int Dim, std::size_t NumVertices) {
    const auto& geoms = mf_geoms();
    MfNumbering n;
    n.mNumVertices = NumVertices;
    for (std::size_t v = 0; v < NumVertices; ++v)
        n.mNodeKeys.push_back({static_cast<std::int64_t>(v)});
    if (Dim >= 2) {
        for (const MfElement& el : rElements)
            for (const auto& e : geoms[static_cast<std::size_t>(el.mGeom)].mEdges) {
                MfKey key = mf_key({el.mVertices[static_cast<std::size_t>(e[0])],
                                    el.mVertices[static_cast<std::size_t>(e[1])]});
                if (n.mEdges.emplace(key, n.mNodeKeys.size()).second)
                    n.mNodeKeys.push_back(std::move(key));
            }
    }
    if (Dim == 3) {
        // Faces are numbered over every face; only the quad ones carry a dof.
        std::unordered_map<MfKey, bool, MfKeyHash> seen;
        std::vector<MfKey> quads;
        for (const MfElement& el : rElements)
            for (const auto& f : geoms[static_cast<std::size_t>(el.mGeom)].mFaces) {
                std::vector<std::int64_t> v;
                for (int k : f)
                    v.push_back(el.mVertices[static_cast<std::size_t>(k)]);
                MfKey key = mf_key(std::move(v));
                if (seen.emplace(key, true).second && key.size() == 4)
                    quads.push_back(std::move(key));
            }
        for (MfKey& key : quads) {
            n.mFaces.emplace(key, n.mNodeKeys.size());
            n.mNodeKeys.push_back(std::move(key));
        }
    }
    for (const MfElement& el : rElements) {
        const MfGeom& g = geoms[static_cast<std::size_t>(el.mGeom)];
        if (g.mInterior && g.mDim == Dim) {
            n.mInterior.push_back(static_cast<std::int64_t>(n.mNodeKeys.size()));
            n.mNodeKeys.push_back(mf_key(el.mVertices));
        } else {
            n.mInterior.push_back(-1);
        }
    }
    return n;
}

// The dof of one non-corner node of a cell: `rSlotVertices` are the global
// vertices it sits between; `Element` is the element index when the cell is an
// element (for its interior), else -1. Returns -1 when the entity is unknown.
std::int64_t mf_slot_dof(const MfNumbering& rN, int Dim, const std::vector<std::int64_t>& rSlot,
                         std::size_t NumCorners, std::int64_t Element) {
    if (rSlot.size() == NumCorners && Element >= 0 && (Dim != 3 || rSlot.size() != 4) &&
        (Dim != 2 || rSlot.size() != 2)) {
        return rN.mInterior[static_cast<std::size_t>(Element)];
    }
    const MfKey key = mf_key(rSlot);
    if (key.size() == 2 && Dim >= 2) {
        const auto it = rN.mEdges.find(key);
        return it == rN.mEdges.end() ? -1 : static_cast<std::int64_t>(it->second);
    }
    if (key.size() == 4 && Dim == 3) {
        const auto it = rN.mFaces.find(key);
        return it == rN.mFaces.end() ? -1 : static_cast<std::int64_t>(it->second);
    }
    if (Element >= 0)
        return rN.mInterior[static_cast<std::size_t>(Element)];
    return -1;
}

// ---------------------------------------------------------------------------
// Tokenizer: whitespace-separated tokens; a line whose first non-blank
// character is `#` is a comment; `"..."` (with backslash escapes) is one token.
// ---------------------------------------------------------------------------

struct MfToken {
    std::string mText;
    std::size_t mLine;
    bool mQuoted;
};

class MfLexer {
public:
    MfLexer(const std::string& rWhat, const std::string& rText) : mWhat(rWhat) {
        std::size_t pos = 0;
        std::size_t line_no = 0;
        bool header_seen = false;
        while (pos < rText.size()) {
            std::size_t eol = rText.find('\n', pos);
            if (eol == std::string::npos)
                eol = rText.size();
            std::string_view line(rText.data() + pos, eol - pos);
            pos = eol + 1;
            ++line_no;
            if (!line.empty() && line.back() == '\r')
                line.remove_suffix(1);
            const std::size_t first = line.find_first_not_of(" \t");
            if (first == std::string_view::npos || line[first] == '#')
                continue;
            if (!header_seen) {
                const std::size_t last = line.find_last_not_of(" \t");
                mHeader = std::string(line.substr(first, last - first + 1));
                mHeaderLine = line_no;
                header_seen = true;
                continue;
            }
            // A few header-like lines are kept whole.
            std::string_view rest = line.substr(first);
            if (rest.rfind("FiniteElementCollection:", 0) == 0 || rest.rfind("VDim:", 0) == 0 ||
                rest.rfind("Ordering:", 0) == 0) {
                const std::size_t colon = rest.find(':');
                mTokens.push_back({std::string(rest.substr(0, colon + 1)), line_no, false});
                std::string_view value = rest.substr(colon + 1);
                const std::size_t b = value.find_first_not_of(" \t");
                const std::size_t e = value.find_last_not_of(" \t");
                mTokens.push_back({b == std::string_view::npos
                                       ? std::string()
                                       : std::string(value.substr(b, e - b + 1)),
                                   line_no, true});
                continue;
            }
            std::size_t i = first;
            while (i < line.size()) {
                while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
                    ++i;
                if (i >= line.size())
                    break;
                if (line[i] == '"') {
                    std::string text;
                    ++i;
                    while (i < line.size() && line[i] != '"') {
                        if (line[i] == '\\' && i + 1 < line.size())
                            ++i;
                        text.push_back(line[i++]);
                    }
                    if (i >= line.size())
                        Fail("unterminated quoted name", line_no);
                    ++i;
                    mTokens.push_back({std::move(text), line_no, true});
                    continue;
                }
                const std::size_t start = i;
                while (i < line.size() && line[i] != ' ' && line[i] != '\t')
                    ++i;
                mTokens.push_back({std::string(line.substr(start, i - start)), line_no, false});
            }
        }
        mEndLine = line_no;
    }

    [[noreturn]] void Fail(const std::string& rWhy, std::size_t Line) const {
        throw ReadError(mWhat + ": " + rWhy + " (line " + std::to_string(Line) + ")");
    }

    bool AtEnd() const { return mPos >= mTokens.size(); }
    const MfToken& Peek() const { return mTokens[mPos]; }
    std::size_t Line() const { return AtEnd() ? mEndLine : mTokens[mPos].mLine; }

    const MfToken& Next(const char* pExpected) {
        if (AtEnd())
            Fail(std::string("the file ends where ") + pExpected + " was expected", mEndLine);
        return mTokens[mPos++];
    }

    std::int64_t Int(const char* pExpected) {
        const MfToken& t = Next(pExpected);
        std::int64_t value = 0;
        if (!ParseInt(t.mText, value))
            Fail(std::string("expected ") + pExpected + ", found '" + t.mText + "'", t.mLine);
        return value;
    }

    double Real(const char* pExpected) {
        const MfToken& t = Next(pExpected);
        double value = 0;
        if (!ParseReal(t.mText, value))
            Fail(std::string("expected ") + pExpected + ", found '" + t.mText + "'", t.mLine);
        return value;
    }

    // Every number up to the next word (or the end); for a grid function's values.
    std::vector<double> Reals() {
        std::vector<double> out;
        double value = 0;
        while (!AtEnd() && ParseReal(mTokens[mPos].mText, value)) {
            out.push_back(value);
            ++mPos;
        }
        return out;
    }

    static bool ParseInt(const std::string& rText, std::int64_t& rValue) {
        if (rText.empty())
            return false;
        std::size_t i = rText[0] == '-' || rText[0] == '+' ? 1 : 0;
        if (i >= rText.size())
            return false;
        std::int64_t v = 0;
        for (; i < rText.size(); ++i) {
            if (rText[i] < '0' || rText[i] > '9')
                return false;
            v = v * 10 + (rText[i] - '0');
        }
        rValue = rText[0] == '-' ? -v : v;
        return true;
    }

    static bool ParseReal(const std::string& rText, double& rValue) {
        if (rText.empty())
            return false;
        const char c = rText[0];
        if (!((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.'))
            return false;
        const char* end = nullptr;
        rValue = detail::parse_double(rText.c_str(), end);
        return end == rText.c_str() + rText.size();
    }

    const std::string& Header() const { return mHeader; }
    std::size_t HeaderLine() const { return mHeaderLine; }

private:
    std::string mWhat;
    std::vector<MfToken> mTokens;
    std::size_t mPos = 0;
    std::string mHeader;
    std::size_t mHeaderLine = 0;
    std::size_t mEndLine = 0;
};

std::string mf_read_text(const std::string& rPath, const char* pWhat) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError(std::string(pWhat) + ": cannot open " + rPath);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

// A finite element space, as a `FiniteElementSpace` header names it.
struct MfSpace {
    enum Kind { H1, H1Other, L2T1, L2, Other } mKind = Other;
    int mOrder = -1;
    std::string mCollection;
    int mVDim = 1;
    int mOrdering = 0;  // 0 byNODES, 1 byVDIM
};

// `H1_2D_P2`, `H1@GL_3D_P1`, `Linear`, `Quadratic`, `L2_T1_2D_P1`, `L2_3D_P0`...
MfSpace mf_classify(const std::string& rName) {
    MfSpace s;
    s.mCollection = rName;
    auto order_after_P = [&]() {
        const std::size_t p = rName.rfind("_P");
        std::int64_t v = -1;
        if (p != std::string::npos && MfLexer::ParseInt(rName.substr(p + 2), v))
            return static_cast<int>(v);
        return -1;
    };
    if (rName == "Linear") {
        s.mKind = MfSpace::H1;
        s.mOrder = 1;
    } else if (rName == "Quadratic") {
        s.mKind = MfSpace::H1;
        s.mOrder = 2;
    } else if (rName == "QuadraticPos" || rName == "Cubic") {
        s.mKind = MfSpace::H1Other;
        s.mOrder = rName == "Cubic" ? 3 : 2;
    } else if (rName.rfind("H1_", 0) == 0 || rName.rfind("H1@", 0) == 0) {
        s.mOrder = order_after_P();
        // Gauss-Lobatto (the default) and equispaced points coincide up to order 2.
        s.mKind = s.mOrder >= 1 && s.mOrder <= 2 ? MfSpace::H1 : MfSpace::H1Other;
    } else if (rName.rfind("H1Pos_", 0) == 0 || rName.rfind("H1Ser_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = s.mOrder == 1 ? MfSpace::H1 : MfSpace::H1Other;
    } else if (rName.rfind("L2_T1_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = MfSpace::L2T1;
    } else if (rName.rfind("L2_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = MfSpace::L2;
    }
    return s;
}

MfSpace mf_read_space(MfLexer& rLex) {
    const MfToken& fes = rLex.Next("FiniteElementSpace");
    if (fes.mText != "FiniteElementSpace")
        rLex.Fail("expected FiniteElementSpace, found '" + fes.mText +
                      "' (NURBS and variable-order spaces are not supported)",
                  fes.mLine);
    MfSpace s;
    bool have_fec = false;
    while (!rLex.AtEnd()) {
        const std::string& key = rLex.Peek().mText;
        if (key == "FiniteElementCollection:") {
            rLex.Next("a collection");
            const std::string name = rLex.Next("a collection name").mText;
            const int vdim = s.mVDim, ordering = s.mOrdering;
            s = mf_classify(name);
            s.mVDim = vdim;
            s.mOrdering = ordering;
            have_fec = true;
        } else if (key == "VDim:") {
            rLex.Next("VDim");
            std::int64_t v = 0;
            const MfToken& t = rLex.Next("a VDim");
            if (!MfLexer::ParseInt(t.mText, v) || v < 1)
                rLex.Fail("bad VDim '" + t.mText + "'", t.mLine);
            s.mVDim = static_cast<int>(v);
        } else if (key == "Ordering:") {
            rLex.Next("Ordering");
            std::int64_t v = 0;
            const MfToken& t = rLex.Next("an Ordering");
            if (!MfLexer::ParseInt(t.mText, v) || (v != 0 && v != 1))
                rLex.Fail("bad Ordering '" + t.mText + "'", t.mLine);
            s.mOrdering = static_cast<int>(v);
        } else {
            break;
        }
    }
    if (!have_fec)
        rLex.Fail("a FiniteElementSpace without a FiniteElementCollection", fes.mLine);
    return s;
}

// Value `Component` of dof `Dof` in a grid function vector.
double mf_value(const std::vector<double>& rValues, const MfSpace& rSpace, std::size_t NumDofs,
                std::size_t Dof, std::size_t Component) {
    return rSpace.mOrdering == 0
               ? rValues[Component * NumDofs + Dof]
               : rValues[Dof * static_cast<std::size_t>(rSpace.mVDim) + Component];
}

struct MfAttributeSet {
    std::string mName;
    std::vector<std::int64_t> mAttributes;
};

struct MfFile {
    int mDim = -1;
    std::vector<MfElement> mElements;
    std::vector<MfElement> mBoundary;
    std::vector<MfAttributeSet> mSets, mBdrSets;
    std::size_t mNumVertices = 0;
    int mSpaceDim = 0;
    std::vector<double> mCoords;  // from `vertices`
    bool mHasNodes = false;
    MfSpace mNodesSpace;
    std::vector<double> mNodes;
};

std::vector<MfElement> mf_read_elements(MfLexer& rLex, const char* pWhat) {
    const std::int64_t n = rLex.Int("an element count");
    if (n < 0)
        rLex.Fail(std::string("negative ") + pWhat + " count", rLex.Line());
    std::vector<MfElement> out;
    out.reserve(static_cast<std::size_t>(n));
    for (std::int64_t k = 0; k < n; ++k) {
        const std::size_t line = rLex.Line();
        MfElement el;
        el.mAttribute = rLex.Int("an attribute");
        const std::int64_t geom = rLex.Int("a geometry type");
        if (geom < 0 || geom > 7)
            rLex.Fail("unknown geometry type " + std::to_string(geom), line);
        el.mGeom = static_cast<int>(geom);
        el.mLine = line;
        const std::size_t nv = mf_geoms()[static_cast<std::size_t>(geom)].mNumVertices;
        for (std::size_t j = 0; j < nv; ++j) {
            const std::int64_t v = rLex.Int("a vertex index");
            if (v < 0)
                rLex.Fail("negative vertex index", line);
            el.mVertices.push_back(v);
        }
        out.push_back(std::move(el));
    }
    return out;
}

std::vector<MfAttributeSet> mf_read_sets(MfLexer& rLex) {
    const std::int64_t n = rLex.Int("an attribute set count");
    std::vector<MfAttributeSet> out;
    for (std::int64_t k = 0; k < n; ++k) {
        MfAttributeSet s;
        s.mName = rLex.Next("an attribute set name").mText;
        const std::int64_t size = rLex.Int("an attribute set size");
        for (std::int64_t j = 0; j < size; ++j)
            s.mAttributes.push_back(rLex.Int("an attribute"));
        std::sort(s.mAttributes.begin(), s.mAttributes.end());
        s.mAttributes.erase(std::unique(s.mAttributes.begin(), s.mAttributes.end()),
                            s.mAttributes.end());
        out.push_back(std::move(s));
    }
    return out;
}

MfFile mf_parse(const std::string& rPath) {
    const std::string what = "MFEM mesh";
    MfLexer lex(what, mf_read_text(rPath, "MFEM mesh"));
    const std::string& header = lex.Header();
    if (header.rfind("MFEM NC mesh", 0) == 0)
        lex.Fail("non-conforming meshes ('" + header + "') are not supported", lex.HeaderLine());
    if (header.rfind("MFEM NURBS", 0) == 0 || header.rfind("MFEM INLINE", 0) == 0)
        lex.Fail("'" + header + "' meshes are not supported", lex.HeaderLine());
    if (header != "MFEM mesh v1.0" && header != "MFEM mesh v1.1" && header != "MFEM mesh v1.2" &&
        header != "MFEM mesh v1.3")
        lex.Fail("not an MFEM mesh (the first line is '" + header + "')", lex.HeaderLine());
    MfFile f;
    bool saw_vertices = false;
    while (!lex.AtEnd()) {
        const MfToken& t = lex.Next("a section");
        if (t.mText == "dimension") {
            const std::int64_t d = lex.Int("a dimension");
            if (d < 1 || d > 3)
                lex.Fail("dimension " + std::to_string(d) + " (1, 2 or 3)", t.mLine);
            f.mDim = static_cast<int>(d);
        } else if (t.mText == "elements") {
            f.mElements = mf_read_elements(lex, "element");
        } else if (t.mText == "boundary") {
            f.mBoundary = mf_read_elements(lex, "boundary element");
        } else if (t.mText == "attribute_sets") {
            f.mSets = mf_read_sets(lex);
        } else if (t.mText == "bdr_attribute_sets") {
            f.mBdrSets = mf_read_sets(lex);
        } else if (t.mText == "vertices") {
            saw_vertices = true;
            const std::int64_t nv = lex.Int("a vertex count");
            if (nv < 0)
                lex.Fail("negative vertex count", t.mLine);
            f.mNumVertices = static_cast<std::size_t>(nv);
            if (!lex.AtEnd() && lex.Peek().mText == "nodes") {
                lex.Next("nodes");
                f.mHasNodes = true;
                f.mNodesSpace = mf_read_space(lex);
                f.mNodes = lex.Reals();
                f.mSpaceDim = f.mNodesSpace.mVDim;
            } else {
                const std::int64_t sd = lex.Int("a space dimension");
                if (sd < 1 || sd > 3)
                    lex.Fail("space dimension " + std::to_string(sd) + " (1, 2 or 3)", t.mLine);
                f.mSpaceDim = static_cast<int>(sd);
                f.mCoords.reserve(f.mNumVertices * static_cast<std::size_t>(sd));
                for (std::size_t k = 0; k < f.mNumVertices * static_cast<std::size_t>(sd); ++k)
                    f.mCoords.push_back(lex.Real("a coordinate"));
            }
        } else if (t.mText == "vertex_parents" || t.mText == "coarse_elements") {
            lex.Fail("non-conforming meshes (a '" + t.mText + "' section) are not supported",
                     t.mLine);
        } else if (t.mText == "mfem_serial_mesh_end") {
            log::warn(
                "MFEM mesh: '{}' is one rank of a parallel mesh; reading its local part and "
                "ignoring the communication groups",
                rPath);
            break;
        } else if (t.mText == "mfem_mesh_end") {
            break;
        } else {
            lex.Fail("unexpected '" + t.mText + "'", t.mLine);
        }
    }
    if (f.mDim < 0)
        throw ReadError("MFEM mesh: no dimension section in " + rPath);
    if (!saw_vertices)
        throw ReadError("MFEM mesh: no vertices section in " + rPath);
    for (const auto* list : {&f.mElements, &f.mBoundary})
        for (const MfElement& el : *list)
            for (std::int64_t v : el.mVertices)
                if (static_cast<std::size_t>(v) >= f.mNumVertices)
                    lex.Fail("vertex " + std::to_string(v) + " out of range (" +
                                 std::to_string(f.mNumVertices) + " vertices)",
                             el.mLine);
    return f;
}

struct MfGridData {
    std::string mName;
    MfSpace mSpace;
    std::vector<double> mValues;
};

MfGridData mf_parse_gf(const MfemGridFunction& rGf) {
    MfLexer lex("MFEM grid function", "\n" + mf_read_text(rGf.mPath, "MFEM grid function"));
    // The lexer took the first line as a header; a grid function has none, so
    // it was prefixed with an empty line and the header is "FiniteElementSpace".
    MfGridData g;
    g.mName = rGf.mName;
    if (lex.Header() != "FiniteElementSpace")
        throw ReadError("MFEM grid function: '" + rGf.mPath +
                        "' does not start with FiniteElementSpace (NURBS and variable-order "
                        "spaces are not supported)");
    // Re-parse the rest of the header from the token stream.
    MfSpace s;
    bool have_fec = false;
    while (!lex.AtEnd()) {
        const std::string& key = lex.Peek().mText;
        if (key == "FiniteElementCollection:") {
            lex.Next("a collection");
            const std::string name = lex.Next("a collection name").mText;
            const int vdim = s.mVDim, ordering = s.mOrdering;
            s = mf_classify(name);
            s.mVDim = vdim;
            s.mOrdering = ordering;
            have_fec = true;
        } else if (key == "VDim:") {
            lex.Next("VDim");
            std::int64_t v = 0;
            const MfToken& t = lex.Next("a VDim");
            if (!MfLexer::ParseInt(t.mText, v) || v < 1)
                lex.Fail("bad VDim '" + t.mText + "'", t.mLine);
            s.mVDim = static_cast<int>(v);
        } else if (key == "Ordering:") {
            lex.Next("Ordering");
            std::int64_t v = 0;
            const MfToken& t = lex.Next("an Ordering");
            if (!MfLexer::ParseInt(t.mText, v) || (v != 0 && v != 1))
                lex.Fail("bad Ordering '" + t.mText + "'", t.mLine);
            s.mOrdering = static_cast<int>(v);
        } else {
            break;
        }
    }
    if (!have_fec)
        throw ReadError("MFEM grid function: '" + rGf.mPath + "' names no FiniteElementCollection");
    g.mSpace = s;
    g.mValues = lex.Reals();
    if (!lex.AtEnd())
        lex.Fail("unexpected '" + lex.Peek().mText + "'", lex.Line());
    return g;
}

NDArray mf_entries(const std::vector<std::int64_t>& rIds) {
    NDArray a(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

Mesh read_mfem(const std::string& rPath) {
    return read_mfem(rPath, {});
}

Mesh read_mfem(const std::string& rPath, const std::vector<MfemGridFunction>& rGridFunctions) {
    MfFile f = mf_parse(rPath);
    const auto& geoms = mf_geoms();
    const int dim = f.mDim;
    const std::size_t nv = f.mNumVertices;

    // --- the nodes space ---------------------------------------------------------
    enum class Coords { Vertices, H1, Corners, Discontinuous } coords = Coords::Vertices;
    int mesh_order = 1;
    std::size_t node_dofs = 0;
    const int sdim = f.mSpaceDim;
    if (f.mHasNodes) {
        const MfSpace& s = f.mNodesSpace;
        if (f.mNodes.size() % static_cast<std::size_t>(s.mVDim) != 0)
            throw ReadError("MFEM mesh: the nodes of " + rPath + " hold " +
                            std::to_string(f.mNodes.size()) + " values, not a multiple of VDim " +
                            std::to_string(s.mVDim));
        node_dofs = f.mNodes.size() / static_cast<std::size_t>(s.mVDim);
        if (s.mKind == MfSpace::H1) {
            coords = Coords::H1;
            mesh_order = s.mOrder;
        } else if (s.mKind == MfSpace::H1Other) {
            coords = Coords::Corners;
            log::warn(
                "MFEM mesh: nodes in '{}' (order {}) are read at the vertices only; meshio++ "
                "keeps curved cells up to order 2",
                s.mCollection, s.mOrder);
        } else if (s.mKind == MfSpace::L2T1 && s.mOrder == 1) {
            coords = Coords::Discontinuous;
        } else {
            throw ReadError("MFEM mesh: nodes in '" + s.mCollection + "' are not supported");
        }
        if (coords != Coords::Discontinuous && node_dofs < nv)
            throw ReadError("MFEM mesh: the nodes of " + rPath + " have " +
                            std::to_string(node_dofs) + " dofs for " + std::to_string(nv) +
                            " vertices");
    }
    bool has_pyramid = false;
    for (const MfElement& el : f.mElements)
        has_pyramid = has_pyramid || el.mGeom == 7;

    // --- grid functions -------------------------------------------------------------
    std::vector<MfGridData> gfs;
    for (const MfemGridFunction& g : rGridFunctions) {
        MfGridData data = mf_parse_gf(g);
        const MfSpace& s = data.mSpace;
        const bool h1 = s.mKind == MfSpace::H1;
        const bool l2p0 = (s.mKind == MfSpace::L2 || s.mKind == MfSpace::L2T1) && s.mOrder == 0;
        if (!h1 && !l2p0) {
            log::warn("MFEM grid function '{}': the '{}' space is not supported; skipped", g.mPath,
                      s.mCollection);
            continue;
        }
        if (h1 && s.mOrder == 2 && (coords == Coords::Discontinuous || has_pyramid)) {
            log::warn(
                "MFEM grid function '{}': an order-2 field on this mesh is not supported; "
                "skipped",
                g.mPath);
            continue;
        }
        gfs.push_back(std::move(data));
    }
    int order = coords == Coords::Discontinuous ? 1 : mesh_order;
    for (const MfGridData& g : gfs)
        if (g.mSpace.mKind == MfSpace::H1)
            order = std::max(order, g.mSpace.mOrder);
    if (order == 2 && has_pyramid) {
        log::warn(
            "MFEM mesh: '{}' has pyramids, which meshio++ holds at order 1 only; reading "
            "the vertices",
            rPath);
        order = 1;
        if (coords == Coords::H1 && mesh_order == 2)
            coords = Coords::Corners;
    }
    if (order == 2)
        for (const MfElement& el : f.mElements)
            if (el.mGeom == 0 || geoms[static_cast<std::size_t>(el.mGeom)].mDim != dim)
                throw ReadError("MFEM mesh: element of geometry " + std::to_string(el.mGeom) +
                                " in a " + std::to_string(dim) + "-D mesh (line " +
                                std::to_string(el.mLine) + ")");

    // --- numbering and points ---------------------------------------------------
    MfNumbering numbering;
    if (order == 2) {
        numbering = mf_number(f.mElements, dim, nv);
    } else {
        numbering.mNumVertices = nv;
        for (std::size_t v = 0; v < nv; ++v)
            numbering.mNodeKeys.push_back({static_cast<std::int64_t>(v)});
    }
    if (coords == Coords::H1 && mesh_order == 2 && node_dofs != numbering.Size())
        throw ReadError("MFEM mesh: the order-2 nodes of " + rPath + " have " +
                        std::to_string(node_dofs) + " dofs; the mesh numbers " +
                        std::to_string(numbering.Size()));
    const std::size_t pdim = static_cast<std::size_t>(sdim);
    // Vertex coordinates first.
    std::vector<double> vxyz(nv * pdim, 0.0);
    if (coords == Coords::Vertices) {
        vxyz = f.mCoords;
    } else if (coords == Coords::H1 || coords == Coords::Corners) {
        for (std::size_t v = 0; v < nv; ++v)
            for (std::size_t c = 0; c < pdim; ++c)
                vxyz[v * pdim + c] = mf_value(f.mNodes, f.mNodesSpace, node_dofs, v, c);
    }

    Mesh mesh;
    // In the discontinuous case every element (and boundary element) gets its
    // own points; `point_vertex` records the MFEM vertex each point stands for.
    std::vector<std::int64_t> point_vertex;
    std::vector<MfKey> point_keys;  // continuous case: the dof keys
    NDArray points;
    // Per element: the point of each MFEM local vertex (discontinuous case).
    std::vector<std::vector<std::int64_t>> dg_points;
    if (coords == Coords::Discontinuous) {
        static const int quad_lex[4] = {0, 1, 3, 2};
        static const int hex_lex[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        std::vector<double> xyz;
        std::size_t offset = 0;
        for (const MfElement& el : f.mElements) {
            const MfGeom& g = geoms[static_cast<std::size_t>(el.mGeom)];
            if (el.mGeom != 1 && el.mGeom != 2 && el.mGeom != 3 && el.mGeom != 4 && el.mGeom != 5)
                throw ReadError("MFEM mesh: discontinuous nodes on geometry " +
                                std::to_string(el.mGeom) + " are not supported");
            std::vector<std::int64_t> pts;
            for (std::size_t j = 0; j < g.mNumVertices; ++j) {
                const std::size_t local =
                    el.mGeom == 3 ? static_cast<std::size_t>(quad_lex[j])
                                  : (el.mGeom == 5 ? static_cast<std::size_t>(hex_lex[j]) : j);
                if (offset + local >= node_dofs)
                    throw ReadError("MFEM mesh: the discontinuous nodes of " + rPath +
                                    " are too few for its elements");
                pts.push_back(static_cast<std::int64_t>(point_vertex.size()));
                point_vertex.push_back(el.mVertices[j]);
                for (std::size_t c = 0; c < pdim; ++c)
                    xyz.push_back(mf_value(f.mNodes, f.mNodesSpace, node_dofs, offset + local, c));
            }
            offset += g.mNumVertices;
            dg_points.push_back(std::move(pts));
        }
        if (offset != node_dofs)
            throw ReadError("MFEM mesh: the discontinuous nodes of " + rPath + " have " +
                            std::to_string(node_dofs) + " dofs; the elements need " +
                            std::to_string(offset));
        points = NDArray(DType::Float64, {point_vertex.size(), pdim});
        std::copy(xyz.begin(), xyz.end(), points.As<double>());
    } else {
        const std::size_t n = numbering.Size();
        points = NDArray(DType::Float64, {n, pdim});
        double* p = points.As<double>();
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t c = 0; c < pdim; ++c) {
                if (coords == Coords::H1 && mesh_order == 2) {
                    p[k * pdim + c] = mf_value(f.mNodes, f.mNodesSpace, node_dofs, k, c);
                } else {
                    double sum = 0;
                    for (std::int64_t v : numbering.mNodeKeys[k])
                        sum += vxyz[static_cast<std::size_t>(v) * pdim + c];
                    p[k * pdim + c] = sum / static_cast<double>(numbering.mNodeKeys[k].size());
                }
            }
        }
        point_keys = numbering.mNodeKeys;
    }
    const std::size_t npts =
        coords == Coords::Discontinuous ? point_vertex.size() : numbering.Size();
    mesh.AssignPoints(std::move(points));

    // --- cells --------------------------------------------------------------------
    // Elements then boundary elements, one block per type in order of first
    // appearance within each group.
    struct Cell {
        std::size_t mSource;  // index into elements or boundary
        bool mBoundary;
    };
    std::vector<std::pair<std::string, std::vector<Cell>>> blocks;
    auto type_of = [&](const MfElement& el) {
        const MfGeom& g = geoms[static_cast<std::size_t>(el.mGeom)];
        return std::string(order == 2 ? g.mQuadratic : g.mLinear);
    };
    auto add_group = [&](const std::vector<MfElement>& rList, bool Boundary) {
        std::map<std::string, std::size_t> index;
        for (std::size_t k = 0; k < rList.size(); ++k) {
            const std::string t = type_of(rList[k]);
            auto [it, fresh] = index.emplace(t, blocks.size());
            if (fresh)
                blocks.push_back({t, {}});
            blocks[it->second].second.push_back({k, Boundary});
        }
    };
    add_group(f.mElements, false);

    // Discontinuous case: a boundary element takes the points of an element
    // that holds all its vertices.
    std::unordered_map<std::int64_t, std::vector<std::size_t>> vertex_elements;
    if (coords == Coords::Discontinuous)
        for (std::size_t e = 0; e < f.mElements.size(); ++e)
            for (std::int64_t v : f.mElements[e].mVertices)
                vertex_elements[v].push_back(e);
    std::vector<MfElement> boundary;
    std::vector<std::vector<std::int64_t>> boundary_dg;
    std::size_t orphan_boundary = 0;
    for (const MfElement& b : f.mBoundary) {
        if (coords == Coords::Discontinuous) {
            std::vector<std::int64_t> pts;
            const auto it = vertex_elements.find(b.mVertices[0]);
            if (it != vertex_elements.end()) {
                for (std::size_t e : it->second) {
                    const MfElement& el = f.mElements[e];
                    pts.clear();
                    for (std::int64_t v : b.mVertices) {
                        const auto pos = std::find(el.mVertices.begin(), el.mVertices.end(), v);
                        if (pos == el.mVertices.end())
                            break;
                        pts.push_back(
                            dg_points[e][static_cast<std::size_t>(pos - el.mVertices.begin())]);
                    }
                    if (pts.size() == b.mVertices.size())
                        break;
                }
            }
            if (pts.size() != b.mVertices.size()) {
                ++orphan_boundary;
                continue;
            }
            boundary_dg.push_back(std::move(pts));
        }
        boundary.push_back(b);
    }
    if (orphan_boundary)
        log::warn("MFEM mesh: {} boundary element(s) lie on no element; dropped", orphan_boundary);
    add_group(boundary, true);

    std::size_t unresolved = 0;
    std::vector<NDArray> attr_blocks;
    std::vector<std::int64_t> cell_attr;
    std::vector<bool> cell_is_boundary;
    std::vector<std::size_t> element_cell(f.mElements.size());
    std::size_t global = 0;
    for (const auto& [type, members] : blocks) {
        const std::size_t k =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(type)));
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray attrs(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const Cell& cell = members[r];
            const MfElement& el =
                cell.mBoundary ? boundary[cell.mSource] : f.mElements[cell.mSource];
            const MfGeom& g = geoms[static_cast<std::size_t>(el.mGeom)];
            // Corners in meshio++ order.
            std::vector<std::int64_t> corner(g.mNumVertices);
            for (std::size_t j = 0; j < g.mNumVertices; ++j) {
                const std::size_t src = mf_vtk_corner(el.mGeom, j);
                corner[j] = el.mVertices[src];
                if (coords == Coords::Discontinuous) {
                    const auto& own =
                        cell.mBoundary ? boundary_dg[cell.mSource] : dg_points[cell.mSource];
                    c[r * k + j] = own[src];
                } else {
                    c[r * k + j] = el.mVertices[src];
                }
            }
            const auto& slots = mf_slots(type);
            const std::int64_t element =
                cell.mBoundary ? -1 : static_cast<std::int64_t>(cell.mSource);
            for (std::size_t s = 0; s < slots.size(); ++s) {
                std::vector<std::int64_t> sv;
                for (int q : slots[s])
                    sv.push_back(corner[static_cast<std::size_t>(q)]);
                std::int64_t dof = mf_slot_dof(numbering, dim, sv, g.mNumVertices, element);
                if (dof < 0) {
                    ++unresolved;
                    dof = sv[0];
                }
                c[r * k + g.mNumVertices + s] = dof;
            }
            attrs.As<std::int64_t>()[r] = el.mAttribute;
            cell_attr.push_back(el.mAttribute);
            cell_is_boundary.push_back(cell.mBoundary);
            if (!cell.mBoundary)
                element_cell[cell.mSource] = global;
            ++global;
        }
        mesh.AddCellBlock(type, std::move(conn));
        attr_blocks.push_back(std::move(attrs));
    }
    if (unresolved)
        log::warn("MFEM mesh: {} boundary node(s) lie on no element edge or face", unresolved);
    if (!attr_blocks.empty())
        mesh.AddCellData("mfem:attribute", std::move(attr_blocks));

    // --- grid functions as data -------------------------------------------------
    for (const MfGridData& g : gfs) {
        const MfSpace& s = g.mSpace;
        const std::size_t vdim = static_cast<std::size_t>(s.mVDim);
        if (g.mValues.size() % vdim != 0)
            throw ReadError("MFEM grid function '" + g.mName +
                            "': " + std::to_string(g.mValues.size()) +
                            " values, not a multiple of VDim " + std::to_string(vdim));
        const std::size_t ndofs = g.mValues.size() / vdim;
        if (s.mKind == MfSpace::H1) {
            const std::size_t expect = s.mOrder == 2 ? numbering.Size() : nv;
            if (ndofs != expect)
                throw ReadError("MFEM grid function '" + g.mName + "' has " +
                                std::to_string(ndofs) + " dofs; the mesh has " +
                                std::to_string(expect) + " at order " + std::to_string(s.mOrder));
            NDArray data =
                vdim == 1 ? NDArray(DType::Float64, {npts}) : NDArray(DType::Float64, {npts, vdim});
            double* d = data.As<double>();
            for (std::size_t p = 0; p < npts; ++p) {
                for (std::size_t c = 0; c < vdim; ++c) {
                    if (coords == Coords::Discontinuous) {
                        d[p * vdim + c] = mf_value(g.mValues, s, ndofs,
                                                   static_cast<std::size_t>(point_vertex[p]), c);
                    } else if (s.mOrder == 2) {
                        d[p * vdim + c] = mf_value(g.mValues, s, ndofs, p, c);
                    } else {
                        double sum = 0;
                        for (std::int64_t v : point_keys[p])
                            sum += mf_value(g.mValues, s, ndofs, static_cast<std::size_t>(v), c);
                        d[p * vdim + c] = sum / static_cast<double>(point_keys[p].size());
                    }
                }
            }
            mesh.AddPointData(g.mName, std::move(data));
        } else {
            if (ndofs != f.mElements.size())
                throw ReadError("MFEM grid function '" + g.mName + "' has " +
                                std::to_string(ndofs) + " dofs for " +
                                std::to_string(f.mElements.size()) + " elements");
            std::vector<double> per_cell(global * vdim, std::numeric_limits<double>::quiet_NaN());
            for (std::size_t e = 0; e < f.mElements.size(); ++e)
                for (std::size_t c = 0; c < vdim; ++c)
                    per_cell[element_cell[e] * vdim + c] = mf_value(g.mValues, s, ndofs, e, c);
            std::vector<NDArray> out;
            std::size_t start = 0;
            for (const auto& [type, members] : blocks) {
                NDArray a = vdim == 1 ? NDArray(DType::Float64, {members.size()})
                                      : NDArray(DType::Float64, {members.size(), vdim});
                std::copy(
                    per_cell.begin() + static_cast<std::ptrdiff_t>(start * vdim),
                    per_cell.begin() + static_cast<std::ptrdiff_t>((start + members.size()) * vdim),
                    a.As<double>());
                start += members.size();
                out.push_back(std::move(a));
            }
            mesh.AddCellData(g.mName, std::move(out));
        }
    }

    // --- regions ----------------------------------------------------------------
    std::map<std::pair<bool, std::int64_t>, std::vector<std::int64_t>> by_attr;
    for (std::size_t g = 0; g < cell_attr.size(); ++g)
        by_attr[{cell_is_boundary[g], cell_attr[g]}].push_back(static_cast<std::int64_t>(g));
    for (const auto& [key, ids] : by_attr) {
        const bool bdr = key.first;
        mesh.AddRegion(Region((bdr ? "boundary_" : "attribute_") + std::to_string(key.second),
                              RegionKind::Cell, bdr ? dim - 1 : dim, key.second, mf_entries(ids)));
    }
    auto add_sets = [&](const std::vector<MfAttributeSet>& rSets, bool Bdr) {
        for (const MfAttributeSet& s : rSets) {
            std::vector<std::int64_t> ids;
            for (std::int64_t a : s.mAttributes) {
                const auto it = by_attr.find({Bdr, a});
                if (it != by_attr.end())
                    ids.insert(ids.end(), it->second.begin(), it->second.end());
            }
            mesh.AddRegion(
                Region(s.mName, RegionKind::Cell, Bdr ? dim - 1 : dim, -1, mf_entries(ids)));
        }
    };
    add_sets(f.mSets, false);
    add_sets(f.mBdrSets, true);
    return mesh;
}

// ===========================================================================
// Writer
// ===========================================================================

namespace {

struct MfOutCell {
    int mGeom;
    std::int64_t mAttribute;
    std::size_t mCell;                 // global meshio++ cell (-1 for side facets)
    std::vector<std::int64_t> mNodes;  // meshio++ order, point indices
    std::string mType;                 // meshio++ type of mNodes
};

void mf_append_real(std::string& rOut, double Value) {
    char buf[40];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", Value);
    rOut += buf;
}

// The weights that express one output dof from mesh points.
using MfWeights = std::vector<std::pair<std::int64_t, double>>;

// The position (as point weights) of the node of a quadratic cell that sits
// between `rSlot` corners, where the cell has no such node itself.
MfWeights mf_centre(const MfOutCell& rCell, const std::vector<int>& rSlot) {
    MfWeights w;
    const std::string& t = rCell.mType;
    const std::size_t n = rSlot.size();
    // Serendipity completion: a quad8 face centre and a hex20 body centre.
    auto mid_of = [&](int a, int b) -> std::int64_t {
        const auto& slots =
            mf_slots(t == "quad8" ? "quad9"
                                  : (t == "hexahedron20" ? "hexahedron27"
                                                         : (t == "wedge15" ? "wedge18" : t)));
        const int nc =
            static_cast<int>(mf_geoms()[static_cast<std::size_t>(rCell.mGeom)].mNumVertices);
        for (std::size_t s = 0; s < slots.size(); ++s)
            if (slots[s].size() == 2 &&
                ((slots[s][0] == a && slots[s][1] == b) || (slots[s][0] == b && slots[s][1] == a)))
                return static_cast<std::int64_t>(nc) + static_cast<std::int64_t>(s);
        return -1;
    };
    const bool serendipity = t == "quad8" || t == "hexahedron20" || t == "wedge15";
    if (serendipity && n == 4) {
        for (std::size_t k = 0; k < 4; ++k) {
            w.emplace_back(rCell.mNodes[static_cast<std::size_t>(rSlot[k])], -0.25);
            const std::int64_t m = mid_of(rSlot[k], rSlot[(k + 1) % 4]);
            w.emplace_back(rCell.mNodes[static_cast<std::size_t>(m)], 0.5);
        }
        return w;
    }
    if (t == "hexahedron20" && n == 8) {
        for (std::size_t k = 0; k < 8; ++k)
            w.emplace_back(rCell.mNodes[k], -0.25);
        for (std::size_t k = 8; k < 20; ++k)
            w.emplace_back(rCell.mNodes[k], 0.25);
        return w;
    }
    for (int k : rSlot)
        w.emplace_back(rCell.mNodes[static_cast<std::size_t>(k)], 1.0 / static_cast<double>(n));
    return w;
}

std::string mf_sanitise(const std::string& rName) {
    std::string out;
    for (char c : rName) {
        const bool keep = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                          (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
        out.push_back(keep ? c : '_');
    }
    return out.empty() ? "data" : out;
}

std::string mf_quote(const std::string& rName) {
    std::string out = "\"";
    for (char c : rName) {
        if (c == '"' || c == '\\')
            out.push_back('\\');
        out.push_back(c == '\n' || c == '\r' ? ' ' : c);
    }
    return out + "\"";
}

bool mf_generated_name(const std::string& rName, const char* pPrefix, std::int64_t& rId) {
    const std::string prefix(pPrefix);
    if (rName.rfind(prefix, 0) != 0)
        return false;
    return MfLexer::ParseInt(rName.substr(prefix.size()), rId) && rId > 0 &&
           rName.substr(prefix.size()) == std::to_string(rId);
}

}  // namespace

void write_mfem(const std::string& rPath, const Mesh& rMesh) {
    write_mfem(rPath, rMesh, false);
}

void write_mfem(const std::string& rPath, const Mesh& rMesh, bool GridFunctions) {
    const std::size_t pdim = rMesh.PointDim();
    if (pdim < 1 || pdim > 3)
        throw WriteError("MFEM mesh writer: points of dimension " + std::to_string(pdim) +
                         " (1, 2 or 3)");
    const auto& geoms = mf_geoms();

    // --- classify cells --------------------------------------------------------
    int dim = -1;
    std::vector<int> block_geom(rMesh.NumCellBlocks(), -1);
    std::set<std::string> dropped;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const int g = cb.IsRagged() ? -1 : mf_geom_of_type(type);
        const bool known =
            g >= 0 && (type == geoms[static_cast<std::size_t>(g)].mLinear ||
                       type == geoms[static_cast<std::size_t>(g)].mQuadratic || type == "quad8" ||
                       type == "hexahedron20" || type == "wedge15" || type == "pyramid13" ||
                       type == "pyramid14");
        block_geom[b] = known ? g : -1;
        if (known && cb.NumCells())
            dim = std::max(dim, geoms[static_cast<std::size_t>(g)].mDim);
    }
    if (dim < 1)
        throw WriteError("MFEM mesh writer: no cells MFEM can hold as elements");
    if (pdim < static_cast<std::size_t>(dim))
        throw WriteError("MFEM mesh writer: " + std::to_string(pdim) + "-D points for " +
                         std::to_string(dim) + "-D cells");
    // Types MFEM names but not at the right dimension, and unknown types.
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const int g = block_geom[b];
        if (g < 0 || (geoms[static_cast<std::size_t>(g)].mDim != dim &&
                      geoms[static_cast<std::size_t>(g)].mDim != dim - 1)) {
            if (rMesh.Cells(b).NumCells())
                dropped.insert(std::string(rMesh.Cells(b).Type()));
            block_geom[b] = -1;
        }
    }

    // mfem:attribute, then cell regions, then 1.
    const bool has_attr = rMesh.HasCellData("mfem:attribute");
    std::vector<std::size_t> block_start(rMesh.NumCellBlocks() + 1, 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        block_start[b + 1] = block_start[b] + rMesh.Cells(b).NumCells();
    const std::size_t ncells = block_start.back();
    std::vector<int> cell_dim(ncells, -1);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        if (block_geom[b] >= 0)
            for (std::size_t g = block_start[b]; g < block_start[b + 1]; ++g)
                cell_dim[g] = geoms[static_cast<std::size_t>(block_geom[b])].mDim;
    std::vector<std::int64_t> attr(ncells, 0);
    if (has_attr) {
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& a = rMesh.CellData("mfem:attribute", b);
            for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r)
                attr[block_start[b] + r] = detail::read_int(a, r);
        }
    }
    // Regions: tags become attributes where mfem:attribute is absent; other
    // named regions become attribute sets.
    std::set<std::int64_t> used_attr[2];
    for (std::size_t g = 0; g < ncells; ++g)
        if (attr[g] > 0 && cell_dim[g] >= 0)
            used_attr[cell_dim[g] == dim ? 0 : 1].insert(attr[g]);
    struct NamedSet {
        std::string mName;
        std::vector<std::int64_t> mCells;
    };
    std::vector<NamedSet> named[2];
    std::int64_t next_attr[2] = {1, 1};
    auto fresh_attr = [&](int Which) {
        while (used_attr[Which].count(next_attr[Which]))
            ++next_attr[Which];
        used_attr[Which].insert(next_attr[Which]);
        return next_attr[Which];
    };
    struct SideSet {
        std::int64_t mAttribute;  // 0: a fresh one
        std::vector<std::pair<std::int64_t, std::int64_t>> mFacets;
        std::ptrdiff_t mSet;  // its entry in named[1], or -1
    };
    std::vector<SideSet> sides;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const std::int64_t* e = reg.Entries();
        if (reg.mKind == RegionKind::Point)
            continue;
        std::int64_t id = 0;
        const bool generated = mf_generated_name(reg.mName, "attribute_", id) ||
                               mf_generated_name(reg.mName, "boundary_", id);
        if (reg.mKind == RegionKind::Side) {
            SideSet side{generated ? id : reg.mTag, {}, -1};
            for (std::size_t j = 0; j < reg.NumEntries(); ++j)
                side.mFacets.emplace_back(e[2 * j], e[2 * j + 1]);
            if (side.mAttribute <= 0 || used_attr[1].count(side.mAttribute))
                side.mAttribute = 0;
            else
                used_attr[1].insert(side.mAttribute);
            if (!generated) {
                side.mSet = static_cast<std::ptrdiff_t>(named[1].size());
                named[1].push_back({reg.mName, {}});  // filled once the attribute is known
            }
            sides.push_back(std::move(side));
            continue;
        }
        std::vector<std::int64_t> cells(e, e + reg.NumEntries());
        if (!has_attr) {
            // The region's cells of each kind take its tag (or a fresh id) when
            // no earlier region gave them one.
            for (int which = 0; which < 2; ++which) {
                bool any = false;
                for (std::int64_t c : cells) {
                    const std::size_t g = static_cast<std::size_t>(c);
                    if (g < ncells && cell_dim[g] == (which == 0 ? dim : dim - 1) && attr[g] == 0)
                        any = true;
                }
                if (!any)
                    continue;
                std::int64_t a = generated ? id : reg.mTag;
                if (a <= 0 || used_attr[which].count(a))
                    a = fresh_attr(which);
                else
                    used_attr[which].insert(a);
                for (std::int64_t c : cells) {
                    const std::size_t g = static_cast<std::size_t>(c);
                    if (g < ncells && cell_dim[g] == (which == 0 ? dim : dim - 1) && attr[g] == 0)
                        attr[g] = a;
                }
            }
        }
        if (!generated) {
            std::vector<std::int64_t> el, bd;
            for (std::int64_t c : cells) {
                const std::size_t g = static_cast<std::size_t>(c);
                if (g < ncells && cell_dim[g] == dim)
                    el.push_back(c);
                else if (g < ncells && cell_dim[g] == dim - 1)
                    bd.push_back(c);
            }
            if (!el.empty())
                named[0].push_back({reg.mName, el});
            if (!bd.empty())
                named[1].push_back({reg.mName, bd});
        }
    }
    for (std::size_t g = 0; g < ncells; ++g)
        if (cell_dim[g] >= 0 && attr[g] <= 0) {
            attr[g] = 1;
            used_attr[cell_dim[g] == dim ? 0 : 1].insert(1);
        }

    // --- output cells -------------------------------------------------------------
    std::vector<MfOutCell> elements, boundary;
    bool quadratic = false, has_pyramid = false;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const int g = block_geom[b];
        if (g < 0)
            continue;
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        quadratic = quadratic || k > geoms[static_cast<std::size_t>(g)].mNumVertices;
        has_pyramid = has_pyramid || g == 7;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            MfOutCell c{g, attr[block_start[b] + r], block_start[b] + r, {}, type};
            for (std::size_t j = 0; j < k; ++j)
                c.mNodes.push_back(detail::read_int(conn, r * k + j));
            (geoms[static_cast<std::size_t>(g)].mDim == dim ? elements : boundary)
                .push_back(std::move(c));
        }
    }
    // Side regions: each facet one more boundary element.
    std::size_t bad_facets = 0;
    for (SideSet& side : sides) {
        if (side.mAttribute <= 0)
            side.mAttribute = fresh_attr(1);
        if (side.mSet >= 0)
            named[1][static_cast<std::size_t>(side.mSet)].mCells = {-side.mAttribute};
        for (const auto& [cell, facet] : side.mFacets) {
            CellType ft = CellType::Custom;
            std::vector<std::int64_t> nodes;
            if (!detail::facet_nodes(rMesh, cell, facet, ft, nodes)) {
                ++bad_facets;
                continue;
            }
            const std::string type = cell_type_name(ft);
            const int g = mf_geom_of_type(type);
            if (g < 0 || geoms[static_cast<std::size_t>(g)].mDim != dim - 1) {
                ++bad_facets;
                continue;
            }
            quadratic = quadratic || nodes.size() > geoms[static_cast<std::size_t>(g)].mNumVertices;
            boundary.push_back({g, side.mAttribute, static_cast<std::size_t>(-1), nodes, type});
        }
    }
    if (quadratic && has_pyramid) {
        log::warn("MFEM mesh writer: MFEM order-2 meshes hold no pyramids; writing corners only");
        detail::provenance_note("high-order-dropped",
                                "a mesh with pyramids is written with linear MFEM cells");
        quadratic = false;
    }
    for (const std::string& t : dropped) {
        log::warn("MFEM mesh writer: '{}' cells are neither elements nor boundary; dropped", t);
        detail::provenance_note("cells-dropped", "'" + t + "' cells have no MFEM equivalent here");
    }
    std::size_t point_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mKind == RegionKind::Point)
            ++point_regions;
    if (point_regions) {
        log::warn("MFEM mesh writer: MFEM has no node sets; {} point region(s) dropped",
                  point_regions);
        detail::provenance_note("regions-dropped", std::to_string(point_regions) +
                                                       " point region(s) have no MFEM equivalent");
    }
    if (bad_facets) {
        log::warn("MFEM mesh writer: {} side region entr(ies) name no facet and were dropped",
                  bad_facets);
        detail::provenance_note("regions-dropped",
                                std::to_string(bad_facets) + " side region entries name no facet");
    }

    // --- vertices ----------------------------------------------------------------
    // Linear: every point is a vertex. Quadratic: the corners only, renumbered.
    const std::size_t npts = rMesh.NumPoints();
    std::vector<std::int64_t> vertex_of(npts, -1);
    std::vector<std::int64_t> vertex_point;
    auto corners_of = [&](const MfOutCell& c) {
        return geoms[static_cast<std::size_t>(c.mGeom)].mNumVertices;
    };
    if (!quadratic) {
        bool all_corners = true;
        for (const auto* list : {&elements, &boundary})
            for (const MfOutCell& c : *list)
                if (c.mNodes.size() > corners_of(c))
                    all_corners = false;
        if (all_corners) {
            for (std::size_t p = 0; p < npts; ++p) {
                vertex_of[p] = static_cast<std::int64_t>(p);
                vertex_point.push_back(static_cast<std::int64_t>(p));
            }
        }
    }
    if (vertex_point.empty()) {
        std::vector<bool> is_corner(npts, false);
        for (const auto* list : {&elements, &boundary})
            for (const MfOutCell& c : *list)
                for (std::size_t j = 0; j < corners_of(c); ++j)
                    is_corner[static_cast<std::size_t>(c.mNodes[j])] = true;
        for (std::size_t p = 0; p < npts; ++p)
            if (is_corner[p]) {
                vertex_of[p] = static_cast<std::int64_t>(vertex_point.size());
                vertex_point.push_back(static_cast<std::int64_t>(p));
            }
    }
    // MFEM vertex lists.
    auto mfem_vertices = [&](const MfOutCell& c) {
        std::vector<std::int64_t> v(corners_of(c));
        for (std::size_t j = 0; j < v.size(); ++j)
            v[mf_vtk_corner(c.mGeom, j)] = vertex_of[static_cast<std::size_t>(c.mNodes[j])];
        return v;
    };

    // --- attribute sets -------------------------------------------------------------
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> sets[2];
    for (int which = 0; which < 2; ++which) {
        for (const NamedSet& s : named[which]) {
            std::set<std::int64_t> attrs;
            for (std::int64_t c : s.mCells) {
                if (c < 0)
                    attrs.insert(-c);
                else
                    attrs.insert(attr[static_cast<std::size_t>(c)]);
            }
            if (attrs.empty())
                continue;
            sets[which].emplace_back(s.mName,
                                     std::vector<std::int64_t>(attrs.begin(), attrs.end()));
        }
    }
    const bool v13 = !sets[0].empty() || !sets[1].empty();

    // --- numbering and dof weights (quadratic) -------------------------------
    std::vector<MfElement> numbered;
    for (const MfOutCell& c : elements)
        numbered.push_back({c.mAttribute, c.mGeom, mfem_vertices(c), 0});
    const std::size_t nv = vertex_point.size();
    MfNumbering numbering;
    std::vector<MfWeights> dof_weights;
    if (quadratic) {
        numbering = mf_number(numbered, dim, nv);
        dof_weights.resize(numbering.Size());
        std::vector<bool> exact(numbering.Size(), false);
        for (std::size_t v = 0; v < nv; ++v) {
            dof_weights[v] = {{vertex_point[v], 1.0}};
            exact[v] = true;
        }
        auto visit = [&](const MfOutCell& c, std::int64_t Element) {
            // The quadratic type whose slots this cell's nodes follow.
            std::string full = geoms[static_cast<std::size_t>(c.mGeom)].mQuadratic;
            const std::size_t nc = corners_of(c);
            const auto& slots = mf_slots(full);
            std::vector<std::int64_t> corner(nc);
            for (std::size_t j = 0; j < nc; ++j)
                corner[j] = vertex_of[static_cast<std::size_t>(c.mNodes[j])];
            for (std::size_t s = 0; s < slots.size(); ++s) {
                std::vector<std::int64_t> sv;
                for (int q : slots[s])
                    sv.push_back(corner[static_cast<std::size_t>(q)]);
                const std::int64_t dof = mf_slot_dof(numbering, dim, sv, nc, Element);
                if (dof < 0 || exact[static_cast<std::size_t>(dof)])
                    continue;
                const std::size_t slot_node = nc + s;
                if (slot_node < c.mNodes.size()) {
                    dof_weights[static_cast<std::size_t>(dof)] = {{c.mNodes[slot_node], 1.0}};
                    exact[static_cast<std::size_t>(dof)] = true;
                } else if (dof_weights[static_cast<std::size_t>(dof)].empty() ||
                           c.mNodes.size() > nc) {
                    dof_weights[static_cast<std::size_t>(dof)] = mf_centre(c, slots[s]);
                    if (c.mNodes.size() > nc)
                        exact[static_cast<std::size_t>(dof)] = true;
                }
            }
        };
        for (std::size_t e = 0; e < elements.size(); ++e)
            visit(elements[e], static_cast<std::int64_t>(e));
        for (const MfOutCell& c : boundary)
            visit(c, -1);
        // Anything no cell placed: the average of its vertices.
        for (std::size_t d = 0; d < numbering.Size(); ++d)
            if (dof_weights[d].empty())
                for (std::int64_t v : numbering.mNodeKeys[d])
                    dof_weights[d].emplace_back(
                        vertex_point[static_cast<std::size_t>(v)],
                        1.0 / static_cast<double>(numbering.mNodeKeys[d].size()));
    }

    // --- data --------------------------------------------------------------------
    const std::size_t data_arrays =
        rMesh.NumPointData() + rMesh.NumCellData() - (has_attr ? 1 : 0) + rMesh.NumFieldData();
    if (!GridFunctions && data_arrays) {
        log::warn("MFEM mesh writer: data arrays are dropped (write grid functions to keep them)");
        detail::provenance_note("data-dropped", "an MFEM mesh holds no data arrays");
    } else if (GridFunctions && rMesh.NumFieldData()) {
        log::warn("MFEM mesh writer: field data has no grid function; dropped");
        detail::provenance_note("data-dropped", "field data has no MFEM grid function");
    }

    // --- mesh file ---------------------------------------------------------------
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    std::string out = v13 ? "MFEM mesh v1.3\n" : "MFEM mesh v1.0\n";
    out += detail::provenance_render_lines(detail::SlotTier::Block, "# ");
    out += "\ndimension\n" + std::to_string(dim) + "\n\nelements\n" +
           std::to_string(elements.size()) + "\n";
    auto append_cells = [&](const std::vector<MfOutCell>& rList) {
        for (const MfOutCell& c : rList) {
            out += std::to_string(c.mAttribute) + " " + std::to_string(c.mGeom);
            for (std::int64_t v : mfem_vertices(c))
                out += " " + std::to_string(v);
            out += '\n';
        }
    };
    auto append_sets =
        [&](const std::vector<std::pair<std::string, std::vector<std::int64_t>>>& rSets) {
            out += std::to_string(rSets.size()) + "\n";
            for (const auto& [name, attrs] : rSets) {
                out += mf_quote(name) + " " + std::to_string(attrs.size());
                for (std::int64_t a : attrs)
                    out += " " + std::to_string(a);
                out += '\n';
            }
        };
    append_cells(elements);
    if (v13) {
        out += "\nattribute_sets\n";
        append_sets(sets[0]);
    }
    out += "\nboundary\n" + std::to_string(boundary.size()) + "\n";
    append_cells(boundary);
    if (v13) {
        out += "\nbdr_attribute_sets\n";
        append_sets(sets[1]);
    }
    out += "\nvertices\n" + std::to_string(nv) + "\n";
    const NDArray& pts = rMesh.Points();
    auto evaluate = [&](const MfWeights& rW, const NDArray& rData, std::size_t Cols,
                        std::size_t C) {
        double sum = 0;
        for (const auto& [p, w] : rW)
            sum += w * detail::read_double(rData, static_cast<std::size_t>(p) * Cols + C);
        return sum;
    };
    if (!quadratic) {
        out += std::to_string(pdim) + "\n";
        for (std::size_t v = 0; v < nv; ++v) {
            for (std::size_t c = 0; c < pdim; ++c) {
                if (c)
                    out += ' ';
                mf_append_real(out, detail::read_double(
                                        pts, static_cast<std::size_t>(vertex_point[v]) * pdim + c));
            }
            out += '\n';
        }
    } else {
        out += "\nnodes\nFiniteElementSpace\nFiniteElementCollection: H1_" + std::to_string(dim) +
               "D_P2\nVDim: " + std::to_string(pdim) + "\nOrdering: 1\n\n";
        for (std::size_t d = 0; d < numbering.Size(); ++d) {
            for (std::size_t c = 0; c < pdim; ++c) {
                if (c)
                    out += ' ';
                mf_append_real(out, evaluate(dof_weights[d], pts, pdim, c));
            }
            out += '\n';
        }
    }
    if (v13)
        out += "\nmfem_mesh_end\n";
    f << out;
    if (!f)
        throw WriteError("MFEM mesh writer: failed writing " + rPath);
    if (!GridFunctions)
        return;

    // --- grid functions ------------------------------------------------------------
    namespace fs = std::filesystem;
    const fs::path mesh_path(rPath);
    const std::string stem = (mesh_path.parent_path() / mesh_path.stem()).string();
    const int order = quadratic ? 2 : 1;
    const std::size_t ndofs = quadratic ? numbering.Size() : nv;
    auto open_gf = [&](const std::string& rName) {
        const std::string path = stem + "." + mf_sanitise(rName) + ".gf";
        auto g = detail::make_classic_ofstream(path, std::ios::binary);
        if (!g)
            throw WriteError("Could not open file for writing: " + path);
        return g;
    };
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const std::size_t cols = a.Shape().size() > 1 ? detail::cols(a) : 1;
        std::string text = "FiniteElementSpace\nFiniteElementCollection: H1_" +
                           std::to_string(dim) + "D_P" + std::to_string(order) +
                           "\nVDim: " + std::to_string(cols) + "\nOrdering: 1\n\n";
        for (std::size_t d = 0; d < ndofs; ++d) {
            for (std::size_t c = 0; c < cols; ++c) {
                if (c)
                    text += ' ';
                if (quadratic)
                    mf_append_real(text, evaluate(dof_weights[d], a, cols, c));
                else
                    mf_append_real(text,
                                   detail::read_double(
                                       a, static_cast<std::size_t>(vertex_point[d]) * cols + c));
            }
            text += '\n';
        }
        auto g = open_gf(name);
        g << text;
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name == "mfem:attribute")
            continue;
        std::size_t cols = 1;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& a = rMesh.CellData(name, b);
            if (a.Shape().size() > 1)
                cols = detail::cols(a);
        }
        std::string text = "FiniteElementSpace\nFiniteElementCollection: L2_" +
                           std::to_string(dim) + "D_P0\nVDim: " + std::to_string(cols) +
                           "\nOrdering: 1\n\n";
        for (const MfOutCell& c : elements) {
            std::size_t b = 0;
            while (c.mCell >= block_start[b + 1])
                ++b;
            const NDArray& a = rMesh.CellData(name, b);
            for (std::size_t k = 0; k < cols; ++k) {
                if (k)
                    text += ' ';
                mf_append_real(text, detail::read_double(a, (c.mCell - block_start[b]) * cols + k));
            }
            text += '\n';
        }
        auto g = open_gf(name);
        g << text;
    }
}

}  // namespace meshioplusplus
