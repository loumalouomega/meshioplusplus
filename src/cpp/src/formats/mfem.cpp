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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
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
#include "lagrange_common.hpp"

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
    // VTK Lagrange cells of any order: written as order-p H1 nodes.
    static const std::pair<std::string_view, int> kLagrange[] = {
        {"VTK_LAGRANGE_CURVE", 1},         {"VTK_LAGRANGE_TRIANGLE", 2},
        {"VTK_LAGRANGE_QUADRILATERAL", 3}, {"VTK_LAGRANGE_TETRAHEDRON", 4},
        {"VTK_LAGRANGE_HEXAHEDRON", 5},    {"VTK_LAGRANGE_WEDGE", 6}};
    for (const auto& [name, geom] : kLagrange)
        if (Type == name)
            return geom;
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
    enum Kind { H1, H1Other, L2T1, L2, Nurbs, Other } mKind = Other;
    // An H1 space's nodes: Gauss-Lobatto (the default), equispaced (`H1@U`), or
    // the legacy `Cubic` collection (equispaced, its own hexahedron interior).
    enum Points { Gll, Uniform, Cubic } mPoints = Gll;
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
    } else if (rName == "Cubic") {
        s.mKind = MfSpace::H1;
        s.mOrder = 3;
        s.mPoints = MfSpace::Cubic;
    } else if (rName == "QuadraticPos") {
        s.mKind = MfSpace::H1Other;
        s.mOrder = 2;
    } else if (rName.rfind("H1_", 0) == 0 || rName.rfind("H1@", 0) == 0) {
        s.mOrder = order_after_P();
        // Gauss-Lobatto (the default) and equispaced (`H1@U`) points are read at
        // any order; they coincide up to order 2, where any closed basis is.
        const char basis = rName[2] == '@' && rName.size() > 3 ? rName[3] : 'G';
        s.mPoints = basis == 'U' ? MfSpace::Uniform : MfSpace::Gll;
        const bool nodal = basis == 'G' || basis == 'U' || s.mOrder <= 2;
        s.mKind = s.mOrder >= 1 && nodal ? MfSpace::H1 : MfSpace::H1Other;
    } else if (rName.rfind("H1Pos_", 0) == 0 || rName.rfind("H1Ser_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = s.mOrder == 1 ? MfSpace::H1 : MfSpace::H1Other;
    } else if (rName.rfind("NURBS", 0) == 0) {
        // NURBS<p>, or NURBS alone for the orders of the mesh's knot vectors
        std::int64_t v = -1;
        if (rName.size() == 5 || MfLexer::ParseInt(rName.substr(5), v)) {
            s.mKind = MfSpace::Nurbs;
            s.mOrder = static_cast<int>(v);
        }
    } else if (rName.rfind("L2_T1_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = MfSpace::L2T1;
    } else if (rName.rfind("L2_", 0) == 0) {
        s.mOrder = order_after_P();
        s.mKind = MfSpace::L2;
    }
    return s;
}

// The collection, VDim and Ordering lines of a FiniteElementSpace header.
MfSpace mf_read_space_body(MfLexer& rLex, std::size_t Line);

MfSpace mf_read_space(MfLexer& rLex) {
    const MfToken& fes = rLex.Next("FiniteElementSpace");
    if (fes.mText != "FiniteElementSpace")
        rLex.Fail("expected FiniteElementSpace, found '" + fes.mText +
                      "' (NURBS and variable-order spaces are not supported)",
                  fes.mLine);
    return mf_read_space_body(rLex, fes.mLine);
}

// The `End: MFEM FiniteElementSpace v1.0` closing the versioned header MFEM
// writes for NURBS spaces (its variable-order element lists are not read).
void mf_read_space_end(MfLexer& rLex) {
    const MfToken& end = rLex.Next("'End:'");
    if (end.mText != "End:")
        rLex.Fail("'" + end.mText + "' in a versioned FiniteElementSpace (not supported)",
                  end.mLine);
    for (const char* word : {"MFEM", "FiniteElementSpace", "v1.0"}) {
        const MfToken& t = rLex.Next(word);
        if (t.mText != word)
            rLex.Fail(std::string("expected '") + word + "', found '" + t.mText + "'", t.mLine);
    }
}

MfSpace mf_read_space_body(MfLexer& rLex, std::size_t Line) {
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
        rLex.Fail("a FiniteElementSpace without a FiniteElementCollection", Line);
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

struct MfNurbs;

struct MfFile {
    std::shared_ptr<const MfNurbs> mpNurbs;  // a NURBS mesh: everything is there
    bool mNonConforming = false;             // read from an `MFEM NC mesh`: its leaves
    // One rank of a parallel mesh (ParMesh::Print): its rank (group 0's only
    // member), and per communication group its ranks and shared vertices
    // (local ids, in the order every rank of the group lists them).
    bool mParallel = false;
    std::int64_t mRank = 0;
    std::vector<std::vector<std::int64_t>> mGroups;
    std::vector<std::vector<std::int64_t>> mGroupVertices;
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

// An `MFEM NC mesh` (ncmesh.cpp's NCMesh::Print): the refinement tree, read as
// its leaf elements. Rows are `rank attribute geometry ref_type` then the
// vertices (a leaf) or the children (refined); `-1` after the attribute marks an
// unused slot. Top-level vertices are the `coordinates`; every other vertex is
// placed between its two `vertex_parents` (at the v1.1 scale, else halfway).
// Only the leaves of the file's own rank are kept.
MfFile mf_parse_nc(MfLexer& rLex, const std::string& rPath, bool Scaled) {
    const auto& geoms = mf_geoms();
    MfFile f;
    f.mNonConforming = true;
    struct NcElement {
        std::int64_t mRank = 0, mAttribute = 0;
        int mGeom = -1, mRefType = 0;
        std::vector<std::int64_t> mIds;
        std::size_t mLine = 0;
    };
    std::vector<NcElement> elements;
    std::map<std::int64_t, std::pair<std::pair<std::int64_t, std::int64_t>, double>> parents;
    std::vector<double> top;  // top-level coordinates, 3 per node
    std::size_t top_count = 0;
    std::int64_t my_rank = 0;
    int sdim = 0;
    bool have_coordinates = false;
    while (!rLex.AtEnd()) {
        const MfToken& t = rLex.Next("a section");
        if (t.mText == "dimension") {
            const std::int64_t d = rLex.Int("a dimension");
            if (d < 1 || d > 3)
                rLex.Fail("dimension " + std::to_string(d) + " (1, 2 or 3)", t.mLine);
            f.mDim = static_cast<int>(d);
        } else if (t.mText == "rank") {
            my_rank = rLex.Int("a rank");
        } else if (t.mText == "sfc_version") {
            rLex.Int("an SFC version");
        } else if (t.mText == "elements") {
            const std::int64_t n = rLex.Int("an element count");
            if (n < 0)
                rLex.Fail("negative element count", t.mLine);
            elements.resize(static_cast<std::size_t>(n));
            for (NcElement& el : elements) {
                el.mLine = rLex.Line();
                el.mRank = rLex.Int("a rank");
                el.mAttribute = rLex.Int("an attribute");
                const std::int64_t geom = rLex.Int("a geometry type");
                if (geom == -1)
                    continue;  // an unused slot
                if (geom < 1 || geom > 7)
                    rLex.Fail("unknown geometry type " + std::to_string(geom), el.mLine);
                el.mGeom = static_cast<int>(geom);
                el.mRefType = static_cast<int>(rLex.Int("a refinement type"));
                // The rest of the row: the vertices, or as many children as it lists.
                while (!rLex.AtEnd() && rLex.Peek().mLine == el.mLine)
                    el.mIds.push_back(rLex.Int("a node or child"));
                if (el.mRefType == 0 &&
                    el.mIds.size() != geoms[static_cast<std::size_t>(geom)].mNumVertices)
                    rLex.Fail(
                        "a leaf element lists " + std::to_string(el.mIds.size()) + " vertices",
                        el.mLine);
            }
        } else if (t.mText == "boundary") {
            f.mBoundary = mf_read_elements(rLex, "boundary element");
        } else if (t.mText == "vertex_parents") {
            const std::int64_t n = rLex.Int("a vertex count");
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t id = rLex.Int("a vertex");
                const std::int64_t p1 = rLex.Int("a parent");
                const std::int64_t p2 = rLex.Int("a parent");
                const double scale = Scaled ? rLex.Real("a scale") : 0.5;
                parents[id] = {{p1, p2}, scale};
            }
        } else if (t.mText == "root_state") {
            const std::int64_t n = rLex.Int("a root count");
            for (std::int64_t k = 0; k < n; ++k)
                rLex.Int("a root state");
        } else if (t.mText == "coordinates") {
            const std::int64_t n = rLex.Int("a vertex count");
            if (n < 0)
                rLex.Fail("negative vertex count", t.mLine);
            top_count = static_cast<std::size_t>(n);
            if (n > 0) {
                const std::int64_t sd = rLex.Int("a space dimension");
                if (sd < 1 || sd > 3)
                    rLex.Fail("space dimension " + std::to_string(sd) + " (1, 2 or 3)", t.mLine);
                sdim = static_cast<int>(sd);
                top.assign(top_count * 3, 0.0);
                for (std::size_t k = 0; k < top_count; ++k)
                    for (int c = 0; c < sdim; ++c)
                        top[3 * k + static_cast<std::size_t>(c)] = rLex.Real("a coordinate");
            }
            have_coordinates = true;
        } else if (t.mText == "nodes") {
            rLex.Fail("curved non-conforming meshes (a 'nodes' section) are not supported",
                      t.mLine);
        } else if (t.mText == "mfem_mesh_end" || t.mText == "mfem_serial_mesh_end") {
            break;
        } else {
            rLex.Fail("unexpected '" + t.mText + "'", t.mLine);
        }
    }
    if (f.mDim < 0)
        throw ReadError("MFEM mesh: no dimension section in " + rPath);
    if (!have_coordinates)
        throw ReadError("MFEM mesh: the non-conforming mesh " + rPath +
                        " has no top-level coordinates");

    // Leaves, depth first from the roots (elements no other one lists as a child).
    std::vector<bool> is_child(elements.size(), false);
    for (const NcElement& el : elements)
        if (el.mGeom > 0 && el.mRefType != 0)
            for (std::int64_t c : el.mIds) {
                if (c < 0 || static_cast<std::size_t>(c) >= elements.size())
                    throw ReadError("MFEM mesh: child element " + std::to_string(c) +
                                    " out of range (line " + std::to_string(el.mLine) + ")");
                is_child[static_cast<std::size_t>(c)] = true;
            }
    std::vector<std::size_t> leaves;
    std::size_t ghosts = 0;
    std::vector<std::size_t> stack;
    for (std::size_t r = elements.size(); r-- > 0;)
        if (elements[r].mGeom > 0 && !is_child[r])
            stack.push_back(r);
    std::vector<bool> seen(elements.size(), false);
    while (!stack.empty()) {
        const std::size_t e = stack.back();
        stack.pop_back();
        if (seen[e])
            throw ReadError("MFEM mesh: element " + std::to_string(e) +
                            " is reached twice in the refinement tree");
        seen[e] = true;
        const NcElement& el = elements[e];
        if (el.mRefType == 0) {
            if (el.mRank == my_rank)
                leaves.push_back(e);
            else
                ++ghosts;
            continue;
        }
        for (std::size_t k = el.mIds.size(); k-- > 0;)
            stack.push_back(static_cast<std::size_t>(el.mIds[k]));
    }
    if (ghosts)
        log::warn("MFEM mesh: {} ghost element(s) of other ranks in {} dropped", ghosts, rPath);

    // Vertices: every node a leaf or boundary element uses, by node id.
    std::set<std::int64_t> used;
    for (std::size_t e : leaves)
        used.insert(elements[e].mIds.begin(), elements[e].mIds.end());
    for (const MfElement& b : f.mBoundary)
        used.insert(b.mVertices.begin(), b.mVertices.end());
    std::map<std::int64_t, std::int64_t> index;
    for (std::int64_t id : used)
        index.emplace(id, static_cast<std::int64_t>(index.size()));
    std::map<std::int64_t, std::array<double, 3>> pos;
    std::set<std::int64_t> visiting;
    std::function<std::array<double, 3>(std::int64_t)> position = [&](std::int64_t Id) {
        const auto found = pos.find(Id);
        if (found != pos.end())
            return found->second;
        std::array<double, 3> x = {0, 0, 0};
        const auto par = parents.find(Id);
        if (par == parents.end()) {
            if (Id < 0 || static_cast<std::size_t>(Id) >= top_count)
                throw ReadError("MFEM mesh: vertex " + std::to_string(Id) +
                                " has neither coordinates nor parents");
            for (std::size_t c = 0; c < 3; ++c)
                x[c] = top[3 * static_cast<std::size_t>(Id) + c];
        } else {
            if (!visiting.insert(Id).second)
                throw ReadError("MFEM mesh: cyclic vertex parents at vertex " + std::to_string(Id));
            const auto a = position(par->second.first.first);
            const auto b = position(par->second.first.second);
            const double s = par->second.second;
            for (std::size_t c = 0; c < 3; ++c)
                x[c] = (1.0 - s) * a[c] + s * b[c];
        }
        pos.emplace(Id, x);
        return x;
    };
    f.mSpaceDim = sdim > 0 ? sdim : f.mDim;
    f.mNumVertices = index.size();
    for (const auto& [id, k] : index) {
        const auto x = position(id);
        for (int c = 0; c < f.mSpaceDim; ++c)
            f.mCoords.push_back(x[static_cast<std::size_t>(c)]);
    }
    for (std::size_t e : leaves) {
        const NcElement& el = elements[e];
        MfElement out{el.mAttribute, el.mGeom, {}, el.mLine};
        for (std::int64_t id : el.mIds)
            out.mVertices.push_back(index.at(id));
        f.mElements.push_back(std::move(out));
    }
    for (MfElement& b : f.mBoundary)
        for (std::int64_t& v : b.mVertices)
            v = index.at(v);
    log::warn(
        "MFEM mesh: the non-conforming mesh {} is read as its {} leaf element(s); hanging "
        "nodes are left unconstrained",
        rPath, leaves.size());
    return f;
}

// A parallel mesh's `communication_groups` (GroupTopology::Save) and shared
// entities (ParMesh::ParPrint), after `mfem_serial_mesh_end`. Only the shared
// vertices matter: shared edges and faces follow from them.
void mf_read_groups(MfLexer& rLex, MfFile& rF) {
    rF.mParallel = true;
    const MfToken& head = rLex.Next("communication_groups");
    if (head.mText != "communication_groups")
        rLex.Fail("expected communication_groups, found '" + head.mText + "'", head.mLine);
    if (rLex.Next("number_of_groups").mText != "number_of_groups")
        rLex.Fail("expected number_of_groups", head.mLine);
    const std::int64_t ngroups = rLex.Int("a group count");
    if (ngroups < 1)
        rLex.Fail("a parallel mesh needs at least one communication group", head.mLine);
    for (std::int64_t g = 0; g < ngroups; ++g) {
        const std::int64_t size = rLex.Int("a group size");
        if (size < 1)
            rLex.Fail("an empty communication group", rLex.Line());
        std::vector<std::int64_t> ranks;
        for (std::int64_t k = 0; k < size; ++k)
            ranks.push_back(rLex.Int("a rank"));
        std::sort(ranks.begin(), ranks.end());
        rF.mGroups.push_back(std::move(ranks));
    }
    if (rF.mGroups[0].size() != 1)
        rLex.Fail("communication group 0 must hold this rank alone", head.mLine);
    rF.mRank = rF.mGroups[0][0];
    rF.mGroupVertices.assign(static_cast<std::size_t>(ngroups), {});
    std::int64_t group = 0;
    while (!rLex.AtEnd()) {
        const MfToken& t = rLex.Next("a shared-entity section");
        if (t.mText == "mfem_mesh_end")
            break;
        if (t.mText == "total_shared_vertices" || t.mText == "total_shared_edges" ||
            t.mText == "total_shared_faces") {
            rLex.Int("a count");
        } else if (t.mText == "shared_vertices") {
            if (++group >= ngroups)
                rLex.Fail("more shared-vertex groups than communication groups", t.mLine);
            const std::int64_t n = rLex.Int("a vertex count");
            for (std::int64_t k = 0; k < n; ++k)
                rF.mGroupVertices[static_cast<std::size_t>(group)].push_back(
                    rLex.Int("a shared vertex"));
        } else if (t.mText == "shared_edges") {
            const std::int64_t n = rLex.Int("an edge count");
            for (std::int64_t k = 0; k < 2 * n; ++k)
                rLex.Int("an edge vertex");
        } else if (t.mText == "shared_faces") {
            const std::int64_t n = rLex.Int("a face count");
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t geom = rLex.Int("a face geometry");
                for (int v = 0; v < (geom == 3 ? 4 : 3); ++v)
                    rLex.Int("a face vertex");
            }
        } else {
            rLex.Fail("unexpected '" + t.mText + "' in the communication groups", t.mLine);
        }
    }
}

// --- NURBS meshes ----------------------------------------------------------------
//
// `MFEM NURBS mesh v1.0`/`v1.1`: the patch topology, its knot vectors and
// control points, numbered the way MFEM's NURBSExtension numbers them
// (mesh/nurbs.cpp, BSD-3-Clause), so the knot-span elements, their vertices
// and every patch's control points are MFEM's. The global form keeps every
// control point once: the topological vertices, then the interior points of
// every edge, face and patch; each patch reaches its own through
// NURBSPatchMap (the patch's vertices, its edges and faces with their
// orientations, its interior block). The `patches` form lists each patch's
// control points itself, merged into the same numbering. Twin of _nurbs.py.

constexpr int kMfQuadEdges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
constexpr int kMfHexEdges[12][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {4, 5}, {5, 6},
                                    {7, 6}, {4, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
constexpr int kMfHexFaces[6][4] = {{3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5},
                                   {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};

int mf_nurbs_geom(int Dim) {
    return Dim == 1 ? 1 : (Dim == 2 ? 3 : 5);
}
int mf_nurbs_bdr_geom(int Dim) {
    return Dim == 1 ? 0 : (Dim == 2 ? 1 : 3);
}
std::int64_t mf_flip_sign(std::int64_t K) {
    return -1 - K;
}
std::int64_t mf_unsign(std::int64_t K) {
    return K >= 0 ? K : -1 - K;
}

struct MfKnotVector {
    int mOrder = 0;
    std::size_t mNumCp = 0;
    std::vector<double> mKnots;
    std::vector<std::size_t> mSpans;  // the knot spans of nonzero length: the elements

    void Finish() {
        mNumCp = mKnots.size() - static_cast<std::size_t>(mOrder) - 1;
        mSpans.clear();
        for (std::size_t s = static_cast<std::size_t>(mOrder); s < mNumCp; ++s)
            if (mKnots[s + 1] > mKnots[s])
                mSpans.push_back(s);
    }
    std::size_t NumElements() const { return mSpans.size(); }
    MfKnotVector Flipped() const {
        MfKnotVector k = *this;
        const double apb = mKnots.front() + mKnots.back();
        for (std::size_t i = 0; i < mKnots.size(); ++i)
            k.mKnots[i] = apb - mKnots[mKnots.size() - 1 - i];
        k.Finish();
        return k;
    }
    // The order+1 B-spline values at reference R in [0, 1] of knot span Span
    // (Cox-de Boor), for control points Span-order..Span.
    void Basis(std::size_t Span, double R, std::vector<double>& rN) const {
        const std::size_t p = static_cast<std::size_t>(mOrder);
        const double u = mKnots[Span] + R * (mKnots[Span + 1] - mKnots[Span]);
        rN.assign(p + 1, 0.0);
        rN[0] = 1.0;
        std::vector<double> left(p + 1, 0.0), right(p + 1, 0.0);
        for (std::size_t j = 1; j <= p; ++j) {
            left[j] = u - mKnots[Span + 1 - j];
            right[j] = mKnots[Span + j] - u;
            double saved = 0.0;
            for (std::size_t q = 0; q < j; ++q) {
                const double t = rN[q] / (right[q + 1] + left[j - q]);
                rN[q] = saved + right[q + 1] * t;
                saved = left[j - q] * t;
            }
            rN[j] = saved;
        }
    }
};

MfKnotVector mf_read_knot(MfLexer& rLex) {
    const std::size_t line = rLex.Line();
    MfKnotVector k;
    const std::int64_t order = rLex.Int("a knot vector order");
    const std::int64_t ncp = rLex.Int("a control point count");
    if (order < 0 || ncp < order + 1)
        rLex.Fail("knot vector of order " + std::to_string(order) + " with " + std::to_string(ncp) +
                      " control points",
                  line);
    k.mOrder = static_cast<int>(order);
    for (std::int64_t i = 0; i < ncp + order + 1; ++i)
        k.mKnots.push_back(rLex.Real("a knot"));
    for (std::size_t i = 1; i < k.mKnots.size(); ++i)
        if (k.mKnots[i] < k.mKnots[i - 1])
            rLex.Fail("knots out of order", line);
    k.Finish();
    return k;
}

// Mesh::GetQuadOrientation.
int mf_quad_orientation(const std::int64_t* pBase, const std::int64_t* pTest) {
    int i = 0;
    while (i < 3 && pTest[i] != pBase[0])
        ++i;
    return pTest[(i + 1) % 4] == pBase[1] ? 2 * i : 2 * i + 1;
}

std::int64_t mf_or1d(std::int64_t N, std::int64_t BigN, int Or) {
    return Or > 0 ? N : BigN - 1 - N;
}

std::int64_t mf_or2d(std::int64_t N1, std::int64_t N2, std::int64_t B1, std::int64_t B2, int Or) {
    switch (Or) {
        case 0:
            return N1 + N2 * B1;
        case 1:
            return N2 + N1 * B2;
        case 2:
            return N2 + (B1 - 1 - N1) * B2;
        case 3:
            return (B1 - 1 - N1) + N2 * B1;
        case 4:
            return (B1 - 1 - N1) + (B2 - 1 - N2) * B1;
        case 5:
            return (B2 - 1 - N2) + (B1 - 1 - N1) * B2;
        case 6:
            return (B2 - 1 - N2) + N1 * B2;
        default:
            return N1 + (B2 - 1 - N2) * B1;
    }
}

int mf_f(std::int64_t N, std::int64_t BigN) {
    return N < 0 ? 0 : (N >= BigN ? 2 : 1);
}

// NURBSPatchMap: patch lattice index -> global number, of the mesh vertices or
// of the control points.
struct MfPatchMap {
    std::vector<std::int64_t> mVerts, mEdges, mFaces;
    std::vector<int> mOEdge, mOFace;
    std::int64_t mPOffset = 0;
    std::int64_t mN[3] = {0, 0, 0};  // interior counts per direction
    int mOPatch = 0;

    std::int64_t Ec(int E, std::int64_t N, std::int64_t BigN, int S = 1) const {
        return mEdges[static_cast<std::size_t>(E)] +
               mf_or1d(N, BigN, S * mOEdge[static_cast<std::size_t>(E)]);
    }
    std::int64_t Fc(int F, std::int64_t M, std::int64_t N, std::int64_t BigM,
                    std::int64_t BigN) const {
        return mFaces[static_cast<std::size_t>(F)] +
               mf_or2d(M, N, BigM, BigN, mOFace[static_cast<std::size_t>(F)]);
    }
    std::int64_t operator()(std::int64_t I) const {
        const std::int64_t i1 = I - 1;
        switch (mf_f(i1, mN[0])) {
            case 0:
                return mVerts[0];
            case 2:
                return mVerts[1];
            default:
                return mPOffset + mf_or1d(i1, mN[0], mOPatch);
        }
    }
    std::int64_t operator()(std::int64_t I, std::int64_t J) const {
        const std::int64_t i1 = I - 1, j1 = J - 1, bi = mN[0], bj = mN[1];
        switch (3 * mf_f(j1, bj) + mf_f(i1, bi)) {
            case 0:
                return mVerts[0];
            case 1:
                return Ec(0, i1, bi);
            case 2:
                return mVerts[1];
            case 3:
                return Ec(3, j1, bj, -1);
            case 4:
                return mPOffset + mf_or2d(i1, j1, bi, bj, mOPatch);
            case 5:
                return Ec(1, j1, bj);
            case 6:
                return mVerts[3];
            case 7:
                return Ec(2, i1, bi, -1);
            default:
                return mVerts[2];
        }
    }
    std::int64_t operator()(std::int64_t I, std::int64_t J, std::int64_t K) const {
        const std::int64_t i1 = I - 1, j1 = J - 1, k1 = K - 1;
        const std::int64_t bi = mN[0], bj = mN[1], bk = mN[2];
        switch (3 * (3 * mf_f(k1, bk) + mf_f(j1, bj)) + mf_f(i1, bi)) {
            case 0:
                return mVerts[0];
            case 1:
                return Ec(0, i1, bi);
            case 2:
                return mVerts[1];
            case 3:
                return Ec(3, j1, bj);
            case 4:
                return Fc(0, i1, bj - 1 - j1, bi, bj);
            case 5:
                return Ec(1, j1, bj);
            case 6:
                return mVerts[3];
            case 7:
                return Ec(2, i1, bi);
            case 8:
                return mVerts[2];
            case 9:
                return Ec(8, k1, bk);
            case 10:
                return Fc(1, i1, k1, bi, bk);
            case 11:
                return Ec(9, k1, bk);
            case 12:
                return Fc(4, bj - 1 - j1, k1, bj, bk);
            case 13:
                return mPOffset + bi * (bj * k1 + j1) + i1;
            case 14:
                return Fc(2, j1, k1, bj, bk);
            case 15:
                return Ec(11, k1, bk);
            case 16:
                return Fc(3, bi - 1 - i1, k1, bi, bk);
            case 17:
                return Ec(10, k1, bk);
            case 18:
                return mVerts[4];
            case 19:
                return Ec(4, i1, bi);
            case 20:
                return mVerts[5];
            case 21:
                return Ec(7, j1, bj);
            case 22:
                return Fc(5, i1, j1, bi, bj);
            case 23:
                return Ec(5, j1, bj);
            case 24:
                return mVerts[7];
            case 25:
                return Ec(6, i1, bi);
            default:
                return mVerts[6];
        }
    }
    std::int64_t At(const std::int64_t* pIdx, int Dim) const {
        return Dim == 1
                   ? (*this)(pIdx[0])
                   : (Dim == 2 ? (*this)(pIdx[0], pIdx[1]) : (*this)(pIdx[0], pIdx[1], pIdx[2]));
    }
};

// A knot-span element: its patch and span indices.
struct MfNurbsElement {
    MfElement mElement;
    std::size_t mPatch = 0;
    std::int64_t mSpan[3] = {0, 0, 0};
};

struct MfNurbsPatchData {
    std::vector<MfKnotVector> mKvs;
    int mSpaceDim = 0;
    std::vector<double> mCps;  // homogeneous, (SpaceDim + 1) per point
};

struct MfNurbs {
    std::string mPath;
    int mDim = -1;
    std::vector<MfElement> mPatches, mBPatches;
    std::vector<std::array<std::int64_t, 3>> mEdgeRows;  // (kv, v0, v1) as in the file
    std::size_t mNumTopoVertices = 0;
    std::vector<MfKnotVector> mKvs;
    bool mPatchForm = false;
    std::vector<MfNurbsPatchData> mPatchData;
    std::vector<double> mWeights;
    bool mHasNodes = false;
    MfSpace mSpace;
    std::vector<double> mNodes;

    // topology
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> mEdgeOf;
    std::vector<std::int64_t> mEdgeToUkv;
    std::vector<std::vector<std::size_t>> mPEdges;
    std::vector<std::vector<int>> mPEdgeOr;
    std::vector<std::vector<std::size_t>> mPFaces;
    std::vector<std::vector<int>> mPFaceOr;
    std::vector<std::array<std::int64_t, 4>> mFaceVerts;
    std::map<std::array<std::int64_t, 4>, std::size_t> mFaceIndex;
    std::vector<std::vector<MfKnotVector>> mCompr;
    // numbering: edge, face and patch offsets and the total, of the mesh
    // vertices [0] and of the control points [1]
    std::vector<std::int64_t> mEOff[2], mFOff[2], mPOff[2];
    std::int64_t mTotal[2] = {0, 0};

    std::int64_t NumVertices() const { return mTotal[0]; }
    std::int64_t NumDofs() const { return mTotal[1]; }
    [[noreturn]] void Fail(const std::string& rWhy) const {
        throw ReadError("MFEM NURBS mesh: " + rWhy + " (" + mPath + ")");
    }

    const MfKnotVector& KvOfEdge(std::size_t E) const {
        return mKvs[static_cast<std::size_t>(mf_unsign(mEdgeToUkv[E]))];
    }
    int KvSign(std::size_t E) const { return mEdgeToUkv[E] >= 0 ? 1 : -1; }
    std::size_t EdgeOf(std::int64_t A, std::int64_t B) const {
        const auto it = mEdgeOf.find({std::min(A, B), std::max(A, B)});
        if (it == mEdgeOf.end())
            Fail("the edge " + std::to_string(A) + "-" + std::to_string(B) +
                 " is not in the edges section");
        return it->second;
    }
    std::vector<std::size_t> DirEdges(std::size_t P) const {
        if (mDim == 1)
            return {P};
        const auto& es = mPEdges[P];
        if (mDim == 2)
            return {es[0], es[1]};
        return {es[0], es[3], es[8]};
    }
    // NURBSExtension::CheckKVDirection
    std::vector<int> KvDir(std::size_t P) const {
        if (mDim == 1)
            return {KvSign(P)};
        const auto& pv = mPatches[P].mVertices;
        std::vector<std::pair<std::int64_t, std::int64_t>> pairs = {{pv[0], pv[1]}, {pv[0], pv[3]}};
        if (mDim == 3)
            pairs.push_back({pv[0], pv[4]});
        std::vector<int> kvdir(static_cast<std::size_t>(mDim), 0);
        for (std::size_t e : mPEdges[P]) {
            const std::int64_t a = mEdgeRows[e][1], b = mEdgeRows[e][2];
            const int ks = KvSign(e);
            for (std::size_t d = 0; d < pairs.size(); ++d) {
                if (a == pairs[d].first && b == pairs[d].second)
                    kvdir[d] = ks;
                else if (a == pairs[d].second && b == pairs[d].first)
                    kvdir[d] = -ks;
            }
        }
        return kvdir;
    }

    void Topology();
    void BoundaryPatches();
    void KnotVectors();
    void Offsets();
    MfPatchMap PatchMap(std::size_t P, bool Space) const;
    MfPatchMap BdrPatchMap(std::size_t B, std::vector<int>& rOkv,
                           std::vector<const MfKnotVector*>& rKvs) const;
    std::vector<MfNurbsElement> Elements() const;
    std::vector<MfElement> Boundary() const;
    void ControlPoints(std::vector<double>& rWeights, std::vector<double>& rXyz,
                       int& rSpaceDim) const;
    void Evaluate(const std::vector<double>& rWeights, const std::vector<double>& rTable,
                  std::size_t NumComp, const MfNurbsElement& rEl,
                  const std::vector<std::array<double, 3>>& rRefs, std::vector<double>& rOut) const;
};

void MfNurbs::Topology() {
    const std::size_t np = mPatches.size();
    if (mDim == 1) {
        if (mEdgeRows.size() != np)
            Fail(std::to_string(mEdgeRows.size()) + " edges for " + std::to_string(np) +
                 " patches");
        for (const auto& r : mEdgeRows)
            mEdgeToUkv.push_back(r[1] <= r[2] ? r[0] : mf_flip_sign(r[0]));
    } else {
        for (std::size_t e = 0; e < mEdgeRows.size(); ++e) {
            const auto& r = mEdgeRows[e];
            mEdgeOf[{std::min(r[1], r[2]), std::max(r[1], r[2])}] = e;
            mEdgeToUkv.push_back(r[1] <= r[2] ? r[0] : mf_flip_sign(r[0]));
        }
    }
    mPEdges.assign(np, {});
    mPEdgeOr.assign(np, {});
    if (mDim > 1) {
        for (std::size_t p = 0; p < np; ++p) {
            const auto& v = mPatches[p].mVertices;
            const std::size_t ne = mDim == 2 ? 4 : 12;
            for (std::size_t j = 0; j < ne; ++j) {
                const int a = mDim == 2 ? kMfQuadEdges[j][0] : kMfHexEdges[j][0];
                const int b = mDim == 2 ? kMfQuadEdges[j][1] : kMfHexEdges[j][1];
                const std::int64_t va = v[static_cast<std::size_t>(a)];
                const std::int64_t vb = v[static_cast<std::size_t>(b)];
                mPEdges[p].push_back(EdgeOf(va, vb));
                mPEdgeOr[p].push_back(va < vb ? 1 : -1);
            }
        }
    }
    mPFaces.assign(np, {});
    mPFaceOr.assign(np, {});
    if (mDim == 3) {
        for (std::size_t p = 0; p < np; ++p) {
            const auto& v = mPatches[p].mVertices;
            for (const auto& fv : kMfHexFaces) {
                std::array<std::int64_t, 4> verts;
                for (std::size_t c = 0; c < 4; ++c)
                    verts[c] = v[static_cast<std::size_t>(fv[c])];
                std::array<std::int64_t, 4> key = verts;
                std::sort(key.begin(), key.end());
                const auto [it, fresh] = mFaceIndex.emplace(key, mFaceVerts.size());
                if (fresh) {
                    mFaceVerts.push_back(verts);
                    mPFaceOr[p].push_back(0);
                } else {
                    mPFaceOr[p].push_back(
                        mf_quad_orientation(mFaceVerts[it->second].data(), verts.data()));
                }
                mPFaces[p].push_back(it->second);
            }
        }
    }
    BoundaryPatches();
    KnotVectors();
}

// FinalizeTopology and CheckBdrElementOrientation: without a boundary
// section, the faces of one patch become the boundary (attribute 1); a
// boundary patch is turned to run as the face does in the first patch holding
// it.
void MfNurbs::BoundaryPatches() {
    using Face = std::vector<std::int64_t>;
    std::map<Face, std::pair<Face, int>> stored;  // sorted -> (first patch's order, count)
    std::vector<Face> order;
    for (const MfElement& el : mPatches) {
        const auto& v = el.mVertices;
        std::vector<Face> faces;
        if (mDim == 1) {
            faces = {{v[0]}, {v[1]}};
        } else if (mDim == 2) {
            for (const auto& e : kMfQuadEdges)
                faces.push_back(
                    {v[static_cast<std::size_t>(e[0])], v[static_cast<std::size_t>(e[1])]});
        } else {
            for (const auto& fv : kMfHexFaces) {
                Face f;
                for (int c : fv)
                    f.push_back(v[static_cast<std::size_t>(c)]);
                faces.push_back(f);
            }
        }
        for (const Face& f : faces) {
            Face key = f;
            std::sort(key.begin(), key.end());
            const auto [it, fresh] = stored.emplace(key, std::make_pair(f, 0));
            if (fresh)
                order.push_back(key);
            ++it->second.second;
        }
    }
    if (mBPatches.empty()) {
        if (mDim == 1)
            std::sort(order.begin(), order.end());  // 1-D faces are the vertices
        else if (mDim == 2)
            std::sort(order.begin(), order.end(), [&](const Face& rA, const Face& rB) {
                return EdgeOf(rA[0], rA[1]) < EdgeOf(rB[0], rB[1]);
            });
        for (const Face& k : order) {
            const auto& [verts, count] = stored.at(k);
            if (count != 1)
                continue;
            MfElement b;
            b.mAttribute = 1;
            b.mGeom = mf_nurbs_bdr_geom(mDim);
            b.mVertices = verts;
            b.mLine = 0;
            mBPatches.push_back(std::move(b));
        }
        return;
    }
    for (MfElement& b : mBPatches) {
        Face key = b.mVertices;
        std::sort(key.begin(), key.end());
        const auto it = stored.find(key);
        if (it == stored.end() || it->second.second != 1)
            continue;
        const Face& fv = it->second.first;
        auto& bv = b.mVertices;
        if (mDim == 2 && bv[0] != fv[0])
            std::swap(bv[0], bv[1]);
        else if (mDim == 3 && mf_quad_orientation(fv.data(), bv.data()) % 2)
            std::swap(bv[0], bv[2]);
    }
}

void MfNurbs::KnotVectors() {
    const std::size_t np = mPatches.size();
    if (mPatchForm) {
        std::int64_t nkv = 0;
        for (std::int64_t k : mEdgeToUkv)
            nkv = std::max(nkv, mf_unsign(k) + 1);
        std::vector<bool> have(static_cast<std::size_t>(nkv), false);
        mKvs.assign(static_cast<std::size_t>(nkv), MfKnotVector{});
        for (std::size_t p = 0; p < np; ++p) {
            const std::vector<int> kvdir = KvDir(p);
            const std::vector<std::size_t> dirs = DirEdges(p);
            for (std::size_t d = 0; d < dirs.size(); ++d) {
                const auto k = static_cast<std::size_t>(mf_unsign(mEdgeToUkv[dirs[d]]));
                if (have[k])
                    continue;
                const MfKnotVector& kv = mPatchData[p].mKvs[d];
                mKvs[k] = kvdir[d] == -1 ? kv.Flipped() : kv;
                have[k] = true;
            }
        }
        for (bool h : have)
            if (!h)
                Fail("a knot vector no patch defines");
    }
    for (std::int64_t k : mEdgeToUkv)
        if (static_cast<std::size_t>(mf_unsign(k)) >= mKvs.size())
            Fail("knot vector " + std::to_string(mf_unsign(k)) + " is not defined");
    mCompr.assign(np, {});
    for (std::size_t p = 0; p < np; ++p) {
        const std::vector<int> kvdir = KvDir(p);
        const std::vector<std::size_t> dirs = DirEdges(p);
        for (std::size_t d = 0; d < dirs.size(); ++d) {
            const MfKnotVector& kv = KvOfEdge(dirs[d]);
            mCompr[p].push_back(kvdir[d] == -1 ? kv.Flipped() : kv);
        }
    }
}

// NURBSExtension::GenerateOffsets
void MfNurbs::Offsets() {
    for (int kind = 0; kind < 2; ++kind) {
        const auto count = [kind](const MfKnotVector& rKv) {
            return kind == 0 ? static_cast<std::int64_t>(rKv.NumElements()) - 1
                             : static_cast<std::int64_t>(rKv.mNumCp) - 2;
        };
        std::int64_t n = static_cast<std::int64_t>(mNumTopoVertices);
        mEOff[kind].clear();
        mFOff[kind].clear();
        mPOff[kind].clear();
        if (mDim > 1)
            for (std::size_t e = 0; e < mEdgeRows.size(); ++e) {
                mEOff[kind].push_back(n);
                n += count(KvOfEdge(e));
            }
        for (const auto& fv : mFaceVerts) {
            mFOff[kind].push_back(n);
            n += count(KvOfEdge(EdgeOf(fv[0], fv[1]))) * count(KvOfEdge(EdgeOf(fv[1], fv[2])));
        }
        for (std::size_t p = 0; p < mPatches.size(); ++p) {
            mPOff[kind].push_back(n);
            std::int64_t size = 1;
            for (std::size_t e : DirEdges(p))
                size *= count(KvOfEdge(e));
            n += size;
        }
        mTotal[kind] = n;
    }
}

MfPatchMap MfNurbs::PatchMap(std::size_t P, bool Space) const {
    const int kind = Space ? 1 : 0;
    MfPatchMap m;
    m.mVerts = mPatches[P].mVertices;
    for (std::size_t d = 0; d < mCompr[P].size(); ++d) {
        const MfKnotVector& kv = mCompr[P][d];
        m.mN[d] = Space ? static_cast<std::int64_t>(kv.mNumCp) - 2
                        : static_cast<std::int64_t>(kv.NumElements()) - 1;
    }
    if (mDim > 1) {
        for (std::size_t e : mPEdges[P])
            m.mEdges.push_back(mEOff[kind][e]);
        m.mOEdge = mPEdgeOr[P];
    }
    if (mDim == 3) {
        for (std::size_t f : mPFaces[P])
            m.mFaces.push_back(mFOff[kind][f]);
        m.mOFace = mPFaceOr[P];
    }
    m.mPOffset = mPOff[kind][P];
    m.mOPatch = 0;
    return m;
}

// SetBdrPatchVertexMap: the map, the knot-vector orientations and the knot
// vectors of boundary patch B.
MfPatchMap MfNurbs::BdrPatchMap(std::size_t B, std::vector<int>& rOkv,
                                std::vector<const MfKnotVector*>& rKvs) const {
    MfPatchMap m;
    m.mVerts = mBPatches[B].mVertices;
    rOkv.clear();
    rKvs.clear();
    if (mDim == 1) {
        rOkv.push_back(1);
        return m;
    }
    const auto& v = m.mVerts;
    if (mDim == 2) {
        const std::size_t e = EdgeOf(v[0], v[1]);
        const int o = v[0] < v[1] ? 1 : -1;
        rKvs.push_back(&KvOfEdge(e));
        rOkv.push_back(mEdgeToUkv[e] >= 0 ? o : -o);
        m.mPOffset = mEOff[0][e];
        m.mN[0] = static_cast<std::int64_t>(KvOfEdge(e).NumElements()) - 1;
        m.mOPatch = o;
        return m;
    }
    std::vector<std::size_t> es;
    for (const auto& ed : kMfQuadEdges) {
        const std::int64_t a = v[static_cast<std::size_t>(ed[0])];
        const std::int64_t c = v[static_cast<std::size_t>(ed[1])];
        es.push_back(EdgeOf(a, c));
        m.mOEdge.push_back(a < c ? 1 : -1);
        m.mEdges.push_back(mEOff[0][es.back()]);
    }
    for (std::size_t d = 0; d < 2; ++d) {
        rKvs.push_back(&KvOfEdge(es[d]));
        rOkv.push_back(mEdgeToUkv[es[d]] >= 0 ? m.mOEdge[d] : -m.mOEdge[d]);
        m.mN[d] = static_cast<std::int64_t>(rKvs.back()->NumElements()) - 1;
    }
    std::array<std::int64_t, 4> key = {v[0], v[1], v[2], v[3]};
    std::sort(key.begin(), key.end());
    const auto it = mFaceIndex.find(key);
    if (it == mFaceIndex.end())
        Fail("a boundary patch is not a face of the mesh");
    m.mOPatch = mf_quad_orientation(mFaceVerts[it->second].data(), v.data());
    m.mPOffset = mFOff[0][it->second];
    return m;
}

std::vector<MfNurbsElement> MfNurbs::Elements() const {
    std::vector<MfNurbsElement> out;
    const int geom = mf_nurbs_geom(mDim);
    for (std::size_t p = 0; p < mPatches.size(); ++p) {
        const MfPatchMap m = PatchMap(p, false);
        const auto& kvs = mCompr[p];
        const std::int64_t nx = static_cast<std::int64_t>(kvs[0].NumElements());
        const std::int64_t ny = mDim > 1 ? static_cast<std::int64_t>(kvs[1].NumElements()) : 1;
        const std::int64_t nz = mDim > 2 ? static_cast<std::int64_t>(kvs[2].NumElements()) : 1;
        for (std::int64_t k = 0; k < nz; ++k)
            for (std::int64_t j = 0; j < ny; ++j)
                for (std::int64_t i = 0; i < nx; ++i) {
                    MfNurbsElement el;
                    el.mPatch = p;
                    el.mSpan[0] = i;
                    el.mSpan[1] = j;
                    el.mSpan[2] = k;
                    el.mElement.mAttribute = mPatches[p].mAttribute;
                    el.mElement.mGeom = geom;
                    el.mElement.mLine = mPatches[p].mLine;
                    auto& v = el.mElement.mVertices;
                    if (mDim == 1) {
                        v = {m(i), m(i + 1)};
                    } else if (mDim == 2) {
                        v = {m(i, j), m(i + 1, j), m(i + 1, j + 1), m(i, j + 1)};
                    } else {
                        v = {m(i, j, k),
                             m(i + 1, j, k),
                             m(i + 1, j + 1, k),
                             m(i, j + 1, k),
                             m(i, j, k + 1),
                             m(i + 1, j, k + 1),
                             m(i + 1, j + 1, k + 1),
                             m(i, j + 1, k + 1)};
                    }
                    out.push_back(std::move(el));
                }
    }
    return out;
}

std::vector<MfElement> MfNurbs::Boundary() const {
    std::vector<MfElement> out;
    const int geom = mf_nurbs_bdr_geom(mDim);
    std::vector<int> okv;
    std::vector<const MfKnotVector*> kvs;
    for (std::size_t b = 0; b < mBPatches.size(); ++b) {
        const MfPatchMap m = BdrPatchMap(b, okv, kvs);
        MfElement el;
        el.mAttribute = mBPatches[b].mAttribute;
        el.mGeom = geom;
        el.mLine = mBPatches[b].mLine;
        if (mDim == 1) {
            el.mVertices = {m(0)};
            out.push_back(el);
        } else if (mDim == 2) {
            const std::int64_t nx = static_cast<std::int64_t>(kvs[0]->NumElements());
            for (std::int64_t i = 0; i < nx; ++i) {
                const std::int64_t i_ = okv[0] >= 0 ? i : nx - 1 - i;
                el.mVertices = {m(i_), m(i_ + 1)};
                out.push_back(el);
            }
        } else {
            const std::int64_t nx = static_cast<std::int64_t>(kvs[0]->NumElements());
            const std::int64_t ny = static_cast<std::int64_t>(kvs[1]->NumElements());
            for (std::int64_t j = 0; j < ny; ++j) {
                const std::int64_t j_ = okv[1] >= 0 ? j : ny - 1 - j;
                for (std::int64_t i = 0; i < nx; ++i) {
                    const std::int64_t i_ = okv[0] >= 0 ? i : nx - 1 - i;
                    el.mVertices = {m(i_, j_), m(i_ + 1, j_), m(i_ + 1, j_ + 1), m(i_, j_ + 1)};
                    out.push_back(el);
                }
            }
        }
    }
    return out;
}

// The weights and Cartesian control points (NumDofs x SpaceDim) in MFEM's
// global numbering.
void MfNurbs::ControlPoints(std::vector<double>& rWeights, std::vector<double>& rXyz,
                            int& rSpaceDim) const {
    const auto ndofs = static_cast<std::size_t>(NumDofs());
    if (mPatchForm) {
        rSpaceDim = mPatchData[0].mSpaceDim;
        const auto sd = static_cast<std::size_t>(rSpaceDim);
        rWeights.assign(ndofs, std::numeric_limits<double>::quiet_NaN());
        rXyz.assign(ndofs * sd, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t p = 0; p < mPatchData.size(); ++p) {
            const MfNurbsPatchData& pd = mPatchData[p];
            if (pd.mSpaceDim != rSpaceDim)
                Fail("patches of different dimensions");
            std::int64_t n[3] = {1, 1, 1};
            for (std::size_t d = 0; d < pd.mKvs.size(); ++d) {
                if (pd.mKvs[d].mNumCp != mCompr[p][d].mNumCp)
                    Fail("a patch disagrees with its knot vectors");
                n[d] = static_cast<std::int64_t>(pd.mKvs[d].mNumCp);
            }
            const MfPatchMap m = PatchMap(p, true);
            std::size_t t = 0;
            for (std::int64_t k = 0; k < n[2]; ++k)
                for (std::int64_t j = 0; j < n[1]; ++j)
                    for (std::int64_t i = 0; i < n[0]; ++i, ++t) {
                        const std::int64_t idx[3] = {i, j, k};
                        const auto g = static_cast<std::size_t>(m.At(idx, mDim));
                        const double* cp = pd.mCps.data() + t * (sd + 1);
                        rWeights[g] = cp[sd];
                        for (std::size_t c = 0; c < sd; ++c)
                            rXyz[g * sd + c] = cp[c] / cp[sd];
                    }
        }
        return;
    }
    if (!mHasNodes)
        Fail("no control points");
    const auto vdim = static_cast<std::size_t>(mSpace.mVDim);
    if (mNodes.size() != ndofs * vdim)
        Fail(std::to_string(mNodes.size()) + " node values for " + std::to_string(ndofs) +
             " control points of dimension " + std::to_string(vdim));
    rSpaceDim = mSpace.mVDim;
    rXyz.resize(ndofs * vdim);
    for (std::size_t g = 0; g < ndofs; ++g)
        for (std::size_t c = 0; c < vdim; ++c)
            rXyz[g * vdim + c] = mf_value(mNodes, mSpace, ndofs, g, c);
    rWeights = mWeights.empty() ? std::vector<double>(ndofs, 1.0) : mWeights;
}

// rTable (NumDofs x NumComp) at the reference points of a knot-span element:
// the rational (NURBS) interpolant, NumComp values per point into rOut.
void MfNurbs::Evaluate(const std::vector<double>& rWeights, const std::vector<double>& rTable,
                       std::size_t NumComp, const MfNurbsElement& rEl,
                       const std::vector<std::array<double, 3>>& rRefs,
                       std::vector<double>& rOut) const {
    const auto& kvs = mCompr[rEl.mPatch];
    const MfPatchMap m = PatchMap(rEl.mPatch, true);
    const std::size_t dim = kvs.size();
    rOut.assign(rRefs.size() * NumComp, 0.0);
    std::vector<double> basis[3];
    std::int64_t first[3] = {0, 0, 0};
    std::int64_t n[3] = {1, 1, 1};
    std::vector<double> num(NumComp);
    for (std::size_t t = 0; t < rRefs.size(); ++t) {
        for (std::size_t d = 0; d < dim; ++d) {
            const std::size_t s = kvs[d].mSpans[static_cast<std::size_t>(rEl.mSpan[d])];
            first[d] = static_cast<std::int64_t>(s) - kvs[d].mOrder;
            n[d] = kvs[d].mOrder + 1;
            kvs[d].Basis(s, rRefs[t][d], basis[d]);
        }
        std::fill(num.begin(), num.end(), 0.0);
        double den = 0.0;
        for (std::int64_t c = 0; c < n[2]; ++c)
            for (std::int64_t b = 0; b < n[1]; ++b)
                for (std::int64_t a = 0; a < n[0]; ++a) {
                    const std::int64_t ab[3] = {a, b, c};
                    double w = 1.0;
                    std::int64_t cp[3] = {0, 0, 0};
                    for (std::size_t d = 0; d < dim; ++d) {
                        w *= basis[d][static_cast<std::size_t>(ab[d])];
                        cp[d] = first[d] + ab[d];
                    }
                    const auto g = static_cast<std::size_t>(m.At(cp, static_cast<int>(dim)));
                    const double nw = w * rWeights[g];
                    for (std::size_t k = 0; k < NumComp; ++k)
                        num[k] += nw * rTable[g * NumComp + k];
                    den += nw;
                }
        for (std::size_t k = 0; k < NumComp; ++k)
            rOut[t * NumComp + k] = num[k] / den;
    }
}

// Skips v1.1's refinement and spacing sections after the knot vectors.
void mf_nurbs_skip_spacing(MfLexer& rLex) {
    if (!rLex.AtEnd() && rLex.Peek().mText == "refinements") {
        rLex.Next("refinements");
        rLex.Reals();
    }
    if (!rLex.AtEnd() && rLex.Peek().mText == "knotvector_refinements") {
        rLex.Next("knotvector_refinements");
        rLex.Reals();
    }
    if (rLex.AtEnd() || rLex.Peek().mText != "spacing")
        return;
    rLex.Next("spacing");
    const std::int64_t n = rLex.Int("a spacing count");
    for (std::int64_t k = 0; k < n; ++k) {
        rLex.Int("a knot vector");
        rLex.Int("a spacing type");
        const std::int64_t ni = rLex.Int("a parameter count");
        const std::int64_t nr = rLex.Int("a parameter count");
        for (std::int64_t i = 0; i < ni; ++i)
            rLex.Int("a parameter");
        for (std::int64_t i = 0; i < nr; ++i)
            rLex.Real("a parameter");
    }
}

MfNurbs mf_parse_nurbs(MfLexer& rLex, const std::string& rPath, bool V11) {
    MfNurbs n;
    n.mPath = rPath;
    const auto section = [&](const char* pName) {
        const MfToken& t = rLex.Next(pName);
        if (t.mText != pName)
            rLex.Fail(std::string("expected '") + pName + "', found '" + t.mText + "'", t.mLine);
    };
    section("dimension");
    const std::int64_t d = rLex.Int("a dimension");
    if (d < 1 || d > 3)
        rLex.Fail("dimension " + std::to_string(d) + " (1, 2 or 3)", rLex.Line());
    n.mDim = static_cast<int>(d);
    section("elements");
    n.mPatches = mf_read_elements(rLex, "patch");
    section("boundary");
    n.mBPatches = mf_read_elements(rLex, "boundary patch");
    for (const MfElement& el : n.mPatches)
        if (el.mGeom != mf_nurbs_geom(n.mDim))
            rLex.Fail("a patch of geometry " + std::to_string(el.mGeom) + " in a " +
                          std::to_string(n.mDim) + "-D mesh",
                      el.mLine);
    for (const MfElement& el : n.mBPatches)
        if (el.mGeom != mf_nurbs_bdr_geom(n.mDim))
            rLex.Fail("a boundary patch of geometry " + std::to_string(el.mGeom), el.mLine);
    section("edges");
    const std::int64_t ne = rLex.Int("an edge count");
    for (std::int64_t e = 0; e < ne; ++e) {
        std::array<std::int64_t, 3> row;
        row[0] = rLex.Int("a knot vector");
        row[1] = rLex.Int("a vertex");
        row[2] = rLex.Int("a vertex");
        n.mEdgeRows.push_back(row);
    }
    section("vertices");
    const std::int64_t nv = rLex.Int("a vertex count");
    if (nv < 0)
        rLex.Fail("negative vertex count", rLex.Line());
    n.mNumTopoVertices = static_cast<std::size_t>(nv);
    for (const auto* list : {&n.mPatches, &n.mBPatches})
        for (const MfElement& el : *list)
            for (std::int64_t v : el.mVertices)
                if (v >= nv)
                    rLex.Fail("vertex out of range (" + std::to_string(nv) + " vertices)",
                              el.mLine);
    for (const auto& r : n.mEdgeRows)
        if (r[0] < 0 || r[1] < 0 || r[2] < 0 || r[1] >= nv || r[2] >= nv)
            rLex.Fail("bad edge row", rLex.Line());
    if (n.mEdgeRows.empty() && n.mDim > 1)
        throw ReadError("MFEM NURBS mesh: " + rPath +
                        " has no edges section (edges MFEM derives itself are not supported)");
    const MfToken& t = rLex.Next("'knotvectors' or 'patches'");
    if (t.mText == "knotvectors") {
        const std::int64_t nk = rLex.Int("a count");
        for (std::int64_t k = 0; k < nk; ++k)
            n.mKvs.push_back(mf_read_knot(rLex));
        if (V11)
            mf_nurbs_skip_spacing(rLex);
    } else if (t.mText == "patches") {
        n.mPatchForm = true;
        for (std::size_t p = 0; p < n.mPatches.size(); ++p) {
            MfNurbsPatchData pd;
            section("knotvectors");
            const std::int64_t nk = rLex.Int("a count");
            for (std::int64_t k = 0; k < nk; ++k)
                pd.mKvs.push_back(mf_read_knot(rLex));
            if (static_cast<int>(pd.mKvs.size()) != n.mDim)
                rLex.Fail("a patch with " + std::to_string(pd.mKvs.size()) + " knot vectors",
                          t.mLine);
            section("dimension");
            const std::int64_t sd = rLex.Int("a space dimension");
            if (sd < 1 || sd > 3)
                rLex.Fail("space dimension " + std::to_string(sd), rLex.Line());
            pd.mSpaceDim = static_cast<int>(sd);
            const MfToken& cpt = rLex.Next("control points");
            if (cpt.mText != "controlpoints" && cpt.mText != "controlpoints_homogeneous" &&
                cpt.mText != "controlpoints_cartesian")
                rLex.Fail("expected control points, found '" + cpt.mText + "'", cpt.mLine);
            const bool cartesian = cpt.mText == "controlpoints_cartesian";
            std::size_t count = 1;
            for (const MfKnotVector& kv : pd.mKvs)
                count *= kv.mNumCp;
            const auto w = static_cast<std::size_t>(sd) + 1;
            pd.mCps.reserve(count * w);
            for (std::size_t k = 0; k < count * w; ++k)
                pd.mCps.push_back(rLex.Real("a coordinate"));
            if (cartesian)
                for (std::size_t k = 0; k < count; ++k)
                    for (std::size_t c = 0; c + 1 < w; ++c)
                        pd.mCps[k * w + c] *= pd.mCps[k * w + w - 1];
            n.mPatchData.push_back(std::move(pd));
        }
    } else {
        rLex.Fail("expected 'knotvectors' or 'patches', found '" + t.mText + "'", t.mLine);
    }
    n.Topology();
    n.Offsets();
    while (!rLex.AtEnd()) {
        const MfToken& s = rLex.Next("a section");
        if (s.mText == "mesh_elements" || s.mText == "periodic") {
            throw ReadError("MFEM NURBS mesh: the '" + s.mText + "' section of " + rPath +
                            " is not supported");
        } else if (s.mText == "weights") {
            if (n.mPatchForm)
                rLex.Fail("weights in a mesh whose patches carry their own", s.mLine);
            for (std::int64_t k = 0; k < n.NumDofs(); ++k)
                n.mWeights.push_back(rLex.Real("a weight"));
        } else if (s.mText == "unitweights" || s.mText == "autoweights") {
            n.mWeights.assign(static_cast<std::size_t>(n.NumDofs()), 1.0);
        } else if (s.mText == "FiniteElementSpace" || (s.mText == "MFEM" && !rLex.AtEnd() &&
                                                       rLex.Peek().mText == "FiniteElementSpace")) {
            const bool versioned = s.mText == "MFEM";
            if (versioned) {
                rLex.Next("FiniteElementSpace");
                const MfToken& v = rLex.Next("a version");
                if (v.mText != "v1.0")
                    rLex.Fail("FiniteElementSpace version '" + v.mText + "'", v.mLine);
            }
            n.mSpace = mf_read_space_body(rLex, s.mLine);
            if (n.mSpace.mKind != MfSpace::Nurbs)
                rLex.Fail("NURBS mesh nodes outside a NURBS space", s.mLine);
            if (versioned)
                mf_read_space_end(rLex);
            n.mNodes = rLex.Reals();
            n.mHasNodes = true;
        } else if (s.mText == "mfem_mesh_end") {
            break;
        } else {
            rLex.Fail("unexpected '" + s.mText + "'", s.mLine);
        }
    }
    return n;
}

MfFile mf_parse(const std::string& rPath) {
    const std::string what = "MFEM mesh";
    MfLexer lex(what, mf_read_text(rPath, "MFEM mesh"));
    const std::string& header = lex.Header();
    if (header == "MFEM NC mesh v1.0" || header == "MFEM NC mesh v1.1")
        return mf_parse_nc(lex, rPath, header == "MFEM NC mesh v1.1");
    if (header.rfind("MFEM NC mesh", 0) == 0)
        lex.Fail("non-conforming mesh version '" + header + "' is not supported", lex.HeaderLine());
    if (header == "MFEM NURBS mesh v1.0" || header == "MFEM NURBS mesh v1.1") {
        MfFile f;
        auto nurbs =
            std::make_shared<MfNurbs>(mf_parse_nurbs(lex, rPath, header == "MFEM NURBS mesh v1.1"));
        f.mDim = nurbs->mDim;
        f.mpNurbs = std::move(nurbs);
        return f;
    }
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
            // The legacy non-conforming layout: the leaf mesh, then how it
            // refines; read as its leaves.
            f.mNonConforming = true;
            const std::int64_t n = lex.Int("a count");
            for (std::int64_t k = 0; k < n; ++k) {
                const std::size_t line = lex.Line();
                while (!lex.AtEnd() && lex.Line() == line)
                    lex.Next("a value");
            }
        } else if (t.mText == "mfem_serial_mesh_end") {
            mf_read_groups(lex, f);
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
    // it was prefixed with an empty line and the header is "FiniteElementSpace"
    // (or the versioned form MFEM writes for NURBS spaces).
    MfGridData g;
    g.mName = rGf.mName;
    const bool versioned = lex.Header() == "MFEM FiniteElementSpace v1.0";
    if (lex.Header() != "FiniteElementSpace" && !versioned)
        throw ReadError("MFEM grid function: '" + rGf.mPath +
                        "' does not start with FiniteElementSpace (variable-order spaces are "
                        "not supported)");
    g.mSpace = mf_read_space_body(lex, 1);
    if (versioned)
        mf_read_space_end(lex);
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

// Cell regions from the element and boundary attributes and the attribute
// sets: `rCellAttr`/`rCellIsBoundary` per global cell.
void mf_add_regions(Mesh& rMesh, const MfFile& rF, const std::vector<std::int64_t>& rCellAttr,
                    const std::vector<bool>& rCellIsBoundary) {
    const int dim = rF.mDim;
    std::map<std::pair<bool, std::int64_t>, std::vector<std::int64_t>> by_attr;
    for (std::size_t g = 0; g < rCellAttr.size(); ++g)
        by_attr[{rCellIsBoundary[g], rCellAttr[g]}].push_back(static_cast<std::int64_t>(g));
    for (const auto& [key, ids] : by_attr) {
        const bool bdr = key.first;
        rMesh.AddRegion(Region((bdr ? "boundary_" : "attribute_") + std::to_string(key.second),
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
            rMesh.AddRegion(
                Region(s.mName, RegionKind::Cell, Bdr ? dim - 1 : dim, -1, mf_entries(ids)));
        }
    };
    add_sets(rF.mSets, false);
    add_sets(rF.mBdrSets, true);
}

// ---------------------------------------------------------------------------
// Order 3 and up: VTK Lagrange cells.
//
// A nodal H1 space's degrees of freedom are MFEM's: every vertex, then every
// edge's (q-1), every face's, and every element's interior ones, entities in
// order of first appearance (FiniteElementSpace::Construct). Rather than
// replaying MFEM's orientation tables, each element's dofs are placed at their
// reference positions -- an edge's from its lower global vertex (the edge's
// own direction), a face's in the frame of the vertex order its first element
// gave it, the interior in the element's own order (fem/fe/fe_h1.cpp) -- and
// matched by position to one canonical node set per (geometry, order, point
// family). One matrix per such set then interpolates to the equispaced nodes
// of the VTK Lagrange cell, which share their edge and face nodes.
// ---------------------------------------------------------------------------

using MfPos = std::array<double, 3>;

lagrange::Shape mf_shape(int Geom) {
    switch (Geom) {
        case 1:
            return lagrange::Shape::Line;
        case 2:
            return lagrange::Shape::Triangle;
        case 3:
            return lagrange::Shape::Quad;
        case 4:
            return lagrange::Shape::Tetra;
        case 5:
            return lagrange::Shape::Hexahedron;
        case 6:
            return lagrange::Shape::Wedge;
        default:
            throw ReadError("MFEM mesh: geometry " + std::to_string(Geom) +
                            " has no Lagrange cell");
    }
}

const char* mf_lagrange_type(int Geom) {
    static const char* const kTypes[] = {"vertex",
                                         "VTK_LAGRANGE_CURVE",
                                         "VTK_LAGRANGE_TRIANGLE",
                                         "VTK_LAGRANGE_QUADRILATERAL",
                                         "VTK_LAGRANGE_TETRAHEDRON",
                                         "VTK_LAGRANGE_HEXAHEDRON",
                                         "VTK_LAGRANGE_WEDGE"};
    return kTypes[Geom];
}

MfPos mf_ref_vertex(int Geom, std::size_t K) {
    static const MfPos kLine[2] = {{0, 0, 0}, {1, 0, 0}};
    static const MfPos kTri[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    static const MfPos kQuad[4] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    static const MfPos kTet[4] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    static const MfPos kHex[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                  {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    static const MfPos kWedge[6] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0},
                                    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    switch (Geom) {
        case 1:
            return kLine[K];
        case 2:
            return kTri[K];
        case 3:
            return kQuad[K];
        case 4:
            return kTet[K];
        case 5:
            return kHex[K];
        default:
            return kWedge[K];
    }
}

MfPos mf_mix(const std::vector<std::pair<double, MfPos>>& rTerms) {
    MfPos out = {0, 0, 0};
    for (const auto& [w, p] : rTerms)
        for (int d = 0; d < 3; ++d)
            out[static_cast<std::size_t>(d)] += w * p[static_cast<std::size_t>(d)];
    return out;
}

// The mesh's edges and (3-D) faces, numbered by first appearance, and each
// face's vertices in the order its first element lists them.
struct MfEntities {
    std::unordered_map<MfKey, std::size_t, MfKeyHash> mEdges;
    std::unordered_map<MfKey, std::size_t, MfKeyHash> mFaces;
    std::vector<std::vector<std::int64_t>> mFaceVertices;
};

MfEntities mf_entities(const std::vector<MfElement>& rElements, int Dim) {
    const auto& geoms = mf_geoms();
    MfEntities ent;
    for (const MfElement& el : rElements) {
        const MfGeom& g = geoms[static_cast<std::size_t>(el.mGeom)];
        if (Dim >= 2)
            for (const auto& e : g.mEdges)
                ent.mEdges.emplace(mf_key({el.mVertices[static_cast<std::size_t>(e[0])],
                                           el.mVertices[static_cast<std::size_t>(e[1])]}),
                                   ent.mEdges.size());
        if (Dim == 3)
            for (const auto& f : g.mFaces) {
                std::vector<std::int64_t> v;
                for (int k : f)
                    v.push_back(el.mVertices[static_cast<std::size_t>(k)]);
                if (ent.mFaces.emplace(mf_key(v), ent.mFaces.size()).second)
                    ent.mFaceVertices.push_back(std::move(v));
            }
    }
    return ent;
}

// The dof layout of one nodal space over a mesh.
struct MfDofs {
    int mOrder = 1;
    MfSpace::Points mPoints = MfSpace::Gll;
    std::vector<double> mCp;  // the 1-D points
    std::size_t mEdgeBase = 0;
    std::vector<std::size_t> mFaceOffset;      // per face
    std::vector<std::size_t> mInteriorOffset;  // per element
    std::size_t mSize = 0;
};

std::size_t mf_interior_count(int Geom, int Q) {
    const std::size_t q = static_cast<std::size_t>(Q);
    switch (Geom) {
        case 1:
            return q - 1;
        case 2:
            return (q - 1) * (q - 2) / 2;
        case 3:
            return (q - 1) * (q - 1);
        case 4:
            return q < 3 ? 0 : (q - 1) * (q - 2) * (q - 3) / 6;
        case 5:
            return (q - 1) * (q - 1) * (q - 1);
        case 6:
            return (q - 1) * (q - 2) / 2 * (q - 1);
        default:
            return 0;
    }
}

MfDofs mf_dofs(const MfFile& rF, const MfEntities& rEnt, int Q, MfSpace::Points Points) {
    const auto& geoms = mf_geoms();
    MfDofs d;
    d.mOrder = Q;
    d.mPoints = Points;
    d.mCp = Points == MfSpace::Gll ? lagrange::gll_points(Q) : lagrange::uniform_points(Q);
    const std::size_t q = static_cast<std::size_t>(Q);
    std::size_t next = rF.mNumVertices;
    d.mEdgeBase = next;
    next += rEnt.mEdges.size() * (q - 1);
    for (const auto& fv : rEnt.mFaceVertices) {
        d.mFaceOffset.push_back(next);
        next += fv.size() == 3 ? (q - 1) * (q - 2) / 2 : (q - 1) * (q - 1);
    }
    for (const MfElement& el : rF.mElements) {
        d.mInteriorOffset.push_back(next);
        if (geoms[static_cast<std::size_t>(el.mGeom)].mDim == rF.mDim)
            next += mf_interior_count(el.mGeom, Q);
    }
    d.mSize = next;
    return d;
}

// An element's dofs of a space and their reference positions. `pEnt` null:
// the reference element itself (vertices 0..n-1), for the canonical set.
void mf_element_dofs(const MfElement& rEl, int Dim, const MfEntities* pEnt, const MfDofs& rD,
                     std::size_t Element, std::vector<std::pair<std::size_t, MfPos>>& rOut) {
    rOut.clear();
    const auto& geoms = mf_geoms();
    const MfGeom& g = geoms[static_cast<std::size_t>(rEl.mGeom)];
    const int q = rD.mOrder;
    const auto& cp = rD.mCp;
    const auto cpq = [&](int i) { return cp[static_cast<std::size_t>(i)]; };
    const auto rv = [&](std::size_t k) { return mf_ref_vertex(rEl.mGeom, k); };
    for (std::size_t j = 0; j < g.mNumVertices; ++j)
        rOut.emplace_back(pEnt ? static_cast<std::size_t>(rEl.mVertices[j]) : j, rv(j));
    if (q < 2)
        return;
    std::size_t local = g.mNumVertices;  // canonical numbering when pEnt is null
    const auto local_of = [&](std::int64_t V) {
        return static_cast<std::size_t>(std::find(rEl.mVertices.begin(), rEl.mVertices.end(), V) -
                                        rEl.mVertices.begin());
    };
    if (Dim >= 2 && g.mDim >= 2) {
        for (const auto& e : g.mEdges) {
            const std::size_t a = static_cast<std::size_t>(e[0]),
                              b = static_cast<std::size_t>(e[1]);
            const bool a_low = rEl.mVertices[a] < rEl.mVertices[b];
            const std::size_t lo = a_low ? a : b, hi = a_low ? b : a;
            std::size_t base = 0;
            if (pEnt)
                base =
                    rD.mEdgeBase + pEnt->mEdges.at(mf_key({rEl.mVertices[a], rEl.mVertices[b]})) *
                                       static_cast<std::size_t>(q - 1);
            for (int k = 1; k < q; ++k) {
                const double t = cpq(k);
                rOut.emplace_back(pEnt ? base + static_cast<std::size_t>(k - 1) : local++,
                                  mf_mix({{1 - t, rv(lo)}, {t, rv(hi)}}));
            }
        }
    }
    if (Dim == 3 && g.mDim == 3) {
        for (const auto& f : g.mFaces) {
            std::vector<std::int64_t> v;
            for (int k : f)
                v.push_back(rEl.mVertices[static_cast<std::size_t>(k)]);
            std::size_t base = 0;
            std::vector<MfPos> r;
            if (pEnt) {
                const std::size_t idx = pEnt->mFaces.at(mf_key(v));
                base = rD.mFaceOffset[idx];
                for (std::int64_t gv : pEnt->mFaceVertices[idx])
                    r.push_back(rv(local_of(gv)));
            } else {
                for (int k : f)
                    r.push_back(rv(static_cast<std::size_t>(k)));
            }
            std::size_t o = 0;
            if (r.size() == 3) {
                for (int j = 1; j < q; ++j)
                    for (int i = 1; i + j < q; ++i) {
                        const double w = cpq(i) + cpq(j) + cpq(q - i - j);
                        const double x = cpq(i) / w, y = cpq(j) / w;
                        rOut.emplace_back(pEnt ? base + o++ : local++,
                                          mf_mix({{1 - x - y, r[0]}, {x, r[1]}, {y, r[2]}}));
                    }
            } else {
                for (int j = 1; j < q; ++j)
                    for (int i = 1; i < q; ++i) {
                        const double x = cpq(i), y = cpq(j);
                        rOut.emplace_back(pEnt ? base + o++ : local++,
                                          mf_mix({{(1 - x) * (1 - y), r[0]},
                                                  {x * (1 - y), r[1]},
                                                  {x * y, r[2]},
                                                  {(1 - x) * y, r[3]}}));
                    }
            }
        }
    }
    if (g.mDim != Dim)
        return;
    std::size_t o = pEnt ? rD.mInteriorOffset[Element] : 0;
    const auto add = [&](const MfPos& rP) { rOut.emplace_back(pEnt ? o++ : local++, rP); };
    switch (rEl.mGeom) {
        case 1:
            for (int i = 1; i < q; ++i)
                add({cpq(i), 0, 0});
            break;
        case 2:
        case 6: {
            std::vector<std::pair<double, double>> tri;
            for (int j = 1; j < q; ++j)
                for (int i = 1; i + j < q; ++i) {
                    const double w = cpq(i) + cpq(j) + cpq(q - i - j);
                    tri.emplace_back(cpq(i) / w, cpq(j) / w);
                }
            if (rEl.mGeom == 2) {
                for (const auto& [x, y] : tri)
                    add({x, y, 0});
            } else {
                for (int k = 1; k < q; ++k)
                    for (const auto& [x, y] : tri)
                        add({x, y, cpq(k)});
            }
            break;
        }
        case 3:
            for (int j = 1; j < q; ++j)
                for (int i = 1; i < q; ++i)
                    add({cpq(i), cpq(j), 0});
            break;
        case 4:
            for (int k = 1; k < q; ++k)
                for (int j = 1; j + k < q; ++j)
                    for (int i = 1; i + j + k < q; ++i) {
                        const double w = cpq(i) + cpq(j) + cpq(k) + cpq(q - i - j - k);
                        add({cpq(i) / w, cpq(j) / w, cpq(k) / w});
                    }
            break;
        case 5:
            if (rD.mPoints == MfSpace::Cubic) {
                // LagrangeHexFiniteElement(3): the eight interior nodes go
                // round each layer like the corners do.
                static const int kCubic[8][3] = {{1, 1, 1}, {2, 1, 1}, {2, 2, 1}, {1, 2, 1},
                                                 {1, 1, 2}, {2, 1, 2}, {2, 2, 2}, {1, 2, 2}};
                for (const auto& c : kCubic)
                    add({cpq(c[0]), cpq(c[1]), cpq(c[2])});
            } else {
                for (int k = 1; k < q; ++k)
                    for (int j = 1; j < q; ++j)
                        for (int i = 1; i < q; ++i)
                            add({cpq(i), cpq(j), cpq(k)});
            }
            break;
        default:
            break;
    }
}

// Quantised reference positions -> canonical node, for the matching.
struct MfPosIndex {
    std::map<std::array<std::int64_t, 3>, std::size_t> mMap;
    static std::array<std::int64_t, 3> Key(const MfPos& rP) {
        return {std::llround(rP[0] * 1e7), std::llround(rP[1] * 1e7), std::llround(rP[2] * 1e7)};
    }
    void Add(const MfPos& rP, std::size_t Index) { mMap.emplace(Key(rP), Index); }
    std::size_t Find(const MfPos& rP) const {
        const auto key = Key(rP);
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = mMap.find({key[0] + dx, key[1] + dy, key[2] + dz});
                    if (it != mMap.end())
                        return it->second;
                }
        return static_cast<std::size_t>(-1);
    }
};

// One (geometry, space) pair's canonical nodes and its matrix to the VTK
// Lagrange nodes of the output order.
struct MfInterp {
    MfPosIndex mIndex;
    std::size_t mNodes = 0;
    std::vector<double> mMatrix;  // VTK nodes x canonical nodes
};

// The canonical nodes of a (geometry, space) pair, in the reference element's
// own dof order, and their position index.
std::vector<MfPos> mf_canonical(int Geom, int Dim, const MfDofs& rD, MfPosIndex& rIndex) {
    MfElement ref;
    ref.mGeom = Geom;
    for (std::size_t k = 0; k < mf_geoms()[static_cast<std::size_t>(Geom)].mNumVertices; ++k)
        ref.mVertices.push_back(static_cast<std::int64_t>(k));
    std::vector<std::pair<std::size_t, MfPos>> nodes;
    mf_element_dofs(ref, Dim, nullptr, rD, 0, nodes);
    std::vector<MfPos> pos(nodes.size());
    for (const auto& [k, p] : nodes) {
        pos[k] = p;
        rIndex.Add(p, k);
    }
    return pos;
}

// The equispaced VTK Lagrange nodes of order `Order` on a geometry, as
// reference positions, in VTK order.
std::vector<MfPos> mf_vtk_positions(int Geom, int Order) {
    std::vector<MfPos> out;
    for (const auto& ijk : lagrange::vtk_lattice(mf_shape(Geom), Order))
        out.push_back({static_cast<double>(ijk[0]) / Order, static_cast<double>(ijk[1]) / Order,
                       static_cast<double>(ijk[2]) / Order});
    return out;
}

// The Lagrange order of a written cell by its node count: 1 for its corners,
// 2 for the complete quadratic types (VTK order-2 Lagrange node order), p for
// a VTK Lagrange cell; -1 for a serendipity cell.
int mf_cell_order(int Geom, std::size_t NumNodes) {
    if (Geom <= 0 || Geom > 6)
        return 1;
    for (int q = 1; q <= 64; ++q) {
        const std::size_t n = lagrange::num_nodes(mf_shape(Geom), q);
        if (n == NumNodes)
            return q;
        if (n > NumNodes)
            break;
    }
    return -1;
}

const MfInterp& mf_interp(std::map<std::pair<int, int>, MfInterp>& rCache, int Geom, int Dim,
                          const MfDofs& rD, int Order, int Slot) {
    const auto key = std::make_pair(Geom, Slot);
    auto it = rCache.find(key);
    if (it != rCache.end())
        return it->second;
    MfInterp in;
    const std::vector<MfPos> pos = mf_canonical(Geom, Dim, rD, in.mIndex);
    in.mNodes = pos.size();
    const lagrange::Shape shape = mf_shape(Geom);
    in.mMatrix =
        lagrange::interpolation_matrix(shape, rD.mOrder, pos, mf_vtk_positions(Geom, Order));
    return rCache.emplace(key, std::move(in)).first->second;
}

// A nodal field over one part: its space layout and a value accessor.
struct MfField {
    MfDofs mDofs;
    std::size_t mComponents = 1;
    std::function<double(std::size_t, std::size_t)> mValue;  // (dof, component)
};

// One part of the mesh to read: a whole serial file, or one rank of a parallel
// mesh, with its local vertices' output numbers.
struct MfPart {
    const MfFile* mpF = nullptr;
    std::vector<const MfGridData*> mGfs;  // the grid functions, in the caller's order
    bool mNodesH1 = false;                // coordinates from H1 nodes, else the vertices
    std::vector<double> mVertexXyz;       // local vertices, SpaceDim each
    std::vector<std::int64_t> mGlobal;    // local vertex -> output vertex
    // A part whose values it computes itself (a NURBS mesh): per element, every
    // field's values at the cell's lattice nodes (mEvalComponents each).
    std::function<void(std::size_t, std::vector<std::vector<double>>&)> mEval;
    std::vector<std::size_t> mEvalComponents;
};

// The output type of a cell: VTK Lagrange at order 3 and up, else the linear or
// complete quadratic type (their node order is VTK's order-1 and order-2
// Lagrange order).
std::string mf_cell_type(int Geom, int Order) {
    if (Geom == 0 || Order <= 1)
        return mf_geoms()[static_cast<std::size_t>(Geom)].mLinear;
    if (Order == 2)
        return mf_geoms()[static_cast<std::size_t>(Geom)].mQuadratic;
    return mf_lagrange_type(Geom);
}

// Reads the parts into one mesh of order-`Order` cells: every cell's nodes
// keyed by the output vertices they combine, so parts share their interface
// nodes; values from each part's own dof numbering; `partition:part` (the
// rank) on every cell when `Labels`.
Mesh mf_read_parts(const std::vector<MfPart>& rParts, std::size_t NumVertices, int Order,
                   bool Labels, const std::string& rPath) {
    const auto& geoms = mf_geoms();
    const MfFile& first = *rParts[0].mpF;
    const int dim = first.mDim;
    const std::size_t sdim = static_cast<std::size_t>(first.mSpaceDim);
    const std::size_t ngfs = rParts[0].mGfs.size();

    // Per part: its entities and fields (the coordinates, then each H1 grid
    // function); the L2 ones become cell data.
    std::vector<MfEntities> ents(rParts.size());
    std::vector<std::vector<MfField>> fields(rParts.size());
    std::vector<std::size_t> point_gfs;  // indices into mGfs
    for (std::size_t g = 0; g < ngfs; ++g)
        if (rParts[0].mGfs[g]->mSpace.mKind == MfSpace::H1 ||
            rParts[0].mGfs[g]->mSpace.mKind == MfSpace::Nurbs)
            point_gfs.push_back(g);
    for (std::size_t q = 0; q < rParts.size(); ++q) {
        const MfPart& part = rParts[q];
        const MfFile& f = *part.mpF;
        if (f.mDim != dim || static_cast<std::size_t>(f.mSpaceDim) != sdim)
            throw ReadError("MFEM mesh: the parts of " + rPath + " disagree on the dimension");
        for (const MfElement& el : f.mElements)
            if (el.mGeom == 0 || geoms[static_cast<std::size_t>(el.mGeom)].mDim != dim ||
                (el.mGeom == 7 && Order > 1))
                throw ReadError("MFEM mesh: element of geometry " + std::to_string(el.mGeom) +
                                " in a " + std::to_string(dim) + "-D mesh of order " +
                                std::to_string(Order) + " (line " + std::to_string(el.mLine) + ")");
        ents[q] = mf_entities(f.mElements, dim);
        if (part.mEval) {  // values come from the part itself
            for (std::size_t c : part.mEvalComponents) {
                MfField fld;
                fld.mComponents = c;
                fields[q].push_back(std::move(fld));
            }
            continue;
        }
        MfField xyz;
        xyz.mComponents = sdim;
        if (part.mNodesH1) {
            const MfSpace& s = f.mNodesSpace;
            xyz.mDofs = mf_dofs(f, ents[q], s.mOrder, s.mPoints);
            const std::size_t ndofs = f.mNodes.size() / static_cast<std::size_t>(s.mVDim);
            if (ndofs != xyz.mDofs.mSize)
                throw ReadError("MFEM mesh: the order-" + std::to_string(s.mOrder) + " nodes of " +
                                rPath + " have " + std::to_string(ndofs) +
                                " dofs; the mesh numbers " + std::to_string(xyz.mDofs.mSize));
            xyz.mValue = [&f, ndofs](std::size_t Dof, std::size_t C) {
                return mf_value(f.mNodes, f.mNodesSpace, ndofs, Dof, C);
            };
        } else {
            xyz.mDofs = mf_dofs(f, ents[q], 1, MfSpace::Gll);
            xyz.mValue = [&part, sdim](std::size_t Dof, std::size_t C) {
                return part.mVertexXyz[Dof * sdim + C];
            };
        }
        fields[q].push_back(std::move(xyz));
        for (std::size_t g : point_gfs) {
            const MfGridData& gd = *part.mGfs[g];
            const MfSpace& s = gd.mSpace;
            const std::size_t vdim = static_cast<std::size_t>(s.mVDim);
            if (gd.mValues.size() % vdim != 0)
                throw ReadError("MFEM grid function '" + gd.mName +
                                "': " + std::to_string(gd.mValues.size()) +
                                " values, not a multiple of VDim " + std::to_string(vdim));
            MfField field;
            field.mDofs = mf_dofs(f, ents[q], s.mOrder, s.mPoints);
            field.mComponents = vdim;
            const std::size_t ndofs = gd.mValues.size() / vdim;
            if (ndofs != field.mDofs.mSize)
                throw ReadError("MFEM grid function '" + gd.mName + "' has " +
                                std::to_string(ndofs) + " dofs; the mesh has " +
                                std::to_string(field.mDofs.mSize) + " at order " +
                                std::to_string(s.mOrder));
            field.mValue = [&gd, ndofs](std::size_t Dof, std::size_t C) {
                return mf_value(gd.mValues, gd.mSpace, ndofs, Dof, C);
            };
            fields[q].push_back(std::move(field));
        }
    }
    const std::size_t nfields = 1 + point_gfs.size();

    // Output points: the vertices first, then every other VTK node by the
    // output vertices and corner weights that place it.
    using NodeKey = std::vector<std::pair<std::int64_t, std::int64_t>>;
    std::map<NodeKey, std::size_t> node_of;
    const std::int64_t p3 = static_cast<std::int64_t>(Order) * Order * Order;
    for (std::size_t v = 0; v < NumVertices; ++v)
        node_of.emplace(NodeKey{{static_cast<std::int64_t>(v), p3}}, v);
    std::size_t npoints = NumVertices;
    const auto cell_nodes = [&](const MfPart& rPart, const MfElement& rEl,
                                std::vector<std::size_t>& rIds) {
        rIds.clear();
        const auto global = [&](std::int64_t V) {
            return rPart.mGlobal[static_cast<std::size_t>(V)];
        };
        if (rEl.mGeom == 0 || rEl.mGeom == 7) {  // a point; a (linear) pyramid
            for (std::int64_t v : rEl.mVertices)
                rIds.push_back(static_cast<std::size_t>(global(v)));
            return;
        }
        const lagrange::Shape shape = mf_shape(rEl.mGeom);
        for (const auto& ijk : lagrange::vtk_lattice(shape, Order)) {
            NodeKey key;
            for (const auto& [corner, w] : lagrange::lattice_weights(shape, Order, ijk))
                key.emplace_back(global(rEl.mVertices[static_cast<std::size_t>(corner)]), w);
            std::sort(key.begin(), key.end());
            const auto [it, fresh] = node_of.emplace(std::move(key), npoints);
            if (fresh)
                ++npoints;
            rIds.push_back(it->second);
        }
    };
    std::vector<std::vector<std::vector<std::size_t>>> element_nodes(rParts.size()),
        boundary_nodes(rParts.size());
    for (std::size_t q = 0; q < rParts.size(); ++q) {
        const MfFile& f = *rParts[q].mpF;
        element_nodes[q].resize(f.mElements.size());
        for (std::size_t e = 0; e < f.mElements.size(); ++e)
            cell_nodes(rParts[q], f.mElements[e], element_nodes[q][e]);
        boundary_nodes[q].resize(f.mBoundary.size());
        for (std::size_t b = 0; b < f.mBoundary.size(); ++b)
            cell_nodes(rParts[q], f.mBoundary[b], boundary_nodes[q][b]);
    }

    // Values at the nodes, element by element (a node shared by several takes
    // the first's; a conforming mesh gives the same value).
    std::vector<std::vector<double>> values(nfields);
    values[0].assign(npoints * sdim, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t k = 1; k < nfields; ++k)
        values[k].assign(npoints * fields[0][k].mComponents,
                         std::numeric_limits<double>::quiet_NaN());
    std::vector<bool> known(npoints, false);
    std::vector<std::pair<std::size_t, MfPos>> dofs;
    std::vector<std::size_t> perm;
    for (std::size_t q = 0; q < rParts.size(); ++q) {
        const MfFile& f = *rParts[q].mpF;
        std::vector<std::map<std::pair<int, int>, MfInterp>> caches(nfields);
        for (std::size_t e = 0; e < f.mElements.size(); ++e) {
            const MfElement& el = f.mElements[e];
            const auto& ids = element_nodes[q][e];
            bool any = false;
            for (std::size_t id : ids)
                any = any || !known[id];
            if (!any)
                continue;
            if (rParts[q].mEval) {
                std::vector<std::vector<double>> vals;
                rParts[q].mEval(e, vals);
                for (std::size_t k = 0; k < nfields; ++k) {
                    const std::size_t nc = fields[q][k].mComponents;
                    for (std::size_t t = 0; t < ids.size(); ++t)
                        if (!known[ids[t]])
                            for (std::size_t c = 0; c < nc; ++c)
                                values[k][ids[t] * nc + c] = vals[k][t * nc + c];
                }
                for (std::size_t id : ids)
                    known[id] = true;
                continue;
            }
            if (el.mGeom == 7) {  // a linear pyramid: its vertices
                for (std::size_t k = 0; k < nfields; ++k) {
                    const MfField& fld = fields[q][k];
                    for (std::size_t j = 0; j < ids.size(); ++j)
                        for (std::size_t c = 0; c < fld.mComponents; ++c)
                            values[k][ids[j] * fld.mComponents + c] =
                                fld.mValue(static_cast<std::size_t>(el.mVertices[j]), c);
                }
                for (std::size_t id : ids)
                    known[id] = true;
                continue;
            }
            for (std::size_t k = 0; k < nfields; ++k) {
                const MfField& fld = fields[q][k];
                const MfInterp& in = mf_interp(caches[k], el.mGeom, dim, fld.mDofs, Order, 0);
                mf_element_dofs(el, dim, &ents[q], fld.mDofs, e, dofs);
                perm.assign(in.mNodes, static_cast<std::size_t>(-1));
                for (const auto& [dof, pos] : dofs) {
                    const std::size_t c = in.mIndex.Find(pos);
                    if (c >= in.mNodes)
                        throw ReadError("MFEM mesh: a degree of freedom of element " +
                                        std::to_string(e) + " matches no node of its element");
                    perm[c] = dof;
                }
                const std::size_t nc = fld.mComponents;
                std::vector<double> u(in.mNodes * nc);
                for (std::size_t m = 0; m < in.mNodes; ++m)
                    for (std::size_t c = 0; c < nc; ++c)
                        u[m * nc + c] = fld.mValue(perm[m], c);
                for (std::size_t t = 0; t < ids.size(); ++t) {
                    if (known[ids[t]])
                        continue;
                    const double* row = in.mMatrix.data() + t * in.mNodes;
                    for (std::size_t c = 0; c < nc; ++c) {
                        double v = 0.0;
                        for (std::size_t m = 0; m < in.mNodes; ++m)
                            v += row[m] * u[m * nc + c];
                        values[k][ids[t] * nc + c] = v;
                    }
                }
            }
            for (std::size_t id : ids)
                known[id] = true;
        }
    }
    // Vertices no element holds keep their own coordinates; a boundary node on
    // no element face is placed from its corners.
    for (const MfPart& part : rParts)
        for (std::size_t v = 0; v < part.mGlobal.size() && !part.mEval; ++v) {
            const auto g = static_cast<std::size_t>(part.mGlobal[v]);
            if (known[g])
                continue;
            for (std::size_t c = 0; c < sdim; ++c)
                values[0][g * sdim + c] = part.mVertexXyz[v * sdim + c];
            known[g] = true;
        }
    std::size_t orphans = 0;
    for (std::size_t q = 0; q < rParts.size(); ++q) {
        const MfPart& part = rParts[q];
        for (std::size_t b = 0; b < part.mpF->mBoundary.size() && !part.mEval; ++b) {
            const MfElement& el = part.mpF->mBoundary[b];
            if (el.mGeom == 0 || el.mGeom == 7)
                continue;
            const lagrange::Shape shape = mf_shape(el.mGeom);
            const auto lattice = lagrange::vtk_lattice(shape, Order);
            for (std::size_t t = 0; t < lattice.size(); ++t) {
                const std::size_t id = boundary_nodes[q][b][t];
                if (known[id])
                    continue;
                ++orphans;
                for (const auto& [corner, w] :
                     lagrange::lattice_weights(shape, Order, lattice[t])) {
                    const auto v =
                        static_cast<std::size_t>(el.mVertices[static_cast<std::size_t>(corner)]);
                    for (std::size_t c = 0; c < sdim; ++c) {
                        double& x = values[0][id * sdim + c];
                        x = (std::isnan(x) ? 0.0 : x) + static_cast<double>(w) /
                                                            static_cast<double>(p3) *
                                                            part.mVertexXyz[v * sdim + c];
                    }
                }
                known[id] = true;
            }
        }
    }
    if (orphans)
        log::warn(
            "MFEM mesh: {} boundary node(s) lie on no element face; placed from their "
            "corners",
            orphans);

    Mesh mesh;
    NDArray points(DType::Float64, {npoints, sdim});
    std::copy(values[0].begin(), values[0].end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // Cells: elements then boundary elements (every part's), a block per type
    // in order of first appearance within each group.
    struct Cell {
        std::size_t mPart;
        std::size_t mSource;
        bool mBoundary;
    };
    std::vector<std::pair<std::string, std::vector<Cell>>> blocks;
    for (const bool boundary : {false, true}) {
        std::map<std::string, std::size_t> index;
        for (std::size_t q = 0; q < rParts.size(); ++q) {
            const auto& list = boundary ? rParts[q].mpF->mBoundary : rParts[q].mpF->mElements;
            for (std::size_t k = 0; k < list.size(); ++k) {
                const std::string t = mf_cell_type(list[k].mGeom, Order);
                auto [it, fresh] = index.emplace(t, blocks.size());
                if (fresh)
                    blocks.push_back({t, {}});
                blocks[it->second].second.push_back({q, k, boundary});
            }
        }
    }
    std::vector<NDArray> attr_blocks, part_blocks;
    std::vector<std::int64_t> cell_attr;
    std::vector<bool> cell_is_boundary;
    std::vector<std::vector<std::size_t>> element_cell(rParts.size());
    for (std::size_t q = 0; q < rParts.size(); ++q)
        element_cell[q].resize(rParts[q].mpF->mElements.size());
    std::vector<std::size_t> block_sizes;
    std::size_t global = 0;
    for (const auto& [type, members] : blocks) {
        const Cell& c0 = members[0];
        const std::size_t k =
            (c0.mBoundary ? boundary_nodes : element_nodes)[c0.mPart][c0.mSource].size();
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray attrs(DType::Int64, {members.size()});
        NDArray parts(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const Cell& cell = members[r];
            const MfFile& f = *rParts[cell.mPart].mpF;
            const MfElement& el =
                cell.mBoundary ? f.mBoundary[cell.mSource] : f.mElements[cell.mSource];
            const auto& ids =
                (cell.mBoundary ? boundary_nodes : element_nodes)[cell.mPart][cell.mSource];
            for (std::size_t j = 0; j < k; ++j)
                c[r * k + j] = static_cast<std::int64_t>(ids[j]);
            attrs.As<std::int64_t>()[r] = el.mAttribute;
            parts.As<std::int64_t>()[r] = f.mRank;
            cell_attr.push_back(el.mAttribute);
            cell_is_boundary.push_back(cell.mBoundary);
            if (!cell.mBoundary)
                element_cell[cell.mPart][cell.mSource] = global;
            ++global;
        }
        mesh.AddCellBlock(type, std::move(conn));
        attr_blocks.push_back(std::move(attrs));
        part_blocks.push_back(std::move(parts));
        block_sizes.push_back(members.size());
    }
    if (!attr_blocks.empty())
        mesh.AddCellData("mfem:attribute", std::move(attr_blocks));
    if (Labels && !part_blocks.empty())
        mesh.AddCellData("partition:part", std::move(part_blocks));

    for (std::size_t k = 1; k < nfields; ++k) {
        const std::size_t nc = fields[0][k].mComponents;
        NDArray data =
            nc == 1 ? NDArray(DType::Float64, {npoints}) : NDArray(DType::Float64, {npoints, nc});
        std::copy(values[k].begin(), values[k].end(), data.As<double>());
        mesh.AddPointData(rParts[0].mGfs[point_gfs[k - 1]]->mName, std::move(data));
    }
    // Element-wise grid functions: each part's values on its elements.
    for (std::size_t g = 0; g < ngfs; ++g) {
        if (rParts[0].mGfs[g]->mSpace.mKind == MfSpace::H1 ||
            rParts[0].mGfs[g]->mSpace.mKind == MfSpace::Nurbs)
            continue;
        const std::size_t vdim = static_cast<std::size_t>(rParts[0].mGfs[g]->mSpace.mVDim);
        std::vector<double> per_cell(global * vdim, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t q = 0; q < rParts.size(); ++q) {
            const MfGridData& gd = *rParts[q].mGfs[g];
            const std::size_t ne = rParts[q].mpF->mElements.size();
            const std::size_t ndofs = gd.mValues.size() / vdim;
            if (ndofs != ne)
                throw ReadError("MFEM grid function '" + gd.mName + "' has " +
                                std::to_string(ndofs) + " dofs for " + std::to_string(ne) +
                                " elements");
            for (std::size_t e = 0; e < ne; ++e)
                for (std::size_t c = 0; c < vdim; ++c)
                    per_cell[element_cell[q][e] * vdim + c] =
                        mf_value(gd.mValues, gd.mSpace, ndofs, e, c);
        }
        std::vector<NDArray> out;
        std::size_t start = 0;
        for (std::size_t n : block_sizes) {
            NDArray a =
                vdim == 1 ? NDArray(DType::Float64, {n}) : NDArray(DType::Float64, {n, vdim});
            std::copy(per_cell.begin() + static_cast<std::ptrdiff_t>(start * vdim),
                      per_cell.begin() + static_cast<std::ptrdiff_t>((start + n) * vdim),
                      a.As<double>());
            start += n;
            out.push_back(std::move(a));
        }
        mesh.AddCellData(rParts[0].mGfs[g]->mName, std::move(out));
    }
    mf_add_regions(mesh, first, cell_attr, cell_is_boundary);
    return mesh;
}

// The order-3-and-up read of one serial file.
Mesh mf_read_high_order(const MfFile& rF, const std::vector<MfGridData>& rGfs, bool NodesH1,
                        const std::vector<double>& rVertexXyz, int Order,
                        const std::string& rPath) {
    MfPart part;
    part.mpF = &rF;
    for (const MfGridData& g : rGfs)
        part.mGfs.push_back(&g);
    part.mNodesH1 = NodesH1;
    part.mVertexXyz = rVertexXyz;
    part.mGlobal.resize(rF.mNumVertices);
    for (std::size_t v = 0; v < rF.mNumVertices; ++v)
        part.mGlobal[v] = static_cast<std::int64_t>(v);
    return mf_read_parts({part}, rF.mNumVertices, Order, false, rPath);
}

// A NURBS mesh as its knot-span elements: VTK Lagrange cells (or linear and
// quadratic ones) of the highest knot-vector order, their nodes the rational
// patch geometry at the cell's lattice; NURBS grid functions on the mesh's own
// space the same way, element-wise ones as cell data.
Mesh mf_read_nurbs(const MfNurbs& rN, const std::vector<MfemGridFunction>& rGridFunctions,
                   const std::string& rPath) {
    const std::vector<MfNurbsElement> elements = rN.Elements();
    std::vector<double> weights, xyz;
    int sdim = 0;
    rN.ControlPoints(weights, xyz, sdim);
    int order = 1;
    for (const auto& row : rN.mCompr)
        for (const MfKnotVector& kv : row)
            order = std::max(order, kv.mOrder);
    std::vector<MfGridData> gfs;
    std::vector<std::vector<double>> tables;  // NURBS grid functions, dof-major
    std::vector<std::size_t> components = {static_cast<std::size_t>(sdim)};
    const auto ndofs = static_cast<std::size_t>(rN.NumDofs());
    for (const MfemGridFunction& gf : rGridFunctions) {
        MfGridData g = mf_parse_gf(gf);
        const auto vdim = static_cast<std::size_t>(g.mSpace.mVDim);
        if (g.mSpace.mKind == MfSpace::Nurbs && g.mValues.size() == ndofs * vdim) {
            std::vector<double> t(ndofs * vdim);
            for (std::size_t d = 0; d < ndofs; ++d)
                for (std::size_t c = 0; c < vdim; ++c)
                    t[d * vdim + c] = mf_value(g.mValues, g.mSpace, ndofs, d, c);
            tables.push_back(std::move(t));
            components.push_back(vdim);
            gfs.push_back(std::move(g));
        } else if ((g.mSpace.mKind == MfSpace::L2 || g.mSpace.mKind == MfSpace::L2T1) &&
                   g.mSpace.mOrder == 0) {
            gfs.push_back(std::move(g));
        } else {
            log::warn(
                "MFEM grid function '{}': the '{}' space is not the NURBS mesh's own nor "
                "element-wise; skipped",
                gf.mPath, g.mSpace.mCollection);
        }
    }
    MfFile f;
    f.mDim = rN.mDim;
    f.mSpaceDim = sdim;
    f.mNumVertices = static_cast<std::size_t>(rN.NumVertices());
    for (const MfNurbsElement& el : elements)
        f.mElements.push_back(el.mElement);
    f.mBoundary = rN.Boundary();

    const lagrange::Shape shape = mf_shape(mf_nurbs_geom(rN.mDim));
    std::vector<std::array<double, 3>> refs;
    for (const auto& ijk : lagrange::vtk_lattice(shape, order))
        refs.push_back({static_cast<double>(ijk[0]) / order, static_cast<double>(ijk[1]) / order,
                        static_cast<double>(ijk[2]) / order});
    std::vector<std::vector<double>> nurbs_tables;  // the coordinates, then each field
    nurbs_tables.push_back(std::move(xyz));
    for (auto& t : tables)
        nurbs_tables.push_back(std::move(t));

    MfPart part;
    part.mpF = &f;
    for (const MfGridData& g : gfs)
        part.mGfs.push_back(&g);
    part.mGlobal.resize(f.mNumVertices);
    for (std::size_t v = 0; v < f.mNumVertices; ++v)
        part.mGlobal[v] = static_cast<std::int64_t>(v);
    part.mEvalComponents = components;
    part.mEval = [&](std::size_t E, std::vector<std::vector<double>>& rOut) {
        rOut.resize(nurbs_tables.size());
        for (std::size_t k = 0; k < nurbs_tables.size(); ++k)
            rN.Evaluate(weights, nurbs_tables[k], components[k], elements[E], refs, rOut[k]);
    };
    return mf_read_parts({part}, f.mNumVertices, order, false, rPath);
}

// `<prefix>.NNNNNN`: the rank a parallel-mesh file name carries, or -1.
std::int64_t mf_rank_suffix(const std::string& rPath, std::string& rPrefix) {
    const std::size_t dot = rPath.rfind('.');
    if (dot == std::string::npos || rPath.size() - dot != 7)
        return -1;
    std::int64_t rank = 0;
    for (std::size_t k = dot + 1; k < rPath.size(); ++k) {
        if (rPath[k] < '0' || rPath[k] > '9')
            return -1;
        rank = 10 * rank + (rPath[k] - '0');
    }
    rPrefix = rPath.substr(0, dot + 1);
    return rank;
}

std::string mf_rank_path(const std::string& rPrefix, std::int64_t Rank) {
    std::string digits = std::to_string(Rank);
    if (digits.size() < 6)
        digits.insert(0, 6 - digits.size(), '0');
    return rPrefix + digits;
}

// Whether `rPath` is `<prefix>.NNNNNN` beside rank files 000000 and 000001:
// one rank of a parallel mesh as ParMesh::Save writes it.
bool mf_rank_siblings(const std::string& rPath) {
    std::string prefix;
    if (mf_rank_suffix(rPath, prefix) < 0)
        return false;
    std::error_code ec;
    return std::filesystem::exists(mf_rank_path(prefix, 0), ec) &&
           std::filesystem::exists(mf_rank_path(prefix, 1), ec);
}

// A parallel mesh: every rank file beside `rPath` (or the one piece asked for).
// With communication groups (ParMesh::ParPrint) the ranks' shared vertices are
// merged by group; without them (ParMesh::Save, the files GLVis reads) every
// rank is a serial mesh whose boundary also lists its faces on other ranks:
// boundary vertices at identical coordinates are merged, and a boundary face
// two ranks list is dropped.
Mesh mf_read_parallel(const std::string& rPath, MfFile First,
                      const std::vector<MfemGridFunction>& rGridFunctions,
                      const ReadOptions& rOptions) {
    namespace fs = std::filesystem;
    std::string prefix;
    std::vector<std::string> paths;
    if (mf_rank_suffix(rPath, prefix) >= 0) {
        std::error_code ec;
        for (std::int64_t r = 0; fs::exists(mf_rank_path(prefix, r), ec); ++r)
            paths.push_back(mf_rank_path(prefix, r));
    }
    if (paths.empty()) {
        log::warn(
            "MFEM mesh: '{}' is one rank of a parallel mesh, not named <prefix>.NNNNNN "
            "beside its siblings; reading that rank alone",
            rPath);
        paths = {rPath};
    }
    if (paths.size() == 1 && !First.mParallel)
        paths = {rPath};
    std::vector<std::size_t> selected;
    if (rOptions.mPieceSet) {
        if (rOptions.mPiece < 0 || static_cast<std::size_t>(rOptions.mPiece) >= paths.size())
            throw ReadError("MFEM mesh: piece " + std::to_string(rOptions.mPiece) +
                            " is out of range: " + rPath + " has " + std::to_string(paths.size()) +
                            " ranks");
        selected = {static_cast<std::size_t>(rOptions.mPiece)};
    } else {
        for (std::size_t r = 0; r < paths.size(); ++r)
            selected.push_back(r);
    }
    std::vector<MfFile> files;
    for (std::size_t r : selected) {
        files.push_back(paths[r] == rPath ? First : mf_parse(paths[r]));
        MfFile& f = files.back();
        if (f.mNonConforming)
            throw ReadError("MFEM mesh: " + paths[r] +
                            " is a non-conforming rank; parallel non-conforming meshes are not "
                            "read");
        if (f.mParallel && paths.size() > 1 && f.mRank != static_cast<std::int64_t>(r))
            throw ReadError("MFEM mesh: " + paths[r] + " holds rank " + std::to_string(f.mRank));
        if (!f.mParallel)
            f.mRank = static_cast<std::int64_t>(r);
    }

    bool groups = true;
    for (const MfFile& f : files)
        groups = groups && f.mParallel;
    std::vector<MfPart> parts(files.size());
    bool pyramid = false, corners = false;
    int order = 1;
    for (std::size_t q = 0; q < files.size(); ++q) {
        const MfFile& f = files[q];
        MfPart& part = parts[q];
        part.mpF = &f;
        const std::size_t sdim = static_cast<std::size_t>(f.mSpaceDim);
        if (f.mHasNodes) {
            const MfSpace& s = f.mNodesSpace;
            if (s.mKind == MfSpace::L2T1 || s.mKind == MfSpace::L2 || s.mKind == MfSpace::Other)
                throw ReadError("MFEM mesh: nodes in '" + s.mCollection +
                                "' are not supported in a parallel mesh");
            const std::size_t ndofs = f.mNodes.size() / static_cast<std::size_t>(s.mVDim);
            if (ndofs < f.mNumVertices)
                throw ReadError("MFEM mesh: the nodes of " + paths[selected[q]] + " have " +
                                std::to_string(ndofs) + " dofs for " +
                                std::to_string(f.mNumVertices) + " vertices");
            part.mVertexXyz.resize(f.mNumVertices * sdim);
            for (std::size_t v = 0; v < f.mNumVertices; ++v)
                for (std::size_t c = 0; c < sdim; ++c)
                    part.mVertexXyz[v * sdim + c] = mf_value(f.mNodes, s, ndofs, v, c);
            part.mNodesH1 = s.mKind == MfSpace::H1;
            corners = corners || s.mKind == MfSpace::H1Other;
            if (part.mNodesH1)
                order = std::max(order, s.mOrder);
        } else {
            part.mVertexXyz = f.mCoords;
        }
        for (const MfElement& el : f.mElements)
            pyramid = pyramid || el.mGeom == 7;
    }

    // Output vertices. With communication groups (ParPrint), a shared vertex by
    // (its group's ranks, its place in the group's list: every rank of a group
    // lists its vertices in the same order). Without them (ParMesh::Save), a
    // boundary vertex (any vertex of a 1-D mesh) is merged with the one at the
    // same position, to 1e-9 of the mesh's extent: ranks compute curved nodes
    // from their own elements, so they agree to rounding only. Every other
    // vertex is its own.
    std::size_t sdim_all = static_cast<std::size_t>(files[0].mSpaceDim);
    double extent = 0.0;
    {
        std::vector<double> lo(sdim_all, std::numeric_limits<double>::max()),
            hi(sdim_all, std::numeric_limits<double>::lowest());
        for (const MfPart& part : parts)
            for (std::size_t k = 0; k < part.mVertexXyz.size(); ++k) {
                lo[k % sdim_all] = std::min(lo[k % sdim_all], part.mVertexXyz[k]);
                hi[k % sdim_all] = std::max(hi[k % sdim_all], part.mVertexXyz[k]);
            }
        for (std::size_t c = 0; c < sdim_all; ++c)
            extent = std::max(extent, hi[c] - lo[c]);
    }
    const double tol = std::max(extent, 1.0) * 1e-9;
    std::map<std::vector<std::int64_t>, std::vector<std::int64_t>> buckets;
    std::vector<std::vector<double>> global_xyz;
    std::map<std::pair<std::vector<std::int64_t>, std::size_t>, std::int64_t> by_group;
    std::int64_t nglobal = 0;
    for (std::size_t q = 0; q < files.size(); ++q) {
        const MfFile& f = files[q];
        MfPart& part = parts[q];
        part.mGlobal.assign(f.mNumVertices, -1);
        const std::size_t sdim = static_cast<std::size_t>(f.mSpaceDim);
        std::vector<bool> candidate(f.mNumVertices, f.mDim == 1 && !groups);
        if (groups) {
            for (std::size_t g = 1; g < f.mGroupVertices.size(); ++g)
                for (std::size_t k = 0; k < f.mGroupVertices[g].size(); ++k) {
                    const std::int64_t v = f.mGroupVertices[g][k];
                    if (v < 0 || static_cast<std::size_t>(v) >= f.mNumVertices)
                        throw ReadError("MFEM mesh: shared vertex " + std::to_string(v) +
                                        " out of range in " + paths[selected[q]]);
                    const auto [it, fresh] =
                        by_group.emplace(std::make_pair(f.mGroups[g], k), nglobal);
                    if (fresh) {
                        ++nglobal;
                        global_xyz.emplace_back();
                    }
                    part.mGlobal[static_cast<std::size_t>(v)] = it->second;
                }
        } else {
            for (const MfElement& b : f.mBoundary)
                for (std::int64_t v : b.mVertices)
                    candidate[static_cast<std::size_t>(v)] = true;
        }
        for (std::size_t v = 0; v < f.mNumVertices; ++v) {
            if (!candidate[v])
                continue;
            const double* x = part.mVertexXyz.data() + v * sdim;
            std::vector<std::int64_t> cell(sdim);
            for (std::size_t c = 0; c < sdim; ++c)
                cell[c] = static_cast<std::int64_t>(std::floor(x[c] / tol));
            // Search the 3^sdim neighbouring buckets for a vertex within tol.
            std::int64_t match = -1;
            std::vector<std::int64_t> probe(sdim);
            const std::size_t combos = sdim == 1 ? 3 : (sdim == 2 ? 9 : 27);
            for (std::size_t m = 0; m < combos && match < 0; ++m) {
                std::size_t code = m;
                for (std::size_t c = 0; c < sdim; ++c) {
                    probe[c] = cell[c] + static_cast<std::int64_t>(code % 3) - 1;
                    code /= 3;
                }
                const auto it = buckets.find(probe);
                if (it == buckets.end())
                    continue;
                for (std::int64_t g : it->second) {
                    bool close = true;
                    for (std::size_t c = 0; c < sdim; ++c)
                        close = close &&
                                std::fabs(global_xyz[static_cast<std::size_t>(g)][c] - x[c]) <= tol;
                    if (close) {
                        match = g;
                        break;
                    }
                }
            }
            if (match < 0) {
                match = nglobal++;
                global_xyz.emplace_back(x, x + sdim);
                buckets[cell].push_back(match);
            }
            part.mGlobal[v] = match;
        }
        for (std::int64_t& g : part.mGlobal)
            if (g < 0) {
                g = nglobal++;
                global_xyz.emplace_back();
            }
    }
    // Without groups, a boundary face two ranks list is their interface.
    if (!groups && files.size() > 1) {
        std::map<std::vector<std::int64_t>, std::set<std::size_t>> owners;
        const auto key_of = [&](std::size_t Q, const MfElement& rB) {
            std::vector<std::int64_t> key;
            for (std::int64_t v : rB.mVertices)
                key.push_back(parts[Q].mGlobal[static_cast<std::size_t>(v)]);
            std::sort(key.begin(), key.end());
            return key;
        };
        for (std::size_t q = 0; q < files.size(); ++q)
            for (const MfElement& b : files[q].mBoundary)
                owners[key_of(q, b)].insert(q);
        std::size_t dropped = 0;
        for (std::size_t q = 0; q < files.size(); ++q) {
            std::vector<MfElement> kept;
            for (MfElement& b : files[q].mBoundary) {
                if (owners[key_of(q, b)].size() > 1)
                    ++dropped;
                else
                    kept.push_back(std::move(b));
            }
            files[q].mBoundary = std::move(kept);
        }
        if (dropped)
            log::debug("MFEM mesh: {} rank-interface boundary element(s) of {} dropped", dropped,
                       rPath);
    }
    if (corners)
        log::warn("MFEM mesh: nodes of {} that are not nodal H1 are read at the vertices only",
                  rPath);

    // Grid functions: each rank its own `<name>.NNNNNN`.
    std::vector<std::vector<MfGridData>> gdata(files.size());
    std::vector<std::string> gf_names;
    for (const MfemGridFunction& g : rGridFunctions) {
        std::string gprefix;
        if (mf_rank_suffix(g.mPath, gprefix) < 0)
            throw ReadError("MFEM grid function '" + g.mPath +
                            "': a parallel mesh's grid function is named by one of its rank "
                            "files, <name>.NNNNNN");
        std::vector<MfGridData> per_rank;
        for (std::size_t q = 0; q < files.size(); ++q)
            per_rank.push_back(mf_parse_gf({g.mName, mf_rank_path(gprefix, files[q].mRank)}));
        const MfSpace& s = per_rank[0].mSpace;
        const bool h1 = s.mKind == MfSpace::H1;
        const bool l2p0 = (s.mKind == MfSpace::L2 || s.mKind == MfSpace::L2T1) && s.mOrder == 0;
        if (!h1 && !l2p0) {
            log::warn("MFEM grid function '{}': the '{}' space is not supported; skipped", g.mPath,
                      s.mCollection);
            continue;
        }
        if (h1 && s.mOrder >= 2 && pyramid) {
            log::warn(
                "MFEM grid function '{}': an order-{} field on this mesh is not supported; "
                "skipped",
                g.mPath, s.mOrder);
            continue;
        }
        if (h1)
            order = std::max(order, s.mOrder);
        for (std::size_t q = 0; q < files.size(); ++q)
            gdata[q].push_back(std::move(per_rank[q]));
    }
    if (pyramid && order >= 2) {
        log::warn(
            "MFEM mesh: '{}' has pyramids, which meshio++ holds at order 1 only; reading "
            "the vertices",
            rPath);
        order = 1;
        for (MfPart& part : parts)
            part.mNodesH1 = false;
    }
    for (std::size_t q = 0; q < files.size(); ++q)
        for (const MfGridData& g : gdata[q])
            parts[q].mGfs.push_back(&g);
    return mf_read_parts(parts, static_cast<std::size_t>(nglobal), order, true, rPath);
}

}  // namespace

Mesh read_mfem(const std::string& rPath) {
    return read_mfem(rPath, {});
}

Mesh read_mfem(const std::string& rPath, const std::vector<MfemGridFunction>& rGridFunctions) {
    return read_mfem(rPath, rGridFunctions, ReadOptions{});
}

Mesh read_mfem(const std::string& rPath, const std::vector<MfemGridFunction>& rGridFunctions,
               const ReadOptions& rOptions) {
    MfFile f = mf_parse(rPath);
    if (f.mpNurbs) {
        if (rOptions.mPieceSet && rOptions.mPiece != 0)
            throw ReadError("MFEM mesh: " + rPath + " is not parallel; its only piece is 0");
        return mf_read_nurbs(*f.mpNurbs, rGridFunctions, rPath);
    }
    if (f.mParallel || mf_rank_siblings(rPath))
        return mf_read_parallel(rPath, std::move(f), rGridFunctions, rOptions);
    if (rOptions.mPieceSet && rOptions.mPiece != 0)
        throw ReadError("MFEM mesh: " + rPath + " is not parallel; its only piece is 0");
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
                "MFEM mesh: nodes in '{}' (order {}) are read at the vertices only: only "
                "nodal H1 spaces (Gauss-Lobatto or equispaced) are read as curved cells",
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
    if (f.mNonConforming && !rGridFunctions.empty())
        log::warn(
            "MFEM mesh: grid functions on the non-conforming mesh {} follow MFEM's "
            "space-filling-curve numbering of its leaves, which is not read; skipped",
            rPath);
    const std::vector<MfemGridFunction> no_gfs;
    const std::vector<MfemGridFunction>& gf_list = f.mNonConforming ? no_gfs : rGridFunctions;
    for (const MfemGridFunction& g : gf_list) {
        MfGridData data = mf_parse_gf(g);
        const MfSpace& s = data.mSpace;
        const bool h1 = s.mKind == MfSpace::H1;
        const bool l2p0 = (s.mKind == MfSpace::L2 || s.mKind == MfSpace::L2T1) && s.mOrder == 0;
        if (!h1 && !l2p0) {
            log::warn("MFEM grid function '{}': the '{}' space is not supported; skipped", g.mPath,
                      s.mCollection);
            continue;
        }
        if (h1 && s.mOrder >= 2 && (coords == Coords::Discontinuous || has_pyramid)) {
            log::warn(
                "MFEM grid function '{}': an order-{} field on this mesh is not supported; "
                "skipped",
                g.mPath, s.mOrder);
            continue;
        }
        gfs.push_back(std::move(data));
    }
    int order = coords == Coords::Discontinuous ? 1 : mesh_order;
    for (const MfGridData& g : gfs)
        if (g.mSpace.mKind == MfSpace::H1)
            order = std::max(order, g.mSpace.mOrder);
    if (order >= 2 && has_pyramid) {
        log::warn(
            "MFEM mesh: '{}' has pyramids, which meshio++ holds at order 1 only; reading "
            "the vertices",
            rPath);
        order = 1;
        if (coords == Coords::H1 && mesh_order >= 2)
            coords = Coords::Corners;
    }
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
    if (order >= 3)
        return mf_read_high_order(f, gfs, coords == Coords::H1, vxyz, order, rPath);
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
    mf_add_regions(mesh, f, cell_attr, cell_is_boundary);
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
                       type == "pyramid14" || type.rfind("VTK_LAGRANGE_", 0) == 0);
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
    // VTK Lagrange cells: the whole mesh is written with order-p H1 nodes, p
    // the highest order among its cells; linear and complete quadratic cells
    // are evaluated at those nodes too, serendipity cells from their corners.
    int high = 0;
    std::size_t serendipity = 0;
    for (const auto* list : {&elements, &boundary})
        for (const MfOutCell& c : *list)
            if (c.mType.rfind("VTK_LAGRANGE_", 0) == 0)
                high = std::max(high, mf_cell_order(c.mGeom, c.mNodes.size()));
    if (high > 0) {
        for (const auto* list : {&elements, &boundary})
            for (const MfOutCell& c : *list) {
                const int q = mf_cell_order(c.mGeom, c.mNodes.size());
                if (q < 0)
                    ++serendipity;
                high = std::max(high, q);
            }
        quadratic = false;
    }
    if ((quadratic || high > 1) && has_pyramid) {
        log::warn(
            "MFEM mesh writer: MFEM order-{} meshes hold no pyramids here; writing "
            "corners only",
            high > 1 ? high : 2);
        detail::provenance_note("high-order-dropped",
                                "a mesh with pyramids is written with linear MFEM cells");
        quadratic = false;
        high = 0;
    }
    if (high == 1)
        high = 0;  // order-1 Lagrange cells are their corners
    if (high > 0 && serendipity) {
        log::warn(
            "MFEM mesh writer: {} serendipity cell(s) in an order-{} mesh are placed "
            "from their corners",
            serendipity, high);
        detail::provenance_note("high-order-dropped",
                                "serendipity cells in a Lagrange mesh are written from corners");
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
    if (!quadratic && high == 0) {
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

    // Order-p H1 dofs (VTK Lagrange meshes): each element evaluated at its dofs'
    // positions, a dof shared by several taking the first's value.
    MfFile hf;
    MfEntities hent;
    MfDofs hdofs;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> high_map;  // per element
    std::map<std::pair<int, int>, std::vector<double>> matrices;             // (geom, cell order)
    std::vector<const std::vector<double>*> high_matrix;                     // per element
    if (high > 0) {
        hf.mDim = dim;
        hf.mNumVertices = nv;
        hf.mElements = numbered;
        hent = mf_entities(numbered, dim);
        hdofs = mf_dofs(hf, hent, high, MfSpace::Gll);
        std::map<int, std::pair<std::vector<MfPos>, MfPosIndex>> canon;
        std::vector<bool> taken(hdofs.mSize, false);
        std::vector<std::pair<std::size_t, MfPos>> dofs;
        high_map.resize(elements.size());
        for (std::size_t e = 0; e < elements.size(); ++e) {
            const MfOutCell& c = elements[e];
            auto cit = canon.find(c.mGeom);
            if (cit == canon.end()) {
                MfPosIndex index;
                std::vector<MfPos> pos = mf_canonical(c.mGeom, dim, hdofs, index);
                cit =
                    canon.emplace(c.mGeom, std::make_pair(std::move(pos), std::move(index))).first;
            }
            const int q = std::max(1, mf_cell_order(c.mGeom, c.mNodes.size()));
            const auto mkey = std::make_pair(c.mGeom, q);
            if (!matrices.count(mkey))
                matrices[mkey] = lagrange::interpolation_matrix(
                    mf_shape(c.mGeom), q, mf_vtk_positions(c.mGeom, q), cit->second.first);
            mf_element_dofs(numbered[e], dim, &hent, hdofs, e, dofs);
            for (const auto& [dof, pos] : dofs) {
                if (taken[dof])
                    continue;
                const std::size_t k = cit->second.second.Find(pos);
                if (k >= cit->second.first.size())
                    throw WriteError("MFEM mesh writer: a node of element " + std::to_string(e) +
                                     " matches no degree of freedom");
                taken[dof] = true;
                high_map[e].emplace_back(dof, k);
            }
        }
        high_matrix.resize(elements.size());
        for (std::size_t e = 0; e < elements.size(); ++e) {
            const MfOutCell& c = elements[e];
            const int q = std::max(1, mf_cell_order(c.mGeom, c.mNodes.size()));
            high_matrix[e] = &matrices.at(std::make_pair(c.mGeom, q));
        }
    }
    // The order-p dof values of one point array (row-major, dofs x Cols).
    auto high_values = [&](const NDArray& rData, std::size_t Cols) {
        std::vector<double> out(hdofs.mSize * Cols, 0.0);
        for (std::size_t e = 0; e < elements.size(); ++e) {
            const MfOutCell& c = elements[e];
            const std::size_t n = c.mNodes.size();
            const int q = std::max(1, mf_cell_order(c.mGeom, n));
            const std::size_t used = q == 1 ? corners_of(c) : n;
            for (const auto& [dof, k] : high_map[e]) {
                const double* row = high_matrix[e]->data() + k * used;
                for (std::size_t col = 0; col < Cols; ++col) {
                    double v = 0.0;
                    for (std::size_t m = 0; m < used; ++m)
                        v +=
                            row[m] * detail::read_double(
                                         rData, static_cast<std::size_t>(c.mNodes[m]) * Cols + col);
                    out[dof * Cols + col] = v;
                }
            }
        }
        return out;
    };

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
    if (high > 0) {
        out += "\nnodes\nFiniteElementSpace\nFiniteElementCollection: H1_" + std::to_string(dim) +
               "D_P" + std::to_string(high) + "\nVDim: " + std::to_string(pdim) +
               "\nOrdering: 1\n\n";
        const std::vector<double> xyz = high_values(pts, pdim);
        for (std::size_t d = 0; d < hdofs.mSize; ++d) {
            for (std::size_t c = 0; c < pdim; ++c) {
                if (c)
                    out += ' ';
                mf_append_real(out, xyz[d * pdim + c]);
            }
            out += '\n';
        }
    } else if (!quadratic) {
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
    const int order = high > 0 ? high : (quadratic ? 2 : 1);
    const std::size_t ndofs = high > 0 ? hdofs.mSize : (quadratic ? numbering.Size() : nv);
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
        const std::vector<double> hv = high > 0 ? high_values(a, cols) : std::vector<double>();
        for (std::size_t d = 0; d < ndofs; ++d) {
            for (std::size_t c = 0; c < cols; ++c) {
                if (c)
                    text += ' ';
                if (high > 0)
                    mf_append_real(text, hv[d * cols + c]);
                else if (quadratic)
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
