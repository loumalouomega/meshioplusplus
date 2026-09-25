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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// COMSOL element type <-> meshio type.
const std::unordered_map<std::string, std::string>& comsol_types() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vtx", "vertex"},     {"edg", "line"},          {"tri", "triangle"},
        {"quad", "quad"},      {"tet", "tetra"},         {"prism", "wedge"},
        {"pyr", "pyramid"},    {"hex", "hexahedron"},    {"edg2", "line3"},
        {"tri2", "triangle6"}, {"quad2", "quad9"},       {"tet2", "tetra10"},
        {"prism2", "wedge18"}, {"hex2", "hexahedron27"}, {"pyr2", "pyramid14"}};
    return m;
}

std::string meshio_to_comsol(const std::string& rT) {
    for (const auto& [c, m] : comsol_types())
        if (m == rT)
            return c;
    return "";
}

// ---------------------------------------------------------------------------
// Sources: the text and the binary serialisation hold the same sequence of
// integers, doubles and strings.
// ---------------------------------------------------------------------------

void comsol_append_utf8(std::string& rOut, std::uint32_t Cp) {
    if (Cp < 0x80) {
        rOut += static_cast<char>(Cp);
    } else if (Cp < 0x800) {
        rOut += static_cast<char>(0xC0 | (Cp >> 6));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    } else if (Cp < 0x10000) {
        rOut += static_cast<char>(0xE0 | (Cp >> 12));
        rOut += static_cast<char>(0x80 | ((Cp >> 6) & 0x3F));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    } else {
        rOut += static_cast<char>(0xF0 | (Cp >> 18));
        rOut += static_cast<char>(0x80 | ((Cp >> 12) & 0x3F));
        rOut += static_cast<char>(0x80 | ((Cp >> 6) & 0x3F));
        rOut += static_cast<char>(0x80 | (Cp & 0x3F));
    }
}

// The code points of a UTF-8 string (a stray byte counts as itself).
std::vector<std::uint32_t> comsol_code_points(const std::string& rS) {
    std::vector<std::uint32_t> out;
    for (std::size_t i = 0; i < rS.size();) {
        const unsigned char c = static_cast<unsigned char>(rS[i]);
        int extra = c >= 0xF0 ? 3 : c >= 0xE0 ? 2 : c >= 0xC0 ? 1 : 0;
        if (i + static_cast<std::size_t>(extra) >= rS.size() + (extra ? 0 : 1))
            extra = 0;
        std::uint32_t cp = extra == 0 ? c : (c & (0x3F >> extra));
        for (int k = 1; k <= extra; ++k)
            cp = (cp << 6) | (static_cast<unsigned char>(rS[i + k]) & 0x3F);
        out.push_back(cp);
        i += static_cast<std::size_t>(extra) + 1;
    }
    return out;
}

class ComsolSource {
public:
    virtual ~ComsolSource() = default;
    virtual std::int64_t Int() = 0;
    virtual double Real() = 0;
    virtual std::string String() = 0;
    virtual bool AtEnd() = 0;
    virtual std::size_t Mark() const = 0;
    virtual void Reset(std::size_t Pos) = 0;
    // Skip one value of a parameter record (a double).
    virtual void SkipValue() = 0;
    // Whether the next value reads as an integer (text: no '.', 'e').
    virtual bool NextIsInteger() = 0;
    // Whether the next value is a type-name string (text: a length, then a letter).
    virtual bool NextIsName() = 0;
    // The input's size in bytes: the bound on any count it declares.
    virtual std::size_t Size() const = 0;

    // A count read from the file, checked before it sizes anything: every
    // entry it counts takes at least one byte of the input.
    std::int64_t Count(const char* pWhat) {
        return static_cast<std::int64_t>(detail::checked_count(Int(), Size(), "COMSOL", pWhat));
    }
};

class ComsolText : public ComsolSource {
public:
    explicit ComsolText(std::string Text) : mText(std::move(Text)) {}
    std::size_t Size() const override { return mText.size(); }

    std::int64_t Int() override {
        const std::string t = Token();
        char* end = nullptr;
        const long long v = std::strtoll(t.c_str(), &end, 10);
        if (end == t.c_str() || *end != '\0')
            throw ReadError("mphtxt: expected an integer, found '" + t + "'");
        return v;
    }
    double Real() override {
        const std::string t = Token();
        const char* end = nullptr;
        const double v = detail::parse_double(t.c_str(), end);
        if (end == t.c_str() || *end != '\0')
            throw ReadError("mphtxt: expected a number, found '" + t + "'");
        return v;
    }
    // A length, one blank, then exactly that many characters (which may
    // include blanks and '#').
    std::string String() override {
        const std::int64_t n = Int();
        if (n < 0)
            throw ReadError("mphtxt: negative string length");
        if (n == 0)
            return "";
        if (mPos < mText.size() && (mText[mPos] == ' ' || mText[mPos] == '\t'))
            ++mPos;
        std::string out;
        for (std::int64_t k = 0; k < n; ++k) {
            if (mPos >= mText.size())
                throw ReadError("mphtxt: unexpected end of file in a string");
            const unsigned char c = static_cast<unsigned char>(mText[mPos]);
            std::size_t len = c >= 0xF0 ? 4 : c >= 0xE0 ? 3 : c >= 0xC0 ? 2 : 1;
            len = std::min(len, mText.size() - mPos);
            out.append(mText, mPos, len);
            mPos += len;
        }
        return out;
    }
    bool AtEnd() override {
        SkipBlank();
        return mPos >= mText.size();
    }
    std::size_t Mark() const override { return mPos; }
    void Reset(std::size_t Pos) override { mPos = Pos; }
    void SkipValue() override { Token(); }
    bool NextIsInteger() override {
        const std::size_t keep = mPos;
        if (AtEnd())
            return false;
        const std::string t = Token();
        mPos = keep;
        if (t.empty())
            return false;
        for (std::size_t k = (t[0] == '-' || t[0] == '+') ? 1 : 0; k < t.size(); ++k)
            if (t[k] < '0' || t[k] > '9')
                return false;
        return true;
    }
    bool NextIsName() override {
        const std::size_t keep = mPos;
        bool ok = NextIsInteger();
        if (ok) {
            Token();
            SkipBlank();
            ok = mPos < mText.size() && ((mText[mPos] >= 'a' && mText[mPos] <= 'z') ||
                                         (mText[mPos] >= 'A' && mText[mPos] <= 'Z'));
        }
        mPos = keep;
        return ok;
    }

private:
    void SkipBlank() {
        while (mPos < mText.size()) {
            const char c = mText[mPos];
            if (c == '#') {
                while (mPos < mText.size() && mText[mPos] != '\n')
                    ++mPos;
            } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                ++mPos;
            } else {
                break;
            }
        }
    }
    std::string Token() {
        SkipBlank();
        if (mPos >= mText.size())
            throw ReadError("mphtxt: unexpected end of file");
        const std::size_t b = mPos;
        while (mPos < mText.size()) {
            const char c = mText[mPos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '#')
                break;
            ++mPos;
        }
        return mText.substr(b, mPos - b);
    }

    std::string mText;
    std::size_t mPos = 0;
};

class ComsolBinary : public ComsolSource {
public:
    explicit ComsolBinary(std::string Bytes) : mBytes(std::move(Bytes)) {}
    std::size_t Size() const override { return mBytes.size(); }

    std::int64_t Int() override {
        Need(4);
        const unsigned char* p = reinterpret_cast<const unsigned char*>(mBytes.data() + mPos);
        const std::uint32_t u =
            static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
            (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
        mPos += 4;
        return static_cast<std::int32_t>(u);
    }
    double Real() override {
        Need(8);
        std::uint64_t u = 0;
        for (int k = 7; k >= 0; --k)
            u = (u << 8) | static_cast<unsigned char>(mBytes[mPos + static_cast<std::size_t>(k)]);
        mPos += 8;
        double v;
        std::memcpy(&v, &u, sizeof(v));
        return v;
    }
    std::string String() override {
        const std::int64_t n = Int();
        if (n < 0 || static_cast<std::uint64_t>(n) > (mBytes.size() - mPos) / 4)
            throw ReadError("mphbin: invalid string length " + std::to_string(n));
        std::string out;
        for (std::int64_t k = 0; k < n; ++k)
            comsol_append_utf8(out, static_cast<std::uint32_t>(Int()));
        return out;
    }
    bool AtEnd() override { return mPos >= mBytes.size(); }
    std::size_t Mark() const override { return mPos; }
    void Reset(std::size_t Pos) override { mPos = Pos; }
    void SkipValue() override {
        Need(8);
        mPos += 8;
    }
    bool NextIsInteger() override { return mPos + 4 <= mBytes.size(); }
    bool NextIsName() override {
        if (mPos + 8 > mBytes.size())
            return false;
        const std::size_t keep = mPos;
        const std::int64_t n = Int();
        const std::int64_t c = Int();
        mPos = keep;
        return n > 0 && n < 64 && ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'));
    }

private:
    void Need(std::size_t N) const {
        if (mPos + N > mBytes.size())
            throw ReadError("mphbin: unexpected end of file");
    }
    std::string mBytes;
    std::size_t mPos = 0;
};

// One Mesh object as read.
struct ComsolBlock {
    std::string mType;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;  // 0-based within its object, meshio order
    std::vector<std::int64_t> mGeom;
};
struct ComsolMesh {
    std::string mTag;
    std::size_t mSdim = 0;
    std::vector<double> mPoints;
    std::vector<ComsolBlock> mBlocks;
};
struct ComsolSelection {
    std::string mLabel;
    std::string mMeshTag;
    int mDim = 0;
    std::vector<std::int64_t> mEntities;
};

// Whether the next values end the file or open an object: 0 0 1 and a class name.
bool comsol_object_or_end(ComsolSource& rIn) {
    if (rIn.AtEnd())
        return true;
    for (int expected : {0, 0, 1})
        if (!rIn.NextIsInteger() || rIn.Int() != expected)
            return false;
    return rIn.NextIsName();
}

// Whether the rest of an element-type record of a version < 4 Mesh fits once
// the parameter values are skipped: `ne` entity indices, then up/down pairs,
// then the next type name, the next object or the end.
bool comsol_tail_fits(ComsolSource& rIn, std::int64_t Ne, bool Last) {
    const std::size_t keep = rIn.Mark();
    bool ok = false;
    try {
        const std::int64_t ngeom = rIn.NextIsInteger() ? rIn.Int() : -1;
        if (ngeom == Ne || ngeom == 0) {
            ok = true;
            for (std::int64_t k = 0; k < ngeom && ok; ++k)
                ok = rIn.NextIsInteger() && rIn.Int() >= 0;
            if (ok && rIn.NextIsInteger()) {
                const std::int64_t nud = rIn.Int();
                // Two integers per pair, each at least a byte of the input.
                ok = nud >= 0 && static_cast<std::uint64_t>(nud) <= rIn.Size() / 2;
                for (std::int64_t k = 0; ok && k < 2 * nud; ++k)
                    ok = rIn.NextIsInteger() && (rIn.Int(), true);
                if (ok)
                    ok = Last ? comsol_object_or_end(rIn) : rIn.NextIsName();
            } else {
                ok = false;
            }
        }
    } catch (const ReadError&) {
        ok = false;
    }
    rIn.Reset(keep);
    return ok;
}

ComsolMesh comsol_read_mesh(ComsolSource& rIn, const char* pFormat) {
    ComsolMesh m;
    const std::int64_t version = rIn.Int();
    const std::int64_t sdim = rIn.Int();
    const std::int64_t np = rIn.Int();
    const std::int64_t lowest = rIn.Int();
    if (sdim < 1 || sdim > 3 || np < 0 ||
        static_cast<std::uint64_t>(np) > rIn.Size() / static_cast<std::uint64_t>(sdim))
        throw ReadError(std::string(pFormat) + ": invalid Mesh header (sdim " +
                        std::to_string(sdim) + ", " + std::to_string(np) + " vertices)");
    m.mSdim = static_cast<std::size_t>(sdim);
    m.mPoints.resize(static_cast<std::size_t>(np * sdim));
    for (double& x : m.mPoints)
        x = rIn.Real();
    const std::int64_t ntypes = rIn.Int();
    for (std::int64_t t = 0; t < ntypes; ++t) {
        const std::string ctype = rIn.String();
        auto it = comsol_types().find(ctype);
        if (it == comsol_types().end())
            throw ReadError(std::string(pFormat) + ": unknown element type '" + ctype + "'");
        ComsolBlock b;
        b.mType = it->second;
        const std::int64_t nep = rIn.Int();
        const std::int64_t ne = rIn.Int();
        const int expected = cell_type_num_nodes(cell_type_from_name(b.mType));
        if (nep != expected || ne < 0 ||
            static_cast<std::uint64_t>(ne) > rIn.Size() / static_cast<std::uint64_t>(nep))
            throw ReadError(std::string(pFormat) + ": '" + ctype + "' elements with " +
                            std::to_string(nep) + " nodes");
        b.mNodes = static_cast<std::size_t>(nep);
        std::vector<std::int64_t> raw(static_cast<std::size_t>(ne * nep));
        for (std::int64_t& v : raw) {
            v = rIn.Int() - lowest;
            if (v < 0 || v >= np)
                throw ReadError(std::string(pFormat) + ": '" + ctype + "' element names vertex " +
                                std::to_string(v + lowest));
        }
        const detail::NodeOrder* order = detail::node_order("mphtxt", b.mType);
        b.mConn.resize(raw.size());
        for (std::size_t r = 0; r < static_cast<std::size_t>(ne); ++r)
            for (std::size_t j = 0; j < b.mNodes; ++j)
                b.mConn[r * b.mNodes + j] =
                    raw[r * b.mNodes + (order ? static_cast<std::size_t>(order->mToMeshio[j]) : j)];
        if (version < 4) {
            // Parameter records of npp values each, one value per parameter
            // dimension: find how many by where the rest of the record fits.
            const std::int64_t npp = rIn.Int();
            const std::int64_t npar = rIn.Int();
            // npp * npar * 3 values follow at most, each a byte or more.
            if (npar < 0 || npp < 0 ||
                (npp > 0 && static_cast<std::uint64_t>(npar) >
                                rIn.Size() / 3 / static_cast<std::uint64_t>(npp)))
                throw ReadError(std::string(pFormat) + ": implausible parameter count");
            const std::size_t start = rIn.Mark();
            bool found = false;
            for (int k = 1; k <= 3 && !found; ++k) {
                rIn.Reset(start);
                try {
                    for (std::int64_t v = 0; v < npp * npar * k; ++v)
                        rIn.SkipValue();
                    found = comsol_tail_fits(rIn, ne, t + 1 == ntypes);
                } catch (const ReadError&) {
                    found = false;
                }
            }
            if (!found)
                throw ReadError(std::string(pFormat) + ": cannot delimit the parameters of the '" +
                                ctype + "' elements");
        }
        // No entity indices at all is allowed: every element then gets
        // COMSOL's default, domain 1 or entity 0 below the space dimension.
        const std::int64_t ngeom = rIn.Int();
        if (ngeom != ne && ngeom != 0)
            throw ReadError(std::string(pFormat) + ": " + std::to_string(ngeom) +
                            " entity indices for " + std::to_string(ne) + " '" + ctype +
                            "' elements");
        const int dim = cell_type_dimension(cell_type_from_name(b.mType));
        b.mGeom.assign(static_cast<std::size_t>(ne), dim == sdim ? 1 : 0);
        for (std::int64_t k = 0; k < ngeom; ++k)
            b.mGeom[static_cast<std::size_t>(k)] = rIn.Int();
        if (version < 4) {
            const std::int64_t nud = rIn.Int();
            for (std::int64_t k = 0; k < 2 * nud; ++k)
                rIn.Int();
        }
        m.mBlocks.push_back(std::move(b));
    }
    return m;
}

Mesh comsol_read(ComsolSource& rIn, const char* pFormat) {
    const std::int64_t major = rIn.Int();
    const std::int64_t minor = rIn.Int();
    if (major != 0 || minor != 1)
        throw ReadError(std::string(pFormat) + ": unsupported file version " +
                        std::to_string(major) + "." + std::to_string(minor));
    std::vector<std::string> tags(static_cast<std::size_t>(rIn.Count("tag")));
    for (std::string& t : tags)
        t = rIn.String();
    const std::int64_t ntypes = rIn.Int();
    for (std::int64_t k = 0; k < ntypes; ++k)
        rIn.String();

    std::vector<ComsolMesh> meshes;
    std::vector<ComsolSelection> selections;
    for (std::int64_t obj = 0; obj < ntypes; ++obj) {
        if (rIn.AtEnd())
            break;
        rIn.Int();
        rIn.Int();
        rIn.Int();
        const std::string cls = rIn.String();
        const std::string tag = static_cast<std::size_t>(obj) < tags.size()
                                    ? tags[static_cast<std::size_t>(obj)]
                                    : std::string();
        if (cls == "Mesh") {
            meshes.push_back(comsol_read_mesh(rIn, pFormat));
            meshes.back().mTag = tag;
        } else if (cls == "Selection") {
            ComsolSelection s;
            rIn.Int();  // version
            s.mLabel = rIn.String();
            s.mMeshTag = rIn.String();
            s.mDim = static_cast<int>(rIn.Int());
            s.mEntities.resize(static_cast<std::size_t>(rIn.Count("entity")));
            for (std::int64_t& e : s.mEntities)
                e = rIn.Int();
            selections.push_back(std::move(s));
        } else {
            log::warn(
                "{}: object {} is a '{}', which meshio++ does not read; the objects "
                "after it are skipped too",
                pFormat, obj, cls);
            break;
        }
    }
    if (meshes.empty())
        throw ReadError(std::string(pFormat) + ": the file holds no Mesh object");

    // All Mesh objects in one mesh, their points one after the other.
    const std::size_t sdim = meshes.front().mSdim;
    std::size_t npoints = 0;
    for (const ComsolMesh& m : meshes) {
        if (m.mSdim != sdim)
            throw ReadError(std::string(pFormat) + ": Mesh objects of different dimensions");
        npoints += m.mPoints.size() / sdim;
    }
    Mesh mesh;
    NDArray pts(DType::Float64, {npoints, sdim});
    double* pp = pts.As<double>();
    for (const ComsolMesh& m : meshes)
        pp = std::copy(m.mPoints.begin(), m.mPoints.end(), pp);
    mesh.AssignPoints(std::move(pts));

    std::vector<NDArray> geom;
    std::vector<std::size_t> first_block;  // per object
    std::vector<int> block_dims;
    std::size_t offset = 0;
    for (const ComsolMesh& m : meshes) {
        first_block.push_back(geom.size());
        for (const ComsolBlock& b : m.mBlocks) {
            const std::size_t ne = b.mGeom.size();
            NDArray conn(DType::Int64, {ne, b.mNodes});
            std::int64_t* c = conn.As<std::int64_t>();
            for (std::size_t k = 0; k < b.mConn.size(); ++k)
                c[k] = b.mConn[k] + static_cast<std::int64_t>(offset);
            mesh.AddCellBlock(b.mType, std::move(conn));
            NDArray g(DType::Int64, {ne});
            std::copy(b.mGeom.begin(), b.mGeom.end(), g.As<std::int64_t>());
            geom.push_back(std::move(g));
            block_dims.push_back(cell_type_dimension(cell_type_from_name(b.mType)));
        }
        offset += m.mPoints.size() / sdim;
    }
    if (!geom.empty())
        mesh.AddCellData("mphtxt:geom", std::vector<NDArray>(geom.begin(), geom.end()));

    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    auto add = [&](const std::string& rName, int Dim, const std::vector<std::int64_t>& rIds) {
        NDArray e(DType::Int64, {rIds.size()});
        std::copy(rIds.begin(), rIds.end(), e.As<std::int64_t>());
        mesh.AddRegion(Region(rName, RegionKind::Cell, Dim, -1, std::move(e)));
    };
    std::set<std::string> taken;
    auto unique_name = [&](const std::string& rName) {
        std::string name = rName;
        for (int k = 2; taken.count(name); ++k)
            name = rName + " (" + std::to_string(k) + ")";
        if (name != rName)
            log::warn("{}: a second region named '{}' is read as '{}'", pFormat, rName, name);
        taken.insert(name);
        return name;
    };
    // With several Mesh objects, each one's cells are a region named by its tag.
    if (meshes.size() > 1)
        for (std::size_t o = 0; o < meshes.size(); ++o) {
            std::vector<std::int64_t> ids;
            const std::size_t end = first_block[o] + meshes[o].mBlocks.size();
            for (std::size_t b = first_block[o]; b < end; ++b)
                for (std::int64_t r = 0; r < bases[b + 1] - bases[b]; ++r)
                    ids.push_back(bases[b] + r);
            add(unique_name(meshes[o].mTag), -1, ids);
        }
    // A Selection: the elements of its dimension whose entity it lists.
    for (const ComsolSelection& s : selections) {
        std::size_t o = 0;
        while (o < meshes.size() && meshes[o].mTag != s.mMeshTag)
            ++o;
        if (o == meshes.size()) {
            log::warn("{}: selection '{}' refers to '{}', which is not a Mesh here; skipped",
                      pFormat, s.mLabel, s.mMeshTag);
            continue;
        }
        const std::set<std::int64_t> wanted(s.mEntities.begin(), s.mEntities.end());
        std::vector<std::int64_t> ids;
        const std::size_t end = first_block[o] + meshes[o].mBlocks.size();
        for (std::size_t b = first_block[o]; b < end; ++b) {
            if (block_dims[b] != s.mDim)
                continue;
            const std::int64_t* g = geom[b].As<std::int64_t>();
            for (std::int64_t r = 0; r < bases[b + 1] - bases[b]; ++r)
                if (wanted.count(g[r]))
                    ids.push_back(bases[b] + r);
        }
        add(unique_name(s.mLabel), s.mDim, ids);
    }
    return mesh;
}

std::string comsol_slurp(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// ---------------------------------------------------------------------------
// Writing: one serialiser over two sinks.
// ---------------------------------------------------------------------------

class ComsolSink {
public:
    virtual ~ComsolSink() = default;
    virtual void Int(std::int64_t V, const char* pComment = nullptr) = 0;
    virtual void Ints(const std::int64_t* pV, std::size_t N) = 0;  // one record
    virtual void Real(const double* pV, std::size_t N) = 0;        // one record
    virtual void String(const std::string& rS, const char* pComment = nullptr) = 0;
    virtual void Comment(const std::string& rText) = 0;
    virtual void Blank() = 0;
};

class ComsolTextSink : public ComsolSink {
public:
    explicit ComsolTextSink(std::string& rOut) : mOut(rOut) {}
    void Int(std::int64_t V, const char* pComment) override {
        mOut += std::to_string(V);
        End(pComment);
    }
    void Ints(const std::int64_t* pV, std::size_t N) override {
        for (std::size_t k = 0; k < N; ++k) {
            mOut += std::to_string(pV[k]);
            mOut += k + 1 == N ? '\n' : ' ';
        }
    }
    void Real(const double* pV, std::size_t N) override {
        char buf[40];
        for (std::size_t k = 0; k < N; ++k) {
            detail::snprintf_c(buf, sizeof(buf), "%.17g", pV[k]);
            mOut += buf;
            mOut += k + 1 == N ? '\n' : ' ';
        }
    }
    void String(const std::string& rS, const char* pComment) override {
        mOut += std::to_string(comsol_code_points(rS).size());
        mOut += ' ';
        mOut += rS;
        End(pComment);
    }
    void Comment(const std::string& rText) override { mOut += "# " + rText + "\n"; }
    void Blank() override { mOut += "\n"; }

private:
    void End(const char* pComment) {
        if (pComment) {
            mOut += " # ";
            mOut += pComment;
        }
        mOut += '\n';
    }
    std::string& mOut;
};

class ComsolBinarySink : public ComsolSink {
public:
    explicit ComsolBinarySink(std::string& rOut) : mOut(rOut) {}
    void Int(std::int64_t V, const char*) override { Put32(V); }
    void Ints(const std::int64_t* pV, std::size_t N) override {
        for (std::size_t k = 0; k < N; ++k)
            Put32(pV[k]);
    }
    void Real(const double* pV, std::size_t N) override {
        for (std::size_t k = 0; k < N; ++k) {
            std::uint64_t u;
            std::memcpy(&u, &pV[k], sizeof(u));
            for (int b = 0; b < 8; ++b)
                mOut += static_cast<char>((u >> (8 * b)) & 0xFF);
        }
    }
    void String(const std::string& rS, const char*) override {
        const std::vector<std::uint32_t> cps = comsol_code_points(rS);
        Put32(static_cast<std::int64_t>(cps.size()));
        for (std::uint32_t cp : cps)
            Put32(static_cast<std::int64_t>(cp));
    }
    void Comment(const std::string&) override {}
    void Blank() override {}

private:
    void Put32(std::int64_t V) {
        if (V < INT32_MIN || V > INT32_MAX)
            throw WriteError("mphbin: " + std::to_string(V) + " does not fit a 32-bit integer");
        const std::uint32_t u = static_cast<std::uint32_t>(static_cast<std::int32_t>(V));
        for (int b = 0; b < 4; ++b)
            mOut += static_cast<char>((u >> (8 * b)) & 0xFF);
    }
    std::string& mOut;
};

struct ComsolSelectionOut {
    std::string mLabel;
    int mDim = 0;
    std::vector<std::int64_t> mEntities;
};

void comsol_write(ComsolSink& rOut, const Mesh& rMesh, const char* pFormat) {
    const std::size_t sdim = rMesh.PointDim();
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<std::string> ctypes;
    std::vector<int> dims;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string c = cb.IsRagged() ? std::string() : meshio_to_comsol(cb.Type());
        if (c.empty())
            throw WriteError(std::string(pFormat) + ": unsupported cell type " + cb.Type());
        ctypes.push_back(c);
        dims.push_back(cell_type_dimension(cell_type_from_name(cb.Type())));
    }
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<int> cell_dim(ncells);
    for (std::size_t b = 0; b < nblocks; ++b)
        for (std::int64_t c = bases[b]; c < bases[b + 1]; ++c)
            cell_dim[static_cast<std::size_t>(c)] = dims[b];

    // Cell regions, validated.
    std::vector<const Region*> regions;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind != RegionKind::Cell) {
            log::warn(
                "{}: {} region '{}' dropped; COMSOL selections are written for cell "
                "regions only",
                pFormat, reg.mKind == RegionKind::Point ? "point" : "side", reg.mName);
            continue;
        }
        const std::int64_t* e = reg.Entries();
        for (std::size_t j = 0; j < reg.NumEntries(); ++j)
            if (e[j] < 0 || static_cast<std::size_t>(e[j]) >= ncells)
                throw WriteError(std::string(pFormat) + ": cell region '" + reg.mName +
                                 "' names cell " + std::to_string(e[j]) + " of " +
                                 std::to_string(ncells));
        regions.push_back(&reg);
    }

    // Geometric entity of every cell: mphtxt:geom when the mesh has it; else
    // per dimension, the pairwise-disjoint single-dimension cell regions in
    // order, then one entity for the cells in none. Domains (dimension sdim)
    // count from 1, lower dimensions from 0.
    std::vector<std::int64_t> entity(ncells, 0);
    if (rMesh.HasCellData("mphtxt:geom")) {
        for (std::size_t b = 0; b < nblocks; ++b) {
            const NDArray& g = rMesh.CellData("mphtxt:geom", b);
            for (std::int64_t c = bases[b]; c < bases[b + 1]; ++c)
                entity[static_cast<std::size_t>(c)] =
                    detail::read_int(g, static_cast<std::size_t>(c - bases[b]));
        }
    } else {
        std::vector<char> assigned(ncells, 0);
        std::map<int, std::int64_t> next;
        auto base = [&](int d) { return d == static_cast<int>(sdim) ? 1 : 0; };
        for (const Region* reg : regions) {
            const std::int64_t* e = reg->Entries();
            if (reg->NumEntries() == 0)
                continue;
            const int d = cell_dim[static_cast<std::size_t>(e[0])];
            bool ok = true;
            for (std::size_t j = 0; j < reg->NumEntries() && ok; ++j)
                ok = cell_dim[static_cast<std::size_t>(e[j])] == d &&
                     !assigned[static_cast<std::size_t>(e[j])];
            if (!ok)
                continue;
            auto it = next.try_emplace(d, base(d)).first;
            for (std::size_t j = 0; j < reg->NumEntries(); ++j) {
                entity[static_cast<std::size_t>(e[j])] = it->second;
                assigned[static_cast<std::size_t>(e[j])] = 1;
            }
            ++it->second;
        }
        for (std::size_t c = 0; c < ncells; ++c)
            if (!assigned[c]) {
                auto it = next.find(cell_dim[c]);
                entity[c] = it == next.end() ? base(cell_dim[c]) : it->second;
            }
    }

    // A cell region becomes a Selection when it is exactly a union of whole
    // entities of one dimension.
    std::vector<ComsolSelectionOut> selections;
    for (const Region* reg : regions) {
        const std::int64_t* e = reg->Entries();
        ComsolSelectionOut s;
        s.mLabel = reg->mName;
        if (reg->NumEntries() == 0) {
            s.mDim = reg->mDim >= 0 ? reg->mDim : static_cast<int>(sdim);
            selections.push_back(std::move(s));
            continue;
        }
        const int d = cell_dim[static_cast<std::size_t>(e[0])];
        std::set<std::int64_t> ents;
        bool ok = true;
        for (std::size_t j = 0; j < reg->NumEntries() && ok; ++j) {
            ok = cell_dim[static_cast<std::size_t>(e[j])] == d;
            ents.insert(entity[static_cast<std::size_t>(e[j])]);
        }
        std::size_t covered = 0;
        for (std::size_t c = 0; c < ncells && ok; ++c)
            if (cell_dim[c] == d && ents.count(entity[c]))
                ++covered;
        if (!ok || covered != reg->NumEntries()) {
            log::warn(
                "{}: cell region '{}' is not a union of whole geometric entities of one "
                "dimension; dropped",
                pFormat, reg->mName);
            continue;
        }
        s.mDim = d;
        s.mEntities.assign(ents.begin(), ents.end());
        selections.push_back(std::move(s));
    }

    const std::size_t nobjects = 1 + selections.size();
    rOut.Int(0);
    rOut.Int(1, "version");
    rOut.Int(static_cast<std::int64_t>(nobjects), "number of tags");
    rOut.String("mesh1");
    for (std::size_t k = 1; k < nobjects; ++k)
        rOut.String("mesh1_sel" + std::to_string(k));
    rOut.Int(static_cast<std::int64_t>(nobjects), "number of types");
    for (std::size_t k = 0; k < nobjects; ++k)
        rOut.String("obj");
    rOut.Blank();
    rOut.Comment("--------- Object 0 ----------");
    rOut.Int(0);
    rOut.Int(0);
    rOut.Int(1);
    rOut.String("Mesh", "class");
    rOut.Int(4, "version");
    rOut.Int(static_cast<std::int64_t>(sdim), "sdim");
    rOut.Int(static_cast<std::int64_t>(rMesh.NumPoints()), "number of mesh vertices");
    rOut.Int(0, "lowest mesh vertex index");
    rOut.Blank();
    rOut.Comment("Mesh vertex coordinates");
    const NDArray& points = rMesh.Points();
    std::vector<double> row(sdim);
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        for (std::size_t c = 0; c < sdim; ++c)
            row[c] = detail::read_double(points, i * sdim + c);
        rOut.Real(row.data(), sdim);
    }
    rOut.Blank();
    rOut.Int(static_cast<std::int64_t>(nblocks), "number of element types");
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        rOut.Blank();
        rOut.Comment("Type #" + std::to_string(b));
        rOut.String(ctypes[b], "type name");
        const NDArray& conn = cb.Conn();
        const std::size_t nn = detail::cols(conn);
        const std::size_t ne = cb.NumCells();
        rOut.Int(static_cast<std::int64_t>(nn), "number of vertices per element");
        rOut.Int(static_cast<std::int64_t>(ne), "number of elements");
        rOut.Comment("Elements");
        const detail::NodeOrder* order = detail::node_order("mphtxt", cb.Type());
        std::vector<std::int64_t> nodes(nn);
        for (std::size_t r = 0; r < ne; ++r) {
            for (std::size_t j = 0; j < nn; ++j)
                nodes[j] = detail::read_int(
                    conn, r * nn + (order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j));
            rOut.Ints(nodes.data(), nn);
        }
        rOut.Blank();
        rOut.Int(static_cast<std::int64_t>(ne), "number of geometric entity indices");
        rOut.Comment("Geometric entity indices");
        for (std::int64_t c = bases[b]; c < bases[b + 1]; ++c)
            rOut.Ints(&entity[static_cast<std::size_t>(c)], 1);
    }
    for (std::size_t k = 0; k < selections.size(); ++k) {
        const ComsolSelectionOut& s = selections[k];
        rOut.Blank();
        rOut.Comment("--------- Object " + std::to_string(k + 1) + " ----------");
        rOut.Int(0);
        rOut.Int(0);
        rOut.Int(1);
        rOut.String("Selection", "class");
        rOut.Int(0, "Version");
        rOut.String(s.mLabel, "Label");
        rOut.String("mesh1", "Geometry/mesh tag");
        rOut.Int(s.mDim, "Dimension");
        rOut.Int(static_cast<std::int64_t>(s.mEntities.size()), "Number of entities");
        rOut.Comment("Entities");
        for (std::int64_t e : s.mEntities)
            rOut.Ints(&e, 1);
    }
}

void comsol_check_data(const Mesh& rMesh, const char* pFormat) {
    std::size_t dropped = rMesh.PointDataNames().size() + rMesh.FieldDataNames().size();
    for (const std::string& n : rMesh.CellDataNames())
        dropped += n != "mphtxt:geom";
    if (dropped > 0)
        log::warn("{}: a COMSOL mesh holds no data arrays; {} array(s) dropped", pFormat, dropped);
}

}  // namespace

Mesh read_mphtxt(const std::string& rPath) {
    ComsolText in(comsol_slurp(rPath));
    return comsol_read(in, "mphtxt");
}

Mesh read_mphbin(const std::string& rPath) {
    ComsolBinary in(comsol_slurp(rPath));
    return comsol_read(in, "mphbin");
}

void write_mphtxt(const std::string& rPath, const Mesh& rMesh) {
    std::string out = detail::provenance_render_lines(detail::SlotTier::Block, "# ");
    ComsolTextSink sink(out);
    comsol_write(sink, rMesh, "mphtxt");
    comsol_check_data(rMesh, "mphtxt");
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    f << out;
    if (!f)
        throw WriteError("mphtxt: failed writing " + rPath);
}

void write_mphbin(const std::string& rPath, const Mesh& rMesh) {
    // No comment slot: consumes the pending record, and fails under
    // Mode::Required like the other slotless formats.
    detail::provenance_lines(detail::SlotTier::None);
    std::string out;
    ComsolBinarySink sink(out);
    comsol_write(sink, rMesh, "mphbin");
    comsol_check_data(rMesh, "mphbin");
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    f << out;
    if (!f)
        throw WriteError("mphbin: failed writing " + rPath);
}

}  // namespace meshioplusplus
