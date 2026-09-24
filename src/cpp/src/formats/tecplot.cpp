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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/detail/binary_stream.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

namespace meshioplusplus {

namespace {

std::string tecplot_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string tecplot_strip(const std::string& rS) {
    std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}
std::vector<std::string> tecplot_tokens(const std::string& rS) {
    std::vector<std::string> out;
    auto iss = detail::make_classic_istringstream(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}
bool is_float_token(const std::string& rS) {
    if (rS.empty())
        return false;
    const char* endp = nullptr;
    detail::parse_double(rS.c_str(), endp);
    return endp == rS.c_str() + rS.size();
}

/// Splits a joined header (e.g. `ZONE T = "VARLOCATION" N=4 ...`) into tokens,
/// keeping a quoted `"..."` span or a parenthesized `(...)` span as one
/// token (quotes/parens included) rather than letting their commas, `=` or
/// spaces fall through to the general splitter -- a title can legitimately
/// contain any of those (a real adversarial fixture: `T = "VARLOCATION"`,
/// where a naive substring search for the VARLOCATION *field* would instead
/// find this title). `=` becomes its own one-character token; everything
/// else splits on whitespace/commas exactly like tecplot_tokens.
std::vector<std::string> tecplot_header_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::size_t i = 0;
    const std::size_t n = rS.size();
    while (i < n) {
        const char c = rS[i];
        if (c == ' ' || c == '\t' || c == ',') {
            ++i;
        } else if (c == '=') {
            out.emplace_back("=");
            ++i;
        } else if (c == '"') {
            std::size_t j = rS.find('"', i + 1);
            if (j == std::string::npos)
                j = n - 1;
            out.push_back(rS.substr(i, j - i + 1));
            i = j + 1;
        } else if (c == '(') {
            std::size_t j = rS.find(')', i + 1);
            if (j == std::string::npos)
                j = n - 1;
            out.push_back(rS.substr(i, j - i + 1));
            i = j + 1;
        } else {
            std::size_t j = i;
            while (j < n && rS[j] != ' ' && rS[j] != '\t' && rS[j] != ',' && rS[j] != '=' &&
                  rS[j] != '"' && rS[j] != '(')
                ++j;
            out.push_back(rS.substr(i, j - i));
            i = j;
        }
    }
    return out;
}

/// Strips the surrounding quotes of a token tecplot_header_tokens quoted
/// (a no-op on a token that was never quoted).
std::string tecplot_unquote(const std::string& rTok) {
    if (rTok.size() >= 2 && rTok.front() == '"' && rTok.back() == '"')
        return rTok.substr(1, rTok.size() - 2);
    return rTok;
}

/// Splits `rInner` on commas that are not inside a `[...]` range, so
/// `"[1,2]=1,[4-6]=2"` splits into `{"[1,2]=1", "[4-6]=2"}` rather than
/// breaking the first range apart. Used for VARLOCATION/VARSHARELIST's
/// `([range]=target, ...)` value and for PASSIVEVARLIST's `([range])`.
std::vector<std::string> tecplot_split_entries(const std::string& rInner) {
    std::vector<std::string> out;
    std::string cur;
    int depth = 0;
    for (char c : rInner) {
        if (c == '[')
            ++depth;
        else if (c == ']')
            --depth;
        if (c == ',' && depth == 0) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty())
        out.push_back(cur);
    return out;
}

/// Expands a `[a,b,c-d,...]` (or bare, bracket-less) range spec into 0-based
/// indices. Tecplot variable/zone numbers are 1-based in the file.
std::vector<std::size_t> tecplot_parse_ranges(std::string rSpec) {
    if (!rSpec.empty() && rSpec.front() == '[')
        rSpec = rSpec.substr(1);
    if (!rSpec.empty() && rSpec.back() == ']')
        rSpec.pop_back();
    std::vector<std::size_t> out;
    for (const std::string& part : tecplot_split_entries(rSpec)) {
        const std::size_t dash = part.find('-');
        if (dash == std::string::npos) {
            if (!part.empty())
                out.push_back(static_cast<std::size_t>(std::stoul(part) - 1));
        } else {
            const int a = std::stoi(part.substr(0, dash));
            const int b = std::stoi(part.substr(dash + 1));
            for (int k = a; k <= b; ++k)
                out.push_back(static_cast<std::size_t>(k - 1));
        }
    }
    return out;
}

std::string tecplot_to_meshio(const std::string& rZ) {
    std::string u = tecplot_upper(rZ);
    if (u == "LINESEG" || u == "FELINESEG")
        return "line";
    if (u == "TRIANGLE" || u == "FETRIANGLE")
        return "triangle";
    if (u == "QUADRILATERAL" || u == "FEQUADRILATERAL")
        return "quad";
    if (u == "TETRAHEDRON" || u == "FETETRAHEDRON")
        return "tetra";
    if (u == "BRICK" || u == "FEBRICK")
        return "hexahedron";
    return "";
}
std::string meshio_to_tecplot(const std::string& rM) {
    if (rM == "line")
        return "FELINESEG";
    if (rM == "triangle")
        return "FETRIANGLE";
    if (rM == "quad")
        return "FEQUADRILATERAL";
    if (rM == "tetra")
        return "FETETRAHEDRON";
    if (rM == "pyramid" || rM == "wedge" || rM == "hexahedron")
        return "FEBRICK";
    if (rM.rfind("polygon", 0) == 0)
        return "FEPOLYGON";
    if (rM.rfind("polyhedron", 0) == 0)
        return "FEPOLYHEDRON";
    return "";
}
const std::vector<int>& tecplot_order(const std::string& rM) {
    static const std::map<std::string, std::vector<int>> o = {
        {"line", {0, 1}},
        {"triangle", {0, 1, 2}},
        {"quad", {0, 1, 2, 3}},
        {"tetra", {0, 1, 2, 3}},
        {"pyramid", {0, 1, 2, 3, 4, 4, 4, 4}},
        {"wedge", {0, 1, 4, 3, 2, 2, 5, 5}},
        {"hexahedron", {0, 1, 2, 3, 4, 5, 6, 7}},
    };
    static const std::vector<int> empty;
    auto it = o.find(rM);
    return it == o.end() ? empty : it->second;
}

}  // namespace

namespace {

// The zone model both the ASCII and the binary (.plt) reader decode into: a
// zone's header fields normalised (0-based sharing, ordered I/J/K, per-variable
// location), so the sharing resolution, the ordered connectivity, the
// timeline and the mesh assembly below are one code path for both encodings.
// read_tecplot and read_tecplot_metadata start from the same scan, so which
// zones a step holds and which zones the metadata lists never drift.
struct TecplotZone {
    std::string mTitle;                 // unquoted; empty if not given
    std::string mTypeName = "ORDERED";  // ORDERED, FELINESEG, ..., FEPOLYHEDRON
    bool mOrdered = false;
    std::size_t mI = 1, mJ = 1, mK = 1;
    std::size_t mNumNodes = 0;
    std::size_t mNumCells = 0;
    std::vector<int> mCellCentered;  // per variable
    // Variable index -> the 0-based zone it is shared from (VARSHARELIST):
    // that variable has no data of its own here.
    std::map<std::size_t, std::size_t> mVarShareZone;
    // Variables with no data anywhere for this zone (PASSIVEVARLIST).
    std::set<std::size_t> mPassiveVars;
    // 0-based zone this zone's connectivity is shared from, or -1.
    long long mConnShareZone = -1;
    bool mHasSolutionTime = false;
    double mSolutionTime = 0.0;
    bool mHasStrandId = false;
    int mStrandId = 0;

    // ASCII: BLOCK (else POINT) packing and the first data line.
    bool mBlock = true;
    std::size_t mDataStart = 0;
    // Binary: where each owned variable's values start and their data format
    // (1 float, 2 double, 3 int32, 4 int16, 5 byte), and the connectivity.
    std::map<std::size_t, std::pair<std::size_t, int>> mVarData;
    std::size_t mConnOffset = 0;
    bool mHasConn = false;
    int mRawFaceNeighbors = 0;
    int mMiscFaceNeighbors = 0;
    int mFaceNeighborMode = 0;
    // Face-based (FEPOLYGON/FEPOLYHEDRON) zones: the face map's sizes. The
    // boundary counts are as the encoding stores them (a .plt counts one
    // extra boundary face for "no neighbour" when there is any).
    std::size_t mNumFaces = 0;
    std::size_t mTotalFaceNodes = 0;
    std::size_t mNumBoundaryFaces = 0;
    std::size_t mNumBoundaryConns = 0;

    bool IsPoly() const { return mTypeName == "FEPOLYGON" || mTypeName == "FEPOLYHEDRON"; }
    bool IsPolyhedron() const { return mTypeName == "FEPOLYHEDRON"; }

    bool Owns(std::size_t Var) const {
        return !mVarShareZone.count(Var) && !mPassiveVars.count(Var);
    }
    std::size_t DataLength(std::size_t Var) const {
        return mCellCentered[Var] ? mNumCells : mNumNodes;
    }
    // Node and cell counts of an ordered zone: lines, quads or hexahedra over
    // the dimensions longer than one; a single point is one vertex.
    void FinishOrdered() {
        mNumNodes = mI * mJ * mK;
        std::size_t cells = 1;
        bool any = false;
        for (std::size_t d : {mI, mJ, mK}) {
            if (d > 1) {
                cells *= d - 1;
                any = true;
            }
        }
        mNumCells = any ? cells : 1;
    }
    std::size_t OrderedDims() const {
        return static_cast<std::size_t>(mI > 1) + static_cast<std::size_t>(mJ > 1) +
               static_cast<std::size_t>(mK > 1);
    }
};

std::string tecplot_zone_meshio_type(const TecplotZone& rZ) {
    if (rZ.mOrdered) {
        static const char* const kTypes[] = {"vertex", "line", "quad", "hexahedron"};
        return kTypes[rZ.OrderedDims()];
    }
    if (rZ.mTypeName == "FEPOLYGON")
        return "polygon";
    if (rZ.mTypeName == "FEPOLYHEDRON")
        return "polyhedron";
    const std::string mtype = tecplot_to_meshio(rZ.mTypeName);
    if (mtype.empty())
        throw ReadError("Tecplot: unsupported zone type " + rZ.mTypeName);
    return mtype;
}

std::size_t tecplot_nodes_per_cell(const std::string& rMeshioType) {
    if (rMeshioType == "vertex")
        return 1;
    if (rMeshioType == "line")
        return 2;
    if (rMeshioType == "triangle")
        return 3;
    if (rMeshioType == "quad" || rMeshioType == "tetra")
        return 4;
    return 8;
}

/// The cells of an ordered zone over its node index `i + I*(j + J*k)`, in VTK
/// corner order and `i` fastest (the order of its cell-centred values). Unit
/// dimensions are dropped first, so a J- or JK-ordered zone meshes like an I-
/// or IJ-ordered one.
NDArray tecplot_ordered_connectivity(const TecplotZone& rZ) {
    const std::size_t dims[3] = {rZ.mI, rZ.mJ, rZ.mK};
    const std::size_t stride[3] = {1, rZ.mI, rZ.mI * rZ.mJ};
    std::vector<std::size_t> n, s;
    for (int a = 0; a < 3; ++a) {
        if (dims[a] > 1) {
            n.push_back(dims[a]);
            s.push_back(stride[a]);
        }
    }
    if (n.empty()) {
        NDArray conn(DType::Int64, {1, 1});
        conn.As<std::int64_t>()[0] = 0;
        return conn;
    }
    const std::size_t npc = std::size_t{1} << n.size();
    NDArray conn(DType::Int64, {rZ.mNumCells, npc});
    std::int64_t* cp = conn.As<std::int64_t>();
    const std::size_t ni = n[0] - 1;
    const std::size_t nj = n.size() > 1 ? n[1] - 1 : 1;
    const std::size_t nk = n.size() > 2 ? n[2] - 1 : 1;
    const std::size_t sj = n.size() > 1 ? s[1] : 0;
    const std::size_t sk = n.size() > 2 ? s[2] : 0;
    std::size_t c = 0;
    for (std::size_t k = 0; k < nk; ++k) {
        for (std::size_t j = 0; j < nj; ++j) {
            for (std::size_t i = 0; i < ni; ++i, ++c) {
                const std::size_t b = i * s[0] + j * sj + k * sk;
                std::int64_t* row = cp + c * npc;
                if (n.size() == 1) {
                    row[0] = static_cast<std::int64_t>(b);
                    row[1] = static_cast<std::int64_t>(b + s[0]);
                    continue;
                }
                const std::size_t q[4] = {b, b + s[0], b + s[0] + sj, b + sj};
                for (int m = 0; m < 4; ++m)
                    row[m] = static_cast<std::int64_t>(q[m]);
                if (n.size() == 3)
                    for (int m = 0; m < 4; ++m)
                        row[4 + m] = static_cast<std::int64_t>(q[m] + sk);
            }
        }
    }
    return conn;
}

/// The face map of an FEPOLYGON/FEPOLYHEDRON zone: every face's nodes and the
/// two elements it separates, all 0-based. A neighbour in another zone (a
/// boundary connection) or none at all is -1: only this zone's cells are
/// assembled from it, and a face with no neighbour on one side still bounds
/// the cell on the other.
struct TecplotFaceMap {
    std::vector<std::int64_t> mStart;  // NumFaces + 1 offsets into mNodes
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mLeft;
    std::vector<std::int64_t> mRight;
};

/// Checks and stores a zone's face map from its raw arrays: `rCounts` is empty
/// for a polygonal zone (every face is an edge), node and element numbers are
/// 1-based in ASCII (0 = no neighbour) and 0-based in a .plt (-1 = none); a
/// negative element names a boundary connection, i.e. another zone.
TecplotFaceMap tecplot_face_map(const TecplotZone& rZ, std::size_t ZoneIdx,
                                const std::vector<std::int64_t>& rCounts,
                                std::vector<std::int64_t> Nodes, std::vector<std::int64_t> Left,
                                std::vector<std::int64_t> Right, bool OneBased) {
    const std::string where = "Tecplot: zone " + std::to_string(ZoneIdx + 1);
    const std::int64_t base = OneBased ? 1 : 0;
    TecplotFaceMap fm;
    fm.mStart.reserve(rZ.mNumFaces + 1);
    fm.mStart.push_back(0);
    for (std::size_t f = 0; f < rZ.mNumFaces; ++f) {
        const std::int64_t n = rCounts.empty() ? 2 : rCounts[f];
        if (n < 2)
            throw ReadError(where + ": face " + std::to_string(f + 1) + " has " +
                            std::to_string(n) + " nodes");
        fm.mStart.push_back(fm.mStart.back() + n);
    }
    if (static_cast<std::size_t>(fm.mStart.back()) != Nodes.size())
        throw ReadError(where + ": the face node counts add up to " +
                        std::to_string(fm.mStart.back()) + ", not TOTALNUMFACENODES " +
                        std::to_string(Nodes.size()));
    const auto num_nodes = static_cast<std::int64_t>(rZ.mNumNodes);
    for (std::int64_t& v : Nodes) {
        v -= base;
        if (v < 0 || v >= num_nodes)
            throw ReadError(where + ": face node " + std::to_string(v + base) + " is out of range");
    }
    const auto num_cells = static_cast<std::int64_t>(rZ.mNumCells);
    for (std::vector<std::int64_t>* pSide : {&Left, &Right}) {
        for (std::int64_t& e : *pSide) {
            if (e < 0) {
                e = -1;  // a boundary connection: the neighbour is in another zone
                continue;
            }
            e -= base;
            if (e >= num_cells)
                throw ReadError(where + ": face neighbour " + std::to_string(e + base) +
                                " is out of range");
        }
    }
    fm.mNodes = std::move(Nodes);
    fm.mLeft = std::move(Left);
    fm.mRight = std::move(Right);
    return fm;
}

/// How many integers a face-based zone's face map holds after its data: the
/// node count per face (polyhedra only), the face nodes, the left and the
/// right elements, then the boundary connections (ASCII: a count per
/// connected boundary face, then that many elements and as many zones).
std::size_t tecplot_ascii_face_map_tokens(const TecplotZone& rZ) {
    std::size_t n = (rZ.IsPolyhedron() ? rZ.mNumFaces : 0) + rZ.mTotalFaceNodes + 2 * rZ.mNumFaces;
    if (rZ.mNumBoundaryFaces > 0)
        n += rZ.mNumBoundaryFaces + 2 * rZ.mNumBoundaryConns;
    return n;
}

/// A zone's own data -- the variables it neither shares nor leaves passive,
/// and its own FE connectivity (0-based) or face map -- from wherever the
/// encoding keeps them.
class TecplotSource {
public:
    virtual ~TecplotSource() = default;
    /// Fills `rCols[v]` for every owned variable `v` and, for an FE zone that
    /// does not share its connectivity, `rConn` (cell-based zones) or `rFaces`
    /// (face-based zones).
    virtual void OwnData(std::size_t ZoneIdx, std::vector<std::vector<double>>& rCols,
                         NDArray& rConn, TecplotFaceMap& rFaces) const = 0;
};

// --- ASCII --------------------------------------------------------------------

/// The zone type, packing and cell-centred variables from its header fields.
/// `F=` is the old form (FEPOINT/FEBLOCK with `ET=`, or POINT/BLOCK for an
/// ordered zone); otherwise ZONETYPE (default ORDERED) and DATAPACKING
/// (default BLOCK, as the Data Format Guide specifies).
void tecplot_ascii_zone_kind(TecplotZone& rZ, const std::map<std::string, std::string>& rFields,
                             const std::string& rVarloc, std::size_t NumVariables) {
    auto field = [&](const char* pKey, const char* pDefault) {
        const auto it = rFields.find(pKey);
        return tecplot_upper(it != rFields.end() ? it->second : std::string(pDefault));
    };
    const std::string f = field("F", "");
    if (!f.empty()) {
        if (f == "FEPOINT" || f == "FEBLOCK") {
            const std::string et = field("ET", "");
            rZ.mTypeName = et.rfind("FE", 0) == 0 ? et : "FE" + et;
            rZ.mBlock = f == "FEBLOCK";
        } else {
            rZ.mTypeName = "ORDERED";
            rZ.mBlock = f == "BLOCK";
        }
    } else {
        rZ.mTypeName = field("ZONETYPE", "ORDERED");
        rZ.mBlock = field("DATAPACKING", "BLOCK") == "BLOCK";
    }
    rZ.mOrdered = rZ.mTypeName == "ORDERED";

    rZ.mCellCentered.assign(NumVariables, 0);
    if (!rZ.mBlock)
        return;  // POINT packing is nodal only
    if (rFields.count("NV") && !rZ.mOrdered) {
        const int nv = std::stoi(rFields.at("NV"));
        for (std::size_t k = static_cast<std::size_t>(nv); k < NumVariables; ++k)
            rZ.mCellCentered[k] = 1;
    } else if (!rVarloc.empty()) {
        const std::string vc = rVarloc.substr(1, rVarloc.size() - 2);  // strip outer ()
        for (const std::string& entry : tecplot_split_entries(vc)) {
            const std::size_t eq = entry.find('=');
            if (eq == std::string::npos)
                continue;
            if (tecplot_upper(entry.substr(eq + 1)) != "CELLCENTERED")
                continue;
            for (std::size_t idx : tecplot_parse_ranges(entry.substr(0, eq)))
                if (idx < NumVariables)
                    rZ.mCellCentered[idx] = 1;
        }
    }
}

// How many numeric tokens this zone's data block holds: BLOCK is one run per
// *owned* variable (cell-centred ones NumCells long, the rest NumNodes long);
// a shared or passive variable has no data here. POINT is NumNodes rows of
// every variable.
std::size_t tecplot_ascii_token_count(const TecplotZone& rZ, std::size_t NumVariables) {
    if (!rZ.mBlock)
        return rZ.mNumNodes * NumVariables;
    std::size_t total = 0;
    for (std::size_t k = 0; k < NumVariables; ++k)
        if (rZ.Owns(k))
            total += rZ.DataLength(k);
    return total;
}

// One pass over every line, splitting VARIABLES from every ZONE header and
// locating each zone's data section by counting off its own token budget plus
// its connectivity lines -- Tecplot ASCII has no fixed tokens-per-line
// convention, so this is the only reliable way to find where one zone's data
// ends and the next one's header begins. Anything else (TEXT, GEOMETRY,
// DATASETAUXDATA, face-neighbour lines) is skipped.
std::vector<TecplotZone> tecplot_scan_zones(const std::vector<std::string>& rLines,
                                            std::vector<std::string>& rVariables) {
    std::vector<TecplotZone> zones;
    std::size_t i = 0;
    while (i < rLines.size()) {
        const std::string u = tecplot_upper(rLines[i]);
        if (u.rfind("VARIABLES", 0) == 0) {
            std::string joined = rLines[i];
            while (i + 1 < rLines.size() && rLines[i + 1][0] == '"')
                joined += " " + rLines[++i];
            rVariables.clear();
            const std::string rhs = joined.substr(joined.find('=') + 1);
            std::size_t p = 0;
            while (p < rhs.size()) {
                if (rhs[p] == '"') {
                    std::size_t q = rhs.find('"', p + 1);
                    if (q == std::string::npos)
                        q = rhs.size();
                    rVariables.push_back(rhs.substr(p + 1, q - p - 1));
                    p = q + 1;
                } else if (std::isspace((unsigned char)rhs[p]) || rhs[p] == ',') {
                    ++p;
                } else {
                    std::size_t q = p;
                    while (q < rhs.size() && !std::isspace((unsigned char)rhs[q]) && rhs[q] != ',')
                        ++q;
                    rVariables.push_back(rhs.substr(p, q - p));
                    p = q;
                }
            }
            ++i;
            continue;
        }
        if (u.rfind("ZONE", 0) != 0) {
            ++i;
            continue;
        }

        TecplotZone z;
        std::string joined = rLines[i];
        while (i + 1 < rLines.size() && !is_float_token(tecplot_tokens(rLines[i + 1])[0]))
            joined += " " + rLines[++i];
        z.mDataStart = i + 1;

        // tk[0] is the "ZONE" keyword itself; the rest is KEY = VALUE triples,
        // where a quoted or parenthesized VALUE is already one token (so a
        // title like T = "VARLOCATION" cannot be mistaken for the field of
        // the same name).
        std::map<std::string, std::string> fields;
        std::string varloc;
        const std::vector<std::string> tk = tecplot_header_tokens(joined);
        for (std::size_t k = 1; k + 2 < tk.size();) {
            if (tk[k + 1] != "=") {
                ++k;  // resync: not a KEY = VALUE triple after all
                continue;
            }
            const std::string key = tecplot_upper(tk[k]);
            const std::string& val = tk[k + 2];
            if (key == "NODES" || key == "N" || key == "ELEMENTS" || key == "E" ||
                key == "DATAPACKING" || key == "ZONETYPE" || key == "F" || key == "ET" ||
                key == "NV" || key == "I" || key == "J" || key == "K" || key == "FACES" ||
                key == "TOTALNUMFACENODES" || key == "NUMCONNECTEDBOUNDARYFACES" ||
                key == "TOTALNUMBOUNDARYCONNECTIONS") {
                fields[key] = val;
            } else if (key == "T") {
                z.mTitle = tecplot_unquote(val);
            } else if (key == "SOLUTIONTIME") {
                z.mSolutionTime = detail::parse_double(val);
                z.mHasSolutionTime = true;
            } else if (key == "STRANDID") {
                z.mStrandId = std::stoi(val);
                z.mHasStrandId = true;
            } else if (key == "VARLOCATION") {
                varloc = val;
                varloc.erase(std::remove(varloc.begin(), varloc.end(), ' '), varloc.end());
            } else if (key == "VARSHARELIST") {
                const std::string inner = val.substr(1, val.size() - 2);  // strip outer ()
                for (const std::string& entry : tecplot_split_entries(inner)) {
                    const std::size_t eq = entry.find('=');
                    // No zone number: shared from the previous zone.
                    const long long src = eq == std::string::npos
                                              ? static_cast<long long>(zones.size()) - 1
                                              : std::stoll(entry.substr(eq + 1)) - 1;
                    if (src < 0)
                        throw ReadError("Tecplot: VARSHARELIST names no earlier zone");
                    for (std::size_t idx : tecplot_parse_ranges(entry.substr(0, eq)))
                        z.mVarShareZone[idx] = static_cast<std::size_t>(src);
                }
            } else if (key == "PASSIVEVARLIST") {
                const std::string inner = val.substr(1, val.size() - 2);  // strip outer ()
                for (std::size_t idx : tecplot_parse_ranges(inner))
                    z.mPassiveVars.insert(idx);
            } else if (key == "CONNECTIVITYSHAREZONE") {
                z.mConnShareZone = std::stoll(val) - 1;
            }
            k += 3;
        }
        if (rVariables.empty())
            throw ReadError("Tecplot: no VARIABLES");
        tecplot_ascii_zone_kind(z, fields, varloc, rVariables.size());
        try {
            if (z.mOrdered) {
                z.mI = fields.count("I") ? std::stoull(fields.at("I")) : 1;
                z.mJ = fields.count("J") ? std::stoull(fields.at("J")) : 1;
                z.mK = fields.count("K") ? std::stoull(fields.at("K")) : 1;
                z.FinishOrdered();
            } else {
                auto get = [&](const char* pA, const char* pB) {
                    auto it = fields.find(pA);
                    if (it == fields.end())
                        it = fields.find(pB);
                    return it == fields.end() ? std::string() : it->second;
                };
                z.mNumNodes = std::stoull(get("NODES", "N"));
                z.mNumCells = std::stoull(get("ELEMENTS", "E"));
            }
        } catch (const std::logic_error&) {
            throw ReadError(z.mOrdered ? "Tecplot: bad I/J/K in an ordered zone"
                                       : "Tecplot: an FE zone needs NODES and ELEMENTS");
        }
        if (z.IsPoly()) {
            auto count = [&](const char* pKey, bool Required, std::size_t Default) {
                const auto it = fields.find(pKey);
                if (it == fields.end()) {
                    if (Required)
                        throw ReadError("Tecplot: a " + z.mTypeName + " zone needs " + pKey);
                    return Default;
                }
                try {
                    return static_cast<std::size_t>(std::stoull(it->second));
                } catch (const std::logic_error&) {
                    throw ReadError(std::string("Tecplot: bad ") + pKey);
                }
            };
            z.mNumFaces = count("FACES", true, 0);
            z.mTotalFaceNodes = count("TOTALNUMFACENODES", z.IsPolyhedron(), 2 * z.mNumFaces);
            z.mNumBoundaryFaces = count("NUMCONNECTEDBOUNDARYFACES", false, 0);
            z.mNumBoundaryConns = count("TOTALNUMBOUNDARYCONNECTIONS", false, 0);
            if (!z.mBlock)
                throw ReadError("Tecplot: a " + z.mTypeName + " zone must be DATAPACKING=BLOCK");
        }

        const std::size_t want = tecplot_ascii_token_count(z, rVariables.size());
        std::size_t li = z.mDataStart, got = 0;
        while (got < want && li < rLines.size()) {
            got += tecplot_tokens(rLines[li]).size();
            ++li;
        }
        if (!z.mOrdered && z.mConnShareZone < 0) {
            if (z.IsPoly()) {
                // The face map is a token stream, not one line per anything.
                const std::size_t face_tokens = tecplot_ascii_face_map_tokens(z);
                std::size_t seen = 0;
                while (seen < face_tokens && li < rLines.size())
                    seen += tecplot_tokens(rLines[li++]).size();
            } else {
                li += z.mNumCells;  // one connectivity line per cell -- none if shared
            }
        }
        zones.push_back(std::move(z));
        i = li;
    }
    if (rVariables.empty())
        throw ReadError("Tecplot: no VARIABLES");
    if (zones.empty())
        throw ReadError("Tecplot: no ZONE");
    return zones;
}

class TecplotAsciiSource final : public TecplotSource {
public:
    TecplotAsciiSource(const std::vector<std::string>& rLines,
                       const std::vector<TecplotZone>& rZones, std::size_t NumVariables)
        : mrLines(rLines), mrZones(rZones), mNumVariables(NumVariables) {}

    void OwnData(std::size_t ZoneIdx, std::vector<std::vector<double>>& rCols, NDArray& rConn,
                 TecplotFaceMap& rFaces) const override {
        const TecplotZone& z = mrZones[ZoneIdx];
        const std::size_t nv = mNumVariables;
        const std::size_t want = tecplot_ascii_token_count(z, nv);
        std::vector<double> flat;
        flat.reserve(want);
        std::size_t li = z.mDataStart;
        while (flat.size() < want && li < mrLines.size()) {
            for (const auto& t : tecplot_tokens(mrLines[li]))
                flat.push_back(detail::parse_double(t));
            ++li;
        }
        if (flat.size() < want)
            throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                            " has fewer values than its header announces");
        if (z.mBlock) {
            std::size_t off = 0;
            for (std::size_t k = 0; k < nv; ++k) {
                if (!z.Owns(k))
                    continue;
                const std::size_t n = z.DataLength(k);
                rCols[k].assign(flat.begin() + static_cast<std::ptrdiff_t>(off),
                                flat.begin() + static_cast<std::ptrdiff_t>(off + n));
                off += n;
            }
        } else {
            for (std::size_t k = 0; k < nv; ++k) {
                if (!z.Owns(k))
                    continue;
                rCols[k].resize(z.mNumNodes);
                for (std::size_t r = 0; r < z.mNumNodes; ++r)
                    rCols[k][r] = flat[r * nv + k];
            }
        }
        if (z.mOrdered || z.mConnShareZone >= 0)
            return;
        if (z.IsPoly()) {
            const std::size_t want_ints = tecplot_ascii_face_map_tokens(z);
            std::vector<std::int64_t> ints;
            ints.reserve(want_ints);
            while (ints.size() < want_ints && li < mrLines.size())
                for (const auto& t : tecplot_tokens(mrLines[li++]))
                    ints.push_back(std::strtoll(t.c_str(), nullptr, 10));
            if (ints.size() < want_ints)
                throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                                " face map is truncated");
            auto take = [&, pos = std::size_t{0}](std::size_t N) mutable {
                std::vector<std::int64_t> out(ints.begin() + static_cast<std::ptrdiff_t>(pos),
                                              ints.begin() + static_cast<std::ptrdiff_t>(pos + N));
                pos += N;
                return out;
            };
            const std::vector<std::int64_t> counts =
                z.IsPolyhedron() ? take(z.mNumFaces) : std::vector<std::int64_t>{};
            std::vector<std::int64_t> nodes = take(z.mTotalFaceNodes);
            std::vector<std::int64_t> left = take(z.mNumFaces);
            std::vector<std::int64_t> right = take(z.mNumFaces);
            rFaces = tecplot_face_map(z, ZoneIdx, counts, std::move(nodes), std::move(left),
                                      std::move(right), true);
            return;
        }
        const std::size_t nn = tecplot_nodes_per_cell(tecplot_zone_meshio_type(z));
        rConn = NDArray(DType::Int64, {z.mNumCells, nn});
        std::int64_t* cp = rConn.As<std::int64_t>();
        for (std::size_t c = 0; c < z.mNumCells; ++c) {
            if (li >= mrLines.size())
                throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                                " connectivity is truncated");
            const auto t = tecplot_tokens(mrLines[li++]);
            if (t.size() < nn)
                throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                                " has a short connectivity line");
            for (std::size_t j = 0; j < nn; ++j)
                cp[c * nn + j] = std::strtoll(t[j].c_str(), nullptr, 10) - 1;
        }
    }

private:
    const std::vector<std::string>& mrLines;
    const std::vector<TecplotZone>& mrZones;
    std::size_t mNumVariables;
};

// --- binary (.plt) --------------------------------------------------------------
//
// Appendix A of the Data Format Guide ("Binary Data File Format"), #!TDV112.
// Facts checked against files TecIO (the library preplot is built on) wrote,
// in both byte orders:
// * every header string is int32 per character, 0-terminated;
// * a zone record's strand is 0-based (-1 static), its "not used" word -1;
// * a geometry record carries coordinate system, scope and draw order before
//   its anchor (the order TecIO's own reader expects for V112);
// * the data section repeats the 299.0 marker per zone, then per-variable data
//   formats, the passive and sharing lists, the connectivity share zone and
//   one min/max pair per owned variable before the values;
// * a cell-centred variable of an ordered zone is stored over the node
//   dimensions with the last one longer than 1 shortened by one (note 5), so
//   the other directions carry a ghost value at their end: I-1 values for an
//   I-ordered zone, I*(J-1) for IJ, I*J*(K-1) for IJK;
// * FE connectivity is zero-based int32.

constexpr float kTecplotZoneMarker = 299.0F;
constexpr float kTecplotGeometryMarker = 399.0F;
constexpr float kTecplotTextMarker = 499.0F;
constexpr float kTecplotLabelMarker = 599.0F;
constexpr float kTecplotUserRecMarker = 699.0F;
constexpr float kTecplotDataSetAuxMarker = 799.0F;
constexpr float kTecplotVarAuxMarker = 899.0F;
constexpr float kTecplotEohMarker = 357.0F;

bool tecplot_is_plt(const char* pData, std::size_t Size) {
    return Size >= 5 && std::memcmp(pData, "#!TDV", 5) == 0;
}

std::string tecplot_plt_string(detail::ByteCursor& rCur) {
    std::string s;
    for (std::int32_t c = rCur.I32(); c != 0; c = rCur.I32())
        s += static_cast<char>(c);
    return s;
}

void tecplot_plt_skip_geometry(detail::ByteCursor& rCur) {
    const std::int32_t coord_sys = rCur.I32();  // 4 = Grid3D: polylines carry Z too
    rCur.Skip(2 * 4);                           // scope, draw order
    rCur.Skip(3 * 8);                           // anchor
    rCur.Skip(4 * 4);                           // zone, color, fill color, is filled
    const std::int32_t gtype = rCur.I32();
    rCur.Skip(4);              // line pattern
    rCur.Skip(2 * 8);          // pattern length, line thickness
    rCur.Skip(3 * 4);          // ellipse points, arrowhead style and attachment
    rCur.Skip(2 * 8);          // arrowhead size and angle
    tecplot_plt_string(rCur);  // macro function command
    const std::size_t width = rCur.I32() == 2 ? 8 : 4;
    rCur.Skip(4);  // clipping
    if (gtype == 0) {
        const std::int32_t lines = rCur.I32();
        for (std::int32_t l = 0; l < lines; ++l) {
            const std::size_t n = static_cast<std::size_t>(rCur.I32());
            rCur.Skip((coord_sys == 4 ? 3 : 2) * n * width);
        }
    } else if (gtype == 1 || gtype == 4) {
        rCur.Skip(2 * width);  // rectangle / ellipse
    } else if (gtype == 2 || gtype == 3) {
        rCur.Skip(width);  // square / circle
    } else {
        throw ReadError("Tecplot .plt: unknown geometry type " + std::to_string(gtype));
    }
}

void tecplot_plt_skip_text(detail::ByteCursor& rCur) {
    rCur.Skip(2 * 4);          // coordinate system, scope
    rCur.Skip(3 * 8);          // anchor
    rCur.Skip(2 * 4);          // font, height units
    rCur.Skip(8);              // height
    rCur.Skip(4);              // box type
    rCur.Skip(2 * 8);          // box margin, line width
    rCur.Skip(2 * 4);          // box outline, fill colors
    rCur.Skip(2 * 8);          // angle, line spacing
    rCur.Skip(3 * 4);          // anchor, zone, color
    tecplot_plt_string(rCur);  // macro function command
    rCur.Skip(4);              // clipping
    tecplot_plt_string(rCur);  // the text
}

// The header section: variables and every zone record, up to the 357.0 marker.
std::vector<TecplotZone> tecplot_plt_header(detail::ByteCursor& rCur,
                                            std::vector<std::string>& rVariables) {
    rCur.Skip(4);              // FileType: full, grid or solution -- all read the same way
    tecplot_plt_string(rCur);  // title
    const std::int32_t nvar = rCur.I32();
    if (nvar <= 0)
        throw ReadError("Tecplot .plt: bad variable count " + std::to_string(nvar));
    for (std::int32_t v = 0; v < nvar; ++v)
        rVariables.push_back(tecplot_plt_string(rCur));
    const std::size_t nv = rVariables.size();

    static const char* const kTypes[] = {"ORDERED",         "FELINESEG",     "FETRIANGLE",
                                         "FEQUADRILATERAL", "FETETRAHEDRON", "FEBRICK",
                                         "FEPOLYGON",       "FEPOLYHEDRON"};
    std::vector<TecplotZone> zones;
    for (;;) {
        const float marker = rCur.F32();
        if (marker == kTecplotZoneMarker) {
            TecplotZone z;
            z.mTitle = tecplot_plt_string(rCur);
            rCur.Skip(4);  // parent zone
            const std::int32_t strand = rCur.I32();
            const double time = rCur.F64();
            rCur.Skip(4);  // not used (-1)
            const std::int32_t ztype = rCur.I32();
            if (ztype < 0 || ztype > 7)
                throw ReadError("Tecplot .plt: unknown zone type " + std::to_string(ztype));
            z.mTypeName = kTypes[ztype];
            z.mOrdered = ztype == 0;
            z.mCellCentered.assign(nv, 0);
            if (rCur.I32() == 1)
                for (std::size_t v = 0; v < nv; ++v)
                    z.mCellCentered[v] = rCur.I32() == 1 ? 1 : 0;
            z.mRawFaceNeighbors = rCur.I32();
            z.mMiscFaceNeighbors = rCur.I32();
            if (z.mMiscFaceNeighbors != 0) {
                z.mFaceNeighborMode = rCur.I32();
                if (!z.mOrdered)
                    rCur.Skip(4);  // FE face neighbours completely specified
            }
            if (z.mOrdered) {
                const std::int32_t I = rCur.I32(), J = rCur.I32(), K = rCur.I32();
                if (I < 1 || J < 1 || K < 1)
                    throw ReadError("Tecplot .plt: bad I/J/K in an ordered zone");
                z.mI = static_cast<std::size_t>(I);
                z.mJ = static_cast<std::size_t>(J);
                z.mK = static_cast<std::size_t>(K);
                z.FinishOrdered();
            } else {
                const std::int32_t pts = rCur.I32();
                if (ztype == 6 || ztype == 7) {
                    // faces, face nodes, boundary faces (+1 when any) and connections
                    std::int32_t counts[4];
                    for (std::int32_t& c : counts) {
                        c = rCur.I32();
                        if (c < 0)
                            throw ReadError("Tecplot .plt: negative face map size");
                    }
                    z.mNumFaces = static_cast<std::size_t>(counts[0]);
                    z.mTotalFaceNodes = static_cast<std::size_t>(counts[1]);
                    z.mNumBoundaryFaces = static_cast<std::size_t>(counts[2]);
                    z.mNumBoundaryConns = static_cast<std::size_t>(counts[3]);
                }
                const std::int32_t elems = rCur.I32();
                if (pts < 0 || elems < 0)
                    throw ReadError("Tecplot .plt: negative node or element count");
                z.mNumNodes = static_cast<std::size_t>(pts);
                z.mNumCells = static_cast<std::size_t>(elems);
                rCur.Skip(3 * 4);  // I/J/K cell dimensions, unused
            }
            while (rCur.I32() == 1) {  // auxiliary name/value pairs
                tecplot_plt_string(rCur);
                rCur.Skip(4);
                tecplot_plt_string(rCur);
            }
            // The file's strands are 0-based (-1 static); only the grouping by
            // strand matters, so the ASCII STRANDID's 1-based numbering is moot.
            z.mHasSolutionTime = strand != -1 || time != 0.0;
            z.mSolutionTime = time;
            z.mHasStrandId = strand >= 0;
            z.mStrandId = strand;
            zones.push_back(std::move(z));
        } else if (marker == kTecplotGeometryMarker) {
            tecplot_plt_skip_geometry(rCur);
        } else if (marker == kTecplotTextMarker) {
            tecplot_plt_skip_text(rCur);
        } else if (marker == kTecplotLabelMarker) {
            const std::int32_t n = rCur.I32();
            for (std::int32_t l = 0; l < n; ++l)
                tecplot_plt_string(rCur);
        } else if (marker == kTecplotUserRecMarker) {
            tecplot_plt_string(rCur);
        } else if (marker == kTecplotDataSetAuxMarker) {
            tecplot_plt_string(rCur);
            rCur.Skip(4);
            tecplot_plt_string(rCur);
        } else if (marker == kTecplotVarAuxMarker) {
            rCur.Skip(4);
            tecplot_plt_string(rCur);
            rCur.Skip(4);
            tecplot_plt_string(rCur);
        } else if (marker == kTecplotEohMarker) {
            break;
        } else {
            throw ReadError("Tecplot .plt: unexpected header marker " + std::to_string(marker) +
                            " at byte " + std::to_string(rCur.Offset() - 4));
        }
    }
    if (zones.empty())
        throw ReadError("Tecplot: no ZONE");
    return zones;
}

std::size_t tecplot_plt_format_width(int Format) {
    switch (Format) {
        case 1:
        case 3:
            return 4;
        case 2:
            return 8;
        case 4:
            return 2;
        case 5:
            return 1;
        default:
            return 0;
    }
}

/// The stored dimensions (and length) of an ordered zone's cell-centred
/// variable: the node dimensions with the last one longer than 1 shortened.
std::size_t tecplot_plt_ordered_cc_stored(const TecplotZone& rZ, std::size_t* pStored) {
    std::size_t dims[3] = {rZ.mI, rZ.mJ, rZ.mK};
    for (int a = 2; a >= 0; --a) {
        if (dims[a] > 1) {
            --dims[a];
            break;
        }
    }
    for (int a = 0; a < 3; ++a)
        pStored[a] = dims[a];
    return dims[0] * dims[1] * dims[2];
}

void tecplot_plt_skip_face_neighbors(detail::ByteCursor& rCur, const TecplotZone& rZ) {
    for (int c = 0; c < rZ.mMiscFaceNeighbors; ++c) {
        switch (rZ.mFaceNeighborMode) {
            case 0:  // local one-to-one: cz, fz, cz
                rCur.Skip(3 * 4);
                break;
            case 2:  // global one-to-one: cz, fz, ZZ, CZ
                rCur.Skip(4 * 4);
                break;
            case 1:
            case 3: {  // one-to-many: cz, fz, oz, nz, then nz cells (or zone/cell pairs)
                rCur.Skip(3 * 4);
                const std::int32_t nz = rCur.I32();
                if (nz < 0)
                    throw ReadError("Tecplot .plt: bad face-neighbour count");
                rCur.Skip(static_cast<std::size_t>(nz) * (rZ.mFaceNeighborMode == 1 ? 4 : 8));
                break;
            }
            default:
                throw ReadError("Tecplot .plt: unknown face-neighbour mode " +
                                std::to_string(rZ.mFaceNeighborMode));
        }
    }
}

// Walks every zone's data section, recording where each owned variable and the
// connectivity start, and the sharing a binary file keeps here, not in the
// header.
void tecplot_plt_scan_data(detail::ByteCursor& rCur, const std::vector<std::string>& rVariables,
                           std::vector<TecplotZone>& rZones) {
    const std::size_t nv = rVariables.size();
    for (std::size_t zi = 0; zi < rZones.size(); ++zi) {
        TecplotZone& z = rZones[zi];
        const float marker = rCur.F32();
        if (marker != kTecplotZoneMarker)
            throw ReadError("Tecplot .plt: expected the data marker of zone " +
                            std::to_string(zi + 1) + ", got " + std::to_string(marker));
        std::vector<int> formats(nv);
        for (std::size_t v = 0; v < nv; ++v)
            formats[v] = rCur.I32();
        if (rCur.I32() != 0)
            for (std::size_t v = 0; v < nv; ++v)
                if (rCur.I32() != 0)
                    z.mPassiveVars.insert(v);
        if (rCur.I32() != 0) {
            for (std::size_t v = 0; v < nv; ++v) {
                const std::int32_t src = rCur.I32();
                if (src >= 0)
                    z.mVarShareZone[v] = static_cast<std::size_t>(src);
            }
        }
        z.mConnShareZone = rCur.I32();
        std::size_t owned = 0;
        for (std::size_t v = 0; v < nv; ++v)
            owned += z.Owns(v) ? 1 : 0;
        rCur.Skip(owned * 2 * 8);  // min/max pairs
        for (std::size_t v = 0; v < nv; ++v) {
            if (!z.Owns(v))
                continue;
            if (formats[v] == 6)
                throw ReadError("Tecplot .plt: variable '" + rVariables[v] +
                                "' is BIT-packed, which is not supported");
            const std::size_t width = tecplot_plt_format_width(formats[v]);
            if (width == 0)
                throw ReadError("Tecplot .plt: unknown data format " + std::to_string(formats[v]));
            std::size_t n = z.DataLength(v);
            if (z.mOrdered && z.mCellCentered[v]) {
                std::size_t stored[3];
                n = tecplot_plt_ordered_cc_stored(z, stored);
            }
            z.mVarData[v] = {rCur.Offset(), formats[v]};
            rCur.Skip(n * width);
        }
        if (z.mOrdered) {
            if (z.mConnShareZone < 0 && z.mMiscFaceNeighbors != 0)
                tecplot_plt_skip_face_neighbors(rCur, z);
            continue;
        }
        const std::string mtype = tecplot_zone_meshio_type(z);
        if (z.IsPoly()) {
            if (z.mConnShareZone < 0) {
                // Face node offsets (polyhedra only), face nodes, left and right
                // elements, then the boundary connections: offsets over the
                // (already +1) boundary face count plus one, elements and zones.
                z.mHasConn = true;
                z.mConnOffset = rCur.Offset();
                std::size_t n =
                    (z.IsPolyhedron() ? z.mNumFaces + 1 : 0) + z.mTotalFaceNodes + 2 * z.mNumFaces;
                if (z.mNumBoundaryFaces > 0)
                    n += z.mNumBoundaryFaces + 1 + 2 * z.mNumBoundaryConns;
                rCur.Skip(n * 4);
            }
            continue;
        }
        if (z.mConnShareZone < 0) {
            static const std::map<std::string, std::size_t> kFaces = {
                {"line", 0}, {"triangle", 3}, {"quad", 4}, {"tetra", 4}, {"hexahedron", 6}};
            z.mHasConn = true;
            z.mConnOffset = rCur.Offset();
            rCur.Skip(tecplot_nodes_per_cell(mtype) * z.mNumCells * 4);
            if (z.mRawFaceNeighbors != 0)
                rCur.Skip(kFaces.at(mtype) * z.mNumCells * 4);
            if (z.mMiscFaceNeighbors != 0)
                tecplot_plt_skip_face_neighbors(rCur, z);
        }
    }
}

class TecplotPltSource final : public TecplotSource {
public:
    TecplotPltSource(const char* pData, std::size_t Size, bool BigEndian,
                     const std::vector<TecplotZone>& rZones)
        : mpData(pData), mSize(Size), mBigEndian(BigEndian), mrZones(rZones) {}

    void OwnData(std::size_t ZoneIdx, std::vector<std::vector<double>>& rCols, NDArray& rConn,
                 TecplotFaceMap& rFaces) const override {
        const TecplotZone& z = mrZones[ZoneIdx];
        detail::ByteCursor cur(mpData, mSize, mBigEndian, "Tecplot .plt");
        const bool swap = mBigEndian != (std::endian::native == std::endian::big);
        for (const auto& [v, where] : z.mVarData) {
            const auto [offset, format] = where;
            std::size_t stored[3] = {0, 0, 0};
            const bool ghosts = z.mOrdered && z.mCellCentered[v];
            const std::size_t n =
                ghosts ? tecplot_plt_ordered_cc_stored(z, stored) : z.DataLength(v);
            std::vector<double> raw(n);
            cur.Seek(offset);
            for (std::size_t r = 0; r < n; ++r) {
                switch (format) {
                    case 1:
                        raw[r] = cur.F32();
                        break;
                    case 2:
                        raw[r] = cur.F64();
                        break;
                    case 3:
                        raw[r] = cur.I32();
                        break;
                    case 4: {
                        std::uint16_t u;
                        std::memcpy(&u, cur.Bytes(2).data(), 2);
                        if (swap)
                            u = detail::bswap16(u);
                        raw[r] = static_cast<std::int16_t>(u);
                        break;
                    }
                    default:
                        raw[r] = static_cast<unsigned char>(cur.Bytes(1)[0]);
                        break;
                }
            }
            if (!ghosts) {
                rCols[v] = std::move(raw);
                continue;
            }
            // Keep i < I-1, j < J-1, k < K-1 (index 0 of a unit dimension).
            const std::size_t keep_i = z.mI > 1 ? z.mI - 1 : 1;
            const std::size_t keep_j = z.mJ > 1 ? z.mJ - 1 : 1;
            const std::size_t keep_k = z.mK > 1 ? z.mK - 1 : 1;
            std::vector<double>& col = rCols[v];
            col.clear();
            col.reserve(z.mNumCells);
            for (std::size_t k = 0; k < keep_k; ++k)
                for (std::size_t j = 0; j < keep_j; ++j)
                    for (std::size_t i = 0; i < keep_i; ++i)
                        col.push_back(raw[i + stored[0] * (j + stored[1] * k)]);
        }
        if (!z.mHasConn)
            return;
        cur.Seek(z.mConnOffset);
        if (z.IsPoly()) {
            auto ints = [&](std::size_t N) {
                std::vector<std::int64_t> out(N);
                for (std::int64_t& v : out)
                    v = cur.I32();
                return out;
            };
            std::vector<std::int64_t> counts;
            if (z.IsPolyhedron()) {
                const std::vector<std::int64_t> offsets = ints(z.mNumFaces + 1);
                if (offsets[0] != 0)
                    throw ReadError("Tecplot .plt: zone " + std::to_string(ZoneIdx + 1) +
                                    " face node offsets do not start at 0");
                counts.resize(z.mNumFaces);
                for (std::size_t f = 0; f < z.mNumFaces; ++f)
                    counts[f] = offsets[f + 1] - offsets[f];
            }
            std::vector<std::int64_t> nodes = ints(z.mTotalFaceNodes);
            std::vector<std::int64_t> left = ints(z.mNumFaces);
            std::vector<std::int64_t> right = ints(z.mNumFaces);
            rFaces = tecplot_face_map(z, ZoneIdx, counts, std::move(nodes), std::move(left),
                                      std::move(right), false);
            return;
        }
        const std::size_t nn = tecplot_nodes_per_cell(tecplot_zone_meshio_type(z));
        rConn = NDArray(DType::Int64, {z.mNumCells, nn});
        std::int64_t* cp = rConn.As<std::int64_t>();
        for (std::size_t r = 0; r < z.mNumCells * nn; ++r)
            cp[r] = cur.I32();
    }

private:
    const char* mpData;
    std::size_t mSize;
    bool mBigEndian;
    const std::vector<TecplotZone>& mrZones;
};

// --- shared -------------------------------------------------------------------------

// One entry per step, each the zone indices (file order preserved) that make
// up that step's mesh. A non-transient file (the first zone carries no
// SOLUTIONTIME) is a single "step" listing every zone -- the several-static-
// parts case. A transient file's steps are its distinct SOLUTIONTIME values,
// ascending, each grouping every zone at that time sharing zones[0]'s strand
// when both carry one; a zone with no SOLUTIONTIME in an otherwise transient
// file belongs to no step and is dropped.
std::vector<std::vector<std::size_t>> tecplot_timeline(const std::vector<TecplotZone>& rZones) {
    if (!rZones[0].mHasSolutionTime) {
        std::vector<std::size_t> all(rZones.size());
        for (std::size_t k = 0; k < rZones.size(); ++k)
            all[k] = k;
        return {std::move(all)};
    }
    std::map<double, std::vector<std::size_t>> by_time;
    for (std::size_t k = 0; k < rZones.size(); ++k) {
        if (!rZones[k].mHasSolutionTime)
            continue;
        if (rZones[0].mHasStrandId && rZones[k].mHasStrandId &&
            rZones[k].mStrandId != rZones[0].mStrandId)
            continue;
        by_time[rZones[k].mSolutionTime].push_back(k);
    }
    std::vector<std::vector<std::size_t>> steps;
    steps.reserve(by_time.size());
    for (auto& [time, idxs] : by_time)  // std::map is sorted by key (time)
        steps.push_back(std::move(idxs));
    return steps;
}

/// One cell block of a zone. A cell-based or ordered zone is one piece over
/// all its cells; a face-based zone is one ragged `polygon` piece, or one
/// `polyhedron<N>` piece per distinct node count (the naming every other
/// polyhedral reader uses), each listing which of the zone's cells it holds.
struct TecplotPiece {
    std::string mType;
    NDArray mConn;                                              // rectangular pieces
    std::vector<std::vector<std::int64_t>> mRows;               // polygon
    std::vector<std::vector<std::vector<std::int64_t>>> mPoly;  // polyhedra
    std::vector<std::size_t> mCells;  // the zone's cells, in order; empty = all of them

    std::size_t NumCells(std::size_t ZoneCells) const {
        return mCells.empty() ? ZoneCells : mCells.size();
    }
    std::size_t Cell(std::size_t Row) const { return mCells.empty() ? Row : mCells[Row]; }
};

/// The ring of a polygonal element from its edges, directed so the element
/// lies on their left (Tecplot's convention: walking a face from its first
/// node to its second, the left element is on the left). The ring starts at
/// the first edge and follows it; when most edges disagree with that
/// direction (an inconsistently wound file) it is reversed. Empty when the
/// edges do not close one loop.
std::vector<std::int64_t> tecplot_polygon_ring(
    const std::vector<std::pair<std::int64_t, std::int64_t>>& rEdges) {
    if (rEdges.size() < 3)
        return {};
    std::map<std::int64_t, std::vector<std::int64_t>> adj;
    std::set<std::pair<std::int64_t, std::int64_t>> directed;
    for (const auto& [a, b] : rEdges) {
        adj[a].push_back(b);
        adj[b].push_back(a);
        directed.insert({a, b});
    }
    for (const auto& [node, nbrs] : adj)
        if (nbrs.size() != 2)
            return {};
    const std::int64_t start = rEdges[0].first;
    std::vector<std::int64_t> ring = {start, rEdges[0].second};
    while (ring.size() < rEdges.size()) {
        const auto& nb = adj[ring.back()];
        const std::int64_t next = nb[0] != ring[ring.size() - 2] ? nb[0] : nb[1];
        if (next == start)
            return {};
        ring.push_back(next);
    }
    const auto& last = adj[ring.back()];
    if (last[0] != start && last[1] != start)
        return {};
    std::size_t agree = 0;
    for (std::size_t i = 0; i < ring.size(); ++i)
        agree += directed.count({ring[i], ring[(i + 1) % ring.size()]});
    if (2 * agree < ring.size())
        std::reverse(ring.begin() + 1, ring.end());
    return ring;
}

/// A face-based zone's cells from its face map. A polyhedral face is wound
/// with its right-hand normal pointing at the right element, so it is kept
/// as-is for the left element (outward) and reversed for the right one.
std::vector<TecplotPiece> tecplot_face_pieces(const TecplotZone& rZ, std::size_t ZoneIdx,
                                              const TecplotFaceMap& rFm) {
    const std::string where = "Tecplot: zone " + std::to_string(ZoneIdx + 1);
    const std::size_t ncells = rZ.mNumCells;
    std::vector<TecplotPiece> pieces;
    if (!rZ.IsPolyhedron()) {
        std::vector<std::vector<std::pair<std::int64_t, std::int64_t>>> edges(ncells);
        for (std::size_t f = 0; f + 1 < rFm.mStart.size(); ++f) {
            const std::int64_t a = rFm.mNodes[static_cast<std::size_t>(rFm.mStart[f])];
            const std::int64_t b = rFm.mNodes[static_cast<std::size_t>(rFm.mStart[f]) + 1];
            if (rFm.mLeft[f] >= 0)
                edges[static_cast<std::size_t>(rFm.mLeft[f])].push_back({a, b});
            if (rFm.mRight[f] >= 0)
                edges[static_cast<std::size_t>(rFm.mRight[f])].push_back({b, a});
        }
        TecplotPiece piece;
        piece.mType = "polygon";
        piece.mRows.reserve(ncells);
        for (std::size_t c = 0; c < ncells; ++c) {
            std::vector<std::int64_t> ring = tecplot_polygon_ring(edges[c]);
            if (ring.empty())
                throw ReadError(where + ": element " + std::to_string(c + 1) +
                                " is not one closed polygon");
            piece.mRows.push_back(std::move(ring));
        }
        pieces.push_back(std::move(piece));
        return pieces;
    }
    std::vector<std::vector<std::vector<std::int64_t>>> faces(ncells);
    for (std::size_t f = 0; f + 1 < rFm.mStart.size(); ++f) {
        const auto b = rFm.mNodes.begin() + static_cast<std::ptrdiff_t>(rFm.mStart[f]);
        const auto e = rFm.mNodes.begin() + static_cast<std::ptrdiff_t>(rFm.mStart[f + 1]);
        if (rFm.mLeft[f] >= 0)
            faces[static_cast<std::size_t>(rFm.mLeft[f])].emplace_back(b, e);
        if (rFm.mRight[f] >= 0)
            faces[static_cast<std::size_t>(rFm.mRight[f])].emplace_back(
                std::make_reverse_iterator(e), std::make_reverse_iterator(b));
    }
    std::map<std::size_t, std::size_t> piece_of;  // node count -> piece
    for (std::size_t c = 0; c < ncells; ++c) {
        if (faces[c].size() < 4)
            throw ReadError(where + ": element " + std::to_string(c + 1) + " has " +
                            std::to_string(faces[c].size()) + " faces");
        std::set<std::int64_t> distinct;
        for (const auto& fc : faces[c])
            distinct.insert(fc.begin(), fc.end());
        const auto [it, fresh] = piece_of.emplace(distinct.size(), pieces.size());
        if (fresh) {
            pieces.emplace_back();
            pieces.back().mType = "polyhedron" + std::to_string(distinct.size());
        }
        pieces[it->second].mPoly.push_back(std::move(faces[c]));
        pieces[it->second].mCells.push_back(c);
    }
    if (pieces.size() == 1)
        pieces[0].mCells.clear();  // every cell, in order
    return pieces;
}

/// One zone's resolved data columns, connectivity and topology -- shared by
/// every step that references it (a zone can be VARSHARELIST'd or
/// CONNECTIVITYSHAREZONE'd from more than one later zone), so it is decoded
/// once and cached.
struct TecplotDecodedZone {
    std::vector<std::vector<double>> mCols;  // per variable; length NumNodes or NumCells
    std::vector<int> mCellCentered;
    std::string mMeshioType;
    NDArray mConn;                      // (NumCells, NodesPerCell) int64, 0-based
    TecplotFaceMap mFaces;              // face-based zones
    std::vector<TecplotPiece> mPieces;  // face-based zones; else one block from mConn
    std::size_t mNumNodes = 0;
    std::size_t mNumCells = 0;
};

class TecplotDecoder {
public:
    TecplotDecoder(const std::vector<TecplotZone>& rZones, std::size_t NumVariables,
                   const TecplotSource& rSource)
        : mrZones(rZones), mNumVariables(NumVariables), mrSource(rSource) {}

    /// Decodes zone `ZoneIdx`, resolving VARSHARELIST/CONNECTIVITYSHAREZONE by
    /// recursively decoding (and caching) whichever zone they name.
    const TecplotDecodedZone& Zone(std::size_t ZoneIdx, std::size_t Depth = 0) {
        const auto found = mCache.find(ZoneIdx);
        if (found != mCache.end())
            return found->second;
        if (Depth > mrZones.size())
            throw ReadError("Tecplot: circular VARSHARELIST/CONNECTIVITYSHAREZONE");
        const TecplotZone& z = mrZones[ZoneIdx];
        TecplotDecodedZone result;
        result.mMeshioType = tecplot_zone_meshio_type(z);
        result.mCols.resize(mNumVariables);
        mrSource.OwnData(ZoneIdx, result.mCols, result.mConn, result.mFaces);
        for (std::size_t k = 0; k < mNumVariables; ++k) {
            const auto share = z.mVarShareZone.find(k);
            if (share != z.mVarShareZone.end()) {
                CheckSource(share->second, ZoneIdx);
                result.mCols[k] = Zone(share->second, Depth + 1).mCols[k];
            } else if (z.mPassiveVars.count(k)) {
                result.mCols[k].assign(z.DataLength(k), std::numeric_limits<double>::quiet_NaN());
            }
        }
        if (z.mOrdered) {
            result.mConn = tecplot_ordered_connectivity(z);
        } else if (z.mConnShareZone >= 0) {
            const auto src = static_cast<std::size_t>(z.mConnShareZone);
            CheckSource(src, ZoneIdx);
            const TecplotDecodedZone& from = Zone(src, Depth + 1);
            if (mrZones[src].mTypeName != z.mTypeName || mrZones[src].mNumCells != z.mNumCells ||
                (z.IsPoly() && mrZones[src].mNumNodes != z.mNumNodes))
                throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                                " shares the connectivity of a different zone type or size");
            result.mConn = from.mConn;
            result.mFaces = from.mFaces;
        }
        if (z.IsPoly())
            result.mPieces = tecplot_face_pieces(z, ZoneIdx, result.mFaces);
        result.mCellCentered = z.mCellCentered;
        result.mNumNodes = z.mNumNodes;
        result.mNumCells = z.mNumCells;
        const auto [inserted, _] = mCache.emplace(ZoneIdx, std::move(result));
        return inserted->second;
    }

private:
    void CheckSource(std::size_t Src, std::size_t ZoneIdx) const {
        if (Src >= mrZones.size() || Src == ZoneIdx)
            throw ReadError("Tecplot: zone " + std::to_string(ZoneIdx + 1) +
                            " shares from bad zone " + std::to_string(Src + 1));
    }

    const std::vector<TecplotZone>& mrZones;
    std::size_t mNumVariables;
    const TecplotSource& mrSource;
    std::map<std::size_t, TecplotDecodedZone> mCache;
};

/// The zone whose points zone `ZoneIdx` literally reuses: itself, unless every
/// one of its nodal variables (coordinates and fields alike) is shared from
/// one earlier zone, in which case the chain is followed to that zone's own
/// origin. A zone that shares only its coordinates but carries its own nodal
/// values keeps its own copy of the points, so no value is lost. Tecplot
/// itself requires VARSHARELIST partners to agree on NODES, so this only
/// matters for zones that concatenate several element-type zones over one
/// shared point cloud (a common hybrid-mesh export shape), not for the
/// zones-as-a-timeline case, where sharing spans different steps.
std::size_t tecplot_points_origin_zone(std::size_t ZoneIdx, const std::vector<TecplotZone>& rZones,
                                       int Xi) {
    std::set<std::size_t> seen;
    std::size_t cur = ZoneIdx;
    while (seen.insert(cur).second) {
        const TecplotZone& z = rZones[cur];
        const auto it_x = z.mVarShareZone.find(static_cast<std::size_t>(Xi));
        if (it_x == z.mVarShareZone.end())
            return cur;
        for (std::size_t v = 0; v < z.mCellCentered.size(); ++v) {
            if (z.mCellCentered[v])
                continue;
            const auto it = z.mVarShareZone.find(v);
            if (it == z.mVarShareZone.end() || it->second != it_x->second)
                return cur;
        }
        cur = it_x->second;
        if (cur >= rZones.size())
            return ZoneIdx;
    }
    return ZoneIdx;  // cycle guard: fall back to owning its own points
}

/// The 0-based (X, Y, Z) variable indices, Z absent (-1) for a 2D file. The
/// same three variables for every zone in the file: VARIABLES is per-file.
void tecplot_xyz_indices(const std::vector<std::string>& rVariables, int& rXi, int& rYi, int& rZi) {
    rXi = rYi = rZi = -1;
    for (std::size_t k = 0; k < rVariables.size(); ++k) {
        const std::string v = tecplot_upper(rVariables[k]);
        if (v == "X")
            rXi = static_cast<int>(k);
        else if (v == "Y")
            rYi = static_cast<int>(k);
        else if (v == "Z")
            rZi = static_cast<int>(k);
    }
}

/// Per zone of a step: its point offset and whether it owns its points (a zone
/// reusing an earlier zone's points of the same step contributes none).
/// Returns the step's point count.
std::size_t tecplot_point_layout(const std::vector<std::size_t>& rZoneIdxs,
                                 const std::vector<TecplotZone>& rZones, int Xi,
                                 std::vector<std::size_t>& rOffset, std::vector<bool>& rOwns) {
    std::map<std::size_t, std::size_t> idx_to_pos;
    for (std::size_t k = 0; k < rZoneIdxs.size(); ++k)
        idx_to_pos.emplace(rZoneIdxs[k], k);
    rOffset.assign(rZoneIdxs.size(), 0);
    rOwns.assign(rZoneIdxs.size(), true);
    std::size_t total = 0;
    for (std::size_t k = 0; k < rZoneIdxs.size(); ++k) {
        const std::size_t origin =
            Xi < 0 ? rZoneIdxs[k] : tecplot_points_origin_zone(rZoneIdxs[k], rZones, Xi);
        const auto it = idx_to_pos.find(origin);
        if (origin != rZoneIdxs[k] && it != idx_to_pos.end() && it->second < k) {
            rOffset[k] = rOffset[it->second];
            rOwns[k] = false;
        } else {
            rOffset[k] = total;
            total += rZones[rZoneIdxs[k]].mNumNodes;
        }
    }
    return total;
}

/// Builds one step's Mesh by concatenating every zone in `rZoneIdxs`: one
/// cell block per zone (never merged by type -- each zone is its own named
/// part), points offset per zone (zones are not welded unless one reuses
/// another's points outright), a variable as cell data where a zone stores it
/// cell-centred and as point data where it stores it at the nodes (both, NaN
/// where absent, when zones disagree), `tecplot:zone` naming each cell's
/// zone, and one Cell region per zone (named by its own title when it has
/// one, de-duplicated, else `zone_<i>`).
Mesh tecplot_build_step_mesh(const std::vector<std::size_t>& rZoneIdxs,
                             const std::vector<TecplotZone>& rZones,
                             const std::vector<std::string>& rVariables,
                             const TecplotSource& rSource) {
    int xi = -1, yi = -1, zi = -1;
    tecplot_xyz_indices(rVariables, xi, yi, zi);
    if (xi < 0)
        throw ReadError("Tecplot: variable 'X' not found");
    if (yi < 0)
        throw ReadError("Tecplot: variable 'Y' not found");
    const std::size_t ndim = (zi >= 0) ? 3 : 2;

    TecplotDecoder decoder(rZones, rVariables.size(), rSource);
    std::vector<const TecplotDecodedZone*> decoded;
    decoded.reserve(rZoneIdxs.size());
    for (std::size_t idx : rZoneIdxs)
        decoded.push_back(&decoder.Zone(idx));

    std::vector<std::size_t> point_offset;
    std::vector<bool> owns_points;
    const std::size_t total_points =
        tecplot_point_layout(rZoneIdxs, rZones, xi, point_offset, owns_points);

    Mesh mesh;
    NDArray pts(DType::Float64, {total_points, ndim});
    double* pp = pts.As<double>();
    for (std::size_t k = 0; k < decoded.size(); ++k) {
        if (!owns_points[k])
            continue;
        const TecplotDecodedZone& d = *decoded[k];
        const std::size_t poff = point_offset[k];
        for (std::size_t r = 0; r < d.mNumNodes; ++r) {
            pp[(poff + r) * ndim + 0] = d.mCols[static_cast<std::size_t>(xi)][r];
            pp[(poff + r) * ndim + 1] = d.mCols[static_cast<std::size_t>(yi)][r];
            if (zi >= 0)
                pp[(poff + r) * ndim + 2] = d.mCols[static_cast<std::size_t>(zi)][r];
        }
    }
    mesh.AssignPoints(std::move(pts));

    // The step's cell blocks: one per zone, or one per piece of a face-based
    // zone, in zone order.
    struct BlockRef {
        std::size_t mZone;            // position in rZoneIdxs
        const TecplotPiece* mpPiece;  // nullptr: the zone's whole mConn
        std::size_t mNumCells;
    };
    std::vector<BlockRef> refs;
    for (std::size_t k = 0; k < decoded.size(); ++k) {
        const TecplotDecodedZone& d = *decoded[k];
        if (d.mPieces.empty()) {
            refs.push_back({k, nullptr, d.mNumCells});
            continue;
        }
        for (const TecplotPiece& piece : d.mPieces)
            refs.push_back({k, &piece, piece.NumCells(d.mNumCells)});
    }
    auto zone_cell = [](const BlockRef& rRef, std::size_t Row) {
        return rRef.mpPiece ? rRef.mpPiece->Cell(Row) : Row;
    };

    for (const BlockRef& ref : refs) {
        const TecplotDecodedZone& d = *decoded[ref.mZone];
        const auto off = static_cast<std::int64_t>(point_offset[ref.mZone]);
        if (!ref.mpPiece) {
            NDArray conn = d.mConn;  // deep copy: offset in place below
            std::int64_t* cp = conn.As<std::int64_t>();
            for (std::size_t j = 0; j < conn.Size(); ++j)
                cp[j] += off;
            mesh.AddCellBlock(d.mMeshioType, std::move(conn));
        } else if (!ref.mpPiece->mRows.empty()) {
            std::vector<std::vector<std::int64_t>> rows = ref.mpPiece->mRows;
            for (auto& row : rows)
                for (std::int64_t& v : row)
                    v += off;
            mesh.AddPolygonBlock(ref.mpPiece->mType, std::move(rows));
        } else {
            std::vector<std::vector<std::vector<std::int64_t>>> cells = ref.mpPiece->mPoly;
            for (auto& cell : cells)
                for (auto& face : cell)
                    for (std::int64_t& v : face)
                        v += off;
            mesh.AddPolyhedronBlock(ref.mpPiece->mType, std::move(cells));
        }
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < rVariables.size(); ++k) {
        if (static_cast<int>(k) == xi || static_cast<int>(k) == yi || static_cast<int>(k) == zi)
            continue;
        bool any_cc = false, all_cc = true;
        for (const TecplotDecodedZone* d : decoded) {
            any_cc = any_cc || d->mCellCentered[k];
            all_cc = all_cc && d->mCellCentered[k];
        }
        if (any_cc) {
            std::vector<NDArray> blk;
            blk.reserve(refs.size());
            for (const BlockRef& ref : refs) {
                const TecplotDecodedZone& d = *decoded[ref.mZone];
                NDArray arr(DType::Float64, {ref.mNumCells});
                double* ap = arr.As<double>();
                if (d.mCellCentered[k])
                    for (std::size_t r = 0; r < ref.mNumCells; ++r)
                        ap[r] = d.mCols[k][zone_cell(ref, r)];
                else
                    std::fill(ap, ap + ref.mNumCells, nan);
                blk.push_back(std::move(arr));
            }
            mesh.AddCellData(rVariables[k], std::move(blk));
        }
        if (!all_cc) {
            NDArray arr(DType::Float64, {total_points});
            double* ap = arr.As<double>();
            std::fill(ap, ap + total_points, nan);
            for (std::size_t z = 0; z < decoded.size(); ++z)
                if (owns_points[z] && !decoded[z]->mCellCentered[k])
                    std::memcpy(ap + point_offset[z], decoded[z]->mCols[k].data(),
                                decoded[z]->mNumNodes * sizeof(double));
            mesh.AddPointData(rVariables[k], std::move(arr));
        }
    }

    std::vector<NDArray> zone_ids;
    zone_ids.reserve(refs.size());
    for (const BlockRef& ref : refs) {
        NDArray zn(DType::Int64, {ref.mNumCells});
        std::int64_t* zp = zn.As<std::int64_t>();
        std::fill(zp, zp + ref.mNumCells, static_cast<std::int64_t>(rZoneIdxs[ref.mZone]));
        zone_ids.push_back(std::move(zn));
    }
    mesh.AddCellData("tecplot:zone", std::move(zone_ids));

    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    std::vector<std::string> used_names;
    for (std::size_t k = 0; k < rZoneIdxs.size(); ++k) {
        std::string name = rZones[rZoneIdxs[k]].mTitle;
        if (name.empty())
            name = "zone_" + std::to_string(k);
        std::string unique = name;
        for (int suffix = 2;
             std::find(used_names.begin(), used_names.end(), unique) != used_names.end(); ++suffix)
            unique = name + "_" + std::to_string(suffix);
        used_names.push_back(unique);

        NDArray entries(DType::Int64, {decoded[k]->mNumCells});
        std::int64_t* ep = entries.As<std::int64_t>();
        std::size_t n = 0;
        for (std::size_t b = 0; b < refs.size(); ++b)
            if (refs[b].mZone == k)
                for (std::size_t r = 0; r < refs[b].mNumCells; ++r)
                    ep[n++] = detail::block_row_to_global(bases, b, static_cast<std::int64_t>(r));
        mesh.AddRegion(Region(unique, RegionKind::Cell, -1, static_cast<std::int64_t>(rZoneIdxs[k]),
                              std::move(entries)));
    }
    return mesh;
}

/// A parsed file: its variables, zones and data source, ASCII or binary.
struct TecplotFile {
    std::vector<std::string> mVariables;
    std::vector<TecplotZone> mZones;
    std::vector<std::string> mLines;           // ASCII
    std::unique_ptr<detail::FileSource> mRaw;  // binary
    std::unique_ptr<TecplotSource> mSource;
};

void tecplot_open(const std::string& rPath, const ReadOptions& rOptions, TecplotFile& rFile) {
    bool binary = false;
    {
        auto probe = detail::make_classic_ifstream(rPath, std::ios::binary);
        if (!probe)
            throw ReadError("Could not open file: " + rPath);
        char head[5] = {0, 0, 0, 0, 0};
        probe.read(head, 5);
        binary = tecplot_is_plt(head, static_cast<std::size_t>(probe.gcount()));
    }
    if (!binary) {
        auto in = detail::make_classic_ifstream(rPath);
        if (!in)
            throw ReadError("Could not open file: " + rPath);
        std::string l;
        while (std::getline(in, l)) {
            std::string s = tecplot_strip(l);
            if (s.empty() || s[0] == '#')
                continue;
            rFile.mLines.push_back(std::move(s));
        }
        rFile.mZones = tecplot_scan_zones(rFile.mLines, rFile.mVariables);
        rFile.mSource = std::make_unique<TecplotAsciiSource>(rFile.mLines, rFile.mZones,
                                                             rFile.mVariables.size());
        return;
    }
    rFile.mRaw = std::make_unique<detail::FileSource>(rPath, rOptions.mMmap);
    const char* data = rFile.mRaw->Data();
    const std::size_t size = rFile.mRaw->Size();
    const std::string version = size >= 8 ? std::string(data + 5, 3) : std::string();
    if (version != "112")
        throw ReadError("Tecplot .plt: version '" + version +
                        "' is not supported (only #!TDV112, written by Tecplot 360 2009 and "
                        "later, is read)");
    // The int32 1 after the magic fixes the byte order.
    bool big_endian = false;
    {
        detail::ByteCursor little(data, size, false, "Tecplot .plt");
        little.Seek(8);
        if (little.I32() != 1) {
            detail::ByteCursor big(data, size, true, "Tecplot .plt");
            big.Seek(8);
            if (big.I32() != 1)
                throw ReadError("Tecplot .plt: bad byte-order word");
            big_endian = true;
        }
    }
    detail::ByteCursor cur(data, size, big_endian, "Tecplot .plt");
    cur.Seek(12);
    rFile.mZones = tecplot_plt_header(cur, rFile.mVariables);
    tecplot_plt_scan_data(cur, rFile.mVariables, rFile.mZones);
    rFile.mSource = std::make_unique<TecplotPltSource>(data, size, big_endian, rFile.mZones);
}

}  // namespace

MeshMetadata read_tecplot_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    TecplotFile file;
    tecplot_open(rPath, rOptions, file);
    const std::vector<std::vector<std::size_t>> timeline = tecplot_timeline(file.mZones);

    MeshMetadata meta;
    meta.mFormat = "tecplot";
    const std::vector<std::size_t>& first_step = timeline[0];
    int xi = -1, yi = -1, zi = -1;
    tecplot_xyz_indices(file.mVariables, xi, yi, zi);
    std::vector<std::size_t> offset;
    std::vector<bool> owns;
    meta.mNumPoints = tecplot_point_layout(first_step, file.mZones, xi, offset, owns);
    TecplotDecoder decoder(file.mZones, file.mVariables.size(), *file.mSource);
    for (std::size_t idx : first_step) {
        const TecplotZone& z = file.mZones[idx];
        if (z.IsPolyhedron()) {  // its blocks depend on its cells' node counts
            const TecplotDecodedZone& d = decoder.Zone(idx);
            for (const TecplotPiece& piece : d.mPieces) {
                CellBlockInfo block;
                block.mType = piece.mType;
                block.mNumCells = piece.NumCells(d.mNumCells);
                meta.mCellBlocks.push_back(std::move(block));
            }
            continue;
        }
        CellBlockInfo block;
        block.mType = tecplot_zone_meshio_type(z);
        block.mNumCells = z.mNumCells;
        meta.mCellBlocks.push_back(std::move(block));
    }
    meta.mPointDim = 0;  // not knowable without decoding X/Y/Z columns
    if (file.mZones[0].mHasSolutionTime)
        for (const std::vector<std::size_t>& step : timeline)
            meta.mTimeValues.push_back(file.mZones[step[0]].mSolutionTime);
    return meta;
}

Mesh read_tecplot(const std::string& rPath, const ReadOptions& rOptions) {
    TecplotFile file;
    tecplot_open(rPath, rOptions, file);
    const std::vector<std::vector<std::size_t>> timeline = tecplot_timeline(file.mZones);
    const std::size_t step = rOptions.ResolveTimeStep(timeline.size());
    return tecplot_build_step_mesh(timeline[step], file.mZones, file.mVariables, *file.mSource);
}

Mesh read_tecplot(const std::string& rPath) {
    return read_tecplot(rPath, ReadOptions{});
}

namespace {

/// The face map of a polygon or polyhedron block for a face-based zone:
/// faces in first-seen order (cells in order, each cell's faces in order),
/// a face shared by two cells written once. A cell's face is outward, so
/// that cell is its left element (right-hand normal towards the right one);
/// the second cell to use it becomes the right element. Elements 1-based,
/// 0 = none. Twin of `_tecplot._face_map`.
struct TecplotWriteFaces {
    std::vector<std::vector<std::int64_t>> mFaces;
    std::vector<std::int64_t> mLeft;
    std::vector<std::int64_t> mRight;
};

template <class TCellView>
TecplotWriteFaces tecplot_write_faces(const TCellView& rCb) {
    TecplotWriteFaces out;
    std::map<std::vector<std::int64_t>, std::size_t> seen;
    auto add = [&](std::vector<std::int64_t> Face, std::int64_t Cell) {
        std::vector<std::int64_t> key = Face;
        std::sort(key.begin(), key.end());
        const auto it = seen.find(key);
        if (it != seen.end() && out.mRight[it->second] == 0) {
            out.mRight[it->second] = Cell;
            return;
        }
        seen[key] = out.mFaces.size();
        out.mFaces.push_back(std::move(Face));
        out.mLeft.push_back(Cell);
        out.mRight.push_back(0);
    };
    const std::size_t n = rCb.NumCells();
    if (rCb.IsPolyhedron()) {
        for (std::size_t c = 0; c < n; ++c)
            for (std::size_t f = 0; f < rCb.NumFaces(c); ++f) {
                const auto [ptr, size] = rCb.Face(c, f);
                add(std::vector<std::int64_t>(ptr, ptr + size), static_cast<std::int64_t>(c + 1));
            }
        return out;
    }
    const NDArray* pConn = rCb.IsRagged() ? nullptr : &rCb.Conn();
    const std::size_t k = rCb.IsRagged() ? 0 : rCb.NodesPerCell();
    for (std::size_t c = 0; c < n; ++c) {
        std::vector<std::int64_t> row;
        if (pConn) {
            for (std::size_t j = 0; j < k; ++j)
                row.push_back(detail::read_int(*pConn, c * k + j));
        } else {
            row.assign(rCb.Row(c), rCb.Row(c) + rCb.RowSize(c));
        }
        for (std::size_t j = 0; j < row.size(); ++j)
            add({row[j], row[(j + 1) % row.size()]}, static_cast<std::int64_t>(c + 1));
    }
    return out;
}

}  // namespace

// Writes one Tecplot ZONE per cell block (no single-type restriction). Zone
// 1 carries the coordinates and every nodal field; later zones reuse them
// through VARSHARELIST rather than duplicating the (unwelded, shared) point
// array. Cell-centred variables are per-block: a zone that owns no values
// for a given variable declares it PASSIVEVARLIST instead of writing zeros.
void write_tecplot(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::size_t> blocks;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const auto cb = rMesh.Cells(i);
        if (!meshio_to_tecplot(cb.Type()).empty())
            blocks.push_back(i);
        else
            log::warn("tecplot: skipping unsupported cell type '{}'", cb.Type());
    }
    if (blocks.empty())
        throw WriteError("Tecplot writer found no supported cell blocks");

    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t dim = rMesh.PointDim();
    const std::size_t num_nodes = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();

    // Shared (nodal) variables: X, Y, (Z), then every point-data variable.
    std::vector<std::string> variables = {"X", "Y"};
    if (dim == 3)
        variables.push_back("Z");
    std::vector<std::vector<double>> shared_data;
    auto push_point_col = [&](std::size_t comp) {
        std::vector<double> col(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            col[r] = detail::read_double(points, r * dim + comp);
        shared_data.push_back(std::move(col));
    };
    push_point_col(0);
    push_point_col(1);
    if (dim == 3)
        push_point_col(2);
    for (const auto& k : rMesh.PointDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        const NDArray& v = rMesh.PointData(k);
        std::size_t ncomp = v.Shape().size() >= 2 ? v.Shape()[1] : 1;
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            std::vector<double> col(num_nodes);
            for (std::size_t r = 0; r < num_nodes; ++r)
                col[r] = detail::read_double(v, r * ncomp + c);
            shared_data.push_back(std::move(col));
        }
    }
    const std::size_t num_shared = variables.size();

    // Cell-centred variables (component-expanded); `mKey`/`mComp`/`mNumComp`
    // let each zone look its own values up (or find them absent) later.
    struct CellVarComp {
        std::string mKey;
        std::size_t mComp;
        std::size_t mNumComp;
    };
    std::vector<CellVarComp> cell_vars;
    for (const auto& k : rMesh.CellDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        if (rMesh.CellDataNumBlocks(k) == 0)
            continue;
        std::size_t ncomp = 1;
        for (std::size_t b : blocks) {
            if (b < rMesh.CellDataNumBlocks(k)) {
                const NDArray& vv = rMesh.CellData(k, b);
                ncomp = vv.Shape().size() >= 2 ? vv.Shape()[1] : 1;
                break;
            }
        }
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            cell_vars.push_back({k, c, ncomp});
        }
    }

    os << "TITLE = \"" << detail::provenance_lines(detail::SlotTier::SingleLine)[0] << "\"\n";
    os << "VARIABLES = ";
    for (std::size_t k = 0; k < variables.size(); ++k)
        os << (k ? ", " : "") << "\"" << variables[k] << "\"";
    os << "\n";

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    char buf[40];
    auto write_column = [&](const std::vector<double>& col) {
        for (std::size_t i = 0; i < col.size(); ++i) {
            detail::snprintf_c(buf, sizeof(buf), "%.17g", col[i]);
            os << buf << ((i + 1) % 20 == 0 || i + 1 == col.size() ? '\n' : ' ');
        }
        if (col.empty())
            os << "\n";
    };

    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        const std::size_t block = blocks[bi];
        const auto cb = rMesh.Cells(block);
        const std::string ztype = meshio_to_tecplot(cb.Type());
        const std::vector<int>& order = tecplot_order(cb.Type());
        const std::size_t num_cells = cb.NumCells();

        // Zone title: an exactly-matching Cell region's name, else block_<i>.
        std::string title = "block_" + std::to_string(block);
        for (std::size_t ri = 0; ri < rMesh.NumRegions(); ++ri) {
            const meshioplusplus::Region& reg = rMesh.Region(ri);
            if (reg.mKind != RegionKind::Cell || reg.NumEntries() == 0)
                continue;
            if (reg.NumEntries() !=
                static_cast<std::size_t>(bases[block + 1] - bases[block]))
                continue;
            const std::int64_t* e = reg.Entries();
            if (e[0] == bases[block] && e[reg.NumEntries() - 1] == bases[block + 1] - 1) {
                title = reg.mName;
                break;
            }
        }

        // Which cell variables this zone owns vs. leaves passive.
        std::vector<std::size_t> present, passive;
        for (std::size_t j = 0; j < cell_vars.size(); ++j) {
            if (block < rMesh.CellDataNumBlocks(cell_vars[j].mKey))
                present.push_back(j);
            else
                passive.push_back(j);
        }

        const bool face_based = ztype == "FEPOLYGON" || ztype == "FEPOLYHEDRON";
        TecplotWriteFaces faces;
        std::size_t total_face_nodes = 0;
        if (face_based) {
            faces = tecplot_write_faces(cb);
            for (const auto& f : faces.mFaces)
                total_face_nodes += f.size();
        }
        os << "ZONE T = \"" << title << "\", NODES = " << num_nodes
           << ", ELEMENTS = " << num_cells << ",\n";
        if (face_based)
            os << "FACES = " << faces.mFaces.size() << ", TOTALNUMFACENODES = " << total_face_nodes
               << ",\nNUMCONNECTEDBOUNDARYFACES = 0, TOTALNUMBOUNDARYCONNECTIONS = 0,\n";
        os << "DATAPACKING = BLOCK, ZONETYPE = " << ztype;
        if (bi > 0)
            os << ",\nVARSHARELIST = ([1-" << num_shared << "] = 1)";
        if (!present.empty()) {
            os << ",\nVARLOCATION = ([";
            for (std::size_t i = 0; i < present.size(); ++i)
                os << (i ? "," : "") << (num_shared + present[i] + 1);
            os << "] = CELLCENTERED)";
        }
        if (!passive.empty()) {
            os << ",\nPASSIVEVARLIST = ([";
            for (std::size_t i = 0; i < passive.size(); ++i)
                os << (i ? "," : "") << (num_shared + passive[i] + 1);
            os << "])";
        }
        os << "\n";

        if (bi == 0)
            for (const auto& col : shared_data)
                write_column(col);
        for (std::size_t j : present) {
            const CellVarComp& cv = cell_vars[j];
            const NDArray& vv = rMesh.CellData(cv.mKey, block);
            std::vector<double> col(vv.Shape()[0]);
            for (std::size_t r = 0; r < vv.Shape()[0]; ++r)
                col[r] = detail::read_double(vv, r * cv.mNumComp + cv.mComp);
            write_column(col);
        }

        if (face_based) {
            auto write_ints = [&](const std::vector<std::int64_t>& rInts) {
                for (std::size_t i = 0; i < rInts.size(); ++i)
                    os << rInts[i] << ((i + 1) % 20 == 0 || i + 1 == rInts.size() ? '\n' : ' ');
                if (rInts.empty())
                    os << "\n";
            };
            if (ztype == "FEPOLYHEDRON") {
                std::vector<std::int64_t> counts;
                for (const auto& f : faces.mFaces)
                    counts.push_back(static_cast<std::int64_t>(f.size()));
                write_ints(counts);
            }
            for (const auto& f : faces.mFaces)
                for (std::size_t j = 0; j < f.size(); ++j)
                    os << (f[j] + 1) << (j + 1 == f.size() ? '\n' : ' ');
            write_ints(faces.mLeft);
            write_ints(faces.mRight);
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < num_cells; ++r) {
            for (std::size_t j = 0; j < order.size(); ++j) {
                std::size_t src = static_cast<std::size_t>(order[j]);
                if (src >= k)
                    src = k - 1;
                os << (detail::read_int(conn, r * k + src) + 1)
                   << (j + 1 == order.size() ? '\n' : ' ');
            }
        }
    }
}

}  // namespace meshioplusplus
