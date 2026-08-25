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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

// (abaqus type, meshio type) in source order; the meshio->abaqus inverse keeps
// the last entry per meshio type (matching the Python dict comprehension).
const std::vector<std::pair<std::string, std::string>>& type_table() {
    static const std::vector<std::pair<std::string, std::string>> t = {
        {"T2D2", "line"},
        {"T2D2H", "line"},
        {"T2D3", "line3"},
        {"T2D3H", "line3"},
        {"T3D2", "line"},
        {"T3D2H", "line"},
        {"T3D3", "line3"},
        {"T3D3H", "line3"},
        {"B21", "line"},
        {"B21H", "line"},
        {"B22", "line3"},
        {"B22H", "line3"},
        {"B31", "line"},
        {"B31H", "line"},
        {"B32", "line3"},
        {"B32H", "line3"},
        {"B33", "line3"},
        {"B33H", "line3"},
        {"CPS4", "quad"},
        {"CPS4R", "quad"},
        {"S4", "quad"},
        {"S4R", "quad"},
        {"S4RS", "quad"},
        {"S4RSW", "quad"},
        {"S4R5", "quad"},
        {"S8R", "quad8"},
        {"S8R5", "quad8"},
        {"S9R5", "quad9"},
        {"CPS3", "triangle"},
        {"STRI3", "triangle"},
        {"S3", "triangle"},
        {"S3R", "triangle"},
        {"S3RS", "triangle"},
        {"R3D3", "triangle"},
        {"STRI65", "triangle6"},
        {"C3D8", "hexahedron"},
        {"C3D8H", "hexahedron"},
        {"C3D8I", "hexahedron"},
        {"C3D8IH", "hexahedron"},
        {"C3D8R", "hexahedron"},
        {"C3D8RH", "hexahedron"},
        {"C3D20", "hexahedron20"},
        {"C3D20H", "hexahedron20"},
        {"C3D20R", "hexahedron20"},
        {"C3D20RH", "hexahedron20"},
        {"C3D4", "tetra"},
        {"C3D4H", "tetra4"},
        {"C3D10", "tetra10"},
        {"C3D10H", "tetra10"},
        {"C3D10I", "tetra10"},
        {"C3D10M", "tetra10"},
        {"C3D10MH", "tetra10"},
        {"C3D6", "wedge"},
        {"C3D15", "wedge15"},
        {"CAX4P", "quad"},
        {"CPE6", "triangle6"},
    };
    return t;
}

const std::unordered_map<std::string, std::string>& abaqus_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> r;
        for (const auto& kv : type_table())
            r[kv.first] = kv.second;
        return r;
    }();
    return m;
}

const std::unordered_map<std::string, std::string>& meshio_to_abaqus() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> r;
        for (const auto& kv : type_table())
            r[kv.second] = kv.first;  // last wins
        return r;
    }();
    return m;
}

std::string abaqus_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}
std::string abaqus_trim(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}
std::vector<std::string> split(const std::string& rS, char sep) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream iss(rS);
    while (std::getline(iss, cur, sep))
        out.push_back(abaqus_trim(cur));
    return out;
}

// --- parameter parsing --------------------------------------------------------

// The comma-separated `KEY` / `KEY=VALUE` parameters on a `*KEYWORD` line, keys
// upper-cased (the Python reference's `get_param_map`). A bare key maps to an
// empty value, which is how `GENERATE` is detected.
//
// The leading `*` on the keyword token is deliberately **kept**, exactly as the
// Python reference does: `*ELSET, ELSET=solid` would otherwise put the bare
// keyword and the real parameter under the same key, and the first one in wins.
std::unordered_map<std::string, std::string> abq_param_map(const std::string& rLine) {
    std::unordered_map<std::string, std::string> out;
    for (const std::string& word : split(rLine, ',')) {
        const std::size_t eq = word.find('=');
        if (eq == std::string::npos)
            out.insert_or_assign(abaqus_upper(abaqus_trim(word)), std::string());
        else
            out.insert_or_assign(abaqus_upper(abaqus_trim(word.substr(0, eq))),
                                 abaqus_trim(word.substr(eq + 1)));
    }
    return out;
}

/// A parameter's value, or an empty string when it is absent or bare.
std::string abq_param(const std::unordered_map<std::string, std::string>& rParams,
                      const std::string& rKey) {
    auto it = rParams.find(rKey);
    return it == rParams.end() ? std::string() : it->second;
}

// --- side-set face identifiers ------------------------------------------------

/**
 * @brief Abaqus face identifier (`S1`..`S6`) -> the meshio++ local facet index
 * of `detail/cell_faces.hpp` / `cell_edges.hpp`.
 *
 * The two numberings genuinely differ, so this table is not the identity. It is
 * derived by matching node sets: Abaqus C3D8 `S1` is the 1-2-3-4 face, i.e.
 * local nodes {0,1,2,3}, which is meshio++'s face 4 (`{0,3,2,1}`) — same face,
 * different slot and winding. Getting this wrong yields a plausible-looking
 * side set pointing at the wrong faces, so each row is spelled out.
 *
 * Shell elements use `SPOS`/`SNEG` for their two sides rather than a facet;
 * there is no facet to name, so they map to 0 and 1 respectively and are
 * documented as such in doc/regions.md.
 *
 * @param rCellType The meshio++ cell type of the element.
 * @param rFace The Abaqus face identifier, upper-cased (e.g. `"S3"`).
 * @return The local facet index, or -1 when the pair has no mapping.
 */
int abq_face_index(const std::string& rCellType, const std::string& rFace) {
    if (rFace == "SPOS")
        return 0;
    if (rFace == "SNEG")
        return 1;
    if (rFace.size() < 2 || rFace[0] != 'S')
        return -1;
    int n = 0;
    for (std::size_t k = 1; k < rFace.size(); ++k) {
        if (!std::isdigit(static_cast<unsigned char>(rFace[k])))
            return -1;
        n = n * 10 + (rFace[k] - '0');
    }
    if (n < 1)
        return -1;
    const int s = n - 1;  // 0-based Abaqus face number

    // tetra: Abaqus S1=1-2-3, S2=1-2-4, S3=2-3-4, S4=1-3-4
    //        meshio++ 0={0,1,3} 1={1,2,3} 2={2,0,3} 3={0,2,1}
    static const int tetra[4] = {3, 0, 1, 2};
    // hexahedron: Abaqus S1=1-2-3-4, S2=5-8-7-6, S3=1-5-6-2,
    //                    S4=2-6-7-3,  S5=3-7-8-4, S6=4-8-5-1
    //             meshio++ 0={0,4,7,3} 1={1,2,6,5} 2={0,1,5,4}
    //                      3={3,7,6,2} 4={0,3,2,1} 5={4,5,6,7}
    static const int hexa[6] = {4, 5, 2, 1, 3, 0};
    // wedge: Abaqus S1=1-2-3, S2=4-5-6, S3=1-2-5-4, S4=2-3-6-5, S5=3-1-4-6
    //        meshio++ 0={0,2,1} 1={3,4,5} 2={0,1,4,3} 3={1,2,5,4} 4={2,0,3,5}
    static const int wedge[5] = {0, 1, 2, 3, 4};
    // 2-D elements: Abaqus numbers the edges 1-2, 2-3, ... in the same order
    // detail/cell_edges.hpp does.
    static const int identity[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

    const int* table = nullptr;
    int count = 0;
    if (rCellType == "tetra" || rCellType == "tetra10") {
        table = tetra;
        count = 4;
    } else if (rCellType == "hexahedron" || rCellType == "hexahedron20") {
        table = hexa;
        count = 6;
    } else if (rCellType == "wedge" || rCellType == "wedge15") {
        table = wedge;
        count = 5;
    } else if (rCellType == "triangle" || rCellType == "triangle6") {
        table = identity;
        count = 3;
    } else if (rCellType == "quad" || rCellType == "quad8" || rCellType == "quad9") {
        table = identity;
        count = 4;
    }
    if (table == nullptr || s >= count)
        return -1;
    return table[s];
}

/// The inverse of `abq_face_index`, for the writer.
std::string abq_face_name(const std::string& rCellType, std::int64_t Facet) {
    for (int n = 1; n <= 6; ++n) {
        const std::string face = "S" + std::to_string(n);
        if (abq_face_index(rCellType, face) == static_cast<int>(Facet))
            return face;
    }
    return {};
}

// --- reader state -------------------------------------------------------------

/// One `*ELEMENT` block being accumulated.
struct AbqBlock {
    std::string mType;  // meshio type
    std::size_t mNodesPerCell = 0;
    std::vector<std::int64_t> mConn;        // row-major, 0-based point indices
    std::vector<std::int64_t> mElementIds;  // the file's element ids, in row order
};

/// Everything one `.inp` file (and everything it `*INCLUDE`s) contributes.
struct AbqFile {
    std::vector<std::vector<double>> mPoints;
    std::unordered_map<std::int64_t, std::int64_t> mPointIds;  // file id -> point index
    std::vector<AbqBlock> mBlocks;
    // Named groups, kept as *file ids* until the very end: an ELSET may name
    // another ELSET declared later, and a SURFACE may name an ELSET, so the
    // resolution to indices has to wait until everything is known.
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> mNodeSets;
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> mElemSets;
    // (name, [(element id or elset name, face identifier)])
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>> mSurfaces;
};

void abq_read_file(const std::string& rPath, AbqFile& rOut, int Depth);

/// The data lines following a keyword, up to the next `*` line.
std::vector<std::string> abq_data_lines(const std::vector<std::string>& rLines, std::size_t& rI) {
    std::vector<std::string> out;
    while (rI < rLines.size() && (rLines[rI].empty() || rLines[rI][0] != '*')) {
        const std::string row = abaqus_trim(rLines[rI]);
        ++rI;
        if (!row.empty())
            out.push_back(row);
    }
    return out;
}

/**
 * @brief Read a `*NSET` / `*ELSET` body: numeric ids and/or referenced set names.
 *
 * `GENERATE` turns a `first, last, step` triple into the explicit range, which
 * is what Abaqus means by it.
 */
void abq_read_set(const std::vector<std::string>& rRows, bool Generate,
                  std::vector<std::int64_t>& rIds, std::vector<std::string>& rNames) {
    for (const std::string& row : rRows) {
        for (const std::string& tok : split(row, ',')) {
            if (tok.empty())
                continue;
            const bool numeric =
                std::isdigit(static_cast<unsigned char>(tok[0])) || tok[0] == '-' || tok[0] == '+';
            if (numeric)
                rIds.push_back(std::strtoll(tok.c_str(), nullptr, 10));
            else
                rNames.push_back(tok);
        }
    }
    if (Generate) {
        if (rIds.size() < 3)
            throw ReadError("Abaqus: GENERATE needs first, last, step");
        const std::int64_t first = rIds[0], last = rIds[1];
        const std::int64_t step = rIds[2] == 0 ? 1 : rIds[2];
        std::vector<std::int64_t> gen;
        for (std::int64_t v = first; step > 0 ? v <= last : v >= last; v += step)
            gen.push_back(v);
        rIds = std::move(gen);
    }
}

void abq_read_lines(const std::vector<std::string>& rLines, const std::string& rPath, AbqFile& rOut,
                    int Depth) {
    const auto& a2m = abaqus_to_meshio();
    std::size_t i = 0;
    while (i < rLines.size()) {
        const std::string& line = rLines[i];
        if (line.rfind("**", 0) == 0) {  // comment
            ++i;
            continue;
        }
        std::string kw = abaqus_upper(abaqus_trim(split(line, ',')[0]));
        if (!kw.empty() && kw[0] == '*')
            kw = kw.substr(1);
        const std::unordered_map<std::string, std::string> params = abq_param_map(line);

        if (kw == "NODE") {
            ++i;
            for (const std::string& row : abq_data_lines(rLines, i)) {
                const std::vector<std::string> tok = split(row, ',');
                const std::int64_t id = std::strtoll(tok[0].c_str(), nullptr, 10);
                rOut.mPointIds[id] = static_cast<std::int64_t>(rOut.mPoints.size());
                std::vector<double> c;
                for (std::size_t k = 1; k < tok.size(); ++k)
                    if (!tok[k].empty())
                        c.push_back(std::strtod(tok[k].c_str(), nullptr));
                rOut.mPoints.push_back(std::move(c));
            }
        } else if (kw == "ELEMENT") {
            const std::string etype = abq_param(params, "TYPE");
            if (etype.empty())
                throw ReadError("Abaqus ELEMENT without TYPE");
            auto it = a2m.find(abaqus_upper(etype));
            if (it == a2m.end())
                it = a2m.find(etype);  // types are case-sensitive in some files
            if (it == a2m.end())
                throw ReadError("Abaqus element type not supported: " + etype);
            const std::string mtype = it->second;
            const int n = num_nodes_per_cell().count(mtype) ? num_nodes_per_cell().at(mtype) : 0;
            if (n == 0)
                throw ReadError("Abaqus: unknown node count for " + mtype);

            ++i;
            std::vector<std::int64_t> vals;
            for (const std::string& row : abq_data_lines(rLines, i))
                for (const std::string& t : split(row, ','))
                    if (!t.empty())
                        vals.push_back(std::strtoll(t.c_str(), nullptr, 10));

            const std::size_t stride = static_cast<std::size_t>(n) + 1;
            if (vals.size() % stride != 0)
                throw ReadError("Abaqus: bad element data");
            AbqBlock block;
            block.mType = mtype;
            block.mNodesPerCell = static_cast<std::size_t>(n);
            const std::size_t ncells = vals.size() / stride;
            block.mConn.reserve(ncells * static_cast<std::size_t>(n));
            block.mElementIds.reserve(ncells);
            for (std::size_t r = 0; r < ncells; ++r) {
                block.mElementIds.push_back(vals[r * stride]);
                for (int j = 0; j < n; ++j) {
                    const std::int64_t node = vals[r * stride + 1 + j];
                    auto pit = rOut.mPointIds.find(node);
                    if (pit == rOut.mPointIds.end())
                        throw ReadError("Abaqus: unknown node id");
                    block.mConn.push_back(pit->second);
                }
            }
            // `*ELEMENT, ELSET=name` declares a set covering exactly this block.
            const std::string elset = abq_param(params, "ELSET");
            if (!elset.empty())
                rOut.mElemSets.emplace_back(elset, block.mElementIds);
            rOut.mBlocks.push_back(std::move(block));
        } else if (kw == "NSET" || kw == "ELSET") {
            const bool is_node = (kw == "NSET");
            const std::string name = abq_param(params, is_node ? "NSET" : "ELSET");
            if (name.empty())
                throw ReadError("Abaqus " + kw + " without a name");
            ++i;
            const std::vector<std::string> rows = abq_data_lines(rLines, i);
            std::vector<std::int64_t> ids;
            std::vector<std::string> refs;
            abq_read_set(rows, params.count("GENERATE") > 0, ids, refs);
            // A set that names other sets is expanded once every set is known.
            for (const std::string& ref : refs) {
                const auto& src = is_node ? rOut.mNodeSets : rOut.mElemSets;
                bool found = false;
                for (const auto& kv : src)
                    if (kv.first == ref) {
                        ids.insert(ids.end(), kv.second.begin(), kv.second.end());
                        found = true;
                    }
                if (!found)
                    throw ReadError("Abaqus: unknown " + kw + " '" + ref + "'");
            }
            (is_node ? rOut.mNodeSets : rOut.mElemSets).emplace_back(name, std::move(ids));
        } else if (kw == "SURFACE") {
            // `*SURFACE, NAME=..., TYPE=ELEMENT`: each data row is
            // `<element id | elset name>, <face identifier>`.
            const std::string name = abq_param(params, "NAME");
            const std::string type = abaqus_upper(abq_param(params, "TYPE"));
            ++i;
            const std::vector<std::string> rows = abq_data_lines(rLines, i);
            if (name.empty() || (!type.empty() && type != "ELEMENT"))
                continue;  // node-based surfaces have no facets: skip, don't fail
            std::vector<std::pair<std::string, std::string>> members;
            for (const std::string& row : rows) {
                const std::vector<std::string> tok = split(row, ',');
                if (tok.size() < 2 || tok[0].empty() || tok[1].empty())
                    continue;
                members.emplace_back(tok[0], abaqus_upper(tok[1]));
            }
            rOut.mSurfaces.emplace_back(name, std::move(members));
        } else if (kw == "INCLUDE") {
            if (Depth > 8)
                throw ReadError("Abaqus: *INCLUDE nested too deeply");
            const std::size_t eq = line.rfind('=');
            if (eq == std::string::npos)
                throw ReadError("Abaqus: *INCLUDE without INPUT=");
            std::string inc = abaqus_trim(line.substr(eq + 1));
            std::filesystem::path p(inc);
            if (!std::filesystem::exists(p))
                p = std::filesystem::path(rPath).parent_path() / inc;
            abq_read_file(p.string(), rOut, Depth + 1);
            ++i;
        } else {
            // There are far too many Abaqus keywords to enumerate; skip the
            // keyword line and its data lines.
            ++i;
            abq_data_lines(rLines, i);
        }
    }
}

void abq_read_file(const std::string& rPath, AbqFile& rOut, int Depth) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        if (!l.empty() && l.back() == '\r')
            l.pop_back();
        lines.push_back(l);
    }
    abq_read_lines(lines, rPath, rOut, Depth);
}

}  // namespace

Mesh read_abaqus(const std::string& rPath) {
    AbqFile file;
    abq_read_file(rPath, file, 0);

    Mesh mesh;

    // --- points -------------------------------------------------------------
    std::size_t dim = 3;
    if (!file.mPoints.empty()) {
        dim = file.mPoints[0].size();
        if (dim == 0)
            dim = 3;
    }
    NDArray points(DType::Float64, {file.mPoints.size(), dim});
    double* pp = points.As<double>();
    for (std::size_t r = 0; r < file.mPoints.size(); ++r)
        for (std::size_t c = 0; c < dim; ++c)
            pp[r * dim + c] = (c < file.mPoints[r].size()) ? file.mPoints[r][c] : 0.0;
    mesh.AssignPoints(std::move(points));

    // --- cells, plus element id -> global (block-major) cell index -----------
    std::unordered_map<std::int64_t, std::int64_t> elem_index;
    std::vector<std::string> block_types;
    std::int64_t global = 0;
    for (const AbqBlock& block : file.mBlocks) {
        const std::size_t ncells = block.mElementIds.size();
        NDArray data(DType::Int64, {ncells, block.mNodesPerCell});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t k = 0; k < block.mConn.size(); ++k)
            dp[k] = block.mConn[k];
        mesh.AddCellBlock(block.mType, std::move(data));
        block_types.push_back(block.mType);
        for (std::size_t c = 0; c < ncells; ++c)
            elem_index[block.mElementIds[c]] = global + static_cast<std::int64_t>(c);
        global += static_cast<std::int64_t>(ncells);
    }

    // --- named groups -------------------------------------------------------
    // Ids are resolved to indices only now, because a set may reference another
    // set (or a surface an element set) declared later in the file.
    for (const auto& [name, ids] : file.mNodeSets) {
        std::vector<std::int64_t> idx;
        idx.reserve(ids.size());
        for (std::int64_t id : ids) {
            auto it = file.mPointIds.find(id);
            if (it != file.mPointIds.end())
                idx.push_back(it->second);
        }
        NDArray entries = NDArray::Uninit(DType::Int64, {idx.size()});
        for (std::size_t k = 0; k < idx.size(); ++k)
            entries.As<std::int64_t>()[k] = idx[k];
        mesh.AddRegion(Region(name, RegionKind::Point, std::move(entries)));
    }

    std::unordered_map<std::string, std::vector<std::int64_t>> elset_cells;
    for (const auto& [name, ids] : file.mElemSets) {
        std::vector<std::int64_t> idx;
        idx.reserve(ids.size());
        for (std::int64_t id : ids) {
            auto it = elem_index.find(id);
            if (it != elem_index.end())
                idx.push_back(it->second);
        }
        NDArray entries = NDArray::Uninit(DType::Int64, {idx.size()});
        for (std::size_t k = 0; k < idx.size(); ++k)
            entries.As<std::int64_t>()[k] = idx[k];
        elset_cells[name] = idx;
        mesh.AddRegion(Region(name, RegionKind::Cell, std::move(entries)));
    }

    // Block-major cell index -> its block's meshio type, for the face mapping.
    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    auto type_of = [&](std::int64_t g) -> std::string {
        const auto [b, row] = detail::global_to_block_row(bases, g);
        (void)row;
        return b == static_cast<std::size_t>(-1) ? std::string() : block_types[b];
    };

    for (const auto& [name, members] : file.mSurfaces) {
        std::vector<std::int64_t> pairs;
        for (const auto& [who, face] : members) {
            std::vector<std::int64_t> cells;
            auto sit = elset_cells.find(who);
            if (sit != elset_cells.end()) {
                cells = sit->second;
            } else {
                const std::int64_t id = std::strtoll(who.c_str(), nullptr, 10);
                auto eit = elem_index.find(id);
                if (eit != elem_index.end())
                    cells.push_back(eit->second);
            }
            for (std::int64_t g : cells) {
                const int facet = abq_face_index(type_of(g), face);
                if (facet < 0)
                    continue;  // an identifier this element type has no facet for
                pairs.push_back(g);
                pairs.push_back(facet);
            }
        }
        NDArray entries = NDArray::Uninit(DType::Int64, {pairs.size() / 2, 2});
        for (std::size_t k = 0; k < pairs.size(); ++k)
            entries.As<std::int64_t>()[k] = pairs[k];
        mesh.AddRegion(Region(name, RegionKind::Side, std::move(entries)));
    }

    return mesh;
}

void write_abaqus(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t n = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;

    os << "*HEADING\n";
    os << "Abaqus DataFile Version 6.14\n";
    os << detail::provenance_render_lines(detail::SlotTier::Block, "");
    os << "*NODE\n";
    {
        // Format node rows in parallel (snprintf per row, bytes unchanged),
        // then stream sequentially.
        std::vector<std::string> rows(n);
        parallel_for(n, [&](std::size_t i) {
            char buf[48];
            std::string& row = rows[i];
            row = std::to_string(i + 1);
            for (std::size_t c = 0; c < dim; ++c) {
                std::snprintf(buf, sizeof(buf), ", %.16e",
                              detail::read_double(points, i * dim + c));
                row += buf;
            }
            row += '\n';
        });
        for (const auto& row : rows)
            os << row;
    }

    const auto& m2a = meshio_to_abaqus();
    std::size_t eid = 0;
    for (const auto cb : rMesh.CellRange()) {
        auto it = m2a.find(cb.Type());
        if (it == m2a.end())
            throw WriteError("Abaqus writer: unsupported cell type " + cb.Type());
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        os << "*ELEMENT, TYPE=" << it->second << "\n";
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            os << (++eid);
            for (std::size_t j = 0; j < k; ++j)
                os << "," << (detail::read_int(conn, r * k + j) + 1);
            os << "\n";
        }
    }

    // --- named groups -------------------------------------------------------
    // Ids are 1-based on the way out, and wrapped at 8 per line, matching the
    // Python reference writer so a set-carrying mesh produces the same file
    // whichever path served it.
    constexpr std::size_t per_line = 8;
    auto write_ids = [&](const std::vector<std::int64_t>& rIds) {
        for (std::size_t k = 0; k < rIds.size(); ++k) {
            os << (rIds[k] + 1);
            const bool last = (k + 1 == rIds.size());
            os << (last ? "\n" : ((k + 1) % per_line == 0 ? ",\n" : ","));
        }
    };

    // An empty group is still declared: the *name* is information, and it is
    // what a round-trip through the Python reference writer has always
    // produced (its per-block emission left the key behind even when nothing
    // landed in it). Same rule as detail/region_remap.hpp's empty-region carry.
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        os << "*ELSET, ELSET=" << r.mName << "\n";
        write_ids(std::vector<std::int64_t>(r.Entries(), r.Entries() + r.NumEntries()));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Point)
            continue;
        os << "*NSET, NSET=" << r.mName << "\n";
        write_ids(std::vector<std::int64_t>(r.Entries(), r.Entries() + r.NumEntries()));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Side)
            continue;
        os << "*SURFACE, NAME=" << r.mName << ", TYPE=ELEMENT\n";
        const std::int64_t* e = r.Entries();
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            const auto [b, row] = detail::global_to_block_row(bases, e[k * 2]);
            (void)row;
            if (b == static_cast<std::size_t>(-1))
                continue;
            const std::string face =
                abq_face_name(std::string(rMesh.Cells(b).Type()), e[k * 2 + 1]);
            if (face.empty())
                continue;  // no Abaqus identifier for this type/facet pair
            os << (e[k * 2] + 1) << ", " << face << "\n";
        }
    }
}

}  // namespace meshioplusplus
