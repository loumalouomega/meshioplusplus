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
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/fast_number.hpp"

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
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}
bool is_float_token(const std::string& rS) {
    if (rS.empty())
        return false;
    char* endp = nullptr;
    std::strtod(rS.c_str(), &endp);
    return endp == rS.c_str() + rS.size();
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

// One ZONE header's parsed fields, its data section's line range, and its
// transient identity (SOLUTIONTIME/STRANDID). Shared by read_tecplot (which
// decodes exactly one zone's data) and read_tecplot_metadata (which decodes
// none): both start from the same tecplot_scan_zones pass, so which zone the
// data-decoding step picks and which zones the metadata's timeline lists can
// never drift against each other.
struct TecplotZoneHeader {
    std::map<std::string, std::string> mFields;  // NODES/N/ELEMENTS/E/DATAPACKING/ZONETYPE/F/ET/NV
    std::string mVarloc;
    bool mHasSolutionTime = false;
    double mSolutionTime = 0.0;
    bool mHasStrandId = false;
    int mStrandId = 0;
    std::size_t mDataStart = 0;
    std::size_t mNumNodes = 0;
    std::size_t mNumCells = 0;
};

std::string tecplot_zone_field(const TecplotZoneHeader& rZ, const char* pA, const char* pB) {
    auto it = rZ.mFields.find(pA);
    if (it != rZ.mFields.end())
        return it->second;
    it = rZ.mFields.find(pB);
    return it != rZ.mFields.end() ? it->second : std::string();
}

// The zone's element type and, for FEBLOCK data, which variables are
// cell-centered -- the same derivation whether or not this zone's data is
// ever decoded.
void tecplot_zone_format(const TecplotZoneHeader& rZ, std::size_t NumVariables, bool& rFeblock,
                         std::string& rZtype, std::vector<int>& rCellCentered) {
    std::string fmt;
    if (rZ.mFields.count("F")) {
        fmt = tecplot_upper(rZ.mFields.at("F"));
        rZtype = rZ.mFields.count("ET") ? rZ.mFields.at("ET") : "";
    } else {
        fmt = "FE" + tecplot_upper(tecplot_zone_field(rZ, "DATAPACKING", ""));
        rZtype = tecplot_zone_field(rZ, "ZONETYPE", "");
    }
    rFeblock = (fmt == "FEBLOCK");

    rCellCentered.assign(NumVariables, 0);
    if (!rFeblock)
        return;
    if (rZ.mFields.count("NV")) {
        int nv = std::stoi(rZ.mFields.at("NV"));
        for (std::size_t k = static_cast<std::size_t>(nv); k < NumVariables; ++k)
            rCellCentered[k] = 1;
    } else if (!rZ.mVarloc.empty()) {
        std::string vc = rZ.mVarloc.substr(1, rZ.mVarloc.size() - 2);  // strip ()
        std::vector<std::string> entries;
        std::string cur;
        for (char c : vc) {
            if (c == ',') {
                entries.push_back(cur);
                cur.clear();
            } else {
                cur += c;
            }
        }
        if (!cur.empty())
            entries.push_back(cur);
        for (const auto& entry : entries) {
            std::size_t eq = entry.find('=');
            if (eq == std::string::npos)
                continue;
            std::string rng = entry.substr(0, eq), loc = tecplot_upper(entry.substr(eq + 1));
            if (loc != "CELLCENTERED")
                continue;
            rng = rng.substr(1, rng.size() - 2);  // strip []
            std::size_t dash = rng.find('-');
            if (dash == std::string::npos) {
                rCellCentered[static_cast<std::size_t>(std::stoi(rng) - 1)] = 1;
            } else {
                int a = std::stoi(rng.substr(0, dash)), b = std::stoi(rng.substr(dash + 1));
                for (int k = a; k <= b; ++k)
                    rCellCentered[static_cast<std::size_t>(k - 1)] = 1;
            }
        }
    }
}

// How many numeric tokens this zone's data block holds, in file order --
// FEBLOCK is one run per variable (cell-centered ones NumCells long, the
// rest NumNodes long); POINT/FEPOINT is NumNodes rows of NumVariables each.
std::size_t tecplot_zone_data_token_count(const TecplotZoneHeader& rZ, std::size_t NumVariables,
                                          bool Feblock, const std::vector<int>& rCellCentered) {
    if (!Feblock)
        return rZ.mNumNodes * NumVariables;
    std::size_t total = 0;
    for (std::size_t k = 0; k < NumVariables; ++k)
        total += rCellCentered[k] ? rZ.mNumCells : rZ.mNumNodes;
    return total;
}

// One pass over every line, splitting VARIABLES from every ZONE header (not
// just the first) and locating each zone's data section by actually counting
// off its own token budget plus its connectivity lines -- Tecplot ASCII has
// no fixed tokens-per-line convention, so this is the only reliable way to
// find where one zone's data ends and the next one's header begins.
std::vector<TecplotZoneHeader> tecplot_scan_zones(const std::vector<std::string>& rLines,
                                                  std::vector<std::string>& rVariables) {
    std::vector<TecplotZoneHeader> zones;
    std::size_t i = 0;
    for (; i < rLines.size(); ++i) {
        std::string u = tecplot_upper(rLines[i]);
        if (u.rfind("VARIABLES", 0) == 0) {
            std::string joined = rLines[i];
            while (i + 1 < rLines.size() && tecplot_strip(rLines[i + 1])[0] == '"')
                joined += " " + rLines[++i];
            std::string rhs = joined.substr(joined.find('=') + 1);
            std::size_t p = 0;
            while (p < rhs.size()) {
                if (rhs[p] == '"') {
                    std::size_t q = rhs.find('"', p + 1);
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
            continue;
        }
        if (u.rfind("ZONE", 0) != 0)
            continue;

        TecplotZoneHeader z;
        std::string joined = rLines[i];
        while (i + 1 < rLines.size() && !is_float_token(tecplot_tokens(rLines[i + 1])[0]))
            joined += " " + rLines[++i];
        z.mDataStart = i + 1;

        std::string ju = joined;
        std::size_t vp = tecplot_upper(ju).find("VARLOCATION");
        if (vp != std::string::npos) {
            std::size_t p1 = ju.find('(', vp), p2 = ju.find(')', p1);
            z.mVarloc = ju.substr(p1, p2 - p1 + 1);
            z.mVarloc.erase(std::remove(z.mVarloc.begin(), z.mVarloc.end(), ' '), z.mVarloc.end());
            ju = ju.substr(0, vp) + ju.substr(p2 + 1);
        }
        std::string body = ju.substr(4);
        for (auto& c : body)
            if (c == ',' || c == '=')
                c = ' ';
        auto tk = tecplot_tokens(body);
        for (std::size_t k = 0; k + 1 < tk.size(); ++k) {
            std::string key = tecplot_upper(tk[k]);
            if (key == "NODES" || key == "N" || key == "ELEMENTS" || key == "E" ||
                key == "DATAPACKING" || key == "ZONETYPE" || key == "F" || key == "ET" ||
                key == "NV") {
                z.mFields[key] = tk[k + 1];
            } else if (key == "SOLUTIONTIME" || key == "STRANDID") {
                if (key == "SOLUTIONTIME") {
                    z.mSolutionTime = std::strtod(tk[k + 1].c_str(), nullptr);
                    z.mHasSolutionTime = true;
                } else {
                    z.mStrandId = std::stoi(tk[k + 1]);
                    z.mHasStrandId = true;
                }
            }
        }
        z.mNumNodes = std::stoull(tecplot_zone_field(z, "NODES", "N"));
        z.mNumCells = std::stoull(tecplot_zone_field(z, "ELEMENTS", "E"));

        bool feblock = false;
        std::string ztype;
        std::vector<int> cell_centered;
        tecplot_zone_format(z, rVariables.size(), feblock, ztype, cell_centered);
        const std::size_t want =
            tecplot_zone_data_token_count(z, rVariables.size(), feblock, cell_centered);

        std::size_t li = z.mDataStart, got = 0;
        while (got < want && li < rLines.size()) {
            got += tecplot_tokens(rLines[li]).size();
            ++li;
        }
        li += z.mNumCells;  // one connectivity line per cell
        zones.push_back(z);
        i = li - 1;  // the for-loop's ++i resumes scanning right after
    }
    if (rVariables.empty())
        throw ReadError("Tecplot: no VARIABLES");
    if (zones.empty())
        throw ReadError("Tecplot: no ZONE");
    return zones;
}

// The zones read_tecplot/read_tecplot_metadata treat as one timeline: those
// sharing rZones[0]'s STRANDID when any zone carries one (SOLUTIONTIME with
// no STRANDID at all groups every zone together), sorted by SOLUTIONTIME.
// Zones with no SOLUTIONTIME at all are not a timeline -- multiple such
// zones is the "several static zones" case roadmap §7 owns, not this one;
// only the first is read here, with a warning if there is more than one.
std::vector<std::size_t> tecplot_timeline(const std::vector<TecplotZoneHeader>& rZones) {
    if (!rZones[0].mHasSolutionTime) {
        if (rZones.size() > 1)
            log::warn(
                "Tecplot: {} zones with no SOLUTIONTIME; reading the first only (multiple "
                "non-transient zones are not yet concatenated)",
                rZones.size());
        return {0};
    }
    std::vector<std::size_t> idx;
    for (std::size_t k = 0; k < rZones.size(); ++k) {
        if (!rZones[k].mHasSolutionTime)
            continue;
        if (rZones[0].mHasStrandId && rZones[k].mHasStrandId &&
            rZones[k].mStrandId != rZones[0].mStrandId)
            continue;
        idx.push_back(k);
    }
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        return rZones[a].mSolutionTime < rZones[b].mSolutionTime;
    });
    return idx;
}

}  // namespace

MeshMetadata read_tecplot_metadata(const std::string& rPath, const ReadOptions& /*rOptions*/) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        std::string s = tecplot_strip(l);
        if (s.empty() || s[0] == '#')
            continue;
        lines.push_back(s);
    }

    std::vector<std::string> variables;
    const std::vector<TecplotZoneHeader> zones = tecplot_scan_zones(lines, variables);
    const std::vector<std::size_t> timeline = tecplot_timeline(zones);

    MeshMetadata meta;
    meta.mFormat = "tecplot";
    const TecplotZoneHeader& first = zones[timeline[0]];
    meta.mNumPoints = first.mNumNodes;
    meta.mPointDim = 0;  // not knowable without decoding X/Y/Z columns
    CellBlockInfo block;
    bool feblock = false;
    std::string ztype;
    std::vector<int> cell_centered;
    tecplot_zone_format(first, variables.size(), feblock, ztype, cell_centered);
    block.mType = tecplot_to_meshio(ztype);
    block.mNumCells = first.mNumCells;
    meta.mCellBlocks.push_back(std::move(block));
    if (first.mHasSolutionTime)
        for (std::size_t idx : timeline)
            meta.mTimeValues.push_back(zones[idx].mSolutionTime);
    return meta;
}

Mesh read_tecplot(const std::string& rPath, const ReadOptions& rOptions) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l)) {
        std::string s = tecplot_strip(l);
        if (s.empty() || s[0] == '#')
            continue;
        lines.push_back(s);
    }

    std::vector<std::string> variables;
    const std::vector<TecplotZoneHeader> zones = tecplot_scan_zones(lines, variables);
    const std::vector<std::size_t> timeline = tecplot_timeline(zones);
    const std::size_t step = rOptions.ResolveTimeStep(timeline.size());
    const TecplotZoneHeader& zone_hdr = zones[timeline[step]];
    const std::size_t data_start = zone_hdr.mDataStart;
    const std::size_t num_nodes = zone_hdr.mNumNodes;
    const std::size_t num_cells = zone_hdr.mNumCells;

    bool feblock = false;
    std::string ztype;
    std::vector<int> cell_centered;
    tecplot_zone_format(zone_hdr, variables.size(), feblock, ztype, cell_centered);

    // Read data values.
    std::vector<std::size_t> ndata(variables.size());
    std::size_t total = 0;
    for (std::size_t k = 0; k < variables.size(); ++k) {
        ndata[k] = cell_centered[k] ? num_cells : num_nodes;
        total += ndata[k];
    }
    std::size_t want = feblock ? total : num_nodes * variables.size();

    std::vector<double> flat;
    flat.reserve(want);
    std::size_t li = data_start;
    while (flat.size() < want && li < lines.size()) {
        for (const auto& t : tecplot_tokens(lines[li]))
            flat.push_back(std::strtod(t.c_str(), nullptr));
        ++li;
    }

    // Per-variable columns.
    std::vector<std::vector<double>> cols(variables.size());
    if (feblock) {
        std::size_t off = 0;
        for (std::size_t k = 0; k < variables.size(); ++k) {
            cols[k].assign(flat.begin() + off, flat.begin() + off + ndata[k]);
            off += ndata[k];
        }
    } else {
        std::size_t nv = variables.size();
        for (std::size_t k = 0; k < nv; ++k)
            cols[k].resize(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            for (std::size_t k = 0; k < nv; ++k)
                cols[k][r] = flat[r * nv + k];
    }

    // Cells.
    std::string mtype = tecplot_to_meshio(ztype);
    if (mtype.empty())
        throw ReadError("Tecplot: unsupported zone type " + ztype);
    std::size_t nn;
    if (mtype == "line")
        nn = 2;
    else if (mtype == "triangle")
        nn = 3;
    else if (mtype == "quad" || mtype == "tetra")
        nn = 4;
    else
        nn = 8;
    NDArray celldata(DType::Int64, {num_cells, nn});
    std::int64_t* cp = celldata.As<std::int64_t>();
    for (std::size_t c = 0; c < num_cells; ++c) {
        auto t = tecplot_tokens(lines.at(li++));
        for (std::size_t j = 0; j < nn; ++j)
            cp[c * nn + j] = std::strtoll(t[j].c_str(), nullptr, 10) - 1;
    }

    // Assemble.
    Mesh mesh;
    int xi = -1, yi = -1, zi = -1;
    for (std::size_t k = 0; k < variables.size(); ++k) {
        std::string v = tecplot_upper(variables[k]);
        if (v == "X")
            xi = (int)k;
        else if (v == "Y")
            yi = (int)k;
        else if (v == "Z")
            zi = (int)k;
    }
    std::size_t ndim = (zi >= 0) ? 3 : 2;
    NDArray pts(DType::Float64, {num_nodes, ndim});
    double* pp = pts.As<double>();
    for (std::size_t r = 0; r < num_nodes; ++r) {
        pp[r * ndim + 0] = cols[xi][r];
        pp[r * ndim + 1] = cols[yi][r];
        if (zi >= 0)
            pp[r * ndim + 2] = cols[zi][r];
    }
    mesh.AssignPoints(std::move(pts));
    for (std::size_t k = 0; k < variables.size(); ++k) {
        if ((int)k == xi || (int)k == yi || (int)k == zi)
            continue;
        NDArray arr(DType::Float64, {cols[k].size()});
        std::memcpy(arr.Data(), cols[k].data(), cols[k].size() * sizeof(double));
        if (cell_centered[k]) {
            std::vector<NDArray> blk;
            blk.push_back(std::move(arr));
            mesh.AddCellData(variables[k], std::move(blk));
        } else {
            mesh.AddPointData(variables[k], std::move(arr));
        }
    }
    mesh.AddCellBlock(mtype, std::move(celldata));
    return mesh;
}

Mesh read_tecplot(const std::string& rPath) {
    return read_tecplot(rPath, ReadOptions{});
}

void write_tecplot(const std::string& rPath, const Mesh& rMesh) {
    // Gather supported cell blocks; require a single unique type.
    std::vector<std::size_t> blocks;
    std::set<std::string> types;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const auto cb = rMesh.Cells(i);
        if (!meshio_to_tecplot(cb.Type()).empty()) {
            blocks.push_back(i);
            types.insert(cb.Type());
        }
    }
    if (types.size() != 1)
        throw WriteError("C++ Tecplot writer supports a single cell type");
    std::string mtype = *types.begin();
    std::string ztype = meshio_to_tecplot(mtype);
    const std::vector<int>& order = tecplot_order(mtype);

    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t dim = rMesh.PointDim();
    const std::size_t num_nodes = rMesh.NumPoints();
    std::size_t num_cells = 0;
    for (std::size_t b : blocks)
        num_cells += rMesh.Cells(b).NumCells();

    // Variables + data columns.
    const NDArray& points = rMesh.Points();
    std::vector<std::string> variables = {"X", "Y"};
    std::vector<std::vector<double>> data;
    auto push_point_col = [&](std::size_t comp) {
        std::vector<double> col(num_nodes);
        for (std::size_t r = 0; r < num_nodes; ++r)
            col[r] = detail::read_double(points, r * dim + comp);
        data.push_back(std::move(col));
    };
    push_point_col(0);
    push_point_col(1);
    int varrange0 = 3, varrange1 = 0;
    if (dim == 3) {
        variables.push_back("Z");
        push_point_col(2);
        varrange0 += 1;
    }

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
            data.push_back(std::move(col));
            varrange0 += 1;
        }
    }
    bool have_cell_data = false;
    varrange1 = varrange0 - 1;
    for (const auto& k : rMesh.CellDataNames()) {
        std::string ku = tecplot_upper(k);
        if (ku == "X" || ku == "Y" || ku == "Z")
            continue;
        if (rMesh.CellDataNumBlocks(k) == 0)
            continue;
        // concatenate the (single-type) blocks
        const NDArray& first = rMesh.CellData(k, 0);
        std::size_t ncomp = first.Shape().size() >= 2 ? first.Shape()[1] : 1;
        for (std::size_t c = 0; c < ncomp; ++c) {
            variables.push_back(ncomp == 1 ? k : k + "_" + std::to_string(c));
            std::vector<double> col;
            for (std::size_t b : blocks) {
                const NDArray& vv = rMesh.CellData(k, b);
                for (std::size_t r = 0; r < vv.Shape()[0]; ++r)
                    col.push_back(detail::read_double(vv, r * ncomp + c));
            }
            data.push_back(std::move(col));
            varrange1 += 1;
            have_cell_data = true;
        }
    }

    os << "TITLE = \"" << detail::provenance_lines(detail::SlotTier::SingleLine)[0] << "\"\n";
    os << "VARIABLES = ";
    for (std::size_t k = 0; k < variables.size(); ++k)
        os << (k ? ", " : "") << "\"" << variables[k] << "\"";
    os << "\n";
    os << "ZONE NODES = " << num_nodes << ", ELEMENTS = " << num_cells << ",\n";
    os << "DATAPACKING = BLOCK, ZONETYPE = " << ztype;
    if (have_cell_data && varrange0 <= varrange1) {
        os << ",\n";
        std::string r = (varrange0 == varrange1)
                            ? std::to_string(varrange0)
                            : std::to_string(varrange0) + "-" + std::to_string(varrange1);
        os << "VARLOCATION = ([" << r << "] = CELLCENTERED)\n";
    } else {
        os << "\n";
    }

    char buf[40];
    for (const auto& col : data) {
        for (std::size_t i = 0; i < col.size(); ++i) {
            detail::snprintf_c(buf, sizeof(buf), "%.17g", col[i]);
            os << buf << ((i + 1) % 20 == 0 || i + 1 == col.size() ? '\n' : ' ');
        }
        if (col.empty())
            os << "\n";
    }

    for (std::size_t b : blocks) {
        const auto cb = rMesh.Cells(b);
        const NDArray& conn = cb.Conn();
        std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
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
