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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/log.hpp"
#include "../detail/text_cursor.hpp"

namespace meshioplusplus {

namespace {

std::string desc3_to_meshio(int d) {
    static const std::unordered_map<int, std::string> m = {
        {2, "vertex"},        {3, "line"},    {4, "line3"},    {5, "triangle"},
        {6, "triangle6"},     {7, "quad"},    {8, "quad8"},    {10, "tetra"},
        {11, "tetra10"},      {12, "wedge"},  {13, "wedge15"}, {15, "hexahedron"},
        {16, "hexahedron20"}, {17, "pyramid"}};
    auto it = m.find(d);
    return it == m.end() ? std::string() : it->second;
}

// meshio type -> (desc1, desc2, desc3)
bool meshio_to_desc(const std::string& rT, std::array<int, 3>& rOut) {
    static const std::unordered_map<std::string, std::array<int, 3>> m = {
        {"vertex", {1, 1, 2}},           {"line", {2, 2, 3}},       {"line3", {2, 3, 4}},
        {"triangle", {3, 7, 5}},         {"triangle6", {3, 7, 6}},  {"quad", {4, 202, 7}},
        {"quad8", {4, 303, 8}},          {"tetra", {5, 4, 10}},     {"tetra10", {5, 15, 11}},
        {"wedge", {6, 207, 12}},         {"wedge15", {6, 307, 13}}, {"hexahedron", {7, 2202, 15}},
        {"hexahedron20", {7, 3303, 16}}, {"pyramid", {8, 4202, 17}}};
    auto it = m.find(rT);
    if (it == m.end())
        return false;
    rOut = it->second;
    return true;
}

bool contains(const std::string& rHay, const char* pNeedle) {
    return rHay.find(pNeedle) != std::string::npos;
}

long long leading_int(const std::string& rLine) {
    detail::TextStream iss(rLine);
    long long v = 0;
    iss >> v;
    return v;
}

}  // namespace

Mesh read_flux(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);

    long long dim = 0, nel = 0, nnod = 0;
    std::size_t di = lines.size(), ci = lines.size();
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string& L = lines[i];
        if (contains(L, "NOMBRE DE DIMENSIONS"))
            dim = leading_int(L);
        else if (contains(L, "D'ELEMENTS") && !contains(L, "VOLUMIQUES") &&
                 !contains(L, "SURFACIQUES") && !contains(L, "LINEIQUES") &&
                 !contains(L, "PONCTUELS") && !contains(L, "MACRO"))
            nel = leading_int(L);
        else if (contains(L, "NOMBRE DE POINTS") && !contains(L, "INTEGRATION"))
            nnod = leading_int(L);
        else if (contains(L, "DESCRIPTEUR DE TOPOLOGIE"))
            di = i;
        else if (contains(L, "COORDONNEES DES NOEUDS"))
            ci = i;
    }
    if (di >= lines.size() || ci >= lines.size())
        throw ReadError("pf3: missing element/coordinate section");

    // element tokens
    std::vector<std::string> etok;
    for (std::size_t i = di + 1; i < ci; ++i) {
        detail::TextStream iss(lines[i]);
        std::string w;
        while (iss >> w)
            etok.push_back(w);
    }

    struct Group {
        std::string mType;
        std::vector<std::vector<std::int64_t>> mRows;
        std::vector<std::int64_t> mRef;
    };
    std::vector<Group> groups;
    std::unordered_map<std::string, std::size_t> gindex;
    std::size_t pos = 0;
    for (long long e = 0; e < nel; ++e) {
        // A file cut short (FLUX excerpts are) keeps the elements it has.
        if (pos + 12 > etok.size()) {
            log::warn("pf3: {} elements declared, {} present; reading those", nel, e);
            break;
        }
        long long ref = std::strtoll(etok[pos + 3].c_str(), nullptr, 10);
        int desc3 = std::atoi(etok[pos + 6].c_str());
        int lnn = std::atoi(etok[pos + 7].c_str());
        pos += 12;
        if (lnn < 0 || pos + static_cast<std::size_t>(lnn) > etok.size())
            throw ReadError("pf3: truncated element connectivity");
        std::string mtype = desc3_to_meshio(desc3);
        if (mtype.empty())
            throw ReadError("pf3: unknown element descriptor " + std::to_string(desc3));
        std::vector<std::int64_t> nodes(lnn);
        for (int j = 0; j < lnn; ++j)
            nodes[j] = std::strtoll(etok[pos + j].c_str(), nullptr, 10);
        pos += lnn;
        auto it = gindex.find(mtype);
        if (it == gindex.end()) {
            gindex[mtype] = groups.size();
            groups.push_back({mtype, {}, {}});
            it = gindex.find(mtype);
        }
        if (!groups[it->second].mRows.empty() &&
            groups[it->second].mRows.front().size() != nodes.size())
            throw ReadError("pf3: '" + mtype + "' elements with different node counts");
        groups[it->second].mRows.push_back(std::move(nodes));
        groups[it->second].mRef.push_back(ref);
    }

    // `id x1 .. x_dim` rows up to the `==== DECOUPAGE TERMINE` trailer.
    std::vector<std::string> ctok;
    for (std::size_t i = ci + 1; i < lines.size(); ++i) {
        detail::TextStream iss(lines[i]);
        std::string w;
        while (iss >> w)
            ctok.push_back(w);
    }
    std::vector<long long> ids;
    std::vector<double> coords;
    const std::size_t udim = static_cast<std::size_t>(dim < 0 ? 0 : dim);
    for (std::size_t cp = 0; cp + udim < ctok.size() && ctok[cp][0] != '=';) {
        ids.push_back(std::strtoll(ctok[cp].c_str(), nullptr, 10));
        for (std::size_t j = 0; j < udim; ++j)
            coords.push_back(detail::parse_double(ctok[cp + 1 + j]));
        cp += 1 + udim;
    }
    if (static_cast<long long>(ids.size()) != nnod)
        log::warn("pf3: {} points declared, {} present", nnod, ids.size());

    Mesh mesh;
    NDArray pts(DType::Float64, {ids.size(), udim});
    std::copy(coords.begin(), coords.end(), pts.As<double>());
    mesh.AssignPoints(std::move(pts));

    // Node ids are normally 1..n in row order; otherwise number the rows.
    bool sequential = true;
    for (std::size_t r = 0; r < ids.size() && sequential; ++r)
        sequential = ids[r] == static_cast<long long>(r) + 1;
    std::unordered_map<long long, std::int64_t> row_of;
    if (!sequential)
        for (std::size_t r = 0; r < ids.size(); ++r)
            row_of.emplace(ids[r], static_cast<std::int64_t>(r));
    auto to_row = [&](std::int64_t Id) -> std::int64_t {
        if (sequential) {
            if (Id < 1 || Id > static_cast<std::int64_t>(ids.size()))
                throw ReadError("pf3: element references an undefined node");
            return Id - 1;
        }
        auto it = row_of.find(Id);
        if (it == row_of.end())
            throw ReadError("pf3: element references undefined node " + std::to_string(Id));
        return it->second;
    };

    std::vector<NDArray> refs;
    for (auto& g : groups) {
        std::size_t ne = g.mRows.size();
        std::size_t k = ne ? g.mRows[0].size() : 0;
        const detail::NodeOrder* order = detail::node_order("flux", g.mType);
        if (order && order->mToMeshio.size() != k)
            order = nullptr;
        NDArray data(DType::Int64, {ne, k});
        std::int64_t* pData = data.As<std::int64_t>();
        for (std::size_t r = 0; r < ne; ++r)
            for (std::size_t j = 0; j < k; ++j)
                pData[r * k + j] =
                    to_row(g.mRows[r][order ? static_cast<std::size_t>(order->mToMeshio[j]) : j]);
        mesh.AddCellBlock(g.mType, std::move(data));
        NDArray rf(DType::Int64, {ne});
        for (std::size_t r = 0; r < ne; ++r)
            rf.As<std::int64_t>()[r] = g.mRef[r];
        refs.push_back(std::move(rf));
    }
    if (!refs.empty())
        mesh.AddCellData("pf3:ref", std::move(refs));
    return mesh;
}

void write_flux(const std::string& rPath, const Mesh& rMesh) {
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);

    const int dim = static_cast<int>(rMesh.PointDim());
    long long counts[4] = {0, 0, 0, 0};  // by topological dim
    std::vector<std::size_t> blocks;
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        std::array<int, 3> d;
        if (!meshio_to_desc(cb.Type(), d))
            throw WriteError("pf3: unsupported cell type " + cb.Type());
        auto it = topological_dimension().find(cb.Type());
        int td = it == topological_dimension().end() ? 3 : it->second;
        counts[td] += static_cast<long long>(cb.NumCells());
        blocks.push_back(k);
    }
    long long nel = 0;
    for (auto k : blocks)
        nel += static_cast<long long>(rMesh.Cells(k).NumCells());

    const bool has_ref = rMesh.HasCellData("pf3:ref");

    char buf[128];
    f << detail::provenance_render_lines(detail::SlotTier::Block, " ");
    auto hdr = [&](long long v, const char* label) {
        std::snprintf(buf, sizeof(buf), "%8lld           %s\n", v, label);
        f << buf;
    };
    hdr(dim, "NOMBRE DE DIMENSIONS DU DECOUPAGE");
    hdr(nel, "NOMBRE  D'ELEMENTS");
    hdr(counts[3], "NOMBRE  D'ELEMENTS VOLUMIQUES");
    hdr(counts[2], "NOMBRE  D'ELEMENTS SURFACIQUES");
    hdr(counts[1], "NOMBRE  D'ELEMENTS LINEIQUES");
    hdr(counts[0], "NOMBRE  D'ELEMENTS PONCTUELS");
    hdr(0, "NOMBRE DE MACRO-ELEMENTS");
    hdr(static_cast<long long>(rMesh.NumPoints()), "NOMBRE DE POINTS");
    hdr(1, "NOMBRE DE REGIONS");
    hdr(0, "NOMBRE DE REGIONS VOLUMIQUES");
    hdr(0, "NOMBRE DE REGIONS SURFACIQUES");
    hdr(0, "NOMBRE DE REGIONS LINEIQUES");
    hdr(0, "NOMBRE DE REGIONS PONCTUELLES");
    hdr(0, "NOMBRE DE REGIONS MACRO-ELEMENTAIRES");
    hdr(20, "NOMBRE DE NOEUDS DANS 1 ELEMENT (MAX)");
    hdr(20, "NOMBRE DE POINTS D'INTEGRATION / ELEMENT (MAX)");
    f << " NOMS DES REGIONS\n";
    f << " DESCRIPTEUR DE TOPOLOGIE DES ELEMENTS\n";

    long long eid = 0;
    for (auto k : blocks) {
        const auto cb = rMesh.Cells(k);
        std::array<int, 3> d;
        meshio_to_desc(cb.Type(), d);
        const NDArray& conn = cb.Conn();
        int lnn = static_cast<int>(detail::cols(conn));
        const detail::NodeOrder* order = detail::node_order("flux", cb.Type());
        if (order && order->mFromMeshio.size() != static_cast<std::size_t>(lnn))
            order = nullptr;
        const NDArray* ref = (has_ref && k < rMesh.CellDataNumBlocks("pf3:ref"))
                                 ? &rMesh.CellData("pf3:ref", k)
                                 : nullptr;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            ++eid;
            long long rv = ref ? detail::read_int(*ref, r) : 0;
            std::snprintf(buf, sizeof(buf), "%8lld%8d%8d%8lld%8d%8d%8d%8d%8d%8d%8d%8d\n", eid, d[0],
                          d[1], rv, lnn, 0, d[2], lnn, 0, 0, 0, 0);
            f << buf;
            for (int j = 0; j < lnn; ++j) {
                const std::size_t src =
                    order ? static_cast<std::size_t>(order->mFromMeshio[j]) : std::size_t(j);
                std::snprintf(buf, sizeof(buf), "%8lld",
                              static_cast<long long>(detail::read_int(conn, r * lnn + src) + 1));
                f << buf;
            }
            f << "\n";
        }
    }

    f << " COORDONNEES DES NOEUDS\n";
    const NDArray& points = rMesh.Points();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        std::snprintf(buf, sizeof(buf), "%8zu", i + 1);
        f << buf;
        for (int j = 0; j < dim; ++j) {
            detail::snprintf_c(buf, sizeof(buf), " %.16g",
                               detail::read_double(points, i * dim + j));
            f << buf;
        }
        f << "\n";
    }
}

}  // namespace meshioplusplus
