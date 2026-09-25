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
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

namespace meshioplusplus {

namespace {

int su2_numnodes(int t) {
    switch (t) {
        case 3:
            return 2;  // line
        case 5:
            return 3;  // triangle
        case 9:
            return 4;  // quad
        case 10:
            return 4;  // tetra
        case 12:
            return 8;  // hexahedron
        case 13:
            return 6;  // wedge
        case 14:
            return 5;  // pyramid
        default:
            return 0;
    }
}
std::string su2_to_meshio(int t) {
    switch (t) {
        case 3:
            return "line";
        case 5:
            return "triangle";
        case 9:
            return "quad";
        case 10:
            return "tetra";
        case 12:
            return "hexahedron";
        case 13:
            return "wedge";
        case 14:
            return "pyramid";
        default:
            return "";
    }
}
int meshio_to_su2(const std::string& rT) {
    if (rT == "line")
        return 3;
    if (rT == "triangle")
        return 5;
    if (rT == "quad")
        return 9;
    if (rT == "tetra")
        return 10;
    if (rT == "hexahedron")
        return 12;
    if (rT == "wedge")
        return 13;
    if (rT == "pyramid")
        return 14;
    return -1;
}

std::string su2_strip(const std::string& rS) {
    std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}
std::vector<std::string> su2_tokens(const std::string& rS) {
    std::vector<std::string> out;
    auto iss = detail::make_classic_istringstream(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

struct Blk {
    std::string mType;
    int mN = 0;
    std::vector<std::int64_t> mConn;
    std::vector<std::int32_t> mTag;
    std::vector<std::int32_t> mZone;  // which IZONE this row came from (0 for a single-zone file)
    std::size_t mCount = 0;
};

// Parse `count` element lines (each "vtk_type n0 n1 ... [extra]") into type-
// grouped blocks (sorted by vtk type code, matching numpy.unique), all with
// the given tag and zone.
void read_elem_block(const std::vector<std::string>& rLines, std::size_t& rLi, std::size_t count,
                     std::int32_t tag, std::int32_t zone, std::vector<Blk>& rOut) {
    std::vector<std::pair<int, std::vector<std::int64_t>>> elems;
    std::set<int> types;
    for (std::size_t e = 0; e < count; ++e) {
        auto t = su2_tokens(rLines.at(rLi++));
        detail::need_tokens(t, 1, "SU2");
        int vt = std::stoi(t[0]);
        int nn = su2_numnodes(vt);
        if (nn == 0)
            throw ReadError("SU2: unsupported element type " + t[0]);
        detail::need_tokens(t, 1 + static_cast<std::size_t>(nn), "SU2");
        std::vector<std::int64_t> nodes(nn);
        for (int j = 0; j < nn; ++j)
            nodes[j] = std::strtoll(t[1 + j].c_str(), nullptr, 10);
        elems.emplace_back(vt, std::move(nodes));
        types.insert(vt);
    }
    for (int vt : types) {  // std::set is sorted
        Blk b;
        b.mType = su2_to_meshio(vt);
        b.mN = su2_numnodes(vt);
        for (auto& e : elems) {
            if (e.first != vt)
                continue;
            b.mConn.insert(b.mConn.end(), e.second.begin(), e.second.end());
            b.mTag.push_back(tag);
            b.mZone.push_back(zone);
            ++b.mCount;
        }
        rOut.push_back(std::move(b));
    }
}

/// One zone's own NDIME/NPOIN/NELEM/NMARK/MARKER_* body, everything a
/// standalone (single-zone) SU2 file also has.
struct Su2ZoneBody {
    int mDim = 0;
    NDArray mPoints;
    std::vector<Blk> mBlocks;
    // marker tag -> its string name, only for markers whose MARKER_TAG was
    // non-numeric (a plain integer tag needs no name to preserve).
    std::map<std::int32_t, std::string> mMarkerNames;
};

/// Parses one zone body starting at `*pLi`, advancing it past everything
/// consumed. Stops at the next "IZONE=" (the next zone) or end of file --
/// shared by the single-zone and multizone read paths, so both go through
/// exactly the same per-record logic.
Su2ZoneBody read_su2_zone_body(const std::vector<std::string>& rLines, std::size_t& rLi,
                               std::int32_t zone) {
    Su2ZoneBody zoneBody;
    std::int32_t next_tag_id = 0;
    std::int32_t current_tag = 0;
    std::string current_tag_name;  // empty when the marker's own tag is numeric

    while (rLi < rLines.size()) {
        std::string line = su2_strip(rLines[rLi]);
        if (line.empty() || line[0] == '%') {
            ++rLi;
            continue;
        }
        std::size_t eq = line.find('=');
        if (eq == std::string::npos) {
            ++rLi;
            continue;
        }
        std::string name = su2_strip(line.substr(0, eq));
        if (name == "IZONE")
            break;  // the next zone: let the caller consume it
        std::string rest = su2_strip(line.substr(eq + 1));
        ++rLi;

        if (name == "NDIME") {
            zoneBody.mDim = std::stoi(rest);
            if (zoneBody.mDim != 2 && zoneBody.mDim != 3)
                throw ReadError("SU2: invalid NDIME");
        } else if (name == "NPOIN") {
            const auto npoin_tok = su2_tokens(rest);
            detail::need_tokens(npoin_tok, 1, "SU2");
            if (zoneBody.mDim == 0)
                throw ReadError("SU2: NPOIN before NDIME");
            // One point per line: the count cannot exceed the lines left.
            const std::size_t npoin = detail::checked_count(std::stoll(npoin_tok[0]),
                                                            rLines.size() - rLi, "SU2", "point");
            NDArray pts(DType::Float64, {npoin, static_cast<std::size_t>(zoneBody.mDim)});
            double* pp = pts.As<double>();
            for (std::size_t i = 0; i < npoin; ++i) {
                auto t = su2_tokens(rLines.at(rLi++));
                detail::need_tokens(t, static_cast<std::size_t>(zoneBody.mDim), "SU2");
                for (int c = 0; c < zoneBody.mDim; ++c)
                    pp[i * static_cast<std::size_t>(zoneBody.mDim) + static_cast<std::size_t>(c)] =
                        detail::parse_double(t[static_cast<std::size_t>(c)]);
            }
            zoneBody.mPoints = std::move(pts);
        } else if (name == "NELEM") {
            const std::size_t ne =
                detail::checked_count(std::stoll(rest), rLines.size() - rLi, "SU2", "element");
            read_elem_block(rLines, rLi, ne, 0, zone, zoneBody.mBlocks);
        } else if (name == "NMARK") {
            // handled implicitly via MARKER_TAG/MARKER_ELEMS
        } else if (name == "MARKER_TAG") {
            try {
                std::size_t pos;
                int v = std::stoi(rest, &pos);
                if (pos == rest.size()) {
                    current_tag = v;
                    current_tag_name.clear();
                } else {
                    current_tag = ++next_tag_id;
                    current_tag_name = rest;
                }
            } catch (...) {
                current_tag = ++next_tag_id;
                current_tag_name = rest;
            }
            if (!current_tag_name.empty())
                zoneBody.mMarkerNames[current_tag] = current_tag_name;
        } else if (name == "MARKER_ELEMS") {
            const std::size_t ne =
                detail::checked_count(std::stoll(rest), rLines.size() - rLi, "SU2", "element");
            read_elem_block(rLines, rLi, ne, current_tag, zone, zoneBody.mBlocks);
        }
    }
    return zoneBody;
}

/// Merges same-type blocks (across every zone, volume or boundary alike)
/// into one block per type -- a straight generalization of the single-zone
/// reader's own "merge same-type boundary blocks" pass, which becomes a
/// no-op there (a single NELEM/MARKER_ELEMS call already groups by type).
void su2_merge_by_type(std::vector<Blk>& rBlocks) {
    std::vector<std::string> types;
    for (const Blk& b : rBlocks)
        if (std::find(types.begin(), types.end(), b.mType) == types.end())
            types.push_back(b.mType);
    for (const std::string& ty : types) {
        int first = -1;
        for (std::size_t i = 0; i < rBlocks.size(); ++i) {
            if (rBlocks[i].mType != ty)
                continue;
            if (first < 0) {
                first = static_cast<int>(i);
                continue;
            }
            Blk& dst = rBlocks[static_cast<std::size_t>(first)];
            Blk& src = rBlocks[i];
            dst.mConn.insert(dst.mConn.end(), src.mConn.begin(), src.mConn.end());
            dst.mTag.insert(dst.mTag.end(), src.mTag.begin(), src.mTag.end());
            dst.mZone.insert(dst.mZone.end(), src.mZone.begin(), src.mZone.end());
            dst.mCount += src.mCount;
            src.mCount = 0;  // mark for removal
            src.mConn.clear();
        }
    }
}

}  // namespace

Mesh read_su2(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l))
        lines.push_back(l);

    // NZONE= (if present) is always the first key: a single-file multizone
    // mesh's own header, before the first (implicit) IZONE= 1.
    std::size_t nzone = 1;
    bool multizone = false;
    {
        std::size_t li = 0;
        while (li < lines.size()) {
            std::string line = su2_strip(lines[li]);
            if (line.empty() || line[0] == '%') {
                ++li;
                continue;
            }
            std::size_t eq = line.find('=');
            if (eq == std::string::npos)
                break;
            if (su2_strip(line.substr(0, eq)) == "NZONE") {
                const auto value = su2_tokens(su2_strip(line.substr(eq + 1)));
                detail::need_tokens(value, 1, "SU2");
                // Every zone takes lines of its own: bounded by the file.
                nzone = detail::checked_count(std::stoll(value[0]), lines.size(), "SU2", "zone");
                multizone = true;
            }
            break;
        }
    }

    std::vector<Su2ZoneBody> zones;
    std::size_t li = 0;
    for (std::size_t z = 0; z < nzone; ++z) {
        // Skip blank/comment lines and (for a multizone file) the NZONE=/IZONE=
        // markers themselves; the shared body parser starts at NDIME.
        while (li < lines.size()) {
            std::string line = su2_strip(lines[li]);
            if (line.empty() || line[0] == '%') {
                ++li;
                continue;
            }
            std::size_t eq = line.find('=');
            if (eq == std::string::npos)
                break;
            std::string name = su2_strip(line.substr(0, eq));
            if (name == "NZONE" || name == "IZONE") {
                const auto value = su2_tokens(su2_strip(line.substr(eq + 1)));
                if (name == "IZONE" &&
                    (value.empty() || std::stoll(value[0]) != static_cast<long long>(z + 1)))
                    throw ReadError("SU2: IZONE out of order (expected " + std::to_string(z + 1) +
                                    ")");
                ++li;
                continue;
            }
            break;
        }
        zones.push_back(read_su2_zone_body(lines, li, static_cast<std::int32_t>(z)));
    }
    if (multizone && zones.size() != nzone)
        throw ReadError("SU2: NZONE=" + std::to_string(nzone) + " but found " +
                        std::to_string(zones.size()) + " IZONE section(s)");

    Mesh mesh;
    const int dim = zones.empty() ? 0 : zones.front().mDim;
    {
        // Concatenate zone points, offsetting each zone's connectivity by the
        // running point count -- zones are independent meshes, never welded.
        std::size_t total_points = 0;
        for (const Su2ZoneBody& z : zones) {
            // Also catches an NDIME after NPOIN, which leaves the points sized
            // for the earlier dimension.
            if (z.mDim != dim || (z.mPoints.Shape().size() == 2 &&
                                  z.mPoints.Shape()[1] != static_cast<std::size_t>(dim)))
                throw ReadError("SU2: zones of different dimensions (NDIME)");
            total_points += z.mPoints.Shape().empty() ? 0 : z.mPoints.Shape()[0];
        }
        NDArray pts(DType::Float64, {total_points, static_cast<std::size_t>(dim)});
        double* pp = pts.As<double>();
        std::size_t offset = 0;
        for (Su2ZoneBody& z : zones) {
            const std::size_t n = z.mPoints.Shape().empty() ? 0 : z.mPoints.Shape()[0];
            if (n)
                std::memcpy(pp + offset * static_cast<std::size_t>(dim), z.mPoints.Data(),
                           z.mPoints.Nbytes());
            for (Blk& b : z.mBlocks)
                for (std::int64_t& id : b.mConn)
                    id += static_cast<std::int64_t>(offset);
            offset += n;
        }
        mesh.AssignPoints(std::move(pts));
    }

    std::vector<Blk> blocks;
    for (Su2ZoneBody& z : zones)
        for (Blk& b : z.mBlocks)
            blocks.push_back(std::move(b));
    su2_merge_by_type(blocks);

    std::vector<NDArray> tags;
    std::vector<NDArray> zone_ids;
    std::vector<Blk*> kept;  // in AddCellBlock order, for the region pass below
    for (auto& b : blocks) {
        if (b.mCount == 0)
            continue;  // merged-away or empty
        NDArray data(DType::Int64, {b.mCount, static_cast<std::size_t>(b.mN)});
        std::memcpy(data.Data(), b.mConn.data(), b.mConn.size() * sizeof(std::int64_t));
        mesh.AddCellBlock(b.mType, std::move(data));
        NDArray tg(DType::Int32, {b.mCount});
        std::memcpy(tg.Data(), b.mTag.data(), b.mTag.size() * sizeof(std::int32_t));
        tags.push_back(std::move(tg));
        if (multizone) {
            NDArray zn(DType::Int32, {b.mCount});
            std::memcpy(zn.Data(), b.mZone.data(), b.mZone.size() * sizeof(std::int32_t));
            zone_ids.push_back(std::move(zn));
        }
        kept.push_back(&b);
    }
    mesh.AddCellData("su2:tag", std::move(tags));
    if (multizone)
        mesh.AddCellData("su2:zone", std::move(zone_ids));

    // Regions: one "zone_<i>" per zone (multizone only) and one per named
    // marker -- "zone_<i>/<name>" when multizone, "<name>" for a single-zone
    // file, so a marker's own string name survives instead of collapsing to
    // an anonymous auto-incremented su2:tag. A purely numeric marker gets no
    // region: su2:tag already carries it losslessly.
    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    std::map<std::int32_t, std::vector<std::int64_t>> zone_cells;
    std::map<std::pair<std::int32_t, std::int32_t>, std::vector<std::int64_t>> marker_cells;
    for (std::size_t bi = 0; bi < kept.size(); ++bi) {
        const Blk& b = *kept[bi];
        for (std::size_t r = 0; r < b.mCount; ++r) {
            const std::int64_t global = detail::block_row_to_global(bases, bi, static_cast<std::int64_t>(r));
            zone_cells[b.mZone[r]].push_back(global);
            marker_cells[{b.mZone[r], b.mTag[r]}].push_back(global);
        }
    }
    if (multizone)
        for (const auto& [zone, cells] : zone_cells) {
            NDArray entries(DType::Int64, {cells.size()});
            std::memcpy(entries.Data(), cells.data(), cells.size() * sizeof(std::int64_t));
            // -1: a zone region spans both volume and boundary cells, so it has
            // no single topological dimension to declare.
            mesh.AddRegion(Region("zone_" + std::to_string(zone), RegionKind::Cell, -1,
                                  static_cast<std::int64_t>(zone), std::move(entries)));
        }
    for (const auto& [key, cells] : marker_cells) {
        const auto& [zone, tag] = key;
        const auto it = zones[static_cast<std::size_t>(zone)].mMarkerNames.find(tag);
        if (it == zones[static_cast<std::size_t>(zone)].mMarkerNames.end())
            continue;  // a numeric marker: su2:tag already carries it
        const std::string name =
            multizone ? "zone_" + std::to_string(zone) + "/" + it->second : it->second;
        NDArray entries(DType::Int64, {cells.size()});
        std::memcpy(entries.Data(), cells.data(), cells.size() * sizeof(std::int64_t));
        mesh.AddRegion(Region(name, RegionKind::Cell, dim - 1, tag, std::move(entries)));
    }
    return mesh;
}

namespace {

std::vector<std::string> su2_vtypes(std::size_t dim) {
    return dim == 2 ? std::vector<std::string>{"triangle", "quad"}
                    : std::vector<std::string>{"tetra", "hexahedron", "wedge", "pyramid"};
}
std::vector<std::string> su2_btypes(std::size_t dim) {
    return dim == 2 ? std::vector<std::string>{"line"} : std::vector<std::string>{"triangle", "quad"};
}
bool su2_in(const std::vector<std::string>& rV, const std::string& rT) {
    return std::find(rV.begin(), rV.end(), rT) != rV.end();
}

/// Which int cell_data array drives markers: the first one found, matching
/// the read side and the "first-int-array" convention several other formats
/// share.
std::string su2_tag_key(const Mesh& rMesh) {
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        DType t = rMesh.CellData(name, 0).Dtype();
        if (t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
            t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64)
            return name;
    }
    return "";
}

/// Writes one zone's own NDIME/NPOIN/NELEM/NMARK/MARKER_* body -- everything a
/// standalone (single-zone) file also has -- restricted to the cells named by
/// their (block, row) in `rCells`. `Subset`: when false (a whole-mesh,
/// single-zone write), every point is written in its original order and
/// connectivity is untouched, exactly like a pre-multizone write. When true
/// (one zone of a multizone write), points are remapped to a dense, zone-local
/// 0-based numbering in first-appearance order, so a multizone write never
/// shares a point index between zones (SU2 zones are independent meshes).
/// `rMarkerNames` resolves a marker's SU2 text: `zone_<i>/<name>` (or, for a
/// single-zone write, plain `<name>`) when the region carrying that name
/// exists, else the bare integer tag.
void write_su2_zone_body(std::ostream& rOs, const Mesh& rMesh, std::size_t Dim,
                         const std::vector<std::pair<std::size_t, std::int64_t>>& rCells,
                         const std::string& rTagKey,
                         const std::map<std::int64_t, std::string>& rMarkerNames, bool Subset) {
    const std::vector<std::string> vtypes = su2_vtypes(Dim);
    const std::vector<std::string> btypes = su2_btypes(Dim);
    const NDArray& points = rMesh.Points();

    // Point remap: identity (whole mesh, original order) unless Subset, in
    // which case only the cells in rCells' own points are written, renumbered
    // by first appearance.
    std::unordered_map<std::int64_t, std::int64_t> remap;
    std::vector<std::int64_t> old_ids;
    auto remap_of = [&](std::int64_t old_id) {
        if (!Subset)
            return old_id;
        const auto [it, inserted] = remap.try_emplace(old_id, static_cast<std::int64_t>(old_ids.size()));
        if (inserted)
            old_ids.push_back(old_id);
        return it->second;
    };
    if (Subset) {
        for (const auto& [block, row] : rCells) {
            const auto cb = rMesh.Cells(block);
            const NDArray& conn = cb.Conn();
            const std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
            for (std::size_t j = 0; j < k; ++j)
                remap_of(detail::read_int(conn, static_cast<std::size_t>(row) * k + j));
        }
    } else {
        old_ids.resize(rMesh.NumPoints());
        for (std::size_t i = 0; i < old_ids.size(); ++i)
            old_ids[i] = static_cast<std::int64_t>(i);
    }

    rOs << "NDIME= " << Dim << "\n";
    rOs << "NPOIN= " << old_ids.size() << "\n";
    for (std::int64_t old_id : old_ids) {
        for (std::size_t c = 0; c < Dim; ++c) {
            char buf[64];
            detail::snprintf_c(buf, sizeof(buf), "%.16e",
                               detail::read_double(points, static_cast<std::size_t>(old_id) * Dim + c));
            rOs << buf << (c + 1 == Dim ? '\n' : ' ');
        }
    }

    auto write_cell = [&](std::size_t block, std::int64_t row) {
        const auto cb = rMesh.Cells(block);
        const NDArray& conn = cb.Conn();
        const std::size_t k = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        rOs << meshio_to_su2(cb.Type());
        for (std::size_t j = 0; j < k; ++j)
            rOs << " " << remap_of(detail::read_int(conn, static_cast<std::size_t>(row) * k + j));
        rOs << "\n";
    };
    auto tag_of = [&](std::size_t block, std::int64_t row) -> std::int64_t {
        if (rTagKey.empty())
            return 1;
        return detail::read_int(rMesh.CellData(rTagKey, block), static_cast<std::size_t>(row));
    };

    std::size_t nelem = 0;
    for (const auto& [block, row] : rCells)
        if (su2_in(vtypes, rMesh.Cells(block).Type()))
            ++nelem;
    rOs << "NELEM= " << nelem << "\n";
    for (const auto& [block, row] : rCells)
        if (su2_in(vtypes, rMesh.Cells(block).Type()))
            write_cell(block, row);

    std::map<std::int64_t, std::size_t> tag_counts;
    for (const auto& [block, row] : rCells)
        if (su2_in(btypes, rMesh.Cells(block).Type()))
            ++tag_counts[tag_of(block, row)];

    rOs << "NMARK= " << tag_counts.size() << "\n";
    for (const auto& [tag, count] : tag_counts) {
        const auto name_it = rMarkerNames.find(tag);
        rOs << "MARKER_TAG= " << (name_it != rMarkerNames.end() ? name_it->second : std::to_string(tag))
            << "\n";
        rOs << "MARKER_ELEMS= " << count << "\n";
        for (const auto& [block, row] : rCells)
            if (su2_in(btypes, rMesh.Cells(block).Type()) && tag_of(block, row) == tag)
                write_cell(block, row);
    }
}

}  // namespace

void write_su2(const std::string& rPath, const Mesh& rMesh) {
    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t dim = rMesh.PointDim();
    const std::string tag_key = su2_tag_key(rMesh);
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);

    // Multizone iff su2:zone exists and carries more than one distinct value;
    // a region "zone_<i>/<name>" (or, single-zone, "<name>") supplies a
    // marker's text where the mesh has one -- see write_su2_zone_body.
    std::map<std::int32_t, std::vector<std::pair<std::size_t, std::int64_t>>> by_zone;
    bool has_zone_data = rMesh.HasCellData("su2:zone");
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::int32_t zone = has_zone_data
                                          ? static_cast<std::int32_t>(detail::read_int(
                                                rMesh.CellData("su2:zone", bi), r))
                                          : 0;
            by_zone[zone].emplace_back(bi, static_cast<std::int64_t>(r));
        }
    }
    const bool multizone = by_zone.size() > 1;

    std::map<std::int32_t, std::map<std::int64_t, std::string>> marker_names;  // zone -> tag -> name
    for (std::size_t ri = 0; ri < rMesh.NumRegions(); ++ri) {
        const Region& region = rMesh.Region(ri);
        if (region.mKind != RegionKind::Cell || region.NumEntries() == 0)
            continue;
        const auto [block, row] = detail::global_to_block_row(bases, region.Entries()[0]);
        if (block == static_cast<std::size_t>(-1))
            continue;
        const std::int32_t zone = has_zone_data
                                      ? static_cast<std::int32_t>(
                                            detail::read_int(rMesh.CellData("su2:zone", block),
                                                             static_cast<std::size_t>(row)))
                                      : 0;
        const std::int64_t tag = tag_key.empty()
                                     ? 1
                                     : detail::read_int(rMesh.CellData(tag_key, block),
                                                        static_cast<std::size_t>(row));
        std::string name = region.mName;
        const std::string prefix = "zone_" + std::to_string(zone) + "/";
        if (multizone && name.rfind(prefix, 0) == 0)
            name = name.substr(prefix.size());
        else if (multizone)
            continue;  // a "zone_<i>" region itself, or one from another zone's write
        marker_names[zone][tag] = name;
    }

    if (!multizone) {
        std::vector<std::pair<std::size_t, std::int64_t>> all;
        for (auto& [zone, cells] : by_zone)
            all.insert(all.end(), cells.begin(), cells.end());
        write_su2_zone_body(os, rMesh, dim, all, tag_key,
                            marker_names.count(0) ? marker_names[0]
                                                  : std::map<std::int64_t, std::string>{},
                            /*Subset=*/false);
        return;
    }

    os << "NZONE= " << by_zone.size() << "\n";
    std::int32_t ordinal = 1;
    for (auto& [zone, cells] : by_zone) {
        os << "\nIZONE= " << ordinal++ << "\n";
        write_su2_zone_body(os, rMesh, dim, cells, tag_key,
                            marker_names.count(zone) ? marker_names[zone]
                                                     : std::map<std::int64_t, std::string>{},
                            /*Subset=*/true);
    }
}

}  // namespace meshioplusplus
