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
#include <cstddef>
#include <cstdint>
#include <ios>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/patran.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr std::int64_t kPatMaxId = 99999999;  // I8
constexpr std::size_t kPatMaxName = 12;       // packet 21's (A12) name card

// Packet 02's shape code (IV) and node count -> meshio++ type.
struct PatShape {
    int mShape;
    std::size_t mNodes;
    const char* mType;
    int mComponentType;  // packet 21's entity code for this shape
};

const std::vector<PatShape>& pat_shapes() {
    static const std::vector<PatShape> shapes = {
        {2, 2, "line", 6},           {2, 3, "line3", 6},     {3, 3, "triangle", 7},
        {3, 6, "triangle6", 7},      {3, 7, "triangle7", 7}, {4, 4, "quad", 8},
        {4, 8, "quad8", 8},          {4, 9, "quad9", 8},     {5, 4, "tetra", 9},
        {5, 10, "tetra10", 9},       {6, 5, "pyramid", 10},  {6, 13, "pyramid13", 10},
        {7, 6, "wedge", 11},         {7, 15, "wedge15", 11}, {8, 8, "hexahedron", 12},
        {8, 20, "hexahedron20", 12},
    };
    return shapes;
}

const PatShape* pat_shape(int Shape, std::size_t Nodes) {
    for (const PatShape& s : pat_shapes())
        if (s.mShape == Shape && s.mNodes == Nodes)
            return &s;
    return nullptr;
}

const PatShape* pat_shape_by_type(std::string_view Type) {
    for (const PatShape& s : pat_shapes())
        if (Type == s.mType)
            return &s;
    return nullptr;
}

// Packet 21 entity codes that name an element (6 bar ... 12 hex).
bool pat_is_element_code(std::int64_t Code) {
    return Code >= 6 && Code <= 12;
}

[[noreturn]] void pat_fail(const std::string& rWhat, std::size_t Line) {
    throw ReadError("Patran neutral: " + rWhat + " (line " + std::to_string(Line) + ")");
}

std::vector<std::string_view> pat_lines(const std::string& rText) {
    std::vector<std::string_view> lines;
    std::size_t pos = 0;
    while (pos < rText.size()) {
        std::size_t eol = rText.find('\n', pos);
        if (eol == std::string::npos)
            eol = rText.size();
        std::string_view line(rText.data() + pos, eol - pos);
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        lines.push_back(line);
        pos = eol + 1;
    }
    return lines;
}

std::int64_t pat_int(const std::string& rText, std::size_t Line) {
    return detail::card_to_int(rText, " (line " + std::to_string(Line) + ")", "Patran neutral");
}

double pat_real(const std::string& rText, std::size_t Line) {
    return detail::card_to_real(rText, " (line " + std::to_string(Line) + ")", "Patran neutral");
}

struct PatHeader {
    std::int64_t mIt, mId, mIv, mKc;
    std::int64_t mN[5];
};

PatHeader pat_header(std::string_view Line, std::size_t LineNo) {
    static const std::vector<detail::CardField> layout = detail::parse_fortran_format("(I2,8I8)");
    const std::vector<std::string> f = detail::split_fixed(Line, layout);
    auto at = [&](std::size_t k) { return k < f.size() ? pat_int(f[k], LineNo) : 0; };
    PatHeader h{at(0), at(1), at(2), at(3), {at(4), at(5), at(6), at(7), at(8)}};
    if (h.mKc < 0)
        pat_fail("negative card count " + std::to_string(h.mKc), LineNo);
    return h;
}

// Every `I8` field of the cards `rLines[First, First + Count)`, in order.
std::vector<std::int64_t> pat_int_cards(const std::vector<std::string_view>& rLines,
                                        std::size_t First, std::size_t Count) {
    static const std::vector<detail::CardField> layout = detail::parse_fortran_format("(10I8)");
    std::vector<std::int64_t> out;
    for (std::size_t c = 0; c < Count; ++c) {
        const std::vector<std::string> f = detail::split_fixed(rLines[First + c], layout);
        for (std::size_t k = 0; k < 10; ++k)
            out.push_back(k < f.size() ? pat_int(f[k], First + c + 1) : 0);
    }
    return out;
}

struct PatElement {
    std::int64_t mId;
    const PatShape* mShape;
    std::int64_t mPid;
    std::vector<std::int64_t> mNodes;  // node ids, file order
    std::size_t mLine;
};

struct PatComponent {
    std::int64_t mNumber;
    std::string mName;
    std::vector<std::pair<std::int64_t, std::int64_t>> mEntries;  // (code, id)
    std::size_t mLine;
};

std::string pat_trim(std::string_view Text) {
    const std::size_t b = Text.find_first_not_of(" \t");
    if (b == std::string_view::npos)
        return {};
    const std::size_t e = Text.find_last_not_of(" \t");
    return std::string(Text.substr(b, e - b + 1));
}

}  // namespace

Mesh read_patran(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Patran neutral: cannot open " + rPath);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::vector<std::string_view> lines = pat_lines(text);

    static const std::vector<detail::CardField> xyz_layout =
        detail::parse_fortran_format("(3E16.9)");
    static const std::vector<detail::CardField> elem_layout =
        detail::parse_fortran_format("(4I8,3E16.9)");

    std::vector<std::int64_t> node_ids;
    std::vector<double> coords;
    std::vector<PatElement> elements;
    std::vector<PatComponent> components;
    std::set<std::pair<std::int64_t, std::size_t>> warned_shapes;
    std::set<std::int64_t> warned_packets;
    bool saw_end = false;

    std::size_t i = 0;
    const std::size_t n = lines.size();
    while (i < n) {
        if (pat_trim(lines[i]).empty()) {
            ++i;
            continue;
        }
        const std::size_t head_line = i + 1;
        const PatHeader h = pat_header(lines[i], head_line);
        ++i;
        if (h.mIt == 99) {
            saw_end = true;
            break;
        }
        const std::size_t kc = static_cast<std::size_t>(h.mKc);
        if (i + kc > n)
            pat_fail("packet " + std::to_string(h.mIt) + " announces " + std::to_string(kc) +
                         " data cards but the file ends first",
                     head_line);
        switch (h.mIt) {
            case 1: {
                if (kc < 1)
                    pat_fail("node " + std::to_string(h.mId) + " has no coordinate card",
                             head_line);
                const std::vector<std::string> f = detail::split_fixed(lines[i], xyz_layout);
                node_ids.push_back(h.mId);
                for (std::size_t d = 0; d < 3; ++d)
                    coords.push_back(d < f.size() ? pat_real(f[d], i + 1) : 0.0);
                break;
            }
            case 2: {
                if (kc < 1)
                    pat_fail("element " + std::to_string(h.mId) + " has no data card", head_line);
                const std::vector<std::string> f = detail::split_fixed(lines[i], elem_layout);
                const std::int64_t nodes = f.empty() ? 0 : pat_int(f[0], i + 1);
                const std::int64_t pid = f.size() > 2 ? pat_int(f[2], i + 1) : 0;
                if (nodes < 0)
                    pat_fail("element " + std::to_string(h.mId) + " has a negative node count",
                             i + 1);
                const std::size_t k = static_cast<std::size_t>(nodes);
                const std::size_t node_cards = (k + 9) / 10;
                if (1 + node_cards > kc)
                    pat_fail("element " + std::to_string(h.mId) + " lists " + std::to_string(k) +
                                 " nodes but has " + std::to_string(kc) + " data cards",
                             head_line);
                const PatShape* shape = pat_shape(static_cast<int>(h.mIv), k);
                if (!shape) {
                    if (warned_shapes.emplace(h.mIv, k).second)
                        log::warn(
                            "Patran neutral: skipping elements of shape {} with {} nodes (no "
                            "meshio++ equivalent; first one is element {})",
                            h.mIv, k, h.mId);
                    break;
                }
                std::vector<std::int64_t> ids = pat_int_cards(lines, i + 1, node_cards);
                ids.resize(k);
                elements.push_back({h.mId, shape, pid, std::move(ids), head_line});
                break;
            }
            case 21: {
                if (kc < 1)
                    pat_fail("component " + std::to_string(h.mId) + " has no name card", head_line);
                PatComponent comp{h.mId, pat_trim(lines[i]), {}, head_line};
                const std::vector<std::int64_t> values = pat_int_cards(lines, i + 1, kc - 1);
                const std::size_t count = std::min<std::size_t>(
                    static_cast<std::size_t>(std::max<std::int64_t>(h.mIv, 0)), values.size());
                for (std::size_t k = 0; k + 1 < count; k += 2)
                    comp.mEntries.emplace_back(values[k], values[k + 1]);
                components.push_back(std::move(comp));
                break;
            }
            case 25:
            case 26:
                break;
            default:
                if (warned_packets.insert(h.mIt).second)
                    log::debug("Patran neutral: skipping packet type {}", h.mIt);
                break;
        }
        i += kc;
    }
    if (!saw_end)
        log::warn("Patran neutral: '{}' has no end packet (99); reading it to the end", rPath);

    // --- points -----------------------------------------------------------------
    std::unordered_map<std::int64_t, std::int64_t> node_index;
    node_index.reserve(node_ids.size());
    for (std::size_t p = 0; p < node_ids.size(); ++p)
        if (!node_index.emplace(node_ids[p], static_cast<std::int64_t>(p)).second)
            throw ReadError("Patran neutral: node " + std::to_string(node_ids[p]) +
                            " is defined twice");
    Mesh mesh;
    NDArray points(DType::Float64, {node_ids.size(), 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // --- cells: one block per type, in order of first appearance ----------------
    std::vector<const PatShape*> block_shapes;
    std::map<const PatShape*, std::vector<std::size_t>> by_shape;
    for (std::size_t e = 0; e < elements.size(); ++e) {
        auto [it, fresh] = by_shape.emplace(elements[e].mShape, std::vector<std::size_t>{});
        if (fresh)
            block_shapes.push_back(elements[e].mShape);
        it->second.push_back(e);
    }
    std::unordered_map<std::int64_t, std::int64_t> element_index;
    std::vector<std::int64_t> cell_pid;
    std::vector<int> cell_dim;
    std::vector<NDArray> pid_blocks;
    for (const PatShape* shape : block_shapes) {
        const std::vector<std::size_t>& members = by_shape[shape];
        const std::size_t k = shape->mNodes;
        const detail::NodeOrder* order = detail::node_order("patran", shape->mType);
        const int dim = cell_type_dimension(cell_type_from_name(shape->mType));
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray pids(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        std::int64_t* pp = pids.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const PatElement& el = elements[members[r]];
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mToMeshio[j]) : j;
                const auto it = node_index.find(el.mNodes[src]);
                if (it == node_index.end())
                    pat_fail("element " + std::to_string(el.mId) + " names undefined node " +
                                 std::to_string(el.mNodes[src]),
                             el.mLine);
                c[r * k + j] = it->second;
            }
            if (!element_index.emplace(el.mId, static_cast<std::int64_t>(cell_pid.size())).second)
                pat_fail("element " + std::to_string(el.mId) + " is defined twice", el.mLine);
            pp[r] = el.mPid;
            cell_pid.push_back(el.mPid);
            cell_dim.push_back(dim);
        }
        mesh.AddCellBlock(shape->mType, std::move(conn));
        pid_blocks.push_back(std::move(pids));
    }
    if (!pid_blocks.empty())
        mesh.AddCellData("patran:property", std::move(pid_blocks));

    // --- named components -------------------------------------------------------
    auto entries_of = [](const std::vector<std::int64_t>& rIds) {
        NDArray a(DType::Int64, {rIds.size()});
        std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
        return a;
    };
    std::vector<bool> named(cell_pid.size(), false);
    std::set<std::int64_t> warned_codes;
    for (const PatComponent& comp : components) {
        std::vector<std::int64_t> pts, cls;
        std::set<std::int64_t> seen_pts, seen_cls;
        int dim = -1;
        std::size_t missing = 0;
        for (const auto& [code, id] : comp.mEntries) {
            if (code == 5) {
                const auto it = node_index.find(id);
                if (it == node_index.end())
                    ++missing;
                else if (seen_pts.insert(it->second).second)
                    pts.push_back(it->second);
            } else if (pat_is_element_code(code)) {
                const auto it = element_index.find(id);
                if (it == element_index.end()) {
                    ++missing;
                } else if (seen_cls.insert(it->second).second) {
                    cls.push_back(it->second);
                    named[static_cast<std::size_t>(it->second)] = true;
                    dim = std::max(dim, cell_dim[static_cast<std::size_t>(it->second)]);
                }
            } else if (warned_codes.insert(code).second) {
                log::warn("Patran neutral: component '{}' lists entities of type {}; skipped",
                          comp.mName, code);
            }
        }
        if (missing)
            log::warn("Patran neutral: component '{}' names {} undefined node(s) or element(s)",
                      comp.mName, missing);
        if (!cls.empty() || pts.empty())
            mesh.AddRegion(
                Region(comp.mName, RegionKind::Cell, dim, comp.mNumber, entries_of(cls)));
        if (!pts.empty())
            mesh.AddRegion(
                Region(comp.mName, RegionKind::Point, -1, comp.mNumber, entries_of(pts)));
    }

    // --- property fallback for elements no component names ----------------------
    std::map<std::int64_t, std::vector<std::int64_t>> by_pid;
    std::map<std::int64_t, int> pid_dim;
    for (std::size_t g = 0; g < cell_pid.size(); ++g) {
        if (named[g])
            continue;
        by_pid[cell_pid[g]].push_back(static_cast<std::int64_t>(g));
        auto [it, fresh] = pid_dim.emplace(cell_pid[g], cell_dim[g]);
        if (!fresh)
            it->second = std::max(it->second, cell_dim[g]);
    }
    for (const auto& [pid, ids] : by_pid)
        mesh.AddRegion(Region("property_" + std::to_string(pid), RegionKind::Cell, pid_dim[pid],
                              pid, entries_of(ids)));
    return mesh;
}

namespace {

void pat_append_header(std::string& rOut, std::int64_t It, std::int64_t Id, std::int64_t Iv,
                       std::int64_t Kc, std::int64_t N1 = 0, std::int64_t N2 = 0,
                       std::int64_t N3 = 0, std::int64_t N4 = 0, std::int64_t N5 = 0) {
    char buf[96];
    detail::snprintf_c(
        buf, sizeof(buf), "%2lld%8lld%8lld%8lld%8lld%8lld%8lld%8lld%8lld\n",
        static_cast<long long>(It), static_cast<long long>(Id), static_cast<long long>(Iv),
        static_cast<long long>(Kc), static_cast<long long>(N1), static_cast<long long>(N2),
        static_cast<long long>(N3), static_cast<long long>(N4), static_cast<long long>(N5));
    rOut += buf;
}

void pat_append_ints(std::string& rOut, const std::vector<std::int64_t>& rValues) {
    char buf[16];
    for (std::size_t k = 0; k < rValues.size(); ++k) {
        detail::snprintf_c(buf, sizeof(buf), "%8lld", static_cast<long long>(rValues[k]));
        rOut += buf;
        if (k % 10 == 9 || k + 1 == rValues.size())
            rOut += '\n';
    }
}

void pat_append_real(std::string& rOut, double Value) {
    char buf[32];
    detail::snprintf_c(buf, sizeof(buf), "%16.9E", Value);
    rOut += buf;
}

// A component name: at most 12 characters, unique.
std::string pat_component_name(const std::string& rName, std::set<std::string>& rTaken) {
    std::string clean = rName.empty() ? "COMPONENT" : rName;
    for (char& c : clean)
        if (c == '\n' || c == '\r')
            c = ' ';
    if (clean.size() > kPatMaxName)
        clean.resize(kPatMaxName);
    std::string out = clean;
    for (int k = 1; rTaken.count(out); ++k) {
        const std::string suffix = "_" + std::to_string(k);
        out = clean.substr(0, std::min(clean.size(), kPatMaxName - suffix.size())) + suffix;
    }
    rTaken.insert(out);
    if (out != rName)
        log::warn("Patran neutral: component '{}' is written as '{}'", rName, out);
    return out;
}

}  // namespace

void write_patran(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t pdim = rMesh.PointDim();
    if (pdim > 3)
        throw WriteError("Patran neutral writer: points of dimension " + std::to_string(pdim) +
                         " (at most 3)");
    const std::size_t npts = rMesh.NumPoints();

    // Blocks with a Patran shape; the rest are dropped.
    std::vector<const PatShape*> shapes(rMesh.NumCellBlocks(), nullptr);
    std::vector<std::int64_t> cell_label;  // global cell -> element id, 0 if dropped
    std::size_t written = 0;
    std::set<std::string> dropped_types;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const PatShape* shape = cb.IsRagged() ? nullptr : pat_shape_by_type(cb.Type());
        shapes[b] = shape;
        if (!shape)
            dropped_types.insert(std::string(cb.Type()));
        for (std::size_t r = 0; r < cb.NumCells(); ++r)
            cell_label.push_back(shape ? static_cast<std::int64_t>(++written) : 0);
    }
    if (npts > static_cast<std::size_t>(kPatMaxId) || written > static_cast<std::size_t>(kPatMaxId))
        throw WriteError(
            "Patran neutral writer: more than 99,999,999 nodes or elements do not fit the I8 "
            "id fields");

    // Notes first: the title renders the provenance record.
    for (const std::string& t : dropped_types) {
        log::warn("Patran neutral has no '{}' element; those cells are dropped", t);
        detail::provenance_note("cells-dropped", "Patran neutral has no '" + t + "' element");
    }
    const bool has_pid = rMesh.HasCellData("patran:property");
    std::size_t side_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mKind == RegionKind::Side)
            ++side_regions;
    if (side_regions) {
        log::warn("Patran neutral has no facet components; {} side region(s) dropped",
                  side_regions);
        detail::provenance_note(
            "regions-dropped",
            std::to_string(side_regions) + " side region(s) have no Patran neutral component");
    }
    const std::size_t other_data =
        rMesh.NumPointData() + rMesh.NumFieldData() + rMesh.NumCellData() - (has_pid ? 1 : 0);
    if (other_data) {
        log::warn("Patran neutral holds no data arrays; point, cell and field data dropped");
        detail::provenance_note("data-dropped",
                                "a Patran neutral file holds no data arrays besides the "
                                "element property");
    }

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    std::string out;
    std::string title = detail::provenance_lines(detail::SlotTier::Bounded)[0];
    if (title.size() > 80)
        title.resize(80);
    pat_append_header(out, 25, 0, 0, 1);
    out += title + "\n";

    std::set<std::int64_t> pids;
    if (has_pid) {
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            if (!shapes[b])
                continue;
            const NDArray& a = rMesh.CellData("patran:property", b);
            for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r)
                pids.insert(detail::read_int(a, r));
        }
    } else if (written) {
        pids.insert(1);
    }
    pat_append_header(out, 26, 0, 0, 1, static_cast<std::int64_t>(npts),
                      static_cast<std::int64_t>(written), 0, static_cast<std::int64_t>(pids.size()),
                      0);
    out += "                                        3.0\n";

    const NDArray& points = rMesh.Points();
    for (std::size_t p = 0; p < npts; ++p) {
        pat_append_header(out, 1, static_cast<std::int64_t>(p + 1), 0, 2);
        for (std::size_t d = 0; d < 3; ++d)
            pat_append_real(out, d < pdim ? detail::read_double(points, p * pdim + d) : 0.0);
        out += "\n1G       6       0       0  000000\n";
        if (out.size() > (1u << 20)) {
            f << out;
            out.clear();
        }
    }

    std::size_t global = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const PatShape* shape = shapes[b];
        if (!shape) {
            global += cb.NumCells();
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t k = shape->mNodes;
        const detail::NodeOrder* order = detail::node_order("patran", shape->mType);
        const NDArray* pid = has_pid ? &rMesh.CellData("patran:property", b) : nullptr;
        const std::int64_t kc = 1 + static_cast<std::int64_t>((k + 9) / 10);
        std::vector<std::int64_t> ids(k);
        char buf[64];
        for (std::size_t r = 0; r < cb.NumCells(); ++r, ++global) {
            pat_append_header(out, 2, cell_label[global], shape->mShape, kc);
            detail::snprintf_c(buf, sizeof(buf), "%8zu%8d%8lld%8d", k, 0,
                               static_cast<long long>(pid ? detail::read_int(*pid, r) : 1), 0);
            out += buf;
            for (int t = 0; t < 3; ++t)
                pat_append_real(out, 0.0);
            out += '\n';
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                ids[j] = detail::read_int(conn, r * k + src) + 1;
            }
            pat_append_ints(out, ids);
            if (out.size() > (1u << 20)) {
                f << out;
                out.clear();
            }
        }
    }

    // One component per region name: a point and a cell region of one name merge.
    std::vector<std::string> names;
    std::map<std::string, std::vector<std::pair<std::int64_t, std::int64_t>>> comps;
    std::map<std::string, std::int64_t> tags;
    std::vector<int> cell_code(global, 0);
    {
        std::size_t g = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r, ++g)
                cell_code[g] = shapes[b] ? shapes[b]->mComponentType : 0;
    }
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Side)
            continue;
        auto [it, fresh] =
            comps.emplace(reg.mName, std::vector<std::pair<std::int64_t, std::int64_t>>{});
        if (fresh) {
            names.push_back(reg.mName);
            tags[reg.mName] = reg.mTag;
        }
        const std::int64_t* e = reg.Entries();
        for (std::size_t j = 0; j < reg.NumEntries(); ++j) {
            if (reg.mKind == RegionKind::Point) {
                it->second.emplace_back(5, e[j] + 1);
            } else {
                const std::size_t g = static_cast<std::size_t>(e[j]);
                if (g < cell_label.size() && cell_label[g])
                    it->second.emplace_back(cell_code[g], cell_label[g]);
            }
        }
    }
    std::set<std::int64_t> used_numbers;
    for (const std::string& name : names)
        if (tags[name] > 0)
            used_numbers.insert(tags[name]);
    std::set<std::int64_t> assigned;
    std::set<std::string> taken;
    std::int64_t next = 1;
    for (const std::string& name : names) {
        std::int64_t number = tags[name];
        if (number <= 0 || !assigned.insert(number).second) {
            while (used_numbers.count(next) || assigned.count(next))
                ++next;
            number = next;
            assigned.insert(number);
        }
        const auto& entries = comps[name];
        std::vector<std::int64_t> flat;
        flat.reserve(entries.size() * 2);
        for (const auto& [code, id] : entries) {
            flat.push_back(code);
            flat.push_back(id);
        }
        const std::int64_t iv = static_cast<std::int64_t>(flat.size());
        pat_append_header(out, 21, number, iv, 1 + (iv + 9) / 10);
        out += pat_component_name(name, taken) + "\n";
        pat_append_ints(out, flat);
    }
    pat_append_header(out, 99, 0, 0, 1);
    f << out;
    if (!f)
        throw WriteError("Patran neutral writer: failed writing " + rPath);
}

}  // namespace meshioplusplus
