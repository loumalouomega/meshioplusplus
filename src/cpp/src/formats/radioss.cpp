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
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <ios>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/radioss.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

constexpr int kRadMaxIncludeDepth = 8;

struct RadLine {
    std::string mText;
    std::string mWhere;  // "file:line"
};

std::string rad_trim(std::string_view s) {
    const std::size_t b = s.find_first_not_of(" \t\r");
    if (b == std::string_view::npos)
        return {};
    const std::size_t e = s.find_last_not_of(" \t\r");
    return std::string(s.substr(b, e - b + 1));
}

std::string rad_upper(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// The deck's lines with `#include` files inlined and comments dropped, up to
// `/END` (or `#enddata` in an include).
void rad_collect(const fs::path& rPath, int Depth, std::vector<RadLine>& rOut, bool& rEnded) {
    if (Depth > kRadMaxIncludeDepth)
        throw ReadError("Radioss: #include nested deeper than " +
                        std::to_string(kRadMaxIncludeDepth));
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Radioss: cannot open " + rPath.string());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string label = rPath.filename().string();
    std::size_t pos = 0, number = 0;
    while (pos < text.size() && !rEnded) {
        std::size_t eol = text.find('\n', pos);
        if (eol == std::string::npos)
            eol = text.size();
        std::string line = text.substr(pos, eol - pos);
        pos = eol + 1;
        ++number;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        const std::string lower = rad_upper(line.substr(0, 9));
        if (lower.rfind("#INCLUDE", 0) == 0) {
            std::string name = rad_trim(std::string_view(line).substr(8));
            if (name.size() >= 2 && (name.front() == '"' || name.front() == '\'') &&
                name.back() == name.front())
                name = name.substr(1, name.size() - 2);
            std::replace(name.begin(), name.end(), '\\', '/');
            const fs::path inc =
                fs::path(name).is_absolute() ? fs::path(name) : rPath.parent_path() / name;
            std::error_code ec;
            if (!fs::is_regular_file(inc, ec)) {
                log::warn("Radioss: include file '{}' not found; skipped", name);
                continue;
            }
            bool sub_ended = false;
            rad_collect(inc, Depth + 1, rOut, sub_ended);
            continue;
        }
        if (Depth > 0 && rad_upper(rad_trim(line)).rfind("#ENDDATA", 0) == 0)
            return;
        // Blank lines stay: a blank /PART or group title is still its title line.
        if (!line.empty() && (line[0] == '#' || line[0] == '$'))
            continue;
        if (rad_upper(rad_trim(line)) == "/END") {
            if (Depth == 0)
                rEnded = true;
            if (Depth == 0)
                return;
            continue;
        }
        rOut.push_back({std::move(line), label + ":" + std::to_string(number)});
    }
}

// Fixed fields of `Width` columns, or comma-separated values.
std::vector<std::string> rad_fields(const std::string& rLine, int Width, std::size_t Count) {
    std::vector<std::string> out;
    if (rLine.find(',') != std::string::npos) {
        std::size_t start = 0;
        while (true) {
            const std::size_t k = rLine.find(',', start);
            out.push_back(rad_trim(std::string_view(rLine).substr(
                start, k == std::string::npos ? std::string::npos : k - start)));
            if (k == std::string::npos)
                break;
            start = k + 1;
        }
        return out;
    }
    const std::size_t w = static_cast<std::size_t>(Width);
    for (std::size_t f = 0; f < Count; ++f) {
        const std::size_t at = f * w;
        out.push_back(at < rLine.size() ? rad_trim(std::string_view(rLine).substr(at, w))
                                        : std::string());
    }
    return out;
}

std::int64_t rad_int(const std::string& rText, const RadLine& rLine) {
    return detail::card_to_int(rText, " (" + rLine.mWhere + ")", "Radioss");
}

double rad_real(const std::string& rText, const RadLine& rLine) {
    return detail::card_to_real(rText, " (" + rLine.mWhere + ")", "Radioss");
}

// Element keywords: their group family, node count, lines per element and the
// meshio++ type they read as (before degenerate-brick and triangle checks).
struct RadElementKind {
    const char* mKeyword;
    const char* mFamily;  // the /GR<family> keyword and /GRxxx/<subtype> that list them
    std::size_t mNodes;
    const char* mType;
};

const RadElementKind* rad_element_kind(const std::string& rKeyword) {
    static const RadElementKind kinds[] = {
        {"BRICK", "BRIC", 8, "hexahedron"},     {"PENTA6", "BRIC", 6, "wedge"},
        {"TETRA4", "BRIC", 4, "tetra"},         {"TETRA10", "BRIC", 10, "tetra10"},
        {"BRIC20", "BRIC", 20, "hexahedron20"}, {"SHELL", "SHEL", 4, "quad"},
        {"SH3N", "SH3N", 3, "triangle"},        {"QUAD", "QUAD", 4, "quad"},
        {"TRIA", "TRIA", 3, "triangle"},        {"BEAM", "BEAM", 2, "line"},
        {"TRUSS", "TRUS", 2, "line"},           {"SPRING", "SPRI", 2, "line"},
    };
    for (const RadElementKind& k : kinds)
        if (rKeyword == k.mKeyword)
            return &k;
    return nullptr;
}

// The /GRxxx keyword -> its element family (GRNOD -> "NODE").
const char* rad_group_family(const std::string& rKeyword) {
    static const std::pair<const char*, const char*> groups[] = {
        {"GRNOD", "NODE"},  {"GRBRIC", "BRIC"}, {"GRSHEL", "SHEL"},
        {"GRSH3N", "SH3N"}, {"GRQUAD", "QUAD"}, {"GRTRIA", "TRIA"},
        {"GRBEAM", "BEAM"}, {"GRTRUS", "TRUS"}, {"GRSPRI", "SPRI"},
    };
    for (const auto& [k, f] : groups)
        if (rKeyword == k)
            return f;
    return nullptr;
}

struct RadElement {
    std::string mFamily;
    std::int64_t mId, mPart;
    std::string mType;
    std::vector<std::int64_t> mNodes;  // node ids, meshio++ order
    std::string mWhere;
    bool mTetra = false;  // a /TETRA4 or /TETRA10, reoriented when inverted
};

struct RadPart {
    std::string mTitle;
    std::int64_t mProperty = 0, mMaterial = 0, mSubset = 0;
};

struct RadGroup {
    std::string mKeyword, mSubtype, mTitle;
    std::int64_t mId;
    std::vector<std::int64_t> mIds;
};

struct RadSurface {
    std::int64_t mId;
    std::string mTitle;
    std::vector<std::array<std::int64_t, 4>> mSegments;
};

struct RadSubset {
    std::string mTitle;
    std::vector<std::int64_t> mChildren;
};

NDArray rad_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride = 1) {
    NDArray a(DType::Int64, Stride == 1 ? std::vector<std::size_t>{rIds.size()}
                                        : std::vector<std::size_t>{rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

Mesh read_radioss(const std::string& rPath) {
    std::vector<RadLine> lines;
    bool ended = false;
    {
        auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
        if (!in)
            throw ReadError("Radioss: cannot open " + rPath);
        std::string head(512, '\0');
        in.read(head.data(), static_cast<std::streamsize>(head.size()));
        head.resize(static_cast<std::size_t>(in.gcount()));
        if (rad_upper(head).find("#RADIOSS ENGINE") != std::string::npos)
            throw ReadError("Radioss: " + rPath +
                            " is an engine deck (_0001.rad); read the starter deck (_0000.rad)");
    }
    rad_collect(fs::path(rPath), 0, lines, ended);

    int version = 2019;
    int iw = 10, rw = 20;
    std::vector<std::int64_t> node_ids;
    std::vector<double> coords;
    std::vector<RadElement> elements;
    std::map<std::int64_t, RadPart> parts;
    std::vector<std::int64_t> part_order;
    std::vector<RadGroup> groups;
    std::vector<RadSurface> surfaces;
    std::map<std::int64_t, RadSubset> subsets;
    std::set<std::string> skipped_keywords;
    std::size_t zero_springs = 0, linear_bric20 = 0;

    std::size_t i = 0;
    const std::size_t n = lines.size();
    auto block_end = [&](std::size_t From) {
        std::size_t k = From;
        while (k < n && (lines[k].mText.empty() || lines[k].mText[0] != '/'))
            ++k;
        return k;
    };
    while (i < n) {
        const RadLine& head = lines[i];
        if (head.mText.empty() || head.mText[0] != '/') {
            ++i;
            continue;
        }
        std::vector<std::string> path;
        {
            std::string kw = rad_trim(head.mText);
            std::size_t start = 1;
            while (start <= kw.size()) {
                const std::size_t k = kw.find('/', start);
                path.push_back(rad_upper(rad_trim(std::string_view(kw).substr(
                    start, k == std::string::npos ? std::string::npos : k - start))));
                if (k == std::string::npos)
                    break;
                start = k + 1;
            }
        }
        const std::string& key = path.empty() ? std::string() : path[0];
        const std::size_t body = i + 1;
        const std::size_t end = block_end(body);
        auto last_id = [&](const RadLine& rLine) -> std::int64_t {
            return path.size() >= 2 ? rad_int(path.back(), rLine) : 0;
        };

        if (key == "BEGIN") {
            // run name; Invers Irun; two unit lines.
            if (body + 1 < end) {
                const std::vector<std::string> f = rad_fields(lines[body + 1].mText, 10, 2);
                if (!f.empty() && !f[0].empty())
                    version = static_cast<int>(rad_int(f[0], lines[body + 1]));
            }
            iw = version >= 51 ? 10 : 8;
            rw = version >= 51 ? 20 : 16;
        } else if (key == "NODE") {
            for (std::size_t k = body; k < end; ++k) {
                const RadLine& ln = lines[k];
                std::vector<std::string> f;
                if (ln.mText.find(',') != std::string::npos) {
                    f = rad_fields(ln.mText, iw, 4);
                } else {
                    f.push_back(rad_trim(ln.mText.substr(
                        0, std::min<std::size_t>(ln.mText.size(), static_cast<std::size_t>(iw)))));
                    for (std::size_t d = 0; d < 3; ++d) {
                        const std::size_t at =
                            static_cast<std::size_t>(iw + rw * static_cast<int>(d));
                        f.push_back(at < ln.mText.size()
                                        ? rad_trim(std::string_view(ln.mText).substr(
                                              at, static_cast<std::size_t>(rw)))
                                        : std::string());
                    }
                }
                if (f.empty() || f[0].empty())
                    continue;
                node_ids.push_back(rad_int(f[0], ln));
                for (std::size_t d = 1; d <= 3; ++d)
                    coords.push_back(d < f.size() && !f[d].empty() ? rad_real(f[d], ln) : 0.0);
            }
        } else if (const RadElementKind* kind = rad_element_kind(key)) {
            const std::int64_t part = last_id(head);
            // Every element's integers, gathered from its line(s).
            const std::size_t per_record = key == "TETRA10" ? 11 : kind->mNodes + 1;
            std::vector<std::int64_t> values;
            std::string where;
            auto flush = [&]() {
                RadElement el{kind->mFamily, values[0], part, kind->mType, {}, where};
                std::vector<std::int64_t> nodes(
                    values.begin() + 1,
                    values.begin() + 1 + static_cast<std::ptrdiff_t>(kind->mNodes));
                if (key == "BRICK") {
                    std::array<std::int64_t, 8> a;
                    std::copy(nodes.begin(), nodes.end(), a.begin());
                    detail::CollapsedBrick c = detail::collapse_brick(a);
                    el.mType = c.mType;
                    nodes = std::move(c.mNodes);
                } else if (key == "SHELL" && (nodes[3] == nodes[2] || nodes[3] == 0)) {
                    el.mType = "triangle";
                    nodes.resize(3);
                } else if (key == "BRIC20") {
                    if (std::find(nodes.begin() + 8, nodes.end(), 0) != nodes.end()) {
                        ++linear_bric20;
                        el.mType = "hexahedron";
                        nodes.resize(8);
                    } else {
                        const detail::NodeOrder* order =
                            detail::node_order("radioss", "hexahedron20");
                        std::vector<std::int64_t> ordered(20);
                        for (std::size_t j = 0; j < 20; ++j)
                            ordered[j] =
                                nodes[order ? static_cast<std::size_t>(order->mToMeshio[j]) : j];
                        nodes = std::move(ordered);
                    }
                } else if (key == "SPRING" && nodes[1] == 0) {
                    ++zero_springs;
                    values.clear();
                    return;
                }
                el.mNodes = std::move(nodes);
                el.mTetra = key == "TETRA4" || key == "TETRA10";
                elements.push_back(std::move(el));
                values.clear();
            };
            for (std::size_t k = body; k < end; ++k) {
                const RadLine& ln = lines[k];
                if (rad_trim(ln.mText).empty())
                    continue;
                // Integer fields of the line: BRIC20 and TETRA10 records span
                // several lines; the others are one line each.
                std::size_t take = 10;
                if (key == "TETRA10")
                    take = values.empty() ? 1 : 10;
                else if (key == "BRIC20")
                    take = values.empty() ? 9 : (values.size() == 9 ? 8 : 4);
                else
                    take = per_record;
                const std::vector<std::string> f = rad_fields(ln.mText, iw, take);
                if (values.empty())
                    where = ln.mWhere;
                for (std::size_t j = 0; j < take && values.size() < per_record; ++j)
                    values.push_back(j < f.size() && !f[j].empty() ? rad_int(f[j], ln) : 0);
                if (values.size() >= per_record || (key != "TETRA10" && key != "BRIC20"))
                    flush();
            }
            if (!values.empty())
                throw ReadError("Radioss: /" + key + " element " + std::to_string(values[0]) +
                                " is cut short (" + where + ")");
        } else if (key == "PART") {
            const std::int64_t id = last_id(head);
            RadPart part;
            std::size_t k = body;
            if (k < end)
                part.mTitle = rad_trim(lines[k++].mText);
            if (k < end) {
                const std::vector<std::string> f = rad_fields(lines[k].mText, iw, 3);
                part.mProperty = !f.empty() && !f[0].empty() ? rad_int(f[0], lines[k]) : 0;
                part.mMaterial = f.size() > 1 && !f[1].empty() ? rad_int(f[1], lines[k]) : 0;
                part.mSubset = f.size() > 2 && !f[2].empty() ? rad_int(f[2], lines[k]) : 0;
            }
            if (parts.emplace(id, part).second)
                part_order.push_back(id);
        } else if (rad_group_family(key) && path.size() >= 3) {
            RadGroup g{key, path[1], {}, last_id(head), {}};
            std::size_t k = body;
            if (k < end)
                g.mTitle = rad_trim(lines[k++].mText);
            for (; k < end; ++k)
                for (const std::string& f : rad_fields(lines[k].mText, iw, 10))
                    if (!f.empty())
                        g.mIds.push_back(rad_int(f, lines[k]));
            groups.push_back(std::move(g));
        } else if (key == "SURF" && path.size() >= 3 && path[1] == "SEG") {
            RadSurface s{last_id(head), {}, {}};
            std::size_t k = body;
            if (k < end)
                s.mTitle = rad_trim(lines[k++].mText);
            for (; k < end; ++k) {
                const std::vector<std::string> f = rad_fields(lines[k].mText, iw, 5);
                std::array<std::int64_t, 4> seg{0, 0, 0, 0};
                for (std::size_t j = 0; j < 4; ++j)
                    seg[j] =
                        j + 1 < f.size() && !f[j + 1].empty() ? rad_int(f[j + 1], lines[k]) : 0;
                if (seg[0] || seg[1] || seg[2])
                    s.mSegments.push_back(seg);
            }
            surfaces.push_back(std::move(s));
        } else if (key == "SURF") {
            skipped_keywords.insert("/SURF/" + (path.size() > 1 ? path[1] : std::string()));
        } else if (key == "SUBSET") {
            RadSubset s;
            std::size_t k = body;
            if (k < end)
                s.mTitle = rad_trim(lines[k++].mText);
            for (; k < end; ++k)
                for (const std::string& f : rad_fields(lines[k].mText, iw, 10))
                    if (!f.empty())
                        s.mChildren.push_back(rad_int(f, lines[k]));
            subsets[last_id(head)] = std::move(s);
        } else if (key.rfind("GR", 0) == 0 || key == "TSHELL" || key == "TSH3N" ||
                   key == "SHEL16" || key == "SPHCEL" || key == "RIVET" || key == "XELEM") {
            skipped_keywords.insert("/" + key);
        }
        i = end;
    }
    if (!skipped_keywords.empty()) {
        std::string list;
        for (const std::string& k : skipped_keywords)
            list += (list.empty() ? "" : ", ") + k;
        log::warn("Radioss: mesh keywords not read: {}", list);
    }
    if (zero_springs)
        log::warn("Radioss: {} /SPRING element(s) with a single node skipped", zero_springs);
    if (linear_bric20)
        log::warn(
            "Radioss: {} /BRIC20 element(s) with missing mid-edge nodes read as "
            "hexahedra",
            linear_bric20);

    // --- points -------------------------------------------------------------------
    Mesh mesh;
    std::unordered_map<std::int64_t, std::int64_t> node_index;
    node_index.reserve(node_ids.size());
    for (std::size_t p = 0; p < node_ids.size(); ++p)
        if (!node_index.emplace(node_ids[p], static_cast<std::int64_t>(p)).second)
            throw ReadError("Radioss: node " + std::to_string(node_ids[p]) + " is defined twice");
    NDArray points(DType::Float64, {node_ids.size(), 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    mesh.AddFieldData("radioss:version", [&] {
        NDArray a(DType::Int64, {});
        a.As<std::int64_t>()[0] = version;
        return a;
    }());

    // --- cells: one block per type, in order of first appearance --------------------
    std::vector<std::string> block_types;
    std::map<std::string, std::vector<std::size_t>> by_type;
    for (std::size_t e = 0; e < elements.size(); ++e) {
        auto [it, fresh] = by_type.emplace(elements[e].mType, std::vector<std::size_t>{});
        if (fresh)
            block_types.push_back(elements[e].mType);
        it->second.push_back(e);
    }
    std::map<std::string, std::unordered_map<std::int64_t, std::int64_t>>
        owner;  // family -> id -> cell
    std::map<std::int64_t, std::vector<std::int64_t>> part_cells;
    std::vector<std::int64_t> cell_part;
    std::vector<int> cell_dim;
    std::vector<NDArray> part_blocks, prop_blocks, mat_blocks;
    // Tetrahedra come in either winding (every /TETRA4 of OpenRadioss's INT_25 QA
    // deck is inverted, gmsh writes them positive): an inverted one is mirrored.
    std::size_t reoriented = 0;
    for (RadElement& el : elements) {
        if (!el.mTetra)
            continue;
        std::array<std::array<double, 3>, 4> p{};
        bool known = true;
        for (std::size_t c = 0; c < 4 && known; ++c) {
            const auto it = node_index.find(el.mNodes[c]);
            known = it != node_index.end();
            if (known)
                for (std::size_t d = 0; d < 3; ++d)
                    p[c][d] = coords[3 * static_cast<std::size_t>(it->second) + d];
        }
        if (!known)
            continue;
        double a[3], b[3], c[3];
        for (std::size_t d = 0; d < 3; ++d) {
            a[d] = p[1][d] - p[0][d];
            b[d] = p[2][d] - p[0][d];
            c[d] = p[3][d] - p[0][d];
        }
        const double det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                           a[2] * (b[0] * c[1] - b[1] * c[0]);
        if (det >= 0.0)
            continue;
        std::swap(el.mNodes[1], el.mNodes[2]);
        if (el.mNodes.size() == 10) {  // edges 0-1 <-> 0-2 and 1-3 <-> 2-3
            std::swap(el.mNodes[4], el.mNodes[6]);
            std::swap(el.mNodes[8], el.mNodes[9]);
        }
        ++reoriented;
    }
    if (reoriented)
        log::warn("Radioss: {} inverted tetrahedra reoriented", reoriented);

    for (const std::string& type : block_types) {
        const std::vector<std::size_t>& members = by_type[type];
        const std::size_t k =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(type)));
        const int dim = cell_type_dimension(cell_type_from_name(type));
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray pa(DType::Int64, {members.size()}), pr(DType::Int64, {members.size()}),
            ma(DType::Int64, {members.size()});
        for (std::size_t r = 0; r < members.size(); ++r) {
            const RadElement& el = elements[members[r]];
            for (std::size_t j = 0; j < k; ++j) {
                const auto it = node_index.find(el.mNodes[j]);
                if (it == node_index.end())
                    throw ReadError("Radioss: element " + std::to_string(el.mId) +
                                    " names undefined node " + std::to_string(el.mNodes[j]) + " (" +
                                    el.mWhere + ")");
                conn.As<std::int64_t>()[r * k + j] = it->second;
            }
            const std::int64_t cell = static_cast<std::int64_t>(cell_part.size());
            if (!owner[el.mFamily].emplace(el.mId, cell).second)
                throw ReadError("Radioss: element " + std::to_string(el.mId) +
                                " is defined twice (" + el.mWhere + ")");
            const auto p = parts.find(el.mPart);
            pa.As<std::int64_t>()[r] = el.mPart;
            pr.As<std::int64_t>()[r] = p == parts.end() ? 0 : p->second.mProperty;
            ma.As<std::int64_t>()[r] = p == parts.end() ? 0 : p->second.mMaterial;
            part_cells[el.mPart].push_back(cell);
            cell_part.push_back(el.mPart);
            cell_dim.push_back(dim);
        }
        mesh.AddCellBlock(type, std::move(conn));
        part_blocks.push_back(std::move(pa));
        prop_blocks.push_back(std::move(pr));
        mat_blocks.push_back(std::move(ma));
    }
    if (!block_types.empty()) {
        mesh.AddCellData("radioss:part", std::move(part_blocks));
        mesh.AddCellData("radioss:property", std::move(prop_blocks));
        mesh.AddCellData("radioss:material", std::move(mat_blocks));
    }

    // --- regions ------------------------------------------------------------------
    std::set<std::tuple<int, std::string, std::int64_t>> seen;
    auto add_region = [&](std::string name, RegionKind kind, std::int64_t tag,
                          const std::string& rKeyword, std::vector<std::int64_t> entries) {
        int dim = -1;
        if (kind == RegionKind::Cell)
            for (std::int64_t c : entries)
                dim = std::max(dim, cell_dim[static_cast<std::size_t>(c)]);
        else if (kind == RegionKind::Side)
            for (std::size_t k = 0; k < entries.size(); k += 2)
                dim = std::max(dim, cell_dim[static_cast<std::size_t>(entries[k])] - 1);
        const int kid = static_cast<int>(kind);
        if (seen.count({kid, name, tag}))
            name += " [" + rKeyword + "]";
        seen.emplace(kid, name, tag);
        mesh.AddRegion(
            Region(name, kind, dim, tag, rad_ids(entries, kind == RegionKind::Side ? 2 : 1)));
    };
    // Parts named in elements but never defined still get a region.
    for (const auto& [pid, cells] : part_cells)
        if (!parts.count(pid)) {
            parts.emplace(pid, RadPart{});
            part_order.push_back(pid);
        }
    for (std::int64_t pid : part_order) {
        const RadPart& p = parts[pid];
        const auto it = part_cells.find(pid);
        add_region(p.mTitle.empty() ? "Part " + std::to_string(pid) : p.mTitle, RegionKind::Cell,
                   pid, "PART", it == part_cells.end() ? std::vector<std::int64_t>{} : it->second);
    }
    // A subset and every subset under it.
    auto subset_closure = [&](std::int64_t Root) {
        std::set<std::int64_t> ids;
        std::vector<std::int64_t> stack{Root};
        while (!stack.empty()) {
            const std::int64_t c = stack.back();
            stack.pop_back();
            if (!ids.insert(c).second)
                continue;
            const auto it = subsets.find(c);
            if (it != subsets.end())
                stack.insert(stack.end(), it->second.mChildren.begin(), it->second.mChildren.end());
        }
        return ids;
    };
    for (const auto& [sid, s] : subsets) {
        const std::set<std::int64_t> ids = subset_closure(sid);
        std::vector<std::int64_t> cells;
        for (std::int64_t pid : part_order)
            if (ids.count(parts[pid].mSubset)) {
                const auto it = part_cells.find(pid);
                if (it != part_cells.end())
                    cells.insert(cells.end(), it->second.begin(), it->second.end());
            }
        add_region(s.mTitle.empty() ? "Subset " + std::to_string(sid) : s.mTitle, RegionKind::Cell,
                   sid, "SUBSET", std::move(cells));
    }

    // Groups: resolved recursively (a group can list other groups).
    std::map<std::pair<std::string, std::int64_t>, std::size_t> group_of;
    for (std::size_t g = 0; g < groups.size(); ++g)
        group_of[{groups[g].mKeyword, groups[g].mId}] = g;
    std::set<std::string> skipped_subtypes;
    std::size_t dropped = 0;
    // A PART subtype stands for the parts' nodes in a GRNOD, and for the parts'
    // cells of the group's own family in an element group.
    std::vector<std::string> cell_family(cell_part.size());
    for (const auto& [family, map] : owner)
        for (const auto& [id, cell] : map)
            cell_family[static_cast<std::size_t>(cell)] = family;
    std::vector<std::vector<std::int64_t>> cell_conn;
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
        const auto cb = mesh.Cells(b);
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r)
            cell_conn.emplace_back(conn.As<std::int64_t>() + r * k,
                                   conn.As<std::int64_t>() + (r + 1) * k);
    }
    std::function<std::set<std::int64_t>(std::size_t, std::set<std::size_t>&)> resolve;
    resolve = [&](std::size_t g, std::set<std::size_t>& rVisiting) -> std::set<std::int64_t> {
        std::set<std::int64_t> out;
        if (!rVisiting.insert(g).second)
            return out;
        const RadGroup& grp = groups[g];
        const std::string family = rad_group_family(grp.mKeyword);
        const bool nodes = family == "NODE";
        std::set<std::int64_t> removed;
        for (std::int64_t raw : grp.mIds) {
            const std::int64_t id = raw < 0 ? -raw : raw;
            std::set<std::int64_t> hits;
            if (grp.mSubtype == family) {
                if (nodes) {
                    const auto it = node_index.find(id);
                    if (it != node_index.end())
                        hits.insert(it->second);
                    else
                        ++dropped;
                } else {
                    const auto& map = owner[family];
                    const auto it = map.find(id);
                    if (it != map.end())
                        hits.insert(it->second);
                    else
                        ++dropped;
                }
            } else if (grp.mSubtype == "PART" || grp.mSubtype == "SUBSET") {
                // The parts' nodes for GRNOD, else their cells of the group's family.
                std::vector<std::int64_t> pids{id};
                if (grp.mSubtype == "SUBSET") {
                    pids.clear();
                    for (std::int64_t pid : part_order)
                        if (subset_closure(id).count(parts[pid].mSubset))
                            pids.push_back(pid);
                }
                for (std::int64_t pid : pids) {
                    const auto it = part_cells.find(pid);
                    if (it == part_cells.end())
                        continue;
                    for (std::int64_t c : it->second) {
                        if (nodes)
                            hits.insert(cell_conn[static_cast<std::size_t>(c)].begin(),
                                        cell_conn[static_cast<std::size_t>(c)].end());
                        else if (cell_family[static_cast<std::size_t>(c)] == family)
                            hits.insert(c);
                    }
                }
            } else if (nodes && grp.mSubtype == "SURF") {
                for (const RadSurface& surf : surfaces)
                    if (surf.mId == id)
                        for (const auto& seg : surf.mSegments)
                            for (std::int64_t v : seg) {
                                const auto it = node_index.find(v);
                                if (v && it != node_index.end())
                                    hits.insert(it->second);
                            }
            } else if (nodes && rad_group_family(grp.mSubtype) &&
                       std::string(rad_group_family(grp.mSubtype)) != "NODE") {
                // A GRNOD of element groups: those elements' nodes.
                const auto it = group_of.find({grp.mSubtype, id});
                if (it != group_of.end())
                    for (std::int64_t c : resolve(it->second, rVisiting))
                        hits.insert(cell_conn[static_cast<std::size_t>(c)].begin(),
                                    cell_conn[static_cast<std::size_t>(c)].end());
            } else if (grp.mSubtype == grp.mKeyword) {
                const auto it = group_of.find({grp.mKeyword, id});
                if (it != group_of.end()) {
                    const std::set<std::int64_t> sub = resolve(it->second, rVisiting);
                    hits.insert(sub.begin(), sub.end());
                }
            } else {
                skipped_subtypes.insert("/" + grp.mKeyword + "/" + grp.mSubtype);
                break;
            }
            (raw < 0 ? removed : out).insert(hits.begin(), hits.end());
        }
        for (std::int64_t r : removed)
            out.erase(r);
        rVisiting.erase(g);
        return out;
    };
    for (std::size_t g = 0; g < groups.size(); ++g) {
        std::set<std::size_t> visiting;
        const RadGroup& grp = groups[g];
        const std::string family = rad_group_family(grp.mKeyword);
        const std::set<std::int64_t> ids = resolve(g, visiting);
        std::vector<std::int64_t> entries(ids.begin(), ids.end());
        add_region(grp.mTitle.empty() ? grp.mKeyword + "_" + std::to_string(grp.mId) : grp.mTitle,
                   family == "NODE" ? RegionKind::Point : RegionKind::Cell, grp.mId, grp.mKeyword,
                   std::move(entries));
    }
    for (const std::string& t : skipped_subtypes)
        log::warn(
            "Radioss: {} groups are not resolved (only entity ids, parts and other groups "
            "are); their regions are left empty or partial",
            t);

    // Surfaces -> side regions; a shell segment is the shell's own face.
    if (!surfaces.empty()) {
        detail::FacetIndexOptions options;
        options.mSurfaceSelf = true;
        const detail::FacetIndex facets(mesh, options);
        for (const RadSurface& s : surfaces) {
            std::vector<std::int64_t> entries;
            for (const auto& seg : s.mSegments) {
                std::array<std::int64_t, 4> idx{};
                bool defined = true;
                // 2 nodes in a 2-D analysis, 3 for a triangle (n4 blank or n3).
                const std::size_t count =
                    seg[2] == 0 ? 2 : ((seg[3] == 0 || seg[3] == seg[2]) ? 3 : 4);
                for (std::size_t k = 0; k < count; ++k) {
                    const auto it = node_index.find(seg[k]);
                    defined = defined && it != node_index.end();
                    idx[k] = defined ? it->second : -1;
                }
                const detail::FacetHit* hit = defined ? facets.Find(idx.data(), count) : nullptr;
                if (!hit) {
                    ++dropped;
                    continue;
                }
                entries.push_back(hit->mFirst.mCell);
                entries.push_back(hit->mFirst.mFacet);
            }
            add_region(s.mTitle.empty() ? "SURF_" + std::to_string(s.mId) : s.mTitle,
                       RegionKind::Side, s.mId, "SURF", std::move(entries));
        }
    }
    if (dropped)
        log::warn(
            "Radioss: {} group or surface entries name undefined ids or no cell facet and "
            "were dropped",
            dropped);
    return mesh;
}

}  // namespace meshioplusplus
