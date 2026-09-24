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
// FEBio `.feb` mesh reader (spec 2.5, 3.0, 4.0) and spec-4.0 writer. The tag
// layouts follow FEBio's own XML parsers (FEBioXML/FEBioGeometrySection.cpp,
// FEBioMeshSection.cpp, FEBioMeshSection4.cpp, FEBioMeshDataSection*.cpp); see
// doc/formats/febio.md. Python twin: src/python/meshioplusplus/febio/_febio.py.

// System includes
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/febio.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

// FEBio element type -> meshio++ cell type and node count. Integration-rule
// aliases (`TET10G4`, `HEX8G1`, ...) and the shell formulations (`tri3s`,
// `q4eas`, ...) name the same nodes.
struct FebType {
    const char* mFebio;
    const char* mType;
    std::size_t mNodes;
};

constexpr FebType kFebTypes[] = {
    {"tet4", "tetra", 4},          {"ut4", "tetra", 4},        {"tet10", "tetra10", 10},
    {"penta6", "wedge", 6},        {"penta15", "wedge15", 15}, {"pyra5", "pyramid", 5},
    {"pyra13", "pyramid13", 13},   {"hex8", "hexahedron", 8},  {"hex20", "hexahedron20", 20},
    {"hex27", "hexahedron27", 27}, {"quad4", "quad", 4},       {"q4eas", "quad", 4},
    {"q4ans", "quad", 4},          {"q4s", "quad", 4},         {"quad8", "quad8", 8},
    {"quad9", "quad9", 9},         {"tri3", "triangle", 3},    {"tri3s", "triangle", 3},
    {"tri6", "triangle6", 6},      {"tri7", "triangle7", 7},   {"tri10", "triangle10", 10},
    {"line2", "line", 2},          {"truss2", "line", 2},      {"line3", "line3", 3},
};

// Types FEBio reads that meshio++ has no cell for: the node count, and the
// meshio++ type `mLenient` downgrades to by keeping the leading nodes.
struct FebLossy {
    const char* mFebio;
    std::size_t mNodes;
    const char* mDowngrade;
};

constexpr FebLossy kFebLossy[] = {
    {"tet5", 5, "tetra"}, {"tet15", 15, "tetra10"}, {"tet20", 20, nullptr}};

std::string feb_lower(std::string_view Text) {
    std::string out(Text);
    for (char& c : out)
        if (c >= 'A' && c <= 'Z')
            c = static_cast<char>(c - 'A' + 'a');
    return out;
}

// `TET10G4` -> `tet10`: an integration-rule suffix `g<digits>` after the node count.
std::string feb_base_type(std::string_view Type) {
    std::string t = feb_lower(Type);
    const std::size_t g = t.rfind('g');
    if (g != std::string::npos && g > 0 && g + 1 < t.size() && t[g - 1] >= '0' && t[g - 1] <= '9' &&
        std::all_of(t.begin() + static_cast<std::ptrdiff_t>(g + 1), t.end(),
                    [](char c) { return c >= '0' && c <= '9'; }))
        t.resize(g);
    return t;
}

const FebType* feb_type(std::string_view Type) {
    const std::string t = feb_base_type(Type);
    for (const FebType& e : kFebTypes)
        if (t == e.mFebio)
            return &e;
    return nullptr;
}

// meshio++ cell type -> the FEBio element type the writer uses.
const char* feb_write_type(std::string_view Type) {
    static constexpr std::pair<std::string_view, const char*> kWrite[] = {
        {"line", "line2"},         {"line3", "line3"},        {"triangle", "tri3"},
        {"triangle6", "tri6"},     {"triangle7", "tri7"},     {"quad", "quad4"},
        {"quad8", "quad8"},        {"quad9", "quad9"},        {"tetra", "tet4"},
        {"tetra10", "tet10"},      {"pyramid", "pyra5"},      {"pyramid13", "pyra13"},
        {"wedge", "penta6"},       {"wedge15", "penta15"},    {"hexahedron", "hex8"},
        {"hexahedron20", "hex20"}, {"hexahedron27", "hex27"},
    };
    for (const auto& [type, febio] : kWrite)
        if (Type == type)
            return febio;
    return nullptr;
}

std::size_t feb_num_corners(std::string_view Type) {
    static constexpr std::pair<std::string_view, std::size_t> kFamilies[] = {
        {"vertex", 1}, {"line", 2},    {"triangle", 3}, {"quad", 4},
        {"tetra", 4},  {"pyramid", 5}, {"wedge", 6},    {"hexahedron", 8},
    };
    for (const auto& [prefix, n] : kFamilies)
        if (Type.substr(0, prefix.size()) == prefix)
            return n;
    return 0;
}

[[noreturn]] void feb_fail(const std::string& rWhat) {
    throw ReadError("FEBio .feb: " + rWhat);
}

// Comma- and/or blank-separated values of an element's text.
void feb_split(std::string_view Text, std::vector<std::string_view>& rOut) {
    rOut.clear();
    std::size_t i = 0;
    const auto sep = [](char c) {
        return c == ',' || c == ' ' || c == '\t' || c == '\n' || c == '\r';
    };
    while (i < Text.size()) {
        while (i < Text.size() && sep(Text[i]))
            ++i;
        const std::size_t start = i;
        while (i < Text.size() && !sep(Text[i]))
            ++i;
        if (i > start)
            rOut.push_back(Text.substr(start, i - start));
    }
}

std::optional<std::int64_t> feb_int(std::string_view Token) {
    std::int64_t v = 0;
    const char* end = Token.data() + Token.size();
    const auto [ptr, ec] = std::from_chars(Token.data(), end, v);
    if (ec != std::errc() || ptr != end)
        return std::nullopt;
    return v;
}

std::int64_t feb_need_int(std::string_view Token, const char* pWhat) {
    const auto v = feb_int(Token);
    if (!v)
        feb_fail(std::string("bad ") + pWhat + " '" + std::string(Token) + "'");
    return *v;
}

double feb_need_double(std::string_view Token) {
    const std::string text(Token);
    const char* end = nullptr;
    const double v = detail::parse_double(text.c_str(), end);
    if (end != text.c_str() + text.size())
        feb_fail("bad number '" + text + "'");
    return v;
}

// A list of ids with FEBio 4's `a:b[:step]` ranges.
std::vector<std::int64_t> feb_id_list(std::string_view Text) {
    std::vector<std::string_view> tokens;
    feb_split(Text, tokens);
    std::vector<std::int64_t> out;
    for (std::string_view t : tokens) {
        const std::size_t colon = t.find(':');
        if (colon == std::string_view::npos) {
            out.push_back(feb_need_int(t, "id"));
            continue;
        }
        const std::size_t colon2 = t.find(':', colon + 1);
        const std::int64_t a = feb_need_int(t.substr(0, colon), "range");
        const std::int64_t b = feb_need_int(
            t.substr(colon + 1, colon2 == std::string_view::npos ? std::string_view::npos
                                                                 : colon2 - colon - 1),
            "range");
        const std::int64_t step =
            colon2 == std::string_view::npos ? 1 : feb_need_int(t.substr(colon2 + 1), "range step");
        if (step <= 0)
            feb_fail("bad range '" + std::string(t) + "'");
        for (std::int64_t v = a; v <= b; v += step)
            out.push_back(v);
    }
    return out;
}

// A region under construction: entries in first-seen order, merged by name.
struct FebGroup {
    std::string mName;
    RegionKind mKind;
    int mDim = -1;
    std::int64_t mTag = -1;
    std::vector<std::int64_t> mEntries;
};

// A block of cells made from a Surface, Edge or DiscreteSet, added after the
// <Elements> blocks.
struct FebExtraBlock {
    std::string mType;
    std::vector<std::int64_t> mConn;
    std::string mRegion;
};

struct FebReader {
    std::string mPath;
    bool mLenient = false;
    std::vector<double> mCoords;
    std::unordered_map<std::int64_t, std::int64_t> mNodeIndex;
    // <Elements> blocks.
    std::vector<std::string> mBlockTypes;
    std::vector<std::vector<std::int64_t>> mBlockConn;
    std::unordered_map<std::int64_t, std::int64_t> mElemIndex;  // element id -> global cell
    std::int64_t mNumCells = 0;
    // Ordered member lists by set name, for MeshData's 1-based `lid`.
    std::map<std::string, std::vector<std::int64_t>> mNodeSets, mElemSets;
    std::vector<FebGroup> mGroups;
    std::map<std::pair<int, std::string>, std::size_t> mGroupIndex;
    std::vector<
        std::pair<std::string, std::vector<std::pair<std::string, std::vector<std::int64_t>>>>>
        mSurfaces;
    std::vector<FebExtraBlock> mExtra;
    std::map<std::string, std::int64_t> mMaterialIds;     // material name -> id
    std::map<std::string, std::int64_t> mDomainMaterial;  // domain name -> material id
    std::size_t mSkippedLossy = 0;

    FebGroup& Group(const std::string& rName, RegionKind Kind) {
        const auto key = std::make_pair(static_cast<int>(Kind), rName);
        const auto it = mGroupIndex.find(key);
        if (it != mGroupIndex.end())
            return mGroups[it->second];
        mGroupIndex.emplace(key, mGroups.size());
        mGroups.push_back(FebGroup{rName, Kind, -1, -1, {}});
        return mGroups.back();
    }

    std::int64_t Node(std::int64_t Id) const {
        const auto it = mNodeIndex.find(Id);
        if (it == mNodeIndex.end())
            feb_fail("undefined node " + std::to_string(Id));
        return it->second;
    }

    std::int64_t Element(std::int64_t Id) const {
        const auto it = mElemIndex.find(Id);
        if (it == mElemIndex.end())
            feb_fail("undefined element " + std::to_string(Id));
        return it->second;
    }

    void ReadNodes(const pugi::xml_node& rNodes) {
        std::vector<std::string_view> tokens;
        std::vector<std::int64_t> members;
        for (const pugi::xml_node& n : rNodes.children()) {
            if (n.type() != pugi::node_element)
                continue;
            const std::int64_t id = feb_need_int(n.attribute("id").value(), "node id");
            feb_split(n.child_value(), tokens);
            if (tokens.size() != 3)
                feb_fail("node " + std::to_string(id) + " needs x,y,z");
            if (!mNodeIndex.emplace(id, static_cast<std::int64_t>(mCoords.size() / 3)).second)
                feb_fail("node " + std::to_string(id) + " is defined twice");
            members.push_back(static_cast<std::int64_t>(mCoords.size() / 3));
            for (std::string_view t : tokens)
                mCoords.push_back(feb_need_double(t));
        }
        // A named <Nodes> block is also a node set, in FEBio.
        const std::string name = rNodes.attribute("name").value();
        if (!name.empty()) {
            FebGroup& g = Group(name, RegionKind::Point);
            g.mEntries.insert(g.mEntries.end(), members.begin(), members.end());
            auto& ordered = mNodeSets[name];
            ordered.insert(ordered.end(), members.begin(), members.end());
        }
    }

    void ReadElements(const pugi::xml_node& rElems, std::size_t Index, bool Spec25) {
        const std::string raw_type = rElems.attribute("type").value();
        const FebType* type = feb_type(raw_type);
        std::size_t nodes = type ? type->mNodes : 0;
        std::string meshio_type = type ? type->mType : "";
        if (!type) {
            const std::string base = feb_base_type(raw_type);
            const FebLossy* lossy = nullptr;
            for (const FebLossy& l : kFebLossy)
                if (base == l.mFebio)
                    lossy = &l;
            if (!lossy)
                feb_fail("unknown element type '" + raw_type + "'");
            if (!mLenient || !lossy->mDowngrade)
                feb_fail(
                    "element type '" + raw_type + "' has no meshio++ cell type" +
                    std::string(lossy->mDowngrade ? " (read with lenient to downgrade it)" : ""));
            nodes = lossy->mNodes;
            meshio_type = lossy->mDowngrade;
            ++mSkippedLossy;
        }
        const std::size_t keep =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(meshio_type)));
        const detail::NodeOrder* order = detail::node_order("febio", meshio_type);
        std::string name = rElems.attribute("name").value();
        if (name.empty())
            name = "Part" + std::to_string(Index + 1);
        std::int64_t tag = -1;
        const std::string mat = rElems.attribute("mat").value();
        if (!mat.empty()) {
            const auto it = mMaterialIds.find(mat);
            tag = it != mMaterialIds.end() ? it->second : feb_int(mat).value_or(-1);
        } else if (const auto it = mDomainMaterial.find(name); it != mDomainMaterial.end()) {
            tag = it->second;
        }
        (void)Spec25;

        std::vector<std::int64_t> conn;
        std::vector<std::string_view> tokens;
        std::vector<std::int64_t> row(nodes);
        std::vector<std::int64_t> members;
        for (const pugi::xml_node& e : rElems.children()) {
            if (e.type() != pugi::node_element)
                continue;
            const std::int64_t id = feb_need_int(e.attribute("id").value(), "element id");
            feb_split(e.child_value(), tokens);
            if (tokens.size() != nodes)
                feb_fail("element " + std::to_string(id) + " of type " + raw_type + " needs " +
                         std::to_string(nodes) + " nodes");
            for (std::size_t k = 0; k < nodes; ++k)
                row[k] = Node(feb_need_int(tokens[k], "node id"));
            for (std::size_t k = 0; k < keep; ++k)
                conn.push_back(row[order ? static_cast<std::size_t>(order->mToMeshio[k]) : k]);
            const std::int64_t global = mNumCells++;
            if (!mElemIndex.emplace(id, global).second)
                feb_fail("element " + std::to_string(id) + " is defined twice");
            members.push_back(global);
        }
        mBlockTypes.push_back(meshio_type);
        mBlockConn.push_back(std::move(conn));
        FebGroup& g = Group(name, RegionKind::Cell);
        g.mDim = std::max(g.mDim, cell_type_dimension(cell_type_from_name(meshio_type)));
        if (tag >= 0)
            g.mTag = tag;
        g.mEntries.insert(g.mEntries.end(), members.begin(), members.end());
        auto& ordered = mElemSets[name];
        ordered.insert(ordered.end(), members.begin(), members.end());
    }

    // NodeSet in any spec's form: a comma list (ranges in 4.0), children with
    // an `id`, `<node_list>`, or `<NodeSet node_set=>` including another set.
    void ReadNodeSet(const pugi::xml_node& rSet) {
        const std::string name = rSet.attribute("name").value();
        std::vector<std::int64_t> members;
        const std::string text = rSet.child_value();
        for (std::int64_t id : feb_id_list(text))
            members.push_back(Node(id));
        for (const pugi::xml_node& c : rSet.children()) {
            if (c.type() != pugi::node_element)
                continue;
            const std::string tag = c.name();
            if (tag == "node_list") {
                for (std::int64_t id : feb_id_list(c.child_value()))
                    members.push_back(Node(id));
            } else if (tag == "NodeSet") {
                const auto it = mNodeSets.find(c.attribute("node_set").value());
                if (it == mNodeSets.end())
                    feb_fail(std::string("NodeSet '") + name + "' includes undefined set '" +
                             c.attribute("node_set").value() + "'");
                members.insert(members.end(), it->second.begin(), it->second.end());
            } else {
                members.push_back(Node(feb_need_int(c.attribute("id").value(), "node id")));
            }
        }
        // meshio++'s own MeshData sets (`meshdata:<array>`) carry an array, not a group.
        if (name.rfind("meshdata:", 0) != 0) {
            FebGroup& g = Group(name, RegionKind::Point);
            g.mEntries.insert(g.mEntries.end(), members.begin(), members.end());
        }
        auto& ordered = mNodeSets[name];
        ordered.insert(ordered.end(), members.begin(), members.end());
    }

    void ReadElementSet(const pugi::xml_node& rSet) {
        const std::string name = rSet.attribute("name").value();
        std::vector<std::int64_t> members;
        for (std::int64_t id : feb_id_list(rSet.child_value()))
            members.push_back(Element(id));
        for (const pugi::xml_node& c : rSet.children())
            if (c.type() == pugi::node_element)
                members.push_back(Element(feb_need_int(c.attribute("id").value(), "element id")));
        if (name.rfind("meshdata:", 0) != 0) {
            FebGroup& g = Group(name, RegionKind::Cell);
            g.mEntries.insert(g.mEntries.end(), members.begin(), members.end());
        }
        auto& ordered = mElemSets[name];
        ordered.insert(ordered.end(), members.begin(), members.end());
    }

    // A facet of a Surface, or a line of an Edge: its meshio++ type and nodes.
    std::pair<std::string, std::vector<std::int64_t>> Facet(const pugi::xml_node& rFacet) {
        const FebType* type = feb_type(rFacet.name());
        if (!type || cell_type_dimension(cell_type_from_name(type->mType)) > 2)
            feb_fail(std::string("unknown facet type '") + rFacet.name() + "'");
        std::vector<std::string_view> tokens;
        feb_split(rFacet.child_value(), tokens);
        if (tokens.size() != type->mNodes)
            feb_fail(std::string("a ") + rFacet.name() + " facet needs " +
                     std::to_string(type->mNodes) + " nodes");
        std::vector<std::int64_t> nodes;
        for (std::string_view t : tokens)
            nodes.push_back(Node(feb_need_int(t, "node id")));
        return {type->mType, std::move(nodes)};
    }

    void ReadSurface(const pugi::xml_node& rSurface) {
        std::vector<std::pair<std::string, std::vector<std::int64_t>>> facets;
        for (const pugi::xml_node& f : rSurface.children())
            if (f.type() == pugi::node_element)
                facets.push_back(Facet(f));
        mSurfaces.emplace_back(rSurface.attribute("name").value(), std::move(facets));
    }

    // Edge lines and discrete-set springs become line blocks.
    void ReadLines(const pugi::xml_node& rSet, bool Discrete) {
        const std::string name = rSet.attribute("name").value();
        std::map<std::string, std::size_t> by_type;
        for (const pugi::xml_node& c : rSet.children()) {
            if (c.type() != pugi::node_element)
                continue;
            std::string type = "line";
            std::vector<std::int64_t> nodes;
            if (Discrete) {
                std::vector<std::string_view> tokens;
                feb_split(c.child_value(), tokens);
                if (tokens.size() != 2)
                    feb_fail("a <delem> needs two nodes");
                for (std::string_view t : tokens)
                    nodes.push_back(Node(feb_need_int(t, "node id")));
            } else {
                auto facet = Facet(c);
                type = facet.first;
                nodes = std::move(facet.second);
            }
            auto [it, fresh] = by_type.emplace(type, mExtra.size());
            if (fresh)
                mExtra.push_back(FebExtraBlock{type, {}, name});
            auto& conn = mExtra[it->second].mConn;
            conn.insert(conn.end(), nodes.begin(), nodes.end());
        }
    }
};

// Values per `<e lid=>`/`<node lid=>` child. With no data type (FEBio's
// `type="fiber"`/`var="fiber"` element data is a vec3 without saying so), the
// first entry's value count.
std::size_t feb_components(const std::string& rType, const pugi::xml_node& rData) {
    if (rType.empty()) {
        std::vector<std::string_view> tokens;
        for (const pugi::xml_node& c : rData.children())
            if (c.type() == pugi::node_element) {
                feb_split(c.child_value(), tokens);
                return std::max<std::size_t>(tokens.size(), 1);
            }
        return 1;
    }
    if (rType == "scalar")
        return 1;
    if (rType == "vec2")
        return 2;
    if (rType == "vec3")
        return 3;
    if (rType == "mat3s")
        return 6;
    if (rType == "mat3")
        return 9;
    feb_fail("unknown data type '" + rType + "'");
}

// Fills `rOut` (rows of `Width`) at the set's members from the `lid` children.
void feb_fill(const pugi::xml_node& rData, const std::vector<std::int64_t>& rMembers,
              std::size_t Width, std::vector<double>& rOut, std::int64_t Base) {
    std::vector<std::string_view> tokens;
    for (const pugi::xml_node& c : rData.children()) {
        if (c.type() != pugi::node_element)
            continue;
        const std::int64_t lid = feb_need_int(c.attribute("lid").value(), "lid");
        if (lid < 1 || static_cast<std::size_t>(lid) > rMembers.size())
            feb_fail(std::string("lid ") + std::to_string(lid) + " is outside its set");
        feb_split(c.child_value(), tokens);
        if (tokens.size() != Width)
            feb_fail(std::string("MeshData '") + rData.attribute("name").value() + "' needs " +
                     std::to_string(Width) + " values per entry");
        const auto row =
            static_cast<std::size_t>(rMembers[static_cast<std::size_t>(lid - 1)] - Base);
        for (std::size_t k = 0; k < Width; ++k)
            rOut[row * Width + k] = feb_need_double(tokens[k]);
    }
}

}  // namespace

Mesh read_febio(const std::string& rPath, const ReadOptions& rOptions) {
    pugi::xml_document doc;
    const pugi::xml_parse_result parsed = doc.load_file(rPath.c_str());
    if (!parsed)
        throw ReadError("FEBio .feb: could not parse " + rPath + ": " + parsed.description());
    const pugi::xml_node root = doc.child("febio_spec");
    if (!root)
        throw ReadError("FEBio .feb: " + rPath + " has no <febio_spec> root");
    const std::string version = root.attribute("version").value();
    if (version != "2.5" && version != "3.0" && version != "4.0")
        throw ReadError("FEBio .feb: febio_spec version '" + version +
                        "' is not supported (2.5, 3.0 and 4.0 are)");
    const bool spec25 = version == "2.5";

    FebReader reader;
    reader.mPath = rPath;
    reader.mLenient = rOptions.mLenient;

    for (const pugi::xml_node& m : root.child("Material").children("material")) {
        const std::string name = m.attribute("name").value();
        const auto id = feb_int(m.attribute("id").value());
        if (!name.empty() && id)
            reader.mMaterialIds[name] = *id;
    }
    for (const pugi::xml_node& d : root.child("MeshDomains").children()) {
        const auto it = reader.mMaterialIds.find(d.attribute("mat").value());
        if (it != reader.mMaterialIds.end())
            reader.mDomainMaterial[d.attribute("name").value()] = it->second;
    }

    // The mesh section, possibly pulled from another file with `from=`.
    pugi::xml_document included;
    pugi::xml_node mesh_node =
        spec25 ? root.child("Geometry")
               : (root.child("Mesh") ? root.child("Mesh") : root.child("Geometry"));
    if (!mesh_node)
        throw ReadError("FEBio .feb: " + rPath + " has no " +
                        std::string(spec25 ? "<Geometry>" : "<Mesh>") + " section");
    if (const char* from = mesh_node.attribute("from").value(); *from) {
        fs::path other(from);
        if (other.is_relative())
            other = fs::path(rPath).parent_path() / other;
        if (!included.load_file(other.string().c_str()))
            throw ReadError("FEBio .feb: could not read '" + other.string() + "' named by from=");
        mesh_node = included.child("febio_spec").child(mesh_node.name());
        if (!mesh_node)
            throw ReadError("FEBio .feb: '" + other.string() + "' has no mesh section");
    }

    std::size_t n_blocks = 0;
    for (const pugi::xml_node& c : mesh_node.children()) {
        if (c.type() != pugi::node_element)
            continue;
        const std::string tag = c.name();
        if (tag == "Nodes")
            reader.ReadNodes(c);
        else if (tag == "Elements")
            reader.ReadElements(c, n_blocks++, spec25);
        else if (tag == "NodeSet")
            reader.ReadNodeSet(c);
        else if (tag == "ElementSet")
            reader.ReadElementSet(c);
        else if (tag == "Surface")
            reader.ReadSurface(c);
        else if (tag == "Edge")
            reader.ReadLines(c, false);
        else if (tag == "DiscreteSet")
            reader.ReadLines(c, true);
        else if (tag == "Part" || tag == "Instance")
            throw ReadError(
                "FEBio .feb: the <Part>/<Instance> form is not supported; export the "
                "model from FEBio Studio as a plain mesh");
        else if (tag != "SurfacePair" && tag != "PartList" && tag != "NodeSetPair" &&
                 tag != "NodeSetSet")
            log::warn("FEBio .feb: <{}> in the mesh section is not read", tag);
    }
    if (reader.mSkippedLossy)
        log::warn("FEBio .feb: {} element block(s) downgraded to a meshio++ cell type",
                  reader.mSkippedLossy);

    Mesh mesh;
    NDArray points(DType::Float64, {reader.mCoords.size() / 3, 3});
    std::copy(reader.mCoords.begin(), reader.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    for (std::size_t b = 0; b < reader.mBlockTypes.size(); ++b) {
        const std::size_t k = static_cast<std::size_t>(
            cell_type_num_nodes(cell_type_from_name(reader.mBlockTypes[b])));
        const auto& conn = reader.mBlockConn[b];
        NDArray arr(DType::Int64, {conn.size() / k, k});
        std::copy(conn.begin(), conn.end(), arr.As<std::int64_t>());
        mesh.AddCellBlock(reader.mBlockTypes[b], std::move(arr));
    }
    const std::int64_t n_elements = reader.mNumCells;

    // Surfaces: Side regions when every facet is a face of a solid element,
    // otherwise their own cell blocks, as FEBio also allows.
    std::optional<detail::FacetIndex> faces;
    std::size_t detached = 0;
    std::int64_t extra_base = n_elements;
    std::vector<FebExtraBlock> surface_blocks;
    for (auto& [name, facets] : reader.mSurfaces) {
        if (!faces) {
            detail::FacetIndexOptions options;
            options.mSurfaceEdges = false;
            faces.emplace(mesh, options);
        }
        std::vector<std::int64_t> sides;
        bool all = true;
        for (const auto& [type, nodes] : facets) {
            const std::size_t corners = feb_num_corners(type);
            const detail::FacetHit* hit = faces->Find(nodes.data(), corners);
            if (!hit) {
                all = false;
                break;
            }
            sides.push_back(hit->mFirst.mCell);
            sides.push_back(hit->mFirst.mFacet);
        }
        if (all) {
            FebGroup& g = reader.Group(name, RegionKind::Side);
            g.mDim = 2;
            g.mEntries.insert(g.mEntries.end(), sides.begin(), sides.end());
            continue;
        }
        ++detached;
        std::map<std::string, std::size_t> by_type;
        for (const auto& [type, nodes] : facets) {
            auto [it, fresh] = by_type.emplace(type, surface_blocks.size());
            if (fresh)
                surface_blocks.push_back(FebExtraBlock{type, {}, name});
            auto& conn = surface_blocks[it->second].mConn;
            conn.insert(conn.end(), nodes.begin(), nodes.end());
        }
    }
    if (detached)
        log::warn("FEBio .feb: {} surface(s) do not lie on solid elements; read as cell blocks",
                  detached);
    for (const std::vector<FebExtraBlock>* list : {&surface_blocks, &reader.mExtra}) {
        for (const FebExtraBlock& e : *list) {
            const std::size_t k =
                static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(e.mType)));
            const std::size_t rows = e.mConn.size() / k;
            NDArray arr(DType::Int64, {rows, k});
            std::copy(e.mConn.begin(), e.mConn.end(), arr.As<std::int64_t>());
            mesh.AddCellBlock(e.mType, std::move(arr));
            FebGroup& g = reader.Group(e.mRegion, RegionKind::Cell);
            g.mDim = std::max(g.mDim, cell_type_dimension(cell_type_from_name(e.mType)));
            for (std::size_t r = 0; r < rows; ++r)
                g.mEntries.push_back(extra_base + static_cast<std::int64_t>(r));
            extra_base += static_cast<std::int64_t>(rows);
        }
    }

    // MeshData: NaN outside the set each array is defined on.
    const std::vector<std::int64_t> bases = detail::block_bases(mesh);
    for (const pugi::xml_node& d : root.child("MeshData").children()) {
        if (d.type() != pugi::node_element)
            continue;
        const std::string tag = d.name();
        std::string type = d.attribute("data_type").value();
        if (type.empty())
            type = d.attribute("datatype").value();
        std::string name = d.attribute("name").value();
        if (name.empty())
            name = d.attribute("var").value();
        if (tag == "NodeData") {
            const auto set = reader.mNodeSets.find(d.attribute("node_set").value());
            if (set == reader.mNodeSets.end())
                feb_fail("NodeData '" + name + "' names an undefined node set");
            const std::size_t width = feb_components(type, d);
            std::vector<double> values(mesh.NumPoints() * width,
                                       std::numeric_limits<double>::quiet_NaN());
            feb_fill(d, set->second, width, values, 0);
            NDArray arr = width == 1 ? NDArray(DType::Float64, {mesh.NumPoints()})
                                     : NDArray(DType::Float64, {mesh.NumPoints(), width});
            std::copy(values.begin(), values.end(), arr.As<double>());
            mesh.AddPointData(name, std::move(arr));
        } else if (tag == "ElementData" && d.first_child()) {
            const auto set = reader.mElemSets.find(d.attribute("elem_set").value());
            if (set == reader.mElemSets.end())
                feb_fail("ElementData '" + name + "' names an undefined element set");
            const std::size_t width = feb_components(type, d);
            std::vector<double> values(static_cast<std::size_t>(bases.back()) * width,
                                       std::numeric_limits<double>::quiet_NaN());
            feb_fill(d, set->second, width, values, 0);
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
                const auto rows = static_cast<std::size_t>(bases[b + 1] - bases[b]);
                NDArray arr = width == 1 ? NDArray(DType::Float64, {rows})
                                         : NDArray(DType::Float64, {rows, width});
                std::copy(values.begin() + bases[b] * static_cast<std::int64_t>(width),
                          values.begin() + bases[b + 1] * static_cast<std::int64_t>(width),
                          arr.As<double>());
                blocks.push_back(std::move(arr));
            }
            mesh.AddCellData(name, std::move(blocks));
        } else {
            log::warn("FEBio .feb: <{} name=\"{}\"> in MeshData is not read", tag, name);
        }
    }

    for (FebGroup& g : reader.mGroups) {
        const std::size_t stride = g.mKind == RegionKind::Side ? 2 : 1;
        NDArray entries = stride == 2 ? NDArray(DType::Int64, {g.mEntries.size() / 2, 2})
                                      : NDArray(DType::Int64, {g.mEntries.size()});
        std::copy(g.mEntries.begin(), g.mEntries.end(), entries.As<std::int64_t>());
        mesh.AddRegion(Region(g.mName, g.mKind, g.mDim, g.mTag, std::move(entries)));
    }
    return mesh;
}

namespace {

std::string feb_escape(const std::string& rText) {
    std::string out;
    for (char c : rText) {
        switch (c) {
            case '&':
                out += "&amp;";
                break;
            case '<':
                out += "&lt;";
                break;
            case '>':
                out += "&gt;";
                break;
            case '"':
                out += "&quot;";
                break;
            default:
                out += c;
        }
    }
    return out;
}

void feb_append_int(std::string& rOut, std::int64_t Value) {
    char buf[24];
    const auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), Value);
    rOut.append(buf, ptr);
}

void feb_append_ids(std::string& rOut, const std::vector<std::int64_t>& rIds) {
    for (std::size_t k = 0; k < rIds.size(); ++k) {
        if (k)
            rOut += ',';
        feb_append_int(rOut, rIds[k]);
    }
}

enum class FebBlockKind { Elements, Surface, Edge, Discrete, Dropped };

// Whether every line of block `Line` joins two corners that are adjacent in a
// face of a solid or along the edge of a surface cell.
bool feb_all_on_edges(const Mesh& rMesh, std::size_t Line, const std::vector<int>& rDims) {
    std::set<std::pair<std::int64_t, std::int64_t>> edges;
    const auto add = [&](std::int64_t A, std::int64_t B) {
        edges.emplace(std::min(A, B), std::max(A, B));
    };
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        if (rDims[b] < 2)
            continue;
        const auto cb = rMesh.Cells(b);
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            if (rDims[b] == 3) {
                for (const detail::CellFaceDef& f : detail::cell_faces(type))
                    for (std::size_t c = 0; c < f.mNumCorners; ++c)
                        add(detail::read_int(conn, r * k + f.mNodes[c]),
                            detail::read_int(conn, r * k + f.mNodes[(c + 1) % f.mNumCorners]));
            } else {
                for (const detail::CellEdgeDef& e : detail::cell_edges(type))
                    add(detail::read_int(conn, r * k + e.mNodes[0]),
                        detail::read_int(conn, r * k + e.mNodes[1]));
            }
        }
    }
    const auto cb = rMesh.Cells(Line);
    const NDArray& conn = cb.Conn();
    const std::size_t k = cb.NodesPerCell();
    for (std::size_t r = 0; r < cb.NumCells(); ++r) {
        const std::int64_t a = detail::read_int(conn, r * k);
        const std::int64_t b = detail::read_int(conn, r * k + 1);
        if (!edges.count({std::min(a, b), std::max(a, b)}))
            return false;
    }
    return true;
}

}  // namespace

void write_febio(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t n_blocks = rMesh.NumCellBlocks();
    int max_dim = 0;
    std::vector<int> dims(n_blocks);
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        if (cb.IsRagged() || (type != "vertex" && !feb_write_type(type)))
            throw WriteError("FEBio .feb writer: cell type '" + type +
                             "' has no FEBio element type");
        dims[b] = cell_type_dimension(cell_type_from_name(type));
        if (cb.NumCells())
            max_dim = std::max(max_dim, dims[b]);
    }
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);

    // What each block becomes: faces of solids a Surface, lines next to higher
    // cells an Edge, vertices nothing, everything else elements.
    std::optional<detail::FacetIndex> faces;
    std::vector<FebBlockKind> kinds(n_blocks, FebBlockKind::Elements);
    std::size_t dropped_vertices = 0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        if (dims[b] == 0) {
            kinds[b] = FebBlockKind::Dropped;
            dropped_vertices += cb.NumCells();
        } else if (dims[b] == 1 && max_dim > 1) {
            // Lines along the edges of other cells are an <Edge>; two-node lines
            // between nodes no cell connects are springs, a <DiscreteSet>.
            kinds[b] = std::string(cb.Type()) == "line" && !feb_all_on_edges(rMesh, b, dims)
                           ? FebBlockKind::Discrete
                           : FebBlockKind::Edge;
        } else if (dims[b] == 2 && max_dim == 3 && cb.NumCells() > 0) {
            if (!faces) {
                detail::FacetIndexOptions options;
                options.mSurfaceEdges = false;
                faces.emplace(rMesh, options);
            }
            const NDArray& conn = cb.Conn();
            const std::size_t k = cb.NodesPerCell();
            const std::size_t corners = feb_num_corners(cb.Type());
            std::vector<std::int64_t> row(corners);
            bool all = true;
            for (std::size_t r = 0; r < cb.NumCells() && all; ++r) {
                for (std::size_t c = 0; c < corners; ++c)
                    row[c] = detail::read_int(conn, r * k + c);
                all = faces->Find(row.data(), corners) != nullptr;
            }
            if (all)
                kinds[b] = FebBlockKind::Surface;
        }
    }

    // Block names: the name of a cell region covering exactly that block.
    std::vector<std::string> names(n_blocks);
    std::vector<bool> region_names_block(rMesh.NumRegions(), false);
    std::set<std::string> taken;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        for (std::size_t r = 0; r < rMesh.NumRegions() && names[b].empty(); ++r) {
            const Region& reg = rMesh.Region(r);
            if (reg.mKind != RegionKind::Cell || region_names_block[r] || reg.NumEntries() == 0 ||
                static_cast<std::int64_t>(reg.NumEntries()) != bases[b + 1] - bases[b])
                continue;
            const std::int64_t* e = reg.Entries();
            if (e[0] == bases[b] && e[reg.NumEntries() - 1] == bases[b + 1] - 1 &&
                !taken.count(reg.mName)) {
                names[b] = reg.mName;
                region_names_block[r] = true;
            }
        }
        if (names[b].empty()) {
            const char* stem = kinds[b] == FebBlockKind::Surface    ? "Surface"
                               : kinds[b] == FebBlockKind::Edge     ? "Edge"
                               : kinds[b] == FebBlockKind::Discrete ? "DiscreteSet"
                                                                    : "Part";
            names[b] = stem + std::to_string(b + 1);
            while (taken.count(names[b]))
                names[b] += "_";
        }
        taken.insert(names[b]);
    }

    // Element numbers, 1-based, over the <Elements> blocks only.
    std::vector<std::int64_t> element_no(static_cast<std::size_t>(bases.back()), 0);
    std::int64_t n_elements = 0;
    for (std::size_t b = 0; b < n_blocks; ++b)
        if (kinds[b] == FebBlockKind::Elements)
            for (std::int64_t g = bases[b]; g < bases[b + 1]; ++g)
                element_no[static_cast<std::size_t>(g)] = ++n_elements;

    // --- notes --------------------------------------------------------------------
    if (dropped_vertices) {
        log::warn("FEBio .feb writer: {} vertex cell(s) have no FEBio element and were dropped",
                  dropped_vertices);
        detail::provenance_note("cells-dropped", std::to_string(dropped_vertices) +
                                                     " vertex cell(s) have no FEBio element");
    }
    // MeshData: every point array over the nodes where it is defined, every cell
    // array over the <Elements> cells where it is; each on its own set,
    // `meshdata:<name>`, which the reader recognises and makes no region of.
    struct FebArray {
        std::string mName;
        bool mNodal = true;
        std::string mType;
        std::size_t mWidth = 1;
        std::vector<std::int64_t> mMembers;  // 1-based node or element numbers
        std::vector<double> mValues;         // mMembers.size() * mWidth
    };
    std::vector<FebArray> arrays;
    std::vector<std::string> unwritable;
    auto data_type = [](std::size_t Width) -> const char* {
        switch (Width) {
            case 1:
                return "scalar";
            case 2:
                return "vec2";
            case 3:
                return "vec3";
            case 6:
                return "mat3s";
            case 9:
                return "mat3";
            default:
                return nullptr;
        }
    };
    auto width_of = [](const NDArray& rA) {
        std::size_t w = 1;
        for (std::size_t k = 1; k < rA.Shape().size(); ++k)
            w *= rA.Shape()[k];
        return w;
    };
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const std::size_t w = width_of(a);
        const char* type = data_type(w);
        if (!type || a.Shape().empty()) {
            unwritable.push_back(name);
            continue;
        }
        FebArray arr{name, true, type, w, {}, {}};
        for (std::size_t p = 0; p < rMesh.NumPoints(); ++p) {
            bool defined = true;
            for (std::size_t c = 0; c < w && defined; ++c)
                defined = !std::isnan(detail::read_double(a, p * w + c));
            if (!defined)
                continue;
            arr.mMembers.push_back(static_cast<std::int64_t>(p + 1));
            for (std::size_t c = 0; c < w; ++c)
                arr.mValues.push_back(detail::read_double(a, p * w + c));
        }
        if (!arr.mMembers.empty())
            arrays.push_back(std::move(arr));
    }
    std::size_t off_elements = 0;
    for (const std::string& name : rMesh.CellDataNames()) {
        std::size_t w = 0;
        bool ok = true;
        for (std::size_t b = 0; b < n_blocks && ok; ++b) {
            if (rMesh.Cells(b).NumCells() == 0)
                continue;
            const NDArray& a = rMesh.CellData(name, b);
            if (a.Shape().empty()) {
                ok = false;
                continue;
            }
            const std::size_t wb = width_of(a);
            if (w == 0)
                w = wb;
            else if (w != wb)
                ok = false;
        }
        const char* type = ok ? data_type(w == 0 ? 1 : w) : nullptr;
        if (!type) {
            unwritable.push_back(name);
            continue;
        }
        w = w == 0 ? 1 : w;
        FebArray arr{name, false, type, w, {}, {}};
        for (std::size_t b = 0; b < n_blocks; ++b) {
            const NDArray& a = rMesh.CellData(name, b);
            for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r) {
                bool defined = true;
                for (std::size_t c = 0; c < w && defined; ++c)
                    defined = !std::isnan(detail::read_double(a, r * w + c));
                if (!defined)
                    continue;
                if (kinds[b] != FebBlockKind::Elements) {
                    ++off_elements;
                    continue;
                }
                arr.mMembers.push_back(
                    element_no[static_cast<std::size_t>(bases[b] + static_cast<std::int64_t>(r))]);
                for (std::size_t c = 0; c < w; ++c)
                    arr.mValues.push_back(detail::read_double(a, r * w + c));
            }
        }
        if (!arr.mMembers.empty())
            arrays.push_back(std::move(arr));
    }
    if (!unwritable.empty() || off_elements != 0 || rMesh.NumFieldData() != 0) {
        std::string what;
        for (const std::string& n : unwritable)
            what += (what.empty() ? "" : ", ") + n;
        if (!unwritable.empty())
            log::warn("FEBio .feb writer: arrays FEBio has no data type for are dropped: {}", what);
        if (off_elements != 0)
            log::warn(
                "FEBio .feb writer: {} cell value(s) on surfaces, edges or discrete sets "
                "are dropped",
                off_elements);
        detail::provenance_note("data-dropped",
                                "field data, and arrays that are not scalar, vec2, vec3, mat3s "
                                "or mat3 or lie off the elements, are not written");
    }

    // The provenance comment goes inside the root, as in VTU: the root tag stays
    // within the bytes sniff_format reads.
    std::string out =
        "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<febio_spec version=\"4.0\">\n";
    out += detail::provenance_render_xml_comment(detail::SlotTier::Block);
    out += "\n\t<Module type=\"solid\"/>\n";
    out +=
        "\t<Material>\n\t\t<!-- Placeholder materials: FEBio needs one per domain. Replace them. "
        "-->\n";
    std::int64_t mat_id = 0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        if (kinds[b] != FebBlockKind::Elements)
            continue;
        out += "\t\t<material id=\"";
        feb_append_int(out, ++mat_id);
        out +=
            "\" name=\"" + feb_escape(names[b]) +
            "\" type=\"isotropic elastic\">\n\t\t\t<E>1</E>\n\t\t\t<v>0.3</v>\n\t\t</material>\n";
    }
    out += "\t</Material>\n\t<Mesh>\n\t\t<Nodes>\n";
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    char buf[40];
    for (std::size_t p = 0; p < rMesh.NumPoints(); ++p) {
        out += "\t\t\t<node id=\"";
        feb_append_int(out, static_cast<std::int64_t>(p + 1));
        out += "\">";
        for (std::size_t d = 0; d < 3; ++d) {
            const double v = d < pdim ? detail::read_double(points, p * pdim + d) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), d ? ",%.17g" : "%.17g", v);
            out += buf;
        }
        out += "</node>\n";
    }
    out += "\t\t</Nodes>\n";

    std::vector<std::int64_t> ids;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        if (kinds[b] == FebBlockKind::Dropped)
            continue;
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const char* ftype = feb_write_type(type);
        const detail::NodeOrder* order = detail::node_order("febio", type);
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const bool elements = kinds[b] == FebBlockKind::Elements;
        if (kinds[b] == FebBlockKind::Discrete) {
            out += "\t\t<DiscreteSet name=\"" + feb_escape(names[b]) + "\">\n";
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                out += "\t\t\t<delem>";
                feb_append_int(out, detail::read_int(conn, r * k) + 1);
                out += ',';
                feb_append_int(out, detail::read_int(conn, r * k + 1) + 1);
                out += "</delem>\n";
            }
            out += "\t\t</DiscreteSet>\n";
            continue;
        }
        const char* outer = elements                            ? "Elements"
                            : kinds[b] == FebBlockKind::Surface ? "Surface"
                                                                : "Edge";
        out += std::string("\t\t<") + outer;
        if (elements)
            out += std::string(" type=\"") + ftype + "\"";
        out += " name=\"" + feb_escape(names[b]) + "\">\n";
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::int64_t g = bases[b] + static_cast<std::int64_t>(r);
            out += std::string("\t\t\t<") + (elements ? "elem" : ftype) + " id=\"";
            feb_append_int(out, elements ? element_no[static_cast<std::size_t>(g)]
                                         : static_cast<std::int64_t>(r + 1));
            out += "\">";
            ids.clear();
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                ids.push_back(detail::read_int(conn, r * k + src) + 1);
            }
            feb_append_ids(out, ids);
            out += std::string("</") + (elements ? "elem" : ftype) + ">\n";
        }
        out += std::string("\t\t</") + outer + ">\n";
    }

    // Regions: point -> NodeSet, side -> Surface (or Edge on a surface mesh),
    // cell -> ElementSet over <Elements> cells.
    std::size_t skipped = 0;
    std::vector<std::int64_t> facet;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const std::int64_t* e = reg.Entries();
        const std::size_t n = reg.NumEntries();
        if (region_names_block[r])
            continue;
        if (reg.mKind == RegionKind::Point) {
            if (n == 0) {
                ++skipped;
                continue;
            }
            ids.assign(e, e + n);
            for (std::int64_t& v : ids)
                ++v;
            out += "\t\t<NodeSet name=\"" + feb_escape(reg.mName) + "\">";
            feb_append_ids(out, ids);
            out += "</NodeSet>\n";
        } else if (reg.mKind == RegionKind::Side) {
            std::string body;
            bool edge = false;
            std::size_t written = 0;
            for (std::size_t j = 0; j < n; ++j) {
                CellType ftype = CellType::Custom;
                if (!detail::facet_nodes(rMesh, e[2 * j], e[2 * j + 1], ftype, facet))
                    continue;
                const char* tag = feb_write_type(cell_type_name(ftype));
                edge = cell_type_dimension(ftype) == 1;
                ids.assign(facet.begin(), facet.end());
                for (std::int64_t& v : ids)
                    ++v;
                body += std::string("\t\t\t<") + tag + " id=\"";
                feb_append_int(body, static_cast<std::int64_t>(++written));
                body += "\">";
                feb_append_ids(body, ids);
                body += std::string("</") + tag + ">\n";
            }
            if (written == 0) {
                ++skipped;
                continue;
            }
            const char* outer = edge ? "Edge" : "Surface";
            out += std::string("\t\t<") + outer + " name=\"" + feb_escape(reg.mName) + "\">\n" +
                   body + "\t\t</" + outer + ">\n";
        } else {
            ids.clear();
            for (std::size_t j = 0; j < n; ++j)
                if (e[j] >= 0 && e[j] < bases.back() && element_no[static_cast<std::size_t>(e[j])])
                    ids.push_back(element_no[static_cast<std::size_t>(e[j])]);
            if (ids.empty()) {
                ++skipped;
                continue;
            }
            out += "\t\t<ElementSet name=\"" + feb_escape(reg.mName) + "\">";
            feb_append_ids(out, ids);
            out += "</ElementSet>\n";
        }
    }
    if (skipped)
        log::warn("FEBio .feb writer: {} region(s) with nothing FEBio can hold were dropped",
                  skipped);
    // The sets the MeshData arrays live on.
    for (const FebArray& a : arrays) {
        out += std::string("\t\t<") + (a.mNodal ? "NodeSet" : "ElementSet") + " name=\"" +
               feb_escape("meshdata:" + a.mName) + "\">";
        feb_append_ids(out, a.mMembers);
        out += std::string("</") + (a.mNodal ? "NodeSet" : "ElementSet") + ">\n";
    }

    out += "\t</Mesh>\n\t<MeshDomains>\n";
    for (std::size_t b = 0; b < n_blocks; ++b) {
        if (kinds[b] != FebBlockKind::Elements)
            continue;
        const char* domain = dims[b] == 3   ? "SolidDomain"
                             : dims[b] == 2 ? "ShellDomain"
                                            : "BeamDomain";
        out += std::string("\t\t<") + domain + " name=\"" + feb_escape(names[b]) + "\" mat=\"" +
               feb_escape(names[b]) + "\"/>\n";
    }
    out += "\t</MeshDomains>\n";
    if (!arrays.empty()) {
        out += "\t<MeshData>\n";
        for (const FebArray& a : arrays) {
            const char* tag = a.mNodal ? "NodeData" : "ElementData";
            out += std::string("\t\t<") + tag + " name=\"" + feb_escape(a.mName) + "\" " +
                   (a.mNodal ? "node_set" : "elem_set") + "=\"" +
                   feb_escape("meshdata:" + a.mName) + "\" data_type=\"" + a.mType + "\">\n";
            for (std::size_t m = 0; m < a.mMembers.size(); ++m) {
                out += std::string("\t\t\t<") + (a.mNodal ? "node" : "e") + " lid=\"";
                feb_append_int(out, static_cast<std::int64_t>(m + 1));
                out += "\">";
                for (std::size_t c = 0; c < a.mWidth; ++c) {
                    detail::snprintf_c(buf, sizeof(buf), c ? ",%.17g" : "%.17g",
                                       a.mValues[m * a.mWidth + c]);
                    out += buf;
                }
                out += std::string("</") + (a.mNodal ? "node" : "e") + ">\n";
            }
            out += std::string("\t\t</") + tag + ">\n";
        }
        out += "\t</MeshData>\n";
    }
    out += "</febio_spec>\n";

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    f << out;
}

}  // namespace meshioplusplus
