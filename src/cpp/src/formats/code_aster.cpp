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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/code_aster.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// Code_Aster reads the first 80 columns of every line (lirlig.F90) and nothing else.
constexpr std::size_t kCaColumns = 80;
// "N" or "M" plus seven digits: node and element names are at most 8 characters.
constexpr std::size_t kCaMaxEntities = 9999999;
constexpr std::size_t kCaMaxGroupName = 24;

struct CaTypeSpec {
    const char* mKeyword;
    const char* mType;
    std::size_t mNodes;
};

const std::vector<CaTypeSpec>& ca_types() {
    static const std::vector<CaTypeSpec> types = {
        {"POI1", "vertex", 1},          {"SEG2", "line", 2},
        {"SEG3", "line3", 3},           {"SEG4", "line4", 4},
        {"TRIA3", "triangle", 3},       {"TRIA6", "triangle6", 6},
        {"TRIA7", "triangle7", 7},      {"QUAD4", "quad", 4},
        {"QUAD8", "quad8", 8},          {"QUAD9", "quad9", 9},
        {"TETRA4", "tetra", 4},         {"TETRA10", "tetra10", 10},
        {"PENTA6", "wedge", 6},         {"PENTA15", "wedge15", 15},
        {"PENTA18", "wedge18", 18},     {"PYRAM5", "pyramid", 5},
        {"PYRAM13", "pyramid13", 13},   {"HEXA8", "hexahedron", 8},
        {"HEXA20", "hexahedron20", 20}, {"HEXA27", "hexahedron27", 27},
    };
    return types;
}

const CaTypeSpec* ca_spec_by_keyword(const std::string& rKeyword) {
    for (const CaTypeSpec& s : ca_types())
        if (rKeyword == s.mKeyword)
            return &s;
    return nullptr;
}

const CaTypeSpec* ca_spec_by_type(const std::string& rType) {
    for (const CaTypeSpec& s : ca_types())
        if (rType == s.mType)
            return &s;
    return nullptr;
}

std::string ca_upper(std::string s) {
    // ASCII only, like the name sanitiser: no locale may touch bytes above 0x7F.
    std::transform(s.begin(), s.end(), s.begin(), [](char c) {
        return c >= 'a' && c <= 'z' ? static_cast<char>(c - 'a' + 'A') : c;
    });
    return s;
}

struct CaToken {
    std::string mText;
    std::size_t mLine;
};

[[noreturn]] void ca_fail(const std::string& rWhat, std::size_t Line) {
    throw ReadError("Code_Aster .mail: " + rWhat + " (line " + std::to_string(Line) + ")");
}

// Splits the file into tokens: each line cut at column 80, a `%` comment
// dropped, then split on blanks and commas. `KEY = VALUE`, `KEY= VALUE` and
// `KEY =VALUE` are joined into one `KEY=VALUE` token.
std::vector<CaToken> ca_tokenize(const std::string& rText) {
    std::vector<CaToken> tokens;
    bool warned_long = false;
    std::size_t line_no = 0;
    std::size_t pos = 0;
    while (pos < rText.size()) {
        std::size_t eol = rText.find('\n', pos);
        if (eol == std::string::npos)
            eol = rText.size();
        std::string line = rText.substr(pos, eol - pos);
        pos = eol + 1;
        ++line_no;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.size() > kCaColumns) {
            const std::size_t comment = line.find('%');
            const bool beyond = line.find_first_not_of(" \t", kCaColumns) != std::string::npos &&
                                (comment == std::string::npos || comment >= kCaColumns);
            if (beyond && !warned_long) {
                log::warn(
                    "Code_Aster .mail: line {} is longer than 80 columns; like "
                    "Code_Aster, the reader ignores everything past column 80",
                    line_no);
                warned_long = true;
            }
            line.resize(kCaColumns);
        }
        const std::size_t comment = line.find('%');
        if (comment != std::string::npos)
            line.resize(comment);
        std::vector<std::string> words;
        std::string word;
        for (char c : line) {
            if (c == ' ' || c == '\t' || c == ',') {
                if (!word.empty())
                    words.push_back(std::move(word));
                word.clear();
            } else {
                word.push_back(c);
            }
        }
        if (!word.empty())
            words.push_back(std::move(word));
        for (std::size_t i = 0; i < words.size(); ++i) {
            std::string w = words[i];
            while (i + 1 < words.size() && (w.back() == '=' || words[i + 1].front() == '=')) {
                w += words[++i];
            }
            tokens.push_back({std::move(w), line_no});
        }
    }
    return tokens;
}

// Keywords Code_Aster's reader (lrmast.F90) skips without a word.
const std::set<std::string>& ca_skipped_keywords() {
    static const std::set<std::string> skipped = {"TITRE",    "DUMP",     "DEBUG",    "GROUP_FA",
                                                  "SYS_UNIT", "SYS_COOR", "MACRO_AR", "MACRO_FA",
                                                  "MACRO_EL", "MATERIAU"};
    return skipped;
}

// Something shaped like a block keyword: a letter, then letters, digits and `_`.
bool ca_is_keyword(const std::string& rUpper) {
    if (rUpper.empty() || rUpper[0] < 'A' || rUpper[0] > 'Z')
        return false;
    return std::all_of(rUpper.begin(), rUpper.end(), [](char c) {
        return (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
    });
}

bool ca_is_option(const std::string& rToken) {
    return rToken.find('=') != std::string::npos;
}

bool ca_is_block_end(const std::string& rUpper) {
    return rUpper == "FINSF";
}

bool ca_is_file_end(const std::string& rUpper) {
    return rUpper == "FIN" || rUpper.rfind("FIN(", 0) == 0;
}

// A Fortran real or integer: [sign] digits [. digits] [E|D [sign] digits], with
// at least one digit before the exponent. Nothing else (no `inf`, no hex).
bool ca_is_number(const std::string& rText) {
    std::size_t i = 0;
    const std::size_t n = rText.size();
    auto digits = [&]() {
        const std::size_t start = i;
        while (i < n && rText[i] >= '0' && rText[i] <= '9')
            ++i;
        return i - start;
    };
    if (i < n && (rText[i] == '+' || rText[i] == '-'))
        ++i;
    std::size_t mantissa = digits();
    if (i < n && rText[i] == '.') {
        ++i;
        mantissa += digits();
    }
    if (mantissa == 0)
        return false;
    if (i < n && (rText[i] == 'E' || rText[i] == 'e' || rText[i] == 'D' || rText[i] == 'd')) {
        ++i;
        if (i < n && (rText[i] == '+' || rText[i] == '-'))
            ++i;
        if (digits() == 0)
            return false;
    }
    return i == n;
}

double ca_number(const CaToken& rToken) {
    if (!ca_is_number(rToken.mText))
        ca_fail("expected a coordinate, found '" + rToken.mText + "'", rToken.mLine);
    std::string text = rToken.mText;
    for (char& c : text)
        if (c == 'D' || c == 'd')
            c = 'E';
    const char* end = nullptr;
    return detail::parse_double(text.c_str(), end);
}

struct CaElementBlock {
    const CaTypeSpec* mSpec;
    std::vector<std::string> mNames;
    std::vector<std::string> mNodes;  // mNames.size() * mSpec->mNodes, file order
    std::vector<std::size_t> mLines;
};

struct CaGroup {
    std::string mName;
    bool mCells;
    std::vector<std::pair<std::string, std::size_t>> mMembers;  // name, line
};

}  // namespace

Mesh read_code_aster(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Code_Aster .mail: cannot open " + rPath);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::vector<CaToken> tokens = ca_tokenize(text);

    int point_dim = 0;
    std::vector<double> coords;
    std::vector<std::string> node_names;
    std::vector<CaElementBlock> blocks;
    std::vector<CaGroup> groups;
    std::map<std::string, bool> warned_keywords;
    bool saw_fin = false;

    std::size_t i = 0;
    const std::size_t n = tokens.size();
    // Skips the options of a block header: every `KEY=VALUE` token.
    auto skip_options = [&](std::string* pName) {
        while (i < n && ca_is_option(tokens[i].mText)) {
            const std::string& t = tokens[i].mText;
            const std::size_t eq = t.find('=');
            const std::string key = ca_upper(t.substr(0, eq));
            if (pName && (key == "NOM" || key == "NAME"))
                *pName = t.substr(eq + 1);
            ++i;
        }
    };
    auto require_token = [&](const char* pWhat, std::size_t Line) {
        if (i >= n)
            ca_fail(std::string("the file ends inside a block, while reading ") + pWhat, Line);
    };

    while (i < n) {
        const CaToken& head = tokens[i];
        const std::string kw = ca_upper(head.mText);
        ++i;
        if (ca_is_file_end(kw)) {
            saw_fin = true;
            break;
        }
        if (ca_is_block_end(kw))
            continue;
        if (kw == "COOR_1D" || kw == "COOR_2D" || kw == "COOR_3D") {
            const int dim = kw[5] - '0';
            if (point_dim != 0 && point_dim != dim)
                ca_fail(kw + " after a COOR_" + std::to_string(point_dim) + "D block", head.mLine);
            point_dim = dim;
            skip_options(nullptr);
            while (true) {
                require_token("coordinates", head.mLine);
                if (ca_is_block_end(ca_upper(tokens[i].mText))) {
                    ++i;
                    break;
                }
                if (ca_is_option(tokens[i].mText)) {
                    ++i;
                    continue;
                }
                node_names.push_back(tokens[i].mText);
                const std::size_t line = tokens[i].mLine;
                ++i;
                for (int d = 0; d < dim; ++d) {
                    require_token("coordinates", line);
                    coords.push_back(ca_number(tokens[i]));
                    ++i;
                }
            }
            continue;
        }
        if (const CaTypeSpec* spec = ca_spec_by_keyword(kw)) {
            CaElementBlock block{spec, {}, {}, {}};
            skip_options(nullptr);
            while (true) {
                require_token("elements", head.mLine);
                if (ca_is_block_end(ca_upper(tokens[i].mText))) {
                    ++i;
                    break;
                }
                if (ca_is_option(tokens[i].mText)) {
                    ++i;
                    continue;
                }
                block.mNames.push_back(tokens[i].mText);
                block.mLines.push_back(tokens[i].mLine);
                const std::size_t line = tokens[i].mLine;
                ++i;
                for (std::size_t k = 0; k < spec->mNodes; ++k) {
                    require_token("element nodes", line);
                    if (ca_is_block_end(ca_upper(tokens[i].mText)))
                        ca_fail(std::string(spec->mKeyword) + " element '" + block.mNames.back() +
                                    "' has fewer than " + std::to_string(spec->mNodes) + " nodes",
                                line);
                    block.mNodes.push_back(tokens[i].mText);
                    ++i;
                }
            }
            if (!block.mNames.empty())
                blocks.push_back(std::move(block));
            continue;
        }
        if (kw == "GROUP_NO" || kw == "GROUP_MA") {
            CaGroup group{{}, kw == "GROUP_MA", {}};
            skip_options(&group.mName);
            while (true) {
                require_token("a group", head.mLine);
                if (ca_is_block_end(ca_upper(tokens[i].mText))) {
                    ++i;
                    break;
                }
                if (ca_is_option(tokens[i].mText)) {
                    ++i;
                    continue;
                }
                if (group.mName.empty())
                    group.mName = tokens[i].mText;
                else
                    group.mMembers.emplace_back(tokens[i].mText, tokens[i].mLine);
                ++i;
            }
            if (group.mName.empty())
                ca_fail(kw + " block without a name", head.mLine);
            groups.push_back(std::move(group));
            continue;
        }
        if (ca_is_keyword(kw)) {
            // Code_Aster skips its title, debug switches and inactive keywords;
            // anything else (another cell type, a typo) is skipped with a warning.
            if (!ca_skipped_keywords().count(kw) && !warned_keywords[kw]) {
                log::warn("Code_Aster .mail: skipping the unsupported '{}' block (line {})",
                          head.mText, head.mLine);
                warned_keywords[kw] = true;
            }
            while (i < n && !ca_is_block_end(ca_upper(tokens[i].mText)))
                ++i;
            if (i == n)
                ca_fail("block '" + head.mText + "' has no FINSF", head.mLine);
            ++i;
            continue;
        }
        ca_fail("expected a keyword, found '" + head.mText + "'", head.mLine);
    }
    if (!saw_fin)
        log::warn("Code_Aster .mail: '{}' has no FIN line; reading it to the end", rPath);

    // --- points -----------------------------------------------------------------
    std::unordered_map<std::string, std::int64_t> node_index;
    node_index.reserve(node_names.size());
    for (std::size_t p = 0; p < node_names.size(); ++p)
        if (!node_index.emplace(node_names[p], static_cast<std::int64_t>(p)).second)
            throw ReadError("Code_Aster .mail: node '" + node_names[p] + "' is defined twice");
    Mesh mesh;
    const std::size_t pdim = point_dim == 0 ? 3 : static_cast<std::size_t>(point_dim);
    NDArray points(DType::Float64, {node_names.size(), pdim});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // --- cells ------------------------------------------------------------------
    std::unordered_map<std::string, std::int64_t> element_index;
    std::vector<int> element_dim;
    std::int64_t global = 0;
    for (const CaElementBlock& b : blocks) {
        const std::size_t k = b.mSpec->mNodes;
        const std::size_t rows = b.mNames.size();
        const detail::NodeOrder* order = detail::node_order("code_aster", b.mSpec->mType);
        const int dim = cell_type_dimension(cell_type_from_name(b.mSpec->mType));
        NDArray conn(DType::Int64, {rows, k});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < rows; ++r) {
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mToMeshio[j]) : j;
                const std::string& name = b.mNodes[r * k + src];
                const auto it = node_index.find(name);
                if (it == node_index.end())
                    ca_fail("element '" + b.mNames[r] + "' names undefined node '" + name + "'",
                            b.mLines[r]);
                c[r * k + j] = it->second;
            }
            if (!element_index.emplace(b.mNames[r], global).second)
                ca_fail("element '" + b.mNames[r] + "' is defined twice", b.mLines[r]);
            element_dim.push_back(dim);
            ++global;
        }
        mesh.AddCellBlock(b.mSpec->mType, std::move(conn));
    }

    // --- groups -----------------------------------------------------------------
    // Two blocks with one name and kind are one group, as far as a Region can tell.
    std::map<std::pair<bool, std::string>, std::vector<std::int64_t>> members;
    std::map<std::pair<bool, std::string>, int> dims;
    std::vector<std::pair<bool, std::string>> order;
    for (const CaGroup& g : groups) {
        const auto key = std::make_pair(g.mCells, g.mName);
        auto [it, fresh] = members.emplace(key, std::vector<std::int64_t>{});
        if (fresh) {
            order.push_back(key);
            dims[key] = -1;
        } else {
            log::warn("Code_Aster .mail: {} '{}' is defined twice; merging the two",
                      g.mCells ? "GROUP_MA" : "GROUP_NO", g.mName);
        }
        std::unordered_set<std::int64_t> seen(it->second.begin(), it->second.end());
        bool warned_duplicate = false;
        for (const auto& [name, line] : g.mMembers) {
            const auto& index = g.mCells ? element_index : node_index;
            const auto found = index.find(name);
            if (found == index.end())
                ca_fail(std::string(g.mCells ? "GROUP_MA" : "GROUP_NO") + " '" + g.mName +
                            "' names undefined " + (g.mCells ? "element" : "node") + " '" + name +
                            "'",
                        line);
            if (!seen.insert(found->second).second) {
                if (!warned_duplicate)
                    log::warn("Code_Aster .mail: {} '{}' lists '{}' more than once",
                              g.mCells ? "GROUP_MA" : "GROUP_NO", g.mName, name);
                warned_duplicate = true;
                continue;
            }
            it->second.push_back(found->second);
            if (g.mCells)
                dims[key] =
                    std::max(dims[key], element_dim[static_cast<std::size_t>(found->second)]);
        }
    }
    for (const auto& key : order) {
        const std::vector<std::int64_t>& ids = members[key];
        NDArray entries(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), entries.As<std::int64_t>());
        mesh.AddRegion(Region(key.second, key.first ? RegionKind::Cell : RegionKind::Point,
                              key.first ? dims[key] : -1, -1, std::move(entries)));
    }
    return mesh;
}

namespace {

// Appends `rToken` to the record being written, wrapping to an indented
// continuation line rather than crossing column 80.
void ca_append(std::string& rOut, std::size_t& rColumn, const std::string& rToken) {
    if (rColumn == 0) {
        rOut += rToken;
        rColumn = rToken.size();
        return;
    }
    if (rColumn + 1 + rToken.size() > kCaColumns) {
        rOut += "\n        ";
        rColumn = 8;
    }
    rOut += ' ';
    rOut += rToken;
    rColumn += 1 + rToken.size();
}

void ca_end_record(std::string& rOut, std::size_t& rColumn) {
    rOut += '\n';
    rColumn = 0;
}

// Letters, digits and `_`, at most 24 characters, unique within its kind.
std::string ca_group_name(const std::string& rName, std::map<std::string, int>& rTaken,
                          const char* pKind) {
    std::string clean;
    // An explicit ASCII test, not std::isalnum: macOS's C locale calls some bytes
    // above 0x7F alphanumeric, which would let UTF-8 bytes through.
    for (char c : rName) {
        const bool keep =
            (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_';
        clean.push_back(keep ? c : '_');
    }
    if (clean.empty())
        clean = "GROUP";
    if (clean.size() > kCaMaxGroupName)
        clean.resize(kCaMaxGroupName);
    std::string out = clean;
    for (int k = 1; rTaken.count(out); ++k) {
        const std::string suffix = "_" + std::to_string(k);
        out = clean.substr(0, std::min(clean.size(), kCaMaxGroupName - suffix.size())) + suffix;
    }
    rTaken[out] = 1;
    if (out != rName)
        log::warn("Code_Aster .mail: {} '{}' is written as '{}'", pKind, rName, out);
    return out;
}

}  // namespace

void write_code_aster(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t pdim = rMesh.PointDim();
    if (pdim > 3)
        throw WriteError("Code_Aster .mail writer: points of dimension " + std::to_string(pdim) +
                         " (at most 3)");
    std::size_t num_cells = 0;
    std::vector<const CaTypeSpec*> specs;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const CaTypeSpec* spec = cb.IsRagged() ? nullptr : ca_spec_by_type(std::string(cb.Type()));
        if (!spec)
            throw WriteError("Code_Aster .mail writer: unsupported cell type '" +
                             std::string(cb.Type()) + "'");
        specs.push_back(spec);
        num_cells += cb.NumCells();
    }
    if (rMesh.NumPoints() > kCaMaxEntities || num_cells > kCaMaxEntities)
        throw WriteError(
            "Code_Aster .mail writer: more than 9,999,999 nodes or elements do "
            "not fit 8-character names; write MED instead");

    // Notes first: they are rendered into the provenance block.
    std::size_t side_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mKind == RegionKind::Side)
            ++side_regions;
    if (side_regions) {
        log::warn("Code_Aster .mail has no facet groups; {} side region(s) dropped", side_regions);
        detail::provenance_note(
            "regions-dropped",
            std::to_string(side_regions) + " side region(s) have no Code_Aster .mail group");
    }
    if (rMesh.NumPointData() + rMesh.NumCellData() + rMesh.NumFieldData() > 0) {
        log::warn("Code_Aster .mail holds no data arrays; point, cell and field data dropped");
        detail::provenance_note("data-dropped", "a Code_Aster .mail mesh holds no data arrays");
    }

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    std::string out = detail::provenance_render_lines(detail::SlotTier::Block, "% ");
    std::size_t column = 0;
    char buf[64];

    out += "COOR_" + std::to_string(pdim == 0 ? 3 : pdim) + "D\n";
    const NDArray& points = rMesh.Points();
    const std::size_t npts = rMesh.NumPoints();
    for (std::size_t p = 0; p < npts; ++p) {
        ca_append(out, column, "N" + std::to_string(p + 1));
        for (std::size_t d = 0; d < pdim; ++d) {
            detail::snprintf_c(buf, sizeof(buf), "%.16E",
                               detail::read_double(points, p * pdim + d));
            ca_append(out, column, buf);
        }
        ca_end_record(out, column);
    }
    out += "FINSF\n";
    f << out;
    out.clear();

    std::size_t label = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const CaTypeSpec* spec = specs[b];
        const NDArray& conn = cb.Conn();
        const std::size_t k = spec->mNodes;
        const detail::NodeOrder* order = detail::node_order("code_aster", spec->mType);
        out += spec->mKeyword;
        out += '\n';
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            ca_append(out, column, "M" + std::to_string(++label));
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t src = order ? static_cast<std::size_t>(order->mFromMeshio[j]) : j;
                ca_append(out, column,
                          "N" + std::to_string(detail::read_int(conn, r * k + src) + 1));
            }
            ca_end_record(out, column);
        }
        out += "FINSF\n";
        f << out;
        out.clear();
    }

    std::map<std::string, int> taken_no, taken_ma;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Side)
            continue;
        const bool cells = reg.mKind == RegionKind::Cell;
        const std::string name =
            ca_group_name(reg.mName, cells ? taken_ma : taken_no, cells ? "GROUP_MA" : "GROUP_NO");
        out += cells ? "GROUP_MA NOM=" : "GROUP_NO NOM=";
        out += name;
        out += '\n';
        const std::int64_t* e = reg.Entries();
        for (std::size_t j = 0; j < reg.NumEntries(); ++j)
            ca_append(out, column, (cells ? "M" : "N") + std::to_string(e[j] + 1));
        if (column)
            ca_end_record(out, column);
        out += "FINSF\n";
        f << out;
        out.clear();
    }
    f << "FIN\n";
    if (!f)
        throw WriteError("Code_Aster .mail writer: failed writing " + rPath);
}

}  // namespace meshioplusplus
