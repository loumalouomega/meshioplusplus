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
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <system_error>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/lsdyna.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;
using detail::CardField;
using detail::CardMode;

constexpr int lsd_max_include_depth = 8;
constexpr std::int64_t lsd_max_std_id = 99999999;

enum class LsdFamily { Solid, TShell, Shell, Beam, Discrete, Mass, Node, Part, Segment };

// Meshio++ cell types an LS-DYNA element card can produce, in a fixed order.
enum class LsdType { Vertex, Line, Triangle, Quad, Tetra, Pyramid, Wedge, Hexahedron, Tetra10 };

const char* lsd_type_name(LsdType Type) {
    switch (Type) {
        case LsdType::Vertex:
            return "vertex";
        case LsdType::Line:
            return "line";
        case LsdType::Triangle:
            return "triangle";
        case LsdType::Quad:
            return "quad";
        case LsdType::Tetra:
            return "tetra";
        case LsdType::Pyramid:
            return "pyramid";
        case LsdType::Wedge:
            return "wedge";
        case LsdType::Hexahedron:
            return "hexahedron";
        case LsdType::Tetra10:
            return "tetra10";
    }
    return "";
}

std::size_t lsd_type_nodes(LsdType Type) {
    switch (Type) {
        case LsdType::Vertex:
            return 1;
        case LsdType::Line:
            return 2;
        case LsdType::Triangle:
            return 3;
        case LsdType::Quad:
            return 4;
        case LsdType::Tetra:
            return 4;
        case LsdType::Pyramid:
            return 5;
        case LsdType::Wedge:
            return 6;
        case LsdType::Hexahedron:
            return 8;
        case LsdType::Tetra10:
            return 10;
    }
    return 0;
}

int lsd_family_dim(LsdFamily Family) {
    switch (Family) {
        case LsdFamily::Solid:
        case LsdFamily::TShell:
            return 3;
        case LsdFamily::Shell:
            return 2;
        case LsdFamily::Beam:
        case LsdFamily::Discrete:
            return 1;
        default:
            return 0;
    }
}

const char* lsd_family_name(LsdFamily Family) {
    switch (Family) {
        case LsdFamily::Solid:
            return "solid";
        case LsdFamily::TShell:
            return "tshell";
        case LsdFamily::Shell:
            return "shell";
        case LsdFamily::Beam:
            return "beam";
        case LsdFamily::Discrete:
            return "discrete";
        case LsdFamily::Mass:
            return "mass";
        case LsdFamily::Node:
            return "node";
        case LsdFamily::Part:
            return "part";
        case LsdFamily::Segment:
            return "segment";
    }
    return "";
}

std::string lsd_upper(std::string Text) {
    for (char& c : Text)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return Text;
}

std::string lsd_strip(const std::string& rText, const char* pChars = " \t\r\n\f\v") {
    const std::size_t b = rText.find_first_not_of(pChars);
    if (b == std::string::npos)
        return std::string();
    const std::size_t e = rText.find_last_not_of(pChars);
    return rText.substr(b, e - b + 1);
}

bool lsd_starts_with(const std::string& rText, const char* pPrefix) {
    return rText.rfind(pPrefix, 0) == 0;
}

std::vector<std::string> lsd_split(const std::string& rText, char Sep) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t k = rText.find(Sep, start);
        if (k == std::string::npos) {
            out.push_back(rText.substr(start));
            return out;
        }
        out.push_back(rText.substr(start, k - start));
        start = k + 1;
    }
}

// -- card layouts ---------------------------------------------------------------

const std::vector<CardField>& lsd_layout_node() {
    static const std::vector<CardField> l = {{'i', 8},  {'r', 16}, {'r', 16},
                                             {'r', 16}, {'r', 8},  {'r', 8}};
    return l;
}

const std::vector<CardField>& lsd_layout_element() {
    static const std::vector<CardField> l(10, CardField{'i', 8});
    return l;
}

const std::vector<CardField>& lsd_layout_mass() {
    static const std::vector<CardField> l = {{'i', 8}, {'i', 8}, {'r', 16}, {'i', 8}};
    return l;
}

const std::vector<CardField>& lsd_layout_ids() {
    static const std::vector<CardField> l(8, CardField{'i', 10});
    return l;
}

const std::vector<CardField>& lsd_layout_segment() {
    static const std::vector<CardField> l = {{'i', 10}, {'i', 10}, {'i', 10}, {'i', 10},
                                             {'r', 10}, {'r', 10}, {'r', 10}, {'r', 10}};
    return l;
}

// -- degenerate hexahedra -------------------------------------------------------

// The cell type and meshio++-ordered nodes of an 8-node LS-DYNA solid. There is no
// tetra, pyramid or wedge card: they are hexahedra with repeated nodes, checked
// most-degenerate first. Twin of `_collapse_solid` in lsdyna/_lsdyna.py.
std::pair<LsdType, std::vector<std::int64_t>> lsd_collapse_solid(
    const std::array<std::int64_t, 8>& n) {
    if (n[3] == n[4] && n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {LsdType::Tetra, {n[0], n[1], n[2], n[3]}};
    if (n[2] == n[3] && n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {LsdType::Tetra, {n[0], n[1], n[2], n[4]}};
    if (n[4] == n[5] && n[5] == n[6] && n[6] == n[7])
        return {LsdType::Pyramid, {n[0], n[1], n[2], n[3], n[4]}};
    if (n[2] == n[3] && n[6] == n[7])
        return {LsdType::Wedge, {n[0], n[1], n[2], n[4], n[5], n[6]}};
    if (n[4] == n[5] && n[6] == n[7])
        return {LsdType::Wedge, {n[0], n[4], n[1], n[3], n[6], n[2]}};
    return {LsdType::Hexahedron, std::vector<std::int64_t>(n.begin(), n.end())};
}

// The 8 LS-DYNA nodes of a meshio++ tetra, pyramid, wedge or hexahedron.
std::array<std::int64_t, 8> lsd_expand_solid(LsdType Type, const std::int64_t* pRow) {
    const std::int64_t* r = pRow;
    switch (Type) {
        case LsdType::Tetra:
            return {r[0], r[1], r[2], r[3], r[3], r[3], r[3], r[3]};
        case LsdType::Pyramid:
            return {r[0], r[1], r[2], r[3], r[4], r[4], r[4], r[4]};
        case LsdType::Wedge:
            return {r[0], r[2], r[5], r[3], r[1], r[1], r[4], r[4]};
        default:
            return {r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]};
    }
}

// -- deck state -----------------------------------------------------------------

struct LsdElement {
    LsdFamily mFamily;
    std::int64_t mEid;
    std::int64_t mPid;
    LsdType mType;
    std::array<std::int64_t, 10> mNodes;
    int mGroup;  // one per element keyword, so blocks follow the deck's sections
};

struct LsdSet {
    LsdFamily mFamily;
    std::int64_t mSid;
    std::string mTitle;
    std::vector<std::int64_t> mIds;
    std::vector<std::array<std::int64_t, 4>> mSegments;
};

struct LsdLine {
    std::size_t mLineNo;
    std::string mText;
};

struct LsdKeyword {
    std::string mName;
    std::optional<CardMode> mMode;
    std::string mRest;
};

struct LsdDeck {
    std::unordered_map<std::int64_t, std::int64_t> mNodeIndex;
    std::vector<double> mCoords;  // x, y, z per node
    std::vector<LsdElement> mElements;
    int mGroup = 0;
    std::vector<std::pair<std::int64_t, std::string>> mParts;  // pid, title; first-seen order
    std::unordered_map<std::int64_t, std::size_t> mPartSlot;
    std::vector<LsdSet> mSets;
    std::vector<std::pair<std::string, fs::path>> mIncludeDirs;
    std::unordered_set<std::string> mWarned;
    std::size_t mParamSkips = 0;
};

using LsdBlock = std::vector<LsdLine>;

struct LsdCtx {
    const std::string& mLabel;
    CardMode mMode;
};

std::string lsd_where(std::size_t LineNo, const std::string& rLabel) {
    return " (line " + std::to_string(LineNo) + " of " + rLabel + ")";
}

void lsd_warn_once(LsdDeck& rDeck, const std::string& rKey, const std::string& rMessage) {
    if (rDeck.mWarned.insert(rKey).second)
        log::warn("{}", rMessage);
}

LsdKeyword lsd_parse_keyword(const std::string& rLine) {
    const std::string body = lsd_strip(rLine.substr(1));
    std::size_t end = body.size();
    for (std::size_t k = 0; k < body.size(); ++k) {
        if (body[k] == ' ' || body[k] == '\t' || body[k] == ',') {
            end = k;
            break;
        }
    }
    LsdKeyword kw;
    std::string up = lsd_upper(body.substr(0, end));
    kw.mRest = lsd_strip(body.substr(end), " \t,");
    auto mode_of = [](char c) {
        return c == '+' ? CardMode::Long : (c == '-' ? CardMode::Standard : CardMode::I10);
    };
    if (!up.empty() && (up.back() == '+' || up.back() == '-' || up.back() == '%')) {
        kw.mMode = mode_of(up.back());
        up.pop_back();
    } else if (kw.mRest == "+" || kw.mRest == "-" || kw.mRest == "%") {
        kw.mMode = mode_of(kw.mRest[0]);
        kw.mRest.clear();
    }
    kw.mName = up;
    return kw;
}

// The deck-wide default a `*KEYWORD` card sets (`LONG=Y|S|K`, `I10=Y`).
CardMode lsd_keyword_mode(const std::string& rRest, CardMode Mode) {
    std::string text = rRest;
    std::replace(text.begin(), text.end(), ',', ' ');
    std::size_t pos = 0;
    while (pos < text.size()) {
        const std::size_t b = text.find_first_not_of(" \t\r\n\f\v", pos);
        if (b == std::string::npos)
            break;
        std::size_t e = text.find_first_of(" \t\r\n\f\v", b);
        if (e == std::string::npos)
            e = text.size();
        const std::string token = text.substr(b, e - b);
        pos = e;
        const std::size_t eq = token.find('=');
        const std::string key = lsd_upper(token.substr(0, eq));
        const std::string value = eq == std::string::npos ? "" : lsd_upper(token.substr(eq + 1));
        if (key == "LONG" && (value == "Y" || value == "S" || value == "K"))
            Mode = CardMode::Long;
        else if (key == "I10" && value == "Y" && Mode != CardMode::Long)
            Mode = CardMode::I10;
    }
    return Mode;
}

bool lsd_skip_param(LsdDeck& rDeck, const std::string& rLine) {
    if (rLine.find('&') != std::string::npos) {
        ++rDeck.mParamSkips;
        return true;
    }
    return false;
}

std::size_t lsd_skip_pgp(LsdDeck& rDeck, const std::vector<std::string>& rLines, std::size_t Pos) {
    lsd_warn_once(rDeck, "pgp", "LS-DYNA: skipped a PGP-encrypted block");
    while (Pos < rLines.size() && !lsd_starts_with(rLines[Pos], "-----END PGP"))
        ++Pos;
    return Pos + 1;
}

std::int64_t lsd_int(const std::vector<std::string>& rFields, std::size_t I,
                     const std::string& rWhere) {
    return detail::card_to_int(I < rFields.size() ? rFields[I] : std::string(), rWhere);
}

double lsd_real(const std::vector<std::string>& rFields, std::size_t I, const std::string& rWhere) {
    return detail::card_to_real(I < rFields.size() ? rFields[I] : std::string(), rWhere);
}

void lsd_read_nodes(LsdDeck& rDeck, const LsdBlock& rBlock, const LsdCtx& rCtx) {
    for (const LsdLine& line : rBlock) {
        if (lsd_strip(line.mText).empty() || lsd_skip_param(rDeck, line.mText))
            continue;
        const std::string where = lsd_where(line.mLineNo, rCtx.mLabel);
        const auto f = detail::split_card(line.mText, lsd_layout_node(), rCtx.mMode);
        const std::int64_t nid = lsd_int(f, 0, where);
        const std::int64_t index = static_cast<std::int64_t>(rDeck.mCoords.size() / 3);
        if (!rDeck.mNodeIndex.emplace(nid, index).second)
            throw ReadError("LS-DYNA: duplicate node id " + std::to_string(nid) + where);
        rDeck.mCoords.push_back(lsd_real(f, 1, where));
        rDeck.mCoords.push_back(lsd_real(f, 2, where));
        rDeck.mCoords.push_back(lsd_real(f, 3, where));
    }
}

void lsd_read_elements(LsdDeck& rDeck, const std::string& rKeyword, const LsdBlock& rBlock,
                       const LsdCtx& rCtx) {
    const std::vector<std::string> tokens = lsd_split(rKeyword, '_');  // ELEMENT, KIND, opts...
    const std::string kind = tokens.size() > 1 ? tokens[1] : std::string();
    const std::vector<std::string> opts(tokens.begin() + std::min<std::size_t>(2, tokens.size()),
                                        tokens.end());
    std::size_t extras = 0;
    if (kind == "SOLID" && opts.empty()) {
    } else if (kind == "SOLID" && opts.size() == 1 && opts[0] == "ORTHO") {
        extras = 2;
    } else if (kind == "SOLID" && !opts.empty() &&
               (opts[0] == "TET4TOTET10" || opts[0] == "H8TOH20")) {
        lsd_warn_once(rDeck, rKeyword,
                      "LS-DYNA: *" + rKeyword + " is a conversion directive and is not applied");
        return;
    } else if (kind == "SHELL" && std::all_of(opts.begin(), opts.end(), [](const std::string& o) {
                   return o == "THICKNESS" || o == "BETA" || o == "MCID" || o == "OFFSET";
               })) {
        extras = opts.size();
    } else if ((kind == "TSHELL" || kind == "BEAM" || kind == "DISCRETE" || kind == "MASS") &&
               opts.empty()) {
    } else if (kind == "SOLID" || kind == "SHELL" || kind == "TSHELL" || kind == "BEAM" ||
               kind == "DISCRETE" || kind == "MASS") {
        lsd_warn_once(rDeck, rKeyword, "LS-DYNA: *" + rKeyword + " is not supported; skipped");
        return;
    } else {
        return;
    }
    ++rDeck.mGroup;
    LsdFamily family = LsdFamily::Solid;
    if (kind == "SHELL")
        family = LsdFamily::Shell;
    else if (kind == "TSHELL")
        family = LsdFamily::TShell;
    else if (kind == "BEAM")
        family = LsdFamily::Beam;
    else if (kind == "DISCRETE")
        family = LsdFamily::Discrete;
    else if (kind == "MASS")
        family = LsdFamily::Mass;

    std::size_t j = 0;
    while (j < rBlock.size()) {
        const LsdLine& line = rBlock[j];
        ++j;
        if (lsd_strip(line.mText).empty())
            continue;
        const std::string where = lsd_where(line.mLineNo, rCtx.mLabel);
        if (lsd_skip_param(rDeck, line.mText)) {
            j += extras;
            continue;
        }
        LsdElement el{};
        el.mFamily = family;
        el.mGroup = rDeck.mGroup;
        if (family == LsdFamily::Mass) {
            const auto f = detail::split_card(line.mText, lsd_layout_mass(), rCtx.mMode);
            el.mEid = lsd_int(f, 0, where);
            el.mPid = lsd_int(f, 3, where);
            el.mType = LsdType::Vertex;
            el.mNodes[0] = lsd_int(f, 1, where);
            rDeck.mElements.push_back(el);
            continue;
        }
        auto f = detail::split_card(line.mText, lsd_layout_element(), rCtx.mMode);
        el.mEid = lsd_int(f, 0, where);
        el.mPid = lsd_int(f, 1, where);
        if (family == LsdFamily::Solid || family == LsdFamily::TShell) {
            bool two_line = false;
            if (family == LsdFamily::Solid) {
                two_line = true;
                for (std::size_t k = 2; k < 10 && k < f.size(); ++k)
                    if (!f[k].empty())
                        two_line = false;
            }
            if (two_line) {
                // Two-line layout: "eid pid" then up to ten node fields.
                if (j >= rBlock.size())
                    throw ReadError("LS-DYNA: truncated element card" + where);
                const LsdLine& line2 = rBlock[j];
                ++j;
                const auto g = detail::split_card(line2.mText, lsd_layout_element(), rCtx.mMode);
                std::vector<std::int64_t> nodes;
                for (std::size_t k = 0; k < 10; ++k)
                    nodes.push_back(lsd_int(g, k, where));
                while (!nodes.empty() && nodes.back() == 0)
                    nodes.pop_back();
                if (nodes.size() == 10) {
                    el.mType = LsdType::Tetra10;
                    std::copy(nodes.begin(), nodes.end(), el.mNodes.begin());
                } else if (nodes.size() == 8) {
                    std::array<std::int64_t, 8> n;
                    std::copy(nodes.begin(), nodes.end(), n.begin());
                    auto [type, ordered] = lsd_collapse_solid(n);
                    el.mType = type;
                    std::copy(ordered.begin(), ordered.end(), el.mNodes.begin());
                } else if (nodes.size() == 4) {
                    el.mType = LsdType::Tetra;
                    std::copy(nodes.begin(), nodes.end(), el.mNodes.begin());
                } else {
                    throw ReadError("LS-DYNA: solid element with " + std::to_string(nodes.size()) +
                                    " nodes is not supported" + where);
                }
            } else {
                std::array<std::int64_t, 8> n;
                for (std::size_t k = 0; k < 8; ++k)
                    n[k] = lsd_int(f, 2 + k, where);
                auto [type, ordered] = lsd_collapse_solid(n);
                el.mType = type;
                std::copy(ordered.begin(), ordered.end(), el.mNodes.begin());
            }
        } else if (family == LsdFamily::Shell) {
            std::array<std::int64_t, 8> n;
            for (std::size_t k = 0; k < 8; ++k)
                n[k] = lsd_int(f, 2 + k, where);
            if (n[4] || n[5] || n[6] || n[7])
                lsd_warn_once(rDeck, "shell-mid",
                              "LS-DYNA: higher-order shell nodes (n5..n8) are ignored");
            if (n[3] == 0 || n[3] == n[2]) {
                el.mType = LsdType::Triangle;
                std::copy(n.begin(), n.begin() + 3, el.mNodes.begin());
            } else {
                el.mType = LsdType::Quad;
                std::copy(n.begin(), n.begin() + 4, el.mNodes.begin());
            }
        } else {  // beam, discrete
            el.mType = LsdType::Line;
            el.mNodes[0] = lsd_int(f, 2, where);
            el.mNodes[1] = lsd_int(f, 3, where);
            if (family == LsdFamily::Beam) {
                // An optional second card; element cards have no ".".
                while (j < rBlock.size() && rBlock[j].mText.find('.') != std::string::npos)
                    ++j;
            }
        }
        rDeck.mElements.push_back(el);
        j += extras;
    }
}

void lsd_read_parts(LsdDeck& rDeck, const std::string& rKeyword, const LsdBlock& rBlock,
                    const LsdCtx& rCtx) {
    std::size_t j = 0;
    while (j + 1 < rBlock.size()) {
        const std::string title = lsd_strip(rBlock[j].mText);
        const LsdLine& card = rBlock[j + 1];
        if (lsd_strip(card.mText).empty())
            break;
        const std::string where = lsd_where(card.mLineNo, rCtx.mLabel);
        const auto f = detail::split_card(card.mText, lsd_layout_ids(), rCtx.mMode);
        const std::int64_t pid = lsd_int(f, 0, where);
        auto slot = rDeck.mPartSlot.find(pid);
        if (slot == rDeck.mPartSlot.end()) {
            rDeck.mPartSlot[pid] = rDeck.mParts.size();
            rDeck.mParts.emplace_back(pid, title);
        } else {
            rDeck.mParts[slot->second].second = title;
        }
        j += 2;
        if (rKeyword != "PART")
            break;
    }
}

std::optional<LsdFamily> lsd_set_family(const std::string& rName) {
    if (rName == "NODE")
        return LsdFamily::Node;
    if (rName == "SOLID")
        return LsdFamily::Solid;
    if (rName == "SHELL")
        return LsdFamily::Shell;
    if (rName == "TSHELL")
        return LsdFamily::TShell;
    if (rName == "BEAM")
        return LsdFamily::Beam;
    if (rName == "DISCRETE")
        return LsdFamily::Discrete;
    if (rName == "PART")
        return LsdFamily::Part;
    if (rName == "SEGMENT")
        return LsdFamily::Segment;
    return std::nullopt;
}

void lsd_read_set(LsdDeck& rDeck, const std::string& rKeyword, const LsdBlock& rBlock,
                  const LsdCtx& rCtx) {
    const std::vector<std::string> tokens = lsd_split(rKeyword, '_');  // SET, FAMILY, opts...
    const auto family = lsd_set_family(tokens.size() > 1 ? tokens[1] : std::string());
    if (!family)
        return;
    std::set<std::string> opts;
    for (std::size_t k = 2; k < tokens.size(); ++k)
        opts.insert(tokens[k]);
    for (const std::string& o : opts) {
        if (o != "LIST" && o != "GENERATE" && o != "TITLE") {
            lsd_warn_once(rDeck, rKeyword, "LS-DYNA: *" + rKeyword + " is not supported; skipped");
            return;
        }
    }
    LsdSet set;
    set.mFamily = *family;
    std::size_t j = 0;
    if (opts.count("TITLE")) {
        if (rBlock.empty())
            return;
        set.mTitle = lsd_strip(rBlock[0].mText);
        j = 1;
    }
    if (j >= rBlock.size())
        return;
    {
        const std::string where = lsd_where(rBlock[j].mLineNo, rCtx.mLabel);
        set.mSid =
            lsd_int(detail::split_card(rBlock[j].mText, lsd_layout_ids(), rCtx.mMode), 0, where);
        ++j;
    }
    const bool generate = opts.count("GENERATE") > 0;
    for (; j < rBlock.size(); ++j) {
        const LsdLine& line = rBlock[j];
        if (lsd_strip(line.mText).empty() || lsd_skip_param(rDeck, line.mText))
            continue;
        const std::string where = lsd_where(line.mLineNo, rCtx.mLabel);
        if (*family == LsdFamily::Segment) {
            const auto f = detail::split_card(line.mText, lsd_layout_segment(), rCtx.mMode);
            std::array<std::int64_t, 4> seg;
            for (std::size_t k = 0; k < 4; ++k)
                seg[k] = lsd_int(f, k, where);
            if (seg[3] == 0)
                seg[3] = seg[2];
            set.mSegments.push_back(seg);
            continue;
        }
        const auto fields = detail::split_card(line.mText, lsd_layout_ids(), rCtx.mMode);
        std::vector<std::int64_t> f;
        for (std::size_t k = 0; k < fields.size(); ++k)
            f.push_back(lsd_int(fields, k, where));
        if (generate) {
            for (std::size_t k = 0; k + 1 < f.size(); k += 2)
                if (f[k] > 0 && f[k + 1] >= f[k])
                    for (std::int64_t v = f[k]; v <= f[k + 1]; ++v)
                        set.mIds.push_back(v);
        } else {
            for (std::int64_t v : f)
                if (v != 0)
                    set.mIds.push_back(v);
        }
    }
    rDeck.mSets.push_back(std::move(set));
}

void lsd_read_text(LsdDeck& rDeck, const std::string& rText, const fs::path& rBaseDir,
                   const std::string& rLabel, int Depth, CardMode Mode);

void lsd_read_file(LsdDeck& rDeck, const fs::path& rPath, int Depth, CardMode Mode) {
    if (Depth > lsd_max_include_depth)
        throw ReadError("LS-DYNA: *INCLUDE nested deeper than " +
                        std::to_string(lsd_max_include_depth));
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("LS-DYNA: could not read " + rPath.string());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::error_code ec;
    const fs::path absolute = fs::absolute(rPath, ec);
    lsd_read_text(rDeck, text, (ec ? rPath : absolute).parent_path(), rPath.string(), Depth, Mode);
}

std::optional<fs::path> lsd_find_include(const LsdDeck& rDeck, const std::string& rName,
                                         const fs::path& rBaseDir) {
    std::vector<std::string> spellings = {rName};
    if (rName.find('\\') != std::string::npos) {
        std::string alt = rName;
        std::replace(alt.begin(), alt.end(), '\\', '/');
        spellings.push_back(alt);
    }
    std::error_code ec;
    for (const std::string& sp : spellings) {
        const fs::path p(sp);
        std::vector<fs::path> cands;
        if (p.is_absolute()) {
            cands.push_back(p);
        } else {
            cands.push_back(rBaseDir / p);
            for (const auto& [d, parent] : rDeck.mIncludeDirs) {
                const fs::path dp(d);
                cands.push_back(dp.is_absolute() ? dp / p : parent / dp / p);
                if (!dp.is_absolute())
                    cands.push_back(dp / p);
            }
            cands.push_back(fs::current_path(ec) / p);
        }
        for (const fs::path& c : cands)
            if (fs::is_regular_file(c, ec))
                return c;
    }
    return std::nullopt;
}

void lsd_read_includes(LsdDeck& rDeck, const std::string& rKeyword, const LsdBlock& rBlock,
                       const fs::path& rBaseDir, int Depth, CardMode Mode) {
    std::vector<std::string> names;
    std::optional<std::string> cur;
    for (const LsdLine& line : rBlock) {
        const std::string s = lsd_strip(line.mText);
        if (s.empty())
            continue;
        cur = cur ? *cur + s : s;
        if (!cur->empty() && cur->back() == '+') {
            cur->pop_back();
            continue;
        }
        names.push_back(*cur);
        cur.reset();
        if (rKeyword == "INCLUDE_TRANSFORM") {
            lsd_warn_once(rDeck, "transform",
                          "LS-DYNA: *INCLUDE_TRANSFORM offsets are not applied; the file is "
                          "included untransformed");
            break;
        }
    }
    if (cur && !cur->empty())
        names.push_back(*cur);
    for (const std::string& name : names) {
        const auto path = lsd_find_include(rDeck, name, rBaseDir);
        if (!path) {
            log::warn("{}", "LS-DYNA: include file '" + name + "' not found; skipped");
            continue;
        }
        lsd_read_file(rDeck, *path, Depth + 1, Mode);
    }
}

void lsd_read_text(LsdDeck& rDeck, const std::string& rText, const fs::path& rBaseDir,
                   const std::string& rLabel, int Depth, CardMode Mode) {
    std::vector<std::string> lines = lsd_split(rText, '\n');
    for (std::string& ln : lines)
        if (!ln.empty() && ln.back() == '\r')
            ln.pop_back();
    std::size_t pos = 0;
    const std::size_t n = lines.size();
    while (pos < n) {
        const std::string& line = lines[pos];
        ++pos;
        if (!lsd_starts_with(line, "*")) {
            if (lsd_starts_with(line, "-----BEGIN PGP"))
                pos = lsd_skip_pgp(rDeck, lines, pos);
            continue;
        }
        const LsdKeyword kw = lsd_parse_keyword(line);
        LsdBlock block;
        while (pos < n && !lsd_starts_with(lines[pos], "*")) {
            const std::string& raw = lines[pos];
            ++pos;
            if (lsd_starts_with(raw, "$"))
                continue;
            if (lsd_starts_with(raw, "-----BEGIN PGP")) {
                pos = lsd_skip_pgp(rDeck, lines, pos);
                continue;
            }
            block.push_back({pos, raw});
        }
        const std::string& name = kw.mName;
        if (name == "END")
            return;
        if (name == "KEYWORD") {
            Mode = lsd_keyword_mode(kw.mRest, Mode);
            continue;
        }
        const LsdCtx ctx{rLabel, kw.mMode ? *kw.mMode : Mode};
        if (name == "NODE")
            lsd_read_nodes(rDeck, block, ctx);
        else if (lsd_starts_with(name, "ELEMENT_"))
            lsd_read_elements(rDeck, name, block, ctx);
        else if (name == "PART" || name == "PART_CONTACT" || name == "PART_COMPOSITE" ||
                 name == "PART_INERTIA")
            lsd_read_parts(rDeck, name, block, ctx);
        else if (lsd_starts_with(name, "SET_"))
            lsd_read_set(rDeck, name, block, ctx);
        else if (name == "INCLUDE" || name == "INCLUDE_NO_TRANSFORM" || name == "INCLUDE_TRANSFORM")
            lsd_read_includes(rDeck, name, block, rBaseDir, Depth, Mode);
        else if (name == "INCLUDE_PATH" || name == "INCLUDE_PATH_RELATIVE") {
            for (const LsdLine& l : block) {
                const std::string d = lsd_strip(l.mText);
                if (!d.empty())
                    rDeck.mIncludeDirs.emplace_back(d, rBaseDir);
            }
        }
    }
}

// -- resolving ------------------------------------------------------------------

NDArray lsd_entries(const std::vector<std::int64_t>& rValues, std::size_t Stride) {
    NDArray out = Stride == 2 ? NDArray::Uninit(DType::Int64, {rValues.size() / 2, 2})
                              : NDArray::Uninit(DType::Int64, {rValues.size()});
    std::int64_t* p = out.As<std::int64_t>();
    for (std::size_t k = 0; k < rValues.size(); ++k)
        p[k] = rValues[k];
    return out;
}

struct LsdFaceKey {
    std::array<std::int64_t, 4> mNodes;
    bool operator==(const LsdFaceKey& rOther) const { return mNodes == rOther.mNodes; }
};

struct LsdFaceKeyHash {
    std::size_t operator()(const LsdFaceKey& rKey) const {
        std::size_t h = 1469598103934665603ULL;
        for (std::int64_t v : rKey.mNodes)
            h = (h ^ static_cast<std::size_t>(v + 2)) * 1099511628211ULL;
        return h;
    }
};

// Sorted corner nodes, padded with -1 for a triangle.
LsdFaceKey lsd_face_key(std::array<std::int64_t, 4> Nodes, std::size_t Count) {
    std::sort(Nodes.begin(), Nodes.begin() + static_cast<std::ptrdiff_t>(Count));
    for (std::size_t k = Count; k < 4; ++k)
        Nodes[k] = -1;
    return LsdFaceKey{Nodes};
}

using LsdFaceMap =
    std::unordered_map<LsdFaceKey, std::pair<std::int64_t, std::int64_t>, LsdFaceKeyHash>;

// Corner-node key -> (global cell, facet) for every face of every cell. An interior
// face is shared by two cells; the lowest cell index wins, and a shell element's own
// face is facet 0.
LsdFaceMap lsd_face_map(const Mesh& rMesh) {
    LsdFaceMap out;
    std::int64_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const auto& faces = detail::cell_faces(cell_type_from_name(type));
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const std::int64_t g = base + static_cast<std::int64_t>(r);
            if (!faces.empty()) {
                for (std::size_t f = 0; f < faces.size(); ++f) {
                    std::array<std::int64_t, 4> nodes{-1, -1, -1, -1};
                    for (std::size_t c = 0; c < faces[f].mNumCorners; ++c)
                        nodes[c] = detail::read_int(conn, r * k + faces[f].mNodes[c]);
                    out.emplace(lsd_face_key(nodes, faces[f].mNumCorners),
                                std::make_pair(g, static_cast<std::int64_t>(f)));
                }
            } else if (type == "triangle" || type == "quad") {
                std::array<std::int64_t, 4> nodes{-1, -1, -1, -1};
                for (std::size_t c = 0; c < k; ++c)
                    nodes[c] = detail::read_int(conn, r * k + c);
                out.emplace(lsd_face_key(nodes, k), std::make_pair(g, std::int64_t{0}));
            }
        }
        base += static_cast<std::int64_t>(cb.NumCells());
    }
    return out;
}

Mesh lsd_build_mesh(LsdDeck& rDeck) {
    if (rDeck.mParamSkips)
        log::warn("LS-DYNA: {} mesh card(s) use *PARAMETER values (&name) and were skipped",
                  rDeck.mParamSkips);
    Mesh mesh;
    const std::size_t npts = rDeck.mCoords.size() / 3;
    NDArray points(DType::Float64, {npts, std::size_t{3}});
    std::copy(rDeck.mCoords.begin(), rDeck.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // One block per (element keyword, cell type), so the blocks follow the deck's
    // sections; a keyword whose first type is the previous block's extends it.
    std::vector<LsdType> order;
    std::vector<std::vector<std::size_t>> rows;
    std::map<std::pair<int, int>, std::size_t> block_of;
    for (std::size_t e = 0; e < rDeck.mElements.size(); ++e) {
        const LsdElement& el = rDeck.mElements[e];
        const std::pair<int, int> key{el.mGroup, static_cast<int>(el.mType)};
        auto it = block_of.find(key);
        if (it == block_of.end()) {
            std::size_t b;
            if (!order.empty() && order.back() == el.mType) {
                b = order.size() - 1;
            } else {
                order.push_back(el.mType);
                rows.emplace_back();
                b = order.size() - 1;
            }
            it = block_of.emplace(key, b).first;
        }
        rows[it->second].push_back(e);
    }
    constexpr std::size_t nfam = 6;
    std::array<std::unordered_map<std::int64_t, std::int64_t>, nfam> owner;
    std::vector<std::int64_t> cell_pid;
    std::vector<LsdFamily> cell_family;
    for (std::size_t b = 0; b < order.size(); ++b) {
        const std::size_t count = lsd_type_nodes(order[b]);
        NDArray conn(DType::Int64, {rows[b].size(), count});
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < rows[b].size(); ++r) {
            const LsdElement& el = rDeck.mElements[rows[b][r]];
            for (std::size_t c = 0; c < count; ++c) {
                const auto it = rDeck.mNodeIndex.find(el.mNodes[c]);
                if (it == rDeck.mNodeIndex.end())
                    throw ReadError("LS-DYNA: element " + std::to_string(el.mEid) +
                                    " references undefined node " + std::to_string(el.mNodes[c]));
                cp[r * count + c] = it->second;
            }
            auto& map = owner[static_cast<std::size_t>(el.mFamily)];
            if (!map.emplace(el.mEid, static_cast<std::int64_t>(cell_pid.size())).second)
                throw ReadError("LS-DYNA: duplicate " + std::string(lsd_family_name(el.mFamily)) +
                                " element id " + std::to_string(el.mEid));
            cell_pid.push_back(el.mPid);
            cell_family.push_back(el.mFamily);
        }
        mesh.AddCellBlock(lsd_type_name(order[b]), std::move(conn));
    }

    std::unordered_map<std::int64_t, std::vector<std::int64_t>> by_pid;
    for (std::size_t g = 0; g < cell_pid.size(); ++g)
        by_pid[cell_pid[g]].push_back(static_cast<std::int64_t>(g));

    using SeenKey = std::tuple<int, std::string, int, std::int64_t>;
    std::set<SeenKey> seen;
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> part_cells;
    for (const auto& [pid, title] : rDeck.mParts) {
        const std::vector<std::int64_t>& members = by_pid[pid];
        int dim = -1;
        for (std::int64_t g : members)
            dim = std::max(dim, lsd_family_dim(cell_family[static_cast<std::size_t>(g)]));
        part_cells[pid] = members;
        const std::string name = title.empty() ? "Part " + std::to_string(pid) : title;
        seen.emplace(static_cast<int>(RegionKind::Cell), name, dim, pid);
        mesh.AddRegion(Region(name, RegionKind::Cell, dim, pid, lsd_entries(members, 1)));
    }

    std::size_t dropped = 0;
    std::optional<LsdFaceMap> face_map;
    for (const LsdSet& set : rDeck.mSets) {
        const std::string family = lsd_family_name(set.mFamily);
        std::string name = set.mTitle.empty()
                               ? lsd_upper(family) + " set " + std::to_string(set.mSid)
                               : set.mTitle;
        std::vector<std::int64_t> entries;
        RegionKind kind = RegionKind::Cell;
        if (set.mFamily == LsdFamily::Node) {
            kind = RegionKind::Point;
            for (std::int64_t id : set.mIds) {
                const auto it = rDeck.mNodeIndex.find(id);
                if (it != rDeck.mNodeIndex.end())
                    entries.push_back(it->second);
            }
            dropped += set.mIds.size() - entries.size();
        } else if (set.mFamily == LsdFamily::Part) {
            for (std::int64_t pid : set.mIds) {
                const auto it = part_cells.find(pid);
                if (it != part_cells.end())
                    entries.insert(entries.end(), it->second.begin(), it->second.end());
            }
        } else if (set.mFamily == LsdFamily::Segment) {
            kind = RegionKind::Side;
            if (!face_map)
                face_map = lsd_face_map(mesh);
            for (const auto& seg : set.mSegments) {
                std::array<std::int64_t, 4> idx;
                bool defined = true;
                for (std::size_t k = 0; k < 4; ++k) {
                    const auto it = rDeck.mNodeIndex.find(seg[k]);
                    idx[k] = it == rDeck.mNodeIndex.end() ? -1 : it->second;
                    defined = defined && it != rDeck.mNodeIndex.end();
                }
                const auto hit = defined
                                     ? face_map->find(lsd_face_key(idx, idx[3] == idx[2] ? 3 : 4))
                                     : face_map->end();
                if (!defined || hit == face_map->end()) {
                    ++dropped;
                } else {
                    entries.push_back(hit->second.first);
                    entries.push_back(hit->second.second);
                }
            }
        } else {
            const auto& map = owner[static_cast<std::size_t>(set.mFamily)];
            std::size_t found = 0;
            for (std::int64_t id : set.mIds) {
                const auto it = map.find(id);
                if (it != map.end()) {
                    entries.push_back(it->second);
                    ++found;
                }
            }
            dropped += set.mIds.size() - found;
        }
        const int kind_id = static_cast<int>(kind);
        if (seen.count({kind_id, name, -1, set.mSid})) {
            name += " [" + family + "]";
            if (seen.count({kind_id, name, -1, set.mSid}))
                continue;
        }
        seen.emplace(kind_id, name, -1, set.mSid);
        mesh.AddRegion(Region(name, kind, -1, set.mSid,
                              lsd_entries(entries, kind == RegionKind::Side ? 2 : 1)));
    }
    if (dropped)
        log::warn("LS-DYNA: {} set entries refer to undefined ids and were dropped", dropped);
    return mesh;
}

// -- writing --------------------------------------------------------------------

std::string lsd_title_line(const std::string& rName) {
    std::string t = rName;
    std::replace(t.begin(), t.end(), '\r', ' ');
    std::replace(t.begin(), t.end(), '\n', ' ');
    return (!t.empty() && (t[0] == '*' || t[0] == '$')) ? " " + t : t;
}

void lsd_put(std::string& rOut, const std::string& rText, std::size_t Width) {
    if (rText.size() < Width)
        rOut.append(Width - rText.size(), ' ');
    rOut += rText;
}

void lsd_put_int(std::string& rOut, std::int64_t Value, std::size_t Width) {
    lsd_put(rOut, std::to_string(Value), Width);
}

std::optional<LsdFamily> lsd_type_family(const std::string& rType) {
    if (rType == "vertex")
        return LsdFamily::Mass;
    if (rType == "line")
        return LsdFamily::Beam;
    if (rType == "triangle" || rType == "quad")
        return LsdFamily::Shell;
    if (rType == "tetra" || rType == "pyramid" || rType == "wedge" || rType == "hexahedron" ||
        rType == "tetra10")
        return LsdFamily::Solid;
    return std::nullopt;
}

LsdType lsd_type_of(const std::string& rType) {
    if (rType == "vertex")
        return LsdType::Vertex;
    if (rType == "line")
        return LsdType::Line;
    if (rType == "triangle")
        return LsdType::Triangle;
    if (rType == "quad")
        return LsdType::Quad;
    if (rType == "tetra")
        return LsdType::Tetra;
    if (rType == "pyramid")
        return LsdType::Pyramid;
    if (rType == "wedge")
        return LsdType::Wedge;
    if (rType == "tetra10")
        return LsdType::Tetra10;
    return LsdType::Hexahedron;
}

const char* lsd_family_keyword(LsdFamily Family) {
    switch (Family) {
        case LsdFamily::Solid:
            return "ELEMENT_SOLID";
        case LsdFamily::Shell:
            return "ELEMENT_SHELL";
        case LsdFamily::Beam:
            return "ELEMENT_BEAM";
        default:
            return "ELEMENT_MASS";
    }
}

struct LsdWritePart {
    std::int64_t mPid;
    std::string mName;
    const Region* mRegion;
};

// A cell region with a topological dimension is a part when none of its cells is
// already claimed; every other cell region is written as a set. Cells that no part
// claims form one part per cell block. Twin of `_assign_parts`.
void lsd_assign_parts(const std::vector<const Region*>& rRegions, const Mesh& rMesh,
                      const std::vector<std::size_t>& rBlockOfCell,
                      std::vector<std::int64_t>& rPidOf, std::vector<LsdWritePart>& rParts,
                      std::vector<const Region*>& rAsSets) {
    const std::size_t ncells = rBlockOfCell.size();
    rPidOf.assign(ncells, 0);
    std::set<std::int64_t> used;
    std::vector<LsdWritePart> parts;
    std::vector<char> claimed(ncells, 0);
    for (const Region* r : rRegions) {
        if (r->mKind != RegionKind::Cell || r->mDim < 0)
            continue;
        const std::int64_t* e = r->Entries();
        bool overlaps = false;
        for (std::size_t k = 0; k < r->NumEntries(); ++k)
            if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < ncells && claimed[e[k]])
                overlaps = true;
        if (overlaps) {
            rAsSets.push_back(r);
            continue;
        }
        for (std::size_t k = 0; k < r->NumEntries(); ++k)
            if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < ncells)
                claimed[e[k]] = 1;
        const std::int64_t want = (r->mTag > 0 && !used.count(r->mTag)) ? r->mTag : 0;
        parts.push_back({want, r->mName, r});
        if (want)
            used.insert(want);
    }
    std::int64_t nxt = used.empty() ? 0 : *used.rbegin();
    for (LsdWritePart& part : parts) {
        if (!part.mPid)
            part.mPid = ++nxt;
        const std::int64_t* e = part.mRegion->Entries();
        for (std::size_t k = 0; k < part.mRegion->NumEntries(); ++k)
            if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < ncells)
                rPidOf[e[k]] = part.mPid;
    }
    rParts = parts;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        bool any = false;
        const std::int64_t pid = nxt + 1;
        for (std::size_t g = 0; g < ncells; ++g) {
            if (rBlockOfCell[g] == b && !rPidOf[g]) {
                rPidOf[g] = pid;
                any = true;
            }
        }
        if (!any)
            continue;
        ++nxt;
        rParts.push_back({pid, std::string(rMesh.Cells(b).Type()), nullptr});
    }
}

void lsd_write_ids(std::string& rOut, const std::vector<std::int64_t>& rIds) {
    for (std::size_t k = 0; k < rIds.size(); k += 8) {
        for (std::size_t j = k; j < std::min(k + 8, rIds.size()); ++j)
            lsd_put_int(rOut, rIds[j], 10);
        rOut += '\n';
    }
}

// Corner nodes of one facet of one cell (a triangle repeats its last node), or
// nothing when the cell has no such facet.
std::optional<std::array<std::int64_t, 4>> lsd_segment_nodes(const Mesh& rMesh, std::int64_t G,
                                                             std::int64_t Facet) {
    std::int64_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t count = static_cast<std::int64_t>(cb.NumCells());
        if (G < base + count) {
            const std::string type(cb.Type());
            const NDArray& conn = cb.Conn();
            const std::size_t k = cb.NodesPerCell();
            const std::size_t r = static_cast<std::size_t>(G - base);
            std::vector<std::int64_t> nodes;
            const auto& faces = detail::cell_faces(cell_type_from_name(type));
            if (!faces.empty()) {
                if (Facet < 0 || static_cast<std::size_t>(Facet) >= faces.size())
                    return std::nullopt;
                const auto& face = faces[static_cast<std::size_t>(Facet)];
                for (std::size_t c = 0; c < face.mNumCorners; ++c)
                    nodes.push_back(detail::read_int(conn, r * k + face.mNodes[c]));
            } else if (type == "triangle" || type == "quad") {
                for (std::size_t c = 0; c < k; ++c)
                    nodes.push_back(detail::read_int(conn, r * k + c));
            } else {
                return std::nullopt;
            }
            if (nodes.size() == 3)
                nodes.push_back(nodes.back());
            return std::array<std::int64_t, 4>{nodes[0], nodes[1], nodes[2], nodes[3]};
        }
        base += count;
    }
    return std::nullopt;
}

}  // namespace

Mesh read_lsdyna(const std::string& rPath) {
    LsdDeck deck;
    lsd_read_file(deck, fs::path(rPath), 0, CardMode::Standard);
    return lsd_build_mesh(deck);
}

void write_lsdyna(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t npts = rMesh.NumPoints();
    if (static_cast<std::int64_t>(npts) > lsd_max_std_id)
        throw WriteError("LS-DYNA writer: too many nodes for 8-column ids");
    std::vector<std::size_t> block_of_cell;
    std::vector<LsdFamily> cell_family;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const auto family = lsd_type_family(type);
        if (!family)
            throw WriteError("LS-DYNA writer: unsupported cell type '" + type + "'");
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            block_of_cell.push_back(b);
            cell_family.push_back(*family);
        }
    }
    if (static_cast<std::int64_t>(block_of_cell.size()) > lsd_max_std_id)
        throw WriteError("LS-DYNA writer: too many elements for 8-column ids");

    std::vector<const Region*> regions;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        regions.push_back(&rMesh.Region(i));
    std::vector<std::int64_t> pid_of;
    std::vector<LsdWritePart> parts;
    std::vector<const Region*> as_sets;
    lsd_assign_parts(regions, rMesh, block_of_cell, pid_of, parts, as_sets);

    auto os = detail::make_classic_ofstream(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os << "*KEYWORD\n";
    os << detail::provenance_render_lines(detail::SlotTier::Block, "$ ");
    os << "*NODE\n";
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
        std::vector<std::string> rows(npts);
        parallel_for(npts, [&](std::size_t i) {
            std::string& row = rows[i];
            lsd_put_int(row, static_cast<std::int64_t>(i + 1), 8);
            for (std::size_t c = 0; c < 3; ++c)
                lsd_put(row,
                        c < dim ? detail::format_real16(detail::read_double(points, i * dim + c))
                                : std::string("0.0"),
                        16);
            row += '\n';
        });
        for (const std::string& row : rows)
            os << row;
    }
    {
        std::size_t g = 0;
        std::string out;
        for (const auto cb : rMesh.CellRange()) {
            const std::string type(cb.Type());
            const LsdFamily family = *lsd_type_family(type);
            const LsdType ltype = lsd_type_of(type);
            const bool two_line = ltype == LsdType::Tetra10;
            const NDArray& conn = cb.Conn();
            const std::size_t k = cb.NodesPerCell();
            out += std::string("*") + lsd_family_keyword(family) + "\n";
            for (std::size_t r = 0; r < cb.NumCells(); ++r, ++g) {
                std::vector<std::int64_t> row(k);
                for (std::size_t c = 0; c < k; ++c)
                    row[c] = detail::read_int(conn, r * k + c);
                const std::int64_t eid = static_cast<std::int64_t>(g + 1);
                const std::int64_t pid = pid_of[g];
                if (family == LsdFamily::Mass) {
                    lsd_put_int(out, eid, 8);
                    lsd_put_int(out, row[0] + 1, 8);
                    lsd_put(out, "0.0", 16);
                    lsd_put_int(out, pid, 8);
                    out += '\n';
                    continue;
                }
                std::vector<std::int64_t> nodes;
                if (family == LsdFamily::Solid) {
                    if (two_line) {
                        nodes = row;
                    } else {
                        const auto n = lsd_expand_solid(ltype, row.data());
                        nodes.assign(n.begin(), n.end());
                    }
                } else if (ltype == LsdType::Triangle) {
                    nodes = {row[0], row[1], row[2], row[2]};
                } else {
                    nodes = row;
                }
                lsd_put_int(out, eid, 8);
                lsd_put_int(out, pid, 8);
                if (two_line)
                    out += '\n';
                for (std::int64_t v : nodes)
                    lsd_put_int(out, v + 1, 8);
                out += '\n';
            }
            os << out;
            out.clear();
        }
    }
    for (const LsdWritePart& part : parts) {
        std::string out = "*PART\n" + lsd_title_line(part.mName) + "\n";
        lsd_put_int(out, part.mPid, 10);
        lsd_put_int(out, 1, 10);
        lsd_put_int(out, 1, 10);
        out += '\n';
        os << out;
    }

    // -- sets -----------------------------------------------------------------
    std::map<LsdFamily, std::set<std::int64_t>> used;
    auto sid_for = [&](LsdFamily family, std::int64_t tag) {
        std::set<std::int64_t>& taken = used[family];
        if (tag > 0 && !taken.count(tag)) {
            taken.insert(tag);
            return tag;
        }
        std::int64_t next = taken.empty() ? 1 : *taken.rbegin() + 1;
        while (taken.count(next))
            ++next;
        taken.insert(next);
        return next;
    };
    std::vector<const Region*> node_regions, side_regions, plain;
    for (const Region* r : regions) {
        if (r->mKind == RegionKind::Point)
            node_regions.push_back(r);
        else if (r->mKind == RegionKind::Side)
            side_regions.push_back(r);
        else if (r->mDim < 0)
            plain.push_back(r);
    }
    for (const Region* r : as_sets)
        if (r->mDim >= 0)
            plain.push_back(r);
    for (const Region* r : node_regions) {
        std::string out = "*SET_NODE_LIST_TITLE\n" + lsd_title_line(r->mName) + "\n";
        lsd_put_int(out, sid_for(LsdFamily::Node, r->mTag), 10);
        out += '\n';
        std::vector<std::int64_t> ids;
        for (std::size_t k = 0; k < r->NumEntries(); ++k)
            ids.push_back(r->Entries()[k] + 1);
        lsd_write_ids(out, ids);
        os << out;
    }
    for (const Region* r : plain) {
        std::map<LsdFamily, std::vector<std::int64_t>> by_family;
        for (std::size_t k = 0; k < r->NumEntries(); ++k) {
            const std::int64_t g = r->Entries()[k];
            if (g >= 0 && static_cast<std::size_t>(g) < cell_family.size())
                by_family[cell_family[static_cast<std::size_t>(g)]].push_back(g + 1);
        }
        if (by_family.empty())
            by_family.emplace(LsdFamily::Solid, std::vector<std::int64_t>{});
        if (by_family.size() > 1)
            log::warn("LS-DYNA writer: region '{}' spans several element families", r->mName);
        for (LsdFamily family : {LsdFamily::Solid, LsdFamily::Shell, LsdFamily::Beam}) {
            const auto it = by_family.find(family);
            if (it == by_family.end())
                continue;
            std::string out = "*SET_" + lsd_upper(lsd_family_name(family)) + "_LIST_TITLE\n" +
                              lsd_title_line(r->mName) + "\n";
            lsd_put_int(out, sid_for(family, r->mTag), 10);
            out += '\n';
            lsd_write_ids(out, it->second);
            os << out;
        }
    }
    for (const Region* r : side_regions) {
        std::string out = "*SET_SEGMENT_TITLE\n" + lsd_title_line(r->mName) + "\n";
        lsd_put_int(out, sid_for(LsdFamily::Segment, r->mTag), 10);
        out += '\n';
        for (std::size_t k = 0; k < r->NumEntries(); ++k) {
            const auto nodes =
                lsd_segment_nodes(rMesh, r->Entries()[2 * k], r->Entries()[2 * k + 1]);
            if (!nodes)
                continue;
            for (std::int64_t v : *nodes)
                lsd_put_int(out, v + 1, 10);
            out += '\n';
        }
        os << out;
    }
    os << "*END\n";
}

}  // namespace meshioplusplus
