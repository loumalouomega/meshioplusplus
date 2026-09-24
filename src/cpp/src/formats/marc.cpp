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
// MSC Marc input decks (.dat) and formatted post files (.t19). The layouts follow
// Marc Volume C (program input) and Volume D (PLDUMP2000, the post file); see
// doc/formats/marc.md. Python twin: src/python/meshioplusplus/marc/_marc.py.

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/marc.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kMarcDat = "Marc .dat";
constexpr const char* kMarcT19 = "Marc .t19";

[[noreturn]] void marc_fail(const char* pLabel, const std::string& rWhat) {
    throw ReadError(std::string(pLabel) + ": " + rWhat);
}

struct MarcType {
    const char* mCell;
    std::size_t mNodes;
};

// Marc element type -> (meshio++ cell type, the element's node count), from
// Volume B. Types whose topology the manuals were not checked for are left out:
// their elements are skipped with a warning. Every type lists its nodes in
// meshio++'s order (corners with the face 1-2-3(-4) normal pointing into the
// element, then mid-edge nodes bottom ring, top ring, verticals), so no
// permutation applies.
const MarcType* marc_type(std::int64_t Type) {
    static const std::unordered_map<std::int64_t, MarcType> kTypes = [] {
        std::unordered_map<std::int64_t, MarcType> m;
        // Plain types, and (v16.12.0) the Herrmann (mixed) types whose pressure
        // sits at the corners, the rebar and the composite types, all with the
        // node lists of their plain twins.
        for (int t : {3, 10, 11, 18, 75, 139, 140, 143, 144, 145, 147, 151, 152})
            m[t] = {"quad", 4};
        for (int t : {2, 6, 138, 158, 201})
            m[t] = {"triangle", 3};
        for (int t :
             {22, 26, 27, 28, 30, 32, 33, 46, 48, 53, 54, 55, 58, 59, 63, 66, 142, 148, 153, 154})
            m[t] = {"quad8", 8};
        for (int t : {124, 125, 126, 128, 129, 200})
            m[t] = {"triangle6", 6};
        for (int t : {7, 43, 117, 123, 146, 149})
            m[t] = {"hexahedron", 8};
        for (int t : {21, 23, 35, 44, 57, 61, 150})
            m[t] = {"hexahedron20", 20};
        for (int t : {134, 135})
            m[t] = {"tetra", 4};
        for (int t : {127, 130, 133})
            m[t] = {"tetra10", 10};
        for (int t : {136, 137})
            m[t] = {"wedge", 6};
        for (int t : {9, 31, 52, 98, 165, 166, 167})
            m[t] = {"line", 2};
        for (int t : {64, 168, 169, 170})
            m[t] = {"line3", 3};
        // Elements with more nodes than their cell keeps: the leading
        // geometric nodes are the cell, the rest (a Herrmann pressure node, a
        // centroid bubble node, generalized plane strain nodes) are dropped.
        for (int t : {80, 82, 83, 118, 119})
            m[t] = {"quad", 5};
        m[81] = {"quad", 7};
        for (int t : {34, 47, 60})
            m[t] = {"quad8", 10};
        for (int t : {84, 120})
            m[t] = {"hexahedron", 9};
        for (int t : {155, 156})
            m[t] = {"triangle", 4};
        m[157] = {"tetra", 5};
        return m;
    }();
    const auto it = kTypes.find(Type);
    return it == kTypes.end() ? nullptr : &it->second;
}

// An element's first line holds its number, type and 14 nodes; more continue.
constexpr std::size_t kMarcFirstLineNodes = 14;

std::string marc_lower(std::string_view Text) {
    std::string out(Text);
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

std::string marc_trim(std::string_view Text) {
    const auto blank = [](char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };
    std::size_t a = 0, b = Text.size();
    while (a < b && blank(Text[a]))
        ++a;
    while (b > a && blank(Text[b - 1]))
        --b;
    return std::string(Text.substr(a, b - a));
}

std::string marc_rstrip(std::string_view Text) {
    std::size_t b = Text.size();
    while (b > 0 && (Text[b - 1] == ' ' || Text[b - 1] == '\t' || Text[b - 1] == '\r' ||
                     Text[b - 1] == '\n'))
        --b;
    return std::string(Text.substr(0, b));
}

// Words split on blanks and commas.
std::vector<std::string> marc_words(std::string_view Text) {
    std::vector<std::string> out;
    std::string word;
    for (char c : Text) {
        if (c == ' ' || c == '\t' || c == ',' || c == '\r' || c == '\n') {
            if (!word.empty())
                out.push_back(std::move(word));
            word.clear();
        } else {
            word += c;
        }
    }
    if (!word.empty())
        out.push_back(std::move(word));
    return out;
}

bool marc_is_comment(std::string_view Line) {
    const std::string s = marc_trim(Line);
    return s.empty() || s[0] == '$';
}

// A data line starts (after blanks) with a number; a keyword with a letter.
bool marc_is_data(std::string_view Line) {
    const std::size_t k = Line.find_first_not_of(" \t");
    if (k == std::string_view::npos || Line[k] == '\r' || Line[k] == '\n')
        return false;
    const char c = Line[k];
    return std::isdigit(static_cast<unsigned char>(c)) || c == '+' || c == '-' || c == '.';
}

std::vector<std::string> marc_lines(const std::string& rPath, const char* pLabel) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        marc_fail(pLabel, "cannot open " + rPath);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<std::string> lines;
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::size_t eol = text.find('\n', pos);
        if (eol == std::string::npos)
            eol = text.size();
        std::string line = text.substr(pos, eol - pos);
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(std::move(line));
        pos = eol + 1;
    }
    return lines;
}

// An `INCLUDE` option line (the keyword at the start of the line, then the
// file name after a blank or comma): the file it names, else empty.
std::string marc_include_target(const std::string& rLine) {
    if (rLine.size() < 7 || marc_lower(rLine.substr(0, 7)) != "include")
        return "";
    std::size_t k = 7;
    if (k < rLine.size() && rLine[k] != ' ' && rLine[k] != '\t' && rLine[k] != ',')
        return "";
    while (k < rLine.size() && (rLine[k] == ' ' || rLine[k] == '\t' || rLine[k] == ','))
        ++k;
    std::string name = marc_trim(rLine.substr(k));
    if (name.size() >= 2 && (name.front() == '"' || name.front() == '\'') &&
        name.back() == name.front())
        name = name.substr(1, name.size() - 2);
    return name;
}

// The deck's lines with every `INCLUDE` replaced by the lines of the file it
// names (relative to the including file), recursively.
std::vector<std::string> marc_deck_lines(const std::string& rPath, int Depth = 0) {
    if (Depth > 16)
        marc_fail(kMarcDat, "INCLUDE files nest more than 16 deep (a cycle?) at " + rPath);
    std::vector<std::string> out;
    for (std::string& line : marc_lines(rPath, kMarcDat)) {
        const std::string target = marc_include_target(line);
        if (target.empty()) {
            out.push_back(std::move(line));
            continue;
        }
        std::filesystem::path file(target);
        if (file.is_relative())
            file = std::filesystem::path(rPath).parent_path() / file;
        std::error_code ec;
        if (!std::filesystem::is_regular_file(file, ec))
            marc_fail(kMarcDat, "INCLUDE names " + file.string() + ", which does not exist");
        std::vector<std::string> inner = marc_deck_lines(file.string(), Depth + 1);
        out.insert(out.end(), std::make_move_iterator(inner.begin()),
                   std::make_move_iterator(inner.end()));
    }
    return out;
}

std::int64_t marc_int(const std::string& rText, const char* pLabel, const std::string& rWhere) {
    return detail::card_to_int(marc_trim(rText), " (" + rWhere + ")", pLabel);
}

double marc_real(const std::string& rText, const char* pLabel, const std::string& rWhere) {
    return detail::card_to_real(marc_trim(rText), " (" + rWhere + ")", pLabel);
}

// The fields of a deck's data line: comma-separated (free format), or fixed
// columns of Width up to the last non-blank one.
std::vector<std::string> marc_fields(std::string_view Line, std::size_t Width) {
    std::vector<std::string> out;
    if (Line.find(',') != std::string_view::npos) {
        std::size_t pos = 0;
        for (;;) {
            const std::size_t comma = Line.find(',', pos);
            out.push_back(marc_trim(Line.substr(
                pos, comma == std::string_view::npos ? std::string_view::npos : comma - pos)));
            if (comma == std::string_view::npos)
                break;
            pos = comma + 1;
        }
        while (!out.empty() && out.back().empty())
            out.pop_back();  // a lone item is followed by a comma
        return out;
    }
    const std::string body = marc_rstrip(Line);
    for (std::size_t k = 0; k < body.size(); k += Width)
        out.push_back(marc_trim(std::string_view(body).substr(k, Width)));
    return out;
}

// The first word of a line that can open a deck (a Marc parameter).
bool marc_is_parameter(const std::string& rWord) {
    static const std::unordered_set<std::string> kParameters = {
        "title",      "sizing",  "elements", "extended",   "version",       "table",
        "processor",  "alloc",   "setname",  "dist",       "large",         "update",
        "finite",     "all",     "no",       "state",      "heat",          "coupled",
        "harmonic",   "buckle",  "dynamic",  "fluid",      "electrostatic", "magnetostatic",
        "joule",      "bearing", "shell",    "print",      "lumping",       "assumed",
        "constant",   "feature", "follow",   "plasticity", "rezoning",      "scale",
        "structural", "thermal", "end"};
    return kParameters.count(rWord) != 0;
}

// -- the input deck -------------------------------------------------------------

struct MarcElement {
    std::int64_t mId = 0, mType = 0;
    std::vector<std::int64_t> mNodes;
};

struct MarcSet {
    std::string mName;
    bool mElements = true;  // else a node set
    std::vector<std::string> mTokens;
    std::vector<std::int64_t> mMembers;
    // Edge (1) and face (2) sets: mMembers are elements, mSides their Marc
    // edge or face numbers (Volume A's numbering, not mapped to facets).
    int mSideKind = 0;
    std::vector<std::int64_t> mSides;
};

struct MarcDeck {
    bool mExtended = false;
    std::int64_t mNcoord = 3;
    std::unordered_map<std::int64_t, std::array<double, 3>> mNodes;  // the last definition wins
    std::vector<std::int64_t> mNodeOrder;
    std::vector<MarcElement> mElements;
    std::vector<MarcSet> mSets;

    std::size_t IntWidth() const { return mExtended ? 10 : 5; }
    std::size_t RealWidth() const { return mExtended ? 20 : 10; }
};

std::size_t marc_next_data(const std::vector<std::string>& rLines, std::size_t I) {
    while (I < rLines.size() && marc_is_comment(rLines[I]))
        ++I;
    return I;
}

std::string marc_where(std::size_t I) {
    return "line " + std::to_string(I + 1);
}

std::size_t marc_connectivity(MarcDeck& rDeck, const std::vector<std::string>& rLines,
                              std::size_t I) {
    I = marc_next_data(rLines, I);
    if (I < rLines.size() && marc_is_data(rLines[I]))
        ++I;  // the header line: element count, unit, print flag ...
    const std::size_t width = rDeck.IntWidth();
    for (;;) {
        I = marc_next_data(rLines, I);
        if (I >= rLines.size() || !marc_is_data(rLines[I]))
            return I;
        std::string where = marc_where(I);
        const auto f = marc_fields(rLines[I], width);
        if (f.size() < 2)
            marc_fail(kMarcDat, "an element needs a number and a type (" + where + ")");
        MarcElement el;
        el.mId = marc_int(f[0], kMarcDat, where);
        el.mType = marc_int(f[1], kMarcDat, where);
        for (std::size_t k = 2; k < f.size(); ++k)
            el.mNodes.push_back(marc_int(f[k], kMarcDat, where));
        ++I;
        const MarcType* known = marc_type(el.mType);
        if (!known) {
            if (f.size() >= 2 + kMarcFirstLineNodes)
                marc_fail(kMarcDat, "element " + std::to_string(el.mId) + " is of type " +
                                        std::to_string(el.mType) +
                                        ", which meshio++ does not know the node count of (its "
                                        "nodes continue on further lines)");
            rDeck.mElements.push_back(std::move(el));
            continue;
        }
        while (el.mNodes.size() < known->mNodes) {
            I = marc_next_data(rLines, I);
            if (I >= rLines.size() || !marc_is_data(rLines[I]))
                marc_fail(kMarcDat, "element " + std::to_string(el.mId) + " lists " +
                                        std::to_string(el.mNodes.size()) + " of its " +
                                        std::to_string(known->mNodes) + " nodes");
            where = marc_where(I);
            for (const std::string& v : marc_fields(rLines[I], width))
                el.mNodes.push_back(marc_int(v, kMarcDat, where));
            ++I;
        }
        el.mNodes.resize(known->mNodes);
        rDeck.mElements.push_back(std::move(el));
    }
}

std::vector<std::string> marc_slices(std::string_view Body, std::size_t Width) {
    std::vector<std::string> out;
    const std::string body = marc_rstrip(Body);
    for (std::size_t k = 0; k < body.size(); k += Width)
        out.push_back(body.substr(k, Width));
    return out;
}

std::size_t marc_coordinates(MarcDeck& rDeck, const std::vector<std::string>& rLines,
                             std::size_t I) {
    I = marc_next_data(rLines, I);
    if (I < rLines.size() && marc_is_data(rLines[I])) {
        const auto header = marc_fields(rLines[I], rDeck.IntWidth());
        if (!header.empty() && !header[0].empty())
            rDeck.mNcoord = std::max<std::int64_t>(1, marc_int(header[0], kMarcDat, marc_where(I)));
        ++I;
    }
    const auto ncoord = static_cast<std::size_t>(rDeck.mNcoord);
    const std::size_t iw = rDeck.IntWidth(), rw = rDeck.RealWidth();
    for (;;) {
        I = marc_next_data(rLines, I);
        if (I >= rLines.size() || !marc_is_data(rLines[I]))
            return I;
        const std::string where = marc_where(I);
        const std::string& line = rLines[I];
        std::int64_t ident = 0;
        std::vector<std::string> values;
        if (line.find(',') != std::string::npos) {
            const auto f = marc_fields(line, iw);
            ident = marc_int(f[0], kMarcDat, where);
            values.assign(f.begin() + 1, f.end());
        } else {
            ident = marc_int(line.substr(0, std::min(iw, line.size())), kMarcDat, where);
            if (line.size() > iw)
                values = marc_slices(std::string_view(line).substr(iw), rw);
        }
        ++I;
        while (values.size() < ncoord) {  // continuation lines, six reals each
            I = marc_next_data(rLines, I);
            if (I >= rLines.size() || !marc_is_data(rLines[I]))
                break;
            const std::string& more = rLines[I];
            const auto extra =
                more.find(',') != std::string::npos ? marc_fields(more, rw) : marc_slices(more, rw);
            values.insert(values.end(), extra.begin(), extra.end());
            ++I;
        }
        std::array<double, 3> xyz{0.0, 0.0, 0.0};
        for (std::size_t d = 0; d < ncoord && d < values.size(); ++d) {
            const double v = marc_real(values[d], kMarcDat, where);
            if (d < 3)
                xyz[d] = v;
        }
        if (!rDeck.mNodes.count(ident))
            rDeck.mNodeOrder.push_back(ident);
        rDeck.mNodes[ident] = xyz;
    }
}

bool marc_all_digits(const std::string& rText) {
    return !rText.empty() && std::all_of(rText.begin(), rText.end(), [](char c) {
        return std::isdigit(static_cast<unsigned char>(c));
    });
}

// The items of a set's data line: words and numbers (touching fixed-width
// integers split apart); `c` or `continue` last means more lines follow.
std::vector<std::string> marc_set_tokens(std::string_view Text, std::size_t Width) {
    std::vector<std::string> out;
    for (const std::string& tok : marc_words(Text)) {
        if (marc_all_digits(tok) && tok.size() > Width && tok.size() % Width == 0) {
            for (std::size_t k = 0; k < tok.size(); k += Width)
                out.push_back(tok.substr(k, Width));
        } else {
            out.push_back(marc_lower(tok));
        }
    }
    return out;
}

bool marc_is_integer(const std::string& rText) {
    std::size_t k = 0;
    while (k < rText.size() && (rText[k] == '+' || rText[k] == '-'))
        ++k;
    return k < rText.size() && marc_all_digits(rText.substr(k));
}

std::size_t marc_define(MarcDeck& rDeck, const std::vector<std::string>& rLines, std::size_t I) {
    const std::string where = marc_where(I);
    const auto raw = marc_words(rLines[I]);
    const std::string kind = raw.size() > 1 ? marc_lower(raw[1]) : "";
    std::size_t k = 2;
    if (raw.size() > k && (marc_lower(raw[k]) == "set" || marc_lower(raw[k]) == "oset"))
        ++k;
    const std::string name = raw.size() > k ? raw[k] : "";
    std::unordered_set<std::string> earlier;
    for (const MarcSet& s : rDeck.mSets)
        earlier.insert(marc_lower(s.mName));
    ++I;
    std::vector<std::string> tokens;
    for (;;) {
        I = marc_next_data(rLines, I);
        if (I >= rLines.size())
            break;
        const std::string& line = rLines[I];
        const auto items = marc_set_tokens(line, rDeck.IntWidth());
        // The first data line starts with a number or an earlier set's name;
        // later ones follow a line that ended in C (continue).
        if (items.empty() || (tokens.empty() && !marc_is_data(line) && !earlier.count(items[0])))
            break;
        tokens.insert(tokens.end(), items.begin(), items.end());
        ++I;
        if (tokens.back() != "c" && tokens.back() != "continue")
            break;
        tokens.pop_back();
    }
    if (kind == "edge" || kind == "face") {
        // `elem:number` members.
        MarcSet s;
        s.mName = name;
        s.mSideKind = kind == "edge" ? 1 : 2;
        for (const std::string& t : tokens) {
            const std::size_t colon = t.find(':');
            if (colon == std::string::npos || !marc_is_integer(t.substr(0, colon)) ||
                !marc_is_integer(t.substr(colon + 1))) {
                log::warn("{}: DEFINE {} SET '{}': '{}' is not an element:number pair ({})",
                          kMarcDat, kind == "edge" ? "EDGE" : "FACE", name, t, where);
                continue;
            }
            s.mMembers.push_back(marc_int(t.substr(0, colon), kMarcDat, where));
            s.mSides.push_back(marc_int(t.substr(colon + 1), kMarcDat, where));
        }
        rDeck.mSets.push_back(std::move(s));
    } else if (kind == "element" || kind == "elsq" || kind == "node" || kind == "ndsq") {
        MarcSet s;
        s.mName = name;
        s.mElements = kind == "element" || kind == "elsq";
        s.mTokens = std::move(tokens);
        rDeck.mSets.push_back(std::move(s));
    } else {
        log::warn(
            "{}: DEFINE {} SET '{}' is not read ({})", kMarcDat,
            [&] {
                std::string up = kind;
                for (char& c : up)
                    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                return up;
            }(),
            name, where);
    }
    return I;
}

MarcDeck marc_parse_deck(const std::vector<std::string>& rLines) {
    MarcDeck deck;
    std::size_t i = 0;
    bool in_parameters = true;
    while (i < rLines.size()) {
        const std::string& line = rLines[i];
        if (marc_is_comment(line) || marc_is_data(line) || line.empty() || line[0] == ' ' ||
            line[0] == '\t') {
            ++i;
            continue;
        }
        const auto words = marc_words(marc_lower(line));
        const std::string key = words.empty() ? "" : words[0];
        const bool end_option = key == "end" && words.size() > 1 && words[1] == "option";
        if (in_parameters) {
            if (key == "extended") {
                deck.mExtended = true;
            } else if (key == "end" && !end_option) {
                in_parameters = false;
                ++i;
                continue;
            }
            if (key != "connectivity" && key != "coordinates" && key != "define") {
                ++i;
                continue;
            }
            in_parameters = false;  // a deck without END: the model starts here
        }
        if (end_option)
            break;
        if (key == "connectivity")
            i = marc_connectivity(deck, rLines, i + 1);
        else if (key == "coordinates")
            i = marc_coordinates(deck, rLines, i + 1);
        else if (key == "define")
            i = marc_define(deck, rLines, i);
        else
            ++i;
    }
    return deck;
}

// The members of a set: numbers, `a TO b [BY c]` ranges and other sets' names,
// combined left to right by AND (the default), EXCEPT and INTERSECT.
std::vector<std::int64_t> marc_expand(
    const std::vector<std::string>& rTokens, const std::string& rName,
    const std::unordered_map<std::string, std::vector<std::int64_t>>& rKnown, const char* pLabel) {
    std::vector<std::int64_t> out;
    std::string op = "and";
    const auto combine = [&](const std::vector<std::int64_t>& rItems) {
        if (op == "and") {
            std::unordered_set<std::int64_t> have(out.begin(), out.end());
            for (std::int64_t v : rItems)
                if (have.insert(v).second)
                    out.push_back(v);
        } else if (op == "except") {
            const std::unordered_set<std::int64_t> drop(rItems.begin(), rItems.end());
            out.erase(std::remove_if(out.begin(), out.end(),
                                     [&](std::int64_t v) { return drop.count(v) != 0; }),
                      out.end());
        } else {
            const std::unordered_set<std::int64_t> keep(rItems.begin(), rItems.end());
            out.erase(std::remove_if(out.begin(), out.end(),
                                     [&](std::int64_t v) { return keep.count(v) == 0; }),
                      out.end());
        }
    };
    const auto number = [&](const std::string& rText) {
        return detail::card_to_int(rText, " (set '" + rName + "')", pLabel);
    };
    std::size_t k = 0;
    while (k < rTokens.size()) {
        const std::string& tok = rTokens[k];
        if (tok == "and" || tok == "except" || tok == "intersect") {
            op = tok;
            ++k;
            continue;
        }
        if (marc_is_integer(tok)) {
            const std::int64_t start = number(tok);
            ++k;
            std::vector<std::int64_t> items;
            if (k < rTokens.size() && (rTokens[k] == "to" || rTokens[k] == "through")) {
                if (k + 1 >= rTokens.size())
                    marc_fail(pLabel, "set '" + rName + "' ends in a range with no end");
                const std::int64_t stop = number(rTokens[k + 1]);
                k += 2;
                std::int64_t step = 1;
                if (k < rTokens.size() && rTokens[k] == "by") {
                    if (k + 1 >= rTokens.size())
                        marc_fail(pLabel, "set '" + rName + "' ends with BY and no step");
                    step = std::abs(number(rTokens[k + 1]));
                    if (step == 0)
                        step = 1;
                    k += 2;
                }
                if (stop >= start)
                    for (std::int64_t v = start; v <= stop; v += step)
                        items.push_back(v);
                else
                    for (std::int64_t v = start; v >= stop; v -= step)
                        items.push_back(v);
            } else {
                items.push_back(start);
            }
            combine(items);
            op = "and";
            continue;
        }
        if (const auto it = rKnown.find(tok); it != rKnown.end()) {
            combine(it->second);
            op = "and";
            ++k;
            continue;
        }
        marc_fail(pLabel, "set '" + rName + "' names '" + tok + "', which is not an earlier set");
    }
    return out;
}

struct MarcLocation {
    std::int64_t mBlock = -1, mRow = -1;
};

// The mesh of nodes, elements and sets (members expanded); each element's
// (block, row) in rLocations (-1 when it has no cell).
Mesh marc_build(const char* pLabel, const std::vector<std::int64_t>& rNodeIds,
                const std::vector<std::array<double, 3>>& rCoords,
                const std::vector<MarcElement>& rElements, const std::vector<MarcSet>& rSets,
                std::vector<MarcLocation>& rLocations) {
    std::unordered_map<std::int64_t, std::int64_t> index;
    for (std::size_t k = 0; k < rNodeIds.size(); ++k)
        index.emplace(rNodeIds[k], static_cast<std::int64_t>(k));
    struct Block {
        std::string mType;
        std::vector<std::int64_t> mConn, mId, mMarcType;
    };
    std::vector<Block> blocks;
    std::map<std::string, std::size_t> block_of;
    std::map<std::int64_t, std::size_t> skipped;
    std::size_t degenerate20 = 0;
    std::unordered_map<std::int64_t, MarcLocation> element_cell;
    std::unordered_set<std::int64_t> all_elements;
    rLocations.assign(rElements.size(), MarcLocation{});
    for (std::size_t e = 0; e < rElements.size(); ++e) {
        const MarcElement& el = rElements[e];
        all_elements.insert(el.mId);
        const MarcType* known = marc_type(el.mType);
        if (!known) {
            ++skipped[el.mType];
            continue;
        }
        std::string cell = known->mCell;
        std::vector<std::int64_t> nodes = el.mNodes;
        nodes.resize(std::min(nodes.size(), static_cast<std::size_t>(
                                                cell_type_num_nodes(cell_type_from_name(cell)))));
        // The 3-node rebar lines list their middle node second.
        if (el.mType >= 168 && el.mType <= 170 && nodes.size() == 3)
            std::swap(nodes[1], nodes[2]);
        for (std::int64_t v : nodes)
            if (!index.count(v))
                marc_fail(pLabel, "element " + std::to_string(el.mId) + " names undefined node " +
                                      std::to_string(v));
        if (cell == "hexahedron") {
            std::array<std::int64_t, 8> a{};
            std::copy(nodes.begin(), nodes.begin() + 8, a.begin());
            detail::CollapsedBrick c = detail::collapse_brick(a);
            cell = c.mType;
            nodes = std::move(c.mNodes);
        } else if (cell == "hexahedron20" &&
                   std::set<std::int64_t>(nodes.begin(), nodes.begin() + 8).size() < 8) {
            ++degenerate20;
        }
        const auto [it, fresh] = block_of.emplace(cell, blocks.size());
        if (fresh)
            blocks.push_back(Block{cell, {}, {}, {}});
        Block& b = blocks[it->second];
        const MarcLocation loc{static_cast<std::int64_t>(it->second),
                               static_cast<std::int64_t>(b.mId.size())};
        rLocations[e] = loc;
        element_cell.emplace(el.mId, loc);
        for (std::int64_t v : nodes)
            b.mConn.push_back(index.at(v));
        b.mId.push_back(el.mId);
        b.mMarcType.push_back(el.mType);
    }
    if (!skipped.empty()) {
        std::string listed;
        for (const auto& [type, count] : skipped)
            listed += (listed.empty() ? "" : ", ") + std::to_string(type) + " (" +
                      std::to_string(count) + ")";
        log::warn("{}: elements of types meshio++ has no cell for were skipped: {}", pLabel,
                  listed);
    }
    if (degenerate20)
        log::warn("{}: {} 20-node brick(s) with repeated corner nodes kept as hexahedron20", pLabel,
                  degenerate20);

    Mesh mesh;
    NDArray points(DType::Float64, {rCoords.size(), 3});
    double* p = points.As<double>();
    for (std::size_t k = 0; k < rCoords.size(); ++k)
        std::copy(rCoords[k].begin(), rCoords[k].end(), p + k * 3);
    mesh.AssignPoints(std::move(points));
    const auto column = [](const std::vector<std::int64_t>& rValues) {
        NDArray a(DType::Int64, {rValues.size()});
        std::copy(rValues.begin(), rValues.end(), a.As<std::int64_t>());
        return a;
    };
    std::vector<NDArray> ids, types;
    std::vector<std::int64_t> starts;
    std::vector<int> dims;
    std::int64_t start = 0;
    for (Block& b : blocks) {
        const std::size_t rows = b.mId.size();
        NDArray conn(DType::Int64, {rows, b.mConn.size() / std::max<std::size_t>(rows, 1)});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
        ids.push_back(column(b.mId));
        types.push_back(column(b.mMarcType));
        starts.push_back(start);
        start += static_cast<std::int64_t>(rows);
        dims.push_back(cell_type_dimension(cell_type_from_name(b.mType)));
    }
    if (!blocks.empty()) {
        mesh.AddCellData("marc:element", std::move(ids));
        mesh.AddCellData("marc:type", std::move(types));
    }
    // One region per (kind, name): a later set replaces an earlier one.
    std::map<std::pair<int, std::string>, meshioplusplus::Region> regions;
    for (const MarcSet& s : rSets) {
        if (s.mSideKind) {
            // (cell or -1, Marc edge/face number): Marc numbers an element's
            // edges and faces its own way (Volume A), not mapped to facets.
            NDArray a(DType::Int64, {s.mMembers.size(), 2});
            for (std::size_t k = 0; k < s.mMembers.size(); ++k) {
                const auto it = element_cell.find(s.mMembers[k]);
                a.As<std::int64_t>()[2 * k] =
                    it == element_cell.end()
                        ? -1
                        : starts[static_cast<std::size_t>(it->second.mBlock)] + it->second.mRow;
                a.As<std::int64_t>()[2 * k + 1] = k < s.mSides.size() ? s.mSides[k] : -1;
            }
            mesh.AddFieldData(
                std::string(s.mSideKind == 1 ? "marc:edge_set:" : "marc:face_set:") + s.mName,
                std::move(a));
            continue;
        }
        std::vector<std::int64_t> entries;
        int dim = -1;
        if (s.mElements) {
            for (std::int64_t m : s.mMembers) {
                const auto it = element_cell.find(m);
                if (it == element_cell.end()) {
                    if (!all_elements.count(m))
                        marc_fail(pLabel, "element set '" + s.mName + "' names undefined element " +
                                              std::to_string(m));
                    continue;
                }
                const auto b = static_cast<std::size_t>(it->second.mBlock);
                entries.push_back(starts[b] + it->second.mRow);
                dim = std::max(dim, dims[b]);
            }
        } else {
            for (std::int64_t m : s.mMembers) {
                const auto it = index.find(m);
                if (it == index.end())
                    marc_fail(pLabel, "node set '" + s.mName + "' names undefined node " +
                                          std::to_string(m));
                entries.push_back(it->second);
            }
        }
        const RegionKind kind = s.mElements ? RegionKind::Cell : RegionKind::Point;
        const auto key = std::make_pair(static_cast<int>(kind), s.mName);
        if (regions.count(key))
            log::warn("{}: a second {} set '{}' replaces the first", pLabel,
                      s.mElements ? "element" : "node", s.mName);
        regions.insert_or_assign(key,
                                 meshioplusplus::Region(s.mName, kind, dim, -1, column(entries)));
    }
    for (auto& [key, region] : regions)
        mesh.AddRegion(std::move(region));
    return mesh;
}

// -- the formatted post file ------------------------------------------------------

constexpr std::size_t kMarcW = 13;  // the post file's column width (i13, e13.6)

std::vector<std::int64_t> marc_ints_of(const std::string& rLine, const std::string& rWhere) {
    std::vector<std::int64_t> out;
    const std::string body = marc_rstrip(rLine);
    for (std::size_t k = 0; k < body.size(); k += kMarcW) {
        const std::string field = marc_trim(std::string_view(body).substr(k, kMarcW));
        if (!field.empty())
            out.push_back(marc_int(field, kMarcT19, rWhere));
    }
    return out;
}

std::vector<double> marc_reals_of(std::string_view Line, const std::string& rWhere) {
    std::vector<double> out;
    const std::string body = marc_rstrip(Line);
    for (std::size_t k = 0; k < body.size(); k += kMarcW) {
        const std::string field = marc_trim(std::string_view(body).substr(k, kMarcW));
        if (!field.empty())
            out.push_back(marc_real(field, kMarcT19, rWhere));
    }
    return out;
}

// A block: its number and its body's lines [mBegin, mEnd); increment markers
// (****, ----, ++++) carry number 0 and their own kind.
struct MarcBlock {
    char mKind = 'b';  // 'b' block, '*' ****, '-' ----, '+' ++++
    std::int64_t mNumber = 0;
    std::size_t mBegin = 0, mEnd = 0;
};

// Reads a block's lines as Fortran records: each record starts a line.
class MarcRecords {
public:
    MarcRecords(const std::vector<std::string>& rLines, const MarcBlock& rBlock)
        : mrLines(rLines), mI(rBlock.mBegin), mEnd(rBlock.mEnd) {}

    const std::string& Line() {
        if (mI >= mEnd)
            marc_fail(kMarcT19, "a block ends early (line " + std::to_string(mI + 1) + ")");
        return mrLines[mI++];
    }
    std::string Where() const { return "line " + std::to_string(mI + 1); }
    std::vector<std::int64_t> Ints(std::size_t Count) {
        std::vector<std::int64_t> out;
        while (out.size() < Count) {
            const std::string where = Where();
            const auto more = marc_ints_of(Line(), where);
            out.insert(out.end(), more.begin(), more.end());
        }
        out.resize(Count);
        return out;
    }
    std::vector<double> Reals(std::size_t Count) {
        std::vector<double> out;
        while (out.size() < Count) {
            const std::string where = Where();
            const auto more = marc_reals_of(Line(), where);
            out.insert(out.end(), more.begin(), more.end());
        }
        out.resize(Count);
        return out;
    }

private:
    const std::vector<std::string>& mrLines;
    std::size_t mI, mEnd;
};

struct MarcIncrement {
    double mTime = 0.0;
    std::int64_t mInc = 0, mIncsub = 0, mJantyp = 0;
    bool mNewModel = false;  // the increment remeshes (newmo)
};

class MarcPost {
public:
    explicit MarcPost(const std::string& rPath) : mLines(marc_lines(rPath, kMarcT19)) {
        if (mLines.empty() || mLines[0].rfind("=beg=501", 0) != 0)
            marc_fail(kMarcT19, "not a Marc formatted post file (no =beg=501 title block)");
        std::vector<MarcBlock> blocks;
        std::size_t i = 0;
        while (i < mLines.size()) {
            const std::string s = marc_trim(mLines[i]);
            if (s.rfind("=beg=", 0) == 0) {
                MarcBlock b;
                b.mNumber = marc_int(s.substr(5, 5), kMarcT19,
                                     "block header, line " + std::to_string(i + 1));
                std::size_t j = i + 1;
                while (j < mLines.size() && mLines[j].rfind("=end=", 0) != 0)
                    ++j;
                if (j >= mLines.size())
                    marc_fail(kMarcT19, "block " + std::to_string(b.mNumber) +
                                            " has no =end= (line " + std::to_string(i + 1) + ")");
                b.mBegin = i + 1;
                b.mEnd = j;
                blocks.push_back(b);
                i = j + 1;
                continue;
            }
            if (s == "****" || s == "----" || s == "++++")
                blocks.push_back(MarcBlock{s[0], 0, i, i});
            ++i;
        }
        // The model blocks come before the first increment.
        bool inside = false;
        for (const MarcBlock& b : blocks) {
            if (b.mKind == '*') {
                mIncrements.emplace_back();
                inside = true;
            } else if (b.mKind == '-' || b.mKind == '+') {
                inside = false;
            } else if (!inside) {
                mModel.push_back(b);
            } else {
                mIncrements.back().push_back(b);
            }
        }
    }

    MarcRecords Reader(const MarcBlock& rBlock) const { return MarcRecords(mLines, rBlock); }

    // The model header (block 502) and element post codes (506) of `pBlocks`,
    // else of the model before the first increment.
    void Header(std::vector<std::int64_t>& rLm,
                std::vector<std::pair<std::int64_t, std::string>>& rCodes,
                const std::vector<MarcBlock>* pBlocks = nullptr) const {
        rLm.assign(30, 0);
        rCodes.clear();
        for (const MarcBlock& b : pBlocks ? *pBlocks : mModel) {
            const std::int64_t family = b.mNumber / 100;
            if (family == 502) {
                const auto values = Reader(b).Ints(30);
                std::copy(values.begin(), values.end(), rLm.begin());
            } else if (family == 506) {
                MarcRecords r = Reader(b);
                for (std::int64_t k = 0; k < rLm[0]; ++k) {
                    const std::string where = r.Where();
                    const std::string& line = r.Line();
                    const std::int64_t code =
                        marc_int(line.substr(0, std::min(kMarcW, line.size())), kMarcT19, where);
                    const std::string label =
                        line.size() > kMarcW ? marc_trim(std::string_view(line).substr(kMarcW, 24))
                                             : "";
                    rCodes.emplace_back(code, label);
                }
            }
        }
    }

    void Model(const std::vector<std::int64_t>& rLm, std::vector<std::int64_t>& rNodeIds,
               std::vector<std::array<double, 3>>& rCoords, std::vector<MarcElement>& rElements,
               std::vector<MarcSet>& rSets, const std::vector<MarcBlock>* pBlocks = nullptr) const {
        const std::int64_t numnp = rLm[1], numel = rLm[2], ncrd = rLm[8], nnodmx = rLm[9];
        const std::int64_t postrv = rLm[13];
        for (const MarcBlock& b : pBlocks ? *pBlocks : mModel) {
            const std::int64_t family = b.mNumber / 100;
            MarcRecords r = Reader(b);
            if (family == 507) {
                for (std::int64_t e = 0; e < numel; ++e) {
                    const auto record = r.Ints(static_cast<std::size_t>(3 + nnodmx));
                    MarcElement el;
                    el.mId = record[0];
                    el.mType = record[1];
                    const auto nnod =
                        static_cast<std::size_t>(std::max<std::int64_t>(record[2], 0));
                    for (std::size_t k = 0; k < nnod && 3 + k < record.size(); ++k)
                        el.mNodes.push_back(record[3 + k]);
                    rElements.push_back(std::move(el));
                }
            } else if (family == 508) {
                for (std::int64_t k = 0; k < numnp; ++k) {
                    const std::string where = r.Where();
                    const std::string& line = r.Line();
                    const std::int64_t ident =
                        marc_int(line.substr(0, std::min(kMarcW, line.size())), kMarcT19, where);
                    std::vector<double> values =
                        line.size() > kMarcW
                            ? marc_reals_of(std::string_view(line).substr(kMarcW), where)
                            : std::vector<double>{};
                    if (values.size() > 5)
                        values.resize(5);
                    if (ncrd > 5) {
                        const auto more = r.Reals(static_cast<std::size_t>(ncrd - 5));
                        values.insert(values.end(), more.begin(), more.end());
                    }
                    std::array<double, 3> xyz{0.0, 0.0, 0.0};
                    for (std::size_t d = 0;
                         d < 3 && d < static_cast<std::size_t>(ncrd) && d < values.size(); ++d)
                        xyz[d] = values[d];
                    rNodeIds.push_back(ident);
                    rCoords.push_back(xyz);
                }
            } else if (family == 513) {
                std::int64_t count = 0;
                std::size_t width = 12;
                if (b.mNumber == 51301 || postrv > 10) {
                    count = r.Ints(1)[0];
                    width = 32;
                } else {
                    count = rLm[15];
                }
                for (std::int64_t s = 0; s < count; ++s) {
                    const std::string& line = r.Line();
                    const std::string name =
                        marc_trim(std::string_view(line).substr(0, std::min(width, line.size())));
                    const auto head = r.Ints(2);
                    const std::int64_t isetn = head[0], isett = head[1];
                    std::vector<std::int64_t> members, sides;
                    if (isetn > 0)
                        members = r.Ints(static_cast<std::size_t>(isetn));
                    const bool side_set = isett == 12 || isett == 13 || isett == 18 || isett == 19;
                    if (side_set && isetn > 0)
                        sides = r.Ints(static_cast<std::size_t>(isetn));  // edge or face numbers
                    if (isett == 0 || isett == 1 || side_set) {
                        MarcSet set;
                        set.mName = name;
                        set.mElements = isett != 1;
                        set.mSideKind = side_set ? (isett == 12 ? 1 : 2) : 0;
                        set.mMembers = std::move(members);
                        set.mSides = std::move(sides);
                        rSets.push_back(std::move(set));
                    } else {
                        log::warn("{}: set '{}' of type {} is not read", kMarcT19, name, isett);
                    }
                }
            }
        }
    }

    MarcIncrement Info(const std::vector<MarcBlock>& rBlocks) const {
        std::vector<std::int64_t> lm(12, 0);
        std::vector<double> xlm(6, 0.0);
        for (const MarcBlock& b : rBlocks) {
            const std::int64_t family = b.mNumber / 100;
            if (family == 517) {
                lm = Reader(b).Ints(12);
            } else if (b.mNumber == 51800) {
                xlm = Reader(b).Reals(6);
            } else if (b.mNumber == 51801) {
                MarcRecords r = Reader(b);
                const std::int64_t nw = r.Ints(1)[0];
                xlm = r.Reals(static_cast<std::size_t>(std::max<std::int64_t>(nw, 0)));
                xlm.resize(xlm.size() + 6, 0.0);
            }
        }
        const std::int64_t ihresp = lm[6];
        MarcIncrement out;
        out.mTime = ihresp >= 1 && ihresp <= 4 ? xlm[1] : xlm[0];
        out.mInc = lm[1];
        out.mIncsub = lm[2];
        out.mJantyp = lm[3];
        out.mNewModel = lm[0] != 0;
        return out;
    }

    std::vector<double> Times() const {
        std::vector<double> out;
        for (const auto& inc : mIncrements)
            out.push_back(Info(inc).mTime);
        return out;
    }

    const std::vector<std::vector<MarcBlock>>& Increments() const { return mIncrements; }

private:
    std::vector<std::string> mLines;
    std::vector<MarcBlock> mModel;
    std::vector<std::vector<MarcBlock>> mIncrements;
};

// Element post codes (Volume C, Table 3-3) that open a symmetric tensor written
// as six codes, components 11 22 33 12 23 31 (xx yy zz xy yz zx).
const char* marc_tensor_name(std::int64_t Code) {
    static const std::map<std::int64_t, const char*> kNames = {
        {301, "Total Strain"},
        {311, "Stress"},
        {321, "Plastic Strain"},
        {331, "Creep Strain"},
        {341, "Cauchy Stress"},
        {351, "Real Harmonic Stress"},
        {361, "Imaginary Harmonic Stress"},
        {371, "Thermal Strain"},
        {381, "Cracking Strain"},
        {391, "Stress in Preferred System"},
        {401, "Elastic Strain"},
        {411, "Global Stress"},
        {421, "Global Elastic Strain"},
        {431, "Global Plastic Strain"},
        {441, "Global Creep Strain"},
        {461, "Elastic Strain in Preferred System"},
        {541, "Phase Transformation Strain"},
    };
    const auto it = kNames.find(Code);
    return it == kNames.end() ? nullptr : it->second;
}

std::string marc_scalar_name(std::int64_t Code) {
    static const std::map<std::int64_t, const char*> kNames = {
        {7, "Equivalent Plastic Strain"},
        {8, "Equivalent Creep Strain"},
        {9, "Temperature"},
        {17, "Equivalent Von Mises Stress"},
        {18, "Mean Normal Stress"},
        {20, "Thickness"},
        {47, "Equivalent Cauchy Stress"},
        {48, "Strain Energy Density"},
        {127, "Equivalent Elastic Strain"},
    };
    const auto it = kNames.find(Code);
    return it == kNames.end() ? "post code " + std::to_string(Code) : it->second;
}

// An element result: its name, its first post code column and its width (1 or 6).
struct MarcColumn {
    std::string mName;
    std::size_t mFirst = 0, mWidth = 1;
};

std::vector<MarcColumn> marc_columns(
    const std::vector<std::pair<std::int64_t, std::string>>& rCodes) {
    std::vector<MarcColumn> out;
    std::size_t k = 0;
    while (k < rCodes.size()) {
        const auto& [code, label] = rCodes[k];
        const std::int64_t layer = code / 1000, base = code % 1000;
        const std::string suffix = layer ? "@layer" + std::to_string(layer) : "";
        if (const char* tensor = marc_tensor_name(base); tensor && k + 5 < rCodes.size()) {
            bool consecutive = true;
            for (std::size_t j = 0; j < 6; ++j)
                consecutive =
                    consecutive && rCodes[k + j].first == code + static_cast<std::int64_t>(j);
            if (consecutive) {
                out.push_back({(label.empty() ? std::string(tensor) : label) + suffix, k, 6});
                k += 6;
                continue;
            }
        }
        out.push_back({(label.empty() ? marc_scalar_name(base) : label) + suffix, k, 1});
        ++k;
    }
    return out;
}

void marc_scalar_field(Mesh& rMesh, const char* pName, DType Type, double Value) {
    NDArray a(Type, {1});
    if (Type == DType::Float64)
        a.As<double>()[0] = Value;
    else
        a.As<std::int64_t>()[0] = static_cast<std::int64_t>(Value);
    rMesh.AddFieldData(pName, std::move(a));
}

Mesh marc_read_t19(const std::string& rPath, const ReadOptions& rOptions) {
    const MarcPost post(rPath);
    const std::size_t n = post.Increments().size();
    const std::size_t index = n ? rOptions.ResolveTimeStep(n) : 0;
    // An increment that remeshes (newmo, BLOCK 517) repeats the model blocks
    // 502 to 514 (BLOCK 519): a step's mesh is the latest model at or before it.
    std::vector<MarcBlock> remeshed;
    for (std::size_t k = 0; n && k <= index; ++k) {
        const auto& inc = post.Increments()[k];
        if (!post.Info(inc).mNewModel)
            continue;
        remeshed.clear();
        for (const MarcBlock& b : inc)
            if (b.mNumber / 100 >= 502 && b.mNumber / 100 <= 514)
                remeshed.push_back(b);
    }
    const std::vector<MarcBlock>* model = remeshed.empty() ? nullptr : &remeshed;
    std::vector<std::int64_t> lm;
    std::vector<std::pair<std::int64_t, std::string>> codes;
    post.Header(lm, codes);
    if (model && std::any_of(remeshed.begin(), remeshed.end(),
                             [](const MarcBlock& rB) { return rB.mNumber / 100 == 502; })) {
        // A remeshed model repeats the header; the post codes stay the first's
        // when it does not repeat them.
        std::vector<std::int64_t> lm2;
        std::vector<std::pair<std::int64_t, std::string>> codes2;
        post.Header(lm2, codes2, model);
        lm = lm2;
        if (!codes2.empty())
            codes = codes2;
    }
    const std::int64_t npost = lm[0], numnp = lm[1], numel = lm[2];
    const std::int64_t nstres = std::max<std::int64_t>(lm[4], 1);
    std::vector<std::int64_t> node_ids;
    std::vector<std::array<double, 3>> coords;
    std::vector<MarcElement> elements;
    std::vector<MarcSet> sets;
    post.Model(lm, node_ids, coords, elements, sets, model);
    std::vector<MarcLocation> locs;
    Mesh mesh = marc_build(kMarcT19, node_ids, coords, elements, sets, locs);
    if (n == 0) {
        if (rOptions.mTimeStep != 0 && rOptions.mTimeStep != -1)
            marc_fail(kMarcT19, "time step " + std::to_string(rOptions.mTimeStep) +
                                    " is out of range: the file has no increments");
        return mesh;
    }
    const auto& blocks = post.Increments()[index];
    const MarcIncrement info = post.Info(blocks);
    marc_scalar_field(mesh, "meshio:time", DType::Float64, info.mTime);
    marc_scalar_field(mesh, "marc:increment", DType::Int64, static_cast<double>(info.mInc));
    marc_scalar_field(mesh, "marc:subincrement", DType::Int64, static_cast<double>(info.mIncsub));
    if (rOptions.mPointsOnly)
        return mesh;
    const auto n_nodes = static_cast<std::size_t>(std::max<std::int64_t>(numnp, 0));
    for (const MarcBlock& b : blocks) {
        const std::int64_t family = b.mNumber / 100;
        MarcRecords r = post.Reader(b);
        if (family == 524) {
            const std::int64_t nnqnod = r.Ints(2)[0];
            for (std::int64_t q = 0; q < nnqnod; ++q) {
                const std::string& line = r.Line();
                const std::string name = marc_trim(
                    std::string_view(line).substr(0, std::min<std::size_t>(48, line.size())));
                const auto ivec = r.Ints(12);
                if (ivec[6] != -1)
                    continue;
                const auto ncomp = static_cast<std::size_t>(std::max<std::int64_t>(ivec[3], 1));
                std::vector<double> real = r.Reals(n_nodes * ncomp), imag;
                if (ivec[5] == 4 || ivec[5] == 5)
                    imag = r.Reals(n_nodes * ncomp);
                for (int part = 0; part < 2; ++part) {
                    const std::vector<double>& data = part ? imag : real;
                    const std::string key = part ? name + "@imag" : name;
                    if ((part && imag.empty()) || !rOptions.WantsArray(key))
                        continue;
                    NDArray a = ncomp == 1 ? NDArray(DType::Float64, {n_nodes})
                                           : NDArray(DType::Float64, {n_nodes, ncomp});
                    std::copy(data.begin(), data.end(), a.As<double>());
                    mesh.AddPointData(key, std::move(a));
                }
            }
        } else if (family == 523 && info.mJantyp > 100 && npost > 0 && numel > 0) {
            const auto np = static_cast<std::size_t>(npost), ns = static_cast<std::size_t>(nstres);
            std::vector<double> values;
            values.reserve(static_cast<std::size_t>(numel) * ns * np);
            for (std::int64_t e = 0; e < numel; ++e)
                for (std::size_t p = 0; p < ns; ++p) {
                    const auto record = r.Reals(np);
                    values.insert(values.end(), record.begin(), record.end());
                }
            for (const MarcColumn& c : marc_columns(codes)) {
                if (!rOptions.WantsArray(c.mName))
                    continue;
                // With several integration points, flattened point-major so that
                // any writer holds it; its (points, components) is the layout.
                const std::size_t width = ns * c.mWidth;
                std::vector<NDArray> per_block;
                for (std::size_t bi = 0; bi < mesh.NumCellBlocks(); ++bi) {
                    const std::size_t rows = mesh.Cells(bi).NumCells();
                    NDArray a = width > 1 ? NDArray(DType::Float64, {rows, width})
                                          : NDArray(DType::Float64, {rows});
                    std::fill(a.As<double>(), a.As<double>() + a.Size(),
                              std::numeric_limits<double>::quiet_NaN());
                    per_block.push_back(std::move(a));
                }
                for (std::size_t e = 0; e < locs.size() && e < static_cast<std::size_t>(numel);
                     ++e) {
                    if (locs[e].mBlock < 0)
                        continue;
                    double* target =
                        per_block[static_cast<std::size_t>(locs[e].mBlock)].As<double>() +
                        static_cast<std::size_t>(locs[e].mRow) * width;
                    for (std::size_t p = 0; p < ns; ++p)
                        for (std::size_t j = 0; j < c.mWidth; ++j)
                            target[p * c.mWidth + j] = values[(e * ns + p) * np + c.mFirst + j];
                }
                if (ns > 1) {
                    NDArray layout(DType::Int64, {2});
                    layout.As<std::int64_t>()[0] = static_cast<std::int64_t>(ns);
                    layout.As<std::int64_t>()[1] = static_cast<std::int64_t>(c.mWidth);
                    mesh.AddFieldData("marc:layout:" + c.mName, std::move(layout));
                }
                mesh.AddCellData(c.mName, std::move(per_block));
            }
        }
    }
    return mesh;
}

}  // namespace

bool is_marc_deck(std::string_view Head) {
    bool first = true;
    std::size_t pos = 0;
    while (pos < Head.size()) {
        std::size_t eol = Head.find('\n', pos);
        if (eol == std::string_view::npos)
            eol = Head.size();
        const std::string stripped = marc_trim(Head.substr(pos, eol - pos));
        pos = eol + 1;
        if (stripped.empty() || stripped[0] == '$')
            continue;
        const auto words = marc_words(marc_lower(stripped));
        const std::string word = words.empty() ? "" : words[0];
        if (first) {
            if (stripped.find('=') != std::string::npos || !marc_is_parameter(word))
                return false;
            first = false;
            if (word == "end")
                return true;
            continue;
        }
        if (word == "end" || word == "connectivity" || word == "coordinates")
            return true;
    }
    return false;
}

Mesh read_marc(const std::string& rPath) {
    const std::vector<std::string> lines = marc_deck_lines(rPath);
    std::string head;
    for (const std::string& line : lines) {
        if (head.size() >= 65536)
            break;
        head += line;
        head += '\n';
    }
    if (!is_marc_deck(std::string_view(head).substr(0, std::min<std::size_t>(head.size(), 65536))))
        marc_fail(kMarcDat, "not a Marc input deck (no Marc parameter opens the file)");
    MarcDeck deck = marc_parse_deck(lines);
    if (deck.mNodes.empty() && deck.mElements.empty())
        marc_fail(kMarcDat, "no COORDINATES or CONNECTIVITY found");
    std::map<std::pair<bool, std::string>, std::vector<std::int64_t>> known;
    for (MarcSet& s : deck.mSets) {
        if (s.mSideKind)
            continue;  // members already read
        std::unordered_map<std::string, std::vector<std::int64_t>> refs;
        for (const auto& [key, members] : known)
            if (key.first == s.mElements)
                refs.emplace(key.second, members);
        s.mMembers = marc_expand(s.mTokens, s.mName, refs, kMarcDat);
        known[{s.mElements, marc_lower(s.mName)}] = s.mMembers;
    }
    std::vector<std::array<double, 3>> coords;
    coords.reserve(deck.mNodeOrder.size());
    for (std::int64_t id : deck.mNodeOrder)
        coords.push_back(deck.mNodes.at(id));
    std::vector<MarcLocation> locs;
    return marc_build(kMarcDat, deck.mNodeOrder, coords, deck.mElements, deck.mSets, locs);
}

Mesh read_marc_t19(const std::string& rPath, const ReadOptions& rOptions) {
    return marc_read_t19(rPath, rOptions);
}

MeshMetadata read_marc_t19_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    ReadOptions options = rOptions;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(marc_read_t19(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "marc_t19";
    meta.mTimeValues = MarcPost(rPath).Times();
    return meta;
}

}  // namespace meshioplusplus
