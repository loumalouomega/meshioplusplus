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
// Ansys MAPDL coded database (.cdb) reader and writer. The block layouts follow
// MAPDL's CDWRITE documentation; the element categories and degenerate-shape
// rules follow the open readers pymapdl-reader and mapdl-archive (MIT). See
// doc/formats/ansysinp.md. Python twin: src/python/meshioplusplus/ansysInp/_ansysInp.py.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// ---- element categories ------------------------------------------------------
// How an element routine number's nodes become a cell (pymapdl-reader's
// ETYPE_MAP): a point, a line (quadratic when its third node is set), a shell or
// plane (triangle when K == L), a degenerate brick, a native tetrahedron, or a
// line whose extra nodes are orientation nodes.
enum class AnsCategory { Skip, Point, Line, Shell, Brick, Tet, LinearLine };

AnsCategory ans_category(int Routine) {
    static const std::unordered_map<int, AnsCategory> kTable = [] {
        std::unordered_map<int, AnsCategory> m;
        for (int n : {7, 21, 71, 175})
            m[n] = AnsCategory::Point;
        for (int n :
             {1,   3,   4,   8,   10,  11,  12,  14,  16,  17,  18,  20,  23,  24,  31,  32,  33,
              34,  37,  38,  39,  40,  44,  59,  60,  61,  66,  68,  116, 126, 129, 151, 153, 156,
              161, 169, 171, 172, 176, 177, 178, 180, 189, 208, 209, 250, 251, 280, 288, 289, 290})
            m[n] = AnsCategory::Line;
        for (int n :
             {2,   13,  22,  25,  28,  29,  35,  41,  42,  43,  51,  53,  54,  55,  57,  63,
              67,  75,  77,  78,  79,  81,  82,  83,  88,  91,  93,  99,  106, 115, 118, 121,
              130, 131, 132, 136, 143, 152, 154, 155, 157, 163, 170, 173, 174, 181, 182, 183,
              212, 213, 218, 219, 222, 223, 230, 233, 238, 252, 281, 282, 283, 292, 293})
            m[n] = AnsCategory::Shell;
        for (int n : {5,   30,  45,  46,  62,  64,  65,  69,  70,  80,  89,  90,  95,  96,
                      97,  100, 101, 102, 103, 104, 105, 107, 108, 117, 120, 122, 164, 185,
                      186, 190, 192, 215, 220, 226, 231, 236, 239, 272, 273, 278, 279})
            m[n] = AnsCategory::Brick;
        for (int n : {87, 92, 98, 119, 123, 140, 168, 187, 221, 227, 232, 237, 240, 285, 291})
            m[n] = AnsCategory::Tet;
        for (int n : {188, 214, 216, 217})
            m[n] = AnsCategory::LinearLine;
        return m;
    }();
    const auto it = kTable.find(Routine);
    return it == kTable.end() ? AnsCategory::Skip : it->second;
}

// MESH200's shape comes from KEYOPT(1) (pymapdl-reader's MESH200_MAP).
AnsCategory ans_mesh200_category(int Keyopt1) {
    if (Keyopt1 >= 0 && Keyopt1 <= 3)
        return AnsCategory::Line;
    if (Keyopt1 >= 4 && Keyopt1 <= 7)
        return AnsCategory::Shell;
    if (Keyopt1 == 8 || Keyopt1 == 9)
        return AnsCategory::Tet;
    if (Keyopt1 == 10 || Keyopt1 == 11)
        return AnsCategory::Brick;
    return AnsCategory::Skip;
}

// A cell resolved from one element row: the meshio++ type and, per node of that
// type, the element's slot (a slot past the row's end is a missing midside).
struct AnsShape {
    const char* mType = nullptr;
    std::vector<int> mSlots;
};

// Slot maps of the degenerate forms, from the ANSYS brick numbering (corners
// I..P = 0..7; mid-edges Q R S T = 8..11 on IJ JK KL LI, U V W X = 12..15 on MN
// NO OP PM, Y Z A B = 16..19 on IM JN KO LP) into meshio++'s (VTK) orders.
AnsShape ans_resolve(AnsCategory Category, const std::vector<std::int64_t>& rNodes) {
    const std::size_t n = rNodes.size();
    const auto at = [&](std::size_t k) { return k < n ? rNodes[k] : 0; };
    const auto slots = [](std::initializer_list<int> l) { return std::vector<int>(l); };
    switch (Category) {
        case AnsCategory::Point:
            if (n >= 1)
                return {"vertex", slots({0})};
            break;
        case AnsCategory::Line:
            if (n >= 3 && at(2) > 0)
                return {"line3", slots({0, 1, 2})};
            if (n >= 2)
                return {"line", slots({0, 1})};
            break;
        case AnsCategory::LinearLine:
            if (n >= 2)
                return {"line", slots({0, 1})};
            break;
        case AnsCategory::Shell:
            if (n == 3)
                return {"triangle", slots({0, 1, 2})};
            if (n == 6)
                return {"triangle6", slots({0, 1, 2, 3, 4, 5})};
            if (n > 5) {  // 8-node (5 is a quad plus an orientation node); absent
                          // trailing midsides are missing ones
                if (at(2) == at(3))
                    return {"triangle6", slots({0, 1, 2, 4, 5, 7})};
                return {"quad8", slots({0, 1, 2, 3, 4, 5, 6, 7})};
            }
            if (n >= 4) {
                if (at(2) == at(3))
                    return {"triangle", slots({0, 1, 2})};
                return {"quad", slots({0, 1, 2, 3})};
            }
            break;
        case AnsCategory::Tet:
            if (n > 4)
                return {"tetra10", slots({0, 1, 2, 3, 4, 5, 6, 7, 8, 9})};
            if (n >= 4)
                return {"tetra", slots({0, 1, 2, 3})};
            break;
        case AnsCategory::Brick: {
            if (n < 8)
                break;
            const bool quad = n > 8;
            if (at(6) != at(7))  // hexahedron
                return quad ? AnsShape{"hexahedron20",
                                       slots({0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                              10, 11, 12, 13, 14, 15, 16, 17, 18, 19})}
                            : AnsShape{"hexahedron", slots({0, 1, 2, 3, 4, 5, 6, 7})};
            if (at(5) != at(6))  // wedge: K == L, O == P
                return quad ? AnsShape{"wedge15",
                                       slots({0, 1, 2, 4, 5, 6, 8, 9, 11, 12, 13, 15, 16, 17, 18})}
                            : AnsShape{"wedge", slots({0, 1, 2, 4, 5, 6})};
            if (at(2) != at(3))  // pyramid: M == N == O == P
                return quad ? AnsShape{"pyramid13",
                                       slots({0, 1, 2, 3, 4, 8, 9, 10, 11, 16, 17, 18, 19})}
                            : AnsShape{"pyramid", slots({0, 1, 2, 3, 4})};
            // tetrahedron: K == L, M == N == O == P
            return quad ? AnsShape{"tetra10", slots({0, 1, 2, 4, 8, 9, 11, 16, 17, 18})}
                        : AnsShape{"tetra", slots({0, 1, 2, 4})};
        }
        case AnsCategory::Skip:
            break;
    }
    return {};
}

// Corner pairs of each mid-edge slot of a meshio++ type (for missing midside
// nodes, written as node 0), indexed by position in the type's node list.
std::vector<std::pair<int, int>> ans_midside_edges(std::string_view Type) {
    if (Type == "line3")
        return {{0, 1}};
    if (Type == "triangle6")
        return {{0, 1}, {1, 2}, {2, 0}};
    if (Type == "quad8")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    if (Type == "tetra10")
        return {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
    if (Type == "pyramid13")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    if (Type == "wedge15")
        return {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    if (Type == "hexahedron20")
        return {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    return {};
}

std::size_t ans_num_corners(std::string_view Type) {
    static const std::pair<std::string_view, std::size_t> kCorners[] = {
        {"vertex", 1}, {"line", 2},    {"triangle", 3}, {"quad", 4},
        {"tetra", 4},  {"pyramid", 5}, {"wedge", 6},    {"hexahedron", 8},
    };
    for (const auto& [prefix, c] : kCorners)
        if (Type.substr(0, prefix.size()) == prefix)
            return c;
    return 0;
}

std::string ans_upper(std::string_view Text) {
    std::string out(Text);
    for (char& c : out)
        if (c >= 'a' && c <= 'z')
            c = static_cast<char>(c - 'a' + 'A');
    return out;
}

std::string ans_strip(std::string_view Text) {
    const std::size_t a = Text.find_first_not_of(" \t\r\n");
    if (a == std::string_view::npos)
        return {};
    const std::size_t b = Text.find_last_not_of(" \t\r\n");
    return std::string(Text.substr(a, b - a + 1));
}

// Comma-separated fields of a command line, stripped (`ET, 4, 186`).
std::vector<std::string> ans_commas(std::string_view Line) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = Line.find(',', start);
        out.push_back(ans_strip(Line.substr(
            start, comma == std::string_view::npos ? std::string_view::npos : comma - start)));
        if (comma == std::string_view::npos)
            return out;
        start = comma + 1;
    }
}

std::optional<std::int64_t> ans_int(const std::string& rText) {
    if (rText.empty())
        return std::nullopt;
    const char* end = nullptr;
    const double v = detail::parse_double(rText.c_str(), end);
    if (end != rText.c_str() + rText.size())
        return std::nullopt;
    return static_cast<std::int64_t>(v);
}

// An element type given by number (`186`) or name (`SOLID186`): the routine.
int ans_routine(const std::string& rText) {
    if (const auto v = ans_int(rText))
        return static_cast<int>(*v);
    std::size_t k = rText.size();
    while (k > 0 && rText[k - 1] >= '0' && rText[k - 1] <= '9')
        --k;
    if (k == rText.size())
        return -1;
    return static_cast<int>(ans_int(rText.substr(k)).value_or(-1));
}

[[noreturn]] void ans_fail(std::size_t Line, const std::string& rWhat) {
    throw ReadError("Ansys .cdb: line " + std::to_string(Line + 1) + ": " + rWhat);
}

std::int64_t ans_field_int(const std::vector<std::string>& rFields, std::size_t K,
                           std::size_t Line) {
    if (K >= rFields.size() || rFields[K].empty())
        return 0;
    const auto v = ans_int(rFields[K]);
    if (!v)
        ans_fail(Line, "bad integer '" + rFields[K] + "'");
    return *v;
}

double ans_field_real(const std::vector<std::string>& rFields, std::size_t K, std::size_t Line) {
    if (K >= rFields.size() || rFields[K].empty())
        return 0.0;
    const std::string& text = rFields[K];
    const char* end = nullptr;
    const double v = detail::parse_double(text.c_str(), end);
    if (end != text.c_str() + text.size())
        ans_fail(Line, "bad number '" + text + "'");
    return v;
}

bool ans_is_terminator(const std::string& rLine) {
    const std::string s = ans_strip(rLine);
    if (s == "-1")
        return true;
    const std::string up = ans_upper(s);
    return up.rfind("N,", 0) == 0 && up.find("LOC") != std::string::npos;
}

// A command line (`FINISH`, `CMBLOCK,...`) rather than a block's data line.
bool ans_is_command(const std::string& rLine) {
    const std::size_t a = rLine.find_first_not_of(" \t");
    if (a == std::string::npos)
        return false;
    const char c = rLine[a];
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '/' || c == '!';
}

// `KEYOPT,` or any abbreviation of it MAPDL accepts (`KEYO,`, `KEYOP,`).
bool ans_is_keyopt(const std::string& rUpper) {
    const std::size_t comma = rUpper.find(',');
    if (comma == std::string::npos || comma < 4 || comma > 6)
        return false;
    return std::string_view("KEYOPT").substr(0, comma) == std::string_view(rUpper).substr(0, comma);
}

struct AnsElement {
    int mSlot = 0;
    std::int64_t mMat = 0, mReal = 0, mSecnum = 0, mId = 0;
    std::vector<std::int64_t> mNodes;
};

struct AnsComponent {
    std::string mName;
    bool mNodes = false;
    std::vector<std::int64_t> mIds;
};

struct AnsDeck {
    std::map<int, int> mRoutine;                // ET slot -> routine
    std::map<int, std::map<int, int>> mKeyopt;  // ET slot -> keyopt k -> value
    std::vector<std::int64_t> mNodeIds;
    std::vector<double> mCoords;
    bool mDroppedRotations = false;
    std::vector<AnsElement> mElements;
    std::vector<AnsComponent> mComponents;
    std::size_t mNonSolidBlocks = 0;
};

// The format line after a block header, parsed.
std::vector<detail::CardField> ans_format(const std::vector<std::string>& rLines, std::size_t K) {
    if (K >= rLines.size())
        ans_fail(K, "a block header is not followed by its format line");
    try {
        return detail::parse_fortran_format(ans_strip(rLines[K]));
    } catch (const ReadError& e) {
        ans_fail(K, e.what());
    }
}

AnsDeck ans_parse(const std::vector<std::string>& rLines) {
    AnsDeck deck;
    bool saw_block = false;
    std::size_t i = 0;
    const std::size_t n = rLines.size();
    while (i < n) {
        // A command line, less any trailing `!` comment.
        const std::string line =
            ans_strip(std::string_view(rLines[i]).substr(0, rLines[i].find('!')));
        const std::string up = ans_upper(line);
        if (up.rfind("ET,", 0) == 0) {
            const auto p = ans_commas(line);
            if (p.size() >= 3) {
                const auto slot = ans_int(p[1]);
                const int routine = ans_routine(ans_upper(p[2]));
                if (slot && routine > 0)
                    deck.mRoutine[static_cast<int>(*slot)] = routine;
            }
            ++i;
        } else if (ans_is_keyopt(up)) {
            const auto p = ans_commas(line);
            if (p.size() >= 4) {
                const auto slot = ans_int(p[1]), k = ans_int(p[2]), v = ans_int(p[3]);
                if (slot && k && v)
                    deck.mKeyopt[static_cast<int>(*slot)][static_cast<int>(*k)] =
                        static_cast<int>(*v);
            }
            ++i;
        } else if (up.rfind("ETBLOCK", 0) == 0) {
            saw_block = true;
            const auto fields = ans_format(rLines, i + 1);
            i += 2;
            while (i < n && !ans_is_terminator(rLines[i])) {
                const auto f = detail::split_fixed(rLines[i], fields);
                if (f.size() >= 2) {
                    const int slot = static_cast<int>(ans_field_int(f, 0, i));
                    deck.mRoutine[slot] = static_cast<int>(ans_field_int(f, 1, i));
                    for (std::size_t k = 2; k < f.size() && k < 20; ++k)
                        if (const auto v = ans_int(f[k]); v && *v != 0)
                            deck.mKeyopt[slot][static_cast<int>(k - 1)] = static_cast<int>(*v);
                }
                ++i;
            }
            ++i;
        } else if (up.rfind("NBLOCK", 0) == 0) {
            saw_block = true;
            const auto fields = ans_format(rLines, i + 1);
            std::size_t n_int = 0;
            while (n_int < fields.size() && fields[n_int].mKind == 'i')
                ++n_int;
            i += 2;
            while (i < n && !ans_is_terminator(rLines[i])) {
                const auto f = detail::split_fixed(rLines[i], fields);
                if (f.empty() || f[0].empty()) {
                    ++i;
                    continue;
                }
                deck.mNodeIds.push_back(ans_field_int(f, 0, i));
                for (std::size_t d = 0; d < 3; ++d)
                    deck.mCoords.push_back(ans_field_real(f, n_int + d, i));
                for (std::size_t d = 3; n_int + d < f.size(); ++d)
                    if (ans_field_real(f, n_int + d, i) != 0.0)
                        deck.mDroppedRotations = true;
                ++i;
            }
            ++i;
        } else if (up.rfind("EBLOCK", 0) == 0) {
            saw_block = true;
            const auto header = ans_commas(up);
            const bool solid = header.size() > 2 && header[2] == "SOLID";
            const auto fields = ans_format(rLines, i + 1);
            i += 2;
            if (!solid) {
                ++deck.mNonSolidBlocks;
                while (i < n && !ans_is_terminator(rLines[i]))
                    ++i;
                ++i;
                continue;
            }
            while (i < n && !ans_is_terminator(rLines[i])) {
                const auto f = detail::split_fixed(rLines[i], fields);
                if (f.size() < 11) {
                    ++i;
                    continue;
                }
                AnsElement e;
                e.mMat = ans_field_int(f, 0, i);
                e.mSlot = static_cast<int>(ans_field_int(f, 1, i));
                e.mReal = ans_field_int(f, 2, i);
                e.mSecnum = ans_field_int(f, 3, i);
                const auto count =
                    static_cast<std::size_t>(std::max<std::int64_t>(ans_field_int(f, 8, i), 0));
                e.mId = ans_field_int(f, 10, i);
                for (std::size_t k = 11; k < f.size() && e.mNodes.size() < count; ++k)
                    e.mNodes.push_back(ans_field_int(f, k, i));
                ++i;
                while (e.mNodes.size() < count && i < n && !ans_is_terminator(rLines[i])) {
                    const auto more = detail::split_fixed(rLines[i], fields);
                    for (std::size_t k = 0; k < more.size() && e.mNodes.size() < count; ++k)
                        e.mNodes.push_back(ans_field_int(more, k, i));
                    ++i;
                }
                deck.mElements.push_back(std::move(e));
            }
            ++i;
        } else if (up.rfind("CMBLOCK", 0) == 0) {
            saw_block = true;
            const auto header = ans_commas(line);
            if (header.size() < 3)
                ans_fail(i, "a CMBLOCK needs a name and an entity type");
            AnsComponent comp;
            comp.mName = header[1];
            const std::string entity = ans_upper(header[2]);
            const bool known = entity == "NODE" || entity.rfind("ELEM", 0) == 0;
            comp.mNodes = entity == "NODE";
            const std::size_t count =
                header.size() > 3 ? static_cast<std::size_t>(ans_int(header[3]).value_or(0)) : 0;
            const auto fields = ans_format(rLines, i + 1);
            i += 2;
            std::vector<std::int64_t> raw;
            // A short block (a header count too large) ends at the next command.
            while (i < n && raw.size() < count && !ans_is_command(rLines[i])) {
                const auto f = detail::split_fixed(rLines[i], fields);
                if (f.empty())
                    break;
                for (std::size_t k = 0; k < f.size() && raw.size() < count; ++k)
                    if (!f[k].empty())
                        raw.push_back(ans_field_int(f, k, i));
                ++i;
            }
            // A negative value closes a range opened by the value before it.
            for (std::size_t k = 0; k < raw.size(); ++k) {
                if (raw[k] >= 0) {
                    comp.mIds.push_back(raw[k]);
                    continue;
                }
                if (comp.mIds.empty())
                    ans_fail(i, "CMBLOCK '" + comp.mName + "' opens with a range end");
                for (std::int64_t v = comp.mIds.back() + 1; v <= -raw[k]; ++v)
                    comp.mIds.push_back(v);
            }
            if (known)
                deck.mComponents.push_back(std::move(comp));
            else
                log::debug("Ansys .cdb: {} component '{}' skipped", entity, comp.mName);
        } else {
            ++i;
        }
    }
    if (!saw_block)
        throw ReadError("Ansys .cdb: no NBLOCK, EBLOCK or CMBLOCK found");
    return deck;
}

std::vector<std::string> ans_read_lines(const std::string& rPath) {
    auto f = detail::make_classic_ifstream(rPath);
    if (!f)
        throw ReadError("Could not open ansysInp file: " + rPath);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

}  // namespace

Mesh read_ansysinp(const std::string& rPath, const ReadOptions& rOptions, AnsysInfo& rInfo) {
    AnsDeck deck = ans_parse(ans_read_lines(rPath));
    if (deck.mNonSolidBlocks)
        log::warn(
            "Ansys .cdb: {} non-solid EBLOCK(s) (MAPDL writes only the SOLID layout) "
            "were skipped",
            deck.mNonSolidBlocks);
    if (deck.mDroppedRotations)
        log::warn("Ansys .cdb: nodal rotation angles are not kept");

    std::vector<double> coords = std::move(deck.mCoords);
    std::unordered_map<std::int64_t, std::int64_t> node_index;
    for (std::size_t k = 0; k < deck.mNodeIds.size(); ++k)
        node_index.emplace(deck.mNodeIds[k], static_cast<std::int64_t>(k));
    // Missing midside nodes (node 0) become new points at their edge midpoints.
    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> midsides;
    std::size_t n_missing = 0;

    struct Block {
        std::string mType;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mRoutine, mSlot, mMat, mReal, mSecnum;
        std::size_t Rows() const { return mSlot.size(); }
    };
    std::vector<Block> blocks;
    std::map<std::string, std::size_t> block_of;
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> element_loc;
    std::map<int, std::size_t> skipped;  // routine -> count

    for (const AnsElement& e : deck.mElements) {
        const auto rt = deck.mRoutine.find(e.mSlot);
        if (rt == deck.mRoutine.end())
            throw ReadError("Ansys .cdb: element " + std::to_string(e.mId) + " uses element type " +
                            std::to_string(e.mSlot) + ", which no ET or ETBLOCK defines");
        const int routine = rt->second;
        AnsCategory category = ans_category(routine);
        if (routine == 200) {
            const auto& ko = deck.mKeyopt[e.mSlot];
            const auto k1 = ko.find(1);
            category = ans_mesh200_category(k1 == ko.end() ? 0 : k1->second);
        }
        const AnsShape shape = ans_resolve(category, e.mNodes);
        if (!shape.mType) {
            if (!rOptions.mLenient)
                throw ReadError("Ansys .cdb: element " + std::to_string(e.mId) + " (element type " +
                                std::to_string(routine) + ", " + std::to_string(e.mNodes.size()) +
                                " nodes) has no meshio++ cell type; read with lenient to skip it");
            ++skipped[routine];
            continue;
        }
        const auto [it, fresh] = block_of.emplace(shape.mType, blocks.size());
        if (fresh)
            blocks.push_back(Block{shape.mType, {}, {}, {}, {}, {}, {}});
        Block& b = blocks[it->second];
        const std::size_t corners = ans_num_corners(shape.mType);
        const auto edges = ans_midside_edges(shape.mType);
        std::vector<std::int64_t> row;
        for (int slot : shape.mSlots) {
            // A row cut short omits trailing midside nodes: they read as node 0.
            const auto k = static_cast<std::size_t>(slot);
            const std::int64_t id = k < e.mNodes.size() ? e.mNodes[k] : 0;
            if (id == 0 && row.size() >= corners) {
                row.push_back(-1);  // resolved below, once the corners are known
                continue;
            }
            const auto found = node_index.find(id);
            if (found == node_index.end())
                throw ReadError("Ansys .cdb: element " + std::to_string(e.mId) +
                                " names undefined node " + std::to_string(id));
            row.push_back(found->second);
        }
        for (std::size_t k = corners; k < row.size(); ++k) {
            if (row[k] >= 0)
                continue;
            const auto [a, c] = edges[k - corners];
            std::int64_t p = row[static_cast<std::size_t>(a)], q = row[static_cast<std::size_t>(c)];
            const auto key = std::make_pair(std::min(p, q), std::max(p, q));
            auto [mit, added] = midsides.emplace(key, static_cast<std::int64_t>(coords.size() / 3));
            if (added) {
                for (std::size_t d = 0; d < 3; ++d)
                    coords.push_back(0.5 * (coords[static_cast<std::size_t>(p) * 3 + d] +
                                            coords[static_cast<std::size_t>(q) * 3 + d]));
                ++n_missing;
            }
            row[k] = mit->second;
        }
        element_loc[e.mId] = {it->second, b.Rows()};
        b.mConn.insert(b.mConn.end(), row.begin(), row.end());
        b.mRoutine.push_back(routine);
        b.mSlot.push_back(e.mSlot);
        b.mMat.push_back(e.mMat);
        b.mReal.push_back(e.mReal);
        b.mSecnum.push_back(e.mSecnum);
    }
    for (const auto& [routine, count] : skipped)
        log::warn("Ansys .cdb: {} element(s) of type {} skipped (no meshio++ cell type)", count,
                  routine);
    if (n_missing)
        log::warn("Ansys .cdb: {} missing midside node(s) placed at their edge midpoints",
                  n_missing);

    Mesh mesh;
    NDArray points(DType::Float64, {coords.size() / 3, 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    const auto column = [](const std::vector<std::int64_t>& rValues) {
        NDArray a(DType::Int64, {rValues.size()});
        std::copy(rValues.begin(), rValues.end(), a.As<std::int64_t>());
        return a;
    };
    std::vector<NDArray> routine_data, slot_data, mat_data, real_data, secnum_data;
    std::vector<std::int64_t> bases;
    std::int64_t base = 0;
    for (const Block& b : blocks) {
        const std::size_t k = b.mConn.size() / std::max<std::size_t>(b.Rows(), 1);
        NDArray conn(DType::Int64, {b.Rows(), k});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
        routine_data.push_back(column(b.mRoutine));
        slot_data.push_back(column(b.mSlot));
        mat_data.push_back(column(b.mMat));
        real_data.push_back(column(b.mReal));
        secnum_data.push_back(column(b.mSecnum));
        bases.push_back(base);
        base += static_cast<std::int64_t>(b.Rows());
    }
    if (!blocks.empty()) {
        mesh.AddCellData("ansys:element", std::move(routine_data));
        mesh.AddCellData("ansys:type", std::move(slot_data));
        mesh.AddCellData("ansys:mat", std::move(mat_data));
        mesh.AddCellData("ansys:real", std::move(real_data));
        mesh.AddCellData("ansys:secnum", std::move(secnum_data));
    }

    // Components: point and cell regions, and the legacy side channel.
    for (const AnsComponent& c : deck.mComponents) {
        std::vector<std::int64_t> entries;
        if (c.mNodes) {
            for (std::int64_t id : c.mIds)
                if (const auto it = node_index.find(id); it != node_index.end())
                    entries.push_back(it->second);
            rInfo.mPointSets[c.mName] = entries;
        } else {
            std::vector<std::vector<std::int64_t>> per(blocks.size());
            for (std::int64_t id : c.mIds)
                if (const auto it = element_loc.find(id); it != element_loc.end()) {
                    per[it->second.first].push_back(static_cast<std::int64_t>(it->second.second));
                    entries.push_back(bases[it->second.first] +
                                      static_cast<std::int64_t>(it->second.second));
                }
            rInfo.mCellSets[c.mName] = per;
        }
        mesh.AddRegion(Region(c.mName, c.mNodes ? RegionKind::Point : RegionKind::Cell, -1, -1,
                              column(entries)));
    }
    return mesh;
}

Mesh read_ansysinp(const std::string& rPath, AnsysInfo& rInfo) {
    return read_ansysinp(rPath, ReadOptions{}, rInfo);
}

namespace {

// The degenerate or native layout of `Type` under an element of `Category`.
std::optional<std::vector<int>> ans_layout(std::string_view Type, AnsCategory Category) {
    using V = std::vector<int>;
    if (Category == AnsCategory::Brick) {
        if (Type == "hexahedron")
            return V{0, 1, 2, 3, 4, 5, 6, 7};
        if (Type == "hexahedron20")
            return V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
        if (Type == "wedge")
            return V{0, 1, 2, 2, 3, 4, 5, 5};
        if (Type == "wedge15")  // K=L, O=P, S=K, W=O, B=A
            return V{0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 2, 8, 9, 10, 5, 11, 12, 13, 14, 14};
        if (Type == "pyramid")
            return V{0, 1, 2, 3, 4, 4, 4, 4};
        if (Type == "pyramid13")  // M=N=O=P, U=V=W=X=M
            return V{0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, 4, 4, 4, 4, 9, 10, 11, 12};
        if (Type == "tetra")
            return V{0, 1, 2, 2, 3, 3, 3, 3};
        if (Type == "tetra10")  // K=L, M..P, S=K, U..X=M, B=A
            return V{0, 1, 2, 2, 3, 3, 3, 3, 4, 5, 2, 6, 3, 3, 3, 3, 7, 8, 9, 9};
    } else if (Category == AnsCategory::Tet) {
        if (Type == "tetra")
            return V{0, 1, 2, 3};
        if (Type == "tetra10")
            return V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    } else if (Category == AnsCategory::Shell) {
        if (Type == "quad")
            return V{0, 1, 2, 3};
        if (Type == "quad8")
            return V{0, 1, 2, 3, 4, 5, 6, 7};
        if (Type == "triangle")
            return V{0, 1, 2, 2};
        if (Type == "triangle6")  // K=L, the K-L midside = K
            return V{0, 1, 2, 2, 3, 4, 2, 5};
    } else if (Category == AnsCategory::Line || Category == AnsCategory::LinearLine) {
        if (Type == "line")
            return V{0, 1};
        if (Type == "line3" && Category == AnsCategory::Line)
            return V{0, 1, 2};
    } else if (Category == AnsCategory::Point) {
        if (Type == "vertex")
            return V{0};
    }
    return std::nullopt;
}

int ans_default_routine(std::string_view Type) {
    static const std::pair<std::string_view, int> kDefault[] = {
        {"vertex", 21},     {"line", 188},       {"line3", 189},        {"triangle", 181},
        {"triangle6", 281}, {"quad", 181},       {"quad8", 281},        {"tetra", 285},
        {"tetra10", 187},   {"pyramid", 185},    {"pyramid13", 186},    {"wedge", 185},
        {"wedge15", 186},   {"hexahedron", 185}, {"hexahedron20", 186},
    };
    for (const auto& [type, routine] : kDefault)
        if (Type == type)
            return routine;
    return -1;
}

std::string ans_i(std::int64_t Value, int Width) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%*lld", Width, static_cast<long long>(Value));
    return buf;
}

// A per-cell integer from `Name` cell data, or `Default` when absent.
std::int64_t ans_cell_int(const Mesh& rMesh, const std::string& rName, std::size_t Block,
                          std::size_t Row, std::int64_t Default) {
    if (!rMesh.HasCellData(rName))
        return Default;
    const NDArray& a = rMesh.CellData(rName, Block);
    if (a.Size() <= Row)
        return Default;
    return detail::read_int(a, Row);
}

}  // namespace

void write_ansysinp(const std::string& rPath, const Mesh& rMesh, const AnsysInfo& rInfo) {
    const std::size_t n_blocks = rMesh.NumCellBlocks();
    // Every cell's (element routine, ET slot): kept from `ansys:element` and
    // `ansys:type` when they describe a layout the cell's type has, else the
    // default routine for its type in a fresh slot.
    std::map<int, int> slot_routine;  // slot -> routine
    std::map<int, int> routine_slot;  // routine -> slot, for fresh slots
    std::vector<std::vector<std::pair<int, int>>> cell_etype(n_blocks);
    bool dropped_etype = false;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const int fallback = ans_default_routine(type);
        if (cb.IsRagged() || fallback < 0)
            throw WriteError("Ansys .cdb writer: cell type '" + type + "' has no element type");
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            int routine = static_cast<int>(ans_cell_int(rMesh, "ansys:element", b, r, fallback));
            int slot = static_cast<int>(ans_cell_int(rMesh, "ansys:type", b, r, 0));
            if (!ans_layout(type, ans_category(routine)) ||
                (slot > 0 && slot_routine.count(slot) && slot_routine[slot] != routine)) {
                dropped_etype = dropped_etype || routine != fallback;
                routine = fallback;
                slot = 0;
            }
            if (slot > 0)
                slot_routine.emplace(slot, routine);
            cell_etype[b].emplace_back(routine, slot);
        }
    }
    for (auto& per : cell_etype)
        for (auto& [routine, slot] : per)
            if (slot == 0) {
                auto it = routine_slot.find(routine);
                if (it == routine_slot.end()) {
                    int fresh = 1;
                    while (slot_routine.count(fresh))
                        ++fresh;
                    slot_routine[fresh] = routine;
                    it = routine_slot.emplace(routine, fresh).first;
                }
                slot = it->second;
            }
    if (dropped_etype)
        log::warn(
            "Ansys .cdb writer: some ansys:element values do not fit their cells' shape; "
            "the default element type was written instead");

    // Components: cell and point regions, plus legacy side-channel sets not
    // already named by a region.
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> node_comps, elem_comps;
    std::set<std::pair<bool, std::string>> seen;
    const std::vector<std::int64_t> bases = [&] {
        std::vector<std::int64_t> out{0};
        for (const auto cb : rMesh.CellRange())
            out.push_back(out.back() + static_cast<std::int64_t>(cb.NumCells()));
        return out;
    }();
    std::size_t side_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Side) {
            ++side_regions;
            continue;
        }
        const bool nodes = reg.mKind == RegionKind::Point;
        if (!seen.insert({nodes, reg.mName}).second)
            continue;
        std::vector<std::int64_t> ids(reg.Entries(), reg.Entries() + reg.NumEntries());
        for (std::int64_t& v : ids)
            ++v;
        (nodes ? node_comps : elem_comps).emplace_back(reg.mName, std::move(ids));
    }
    for (const auto& [name, idx] : rInfo.mPointSets)
        if (seen.insert({true, name}).second) {
            std::vector<std::int64_t> ids;
            for (std::int64_t v : idx)
                ids.push_back(v + 1);
            node_comps.emplace_back(name, std::move(ids));
        }
    for (const auto& [name, per] : rInfo.mCellSets)
        if (seen.insert({false, name}).second) {
            std::vector<std::int64_t> ids;
            for (std::size_t b = 0; b < per.size() && b < n_blocks; ++b)
                for (std::int64_t v : per[b])
                    ids.push_back(bases[b] + v + 1);
            std::sort(ids.begin(), ids.end());
            elem_comps.emplace_back(name, std::move(ids));
        }
    if (side_regions) {
        log::warn(
            "Ansys .cdb writer: {} side region(s) have no component equivalent and were "
            "dropped",
            side_regions);
        detail::provenance_note("regions-dropped",
                                std::to_string(side_regions) + " side region(s) dropped");
    }

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open ansysInp file for writing: " + rPath);
    std::string out = detail::provenance_render_lines(detail::SlotTier::Block, "! ");
    out += "/PREP7\n";
    for (const auto& [slot, routine] : slot_routine)
        out += "ET," + std::to_string(slot) + "," + std::to_string(routine) + "\n";

    const NDArray& points = rMesh.Points();
    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    out += "NBLOCK,6,SOLID," + ans_i(static_cast<std::int64_t>(npts), 9) + "," +
           ans_i(static_cast<std::int64_t>(npts), 9) + "\n(3i9,6e21.13e3)\n";
    char buf[48];
    for (std::size_t p = 0; p < npts; ++p) {
        out += ans_i(static_cast<std::int64_t>(p + 1), 9) + ans_i(0, 9) + ans_i(0, 9);
        for (std::size_t d = 0; d < 3; ++d) {
            const double v = d < pdim ? detail::read_double(points, p * pdim + d) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), "%21.13E", v);
            out += buf;
        }
        out += '\n';
    }
    out += "N,R5.3,LOC,       -1,\n";

    const std::int64_t n_cells = bases.back();
    out += "EBLOCK,19,SOLID," + ans_i(n_cells, 10) + "," + ans_i(n_cells, 10) + "\n(19i10)\n";
    std::int64_t element = 0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const auto [routine, slot] = cell_etype[b][r];
            const std::vector<int> layout = *ans_layout(type, ans_category(routine));
            std::vector<std::int64_t> nodes;
            for (int from : layout)
                nodes.push_back(detail::read_int(conn, r * k + static_cast<std::size_t>(from)) + 1);
            const std::int64_t head[11] = {ans_cell_int(rMesh, "ansys:mat", b, r, 1),
                                           slot,
                                           ans_cell_int(rMesh, "ansys:real", b, r, 1),
                                           ans_cell_int(rMesh, "ansys:secnum", b, r, 1),
                                           0,
                                           0,
                                           0,
                                           0,
                                           static_cast<std::int64_t>(nodes.size()),
                                           0,
                                           ++element};
            for (std::int64_t v : head)
                out += ans_i(v, 10);
            for (std::size_t c = 0; c < nodes.size(); ++c) {
                if (c == 8)
                    out += '\n';
                out += ans_i(nodes[c], 10);
            }
            out += '\n';
        }
    }
    out += "        -1\n";

    // Components, with runs of consecutive ids written as `first, -last`.
    const auto write_component = [&](const std::string& rName, const char* pEntity,
                                     std::vector<std::int64_t> ids) {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        std::vector<std::int64_t> packed;
        for (std::size_t k = 0; k < ids.size();) {
            std::size_t j = k;
            while (j + 1 < ids.size() && ids[j + 1] == ids[j] + 1)
                ++j;
            packed.push_back(ids[k]);
            if (j > k)
                packed.push_back(-ids[j]);
            k = j + 1;
        }
        out += "CMBLOCK," + rName + "," + pEntity + "," +
               ans_i(static_cast<std::int64_t>(packed.size()), 8) + "\n(8i10)\n";
        for (std::size_t c = 0; c < packed.size(); ++c) {
            out += ans_i(packed[c], 10);
            if (c % 8 == 7 || c + 1 == packed.size())
                out += '\n';
        }
    };
    for (auto& [name, ids] : elem_comps)
        write_component(name, "ELEM", ids);
    for (auto& [name, ids] : node_comps)
        write_component(name, "NODE", ids);
    out += "FINISH\n";
    f << out;
}

}  // namespace meshioplusplus
