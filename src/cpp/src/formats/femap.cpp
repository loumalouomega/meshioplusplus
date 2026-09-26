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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

// Femap topology code -> meshio++ type and, per meshio++ node, the slot of the
// 20-slot element record it comes from. The slots follow Femap's degenerate
// brick: corners 0-3 bottom, 4-7 top; mid-edges 8-11 bottom, 12-15 vertical,
// 16-19 top. Pinned against real Femap 8.2 files (FrontISTR's examples).
struct FnTopology {
    int mCode;
    const char* mType;
    std::vector<int> mSlots;
    int mDefaultType;  // Femap element type written for this topology
};

const std::vector<FnTopology>& fn_topologies() {
    static const std::vector<FnTopology> t = {
        {0, "line", {0, 1}, 1},
        {1, "line3", {0, 1, 2}, 1},
        {2, "triangle", {0, 1, 2}, 17},
        {3, "triangle6", {0, 1, 2, 4, 5, 6}, 18},
        {4, "quad", {0, 1, 2, 3}, 17},
        {5, "quad8", {0, 1, 2, 3, 4, 5, 6, 7}, 18},
        {6, "tetra", {0, 1, 2, 4}, 25},
        {7, "wedge", {0, 1, 2, 4, 5, 6}, 25},
        {8, "hexahedron", {0, 1, 2, 3, 4, 5, 6, 7}, 25},
        {9, "vertex", {0}, 27},
        {10, "tetra10", {0, 1, 2, 4, 8, 9, 10, 12, 13, 14}, 26},
        {11, "wedge15", {0, 1, 2, 4, 5, 6, 8, 9, 10, 16, 17, 18, 12, 13, 14}, 26},
        {12,
         "hexahedron20",
         {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15},
         26},
        {14, "pyramid", {0, 1, 2, 3, 4}, 25},
    };
    return t;
}

const FnTopology* fn_topology(std::int64_t Code) {
    for (const FnTopology& t : fn_topologies())
        if (t.mCode == Code)
            return &t;
    return nullptr;
}

const FnTopology* fn_topology_of_type(std::string_view Type) {
    for (const FnTopology& t : fn_topologies())
        if (Type == t.mType)
            return &t;
    return nullptr;
}

const char* fn_topology_name(std::int64_t Code) {
    switch (Code) {
        case 13:
            return "rigid";
        case 15:
            return "multi-list";
        case 16:
            return "contact";
        case 17:
            return "weld";
        case 18:
            return "rigid";
        case 19:
            return "pyramid13";
        default:
            return "unknown";
    }
}

std::string fn_trim(std::string_view Text) {
    const std::size_t b = Text.find_first_not_of(" \t");
    if (b == std::string_view::npos)
        return {};
    const std::size_t e = Text.find_last_not_of(" \t");
    return std::string(Text.substr(b, e - b + 1));
}

// The comma-separated fields of a record line; a trailing comma adds none.
std::vector<std::string> fn_fields(std::string_view Line) {
    std::vector<std::string> out;
    std::size_t pos = 0;
    while (pos <= Line.size()) {
        std::size_t comma = Line.find(',', pos);
        if (comma == std::string_view::npos)
            comma = Line.size();
        out.push_back(fn_trim(Line.substr(pos, comma - pos)));
        pos = comma + 1;
    }
    while (!out.empty() && out.back().empty())
        out.pop_back();
    return out;
}

bool fn_parse_int(const std::string& rText, std::int64_t& rValue) {
    if (rText.empty())
        return false;
    std::size_t i = rText[0] == '-' || rText[0] == '+' ? 1 : 0;
    if (i >= rText.size())
        return false;
    std::int64_t v = 0;
    for (; i < rText.size(); ++i) {
        if (rText[i] < '0' || rText[i] > '9')
            return false;
        if (v > (std::numeric_limits<std::int64_t>::max() - (rText[i] - '0')) / 10)
            return false;  // more digits than an int64 holds: not an id
        v = v * 10 + (rText[i] - '0');
    }
    rValue = rText[0] == '-' ? -v : v;
    return true;
}

bool fn_parse_real(const std::string& rText, double& rValue) {
    if (rText.empty())
        return false;
    const char* end = nullptr;
    rValue = detail::parse_double(rText.c_str(), end);
    return end == rText.c_str() + rText.size();
}

struct FnBlock {
    std::int64_t mId;
    std::size_t mFirstLine;  // 1-based line of the first record line
    std::vector<std::string_view> mLines;
};

// A cursor over one block's record lines.
class FnCursor {
public:
    explicit FnCursor(const FnBlock& rBlock) : mBlock(rBlock) {}

    bool AtEnd() const { return mPos >= mBlock.mLines.size(); }
    std::size_t Remaining() const { return mBlock.mLines.size() - mPos; }
    std::size_t Line() const { return mBlock.mFirstLine + mPos; }
    std::string_view Peek(std::size_t Ahead = 0) const { return mBlock.mLines[mPos + Ahead]; }

    std::string_view Next(const char* pWhat) {
        if (AtEnd())
            Fail(std::string("block ") + std::to_string(mBlock.mId) + " ends inside " + pWhat);
        return mBlock.mLines[mPos++];
    }

    std::vector<std::string> Fields(const char* pWhat) { return fn_fields(Next(pWhat)); }

    [[noreturn]] void Fail(const std::string& rWhy) const {
        throw ReadError("Femap neutral: " + rWhy + " (line " + std::to_string(Line()) + ")");
    }

    std::int64_t Int(const std::vector<std::string>& rF, std::size_t K, const char* pWhat) const {
        std::int64_t v = 0;
        if (K >= rF.size() || !fn_parse_int(rF[K], v))
            Fail(std::string("bad ") + pWhat +
                 (K < rF.size() ? " '" + rF[K] + "'" : std::string()));
        return v;
    }

    double Real(const std::vector<std::string>& rF, std::size_t K, const char* pWhat) const {
        double v = 0;
        if (K >= rF.size() || !fn_parse_real(rF[K], v))
            Fail(std::string("bad ") + pWhat +
                 (K < rF.size() ? " '" + rF[K] + "'" : std::string()));
        return v;
    }

private:
    const FnBlock& mBlock;
    std::size_t mPos = 0;
};

std::string fn_title(std::string_view Line) {
    std::string t = fn_trim(Line);
    return t == "<NULL>" ? std::string() : t;
}

struct FnElement {
    std::int64_t mId, mProperty, mType;
    const FnTopology* mTopology;
    std::array<std::int64_t, 20> mSlots;
    std::size_t mLine;
};

struct FnGroup {
    std::int64_t mId;
    std::string mTitle;
    std::vector<std::int64_t> mNodes, mElements;
};

struct FnSet {
    std::int64_t mId;
    std::string mTitle;
    double mValue;
};

struct FnVector {
    std::int64_t mSet, mId;
    std::string mTitle;
    std::int64_t mEntity;  // 7 nodal, 8 elemental
    std::vector<std::pair<std::int64_t, double>> mValues;
};

struct FnFile {
    double mVersion = 0;
    std::vector<std::int64_t> mNodeIds;
    std::vector<double> mCoords;
    std::vector<FnElement> mElements;
    std::map<std::int64_t, std::string> mProperties;
    std::vector<FnGroup> mGroups;
    std::vector<FnSet> mSets;
    std::vector<FnVector> mVectors;
};

std::vector<FnBlock> fn_blocks(std::string_view rText, std::vector<std::string_view>& rLines) {
    std::size_t pos = 0;
    while (pos < rText.size()) {
        std::size_t eol = rText.find('\n', pos);
        if (eol == std::string::npos)
            eol = rText.size();
        std::string_view line(rText.data() + pos, eol - pos);
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        rLines.push_back(line);
        pos = eol + 1;
    }
    std::vector<FnBlock> blocks;
    std::size_t i = 0;
    const std::size_t n = rLines.size();
    while (i < n) {
        if (fn_trim(rLines[i]) != "-1") {
            ++i;  // anything before the first block (MYSTRAN writes a line there)
            continue;
        }
        if (i + 1 >= n)
            break;
        std::int64_t id = 0;
        if (!fn_parse_int(fn_trim(rLines[i + 1]), id))
            throw ReadError("Femap neutral: expected a block id after '-1', found '" +
                            fn_trim(rLines[i + 1]) + "' (line " + std::to_string(i + 2) + ")");
        FnBlock block{id, i + 3, {}};
        std::size_t j = i + 2;
        while (j < n && fn_trim(rLines[j]) != "-1") {
            if (!fn_trim(rLines[j]).empty())
                block.mLines.push_back(rLines[j]);
            ++j;
        }
        if (j >= n)
            log::warn("Femap neutral: block {} (line {}) has no closing -1", id, i + 1);
        blocks.push_back(std::move(block));
        i = j + 1;
    }
    return blocks;
}

void fn_read_nodes(const FnBlock& rBlock, FnFile& rFile) {
    FnCursor c(rBlock);
    while (!c.AtEnd()) {
        const std::vector<std::string> f = c.Fields("a node");
        if (f.size() < 14)
            c.Fail("a node record with " + std::to_string(f.size()) +
                   " fields (x, y, z are 11-13)");
        rFile.mNodeIds.push_back(c.Int(f, 0, "node id"));
        for (std::size_t k = 11; k < 14; ++k)
            rFile.mCoords.push_back(c.Real(f, k, "coordinate"));
    }
}

// Skips one node list of an element record: lines up to one whose first field is -1.
void fn_skip_list(FnCursor& rC) {
    while (true) {
        const std::vector<std::string> f = rC.Fields("an element node list");
        if (!f.empty() && f[0] == "-1")
            return;
    }
}

void fn_read_elements(const FnBlock& rBlock, FnFile& rFile, std::set<std::int64_t>& rWarned) {
    FnCursor c(rBlock);
    while (!c.AtEnd()) {
        const std::size_t line = c.Line();
        const std::vector<std::string> head = c.Fields("an element");
        FnElement el{};
        el.mId = c.Int(head, 0, "element id");
        el.mProperty = c.Int(head, 2, "element property");
        el.mType = c.Int(head, 3, "element type");
        const std::int64_t topology = c.Int(head, 4, "element topology");
        el.mLine = line;
        for (int part = 0; part < 2; ++part) {
            const std::vector<std::string> f = c.Fields("an element's nodes");
            for (std::size_t k = 0; k < 10; ++k) {
                std::int64_t v = 0;
                if (k < f.size() && !fn_parse_int(f[k], v))
                    c.Fail("bad node id '" + f[k] + "'");
                el.mSlots[static_cast<std::size_t>(part) * 10 + k] = v;
            }
        }
        for (int k = 0; k < 3; ++k)
            c.Next("an element record");  // orientation, offsets
        const std::vector<std::string> last = c.Fields("an element record");
        // From 4.5, each non-zero list flag (fields 12-15) is followed by a node list.
        for (std::size_t k = 12; k < 16 && k < last.size(); ++k) {
            std::int64_t flag = 0;
            if (fn_parse_int(last[k], flag) && flag != 0 && rFile.mVersion >= 4.5 - 1e-9)
                fn_skip_list(c);
        }
        el.mTopology = fn_topology(topology);
        if (!el.mTopology) {
            if (rWarned.insert(topology).second)
                log::warn("Femap neutral: skipping {} elements (topology {}; first one is {})",
                          fn_topology_name(topology), topology, el.mId);
            continue;
        }
        rFile.mElements.push_back(el);
    }
}

void fn_read_properties(const FnBlock& rBlock, FnFile& rFile) {
    FnCursor c(rBlock);
    while (!c.AtEnd()) {
        const std::vector<std::string> head = c.Fields("a property");
        const std::int64_t id = c.Int(head, 0, "property id");
        rFile.mProperties[id] = fn_title(c.Next("a property title"));
        c.Next("property flags");
        auto skip_counted = [&](std::size_t PerLine, const char* pWhat) {
            const std::vector<std::string> f = c.Fields(pWhat);
            const std::int64_t count = c.Int(f, 0, pWhat);
            if (count < 0)
                c.Fail(std::string("negative ") + pWhat);
            for (std::int64_t k = 0; k < (count + static_cast<std::int64_t>(PerLine) - 1) /
                                             static_cast<std::int64_t>(PerLine);
                 ++k)
                c.Next(pWhat);
            return count;
        };
        skip_counted(8, "laminate count");
        const std::int64_t values = skip_counted(5, "property value count");
        // Femap 2401 follows the values with as many integers, five to a line
        // (function references): a repeat of the value count before a line of
        // several integers.
        if (c.Remaining() >= 2 && values > 0) {
            const std::vector<std::string> f = fn_fields(c.Peek());
            const std::vector<std::string> next = fn_fields(c.Peek(1));
            std::int64_t count = 0, v = 0;
            bool ints = next.size() > 1;
            for (const std::string& t : next)
                ints = ints && fn_parse_int(t, v);
            if (f.size() == 1 && fn_parse_int(f[0], count) && count == values && ints) {
                c.Next("a function count");
                for (std::int64_t k = 0; k < (values + 4) / 5; ++k)
                    c.Next("a function reference");
            }
        }
        // Outline counts (6.0 and 8.1 on): a line with one integer, then that many lines.
        while (!c.AtEnd()) {
            const std::vector<std::string> f = fn_fields(c.Peek());
            std::int64_t count = 0;
            if (f.size() != 1 || !fn_parse_int(f[0], count) || count < 0)
                break;
            c.Next("an outline count");
            for (std::int64_t k = 0; k < count; ++k)
                c.Next("an outline point");
        }
    }
}

// 408 groups, in the layout of Femap 5-9 files; a record that does not fit it
// ends the block's reading with a warning, keeping the groups read so far.
void fn_read_groups(const FnBlock& rBlock, FnFile& rFile) {
    FnCursor c(rBlock);
    try {
        while (!c.AtEnd()) {
            const std::vector<std::string> head = c.Fields("a group");
            FnGroup g{c.Int(head, 0, "group id"), fn_title(c.Next("a group title")), {}, {}};
            // layers, coordinate clipping, plane clipping and six clipping planes
            for (int k = 0; k < 3 + 18; ++k)
                c.Next("a group record");
            c.Fields("a rule count");
            // Rules: type, then start,stop,inc,include entries up to -1,-1,-1,-1;
            // the rule list ends with a lone -1.
            while (true) {
                const std::vector<std::string> f = c.Fields("a group rule");
                if (c.Int(f, 0, "rule type") == -1)
                    break;
                while (true) {
                    const std::vector<std::string> e = c.Fields("a group rule entry");
                    if (c.Int(e, 0, "rule entry") == -1)
                        break;
                }
            }
            c.Fields("a list count");
            // Lists: type (7 nodes, 8 elements), then one id per line up to -1;
            // the list of lists ends with a lone -1.
            while (true) {
                const std::vector<std::string> f = c.Fields("a group list");
                const std::int64_t type = c.Int(f, 0, "list type");
                if (type == -1)
                    break;
                while (true) {
                    const std::vector<std::string> e = c.Fields("a group list entry");
                    const std::int64_t id = c.Int(e, 0, "list entry");
                    if (id == -1)
                        break;
                    if (type == 7)
                        g.mNodes.push_back(id);
                    else if (type == 8)
                        g.mElements.push_back(id);
                }
            }
            rFile.mGroups.push_back(std::move(g));
        }
    } catch (const ReadError& e) {
        log::warn("Femap neutral: groups (block 408) past the ones read are skipped: {}", e.what());
    }
}

bool fn_is_int_line(std::string_view Line, std::size_t Min, std::size_t Max, bool Positive) {
    const std::vector<std::string> f = fn_fields(Line);
    if (f.size() < Min || f.size() > Max)
        return false;
    for (const std::string& s : f) {
        std::int64_t v = 0;
        if (!fn_parse_int(s, v))
            return false;
        if (Positive && v <= 0)
            return false;
    }
    return true;
}

bool fn_is_real_line(std::string_view Line) {
    const std::vector<std::string> f = fn_fields(Line);
    double v = 0;
    return f.size() == 1 && fn_parse_real(f[0], v);
}

// Whether an output-set record starts at `Ahead`: its id, a title, the
// program/analysis line (2-4 integers), the value, the note-line count.
bool fn_set_starts(const FnCursor& rC, std::size_t Ahead) {
    if (rC.Remaining() < Ahead + 5)
        return false;
    return fn_is_int_line(rC.Peek(Ahead), 1, 1, true) &&
           fn_is_int_line(rC.Peek(Ahead + 2), 2, 4, false) && fn_is_real_line(rC.Peek(Ahead + 3)) &&
           fn_is_int_line(rC.Peek(Ahead + 4), 1, 1, false);
}

void fn_read_sets(const FnBlock& rBlock, FnFile& rFile) {
    FnCursor c(rBlock);
    while (!c.AtEnd()) {
        if (!fn_set_starts(c, 0))
            c.Fail("expected an output set record");
        FnSet s;
        s.mId = c.Int(c.Fields("an output set"), 0, "output set id");
        s.mTitle = fn_title(c.Next("an output set title"));
        c.Next("an output set record");
        s.mValue = c.Real(c.Fields("an output set value"), 0, "output set value");
        const std::int64_t notes = c.Int(c.Fields("a note count"), 0, "note count");
        for (std::int64_t k = 0; k < notes; ++k)
            c.Next("an output set note");
        // The lines after the notes depend on the version (none up to 9, one in
        // 11.0, three from 11.2, seven in 2020.1): skip to the next record.
        while (!c.AtEnd() && !fn_set_starts(c, 0))
            c.Next("an output set record");
        rFile.mSets.push_back(std::move(s));
    }
}

void fn_read_vectors(const FnBlock& rBlock, FnFile& rFile, bool Ranges) {
    FnCursor c(rBlock);
    while (!c.AtEnd()) {
        const std::vector<std::string> head = c.Fields("an output vector");
        FnVector v;
        v.mSet = c.Int(head, 0, "output set id");
        v.mId = c.Int(head, 1, "output vector id");
        v.mTitle = fn_title(c.Next("an output vector title"));
        c.Next("an output vector range");
        c.Next("output vector components");
        c.Next("output vector components");
        std::vector<std::string> f = c.Fields("an output vector record");
        if (f.size() == 1)  // double-sided contour flag (10.0 on)
            f = c.Fields("an output vector record");
        v.mEntity = c.Int(f, 3, "output vector entity type");
        c.Next("an output vector record");
        // Data: `id,value` records (451, and some of 1051) or `start,end,value...`
        // ranges (1051), up to a line whose first field is -1.
        while (true) {
            std::vector<std::string> d = c.Fields("output vector data");
            if (!d.empty() && d[0] == "-1")
                break;
            // A 1051 record is a range when its second field is an integer (the
            // range end); Femap 11 also writes plain `id,value` records there.
            std::int64_t unused = 0;
            if (!Ranges || d.size() < 2 || !fn_parse_int(d[1], unused)) {
                v.mValues.emplace_back(c.Int(d, 0, "entity id"), c.Real(d, 1, "value"));
                continue;
            }
            const std::int64_t start = c.Int(d, 0, "range start");
            const std::int64_t stop = c.Int(d, 1, "range end");
            if (stop < start)
                c.Fail("range " + std::to_string(start) + ".." + std::to_string(stop));
            std::size_t k = 2;
            for (std::int64_t id = start; id <= stop; ++id) {
                if (k >= d.size()) {
                    d = c.Fields("output vector data");
                    k = 0;
                }
                v.mValues.emplace_back(id, c.Real(d, k++, "value"));
            }
        }
        rFile.mVectors.push_back(std::move(v));
    }
}

FnFile fn_parse(const std::string& rPath) {
    const detail::FileSource text_source =
        detail::open_source(rPath, "Femap neutral: cannot open " + rPath);
    const std::string_view text = text_source.View();
    std::vector<std::string_view> lines;
    const std::vector<FnBlock> blocks = fn_blocks(text, lines);
    FnFile f;
    for (const FnBlock& b : blocks)
        if (b.mId == 100 && b.mLines.size() >= 2) {
            FnCursor c(b);
            c.Next("a title");
            f.mVersion = c.Real(c.Fields("a version"), 0, "version");
        }
    std::set<std::int64_t> warned_topologies;
    for (const FnBlock& b : blocks) {
        switch (b.mId) {
            case 403:
                fn_read_nodes(b, f);
                break;
            case 404:
                fn_read_elements(b, f, warned_topologies);
                break;
            case 402:
                fn_read_properties(b, f);
                break;
            case 408:
                fn_read_groups(b, f);
                break;
            case 450:
                fn_read_sets(b, f);
                break;
            case 451:
                fn_read_vectors(b, f, false);
                break;
            case 1051:
                fn_read_vectors(b, f, true);
                break;
            default:
                break;
        }
    }
    return f;
}

NDArray fn_entries(const std::vector<std::int64_t>& rIds) {
    NDArray a(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

NDArray fn_scalar(DType Type, double Value) {
    NDArray a(Type, {});
    detail::write_double(a, 0, Value);
    return a;
}

}  // namespace

Mesh read_femap(const std::string& rPath, const ReadOptions& rOpts) {
    const FnFile f = fn_parse(rPath);
    if (f.mNodeIds.empty())
        throw ReadError("Femap neutral: '" + rPath +
                        "' holds no nodes (block 403): a geometry-only or results-only file has no "
                        "mesh to read");

    // --- points -----------------------------------------------------------------
    std::unordered_map<std::int64_t, std::int64_t> node_index;
    for (std::size_t p = 0; p < f.mNodeIds.size(); ++p)
        if (!node_index.emplace(f.mNodeIds[p], static_cast<std::int64_t>(p)).second)
            throw ReadError("Femap neutral: node " + std::to_string(f.mNodeIds[p]) +
                            " is defined twice");
    Mesh mesh;
    const std::size_t npts = f.mNodeIds.size();
    NDArray points(DType::Float64, {npts, 3});
    std::copy(f.mCoords.begin(), f.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // --- cells ----------------------------------------------------------------------
    std::vector<const FnTopology*> order;
    std::map<const FnTopology*, std::vector<std::size_t>> by_topology;
    for (std::size_t e = 0; e < f.mElements.size(); ++e) {
        auto [it, fresh] =
            by_topology.emplace(f.mElements[e].mTopology, std::vector<std::size_t>{});
        if (fresh)
            order.push_back(f.mElements[e].mTopology);
        it->second.push_back(e);
    }
    std::unordered_map<std::int64_t, std::int64_t> element_index;
    std::vector<std::int64_t> cell_property;
    std::vector<int> cell_dim;
    std::vector<NDArray> prop_blocks, type_blocks;
    std::vector<std::size_t> block_start{0};
    for (const FnTopology* t : order) {
        const std::vector<std::size_t>& members = by_topology[t];
        const std::size_t k = t->mSlots.size();
        const int dim = cell_type_dimension(cell_type_from_name(t->mType));
        NDArray conn(DType::Int64, {members.size(), k});
        NDArray props(DType::Int64, {members.size()});
        NDArray types(DType::Int64, {members.size()});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < members.size(); ++r) {
            const FnElement& el = f.mElements[members[r]];
            for (std::size_t j = 0; j < k; ++j) {
                const std::int64_t id = el.mSlots[static_cast<std::size_t>(t->mSlots[j])];
                const auto it = node_index.find(id);
                if (it == node_index.end())
                    throw ReadError("Femap neutral: element " + std::to_string(el.mId) +
                                    (id == 0
                                         ? " has no node in slot " + std::to_string(t->mSlots[j]) +
                                               " its topology needs"
                                         : " names undefined node " + std::to_string(id)) +
                                    " (line " + std::to_string(el.mLine) + ")");
                c[r * k + j] = it->second;
            }
            if (!element_index.emplace(el.mId, static_cast<std::int64_t>(cell_property.size()))
                     .second)
                throw ReadError("Femap neutral: element " + std::to_string(el.mId) +
                                " is defined twice (line " + std::to_string(el.mLine) + ")");
            props.As<std::int64_t>()[r] = el.mProperty;
            types.As<std::int64_t>()[r] = el.mType;
            cell_property.push_back(el.mProperty);
            cell_dim.push_back(dim);
        }
        mesh.AddCellBlock(t->mType, std::move(conn));
        prop_blocks.push_back(std::move(props));
        type_blocks.push_back(std::move(types));
        block_start.push_back(block_start.back() + members.size());
    }
    if (!prop_blocks.empty()) {
        mesh.AddCellData("femap:property", std::move(prop_blocks));
        mesh.AddCellData("femap:type", std::move(type_blocks));
    }
    const std::size_t ncells = cell_property.size();

    // --- regions: properties, then groups -------------------------------------
    std::map<std::int64_t, std::vector<std::int64_t>> by_property;
    std::map<std::int64_t, int> property_dim;
    for (std::size_t g = 0; g < ncells; ++g) {
        by_property[cell_property[g]].push_back(static_cast<std::int64_t>(g));
        auto [it, fresh] = property_dim.emplace(cell_property[g], cell_dim[g]);
        if (!fresh)
            it->second = std::max(it->second, cell_dim[g]);
    }
    for (const auto& [pid, ids] : by_property) {
        const auto t = f.mProperties.find(pid);
        const std::string name = t != f.mProperties.end() && !t->second.empty()
                                     ? t->second
                                     : "property_" + std::to_string(pid);
        mesh.AddRegion(Region(name, RegionKind::Cell, property_dim[pid], pid, fn_entries(ids)));
    }
    std::size_t missing = 0;
    for (const FnGroup& g : f.mGroups) {
        const std::string name = g.mTitle.empty() ? "group_" + std::to_string(g.mId) : g.mTitle;
        std::vector<std::int64_t> pts, cls;
        int dim = -1;
        for (std::int64_t id : g.mNodes) {
            const auto it = node_index.find(id);
            if (it == node_index.end())
                ++missing;
            else
                pts.push_back(it->second);
        }
        for (std::int64_t id : g.mElements) {
            const auto it = element_index.find(id);
            if (it == element_index.end()) {
                ++missing;
            } else {
                cls.push_back(it->second);
                dim = std::max(dim, cell_dim[static_cast<std::size_t>(it->second)]);
            }
        }
        if (!cls.empty() || pts.empty())
            mesh.AddRegion(Region(name, RegionKind::Cell, dim, g.mId, fn_entries(cls)));
        if (!pts.empty())
            mesh.AddRegion(Region(name, RegionKind::Point, -1, g.mId, fn_entries(pts)));
    }
    if (missing)
        log::warn("Femap neutral: groups name {} undefined or skipped node(s) or element(s)",
                  missing);

    // --- the selected output set -------------------------------------------------
    if (f.mSets.empty()) {
        if (!f.mVectors.empty())
            log::warn(
                "Femap neutral: output vectors without an output set (block 450) are "
                "skipped");
        rOpts.ResolveTimeStep(0);
        return mesh;
    }
    const FnSet& set = f.mSets[rOpts.ResolveTimeStep(f.mSets.size())];
    mesh.AddFieldData(kSequenceTimeKey, fn_scalar(DType::Float64, set.mValue));
    mesh.AddFieldData("femap:set", fn_scalar(DType::Int64, static_cast<double>(set.mId)));
    if (!rOpts.WantsAnyData())
        return mesh;
    std::set<std::string> used;
    std::size_t skipped = 0;
    for (const FnVector& v : f.mVectors) {
        if (v.mSet != set.mId)
            continue;
        if (v.mEntity != 7 && v.mEntity != 8) {
            ++skipped;
            continue;
        }
        std::string name = v.mTitle.empty() ? "vector_" + std::to_string(v.mId) : v.mTitle;
        if (used.count(name))
            name += " (" + std::to_string(v.mId) + ")";
        used.insert(name);
        if (!rOpts.WantsArray(name))
            continue;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        if (v.mEntity == 7) {
            NDArray data(DType::Float64, {npts});
            double* d = data.As<double>();
            std::fill(d, d + npts, nan);
            for (const auto& [id, value] : v.mValues) {
                const auto it = node_index.find(id);
                if (it != node_index.end())
                    d[it->second] = value;
            }
            mesh.AddPointData(name, std::move(data));
        } else {
            std::vector<double> per_cell(ncells, nan);
            for (const auto& [id, value] : v.mValues) {
                const auto it = element_index.find(id);
                if (it != element_index.end())
                    per_cell[static_cast<std::size_t>(it->second)] = value;
            }
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b + 1 < block_start.size(); ++b) {
                NDArray a(DType::Float64, {block_start[b + 1] - block_start[b]});
                std::copy(per_cell.begin() + static_cast<std::ptrdiff_t>(block_start[b]),
                          per_cell.begin() + static_cast<std::ptrdiff_t>(block_start[b + 1]),
                          a.As<double>());
                blocks.push_back(std::move(a));
            }
            mesh.AddCellData(name, std::move(blocks));
        }
    }
    if (skipped)
        log::warn("Femap neutral: {} output vector(s) on neither nodes nor elements skipped",
                  skipped);
    return mesh;
}

std::vector<double> femap_time_values(const std::string& rPath) {
    const FnFile f = fn_parse(rPath);
    std::vector<double> out;
    for (const FnSet& s : f.mSets)
        out.push_back(s.mValue);
    return out;
}

MeshMetadata read_femap_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    MeshMetadata meta = metadata_from_mesh(read_femap(rPath, ReadOptions{}));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "femap";
    meta.mTimeValues = femap_time_values(rPath);
    return meta;
}

// ===========================================================================
// Writer (the Femap 8.2 layout)
// ===========================================================================

namespace {

void fn_append_real(std::string& rOut, double Value) {
    char buf[40];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", Value);
    rOut += buf;
}

std::string fn_clean_title(const std::string& rTitle) {
    std::string t = rTitle;
    for (char& ch : t)
        if (ch == '\n' || ch == '\r')
            ch = ' ';
    return t.empty() ? "<NULL>" : t;
}

void fn_block_open(std::string& rOut, int Id) {
    char buf[32];
    detail::snprintf_c(buf, sizeof(buf), "   -1\n%6d\n", Id);
    rOut += buf;
}

void fn_block_close(std::string& rOut) {
    rOut += "   -1\n";
}

// Lines of `Count` zeros, `PerLine` to a line (Femap 8.2's fixed laminate and
// value arrays).
void fn_zero_lines(std::string& rOut, int Count, int PerLine, const char* pZero) {
    for (int k = 0; k < Count; k += PerLine) {
        for (int j = k; j < std::min(Count, k + PerLine); ++j) {
            rOut += pZero;
            rOut += ',';
        }
        rOut += '\n';
    }
}

}  // namespace

void write_femap(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t pdim = rMesh.PointDim();
    if (pdim > 3)
        throw WriteError("Femap neutral writer: points of dimension " + std::to_string(pdim) +
                         " (at most 3)");
    const std::size_t npts = rMesh.NumPoints();

    std::vector<const FnTopology*> tops(rMesh.NumCellBlocks(), nullptr);
    std::vector<std::size_t> block_start{0};
    std::set<std::string> dropped_types;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        tops[b] = cb.IsRagged() ? nullptr : fn_topology_of_type(cb.Type());
        if (!tops[b] && cb.NumCells())
            dropped_types.insert(std::string(cb.Type()));
        block_start.push_back(block_start.back() + cb.NumCells());
    }
    const std::size_t ncells = block_start.back();
    const bool has_prop = rMesh.HasCellData("femap:property");
    const bool has_type = rMesh.HasCellData("femap:type");
    std::vector<std::int64_t> prop(ncells, 1), etype(ncells, 0), label(ncells, 0);
    std::int64_t written = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r) {
            const std::size_t g = block_start[b] + r;
            if (has_prop)
                prop[g] = detail::read_int(rMesh.CellData("femap:property", b), r);
            etype[g] = has_type ? detail::read_int(rMesh.CellData("femap:type", b), r)
                                : (tops[b] ? tops[b]->mDefaultType : 0);
            if (tops[b])
                label[g] = ++written;
        }
    }

    // Notes first: the title renders the provenance record.
    for (const std::string& t : dropped_types) {
        log::warn("Femap neutral has no '{}' topology; those cells are dropped", t);
        detail::provenance_note("cells-dropped", "Femap neutral has no '" + t + "' topology");
    }
    std::size_t side_regions = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mKind == RegionKind::Side)
            ++side_regions;
    if (side_regions) {
        log::warn("Femap neutral has no facet groups; {} side region(s) dropped", side_regions);
        detail::provenance_note("regions-dropped", std::to_string(side_regions) +
                                                       " side region(s) have no Femap group");
    }
    // Results: one output set (450) of vectors (451), a point array per
    // component as nodal vectors, a cell array per component as elemental ones.
    struct FnOutVector {
        std::string mTitle;
        int mEntity;  // 7 nodes, 8 elements
        std::vector<std::pair<std::int64_t, double>> mValues;
    };
    std::vector<FnOutVector> vectors;
    std::vector<std::string> unwritable;
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (a.Shape().size() > 2 || a.Shape().empty() || a.Shape()[0] != npts) {
            unwritable.push_back(name);
            continue;
        }
        const std::size_t nc = a.Shape().size() == 2 ? a.Shape()[1] : 1;
        for (std::size_t c = 0; c < nc; ++c) {
            FnOutVector v{nc == 1 ? name : name + "_" + std::to_string(c), 7, {}};
            for (std::size_t p = 0; p < npts; ++p) {
                const double x = detail::read_double(a, p * nc + c);
                if (!std::isnan(x))
                    v.mValues.emplace_back(static_cast<std::int64_t>(p + 1), x);
            }
            vectors.push_back(std::move(v));
        }
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name.rfind("femap:", 0) == 0)
            continue;
        std::size_t nc = 0;
        bool ok = true;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks() && ok; ++b) {
            const NDArray& a = rMesh.CellData(name, b);
            const std::size_t c = a.Shape().size() == 2 ? a.Shape()[1] : 1;
            ok = a.Shape().size() >= 1 && a.Shape().size() <= 2 &&
                 a.Shape()[0] == rMesh.Cells(b).NumCells() && (nc == 0 || c == nc);
            nc = c;
        }
        if (!ok || nc == 0) {
            unwritable.push_back(name);
            continue;
        }
        for (std::size_t c = 0; c < nc; ++c) {
            FnOutVector v{nc == 1 ? name : name + "_" + std::to_string(c), 8, {}};
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                const NDArray& a = rMesh.CellData(name, b);
                for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r) {
                    const std::size_t g = block_start[b] + r;
                    const double x = detail::read_double(a, r * nc + c);
                    if (label[g] && !std::isnan(x))
                        v.mValues.emplace_back(label[g], x);
                }
            }
            vectors.push_back(std::move(v));
        }
    }
    std::int64_t set_id = 1;
    double set_value = 0.0;
    for (const std::string& name : rMesh.FieldDataNames()) {
        const NDArray& a = rMesh.FieldData(name);
        if (name == "femap:set" && a.Size() == 1)
            set_id = std::max<std::int64_t>(1, detail::read_int(a, 0));
        else if (name == kSequenceTimeKey && a.Size() == 1)
            set_value = detail::read_double(a, 0);
        else
            unwritable.push_back(name);
    }
    if (!unwritable.empty()) {
        std::string list;
        for (const std::string& n : unwritable)
            list += (list.empty() ? "" : ", ") + n;
        log::warn("Femap neutral writer: arrays with no Femap output vector dropped: {}", list);
        detail::provenance_note("data-dropped", "arrays with no Femap output vector: " + list);
    }

    // Properties: the cells of each, and a title from the cell region the reader
    // made of it (same tag, same cells).
    std::map<std::int64_t, std::vector<std::int64_t>> prop_cells;
    std::map<std::int64_t, std::int64_t> prop_type;
    for (std::size_t g = 0; g < ncells; ++g)
        if (label[g]) {
            prop_cells[prop[g]].push_back(static_cast<std::int64_t>(g));
            prop_type.emplace(prop[g], etype[g]);
        }
    std::map<std::int64_t, std::string> prop_title;
    std::set<std::size_t> property_regions;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const auto it = prop_cells.find(reg.mTag);
        if (reg.mKind != RegionKind::Cell || it == prop_cells.end() || prop_title.count(reg.mTag))
            continue;
        const std::int64_t* e = reg.Entries();
        if (std::vector<std::int64_t>(e, e + reg.NumEntries()) == it->second) {
            prop_title[reg.mTag] = reg.mName;
            property_regions.insert(r);
        }
    }

    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("Could not open file for writing: " + rPath);
    std::string out;
    fn_block_open(out, 100);
    out += fn_clean_title(detail::provenance_lines(detail::SlotTier::SingleLine)[0]) + "\n8.2,\n";
    fn_block_close(out);

    if (!prop_cells.empty()) {
        fn_block_open(out, 402);
        for (const auto& [pid, cells] : prop_cells) {
            const auto t = prop_title.find(pid);
            out += std::to_string(pid) + ",110,0," + std::to_string(prop_type[pid]) + ",1,0,\n";
            out += fn_clean_title(t != prop_title.end() ? t->second
                                                        : "property_" + std::to_string(pid)) +
                   "\n0,0,0,0,\n90,\n";
            fn_zero_lines(out, 90, 8, "0");
            out += "190,\n";
            fn_zero_lines(out, 190, 5, "0.");
            out += "0,\n0,\n";
        }
        fn_block_close(out);
    }

    fn_block_open(out, 403);
    const NDArray& points = rMesh.Points();
    for (std::size_t p = 0; p < npts; ++p) {
        out += std::to_string(p + 1) + ",0,0,1,46,0,0,0,0,0,0,";
        for (std::size_t d = 0; d < 3; ++d) {
            fn_append_real(out, d < pdim ? detail::read_double(points, p * pdim + d) : 0.0);
            out += ',';
        }
        out += "0,\n";
    }
    fn_block_close(out);
    f << out;
    out.clear();

    if (written) {
        fn_block_open(out, 404);
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const FnTopology* t = tops[b];
            if (!t)
                continue;
            const auto cb = rMesh.Cells(b);
            const NDArray& conn = cb.Conn();
            const std::size_t k = t->mSlots.size();
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::size_t g = block_start[b] + r;
                std::array<std::int64_t, 20> slots{};
                for (std::size_t j = 0; j < k; ++j)
                    slots[static_cast<std::size_t>(t->mSlots[j])] =
                        detail::read_int(conn, r * k + j) + 1;
                out += std::to_string(label[g]) + ",124," + std::to_string(prop[g]) + "," +
                       std::to_string(etype[g]) + "," + std::to_string(t->mCode) +
                       ",1,0,0,0,0,0,0,0,\n";
                for (std::size_t s = 0; s < 20; ++s) {
                    out += std::to_string(slots[s]) + ",";
                    if (s == 9 || s == 19)
                        out += '\n';
                }
                out += "0.,0.,0.,\n0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n";
            }
            if (out.size() > (1u << 20)) {
                f << out;
                out.clear();
            }
        }
        fn_block_close(out);
    }

    // Groups: one per region name (a point and a cell region of one name merge),
    // the property regions excepted.
    std::vector<std::string> names;
    std::map<std::string, std::pair<std::vector<std::int64_t>, std::vector<std::int64_t>>> groups;
    std::map<std::string, std::int64_t> tags;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Side || property_regions.count(r))
            continue;
        auto [it, fresh] = groups.emplace(
            reg.mName, std::make_pair(std::vector<std::int64_t>{}, std::vector<std::int64_t>{}));
        if (fresh) {
            names.push_back(reg.mName);
            tags[reg.mName] = reg.mTag;
        }
        const std::int64_t* e = reg.Entries();
        for (std::size_t j = 0; j < reg.NumEntries(); ++j) {
            if (reg.mKind == RegionKind::Point) {
                it->second.first.push_back(e[j] + 1);
            } else {
                const std::size_t g = static_cast<std::size_t>(e[j]);
                if (g < ncells && label[g])
                    it->second.second.push_back(label[g]);
            }
        }
    }
    if (!names.empty()) {
        std::set<std::int64_t> used;
        for (const std::string& n : names)
            if (tags[n] > 0)
                used.insert(tags[n]);
        std::set<std::int64_t> assigned;
        std::int64_t next = 1;
        fn_block_open(out, 408);
        for (const std::string& name : names) {
            std::int64_t id = tags[name];
            if (id <= 0 || assigned.count(id)) {
                while (used.count(next) || assigned.count(next))
                    ++next;
                id = next;
            }
            assigned.insert(id);
            out += std::to_string(id) + ",0,0,\n" + fn_clean_title(name) + "\n";
            out += "0,0,0,\n0,0,0,0,0.,0.,\n0,0,\n";
            for (int k = 0; k < 6; ++k)
                out += "0,0,\n0.,0.,0.,\n0.,0.,0.,\n";
            out += "91,\n-1,\n23,\n";
            const auto& [nodes, elements] = groups[name];
            if (!nodes.empty()) {
                out += "7,\n";
                for (std::int64_t id2 : nodes)
                    out += std::to_string(id2) + ",\n";
                out += "-1,\n";
            }
            if (!elements.empty()) {
                out += "8,\n";
                for (std::int64_t id2 : elements)
                    out += std::to_string(id2) + ",\n";
                out += "-1,\n";
            }
            out += "-1,\n";
        }
        fn_block_close(out);
    }

    if (!vectors.empty()) {
        fn_block_open(out, 450);
        out += std::to_string(set_id) + ",\nmeshio++\n0,1,\n";
        fn_append_real(out, set_value);
        out += ",\n0,\n";
        fn_block_close(out);
        fn_block_open(out, 451);
        for (std::size_t k = 0; k < vectors.size(); ++k) {
            const FnOutVector& v = vectors[k];
            double lo = 0.0, hi = 0.0, absmax = 0.0;
            std::int64_t id_lo = 0, id_hi = 0;
            for (std::size_t j = 0; j < v.mValues.size(); ++j) {
                const auto& [id, x] = v.mValues[j];
                if (j == 0 || x < lo) {
                    lo = x;
                    id_lo = id;
                }
                if (j == 0 || x > hi) {
                    hi = x;
                    id_hi = id;
                }
                absmax = std::max(absmax, std::fabs(x));
            }
            out += std::to_string(set_id) + "," + std::to_string(k + 1) + ",1,\n" +
                   fn_clean_title(v.mTitle) + "\n";
            fn_append_real(out, lo);
            out += ',';
            fn_append_real(out, hi);
            out += ',';
            fn_append_real(out, absmax);
            out += ",\n0,0,0,0,0,0,0,0,0,0,\n0,0,0,0,0,0,0,0,0,0,\n";
            out += std::to_string(id_lo) + "," + std::to_string(id_hi) + ",0," +
                   std::to_string(v.mEntity) + ",\n0,0,1,\n";
            for (const auto& [id, x] : v.mValues) {
                out += std::to_string(id) + ',';
                fn_append_real(out, x);
                out += ",\n";
            }
            out += "-1,0.,\n";
            if (out.size() > (1u << 20)) {
                f << out;
                out.clear();
            }
        }
        fn_block_close(out);
    }
    f << out;
    if (!f)
        throw WriteError("Femap neutral writer: failed writing " + rPath);
}

}  // namespace meshioplusplus
