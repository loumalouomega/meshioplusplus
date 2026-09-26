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
#include "meshioplusplus/detail/ansys_model.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

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
    return detail::checked_integer<std::int64_t>(v, "Ansys .cdb");
}

// The capacity a block's header count asks for, capped at one entry per line
// left: a wrong count only sizes the reservation.
std::size_t ans_count_hint(const std::vector<std::string>& rHeader, std::size_t Field,
                           std::size_t LinesLeft) {
    if (rHeader.size() <= Field)
        return 0;
    const std::int64_t v = ans_int(rHeader[Field]).value_or(0);
    return v > 0 ? std::min(static_cast<std::size_t>(v), LinesLeft) : 0;
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

struct AnsDeck {
    detail::AnsysModel mModel;
    bool mDroppedRotations = false;
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
                    deck.mModel.mRoutine[static_cast<int>(*slot)] = routine;
            }
            ++i;
        } else if (ans_is_keyopt(up)) {
            const auto p = ans_commas(line);
            if (p.size() >= 4) {
                const auto slot = ans_int(p[1]), k = ans_int(p[2]), v = ans_int(p[3]);
                if (slot && k && v)
                    deck.mModel.mKeyopt[static_cast<int>(*slot)][static_cast<int>(*k)] =
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
                    deck.mModel.mRoutine[slot] = static_cast<int>(ans_field_int(f, 1, i));
                    for (std::size_t k = 2; k < f.size() && k < 20; ++k)
                        if (const auto v = ans_int(f[k]); v && *v != 0)
                            deck.mModel.mKeyopt[slot][static_cast<int>(k - 1)] =
                                static_cast<int>(*v);
                }
                ++i;
            }
            ++i;
        } else if (up.rfind("NBLOCK", 0) == 0) {
            saw_block = true;
            // NBLOCK,<fields>,SOLID,<largest id>,<count>
            const std::size_t hint = ans_count_hint(ans_commas(up), 4, n - i);
            deck.mModel.mNodeIds.reserve(deck.mModel.mNodeIds.size() + hint);
            deck.mModel.mCoords.reserve(deck.mModel.mCoords.size() + 3 * hint);
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
                deck.mModel.mNodeIds.push_back(ans_field_int(f, 0, i));
                for (std::size_t d = 0; d < 3; ++d)
                    deck.mModel.mCoords.push_back(ans_field_real(f, n_int + d, i));
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
            // EBLOCK,<fields>,SOLID,<largest id>,<count>
            if (solid)
                deck.mModel.mElements.reserve(deck.mModel.mElements.size() +
                                              ans_count_hint(header, 4, n - i));
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
                detail::AnsysElement e;
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
                deck.mModel.mElements.push_back(std::move(e));
            }
            ++i;
        } else if (up.rfind("CMBLOCK", 0) == 0) {
            saw_block = true;
            const auto header = ans_commas(line);
            if (header.size() < 3)
                ans_fail(i, "a CMBLOCK needs a name and an entity type");
            detail::AnsysComponent comp;
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
                deck.mModel.mComponents.push_back(std::move(comp));
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

    std::unordered_map<std::int64_t, std::int64_t> node_index;
    return detail::ansys_build_mesh(deck.mModel, rOptions.mLenient, "Ansys .cdb", rInfo,
                                    node_index);
}

Mesh read_ansysinp(const std::string& rPath, AnsysInfo& rInfo) {
    return read_ansysinp(rPath, ReadOptions{}, rInfo);
}

namespace {

// The degenerate or native layout of `Type` under an element of `Category`.
std::optional<std::vector<int>> ans_layout(std::string_view Type, detail::AnsysCategory Category) {
    using V = std::vector<int>;
    if (Category == detail::AnsysCategory::Brick) {
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
    } else if (Category == detail::AnsysCategory::Tet) {
        if (Type == "tetra")
            return V{0, 1, 2, 3};
        if (Type == "tetra10")
            return V{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    } else if (Category == detail::AnsysCategory::Shell) {
        if (Type == "quad")
            return V{0, 1, 2, 3};
        if (Type == "quad8")
            return V{0, 1, 2, 3, 4, 5, 6, 7};
        if (Type == "triangle")
            return V{0, 1, 2, 2};
        if (Type == "triangle6")  // K=L, the K-L midside = K
            return V{0, 1, 2, 2, 3, 4, 2, 5};
    } else if (Category == detail::AnsysCategory::Line ||
               Category == detail::AnsysCategory::LinearLine) {
        if (Type == "line")
            return V{0, 1};
        if (Type == "line3" && Category == detail::AnsysCategory::Line)
            return V{0, 1, 2};
    } else if (Category == detail::AnsysCategory::Point) {
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

// Block `Block` of the `Name` cell data, or null when the mesh has none.
const NDArray* ans_cell_column(const Mesh& rMesh, const std::string& rName, std::size_t Block) {
    return rMesh.HasCellData(rName) ? &rMesh.CellData(rName, Block) : nullptr;
}

// Row `Row` of a cell-data column, or `Default` when absent.
std::int64_t ans_column_int(const NDArray* pColumn, std::size_t Row, std::int64_t Default) {
    if (!pColumn || pColumn->Size() <= Row)
        return Default;
    return detail::read_int(*pColumn, Row);
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
        const NDArray* elements = ans_cell_column(rMesh, "ansys:element", b);
        const NDArray* types = ans_cell_column(rMesh, "ansys:type", b);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            int routine = static_cast<int>(ans_column_int(elements, r, fallback));
            int slot = static_cast<int>(ans_column_int(types, r, 0));
            if (!ans_layout(type, detail::ansys_category(routine)) ||
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
    // Rows are formatted in parallel, then joined in order (bytes unchanged).
    std::vector<std::string> rows(npts);
    parallel_for(npts, [&](std::size_t p) {
        char buf[48];
        std::string& row = rows[p];
        row = ans_i(static_cast<std::int64_t>(p + 1), 9) + ans_i(0, 9) + ans_i(0, 9);
        for (std::size_t d = 0; d < 3; ++d) {
            const double v = d < pdim ? detail::read_double(points, p * pdim + d) : 0.0;
            detail::snprintf_c(buf, sizeof(buf), "%21.13E", v);
            row += buf;
        }
        row += '\n';
    });
    for (const std::string& row : rows)
        out += row;
    out += "N,R5.3,LOC,       -1,\n";

    const std::int64_t n_cells = bases.back();
    out += "EBLOCK,19,SOLID," + ans_i(n_cells, 10) + "," + ans_i(n_cells, 10) + "\n(19i10)\n";
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        const NDArray* mat = ans_cell_column(rMesh, "ansys:mat", b);
        const NDArray* real = ans_cell_column(rMesh, "ansys:real", b);
        const NDArray* secnum = ans_cell_column(rMesh, "ansys:secnum", b);
        rows.assign(cb.NumCells(), std::string());
        parallel_for(cb.NumCells(), [&](std::size_t r) {
            const auto [routine, slot] = cell_etype[b][r];
            const std::vector<int> layout = *ans_layout(type, detail::ansys_category(routine));
            const std::int64_t head[11] = {ans_column_int(mat, r, 1),
                                           slot,
                                           ans_column_int(real, r, 1),
                                           ans_column_int(secnum, r, 1),
                                           0,
                                           0,
                                           0,
                                           0,
                                           static_cast<std::int64_t>(layout.size()),
                                           0,
                                           bases[b] + static_cast<std::int64_t>(r) + 1};
            std::string& row = rows[r];
            for (std::int64_t v : head)
                row += ans_i(v, 10);
            for (std::size_t c = 0; c < layout.size(); ++c) {
                if (c == 8)
                    row += '\n';
                row += ans_i(
                    detail::read_int(conn, r * k + static_cast<std::size_t>(layout[c])) + 1, 10);
            }
            row += '\n';
        });
        for (const std::string& row : rows)
            out += row;
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
