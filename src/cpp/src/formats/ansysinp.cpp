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
#include <charconv>
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
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/keyword_card_view.hpp"
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
#include "../detail/row_writer.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

std::string ans_upper(std::string_view Text) {
    std::string out(Text);
    for (char& c : out)
        if (c >= 'a' && c <= 'z')
            c = static_cast<char>(c - 'a' + 'A');
    return out;
}

std::string_view ans_strip_view(std::string_view Text) {
    const std::size_t a = Text.find_first_not_of(" \t\r\n");
    if (a == std::string_view::npos)
        return {};
    const std::size_t b = Text.find_last_not_of(" \t\r\n");
    return Text.substr(a, b - a + 1);
}

std::string ans_strip(std::string_view Text) {
    return std::string(ans_strip_view(Text));
}

// Comma-separated fields of a command line, stripped (`ET, 4, 186`).
std::vector<std::string_view> ans_commas(std::string_view Line) {
    std::vector<std::string_view> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t comma = Line.find(',', start);
        out.push_back(ans_strip_view(Line.substr(
            start, comma == std::string_view::npos ? std::string_view::npos : comma - start)));
        if (comma == std::string_view::npos)
            return out;
        start = comma + 1;
    }
}

std::optional<std::int64_t> ans_int(std::string_view rText) {
    if (rText.empty())
        return std::nullopt;
    double v = 0.0;
    if (!detail::parse_double_token(rText, v))
        return std::nullopt;
    return detail::checked_integer<std::int64_t>(v, "Ansys .cdb");
}

// The capacity a block's header count asks for, capped at one entry per line
// left: a wrong count only sizes the reservation.
std::size_t ans_count_hint(const std::vector<std::string_view>& rHeader, std::size_t Field,
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

std::int64_t ans_field_int(const std::vector<std::string_view>& rFields, std::size_t K,
                           std::size_t Line) {
    if (K >= rFields.size() || rFields[K].empty())
        return 0;
    const auto v = ans_int(rFields[K]);
    if (!v)
        ans_fail(Line, "bad integer '" + std::string(rFields[K]) + "'");
    return *v;
}

double ans_field_real(const std::vector<std::string_view>& rFields, std::size_t K,
                      std::size_t Line) {
    if (K >= rFields.size() || rFields[K].empty())
        return 0.0;
    const std::string_view text = rFields[K];
    double v = 0.0;
    if (!detail::parse_double_token(text, v))
        ans_fail(Line, "bad number '" + std::string(text) + "'");
    return v;
}

bool ans_is_terminator(std::string_view rLine) {
    const std::string s = ans_strip(rLine);
    if (s == "-1")
        return true;
    const std::string up = ans_upper(s);
    return up.rfind("N,", 0) == 0 && up.find("LOC") != std::string::npos;
}

// A command line (`FINISH`, `CMBLOCK,...`) rather than a block's data line.
bool ans_is_command(std::string_view rLine) {
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
std::vector<detail::CardField> ans_format(const std::vector<std::string_view>& rLines,
                                          std::size_t K) {
    if (K >= rLines.size())
        ans_fail(K, "a block header is not followed by its format line");
    try {
        return detail::parse_fortran_format(ans_strip(rLines[K]));
    } catch (const ReadError& e) {
        ans_fail(K, e.what());
    }
}

AnsDeck ans_parse(const std::vector<std::string_view>& rLines) {
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
                const auto f = detail::split_fixed_view(rLines[i], fields);
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
                const auto f = detail::split_fixed_view(rLines[i], fields);
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
                const auto f = detail::split_fixed_view(rLines[i], fields);
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
                    const auto more = detail::split_fixed_view(rLines[i], fields);
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
                const auto f = detail::split_fixed_view(rLines[i], fields);
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

std::vector<std::string_view> ans_read_lines(std::string_view Text) {
    auto lines = detail::split_lines(Text);
    for (auto& line : lines) {
        if (!line.empty() && line.back() == '\r')
            line.remove_suffix(1);
    }
    return lines;
}

}  // namespace

Mesh read_ansysinp(const std::string& rPath, const ReadOptions& rOptions, AnsysInfo& rInfo) {
    const detail::FileSource source =
        detail::open_source(rPath, "Could not open ansysInp file: " + rPath);
    AnsDeck deck = ans_parse(ans_read_lines(source.View()));
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
    std::vector<int> layout = detail::ansys_cell_layout(Type, Category);
    if (layout.empty())
        return std::nullopt;
    return layout;
}

std::string ans_i(std::int64_t Value, int Width) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%*lld", Width, static_cast<long long>(Value));
    return buf;
}

// `ans_i` appended in place: right-aligned in `Width` columns, never cut.
void ans_append_i(std::string& rOut, std::int64_t Value, int Width) {
    char buf[24];
    const auto r = std::to_chars(buf, buf + sizeof buf, Value);
    const auto len = static_cast<std::size_t>(r.ptr - buf);
    if (len < static_cast<std::size_t>(Width))
        rOut.append(static_cast<std::size_t>(Width) - len, ' ');
    rOut.append(buf, len);
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
    detail::AnsysElementTypes element_types =
        detail::ansys_element_types(rMesh, "Ansys .cdb writer");
    const std::map<int, int>& slot_routine = element_types.mSlotRoutine;
    const std::vector<std::vector<std::pair<int, int>>>& cell_etype = element_types.mCellEtype;

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
    // Rows formatted in parallel chunks, appended in order (row_writer.hpp).
    {
        const detail::DoubleView pv(points);
        const detail::CNumber num;
        detail::append_row_chunks(
            out, npts, [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t p = First; p < Last; ++p) {
                    ans_append_i(rBuf, static_cast<std::int64_t>(p + 1), 9);
                    ans_append_i(rBuf, 0, 9);
                    ans_append_i(rBuf, 0, 9);
                    for (std::size_t d = 0; d < 3; ++d)
                        num.Append(rBuf, "%21.13E", d < pdim ? pv[p * pdim + d] : 0.0);
                    rBuf += '\n';
                }
            });
    }
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
        // One layout per element category, not a vector per row.
        std::array<std::optional<std::vector<int>>, 8> layouts;
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            const auto category = detail::ansys_category(cell_etype[b][r].first);
            auto& slot = layouts[static_cast<std::size_t>(category)];
            if (!slot)
                slot = *ans_layout(type, category);
        }
        const detail::Int64View cv(conn);
        detail::append_row_chunks(
            out, cb.NumCells(), [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t r = First; r < Last; ++r) {
                    const auto [routine, slot] = cell_etype[b][r];
                    const std::vector<int>& layout =
                        *layouts[static_cast<std::size_t>(detail::ansys_category(routine))];
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
                    for (std::int64_t v : head)
                        ans_append_i(rBuf, v, 10);
                    for (std::size_t c = 0; c < layout.size(); ++c) {
                        if (c == 8)
                            rBuf += '\n';
                        ans_append_i(rBuf, cv[r * k + static_cast<std::size_t>(layout[c])] + 1, 10);
                    }
                    rBuf += '\n';
                }
            });
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
