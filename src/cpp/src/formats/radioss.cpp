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
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

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
    const detail::FileSource text_source =
        detail::open_source(rPath.string(), "Radioss: cannot open " + rPath.string());
    const std::string_view text = text_source.View();
    const std::string label = rPath.filename().string();
    std::size_t pos = 0, number = 0;
    while (pos < text.size() && !rEnded) {
        std::size_t eol = text.find('\n', pos);
        if (eol == std::string::npos)
            eol = text.size();
        std::string line(text.substr(pos, eol - pos));
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
    std::int64_t mId = 0;
    std::string mTitle;
    std::vector<std::array<std::int64_t, 4>> mSegments;
    // `/SURF/<mKind>[/<mMode>]/id`: SEG (the segments above), PART, SUBSET, MAT,
    // PROP, GRBRIC, GRSHEL, GRSH3N, GRTRIA or SURF, with EXT, ALL or FREE.
    std::string mKind = "SEG", mMode;
    std::vector<std::int64_t> mIds;
};

// `/BOX/RECTA`, `/CYLIN`, `/SPHER` (corners, axis ends or centre, from nodes or
// coordinates) and `/BOX/BOX` (other boxes, a negative id subtracted).
struct RadBox {
    std::string mKind;
    std::int64_t mSkew = 0;
    std::int64_t mNode1 = 0, mNode2 = 0;
    std::array<double, 3> mP1{}, mP2{};
    double mDiameter = 0.0;
    std::vector<std::int64_t> mChildren;
};

// `/SKEW/FIX`: a fixed frame, its origin and the X and Y axes as given.
struct RadSkew {
    std::array<double, 3> mOrigin{}, mX{}, mY{};
};

// The unit axes of a fixed skew (X, Z = X x Y, Y = Z x X), as the rows of a
// rotation into the skew frame; false when X and Y are parallel or zero.
bool rad_skew_axes(const RadSkew& rSkew, std::array<std::array<double, 3>, 3>& rAxes) {
    const auto cross = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return std::array<double, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                                     a[0] * b[1] - a[1] * b[0]};
    };
    const auto unit = [](std::array<double, 3> a) {
        const double n = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if (n > 0.0)
            for (double& x : a)
                x /= n;
        return std::make_pair(a, n > 0.0);
    };
    const auto [x, okx] = unit(rSkew.mX);
    const auto [z, okz] = unit(cross(x, rSkew.mY));
    if (!okx || !okz)
        return false;
    rAxes = {x, cross(z, x), z};
    return true;
}

struct RadSubset {
    std::string mTitle;
    std::vector<std::int64_t> mChildren;
};

// A /BEGIN length unit in metres: an SI prefix and "m" ("mm", "mum", "km"), a
// few imperial names, or a number; 0 when unknown.
double rad_length_unit(const std::string& rUnit) {
    if (rUnit.empty())
        return 0.0;
    static const std::pair<const char*, double> named[] = {
        {"in", 0.0254}, {"ft", 0.3048}, {"yd", 0.9144}, {"mi", 1609.344}};
    for (const auto& [name, metres] : named)
        if (rUnit == name)
            return metres;
    if (rUnit.back() == 'm') {
        static const std::pair<const char*, double> prefixes[] = {
            {"", 1.0},   {"y", 1e-24}, {"z", 1e-21}, {"a", 1e-18}, {"f", 1e-15}, {"p", 1e-12},
            {"n", 1e-9}, {"mu", 1e-6}, {"m", 1e-3},  {"c", 1e-2},  {"d", 1e-1},  {"da", 1e1},
            {"h", 1e2},  {"k", 1e3},   {"M", 1e6},   {"G", 1e9},   {"T", 1e12},  {"P", 1e15},
            {"E", 1e18}, {"Z", 1e21},  {"Y", 1e24}};
        const std::string prefix = rUnit.substr(0, rUnit.size() - 1);
        for (const auto& [name, factor] : prefixes)
            if (prefix == name)
                return factor;
    }
    try {
        return detail::card_to_real(rUnit, "", "Radioss");
    } catch (const ReadError&) {
        return 0.0;
    }
}

// A fixed-width slice of a line (the whole comma field when the line has commas).
std::string rad_slice(const std::string& rLine, std::size_t At, std::size_t Width,
                      std::size_t CommaField) {
    if (rLine.find(',') != std::string::npos) {
        const std::vector<std::string> f = rad_fields(rLine, 1, 0);
        return CommaField < f.size() ? f[CommaField] : std::string();
    }
    return At < rLine.size() ? rad_trim(std::string_view(rLine).substr(At, Width)) : std::string();
}

NDArray rad_ids(const std::vector<std::int64_t>& rIds, std::size_t Stride = 1) {
    NDArray a(DType::Int64, Stride == 1 ? std::vector<std::size_t>{rIds.size()}
                                        : std::vector<std::size_t>{rIds.size() / Stride, Stride});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

// The input version on a `#RADIOSS STARTER` line (`41` in
// `#RADIOSS STARTER      41BAR2V41B`), 0 if none.
int rad_header_version(const std::string& rHead) {
    const std::string up = rad_upper(rHead);
    const std::size_t at = up.find("#RADIOSS STARTER");
    if (at == std::string::npos)
        return 0;
    std::size_t k = at + 16;
    const std::size_t eol = up.find('\n', k);
    while (k < up.size() && k < eol && (up[k] == ' ' || up[k] == '\t'))
        ++k;
    int v = 0;
    std::size_t digits = 0;
    for (; k < up.size() && k < eol && up[k] >= '0' && up[k] <= '9' && digits < 4; ++k, ++digits)
        v = v * 10 + (up[k] - '0');
    return v;
}

// An engine deck (`<run>_0001.rad`): every keyword and the numbers of its
// lines, as `radioss:engine:<keyword>` field data (`/RUN/<name>/1` holds the
// end time, `/ANIM/DT` start and interval, `/TFILE` the history interval; an
// output request such as `/ANIM/ELEM/SIGX` is an empty array).
std::vector<std::pair<std::string, std::vector<double>>> rad_engine_fields(
    const std::string& rPath) {
    std::vector<RadLine> lines;
    bool ended = false;
    rad_collect(fs::path(rPath), 0, lines, ended);
    std::vector<std::pair<std::string, std::vector<double>>> out;
    for (const RadLine& ln : lines) {
        const std::string t = rad_trim(ln.mText);
        if (t.empty())
            continue;
        if (t[0] == '/') {
            out.emplace_back("radioss:engine:" + t.substr(1), std::vector<double>{});
            continue;
        }
        if (out.empty())
            continue;
        auto iss = detail::make_classic_istringstream(t);
        std::string tok;
        while (iss >> tok) {
            const char* e = nullptr;
            const double v = detail::parse_double(tok.c_str(), e);
            if (e == tok.c_str() + tok.size())
                out.back().second.push_back(v);
        }
    }
    return out;
}

void rad_add_engine_fields(Mesh& rMesh, const std::string& rPath) {
    for (const auto& [name, values] : rad_engine_fields(rPath)) {
        std::vector<double> all = values;
        if (rMesh.HasFieldData(name)) {  // a repeated keyword: its numbers appended
            const NDArray& prev = rMesh.FieldData(name);
            std::vector<double> joined(prev.Size());
            for (std::size_t k = 0; k < joined.size(); ++k)
                joined[k] = detail::read_double(prev, k);
            joined.insert(joined.end(), all.begin(), all.end());
            all = std::move(joined);
        }
        NDArray a(DType::Float64, {all.size()});
        std::copy(all.begin(), all.end(), a.As<double>());
        rMesh.AddFieldData(name, std::move(a));
    }
}

// An engine deck read on its own: no mesh, its run controls.
Mesh rad_engine_mesh(const std::string& rPath) {
    Mesh mesh;
    mesh.AssignPoints(NDArray(DType::Float64, {0, 3}));
    rad_add_engine_fields(mesh, rPath);
    return mesh;
}

}  // namespace

Mesh read_radioss(const std::string& rPath) {
    std::vector<RadLine> lines;
    bool ended = false;
    int header_version = 0;
    {
        auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
        if (!in)
            throw ReadError("Radioss: cannot open " + rPath);
        std::string head(512, '\0');
        in.read(head.data(), static_cast<std::streamsize>(head.size()));
        head.resize(static_cast<std::size_t>(in.gcount()));
        if (rad_upper(head).find("#RADIOSS ENGINE") != std::string::npos)
            return rad_engine_mesh(rPath);
        header_version = rad_header_version(head);
    }
    rad_collect(fs::path(rPath), 0, lines, ended);

    int version = 2019;
    int iw = 10, rw = 20;
    // Before input version 5.1 (4.1, 4.4) there is no /BEGIN: the version is
    // on the #RADIOSS STARTER line, fields are 8 and 16 columns wide, and a
    // title is the keyword's last part rather than a line of its own.
    bool titles_in_path = header_version > 0 && header_version < 51;
    if (titles_in_path) {
        version = header_version;
        iw = 8;
        rw = 16;
    }
    std::vector<std::int64_t> node_ids;
    std::vector<double> coords;
    std::vector<RadElement> elements;
    std::map<std::int64_t, RadPart> parts;
    std::vector<std::int64_t> part_order;
    std::vector<RadGroup> groups;
    std::vector<RadSurface> surfaces;
    std::map<std::int64_t, RadSubset> subsets;
    std::map<std::int64_t, RadBox> boxes;
    std::map<std::int64_t, RadSkew> skews;
    std::vector<std::pair<std::string, std::vector<double>>> analytic;  // /SURF/PLANE, /ELLIPS
    std::map<std::int64_t, std::array<double, 3>> ellipsoid_skew;       // surf -> skew, n, -
    double length_scale = 1.0;
    double work_length = 0.0;  // the work length unit in metres, 0 if /BEGIN names none
    // /UNIT/<id>: the local length units in metres, read first (a keyword can
    // name a unit defined further down).
    std::map<std::int64_t, double> unit_length;
    for (std::size_t k = 0; k + 2 < lines.size(); ++k) {
        const std::string t = rad_upper(rad_trim(lines[k].mText));
        if (t.rfind("/UNIT/", 0) != 0)
            continue;
        const std::string id = rad_trim(std::string_view(t).substr(6));
        if (id.empty() || id.find_first_not_of("0123456789") != std::string::npos)
            continue;
        const std::vector<std::string> f = rad_fields(lines[k + 2].mText, 20, 3);
        const double len = f.size() > 1 ? rad_length_unit(f[1]) : 0.0;
        if (len > 0.0)
            unit_length[std::stoll(id)] = len;
        else
            log::warn("Radioss: /UNIT/{} has no length unit meshio++ knows; ignored", id);
    }
    std::set<std::int64_t> warned_units;
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
        std::vector<std::string> path, raw_path;
        {
            std::string kw = rad_trim(head.mText);
            std::size_t start = 1;
            while (start <= kw.size()) {
                const std::size_t k = kw.find('/', start);
                raw_path.push_back(rad_trim(std::string_view(kw).substr(
                    start, k == std::string::npos ? std::string::npos : k - start)));
                path.push_back(rad_upper(raw_path.back()));
                if (k == std::string::npos)
                    break;
                start = k + 1;
            }
        }
        const std::string& key = path.empty() ? std::string() : path[0];
        const std::size_t body = i + 1;
        const std::size_t end = block_end(body);
        // The keyword's id: its first integer field (`/BOX/RECTA/3/1` is box 3 in
        // unit system 1; `/SURF/PART/EXT/12` is surface 12).
        auto last_id = [&](const RadLine& rLine) -> std::int64_t {
            for (std::size_t k = 1; k < path.size(); ++k)
                if (!path[k].empty() &&
                    path[k].find_first_not_of("+-0123456789") == std::string::npos)
                    return rad_int(path[k], rLine);
            return path.size() >= 2 ? rad_int(path.back(), rLine) : 0;
        };
        // The unit system a keyword names: the integer after its option id
        // (`/BOX/RECTA/3/1`), or its only integer when it has none (`/NODE/1`).
        auto unit_of = [&](bool OptionId) -> std::int64_t {
            std::vector<std::int64_t> ints;
            for (std::size_t k = 1; k < path.size(); ++k)
                if (!path[k].empty() &&
                    path[k].find_first_not_of("+-0123456789") == std::string::npos)
                    ints.push_back(rad_int(path[k], head));
            const std::size_t at = OptionId ? 1 : 0;
            return ints.size() > at ? ints[at] : 0;
        };
        // Lengths of this keyword into the work units: its /UNIT's, else /BEGIN's.
        auto scale_of = [&](std::int64_t Unit) {
            if (Unit == 0)
                return length_scale;
            const auto it = unit_length.find(Unit);
            if (it == unit_length.end() || work_length <= 0.0) {
                if (warned_units.insert(Unit).second)
                    log::warn(
                        "Radioss: unit system {} is {}; its lengths are read in the input "
                        "units",
                        Unit,
                        it == unit_length.end() ? "not defined" : "used without /BEGIN units");
                return length_scale;
            }
            return it->second / work_length;
        };
        // A keyword's title: its own line from input version 5.1 on; before,
        // the keyword's last part (`/PART/1/CUIVRE`), the data then starting on
        // the next line.
        auto title_of = [&](std::size_t& rK) -> std::string {
            if (titles_in_path) {
                const std::string& last = raw_path.back();
                return raw_path.size() >= 3 &&
                               last.find_first_not_of("+-0123456789") != std::string::npos
                           ? last
                           : std::string();
            }
            return rK < end ? rad_trim(lines[rK++].mText) : std::string();
        };

        if (key == "BEGIN") {
            // run name; Invers Irun; two unit lines. A deck with /BEGIN has its
            // titles on lines of their own, whatever its input version.
            titles_in_path = false;
            if (body + 1 < end) {
                const std::vector<std::string> f = rad_fields(lines[body + 1].mText, 10, 2);
                if (!f.empty() && !f[0].empty())
                    version = static_cast<int>(rad_int(f[0], lines[body + 1]));
            }
            iw = version >= 51 ? 10 : 8;
            rw = version >= 51 ? 20 : 16;
            // Input and work units (mass, length, time; 20 columns each): the
            // solver works in the work units, so lengths are converted.
            if (body + 2 < end) {
                const std::vector<std::string> in = rad_fields(lines[body + 2].mText, 20, 3);
                const std::vector<std::string> work = body + 3 < end
                                                          ? rad_fields(lines[body + 3].mText, 20, 3)
                                                          : std::vector<std::string>{};
                const std::string li = in.size() > 1 ? in[1] : std::string();
                const std::string lw = work.size() > 1 && !work[1].empty() ? work[1] : li;
                const double fi = rad_length_unit(li), fw = rad_length_unit(lw);
                work_length = fw;
                if (!li.empty() && (fi <= 0.0 || fw <= 0.0))
                    log::warn(
                        "Radioss: unknown length unit '{}' in /BEGIN; lengths are read "
                        "as written",
                        fi <= 0.0 ? li : lw);
                else if (!li.empty())
                    length_scale = fi / fw;
            }
        } else if (key == "BOX" && path.size() >= 3) {
            RadBox b;
            b.mKind = path[1];
            std::size_t k = body;
            title_of(k);
            const double box_scale = scale_of(unit_of(true));
            auto real3 = [&](std::size_t Line) {
                std::array<double, 3> p{};
                if (Line >= end)
                    return p;
                const std::vector<std::string> f = rad_fields(lines[Line].mText, rw, 3);
                for (std::size_t d = 0; d < 3; ++d)
                    p[d] = d < f.size() && !f[d].empty() ? rad_real(f[d], lines[Line]) * box_scale
                                                         : 0.0;
                return p;
            };
            auto int_at = [&](std::size_t Line, std::size_t Field) -> std::int64_t {
                if (Line >= end)
                    return 0;
                const std::string t =
                    rad_slice(lines[Line].mText, Field * static_cast<std::size_t>(iw),
                              static_cast<std::size_t>(iw), Field);
                return t.empty() ? 0 : rad_int(t, lines[Line]);
            };
            auto diameter = [&](std::size_t Line) {
                if (Line >= end)
                    return 0.0;
                const std::string t =
                    rad_slice(lines[Line].mText, 3 * static_cast<std::size_t>(iw),
                              static_cast<std::size_t>(rw), b.mKind == "SPHER" ? 2 : 3);
                return t.empty() ? 0.0 : rad_real(t, lines[Line]) * box_scale;
            };
            if (b.mKind == "RECTA") {
                b.mNode1 = int_at(k, 0);
                b.mNode2 = int_at(k, 1);
                b.mSkew = int_at(k, 2);
                b.mP1 = real3(k + 1);
                b.mP2 = real3(k + 2);
            } else if (b.mKind == "CYLIN") {
                b.mNode1 = int_at(k, 0);
                b.mNode2 = int_at(k, 1);
                b.mDiameter = diameter(k);
                b.mP1 = real3(k + 1);
                b.mP2 = real3(k + 2);
            } else if (b.mKind == "SPHER") {
                b.mNode1 = int_at(k, 0);
                b.mDiameter = diameter(k);
                b.mP1 = real3(k + 1);
            } else if (b.mKind == "BOX") {
                for (; k < end; ++k)
                    for (const std::string& f : rad_fields(lines[k].mText, iw, 10))
                        if (!f.empty())
                            b.mChildren.push_back(rad_int(f, lines[k]));
            } else {
                skipped_keywords.insert("/BOX/" + b.mKind);
            }
            boxes[last_id(head)] = std::move(b);
        } else if (key == "SKEW" && path.size() >= 3 && path[1] == "FIX") {
            // origin (from 5.1), X axis, Y axis
            RadSkew sk;
            std::size_t k = body;
            title_of(k);
            const double sk_scale = scale_of(unit_of(true));
            auto vec = [&](std::size_t Line, double Scale) {
                std::array<double, 3> v{};
                if (Line >= end)
                    return v;
                const std::vector<std::string> f = rad_fields(lines[Line].mText, rw, 3);
                for (std::size_t d = 0; d < 3; ++d)
                    v[d] =
                        d < f.size() && !f[d].empty() ? rad_real(f[d], lines[Line]) * Scale : 0.0;
                return v;
            };
            if (!titles_in_path)  // 4.x skews have no origin line
                sk.mOrigin = vec(k++, sk_scale);
            sk.mX = vec(k, 1.0);
            sk.mY = vec(k + 1, 1.0);
            skews[last_id(head)] = sk;
        } else if (key == "SURF" && path.size() >= 3 &&
                   (path[1] == "PLANE" || path[1] == "ELLIPS")) {
            // Analytical surfaces: no segments, their definition as field data.
            const std::int64_t id = last_id(head);
            std::size_t k = body;
            title_of(k);
            const double sc = scale_of(unit_of(true));
            auto reals = [&](std::size_t Line, std::size_t Count) {
                std::vector<double> v(Count, 0.0);
                if (Line >= end)
                    return v;
                const std::vector<std::string> f = rad_fields(lines[Line].mText, rw, Count);
                for (std::size_t d = 0; d < Count; ++d)
                    v[d] = d < f.size() && !f[d].empty() ? rad_real(f[d], lines[Line]) * sc : 0.0;
                return v;
            };
            if (path[1] == "PLANE") {
                std::vector<double> v = reals(k, 3);
                const std::vector<double> m1 = reals(k + 1, 3);
                v.insert(v.end(), m1.begin(), m1.end());
                analytic.emplace_back("radioss:surf_plane:" + std::to_string(id), std::move(v));
            } else {
                const std::vector<std::string> f =
                    k < end ? rad_fields(lines[k].mText, iw, 2) : std::vector<std::string>{};
                const std::int64_t skew = !f.empty() && !f[0].empty() ? rad_int(f[0], lines[k]) : 0;
                const std::int64_t degree =
                    f.size() > 1 && !f[1].empty() ? rad_int(f[1], lines[k]) : 2;
                std::vector<double> v = {static_cast<double>(degree < 2 ? 2 : degree)};
                const std::vector<double> centre = reals(k + 1, 3), axes = reals(k + 2, 3);
                v.insert(v.end(), centre.begin(), centre.end());
                v.insert(v.end(), axes.begin(), axes.end());
                ellipsoid_skew[id] = {static_cast<double>(skew), 0.0, 0.0};
                analytic.emplace_back("radioss:surf_ellips:" + std::to_string(id), std::move(v));
            }
        } else if (key == "NODE") {
            const double node_scale = scale_of(unit_of(false));
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
                    coords.push_back(d < f.size() && !f[d].empty() ? rad_real(f[d], ln) * node_scale
                                                                   : 0.0);
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
            part.mTitle = title_of(k);
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
            g.mTitle = title_of(k);
            for (; k < end; ++k)
                for (const std::string& f : rad_fields(lines[k].mText, iw, 10))
                    if (!f.empty())
                        g.mIds.push_back(rad_int(f, lines[k]));
            groups.push_back(std::move(g));
        } else if (key == "SURF" && path.size() >= 3 && path[1] == "SEG") {
            RadSurface s;
            s.mId = last_id(head);
            std::size_t k = body;
            s.mTitle = title_of(k);
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
        } else if (key == "SURF" && path.size() >= 3 &&
                   (path[1] == "PART" || path[1] == "SUBSET" || path[1] == "MAT" ||
                    path[1] == "PROP" || path[1] == "GRBRIC" || path[1] == "GRSHEL" ||
                    path[1] == "GRSH3N" || path[1] == "GRTRIA" || path[1] == "SURF" ||
                    path[1] == "BOX" || path[1] == "BOX2")) {
            RadSurface s;
            s.mId = last_id(head);
            s.mKind = path[1];
            if (path.size() >= 4 && (path[2] == "EXT" || path[2] == "ALL" || path[2] == "FREE"))
                s.mMode = path[2];
            std::size_t k = body;
            s.mTitle = title_of(k);
            for (; k < end; ++k)
                for (const std::string& f : rad_fields(lines[k].mText, iw, 10))
                    if (!f.empty())
                        s.mIds.push_back(rad_int(f, lines[k]));
            surfaces.push_back(std::move(s));
        } else if (key == "SURF") {
            skipped_keywords.insert("/SURF/" + (path.size() > 1 ? path[1] : std::string()));
        } else if (key == "SUBSET") {
            RadSubset s;
            std::size_t k = body;
            s.mTitle = title_of(k);
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
    if (length_scale != 1.0)
        mesh.AddFieldData("radioss:length_scale", [&] {
            NDArray a(DType::Float64, {});
            a.As<double>()[0] = length_scale;
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
    // Box membership of a point (skewed boxes are not supported: empty).
    auto point_of = [&](std::int64_t Node, const std::array<double, 3>& rFallback) {
        const auto it = Node ? node_index.find(Node) : node_index.end();
        if (it == node_index.end())
            return rFallback;
        const std::size_t p = static_cast<std::size_t>(it->second);
        return std::array<double, 3>{coords[3 * p], coords[3 * p + 1], coords[3 * p + 2]};
    };
    std::set<std::int64_t> skewed_boxes, missing_boxes;
    std::function<bool(std::int64_t, const double*, int)> in_box;
    in_box = [&](std::int64_t Id, const double* pX, int Depth) -> bool {
        const auto it = boxes.find(Id);
        if (it == boxes.end() || Depth > 16) {
            missing_boxes.insert(Id);
            return false;
        }
        const RadBox& b = it->second;
        if (b.mKind == "BOX") {
            bool in = false;
            for (std::int64_t c : b.mChildren)
                if (c > 0 && in_box(c, pX, Depth + 1))
                    in = true;
            for (std::int64_t c : b.mChildren)
                if (c < 0 && in_box(-c, pX, Depth + 1))
                    in = false;
            return in;
        }
        const std::array<double, 3> p1 = point_of(b.mNode1, b.mP1);
        if (b.mSkew && b.mKind == "RECTA") {
            // Edges along the skew's axes: compare in the skew frame.
            const auto sk = skews.find(b.mSkew);
            std::array<std::array<double, 3>, 3> axes;
            if (sk == skews.end() || !rad_skew_axes(sk->second, axes)) {
                skewed_boxes.insert(Id);
                return false;
            }
            const std::array<double, 3> p2 = point_of(b.mNode2, b.mP2);
            for (const auto& a : axes) {
                const double x = a[0] * pX[0] + a[1] * pX[1] + a[2] * pX[2];
                const double u = a[0] * p1[0] + a[1] * p1[1] + a[2] * p1[2];
                const double v = a[0] * p2[0] + a[1] * p2[1] + a[2] * p2[2];
                if (x < std::min(u, v) || x > std::max(u, v))
                    return false;
            }
            return true;
        }
        if (b.mKind == "SPHER") {
            double d2 = 0.0;
            for (std::size_t d = 0; d < 3; ++d)
                d2 += (pX[d] - p1[d]) * (pX[d] - p1[d]);
            return d2 <= 0.25 * b.mDiameter * b.mDiameter;
        }
        const std::array<double, 3> p2 = point_of(b.mNode2, b.mP2);
        if (b.mKind == "RECTA") {
            for (std::size_t d = 0; d < 3; ++d)
                if (pX[d] < std::min(p1[d], p2[d]) || pX[d] > std::max(p1[d], p2[d]))
                    return false;
            return true;
        }
        if (b.mKind == "CYLIN") {
            double a[3], q[3], aa = 0.0, qa = 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                a[d] = p2[d] - p1[d];
                q[d] = pX[d] - p1[d];
                aa += a[d] * a[d];
                qa += q[d] * a[d];
            }
            if (aa == 0.0 || qa < 0.0 || qa > aa)
                return false;
            double r2 = 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                const double off = q[d] - qa / aa * a[d];
                r2 += off * off;
            }
            return r2 <= 0.25 * b.mDiameter * b.mDiameter;
        }
        return false;
    };
    std::function<std::set<std::int64_t>(std::size_t, std::set<std::size_t>&)> resolve;
    resolve = [&](std::size_t g, std::set<std::size_t>& rVisiting) -> std::set<std::int64_t> {
        std::set<std::int64_t> out;
        if (!rVisiting.insert(g).second)
            return out;
        const RadGroup& grp = groups[g];
        const std::string family = rad_group_family(grp.mKeyword);
        const bool nodes = family == "NODE";
        // GENE: `first last` id ranges; GEN_INCR: `first last step`.
        if (grp.mSubtype == "GENE" || grp.mSubtype == "GEN_INCR") {
            const std::size_t w = grp.mSubtype == "GENE" ? 2 : 3;
            auto take = [&](std::int64_t Id, std::int64_t Entity) {
                for (std::size_t k = 0; k + w <= grp.mIds.size(); k += w) {
                    const std::int64_t first = grp.mIds[k], last = grp.mIds[k + 1];
                    const std::int64_t step = w == 3 ? grp.mIds[k + 2] : 1;
                    if (Id >= first && Id <= last && step > 0 && (Id - first) % step == 0) {
                        out.insert(Entity);
                        return;
                    }
                }
            };
            if (nodes)
                for (const auto& [id, p] : node_index)
                    take(id, p);
            else
                for (const auto& [id, cell] : owner[family])
                    take(id, cell);
            rVisiting.erase(g);
            return out;
        }
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
            } else if (grp.mSubtype == "BOX" || grp.mSubtype == "BOX2") {
                // Nodes inside; elements with all (BOX) or any (BOX2) node inside.
                const bool any = grp.mSubtype == "BOX2";
                if (nodes) {
                    for (std::size_t p = 0; p < node_ids.size(); ++p)
                        if (in_box(id, &coords[3 * p], 0))
                            hits.insert(static_cast<std::int64_t>(p));
                } else {
                    for (const auto& [eid, cell] : owner[family]) {
                        const auto& conn = cell_conn[static_cast<std::size_t>(cell)];
                        std::size_t inside = 0;
                        for (std::int64_t p : conn)
                            inside +=
                                in_box(id, &coords[3 * static_cast<std::size_t>(p)], 0) ? 1 : 0;
                        if (any ? inside > 0 : inside == conn.size())
                            hits.insert(cell);
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
        // Solid faces: (cell, face) with their sorted corners.
        auto solid_faces = [&](std::int64_t Cell) {
            std::vector<std::pair<std::int64_t, std::vector<std::int64_t>>> out;
            CellType type{};
            std::vector<std::int64_t> fnodes;
            for (std::int64_t f = 0; detail::facet_nodes(mesh, Cell, f, type, fnodes); ++f) {
                const std::size_t corners = cell_type_name(type).rfind("triangle", 0) == 0 ? 3 : 4;
                fnodes.resize(std::min(corners, fnodes.size()));
                std::sort(fnodes.begin(), fnodes.end());
                out.emplace_back(f, fnodes);
            }
            return out;
        };
        std::map<std::vector<std::int64_t>, std::size_t> model_faces;  // for FREE
        bool counted = false;
        const std::map<std::string, std::string> group_keyword = {
            {"GRBRIC", "GRBRIC"}, {"GRSHEL", "GRSHEL"}, {"GRSH3N", "GRSH3N"}, {"GRTRIA", "GRTRIA"}};
        std::map<std::int64_t, std::size_t> surface_of;
        for (std::size_t k = 0; k < surfaces.size(); ++k)
            surface_of.emplace(surfaces[k].mId, k);
        std::map<std::size_t, std::set<std::pair<std::int64_t, std::int64_t>>> done;
        std::function<std::set<std::pair<std::int64_t, std::int64_t>>(std::size_t, int)> side_set;
        side_set = [&](std::size_t Index, int Depth) {
            const auto memo = done.find(Index);
            if (memo != done.end())
                return memo->second;
            std::set<std::pair<std::int64_t, std::int64_t>> out;
            const RadSurface& s = surfaces[Index];
            if (s.mKind == "SURF") {
                std::set<std::pair<std::int64_t, std::int64_t>> minus;
                for (std::int64_t raw : s.mIds) {
                    const auto it = surface_of.find(raw < 0 ? -raw : raw);
                    if (it == surface_of.end() || Depth > 16) {
                        ++dropped;
                        continue;
                    }
                    const auto sub = side_set(it->second, Depth + 1);
                    (raw < 0 ? minus : out).insert(sub.begin(), sub.end());
                }
                for (const auto& e : minus)
                    out.erase(e);
            } else if (s.mKind == "BOX" || s.mKind == "BOX2") {
                // Shell faces with all (BOX) or any (BOX2) node in the box;
                // with EXT the model's external solid faces, with ALL every
                // solid face, likewise.
                const bool any = s.mKind == "BOX2";
                const std::int64_t box = s.mIds.empty() ? 0 : s.mIds[0];
                auto inside = [&](const std::vector<std::int64_t>& rPts) {
                    std::size_t n = 0;
                    for (std::int64_t p : rPts)
                        n += in_box(box, &coords[3 * static_cast<std::size_t>(p)], 0) ? 1 : 0;
                    return any ? n > 0 : n == rPts.size();
                };
                if (s.mMode == "EXT" && !counted) {
                    for (std::size_t c = 0; c < cell_dim.size(); ++c)
                        if (cell_dim[c] == 3)
                            for (const auto& [f, key] : solid_faces(static_cast<std::int64_t>(c)))
                                ++model_faces[key];
                    counted = true;
                }
                for (std::size_t c = 0; c < cell_conn.size(); ++c) {
                    const std::string& fam = cell_family[c];
                    if (fam == "SHEL" || fam == "SH3N" || fam == "TRIA") {
                        if (inside(cell_conn[c]))
                            out.emplace(static_cast<std::int64_t>(c), 0);
                    } else if (cell_dim[c] == 3 && !s.mMode.empty()) {
                        for (const auto& [f, key] : solid_faces(static_cast<std::int64_t>(c)))
                            if ((s.mMode == "ALL" || model_faces[key] == 1) && inside(key))
                                out.emplace(static_cast<std::int64_t>(c), f);
                    }
                }
            } else if (s.mKind == "SEG") {
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
                    const detail::FacetHit* hit =
                        defined ? facets.Find(idx.data(), count) : nullptr;
                    if (!hit) {
                        ++dropped;
                        continue;
                    }
                    out.emplace(hit->mFirst.mCell, hit->mFirst.mFacet);
                }
            } else {
                // The cells the surface draws from.
                std::set<std::int64_t> cells;
                if (group_keyword.count(s.mKind)) {
                    for (std::int64_t raw : s.mIds) {
                        const auto it = group_of.find({s.mKind, raw < 0 ? -raw : raw});
                        if (it == group_of.end()) {
                            ++dropped;
                            continue;
                        }
                        std::set<std::size_t> visiting;
                        for (std::int64_t c : resolve(it->second, visiting))
                            cells.insert(c);
                    }
                } else {
                    for (std::int64_t raw : s.mIds) {
                        const std::int64_t id = raw < 0 ? -raw : raw;
                        const std::set<std::int64_t> sub =
                            s.mKind == "SUBSET" ? subset_closure(id) : std::set<std::int64_t>{};
                        for (std::size_t c = 0; c < cell_part.size(); ++c) {
                            const RadPart& part = parts[cell_part[c]];
                            const bool pick = s.mKind == "PART"     ? cell_part[c] == id
                                              : s.mKind == "SUBSET" ? sub.count(part.mSubset) > 0
                                              : s.mKind == "MAT"    ? part.mMaterial == id
                                                                    : part.mProperty == id;
                            if (pick)
                                cells.insert(static_cast<std::int64_t>(c));
                        }
                    }
                }
                // Shells: their own face. Solids: with EXT the faces no other
                // chosen solid shares, with FREE those no solid of the model
                // shares, with ALL every face (GRBRIC without a mode: EXT).
                const std::string mode = s.mKind == "GRBRIC" && s.mMode.empty() ? "EXT" : s.mMode;
                std::map<std::vector<std::int64_t>, std::size_t> chosen;
                std::vector<std::pair<
                    std::int64_t, std::vector<std::pair<std::int64_t, std::vector<std::int64_t>>>>>
                    solids;
                for (std::int64_t c : cells) {
                    const std::string& fam = cell_family[static_cast<std::size_t>(c)];
                    if (fam == "SHEL" || fam == "SH3N" || fam == "TRIA")
                        out.emplace(c, 0);
                    else if (cell_dim[static_cast<std::size_t>(c)] == 3 && !mode.empty()) {
                        solids.emplace_back(c, solid_faces(c));
                        for (const auto& [f, key] : solids.back().second)
                            ++chosen[key];
                    }
                }
                if (mode == "FREE" && !counted) {
                    for (std::size_t c = 0; c < cell_dim.size(); ++c)
                        if (cell_dim[c] == 3)
                            for (const auto& [f, key] : solid_faces(static_cast<std::int64_t>(c)))
                                ++model_faces[key];
                    counted = true;
                }
                for (const auto& [c, faces] : solids)
                    for (const auto& [f, key] : faces) {
                        const std::size_t shared =
                            mode == "FREE" ? model_faces[key] : (mode == "EXT" ? chosen[key] : 1);
                        if (shared == 1)
                            out.emplace(c, f);
                    }
            }
            done[Index] = out;
            return out;
        };
        for (std::size_t index = 0; index < surfaces.size(); ++index) {
            const RadSurface& s = surfaces[index];
            std::vector<std::int64_t> entries;
            for (const auto& [c, f] : side_set(index, 0)) {
                entries.push_back(c);
                entries.push_back(f);
            }
            add_region(s.mTitle.empty() ? "SURF_" + std::to_string(s.mId) : s.mTitle,
                       RegionKind::Side, s.mId, "SURF", std::move(entries));
        }
    }
    for (std::int64_t id : skewed_boxes)
        log::warn("Radioss: /BOX {} names a skew that is not a /SKEW/FIX; it contains nothing", id);
    for (auto& [name, values] : analytic) {
        // An ellipsoid's orientation: its skew's axes (the identity without one).
        const std::int64_t id = std::stoll(name.substr(name.rfind(':') + 1));
        const auto es = name.rfind("radioss:surf_ellips:", 0) == 0 ? ellipsoid_skew.find(id)
                                                                   : ellipsoid_skew.end();
        if (es != ellipsoid_skew.end()) {
            std::array<std::array<double, 3>, 3> axes = {
                std::array<double, 3>{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            const auto skew = static_cast<std::int64_t>(es->second[0]);
            if (skew) {
                const auto sk = skews.find(skew);
                if (sk == skews.end() || !rad_skew_axes(sk->second, axes))
                    log::warn("Radioss: /SURF/ELLIPS {} names a skew that is not a /SKEW/FIX", id);
            }
            for (const auto& a : axes)
                values.insert(values.end(), a.begin(), a.end());
        }
        NDArray a(DType::Float64, {values.size()});
        std::copy(values.begin(), values.end(), a.As<double>());
        mesh.AddFieldData(name, std::move(a));
    }
    for (std::int64_t id : missing_boxes)
        log::warn("Radioss: /BOX {} is not defined; it contains nothing", id);
    if (dropped)
        log::warn(
            "Radioss: {} group or surface entries name undefined ids or no cell facet and "
            "were dropped",
            dropped);
    // A starter deck `<run>_0000.rad`: its engine deck `<run>_0001.rad`'s controls.
    {
        const fs::path starter(rPath);
        const std::string stem = starter.stem().string();
        if (stem.size() > 5 && stem.compare(stem.size() - 5, 5, "_0000") == 0) {
            const fs::path engine =
                starter.parent_path() /
                (stem.substr(0, stem.size() - 5) + "_0001" + starter.extension().string());
            std::error_code ec;
            if (fs::is_regular_file(engine, ec))
                rad_add_engine_fields(mesh, engine.string());
        }
    }
    return mesh;
}

}  // namespace meshioplusplus
