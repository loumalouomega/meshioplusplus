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
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xyz.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

const char* const kXyzChemistry =
    "XYZ: this looks like a chemistry/molecular XYZ file (atom count, a comment line, then "
    "'element x y z' rows), not a point cloud; meshio++ has no reader for it";

bool xyz_is_space(char c) {
    return std::isspace(static_cast<unsigned char>(c)) != 0;
}

std::string xyz_strip(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && xyz_is_space(rS[b]))
        ++b;
    while (e > b && xyz_is_space(rS[e - 1]))
        --e;
    return rS.substr(b, e - b);
}

std::string xyz_lower(std::string s) {
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool xyz_all_digits(const std::string& rS) {
    return !rS.empty() && std::all_of(rS.begin(), rS.end(), [](char c) {
        return std::isdigit(static_cast<unsigned char>(c)) != 0;
    });
}

bool xyz_starts_with(const std::string& rS, const char* pPrefix) {
    return rS.rfind(pPrefix, 0) == 0;
}

std::vector<std::string> xyz_whitespace_split(const std::string& rS) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < rS.size()) {
        while (i < rS.size() && xyz_is_space(rS[i]))
            ++i;
        std::size_t j = i;
        while (j < rS.size() && !xyz_is_space(rS[j]))
            ++j;
        if (j > i)
            out.push_back(rS.substr(i, j - i));
        i = j;
    }
    return out;
}

// `line.split(delimiter)` with each token stripped and one trailing empty token dropped.
std::vector<std::string> xyz_split(const std::string& rLine, const std::string& rDelimiter) {
    if (rDelimiter.empty())
        return xyz_whitespace_split(rLine);
    std::vector<std::string> tokens;
    std::size_t start = 0;
    for (;;) {
        const std::size_t at = rLine.find(rDelimiter, start);
        if (at == std::string::npos) {
            tokens.push_back(xyz_strip(rLine.substr(start)));
            break;
        }
        tokens.push_back(xyz_strip(rLine.substr(start, at - start)));
        start = at + rDelimiter.size();
    }
    if (!tokens.empty() && tokens.back().empty())
        tokens.pop_back();
    return tokens;
}

// An element-symbol-like first token: one to three letters, then optional digits.
bool xyz_is_element(const std::string& rS) {
    std::size_t i = 0;
    while (i < rS.size() && i < 3 && std::isalpha(static_cast<unsigned char>(rS[i])))
        ++i;
    if (i == 0)
        return false;
    for (std::size_t j = i; j < rS.size(); ++j)
        if (!std::isdigit(static_cast<unsigned char>(rS[j])))
            return false;
    return true;
}

bool xyz_is_name(const std::string& rS) {
    if (rS.empty() || std::isdigit(static_cast<unsigned char>(rS[0])))
        return false;
    return std::all_of(rS.begin(), rS.end(), [](char c) {
        return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
    });
}

// Tokens of a header comment: split on whitespace, commas and semicolons.
std::vector<std::string> xyz_name_tokens(const std::string& rBody) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : rBody) {
        if (xyz_is_space(c) || c == ',' || c == ';') {
            if (!cur.empty())
                out.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty())
        out.push_back(cur);
    return out;
}

bool xyz_header_names(const std::vector<std::string>& rComments, std::vector<std::string>& rNames) {
    for (auto it = rComments.rbegin(); it != rComments.rend(); ++it) {
        std::string body;
        if (xyz_starts_with(*it, "//")) {
            body = it->substr(2);
        } else {
            std::size_t k = 0;
            while (k < it->size() && (*it)[k] == '#')
                ++k;
            body = it->substr(k);
        }
        std::vector<std::string> tokens = xyz_name_tokens(xyz_strip(body));
        if (tokens.size() >= 3 && std::all_of(tokens.begin(), tokens.end(), xyz_is_name) &&
            xyz_lower(tokens[0]) == "x" && xyz_lower(tokens[1]) == "y" &&
            xyz_lower(tokens[2]) == "z") {
            rNames = std::move(tokens);
            return true;
        }
    }
    return false;
}

using XyzTable = std::vector<std::vector<double>>;  // rows

struct XyzRange {
    bool mNan = false;
    bool mIntegral = true;
    double mMin = 0;
    double mMax = 0;
};

XyzRange xyz_range(const XyzTable& rTable, const std::vector<std::size_t>& rColumns) {
    XyzRange r;
    bool first = true;
    for (const auto& row : rTable)
        for (std::size_t c : rColumns) {
            const double v = row[c];
            if (std::isnan(v))
                r.mNan = true;
            if (v != std::rint(v))
                r.mIntegral = false;
            if (first || v < r.mMin)
                r.mMin = v;
            if (first || v > r.mMax)
                r.mMax = v;
            first = false;
        }
    return r;
}

bool xyz_byte_valued(const XyzRange& rR) {
    return !rR.mNan && rR.mIntegral && rR.mMin >= 0 && rR.mMax <= 255;
}

bool xyz_unit_valued(const XyzRange& rR) {
    return !rR.mNan && rR.mMin >= 0 && rR.mMax <= 1;
}

bool xyz_looks_like_normals(const XyzTable& rTable) {
    bool all_ok = true, any_unit = false;
    for (const auto& row : rTable) {
        const double norm = std::sqrt((row[3] * row[3] + row[4] * row[4]) + row[5] * row[5]);
        const bool unit = std::fabs(norm - 1.0) < 1e-2;
        any_unit = any_unit || unit;
        all_ok = all_ok && (unit || norm == 0);
    }
    return all_ok && any_unit;
}

std::vector<std::string> xyz_default_columns(const XyzTable& rTable, std::size_t ncols,
                                             const std::string& rSuffix) {
    if (ncols == 3)
        return {"x", "y", "z"};
    if (ncols == 4)
        return {"x", "y", "z", rSuffix == ".pts" ? "intensity" : "scalar"};
    if (ncols == 6) {
        if (rSuffix == ".xyzn")
            return {"x", "y", "z", "nx", "ny", "nz"};
        if (rSuffix == ".xyzrgb")
            return {"x", "y", "z", "r", "g", "b"};
        if (xyz_looks_like_normals(rTable))
            return {"x", "y", "z", "nx", "ny", "nz"};
        const XyzRange r = xyz_range(rTable, {3, 4, 5});
        if (xyz_byte_valued(r) || xyz_unit_valued(r))
            return {"x", "y", "z", "r", "g", "b"};
    }
    if (ncols == 7 && rSuffix == ".pts")
        return {"x", "y", "z", "intensity", "r", "g", "b"};
    throw ReadError("XYZ: cannot tell what the " + std::to_string(ncols) +
                    " columns are; pass columns=[...] (names x y z nx ny nz r g b a, any other "
                    "name is a scalar, '_' skips a column)");
}

NDArray xyz_colours(const XyzTable& rTable, const std::vector<std::size_t>& rColumns) {
    const XyzRange r = xyz_range(rTable, rColumns);
    const bool bytes = xyz_byte_valued(r);
    if (!bytes && !xyz_unit_valued(r))
        throw ReadError("XYZ: colour columns must be bytes (0..255) or unit floats (0..1)");
    NDArray out(DType::UInt8, {rTable.size(), rColumns.size()});
    std::uint8_t* dst = out.As<std::uint8_t>();
    for (std::size_t i = 0; i < rTable.size(); ++i)
        for (std::size_t j = 0; j < rColumns.size(); ++j) {
            const double v = rTable[i][rColumns[j]];
            dst[i * rColumns.size() + j] =
                static_cast<std::uint8_t>(bytes ? v : std::nearbyint(v * 255));
        }
    return out;
}

const std::map<std::string, std::string>& xyz_roles() {
    static const std::map<std::string, std::string> kRoles = {
        {"x", "x"},    {"y", "y"},         {"z", "z"},         {"nx", "nx"},       {"ny", "ny"},
        {"nz", "nz"},  {"normal_x", "nx"}, {"normal_y", "ny"}, {"normal_z", "nz"}, {"r", "r"},
        {"g", "g"},    {"b", "b"},         {"a", "a"},         {"red", "r"},       {"green", "g"},
        {"blue", "b"}, {"alpha", "a"},
    };
    return kRoles;
}

NDArray xyz_columns_to_array(const XyzTable& rTable, const std::vector<std::size_t>& rColumns,
                             bool one_d) {
    NDArray out(DType::Float64, one_d ? std::vector<std::size_t>{rTable.size()}
                                      : std::vector<std::size_t>{rTable.size(), rColumns.size()});
    double* dst = out.As<double>();
    for (std::size_t i = 0; i < rTable.size(); ++i)
        for (std::size_t j = 0; j < rColumns.size(); ++j)
            dst[i * rColumns.size() + j] = rTable[i][rColumns[j]];
    return out;
}

Mesh xyz_to_mesh(const XyzTable& rTable, const std::vector<std::string>& rNames) {
    std::map<std::string, std::size_t> roles;
    std::vector<std::pair<std::string, std::size_t>> scalars;
    for (std::size_t j = 0; j < rNames.size(); ++j) {
        const std::string lower = xyz_lower(rNames[j]);
        if (lower == "_" || lower == "skip" || lower == "ignore")
            continue;
        auto role = xyz_roles().find(lower);
        if (role == xyz_roles().end())
            scalars.emplace_back(rNames[j], j);
        else if (roles.count(role->second))
            throw ReadError("XYZ: column '" + rNames[j] + "' appears twice");
        else
            roles[role->second] = j;
    }
    for (const char* axis : {"x", "y", "z"})
        if (!roles.count(axis))
            throw ReadError("XYZ: the columns must include x, y and z");

    std::map<std::string, NDArray> point_data;
    if (roles.count("nx") || roles.count("ny") || roles.count("nz")) {
        if (!(roles.count("nx") && roles.count("ny") && roles.count("nz")))
            throw ReadError("XYZ: normals need all of nx, ny, nz");
        point_data.emplace("normals", xyz_columns_to_array(
                                          rTable, {roles["nx"], roles["ny"], roles["nz"]}, false));
    }
    if (roles.count("r") || roles.count("g") || roles.count("b")) {
        if (!(roles.count("r") && roles.count("g") && roles.count("b")))
            throw ReadError("XYZ: colours need all of r, g, b");
        std::vector<std::size_t> idx = {roles["r"], roles["g"], roles["b"]};
        if (roles.count("a"))
            idx.push_back(roles["a"]);
        point_data.emplace(roles.count("a") ? "rgba" : "rgb", xyz_colours(rTable, idx));
    } else if (roles.count("a")) {
        throw ReadError("XYZ: an alpha column needs r, g and b");
    }
    for (const auto& [name, j] : scalars) {
        if (point_data.count(name))
            throw ReadError("XYZ: column '" + name + "' appears twice");
        point_data.emplace(name, xyz_columns_to_array(rTable, {j}, true));
    }

    Mesh mesh;
    mesh.AssignPoints(xyz_columns_to_array(rTable, {roles["x"], roles["y"], roles["z"]}, false));
    NDArray cells(DType::Int64, {rTable.size(), std::size_t(1)});
    for (std::size_t i = 0; i < rTable.size(); ++i)
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellBlock("vertex", std::move(cells));
    for (auto& [name, array] : point_data)
        mesh.AddPointData(name, std::move(array));
    return mesh;
}

std::string xyz_suffix(const std::string& rPath) {
    const std::size_t dot = rPath.rfind('.');
    const std::size_t slash = rPath.find_last_of("/\\");
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
        return "";
    return xyz_lower(rPath.substr(dot));
}

std::string xyz_clean(const std::string& rName) {
    std::string out;
    bool in_run = false;
    for (char c : rName) {
        if (xyz_is_space(c) || c == ',' || c == ';') {
            if (!in_run)
                out += '_';
            in_run = true;
        } else {
            out += c;
            in_run = false;
        }
    }
    return out.empty() ? "field" : out;
}

}  // namespace

Mesh read_xyz(const std::string& rPath, const XyzReadOptions& rOptions) {
    const std::string suffix = xyz_suffix(rPath);
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    for (std::string raw; std::getline(in, raw);)
        lines.push_back(xyz_strip(raw));

    if (lines.size() >= 3 && xyz_all_digits(lines[0])) {
        const std::vector<std::string> atom = xyz_whitespace_split(lines[2]);
        if (atom.size() >= 4 && xyz_is_element(atom[0]) && xyz_lower(atom[0]) != "nan" &&
            xyz_lower(atom[0]) != "inf")
            throw ReadError(kXyzChemistry);
    }

    std::vector<std::string> comments, rows;
    std::vector<std::size_t> numbers;
    for (std::size_t k = 0; k < lines.size(); ++k) {
        const std::string& line = lines[k];
        if (line.empty())
            continue;
        if (xyz_starts_with(line, "#") || xyz_starts_with(line, "//")) {
            if (rows.empty())
                comments.push_back(line);
            continue;
        }
        rows.push_back(line);
        numbers.push_back(k + 1);
    }

    bool has_declared = false;
    long long declared = 0;
    if (suffix == ".pts" && !rows.empty() && xyz_all_digits(rows[0])) {
        has_declared = true;
        declared = std::strtoll(rows[0].c_str(), nullptr, 10);
        rows.erase(rows.begin());
        numbers.erase(numbers.begin());
    }

    std::string delimiter = rOptions.mDelimiter;
    if (delimiter.empty() && !rows.empty())
        delimiter = rows[0].find(';') != std::string::npos   ? ";"
                    : rows[0].find(',') != std::string::npos ? ","
                                                             : "";
    if (!delimiter.empty() && xyz_strip(delimiter).empty())
        delimiter.clear();

    if (rows.empty()) {
        Mesh mesh;
        mesh.AssignPoints(NDArray(DType::Float64, {std::size_t(0), std::size_t(3)}));
        mesh.AddCellBlock("vertex", NDArray(DType::Int64, {std::size_t(0), std::size_t(1)}));
        return mesh;
    }

    XyzTable table;
    table.reserve(rows.size());
    std::size_t ncols = 0;
    for (std::size_t r = 0; r < rows.size(); ++r) {
        const std::vector<std::string> tokens = xyz_split(rows[r], delimiter);
        if (r == 0)
            ncols = tokens.size();
        if (tokens.size() != ncols)
            throw ReadError("XYZ: line " + std::to_string(numbers[r]) + ": expected " +
                            std::to_string(ncols) + " columns, found " +
                            std::to_string(tokens.size()));
        std::vector<double> values(ncols);
        for (std::size_t c = 0; c < ncols; ++c) {
            const char* stop = nullptr;
            values[c] = detail::parse_double(tokens[c].c_str(), stop);
            if (tokens[c].empty() || stop == tokens[c].c_str() || *stop != '\0') {
                std::string joined;
                for (std::size_t k = 0; k < tokens.size(); ++k)
                    joined += (k ? " " : "") + tokens[k];
                throw ReadError("XYZ: line " + std::to_string(numbers[r]) + ": '" + joined +
                                "' is not numeric");
            }
        }
        table.push_back(std::move(values));
    }
    if (has_declared && declared != static_cast<long long>(table.size()))
        throw ReadError("XYZ: the .pts header declares " + std::to_string(declared) +
                        " points, found " + std::to_string(table.size()));

    std::vector<std::string> names;
    if (!rOptions.mColumns.empty()) {
        names = rOptions.mColumns;
        if (names.size() != ncols)
            throw ReadError("XYZ: columns= names " + std::to_string(names.size()) +
                            " columns, the file has " + std::to_string(ncols));
    } else {
        if (xyz_header_names(comments, names) && names.size() != ncols) {
            log::warn("XYZ: header names {} columns, the file has {}; ignoring it", names.size(),
                      ncols);
            names.clear();
        }
        if (names.empty())
            names = xyz_default_columns(table, ncols, suffix);
    }
    return xyz_to_mesh(table, names);
}

void write_xyz(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFormat) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    if (dim < 3) {
        log::warn("XYZ requires 3D points; padding with zeros.");
        detail::provenance_note("point-padding", "points padded with zero coordinates to 3D");
    }
    std::string skipped;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() == "vertex")
            continue;
        if (!skipped.empty())
            skipped += ", ";
        skipped += cb.Type();
    }
    if (!skipped.empty()) {
        log::warn("XYZ holds points only. Skipping {} cells.", skipped);
        detail::provenance_note("cells-dropped",
                                "cell block(s) of type " + skipped + " have no XYZ equivalent");
    }
    if (rMesh.NumCellData() > 0)
        detail::provenance_note("data-dropped", "cell data has no XYZ equivalent");

    std::string format;
    if (!rFloatFormat.empty()) {
        std::size_t i = 0;
        while (i < rFloatFormat.size() && std::isdigit(static_cast<unsigned char>(rFloatFormat[i])))
            ++i;
        if (i < rFloatFormat.size() && rFloatFormat[i] == '.') {
            ++i;
            while (i < rFloatFormat.size() &&
                   std::isdigit(static_cast<unsigned char>(rFloatFormat[i])))
                ++i;
        }
        if (i + 1 != rFloatFormat.size() || !std::strchr("eEfFgG", rFloatFormat[i]))
            throw WriteError("XYZ: float_fmt must look like '.16e' or '.9g', got '" + rFloatFormat +
                             "'");
        format = "%" + rFloatFormat;
    }

    struct Column {
        const NDArray* mArray;
        std::size_t mStride;
        std::size_t mOffset;
    };
    std::vector<std::string> names = {"x", "y", "z"};
    std::vector<Column> columns;
    NDArray zeros(DType::Float64, {n, std::size_t(1)});
    for (std::size_t a = 0; a < 3; ++a)
        columns.push_back(a < dim ? Column{&points, dim, a} : Column{&zeros, 1, 0});
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& array = rMesh.PointData(name);
        std::size_t count = 1;
        for (std::size_t d = 1; d < array.Shape().size(); ++d)
            count *= array.Shape()[d];
        const std::size_t width = n == 0 ? 1 : count;  // an empty cloud names one column
        std::vector<std::string> labels;
        if (name == "normals" && width == 3)
            labels = {"nx", "ny", "nz"};
        else if (name == "rgb" && width == 3)
            labels = {"r", "g", "b"};
        else if (name == "rgba" && width == 4)
            labels = {"r", "g", "b", "a"};
        else if (width == 1)
            labels = {xyz_clean(name)};
        else
            for (std::size_t k = 0; k < width; ++k)
                labels.push_back(xyz_clean(name) + "_" + std::to_string(k));
        if (std::any_of(labels.begin(), labels.end(), [&](const std::string& rLabel) {
                return std::find(names.begin(), names.end(), rLabel) != names.end();
            })) {
            // e.g. rgb and rgba together, or point data called "x": the reader would see
            // a repeated role, so fall back to indexed names.
            labels.clear();
            for (std::size_t k = 0; k < width; ++k)
                labels.push_back(xyz_clean(name) + "_" + std::to_string(k));
        }
        for (std::size_t k = 0; k < width; ++k) {
            names.push_back(labels[k]);
            columns.push_back(Column{&array, width, k});
        }
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os << detail::provenance_render_lines(detail::SlotTier::Block, "# ") << "#";
    for (const auto& name : names)
        os << ' ' << name;
    os << '\n';
    char buf[64];
    std::string line;
    for (std::size_t i = 0; i < n; ++i) {
        line.clear();
        for (std::size_t c = 0; c < columns.size(); ++c) {
            const NDArray& array = *columns[c].mArray;
            const std::size_t at = i * columns[c].mStride + columns[c].mOffset;
            if (c)
                line += ' ';
            if (!detail::is_float_dtype(array.Dtype())) {
                line += array.Dtype() == DType::UInt64
                            ? std::to_string(array.As<std::uint64_t>()[at])
                            : std::to_string(detail::read_int(array, at));
                continue;
            }
            const double v = detail::read_double(array, at);
            if (std::isnan(v)) {
                line += "nan";
                continue;
            }
            const char* fmt = !format.empty()                   ? format.c_str()
                              : array.Dtype() == DType::Float32 ? "%.9g"
                                                                : "%.17g";
            detail::snprintf_c(buf, sizeof(buf), fmt, v);
            line += buf;
        }
        line += '\n';
        os << line;
    }
}

}  // namespace meshioplusplus
