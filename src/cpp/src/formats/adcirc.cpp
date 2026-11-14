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
// ADCIRC fort.14 grids and fort.63/fort.64 results. See formats/adcirc.hpp.

// System includes
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/adcirc.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kAdcBoundaryPrefix = "adcirc:boundary:";

/// A file's lines, read one at a time.
struct AdcLines {
    std::vector<std::string> mLines;
    std::size_t mAt = 0;
    std::string mPath;

    explicit AdcLines(const std::string& rPath) : mPath(rPath) {
        const detail::FileSource source =
            detail::open_source(rPath, "meshio++: adcirc: could not open file: " + rPath);
        const std::string_view data = source.View();
        std::size_t b = 0;
        while (b < data.size()) {
            std::size_t e = data.find('\n', b);
            if (e == std::string_view::npos)
                e = data.size();
            std::string line(data.substr(b, e - b));
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            mLines.push_back(std::move(line));
            b = e + 1;
        }
    }

    bool Eof() const { return mAt >= mLines.size(); }

    const std::string& Next(const char* pWhat) {
        if (mAt >= mLines.size())
            throw ReadError("meshio++: adcirc: '" + mPath + "' ends before its " + pWhat);
        return mLines[mAt++];
    }
};

/// The leading numeric tokens of a line (ADCIRC puts labels and `!` comments after
/// the numbers).
std::vector<std::string> adc_numbers(const std::string& rLine) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < rLine.size()) {
        while (i < rLine.size() && (rLine[i] == ' ' || rLine[i] == '\t' || rLine[i] == ','))
            ++i;
        const std::size_t b = i;
        while (i < rLine.size() && rLine[i] != ' ' && rLine[i] != '\t' && rLine[i] != ',')
            ++i;
        if (i == b)
            break;
        std::string tok = rLine.substr(b, i - b);
        const char c = tok[0];
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+' || c == '.'))
            break;
        out.push_back(std::move(tok));
    }
    return out;
}

std::vector<std::string> adc_need(AdcLines& rIn, const char* pWhat, std::size_t Count) {
    std::vector<std::string> t = adc_numbers(rIn.Next(pWhat));
    if (t.size() < Count)
        throw ReadError(std::string("meshio++: adcirc: ") + pWhat + " needs " +
                        std::to_string(Count) + " numbers");
    return t;
}

long long adc_int(const std::string& rTok) {
    long long v = 0;
    const auto [p, ec] = std::from_chars(rTok.data(), rTok.data() + rTok.size(), v);
    if (ec != std::errc() || p != rTok.data() + rTok.size()) {
        // Some writers print integers as reals ("3.0").
        const char* end = nullptr;
        const double d = detail::parse_double(rTok.c_str(), end);
        if (end == rTok.c_str())
            throw ReadError("meshio++: adcirc: '" + rTok + "' is not an integer");
        return detail::checked_integer<long long>(d, "adcirc");
    }
    return v;
}

double adc_real(const std::string& rTok) {
    std::string t = rTok;
    for (char& c : t)
        if (c == 'd' || c == 'D')
            c = 'e';
    const char* end = nullptr;
    const double v = detail::parse_double(t.c_str(), end);
    if (end == t.c_str())
        throw ReadError("meshio++: adcirc: '" + rTok + "' is not a number");
    return v;
}

/// Land boundary types whose node records carry barrier columns (heights,
/// coefficients, connected nodes): read past them, never written.
bool adc_is_barrier(long long Type) {
    const long long t = Type % 10;
    return t == 3 || t == 4 || t == 5;
}

/// The result variable a `fort.NN` file holds.
std::string adc_result_name(const std::string& rPath) {
    const std::string base = std::filesystem::path(rPath).filename().string();
    if (base == "fort.63")
        return "zeta";
    if (base == "fort.64")
        return "velocity";
    if (base == "fort.73")
        return "pressure";
    if (base == "fort.74")
        return "wind_velocity";
    return "adcirc:" + base;
}

struct AdcResult {
    std::string mName;
    std::size_t mWidth = 1;
    std::vector<double> mTimes;
    std::vector<std::size_t> mStarts;  ///< the line of each dataset's header
};

/// A result file's header and the position of every dataset.
AdcResult adc_scan_result(AdcLines& rIn, std::size_t NumPoints) {
    AdcResult r;
    r.mName = adc_result_name(rIn.mPath);
    rIn.Next("title");
    const auto head = adc_need(rIn, "dataset header (NDSETS NP DT NSPOOL IRTYPE)", 5);
    const long long nsets = adc_int(head[0]);
    if (static_cast<std::size_t>(adc_int(head[1])) != NumPoints)
        throw ReadError("meshio++: adcirc: '" + rIn.mPath + "' has " + head[1] +
                        " nodes but the grid has " + std::to_string(NumPoints));
    const long long irtype = adc_int(head[4]);
    if (irtype < 1 || irtype > 3)
        throw ReadError("meshio++: adcirc: '" + rIn.mPath + "' has record type " + head[4]);
    r.mWidth = static_cast<std::size_t>(irtype);
    // NDSETS can be a placeholder (-1 or larger than written): walk what is there.
    (void)nsets;
    while (!rIn.Eof()) {
        const std::string& line = rIn.Next("dataset");
        const auto t = adc_numbers(line);
        if (t.empty())
            continue;
        if (t.size() < 2)
            throw ReadError("meshio++: adcirc: a dataset header needs TIME IT");
        r.mStarts.push_back(rIn.mAt - 1);
        r.mTimes.push_back(adc_real(t[0]));
        const std::size_t rows =
            t.size() >= 4 ? static_cast<std::size_t>(adc_int(t[2])) : NumPoints;
        if (rows > NumPoints)
            throw ReadError("meshio++: adcirc: a sparse dataset lists more nodes than the grid");
        rIn.mAt += rows;
    }
    return r;
}

/// Dataset @p K of a scanned result, as a `(NP[, width])` array.
NDArray adc_read_dataset(AdcLines& rIn, const AdcResult& rR, std::size_t K,
                         const std::unordered_map<long long, std::size_t>& rIndex,
                         std::size_t NumPoints) {
    rIn.mAt = rR.mStarts[K];
    const auto t = adc_numbers(rIn.Next("dataset header"));
    const bool sparse = t.size() >= 4;
    const std::size_t rows = sparse ? static_cast<std::size_t>(adc_int(t[2])) : NumPoints;
    const double fill = sparse ? adc_real(t[3]) : 0.0;
    NDArray out(DType::Float64, rR.mWidth == 1 ? std::vector<std::size_t>{NumPoints}
                                               : std::vector<std::size_t>{NumPoints, rR.mWidth});
    std::fill_n(out.As<double>(), out.Size(), fill);
    for (std::size_t r = 0; r < rows; ++r) {
        const auto v = adc_need(rIn, "dataset value", 1 + rR.mWidth);
        const auto it = rIndex.find(adc_int(v[0]));
        if (it == rIndex.end())
            throw ReadError("meshio++: adcirc: a dataset names node " + v[0] +
                            ", which the grid does not have");
        for (std::size_t c = 0; c < rR.mWidth; ++c)
            out.As<double>()[it->second * rR.mWidth + c] = adc_real(v[1 + c]);
    }
    return out;
}

std::vector<std::string> adc_result_paths(const std::string& rPath,
                                          const AdcircReadOptions& rOptions) {
    if (!rOptions.mResults.empty())
        return rOptions.mResults;
    std::vector<std::string> out;
    const std::filesystem::path p(rPath);
    if (!rOptions.mSiblings || p.filename() != "fort.14")
        return out;
    for (const char* name : {"fort.63", "fort.64"}) {
        const auto sib = p.parent_path() / name;
        std::error_code ec;
        if (std::filesystem::is_regular_file(sib, ec))
            out.push_back(sib.string());
    }
    return out;
}

/// The grid: points, depth, triangles, boundaries; @p rIndex maps ids to rows.
Mesh adc_read_grid(const std::string& rPath, std::unordered_map<long long, std::size_t>& rIndex) {
    AdcLines in(rPath);
    const std::string title = in.Next("title");
    const auto counts = adc_need(in, "element and node counts (NE NP)", 2);
    const std::size_t limit = in.mLines.size();
    const std::size_t ne = detail::checked_count(adc_int(counts[0]), limit, "adcirc", "element");
    const std::size_t np = detail::checked_count(adc_int(counts[1]), limit, "adcirc", "node");
    NDArray pts(DType::Float64, {np, 3});
    NDArray depth(DType::Float64, {np});
    for (std::size_t i = 0; i < np; ++i) {
        const auto t = adc_need(in, "node", 4);
        if (!rIndex.emplace(adc_int(t[0]), i).second)
            throw ReadError("meshio++: adcirc: duplicate node id " + t[0]);
        pts.As<double>()[i * 3] = adc_real(t[1]);
        pts.As<double>()[i * 3 + 1] = adc_real(t[2]);
        pts.As<double>()[i * 3 + 2] = 0.0;
        depth.As<double>()[i] = adc_real(t[3]);
    }
    const auto node = [&](const std::string& rTok) {
        const auto it = rIndex.find(adc_int(rTok));
        if (it == rIndex.end())
            throw ReadError("meshio++: adcirc: node " + rTok + " is not in the grid");
        return static_cast<std::int64_t>(it->second);
    };
    NDArray tri(DType::Int64, {ne, 3});
    for (std::size_t e = 0; e < ne; ++e) {
        const auto t = adc_need(in, "element", 5);
        if (adc_int(t[1]) != 3)
            throw ReadError("meshio++: adcirc: element " + t[0] + " has " + t[1] +
                            " nodes; ADCIRC grids are triangles");
        for (std::size_t k = 0; k < 3; ++k)
            tri.As<std::int64_t>()[e * 3 + k] = node(t[2 + k]);
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    if (ne > 0)
        mesh.AddCellBlock("triangle", std::move(tri));
    mesh.AddPointData("depth", std::move(depth));
    NDArray tbytes(DType::UInt8, {title.size()});
    std::copy(title.begin(), title.end(), reinterpret_cast<char*>(tbytes.Data()));
    mesh.AddFieldData("adcirc:title", std::move(tbytes));

    // Boundaries are optional: a grid can end after its elements.
    const auto skip_blank = [&] {
        while (!in.Eof() && adc_numbers(in.mLines[in.mAt]).empty())
            ++in.mAt;
    };
    skip_blank();
    if (in.Eof())
        return mesh;
    const auto read_strings = [&](const char* pKind, bool Open) {
        const std::size_t count = detail::checked_count(
            adc_int(adc_need(in, Open ? "NOPE" : "NBOU", 1)[0]), limit, "adcirc", "boundary");
        adc_need(in, Open ? "NETA" : "NVEL", 1);
        for (std::size_t k = 0; k < count; ++k) {
            const auto h = adc_need(in, "boundary header", 1);
            const std::size_t n = detail::checked_count(adc_int(h[0]), limit, "adcirc", "node");
            const long long type = h.size() > 1 ? adc_int(h[1]) : 0;
            NDArray ids(DType::Int64, {n});
            for (std::size_t r = 0; r < n; ++r)
                ids.As<std::int64_t>()[r] = node(adc_need(in, "boundary node", 1)[0]);
            const std::string name = std::string(pKind) + "_" + std::to_string(k);
            // A region's entries are a sorted set; the string's order rides beside it.
            mesh.AddFieldData(std::string(kAdcBoundaryPrefix) + name, ids);
            mesh.AddRegion(Region(name, RegionKind::Point, -1, type, std::move(ids)));
        }
    };
    read_strings("open_boundary", true);
    skip_blank();
    if (!in.Eof())
        read_strings("land_boundary", false);
    return mesh;
}

}  // namespace

Mesh read_adcirc(const std::string& rPath, const AdcircReadOptions& rOptions) {
    std::unordered_map<long long, std::size_t> index;
    Mesh mesh = adc_read_grid(rPath, index);
    const std::size_t np = mesh.NumPoints();
    const std::vector<std::string> paths = adc_result_paths(rPath, rOptions);
    std::vector<double> times;
    bool first = true;
    std::size_t step = 0;
    for (const std::string& path : paths) {
        AdcLines in(path);
        const AdcResult r = adc_scan_result(in, np);
        if (first) {
            times = r.mTimes;
            ReadOptions o;
            o.mTimeStep = rOptions.mTimeStep;
            step = o.ResolveTimeStep(times.size());
            first = false;
        } else if (r.mTimes.size() != times.size()) {
            throw ReadError("meshio++: adcirc: '" + path + "' has " +
                            std::to_string(r.mTimes.size()) + " datasets but '" + paths.front() +
                            "' has " + std::to_string(times.size()));
        }
        if (r.mTimes.empty())
            continue;
        mesh.AddPointData(r.mName, adc_read_dataset(in, r, step, index, np));
    }
    if (!times.empty()) {
        NDArray t(DType::Float64, {1});
        t.As<double>()[0] = times[step];
        mesh.AddFieldData("meshio:time", std::move(t));
    }
    return mesh;
}

Mesh read_adcirc(const std::string& rPath, const ReadOptions& rOptions) {
    AdcircReadOptions o;
    o.mTimeStep = rOptions.mTimeStep;
    return read_adcirc(rPath, o);
}

MeshMetadata read_adcirc_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    MeshMetadata meta = metadata_from_mesh(read_adcirc(rPath, rOptions));
    meta.mFellBackToFullRead = true;
    std::unordered_map<long long, std::size_t> index;
    const Mesh grid = adc_read_grid(rPath, index);
    const auto paths = adc_result_paths(rPath, AdcircReadOptions{});
    if (!paths.empty()) {
        AdcLines in(paths.front());
        meta.mTimeValues = adc_scan_result(in, grid.NumPoints()).mTimes;
    }
    return meta;
}

void write_adcirc(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::size_t> tris;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::string type = rMesh.Cells(b).Type();
        if (type != "triangle")
            throw WriteError("meshio++: adcirc: a '" + type +
                             "' block; ADCIRC grids hold triangles only");
        tris.push_back(b);
    }
    std::vector<const Region*> open, land;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Point)
            continue;
        if (r.mName.rfind("open_boundary", 0) == 0)
            open.push_back(&r);
        else if (r.mName.rfind("land_boundary", 0) == 0) {
            if (adc_is_barrier(r.mTag))
                throw WriteError("meshio++: adcirc: land boundary '" + r.mName + "' has type " +
                                 std::to_string(r.mTag) +
                                 ", whose records need barrier columns meshio++ does not keep");
            land.push_back(&r);
        }
    }
    std::string title = "meshio++";
    if (rMesh.HasFieldData("adcirc:title")) {
        const NDArray& t = rMesh.FieldData("adcirc:title");
        title.assign(reinterpret_cast<const char*>(t.Data()), t.Size());
    }
    const std::size_t np = rMesh.NumPoints();
    const NDArray* depth = rMesh.HasPointData("depth") ? &rMesh.PointData("depth") : nullptr;
    std::size_t ne = 0;
    for (std::size_t b : tris)
        ne += rMesh.Cells(b).NumCells();
    std::string out = title + "\n" + std::to_string(ne) + " " + std::to_string(np) + "\n";
    char buf[128];
    const std::size_t pd = rMesh.PointDim();
    for (std::size_t i = 0; i < np; ++i) {
        const double x = detail::read_double(rMesh.Points(), i * pd);
        const double y = pd > 1 ? detail::read_double(rMesh.Points(), i * pd + 1) : 0.0;
        const double d = depth ? detail::read_double(*depth, i) : 0.0;
        detail::snprintf_c(buf, sizeof buf, "%zu %.17g %.17g %.17g\n", i + 1, x, y, d);
        out += buf;
    }
    std::size_t e = 0;
    for (std::size_t b : tris) {
        const auto cb = rMesh.Cells(b);
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            detail::snprintf_c(buf, sizeof buf, "%zu 3 %lld %lld %lld\n", ++e,
                               static_cast<long long>(detail::read_int(cb.Conn(), r * 3)) + 1,
                               static_cast<long long>(detail::read_int(cb.Conn(), r * 3 + 1)) + 1,
                               static_cast<long long>(detail::read_int(cb.Conn(), r * 3 + 2)) + 1);
            out += buf;
        }
    }
    const auto strings = [&](const std::vector<const Region*>& rList, bool Open) {
        std::size_t total = 0;
        for (const Region* r : rList)
            total += r->mEntries.Size();
        out += std::to_string(rList.size()) +
               (Open ? " = Number of open boundaries\n" : " = Number of land boundaries\n");
        out += std::to_string(total) + (Open ? " = Total number of open boundary nodes\n"
                                             : " = Total number of land boundary nodes\n");
        for (std::size_t k = 0; k < rList.size(); ++k) {
            const Region& r = *rList[k];
            out += std::to_string(r.mEntries.Size());
            if (!Open || r.mTag >= 0)
                out += " " + std::to_string(std::max<std::int64_t>(r.mTag, 0));
            out += " = Number of nodes for " + std::string(Open ? "open" : "land") + " boundary " +
                   std::to_string(k + 1) + "\n";
            // The string's own order when the mesh kept it and it is the same set.
            const NDArray* order = nullptr;
            const std::string key = std::string(kAdcBoundaryPrefix) + r.mName;
            if (rMesh.HasFieldData(key) && rMesh.FieldData(key).Size() == r.mEntries.Size()) {
                std::vector<std::int64_t> a, b;
                for (std::size_t i = 0; i < r.mEntries.Size(); ++i) {
                    a.push_back(detail::read_int(rMesh.FieldData(key), i));
                    b.push_back(detail::read_int(r.mEntries, i));
                }
                std::sort(a.begin(), a.end());
                if (a == b)
                    order = &rMesh.FieldData(key);
            }
            const NDArray& ids = order ? *order : r.mEntries;
            for (std::size_t i = 0; i < ids.Size(); ++i)
                out += std::to_string(detail::read_int(ids, i) + 1) + "\n";
        }
    };
    // Boundaries in their number's order (`open_boundary_10` after `_2`).
    const auto by_number = [](const Region* pA, const Region* pB) {
        const auto key = [](const std::string& rName) {
            const std::size_t u = rName.rfind('_');
            const std::string tail = u == std::string::npos ? "" : rName.substr(u + 1);
            const bool digits = !tail.empty() && tail.size() < 18 &&
                                std::all_of(tail.begin(), tail.end(),
                                            [](unsigned char c) { return std::isdigit(c) != 0; });
            return std::make_pair(digits ? std::stoll(tail) : (1LL << 62), rName);
        };
        return key(pA->mName) < key(pB->mName);
    };
    std::sort(open.begin(), open.end(), by_number);
    std::sort(land.begin(), land.end(), by_number);
    strings(open, true);
    strings(land, false);
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: adcirc: could not open file for writing: " + rPath);
    os << out;
    if (!os)
        throw WriteError("meshio++: adcirc: failed writing " + rPath);
}

}  // namespace meshioplusplus
