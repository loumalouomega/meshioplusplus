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
// ASPRS LAS point clouds. See formats/las.hpp.

// System includes
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/las.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

template <class T>
T las_get(std::string_view Data, std::size_t At) {
    if (At + sizeof(T) > Data.size())
        throw ReadError("meshio++: las: the file is truncated");
    std::array<unsigned char, sizeof(T)> b{};
    std::memcpy(b.data(), Data.data() + At, sizeof(T));
    if constexpr (std::endian::native == std::endian::big)
        std::reverse(b.begin(), b.end());
    T v;
    std::memcpy(&v, b.data(), sizeof(T));
    return v;
}

template <class T>
void las_put(std::string& rOut, T Value) {
    std::array<unsigned char, sizeof(T)> b{};
    std::memcpy(b.data(), &Value, sizeof(T));
    if constexpr (std::endian::native == std::endian::big)
        std::reverse(b.begin(), b.end());
    rOut.append(reinterpret_cast<const char*>(b.data()), sizeof(T));
}

/// The standard record length of each point format.
constexpr std::array<std::size_t, 11> kLasRecord = {20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67};

/// One point-data field: its name, dtype and how to take it from a record.
struct LasField {
    const char* mName;
    DType mType;
    std::function<double(const char*)> mGet;
};

std::vector<LasField> las_fields(int Format) {
    const auto u8 = [](std::size_t Off) {
        return [Off](const char* p) {
            return static_cast<double>(static_cast<unsigned char>(p[Off]));
        };
    };
    const auto bits = [](std::size_t Off, int Shift, int Width) {
        return [=](const char* p) {
            return static_cast<double>((static_cast<unsigned char>(p[Off]) >> Shift) &
                                       ((1 << Width) - 1));
        };
    };
    const auto u16 = [](std::size_t Off) {
        return [Off](const char* p) {
            return static_cast<double>(las_get<std::uint16_t>(std::string_view(p, Off + 2), Off));
        };
    };
    std::vector<LasField> f;
    f.push_back({"intensity", DType::UInt16, u16(12)});
    const bool modern = Format >= 6;
    if (!modern) {
        f.push_back({"return_number", DType::UInt8, bits(14, 0, 3)});
        f.push_back({"number_of_returns", DType::UInt8, bits(14, 3, 3)});
        f.push_back({"scan_direction_flag", DType::UInt8, bits(14, 6, 1)});
        f.push_back({"edge_of_flight_line", DType::UInt8, bits(14, 7, 1)});
        f.push_back({"classification", DType::UInt8, bits(15, 0, 5)});
        f.push_back({"synthetic", DType::UInt8, bits(15, 5, 1)});
        f.push_back({"key_point", DType::UInt8, bits(15, 6, 1)});
        f.push_back({"withheld", DType::UInt8, bits(15, 7, 1)});
        f.push_back({"scan_angle_rank", DType::Int8, [](const char* p) {
                         return static_cast<double>(static_cast<std::int8_t>(p[16]));
                     }});
        f.push_back({"user_data", DType::UInt8, u8(17)});
        f.push_back({"point_source_id", DType::UInt16, u16(18)});
    } else {
        f.push_back({"return_number", DType::UInt8, bits(14, 0, 4)});
        f.push_back({"number_of_returns", DType::UInt8, bits(14, 4, 4)});
        f.push_back({"synthetic", DType::UInt8, bits(15, 0, 1)});
        f.push_back({"key_point", DType::UInt8, bits(15, 1, 1)});
        f.push_back({"withheld", DType::UInt8, bits(15, 2, 1)});
        f.push_back({"overlap", DType::UInt8, bits(15, 3, 1)});
        f.push_back({"scanner_channel", DType::UInt8, bits(15, 4, 2)});
        f.push_back({"scan_direction_flag", DType::UInt8, bits(15, 6, 1)});
        f.push_back({"edge_of_flight_line", DType::UInt8, bits(15, 7, 1)});
        f.push_back({"classification", DType::UInt8, u8(16)});
        f.push_back({"user_data", DType::UInt8, u8(17)});
        f.push_back({"scan_angle", DType::Int16, [](const char* p) {
                         return static_cast<double>(
                             las_get<std::int16_t>(std::string_view(p, 20), 18));
                     }});
        f.push_back({"point_source_id", DType::UInt16, u16(20)});
    }
    std::size_t gps = 0, rgb = 0, nir = 0;
    switch (Format) {
        case 1:
        case 4:
            gps = 20;
            break;
        case 2:
            rgb = 20;
            break;
        case 3:
        case 5:
            gps = 20;
            rgb = 28;
            break;
        case 6:
        case 9:
            gps = 22;
            break;
        case 7:
            gps = 22;
            rgb = 30;
            break;
        case 8:
        case 10:
            gps = 22;
            rgb = 30;
            nir = 36;
            break;
        default:
            break;
    }
    if (gps)
        f.push_back({"gps_time", DType::Float64, [gps](const char* p) {
                         return las_get<double>(std::string_view(p, gps + 8), gps);
                     }});
    if (rgb) {
        f.push_back({"red", DType::UInt16, u16(rgb)});
        f.push_back({"green", DType::UInt16, u16(rgb + 2)});
        f.push_back({"blue", DType::UInt16, u16(rgb + 4)});
    }
    if (nir)
        f.push_back({"nir", DType::UInt16, u16(nir)});
    return f;
}

NDArray las_bytes(std::string_view Text) {
    NDArray out(DType::UInt8, {Text.size()});
    if (!Text.empty())
        std::memcpy(out.Data(), Text.data(), Text.size());
    return out;
}

}  // namespace

Mesh read_las(const std::string& rPath, const ReadOptions& rOptions, const LasReadOptions& rLas) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: las: could not open file: " + rPath);
    const std::string_view d = source.View();
    if (d.size() < 227 || d.substr(0, 4) != "LASF")
        throw ReadError("meshio++: las: '" + rPath + "' is not a LAS file (no LASF signature)");
    const int major = static_cast<unsigned char>(d[24]), minor = static_cast<unsigned char>(d[25]);
    if (major != 1 || minor > 4)
        throw ReadError("meshio++: las: version " + std::to_string(major) + "." +
                        std::to_string(minor) + " is not read (1.0 to 1.4)");
    const std::size_t header = las_get<std::uint16_t>(d, 94);
    const std::size_t offset = las_get<std::uint32_t>(d, 96);
    const std::size_t nvlr = las_get<std::uint32_t>(d, 100);
    const int format = static_cast<unsigned char>(d[104]) & 0x3F;  // bits 6-7: LAZ compression
    if (static_cast<unsigned char>(d[104]) & 0xC0)
        throw ReadError("meshio++: las: '" + rPath +
                        "' is LAZ-compressed; meshio++ reads uncompressed LAS (LAZ needs LASzip)");
    if (format > 10)
        throw ReadError("meshio++: las: point format " + std::to_string(format) + " is not read");
    const std::size_t record = las_get<std::uint16_t>(d, 105);
    if (record < kLasRecord[static_cast<std::size_t>(format)])
        throw ReadError("meshio++: las: point record length " + std::to_string(record) +
                        " is shorter than format " + std::to_string(format) + " needs");
    std::uint64_t count = las_get<std::uint32_t>(d, 107);
    if (minor >= 4 && header >= 375 && count == 0)
        count = las_get<std::uint64_t>(d, 247);
    std::array<double, 3> scale{}, origin{};
    for (std::size_t c = 0; c < 3; ++c) {
        scale[c] = las_get<double>(d, 131 + 8 * c);
        origin[c] = las_get<double>(d, 155 + 8 * c);
    }
    if (header > d.size() || offset > d.size() || offset < header)
        throw ReadError("meshio++: las: the header size or point offset is out of range");
    if (count > (d.size() - offset) / record)
        throw ReadError("meshio++: las: the header promises " + std::to_string(count) +
                        " points but the file holds " +
                        std::to_string((d.size() - offset) / record));
    // The VLRs, verbatim.
    std::size_t at = header;
    for (std::size_t k = 0; k < nvlr; ++k) {
        if (at + 54 > offset)
            throw ReadError("meshio++: las: a VLR runs past the point data");
        at += 54 + las_get<std::uint16_t>(d, at + 20);
    }
    if (at > offset)
        throw ReadError("meshio++: las: the VLRs run past the point data");
    const std::string_view vlrs = d.substr(header, at - header);

    const std::size_t every = std::max<std::size_t>(rLas.mEvery, 1);
    const std::size_t n = static_cast<std::size_t>((count + every - 1) / every);
    NDArray pts(DType::Float64, {n, 3});
    const std::vector<LasField> fields = las_fields(format);
    std::vector<NDArray> values;
    std::vector<bool> wanted;
    for (const LasField& f : fields) {
        wanted.push_back(rOptions.WantsAnyData() && rOptions.WantsArray(f.mName));
        values.emplace_back(f.mType, std::vector<std::size_t>{wanted.back() ? n : 0});
    }
    for (std::size_t i = 0; i < n; ++i) {
        const char* p = d.data() + offset + i * every * record;
        for (std::size_t c = 0; c < 3; ++c)
            pts.As<double>()[i * 3 + c] =
                las_get<std::int32_t>(std::string_view(p, 12), 4 * c) * scale[c] + origin[c];
        for (std::size_t k = 0; k < fields.size(); ++k)
            if (wanted[k])
                detail::write_double(values[k], i, fields[k].mGet(p));
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    NDArray cells(DType::Int64, {n, 1});
    for (std::size_t i = 0; i < n; ++i)
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    if (n > 0)
        mesh.AddCellBlock("vertex", std::move(cells));
    for (std::size_t k = 0; k < fields.size(); ++k)
        if (wanted[k])
            mesh.AddPointData(fields[k].mName, std::move(values[k]));
    NDArray s(DType::Float64, {3}), o(DType::Float64, {3}), v(DType::Int64, {2}),
        f(DType::Int64, {1});
    for (std::size_t c = 0; c < 3; ++c) {
        s.As<double>()[c] = scale[c];
        o.As<double>()[c] = origin[c];
    }
    v.As<std::int64_t>()[0] = major;
    v.As<std::int64_t>()[1] = minor;
    f.As<std::int64_t>()[0] = format;
    mesh.AddFieldData("las:scale", std::move(s));
    mesh.AddFieldData("las:offset", std::move(o));
    mesh.AddFieldData("las:version", std::move(v));
    mesh.AddFieldData("las:point_format", std::move(f));
    if (!vlrs.empty())
        mesh.AddFieldData("las:vlrs", las_bytes(vlrs));
    if (minor >= 4 && header >= 375 && las_get<std::uint32_t>(d, 243) > 0)
        log::warn("las: '{}' has extended VLRs, which are not kept", rPath);
    return mesh;
}

Mesh read_las(const std::string& rPath, const ReadOptions& rOptions) {
    return read_las(rPath, rOptions, LasReadOptions{});
}

void write_las(const std::string& rPath, const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints(), pd = rMesh.PointDim();
    const auto has = [&](const char* pName) {
        return rMesh.HasPointData(pName) && rMesh.PointData(pName).Size() == n;
    };
    const bool rgb = has("red") && has("green") && has("blue");
    const int format = rgb ? (has("nir") ? 8 : 7) : 6;
    const std::size_t record = kLasRecord[static_cast<std::size_t>(format)];
    std::array<double, 3> scale = {0.001, 0.001, 0.001}, origin{}, lo{}, hi{};
    const auto coord = [&](std::size_t i, std::size_t c) {
        return c < pd ? detail::read_double(rMesh.Points(), i * pd + c) : 0.0;
    };
    for (std::size_t c = 0; c < 3; ++c) {
        lo[c] = n ? coord(0, c) : 0.0;
        hi[c] = lo[c];
        for (std::size_t i = 1; i < n; ++i) {
            lo[c] = std::min(lo[c], coord(i, c));
            hi[c] = std::max(hi[c], coord(i, c));
        }
        origin[c] = std::floor(lo[c]);
    }
    if (rMesh.HasFieldData("las:scale") && rMesh.FieldData("las:scale").Size() == 3)
        for (std::size_t c = 0; c < 3; ++c)
            scale[c] = detail::read_double(rMesh.FieldData("las:scale"), c);
    if (rMesh.HasFieldData("las:offset") && rMesh.FieldData("las:offset").Size() == 3)
        for (std::size_t c = 0; c < 3; ++c)
            origin[c] = detail::read_double(rMesh.FieldData("las:offset"), c);
    for (const double s : scale)
        if (!(s > 0.0) || !std::isfinite(s))
            throw WriteError("meshio++: las: las:scale must be three positive numbers");
    std::string vlrs;
    std::uint16_t encoding = 16;  // WKT, which point formats 6 to 10 require
    if (rMesh.HasFieldData("las:vlrs")) {
        const NDArray& v = rMesh.FieldData("las:vlrs");
        vlrs.assign(reinterpret_cast<const char*>(v.Data()), v.Size());
    }
    std::size_t nvlr = 0;
    for (std::size_t at = 0; at + 54 <= vlrs.size(); ++nvlr)
        at += 54 + las_get<std::uint16_t>(vlrs, at + 20);
    const auto field = [&](const char* pName, double Default, std::size_t I) {
        return has(pName) ? detail::read_double(rMesh.PointData(pName), I) : Default;
    };
    std::array<std::uint64_t, 15> by_return{};
    std::string points;
    points.reserve(n * record);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t c = 0; c < 3; ++c) {
            const double q = std::round((coord(i, c) - origin[c]) / scale[c]);
            if (!(q >= std::numeric_limits<std::int32_t>::min() &&
                  q <= std::numeric_limits<std::int32_t>::max()))
                throw WriteError(
                    "meshio++: las: a coordinate does not fit the int32 grid of "
                    "las:scale and las:offset");
            las_put(points, static_cast<std::int32_t>(q));
        }
        las_put(points, static_cast<std::uint16_t>(field("intensity", 0, i)));
        const auto ret = static_cast<unsigned>(field("return_number", 1, i)) & 15u;
        const auto nret = static_cast<unsigned>(field("number_of_returns", 1, i)) & 15u;
        if (ret >= 1)
            ++by_return[ret - 1];
        points += static_cast<char>(ret | (nret << 4));
        const auto flag = [&](const char* pName, int Shift, unsigned Mask) {
            return (static_cast<unsigned>(field(pName, 0, i)) & Mask) << Shift;
        };
        points += static_cast<char>(
            flag("synthetic", 0, 1) | flag("key_point", 1, 1) | flag("withheld", 2, 1) |
            flag("overlap", 3, 1) | flag("scanner_channel", 4, 3) |
            flag("scan_direction_flag", 6, 1) | flag("edge_of_flight_line", 7, 1));
        points += static_cast<char>(static_cast<unsigned>(field("classification", 0, i)) & 255u);
        points += static_cast<char>(static_cast<unsigned>(field("user_data", 0, i)) & 255u);
        const double angle = has("scan_angle") ? field("scan_angle", 0, i)
                                               : std::round(field("scan_angle_rank", 0, i) / 0.006);
        las_put(points, static_cast<std::int16_t>(angle));
        las_put(points, static_cast<std::uint16_t>(field("point_source_id", 0, i)));
        las_put(points, field("gps_time", 0, i));
        if (format >= 7) {
            las_put(points, static_cast<std::uint16_t>(field("red", 0, i)));
            las_put(points, static_cast<std::uint16_t>(field("green", 0, i)));
            las_put(points, static_cast<std::uint16_t>(field("blue", 0, i)));
        }
        if (format == 8)
            las_put(points, static_cast<std::uint16_t>(field("nir", 0, i)));
    }
    std::string out = "LASF";
    las_put<std::uint16_t>(out, 0);  // file source id
    las_put<std::uint16_t>(out, encoding);
    out.append(16, '\0');  // GUID
    out += '\1';
    out += '\4';
    std::string ident = "meshio++";
    ident.resize(32, '\0');
    out += ident;                    // system identifier
    out += ident;                    // generating software
    las_put<std::uint16_t>(out, 0);  // creation day
    las_put<std::uint16_t>(out, 0);  // creation year
    las_put<std::uint16_t>(out, 375);
    las_put<std::uint32_t>(out, static_cast<std::uint32_t>(375 + vlrs.size()));
    las_put<std::uint32_t>(out, static_cast<std::uint32_t>(nvlr));
    out += static_cast<char>(format);
    las_put<std::uint16_t>(out, static_cast<std::uint16_t>(record));
    las_put<std::uint32_t>(out, 0);  // legacy count: 0 for formats 6 to 10
    for (int k = 0; k < 5; ++k)
        las_put<std::uint32_t>(out, 0);
    for (const double s : scale)
        las_put(out, s);
    for (const double o : origin)
        las_put(out, o);
    for (std::size_t c = 0; c < 3; ++c) {
        las_put(out, hi[c]);
        las_put(out, lo[c]);
    }
    las_put<std::uint64_t>(out, 0);  // waveform data
    las_put<std::uint64_t>(out, 0);  // first EVLR
    las_put<std::uint32_t>(out, 0);  // EVLR count
    las_put<std::uint64_t>(out, n);
    for (const std::uint64_t c : by_return)
        las_put(out, c);
    out += vlrs;
    out += points;
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: las: could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!os)
        throw WriteError("meshio++: las: failed writing " + rPath);
}

}  // namespace meshioplusplus
