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
// FreeSurfer surfaces, curv overlays and annot parcellations. See
// formats/freesurfer.hpp.

// System includes
#include <bit>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/freesurfer.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

/// Big-endian reads over a file, bounds-checked.
struct FsCursor {
    std::string mData;
    std::size_t mAt = 0;
    std::string mWhat;

    void Need(std::size_t N) const {
        if (mData.size() - mAt < N)
            throw ReadError("meshio++: freesurfer: " + mWhat + " is truncated");
    }
    std::int32_t I32() {
        Need(4);
        std::uint32_t u;
        std::memcpy(&u, mData.data() + mAt, 4);
        mAt += 4;
        if constexpr (std::endian::native == std::endian::little)
            u = detail::bswap32(u);
        return static_cast<std::int32_t>(u);
    }
    float F32() {
        const std::int32_t i = I32();
        return std::bit_cast<float>(i);
    }
    std::uint32_t Magic() {
        Need(3);
        const auto* b = reinterpret_cast<const unsigned char*>(mData.data() + mAt);
        mAt += 3;
        return (std::uint32_t{b[0]} << 16) | (std::uint32_t{b[1]} << 8) | b[2];
    }
};

FsCursor fs_open(const std::string& rPath, const std::string& rWhat) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: freesurfer: could not open file: " + rPath);
    FsCursor c;
    c.mData.assign(source.View());
    c.mWhat = rWhat + " '" + rPath + "'";
    return c;
}

NDArray fs_bytes(std::string_view Text) {
    NDArray out(DType::UInt8, {Text.size()});
    std::memcpy(out.Data(), Text.data(), Text.size());
    return out;
}

std::string fs_overlay_name(const std::string& rPath) {
    std::string base = std::filesystem::path(rPath).filename().string();
    if (base.rfind("lh.", 0) == 0 || base.rfind("rh.", 0) == 0)
        base = base.substr(3);
    return base.empty() ? std::string("curv") : base;
}

void fs_add_curv(Mesh& rMesh, const std::string& rPath) {
    FsCursor c = fs_open(rPath, "curv file");
    if (c.Magic() != 0xFFFFFFu)
        throw ReadError("meshio++: freesurfer: '" + rPath +
                        "' is not a new-format curv file (old-format curv is not read)");
    const std::int32_t nv = c.I32();
    c.I32();  // faces
    const std::int32_t per = c.I32();
    if (nv < 0 || static_cast<std::size_t>(nv) != rMesh.NumPoints())
        throw ReadError("meshio++: freesurfer: '" + rPath + "' has " + std::to_string(nv) +
                        " values but the surface has " + std::to_string(rMesh.NumPoints()) +
                        " vertices");
    if (per != 1)
        throw ReadError("meshio++: freesurfer: '" + rPath + "' has " + std::to_string(per) +
                        " values per vertex; meshio++ reads one");
    NDArray values(DType::Float32, {static_cast<std::size_t>(nv)});
    for (std::int32_t k = 0; k < nv; ++k)
        values.As<float>()[k] = c.F32();
    rMesh.AddPointData(fs_overlay_name(rPath), std::move(values));
}

std::string fs_string(FsCursor& rC) {
    const std::int32_t len = rC.I32();
    if (len < 0)
        throw ReadError("meshio++: freesurfer: a negative string length in " + rC.mWhat);
    rC.Need(static_cast<std::size_t>(len));
    std::string s(rC.mData.data() + rC.mAt, static_cast<std::size_t>(len));
    rC.mAt += static_cast<std::size_t>(len);
    const auto nul = s.find('\0');
    return nul == std::string::npos ? s : s.substr(0, nul);
}

void fs_add_annot(Mesh& rMesh, const std::string& rPath) {
    FsCursor c = fs_open(rPath, "annot file");
    const std::int32_t nv = c.I32();
    if (nv < 0 || static_cast<std::size_t>(nv) != rMesh.NumPoints())
        throw ReadError("meshio++: freesurfer: '" + rPath + "' labels " + std::to_string(nv) +
                        " vertices but the surface has " + std::to_string(rMesh.NumPoints()));
    std::vector<std::int32_t> colour(static_cast<std::size_t>(nv), 0);
    for (std::int32_t k = 0; k < nv; ++k) {
        const std::int32_t v = c.I32();
        const std::int32_t label = c.I32();
        if (v < 0 || v >= nv)
            throw ReadError("meshio++: freesurfer: '" + rPath + "' names vertex " +
                            std::to_string(v));
        colour[static_cast<std::size_t>(v)] = label;
    }
    struct Entry {
        std::int32_t mIndex;
        std::string mName;
        std::int32_t mColour;
    };
    std::vector<Entry> entries;
    if (c.mAt < c.mData.size() && c.I32() == 1) {
        const std::int32_t n = c.I32();
        const auto rgb = [&] {
            const std::int32_t r = c.I32(), g = c.I32(), b = c.I32();
            c.I32();  // transparency
            return r + g * 256 + b * 65536;
        };
        if (n > 0) {
            fs_string(c);  // the colour table's file name
            const std::size_t count =
                detail::checked_count(n, c.mData.size(), "freesurfer", "entry");
            for (std::size_t k = 0; k < count; ++k) {
                const std::string name = fs_string(c);
                entries.push_back({static_cast<std::int32_t>(k), name, rgb()});
            }
        } else {
            if (-n != 2)
                throw ReadError("meshio++: freesurfer: colour table version " + std::to_string(-n) +
                                " is not read");
            c.I32();       // the largest structure index
            fs_string(c);  // the colour table's file name
            const std::size_t count =
                detail::checked_count(c.I32(), c.mData.size(), "freesurfer", "entry");
            for (std::size_t k = 0; k < count; ++k) {
                const std::int32_t index = c.I32();
                const std::string name = fs_string(c);
                entries.push_back({index, name, rgb()});
            }
        }
    }
    std::unordered_map<std::int32_t, std::int32_t> by_colour;
    for (const Entry& e : entries)
        by_colour.emplace(e.mColour, e.mIndex);
    NDArray index(DType::Int32, {static_cast<std::size_t>(nv)});
    for (std::int32_t k = 0; k < nv; ++k) {
        const auto it = by_colour.find(colour[static_cast<std::size_t>(k)]);
        index.As<std::int32_t>()[k] = it == by_colour.end() ? -1 : it->second;
    }
    for (const Entry& e : entries) {
        std::vector<std::int64_t> ids;
        for (std::int32_t k = 0; k < nv; ++k)
            if (index.As<std::int32_t>()[k] == e.mIndex)
                ids.push_back(k);
        NDArray arr(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), arr.As<std::int64_t>());
        rMesh.AddRegion(Region(e.mName, RegionKind::Point, -1, e.mIndex, std::move(arr)));
    }
    rMesh.AddPointData("freesurfer:annot", std::move(index));
}

}  // namespace

Mesh read_freesurfer(const std::string& rPath, const FreesurferReadOptions& rOptions) {
    FsCursor c = fs_open(rPath, "surface");
    const std::uint32_t magic = c.Magic();
    if (magic == 0xFFFFFFu || magic == 0xFFFFFDu)
        throw ReadError("meshio++: freesurfer: '" + rPath +
                        "' is a quad surface; meshio++ reads triangle surfaces");
    if (magic != 0xFFFFFEu)
        throw ReadError("meshio++: freesurfer: '" + rPath + "' is not a FreeSurfer surface");
    const std::size_t end = c.mData.find("\n\n", c.mAt);
    if (end == std::string::npos)
        throw ReadError("meshio++: freesurfer: '" + rPath + "' has no creator line");
    const std::string created = c.mData.substr(c.mAt, end - c.mAt);
    c.mAt = end + 2;
    const std::int32_t nv = c.I32(), nf = c.I32();
    const std::size_t limit = c.mData.size();
    const std::size_t np = detail::checked_count(nv, limit / 12, "freesurfer", "vertex");
    const std::size_t nc = detail::checked_count(nf, limit / 12, "freesurfer", "face");
    NDArray pts(DType::Float64, {np, 3});
    for (std::size_t k = 0; k < np * 3; ++k)
        pts.As<double>()[k] = static_cast<double>(c.F32());
    NDArray tri(DType::Int64, {nc, 3});
    for (std::size_t k = 0; k < nc * 3; ++k) {
        const std::int32_t v = c.I32();
        if (v < 0 || static_cast<std::size_t>(v) >= np)
            throw ReadError("meshio++: freesurfer: a face names vertex " + std::to_string(v));
        tri.As<std::int64_t>()[k] = v;
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    if (nc > 0)
        mesh.AddCellBlock("triangle", std::move(tri));
    mesh.AddFieldData("freesurfer:created_by", fs_bytes(created));
    if (c.mAt < c.mData.size())
        mesh.AddFieldData("freesurfer:tail", fs_bytes(std::string_view(c.mData).substr(c.mAt)));
    for (const std::string& curv : rOptions.mCurv)
        fs_add_curv(mesh, curv);
    if (!rOptions.mAnnot.empty())
        fs_add_annot(mesh, rOptions.mAnnot);
    return mesh;
}

Mesh read_freesurfer(const std::string& rPath) {
    return read_freesurfer(rPath, FreesurferReadOptions{});
}

void write_freesurfer(const std::string& rPath, const Mesh& rMesh) {
    if (rMesh.NumCellBlocks() != 1 || std::string(rMesh.Cells(0).Type()) != "triangle")
        throw WriteError("meshio++: freesurfer: a FreeSurfer surface is one triangle block");
    const auto cb = rMesh.Cells(0);
    std::string out = "\xFF\xFF\xFE";
    std::string created = "created by meshio++";
    if (rMesh.HasFieldData("freesurfer:created_by")) {
        const NDArray& t = rMesh.FieldData("freesurfer:created_by");
        created.assign(reinterpret_cast<const char*>(t.Data()), t.Size());
    }
    out += created + "\n\n";
    const auto put = [&](std::uint32_t u) {
        if constexpr (std::endian::native == std::endian::little)
            u = detail::bswap32(u);
        out.append(reinterpret_cast<const char*>(&u), 4);
    };
    const std::size_t np = rMesh.NumPoints(), pd = rMesh.PointDim();
    put(static_cast<std::uint32_t>(np));
    put(static_cast<std::uint32_t>(cb.NumCells()));
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            put(std::bit_cast<std::uint32_t>(static_cast<float>(
                c < pd ? detail::read_double(rMesh.Points(), i * pd + c) : 0.0)));
    for (std::size_t k = 0; k < cb.NumCells() * 3; ++k)
        put(static_cast<std::uint32_t>(detail::read_int(cb.Conn(), k)));
    if (rMesh.HasFieldData("freesurfer:tail")) {
        const NDArray& t = rMesh.FieldData("freesurfer:tail");
        out.append(reinterpret_cast<const char*>(t.Data()), t.Size());
    }
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: freesurfer: could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!os)
        throw WriteError("meshio++: freesurfer: failed writing " + rPath);
}

}  // namespace meshioplusplus
