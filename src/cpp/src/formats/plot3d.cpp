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
// PLOT3D multi-block structured grids, solutions and function files: the variant
// inferred from the file's size. See formats/plot3d.hpp.

// System includes
#include <array>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/plot3d.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/lattice_cells.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

constexpr std::int64_t kP3dMaxBlocks = 1000000;
constexpr const char* kP3dDimsKey = "plot3d:dims";
constexpr const char* kP3dIblankKey = "plot3d:iblank";

/// A PLOT3D file as one sequential stream: the tokens of an ASCII file, or the
/// bytes of a binary one (a Fortran file's record payloads concatenated).
struct P3dStream {
    bool mAscii = false;
    bool mFortran = false;
    bool mBig = false;  ///< the payload byte order (a Fortran file's markers say it)
    std::string mBytes;
    std::vector<std::string> mTokens;
    std::vector<std::size_t> mRecords;  ///< a Fortran file's record payload sizes

    /// Units in the stream: tokens, or bytes.
    std::size_t Size() const { return mAscii ? mTokens.size() : mBytes.size(); }
    std::size_t IntUnit() const { return mAscii ? 1 : 4; }
    std::size_t RealUnit(int W) const { return mAscii ? 1 : static_cast<std::size_t>(W); }

    /// The integer at unit offset @p At, or nothing when it is not one.
    std::optional<std::int64_t> IntAt(std::size_t At, bool Big) const {
        if (mAscii) {
            if (At >= mTokens.size())
                return std::nullopt;
            const std::string& t = mTokens[At];
            std::size_t k = (t[0] == '-' || t[0] == '+') ? 1 : 0;
            if (k == t.size())
                return std::nullopt;
            std::int64_t v = 0;
            for (; k < t.size(); ++k) {
                if (!std::isdigit(static_cast<unsigned char>(t[k])) || v > (1LL << 40))
                    return std::nullopt;
                v = v * 10 + (t[k] - '0');
            }
            return t[0] == '-' ? -v : v;
        }
        if (At + 4 > mBytes.size())
            return std::nullopt;
        std::uint32_t u;
        std::memcpy(&u, mBytes.data() + At, 4);
        if (Big != (std::endian::native == std::endian::big))
            u = detail::bswap32(u);
        return static_cast<std::int32_t>(u);
    }

    double RealAt(std::size_t At, int W, bool Big) const {
        if (mAscii) {
            std::string t = mTokens[At];
            for (char& c : t)
                if (c == 'd' || c == 'D')
                    c = 'e';
            const char* end = nullptr;
            const double v = detail::parse_double(t.c_str(), end);
            if (end != t.c_str() + t.size())
                throw ReadError("meshio++: plot3d: '" + mTokens[At] + "' is not a number");
            return v;
        }
        const bool swap = Big != (std::endian::native == std::endian::big);
        if (W == 4) {
            std::uint32_t u;
            std::memcpy(&u, mBytes.data() + At, 4);
            if (swap)
                u = detail::bswap32(u);
            return static_cast<double>(std::bit_cast<float>(u));
        }
        std::uint64_t u;
        std::memcpy(&u, mBytes.data() + At, 8);
        if (swap)
            u = detail::bswap64(u);
        return std::bit_cast<double>(u);
    }
};

bool p3d_looks_ascii(std::string_view Data) {
    if (Data.empty())
        return false;
    const std::size_t n = std::min<std::size_t>(Data.size(), 4096);
    for (std::size_t i = 0; i < n; ++i) {
        const unsigned char c = static_cast<unsigned char>(Data[i]);
        if (!(std::isdigit(c) || std::isspace(c) || c == '.' || c == '-' || c == '+' || c == 'e' ||
              c == 'E' || c == 'd' || c == 'D' || c == ','))
            return false;
    }
    return true;
}

P3dStream p3d_load(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: plot3d: could not open file: " + rPath);
    const std::string_view data = source.View();
    P3dStream s;
    if (p3d_looks_ascii(data)) {
        s.mAscii = true;
        std::size_t i = 0;
        while (i < data.size()) {
            while (i < data.size() &&
                   (std::isspace(static_cast<unsigned char>(data[i])) || data[i] == ','))
                ++i;
            const std::size_t b = i;
            while (i < data.size() && !std::isspace(static_cast<unsigned char>(data[i])) &&
                   data[i] != ',')
                ++i;
            if (i > b)
                s.mTokens.emplace_back(data.substr(b, i - b));
        }
        return s;
    }
    if (const auto layout = detail::sniff_fortran_records(data.data(), data.size())) {
        try {
            const auto records =
                detail::fortran_records(data.data(), data.size(), *layout, "plot3d");
            std::size_t covered = 0;
            for (const auto& r : records)
                covered += r.mSize + 2 * static_cast<std::size_t>(layout->mMarkerBytes);
            if (covered == data.size()) {
                s.mFortran = true;
                s.mBig = layout->mBigEndian;
                for (const auto& r : records) {
                    s.mBytes.append(data.data() + r.mOffset, r.mSize);
                    s.mRecords.push_back(r.mSize);
                }
                return s;
            }
        } catch (const ReadError&) {
            // not a Fortran file after all: a raw stream
        }
    }
    s.mBytes.assign(data.data(), data.size());
    return s;
}

/// One way to read a grid file.
struct P3dShape {
    bool mMulti = false;
    int mDim = 3;
    int mW = 8;
    bool mIblank = false;
    bool mBig = false;
    std::vector<std::array<std::int64_t, 3>> mDims;
    std::size_t mHeader = 0;  ///< units before the first block

    std::size_t Points(std::size_t B) const {
        return static_cast<std::size_t>(mDims[B][0] * mDims[B][1] * mDims[B][2]);
    }
    std::string Describe() const {
        return std::string(mMulti ? "multi-block" : "single-block") + ", " + std::to_string(mDim) +
               "-D, " + (mW == 4 ? "single" : "double") + " precision, " +
               (mIblank ? "with" : "without") + " iblank" + (mBig ? ", big-endian" : "");
    }
};

/// The header of a grid (or `.q`/`.f`) file read as @p Multi / @p Dim, with @p Extra
/// integers per block after its dimensions (`nvar` in a function file).
std::optional<std::pair<std::vector<std::array<std::int64_t, 4>>, std::size_t>> p3d_header(
    const P3dStream& rS, bool Multi, int Dim, bool Big, int Extra) {
    std::size_t at = 0;
    std::int64_t n = 1;
    if (Multi) {
        const auto v = rS.IntAt(0, Big);
        if (!v || *v < 1 || *v > kP3dMaxBlocks)
            return std::nullopt;
        n = *v;
        at += rS.IntUnit();
    }
    std::vector<std::array<std::int64_t, 4>> dims(static_cast<std::size_t>(n), {1, 1, 1, 0});
    for (auto& d : dims) {
        for (int c = 0; c < Dim + Extra; ++c) {
            const auto v = rS.IntAt(at, Big);
            if (!v || *v < 1 || *v > (1LL << 31))
                return std::nullopt;
            (c < Dim ? d[static_cast<std::size_t>(c)] : d[3]) = *v;
            at += rS.IntUnit();
        }
    }
    return std::make_pair(dims, at);
}

std::vector<P3dShape> p3d_grid_shapes(const P3dStream& rS, const Plot3dReadOptions& rO) {
    std::vector<P3dShape> out;
    std::vector<bool> bigs = {rS.mBig};
    if (!rS.mAscii && !rS.mFortran)
        bigs = {false, true};
    for (const bool big : bigs)
        for (const bool multi : {true, false})
            for (const int dim : {3, 2}) {
                if ((rO.mMultiblock >= 0 && multi != (rO.mMultiblock == 1)) ||
                    (rO.mDimensions != 0 && dim != rO.mDimensions))
                    continue;
                const auto header = p3d_header(rS, multi, dim, big, 0);
                if (!header)
                    continue;
                for (const int w : rS.mAscii ? std::vector<int>{8} : std::vector<int>{4, 8})
                    for (const bool iblank : {false, true}) {
                        if ((rO.mPrecision != 0 && !rS.mAscii && w != rO.mPrecision) ||
                            (rO.mIblank >= 0 && iblank != (rO.mIblank == 1)))
                            continue;
                        P3dShape s;
                        s.mMulti = multi;
                        s.mDim = dim;
                        s.mW = w;
                        s.mIblank = iblank;
                        s.mBig = big;
                        s.mHeader = header->second;
                        for (const auto& d : header->first)
                            s.mDims.push_back({d[0], d[1], d[2]});
                        // The exact size: a Fortran file's records, else its length.
                        std::size_t total = s.mHeader;
                        std::vector<std::size_t> records;
                        if (multi)
                            records.push_back(4);
                        records.push_back(4 * static_cast<std::size_t>(dim) * s.mDims.size());
                        bool too_big = false;
                        for (std::size_t b = 0; b < s.mDims.size() && !too_big; ++b) {
                            const std::size_t np = s.Points(b);
                            if (np > rS.Size()) {
                                too_big = true;
                                break;
                            }
                            const std::size_t units =
                                np * (static_cast<std::size_t>(dim) * rS.RealUnit(w) +
                                      (iblank ? rS.IntUnit() : 0));
                            total += units;
                            records.push_back(units);
                        }
                        if (too_big || total != rS.Size())
                            continue;
                        if (rS.mFortran && records != rS.mRecords)
                            continue;
                        out.push_back(std::move(s));
                    }
            }
    return out;
}

/// Whether a reading's values look like a grid: finite coordinates of sane magnitude
/// (no subnormals, which a misread precision produces) and small iblank values. Used
/// only to break a tie between readings the size alone cannot tell apart.
bool p3d_plausible(const P3dStream& rS, const P3dShape& rShape) {
    std::size_t at = rShape.mHeader;
    for (std::size_t b = 0; b < rShape.mDims.size(); ++b) {
        const std::size_t np = rShape.Points(b);
        for (std::size_t k = 0; k < np * static_cast<std::size_t>(rShape.mDim);
             ++k, at += rS.RealUnit(rShape.mW)) {
            double v = 0.0;
            try {
                v = rS.RealAt(at, rShape.mW, rShape.mBig);
            } catch (const ReadError&) {
                return false;
            }
            if (!std::isfinite(v) || std::fabs(v) > 1e30 || (v != 0.0 && std::fabs(v) < 1e-30))
                return false;
        }
        if (rShape.mIblank)
            for (std::size_t p = 0; p < np; ++p, at += rS.IntUnit()) {
                const auto v = rS.IntAt(at, rShape.mBig);
                if (!v || *v < -kP3dMaxBlocks || *v > kP3dMaxBlocks)
                    return false;
            }
    }
    return true;
}

P3dShape p3d_pick_grid(const P3dStream& rS, const Plot3dReadOptions& rO, const std::string& rPath) {
    std::vector<P3dShape> shapes = p3d_grid_shapes(rS, rO);
    if (shapes.size() > 1) {
        std::vector<P3dShape> plausible;
        for (const auto& s : shapes)
            if (p3d_plausible(rS, s))
                plausible.push_back(s);
        if (!plausible.empty())
            shapes = std::move(plausible);
    }
    if (shapes.empty())
        throw ReadError("meshio++: plot3d: '" + rPath +
                        "' is not a PLOT3D grid this reader can size (ASCII, Fortran "
                        "unformatted or raw; single or multi-block; 2-D or 3-D; single or "
                        "double precision; with or without iblank)");
    if (shapes.size() > 1) {
        std::string what;
        for (const auto& s : shapes)
            what += "\n  - " + s.Describe();
        throw ReadError("meshio++: plot3d: '" + rPath + "' reads " + std::to_string(shapes.size()) +
                        " ways; pick one with the precision, multiblock, iblank or dimensions "
                        "read option:" +
                        what);
    }
    return shapes.front();
}

/// The `.q` solution: per block four reference values then `dim + 2` variables.
void p3d_add_solution(Mesh& rMesh, const P3dShape& rGrid, const std::string& rPath) {
    const P3dStream s = p3d_load(rPath);
    const std::size_t nvar = static_cast<std::size_t>(rGrid.mDim) + 2;
    std::vector<bool> bigs = {s.mBig};
    if (!s.mAscii && !s.mFortran)
        bigs = {rGrid.mBig, !rGrid.mBig};
    for (const bool big : bigs) {
        const auto header = p3d_header(s, rGrid.mMulti, rGrid.mDim, big, 0);
        if (!header || header->first.size() != rGrid.mDims.size())
            continue;
        bool same = true;
        for (std::size_t b = 0; b < rGrid.mDims.size(); ++b)
            for (int c = 0; c < 3; ++c)
                same = same && header->first[b][static_cast<std::size_t>(c)] ==
                                   rGrid.mDims[b][static_cast<std::size_t>(c)];
        if (!same)
            continue;
        for (const int w : s.mAscii ? std::vector<int>{8} : std::vector<int>{4, 8}) {
            std::size_t total = header->second;
            for (std::size_t b = 0; b < rGrid.mDims.size(); ++b)
                total += (4 + rGrid.Points(b) * nvar) * s.RealUnit(w);
            if (total != s.Size())
                continue;
            std::size_t np_all = 0;
            for (std::size_t b = 0; b < rGrid.mDims.size(); ++b)
                np_all += rGrid.Points(b);
            NDArray rho(DType::Float64, {np_all});
            NDArray mom(DType::Float64, {np_all, 3});
            NDArray energy(DType::Float64, {np_all});
            NDArray ref(DType::Float64, {rGrid.mDims.size(), 4});
            std::size_t at = header->second, first = 0;
            for (std::size_t b = 0; b < rGrid.mDims.size(); ++b) {
                for (std::size_t r = 0; r < 4; ++r, at += s.RealUnit(w))
                    ref.As<double>()[b * 4 + r] = s.RealAt(at, w, big);
                const std::size_t np = rGrid.Points(b);
                for (std::size_t v = 0; v < nvar; ++v)
                    for (std::size_t p = 0; p < np; ++p, at += s.RealUnit(w)) {
                        const double x = s.RealAt(at, w, big);
                        if (v == 0)
                            rho.As<double>()[first + p] = x;
                        else if (v == nvar - 1)
                            energy.As<double>()[first + p] = x;
                        else
                            mom.As<double>()[(first + p) * 3 + (v - 1)] = x;
                    }
                first += np;
            }
            rMesh.AddPointData("density", std::move(rho));
            rMesh.AddPointData("momentum", std::move(mom));
            rMesh.AddPointData("energy", std::move(energy));
            rMesh.AddFieldData("plot3d:q_reference", std::move(ref));
            return;
        }
    }
    throw ReadError("meshio++: plot3d: '" + rPath +
                    "' is not a PLOT3D solution for this grid (block count, dimensions or "
                    "size disagree)");
}

/// The `.f` function file: per block its dimensions and `nvar`, then the values.
void p3d_add_function(Mesh& rMesh, const P3dShape& rGrid, const std::string& rPath) {
    const P3dStream s = p3d_load(rPath);
    std::vector<bool> bigs = {s.mBig};
    if (!s.mAscii && !s.mFortran)
        bigs = {rGrid.mBig, !rGrid.mBig};
    for (const bool big : bigs) {
        const auto header = p3d_header(s, rGrid.mMulti, rGrid.mDim, big, 1);
        if (!header || header->first.size() != rGrid.mDims.size())
            continue;
        bool same = true;
        const std::int64_t nvar = header->first.front()[3];
        for (std::size_t b = 0; b < rGrid.mDims.size(); ++b) {
            for (int c = 0; c < 3; ++c)
                same = same && header->first[b][static_cast<std::size_t>(c)] ==
                                   rGrid.mDims[b][static_cast<std::size_t>(c)];
            same = same && header->first[b][3] == nvar;
        }
        if (!same)
            continue;
        const std::size_t nv = static_cast<std::size_t>(nvar);
        for (const int w : s.mAscii ? std::vector<int>{8} : std::vector<int>{4, 8}) {
            std::size_t total = header->second, np_all = 0;
            for (std::size_t b = 0; b < rGrid.mDims.size(); ++b) {
                total += rGrid.Points(b) * nv * s.RealUnit(w);
                np_all += rGrid.Points(b);
            }
            if (total != s.Size())
                continue;
            NDArray f(DType::Float64, nv == 1 ? std::vector<std::size_t>{np_all}
                                              : std::vector<std::size_t>{np_all, nv});
            std::size_t at = header->second, first = 0;
            for (std::size_t b = 0; b < rGrid.mDims.size(); ++b) {
                const std::size_t np = rGrid.Points(b);
                for (std::size_t v = 0; v < nv; ++v)
                    for (std::size_t p = 0; p < np; ++p, at += s.RealUnit(w))
                        f.As<double>()[(first + p) * nv + v] = s.RealAt(at, w, big);
                first += np;
            }
            rMesh.AddPointData("plot3d:f", std::move(f));
            return;
        }
    }
    throw ReadError("meshio++: plot3d: '" + rPath +
                    "' is not a PLOT3D function file for this grid (block count, dimensions, "
                    "variable count or size disagree)");
}

}  // namespace

Mesh read_plot3d(const std::string& rPath, const Plot3dReadOptions& rOptions) {
    const P3dStream s = p3d_load(rPath);
    const P3dShape shape = p3d_pick_grid(s, rOptions, rPath);
    std::size_t np_all = 0;
    for (std::size_t b = 0; b < shape.mDims.size(); ++b)
        np_all += shape.Points(b);
    NDArray pts(DType::Float64, {np_all, 3});
    NDArray iblank(DType::Int32, {shape.mIblank ? np_all : 0});
    std::size_t at = shape.mHeader, first = 0;
    for (std::size_t b = 0; b < shape.mDims.size(); ++b) {
        const std::size_t np = shape.Points(b);
        for (int c = 0; c < shape.mDim; ++c)
            for (std::size_t p = 0; p < np; ++p, at += s.RealUnit(shape.mW))
                pts.As<double>()[(first + p) * 3 + static_cast<std::size_t>(c)] =
                    s.RealAt(at, shape.mW, shape.mBig);
        if (shape.mIblank)
            for (std::size_t p = 0; p < np; ++p, at += s.IntUnit()) {
                const auto v = s.IntAt(at, shape.mBig);
                if (!v)
                    throw ReadError("meshio++: plot3d: an iblank value is not an integer");
                iblank.As<std::int32_t>()[first + p] = static_cast<std::int32_t>(*v);
            }
        first += np;
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    NDArray dims(DType::Int64, {shape.mDims.size(), 3});
    std::size_t base = 0, cell_base = 0;
    for (std::size_t b = 0; b < shape.mDims.size(); ++b) {
        for (std::size_t c = 0; c < 3; ++c)
            dims.As<std::int64_t>()[b * 3 + c] = shape.mDims[b][c];
        std::string type;
        std::vector<std::int64_t> conn;
        std::size_t npc = 0;
        detail::lattice_cells(shape.mDims[b][0], shape.mDims[b][1], shape.mDims[b][2], type, conn,
                              npc);
        const std::size_t nc = npc ? conn.size() / npc : 0;
        if (nc > 0) {
            NDArray block(DType::Int64, {nc, npc});
            for (std::size_t k = 0; k < conn.size(); ++k)
                block.As<std::int64_t>()[k] = conn[k] + static_cast<std::int64_t>(base);
            mesh.AddCellBlock(type, std::move(block));
            NDArray ids(DType::Int64, {nc});
            for (std::size_t k = 0; k < nc; ++k)
                ids.As<std::int64_t>()[k] = static_cast<std::int64_t>(cell_base + k);
            mesh.AddRegion(Region("block_" + std::to_string(b), RegionKind::Cell, std::move(ids)));
            cell_base += nc;
        }
        base += shape.Points(b);
    }
    mesh.AddFieldData(kP3dDimsKey, std::move(dims));
    if (shape.mIblank)
        mesh.AddPointData(kP3dIblankKey, std::move(iblank));
    if (!rOptions.mSolution.empty())
        p3d_add_solution(mesh, shape, rOptions.mSolution);
    if (!rOptions.mFunction.empty())
        p3d_add_function(mesh, shape, rOptions.mFunction);
    return mesh;
}

Mesh read_plot3d(const std::string& rPath) {
    return read_plot3d(rPath, Plot3dReadOptions{});
}

void write_plot3d(const std::string& rPath, const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kP3dDimsKey))
        throw WriteError(
            "meshio++: plot3d: the mesh is not a set of PLOT3D blocks (it has no "
            "field_data['plot3d:dims']); read one from a PLOT3D file, or write a structured "
            "format such as vts");
    const NDArray& dims = rMesh.FieldData(kP3dDimsKey);
    if (dims.Size() == 0 || dims.Size() % 3 != 0)
        throw WriteError("meshio++: plot3d: field_data['plot3d:dims'] must be (nblocks, 3)");
    const std::size_t nb = dims.Size() / 3;
    std::size_t total = 0;
    std::vector<std::size_t> npts(nb);
    for (std::size_t b = 0; b < nb; ++b) {
        std::int64_t prod = 1;
        for (std::size_t c = 0; c < 3; ++c) {
            const std::int64_t v = detail::read_int(dims, b * 3 + c);
            if (v < 1 || v > std::numeric_limits<std::int32_t>::max())
                throw WriteError("meshio++: plot3d: a block dimension is out of range");
            prod *= v;
        }
        npts[b] = static_cast<std::size_t>(prod);
        total += npts[b];
    }
    if (total != rMesh.NumPoints())
        throw WriteError("meshio++: plot3d: the blocks in field_data['plot3d:dims'] hold " +
                         std::to_string(total) + " points but the mesh has " +
                         std::to_string(rMesh.NumPoints()));
    const NDArray* iblank =
        rMesh.HasPointData(kP3dIblankKey) ? &rMesh.PointData(kP3dIblankKey) : nullptr;
    if (iblank && iblank->Size() != total)
        throw WriteError("meshio++: plot3d: plot3d:iblank needs one value per point");
    const detail::FortranRecordLayout layout{4, false};
    std::string out, rec;
    detail::fortran_put_int32(rec, static_cast<std::int32_t>(nb), false);
    detail::fortran_append_record(out, rec.data(), rec.size(), layout);
    rec.clear();
    for (std::size_t k = 0; k < nb * 3; ++k)
        detail::fortran_put_int32(rec, static_cast<std::int32_t>(detail::read_int(dims, k)), false);
    detail::fortran_append_record(out, rec.data(), rec.size(), layout);
    const NDArray& pts = rMesh.Points();
    const std::size_t pd = rMesh.PointDim();
    std::size_t first = 0;
    for (std::size_t b = 0; b < nb; ++b) {
        rec.clear();
        for (std::size_t c = 0; c < 3; ++c)
            for (std::size_t p = 0; p < npts[b]; ++p)
                detail::fortran_put_float64(
                    rec, c < pd ? detail::read_double(pts, (first + p) * pd + c) : 0.0, false);
        if (iblank)
            for (std::size_t p = 0; p < npts[b]; ++p)
                detail::fortran_put_int32(
                    rec, static_cast<std::int32_t>(detail::read_int(*iblank, first + p)), false);
        detail::fortran_append_record(out, rec.data(), rec.size(), layout);
        first += npts[b];
    }
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: plot3d: could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!os)
        throw WriteError("meshio++: plot3d: failed writing " + rPath);
}

}  // namespace meshioplusplus
