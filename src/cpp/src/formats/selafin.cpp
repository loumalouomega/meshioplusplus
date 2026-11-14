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
// TELEMAC Selafin results. See formats/selafin.hpp.

// System includes
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/selafin.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "../detail/open_source.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kSlfUnitPrefix = "selafin:unit:";

struct SlfRecords {
    std::string mBytes;
    std::vector<detail::FortranRecord> mRecords;
    std::size_t mAt = 0;

    std::string_view Next(const char* pWhat) {
        if (mAt >= mRecords.size())
            throw ReadError(std::string("meshio++: selafin: the file ends before its ") + pWhat);
        const auto& r = mRecords[mAt++];
        return std::string_view(mBytes.data() + r.mOffset, r.mSize);
    }
};

std::int32_t slf_i32(std::string_view Rec, std::size_t K) {
    std::uint32_t u;
    std::memcpy(&u, Rec.data() + K * 4, 4);
    if constexpr (std::endian::native == std::endian::little)
        u = detail::bswap32(u);
    return static_cast<std::int32_t>(u);
}

double slf_real(std::string_view Rec, std::size_t K, std::size_t W) {
    if (W == 4) {
        std::uint32_t u;
        std::memcpy(&u, Rec.data() + K * 4, 4);
        if constexpr (std::endian::native == std::endian::little)
            u = detail::bswap32(u);
        return static_cast<double>(std::bit_cast<float>(u));
    }
    std::uint64_t u;
    std::memcpy(&u, Rec.data() + K * 8, 8);
    if constexpr (std::endian::native == std::endian::little)
        u = detail::bswap64(u);
    return std::bit_cast<double>(u);
}

std::string slf_trim(std::string_view Text) {
    std::size_t e = Text.size();
    while (e > 0 && (Text[e - 1] == ' ' || Text[e - 1] == '\0'))
        --e;
    std::size_t b = 0;
    while (b < e && Text[b] == ' ')
        ++b;
    return std::string(Text.substr(b, e - b));
}

NDArray slf_bytes(const std::string& rText) {
    NDArray out(DType::UInt8, {rText.size()});
    if (!rText.empty())
        std::memcpy(out.Data(), rText.data(), rText.size());
    return out;
}

struct SlfHeader {
    std::string mTitle;
    std::vector<std::string> mNames, mUnits;
    std::int64_t mNplan = 0;
    std::size_t mNelem = 0, mNpoin = 0, mNdp = 0, mW = 4;
    std::vector<std::int64_t> mIkle;
    std::vector<std::int32_t> mIpobo;
    std::vector<double> mX, mY;
    std::size_t mFirstStep = 0;  ///< record index of the first step's time
    std::vector<double> mTimes;
};

SlfRecords slf_load(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: selafin: could not open file: " + rPath);
    SlfRecords s;
    s.mBytes.assign(source.View());
    const detail::FortranRecordLayout layout{4, true};
    s.mRecords = detail::fortran_records(s.mBytes.data(), s.mBytes.size(), layout, "selafin",
                                         /*Lenient=*/true);
    return s;
}

SlfHeader slf_header(SlfRecords& rS) {
    SlfHeader h;
    const std::string_view title = rS.Next("title");
    if (title.size() != 80)
        throw ReadError(
            "meshio++: selafin: the title record is not 80 characters (not a Selafin "
            "file, or not big-endian)");
    h.mTitle = slf_trim(title);
    const std::string_view nbv = rS.Next("variable counts");
    if (nbv.size() != 8)
        throw ReadError("meshio++: selafin: the NBV1 NBV2 record is malformed");
    const std::size_t nvar = detail::checked_count(slf_i32(nbv, 0) + std::int64_t{slf_i32(nbv, 1)},
                                                   rS.mRecords.size(), "selafin", "variable");
    for (std::size_t v = 0; v < nvar; ++v) {
        const std::string_view rec = rS.Next("variable name");
        if (rec.size() != 32)
            throw ReadError("meshio++: selafin: a variable record is not 32 characters");
        h.mNames.push_back(slf_trim(rec.substr(0, 16)));
        h.mUnits.push_back(slf_trim(rec.substr(16, 16)));
    }
    const std::string_view iparam = rS.Next("IPARAM");
    if (iparam.size() != 40)
        throw ReadError("meshio++: selafin: IPARAM is not ten integers");
    h.mNplan = slf_i32(iparam, 6);
    if (slf_i32(iparam, 9) == 1)
        rS.Next("date");
    const std::string_view counts = rS.Next("NELEM NPOIN NDP");
    if (counts.size() != 16)
        throw ReadError("meshio++: selafin: the NELEM NPOIN NDP 1 record is malformed");
    const std::size_t limit = rS.mBytes.size();
    h.mNelem = detail::checked_count(slf_i32(counts, 0), limit, "selafin", "element");
    h.mNpoin = detail::checked_count(slf_i32(counts, 1), limit, "selafin", "point");
    h.mNdp = detail::checked_count(slf_i32(counts, 2), 64, "selafin", "node-per-element");
    const std::string_view ikle = rS.Next("IKLE");
    if (ikle.size() != h.mNelem * h.mNdp * 4)
        throw ReadError("meshio++: selafin: IKLE does not hold NELEM * NDP integers");
    h.mIkle.resize(h.mNelem * h.mNdp);
    for (std::size_t k = 0; k < h.mIkle.size(); ++k) {
        const std::int64_t v = slf_i32(ikle, k);
        if (v < 1 || static_cast<std::size_t>(v) > h.mNpoin)
            throw ReadError("meshio++: selafin: IKLE names point " + std::to_string(v));
        h.mIkle[k] = v - 1;
    }
    const std::string_view ipobo = rS.Next("IPOBO");
    if (ipobo.size() != h.mNpoin * 4)
        throw ReadError("meshio++: selafin: IPOBO does not hold NPOIN integers");
    h.mIpobo.resize(h.mNpoin);
    for (std::size_t k = 0; k < h.mNpoin; ++k)
        h.mIpobo[k] = slf_i32(ipobo, k);
    const std::string_view x = rS.Next("X");
    if (h.mNpoin > 0 && x.size() == h.mNpoin * 8)
        h.mW = 8;
    else if (x.size() != h.mNpoin * 4)
        throw ReadError("meshio++: selafin: X does not hold NPOIN reals");
    const std::string_view y = rS.Next("Y");
    if (y.size() != x.size())
        throw ReadError("meshio++: selafin: Y does not match X");
    h.mX.resize(h.mNpoin);
    h.mY.resize(h.mNpoin);
    for (std::size_t k = 0; k < h.mNpoin; ++k) {
        h.mX[k] = slf_real(x, k, h.mW);
        h.mY[k] = slf_real(y, k, h.mW);
    }
    h.mFirstStep = rS.mAt;
    // Every complete step: a time record and one record per variable.
    for (std::size_t r = rS.mAt; r < rS.mRecords.size(); r += 1 + nvar) {
        if (rS.mRecords[r].mSize != h.mW)
            break;
        bool complete = true;
        for (std::size_t v = 0; v < nvar && complete; ++v)
            complete =
                r + 1 + v < rS.mRecords.size() && rS.mRecords[r + 1 + v].mSize == h.mNpoin * h.mW;
        if (!complete)
            break;
        h.mTimes.push_back(
            slf_real(std::string_view(rS.mBytes.data() + rS.mRecords[r].mOffset, h.mW), 0, h.mW));
    }
    return h;
}

}  // namespace

Mesh read_selafin(const std::string& rPath, const ReadOptions& rOptions) {
    SlfRecords s = slf_load(rPath);
    const SlfHeader h = slf_header(s);
    const std::size_t step = rOptions.ResolveTimeStep(h.mTimes.size());
    NDArray pts(DType::Float64, {h.mNpoin, 3});
    for (std::size_t k = 0; k < h.mNpoin; ++k) {
        pts.As<double>()[k * 3] = h.mX[k];
        pts.As<double>()[k * 3 + 1] = h.mY[k];
        pts.As<double>()[k * 3 + 2] = 0.0;
    }
    std::string type;
    if (h.mNdp == 3)
        type = "triangle";
    else if (h.mNdp == 4)
        type = h.mNplan > 1 ? "tetra" : "quad";
    else if (h.mNdp == 6)
        type = "wedge";
    else
        throw ReadError("meshio++: selafin: elements of " + std::to_string(h.mNdp) +
                        " nodes are not a TELEMAC element");
    std::vector<NDArray> values;
    if (!h.mTimes.empty()) {
        const std::size_t base = h.mFirstStep + step * (1 + h.mNames.size()) + 1;
        for (std::size_t v = 0; v < h.mNames.size(); ++v) {
            const auto& r = s.mRecords[base + v];
            const std::string_view rec(s.mBytes.data() + r.mOffset, r.mSize);
            NDArray a(DType::Float64, {h.mNpoin});
            for (std::size_t k = 0; k < h.mNpoin; ++k)
                a.As<double>()[k] = slf_real(rec, k, h.mW);
            values.push_back(std::move(a));
        }
        // TELEMAC-3D keeps the plane elevations as its first variable.
        if (h.mNplan > 1)
            for (std::size_t v = 0; v < h.mNames.size(); ++v)
                if (h.mNames[v].rfind("ELEVATION Z", 0) == 0 ||
                    h.mNames[v].rfind("COTE Z", 0) == 0) {
                    for (std::size_t k = 0; k < h.mNpoin; ++k)
                        pts.As<double>()[k * 3 + 2] = values[v].As<double>()[k];
                    break;
                }
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    if (h.mNelem > 0) {
        NDArray conn(DType::Int64, {h.mNelem, h.mNdp});
        std::copy(h.mIkle.begin(), h.mIkle.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(type, std::move(conn));
    }
    mesh.AddFieldData("selafin:title", slf_bytes(h.mTitle));
    if (!h.mTimes.empty()) {
        NDArray t(DType::Float64, {1});
        t.As<double>()[0] = h.mTimes[step];
        mesh.AddFieldData("meshio:time", std::move(t));
    }
    if (!rOptions.WantsAnyData())
        return mesh;
    for (std::size_t v = 0; v < values.size(); ++v) {
        if (!rOptions.WantsArray(h.mNames[v]))
            continue;
        mesh.AddPointData(h.mNames[v], std::move(values[v]));
        mesh.AddFieldData(kSlfUnitPrefix + h.mNames[v], slf_bytes(h.mUnits[v]));
    }
    if (rOptions.WantsArray("selafin:ipobo")) {
        NDArray ipobo(DType::Int32, {h.mNpoin});
        std::copy(h.mIpobo.begin(), h.mIpobo.end(), ipobo.As<std::int32_t>());
        mesh.AddPointData("selafin:ipobo", std::move(ipobo));
    }
    return mesh;
}

MeshMetadata read_selafin_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    MeshMetadata meta = metadata_from_mesh(read_selafin(rPath, rOptions));
    meta.mFellBackToFullRead = true;
    SlfRecords s = slf_load(rPath);
    meta.mTimeValues = slf_header(s).mTimes;
    return meta;
}

namespace {

std::string slf_padded(const std::string& rText, std::size_t Width) {
    std::string out = rText.substr(0, Width);
    out.resize(Width, ' ');
    return out;
}

}  // namespace

void write_selafin(const std::string& rPath, const Mesh& rMesh) {
    if (rMesh.NumCellBlocks() != 1 || (std::string(rMesh.Cells(0).Type()) != "triangle" &&
                                       std::string(rMesh.Cells(0).Type()) != "quad"))
        throw WriteError(
            "meshio++: selafin: the writer takes one triangle or quad block (a 2-D "
            "TELEMAC mesh)");
    const auto cb = rMesh.Cells(0);
    const std::size_t npoin = rMesh.NumPoints(), nelem = cb.NumCells(), ndp = cb.NodesPerCell();
    std::vector<std::string> names;
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (name == "selafin:ipobo" || a.Size() != npoin ||
            !(a.Dtype() == DType::Float64 || a.Dtype() == DType::Float32))
            continue;
        names.push_back(name);
    }
    const detail::FortranRecordLayout layout{4, true};
    std::string out, rec;
    const auto flush = [&] {
        detail::fortran_append_record(out, rec.data(), rec.size(), layout);
        rec.clear();
    };
    std::string title = "meshio++";
    if (rMesh.HasFieldData("selafin:title")) {
        const NDArray& t = rMesh.FieldData("selafin:title");
        title.assign(reinterpret_cast<const char*>(t.Data()), t.Size());
    }
    rec = slf_padded(title, 72) + "SERAFIND";
    flush();
    detail::fortran_put_int32(rec, static_cast<std::int32_t>(names.size()), true);
    detail::fortran_put_int32(rec, 0, true);
    flush();
    for (const std::string& name : names) {
        std::string unit;
        const std::string key = kSlfUnitPrefix + name;
        if (rMesh.HasFieldData(key)) {
            const NDArray& u = rMesh.FieldData(key);
            unit.assign(reinterpret_cast<const char*>(u.Data()), u.Size());
        }
        rec = slf_padded(name, 16) + slf_padded(unit, 16);
        flush();
    }
    for (int k = 0; k < 10; ++k)
        detail::fortran_put_int32(rec, k == 0 ? 1 : 0, true);
    flush();
    for (const std::size_t v : {nelem, npoin, ndp, std::size_t{1}})
        detail::fortran_put_int32(rec, static_cast<std::int32_t>(v), true);
    flush();
    const detail::Int64View conn(cb.Conn());
    for (std::size_t k = 0; k < nelem * ndp; ++k)
        detail::fortran_put_int32(rec, static_cast<std::int32_t>(conn[k] + 1), true);
    flush();
    const NDArray* ipobo =
        rMesh.HasPointData("selafin:ipobo") ? &rMesh.PointData("selafin:ipobo") : nullptr;
    for (std::size_t k = 0; k < npoin; ++k)
        detail::fortran_put_int32(rec,
                                  ipobo && ipobo->Size() == npoin
                                      ? static_cast<std::int32_t>(detail::read_int(*ipobo, k))
                                      : 0,
                                  true);
    flush();
    const std::size_t pd = rMesh.PointDim();
    for (std::size_t c = 0; c < 2; ++c) {
        for (std::size_t k = 0; k < npoin; ++k)
            detail::fortran_put_float64(
                rec, c < pd ? detail::read_double(rMesh.Points(), k * pd + c) : 0.0, true);
        flush();
    }
    const double time = rMesh.HasFieldData("meshio:time") && rMesh.FieldData("meshio:time").Size()
                            ? detail::read_double(rMesh.FieldData("meshio:time"), 0)
                            : 0.0;
    detail::fortran_put_float64(rec, time, true);
    flush();
    for (const std::string& name : names) {
        const NDArray& a = rMesh.PointData(name);
        for (std::size_t k = 0; k < npoin; ++k)
            detail::fortran_put_float64(rec, detail::read_double(a, k), true);
        flush();
    }
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: selafin: could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!os)
        throw WriteError("meshio++: selafin: failed writing " + rPath);
}

}  // namespace meshioplusplus
