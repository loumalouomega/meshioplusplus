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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/vtk_common.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/vtk_xml_pieces.hpp"
#include "../detail/text_cursor.hpp"

namespace meshioplusplus {

namespace {

DType dtype_from_vtk_token(std::string t) {
    for (auto& ch : t)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (t == "float")
        return DType::Float32;
    if (t == "double")
        return DType::Float64;
    if (t == "vtktypeint64" || t == "long")
        return DType::Int64;
    if (t == "vtktypeint8" || t == "char")
        return DType::Int8;
    if (t == "vtktypeint16" || t == "short")
        return DType::Int16;
    if (t == "vtktypeint32" || t == "int")
        return DType::Int32;
    if (t == "vtktypeuint8" || t == "unsigned_char")
        return DType::UInt8;
    if (t == "vtktypeuint16" || t == "unsigned_short")
        return DType::UInt16;
    if (t == "vtktypeuint32" || t == "unsigned_int")
        return DType::UInt32;
    if (t == "vtktypeuint64" || t == "unsigned_long")
        return DType::UInt64;
    throw ReadError("VTK data type '" + t + "' not supported by the C++ reader");
}

void store(NDArray& rA, std::size_t i, double d, std::int64_t v) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[i] = static_cast<float>(d);
            break;
        case DType::Float64:
            rA.As<double>()[i] = d;
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[i] = static_cast<std::int8_t>(v);
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[i] = static_cast<std::int16_t>(v);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = v;
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(v);
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[i] = static_cast<std::uint16_t>(v);
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
            break;
    }
}

struct VtkCursor {
    // A view, not a reference to a std::string: the buffer may be a memory
    // mapping rather than an owned string (see detail/file_source.hpp).
    std::string_view mBuf;
    std::size_t mPos = 0;

    explicit VtkCursor(std::string_view b) : mBuf(b) {}

    bool Eof() const { return mPos >= mBuf.size(); }

    std::string ReadLine() {
        std::size_t start = mPos;
        while (mPos < mBuf.size() && mBuf[mPos] != '\n')
            ++mPos;
        std::string line(mBuf.substr(start, mPos - start));
        if (mPos < mBuf.size())
            ++mPos;  // skip '\n'
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    }

    void ConsumeEol() {
        while (mPos < mBuf.size() && mBuf[mPos] != '\n' &&
               std::isspace(static_cast<unsigned char>(mBuf[mPos])))
            ++mPos;
        if (mPos < mBuf.size() && mBuf[mPos] == '\n')
            ++mPos;
    }

    // Read `count` values of dtype `dt`, ascii or big-endian binary.
    NDArray ReadValues(DType dt, std::size_t count, bool is_ascii) {
        const std::size_t isz = dtype_size(dt);
        // A header count is checked against the bytes left before it sizes
        // anything: an ASCII value takes at least one byte, a binary one isz.
        const std::size_t left = mBuf.size() - std::min(mPos, mBuf.size());
        if (count > left / (is_ascii ? 1 : isz))
            throw ReadError("VTK: an array of " + std::to_string(count) +
                            " values is larger than the rest of the file");
        NDArray a = NDArray::Uninit(dt, {count});  // every element written below
        if (is_ascii) {
            const bool flt = detail::is_float_dtype(dt);
            // strtod/strtoll scan for a terminator: a buffered source is a
            // std::string, and a mapping relies on the kernel's zero-filled
            // final page -- which is why FileSource declines page-multiple
            // sized files.
            const char* base = mBuf.data();
            for (std::size_t i = 0; i < count; ++i) {
                char* endp = nullptr;
                if (flt) {
                    const char* fend = nullptr;
                    double x = detail::parse_double(base + mPos, fend);
                    if (fend == base + mPos)
                        throw ReadError("VTK ascii parse error");
                    store(a, i, x, 0);
                    endp = const_cast<char*>(fend);
                } else {
                    long long x = std::strtoll(base + mPos, &endp, 10);
                    if (endp == base + mPos)
                        throw ReadError("VTK ascii parse error");
                    store(a, i, 0.0, static_cast<std::int64_t>(x));
                }
                mPos = static_cast<std::size_t>(endp - base);
            }
        } else {
            if (mPos + count * isz > mBuf.size())
                throw ReadError("VTK binary truncated");
            char* out = reinterpret_cast<char*>(a.Data());
            // Element offsets are i*isz -> byte-swap in parallel (bswap intrinsic).
            const char* src = mBuf.data() + mPos;
            const int w = static_cast<int>(isz);
            parallel_for_bw(
                count, [&](std::size_t i) { detail::bswap_copy(out + i * isz, src + i * isz, w); });
            mPos += count * isz;
        }
        ConsumeEol();
        return a;
    }
};

std::vector<std::string> split(const std::string& rS) {
    std::vector<std::string> out;
    auto iss = detail::make_classic_istringstream(rS);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}

std::string vtk_upper(std::string s) {
    for (auto& c : s)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::int64_t> vtk_to_int64(const NDArray& rA) {
    std::vector<std::int64_t> v(rA.Size());
    std::int64_t* dst = v.data();
    // Hoist the per-element dtype switch out of the loop, then bulk-convert.
    detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        const T* src = rA.As<T>();
        if constexpr (std::is_floating_point_v<T>) {
            for (std::size_t i = 0; i < rA.Size(); ++i)
                dst[i] = detail::checked_integer<std::int64_t>(src[i], "VTK");
        } else {
            parallel_for_bw(rA.Size(),
                            [&](std::size_t i) { dst[i] = static_cast<std::int64_t>(src[i]); });
        }
    });
    return v;
}

}  // namespace

Mesh read_vtk(const std::string& rPath) try {
    // Mapped where that pays, copied otherwise. Function-local: every parsed
    // value is copied into owning mesh storage, so nothing in the returned Mesh
    // points back into this buffer.
    const detail::FileSource source(rPath);
    VtkCursor cur(source.View());

    std::string header = cur.ReadLine();
    if (header.rfind("# vtk DataFile Version", 0) != 0)
        throw ReadError("Illegal VTK header");
    const bool is_v5 = header.find("Version 5") != std::string::npos;
    cur.ReadLine();  // title
    std::string dtype_line = vtk_upper(cur.ReadLine());
    bool is_ascii;
    if (dtype_line.find("ASCII") != std::string::npos)
        is_ascii = true;
    else if (dtype_line.find("BINARY") != std::string::npos)
        is_ascii = false;
    else
        throw ReadError("Unknown VTK data type line: " + dtype_line);

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    // Held alive so reconstruct_cells can read the int64 connectivity buffer
    // directly (VTK 5.1), skipping a to_int64 copy of the whole connectivity.
    NDArray conn_nd;
    const std::int64_t* conn_ptr = nullptr;
    bool conn_owned = false;  // conn_nd owns the int64 connectivity (VTK 5.1)
    std::unordered_map<std::string, NDArray> cell_data_raw;
    std::string active;  // POINT_DATA or CELL_DATA
    std::string dataset;
    std::array<std::int64_t, 3> dims{};
    std::array<double, 3> origin{}, spacing{};
    std::array<NDArray, 3> axes;
    bool have_dims = false, have_origin = false, have_spacing = false;
    std::size_t active_count = 0;
    std::size_t declared_points = 0, declared_cells = 0;
    bool have_point_data = false, have_cell_data = false;
    auto parse_attribute = [&]<class T>(const std::vector<std::string>& rTokens,
                                        std::array<T, 3>& rValues) {
        detail::need_tokens(rTokens, 4, "VTK");
        for (std::size_t k = 0; k < 3; ++k) {
            detail::TextStream in(rTokens[k + 1]);
            std::string extra;
            if (!(in >> rValues[k]) || (in >> extra))
                throw ReadError("VTK: malformed geometry attribute");
        }
    };
    auto add_array = [&](const std::string& rName, NDArray&& rArr) {
        if (active == "POINT_DATA")
            mesh.AddPointData(rName, std::move(rArr));
        else if (active == "CELL_DATA")
            cell_data_raw.insert_or_assign(rName, std::move(rArr));
        else
            mesh.AddFieldData(rName, std::move(rArr));
    };

    while (!cur.Eof()) {
        std::string line = cur.ReadLine();
        if (line.empty())
            continue;
        std::vector<std::string> tok = split(line);
        if (tok.empty())
            continue;
        std::string section = vtk_upper(tok[0]);

        if (section == "DATASET") {
            detail::need_tokens(tok, 2, "VTK");
            dataset = vtk_upper(tok[1]);
            if (dataset != "UNSTRUCTURED_GRID" && dataset != "STRUCTURED_POINTS" &&
                dataset != "STRUCTURED_GRID" && dataset != "RECTILINEAR_GRID")
                throw ReadError("VTK: unsupported DATASET '" + dataset + "'");
        } else if (section == "DIMENSIONS") {
            parse_attribute(tok, dims);
            for (const auto d : dims)
                if (d < 1)
                    throw ReadError("VTK: DIMENSIONS must be positive");
            have_dims = true;
        } else if (section == "ORIGIN") {
            parse_attribute(tok, origin);
            have_origin = true;
        } else if (section == "SPACING" || section == "ASPECT_RATIO") {
            parse_attribute(tok, spacing);
            have_spacing = true;
        } else if (section == "X_COORDINATES" || section == "Y_COORDINATES" ||
                   section == "Z_COORDINATES") {
            detail::need_tokens(tok, 3, "VTK");
            const auto axis = section[0] == 'X' ? 0 : section[0] == 'Y' ? 1 : 2;
            axes[axis] =
                cur.ReadValues(dtype_from_vtk_token(tok[2]), std::stoull(tok[1]), is_ascii);
        } else if (section == "POINTS") {
            detail::need_tokens(tok, 3, "VTK");
            std::size_t n = std::stoull(tok[1]);
            if (n > std::numeric_limits<std::size_t>::max() / 3)
                throw ReadError("VTK: point count overflows");
            DType dt = dtype_from_vtk_token(tok[2]);
            NDArray pts = cur.ReadValues(dt, n * 3, is_ascii);
            pts.Reshape({n, 3});
            mesh.AssignPoints(std::move(pts));
        } else if (section == "CELLS") {
            detail::need_tokens(tok, 3, "VTK");
            if (is_v5) {
                std::size_t num_off = std::stoull(tok[1]);
                std::size_t num_idx = std::stoull(tok[2]);
                std::string l = cur.ReadLine();
                if (vtk_upper(l).rfind("OFFSETS", 0) != 0)
                    throw ReadError("Expected OFFSETS (VTK 5.1 layout)");
                const std::vector<std::string> otok = split(l);
                detail::need_tokens(otok, 2, "VTK");
                DType odt = dtype_from_vtk_token(otok[1]);
                std::vector<std::int64_t> off_all =
                    vtk_to_int64(cur.ReadValues(odt, num_off, is_ascii));
                l = cur.ReadLine();
                if (vtk_upper(l).rfind("CONNECTIVITY", 0) != 0)
                    throw ReadError("Expected CONNECTIVITY");
                const std::vector<std::string> ctok = split(l);
                detail::need_tokens(ctok, 2, "VTK");
                DType cdt = dtype_from_vtk_token(ctok[1]);
                conn_nd = cur.ReadValues(cdt, num_idx, is_ascii);
                if (conn_nd.Dtype() == DType::Int64) {
                    // Already int64 (vtktypeint64) -> read the buffer directly.
                    conn_ptr = conn_nd.As<std::int64_t>();
                    conn_owned = true;
                } else {
                    conn = vtk_to_int64(conn_nd);
                    conn_ptr = conn.data();
                }
                // off_all has a leading 0; end-offsets are the remainder.
                if (!off_all.empty())
                    offsets.assign(off_all.begin() + 1, off_all.end());
            } else {
                // Version 4.2: interleaved [count, nodes...]; int32 values.
                std::size_t num_cells = std::stoull(tok[1]);
                std::size_t total = std::stoull(tok[2]);
                DType dt = is_ascii ? DType::Int64 : DType::Int32;
                std::vector<std::int64_t> raw = vtk_to_int64(cur.ReadValues(dt, total, is_ascii));
                if (num_cells > total)
                    throw ReadError("VTK: more cells than CELLS entries");
                conn.reserve(total - num_cells);
                offsets.reserve(num_cells);
                std::size_t p = 0;
                std::int64_t running = 0;
                for (std::size_t i = 0; i < num_cells; ++i) {
                    if (p >= raw.size())
                        throw ReadError("VTK: CELLS list ends early");
                    std::int64_t n = raw[p++];
                    if (n < 0 || static_cast<std::uint64_t>(n) > raw.size() - p)
                        throw ReadError("VTK: a cell's node count overruns the CELLS list");
                    for (std::int64_t j = 0; j < n; ++j)
                        conn.push_back(raw[p++]);
                    running += n;
                    offsets.push_back(running);
                }
                conn_ptr = conn.data();
            }
        } else if (section == "CELL_TYPES") {
            detail::need_tokens(tok, 2, "VTK");
            std::size_t n = std::stoull(tok[1]);
            DType dt = is_ascii ? DType::Int64 : DType::Int32;
            types = vtk_to_int64(cur.ReadValues(dt, n, is_ascii));
        } else if (section == "POINT_DATA") {
            detail::need_tokens(tok, 2, "VTK");
            active_count = std::stoull(tok[1]);
            declared_points = active_count;
            have_point_data = true;
            active = "POINT_DATA";
        } else if (section == "CELL_DATA") {
            detail::need_tokens(tok, 2, "VTK");
            active_count = std::stoull(tok[1]);
            declared_cells = active_count;
            have_cell_data = true;
            active = "CELL_DATA";
        } else if (section == "SCALARS" || section == "VECTORS" || section == "TENSORS" ||
                   section == "NORMALS" || section == "TEXTURE_COORDINATES") {
            detail::need_tokens(tok, section == "TEXTURE_COORDINATES" ? 4 : 3, "VTK");
            std::size_t components = section == "TENSORS" ? 9
                                     : section == "SCALARS"
                                         ? (tok.size() > 3 ? std::stoull(tok[3]) : 1)
                                         : 3;
            const DType dt = dtype_from_vtk_token(tok[section == "TEXTURE_COORDINATES" ? 3 : 2]);
            if (section == "TEXTURE_COORDINATES")
                components = std::stoull(tok[2]);
            if (active.empty() || components == 0 ||
                active_count > std::numeric_limits<std::size_t>::max() / components)
                throw ReadError("VTK: invalid attribute count");
            if (section == "SCALARS") {
                const auto lookup = split(cur.ReadLine());
                if (lookup.empty() || vtk_upper(lookup[0]) != "LOOKUP_TABLE")
                    throw ReadError("VTK: SCALARS requires LOOKUP_TABLE");
            }
            auto arr = cur.ReadValues(dt, active_count * components, is_ascii);
            if (section == "TENSORS")
                arr.Reshape({active_count, 3, 3});
            else
                arr.Reshape({active_count, components});
            add_array(tok[1], std::move(arr));
        } else if (section == "COLOR_SCALARS" || section == "LOOKUP_TABLE") {
            detail::need_tokens(tok, 3, "VTK");
            const std::size_t count =
                section == "COLOR_SCALARS" ? active_count : std::stoull(tok[2]);
            const std::size_t components = section == "COLOR_SCALARS" ? std::stoull(tok[2]) : 4;
            if (components == 0 || count > std::numeric_limits<std::size_t>::max() / components)
                throw ReadError("VTK: invalid color array count");
            cur.ReadValues(is_ascii ? DType::Float32 : DType::UInt8, count * components, is_ascii);
        } else if (section == "FIELD") {
            detail::need_tokens(tok, 3, "VTK");
            std::size_t k = std::stoull(tok[2]);
            for (std::size_t fi = 0; fi < k; ++fi) {
                std::vector<std::string> ft = split(cur.ReadLine());
                if (!ft.empty() && vtk_upper(ft[0]) == "METADATA") {
                    while (true) {
                        std::string ml = cur.ReadLine();
                        bool blank = true;
                        for (char c : ml)
                            if (!std::isspace(static_cast<unsigned char>(c)))
                                blank = false;
                        if (blank)
                            break;
                    }
                    ft = split(cur.ReadLine());
                }
                detail::need_tokens(ft, 4, "VTK");
                std::string name = ft[0];
                std::size_t ncomp = std::stoull(ft[1]);
                std::size_t ntuples = std::stoull(ft[2]);
                if (ncomp != 0 && ntuples > std::numeric_limits<std::size_t>::max() / ncomp)
                    throw ReadError("VTK: field array size overflows");
                DType dt = dtype_from_vtk_token(ft[3]);
                NDArray arr = cur.ReadValues(dt, ncomp * ntuples, is_ascii);
                if (ncomp != 1)
                    arr.Reshape({ntuples, ncomp});
                add_array(name, std::move(arr));
            }
        } else if (section == "METADATA") {
            while (true) {
                std::string ml = cur.ReadLine();
                bool blank = true;
                for (char c : ml)
                    if (!std::isspace(static_cast<unsigned char>(c)))
                        blank = false;
                if (blank || cur.Eof())
                    break;
            }
        } else {
            throw ReadError("VTK section '" + section + "' not supported by the C++ reader");
        }
    }

    if (dataset.empty())
        throw ReadError("VTK: missing DATASET");
    if (dataset != "UNSTRUCTURED_GRID") {
        if (!have_dims)
            throw ReadError("VTK: missing DIMENSIONS");
        if (!conn.empty() || !offsets.empty() || !types.empty() || conn_owned)
            throw ReadError("VTK: structured datasets must not declare CELLS/CELL_TYPES");
        std::size_t count = 1;
        for (const auto d : dims) {
            if (static_cast<std::uint64_t>(d) > std::numeric_limits<std::size_t>::max() / count)
                throw ReadError("VTK: dimensions overflow");
            count *= d;
        }
        if (dataset == "STRUCTURED_POINTS" || dataset == "RECTILINEAR_GRID") {
            if (dataset == "STRUCTURED_POINTS" && (!have_origin || !have_spacing))
                throw ReadError("VTK: structured points requires ORIGIN and SPACING/ASPECT_RATIO");
            for (std::size_t axis = 0; axis < 3; ++axis)
                if (dataset == "RECTILINEAR_GRID" &&
                    axes[axis].Size() != static_cast<std::size_t>(dims[axis]))
                    throw ReadError("VTK: coordinate count differs from DIMENSIONS");
            if (count > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) /
                            (3 * sizeof(double)))
                throw ReadError("VTK: structured point allocation overflows");
            auto pts = NDArray::Uninit(
                dataset == "RECTILINEAR_GRID" ? axes[0].Dtype() : DType::Float64, {count, 3});
            auto fill_point = [&](std::size_t row) {
                const std::array<std::int64_t, 3> index{
                    {static_cast<std::int64_t>(row % dims[0]),
                     static_cast<std::int64_t>((row / dims[0]) % dims[1]),
                     static_cast<std::int64_t>(row / dims[0] / dims[1])}};
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    double value;
                    if (dataset == "STRUCTURED_POINTS") {
                        // Match the Python reference's linspace endpoints.
                        const double last =
                            origin[axis] + static_cast<double>(dims[axis] - 1) * spacing[axis];
                        const double step = dims[axis] > 1 ? (last - origin[axis]) /
                                                                 static_cast<double>(dims[axis] - 1)
                                                           : 0;
                        value = index[axis] == dims[axis] - 1
                                    ? last
                                    : origin[axis] + static_cast<double>(index[axis]) * step;
                    } else {
                        value = detail::read_double(axes[axis], index[axis]);
                    }
                    const auto integer = detail::is_float_dtype(pts.Dtype())
                                             ? 0
                                             : detail::checked_integer<std::int64_t>(value, "VTK");
                    store(pts, row * 3 + axis, value, integer);
                }
            };
            if (detail::is_float_dtype(pts.Dtype()))
                parallel_for_bw(count, fill_point);
            else
                // Checked integer conversion may throw: keep it outside workers.
                for (std::size_t row = 0; row < count; ++row)
                    fill_point(row);
            mesh.AssignPoints(std::move(pts));
        }
        if (mesh.NumPoints() != count)
            throw ReadError("VTK: POINTS count differs from DIMENSIONS");
        detail::vtk_structured_cells(dims, conn, offsets, types);
        conn_ptr = conn.data();
    }
    if ((have_point_data && declared_points != mesh.NumPoints()) ||
        (have_cell_data && declared_cells != types.size()))
        throw ReadError("VTK: declared data count differs from geometry");
    for (const auto& name : mesh.PointDataNames()) {
        const auto& arr = mesh.PointData(name);
        if (arr.Shape().empty() || arr.Shape()[0] != mesh.NumPoints())
            throw ReadError("VTK: point data count differs from geometry");
    }
    for (const auto& [name, arr] : cell_data_raw)
        if (arr.Shape().empty() || arr.Shape()[0] != types.size())
            throw ReadError("VTK: cell data count differs from geometry");

    // Fast path (zero copy): a single cell type spanning all cells, non-special,
    // with an identity VTK->meshio node order and regular end-offsets
    // (offsets[i] == (i+1)*n) means the owning int64 connectivity NDArray is
    // already the block data -> reshape and move it straight into the cell block
    // instead of gathering a fresh copy.
    bool moved = false;
    if (conn_owned && !types.empty()) {
        const int vt = static_cast<int>(types[0]);
        bool single = true;
        for (std::size_t i = 1; i < types.size(); ++i)
            if (types[i] != types[0]) {
                single = false;
                break;
            }
        const auto& tmap = vtk_to_meshio_type();
        auto it = tmap.find(vt);
        if (single && it != tmap.end() && !is_special_cell(it->second) &&
            vtk_to_meshio_order(vt).empty()) {
            auto nit = num_nodes_per_cell().find(it->second);
            if (nit != num_nodes_per_cell().end()) {
                const std::size_t n = static_cast<std::size_t>(nit->second);
                const std::size_t ncells = types.size();
                bool regular = conn_nd.Size() == ncells * n && offsets.size() == ncells;
                for (std::size_t i = 0; regular && i < ncells; ++i)
                    if (offsets[i] != static_cast<std::int64_t>((i + 1) * n))
                        regular = false;
                if (regular) {
                    conn_nd.Reshape({ncells, n});
                    mesh.AddCellBlock(it->second, std::move(conn_nd));
                    for (auto& kv : cell_data_raw)
                        mesh.AppendCellData(kv.first, std::move(kv.second));
                    moved = true;
                }
            }
        }
    }
    if (!moved) {
        detail::check_vtk_cell_arrays(conn_owned ? conn_nd.Size() : conn.size(), offsets, types,
                                      cell_data_raw);
        detail::reconstruct_cells(conn_ptr, offsets, types, cell_data_raw, mesh);
    }
    return mesh;
} catch (const ReadError&) {
    throw;
} catch (const std::exception& rError) {
    throw ReadError(std::string("VTK: ") + rError.what());
}

}  // namespace meshioplusplus
