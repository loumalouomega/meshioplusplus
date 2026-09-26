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
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <system_error>
#include <type_traits>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// strtoll/strtoull over [pFirst, pLast): an optional sign, decimal digits,
// saturation on overflow and strtoull's modular negation, but never a read
// past the token. The strto* calls it replaces take a NUL-terminated string,
// and ASan's strict_string_checks makes every such call check the whole rest
// of it, so a DataArray of n integers cost O(n^2) under the sanitizers.
// Returns pFirst when there is no digit to convert.
template <class TInt>
const char* vtkxml_parse_decimal(const char* pFirst, const char* pLast, TInt& rValue) {
    const char* p = pFirst;
    const bool negative = p != pLast && *p == '-';
    if (p != pLast && (*p == '+' || *p == '-'))
        ++p;
    if (p == pLast || !std::isdigit(static_cast<unsigned char>(*p)))
        return pFirst;
    std::uint64_t magnitude = 0;
    const auto [end, ec] = std::from_chars(p, pLast, magnitude);
    const bool overflow = ec == std::errc::result_out_of_range;
    if constexpr (std::is_signed_v<TInt>) {
        const std::uint64_t limit =
            static_cast<std::uint64_t>(std::numeric_limits<TInt>::max()) + (negative ? 1 : 0);
        if (overflow || magnitude > limit)
            rValue = negative ? std::numeric_limits<TInt>::min() : std::numeric_limits<TInt>::max();
        else
            rValue = static_cast<TInt>(negative ? std::uint64_t{0} - magnitude : magnitude);
    } else {
        if (overflow)
            rValue = std::numeric_limits<TInt>::max();
        else
            rValue = static_cast<TInt>(negative ? std::uint64_t{0} - magnitude : magnitude);
    }
    return end;
}

}  // namespace

const char* vtu_type_str(DType dt) {
    switch (dt) {
        case DType::Float32:
            return "Float32";
        case DType::Float64:
            return "Float64";
        case DType::Int8:
            return "Int8";
        case DType::Int16:
            return "Int16";
        case DType::Int32:
            return "Int32";
        case DType::Int64:
            return "Int64";
        case DType::UInt8:
            return "UInt8";
        case DType::UInt16:
            return "UInt16";
        case DType::UInt32:
            return "UInt32";
        case DType::UInt64:
            return "UInt64";
    }
    return "Float64";
}

DType dtype_from_vtu(const std::string& rS) {
    if (rS == "Float32")
        return DType::Float32;
    if (rS == "Float64")
        return DType::Float64;
    if (rS == "Int8")
        return DType::Int8;
    if (rS == "Int16")
        return DType::Int16;
    if (rS == "Int32")
        return DType::Int32;
    if (rS == "Int64")
        return DType::Int64;
    if (rS == "UInt8")
        return DType::UInt8;
    if (rS == "UInt16")
        return DType::UInt16;
    if (rS == "UInt32")
        return DType::UInt32;
    if (rS == "UInt64")
        return DType::UInt64;
    throw ReadError("Illegal VTU data type '" + rS + "'");
}

void vtu_ascii_double(std::ostream& rOs, double v) {
    char buf[32];
    detail::snprintf_c(buf, sizeof(buf), "%.11e", v);
    rOs << buf << '\n';
}

void vtu_ascii_ndarray(std::ostream& rOs, const NDArray& rA) {
    const bool flt = is_float_dtype(rA.Dtype());
    const std::size_t n = rA.Size();
    if (rA.Dtype() == DType::UInt64) {
        // read_int would hand a value above INT64_MAX back as a negative one.
        const std::uint64_t* p = rA.As<std::uint64_t>();
        for (std::size_t i = 0; i < n; ++i)
            rOs << p[i] << '\n';
        return;
    }
    if (flt) {
        const DoubleView v(rA);
        for (std::size_t i = 0; i < n; ++i)
            vtu_ascii_double(rOs, v[i]);
    } else {
        const Int64View v(rA);
        for (std::size_t i = 0; i < n; ++i)
            rOs << v[i] << '\n';
    }
}

void vtu_store(NDArray& rA, std::size_t i, double d, std::int64_t v) {
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

NDArray vtu_parse_ascii(const char* pText, DType dt) {
    const bool isflt = is_float_dtype(dt);
    std::vector<double> dv;
    std::vector<std::int64_t> iv;
    const char* p = pText ? pText : "";
    const char* const last = p + std::strlen(p);
    while (*p) {
        while (*p && std::isspace(static_cast<unsigned char>(*p)))
            ++p;
        if (!*p)
            break;
        char* endp = nullptr;
        if (isflt) {
            const char* fend = nullptr;
            double x = detail::parse_double(p, fend);
            if (fend == p)
                break;
            dv.push_back(x);
            endp = const_cast<char*>(fend);
        } else if (dt == DType::UInt64) {
            // Unsigned: a signed parse saturates above INT64_MAX; the bit pattern round-trips
            // through the int64 buffer and vtu_store's cast back.
            std::uint64_t x = 0;
            endp = const_cast<char*>(vtkxml_parse_decimal(p, last, x));
            if (endp == p)
                break;
            iv.push_back(static_cast<std::int64_t>(x));
        } else {
            std::int64_t x = 0;
            endp = const_cast<char*>(vtkxml_parse_decimal(p, last, x));
            if (endp == p)
                break;
            iv.push_back(x);
        }
        p = endp;
    }
    std::size_t n = isflt ? dv.size() : iv.size();
    NDArray a(dt, {n});
    for (std::size_t i = 0; i < n; ++i)
        vtu_store(a, i, isflt ? dv[i] : 0.0, isflt ? 0 : iv[i]);
    return a;
}

std::string vtu_strip(const char* pS) {
    std::string t = pS ? pS : "";
    std::size_t b = 0, e = t.size();
    while (b < e && std::isspace(static_cast<unsigned char>(t[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(t[e - 1])))
        --e;
    return t.substr(b, e - b);
}

NDArray vtu_parse_binary(const std::string& rText, DType dt, VtkCodec codec, std::size_t hsz) {
    std::vector<unsigned char> bytes;
    if (codec == VtkCodec::None)
        bytes = vtu_decode_uncompressed(rText.c_str(), rText.size(), hsz);
    else
        bytes = vtu_decode_blocks(rText.c_str(), rText.size(), hsz, codec);
    std::size_t isz = dtype_size(dt);
    std::size_t n = isz ? bytes.size() / isz : 0;
    NDArray a(dt, {n});
    if (n)
        std::memcpy(a.Data(), bytes.data(), n * isz);
    return a;
}

namespace {

constexpr const char* kVtuGhostName = "vtkGhostType";

}  // namespace

DType vtu_disk_dtype(const std::string& rName, DType Dt) {
    return rName == kVtuGhostName ? DType::UInt8 : Dt;
}

const NDArray& vtu_disk_array(const std::string& rName, const NDArray& rArray, NDArray& rScratch) {
    if (rName != kVtuGhostName || rArray.Dtype() == DType::UInt8)
        return rArray;
    rScratch = NDArray::Uninit(DType::UInt8, rArray.Shape());
    std::uint8_t* out = rScratch.As<std::uint8_t>();
    const Int64View v(rArray);
    parallel_for_bw(rArray.Size(),
                    [&](std::size_t i) { out[i] = static_cast<std::uint8_t>(v[i]); });
    return rScratch;
}

void vtu_write_field_array(std::ostream& rOs, const std::string& rName, const NDArray& rArray,
                           bool Binary, VtkCodec Codec) {
    vtu_write_field_array(rOs, rName, rArray, Binary, Codec, 4);
}

void vtu_write_field_array(std::ostream& rOs, const std::string& rName, const NDArray& rArray,
                           bool Binary, VtkCodec Codec, std::size_t Hsz) {
    const std::vector<std::size_t>& shape = rArray.Shape();
    const std::size_t tuples = shape.empty() ? 1 : shape[0];
    rOs << "<DataArray type=\"" << vtu_type_str(rArray.Dtype()) << "\" Name=\"" << rName
        << "\" NumberOfTuples=\"" << tuples << "\"";
    if (shape.size() >= 2) {
        std::size_t components = 1;
        for (std::size_t d = 1; d < shape.size(); ++d)
            components *= shape[d];
        rOs << " NumberOfComponents=\"" << components << "\"";
    }
    rOs << " format=\"" << (Binary ? "binary" : "ascii") << "\">\n";
    if (Binary)
        rOs << vtu_encode_binary(reinterpret_cast<const unsigned char*>(rArray.Data()),
                                 rArray.Nbytes(), Codec, Hsz)
            << "\n";
    else
        vtu_ascii_ndarray(rOs, rArray);
    rOs << "</DataArray>\n";
}

std::vector<std::int64_t> vtu_to_int64(const NDArray& rA) {
    const Int64View view(rA);
    return std::vector<std::int64_t>(view.Data(), view.Data() + rA.Size());
}

}  // namespace detail
}  // namespace meshioplusplus
