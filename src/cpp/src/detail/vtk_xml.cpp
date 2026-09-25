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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <system_error>
#include <type_traits>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"

namespace meshioplusplus {
namespace detail {

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
    for (std::size_t i = 0; i < n; ++i) {
        if (flt)
            vtu_ascii_double(rOs, read_double(rA, i));
        else
            rOs << read_int(rA, i) << '\n';
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

namespace {

/**
 * @brief Parse one base-10 integer at @p p with `std::strtoll`/`strtoull`'s
 * acceptance -- an optional sign, saturation on overflow, and for an unsigned
 * target a negated magnitude wrapping as `strtoull` does -- but through
 * `std::from_chars`, bounded by @p pEnd. `strto*` finds the end of the string
 * first under instrumentation (ASan's `strict_string_checks`), which made
 * parsing one ASCII `DataArray` quadratic in its length.
 * @return Past the number, or @p p when there is none.
 */
template <class T>
const char* vxml_parse_integer(const char* p, const char* pEnd, T& rOut) {
    const char* q = p;
    bool neg = false;
    if (q < pEnd && (*q == '+' || *q == '-')) {
        neg = *q == '-';
        ++q;
    }
    unsigned long long mag = 0;
    const auto [ptr, ec] = std::from_chars(q, pEnd, mag);
    if (ptr == q)
        return p;
    const bool overflow = ec == std::errc::result_out_of_range;
    if constexpr (std::is_signed_v<T>) {
        constexpr unsigned long long kMax = std::numeric_limits<long long>::max();
        if (!neg)
            rOut = overflow || mag > kMax ? std::numeric_limits<long long>::max()
                                          : static_cast<long long>(mag);
        else  // -(kMax + 1) is exactly the minimum; anything larger saturates to it
            rOut = overflow || mag > kMax ? std::numeric_limits<long long>::min()
                                          : -static_cast<long long>(mag);
    } else {
        rOut = overflow ? std::numeric_limits<unsigned long long>::max() : (neg ? 0ull - mag : mag);
    }
    return ptr;
}

}  // namespace

NDArray vtu_parse_ascii(const char* pText, DType dt) {
    const bool isflt = is_float_dtype(dt);
    std::vector<double> dv;
    std::vector<std::int64_t> iv;
    const char* p = pText ? pText : "";
    const char* const text_end = p + std::strlen(p);
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
            // strtoll saturates above INT64_MAX; the bit pattern round-trips
            // through the int64 buffer and vtu_store's cast back.
            unsigned long long x = 0;
            const char* e = vxml_parse_integer(p, text_end, x);
            if (e == p)
                break;
            iv.push_back(static_cast<std::int64_t>(x));
            endp = const_cast<char*>(e);
        } else {
            long long x = 0;
            const char* e = vxml_parse_integer(p, text_end, x);
            if (e == p)
                break;
            iv.push_back(static_cast<std::int64_t>(x));
            endp = const_cast<char*>(e);
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
    for (std::size_t i = 0; i < rArray.Size(); ++i)
        out[i] = static_cast<std::uint8_t>(read_int(rArray, i));
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
    std::vector<std::int64_t> v(rA.Size());
    for (std::size_t i = 0; i < rA.Size(); ++i)
        v[i] = read_int(rA, i);
    return v;
}

}  // namespace detail
}  // namespace meshioplusplus
