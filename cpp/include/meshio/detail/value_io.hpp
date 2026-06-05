#pragma once
//
// Shared helpers to read scalar values out of an NDArray regardless of its
// dtype, used by the ASCII format writers.

#include <cstddef>
#include <cstdint>

#include "meshio/ndarray.hpp"

namespace meshio {
namespace detail {

inline bool is_float_dtype(DType dt) {
    return dt == DType::Float32 || dt == DType::Float64;
}

inline double read_double(const NDArray& a, std::size_t i) {
    switch (a.dtype()) {
        case DType::Float32: return static_cast<double>(a.as<float>()[i]);
        case DType::Float64: return a.as<double>()[i];
        case DType::Int8: return static_cast<double>(a.as<std::int8_t>()[i]);
        case DType::Int16: return static_cast<double>(a.as<std::int16_t>()[i]);
        case DType::Int32: return static_cast<double>(a.as<std::int32_t>()[i]);
        case DType::Int64: return static_cast<double>(a.as<std::int64_t>()[i]);
        case DType::UInt8: return static_cast<double>(a.as<std::uint8_t>()[i]);
        case DType::UInt16: return static_cast<double>(a.as<std::uint16_t>()[i]);
        case DType::UInt32: return static_cast<double>(a.as<std::uint32_t>()[i]);
        case DType::UInt64: return static_cast<double>(a.as<std::uint64_t>()[i]);
    }
    return 0.0;
}

inline std::int64_t read_int(const NDArray& a, std::size_t i) {
    switch (a.dtype()) {
        case DType::Int8: return a.as<std::int8_t>()[i];
        case DType::Int16: return a.as<std::int16_t>()[i];
        case DType::Int32: return a.as<std::int32_t>()[i];
        case DType::Int64: return a.as<std::int64_t>()[i];
        case DType::UInt8: return a.as<std::uint8_t>()[i];
        case DType::UInt16: return a.as<std::uint16_t>()[i];
        case DType::UInt32: return a.as<std::uint32_t>()[i];
        case DType::UInt64: return static_cast<std::int64_t>(a.as<std::uint64_t>()[i]);
        default: return static_cast<std::int64_t>(read_double(a, i));
    }
}

inline std::size_t rows(const NDArray& a) {
    return a.shape().empty() ? 0 : a.shape()[0];
}

inline std::size_t cols(const NDArray& a) {
    return a.shape().size() >= 2 ? a.shape()[1] : 1;
}

}  // namespace detail
}  // namespace meshio
