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

#pragma once

/**
 * @file detail/typed_view.hpp
 * @brief An array's elements as `int64_t` or `double`, with the dtype switch
 * taken once per array instead of once per element.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * `Int64View` reads every element exactly as `detail::read_int` does, and
 * `DoubleView` exactly as `detail::read_double` does -- a pointer into the
 * array itself when its dtype already is the target, otherwise one converted
 * copy made in parallel through `dispatch_dtype`. A hot loop then indexes a
 * plain pointer. Roadmap §4, "Hoist the dtype switch in operations".
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// One element as `read_int` converts it.
template <class T>
inline std::int64_t typed_view_int(T v) {
    if constexpr (std::is_floating_point_v<T>) {
        const double d = static_cast<double>(v);
        if (!(d == d))
            return 0;
        if (d >= 9223372036854775807.0)
            return std::numeric_limits<std::int64_t>::max();
        if (d < -9223372036854775808.0)
            return std::numeric_limits<std::int64_t>::min();
        return static_cast<std::int64_t>(d);
    } else {
        return static_cast<std::int64_t>(v);
    }
}

/// An array read as `int64_t`, as `read_int` reads each element.
class Int64View {
public:
    explicit Int64View(const NDArray& rA) {
        const std::size_t n = rA.Size();
        if (rA.Dtype() == DType::Int64) {
            mpData = rA.As<std::int64_t>();
            return;
        }
        mOwned.resize(n);
        dispatch_dtype(rA.Dtype(), [&]<class T>() {
            const T* src = rA.As<T>();
            std::int64_t* dst = mOwned.data();
            parallel_for_bw(n, [&](std::size_t i) { dst[i] = typed_view_int<T>(src[i]); });
        });
        mpData = mOwned.data();
    }
    Int64View(const Int64View&) = delete;
    Int64View& operator=(const Int64View&) = delete;

    const std::int64_t* Data() const { return mpData; }
    std::int64_t operator[](std::size_t i) const { return mpData[i]; }

private:
    const std::int64_t* mpData = nullptr;
    std::vector<std::int64_t> mOwned;
};

/// An array read as `double`, as `read_double` reads each element.
class DoubleView {
public:
    explicit DoubleView(const NDArray& rA) {
        const std::size_t n = rA.Size();
        if (rA.Dtype() == DType::Float64) {
            mpData = rA.As<double>();
            return;
        }
        mOwned.resize(n);
        dispatch_dtype(rA.Dtype(), [&]<class T>() {
            const T* src = rA.As<T>();
            double* dst = mOwned.data();
            parallel_for_bw(n, [&](std::size_t i) { dst[i] = static_cast<double>(src[i]); });
        });
        mpData = mOwned.data();
    }
    DoubleView(const DoubleView&) = delete;
    DoubleView& operator=(const DoubleView&) = delete;

    const double* Data() const { return mpData; }
    double operator[](std::size_t i) const { return mpData[i]; }

private:
    const double* mpData = nullptr;
    std::vector<double> mOwned;
};

}  // namespace detail
}  // namespace meshioplusplus
