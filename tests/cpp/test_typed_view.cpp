// License: MIT License, meshio++ default license: LICENSE
// Core-private numeric views preserve the scalar dispatch and caller storage.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/typed_view.hpp"

namespace {
namespace tvd = meshioplusplus::detail;
using meshioplusplus::DType;
using meshioplusplus::NDArray;
constexpr DType tvd_dtypes[] = {DType::Float32, DType::Float64, DType::Int8,  DType::Int16,
                                DType::Int32,   DType::Int64,   DType::UInt8, DType::UInt16,
                                DType::UInt32,  DType::UInt64};
}  // namespace

TEST(TypedViews, AllDtypesMatchScalarReadsWithoutMutatingTheInput) {
    for (const auto dtype : tvd_dtypes) {
        NDArray array(dtype, {9});
        tvd::dispatch_dtype(dtype, [&]<class T>() {
            auto* values = array.As<T>();
            if constexpr (std::is_floating_point_v<T>) {
                const double input[] = {-0.0,
                                        0.0,
                                        1.5,
                                        -2.75,
                                        1e40,
                                        -1e40,
                                        std::numeric_limits<double>::infinity(),
                                        -std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::quiet_NaN()};
                for (std::size_t i = 0; i < 9; ++i)
                    values[i] = static_cast<T>(input[i]);
            } else {
                values[0] = 0;
                values[1] = 1;
                values[2] = 2;
                values[3] = static_cast<T>(-1);
                values[4] = std::numeric_limits<T>::max();
                values[5] = std::numeric_limits<T>::min();
                values[6] = static_cast<T>(std::numeric_limits<T>::max() - 1);
                values[7] = 17;
                values[8] = 23;
            }
        });
        const std::vector<std::byte> before(array.Data(), array.Data() + array.Nbytes());
        const tvd::Int64View ints(array);
        const tvd::DoubleView reals(array);
        if (dtype == DType::Int64)
            EXPECT_EQ(ints.Data(), array.As<std::int64_t>());
        if (dtype == DType::Float64)
            EXPECT_EQ(reals.Data(), array.As<double>());
        for (std::size_t i = 0; i < array.Size(); ++i) {
            EXPECT_EQ(ints[i], tvd::read_int(array, i));
            EXPECT_EQ(std::bit_cast<std::uint64_t>(reals[i]),
                      std::bit_cast<std::uint64_t>(tvd::read_double(array, i)));
        }
        EXPECT_EQ(before, std::vector<std::byte>(array.Data(), array.Data() + array.Nbytes()));
    }
}

TEST(TypedViews, EmptyAndScalarStorageKeepTheirSizes) {
    const NDArray absent;
    const tvd::Int64View absent_ints(absent);
    const tvd::DoubleView absent_reals(absent);
    for (const auto dtype : tvd_dtypes) {
        const NDArray empty(dtype, {0});
        const tvd::Int64View empty_ints(empty);
        const tvd::DoubleView empty_reals(empty);
        NDArray scalar(dtype, {});
        ASSERT_EQ(scalar.Size(), 1u);
        tvd::dispatch_dtype(dtype, [&]<class T>() { scalar.As<T>()[0] = static_cast<T>(7); });
        const tvd::Int64View ints(scalar);
        const tvd::DoubleView reals(scalar);
        EXPECT_EQ(ints[0], 7);
        EXPECT_EQ(reals[0], 7.0);
    }
}
