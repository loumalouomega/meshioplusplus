// License: MIT License, meshio++ default license: LICENSE
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/permas.hpp"

namespace {
namespace tid = meshioplusplus;
constexpr tid::DType tid_dtypes[] = {tid::DType::Float32, tid::DType::Float64, tid::DType::Int8,
                                     tid::DType::Int16,   tid::DType::Int32,   tid::DType::Int64,
                                     tid::DType::UInt8,   tid::DType::UInt16,  tid::DType::UInt32,
                                     tid::DType::UInt64};

tid::NDArray tid_array(tid::DType Dtype, const std::vector<std::size_t>& rShape,
                       const std::vector<double>& rValues) {
    tid::NDArray out(Dtype, rShape);
    tid::detail::dispatch_dtype(Dtype, [&]<class T>() {
        for (std::size_t i = 0; i < rValues.size(); ++i)
            out.As<T>()[i] = static_cast<T>(rValues[i]);
    });
    return out;
}

tid::Mesh tid_mesh(tid::DType RealType, tid::DType IndexType, std::size_t Dim) {
    tid::Mesh mesh;
    mesh.AssignPoints(tid_array(RealType, {3, Dim},
                                Dim == 2 ? std::vector<double>{0, 0, 1, 0, 0, 1}
                                         : std::vector<double>{0, 0, 0, 1, 0, 0, 0, 1, 0}));
    mesh.AddCellBlock("triangle", tid_array(IndexType, {1, 3}, {0, 1, 2}));
    mesh.AddPointData("s", tid_array(RealType, {3}, {1, 2, 3}));
    mesh.AddPointData("v", tid_array(RealType, {3, 2}, {1, 2, 3, 4, 5, 6}));
    mesh.AddCellData("pf3:ref", {tid_array(IndexType, {1}, {7})});
    return mesh;
}

std::string tid_bytes(const std::string& rPath) {
    auto in = tid::detail::make_classic_ifstream(rPath, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
}  // namespace

TEST(TextIoDtypes, FourWritersMatchCanonicalStorageOnAllDtypes) {
    using Writer = void (*)(const std::string&, const tid::Mesh&);
    const std::array<std::pair<const char*, Writer>, 4> writers = {{{"off", tid::write_off},
                                                                    {"ip", tid::write_ip},
                                                                    {"flux", tid::write_flux},
                                                                    {"permas", tid::write_permas}}};
    for (const auto dtype : tid_dtypes) {
        for (const std::size_t dim : {2u, 3u}) {
            const auto input = tid_mesh(dtype, dtype, dim);
            const auto expected = tid_mesh(tid::DType::Float64, tid::DType::Int64, dim);
            for (const auto& [format, writer] : writers) {
                SCOPED_TRACE(format);
                SCOPED_TRACE(static_cast<int>(dtype));
                SCOPED_TRACE(dim);
                const std::string first = mt::temp_path(".data");
                const std::string second = mt::temp_path(".data");
                writer(first, input);
                writer(second, expected);
                EXPECT_EQ(tid_bytes(first), tid_bytes(second));
                std::filesystem::remove(first);
                std::filesystem::remove(second);
            }
        }
    }
}
