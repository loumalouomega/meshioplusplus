// License: MIT License, meshio++ default license: LICENSE
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/febio.hpp"
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/patran.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/z88.hpp"

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

/// The mesh each writer is given: a triangle with point and cell data, plus
/// the integer cell data the format reads per cell (tags, types, properties).
/// Z88 has no linear-triangle element, so it gets a tetrahedron.
tid::Mesh tid_mesh(const std::string& rFormat, tid::DType RealType, tid::DType IndexType,
                   std::size_t Dim) {
    tid::Mesh mesh;
    if (rFormat == "z88") {
        mesh.AssignPoints(tid_array(RealType, {4, 3}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}));
        mesh.AddCellBlock("tetra", tid_array(IndexType, {1, 4}, {0, 1, 2, 3}));
        mesh.AddCellData("z88:type", {tid_array(IndexType, {1}, {17})});
        return mesh;
    }
    mesh.AssignPoints(tid_array(RealType, {3, Dim},
                                Dim == 2 ? std::vector<double>{0, 0, 1, 0, 0, 1}
                                         : std::vector<double>{0, 0, 0, 1, 0, 0, 0, 1, 0}));
    mesh.AddCellBlock("triangle", tid_array(IndexType, {1, 3}, {0, 1, 2}));
    mesh.AddPointData("s", tid_array(RealType, {3}, {1, 2, 3}));
    mesh.AddPointData("v", tid_array(RealType, {3, 2}, {1, 2, 3, 4, 5, 6}));
    mesh.AddCellData("pf3:ref", {tid_array(IndexType, {1}, {7})});
    if (rFormat.rfind("gmsh", 0) == 0) {
        mesh.AddCellData("gmsh:physical", {tid_array(IndexType, {1}, {3})});
        mesh.AddCellData("gmsh:geometrical", {tid_array(IndexType, {1}, {5})});
    } else if (rFormat == "mdpa") {
        mesh.AddCellData("gmsh:physical", {tid_array(IndexType, {1}, {3})});
    } else if (rFormat == "patran") {
        mesh.AddCellData("patran:property", {tid_array(IndexType, {1}, {4})});
    } else if (rFormat == "femap") {
        mesh.AddCellData("femap:property", {tid_array(IndexType, {1}, {4})});
    }
    return mesh;
}

std::string tid_bytes(const std::string& rPath) {
    auto in = tid::detail::make_classic_ifstream(rPath, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

using TidWriter = std::function<void(const std::string&, const tid::Mesh&)>;

const std::vector<std::pair<std::string, TidWriter>>& tid_writers() {
    static const std::vector<std::pair<std::string, TidWriter>> writers = {
        {"off", tid::write_off},
        {"ip", tid::write_ip},
        {"flux", tid::write_flux},
        {"permas", tid::write_permas},
        {"gmsh22-ascii", [](const std::string& rPath,
                            const tid::Mesh& rMesh) { tid::write_gmsh22(rPath, rMesh, false); }},
        {"gmsh22-binary", [](const std::string& rPath,
                             const tid::Mesh& rMesh) { tid::write_gmsh22(rPath, rMesh, true); }},
        {"gmsh41-ascii", [](const std::string& rPath,
                            const tid::Mesh& rMesh) { tid::write_gmsh41(rPath, rMesh, false); }},
        {"gmsh41-binary", [](const std::string& rPath,
                             const tid::Mesh& rMesh) { tid::write_gmsh41(rPath, rMesh, true); }},
        {"febio", tid::write_febio},
        {"z88",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_z88(rPath, rMesh); }},
        {"patran", tid::write_patran},
        {"femap", tid::write_femap},
        {"mdpa",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_mdpa(rPath, rMesh); }},
    };
    return writers;
}
}  // namespace

TEST(TextIoDtypes, HoistedWritersMatchCanonicalStorageOnAllDtypes) {
    for (const auto& [format, writer] : tid_writers()) {
        SCOPED_TRACE(format);
        for (const auto dtype : tid_dtypes) {
            SCOPED_TRACE(static_cast<int>(dtype));
            for (const std::size_t dim : {2u, 3u}) {
                SCOPED_TRACE(dim);
                const auto input = tid_mesh(format, dtype, dtype, dim);
                const auto expected = tid_mesh(format, tid::DType::Float64, tid::DType::Int64, dim);
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

TEST(TextIoDtypes, HoistedWritersLeaveCallerArraysUntouched) {
    for (const auto& [format, writer] : tid_writers()) {
        SCOPED_TRACE(format);
        for (const auto dtype : {tid::DType::Float32, tid::DType::Int32, tid::DType::UInt16}) {
            SCOPED_TRACE(static_cast<int>(dtype));
            const auto input = tid_mesh(format, dtype, dtype, 3);
            const tid::NDArray& points = input.Points();
            const tid::NDArray& conn = input.Cells(0).Conn();
            const std::vector<std::byte> points_before(points.Data(),
                                                       points.Data() + points.Nbytes());
            const std::vector<std::byte> conn_before(conn.Data(), conn.Data() + conn.Nbytes());
            const std::string path = mt::temp_path(".data");
            writer(path, input);
            std::filesystem::remove(path);
            EXPECT_EQ(std::vector<std::byte>(points.Data(), points.Data() + points.Nbytes()),
                      points_before);
            EXPECT_EQ(std::vector<std::byte>(conn.Data(), conn.Data() + conn.Nbytes()),
                      conn_before);
        }
    }
}
