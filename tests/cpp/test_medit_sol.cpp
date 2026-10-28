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
// Medit solutions (.sol / .solb), roadmap 3.2.6, read and written by the C++ core
// directly (tests/python/test_medit_sol.py runs both engines and compares their bytes).

// System includes
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

NDArray ramp(std::vector<std::size_t> Shape, double Offset) {
    NDArray a(DType::Float64, std::move(Shape));
    for (std::size_t i = 0; i < a.Size(); ++i)
        a.As<double>()[i] = Offset + 0.5 * static_cast<double>(i);
    return a;
}

Mesh solution_mesh() {
    Mesh mesh = mt::tri_mesh();
    mesh.AddPointData("metric", ramp({4, 6}, 1.0));
    mesh.AddPointData("size", ramp({4}, 10.0));
    mesh.AddCellData("error", {ramp({2}, 100.0)});
    return mesh;
}

}  // namespace

TEST(MeditSol, RoundTripsVertexAndElementFieldsInBothEncodings) {
    for (const char* ext : {".sol", ".solb"}) {
        const Mesh mesh = solution_mesh();
        const std::string path = mt::temp_path(ext);
        const std::string format = meshioplusplus::resolve_format(path, "");
        EXPECT_EQ(format, "medit_sol");
        meshioplusplus::registry_writers().at(format)(path, mesh);

        // standalone: the vertex fields, in sorted-name order, on a geometry-less mesh
        const Mesh sol = meshioplusplus::registry_readers().at(format)(path);
        EXPECT_EQ(sol.NumPoints(), 4u);
        EXPECT_EQ(sol.PointDim(), 0u);
        EXPECT_EQ(sol.FieldData("medit:dimension").As<std::int64_t>()[0], 3);
        const NDArray& metric = sol.PointData("medit:sol:0");
        ASSERT_EQ(metric.Shape(), (std::vector<std::size_t>{4, 6}));
        EXPECT_EQ(metric.As<double>()[7], 1.0 + 0.5 * 7);
        EXPECT_EQ(sol.PointData("medit:sol:1").Shape(), (std::vector<std::size_t>{4}));

        // attached: the element field lands on its triangles
        Mesh into = mt::tri_mesh();
        meshioplusplus::read_medit_sol_into(path, into);
        EXPECT_EQ(into.CellData("medit:sol:0", 0).As<double>()[1], 100.5);
        EXPECT_EQ(into.PointData("medit:sol:1").As<double>()[3], 11.5);
    }
}

TEST(MeditSol, ElementFieldsFillOtherBlocksWithNaN) {
    // the line's value is NaN: that block is written with no section
    Mesh written = mt::tri_mesh();
    written.AddCellBlock("line", NDArray(DType::Int64, {1, 2}));
    NDArray nan(DType::Float64, {1});
    nan.As<double>()[0] = std::nan("");
    written.AddCellData("error", {ramp({2}, 1.0), std::move(nan)});
    const std::string path = mt::temp_path(".sol");
    meshioplusplus::write_medit_sol(path, written);
    Mesh into = mt::tri_mesh();
    into.AddCellBlock("line", NDArray(DType::Int64, {1, 2}));
    meshioplusplus::read_medit_sol_into(path, into);
    EXPECT_EQ(into.CellData("medit:sol:0", 0).As<double>()[1], 1.5);
    EXPECT_TRUE(std::isnan(into.CellData("medit:sol:0", 1).As<double>()[0]));
}

TEST(MeditSol, RefusesWhatItCannotRepresent) {
    // only integer data: nothing to write
    Mesh labels = mt::tri_mesh();
    NDArray ids(DType::Int64, {4});
    labels.AddPointData("ids", std::move(ids));
    EXPECT_THROW(meshioplusplus::write_medit_sol(mt::temp_path(".sol"), labels),
                 meshioplusplus::WriteError);

    // a solution of another dimension or size than the mesh
    const std::string path = mt::temp_path(".sol");
    meshioplusplus::write_medit_sol(path, solution_mesh());
    Mesh flat = mt::tri_mesh_2d();
    EXPECT_THROW(meshioplusplus::read_medit_sol_into(path, flat), meshioplusplus::ReadError);
    Mesh other = mt::quad_mesh();
    EXPECT_THROW(meshioplusplus::read_medit_sol_into(path, other), meshioplusplus::ReadError);

    // a type no libMeshb solution has, and a count beyond the file
    for (const char* text : {"Dimension 3\nSolAtVertices\n1\n1 9\n1.0\nEnd\n",
                             "Dimension 3\nSolAtVertices\n99999999999\n1 1\n1.0\nEnd\n"}) {
        const std::string bad = mt::temp_path(".sol");
        std::ofstream(bad) << text;
        EXPECT_THROW(meshioplusplus::read_medit_sol(bad), meshioplusplus::ReadError);
    }
}
