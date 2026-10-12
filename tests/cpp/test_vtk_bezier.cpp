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
// VTK Bezier cells (types 75-81), the VTU markers of their rational weights and
// degrees, and the edge order of higher-order hexahedra across VTK file
// versions (VTK 9.1 swapped edges 10 and 11; XML 2.2 and legacy 5.1 say so).

// System includes
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

std::string vb_read_text(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void vb_write_text(const std::string& rPath, const std::string& rText) {
    std::ofstream out(rPath, std::ios::binary);
    out << rText;
}

std::vector<std::int64_t> vb_row(const Mesh& rMesh, std::size_t Block = 0) {
    const auto cb = rMesh.Cells(Block);
    std::vector<std::int64_t> out;
    for (std::size_t i = 0; i < cb.Conn().Size(); ++i)
        out.push_back(meshioplusplus::detail::read_int(cb.Conn(), i));
    return out;
}

// One order-2 hexahedron of `rType` over 27 distinct points.
Mesh vb_hex27(const std::string& rType) {
    std::vector<std::vector<double>> pts;
    std::vector<std::int64_t> row;
    for (int i = 0; i < 27; ++i) {
        pts.push_back(
            {static_cast<double>(i), static_cast<double>(i % 3), static_cast<double>(i % 5)});
        row.push_back(i);
    }
    return mt::make_mesh(std::move(pts), rType, {row});
}

std::vector<std::int64_t> vb_iota(std::size_t N) {
    std::vector<std::int64_t> v(N);
    for (std::size_t i = 0; i < N; ++i)
        v[i] = static_cast<std::int64_t>(i);
    return v;
}

}  // namespace

TEST(VtkBezier, TheSevenTypesRoundTripThroughVtuAndLegacy) {
    const char* types[] = {"VTK_BEZIER_CURVE",         "VTK_BEZIER_TRIANGLE",
                           "VTK_BEZIER_QUADRILATERAL", "VTK_BEZIER_TETRAHEDRON",
                           "VTK_BEZIER_HEXAHEDRON",    "VTK_BEZIER_WEDGE",
                           "VTK_BEZIER_PYRAMID"};
    const std::size_t counts[] = {3, 6, 9, 10, 27, 18, 14};
    for (std::size_t t = 0; t < 7; ++t) {
        std::vector<std::vector<double>> pts;
        for (std::size_t i = 0; i < counts[t]; ++i)
            pts.push_back({static_cast<double>(i), 0.5 * static_cast<double>(i), 1.0});
        const Mesh m = mt::make_mesh(std::move(pts), types[t], {vb_iota(counts[t])});
        const std::string vtu = mt::temp_path(".vtu");
        meshioplusplus::write_vtu(vtu, m, false, false);
        const Mesh a = meshioplusplus::read_vtu(vtu);
        EXPECT_EQ(std::string(a.Cells(0).Type()), types[t]);
        EXPECT_EQ(vb_row(a), vb_iota(counts[t])) << types[t];
        for (const bool v51 : {true, false}) {
            const std::string vtk = mt::temp_path(".vtk");
            meshioplusplus::write_vtk(vtk, m, false, v51);
            const Mesh b = meshioplusplus::read_vtk(vtk);
            EXPECT_EQ(std::string(b.Cells(0).Type()), types[t]);
            EXPECT_EQ(vb_row(b), vb_iota(counts[t])) << types[t] << (v51 ? " 5.1" : " 4.2");
        }
    }
}

TEST(VtkBezier, MarkedWeightsAndDegreesReadUnderTheirVtkNames) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < 9; ++i)
        pts.push_back({static_cast<double>(i % 3), static_cast<double>(i / 3), 0.0});
    Mesh m = mt::make_mesh(std::move(pts), "VTK_BEZIER_QUADRILATERAL", {vb_iota(9)});
    NDArray w(DType::Float64, {9});
    for (std::size_t i = 0; i < 9; ++i)
        w.As<double>()[i] = i == 8 ? 0.5 : 1.0;
    m.AddPointData("RationalWeights", std::move(w));
    NDArray deg(DType::Int64, {1, 3});
    deg.As<std::int64_t>()[0] = deg.As<std::int64_t>()[1] = 2;
    m.AddCellData("HigherOrderDegrees", std::vector<NDArray>{std::move(deg)});
    const std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, m, false, false);
    std::string text = vb_read_text(path);
    EXPECT_NE(text.find("<PointData RationalWeights=\"RationalWeights\">"), std::string::npos);
    EXPECT_NE(text.find("<CellData HigherOrderDegrees=\"HigherOrderDegrees\">"), std::string::npos);

    // A file that marks arrays of other names reads them under the markers'.
    for (const auto& [from, to] :
         {std::pair<std::string, std::string>{"RationalWeights=\"RationalWeights\"",
                                              "RationalWeights=\"w\""},
          {"HigherOrderDegrees=\"HigherOrderDegrees\"", "HigherOrderDegrees=\"deg\""},
          {"Name=\"RationalWeights\"", "Name=\"w\""},
          {"Name=\"HigherOrderDegrees\"", "Name=\"deg\""}})
        text.replace(text.find(from), from.size(), to);
    vb_write_text(path, text);
    const Mesh r = meshioplusplus::read_vtu(path);
    ASSERT_TRUE(r.HasPointData("RationalWeights"));
    EXPECT_FALSE(r.HasPointData("w"));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(r.PointData("RationalWeights"), 8), 0.5);
    ASSERT_TRUE(r.HasCellData("HigherOrderDegrees"));
    EXPECT_EQ(meshioplusplus::detail::read_int(r.CellData("HigherOrderDegrees", 0), 1), 2);
    EXPECT_EQ(meshioplusplus::read_vtu_metadata(path).mPointDataNames,
              std::vector<std::string>{"RationalWeights"});
}

TEST(VtkBezier, HigherOrderHexahedraFollowTheFileVersion) {
    std::vector<std::int64_t> old_order = vb_iota(27);
    std::swap(old_order[18], old_order[19]);
    for (const char* type : {"VTK_LAGRANGE_HEXAHEDRON", "VTK_BEZIER_HEXAHEDRON"}) {
        const Mesh m = vb_hex27(type);
        // VTU: written as 2.2 (the current order); a 0.1 file is converted.
        const std::string vtu = mt::temp_path(".vtu");
        meshioplusplus::write_vtu(vtu, m, false, false);
        std::string text = vb_read_text(vtu);
        ASSERT_NE(text.find("version=\"2.2\""), std::string::npos);
        EXPECT_EQ(vb_row(meshioplusplus::read_vtu(vtu)), vb_iota(27));
        text.replace(text.find("version=\"2.2\""), 13, "version=\"0.1\"");
        vb_write_text(vtu, text);
        EXPECT_EQ(vb_row(meshioplusplus::read_vtu(vtu)), old_order) << type;

        // Legacy: 4.2 holds the old order on disk, 5.1 the current; both read back.
        for (const bool v51 : {false, true}) {
            const std::string vtk = mt::temp_path(".vtk");
            meshioplusplus::write_vtk(vtk, m, false, v51);
            EXPECT_EQ(vb_row(meshioplusplus::read_vtk(vtk)), vb_iota(27)) << type << v51;
        }
    }
    // Other cells keep version 0.1, so their files are byte-identical.
    const std::string vtu = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(vtu, mt::tet_mesh(), false, false);
    EXPECT_NE(vb_read_text(vtu).find("version=\"0.1\""), std::string::npos);
}
