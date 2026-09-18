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

// External includes
#include <gtest/gtest.h>

// System includes
#include <filesystem>
#include <fstream>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/wkt.hpp"

TEST(Su2, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_su2(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_su2(p); };
    mt::roundtrip(w, r, mt::tri_mesh_2d(), ".su2");
    mt::roundtrip(w, r, mt::tet_mesh(), ".su2");
    mt::roundtrip(w, r, mt::hex_mesh(), ".su2");
}

TEST(Flac3d, AsciiAndBinary) {
    for (bool binary : {false, true}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_flac3d(p, m, ".16e", binary);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_flac3d(p); };
        mt::roundtrip(w, r, mt::tet_mesh(), ".f3grid");
        mt::roundtrip(w, r, mt::hex_mesh(), ".f3grid");
    }
}

TEST(Ansys, AsciiAndBinary) {
    for (bool binary : {false, true}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_ansys(p, m, binary);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_ansys(p); };
        mt::roundtrip(w, r, mt::tri_mesh_2d(), ".msh");
        mt::roundtrip(w, r, mt::tet_mesh(), ".msh");
        mt::roundtrip(w, r, mt::hex_mesh(), ".msh");
        mt::roundtrip(w, r, mt::tri_quad_mesh(), ".msh");
    }
}

TEST(Dolfin, TriangleTetra) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_dolfin(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_dolfin(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".xml");
    mt::roundtrip(w, r, mt::tri_mesh_2d(), ".xml");
    mt::roundtrip(w, r, mt::tet_mesh(), ".xml");
}

TEST(Dolfin, PointDataRoundTripsAsADimZeroMeshFunction) {
    // The regression: point_data was dropped outright while cell_data
    // round-tripped through the same sibling-file mechanism. `dim` is the
    // topological dimension of the entities a mesh function is defined on, so
    // vertices are 0 -- which is the whole discriminator on read.
    mt::Mesh in = mt::tri_mesh();
    const std::size_t npts = in.NumPoints();
    meshioplusplus::NDArray pd(meshioplusplus::DType::Float64, {npts});
    for (std::size_t i = 0; i < npts; ++i)
        reinterpret_cast<double*>(pd.Data())[i] = 1.5 + static_cast<double>(i);
    in.AddPointData("temp", std::move(pd));

    const std::string path = mt::temp_path(".xml");
    meshioplusplus::write_dolfin(path, in);

    // The sibling file must actually say dim="0" -- reading our own output back
    // would pass even with a wrong dim, since the reader would then just put it
    // in cell_data and the values would still be there.
    const std::string sibling = std::filesystem::path(path).parent_path().string() + "/" +
                                std::filesystem::path(path).stem().string() + "_temp.xml";
    {
        std::ifstream f(sibling);
        ASSERT_TRUE(f.good()) << "no sibling mesh_function file at " << sibling;
        std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        EXPECT_NE(text.find("dim=\"0\""), std::string::npos) << text;
    }

    mt::Mesh out = meshioplusplus::read_dolfin(path);
    ASSERT_TRUE(out.HasPointData("temp"));
    EXPECT_FALSE(out.HasCellData("temp"));
    const meshioplusplus::NDArray& back = out.PointData("temp");
    ASSERT_EQ(back.Size(), npts);
    for (std::size_t i = 0; i < npts; ++i)
        EXPECT_DOUBLE_EQ(reinterpret_cast<const double*>(back.Data())[i],
                         1.5 + static_cast<double>(i));

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(sibling, ec);
}

TEST(Dolfin, TwoCellDataBlocksBothSurviveInTheSiblingFile) {
    // Regression: write_mesh_function used to reopen the same sibling file
    // "<stem>_<name>.xml" once per contributing cell block (truncating each
    // time), so on a multi-block mesh only the LAST block's values survived.
    // The mesh file itself concatenates every block of the chosen type (here
    // both triangle blocks), so the sibling file must match row for row.
    mt::Mesh in;
    in.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
    in.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));
    in.AddCellBlock("triangle", mt::conn_from({{1, 4, 5}}));

    meshioplusplus::NDArray block0(meshioplusplus::DType::Float64, {2});
    block0.As<double>()[0] = 10.0;
    block0.As<double>()[1] = 11.0;
    meshioplusplus::NDArray block1(meshioplusplus::DType::Float64, {1});
    block1.As<double>()[0] = 20.0;
    std::vector<meshioplusplus::NDArray> blocks;
    blocks.push_back(std::move(block0));
    blocks.push_back(std::move(block1));
    in.AddCellData("region", std::move(blocks));

    const std::string path = mt::temp_path(".xml");
    meshioplusplus::write_dolfin(path, in);

    const std::string sibling = std::filesystem::path(path).parent_path().string() + "/" +
                                std::filesystem::path(path).stem().string() + "_region.xml";
    mt::Mesh out = meshioplusplus::read_dolfin(path);
    ASSERT_TRUE(out.HasCellData("region"));
    ASSERT_EQ(out.CellDataNumBlocks("region"), 1u)
        << "the reader always merges into a single cell block for the file's one cell type";
    const meshioplusplus::NDArray& back = out.CellData("region", 0);
    ASSERT_EQ(back.Size(), 3u) << "the second block's row was overwritten/lost";
    EXPECT_DOUBLE_EQ(reinterpret_cast<const double*>(back.Data())[0], 10.0);
    EXPECT_DOUBLE_EQ(reinterpret_cast<const double*>(back.Data())[1], 11.0);
    EXPECT_DOUBLE_EQ(reinterpret_cast<const double*>(back.Data())[2], 20.0);

    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(sibling, ec);
}

TEST(Dolfin, WarnsAboutTheLegacyFormatAndDiscardedCellTypes) {
    // Both warnings are documented (formats/dolfin.hpp, doc/formats/dolfin.md)
    // and the Python engine already emits them; the C++ writer used to emit
    // neither.
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));

    const std::string path = mt::temp_path(".xml");
    testing::internal::CaptureStderr();
    meshioplusplus::write_dolfin(path, m);
    const std::string err = testing::internal::GetCapturedStderr();
    EXPECT_NE(err.find("legacy format"), std::string::npos) << err;
    EXPECT_NE(err.find("quad"), std::string::npos) << err;

    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Wkt, TriangleGeometry) {
    // WKT (TIN) de-duplicates points, so point order is not preserved; check
    // that the triangle count round-trips.
    mt::Mesh in = mt::tri_mesh();
    std::string path = mt::temp_path(".wkt");
    meshioplusplus::write_wkt(path, in);
    mt::Mesh out = meshioplusplus::read_wkt(path);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    EXPECT_EQ(out.Cells(0).NumCells(), 2u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Tetgen, TetraPair) {
    // TetGen writes a .node/.ele pair sharing a stem.
    mt::Mesh in = mt::tet_mesh();
    std::string node = mt::temp_path(".node");
    meshioplusplus::write_tetgen(node, in);
    mt::Mesh out = meshioplusplus::read_tetgen(node);
    mt::expect_mesh_eq(in, out);
    std::string ele = node.substr(0, node.size() - 5) + ".ele";
    std::error_code ec;
    std::filesystem::remove(node, ec);
    std::filesystem::remove(ele, ec);
}

// Malformed-input paths: the readers must raise ReadError rather than silently
// mis-parse. These exercise the error branches that self round-trips never hit.

TEST(Tetgen, ReadRejectsMalformedNodeHeader) {
    std::string node = mt::temp_path(".node");
    {
        std::ofstream f(node);
        f << "not a valid header\n";
    }
    EXPECT_THROW(meshioplusplus::read_tetgen(node), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(node, ec);
}

TEST(Su2, ReadRejectsInvalidNdime) {
    std::string path = mt::temp_path(".su2");
    {
        std::ofstream f(path);
        f << "NDIME= 9\n";
    }
    EXPECT_THROW(meshioplusplus::read_su2(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
