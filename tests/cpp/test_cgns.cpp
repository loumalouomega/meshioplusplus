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

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include <array>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <hdf5.h>

#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB
#include <cgnslib.h>
#endif
#include "meshioplusplus/operations/stats.hpp"

using meshioplusplus::detail::read_double;  // NOLINT
namespace h5 = meshioplusplus::h5;

namespace {

using P3 = std::array<double, 3>;

P3 cgns_mid(const P3& a, const P3& b) {
    return {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
}

P3 cgns_avg(const std::vector<P3>& rPts) {
    P3 c{0, 0, 0};
    for (const P3& p : rPts)
        for (int k = 0; k < 3; ++k)
            c[k] += p[k];
    for (int k = 0; k < 3; ++k)
        c[k] /= static_cast<double>(rPts.size());
    return c;
}

meshioplusplus::NDArray cgns_points_from(const std::vector<P3>& rPts) {
    meshioplusplus::NDArray out(meshioplusplus::DType::Float64, {rPts.size(), 3});
    double* d = out.As<double>();
    for (std::size_t i = 0; i < rPts.size(); ++i)
        for (int k = 0; k < 3; ++k)
            d[i * 3 + k] = rPts[i][k];
    return out;
}

}  // namespace

TEST(Cgns, RoundTripLinear) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_cgns(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_cgns(p); };
        mt::roundtrip(w, r, mt::line_mesh(), ".cgns");
        mt::roundtrip(w, r, mt::tri_mesh(), ".cgns");
        mt::roundtrip(w, r, mt::quad_mesh(), ".cgns");
        mt::roundtrip(w, r, mt::tet_mesh(), ".cgns");
        mt::roundtrip(w, r, mt::hex_mesh(), ".cgns");
        mt::roundtrip(w, r, mt::wedge_mesh(), ".cgns");
    }
}

TEST(Cgns, RoundTripQuadratic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_cgns(p, m, -1); };
    auto r = [](const std::string& p) { return meshioplusplus::read_cgns(p); };
    mt::roundtrip(w, r, mt::triangle6_mesh(), ".cgns");
    mt::roundtrip(w, r, mt::quad8_mesh(), ".cgns");
    mt::roundtrip(w, r, mt::tet10_mesh(), ".cgns");
    mt::roundtrip(w, r, mt::hex20_mesh(), ".cgns");  // the permutation canary
}

TEST(Cgns, RoundTripMultiBlock) {
    // Unlike MED, CGNS keeps same-type blocks as separate sections -- but
    // block count and order must still round-trip exactly.
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::Mesh m = mt::tri_quad_mesh();  // triangle(2), quad(1), triangle(1)
    meshioplusplus::write_cgns(p, m, -1);
    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);

    ASSERT_EQ(out.NumCellBlocks(), 3u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
    EXPECT_EQ(out.Cells(0).NumCells(), 2u);
    EXPECT_EQ(out.Cells(1).Type(), "quad");
    EXPECT_EQ(out.Cells(1).NumCells(), 1u);
    EXPECT_EQ(out.Cells(2).Type(), "triangle");
    EXPECT_EQ(out.Cells(2).NumCells(), 1u);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, TwoDimensionalPoints) {
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::Mesh m = mt::tri_mesh_2d();  // genuinely (n,2) points
    meshioplusplus::write_cgns(p, m, -1);

    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        h5::Hid base = h5::open_group(f, "Base");
        meshioplusplus::NDArray basedata = h5::read_dataset(base, " data");
        ASSERT_EQ(basedata.Size(), 2u);
        EXPECT_EQ(meshioplusplus::detail::read_int(basedata, 0), 2);  // CellDim
        EXPECT_EQ(meshioplusplus::detail::read_int(basedata, 1), 2);  // PhysDim
        h5::Hid coords = h5::open_group(f, "Base/Zone1/GridCoordinates");
        EXPECT_FALSE(h5::exists(coords, "CoordinateZ"));
        EXPECT_TRUE(h5::exists(coords, "CoordinateX"));
        EXPECT_TRUE(h5::exists(coords, "CoordinateY"));
    }

    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    EXPECT_EQ(out.PointDim(), 2u);
    EXPECT_EQ(out.NumPoints(), m.NumPoints());

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, NodeLayout) {
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, mt::tet_mesh(), -1);

    h5::SilenceErrors silence;
    h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);

    // Root node attributes + reserved datasets.
    EXPECT_EQ(h5::read_attr_string(f, "name"), "HDF5 MotherNode");
    EXPECT_EQ(h5::read_attr_string(f, "label"), "Root Node of HDF5 File");
    EXPECT_EQ(h5::read_attr_string(f, "type"), "MT");
    EXPECT_TRUE(h5::exists(f, " format"));
    EXPECT_TRUE(h5::exists(f, " hdf5version"));
    {
        meshioplusplus::NDArray hv = h5::read_dataset(f, " hdf5version");
        EXPECT_EQ(hv.Size(), 33u);
    }

    // CGNSLibraryVersion.
    h5::Hid cgver = h5::open_group(f, "CGNSLibraryVersion");
    EXPECT_EQ(h5::read_attr_string(cgver, "label"), "CGNSLibraryVersion_t");
    EXPECT_EQ(h5::read_attr_string(cgver, "type"), "R4");
    EXPECT_EQ(h5::read_attr_int(cgver, "flags"), 1);

    // Base.
    h5::Hid base = h5::open_group(f, "Base");
    EXPECT_EQ(h5::read_attr_string(base, "label"), "CGNSBase_t");
    meshioplusplus::NDArray basedata = h5::read_dataset(base, " data");
    ASSERT_EQ(basedata.Size(), 2u);
    EXPECT_EQ(meshioplusplus::detail::read_int(basedata, 0), 3);  // CellDim
    EXPECT_EQ(meshioplusplus::detail::read_int(basedata, 1), 3);  // PhysDim

    // Zone1: shape (3,1) = {NVertex, NCell, 0}.
    h5::Hid zone = h5::open_group(base, "Zone1");
    EXPECT_EQ(h5::read_attr_string(zone, "label"), "Zone_t");
    meshioplusplus::NDArray zdata = h5::read_dataset(zone, " data");
    ASSERT_EQ(zdata.Shape().size(), 2u);
    EXPECT_EQ(zdata.Shape()[0], 3u);
    EXPECT_EQ(zdata.Shape()[1], 1u);
    EXPECT_EQ(meshioplusplus::detail::read_int(zdata, 0), 5);  // tet_mesh: 5 points
    EXPECT_EQ(meshioplusplus::detail::read_int(zdata, 1), 2);  // 2 tetra
    EXPECT_EQ(meshioplusplus::detail::read_int(zdata, 2), 0);

    // ZoneType: 12 bytes, no NUL.
    h5::Hid zt = h5::open_group(zone, "ZoneType");
    meshioplusplus::NDArray ztdata = h5::read_dataset(zt, " data");
    EXPECT_EQ(ztdata.Size(), 12u);

    // A section's "name"/"flags" attributes and a NBR-free structural check.
    h5::Hid sect = h5::open_group(zone, "TETRA_4_1");
    EXPECT_EQ(h5::read_attr_string(sect, "name"), "TETRA_4_1");
    EXPECT_EQ(h5::read_attr_string(sect, "label"), "Elements_t");
    EXPECT_EQ(h5::read_attr_int(sect, "flags"), 1);
    meshioplusplus::NDArray sdata = h5::read_dataset(sect, " data");
    ASSERT_EQ(sdata.Size(), 2u);
    EXPECT_EQ(meshioplusplus::detail::read_int(sdata, 0), 10);  // TETRA_4 code
    h5::Hid rng = h5::open_group(sect, "ElementRange");
    meshioplusplus::NDArray range = h5::read_dataset(rng, " data");
    EXPECT_EQ(meshioplusplus::detail::read_int(range, 0), 1);
    EXPECT_EQ(meshioplusplus::detail::read_int(range, 1), 2);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, FlowSolutionRoundTrip) {
    // point_data at "Vertex", cell_data at "CellCenter" (v9.9.0). A k>1 array
    // becomes k sibling DataArray_t nodes -- CGNS has no NumberOfComponents --
    // and is re-joined on read.
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::Mesh m = mt::tri_mesh();
    const std::size_t np = m.NumPoints();
    const std::size_t nc = m.Cells(0).NumCells();

    meshioplusplus::NDArray temp(meshioplusplus::DType::Float64, {np});
    for (std::size_t i = 0; i < np; ++i)
        temp.As<double>()[i] = static_cast<double>(i) + 0.5;
    m.AddPointData("temperature", std::move(temp));

    meshioplusplus::NDArray vel(meshioplusplus::DType::Float64, {np, 3});
    for (std::size_t i = 0; i < np * 3; ++i)
        vel.As<double>()[i] = static_cast<double>(i) * 2.0;
    m.AddPointData("velocity", std::move(vel));

    meshioplusplus::NDArray mat(meshioplusplus::DType::Float64, {nc});
    for (std::size_t i = 0; i < nc; ++i)
        mat.As<double>()[i] = 7.0 + static_cast<double>(i);
    m.AddCellData("material", {std::move(mat)});

    meshioplusplus::write_cgns(p, m, -1);

    // Raw node layout: the split component nodes really are DataArray_t under
    // a FlowSolution_t whose GridLocation says Vertex/CellCenter.
    {
        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        h5::Hid sol = h5::open_group(f, "Base/Zone1/FlowSolution");
        EXPECT_EQ(h5::read_attr_string(sol, "label"), "FlowSolution_t");
        h5::Hid gl = h5::open_group(sol, "GridLocation");
        EXPECT_EQ(h5::read_attr_string(gl, "label"), "GridLocation_t");
        meshioplusplus::NDArray loc = h5::read_dataset(gl, " data");
        ASSERT_EQ(loc.Size(), 6u);  // "Vertex", no trailing NUL
        for (const char* name : {"velocity_0", "velocity_1", "velocity_2", "temperature"}) {
            ASSERT_TRUE(h5::exists(sol, name)) << name;
            h5::Hid g = h5::open_group(sol, name);
            EXPECT_EQ(h5::read_attr_string(g, "label"), "DataArray_t") << name;
            EXPECT_EQ(h5::read_dataset(g, " data").Size(), np) << name;
        }
        EXPECT_FALSE(h5::exists(sol, "velocity"));  // split, not interleaved

        h5::Hid csol = h5::open_group(f, "Base/Zone1/FlowSolutionCells");
        h5::Hid cgl = h5::open_group(csol, "GridLocation");
        EXPECT_EQ(h5::read_dataset(cgl, " data").Size(), 10u);  // "CellCenter"
        EXPECT_TRUE(h5::exists(csol, "material"));
    }

    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    ASSERT_TRUE(out.HasPointData("velocity"));
    const meshioplusplus::NDArray& v = out.PointData("velocity");
    ASSERT_EQ(v.Shape().size(), 2u);
    EXPECT_EQ(v.Shape()[0], np);
    EXPECT_EQ(v.Shape()[1], 3u);
    for (std::size_t i = 0; i < np * 3; ++i)
        EXPECT_DOUBLE_EQ(read_double(v, i), static_cast<double>(i) * 2.0);

    ASSERT_TRUE(out.HasPointData("temperature"));
    EXPECT_EQ(out.PointData("temperature").Shape().size(), 1u);  // stays 1-D
    ASSERT_TRUE(out.HasCellData("material"));
    ASSERT_EQ(out.CellDataNumBlocks("material"), 1u);
    for (std::size_t i = 0; i < nc; ++i)
        EXPECT_DOUBLE_EQ(read_double(out.CellData("material", 0), i), 7.0 + static_cast<double>(i));

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, CellDataSplitsBackAcrossBlocks) {
    // A zone-wide CellCenter array is concatenated block-major on write and
    // split back by each section's cell count on read.
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::Mesh m = mt::tri_quad_mesh();  // triangle(2), quad(1), triangle(1)
    std::vector<meshioplusplus::NDArray> blocks;
    double v = 0.0;
    for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
        meshioplusplus::NDArray d(meshioplusplus::DType::Float64, {m.Cells(b).NumCells()});
        for (std::size_t c = 0; c < m.Cells(b).NumCells(); ++c)
            d.As<double>()[c] = v++;
        blocks.push_back(std::move(d));
    }
    m.AddCellData("tag", std::move(blocks));

    meshioplusplus::write_cgns(p, m, -1);
    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);

    ASSERT_EQ(out.NumCellBlocks(), 3u);
    ASSERT_EQ(out.CellDataNumBlocks("tag"), 3u);
    double expect = 0.0;
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const meshioplusplus::NDArray& d = out.CellData("tag", b);
        ASSERT_EQ(d.Size(), out.Cells(b).NumCells()) << b;
        for (std::size_t c = 0; c < d.Size(); ++c)
            EXPECT_DOUBLE_EQ(read_double(d, c), expect++);
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, LegacyRead) {
    // The pre-v9.8.0 layout: no node attributes at all, structurally
    // distinguished from the spec layout by the absence of a "label"
    // attribute equal to "CGNSBase_t" anywhere under the root.
    std::string p = mt::temp_path(".cgns");
    {
        h5::SilenceErrors silence;
        h5::Hid f = h5::create_file(p);
        h5::Hid base = h5::create_group(f, "Base");
        h5::Hid zone = h5::create_group(base, "Zone1");
        h5::Hid coords = h5::create_group(zone, "GridCoordinates");
        const std::vector<std::vector<double>> xyz = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int c = 0; c < 3; ++c) {
            h5::Hid g = h5::create_group(coords, c == 0   ? "CoordinateX"
                                                 : c == 1 ? "CoordinateY"
                                                          : "CoordinateZ");
            meshioplusplus::NDArray col(meshioplusplus::DType::Float64, {xyz.size()});
            for (std::size_t i = 0; i < xyz.size(); ++i)
                col.As<double>()[i] = xyz[i][static_cast<std::size_t>(c)];
            h5::write_dataset(g, " data", col);
        }
        h5::Hid elems = h5::create_group(zone, "GridElements");
        h5::Hid rng = h5::create_group(elems, "ElementRange");
        h5::Hid conn = h5::create_group(elems, "ElementConnectivity");
        meshioplusplus::NDArray range(meshioplusplus::DType::Int64, {2});
        range.As<std::int64_t>()[0] = 1;
        range.As<std::int64_t>()[1] = 1;
        h5::write_dataset(rng, " data", range);
        meshioplusplus::NDArray flat(meshioplusplus::DType::Int64, {4});
        std::int64_t ids[4] = {1, 2, 3, 4};
        for (int i = 0; i < 4; ++i)
            flat.As<std::int64_t>()[i] = ids[i];
        h5::write_dataset(conn, " data", flat);
    }

    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    EXPECT_EQ(out.NumPoints(), 4u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "tetra");
    EXPECT_EQ(out.Cells(0).NumCells(), 1u);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// A jagged polygon block used to be refused outright; since v9.21.0 it is an
// NGON_n section -- which is what this step used to assert the opposite of.
TEST(Cgns, PolygonBlockRoundTripsAsNgon) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {1, 4, 2}});
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    EXPECT_EQ(out.NumPoints(), 5u);
    // Grouped by node count, so the quad and the triangle land in two blocks.
    std::size_t n3 = 0, n4 = 0;
    for (const auto cb : out.CellRange()) {
        if (cb.Type() == "polygon3")
            n3 += cb.NumCells();
        if (cb.Type() == "polygon4")
            n4 += cb.NumCells();
    }
    EXPECT_EQ(n3, 1u);
    EXPECT_EQ(n4, 1u);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, UnverifiedOrderingTypeThrows) {
    // tetra20 is a known CGNS type (TETRA_20) but its interior-node order
    // has no verified text source -- must refuse, not guess.
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from(std::vector<std::vector<double>>(20, {0, 0, 0})));
    meshioplusplus::NDArray conn(meshioplusplus::DType::Int64, {1, 20});
    for (int i = 0; i < 20; ++i)
        conn.As<std::int64_t>()[i] = i;
    m.AddCellBlock("tetra20", std::move(conn));
    std::string p = mt::temp_path(".cgns");
    EXPECT_THROW(meshioplusplus::write_cgns(p, m, -1), meshioplusplus::WriteError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, StructuredZoneThrows) {
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, mt::tet_mesh(), -1);
    {
        h5::SilenceErrors silence;
        h5::Hid f = h5::open_file_rw(p);
        h5::Hid zt = h5::open_group(f, "Base/Zone1/ZoneType");
        // Overwrite " data" with "Structured" (10 bytes, same convention).
        H5Ldelete(zt, " data", H5P_DEFAULT);
        meshioplusplus::NDArray sdata(meshioplusplus::DType::Int8, {10});
        const char* s = "Structured";
        for (int i = 0; i < 10; ++i)
            sdata.As<std::int8_t>()[i] = static_cast<std::int8_t>(s[i]);
        h5::write_dataset(zt, " data", sdata);
    }
    EXPECT_THROW(meshioplusplus::read_cgns(p), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, MixedSectionThrows) {
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, mt::tet_mesh(), -1);
    {
        h5::SilenceErrors silence;
        h5::Hid f = h5::open_file_rw(p);
        h5::Hid sect = h5::open_group(f, "Base/Zone1/TETRA_4_1");
        H5Ldelete(sect, " data", H5P_DEFAULT);
        meshioplusplus::NDArray sdata(meshioplusplus::DType::Int32, {2});
        sdata.As<std::int32_t>()[0] = 20;  // MIXED
        sdata.As<std::int32_t>()[1] = 0;
        h5::write_dataset(sect, " data", sdata);
    }
    EXPECT_THROW(meshioplusplus::read_cgns(p), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(Cgns, BadConnectivitySizeThrows) {
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, mt::tet_mesh(), -1);
    {
        h5::SilenceErrors silence;
        h5::Hid f = h5::open_file_rw(p);
        h5::Hid conn = h5::open_group(f, "Base/Zone1/TETRA_4_1/ElementConnectivity");
        H5Ldelete(conn, " data", H5P_DEFAULT);
        meshioplusplus::NDArray flat(meshioplusplus::DType::Int64, {3});  // too short
        for (int i = 0; i < 3; ++i)
            flat.As<std::int64_t>()[i] = i + 1;
        h5::write_dataset(conn, " data", flat);
    }
    EXPECT_THROW(meshioplusplus::read_cgns(p), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// ---- geometric ordering verification --------------------------------------
//
// A round trip through OUR OWN reader and writer cannot prove the CGNS
// ordering is genuinely SIDS-correct -- the self-inverse permutation makes
// read(write(m)) == m even if the permutation itself is wrong. These tests
// instead inspect the RAW file bytes against coordinates placed at their
// true geometric positions, so a wrong permutation shows up as a mismatch
// against an independently-computed expectation.

TEST(CgnsOrdering, Hexahedron27FaceCentresMatchSidsGeometry) {
    const std::vector<P3> corners = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const std::vector<std::array<int, 2>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                                   {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    // meshio/VTK hexahedron27 mid-face order (vtkTriQuadraticHexahedron):
    // 20=(0,4,7,3), 21=(1,2,6,5), 22=(0,1,5,4), 23=(3,7,6,2), 24=bottom, 25=top.
    const std::vector<std::vector<int>> faces = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
                                                 {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

    std::vector<P3> pts = corners;
    for (const auto& e : edges)
        pts.push_back(cgns_mid(corners[e[0]], corners[e[1]]));
    for (const auto& fc : faces) {
        std::vector<P3> fp;
        for (int i : fc)
            fp.push_back(corners[i]);
        pts.push_back(cgns_avg(fp));
    }
    pts.push_back(cgns_avg(corners));  // body center, meshio node 26

    meshioplusplus::Mesh m;
    m.AssignPoints(cgns_points_from(pts));
    meshioplusplus::NDArray conn(meshioplusplus::DType::Int64, {1, 27});
    for (int i = 0; i < 27; ++i)
        conn.As<std::int64_t>()[i] = i;
    m.AddCellBlock("hexahedron27", std::move(conn));

    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    h5::SilenceErrors silence;
    h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    h5::Hid conn_grp = h5::open_group(f, "Base/Zone1/HEXA_27_1/ElementConnectivity");
    meshioplusplus::NDArray raw = h5::read_dataset(conn_grp, " data");
    ASSERT_EQ(raw.Size(), 27u);

    // SIDS face order F1..F6 = (0,3,2,1)(0,1,5,4)(1,2,6,5)(2,3,7,6)(0,4,7,3)
    // (4,5,6,7) => CGNS mid-face column c (0-based, c=20..25) must hold
    // meshio (VTK) node `expect_meshio_index[c-20] + 1` (1-based).
    const std::vector<int> perm = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 16, 17,
                                   18, 19, 12, 13, 14, 15, 24, 22, 21, 23, 20, 25, 26};
    for (int c = 0; c < 27; ++c)
        EXPECT_EQ(meshioplusplus::detail::read_int(raw, static_cast<std::size_t>(c)),
                  perm[static_cast<std::size_t>(c)] + 1)
            << "CGNS column " << c;

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

TEST(CgnsOrdering, Hexahedron20EdgeMidpointsMatchSidsGeometry) {
    // SIDS HEXA_20 mid-edge order is bottom ring, VERTICALS, top ring;
    // meshio's own order (detail::cell_refine_edges) is bottom ring, top
    // ring, verticals -- the two differ, which is exactly what the
    // hexahedron20 permutation table exists to bridge. Build the
    // connectivity in meshio's own order, then check the raw file's mid-edge
    // columns against the true SIDS edge midpoints geometrically (not by
    // array index), independent of cgns.cpp's own permutation table.
    const std::vector<P3> corners = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const std::vector<std::array<int, 2>> sids_edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom ring -> CGNS columns 8-11
        {0, 4}, {1, 5}, {2, 6}, {3, 7},  // verticals   -> CGNS columns 12-15
        {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top ring    -> CGNS columns 16-19
    };

    std::vector<P3> meshio_pts = corners;
    for (const auto& e :
         meshioplusplus::detail::cell_refine_edges(meshioplusplus::CellType::Hexahedron))
        meshio_pts.push_back(cgns_mid(corners[e[0]], corners[e[1]]));
    ASSERT_EQ(meshio_pts.size(), 20u);

    meshioplusplus::Mesh m;
    m.AssignPoints(cgns_points_from(meshio_pts));
    meshioplusplus::NDArray conn(meshioplusplus::DType::Int64, {1, 20});
    for (int i = 0; i < 20; ++i)
        conn.As<std::int64_t>()[i] = i;
    m.AddCellBlock("hexahedron20", std::move(conn));

    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    h5::SilenceErrors silence;
    h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    h5::Hid conn_grp = h5::open_group(f, "Base/Zone1/HEXA_20_1/ElementConnectivity");
    meshioplusplus::NDArray raw = h5::read_dataset(conn_grp, " data");
    ASSERT_EQ(raw.Size(), 20u);

    for (int c = 0; c < 8; ++c)  // corners are always identity
        EXPECT_EQ(meshioplusplus::detail::read_int(raw, static_cast<std::size_t>(c)), c + 1);

    for (std::size_t k = 0; k < sids_edges.size(); ++k) {
        const P3 expect = cgns_mid(corners[sids_edges[k][0]], corners[sids_edges[k][1]]);
        const std::int64_t node1based = meshioplusplus::detail::read_int(raw, 8 + k);
        ASSERT_GE(node1based, 1);
        const P3& got = meshio_pts[static_cast<std::size_t>(node1based - 1)];
        EXPECT_NEAR(got[0], expect[0], 1e-12) << "CGNS column " << (8 + k);
        EXPECT_NEAR(got[1], expect[1], 1e-12) << "CGNS column " << (8 + k);
        EXPECT_NEAR(got[2], expect[2], 1e-12) << "CGNS column " << (8 + k);
    }

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// ---- pure computation: table sanity ----------------------------------------

TEST(CgnsOrdering, PermutationsAreInvolutions) {
    // Every table entry in cgns.cpp's cell type table must be self-inverse,
    // since read_cgns applies the exact same array to undo write_cgns's
    // permutation (shift -1 instead of +1). Extract each permutation as
    // observed in a freshly WRITTEN file (write side only -- this does not
    // exercise read_cgns at all) and check p[p[c]] == c mathematically,
    // which is a necessary (though not sufficient) correctness precondition
    // for reusing one array both ways.
    for (const auto& [type, npc] : std::vector<std::pair<std::string, int>>{
             {"wedge15", 15}, {"wedge18", 18}, {"hexahedron20", 20}, {"hexahedron27", 27}}) {
        meshioplusplus::Mesh m;
        std::vector<std::vector<double>> pts(static_cast<std::size_t>(npc),
                                             std::vector<double>{0, 0, 0});
        m.AssignPoints(mt::points_from(pts));
        meshioplusplus::NDArray conn(meshioplusplus::DType::Int64,
                                     {1, static_cast<std::size_t>(npc)});
        for (int i = 0; i < npc; ++i)
            conn.As<std::int64_t>()[i] = i;  // identity labels: meshio node i has value i
        m.AddCellBlock(type, std::move(conn));

        std::string p = mt::temp_path(".cgns");
        meshioplusplus::write_cgns(p, m, -1);

        h5::SilenceErrors silence;
        h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
        h5::Hid zone = h5::open_group(f, "Base/Zone1");
        std::string section_name;
        for (const auto& child : h5::group_links(zone)) {
            if (!child.empty() && child[0] == ' ')
                continue;  // Zone1's own " data" payload
            h5::Hid g = h5::open_group(zone, child);
            if (h5::has_attr(g, "label") && h5::read_attr_string(g, "label") == "Elements_t") {
                section_name = child;
                break;
            }
        }
        ASSERT_FALSE(section_name.empty()) << type;
        h5::Hid conn_grp = h5::open_group(f, "Base/Zone1/" + section_name + "/ElementConnectivity");
        meshioplusplus::NDArray raw = h5::read_dataset(conn_grp, " data");
        ASSERT_EQ(static_cast<int>(raw.Size()), npc) << type;

        // perm[c] = the meshio node placed at CGNS column c (0-based; file
        // values are 1-based).
        std::vector<int> perm(static_cast<std::size_t>(npc));
        for (int c = 0; c < npc; ++c)
            perm[static_cast<std::size_t>(c)] = static_cast<int>(
                meshioplusplus::detail::read_int(raw, static_cast<std::size_t>(c)) - 1);
        for (int c = 0; c < npc; ++c)
            EXPECT_EQ(perm[static_cast<std::size_t>(perm[static_cast<std::size_t>(c)])], c)
                << type << " index " << c << " is not an involution";

        std::error_code ec;
        std::filesystem::remove(p, ec);
    }
}

// --- the cgnslib (MLL) backend ----------------------------------------------
//
// Two capabilities justify the optional dependency, and both are things the
// raw-HDF5 reader in cgns.cpp fundamentally cannot have. Everything else about
// the backend is covered for free by the suite above, which runs through it
// whenever the flag is on.

TEST(CgnsMll, ReportsWhetherTheBackendIsBuilt) {
    // The predicate must agree with the macro, or every skip below is a lie.
#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB
    EXPECT_TRUE(meshioplusplus::cgns_has_cgnslib());
#else
    EXPECT_FALSE(meshioplusplus::cgns_has_cgnslib());
    // Compiled out, the entry point still exists and throws NAMING the flag --
    // the partition_kahip_parts contract. A link error would break the
    // Python-fallback shim; a silent downgrade would answer a different
    // question than the caller asked.
    EXPECT_THROW(meshioplusplus::read_cgns_mll("nonexistent.cgns"), meshioplusplus::ReadError);
    try {
        meshioplusplus::read_cgns_mll("nonexistent.cgns");
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_CGNSLIB"), std::string::npos);
    }
#endif
}

#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB

TEST(CgnsMll, ReadsAnAdfContainerTheRawHdf5ReaderCannotOpen) {
    // `.cgns` has two on-disk containers. cgns.cpp speaks HDF5 directly, so an
    // ADF file is not merely unimplemented there -- it is unreachable by
    // construction. This is the strongest single argument for the backend, and
    // it is what makes the CGNS project's own example meshes usable.
    if (!std::system("command -v cgnsconvert > /dev/null 2>&1") == 0)
        GTEST_SKIP() << "cgnsconvert not on PATH";

    const mt::Mesh m =
        mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}});
    const std::string hdf5 = mt::temp_path("_mll_hdf5.cgns");
    const std::string adf = mt::temp_path("_mll_adf.cgns");
    meshioplusplus::write_cgns(hdf5, m, -1);

    // `-a` rewrites into the ADF container through cgnslib's own API.
    const std::string cmd = "cgnsconvert -a '" + hdf5 + "' '" + adf + "' > /dev/null 2>&1";
    if (std::system(cmd.c_str()) != 0) {
        std::remove(hdf5.c_str());
        GTEST_SKIP() << "cgnsconvert -a failed (this cgnslib may be built without ADF)";
    }

    // The MLL reads it...
    const mt::Mesh via_mll = meshioplusplus::read_cgns_mll(adf);
    EXPECT_EQ(via_mll.NumPoints(), 4u);
    ASSERT_EQ(via_mll.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(via_mll.Cells(0).Type()), "tetra");

    std::remove(hdf5.c_str());
    std::remove(adf.c_str());
}

TEST(CgnsMll, ReadsNgonNfacePolyhedralSections) {
    // The other capability: NGON_n lists faces, NFACE_n lists each cell as
    // SIGNED face ids (negative = traverse that face reversed, which is how
    // CGNS orients a shared face outward from each of the two cells using it).
    // Written here through cgnslib's own API, so the bytes are the MLL's and
    // the test is not merely reading back our own idea of the encoding.
    const std::string path = mt::temp_path("_mll_ngon.cgns");
    std::remove(path.c_str());

    int fn = 0, B = 0, Z = 0, S = 0;
    ASSERT_EQ(cg_open(path.c_str(), CG_MODE_WRITE, &fn), CG_OK) << cg_get_error();
    ASSERT_EQ(cg_base_write(fn, "Base", 3, 3, &B), CG_OK);

    // A unit cube: 8 points, 6 quad faces, 1 cell.
    cgsize_t zsize[3] = {8, 1, 0};
    ASSERT_EQ(cg_zone_write(fn, B, "Zone", zsize, CGNS_ENUMV(Unstructured), &Z), CG_OK);
    const double x[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    const double y[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    const double z[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    int c = 0;
    ASSERT_EQ(cg_coord_write(fn, B, Z, CGNS_ENUMV(RealDouble), "CoordinateX", x, &c), CG_OK);
    ASSERT_EQ(cg_coord_write(fn, B, Z, CGNS_ENUMV(RealDouble), "CoordinateY", y, &c), CG_OK);
    ASSERT_EQ(cg_coord_write(fn, B, Z, CGNS_ENUMV(RealDouble), "CoordinateZ", z, &c), CG_OK);

    // NGON_n: the cube's six quads, 1-based, outward.
    const cgsize_t faces[24] = {1, 4, 3, 2, 5, 6, 7, 8, 1, 2, 6, 5,
                                3, 4, 8, 7, 1, 5, 8, 4, 2, 3, 7, 6};
    cgsize_t face_off[7] = {0, 4, 8, 12, 16, 20, 24};
    ASSERT_EQ(
        cg_poly_section_write(fn, B, Z, "Faces", CGNS_ENUMV(NGON_n), 1, 6, 0, faces, face_off, &S),
        CG_OK)
        << cg_get_error();

    // NFACE_n: one cell referencing all six faces (ids 1..6).
    const cgsize_t cells[6] = {1, 2, 3, 4, 5, 6};
    cgsize_t cell_off[2] = {0, 6};
    ASSERT_EQ(
        cg_poly_section_write(fn, B, Z, "Cells", CGNS_ENUMV(NFACE_n), 7, 7, 0, cells, cell_off, &S),
        CG_OK)
        << cg_get_error();
    ASSERT_EQ(cg_close(fn), CG_OK);

    const mt::Mesh got = meshioplusplus::read_cgns_mll(path);
    EXPECT_EQ(got.NumPoints(), 8u);
    ASSERT_EQ(got.NumCellBlocks(), 1u);
    const auto cb = got.Cells(0);
    EXPECT_TRUE(cb.IsPolyhedron());
    EXPECT_EQ(std::string(cb.Type()), "polyhedron8");
    ASSERT_EQ(cb.NumCells(), 1u);
    EXPECT_EQ(cb.NumFaces(0), 6u);

    // The geometry must actually be the unit cube -- reading six faces of the
    // right arity proves nothing about whether the node ids landed correctly.
    EXPECT_NEAR(meshioplusplus::compute_stats(got).mUnsignedVolume, 1.0, 1e-12);

    std::remove(path.c_str());
}

#endif  // MESHIOPLUSPLUS_HAS_CGNSLIB


// --------------------------------------------------------------------------
// NGON_n / NFACE_n -- polyhedral cells, hand-rolled in both directions
// --------------------------------------------------------------------------

namespace {

// Unit cube as a polyhedron block, faces wound outward.
meshioplusplus::Mesh cgns_poly_cube() {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1}}));
    m.AddPolyhedronBlock("polyhedron8", {{{
                                             {0, 3, 2, 1},
                                             {4, 5, 6, 7},
                                             {0, 1, 5, 4},
                                             {2, 3, 7, 6},
                                             {1, 2, 6, 5},
                                             {0, 4, 7, 3},
                                         }}});
    return m;
}

// Total volume via the geometric kernel -- the oracle that a face list came
// back with the right node ids AND the right winding. Cell counts cannot see
// either.
double cgns_total_volume(const meshioplusplus::Mesh& rM) {
    double v = 0.0;
    meshioplusplus::detail::CellRings rings;
    std::vector<meshioplusplus::detail::Vec3> coords;
    for (const auto cb : rM.CellRange()) {
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            if (!meshioplusplus::detail::cell_rings(cb, i, rM.Points(), rM.PointDim(), rings,
                                                    coords))
                continue;
            meshioplusplus::detail::orient_rings(rings, coords.data());
            v += meshioplusplus::detail::poly_measure(rings, coords.data()).mVolume;
        }
    }
    return v;
}

}  // namespace

TEST(CgnsPoly, APolyhedronRoundTripsThroughNgonAndNface) {
    const meshioplusplus::Mesh m = cgns_poly_cube();
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    const meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    EXPECT_EQ(out.NumPoints(), 8u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "polyhedron8");
    EXPECT_EQ(out.Cells(0).NumCells(), 1u);
    EXPECT_EQ(out.Cells(0).NumFaces(0), 6u);
    // The oracle: right ids, right winding.
    EXPECT_NEAR(cgns_total_volume(out), 1.0, 1e-12);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// The face pool must NOT come back as cells of its own -- that would double
// every polyhedron's geometry while leaving the cell count looking plausible.
TEST(CgnsPoly, TheSharedNgonFacePoolIsNotEmittedAsCells) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    // Two cubes sharing the plane x = 1, so one face is genuinely internal.
    m.AddPolyhedronBlock("polyhedron8", {{
                                             {
                                                 {0, 3, 2, 1},
                                                 {4, 5, 6, 7},
                                                 {0, 1, 5, 4},
                                                 {2, 3, 7, 6},
                                                 {1, 2, 6, 5},
                                                 {0, 4, 7, 3},
                                             },
                                             {
                                                 {1, 2, 9, 8},
                                                 {5, 11, 10, 6},
                                                 {1, 8, 11, 5},
                                                 {2, 6, 10, 9},
                                                 {8, 9, 10, 11},
                                                 {1, 5, 6, 2},
                                             },
                                         }});
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    const meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    std::size_t ncells = 0, npoly2d = 0;
    for (const auto cb : out.CellRange()) {
        if (std::string(cb.Type()).rfind("polyhedron", 0) == 0)
            ncells += cb.NumCells();
        if (std::string(cb.Type()).rfind("polygon", 0) == 0)
            npoly2d += cb.NumCells();
    }
    EXPECT_EQ(ncells, 2u);
    EXPECT_EQ(npoly2d, 0u) << "the NGON_n face pool was emitted as polygon cells too";
    EXPECT_NEAR(cgns_total_volume(out), 2.0, 1e-12);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// Mixing is the case that makes the block-filtered face table necessary: a
// hexahedron keeps its HEXA_8 section, and its faces must NOT land in the
// NGON_n pool as elements nothing references.
TEST(CgnsPoly, AHexahedronAndAPolyhedronCoexistInOneFile) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {2, 0, 0},
                                    {2, 1, 0},
                                    {2, 1, 1},
                                    {2, 0, 1}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    m.AddPolyhedronBlock("polyhedron8", {{{
                                             {1, 2, 9, 8},
                                             {5, 11, 10, 6},
                                             {1, 8, 11, 5},
                                             {2, 6, 10, 9},
                                             {8, 9, 10, 11},
                                             {1, 5, 6, 2},
                                         }}});
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    const meshioplusplus::Mesh out = meshioplusplus::read_cgns(p);
    bool has_hex = false, has_poly = false;
    for (const auto cb : out.CellRange()) {
        has_hex = has_hex || cb.Type() == "hexahedron";
        has_poly = has_poly || std::string(cb.Type()).rfind("polyhedron", 0) == 0;
    }
    EXPECT_TRUE(has_hex) << "the hexahedron was not kept as a HEXA_8 section";
    EXPECT_TRUE(has_poly);
    EXPECT_NEAR(cgns_total_volume(out), 2.0, 1e-12);

    std::error_code ec;
    std::filesystem::remove(p, ec);
}

// The CGNS 3.x inline-length layout must be refused BY NAME, pointing at the
// cgnslib backend, rather than misread as the 4.0 one.
TEST(CgnsPoly, A3xLayoutSectionIsRefusedByName) {
    const meshioplusplus::Mesh m = cgns_poly_cube();
    std::string p = mt::temp_path(".cgns");
    meshioplusplus::write_cgns(p, m, -1);

    // Delete ElementStartOffset from the NGON section -> a 3.x-shaped file.
    {
        h5::Hid f = h5::open_file_rw(p);
        H5Ldelete(f, "/Base/Zone1/NGON_n_1/ElementStartOffset", H5P_DEFAULT);
    }
    try {
        meshioplusplus::read_cgns(p);
        ADD_FAILURE() << "a 3.x-layout section was read rather than refused";
    } catch (const meshioplusplus::ReadError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("CGNSLIB"), std::string::npos) << msg;
    }
    std::error_code ec;
    std::filesystem::remove(p, ec);
}


// The declared CGNSLibraryVersion is load-bearing for a face-based file, not
// cosmetic: BELOW 4.0 cgnslib reads NGON_n/NFACE_n with the 3.x inline-length
// layout, so a file declaring 3.1 while writing ElementStartOffset arrays makes
// cgnslib splice the offsets into the connectivity and then corrupt its own
// heap. cgnscheck aborted rather than diagnosing it, so nothing but this test
// stands between a regression here and silently unreadable files.
TEST(CgnsPoly, AFaceBasedFileDeclaresCgns40) {
    auto version_of = [](const meshioplusplus::Mesh& m) {
        std::string p = mt::temp_path(".cgns");
        meshioplusplus::write_cgns(p, m, -1);
        float v = 0.0f;
        {
            h5::SilenceErrors silence;
            h5::Hid f(H5Fopen(p.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
            h5::Hid g = h5::open_group(f, "CGNSLibraryVersion");
            meshioplusplus::NDArray d = h5::read_dataset(g, " data");
            v = static_cast<float>(read_double(d, 0));
        }
        std::error_code ec;
        std::filesystem::remove(p, ec);
        return v;
    };
    EXPECT_GE(version_of(cgns_poly_cube()), 4.0f);
    // A file with no face-based section reads identically under either number,
    // so its bytes are deliberately left alone.
    EXPECT_LT(version_of(mt::tet_mesh()), 4.0f);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
