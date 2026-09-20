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

#ifdef MESHIOPLUSPLUS_HAS_HDF5

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <numeric>

#include <hdf5.h>

#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"

namespace h5 = meshioplusplus::h5;
using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

namespace {

using I64Vec = std::vector<std::int64_t>;

struct VtkhdfTempFile {
    std::string mPath;
    explicit VtkhdfTempFile(const std::string& rSuffix = ".vtkhdf")
        : mPath(mt::temp_path(rSuffix)) {}
    ~VtkhdfTempFile() {
        std::error_code ec;
        std::filesystem::remove(mPath, ec);
    }
};

mt::NDArray vtkhdf_test_i64(const I64Vec& rV, std::vector<std::size_t> shape = {}) {
    if (shape.empty())
        shape = {rV.size()};
    mt::NDArray a = mt::NDArray::Uninit(meshioplusplus::DType::Int64, shape);
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

std::vector<double> vtkhdf_test_doubles(const mt::NDArray& rA) {
    std::vector<double> out(rA.Size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = read_double(rA, i);
    return out;
}

I64Vec vtkhdf_test_ints(const mt::NDArray& rA) {
    I64Vec out(rA.Size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = read_int(rA, i);
    return out;
}

// One tetra-pair piece, shifted in x, with a point value and two cell values.
mt::Mesh vtkhdf_test_tets(double Shift = 0.0, double U = 0.0) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0 + Shift, 0, 0},
                                    {1 + Shift, 0, 0},
                                    {0 + Shift, 1, 0},
                                    {0 + Shift, 0, 1},
                                    {1 + Shift, 1, 1}}));
    m.AddCellBlock("tetra", vtkhdf_test_i64({0, 1, 2, 3, 1, 2, 3, 4}, {2, 4}));
    m.AddPointData("u", mt::data_array({U, U, U, U, U}));
    m.AddCellData("p", {mt::data_array({U, U + 0.5})});
    return m;
}

// tetra, cube-shaped polyhedron, tetra: the OpenFOAM shape.
mt::Mesh vtkhdf_test_cube_poly(double Shift = 0.0, double U = 0.0) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0 + Shift, 0, 0},
                                    {1 + Shift, 0, 0},
                                    {1 + Shift, 1, 0},
                                    {0 + Shift, 1, 0},
                                    {0 + Shift, 0, 1},
                                    {1 + Shift, 0, 1},
                                    {1 + Shift, 1, 1},
                                    {0 + Shift, 1, 1},
                                    {0 + Shift, 0, 2}}));
    m.AddCellBlock("tetra", vtkhdf_test_i64({4, 5, 6, 8}, {1, 4}));
    m.AddPolyhedronBlock(
        "polyhedron8",
        {{{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}}});
    m.AddCellBlock("tetra", vtkhdf_test_i64({0, 1, 3, 4}, {1, 4}));
    m.AddPointData("u", mt::data_array(std::vector<double>(9, U)));
    m.AddCellData("p", {mt::data_array({U}), mt::data_array({U + 1}), mt::data_array({U + 2})});
    return m;
}

// polyhedron6, polyhedron4, polyhedron6: one file run whose node counts alternate.
mt::Mesh vtkhdf_test_prisms(double Shift, double U) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0 + Shift, 0, 0},
                                    {1 + Shift, 0, 0},
                                    {0 + Shift, 1, 0},
                                    {0 + Shift, 0, 1},
                                    {1 + Shift, 0, 1},
                                    {0 + Shift, 1, 1},
                                    {3 + Shift, 0, 0},
                                    {4 + Shift, 0, 0},
                                    {3 + Shift, 1, 0},
                                    {3 + Shift, 0, 1}}));
    const std::vector<I64Vec> prism{{0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    const std::vector<I64Vec> tet{{6, 8, 7}, {6, 7, 9}, {7, 8, 9}, {8, 6, 9}};
    m.AddPolyhedronBlock("polyhedron6", {prism});
    m.AddPolyhedronBlock("polyhedron4", {tet});
    m.AddPolyhedronBlock("polyhedron6", {prism});
    m.AddCellData("p", {mt::data_array({U}), mt::data_array({U + 1}), mt::data_array({U + 2})});
    return m;
}

// ---- raw HDF5 building blocks for files this writer does not produce ------------- //
void vtkhdf_raw_type(hid_t Grp, const char* pType, std::int64_t Major = 2, std::int64_t Minor = 0) {
    h5::write_attr_int_array(Grp, "Version", {Major, Minor});
    h5::write_attr_string_fixed(Grp, "Type", pType);
}

// A partitioned UnstructuredGrid (piece-local ids, n+1 offsets per piece) from single-piece
// meshes, laid out the way vtkHDFWriter does it. Each piece is a `vtkhdf_test_tets(Shift, U)`,
// given as the pair and built here one at a time: a KRATOS mesh cannot be copied into a list.
void vtkhdf_write_partitioned(const std::string& rPath,
                              const std::vector<std::pair<double, double>>& rPieces) {
    h5::Hid f = h5::create_file(rPath);
    h5::Hid g = h5::create_group_crt(f, "VTKHDF");
    vtkhdf_raw_type(g, "UnstructuredGrid");
    I64Vec npts, ncells, nconn, conn, offs, types;
    std::vector<double> pts, u, p;
    for (const auto& [shift, value] : rPieces) {
        const mt::Mesh m = vtkhdf_test_tets(shift, value);
        npts.push_back(static_cast<std::int64_t>(m.NumPoints()));
        const auto cb = m.Cells(0);
        ncells.push_back(static_cast<std::int64_t>(cb.NumCells()));
        offs.push_back(0);
        const std::size_t k = meshioplusplus::detail::cols(cb.Conn());
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            for (std::size_t j = 0; j < k; ++j)
                conn.push_back(read_int(cb.Conn(), r * k + j));
            offs.push_back(static_cast<std::int64_t>(conn.size()) -
                           (std::accumulate(nconn.begin(), nconn.end(), std::int64_t{0})));
            types.push_back(10);
        }
        nconn.push_back(static_cast<std::int64_t>(cb.NumCells() * k));
        for (std::size_t i = 0; i < m.NumPoints() * 3; ++i)
            pts.push_back(read_double(m.Points(), i));
        for (double v : vtkhdf_test_doubles(m.PointData("u")))
            u.push_back(v);
        for (double v : vtkhdf_test_doubles(m.CellData("p", 0)))
            p.push_back(v);
    }
    h5::write_dataset(g, "NumberOfPoints", vtkhdf_test_i64(npts));
    h5::write_dataset(g, "NumberOfCells", vtkhdf_test_i64(ncells));
    h5::write_dataset(g, "NumberOfConnectivityIds", vtkhdf_test_i64(nconn));
    h5::write_dataset(g, "Points", mt::data_array(pts, 3));
    h5::write_dataset(g, "Connectivity", vtkhdf_test_i64(conn));
    h5::write_dataset(g, "Offsets", vtkhdf_test_i64(offs));
    mt::NDArray t = mt::NDArray::Uninit(meshioplusplus::DType::UInt8, {types.size()});
    std::fill(t.As<std::uint8_t>(), t.As<std::uint8_t>() + types.size(), std::uint8_t{10});
    h5::write_dataset(g, "Types", t);
    h5::Hid pd = h5::create_group(g, "PointData");
    h5::write_dataset(pd, "u", mt::data_array(u));
    h5::Hid cd = h5::create_group(g, "CellData");
    h5::write_dataset(cd, "p", mt::data_array(p));
    h5::Hid fd = h5::create_group(g, "FieldData");
}

// Geometry written once, fields appended per step: every step's offsets name the same geometry.
void vtkhdf_write_static_once(const std::string& rPath, std::size_t NSteps) {
    vtkhdf_write_partitioned(rPath, {{0.0, 0.0}});
    h5::Hid f(H5Fopen(rPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
    h5::Hid g = h5::open_group(f, "VTKHDF");
    H5Ldelete(g, "PointData", H5P_DEFAULT);
    H5Ldelete(g, "CellData", H5P_DEFAULT);
    std::vector<double> u, p;
    for (std::size_t k = 0; k < NSteps; ++k) {
        for (int i = 0; i < 5; ++i)
            u.push_back(static_cast<double>(k));
        for (int i = 0; i < 2; ++i)
            p.push_back(10.0 * static_cast<double>(k));
    }
    h5::Hid pd = h5::create_group(g, "PointData");
    h5::write_dataset(pd, "u", mt::data_array(u));
    h5::Hid cd = h5::create_group(g, "CellData");
    h5::write_dataset(cd, "p", mt::data_array(p));
    h5::Hid s = h5::create_group(g, "Steps");
    h5::write_attr_int(s, "NSteps", static_cast<std::int64_t>(NSteps));
    std::vector<double> values;
    I64Vec zeros(NSteps, 0), ones(NSteps, 1), pu, pp;
    for (std::size_t k = 0; k < NSteps; ++k) {
        values.push_back(0.5 * static_cast<double>(k));
        pu.push_back(static_cast<std::int64_t>(5 * k));
        pp.push_back(static_cast<std::int64_t>(2 * k));
    }
    h5::write_dataset(s, "Values", mt::data_array(values));
    for (const char* name : {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"})
        h5::write_dataset(s, name, vtkhdf_test_i64(zeros));
    h5::write_dataset(s, "NumberOfParts", vtkhdf_test_i64(ones));
    h5::Hid pdo = h5::create_group(s, "PointDataOffsets");
    h5::write_dataset(pdo, "u", vtkhdf_test_i64(pu));
    h5::Hid cdo = h5::create_group(s, "CellDataOffsets");
    h5::write_dataset(cdo, "p", vtkhdf_test_i64(pp));
}

}  // namespace

// ---------------------------------------------------------------------------
// round trips
// ---------------------------------------------------------------------------
TEST(Vtkhdf, RoundTripsTheSharedFixtures) {
    const std::vector<std::pair<const char*, mt::Mesh (*)()>> meshes = {
        {"line", mt::line_mesh},
        {"tri", mt::tri_mesh},
        {"quad", mt::quad_mesh},
        {"tet", mt::tet_mesh},
        {"hex", mt::hex_mesh},
        {"wedge", mt::wedge_mesh},
        {"triangle6", mt::triangle6_mesh},
        {"quad8", mt::quad8_mesh},
        {"tet10", mt::tet10_mesh},
        {"hex20", mt::hex20_mesh},
        {"wedge15", mt::wedge15_mesh},
        {"pyramid13", mt::pyramid13_mesh},
        {"wedge18", mt::wedge18_mesh},
        {"tri_quad", mt::tri_quad_mesh},
        {"data", mt::data_mesh},
    };
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshioplusplus::write_vtkhdf(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshioplusplus::read_vtkhdf(p); };
        for (const auto& [name, make] : meshes) {
            SCOPED_TRACE(name);
            mt::roundtrip(w, r, make(), ".vtkhdf");
        }
    }
}

TEST(Vtkhdf, GzipOnlyCompressesDatasetsWorthCompressing) {
    mt::Mesh m;
    std::vector<std::vector<double>> pts;
    I64Vec conn;
    for (int i = 0; i < 5000; ++i)
        pts.push_back({static_cast<double>(i), 0, 0});
    for (int i = 0; i + 1 < 5000; ++i) {
        conn.push_back(i);
        conn.push_back(i + 1);
    }
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("line", vtkhdf_test_i64(conn, {4999, 2}));
    VtkhdfTempFile gz, raw;
    meshioplusplus::write_vtkhdf(gz.mPath, m, 4);
    meshioplusplus::write_vtkhdf(raw.mPath, m, -1);
    EXPECT_LT(std::filesystem::file_size(gz.mPath), std::filesystem::file_size(raw.mPath));
    h5::Hid f = h5::open_file_read(gz.mPath);
    h5::Hid g = h5::open_group(f, "VTKHDF");
    h5::Hid d(H5Dopen2(g, "Points", H5P_DEFAULT), H5Dclose);
    h5::Hid dcpl(H5Dget_create_plist(d), H5Pclose);
    EXPECT_EQ(H5Pget_nfilters(dcpl), 1);  // deflate on the big dataset...
    h5::Hid small(H5Dopen2(g, "NumberOfPoints", H5P_DEFAULT), H5Dclose);
    h5::Hid scpl(H5Dget_create_plist(small), H5Pclose);
    EXPECT_EQ(H5Pget_nfilters(scpl), 0);  // ...none on a one-element counter
}

TEST(Vtkhdf, DataDtypesAndShapesSurvive) {
    mt::Mesh m = vtkhdf_test_tets();
    mt::NDArray f32 = mt::NDArray::Uninit(meshioplusplus::DType::Float32, {5});
    for (int i = 0; i < 5; ++i)
        f32.As<float>()[i] = 0.25f * static_cast<float>(i);
    m.AddPointData("f32", std::move(f32));
    m.AddPointData("vec", mt::data_array({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, 3));
    m.AddFieldData("matrix", mt::data_array({1, 2, 3, 4, 5, 6}, 3));
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, m);
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
#if !defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE) && !defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
    // The NATIVE and KRATOS meshes normalize data arrays to Float64/Int64 themselves.
    EXPECT_EQ(back.PointData("f32").Dtype(), meshioplusplus::DType::Float32);
#endif
    EXPECT_EQ(vtkhdf_test_doubles(back.PointData("f32")), vtkhdf_test_doubles(m.PointData("f32")));
    EXPECT_EQ(back.PointData("vec").Shape(), (std::vector<std::size_t>{5, 3}));
    EXPECT_EQ(back.FieldData("matrix").Shape(), (std::vector<std::size_t>{2, 3}));
    EXPECT_EQ(vtkhdf_test_doubles(back.FieldData("matrix")),
              vtkhdf_test_doubles(m.FieldData("matrix")));
}

TEST(Vtkhdf, ArraysAndPointsOnlyNarrowTheRead) {
    mt::Mesh m = vtkhdf_test_tets(0, 3.0);
    m.AddFieldData("g", mt::data_array({1.0}));
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, m);
    meshioplusplus::ReadOptions only_u;
    only_u.mDataArrays = std::vector<std::string>{"u"};
    mt::Mesh a = meshioplusplus::read_vtkhdf(f.mPath, only_u);
    EXPECT_EQ(a.PointDataNames(), std::vector<std::string>{"u"});
    EXPECT_EQ(a.NumCellData(), 0u);
    EXPECT_EQ(a.NumFieldData(), 0u);
    meshioplusplus::ReadOptions geo;
    geo.mPointsOnly = true;
    mt::Mesh b = meshioplusplus::read_vtkhdf(f.mPath, geo);
    EXPECT_EQ(b.NumPointData() + b.NumCellData() + b.NumFieldData(), 0u);
    EXPECT_EQ(b.Cells(0).NumCells(), 2u);
}

// ---------------------------------------------------------------------------
// the on-disk shape: a round trip cannot see any of this
// ---------------------------------------------------------------------------
TEST(Vtkhdf, TypeAttributeIsFixedLengthAscii) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets());
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    h5::Hid a(H5Aopen(g, "Type", H5P_DEFAULT), H5Aclose);
    h5::Hid t(H5Aget_type(a), H5Tclose);
    EXPECT_EQ(H5Tis_variable_str(t), 0);
    EXPECT_EQ(H5Tget_cset(t), H5T_CSET_ASCII);
    EXPECT_EQ(h5::read_attr_string(g, "Type"), "UnstructuredGrid");
}

TEST(Vtkhdf, FileShapeConformance) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_cube_poly(0.0, 1.0));
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    EXPECT_EQ(h5::read_attr_int_array(g, "Version"), (I64Vec{2, 5}));  // polyhedra need 2.5
    const auto n_cells =
        static_cast<std::size_t>(h5::read_dataset(g, "NumberOfCells").As<std::int64_t>()[0]);
    for (const char* name :
         {"Connectivity", "Offsets", "FaceConnectivity", "FaceOffsets", "PolyhedronToFaces",
          "PolyhedronOffsets", "NumberOfPoints", "NumberOfCells"})
        EXPECT_EQ(h5::read_dataset(g, name).Dtype(), meshioplusplus::DType::Int64) << name;
    EXPECT_EQ(h5::read_dataset(g, "Types").Dtype(), meshioplusplus::DType::UInt8);
    const mt::NDArray offsets = h5::read_dataset(g, "Offsets");
    EXPECT_EQ(offsets.Size(), n_cells + 1);
    EXPECT_EQ(vtkhdf_test_ints(offsets)[0], 0);
    EXPECT_EQ(vtkhdf_test_ints(h5::read_dataset(g, "PolyhedronOffsets")), (I64Vec{0, 0, 6, 6}));
    // a polyhedron's Connectivity row is its sorted unique node set
    const I64Vec conn = vtkhdf_test_ints(h5::read_dataset(g, "Connectivity"));
    ASSERT_GE(conn.size(), 12u);
    EXPECT_EQ(I64Vec(conn.begin() + 4, conn.begin() + 12), (I64Vec{0, 1, 2, 3, 4, 5, 6, 7}));
    h5::Hid gcpl(H5Gget_create_plist(g), H5Pclose);
    unsigned crt = 0;
    ASSERT_GE(H5Pget_link_creation_order(gcpl, &crt), 0);
    EXPECT_EQ(crt, unsigned(H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED));
    EXPECT_TRUE(h5::exists(g, "PointData") && h5::exists(g, "CellData") &&
                h5::exists(g, "FieldData"));
}

TEST(Vtkhdf, PlainMeshDeclaresTheOldestCoveringVersion) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets());
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    EXPECT_EQ(h5::read_attr_int_array(g, "Version"), (I64Vec{2, 0}));
}

TEST(Vtkhdf, PinnedVersions) {
    VtkhdfTempFile f;
    using V = meshioplusplus::VtkhdfVersion;
    for (V v : {V{2, 0}, V{1, 0}, V{2, 8}})
        EXPECT_NO_THROW(meshioplusplus::write_vtkhdf(
            f.mPath, vtkhdf_test_tets(), 4, meshioplusplus::VtkhdfType::UnstructuredGrid, v));
    for (V v : {V{2, 9}, V{3, 0}})
        EXPECT_THROW(meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets(), 4,
                                                  meshioplusplus::VtkhdfType::UnstructuredGrid, v),
                     meshioplusplus::WriteError);
    try {
        meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_cube_poly(), 4,
                                     meshioplusplus::VtkhdfType::UnstructuredGrid, V{2, 4});
        FAIL() << "a 2.4 file cannot hold polyhedra";
    } catch (const meshioplusplus::WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("polyhedral cells needs VTKHDF 2.5"),
                  std::string::npos);
    }
    EXPECT_THROW(meshioplusplus::write_vtkhdf(
                     f.mPath, vtkhdf_test_tets(), 4,
                     meshioplusplus::VtkhdfType::PartitionedDataSetCollection, V{2, 0}),
                 meshioplusplus::WriteError);
    EXPECT_THROW(meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets(), 10),
                 meshioplusplus::WriteError);
}

// ---------------------------------------------------------------------------
// reading refusals
// ---------------------------------------------------------------------------
TEST(Vtkhdf, RefusesAnHdf5FileWithoutTheVtkhdfGroup) {
    VtkhdfTempFile f(".hdf");
    {
        h5::Hid file = h5::create_file(f.mPath);
        h5::write_dataset(file, "x", vtkhdf_test_i64({1, 2, 3}));
    }
    try {
        meshioplusplus::read_vtkhdf(f.mPath);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("no /VTKHDF group"), std::string::npos);
    }
}

TEST(Vtkhdf, ReadRejectsANonHdf5File) {
    VtkhdfTempFile f;
    {
        FILE* fp = std::fopen(f.mPath.c_str(), "wb");
        std::fputs("this is not hdf5", fp);
        std::fclose(fp);
    }
    try {
        meshioplusplus::read_vtkhdf(f.mPath);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("cannot open"), std::string::npos);
    }
}

TEST(Vtkhdf, UnsupportedTypesAreNamed) {
    for (const char* type : {"ImageData", "OverlappingAMR", "HyperTreeGrid", "Table"}) {
        VtkhdfTempFile f;
        {
            h5::Hid file = h5::create_file(f.mPath);
            h5::Hid g = h5::create_group(file, "VTKHDF");
            vtkhdf_raw_type(g, type, 2, 8);
        }
        try {
            meshioplusplus::read_vtkhdf(f.mPath);
            FAIL() << type;
        } catch (const meshioplusplus::ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(type), std::string::npos) << type;
        }
    }
}

TEST(Vtkhdf, VersionRanges) {
    auto with_version = [](std::int64_t major, std::int64_t minor) {
        auto f = std::make_shared<VtkhdfTempFile>();
        meshioplusplus::write_vtkhdf(f->mPath, vtkhdf_test_tets());
        h5::Hid file(H5Fopen(f->mPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid g = h5::open_group(file, "VTKHDF");
        H5Adelete(g, "Version");
        h5::write_attr_int_array(g, "Version", {major, minor});
        return f;
    };
    EXPECT_EQ(meshioplusplus::read_vtkhdf(with_version(1, 0)->mPath).Cells(0).NumCells(),
              2u);  // older major reads
    EXPECT_NO_THROW(meshioplusplus::read_vtkhdf(with_version(2, 3)->mPath));
    EXPECT_NO_THROW(
        meshioplusplus::read_vtkhdf(with_version(2, 99)->mPath));  // unknown minor: read, warn
    try {
        meshioplusplus::read_vtkhdf(with_version(3, 0)->mPath);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("unsupported VTKHDF version 3.0"), std::string::npos);
    }
}

TEST(Vtkhdf, PolyhedraAreDecodedWhateverTheDeclaredVersion) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_cube_poly());
    {
        h5::Hid file(H5Fopen(f.mPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid g = h5::open_group(file, "VTKHDF");
        H5Adelete(g, "Version");
        h5::write_attr_int_array(g, "Version", {2, 0});  // a lie: no polyhedra before 2.5
    }
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    ASSERT_EQ(back.NumCellBlocks(), 3u);
    EXPECT_EQ(back.Cells(1).Type(), "polyhedron8");
}

// ---------------------------------------------------------------------------
// polyhedra
// ---------------------------------------------------------------------------
TEST(Vtkhdf, MixedPolyhedronKeepsOrderAndData) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_cube_poly(0.0, 4.0));
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    ASSERT_EQ(back.NumCellBlocks(), 3u);
    EXPECT_EQ(back.Cells(0).Type(), "tetra");
    EXPECT_EQ(back.Cells(1).Type(), "polyhedron8");
    EXPECT_EQ(back.Cells(2).Type(), "tetra");
    EXPECT_EQ(vtkhdf_test_ints(back.Cells(0).Conn()), (I64Vec{4, 5, 6, 8}));
    EXPECT_EQ(vtkhdf_test_ints(back.Cells(2).Conn()), (I64Vec{0, 1, 3, 4}));
    ASSERT_EQ(back.Cells(1).NumFaces(0), 6u);
    const auto face0 = back.Cells(1).Face(0, 0);
    EXPECT_EQ((I64Vec(face0.first, face0.first + face0.second)), (I64Vec{0, 3, 2, 1}));
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 0)), std::vector<double>{4.0});
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 1)), std::vector<double>{5.0});
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 2)), std::vector<double>{6.0});
}

TEST(Vtkhdf, PolyhedraBucketByNodeCountAndKeepCellDataAligned) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_prisms(0.0, 10.0));
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    // the file run [6, 4, 6] buckets to {6: [A, C], 4: [B]}, first-seen order
    ASSERT_EQ(back.NumCellBlocks(), 2u);
    EXPECT_EQ(back.Cells(0).Type(), "polyhedron6");
    EXPECT_EQ(back.Cells(0).NumCells(), 2u);
    EXPECT_EQ(back.Cells(1).Type(), "polyhedron4");
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 0)), (std::vector<double>{10.0, 12.0}));
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 1)), (std::vector<double>{11.0}));
}

// ---------------------------------------------------------------------------
// PolyData
// ---------------------------------------------------------------------------
namespace {
mt::Mesh vtkhdf_test_polydata() {
    mt::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {2, 2, 2}, {3, 3, 3}}));
    m.AddCellBlock("triangle", vtkhdf_test_i64({0, 1, 2}, {1, 3}));
    m.AddCellBlock("vertex", vtkhdf_test_i64({4}, {1, 1}));
    m.AddCellBlock("quad", vtkhdf_test_i64({0, 1, 3, 2}, {1, 4}));
    m.AddCellBlock("line", vtkhdf_test_i64({0, 4}, {1, 2}));
    m.AddCellBlock("polygon", vtkhdf_test_i64({0, 1, 3, 2, 5}, {1, 5}));
    m.AddPointData("u", mt::data_array({0, 1, 2, 3, 4, 5}));
    m.AddCellData("p", {mt::data_array({10}), mt::data_array({20}), mt::data_array({30}),
                        mt::data_array({40}), mt::data_array({50})});
    return m;
}
}  // namespace

TEST(Vtkhdf, PolyDataRegroupsIntoCanonicalOrder) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_polydata(), 4,
                                 meshioplusplus::VtkhdfType::PolyData);
    {
        h5::Hid file = h5::open_file_read(f.mPath);
        h5::Hid g = h5::open_group(file, "VTKHDF");
        EXPECT_EQ(h5::read_attr_string(g, "Type"), "PolyData");
        for (const char* c : {"Vertices", "Lines", "Polygons", "Strips"})
            EXPECT_TRUE(h5::exists(g, c)) << c;
        h5::Hid strips = h5::open_group(g, "Strips");
        EXPECT_EQ(vtkhdf_test_ints(h5::read_dataset(strips, "NumberOfCells")),
                  I64Vec{0});  // present but empty
        h5::Hid polys = h5::open_group(g, "Polygons");
        EXPECT_EQ(vtkhdf_test_ints(h5::read_dataset(polys, "NumberOfCells")), I64Vec{3});
    }
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    std::vector<std::string> types;
    for (std::size_t i = 0; i < back.NumCellBlocks(); ++i)
        types.push_back(back.Cells(i).Type());
    EXPECT_EQ(types, (std::vector<std::string>{"vertex", "line", "triangle", "quad", "polygon"}));
    std::vector<double> tags;
    for (std::size_t i = 0; i < back.NumCellBlocks(); ++i)
        tags.push_back(vtkhdf_test_doubles(back.CellData("p", i))[0]);
    EXPECT_EQ(tags, (std::vector<double>{20, 40, 10, 30, 50}));
}

TEST(Vtkhdf, PolyDataRefusesVolumeCellsByName) {
    VtkhdfTempFile f;
    try {
        meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets(), 4,
                                     meshioplusplus::VtkhdfType::PolyData);
        FAIL();
    } catch (const meshioplusplus::WriteError& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("'tetra'"), std::string::npos);
        EXPECT_NE(what.find("UnstructuredGrid"), std::string::npos);
    }
}

namespace {
void vtkhdf_write_raw_polydata(const std::string& rPath, bool Strips, bool PolyVertex,
                               bool PolyLine) {
    h5::Hid f = h5::create_file(rPath);
    h5::Hid g = h5::create_group_crt(f, "VTKHDF");
    vtkhdf_raw_type(g, "PolyData");
    h5::write_dataset(g, "NumberOfPoints", vtkhdf_test_i64({5}));
    h5::write_dataset(g, "Points", mt::data_array(std::vector<double>(15, 0.0), 3));
    const I64Vec cats[4] = {PolyVertex ? I64Vec{0, 1, 2} : I64Vec{4},
                            PolyLine ? I64Vec{0, 1, 2} : I64Vec{0, 4}, I64Vec{0, 1, 2},
                            Strips ? I64Vec{0, 1, 2, 3} : I64Vec{}};
    const char* names[4] = {"Vertices", "Lines", "Polygons", "Strips"};
    std::vector<double> c;
    for (int i = 0; i < 4; ++i) {
        h5::Hid s = h5::create_group(g, names[i]);
        h5::write_dataset(s, "NumberOfCells", vtkhdf_test_i64({cats[i].empty() ? 0 : 1}));
        h5::write_dataset(s, "NumberOfConnectivityIds",
                          vtkhdf_test_i64({static_cast<std::int64_t>(cats[i].size())}));
        h5::write_dataset(s, "Connectivity", vtkhdf_test_i64(cats[i]));
        h5::write_dataset(
            s, "Offsets",
            vtkhdf_test_i64(cats[i].empty()
                                ? I64Vec{0}
                                : I64Vec{0, static_cast<std::int64_t>(cats[i].size())}));
        if (!cats[i].empty())
            c.push_back(static_cast<double>(c.size()));
    }
    h5::create_group(g, "PointData");
    h5::Hid cd = h5::create_group(g, "CellData");
    h5::write_dataset(cd, "c", mt::data_array(c));
    h5::create_group(g, "FieldData");
}
}  // namespace

TEST(Vtkhdf, PolyDataConstructsWithoutAMeshioTypeAreRefusedOrSkipped) {
    struct Case {
        bool mStrips, mPolyVertex, mPolyLine;
        const char* mWhat;
        std::size_t mKept;
    };
    for (const Case& c :
         {Case{true, false, false, "triangle-strip", 3}, Case{false, true, false, "poly-vertex", 2},
          Case{false, false, true, "poly-line", 2}}) {
        SCOPED_TRACE(c.mWhat);
        VtkhdfTempFile f;
        vtkhdf_write_raw_polydata(f.mPath, c.mStrips, c.mPolyVertex, c.mPolyLine);
        try {
            meshioplusplus::read_vtkhdf(f.mPath);
            FAIL() << "strict read must refuse " << c.mWhat;
        } catch (const meshioplusplus::ReadError& e) {
            EXPECT_NE(std::string(e.what()).find(c.mWhat), std::string::npos);
        }
        meshioplusplus::ReadOptions lenient;
        lenient.mLenient = true;
        mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, lenient);
        std::size_t n = 0, rows = 0;
        for (std::size_t i = 0; i < m.NumCellBlocks(); ++i)
            n += m.Cells(i).NumCells();
        for (std::size_t i = 0; i < m.CellDataNumBlocks("c"); ++i)
            rows += m.CellData("c", i).Size();
        EXPECT_EQ(n, c.mKept);     // the unrepresentable cell is dropped...
        EXPECT_EQ(rows, c.mKept);  // ...with its cell-data row
    }
}

// ---------------------------------------------------------------------------
// composites
// ---------------------------------------------------------------------------
namespace {
mt::Mesh vtkhdf_test_region_mesh() {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 1, 1},
                                    {5, 0, 0},
                                    {6, 0, 0},
                                    {5, 1, 0},
                                    {5, 0, 1}}));
    m.AddCellBlock("tetra", vtkhdf_test_i64({0, 1, 2, 3, 1, 2, 3, 4, 5, 6, 7, 8}, {3, 4}));
    m.AddPointData("u", mt::data_array({0, 1, 2, 3, 4, 5, 6, 7, 8}));
    m.AddCellData("p", {mt::data_array({1, 2, 3})});
    m.AddFieldData("g", mt::data_array({9.0}));
    // NON-alphabetical block order, carried by the tags: blocks are written in (tag, name) order
    m.AddRegion(meshioplusplus::Region("zeta", meshioplusplus::RegionKind::Cell, -1, 0,
                                       vtkhdf_test_i64({2})));
    m.AddRegion(meshioplusplus::Region("alpha", meshioplusplus::RegionKind::Cell, -1, 1,
                                       vtkhdf_test_i64({0, 1})));
    return m;
}
}  // namespace

TEST(Vtkhdf, CompositeRegionsRoundTripInBlockOrder) {
    for (auto type : {meshioplusplus::VtkhdfType::PartitionedDataSetCollection,
                      meshioplusplus::VtkhdfType::MultiBlockDataSet}) {
        VtkhdfTempFile f;
        meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_region_mesh(), 4, type);
        mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
        ASSERT_EQ(back.NumRegions(), 2u);
        // the mesh stores regions by name (alpha, zeta); the tags carry the block order (0, 1)
        EXPECT_EQ(back.Region(0).mName, "alpha");
        EXPECT_EQ(back.Region(0).mTag, 1);
        EXPECT_EQ(back.Region(1).mName, "zeta");
        EXPECT_EQ(back.Region(1).mTag, 0);
        // blocks were read in Assembly order (zeta, then alpha), so cells follow it
        EXPECT_EQ(vtkhdf_test_ints(back.Region(1).mEntries), (I64Vec{0}));
        EXPECT_EQ(vtkhdf_test_ints(back.Region(0).mEntries), (I64Vec{1, 2}));
        EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 0)), (std::vector<double>{3, 1, 2}));
        EXPECT_EQ(vtkhdf_test_doubles(back.FieldData("g")), std::vector<double>{9.0});
        EXPECT_EQ(back.NumPoints(), 4u + 5u);  // each block pruned to its own points
    }
}

TEST(Vtkhdf, CompositeFileShape) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_region_mesh(), 4,
                                 meshioplusplus::VtkhdfType::PartitionedDataSetCollection);
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    EXPECT_EQ(h5::read_attr_string(g, "Type"), "PartitionedDataSetCollection");
    EXPECT_EQ(h5::read_attr_int_array(g, "Version"), (I64Vec{2, 1}));
    // creation order, and nothing but blocks and the Assembly at the root: vtkHDFReader
    // reads any other root group as a block and then fails
    EXPECT_EQ(h5::group_links_crt(g), (std::vector<std::string>{"zeta", "alpha", "Assembly"}));
    h5::Hid zeta = h5::open_group(g, "zeta");
    EXPECT_EQ(h5::read_attr_int(zeta, "Index"), 0);
    EXPECT_EQ(h5::read_attr_string(zeta, "Type"), "UnstructuredGrid");
    EXPECT_TRUE(h5::exists(zeta, "FieldData/g"));  // field_data rides on each block
    h5::Hid asm_ = h5::open_group(g, "Assembly");
    h5::Hid node = h5::open_group(asm_, "zeta");
    EXPECT_TRUE(h5::is_soft_link(node, "zeta"));
    EXPECT_EQ(h5::soft_link_target(node, "zeta"), "/VTKHDF/zeta");
    for (hid_t grp : {hid_t(g), hid_t(asm_)}) {
        h5::Hid gcpl(H5Gget_create_plist(grp), H5Pclose);
        unsigned crt = 0;
        ASSERT_GE(H5Pget_link_creation_order(gcpl, &crt), 0);
        EXPECT_EQ(crt, unsigned(H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED));
    }
}

TEST(Vtkhdf, MultiBlockHasNoIndexAndTopLevelLinks) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_region_mesh(), 4,
                                 meshioplusplus::VtkhdfType::MultiBlockDataSet);
    h5::Hid file = h5::open_file_read(f.mPath);
    h5::Hid g = h5::open_group(file, "VTKHDF");
    h5::Hid zeta = h5::open_group(g, "zeta");
    EXPECT_FALSE(h5::has_attr(zeta, "Index"));
    h5::Hid asm_ = h5::open_group(g, "Assembly");
    EXPECT_TRUE(h5::is_soft_link(asm_, "zeta"));
}

TEST(Vtkhdf, CompositeWithoutRegionsWritesOneBlockPerCellBlock) {
    mt::Mesh m = mt::tri_quad_mesh();
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, m, 4,
                                 meshioplusplus::VtkhdfType::PartitionedDataSetCollection);
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    ASSERT_EQ(back.NumRegions(), m.NumCellBlocks());
    std::size_t covered = 0;
    for (std::size_t i = 0; i < back.NumRegions(); ++i)
        covered += back.Region(i).NumEntries();
    std::size_t cells = 0;
    for (std::size_t i = 0; i < m.NumCellBlocks(); ++i)
        cells += m.Cells(i).NumCells();
    EXPECT_EQ(covered, cells);
}

TEST(Vtkhdf, RegionsThatDoNotPartitionTheCellsFallBackToBlocks) {
    mt::Mesh m = vtkhdf_test_region_mesh();
    m.AddRegion(meshioplusplus::Region("only", meshioplusplus::RegionKind::Cell, -1, 5,
                                       vtkhdf_test_i64({0})));
    // overlapping now: cell 0 sits in both "alpha" and "only" -> one block per cell block
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, m, 4, meshioplusplus::VtkhdfType::MultiBlockDataSet);
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    ASSERT_EQ(back.NumRegions(), 1u);
    EXPECT_EQ(back.Region(0).mName, "block_0");
}

TEST(Vtkhdf, EmptyCompositeIsRefused) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from(std::vector<std::vector<double>>{}));
    VtkhdfTempFile f;
    EXPECT_THROW(
        meshioplusplus::write_vtkhdf(f.mPath, m, 4, meshioplusplus::VtkhdfType::MultiBlockDataSet),
        meshioplusplus::WriteError);
}

TEST(Vtkhdf, TransientCompositesAreRefusedByName) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_region_mesh(), 4,
                                 meshioplusplus::VtkhdfType::PartitionedDataSetCollection);
    {
        h5::Hid file(H5Fopen(f.mPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        h5::Hid z = h5::open_group(file, "VTKHDF/zeta");
        h5::Hid s = h5::create_group(z, "Steps");
        h5::write_attr_int(s, "NSteps", 3);
    }
    try {
        meshioplusplus::read_vtkhdf(f.mPath);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("transient composite"), std::string::npos);
    }
}

TEST(Vtkhdf, CompositePieceSelection) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_region_mesh(), 4,
                                 meshioplusplus::VtkhdfType::PartitionedDataSetCollection);
    meshioplusplus::ReadOptions o;
    o.mPieceSet = true;
    o.mPiece = 1;  // alpha
    mt::Mesh one = meshioplusplus::read_vtkhdf(f.mPath, o);
    EXPECT_EQ(one.Cells(0).NumCells(), 2u);
    EXPECT_EQ(vtkhdf_test_doubles(one.CellData("p", 0)), (std::vector<double>{1, 2}));
    EXPECT_EQ(one.NumRegions(), 0u);
    o.mPiece = -1;
    EXPECT_EQ(vtkhdf_test_doubles(meshioplusplus::read_vtkhdf(f.mPath, o).CellData("p", 0)),
              (std::vector<double>{1, 2}));
    o.mPiece = 2;
    try {
        meshioplusplus::read_vtkhdf(f.mPath, o);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("2 piece"), std::string::npos);
    }
}

// ---------------------------------------------------------------------------
// partitions
// ---------------------------------------------------------------------------
TEST(Vtkhdf, PartitionsMergeWithOneRegionPerPiece) {
    VtkhdfTempFile f;
    vtkhdf_write_partitioned(f.mPath, {{0.0, 1.0}, {5.0, 2.0}});
    mt::Mesh back = meshioplusplus::read_vtkhdf(f.mPath);
    EXPECT_EQ(back.NumPoints(), 10u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);  // adjacent same-type blocks join
    EXPECT_EQ(vtkhdf_test_ints(back.Cells(0).Conn()),
              (I64Vec{0, 1, 2, 3, 1, 2, 3, 4, 5, 6, 7, 8, 6, 7, 8, 9}));
    ASSERT_EQ(back.NumRegions(), 2u);
    EXPECT_EQ(back.Region(0).mName, "piece_0");
    EXPECT_EQ(vtkhdf_test_ints(back.Region(0).mEntries), (I64Vec{0, 1}));
    EXPECT_EQ(vtkhdf_test_ints(back.Region(1).mEntries), (I64Vec{2, 3}));
    EXPECT_EQ(vtkhdf_test_doubles(back.PointData("u")),
              (std::vector<double>{1, 1, 1, 1, 1, 2, 2, 2, 2, 2}));
    EXPECT_EQ(vtkhdf_test_doubles(back.CellData("p", 0)), (std::vector<double>{1, 1.5, 2, 2.5}));
}

TEST(Vtkhdf, PartitionPieceSwitch) {
    VtkhdfTempFile f;
    vtkhdf_write_partitioned(f.mPath, {{0.0, 1.0}, {5.0, 2.0}, {9.0, 3.0}});
    struct Case {
        std::int64_t mPiece;
        double mX, mU;
    };
    for (const Case& c : {Case{0, 0.0, 1.0}, Case{1, 5.0, 2.0}, Case{2, 9.0, 3.0},
                          Case{-1, 9.0, 3.0}, Case{-3, 0.0, 1.0}}) {
        meshioplusplus::ReadOptions o;
        o.mPieceSet = true;
        o.mPiece = c.mPiece;
        mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, o);
        EXPECT_EQ(read_double(m.Points(), 0), c.mX) << c.mPiece;
        EXPECT_EQ(read_double(m.PointData("u"), 0), c.mU) << c.mPiece;
        EXPECT_EQ(vtkhdf_test_ints(m.Cells(0).Conn()),
                  (I64Vec{0, 1, 2, 3, 1, 2, 3, 4}));  // the piece's own ids
        EXPECT_EQ(m.NumRegions(), 0u);
    }
    for (std::int64_t bad : {3, -4, 99}) {
        meshioplusplus::ReadOptions o;
        o.mPieceSet = true;
        o.mPiece = bad;
        try {
            meshioplusplus::read_vtkhdf(f.mPath, o);
            FAIL() << bad;
        } catch (const meshioplusplus::ReadError& e) {
            EXPECT_NE(std::string(e.what()).find("3 piece"), std::string::npos);
        }
    }
}

// ---------------------------------------------------------------------------
// time
// ---------------------------------------------------------------------------
TEST(Vtkhdf, TimeStepSelectionOnStaticGeometryWrittenOnce) {
    VtkhdfTempFile f;
    vtkhdf_write_static_once(f.mPath, 3);
    const std::vector<std::pair<int, double>> steps = {
        {0, 0.0}, {1, 1.0}, {2, 2.0}, {-1, 2.0}, {-3, 0.0}};
    for (const auto& [step, want] : steps) {
        meshioplusplus::ReadOptions o;
        o.mTimeStep = step;
        mt::Mesh m = meshioplusplus::read_vtkhdf(f.mPath, o);
        EXPECT_EQ(m.NumPoints(), 5u) << step;
        EXPECT_EQ(vtkhdf_test_ints(m.Cells(0).Conn()), (I64Vec{0, 1, 2, 3, 1, 2, 3, 4})) << step;
        EXPECT_EQ(vtkhdf_test_doubles(m.PointData("u")), std::vector<double>(5, want)) << step;
        EXPECT_EQ(vtkhdf_test_doubles(m.CellData("p", 0)), std::vector<double>(2, 10 * want))
            << step;
    }
    meshioplusplus::ReadOptions o;
    o.mTimeStep = 1;
    EXPECT_EQ(
        vtkhdf_test_doubles(
            meshioplusplus::read_vtkhdf(f.mPath, o).FieldData(meshioplusplus::kSequenceTimeKey)),
        std::vector<double>{0.5});
    for (int bad : {3, -4}) {
        o.mTimeStep = bad;
        try {
            meshioplusplus::read_vtkhdf(f.mPath, o);
            FAIL() << bad;
        } catch (const meshioplusplus::ReadError& e) {
            EXPECT_NE(std::string(e.what()).find("3 step"), std::string::npos);
        }
    }
}

TEST(Vtkhdf, StaticFileHasOneStep) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_tets(0.0, 7.0));
    meshioplusplus::ReadOptions o;
    o.mTimeStep = -1;
    EXPECT_EQ(read_double(meshioplusplus::read_vtkhdf(f.mPath, o).PointData("u"), 0), 7.0);
    o.mTimeStep = 1;
    try {
        meshioplusplus::read_vtkhdf(f.mPath, o);
        FAIL();
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("1 step"), std::string::npos);
    }
}

TEST(Vtkhdf, MetadataReportsStepsWithoutAFullRead) {
    VtkhdfTempFile f;
    vtkhdf_write_static_once(f.mPath, 4);
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_vtkhdf_metadata(f.mPath, {});
    EXPECT_FALSE(meta.mFellBackToFullRead);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.0, 0.5, 1.0, 1.5}));
    EXPECT_EQ(meta.mNumPoints, 5u);
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "tetra");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 2u);
    EXPECT_EQ(meta.mPointDataNames, std::vector<std::string>{"u"});
    EXPECT_EQ(meta.mCellDataNames, std::vector<std::string>{"p"});
}

TEST(Vtkhdf, MetadataFallsBackForPolyhedraAndStillReportsTime) {
    VtkhdfTempFile f;
    meshioplusplus::write_vtkhdf(f.mPath, vtkhdf_test_cube_poly());
    meshioplusplus::MeshMetadata meta = meshioplusplus::read_vtkhdf_metadata(f.mPath, {});
    EXPECT_TRUE(meta.mFellBackToFullRead);
    EXPECT_TRUE(meta.mTimeValues.empty());
    EXPECT_EQ(meta.mNumPoints, 9u);
}

TEST(Vtkhdf, ResolvePieceCountsFromTheEndAndNamesTheCount) {
    meshioplusplus::ReadOptions o;
    EXPECT_THROW(o.ResolvePiece(3), meshioplusplus::ReadError);  // no piece selected
    o.mPieceSet = true;
    const std::vector<std::pair<std::int64_t, std::size_t>> pieces = {
        {0, 0u}, {2, 2u}, {-1, 2u}, {-3, 0u}};
    for (const auto& [piece, want] : pieces) {
        o.mPiece = piece;
        EXPECT_EQ(o.ResolvePiece(3), want) << piece;
    }
    for (std::int64_t bad : {3, -4}) {
        o.mPiece = bad;
        EXPECT_THROW(o.ResolvePiece(3), meshioplusplus::ReadError) << bad;
    }
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5
