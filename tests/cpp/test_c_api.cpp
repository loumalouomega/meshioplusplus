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

/**
 * @file test_c_api.cpp
 * @brief Tests for the C API (bindings/c/). Compiled into the gtest suite
 *        only when MESHIOPLUSPLUS_BUILD_C_API=ON; written purely against the
 *        public C surface (plus mt:: fixtures for reference meshes), so it
 *        runs identically under every mesh backend.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "meshioplusplus/version.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"

// Project includes
#include "meshioplusplus/meshioplusplus.h"

#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"

#include "brep_solids.hpp"
#include "mesh_fixtures.hpp"
#include "gmsh_fixtures.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/region.hpp"
#endif

namespace {

// Deliberately non-square everywhere (5 points x 3 dims, 2 cells x 4 nodes)
// with asymmetric coordinates so a transposed layout cannot cancel out.
const std::vector<double> kPoints = {
    0.0, 0.0, 0.0,  //
    1.1, 0.2, 0.3,  //
    0.4, 1.2, 0.5,  //
    0.6, 0.7, 1.3,  //
    1.4, 1.5, 1.6,  //
};
const std::vector<std::int64_t> kConn = {0, 1, 2, 3, 1, 2, 3, 4};

// Build the reference tet mesh through the C API.
mio_mesh* build_tet_mesh() {
    mio_mesh* m = mio_mesh_create();
    EXPECT_NE(m, nullptr);
    EXPECT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, kPoints.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, kConn.data()), MIO_OK);
    return m;
}

std::string block_type(const mio_mesh* pMesh, std::int64_t block) {
    char buf[64] = {};
    const std::int64_t n = mio_mesh_cell_block_type(pMesh, block, buf, sizeof(buf));
    EXPECT_GE(n, 0);
    return buf;
}

TEST(CApi, VersionAndBackend) {
    EXPECT_STRNE(mio_version(), "");
    const std::string backend = mio_mesh_backend();
    EXPECT_TRUE(backend == "meshio" || backend == "native" || backend == "kratos") << backend;
}

TEST(CApi, VtkPiecesAndLegacyStructuredWithoutPython) {
    const std::string xml =
        "<VTKFile type='ImageData'><ImageData WholeExtent='0 2 0 1 0 1'>"
        "<Piece Extent='0 1 0 1 0 1'/><Piece Extent='1 2 0 1 0 1'/></ImageData></VTKFile>";
    const std::string legacy =
        "# vtk DataFile Version 4.2\nindependent\nASCII\n"
        "DATASET STRUCTURED_POINTS\nDIMENSIONS 3 2 2\nORIGIN 0 0 0\nSPACING 1 1 1\n";
    for (const std::string format : {"vti", "vtk"}) {
        const auto path = mt::temp_path("." + format);
        {
            auto out = meshioplusplus::detail::make_classic_ofstream(path);
            out << (format == "vti" ? xml : legacy);
        }
        mio_mesh* mesh = mio_read(path.c_str(), format.c_str());
        ASSERT_NE(mesh, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_num_points(mesh), format == "vti" ? 16 : 12);
        const auto target = mt::temp_path(".vtu");
        EXPECT_EQ(mio_write(target.c_str(), mesh, "vtu"), MIO_OK) << mio_last_error();
        mio_mesh_free(mesh);
        std::remove(path.c_str());
        std::remove(target.c_str());
    }
}

TEST(CApi, AppendedVtkXmlReadsWithoutPython) {
    for (const std::string format : {"vtp", "vts", "vtr", "vti"}) {
        const bool surface = format == "vtp";
        const std::size_t n = surface ? 4 : 8;
        const std::string type = surface           ? "PolyData"
                                 : format == "vts" ? "StructuredGrid"
                                 : format == "vtr" ? "RectilinearGrid"
                                                   : "ImageData";
        for (bool base64 : {false, true}) {
            std::string xml =
                "<VTKFile type='" + type + "' byte_order='LittleEndian'><" + type +
                (surface ? "" : " WholeExtent='0 1 0 1 0 1'") + "><Piece " +
                (surface ? "NumberOfPoints='4'" : "Extent='0 1 0 1 0 1'") +
                "><PointData><DataArray type='Int32' Name='tag' format='appended' offset='0'/>"
                "</PointData>";
            if (surface || format == "vts") {
                xml +=
                    "<Points><DataArray type='Float64' NumberOfComponents='3' format='ascii'>"
                    "0 0 0 1 0 0 0 1 0 1 1 0 ";
                if (!surface)
                    xml += "0 0 1 1 0 1 0 1 1 1 1 1 ";
                xml += "</DataArray></Points>";
            }
            if (surface)
                xml +=
                    "<Polys><DataArray type='Int32' Name='connectivity' format='ascii'>0 1 3 "
                    "2</DataArray>"
                    "<DataArray type='Int32' Name='offsets' format='ascii'>4</DataArray></Polys>";
            if (format == "vtr")
                xml +=
                    "<Coordinates><DataArray type='Float64' format='ascii'>0 1</DataArray>"
                    "<DataArray type='Float64' format='ascii'>0 1</DataArray>"
                    "<DataArray type='Float64' format='ascii'>0 1</DataArray></Coordinates>";
            xml += "</Piece></" + type + "><AppendedData encoding='" + (base64 ? "base64" : "raw") +
                   "'>_";
            std::vector<unsigned char> bytes{static_cast<unsigned char>(n * 4), 0, 0, 0};
            for (std::size_t i = 0; i < n; ++i)
                bytes.insert(bytes.end(), {static_cast<unsigned char>(i + 7), 0, 0, 0});
            if (base64)
                xml += meshioplusplus::detail::b64encode(bytes.data(), bytes.size());
            else
                xml.append(reinterpret_cast<const char*>(bytes.data()), bytes.size());
            xml += "</AppendedData></VTKFile>";
            const auto path = mt::temp_path("_appended." + format);
            {
                auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
                out.write(xml.data(), static_cast<std::streamsize>(xml.size()));
            }
            mio_mesh* mesh = mio_read(path.c_str(), format.c_str());
            ASSERT_NE(mesh, nullptr) << mio_last_error();
            EXPECT_EQ(mio_mesh_num_points(mesh), static_cast<std::int64_t>(n));
            EXPECT_EQ(mio_mesh_num_cell_blocks(mesh), 1);
            mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), format.c_str());
            ASSERT_NE(meta, nullptr) << mio_last_error();
            EXPECT_EQ(mio_read_metadata_fell_back(meta), 0);
            mio_read_metadata_free(meta);
            const auto converted = mt::temp_path("_appended_converted.vtu");
            EXPECT_EQ(mio_write(converted.c_str(), mesh, "vtu"), MIO_OK) << mio_last_error();
            mio_mesh_free(mesh);
            std::remove(path.c_str());
            std::remove(converted.c_str());
        }
    }
}

TEST(CApi, Gmsh40ReadsBothBinaryCountWidthsAndConverts) {
    for (int width : {0, 4, 8}) {
        const std::string path = mt::temp_path("_gmsh40.msh");
        const std::string bytes = mt::gmsh40_fixture(width != 0, width);
        {
            std::ofstream os(path, std::ios::binary);
            os.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        }
        mio_mesh* mesh = mio_read(path.c_str(), "gmsh");
        ASSERT_NE(mesh, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_num_points(mesh), 4);
        ASSERT_EQ(mio_mesh_num_cell_blocks(mesh), 3);
        EXPECT_EQ(block_type(mesh, 1), "triangle");
        const void* conn = nullptr;
        mio_dtype dtype;
        ASSERT_EQ(mio_mesh_cell_block_conn(mesh, 1, &conn, &dtype), MIO_OK);
        ASSERT_EQ(dtype, MIO_INT64);
        const auto* ids = static_cast<const std::int64_t*>(conn);
        EXPECT_EQ(std::vector<std::int64_t>(ids, ids + 6),
                  (std::vector<std::int64_t>{1, 0, 2, 1, 2, 3}));
        mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), "gmsh");
        ASSERT_NE(meta, nullptr) << mio_last_error();
        EXPECT_EQ(mio_read_metadata_fell_back(meta), 1);
        EXPECT_EQ(mio_read_metadata_num_regions(meta), 2);
        mio_read_metadata_free(meta);
        const std::string out = mt::temp_path("_gmsh40.vtu");
        ASSERT_EQ(mio_write(out.c_str(), mesh, "vtu"), MIO_OK) << mio_last_error();
        mio_mesh* back = mio_read(out.c_str(), "vtu");
        ASSERT_NE(back, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_num_points(back), 4);
        EXPECT_EQ(mio_mesh_num_cell_blocks(back), 3);
        mio_mesh_free(back);
        mio_mesh_free(mesh);
        std::remove(out.c_str());
        std::remove(path.c_str());
    }
}

// The compile-time macros describe the HEADER; mio_version() describes the
// linked LIBRARY. With a shared build the two genuinely can differ, which is why
// both exist -- but in this build they are the same artifact, so they must agree,
// and that is what catches a version file bumped without its twin.
TEST(CApi, CompileTimeVersionMatchesTheLinkedLibrary) {
    const std::string expected = std::to_string(MIO_VERSION_MAJOR) + "." +
                                 std::to_string(MIO_VERSION_MINOR) + "." +
                                 std::to_string(MIO_VERSION_PATCH);
    EXPECT_EQ(std::string(mio_version()), expected);
    EXPECT_EQ(std::string(MESHIOPLUSPLUS_VERSION_STRING), expected);

    // The feature-detection macros are what consumers actually write.
    EXPECT_TRUE(MIO_VERSION_AT_LEAST(MIO_VERSION_MAJOR, MIO_VERSION_MINOR, MIO_VERSION_PATCH));
    EXPECT_FALSE(MIO_VERSION_BEFORE(MIO_VERSION_MAJOR, MIO_VERSION_MINOR, MIO_VERSION_PATCH));
    EXPECT_TRUE(MIO_VERSION_AT_LEAST(1, 0, 0));
    EXPECT_FALSE(MIO_VERSION_AT_LEAST(MIO_VERSION_MAJOR + 1, 0, 0));
    // Ordering must not break across a component boundary: 9.6.0 is after 9.5.99.
    EXPECT_TRUE(MESHIOPLUSPLUS_VERSION_AT_LEAST(9, 5, 0));
    EXPECT_TRUE(MESHIOPLUSPLUS_VERSION > (9 * 10000 + 5 * 100 + 99));

    // selective refinement landed in 9.5.0; this is the shape a consumer writes.
#if MIO_VERSION_AT_LEAST(9, 5, 0)
    EXPECT_NE(&mio_refine_ex, nullptr);
#endif
}

TEST(CApi, FormatAvailability) {
    EXPECT_EQ(mio_format_readable("vtu"), 1);
    EXPECT_EQ(mio_format_writable("vtu"), 1);
    EXPECT_EQ(mio_format_readable("openfoam"), 1);
    EXPECT_EQ(mio_format_writable("openfoam"), 1);  // writable since v9.20.0
    EXPECT_EQ(mio_format_readable("nonexistent"), 0);
    EXPECT_EQ(mio_format_readable(nullptr), 0);
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    EXPECT_EQ(mio_format_readable("med"), 1);
    EXPECT_EQ(mio_format_writable("med"), 1);
#else
    EXPECT_EQ(mio_format_readable("med"), 0);
#endif
}

TEST(CApi, CellTypeMetadata) {
    EXPECT_STREQ(mio_cell_type_name(MIO_CELL_Tetra10), "tetra10");
    EXPECT_EQ(mio_cell_type_from_name("tetra10"), MIO_CELL_Tetra10);
    EXPECT_EQ(mio_cell_type_from_name("not_a_type"), MIO_CELL_Custom);
    EXPECT_EQ(mio_cell_type_from_name(nullptr), MIO_CELL_Custom);
    EXPECT_EQ(mio_cell_type_num_nodes(MIO_CELL_Hexahedron20), 20);
    EXPECT_EQ(mio_cell_type_num_nodes(MIO_CELL_Polygon), -1);
    EXPECT_EQ(mio_cell_type_dimension(MIO_CELL_Triangle), 2);
    EXPECT_EQ(mio_cell_type_dimension(MIO_CELL_Custom), -1);
    EXPECT_STREQ(mio_cell_type_name(MIO_CELL_Custom), "");
    // The VTK Bezier family (ABI 26): variable size, a fixed dimension.
    EXPECT_STREQ(mio_cell_type_name(MIO_CELL_VtkBezierQuadrilateral), "VTK_BEZIER_QUADRILATERAL");
    EXPECT_EQ(mio_cell_type_from_name("VTK_BEZIER_PYRAMID"), MIO_CELL_VtkBezierPyramid);
    EXPECT_EQ(mio_cell_type_num_nodes(MIO_CELL_VtkBezierHexahedron), -1);
    EXPECT_EQ(mio_cell_type_dimension(MIO_CELL_VtkBezierCurve), 1);
    EXPECT_EQ(mio_cell_type_dimension(MIO_CELL_VtkBezierWedge), 3);
}

TEST(CApi, BuildAndInspect) {
    mio_mesh* m = build_tet_mesh();

    EXPECT_EQ(mio_mesh_num_points(m), 5);
    EXPECT_EQ(mio_mesh_point_dim(m), 3);
    EXPECT_EQ(mio_mesh_num_cell_blocks(m), 1);
    EXPECT_EQ(block_type(m, 0), "tetra");

    std::int64_t num_cells = 0, npc = 0;
    std::int32_t ragged = -1;
    ASSERT_EQ(mio_mesh_cell_block_info(m, 0, &num_cells, &npc, &ragged), MIO_OK);
    EXPECT_EQ(num_cells, 2);
    EXPECT_EQ(npc, 4);
    EXPECT_EQ(ragged, 0);

    const void* pts = nullptr;
    mio_dtype dt = MIO_FLOAT32;
    ASSERT_EQ(mio_mesh_get_points(m, &pts, &dt), MIO_OK);
    ASSERT_EQ(dt, MIO_FLOAT64);
    const double* d = static_cast<const double*>(pts);
    for (std::size_t i = 0; i < kPoints.size(); ++i)
        EXPECT_DOUBLE_EQ(d[i], kPoints[i]);

    const void* conn = nullptr;
    ASSERT_EQ(mio_mesh_cell_block_conn(m, 0, &conn, &dt), MIO_OK);
    ASSERT_EQ(dt, MIO_INT64);
    const std::int64_t* c = static_cast<const std::int64_t*>(conn);
    for (std::size_t i = 0; i < kConn.size(); ++i)
        EXPECT_EQ(c[i], kConn[i]);

    mio_mesh_free(m);
}

TEST(CApi, Int32ConnectivityWidens) {
    const std::vector<std::int32_t> conn32 = {0, 1, 2, 3, 1, 2, 3, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, kPoints.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT32, conn32.data()), MIO_OK);
    const void* conn = nullptr;
    mio_dtype dt = MIO_FLOAT32;
    ASSERT_EQ(mio_mesh_cell_block_conn(m, 0, &conn, &dt), MIO_OK);
    EXPECT_EQ(dt, MIO_INT64);
    const std::int64_t* c = static_cast<const std::int64_t*>(conn);
    for (std::size_t i = 0; i < conn32.size(); ++i)
        EXPECT_EQ(c[i], conn32[i]);
    mio_mesh_free(m);
}

TEST(CApi, NamedDataRoundTrip) {
    mio_mesh* m = build_tet_mesh();

    const std::vector<double> temperature = {1.0, 2.0, 3.0, 4.0, 5.0};
    const std::int64_t shape1[] = {5};
    ASSERT_EQ(mio_mesh_add_point_data(m, "temperature", MIO_FLOAT64, 1, shape1, temperature.data()),
              MIO_OK);
    std::vector<double> velocity(15);
    for (std::size_t i = 0; i < velocity.size(); ++i)
        velocity[i] = 0.5 * static_cast<double>(i);
    const std::int64_t shape2[] = {5, 3};
    ASSERT_EQ(mio_mesh_add_point_data(m, "velocity", MIO_FLOAT64, 2, shape2, velocity.data()),
              MIO_OK);

    const std::vector<double> quality = {0.5, 0.75};
    const std::int64_t shapec[] = {2};
    ASSERT_EQ(mio_mesh_append_cell_data(m, "quality", MIO_FLOAT64, 1, shapec, quality.data()),
              MIO_OK);

    const std::vector<double> gravity = {0.0, 0.0, -9.81};
    const std::int64_t shapef[] = {3};
    ASSERT_EQ(mio_mesh_add_field_data(m, "gravity", MIO_FLOAT64, 1, shapef, gravity.data()),
              MIO_OK);

    // Names come back sorted on every backend.
    EXPECT_EQ(mio_mesh_num_point_data(m), 2);
    char buf[64] = {};
    ASSERT_GE(mio_mesh_point_data_name(m, 0, buf, sizeof(buf)), 0);
    EXPECT_STREQ(buf, "temperature");
    ASSERT_GE(mio_mesh_point_data_name(m, 1, buf, sizeof(buf)), 0);
    EXPECT_STREQ(buf, "velocity");

    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT32;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(m, "velocity", &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 5);
    EXPECT_EQ(shape[1], 3);
    const double* v = static_cast<const double*>(data);
    for (std::size_t i = 0; i < velocity.size(); ++i)
        EXPECT_DOUBLE_EQ(v[i], velocity[i]);

    EXPECT_EQ(mio_mesh_num_cell_data(m), 1);
    EXPECT_EQ(mio_mesh_cell_data_num_blocks(m, "quality"), 1);
    ASSERT_EQ(mio_mesh_get_cell_data(m, "quality", 0, &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(ndim, 1);
    EXPECT_EQ(shape[0], 2);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[1], 0.75);

    EXPECT_EQ(mio_mesh_num_field_data(m), 1);
    ASSERT_GE(mio_mesh_field_data_name(m, 0, buf, sizeof(buf)), 0);
    EXPECT_STREQ(buf, "gravity");
    ASSERT_EQ(mio_mesh_get_field_data(m, "gravity", &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(shape[0], 3);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[2], -9.81);

    mio_mesh_free(m);
}

TEST(CApi, FileRoundTripAndConvert) {
    const std::string vtu = mt::temp_path("_capi.vtu");
    const std::string vtk = mt::temp_path("_capi.vtk");

    mio_mesh* m = build_tet_mesh();
    ASSERT_EQ(mio_write(vtu.c_str(), m, nullptr), MIO_OK);
    mio_mesh_free(m);

    mio_mesh* r = mio_read(vtu.c_str(), nullptr);
    ASSERT_NE(r, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(r), 5);
    EXPECT_EQ(mio_mesh_num_cell_blocks(r), 1);
    EXPECT_EQ(block_type(r, 0), "tetra");
    const void* pts = nullptr;
    mio_dtype dt;
    ASSERT_EQ(mio_mesh_get_points(r, &pts, &dt), MIO_OK);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(pts)[4], 0.2);  // point 1, y
    mio_mesh_free(r);

    ASSERT_EQ(mio_convert(vtu.c_str(), nullptr, vtk.c_str(), nullptr), MIO_OK);
    mio_mesh* r2 = mio_read(vtk.c_str(), "vtk");
    ASSERT_NE(r2, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(r2), 5);
    mio_mesh_free(r2);

    std::remove(vtu.c_str());
    std::remove(vtk.c_str());
}

TEST(CApi, MeshOperations) {
    mio_mesh* m = build_tet_mesh();  // 2 tetra sharing a face

    // extract_surface -> boundary triangles.
    mio_mesh* surf = mio_extract_surface(m, /*record_parent_ids=*/1);
    ASSERT_NE(surf, nullptr) << mio_last_error();
    EXPECT_GT(mio_mesh_num_cell_blocks(surf), 0);
    EXPECT_EQ(block_type(surf, 0), "triangle");
    mio_mesh_free(surf);

    // extract_skin -> also boundary triangles.
    mio_mesh* skin = mio_extract_skin(m, /*linearize=*/0);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(block_type(skin, 0), "triangle");
    mio_mesh_free(skin);

    // attach_quality preserves the cell block(s) and adds cell_data.
    mio_mesh* q = mio_attach_quality(m);
    ASSERT_NE(q, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_cell_blocks(q), 1);
    EXPECT_EQ(block_type(q, 0), "tetra");
    EXPECT_GT(mio_mesh_num_cell_data(q), 0);
    mio_mesh_free(q);

    // quality counts.
    std::int64_t nc = -1, ninv = -1, ndeg = -1;
    ASSERT_EQ(mio_quality_counts(m, &nc, &ninv, &ndeg), MIO_OK) << mio_last_error();
    EXPECT_EQ(nc, 2);
    EXPECT_GE(ninv, 0);
    EXPECT_GE(ndeg, 0);

    // bandwidth + reorder (RCM): result carries a renumbered mesh + node perm.
    EXPECT_GE(mio_compute_bandwidth(m), 0);
    mio_reorder_result* res = mio_reorder(m, "rcm");
    ASSERT_NE(res, nullptr) << mio_last_error();
    const mio_mesh* rmesh = mio_reorder_result_mesh(res);
    ASSERT_NE(rmesh, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(rmesh), 1);
    const void* np = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_reorder_result_node_perm(res, &np, &dt, &n), MIO_OK) << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(n, 5);  // build_tet_mesh has 5 points
    // node permutation is a bijection over [0, n).
    {
        const std::int64_t* p = static_cast<const std::int64_t*>(np);
        std::vector<int> seen(static_cast<std::size_t>(n), 0);
        for (std::int64_t i = 0; i < n; ++i) {
            ASSERT_GE(p[i], 0);
            ASSERT_LT(p[i], n);
            EXPECT_EQ(seen[static_cast<std::size_t>(p[i])], 0);
            seen[static_cast<std::size_t>(p[i])] = 1;
        }
    }
    EXPECT_EQ(mio_reorder_result_num_cell_perms(res), 1);
    // Ownership transfer keeps the mesh alive after freeing the result.
    mio_mesh* taken = mio_reorder_result_take_mesh(res);
    ASSERT_NE(taken, nullptr) << mio_last_error();
    mio_reorder_result_free(res);
    EXPECT_EQ(mio_mesh_num_cell_blocks(taken), 1);
    mio_mesh_free(taken);
    // Unknown method fails cleanly.
    EXPECT_EQ(mio_reorder(m, "bogus"), nullptr);

    mio_mesh_free(m);
}

TEST(CApi, Decimate) {
    // A 4-triangle fan around a centre vertex: the centre collapses into the
    // pinned boundary, leaving 2 triangles.
    const std::array<double, 15> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 0};
    const std::array<std::int64_t, 12> conn = {0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 4, 3, MIO_INT64, conn.data()), MIO_OK);

    mio_decimate_result* res =
        mio_decimate(m, /*target_ratio=*/-1.0, /*target_faces=*/1, /*max_error=*/-1.0,
                     /*placement=*/nullptr, /*preserve_boundary=*/1, /*preserve_features=*/1,
                     /*feature_angle=*/30.0);
    ASSERT_NE(res, nullptr) << mio_last_error();
    const mio_mesh* dm = mio_decimate_result_mesh(res);
    ASSERT_NE(dm, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(dm), 1);
    EXPECT_EQ(block_type(dm, 0), "triangle");
    EXPECT_EQ(mio_decimate_result_faces_removed(res), 2);
    EXPECT_EQ(mio_decimate_result_points_removed(res), 1);
    EXPECT_GE(mio_decimate_result_collapses_rejected(res), 0);
    EXPECT_GE(mio_decimate_result_max_error_applied(res), 0.0);

    // The point map lands every input point on a live output index (the
    // collapsed centre maps to its survivor, not -1).
    const void* pm = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_decimate_result_point_map(res, &pm, &dt, &n), MIO_OK) << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(n, 5);
    {
        const std::int64_t* p = static_cast<const std::int64_t*>(pm);
        for (std::int64_t i = 0; i < n; ++i) {
            EXPECT_GE(p[i], 0);
            EXPECT_LT(p[i], 4);
        }
    }
    ASSERT_EQ(mio_decimate_result_num_cell_maps(res), 1);
    const void* cm = nullptr;
    ASSERT_EQ(mio_decimate_result_cell_map(res, 0, &cm, &dt, &n), MIO_OK) << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(n, 4);
    EXPECT_EQ(mio_decimate_result_cell_map(res, 7, &cm, &dt, &n), MIO_ERR_NOT_FOUND);

    // Ownership transfer keeps the mesh alive after freeing the result.
    mio_mesh* taken = mio_decimate_result_take_mesh(res);
    ASSERT_NE(taken, nullptr) << mio_last_error();
    mio_decimate_result_free(res);
    EXPECT_EQ(mio_mesh_num_cell_blocks(taken), 1);
    mio_mesh_free(taken);

    // Error paths fail cleanly: no criterion, a bad placement, a volume mesh.
    EXPECT_EQ(mio_decimate(m, -1.0, -1, -1.0, nullptr, 1, 1, 30.0), nullptr);
    EXPECT_EQ(mio_decimate(m, -1.0, 1, -1.0, "nearest", 1, 1, 30.0), nullptr);
    mio_mesh_free(m);
    mio_mesh* tet = build_tet_mesh();
    EXPECT_EQ(mio_decimate(tet, 0.5, -1, -1.0, nullptr, 1, 1, 30.0), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("extract_surface"), std::string::npos);
    mio_mesh_free(tet);
}

TEST(CApi, DecimateExFrozenMaskAndNullOpts) {
    // Same 4-triangle fan as Decimate: without frozen, the centre collapses.
    const std::array<double, 15> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 0};
    const std::array<std::int64_t, 12> conn = {0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 4, 3, MIO_INT64, conn.data()), MIO_OK);

    // Unlike mio_smooth_ex, there is no "all defaults" criterion, so NULL
    // opts is a clean error, not mio_decimate's own defaults.
    EXPECT_EQ(mio_decimate_ex(m, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("opts is NULL"), std::string::npos);

    // Freezing the centre vertex (id 4) must prevent its own collapse, even
    // though it is exactly the vertex the unfrozen Decimate test above
    // removes.
    mio_decimate_opts opts;
    mio_decimate_opts_init(&opts);
    opts.target_faces = 1;
    const std::int64_t frozen_id = 4;
    opts.frozen = &frozen_id;
    opts.num_frozen = 1;
    mio_decimate_result* res = mio_decimate_ex(m, &opts);
    ASSERT_NE(res, nullptr) << mio_last_error();
    EXPECT_EQ(mio_decimate_result_points_removed(res), 0);
    mio_decimate_result_free(res);

    // An out-of-range frozen id fails the call by name.
    const std::int64_t bad_id = 999;
    opts.frozen = &bad_id;
    EXPECT_EQ(mio_decimate_ex(m, &opts), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("frozen node id 999"), std::string::npos);

    mio_mesh_free(m);
}

TEST(CApi, DecimateVolume) {
    // A unit cube split into 6 positively-oriented tets sharing the main
    // diagonal 0-6 (the same fixture as test_decimate_volume.cpp).
    const std::array<double, 24> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                        0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::array<std::int64_t, 24> conn = {0, 1, 2, 6, 0, 2, 3, 6, 0, 3, 7, 6,
                                               0, 7, 4, 6, 0, 4, 5, 6, 0, 5, 1, 6};
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 6, 4, MIO_INT64, conn.data()), MIO_OK);

    mio_decimate_volume_result* res =
        mio_decimate_volume(m, /*target_ratio=*/-1.0, /*target_cells=*/1, /*max_error=*/-1.0,
                            /*placement=*/nullptr, /*preserve_boundary=*/0, /*preserve_features=*/0,
                            /*feature_angle=*/30.0);
    ASSERT_NE(res, nullptr) << mio_last_error();
    const mio_mesh* dm = mio_decimate_volume_result_mesh(res);
    ASSERT_NE(dm, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(dm), 1);
    EXPECT_EQ(block_type(dm, 0), "tetra");
    EXPECT_GT(mio_decimate_volume_result_tets_removed(res), 0);
    EXPECT_GE(mio_decimate_volume_result_points_removed(res), 0);
    EXPECT_GE(mio_decimate_volume_result_collapses_rejected(res), 0);
    EXPECT_GE(mio_decimate_volume_result_max_error_applied(res), 0.0);

    const void* pm = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_decimate_volume_result_point_map(res, &pm, &dt, &n), MIO_OK) << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(n, 8);

    ASSERT_EQ(mio_decimate_volume_result_num_cell_maps(res), 1);
    const void* cm = nullptr;
    ASSERT_EQ(mio_decimate_volume_result_cell_map(res, 0, &cm, &dt, &n), MIO_OK)
        << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(n, 6);
    EXPECT_EQ(mio_decimate_volume_result_cell_map(res, 7, &cm, &dt, &n), MIO_ERR_NOT_FOUND);

    mio_mesh* taken = mio_decimate_volume_result_take_mesh(res);
    ASSERT_NE(taken, nullptr) << mio_last_error();
    mio_decimate_volume_result_free(res);
    EXPECT_EQ(mio_mesh_num_cell_blocks(taken), 1);
    mio_mesh_free(taken);

    // Error paths fail cleanly: no criterion, and decimate() itself still
    // refuses a tet mesh, pointing here.
    EXPECT_EQ(mio_decimate_volume(m, -1.0, -1, -1.0, nullptr, 0, 1, 30.0), nullptr);
    EXPECT_EQ(mio_decimate(m, -1.0, 1, -1.0, nullptr, 1, 1, 30.0), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("extract_surface"), std::string::npos);
    mio_mesh_free(m);

    // A non-tet 3D block is refused by name.
    mio_mesh* hex = mio_mesh_create();
    ASSERT_NE(hex, nullptr);
    const std::array<double, 24> hpts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                         0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::array<std::int64_t, 8> hconn = {0, 1, 2, 3, 4, 5, 6, 7};
    ASSERT_EQ(mio_mesh_set_points(hex, MIO_FLOAT64, 8, 3, hpts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(hex, "hexahedron", 1, 8, MIO_INT64, hconn.data()), MIO_OK);
    EXPECT_EQ(mio_decimate_volume(hex, -1.0, 1, -1.0, nullptr, 0, 1, 30.0), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("tet-only"), std::string::npos);
    mio_mesh_free(hex);
}

TEST(CApi, DecimateVolumeExFrozenMaskAndNullOpts) {
    // Same 6-tet cube as DecimateVolume, sharing the main diagonal 0-6.
    const std::array<double, 24> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                        0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::array<std::int64_t, 24> conn = {0, 1, 2, 6, 0, 2, 3, 6, 0, 3, 7, 6,
                                               0, 7, 4, 6, 0, 4, 5, 6, 0, 5, 1, 6};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 6, 4, MIO_INT64, conn.data()), MIO_OK);

    // As mio_decimate_ex: no "all defaults" criterion, so NULL opts errors.
    EXPECT_EQ(mio_decimate_volume_ex(m, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("opts is NULL"), std::string::npos);

    // Freeze every vertex: no collapse can proceed, so nothing is removed
    // rather than the call failing outright (the same "queue exhausted"
    // shape a fully-pinned mesh already has without frozen).
    mio_decimate_volume_opts opts;
    mio_decimate_volume_opts_init(&opts);
    opts.target_cells = 1;
    const std::array<std::int64_t, 8> all_ids = {0, 1, 2, 3, 4, 5, 6, 7};
    opts.frozen = all_ids.data();
    opts.num_frozen = static_cast<std::int64_t>(all_ids.size());
    mio_decimate_volume_result* res = mio_decimate_volume_ex(m, &opts);
    ASSERT_NE(res, nullptr) << mio_last_error();
    EXPECT_EQ(mio_decimate_volume_result_tets_removed(res), 0);
    mio_decimate_volume_result_free(res);

    // An out-of-range frozen id fails the call by name.
    const std::int64_t bad_id = 999;
    opts.frozen = &bad_id;
    opts.num_frozen = 1;
    EXPECT_EQ(mio_decimate_volume_ex(m, &opts), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("frozen node id 999"), std::string::npos);

    mio_mesh_free(m);
}

TEST(CApi, Merge) {
    mio_mesh* a = build_tet_mesh();  // 5 points, 2 tetra
    mio_mesh* b = build_tet_mesh();
    const mio_mesh* inputs[2] = {a, b};

    // Concatenate: point/cell counts sum; source tag present.
    mio_mesh* cat = mio_merge(inputs, 2, /*weld=*/0, 1e-8, /*source_tag=*/1,
                              /*data_policy=*/0, /*drop_dup=*/0);
    ASSERT_NE(cat, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(cat), 10);
    EXPECT_EQ(mio_mesh_num_cell_blocks(cat), 1);
    EXPECT_GT(mio_mesh_num_cell_data(cat), 0);  // source_mesh_id
    mio_mesh_free(cat);

    // Weld two identical meshes -> the 5 coincident points collapse.
    mio_mesh* welded = mio_merge(inputs, 2, /*weld=*/1, 1e-9, /*source_tag=*/1,
                                 /*data_policy=*/0, /*drop_dup=*/0);
    ASSERT_NE(welded, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(welded), 5);
    mio_mesh_free(welded);

    // Error paths.
    EXPECT_EQ(mio_merge(nullptr, 2, 0, 1e-8, 1, 0, 0), nullptr);
    EXPECT_EQ(mio_merge(inputs, 0, 0, 1e-8, 1, 0, 0), nullptr);

    mio_mesh_free(a);
    mio_mesh_free(b);
}

TEST(CApi, Interpolate) {
    // Source: two points on the x-axis carrying a scalar field.
    mio_mesh* src = mio_mesh_create();
    const std::vector<double> spts = {0, 0, 0, 1, 0, 0};
    ASSERT_EQ(mio_mesh_set_points(src, MIO_FLOAT64, 2, 3, spts.data()), MIO_OK);
    const std::vector<double> f = {10.0, 20.0};
    const std::int64_t shape1[] = {2};
    ASSERT_EQ(mio_mesh_add_point_data(src, "f", MIO_FLOAT64, 1, shape1, f.data()), MIO_OK);

    // Target: one point near the second source point.
    mio_mesh* tgt = mio_mesh_create();
    const std::vector<double> tpts = {0.9, 0, 0};
    ASSERT_EQ(mio_mesh_set_points(tgt, MIO_FLOAT64, 1, 3, tpts.data()), MIO_OK);

    // Default arrays (= all source point_data), default method (nearest).
    mio_mesh* out = mio_interpolate(src, tgt, nullptr, nullptr, 0, 0, 0.0, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT32;
    ASSERT_EQ(mio_mesh_get_point_data(out, "f", &data, &dt, nullptr, nullptr), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[0], 20.0);
    mio_mesh_free(out);

    // An explicit name list goes through the char** + count convention.
    const char* names[] = {"f"};
    out = mio_interpolate(src, tgt, "nearest", names, 1, 0, 0.0, "error");
    ASSERT_NE(out, nullptr) << mio_last_error();
    mio_mesh_free(out);

    // Error paths: bad method / unknown array / NULL meshes.
    EXPECT_EQ(mio_interpolate(src, tgt, "bogus", nullptr, 0, 0, 0.0, nullptr), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    const char* bad[] = {"nope"};
    EXPECT_EQ(mio_interpolate(src, tgt, nullptr, bad, 1, 0, 0.0, nullptr), nullptr);
    EXPECT_EQ(mio_interpolate(nullptr, tgt, nullptr, nullptr, 0, 0, 0.0, nullptr), nullptr);
    EXPECT_EQ(mio_interpolate(src, nullptr, nullptr, nullptr, 0, 0, 0.0, nullptr), nullptr);

    mio_mesh_free(src);
    mio_mesh_free(tgt);
}

TEST(CApi, ConservativeInterpolate) {
    // Two identical single-triangle meshes: the full-value-recovery oracle.
    mio_mesh* src = mio_mesh_create();
    const std::vector<double> spts = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    ASSERT_EQ(mio_mesh_set_points(src, MIO_FLOAT64, 3, 3, spts.data()), MIO_OK);
    const std::vector<std::int64_t> sconn = {0, 1, 2};
    ASSERT_EQ(mio_mesh_add_cell_block(src, "triangle", 1, 3, MIO_INT64, sconn.data()), MIO_OK);
    const std::vector<double> f = {7.0};
    const std::int64_t shapec[] = {1};
    ASSERT_EQ(mio_mesh_append_cell_data(src, "f", MIO_FLOAT64, 1, shapec, f.data()), MIO_OK);

    mio_mesh* tgt = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(tgt, MIO_FLOAT64, 3, 3, spts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(tgt, "triangle", 1, 3, MIO_INT64, sconn.data()), MIO_OK);

    mio_mesh* out = mio_conservative_interpolate(src, tgt, nullptr, 0, 0.0, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT32;
    ASSERT_EQ(mio_mesh_get_cell_data(out, "f", 0, &data, &dt, nullptr, nullptr), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[0], 7.0);
    mio_mesh_free(out);

    // An explicit name list goes through the char** + count convention.
    const char* names[] = {"f"};
    out = mio_conservative_interpolate(src, tgt, names, 1, 0.0, "error");
    ASSERT_NE(out, nullptr) << mio_last_error();
    mio_mesh_free(out);

    // Error paths: unknown array / NULL meshes.
    const char* bad[] = {"nope"};
    EXPECT_EQ(mio_conservative_interpolate(src, tgt, bad, 1, 0.0, nullptr), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_conservative_interpolate(nullptr, tgt, nullptr, 0, 0.0, nullptr), nullptr);
    EXPECT_EQ(mio_conservative_interpolate(src, nullptr, nullptr, 0, 0.0, nullptr), nullptr);

    mio_mesh_free(src);
    mio_mesh_free(tgt);
}

TEST(CApi, Diff) {
    mio_mesh* a = build_tet_mesh();
    mio_mesh* b = build_tet_mesh();

    // Identical meshes -> equal, verdict identical.
    int equal = -1;
    ASSERT_EQ(mio_meshes_equal(a, b, 1e-12, 1e-9, 0, &equal), MIO_OK) << mio_last_error();
    EXPECT_EQ(equal, 1);
    mio_diff_result* res = nullptr;
    ASSERT_EQ(mio_diff(a, b, 1e-12, 1e-9, 0, &res), MIO_OK) << mio_last_error();
    ASSERT_NE(res, nullptr);
    EXPECT_EQ(mio_diff_result_verdict(res), MIO_DIFF_IDENTICAL);
    EXPECT_EQ(mio_diff_result_num_block_diffs(res), 1);
    int type_mism = -1, count_mism = -1;
    std::int64_t conn = -1;
    ASSERT_EQ(mio_diff_result_block(res, 0, &type_mism, &count_mism, &conn), MIO_OK);
    EXPECT_EQ(type_mism, 0);
    EXPECT_EQ(count_mism, 0);
    EXPECT_EQ(conn, 0);
    mio_diff_result_free(res);
    mio_mesh_free(b);

    // A mesh with different connectivity -> different.
    const std::vector<std::int64_t> conn2 = {0, 1, 2, 3, 0, 2, 3, 4};
    mio_mesh* c = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(c, MIO_FLOAT64, 5, 3, kPoints.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(c, "tetra", 2, 4, MIO_INT64, conn2.data()), MIO_OK);
    ASSERT_EQ(mio_meshes_equal(a, c, 1e-12, 1e-9, 0, &equal), MIO_OK);
    EXPECT_EQ(equal, 0);
    res = nullptr;
    ASSERT_EQ(mio_diff(a, c, 1e-12, 1e-9, 0, &res), MIO_OK);
    EXPECT_EQ(mio_diff_result_verdict(res), MIO_DIFF_DIFFERENT);
    ASSERT_EQ(mio_diff_result_block(res, 0, nullptr, nullptr, &conn), MIO_OK);
    EXPECT_EQ(conn, 1);
    mio_diff_result_free(res);
    mio_mesh_free(c);

    // NULL mesh is a clean error.
    EXPECT_EQ(mio_diff(nullptr, a, 1e-12, 1e-9, 0, &res), MIO_ERR_INVALID_ARG);

    mio_mesh_free(a);
}

TEST(CApi, SniffFormat) {
    const std::string vtu = mt::temp_path("_capi_sniff.vtu");
    mio_mesh* m = build_tet_mesh();
    ASSERT_EQ(mio_write(vtu.c_str(), m, nullptr), MIO_OK);
    mio_mesh_free(m);

    char buf[32] = {};
    const std::int64_t n = mio_sniff_format(vtu.c_str(), buf, sizeof(buf));
    EXPECT_EQ(std::string(buf), "vtu");
    EXPECT_EQ(n, 3);
    std::remove(vtu.c_str());
}

TEST(CApi, ZeroCopyPointerStability) {
    mio_mesh* m = build_tet_mesh();
    const void* pts1 = nullptr;
    mio_dtype dt;
    ASSERT_EQ(mio_mesh_get_points(m, &pts1, &dt), MIO_OK);
    // Non-mutating traffic must not invalidate or move the borrow.
    (void)mio_mesh_num_cell_blocks(m);
    (void)block_type(m, 0);
    const void* pts2 = nullptr;
    ASSERT_EQ(mio_mesh_get_points(m, &pts2, &dt), MIO_OK);
    EXPECT_EQ(pts1, pts2);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(pts1)[3], 1.1);
    mio_mesh_free(m);
}

TEST(CApi, ErrorPaths) {
    // NULL / invalid arguments.
    EXPECT_EQ(mio_read(nullptr, nullptr), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_write("out.vtu", nullptr, nullptr), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_num_points(nullptr), -1);
    EXPECT_EQ(mio_mesh_set_points(nullptr, MIO_FLOAT64, 1, 3, kPoints.data()), MIO_ERR_INVALID_ARG);

    mio_mesh* m = build_tet_mesh();

    // Wrong dtypes / shapes.
    EXPECT_EQ(mio_mesh_set_points(m, MIO_INT32, 5, 3, kPoints.data()), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 5, MIO_INT64, kConn.data()),
              MIO_ERR_INVALID_ARG);        // tetra has 4 nodes, not 5
    const std::int64_t bad_shape[] = {4};  // 5 points in the mesh
    const std::vector<double> four(4, 1.0);
    EXPECT_EQ(mio_mesh_add_point_data(m, "bad", MIO_FLOAT64, 1, bad_shape, four.data()),
              MIO_ERR_INVALID_ARG);

    // Out-of-range lookups.
    EXPECT_EQ(mio_mesh_cell_block_info(m, 7, nullptr, nullptr, nullptr), MIO_ERR_NOT_FOUND);
    const void* data = nullptr;
    mio_dtype dt;
    EXPECT_EQ(mio_mesh_get_point_data(m, "nope", &data, &dt, nullptr, nullptr), MIO_ERR_NOT_FOUND);
    EXPECT_EQ(mio_mesh_point_data_name(m, 0, nullptr, 0), -1);  // no point data yet

    // Unknown format / extension.
    EXPECT_EQ(mio_write("mesh.not_an_extension", m, nullptr), MIO_ERR_READ);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_write("mesh.vtu", m, "no_such_format"), MIO_ERR_NOT_FOUND);
    // openfoam gained a writer in v9.20.0, so the "resolvable format with no
    // writer" case it used to demonstrate needs a different format. glTF reads
    // since the general reader landed, so a file that does not exist is the failure.
    EXPECT_EQ(mio_read("no_such_mesh.glb", "gltf"), nullptr);

#ifndef MESHIOPLUSPLUS_HAS_HDF5
    // Compiled-out formats name the missing dependency.
    EXPECT_EQ(mio_write("mesh.med", m, nullptr), MIO_ERR_NOT_FOUND);
    EXPECT_NE(std::strstr(mio_last_error(), "HDF5"), nullptr) << mio_last_error();
#endif

    mio_mesh_free(m);
    mio_mesh_free(nullptr);  // NULL-safe
}

TEST(CApi, WritesGltfThroughTheRegistry) {
    mio_mesh* m = build_tet_mesh();
    const std::string glb = std::string(::testing::TempDir()) + "capi_mesh.glb";
    ASSERT_EQ(mio_write(glb.c_str(), m, nullptr), MIO_OK) << mio_last_error();
    std::FILE* f = std::fopen(glb.c_str(), "rb");
    ASSERT_NE(f, nullptr);
    char magic[4] = {0, 0, 0, 0};
    ASSERT_EQ(std::fread(magic, 1, 4, f), 4u);
    std::fclose(f);
    EXPECT_EQ(std::string(magic, 4), "glTF");
    // The suffix picks the container: a .gltf is JSON, with the .bin beside it.
    const std::string gltf = std::string(::testing::TempDir()) + "capi_mesh.gltf";
    ASSERT_EQ(mio_write(gltf.c_str(), m, "gltf"), MIO_OK) << mio_last_error();
    f = std::fopen(gltf.c_str(), "rb");
    ASSERT_NE(f, nullptr);
    char brace = 0;
    ASSERT_EQ(std::fread(&brace, 1, 1, f), 1u);
    std::fclose(f);
    EXPECT_EQ(brace, '{');
    std::remove(glb.c_str());
    std::remove(gltf.c_str());
    std::remove((std::string(::testing::TempDir()) + "capi_mesh.bin").c_str());
    mio_mesh_free(m);
}

TEST(CApi, StringBufferProtocol) {
    mio_mesh* m = build_tet_mesh();
    // Full length is returned even when the buffer is too small ("tetra" = 5).
    char tiny[3] = {'x', 'x', 'x'};
    EXPECT_EQ(mio_mesh_cell_block_type(m, 0, tiny, sizeof(tiny)), 5);
    EXPECT_STREQ(tiny, "te");                                  // truncated + NUL-terminated
    EXPECT_EQ(mio_mesh_cell_block_type(m, 0, nullptr, 0), 5);  // pure length query
    mio_mesh_free(m);
}

/* Ragged (polygon / polyhedron) connectivity.
 *
 * These used to need MED (the one C++ writer that serializes ragged blocks) to
 * construct a ragged mesh at all, and so only ran on an HDF5 build. The setters
 * mean the C API can now build one directly, which is why they are unguarded. */

TEST(CApi, PolygonBlockRoundTripsThroughTheCApi) {
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    const double xyz[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 2, 0.5, 0};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, xyz), MIO_OK);

    // A quad then a triangle -- the point of a jagged block.
    const std::int64_t row_offsets[] = {0, 4, 7};
    const std::int64_t nodes[] = {0, 1, 2, 3, 1, 4, 2};
    ASSERT_EQ(mio_mesh_add_polygon_block(m, "polygon", 2, row_offsets, nodes, 7), MIO_OK)
        << mio_last_error();

    mio_cell_block_info info{};
    ASSERT_EQ(mio_mesh_cell_block_info_ex(m, 0, &info), MIO_OK);
    EXPECT_EQ(info.num_cells, 2);
    EXPECT_EQ(info.nodes_per_cell, 0);
    EXPECT_EQ(info.is_ragged, 1);
    EXPECT_EQ(info.is_polyhedron, 0);
    EXPECT_EQ(info.num_faces, 2);  // 1-level: one face per cell
    EXPECT_EQ(info.num_nodes, 7);

    // The five-argument original must keep agreeing on the fields it shares.
    std::int64_t nc = 0, npc = -1;
    std::int32_t ragged = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(m, 0, &nc, &npc, &ragged), MIO_OK);
    EXPECT_EQ(nc, info.num_cells);
    EXPECT_EQ(npc, info.nodes_per_cell);
    EXPECT_EQ(ragged, info.is_ragged);

    // A ragged block has no rectangular buffer to borrow.
    const void* conn = nullptr;
    mio_dtype dt;
    EXPECT_EQ(mio_mesh_cell_block_conn(m, 0, &conn, &dt), MIO_ERR_UNSUPPORTED);

    mio_poly_conn* pc = mio_poly_conn_create(m, 0);
    ASSERT_NE(pc, nullptr) << mio_last_error();
    mio_poly_conn_shape shape{};
    ASSERT_EQ(mio_poly_conn_get_shape(pc, &shape), MIO_OK);
    EXPECT_EQ(shape.is_polyhedron, 0);
    EXPECT_EQ(shape.num_cells, 2);
    EXPECT_EQ(shape.num_faces, 2);
    EXPECT_EQ(shape.num_nodes, 7);

    std::int64_t n = -1;
    const std::int64_t* out_nodes = mio_poly_conn_nodes(pc, &n);
    ASSERT_NE(out_nodes, nullptr);
    EXPECT_EQ(n, 7);
    EXPECT_TRUE(std::equal(nodes, nodes + 7, out_nodes));
    const std::int64_t* out_rows = mio_poly_conn_face_offsets(pc, &n);
    ASSERT_NE(out_rows, nullptr);
    EXPECT_EQ(n, 3);
    EXPECT_TRUE(std::equal(row_offsets, row_offsets + 3, out_rows));

    // A 1-level block has no cell-offsets array: NULL, not a synthesized
    // identity that a caller could mistake for information.
    n = -1;
    EXPECT_EQ(mio_poly_conn_cell_offsets(pc, &n), nullptr);
    EXPECT_EQ(n, 0);

    mio_poly_conn_free(pc);
    mio_mesh_free(m);
}

TEST(CApi, PolyhedronBlockRoundTripsThroughTheCApi) {
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    const double xyz[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, xyz), MIO_OK);

    // Two cells: a 4-face tetrahedron, then a 3-face sliver (deliberately not
    // the same face count, so cell_offsets carries real information).
    const std::int64_t cell_offsets[] = {0, 4, 7};
    const std::int64_t face_offsets[] = {0, 3, 6, 9, 12, 15, 18, 21};
    const std::int64_t nodes[] = {0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0, 1, 2, 4, 2, 3, 4, 3, 1, 4};
    ASSERT_EQ(
        mio_mesh_add_polyhedron_block(m, "polyhedron", 2, cell_offsets, 7, face_offsets, nodes, 21),
        MIO_OK)
        << mio_last_error();

    mio_cell_block_info info{};
    ASSERT_EQ(mio_mesh_cell_block_info_ex(m, 0, &info), MIO_OK);
    EXPECT_EQ(info.num_cells, 2);
    EXPECT_EQ(info.is_ragged, 1);
    EXPECT_EQ(info.is_polyhedron, 1);
    EXPECT_EQ(info.num_faces, 7);
    EXPECT_EQ(info.num_nodes, 21);

    mio_poly_conn* pc = mio_poly_conn_create(m, 0);
    ASSERT_NE(pc, nullptr) << mio_last_error();
    mio_poly_conn_shape shape{};
    ASSERT_EQ(mio_poly_conn_get_shape(pc, &shape), MIO_OK);
    EXPECT_EQ(shape.is_polyhedron, 1);
    EXPECT_EQ(shape.num_cells, 2);
    EXPECT_EQ(shape.num_faces, 7);
    EXPECT_EQ(shape.num_nodes, 21);

    std::int64_t n = -1;
    const std::int64_t* out_nodes = mio_poly_conn_nodes(pc, &n);
    ASSERT_NE(out_nodes, nullptr);
    ASSERT_EQ(n, 21);
    EXPECT_TRUE(std::equal(nodes, nodes + 21, out_nodes));
    const std::int64_t* out_faces = mio_poly_conn_face_offsets(pc, &n);
    ASSERT_NE(out_faces, nullptr);
    ASSERT_EQ(n, 8);
    EXPECT_TRUE(std::equal(face_offsets, face_offsets + 8, out_faces));
    const std::int64_t* out_cells = mio_poly_conn_cell_offsets(pc, &n);
    ASSERT_NE(out_cells, nullptr);
    ASSERT_EQ(n, 3);
    EXPECT_TRUE(std::equal(cell_offsets, cell_offsets + 3, out_cells));

    // Face f of cell c spans nodes[face_offsets[cell_offsets[c] + f] .. +1) --
    // the documented indexing, checked rather than described.
    const std::int64_t f = out_cells[1] + 0;  // cell 1's first face
    EXPECT_EQ(out_faces[f + 1] - out_faces[f], 3);
    EXPECT_EQ(out_nodes[out_faces[f]], 1);

    mio_poly_conn_free(pc);
    mio_mesh_free(m);
}

TEST(CApi, PolyConnSnapshotSurvivesMutation) {
    // The contract that distinguishes the snapshot from every other getter on
    // this ABI: rule 3 borrows die at the next mutating call, this does not.
    mio_mesh* m = mio_mesh_create();
    const double xyz[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, xyz), MIO_OK);
    const std::int64_t row_offsets[] = {0, 4};
    const std::int64_t nodes[] = {0, 1, 2, 3};
    ASSERT_EQ(mio_mesh_add_polygon_block(m, "polygon", 1, row_offsets, nodes, 4), MIO_OK);

    mio_poly_conn* pc = mio_poly_conn_create(m, 0);
    ASSERT_NE(pc, nullptr) << mio_last_error();
    const std::int64_t* borrowed = mio_poly_conn_nodes(pc, nullptr);
    ASSERT_NE(borrowed, nullptr);

    const std::int64_t tri[] = {0, 1, 2};
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 1, 3, MIO_INT64, tri), MIO_OK);
    EXPECT_EQ(mio_mesh_num_cell_blocks(m), 2);
    // Still readable, and still the values it was created from.
    EXPECT_TRUE(std::equal(nodes, nodes + 4, borrowed));

    mio_poly_conn_free(pc);
    mio_mesh_free(m);
}

TEST(CApi, PolyConnRejectsRectangularBlocks) {
    mio_mesh* m = mio_mesh_create();
    const double xyz[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 3, 3, xyz), MIO_OK);
    const std::int64_t tri[] = {0, 1, 2};
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 1, 3, MIO_INT64, tri), MIO_OK);

    EXPECT_EQ(mio_poly_conn_create(m, 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_poly_conn_create(m, 7), nullptr);  // out of range

    mio_cell_block_info info{};
    ASSERT_EQ(mio_mesh_cell_block_info_ex(m, 0, &info), MIO_OK);
    EXPECT_EQ(info.is_ragged, 0);
    EXPECT_EQ(info.is_polyhedron, 0);
    EXPECT_EQ(info.num_faces, 1);
    EXPECT_EQ(info.num_nodes, 3);
    mio_mesh_free(m);
}

TEST(CApi, RaggedSettersRejectMalformedOffsets) {
    mio_mesh* m = mio_mesh_create();
    const std::int64_t nodes[] = {0, 1, 2, 3};

    const std::int64_t not_zero_based[] = {1, 4};
    EXPECT_EQ(mio_mesh_add_polygon_block(m, "polygon", 1, not_zero_based, nodes, 4),
              MIO_ERR_INVALID_ARG);
    const std::int64_t decreasing[] = {0, 3, 2};
    EXPECT_EQ(mio_mesh_add_polygon_block(m, "polygon", 2, decreasing, nodes, 2),
              MIO_ERR_INVALID_ARG);
    const std::int64_t wrong_total[] = {0, 3};
    EXPECT_EQ(mio_mesh_add_polygon_block(m, "polygon", 1, wrong_total, nodes, 4),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_add_polygon_block(m, "polygon", 1, nullptr, nodes, 4), MIO_ERR_INVALID_ARG);

    const std::int64_t negative[] = {0, -1, 2, 3};
    const std::int64_t ok_offsets[] = {0, 4};
    EXPECT_EQ(mio_mesh_add_polygon_block(m, "polygon", 1, ok_offsets, negative, 4),
              MIO_ERR_INVALID_ARG);

    // The 2-level setter validates both offset arrays independently.
    const std::int64_t cell_offsets[] = {0, 2};
    const std::int64_t bad_faces[] = {0, 3, 5};  // must end at num_nodes == 6
    EXPECT_EQ(
        mio_mesh_add_polyhedron_block(m, "polyhedron", 1, cell_offsets, 2, bad_faces, nodes, 6),
        MIO_ERR_INVALID_ARG);
    const std::int64_t bad_cells[] = {0, 3};  // must end at num_faces == 2
    const std::int64_t ok_faces[] = {0, 3, 6};
    EXPECT_EQ(mio_mesh_add_polyhedron_block(m, "polyhedron", 1, bad_cells, 2, ok_faces, nodes, 6),
              MIO_ERR_INVALID_ARG);

    // Nothing was added by any of the rejected calls.
    EXPECT_EQ(mio_mesh_num_cell_blocks(m), 0);
    mio_mesh_free(m);
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(CApi, RaggedBlockBuiltThroughTheCApiSurvivesAMedRoundTrip) {
    // The end-to-end path the setters exist for: build ragged through C, write
    // it with the one C++ writer that serializes ragged blocks (MED's POG),
    // read it back and compare through the snapshot.
    mio_mesh* m = mio_mesh_create();
    const double xyz[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 2, 0.5, 0};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, xyz), MIO_OK);
    const std::int64_t row_offsets[] = {0, 4, 7};
    const std::int64_t nodes[] = {0, 1, 2, 3, 1, 4, 2};
    ASSERT_EQ(mio_mesh_add_polygon_block(m, "polygon", 2, row_offsets, nodes, 7), MIO_OK);

    const std::string med = mt::temp_path("_capi_ragged.med");
    ASSERT_EQ(mio_write(med.c_str(), m, nullptr), MIO_OK) << mio_last_error();
    mio_mesh_free(m);

    mio_mesh* back = mio_read(med.c_str(), nullptr);
    ASSERT_NE(back, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_num_cell_blocks(back), 1);
    mio_poly_conn* pc = mio_poly_conn_create(back, 0);
    ASSERT_NE(pc, nullptr) << mio_last_error();
    std::int64_t n = -1;
    const std::int64_t* out_nodes = mio_poly_conn_nodes(pc, &n);
    ASSERT_EQ(n, 7);
    EXPECT_TRUE(std::equal(nodes, nodes + 7, out_nodes));
    const std::int64_t* out_rows = mio_poly_conn_face_offsets(pc, &n);
    ASSERT_EQ(n, 3);
    EXPECT_TRUE(std::equal(row_offsets, row_offsets + 3, out_rows));
    mio_poly_conn_free(pc);
    mio_mesh_free(back);
    std::remove(med.c_str());
}
#endif

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(CApi, HdfFormatConvert) {
    const std::string vtu = mt::temp_path("_capi_h5.vtu");
    const std::string med = mt::temp_path("_capi_h5.med");
    mio_mesh* m = build_tet_mesh();
    ASSERT_EQ(mio_write(vtu.c_str(), m, nullptr), MIO_OK);
    mio_mesh_free(m);
    ASSERT_EQ(mio_convert(vtu.c_str(), nullptr, med.c_str(), nullptr), MIO_OK) << mio_last_error();
    mio_mesh* r = mio_read(med.c_str(), nullptr);
    ASSERT_NE(r, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(r), 5);
    mio_mesh_free(r);
    std::remove(vtu.c_str());
    std::remove(med.c_str());
}
#endif

/* --- data operations ---------------------------------------------------- */

// A tet mesh carrying point, cell and field data, for the data operations.
mio_mesh* build_data_mesh() {
    mio_mesh* m = build_tet_mesh();  // 5 points, 2 tetra
    static const std::array<double, 5> temperature = {0.0, 1.0, 2.0, 3.0, 4.0};
    static const std::array<double, 15> velocity = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 2, 0, 0};
    static const std::array<double, 2> mat = {10.0, 20.0};
    static const std::array<double, 3> meta = {1.0, 2.0, 3.0};
    std::int64_t s1[1] = {5};
    std::int64_t s2[2] = {5, 3};
    std::int64_t sc[1] = {2};
    std::int64_t sf[1] = {3};
    EXPECT_EQ(mio_mesh_add_point_data(m, "T", MIO_FLOAT64, 1, s1, temperature.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_add_point_data(m, "v", MIO_FLOAT64, 2, s2, velocity.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_append_cell_data(m, "mat", MIO_FLOAT64, 1, sc, mat.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_add_field_data(m, "meta", MIO_FLOAT64, 1, sf, meta.data()), MIO_OK);
    return m;
}

TEST(CApi, DataDropAndKeep) {
    mio_mesh* m = build_data_mesh();
    const char* names[] = {"T"};
    mio_mesh* dropped = mio_data_drop(m, MIO_DATA_POINT, names, 1, 0);
    ASSERT_NE(dropped, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_point_data(dropped), 1);  // only "v" left
    // Geometry is never modified.
    EXPECT_EQ(mio_mesh_num_points(dropped), 5);
    EXPECT_EQ(mio_mesh_num_cell_blocks(dropped), 1);
    mio_mesh_free(dropped);

    mio_mesh* kept = mio_data_keep(m, MIO_DATA_POINT, names, 1, 0);
    ASSERT_NE(kept, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_point_data(kept), 1);
    mio_mesh_free(kept);
    mio_mesh_free(m);
}

TEST(CApi, SetsDataConversionsAndErrorGuards) {
    mio_mesh* mesh = build_data_mesh();
    const std::int64_t a[] = {0, 1}, b[] = {1, 2};
    ASSERT_EQ(mio_mesh_add_region(mesh, "a", MIO_REGION_POINT, 0, 7, a, 2), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(mesh, "b", MIO_REGION_POINT, -1, -1, b, 2), MIO_OK);
    const char* order[] = {"b", "a"};
    mio_mesh* out = mio_sets_to_data(mesh, MIO_DATA_POINT, nullptr, "-", order, 2);
    ASSERT_NE(out, nullptr) << mio_last_error();
    auto* regions = mio_regions_create(out);
    ASSERT_NE(regions, nullptr);
    EXPECT_EQ(mio_regions_count(regions), 0);
    mio_regions_free(regions);
    const void* data = nullptr;
    mio_dtype dtype{};
    ASSERT_EQ(mio_mesh_get_point_data(out, "b-a", &data, &dtype, nullptr, nullptr), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    const auto* labels = static_cast<const std::int64_t*>(data);
    EXPECT_EQ(labels[0], 1);
    EXPECT_EQ(labels[1], 1);
    EXPECT_EQ(labels[2], 0);
    EXPECT_EQ(labels[4], -1);
    mio_mesh* back = mio_data_to_sets(out, MIO_DATA_POINT, "b-a");
    ASSERT_NE(back, nullptr) << mio_last_error();
    regions = mio_regions_create(back);
    ASSERT_NE(regions, nullptr);
    EXPECT_EQ(mio_regions_count(regions), 3);  // -1 is a tag, not discarded.
    mio_regions_free(regions);
    mio_mesh_free(back);
    mio_mesh_free(out);
    EXPECT_EQ(mio_sets_to_data(nullptr, MIO_DATA_POINT, nullptr, nullptr, nullptr, 0), nullptr);
    EXPECT_EQ(mio_sets_to_data(mesh, MIO_DATA_FIELD, nullptr, nullptr, nullptr, 0), nullptr);
    EXPECT_EQ(mio_sets_to_data(mesh, MIO_DATA_POINT, nullptr, nullptr, order, -1), nullptr);
    EXPECT_EQ(mio_data_to_sets(mesh, MIO_DATA_POINT, nullptr), nullptr);
    EXPECT_EQ(mio_data_to_sets(mesh, MIO_DATA_POINT, "T"), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    mio_mesh_free(mesh);
}

TEST(CApi, DataDropUnknownKeyFails) {
    mio_mesh* m = build_data_mesh();
    const char* names[] = {"nope"};
    EXPECT_EQ(mio_data_drop(m, MIO_DATA_POINT, names, 1, 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    // ignore_missing makes it a no-op instead.
    mio_mesh* ok = mio_data_drop(m, MIO_DATA_POINT, names, 1, 1);
    ASSERT_NE(ok, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_point_data(ok), 2);
    mio_mesh_free(ok);
    mio_mesh_free(m);
}

TEST(CApi, DataRename) {
    mio_mesh* m = build_data_mesh();
    mio_mesh* out = mio_data_rename(m, MIO_DATA_POINT, "T", "temperature");
    ASSERT_NE(out, nullptr) << mio_last_error();
    char buf[64] = {};
    bool found = false;
    for (std::int64_t i = 0; i < mio_mesh_num_point_data(out); ++i) {
        mio_mesh_point_data_name(out, i, buf, sizeof(buf));
        if (std::string(buf) == "temperature")
            found = true;
    }
    EXPECT_TRUE(found);
    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, DataAveraging) {
    mio_mesh* m = build_data_mesh();
    mio_mesh* to_cell = mio_data_point_to_cell(m, nullptr, 0, nullptr);
    ASSERT_NE(to_cell, nullptr) << mio_last_error();
    EXPECT_GT(mio_mesh_num_cell_data(to_cell), 0);
    mio_mesh_free(to_cell);

    const char* names[] = {"mat"};
    mio_mesh* to_point = mio_data_cell_to_point(m, names, 1, MIO_WEIGHT_UNIFORM, nullptr);
    ASSERT_NE(to_point, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(to_point, "mat", &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(shape[0], 5);
    // A mean is not an integer: the output is always Float64.
    EXPECT_EQ(dt, MIO_FLOAT64);
    mio_mesh_free(to_point);

    mio_mesh* weighted = mio_data_cell_to_point(m, names, 1, MIO_WEIGHT_MEASURE, "_w");
    ASSERT_NE(weighted, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_get_point_data(weighted, "mat_w", &data, &dt, &ndim, shape), MIO_OK);
    mio_mesh_free(weighted);
    mio_mesh_free(m);
}

TEST(CApi, DataCalc) {
    mio_mesh* m = build_data_mesh();
    mio_mesh* out = mio_data_calc(m, "norm(v)", MIO_DATA_POINT, "speed", 0);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(out, "speed", &data, &dt, &ndim, shape), MIO_OK);
    ASSERT_EQ(dt, MIO_FLOAT64);
    const double* speed = static_cast<const double*>(data);
    EXPECT_NEAR(speed[0], 1.0, 1e-12);
    EXPECT_NEAR(speed[3], std::sqrt(2.0), 1e-12);
    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, DataCalcRejectsBadExpressions) {
    mio_mesh* m = build_data_mesh();
    EXPECT_EQ(mio_data_calc(m, "log(T)", MIO_DATA_POINT, "o", 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_data_calc(m, "nope + 1", MIO_DATA_POINT, "o", 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    EXPECT_EQ(mio_data_calc(m, "T +", MIO_DATA_POINT, "o", 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    mio_mesh_free(m);
}

TEST(CApi, DataCondition) {
    mio_mesh* m = build_data_mesh();
    const char* names[] = {"T"};
    mio_mesh* out = mio_data_condition(m, MIO_DATA_POINT, names, 1, MIO_COND_NORMALIZE, 0.0, 1.0,
                                       MIO_SCOPE_COMPONENT, MIO_NAN_IGNORE, 0.0, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(out, "T", &data, &dt, &ndim, shape), MIO_OK);
    const double* t = static_cast<const double*>(data);
    // T = {0,1,2,3,4} -> min maps to 0, max to 1.
    EXPECT_NEAR(t[0], 0.0, 1e-12);
    EXPECT_NEAR(t[4], 1.0, 1e-12);
    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, TensorInvariants) {
    mio_mesh* m = build_data_mesh();
    // xx yy zz xy yz zx, repeated once per point (5 points).
    static const std::array<double, 30> s = {1,   2,   3,   0.5, 0.6, 0.7, 1,   2,   3,   0.5,
                                             0.6, 0.7, 1,   2,   3,   0.5, 0.6, 0.7, 1,   2,
                                             3,   0.5, 0.6, 0.7, 1,   2,   3,   0.5, 0.6, 0.7};
    std::int64_t shape[2] = {5, 6};
    ASSERT_EQ(mio_mesh_add_point_data(m, "s", MIO_FLOAT64, 2, shape, s.data()), MIO_OK);

    const char* names[] = {"s"};
    mio_mesh* out = mio_tensor_invariants(
        m, MIO_DATA_POINT, names, 1, MIO_TINV_MISES | MIO_TINV_HYDROSTATIC, nullptr, nullptr, 1);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dt;
    std::int32_t ndim = 0;
    std::int64_t out_shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(out, "s_mises", &data, &dt, &ndim, out_shape), MIO_OK);
    const double* mises = static_cast<const double*>(data);
    const double expect =
        std::sqrt(0.5 * ((1 - 2) * (1 - 2) + (2 - 3) * (2 - 3) + (3 - 1) * (3 - 1) +
                         6.0 * (0.5 * 0.5 + 0.6 * 0.6 + 0.7 * 0.7)));
    EXPECT_NEAR(mises[0], expect, 1e-12);
    ASSERT_EQ(mio_mesh_get_point_data(out, "s_hydrostatic", &data, &dt, &ndim, out_shape), MIO_OK);
    EXPECT_NEAR(static_cast<const double*>(data)[0], 2.0, 1e-12);
    // Only the two requested outputs were computed.
    EXPECT_NE(mio_mesh_get_point_data(out, "s_principal", &data, &dt, &ndim, out_shape), MIO_OK);
    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, TensorInvariantsRejectsFieldLocation) {
    mio_mesh* m = build_data_mesh();
    EXPECT_EQ(mio_tensor_invariants(m, MIO_DATA_FIELD, nullptr, 0, 0, nullptr, nullptr, 1),
              nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    mio_mesh_free(m);
}

TEST(CApi, DataInfoHandle) {
    mio_mesh* m = build_data_mesh();
    mio_data_info* info = mio_data_info_create(m);
    ASSERT_NE(info, nullptr) << mio_last_error();
    const std::int64_t n = mio_data_info_count(info);
    EXPECT_EQ(n, 4);  // T, v, mat, meta

    bool saw_t = false;
    for (std::int64_t i = 0; i < n; ++i) {
        char buf[64] = {};
        const std::int64_t len = mio_data_info_name(info, i, buf, sizeof(buf));
        EXPECT_GE(len, 0);
        mio_data_array_info entry;
        ASSERT_EQ(mio_data_info_entry(info, i, &entry), MIO_OK);
        if (std::string(buf) == "T" && entry.location == MIO_DATA_POINT) {
            saw_t = true;
            EXPECT_EQ(entry.num_entries, 5);
            EXPECT_EQ(entry.num_components, 1);
            EXPECT_EQ(entry.num_finite, 5);
            EXPECT_EQ(entry.num_nan, 0);
            EXPECT_NEAR(entry.min, 0.0, 1e-12);
            EXPECT_NEAR(entry.max, 4.0, 1e-12);
            EXPECT_NEAR(entry.mean, 2.0, 1e-12);
            double cmin = 0, cmax = 0, cmean = 0;
            ASSERT_EQ(mio_data_info_component(info, i, 0, &cmin, &cmax, &cmean), MIO_OK);
            EXPECT_NEAR(cmax, 4.0, 1e-12);
            // Every out pointer is optional.
            EXPECT_EQ(mio_data_info_component(info, i, 0, nullptr, nullptr, nullptr), MIO_OK);
        }
    }
    EXPECT_TRUE(saw_t);

    // Out-of-range indices are rejected, not dereferenced.
    mio_data_array_info entry;
    EXPECT_NE(mio_data_info_entry(info, n, &entry), MIO_OK);
    EXPECT_EQ(mio_data_info_name(info, -1, nullptr, 0), -1);
    EXPECT_NE(mio_data_info_component(info, 0, 999, nullptr, nullptr, nullptr), MIO_OK);

    mio_data_info_free(info);
    mio_data_info_free(nullptr);  // must tolerate NULL
    mio_mesh_free(m);
}

TEST(CApi, DataInfoNameBufferTooSmall) {
    mio_mesh* m = build_data_mesh();
    mio_data_info* info = mio_data_info_create(m);
    ASSERT_NE(info, nullptr);
    // The required length is returned even when the buffer cannot hold it.
    char tiny[2] = {};
    const std::int64_t needed = mio_data_info_name(info, 0, tiny, sizeof(tiny));
    EXPECT_GE(needed, 0);
    EXPECT_EQ(tiny[1], '\0');  // always NUL-terminated
    mio_data_info_free(info);
    mio_mesh_free(m);
}

TEST(CApi, DataNullArgumentsAreRejected) {
    mio_mesh* m = build_data_mesh();
    const char* names[] = {"T"};
    EXPECT_EQ(mio_data_drop(nullptr, MIO_DATA_POINT, names, 1, 0), nullptr);
    EXPECT_EQ(mio_data_rename(m, MIO_DATA_POINT, nullptr, "x"), nullptr);
    EXPECT_EQ(mio_data_calc(m, nullptr, MIO_DATA_POINT, "o", 0), nullptr);
    EXPECT_EQ(mio_data_info_create(nullptr), nullptr);
    EXPECT_EQ(mio_data_info_count(nullptr), -1);
    // A NULL names array with a positive count must be caught, not dereferenced.
    EXPECT_EQ(mio_data_drop(m, MIO_DATA_POINT, nullptr, 3, 0), nullptr);
    EXPECT_STRNE(mio_last_error(), "");
    mio_mesh_free(m);
}

TEST(CApi, DataIntegrateHandle) {
    mio_mesh* m = build_data_mesh();  // cell_data "mat" = {10.0, 20.0} over 2 tets
    const int64_t region_cells[1] = {0};
    ASSERT_EQ(mio_mesh_add_region(m, "one_cell", MIO_REGION_CELL, -1, -1, region_cells, 1), MIO_OK);

    const char* names[] = {"mat"};
    mio_data_integrate* result = mio_data_integrate_create(m, names, 1);
    ASSERT_NE(result, nullptr) << mio_last_error();
    EXPECT_EQ(mio_data_integrate_count(result), 1);

    char buf[64] = {};
    const std::int64_t len = mio_data_integrate_name(result, 0, buf, sizeof(buf));
    EXPECT_GE(len, 0);
    EXPECT_STREQ(buf, "mat");

    mio_field_integral_info entry;
    ASSERT_EQ(mio_data_integrate_entry(result, 0, &entry), MIO_OK);
    EXPECT_EQ(entry.num_components, 1);
    EXPECT_EQ(entry.num_cells, 2);
    EXPECT_EQ(entry.num_skipped, 0);

    double total = 0, mean = 0, domain_measure = 0;
    std::int64_t num_nan = 0;
    ASSERT_EQ(mio_data_integrate_component(result, 0, 0, &total, &mean, &domain_measure, &num_nan),
              MIO_OK);
    EXPECT_GT(total, 0.0);
    EXPECT_GT(domain_measure, 0.0);
    EXPECT_NEAR(mean, total / domain_measure, 1e-9);
    EXPECT_EQ(num_nan, 0);
    // Every out pointer is optional.
    EXPECT_EQ(mio_data_integrate_component(result, 0, 0, nullptr, nullptr, nullptr, nullptr),
              MIO_OK);

    ASSERT_EQ(mio_data_integrate_region_count(result, 0), 1);
    char rbuf[64] = {};
    EXPECT_GE(mio_data_integrate_region_name(result, 0, 0, rbuf, sizeof(rbuf)), 0);
    EXPECT_STREQ(rbuf, "one_cell");
    mio_field_integral_info rentry;
    ASSERT_EQ(mio_data_integrate_region_entry(result, 0, 0, &rentry), MIO_OK);
    EXPECT_EQ(rentry.num_cells, 1);
    double rtotal = 0;
    ASSERT_EQ(
        mio_data_integrate_region_component(result, 0, 0, 0, &rtotal, nullptr, nullptr, nullptr),
        MIO_OK);
    EXPECT_LT(rtotal, total);  // one cell contributes less than both

    // Out-of-range indices are rejected, not dereferenced.
    EXPECT_NE(mio_data_integrate_entry(result, 1, &entry), MIO_OK);
    EXPECT_EQ(mio_data_integrate_name(result, -1, nullptr, 0), -1);
    EXPECT_NE(mio_data_integrate_component(result, 0, 999, nullptr, nullptr, nullptr, nullptr),
              MIO_OK);
    EXPECT_NE(mio_data_integrate_region_entry(result, 0, 1, &rentry), MIO_OK);

    mio_data_integrate_free(result);
    mio_data_integrate_free(nullptr);  // must tolerate NULL

    // A point_data-only name fails, naming the fix.
    const char* point_only[] = {"T"};
    EXPECT_EQ(mio_data_integrate_create(m, point_only, 1), nullptr);
    EXPECT_STRNE(mio_last_error(), "");

    EXPECT_EQ(mio_data_integrate_create(nullptr, nullptr, 0), nullptr);
    EXPECT_EQ(mio_data_integrate_count(nullptr), -1);

    mio_mesh_free(m);
}

}  // namespace

TEST(CApi, GmshPeriodicSideChannel) {
    for (int version : {22, 40, 41})
        for (bool binary : {false, true}) {
            const auto path = mt::temp_path(".msh");
            const auto bytes = mt::gmsh_periodic_fixture(version, binary);
            {
                std::ofstream out(path, std::ios::binary);
                out.write(bytes.data(), bytes.size());
            }
            EXPECT_EQ(mio_read(path.c_str(), "gmsh"), nullptr);
            mio_format_info* info = nullptr;
            mio_mesh* mesh = mio_read_with_info(path.c_str(), "gmsh", nullptr, &info);
            ASSERT_NE(mesh, nullptr) << mio_last_error();
            ASSERT_NE(info, nullptr);
            EXPECT_EQ(mio_gmsh_info_count(info, MIO_GMSH_PERIODIC), 2);
            EXPECT_EQ(mio_gmsh_info_count(info, 123), -1);
            const void* data = nullptr;
            mio_dtype dtype;
            int32_t ndim;
            int64_t shape[MIO_MAX_NDIM]{};
            ASSERT_EQ(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, 1, 2, &data, &dtype, &ndim, shape),
                MIO_OK);
            EXPECT_EQ(dtype, MIO_INT64);
            EXPECT_EQ(ndim, 2);
            EXPECT_EQ(shape[0], 3);
            EXPECT_EQ(shape[1], 2);
            EXPECT_EQ(static_cast<const int64_t*>(data)[1], 1);
            EXPECT_EQ(static_cast<const int64_t*>(data)[3], 3);
            EXPECT_EQ(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, -1, 2, &data, &dtype, &ndim, shape),
                MIO_ERR_INVALID_ARG);
            EXPECT_NE(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, 0, 99, &data, &dtype, &ndim, shape),
                MIO_OK);
            EXPECT_EQ(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, 0, 2, nullptr, &dtype, &ndim, shape),
                MIO_OK);
            ASSERT_EQ(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, 0, 1, &data, &dtype, &ndim, shape),
                MIO_OK);
            EXPECT_EQ(shape[0], 0);
            for (const char* target : {"gmsh", "gmsh22"}) {
                const auto output = mt::temp_path(".msh");
                ASSERT_EQ(mio_write_with_info(output.c_str(), mesh, target, info), MIO_OK)
                    << mio_last_error();
                mio_format_info* info2 = nullptr;
                auto* back = mio_read_with_info(output.c_str(), "gmsh", nullptr, &info2);
                ASSERT_NE(back, nullptr) << mio_last_error();
                ASSERT_EQ(mio_gmsh_info_array(info2, MIO_GMSH_PERIODIC, 1, 2, &data, &dtype, &ndim,
                                              shape),
                          MIO_OK);
                EXPECT_EQ(static_cast<const int64_t*>(data)[3], 3);
                mio_mesh_free(back);
                mio_format_info_free(info2);
                std::remove(output.c_str());
            }
            mio_mesh_free(mesh);  // Info arrays outlive the mesh handle.
            ASSERT_EQ(
                mio_gmsh_info_array(info, MIO_GMSH_PERIODIC, 1, 0, &data, &dtype, &ndim, shape),
                MIO_OK);
            EXPECT_EQ(static_cast<const int32_t*>(data)[1], 12);
            mio_format_info_free(info);
            std::remove(path.c_str());
        }
    EXPECT_EQ(mio_gmsh_info_count(nullptr, MIO_GMSH_PERIODIC), -1);
}

TEST(CApi, NurbsSideChannel) {
    // A quarter annulus, one NURBS patch of degrees (2, 1) (MFEM global form).
    const std::string w = "0.70710678118654757";
    const auto path = mt::temp_path(".mesh");
    {
        std::ofstream out(path, std::ios::binary);
        out << "MFEM NURBS mesh v1.0\ndimension\n2\nelements\n1\n5 3 0 1 2 3\nboundary\n0\n"
               "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\nvertices\n4\nknotvectors\n2\n"
               "2 3 0 0 0 1 1 1\n1 2 0 0 1 1\nweights\n1\n1\n1\n1\n"
            << w << "\n"
            << w
            << "\nFiniteElementSpace\nFiniteElementCollection: NURBS2\nVDim: 2\nOrdering: 1\n"
               "1 0\n0 1\n0 2\n2 0\n1 1\n2 2\n";
    }
    mio_format_info* info = nullptr;
    mio_mesh* mesh = mio_read_with_info(path.c_str(), "mfem", nullptr, &info);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(mio_nurbs_info_count(info, MIO_NURBS_PATCHES), 1);
    EXPECT_EQ(mio_nurbs_info_count(info, MIO_NURBS_FIELDS), 0);
    EXPECT_EQ(mio_nurbs_info_count(info, MIO_NURBS_TOPOLOGY), 1);
    EXPECT_EQ(mio_nurbs_info_count(info, MIO_NURBS_CELL_MAP), 1);
    EXPECT_EQ(mio_nurbs_info_count(info, 99), -1);
    EXPECT_EQ(mio_gmsh_info_count(info, MIO_GMSH_PERIODIC), -1);  // not a gmsh handle
    int64_t v = 0;
    ASSERT_EQ(mio_nurbs_info_int(info, MIO_NURBS_PATCHES, 0, 0, &v), MIO_OK);
    EXPECT_EQ(v, 2);
    ASSERT_EQ(mio_nurbs_info_int(info, MIO_NURBS_PATCHES, 0, 1, &v), MIO_OK);
    EXPECT_EQ(v, 2);
    ASSERT_EQ(mio_nurbs_info_int(info, MIO_NURBS_TOPOLOGY, 0, 1, &v), MIO_OK);
    EXPECT_EQ(v, 4);
    EXPECT_EQ(mio_nurbs_info_int(info, MIO_NURBS_PATCHES, 1, 0, &v), MIO_ERR_INVALID_ARG);
    const void* data = nullptr;
    mio_dtype dtype;
    int32_t ndim;
    int64_t shape[MIO_MAX_NDIM]{};
    ASSERT_EQ(mio_nurbs_info_array(info, MIO_NURBS_PATCHES, 0, 0, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(dtype, MIO_INT32);
    EXPECT_EQ(static_cast<const int32_t*>(data)[0], 2);
    EXPECT_EQ(static_cast<const int32_t*>(data)[1], 1);
    ASSERT_EQ(mio_nurbs_info_array(info, MIO_NURBS_PATCHES, 0, 1, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 6);
    EXPECT_EQ(shape[1], 2);
    ASSERT_EQ(mio_nurbs_info_array(info, MIO_NURBS_PATCHES, 0, 2, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[1], 0.70710678118654757);
    ASSERT_EQ(mio_nurbs_info_array(info, MIO_NURBS_PATCHES, 0, 4, &data, &dtype, &ndim, shape),
              MIO_OK);  // the knot vector of direction 1
    EXPECT_EQ(shape[0], 4);
    EXPECT_NE(mio_nurbs_info_array(info, MIO_NURBS_PATCHES, 0, 5, &data, &dtype, &ndim, shape),
              MIO_OK);  // no direction 2
    ASSERT_EQ(mio_nurbs_info_array(info, MIO_NURBS_CELL_MAP, 0, 1, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(shape[0], 1);
    EXPECT_EQ(shape[1], 3);
    EXPECT_EQ(mio_nurbs_info_string(info, MIO_NURBS_FIELDS, 0, 0, nullptr, 0), -1);

    // Written back as the NURBS mesh, and read back the same.
    const auto output = mt::temp_path(".mesh");
    ASSERT_EQ(mio_write_with_info(output.c_str(), mesh, "mfem", info), MIO_OK) << mio_last_error();
    {
        std::ifstream in(output);
        std::string first;
        std::getline(in, first);
        EXPECT_EQ(first, "MFEM NURBS mesh v1.0");
    }
    mio_format_info* info2 = nullptr;
    mio_mesh* back = mio_read_with_info(output.c_str(), "mfem", nullptr, &info2);
    ASSERT_NE(back, nullptr) << mio_last_error();
    ASSERT_EQ(mio_nurbs_info_array(info2, MIO_NURBS_PATCHES, 0, 1, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(static_cast<const double*>(data)[2], 1.0);  // control point (1, 1)
    EXPECT_EQ(mio_write_with_info(output.c_str(), mesh, "vtu", info), MIO_ERR_INVALID_ARG);
    mio_mesh_free(back);
    mio_format_info_free(info2);
    mio_mesh_free(mesh);
    mio_format_info_free(info);
    std::remove(path.c_str());
    std::remove(output.c_str());
}

TEST(CApi, IgaRefinesAndEvaluatesTheNurbsModel) {
    const std::string w = "0.70710678118654757";
    const auto path = mt::temp_path(".mesh");
    {
        std::ofstream out(path, std::ios::binary);
        out << "MFEM NURBS mesh v1.0\ndimension\n2\nelements\n1\n5 3 0 1 2 3\nboundary\n0\n"
               "edges\n4\n0 0 1\n0 3 2\n1 0 3\n1 1 2\nvertices\n4\nknotvectors\n2\n"
               "2 3 0 0 0 1 1 1\n1 2 0 0 1 1\nweights\n1\n1\n1\n1\n"
            << w << "\n"
            << w
            << "\nFiniteElementSpace\nFiniteElementCollection: NURBS2\nVDim: 2\nOrdering: 1\n"
               "1 0\n0 1\n0 2\n2 0\n1 1\n2 2\n";
    }
    mio_format_info* info = nullptr;
    mio_mesh* mesh = mio_read_with_info(path.c_str(), "mfem", nullptr, &info);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    mio_mesh_free(mesh);

    mio_iga_opts opts;
    mio_iga_opts_init(&opts);
    EXPECT_EQ(opts.subdivisions[1], 1);
    EXPECT_EQ(opts.trim_samples, 16);
    EXPECT_EQ(opts.reserved[0], 0);
    for (int d = 0; d < 3; ++d) {
        opts.subdivisions[d] = 3;
        opts.elevate[d] = 1;
    }
    mio_format_info* refined = nullptr;
    mio_mesh* fine = mio_iga_refine(info, &opts, &refined);
    ASSERT_NE(fine, nullptr) << mio_last_error();
    ASSERT_NE(refined, nullptr);
    char fmt[16];
    EXPECT_EQ(mio_format_info_format(refined, fmt, sizeof fmt), 4);
    EXPECT_STREQ(fmt, "mfem");
    const void* data = nullptr;
    mio_dtype dtype;
    int32_t ndim;
    int64_t shape[MIO_MAX_NDIM]{};
    ASSERT_EQ(mio_nurbs_info_array(refined, MIO_NURBS_PATCHES, 0, 0, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(static_cast<const int32_t*>(data)[0], 3);  // degrees 2 and 1 raised by 1
    EXPECT_EQ(static_cast<const int32_t*>(data)[1], 2);
    // The refined model writes back as an MFEM NURBS mesh through the returned mesh.
    const auto output = mt::temp_path(".mesh");
    ASSERT_EQ(mio_write_with_info(output.c_str(), fine, "mfem", refined), MIO_OK)
        << mio_last_error();
    mio_mesh_free(fine);

    // Evaluate the same point on the original and the refined patch.
    const double uv[4] = {0.3, 0.7, 0.9, 0.2};
    mio_iga_evaluation* a = mio_iga_evaluate(info, 0, uv, 2, 1, 16);
    mio_iga_evaluation* b = mio_iga_evaluate(refined, 0, uv, 2, 1, 16);
    ASSERT_NE(a, nullptr) << mio_last_error();
    ASSERT_NE(b, nullptr) << mio_last_error();
    EXPECT_EQ(mio_iga_evaluation_num_fields(a), 0);
    const void* pa = nullptr;
    const void* pb = nullptr;
    ASSERT_EQ(mio_iga_evaluation_array(a, MIO_IGA_POINTS, 0, &pa, &dtype, &ndim, shape), MIO_OK);
    EXPECT_EQ(dtype, MIO_FLOAT64);
    EXPECT_EQ(shape[0], 2);
    EXPECT_EQ(shape[1], 2);
    ASSERT_EQ(mio_iga_evaluation_array(b, MIO_IGA_POINTS, 0, &pb, nullptr, nullptr, nullptr),
              MIO_OK);
    for (int i = 0; i < 4; ++i)
        EXPECT_NEAR(static_cast<const double*>(pa)[i], static_cast<const double*>(pb)[i], 1e-13);
    ASSERT_EQ(mio_iga_evaluation_array(a, MIO_IGA_DERIVATIVES, 0, &pa, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(ndim, 3);
    EXPECT_EQ(shape[1], 2);
    // a planar patch has no 3-D normals, no trimming, no fields
    EXPECT_EQ(mio_iga_evaluation_array(a, MIO_IGA_NORMALS, 0, &pa, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_iga_evaluation_array(a, MIO_IGA_INSIDE, 0, &pa, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_iga_evaluation_array(a, MIO_IGA_FIELD, 0, &pa, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_iga_evaluation_free(a);
    mio_iga_evaluation_free(b);

    // Errors: out-of-range parameters and patch, a refinement of unequal directions.
    const double bad[2] = {2.0, 0.5};
    EXPECT_EQ(mio_iga_evaluate(info, 0, bad, 1, 0, 16), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("knot range"), std::string::npos);
    EXPECT_EQ(mio_iga_evaluate(info, 3, uv, 1, 0, 16), nullptr);
    opts.subdivisions[0] = 2;
    EXPECT_EQ(mio_iga_refine(info, &opts, nullptr), nullptr);
    EXPECT_EQ(mio_iga_refine(nullptr, nullptr, nullptr), nullptr);
    mio_iga_evaluation_free(nullptr);
    EXPECT_EQ(mio_iga_evaluation_num_fields(nullptr), -1);

    // The sampling alone: one cell per span, the map in the new handle.
    mio_iga_opts_init(&opts);
    mio_format_info* mapped = nullptr;
    mio_mesh* cells = mio_iga_sample(info, &opts, &mapped);
    ASSERT_NE(cells, nullptr) << mio_last_error();
    ASSERT_NE(mapped, nullptr);
    EXPECT_EQ(mio_nurbs_info_count(mapped, MIO_NURBS_PATCHES), 1);
    mio_mesh_free(cells);
    mio_format_info_free(mapped);
    mio_format_info_free(refined);
    mio_format_info_free(info);
    std::remove(path.c_str());
    std::remove(output.c_str());
}

TEST(CApi, MeshCadSurfaceMeshesABrepAndReportsItsFaces) {
    const auto path = mt::temp_path(".cad.json");
    std::ofstream(path, std::ios::binary) << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2,
        "surface": {"is_trimmed": true, "degrees": [1, 1], "knot_vectors": [[0, 2], [0, 2]],
          "control_points": [[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]},
        "boundary_loops": [{"loop_type": "outer", "trimming_curves": [
          {"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 0, 0, 1]], [[1.5, 0, 0, 1]]]}},
          {"trim_index": 1, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[1.5, 0, 0, 1]], [[0, 1.5, 0, 1]]]}},
          {"trim_index": 2, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 1.5, 0, 1]], [[0, 0, 0, 1]]]}}]}]}],
        "edges": []}]})";
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_free(sampled);

    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    EXPECT_EQ(opts.target_size, 0.0);
    EXPECT_EQ(opts.edges, 1);
    EXPECT_EQ(opts.record_face_ids, 1);
    EXPECT_EQ(opts.max_points, 50000000);
    EXPECT_EQ(opts.cells_across_gaps, 0.0);
    // The size is required.
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    EXPECT_EQ(mio_mesh_cad_surface(info, nullptr, nullptr), nullptr);
    opts.target_size = 0.1;
    mio_mesh_cad_surface_report report;
    mio_mesh* mesh = mio_mesh_cad_surface(info, &opts, &report);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_faces, 1);
    EXPECT_EQ(report.num_failed, 0);
    EXPECT_EQ(report.num_misoriented_edges, 0);
    EXPECT_GT(report.num_triangles, 100);
    EXPECT_EQ(report.num_points, mio_mesh_num_points(mesh));
    EXPECT_GT(report.min_angle_deg, 10.0);
    EXPECT_LE(report.max_edge, 0.16);  // up to 1.6 h beside the 45-degree corners
    EXPECT_EQ(report.reserved[0], 0);
    EXPECT_EQ(mio_mesh_cad_surface_face_status(0), MIO_MESH_CAD_OK);
    EXPECT_EQ(mio_mesh_cad_surface_face_status(1), -1);
    EXPECT_EQ(mio_mesh_cad_surface_face_reason(0, nullptr, 0), 0);
    EXPECT_EQ(mio_mesh_cad_surface_warnings(nullptr, 0), 0);
    EXPECT_EQ(opts.min_angle_deg, 0.0);
    EXPECT_EQ(opts.smooth_iterations, 0);
    EXPECT_EQ(opts.smooth_method, MIO_MESH_CAD_SMOOTH_LAPLACIAN);
    EXPECT_EQ(report.num_refinement_points, 0);
    EXPECT_EQ(report.min_input_angle_deg, 180.0);
    EXPECT_EQ(opts.order, 1);
    EXPECT_EQ(report.num_curved_straightened, 0);
    const std::int64_t linear_points = mio_mesh_num_points(mesh);
    mio_mesh_free(mesh);
    // v17.16.0: curved elements, quadratic and cubic; 0 is linear.
    for (const int order : {2, 3, 0}) {
        opts.order = order;
        mesh = mio_mesh_cad_surface(info, &opts, &report);
        ASSERT_NE(mesh, nullptr) << mio_last_error();
        char type[32] = {};
        ASSERT_GT(mio_mesh_cell_block_type(mesh, 0, type, sizeof type), 0);
        EXPECT_STREQ(type, order == 2 ? "triangle6" : order == 3 ? "triangle10" : "triangle");
        if (order == 0)
            EXPECT_EQ(mio_mesh_num_points(mesh), linear_points);
        else
            EXPECT_GT(mio_mesh_num_points(mesh), linear_points);
        EXPECT_EQ(report.num_curved_straightened, 0);
        mio_mesh_free(mesh);
    }
    opts.order = 4;
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    opts.order = 1;
    mesh = mio_mesh_cad_surface(info, &opts, &report);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    mio_mesh_free(mesh);
    // v17.18.0: the size field and pole rings, from the former reserved tail.
    {
        mio_mesh_cad_surface_opts fresh;
        mio_mesh_cad_surface_opts_init(&fresh);
        EXPECT_EQ(fresh.pole_rings, 0);
        EXPECT_EQ(fresh.grading, 0.0);
        EXPECT_EQ(fresh.curvature_deviation, 0.0);
        EXPECT_EQ(fresh.cells_across_gaps, 0.0);
    }
    opts.curvature_deviation = 0.01;
    opts.cells_across_gaps = 2.0;
    opts.grading = 1.3;
    opts.pole_rings = 2;
    mesh = mio_mesh_cad_surface(info, &opts, &report);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_failed, 0);
    mio_mesh_free(mesh);
    opts.grading = 0.5;
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    opts.grading = 0.0;
    opts.pole_rings = 9;
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    opts.pole_rings = 0;
    opts.curvature_deviation = 0.0;
    opts.cells_across_gaps = 0.0;
    // v17.19.0: given sizes through mio_mesh_cad_surface2's extended options.
    {
        mio_mesh_cad_surface_opts_ex ex;
        mio_mesh_cad_surface_opts_ex_init(&ex);
        EXPECT_EQ(ex.metric, MIO_MESH_CAD_METRIC_MEAN);
        EXPECT_EQ(ex.num_face_sizes, 0);
        EXPECT_EQ(ex.background, nullptr);
        mio_mesh_cad_surface_report plain{}, finer{};
        mio_mesh* a = mio_mesh_cad_surface2(info, &opts, nullptr, &plain);
        ASSERT_NE(a, nullptr) << mio_last_error();
        const int64_t face = 0;
        const double size = 0.25 * opts.target_size;
        ex.face_size_faces = &face;
        ex.face_sizes = &size;
        ex.num_face_sizes = 1;
        mio_mesh* b = mio_mesh_cad_surface2(info, &opts, &ex, &finer);
        ASSERT_NE(b, nullptr) << mio_last_error();
        EXPECT_GT(finer.num_triangles, plain.num_triangles);
        EXPECT_EQ(finer.num_failed, 0);
        int64_t limits[4] = {-7, -7, -7, -7};
        EXPECT_GE(mio_mesh_cad_surface_face_size_limits(0, limits, 4), 0);
        EXPECT_EQ(mio_mesh_cad_surface_face_size_limits(99, limits, 4), -1);
        EXPECT_EQ(mio_mesh_cad_surface_num_virtual_edges(), 0);
        EXPECT_EQ(mio_mesh_cad_surface_virtual_edge(0, limits, 4), -1);
        mio_mesh_free(a);
        mio_mesh_free(b);
        const int64_t bad_face = 99;
        ex.face_size_faces = &bad_face;
        EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
        ex.face_size_faces = &face;
        ex.metric = 7;
        EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
        ex.metric = MIO_MESH_CAD_METRIC_MEAN;
        // Virtual edges: none here (three corners, all sharp); refused curved.
        ex.virtual_edges = 1;
        mio_mesh* v = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
        ASSERT_NE(v, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_cad_surface_num_virtual_edges(), 0);
        mio_mesh_free(v);
        // v17.20.0: the reparametrisation, virtual faces and model
        // tolerances; a planar face is never flattened, nothing is a sliver.
        EXPECT_EQ(ex.reparam_condition_number, 0.0);
        EXPECT_EQ(ex.model_tolerances, 0);
        EXPECT_EQ(ex.virtual_faces, 0);
        ex.reparam_condition_number = 4.0;
        ex.model_tolerances = 1;
        ex.virtual_faces = 1;
        mio_mesh* x = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
        ASSERT_NE(x, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_cad_surface_num_reparametrised(), 0);
        EXPECT_EQ(mio_mesh_cad_surface_num_merged_faces(), 0);
        EXPECT_EQ(mio_mesh_cad_surface_num_tolerance_uses(), 0);
        int64_t idx = 0;
        EXPECT_NE(mio_mesh_cad_surface_reparametrised(0, &idx, nullptr, nullptr, nullptr), MIO_OK);
        EXPECT_NE(mio_mesh_cad_surface_merged_face(0, nullptr, nullptr, nullptr, nullptr), MIO_OK);
        EXPECT_NE(mio_mesh_cad_surface_tolerance_use(0, nullptr, nullptr, nullptr, nullptr),
                  MIO_OK);
        EXPECT_EQ(mio_mesh_cad_surface_merged_face_edges(0, limits, 4), -1);
        mio_mesh_free(x);
        ex.reparam_condition_number = 0.5;
        EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("reparam_condition_number"),
                  std::string::npos);
        ex.reparam_condition_number = 0.0;
        ex.model_tolerances = 0;
        ex.virtual_faces = 0;
        opts.order = 2;
        EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("virtual_edges"), std::string::npos);
        opts.order = 1;
        ex.virtual_edges = 0;
        // Roadmap 6.10.2: a metric from a constant Hessian asking for an
        // eighth of target_size along x refines the face.
        EXPECT_EQ(ex.metric_background, nullptr);
        EXPECT_EQ(ex.metric_array, nullptr);
        EXPECT_EQ(ex.metric_error, 0.0);
        mio_mesh* bg = mio_mesh_cad_surface2(info, &opts, nullptr, nullptr);
        ASSERT_NE(bg, nullptr) << mio_last_error();
        const int64_t np = mio_mesh_num_points(bg);
        std::vector<double> hess(static_cast<std::size_t>(9 * np), 0.0);
        const double hx = opts.target_size / 8.0;
        for (int64_t i = 0; i < np; ++i)
            hess[static_cast<std::size_t>(9 * i)] = 8.0 / (hx * hx);
        const int64_t shape[2] = {np, 9};
        ASSERT_EQ(mio_mesh_add_point_data(bg, "h", MIO_FLOAT64, 2, shape, hess.data()), MIO_OK);
        ex.metric_background = bg;
        ex.metric_array = "h";
        ex.metric_error = 1.0;
        mio_mesh_cad_surface_report metric{};
        mio_mesh* m = mio_mesh_cad_surface2(info, &opts, &ex, &metric);
        ASSERT_NE(m, nullptr) << mio_last_error();
        EXPECT_EQ(metric.num_failed, 0);
        EXPECT_GT(metric.num_triangles, plain.num_triangles);
        mio_mesh_free(m);
        ex.metric_error = 0.0;
        EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("metric_error"), std::string::npos);
        ex.metric_background = nullptr;
        ex.metric_array = nullptr;
        mio_mesh_free(bg);
    }
    // v17.15.0: refined to 25 degrees and smoothed. The triangle's two
    // 45-degree corners are reported; every angle is counted.
    opts.min_angle_deg = 25.0;
    opts.smooth_iterations = 2;
    opts.smooth_method = MIO_MESH_CAD_SMOOTH_ODT;
    mesh = mio_mesh_cad_surface(info, &opts, &report);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_failed, 0);
    EXPECT_GE(report.min_angle_deg, 25.0 - 1e-9);
    EXPECT_GT(report.num_smoothed, 0);
    EXPECT_NEAR(report.min_input_angle_deg, 45.0, 1e-9);
    std::int64_t bins[MIO_MESH_CAD_ANGLE_BINS] = {};
    EXPECT_EQ(mio_mesh_cad_surface_angle_histogram(bins, MIO_MESH_CAD_ANGLE_BINS),
              MIO_MESH_CAD_ANGLE_BINS);
    std::int64_t total = 0;
    for (const std::int64_t b : bins)
        total += b;
    EXPECT_EQ(total, 3 * report.num_triangles);
    EXPECT_EQ(mio_mesh_cad_surface_angle_histogram(nullptr, 0), MIO_MESH_CAD_ANGLE_BINS);
    EXPECT_EQ(mio_mesh_cad_surface_angle_histogram(nullptr, 3), -1);
    ASSERT_EQ(mio_mesh_cad_surface_num_corners(), 2);
    for (std::int64_t k = 0; k < 2; ++k) {
        mio_mesh_cad_corner c;
        ASSERT_EQ(mio_mesh_cad_surface_corner(k, &c), MIO_OK) << mio_last_error();
        EXPECT_EQ(c.face, 0);
        EXPECT_EQ(c.pole, 0);
        EXPECT_NEAR(c.angle_deg, 45.0, 1e-9);
        EXPECT_GE(c.point, 0);
    }
    mio_mesh_cad_corner c;
    EXPECT_EQ(mio_mesh_cad_surface_corner(2, &c), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_cad_surface_corner(0, nullptr), MIO_ERR_INVALID_ARG);
    // v17.20.0: corner lopping at 20 degrees lops the two 45-degree corners
    // of the planar face and states the straight-segment bound, 20.7.
    {
        mio_mesh_cad_surface_opts_ex ex;
        mio_mesh_cad_surface_opts_ex_init(&ex);
        EXPECT_EQ(ex.lop_corners, 0);
        ex.lop_corners = 1;
        opts.min_angle_deg = 20.0;
        mio_mesh* l = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
        ASSERT_NE(l, nullptr) << mio_last_error();
        ASSERT_EQ(mio_mesh_cad_surface_num_lopped_corners(), 2);
        for (std::int64_t k = 0; k < 2; ++k) {
            int64_t face = -1, point = -1, forced = -1;
            double angle = 0, lfs0 = 0, radius = 0, turning = -1, bound = 0, min_forced = 0;
            int32_t within = 0, guaranteed = 0;
            ASSERT_EQ(mio_mesh_cad_surface_lopped_corner(k, &face, &point, &angle, &lfs0, &radius,
                                                         &turning, &bound, &forced, &min_forced,
                                                         &within, &guaranteed),
                      MIO_OK)
                << mio_last_error();
            EXPECT_EQ(face, 0);
            EXPECT_GE(point, 0);
            EXPECT_NEAR(angle, 45.0, 1e-9);
            EXPECT_GT(radius, 0.0);
            EXPECT_LE(radius, (std::sqrt(2.0) - 1.0) * lfs0 * (1 + 1e-12));
            EXPECT_NEAR(turning, 0.0, 1e-9);
            EXPECT_NEAR(bound, std::asin(1.0 / (2.0 * std::sqrt(2.0))) * 180.0 / M_PI, 1e-9);
            EXPECT_EQ(within, 1);
            EXPECT_EQ(guaranteed, 1);
            char why[64];
            EXPECT_EQ(mio_mesh_cad_surface_lopped_corner_reason(k, why, sizeof why), 0);
        }
        EXPECT_NE(mio_mesh_cad_surface_lopped_corner(2, nullptr, nullptr, nullptr, nullptr, nullptr,
                                                     nullptr, nullptr, nullptr, nullptr, nullptr,
                                                     nullptr),
                  MIO_OK);
        EXPECT_EQ(mio_mesh_cad_surface_lopped_corner_reason(2, nullptr, 0), -1);
        mio_mesh_free(l);
        opts.min_angle_deg = 25.0;
    }
    // Not implemented yet, or out of range: refused.
    opts.smooth_method = 7;
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    opts.smooth_method = MIO_MESH_CAD_SMOOTH_LAPLACIAN;
    opts.min_angle_deg = 60.0;
    EXPECT_EQ(mio_mesh_cad_surface(info, &opts, nullptr), nullptr);
    mio_mesh_free(mesh);
    mio_format_info_free(info);
}

// Roadmap 6.5.1: trims the file-level `model_tolerance` pairs are reported
// as MIO_MESH_CAD_TOLERANCE_PAIRING, not as a projection. The box without
// edge topology, one face translated by 1e-5, as in test_mesh_cad_reparam.
TEST(CApi, MeshCadSurfaceReportsTheModelTolerancesPairing) {
    meshioplusplus::BrepInfo b = brep_solids::box();
    b.mEdges.clear();
    meshioplusplus::NDArray& r_cp =
        b.mPatches[static_cast<std::size_t>(b.mFaces[0].mPatch)].mControlPoints;
    for (std::size_t i = 0; i < r_cp.Size(); ++i)
        r_cp.As<double>()[i] += 1e-5;
    b.mExtra = R"({"tolerances":{"model_tolerance":0.0001}})";
    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, b);
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_free(sampled);
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.2;
    mio_mesh_cad_surface_opts_ex ex;
    mio_mesh_cad_surface_opts_ex_init(&ex);
    ex.model_tolerances = 1;
    mio_mesh* m = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_cad_surface_num_tolerance_uses(), 4);
    for (int64_t i = 0; i < 4; ++i) {
        int64_t edge = 0;
        int32_t where = -1;
        ASSERT_EQ(mio_mesh_cad_surface_tolerance_use(i, &edge, nullptr, nullptr, &where), MIO_OK);
        EXPECT_EQ(edge, -1);
        EXPECT_EQ(where, MIO_MESH_CAD_TOLERANCE_PAIRING);
    }
    mio_mesh_free(m);
    mio_format_info_free(info);
}

// Roadmap 6.6.3 and 6.8 through the extended options: untangling (refused
// at order 1) and quadrilaterals, quad-dominant and all-quad, on the box.
TEST(CApi, MeshCadSurfaceUntanglesAndMakesQuads) {
    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, brep_solids::box());
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_free(sampled);
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.2;
    opts.min_angle_deg = 25.0;
    mio_mesh_cad_surface_opts_ex ex;
    mio_mesh_cad_surface_opts_ex_init(&ex);
    EXPECT_EQ(ex.untangle, 0);
    EXPECT_EQ(ex.quads, MIO_MESH_CAD_QUADS_NONE);
    EXPECT_EQ(ex.quad_pairs, MIO_MESH_CAD_QUAD_PAIRS_MATCHING);
    // Untangling needs curved elements.
    ex.untangle = 1;
    EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("untangle"), std::string::npos);
    opts.order = 2;
    mio_mesh* curved = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
    ASSERT_NE(curved, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_cad_surface_num_untangled(), 0);  // nothing folds on a box
    mio_mesh_free(curved);
    opts.order = 1;
    ex.untangle = 0;
    const auto block_of = [](const mio_mesh* pM, const char* pType) {
        char buf[32];
        for (int64_t b = 0; b < mio_mesh_num_cell_blocks(pM); ++b)
            if (mio_mesh_cell_block_type(pM, b, buf, sizeof buf) >= 0 && std::string(buf) == pType)
                return b;
        return int64_t{-1};
    };
    for (const int32_t q : {MIO_MESH_CAD_QUADS_DOMINANT, MIO_MESH_CAD_QUADS_ALL}) {
        ex.quads = q;
        mio_mesh_cad_surface_report report{};
        mio_mesh* m = mio_mesh_cad_surface2(info, &opts, &ex, &report);
        ASSERT_NE(m, nullptr) << mio_last_error();
        int64_t quads = -1, left = -1;
        double sj = -1.0;
        ASSERT_EQ(mio_mesh_cad_surface_quad_stats(&quads, &left, &sj), MIO_OK);
        EXPECT_GT(quads, 0) << q;
        EXPECT_GT(sj, 0.3) << q;
        EXPECT_GE(block_of(m, "quad"), 0) << q;
        EXPECT_EQ(report.num_triangles, left) << q;
        EXPECT_EQ(report.num_open_edges, 0) << q;
        if (q == MIO_MESH_CAD_QUADS_ALL) {
            EXPECT_EQ(left, 0);
            EXPECT_EQ(block_of(m, "triangle"), -1);
        }
        mio_mesh_free(m);
    }
    EXPECT_EQ(mio_mesh_cad_surface_quad_stats(nullptr, nullptr, nullptr), MIO_OK);
    // No metric background: no metric angle (roadmap 6.10.2).
    EXPECT_EQ(mio_mesh_cad_surface_min_metric_angle(), 0.0);
    // Greedy pairing (6.8.1) leaves at least as many triangles.
    ex.quads = MIO_MESH_CAD_QUADS_DOMINANT;
    int64_t matched_left = -1, greedy_left = -1;
    mio_mesh* mm = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
    ASSERT_NE(mm, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_cad_surface_quad_stats(nullptr, &matched_left, nullptr), MIO_OK);
    mio_mesh_free(mm);
    ex.quad_pairs = MIO_MESH_CAD_QUAD_PAIRS_GREEDY;
    mm = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
    ASSERT_NE(mm, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_cad_surface_quad_stats(nullptr, &greedy_left, nullptr), MIO_OK);
    mio_mesh_free(mm);
    EXPECT_LE(matched_left, greedy_left);
    ex.quad_pairs = MIO_MESH_CAD_QUAD_PAIRS_MATCHING;
    // All-quad at order 2 writes quad9 (6.8.3).
    ex.quads = MIO_MESH_CAD_QUADS_ALL;
    opts.order = 2;
    mio_mesh* q9 = mio_mesh_cad_surface2(info, &opts, &ex, nullptr);
    ASSERT_NE(q9, nullptr) << mio_last_error();
    EXPECT_GE(block_of(q9, "quad9"), 0);
    EXPECT_GE(block_of(q9, "line3"), 0);
    mio_mesh_free(q9);
    opts.order = 1;
    ex.quads = 7;
    EXPECT_EQ(mio_mesh_cad_surface2(info, &opts, &ex, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("quads"), std::string::npos);
    mio_format_info_free(info);
}

// Roadmap 6.9: layers grown outward from the box's mesh are wedges, all
// valid, under a closed cap.
TEST(CApi, BoundaryLayer) {
    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, brep_solids::box());
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.25;
    mio_mesh* m = mio_mesh_cad_surface(info, &opts, nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();
    mio_boundary_layer_opts bo;
    mio_boundary_layer_opts_init(&bo);
    EXPECT_EQ(bo.num_layers, 3);
    EXPECT_DOUBLE_EQ(bo.growth, 1.2);
    EXPECT_EQ(mio_boundary_layer(m, &bo, nullptr), nullptr);  // no first_height
    EXPECT_NE(std::string(mio_last_error()).find("first_height"), std::string::npos);
    bo.first_height = 0.01;
    mio_boundary_layer_report report{};
    mio_mesh* out = mio_boundary_layer(m, &bo, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_GT(report.num_wedges, 0);
    EXPECT_EQ(report.num_cells, report.num_wedges);
    EXPECT_EQ(report.num_collapsed, 0);
    EXPECT_EQ(report.cap_closed, 1);
    EXPECT_GT(report.min_scaled_jacobian, 0.1);
    mio_mesh_free(out);
    EXPECT_EQ(mio_boundary_layer(m, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
    mio_mesh_free(sampled);
    mio_format_info_free(info);
}

// Roadmap 6.9.1: the fields taken from the opts' and report's reserved tails
// reach the core: the defaults, a layers collapse, smoothed normals and the
// front test switched off, and an out-of-range collapse refused by name.
TEST(CApi, BoundaryLayerCollapseAndFront) {
    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, brep_solids::box());
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.25;
    mio_mesh* m = mio_mesh_cad_surface(info, &opts, nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();
    mio_boundary_layer_opts bo;
    mio_boundary_layer_opts_init(&bo);
    EXPECT_EQ(bo.collapse, MIO_BOUNDARY_LAYER_COLLAPSE_POINTS);
    EXPECT_EQ(bo.normal_smooth_iterations, 0);
    EXPECT_EQ(bo.front_check, 1);
    EXPECT_EQ(bo.reserved_pad, 0);
    EXPECT_EQ(bo.reserved[0], 0);
    bo.first_height = 0.01;
    mio_boundary_layer_report report;
    std::memset(&report, 0xff, sizeof(report));
    mio_mesh* out = mio_boundary_layer(m, &bo, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_fewer_layers, 0);
    EXPECT_EQ(report.num_front_crossings, 0);
    EXPECT_EQ(report.reserved[0], 0);
    EXPECT_EQ(report.reserved[1], 0);
    const int64_t wedges = report.num_wedges;
    mio_mesh_free(out);
    bo.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_LAYERS;
    bo.normal_smooth_iterations = 4;
    bo.front_check = 0;
    out = mio_boundary_layer(m, &bo, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_wedges, wedges);
    EXPECT_EQ(report.num_front_crossings, 0);
    EXPECT_EQ(report.cap_closed, 1);
    mio_mesh_free(out);
    bo.collapse = 7;
    EXPECT_EQ(mio_boundary_layer(m, &bo, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("collapse"), std::string::npos)
        << mio_last_error();
    bo.collapse = MIO_BOUNDARY_LAYER_COLLAPSE_POINTS;
    bo.normal_smooth_iterations = 1001;
    EXPECT_EQ(mio_boundary_layer(m, &bo, nullptr), nullptr);
    mio_mesh_free(m);
    mio_mesh_free(sampled);
    mio_format_info_free(info);
}

// Roadmap 6.7: the box's skin filled with tetrahedra, its triangles kept;
// without Netgen the call fails naming the build flag.
TEST(CApi, Tetrahedralize) {
    const auto block_of = [](const mio_mesh* pM, const char* pType) {
        char buf[32];
        for (int64_t b = 0; b < mio_mesh_num_cell_blocks(pM); ++b)
            if (mio_mesh_cell_block_type(pM, b, buf, sizeof buf) >= 0 && std::string(buf) == pType)
                return b;
        return int64_t{-1};
    };
    mio_tetrahedralize_opts to;
    mio_tetrahedralize_opts_init(&to);
    EXPECT_EQ(to.backend, MIO_TETRAHEDRALIZE_BACKEND_AUTO);
    EXPECT_DOUBLE_EQ(to.grading, 0.3);
    EXPECT_EQ(to.optimize_steps, 3);
    EXPECT_EQ(to.max_cells, 50000000);
    EXPECT_EQ(to.face_array, nullptr);
    EXPECT_EQ(to.max_size, 0.0);
    EXPECT_EQ(mio_has_netgen() != 0, meshioplusplus::tetrahedralize_has_netgen());
    mio_tetrahedralize_opts_init(nullptr);

    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, brep_solids::box());
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.25;
    mio_mesh* m = mio_mesh_cad_surface(info, &opts, nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();

    EXPECT_EQ(mio_tetrahedralize(nullptr, nullptr, nullptr), nullptr);
    to.grading = 2.0;
    EXPECT_EQ(mio_tetrahedralize(m, &to, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("grading"), std::string::npos);
    mio_tetrahedralize_opts_init(&to);
    to.backend = 9;
    EXPECT_EQ(mio_tetrahedralize(m, &to, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("backend"), std::string::npos);
    mio_tetrahedralize_opts_init(&to);

    mio_tetrahedralize_report report;
    std::memset(&report, 0xff, sizeof(report));
    // The native octree backend, always there: the boundary fitted.
    {
        mio_tetrahedralize_opts oct;
        mio_tetrahedralize_opts_init(&oct);
        oct.backend = MIO_TETRAHEDRALIZE_BACKEND_OCTREE;
        mio_mesh* out = mio_tetrahedralize(m, &oct, &report);
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_EQ(report.backend, MIO_TETRAHEDRALIZE_BACKEND_OCTREE);
        EXPECT_EQ(report.boundary_preserved, 0);
        EXPECT_GT(report.hausdorff_to_surface, 0.0);
        EXPECT_GT(report.num_tetrahedra, 0);
        EXPECT_NEAR(report.volume, 1.0, 0.02);
        EXPECT_EQ(report.num_new_points, report.num_points);
        for (int64_t k : report.reserved)
            EXPECT_EQ(k, 0);
        EXPECT_GE(block_of(out, "triangle"), 0);
        mio_mesh_free(out);
    }
    if (!mio_has_netgen()) {
        // AUTO falls back to the octree backend; NETGEN is refused by name.
        mio_mesh* out = mio_tetrahedralize(m, &to, &report);
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_EQ(report.backend, MIO_TETRAHEDRALIZE_BACKEND_OCTREE);
        mio_mesh_free(out);
        to.backend = MIO_TETRAHEDRALIZE_BACKEND_NETGEN;
        EXPECT_EQ(mio_tetrahedralize(m, &to, &report), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_NETGEN"),
                  std::string::npos)
            << mio_last_error();
    } else {
        to.max_size = 0.25;
        mio_mesh* out = mio_tetrahedralize(m, &to, &report);
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_GT(report.num_tetrahedra, 0);
        EXPECT_EQ(report.num_regions, 1);
        EXPECT_EQ(report.boundary_preserved, 1);
        EXPECT_EQ(report.hausdorff_to_surface, 0.0);
        EXPECT_EQ(report.backend, MIO_TETRAHEDRALIZE_BACKEND_NETGEN);
        EXPECT_NEAR(report.volume, 1.0, 1e-9);
        EXPECT_GT(report.min_dihedral_deg, 0.0);
        EXPECT_LT(report.max_dihedral_deg, 180.0);
        EXPECT_EQ(report.num_points, mio_mesh_num_points(out));
        EXPECT_EQ(report.num_points - report.num_new_points, mio_mesh_num_points(m));
        for (int64_t k : report.reserved)
            EXPECT_EQ(k, 0);
        EXPECT_GE(block_of(out, "tetra"), 0);
        EXPECT_GE(block_of(out, "triangle"), 0);
        mio_mesh_free(out);
        // Past max_cells, refused.
        to.max_cells = 2;
        EXPECT_EQ(mio_tetrahedralize(m, &to, nullptr), nullptr);
        // The defaults (NULL opts) mesh too.
        out = mio_tetrahedralize(m, nullptr, nullptr);
        ASSERT_NE(out, nullptr) << mio_last_error();
        mio_mesh_free(out);
    }
    mio_mesh_free(m);
    mio_mesh_free(sampled);
    mio_format_info_free(info);
}

// Roadmap 6.10.1: a mesh of the box projected back onto it moves nothing,
// every point lands on a face or an edge, and none outside its trims.
TEST(CApi, ProjectToCad) {
    const auto path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, meshioplusplus::Mesh{}, brep_solids::box());
    mio_format_info* info = nullptr;
    mio_mesh* sampled = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(sampled, nullptr) << mio_last_error();
    mio_mesh_cad_surface_opts opts;
    mio_mesh_cad_surface_opts_init(&opts);
    opts.target_size = 0.25;
    mio_mesh* m = mio_mesh_cad_surface(info, &opts, nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();
    mio_project_to_cad_opts po;
    mio_project_to_cad_opts_init(&po);
    EXPECT_EQ(po.trim_samples, 16);
    EXPECT_EQ(po.face_array, nullptr);
    mio_project_to_cad_report report{};
    mio_mesh* out = mio_project_to_cad(m, info, &po, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_projected, mio_mesh_num_points(m));
    EXPECT_EQ(report.num_corners, 8);
    EXPECT_GT(report.num_on_edges, 8);
    EXPECT_EQ(report.num_outside_trim, 0);
    EXPECT_LT(report.max_move, 1e-12);
    EXPECT_GE(report.num_on_curves, 0);
    EXPECT_LE(report.num_on_curves, report.num_on_edges);
    EXPECT_EQ(report.reserved[0], 0);
    EXPECT_EQ(mio_project_to_cad_outside_trim(nullptr, 0), 0);
    mio_mesh_free(out);
    // A mesh with no faces named, and a NULL model, are refused by name.
    EXPECT_EQ(mio_project_to_cad(sampled, info, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("cad:face"), std::string::npos);
    EXPECT_EQ(mio_project_to_cad(m, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
    mio_mesh_free(sampled);
    mio_format_info_free(info);
}

TEST(CApi, IgaOnABrepKeepsItsTrimsAndReportsInside) {
    // The same square face as BrepSideChannel with an outer loop of four trims.
    const auto path = mt::temp_path(".cad.json");
    std::ofstream(path, std::ios::binary) << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2,
        "surface": {"is_trimmed": true, "degrees": [1, 1], "knot_vectors": [[0, 2], [0, 2]],
          "control_points": [[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]},
        "boundary_loops": [{"loop_type": "outer", "trimming_curves": [
          {"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 0, 0, 1]], [[1.5, 0, 0, 1]]]}},
          {"trim_index": 1, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[1.5, 0, 0, 1]], [[0, 1.5, 0, 1]]]}},
          {"trim_index": 2, "curve_direction": true, "parameter_curve": {"degree": 1,
           "knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 1.5, 0, 1]], [[0, 0, 0, 1]]]}}]}]}],
        "edges": []}]})";
    mio_format_info* info = nullptr;
    mio_mesh* mesh = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    mio_mesh_free(mesh);
    mio_iga_opts opts;
    mio_iga_opts_init(&opts);
    opts.subdivisions[0] = 2;
    opts.subdivisions[1] = 3;
    opts.elevate[1] = 1;
    mio_format_info* refined = nullptr;
    mio_mesh* cells = mio_iga_refine(info, &opts, &refined);
    ASSERT_NE(cells, nullptr) << mio_last_error();
    ASSERT_NE(refined, nullptr);
    mio_mesh_free(cells);
    // (0.2, 0.2) lies inside the triangle, (1.8, 1.8) beyond its hypotenuse.
    const double uv[4] = {0.2, 0.2, 1.8, 1.8};
    mio_iga_evaluation* ev = mio_iga_evaluate(refined, 0, uv, 2, 1, 16);
    ASSERT_NE(ev, nullptr) << mio_last_error();
    const void* data = nullptr;
    mio_dtype dtype;
    int32_t ndim;
    int64_t shape[MIO_MAX_NDIM]{};
    ASSERT_EQ(mio_iga_evaluation_array(ev, MIO_IGA_INSIDE, 0, &data, &dtype, &ndim, shape), MIO_OK);
    EXPECT_EQ(dtype, MIO_UINT8);
    EXPECT_EQ(static_cast<const uint8_t*>(data)[0], 1);
    EXPECT_EQ(static_cast<const uint8_t*>(data)[1], 0);
    ASSERT_EQ(mio_iga_evaluation_array(ev, MIO_IGA_NORMALS, 0, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_NEAR(std::abs(static_cast<const double*>(data)[2]), 1.0, 1e-13);
    mio_iga_evaluation_free(ev);
    mio_format_info_free(refined);
    mio_format_info_free(info);
    std::remove(path.c_str());
}

// ---------------------------------------------------------------------------
// Selective reads (mio_read_ex) and the opaque file summary (mio_read_metadata)
// ---------------------------------------------------------------------------

TEST(CApi, WriteOptsInitMatchesPlainWrite) {
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    EXPECT_EQ(opts.encoding, MIO_ENCODING_DEFAULT);
    EXPECT_EQ(opts.codec, MIO_CODEC_DEFAULT);
    EXPECT_EQ(opts.float_format, nullptr);
    for (int i = 0; i < 5; ++i)
        EXPECT_EQ(opts.reserved[i], 0) << "reserved must stay zero for ABI growth";
    // Same discipline as mio_read_opts: the tail is the growth budget, so its
    // width is part of the ABI and a field added later must come out of it.
    static_assert(sizeof(mio_write_opts::reserved) == 5 * sizeof(int64_t),
                  "mio_write_opts.reserved width is ABI");
}

TEST(CApi, WriteExHonoursEncodingAndCodec) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string ascii_path = mt::temp_path("_wex_ascii.vtu");
    const std::string binary_path = mt::temp_path("_wex_binary.vtu");

    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.encoding = MIO_ENCODING_ASCII;
    ASSERT_EQ(mio_write_ex(ascii_path.c_str(), m, "vtu", &opts), MIO_OK) << mio_last_error();

    mio_write_opts_init(&opts);
    opts.encoding = MIO_ENCODING_BINARY;
    opts.codec = MIO_CODEC_NONE;
    ASSERT_EQ(mio_write_ex(binary_path.c_str(), m, "vtu", &opts), MIO_OK) << mio_last_error();

    // The two encodings really produced different files, and both read back.
    mio_mesh* a = mio_read(ascii_path.c_str(), "vtu");
    mio_mesh* b = mio_read(binary_path.c_str(), "vtu");
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(mio_mesh_num_points(a), mio_mesh_num_points(m));
    EXPECT_EQ(mio_mesh_num_points(b), mio_mesh_num_points(m));
    mio_mesh_free(a);
    mio_mesh_free(b);
    std::remove(ascii_path.c_str());
    std::remove(binary_path.c_str());
    mio_mesh_free(m);
}

TEST(CApi, WriteExPcdLzfIsReachableAndRejectsIncompatibleOptions) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string path = mt::temp_path("_wex_lzf.pcd");
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.codec = MIO_CODEC_LZF;
    ASSERT_EQ(mio_write_ex(path.c_str(), m, "pcd", &opts), MIO_OK) << mio_last_error();
    mio_mesh* back = mio_read(path.c_str(), "pcd");
    ASSERT_NE(back, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(back), mio_mesh_num_points(m));
    mio_mesh_free(back);
    opts.encoding = MIO_ENCODING_ASCII;
    EXPECT_NE(mio_write_ex(path.c_str(), m, "pcd", &opts), MIO_OK);
    opts.encoding = MIO_ENCODING_BINARY;
    EXPECT_NE(mio_write_ex(path.c_str(), m, "vtu", &opts), MIO_OK);
    mio_mesh_free(m);
    std::remove(path.c_str());
}

TEST(CApi, GltfOptionsExposeColourAxisAndScale) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string path = mt::temp_path("_opts.gltf");
    mio_gltf_opts opts;
    mio_gltf_opts_init(&opts);
    EXPECT_EQ(opts.normals, 1);
    EXPECT_EQ(opts.scale, 1.0);
    const double temperature[] = {0, 1, 2, 3, 4};
    const int64_t shape[] = {5};
    ASSERT_EQ(mio_mesh_add_point_data(m, "temperature", MIO_FLOAT64, 1, shape, temperature),
              MIO_OK);
    opts.color_by = "temperature";
    opts.cmap = "turbo";
    opts.up_axis = 3;
    opts.scale = 0.001;
    opts.split_angle = 45;
    opts.vmin_set = opts.vmax_set = 1;
    opts.vmin = 0;
    opts.vmax = 3;
    ASSERT_EQ(mio_write_gltf(path.c_str(), m, &opts), MIO_OK) << mio_last_error();
    std::ifstream in(path);
    const std::string text{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    EXPECT_NE(text.find("COLOR_0"), std::string::npos);
    EXPECT_NE(text.find("KHR_materials_unlit"), std::string::npos);
    EXPECT_NE(text.find("0.001"), std::string::npos);
    opts.up_axis = 4;
    EXPECT_EQ(mio_write_gltf(path.c_str(), m, &opts), MIO_ERR_INVALID_ARG);
    opts.up_axis = 0;
    opts.scale = -1;
    EXPECT_NE(mio_write_gltf(path.c_str(), m, &opts), MIO_OK);
    EXPECT_EQ(mio_write_gltf(path.c_str(), nullptr, nullptr), MIO_ERR_INVALID_ARG);
    mio_mesh_free(m);
    std::filesystem::remove(path);
    std::filesystem::remove(std::filesystem::path(path).replace_extension(".bin"));
}

TEST(CApi, WriteExRawAppendedWritesVtuAndRefusesOthers) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string path = mt::temp_path("_wex_raw.vtu");
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.encoding = MIO_ENCODING_RAW_APPENDED;
    opts.codec = MIO_CODEC_NONE;
    ASSERT_EQ(mio_write_ex(path.c_str(), m, "vtu", &opts), MIO_OK) << mio_last_error();
    mio_mesh* back = mio_read(path.c_str(), "vtu");
    ASSERT_NE(back, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(back), mio_mesh_num_points(m));
    mio_mesh_free(back);
    std::remove(path.c_str());

    const std::string bad = mt::temp_path("_wex_raw.vtk");
    EXPECT_NE(mio_write_ex(bad.c_str(), m, "vtk", &opts), MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("raw appended"), std::string::npos)
        << mio_last_error();
    std::remove(bad.c_str());
    mio_mesh_free(m);
}

TEST(CApi, ExodusSeriesLifecycleOrExplicitMissingDependency) {
    const auto path = mt::temp_path("_c_series.e");
    EXPECT_EQ(mio_exodus_series_create(nullptr), nullptr);
    EXPECT_EQ(mio_exodus_series_num_steps(nullptr), -1);
    EXPECT_EQ(mio_exodus_series_flush(nullptr), MIO_ERR_INVALID_ARG);
    mio_exodus_series_free(nullptr);
    auto* series = mio_exodus_series_create(path.c_str());
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
    ASSERT_NE(series, nullptr) << mio_last_error();
    auto* mesh = build_tet_mesh();
    ASSERT_NE(mesh, nullptr);
    EXPECT_NE(mio_exodus_series_write_data(series, 0.0, mesh), MIO_OK);
    ASSERT_EQ(mio_exodus_series_write_points_cells(series, mesh), MIO_OK) << mio_last_error();
    ASSERT_EQ(mio_exodus_series_write_data(series, 0.25, mesh), MIO_OK) << mio_last_error();
    ASSERT_EQ(mio_exodus_series_write_data(series, 1.25, mesh), MIO_OK) << mio_last_error();
    EXPECT_EQ(mio_exodus_series_num_steps(series), 2);
    EXPECT_EQ(mio_exodus_series_flush(series), MIO_OK);
    EXPECT_EQ(mio_exodus_series_finalize(series), MIO_OK);
    EXPECT_EQ(mio_exodus_series_finalize(series), MIO_OK);
    EXPECT_EQ(mio_exodus_series_finalized(series), 1);
    EXPECT_NE(mio_exodus_series_write_data(series, 2.0, mesh), MIO_OK);
    mio_exodus_series_free(series);
    mio_mesh_free(mesh);
    std::filesystem::remove(path);
#else
    EXPECT_EQ(series, nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_NETCDF"), std::string::npos);
#endif
}

#ifndef MESHIOPLUSPLUS_HAS_HDF5
TEST(CApi, MedNamedFunctionsFailByNameWithoutHdf5) {
    EXPECT_EQ(mio_med_mesh_count("missing.med"), -1);
    EXPECT_NE(std::string(mio_last_error()).find("HDF5"), std::string::npos);
    EXPECT_EQ(mio_med_read_named("missing.med", "mesh", nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("HDF5"), std::string::npos);
}
#endif

TEST(CApi, WriteExRejectsAnOptionTheFormatCannotHonour) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string path = mt::temp_path("_wex_bad.msh");
    mio_write_opts opts;
    mio_write_opts_init(&opts);
    opts.codec = MIO_CODEC_ZSTD;  // gmsh has no block codec
    // Failing beats silently ignoring: asking for zstd here is a mistake.
    EXPECT_NE(mio_write_ex(path.c_str(), m, "gmsh", &opts), MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("codec"), std::string::npos);
    std::remove(path.c_str());
    mio_mesh_free(m);
}

TEST(CApi, WriteExWithNullOptsIsPlainWrite) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    const std::string path = mt::temp_path("_wex_null.vtu");
    EXPECT_EQ(mio_write_ex(path.c_str(), m, "vtu", nullptr), MIO_OK) << mio_last_error();
    std::remove(path.c_str());
    mio_mesh_free(m);
}

TEST(CApi, ReadOptsInitIsReadEverything) {
    mio_read_opts opts;
    mio_read_opts_init(&opts);
    EXPECT_EQ(opts.points_only, 0);
    EXPECT_EQ(opts.metadata_only, 0);
    EXPECT_EQ(opts.arrays, nullptr);
    EXPECT_EQ(opts.num_arrays, 0);
    EXPECT_EQ(opts.mmap_mode, 0);
    // 0 = the first step, which is the historical behaviour -- so a
    // default-initialized options struct still reads exactly what it always did.
    EXPECT_EQ(opts.time_step, 0);
    // 0 = throw on a construct the reader cannot represent, the historical
    // behaviour, so the defaults still read exactly what they always did.
    EXPECT_EQ(opts.lenient, 0);
    // No piece chosen: a partitioned file merges into one mesh, as it always did.
    EXPECT_EQ(opts.piece, 0);
    EXPECT_EQ(opts.piece_set, 0);
    // 0 = keep the ghost cells (halo) of a partitioned file, as a read always did.
    EXPECT_EQ(opts.drop_ghosts, 0);
    for (std::size_t i = 0; i < sizeof(opts.reserved) / sizeof(opts.reserved[0]); ++i)
        EXPECT_EQ(opts.reserved[i], 0) << "reserved must stay zero for ABI growth";
    // `time_step`, `lenient`, `piece`, `piece_set` and `drop_ghosts` each took one of the
    // six former reserved int64 slots rather than growing the struct, so the tail is
    // still exactly six int64s wide and a caller compiled against an older header
    // passes a correctly-sized object. Stated as the tail's width rather than
    // sizeof(the whole struct), which would be a padding assertion rather than an
    // ABI one.
    static_assert(sizeof(mio_read_opts::time_step) + sizeof(mio_read_opts::lenient) +
                          sizeof(mio_read_opts::piece) + sizeof(mio_read_opts::piece_set) +
                          sizeof(mio_read_opts::drop_ghosts) + sizeof(mio_read_opts::reserved) ==
                      6 * sizeof(std::int64_t),
                  "mio_read_opts grew: that is an ABI break, not additive growth");
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(CApi, ReadExPieceSelectsOnePartitionAndMergesByDefault) {
    // Two cell regions -> a two-block composite VTKHDF; the tags carry the block order.
    mt::Mesh m;
    m.AssignPoints(mt::points_from(
        {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {5, 0, 0}, {6, 0, 0}, {5, 1, 0}, {5, 0, 1}}));
    mt::NDArray conn = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {2, 4});
    const std::int64_t rows[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    std::copy(rows, rows + 8, conn.As<std::int64_t>());
    m.AddCellBlock("tetra", std::move(conn));
    mt::NDArray near_cells = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {1});
    near_cells.As<std::int64_t>()[0] = 0;
    mt::NDArray far_cells = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {1});
    far_cells.As<std::int64_t>()[0] = 1;
    m.AddRegion(meshioplusplus::Region("near", meshioplusplus::RegionKind::Cell, -1, 0,
                                       std::move(near_cells)));
    m.AddRegion(meshioplusplus::Region("far", meshioplusplus::RegionKind::Cell, -1, 1,
                                       std::move(far_cells)));
    const std::string path = mt::temp_path(".vtkhdf");
    meshioplusplus::write_vtkhdf(path, m, 4,
                                 meshioplusplus::VtkhdfType::PartitionedDataSetCollection);

    mio_mesh* whole = mio_read(path.c_str(), "vtkhdf");
    ASSERT_NE(whole, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(whole), 8);

    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.piece = 1;
    opts.piece_set = 1;
    mio_mesh* one = mio_read_ex(path.c_str(), "vtkhdf", &opts);
    ASSERT_NE(one, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(one), 4);  // the "far" block alone

    opts.piece = 2;  // out of range fails the call rather than clamping
    EXPECT_EQ(mio_read_ex(path.c_str(), "vtkhdf", &opts), nullptr);

    // A hand-zeroed struct must NOT ask for piece 0 alone: piece_set is what selects.
    mio_read_opts zeroed{};
    mio_mesh* merged = mio_read_ex(path.c_str(), "vtkhdf", &zeroed);
    ASSERT_NE(merged, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(merged), 8);

    mio_mesh_free(merged);
    mio_mesh_free(one);
    mio_mesh_free(whole);
    std::filesystem::remove(path);
}
#endif

TEST(CApi, ReadExPointsOnlyDropsDataKeepsGeometry) {
    const std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, mt::data_mesh(), /*binary=*/true, /*zlib=*/false);

    mio_mesh* full = mio_read(path.c_str(), "vtu");
    ASSERT_NE(full, nullptr);

    mio_read_opts opts;
    mio_read_opts_init(&opts);
    opts.points_only = 1;
    mio_mesh* bare = mio_read_ex(path.c_str(), "vtu", &opts);
    ASSERT_NE(bare, nullptr) << mio_last_error();

    EXPECT_EQ(mio_mesh_num_point_data(bare), 0);
    EXPECT_EQ(mio_mesh_num_cell_data(bare), 0);
    EXPECT_GT(mio_mesh_num_point_data(full), 0);
    EXPECT_EQ(mio_mesh_num_points(bare), mio_mesh_num_points(full));
    EXPECT_EQ(mio_mesh_num_cell_blocks(bare), mio_mesh_num_cell_blocks(full));

    mio_mesh_free(bare);
    mio_mesh_free(full);
    std::filesystem::remove(path);
}

TEST(CApi, ReadExArraysSubsetAndExplicitNone) {
    const std::string path = mt::temp_path(".vtu");
    const meshioplusplus::Mesh source = mt::data_mesh();
    ASSERT_GE(source.PointDataNames().size(), 2u);
    const std::string keep = source.PointDataNames().front();
    meshioplusplus::write_vtu(path, source, /*binary=*/true, /*zlib=*/false);

    mio_read_opts opts;
    mio_read_opts_init(&opts);
    const char* names[1] = {keep.c_str()};
    opts.arrays = names;
    opts.num_arrays = 1;
    mio_mesh* subset = mio_read_ex(path.c_str(), "vtu", &opts);
    ASSERT_NE(subset, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_point_data(subset), 1);
    mio_mesh_free(subset);

    // Non-NULL pointer with count 0 means "no arrays" -- distinct from NULL.
    mio_read_opts none;
    mio_read_opts_init(&none);
    none.arrays = names;
    none.num_arrays = 0;
    mio_mesh* empty = mio_read_ex(path.c_str(), "vtu", &none);
    ASSERT_NE(empty, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_point_data(empty), 0);
    mio_mesh_free(empty);

    std::filesystem::remove(path);
}

TEST(CApi, ReadExWithNullOptsMatchesPlainRead) {
    const std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, mt::data_mesh(), /*binary=*/true, /*zlib=*/false);

    mio_mesh* a = mio_read(path.c_str(), "vtu");
    mio_mesh* b = mio_read_ex(path.c_str(), "vtu", nullptr);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(mio_mesh_num_points(a), mio_mesh_num_points(b));
    EXPECT_EQ(mio_mesh_num_point_data(a), mio_mesh_num_point_data(b));
    mio_mesh_free(a);
    mio_mesh_free(b);
    std::filesystem::remove(path);
}

TEST(CApi, ReadMetadataReportsShapeAndNames) {
    const std::string path = mt::temp_path(".vtu");
    const meshioplusplus::Mesh source = mt::data_mesh();
    meshioplusplus::write_vtu(path, source, /*binary=*/true, /*zlib=*/false);

    mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), "vtu");
    ASSERT_NE(meta, nullptr) << mio_last_error();

    EXPECT_EQ(mio_read_metadata_num_points(meta), static_cast<int64_t>(source.NumPoints()));
    EXPECT_EQ(mio_read_metadata_num_cell_blocks(meta),
              static_cast<int64_t>(source.NumCellBlocks()));
    EXPECT_GT(mio_read_metadata_num_cells(meta), 0);
    EXPECT_EQ(mio_read_metadata_fell_back(meta), 0) << "vtu has a native metadata path";

    int64_t ncells = 0, npc = 0;
    int ragged = -1;
    ASSERT_EQ(mio_read_metadata_cell_block(meta, 0, &ncells, &npc, &ragged), MIO_OK);
    EXPECT_EQ(ncells, static_cast<int64_t>(source.Cells(0).NumCells()));
    EXPECT_EQ(ragged, 0);

    char buf[64];
    const int64_t n = mio_read_metadata_cell_block_type(meta, 0, buf, sizeof(buf));
    ASSERT_GT(n, 0);
    EXPECT_EQ(std::string(buf), std::string(source.Cells(0).Type()));

    EXPECT_EQ(mio_read_metadata_num_names(meta, MIO_DATA_POINT),
              static_cast<int64_t>(source.PointDataNames().size()));
    ASSERT_EQ(mio_read_metadata_name(meta, MIO_DATA_POINT, 0, buf, sizeof(buf)),
              static_cast<int64_t>(source.PointDataNames().front().size()));
    EXPECT_EQ(std::string(buf), source.PointDataNames().front());

    // A native summary never decodes the coordinates, so it has no bbox.
    EXPECT_EQ(mio_read_metadata_bbox(meta, nullptr, nullptr), MIO_ERR_NOT_FOUND);

    mio_read_metadata_free(meta);
    std::filesystem::remove(path);
}

TEST(CApi, ReadMetadataFallbackFlagsItselfAndHasBBox) {
    const std::string path = mt::temp_path(".stl");
    meshioplusplus::write_stl(path, mt::tri_mesh(), /*binary=*/false, /*skin=*/true);

    mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), "stl");
    ASSERT_NE(meta, nullptr) << mio_last_error();
    EXPECT_EQ(mio_read_metadata_fell_back(meta), 1);

    double lo[3], hi[3];
    EXPECT_EQ(mio_read_metadata_bbox(meta, lo, hi), MIO_OK);
    for (int d = 0; d < 3; ++d)
        EXPECT_LE(lo[d], hi[d]);

    mio_read_metadata_free(meta);
    std::filesystem::remove(path);
}

TEST(CApi, ReadMetadataErrorsAreGuardedNotThrown) {
    EXPECT_EQ(mio_read_metadata_create(nullptr, "vtu"), nullptr);
    EXPECT_EQ(mio_read_metadata_num_points(nullptr), -1);
    EXPECT_EQ(mio_read_metadata_fell_back(nullptr), -1);
    EXPECT_EQ(mio_read_metadata_cell_block(nullptr, 0, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_reader_supports_options(nullptr), -1);
    mio_read_metadata_free(nullptr);  // NULL-safe
    EXPECT_NE(std::string(mio_last_error()), "");
}

TEST(CApi, ReaderSupportsOptions) {
    EXPECT_EQ(mio_reader_supports_options("vtu"), 1);
    EXPECT_EQ(mio_reader_supports_options("stl"), 0);
}

TEST(CApi, ConvertCellsSimplexifyThroughTheAbi) {
    // One unit-cube hexahedron -> 6 tetra.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);

    mio_convert_cells_result* r = mio_convert_cells(m, "simplexify", /*record_parent_ids=*/1);
    ASSERT_NE(r, nullptr);
    const mio_mesh* out = mio_convert_cells_result_mesh(r);
    ASSERT_NE(out, nullptr);
    EXPECT_EQ(mio_mesh_num_points(out), 8);
    ASSERT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "tetra");
    std::int64_t num_cells = 0, nodes_per_cell = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(out, 0, &num_cells, &nodes_per_cell, nullptr), MIO_OK);
    EXPECT_EQ(num_cells, 6);
    EXPECT_EQ(nodes_per_cell, 4);

    // Zero-copy borrows of the index maps.
    const void* data = nullptr;
    mio_dtype dtype = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_convert_cells_result_point_map(r, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    EXPECT_EQ(n, 8);
    EXPECT_EQ(static_cast<const std::int64_t*>(data)[0], 0);

    ASSERT_EQ(mio_convert_cells_result_num_cell_maps(r), 1);
    ASSERT_EQ(mio_convert_cells_result_cell_map(r, 0, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(static_cast<const std::int64_t*>(data)[0], 0);  // parent 0's first child

    mio_mesh* owned = mio_convert_cells_result_take_mesh(r);
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(mio_mesh_num_points(owned), 8);
    mio_mesh_free(owned);
    mio_convert_cells_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, ConvertCellsLinearizeAndElevateRoundTrip) {
    mio_mesh* m = build_tet_mesh();

    mio_convert_cells_result* up = mio_convert_cells(m, "elevate", 0);
    ASSERT_NE(up, nullptr);
    const mio_mesh* quad = mio_convert_cells_result_mesh(up);
    EXPECT_EQ(block_type(quad, 0), "tetra10");
    EXPECT_GT(mio_mesh_num_points(quad), 5);

    mio_convert_cells_result* down = mio_convert_cells(quad, "linearize", 0);
    ASSERT_NE(down, nullptr);
    const mio_mesh* back = mio_convert_cells_result_mesh(down);
    EXPECT_EQ(block_type(back, 0), "tetra");
    EXPECT_EQ(mio_mesh_num_points(back), 5);

    mio_convert_cells_result_free(down);
    mio_convert_cells_result_free(up);
    mio_mesh_free(m);
}

TEST(CApi, ConvertCellsErrorsAreGuardedNotThrown) {
    mio_mesh* m = build_tet_mesh();
    // An unknown mode fails through the status/last-error contract, not by
    // letting the C++ exception escape the ABI.
    EXPECT_EQ(mio_convert_cells(m, "nope", 0), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_convert_cells(nullptr, "linearize", 0), nullptr);
    EXPECT_EQ(mio_convert_cells(m, nullptr, 0), nullptr);
    EXPECT_EQ(mio_convert_cells_result_mesh(nullptr), nullptr);
    EXPECT_EQ(mio_convert_cells_result_num_cell_maps(nullptr), -1);
    EXPECT_EQ(mio_convert_cells_result_point_map(nullptr, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_convert_cells_result_free(nullptr);  // NULL-safe
    mio_mesh_free(m);
}

TEST(CApi, SubdivideHexahedronThroughTheAbi) {
    // One unit-cube hexahedron -> 6 polyhedral children, one per face.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);

    mio_subdivide_result* r = mio_subdivide(m, /*record_parent_ids=*/1);
    ASSERT_NE(r, nullptr);
    const mio_mesh* out = mio_subdivide_result_mesh(r);
    ASSERT_NE(out, nullptr);
    EXPECT_EQ(mio_mesh_num_points(out), 9);  // one new interior point (apex)
    ASSERT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "polyhedron");

    mio_cell_block_info info{};
    ASSERT_EQ(mio_mesh_cell_block_info_ex(out, 0, &info), MIO_OK);
    EXPECT_EQ(info.num_cells, 6);
    EXPECT_EQ(info.is_ragged, 1);
    EXPECT_EQ(info.is_polyhedron, 1);

    // Zero-copy borrow of the cell map -- there is no point map (subdivide
    // never prunes or renumbers an original point).
    const void* data = nullptr;
    mio_dtype dtype = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_subdivide_result_num_cell_maps(r), 1);
    ASSERT_EQ(mio_subdivide_result_cell_map(r, 0, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(static_cast<const std::int64_t*>(data)[0], 0);

    mio_mesh* owned = mio_subdivide_result_take_mesh(r);
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(mio_mesh_num_points(owned), 9);
    mio_mesh_free(owned);
    mio_subdivide_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, SubdivideErrorsAreGuardedNotThrown) {
    mio_mesh* m = build_tet_mesh();
    EXPECT_EQ(mio_subdivide(nullptr, 0), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_subdivide_result_mesh(nullptr), nullptr);
    EXPECT_EQ(mio_subdivide_result_num_cell_maps(nullptr), -1);
    EXPECT_EQ(mio_subdivide_result_cell_map(nullptr, 0, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_subdivide_result_free(nullptr);  // NULL-safe
    mio_mesh_free(m);
}

TEST(CApi, AgglomerateTwoHexesThroughTheAbi) {
    // Two unit hexahedra sharing one face (x=1 plane).
    const std::vector<double> pts = {0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
                                     1, 1, 1, 1, 0, 1, 2, 0, 0, 2, 1, 0, 2, 1, 1, 2, 0, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 12, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 2, 8, MIO_INT64, conn.data()), MIO_OK);

    mio_agglomerate_result* r = mio_agglomerate(m, /*target_group_size=*/2);
    ASSERT_NE(r, nullptr);
    const mio_mesh* out = mio_agglomerate_result_mesh(r);
    ASSERT_NE(out, nullptr);
    ASSERT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "polyhedron");

    mio_cell_block_info info{};
    ASSERT_EQ(mio_mesh_cell_block_info_ex(out, 0, &info), MIO_OK);
    EXPECT_EQ(info.num_cells, 1);
    EXPECT_EQ(info.is_polyhedron, 1);

    // Zero-copy borrow of the FLAT cell map -- both hexes land in the one
    // merged cell.
    const void* data = nullptr;
    mio_dtype dtype = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_agglomerate_result_cell_map(r, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    EXPECT_EQ(n, 2);
    const auto* cm = static_cast<const std::int64_t*>(data);
    EXPECT_EQ(cm[0], cm[1]);

    mio_mesh* owned = mio_agglomerate_result_take_mesh(r);
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(owned), 1);
    mio_mesh_free(owned);
    mio_agglomerate_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, AgglomerateErrorsAreGuardedNotThrown) {
    mio_mesh* m = build_tet_mesh();
    EXPECT_EQ(mio_agglomerate(nullptr, 8), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_agglomerate(m, -1), nullptr);
    EXPECT_EQ(mio_agglomerate_result_mesh(nullptr), nullptr);
    EXPECT_EQ(mio_agglomerate_result_cell_map(nullptr, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_agglomerate_result_free(nullptr);  // NULL-safe
    mio_mesh_free(m);
}

TEST(CApi, SmoothMovesInteriorAndPinsBoundary) {
    // A 3x3 grid of quads: only the centre node is interior, so it is the only
    // one free to move -- which makes both halves of the assertion sharp.
    std::vector<double> pts;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            pts.insert(pts.end(), {static_cast<double>(i), static_cast<double>(j), 0.0});
    pts[4 * 3 + 0] += 0.4;  // push the centre node off-centre
    pts[4 * 3 + 1] -= 0.3;
    const std::vector<std::int64_t> conn = {0, 3, 4, 1, 1, 4, 5, 2, 3, 6, 7, 4, 4, 7, 8, 5};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 9, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 4, 4, MIO_INT64, conn.data()), MIO_OK);

    std::int64_t moved = -1, skipped = -1;
    double max_disp = -1.0;
    mio_mesh* out =
        mio_smooth(m, "laplacian", /*iterations=*/10, /*lambda=*/-1.0, /*mu=*/-0.34,
                   /*fix_boundary=*/1, /*preserve_features=*/1,
                   /*feature_angle=*/30.0, /*guard_inversion=*/1, &moved, &max_disp, &skipped);
    ASSERT_NE(out, nullptr);

    // Geometry only: same counts, same connectivity.
    EXPECT_EQ(mio_mesh_num_points(out), 9);
    ASSERT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "quad");

    EXPECT_EQ(moved, 1);  // exactly the one interior node
    EXPECT_GT(max_disp, 0.0);
    EXPECT_EQ(skipped, 0);

    // The centre node was pulled back toward (1, 1); the corners never moved.
    const void* data = nullptr;
    mio_dtype dtype = MIO_INT64;
    ASSERT_EQ(mio_mesh_get_points(out, &data, &dtype), MIO_OK);
    const double* p = static_cast<const double*>(data);
    // Laplacian converges geometrically: after 10 passes at lambda = 0.5 the
    // residual is the initial 0.4 offset times 0.5^10, i.e. ~3.9e-4.
    EXPECT_NEAR(p[4 * 3 + 0], 1.0, 1e-3);
    EXPECT_NEAR(p[4 * 3 + 1], 1.0, 1e-3);
    EXPECT_EQ(p[0], 0.0);
    EXPECT_EQ(p[8 * 3 + 0], 2.0);

    // A NULL method defaults to taubin rather than failing, and the counter
    // out-params are all individually optional.
    mio_mesh* dflt =
        mio_smooth(m, nullptr, 3, -1.0, -0.34, 1, 1, 30.0, 1, nullptr, nullptr, nullptr);
    EXPECT_NE(dflt, nullptr);
    mio_mesh_free(dflt);

    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, SmoothOdtMovesATetVertexTowardItsCircumcenter) {
    // A single, non-degenerate tet with one free vertex: ODT's closed-form
    // update pulls it toward the tet's own circumcenter -- mio_smooth's
    // method string is passed through to smooth_method_from_name unchanged,
    // so "odt" needs no C API surface of its own beyond the string.
    const std::vector<double> pts = {0, 0, 0, 4, 0, 0, 0, 4, 0, 1, 1, 5};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 1, 4, MIO_INT64, conn.data()), MIO_OK);

    std::int64_t moved = -1, skipped = -1;
    double max_disp = -1.0;
    mio_mesh* out =
        mio_smooth(m, "odt", /*iterations=*/1, /*lambda=*/-1.0, /*mu=*/-0.34,
                   /*fix_boundary=*/0, /*preserve_features=*/0,
                   /*feature_angle=*/30.0, /*guard_inversion=*/1, &moved, &max_disp, &skipped);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(out), 4);
    EXPECT_GT(moved, 0);
    EXPECT_GT(max_disp, 0.0);
    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, SmoothRejectsBadArguments) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 1, 4, MIO_INT64, conn.data()), MIO_OK);

    // Unknown method, lambda out of (0, 1), and a taubin mu that would amplify.
    EXPECT_EQ(mio_smooth(m, "bogus", 1, -1.0, -0.34, 1, 1, 30.0, 1, nullptr, nullptr, nullptr),
              nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_smooth(m, "taubin", 1, 1.5, -0.34, 1, 1, 30.0, 1, nullptr, nullptr, nullptr),
              nullptr);
    EXPECT_EQ(mio_smooth(m, "taubin", 1, 0.4, -0.2, 1, 1, 30.0, 1, nullptr, nullptr, nullptr),
              nullptr);
    // No exception may cross the ABI for a NULL mesh either.
    EXPECT_EQ(
        mio_smooth(nullptr, "taubin", 1, -1.0, -0.34, 1, 1, 30.0, 1, nullptr, nullptr, nullptr),
        nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    mio_mesh_free(m);
}

TEST(CApi, SmoothExFrozenMaskAndNullOptsMatchesDefaults) {
    // Same 3x3 grid as SmoothMovesInteriorAndPinsBoundary: only the centre
    // node is interior, so it is the only one smoothing may move by default.
    std::vector<double> pts;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            pts.insert(pts.end(), {static_cast<double>(i), static_cast<double>(j), 0.0});
    pts[4 * 3 + 0] += 0.4;
    pts[4 * 3 + 1] -= 0.3;
    const std::vector<std::int64_t> conn = {0, 3, 4, 1, 1, 4, 5, 2, 3, 6, 7, 4, 4, 7, 8, 5};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 9, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 4, 4, MIO_INT64, conn.data()), MIO_OK);

    // NULL opts is exactly mio_smooth_opts_init()'s defaults, i.e. plain
    // mio_smooth's own defaults (method NULL -> taubin, 10 iterations) --
    // compare against the plain call rather than a hardcoded convergence
    // target, since taubin's shrink-free two-pass-per-iteration scheme
    // converges at a different rate than laplacian's.
    mio_mesh* plain =
        mio_smooth(m, nullptr, 10, -1.0, -0.34, 1, 1, 30.0, 1, nullptr, nullptr, nullptr);
    ASSERT_NE(plain, nullptr) << mio_last_error();
    mio_mesh* dflt = mio_smooth_ex(m, nullptr, nullptr, nullptr, nullptr);
    ASSERT_NE(dflt, nullptr) << mio_last_error();
    const void* plain_data = nullptr;
    const void* dflt_data = nullptr;
    mio_dtype dtype = MIO_INT64;
    ASSERT_EQ(mio_mesh_get_points(plain, &plain_data, &dtype), MIO_OK);
    ASSERT_EQ(mio_mesh_get_points(dflt, &dflt_data, &dtype), MIO_OK);
    const double* plain_p = static_cast<const double*>(plain_data);
    const double* dflt_p = static_cast<const double*>(dflt_data);
    for (int i = 0; i < 9 * 3; ++i)
        EXPECT_EQ(plain_p[i], dflt_p[i]);
    // And the interior node genuinely moved, so the comparison above is not
    // vacuously true of two unchanged copies.
    EXPECT_NE(dflt_p[4 * 3 + 0], pts[4 * 3 + 0]);
    mio_mesh_free(plain);
    mio_mesh_free(dflt);

    // With fix_boundary off and the interior node explicitly frozen, it must
    // stay exactly where it started -- proving frozen reaches the core
    // independently of the boundary/feature pins.
    mio_smooth_opts opts;
    mio_smooth_opts_init(&opts);
    opts.fix_boundary = 0;
    const std::int64_t frozen_id = 4;
    opts.frozen = &frozen_id;
    opts.num_frozen = 1;
    std::int64_t moved = -1;
    mio_mesh* pinned = mio_smooth_ex(m, &opts, &moved, nullptr, nullptr);
    ASSERT_NE(pinned, nullptr) << mio_last_error();
    const void* pinned_data = nullptr;
    mio_dtype pinned_dtype = MIO_INT64;
    ASSERT_EQ(mio_mesh_get_points(pinned, &pinned_data, &pinned_dtype), MIO_OK);
    const double* pp = static_cast<const double*>(pinned_data);
    EXPECT_EQ(pp[4 * 3 + 0], pts[4 * 3 + 0]);
    EXPECT_EQ(pp[4 * 3 + 1], pts[4 * 3 + 1]);
    mio_mesh_free(pinned);

    // An out-of-range frozen id fails the call by name.
    const std::int64_t bad_id = 999;
    opts.frozen = &bad_id;
    EXPECT_EQ(mio_smooth_ex(m, &opts, nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("frozen node id 999"), std::string::npos);

    mio_mesh_free(m);
}

TEST(CApi, OptimizeVolumeFlipsAndImproves) {
    // The 2-tet fixture whose 2-3 flip raises the worst scaled Jacobian from
    // ~0.066 to ~0.362 (see test_optimize_volume.cpp). relocate off so only the
    // flip acts; every out-param is exercised.
    const std::vector<double> pts = {0.0, 0.0,    0.0,    1.0,    0.0,    0.0,    0.3,    0.9,
                                     0.0, 0.2252, 0.2808, 0.6977, 0.4659, 0.3149, -0.0367};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 0, 1, 2, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, conn.data()), MIO_OK);

    std::int64_t nflips = -1, n23 = -1, n32 = -1, moved = -1, ntets = -1;
    double qbefore = -1.0, qafter = -1.0;
    mio_mesh* out = mio_optimize_volume(m, /*max_iterations=*/5, /*relocate=*/0, /*flip=*/1,
                                        /*preserve_boundary=*/1, /*min_improvement=*/1e-6, &nflips,
                                        &n23, &n32, &moved, &ntets, &qbefore, &qafter);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(n23, 1);
    EXPECT_EQ(n32, 0);
    EXPECT_EQ(nflips, 1);
    EXPECT_EQ(ntets, 3);
    EXPECT_EQ(mio_mesh_num_points(out), 5);  // point set invariant
    EXPECT_GT(qafter, qbefore + 0.2);
    mio_mesh_free(out);

    // A non-tet mesh is refused by name; no exception crosses the ABI.
    mio_mesh* q = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(q, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    const std::vector<std::int64_t> quad = {0, 1, 2, 3};
    ASSERT_EQ(mio_mesh_add_cell_block(q, "quad", 1, 4, MIO_INT64, quad.data()), MIO_OK);
    EXPECT_EQ(mio_optimize_volume(q, 5, 1, 1, 1, 1e-6, nullptr, nullptr, nullptr, nullptr, nullptr,
                                  nullptr, nullptr),
              nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    mio_mesh_free(q);
    mio_mesh_free(m);
}

TEST(CApi, RefineHexIntoEight) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);

    mio_refine_result* r = mio_refine(m, /*levels=*/1, /*record_parent_ids=*/1);
    ASSERT_NE(r, nullptr);
    const mio_mesh* out = mio_refine_result_mesh(r);
    ASSERT_NE(out, nullptr);
    // 8 corners + 12 edge mids + 6 face centres + 1 body.
    EXPECT_EQ(mio_mesh_num_points(out), 27);
    ASSERT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "hexahedron");
    std::int64_t num_cells = 0, nodes_per_cell = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(out, 0, &num_cells, &nodes_per_cell, nullptr), MIO_OK);
    EXPECT_EQ(num_cells, 8);
    EXPECT_EQ(nodes_per_cell, 8);

    // Zero-copy borrows of the index maps.
    const void* data = nullptr;
    mio_dtype dtype = MIO_FLOAT64;
    std::int64_t n = -1;
    ASSERT_EQ(mio_refine_result_point_map(r, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    EXPECT_EQ(n, 8);
    EXPECT_EQ(static_cast<const std::int64_t*>(data)[0], 0);

    ASSERT_EQ(mio_refine_result_num_cell_maps(r), 1);
    ASSERT_EQ(mio_refine_result_cell_map(r, 0, &data, &dtype, &n), MIO_OK);
    EXPECT_EQ(n, 1);
    EXPECT_EQ(static_cast<const std::int64_t*>(data)[0], 0);  // parent 0's first child

    mio_mesh* owned = mio_refine_result_take_mesh(r);
    ASSERT_NE(owned, nullptr);
    EXPECT_EQ(mio_mesh_num_points(owned), 27);
    mio_mesh_free(owned);
    mio_refine_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, RefineExWithDefaultsMatchesMioRefine) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);

    // A NULL options pointer and a zero-initialized struct must both reproduce
    // mio_refine(mesh, 1, 0) exactly -- the append-only-tail contract.
    mio_refine_opts opts;
    mio_refine_opts_init(&opts);
    mio_refine_result* a = mio_refine(m, 1, 0);
    mio_refine_result* b = mio_refine_ex(m, &opts);
    mio_refine_result* c = mio_refine_ex(m, nullptr);
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    ASSERT_NE(c, nullptr);
    EXPECT_EQ(mio_mesh_num_points(mio_refine_result_mesh(a)), 27);
    EXPECT_EQ(mio_mesh_num_points(mio_refine_result_mesh(b)), 27);
    EXPECT_EQ(mio_mesh_num_points(mio_refine_result_mesh(c)), 27);
    mio_refine_result_free(a);
    mio_refine_result_free(b);
    mio_refine_result_free(c);
    mio_mesh_free(m);
}

TEST(CApi, RefineExRecordHierarchyAttachesCellIdAndParentId) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);

    mio_refine_opts opts;
    mio_refine_opts_init(&opts);
    opts.record_hierarchy = 1;
    mio_refine_result* r = mio_refine_ex(m, &opts);
    ASSERT_NE(r, nullptr);
    const mio_mesh* out = mio_refine_result_mesh(r);

    const void* data = nullptr;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[8] = {};
    ASSERT_EQ(mio_mesh_get_cell_data(out, "refine:cell_id", 0, &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(shape[0], 8);  // one hexahedron -> 8 children, uniform refinement
    ASSERT_EQ(mio_mesh_get_cell_data(out, "refine:parent_id", 0, &data, &dt, &ndim, shape), MIO_OK);
    const std::int64_t* parents = static_cast<const std::int64_t*>(data);
    for (int64_t i = 0; i < shape[0]; ++i)
        EXPECT_EQ(parents[i], 0) << "every child of the sole input cell shares parent 0";

    // Also proves the multigrid-stencil fix: RedGreen leaves no hanging nodes,
    // so refine:entity would normally never be attached at all.
    EXPECT_EQ(mio_mesh_get_point_data(out, "refine:entity", &data, &dt, &ndim, shape), MIO_OK);

    mio_refine_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, RefineExSelectsASubsetAndClosesUp) {
    // A 3 x 3 grid of quadrilaterals; refine the middle one only.
    std::vector<double> pts;
    for (int j = 0; j <= 3; ++j)
        for (int i = 0; i <= 3; ++i) {
            pts.push_back(i);
            pts.push_back(j);
            pts.push_back(0);
        }
    std::vector<std::int64_t> conn;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const std::int64_t a = j * 4 + i;
            conn.insert(conn.end(), {a, a + 1, a + 5, a + 4});
        }
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 16, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 9, 4, MIO_INT64, conn.data()), MIO_OK);

    const std::int64_t selected[] = {4};
    mio_refine_opts opts;
    mio_refine_opts_init(&opts);
    opts.cells = selected;
    opts.num_cells = 1;
    opts.record_levels = 1;
    mio_refine_result* r = mio_refine_ex(m, &opts);
    ASSERT_NE(r, nullptr);
    std::int64_t num_cells = 0, npc = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(mio_refine_result_mesh(r), 0, &num_cells, &npc, nullptr),
              MIO_OK);
    // More than the input, far fewer than the 36 a uniform refinement gives.
    EXPECT_GT(num_cells, 9);
    EXPECT_LT(num_cells, 36);
    EXPECT_EQ(npc, 4) << "green quadrilaterals stay quadrilaterals";
    mio_refine_result_free(r);

    // Two selectors at once is an error, reported rather than thrown.
    opts.region = "anything";
    EXPECT_EQ(mio_refine_ex(m, &opts), nullptr);
    EXPECT_NE(mio_last_error(), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, UndoGreenRestoresTheCoarseParentThroughTheAbi) {
    // A 3 x 3 grid of quadrilaterals; refine the middle one only, so the
    // neighbours around it get transitional (green) closures.
    std::vector<double> pts;
    for (int j = 0; j <= 3; ++j)
        for (int i = 0; i <= 3; ++i) {
            pts.push_back(i);
            pts.push_back(j);
            pts.push_back(0);
        }
    std::vector<std::int64_t> conn;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const std::int64_t a = j * 4 + i;
            conn.insert(conn.end(), {a, a + 1, a + 5, a + 4});
        }
    mio_mesh* coarse = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(coarse, MIO_FLOAT64, 16, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(coarse, "quad", 9, 4, MIO_INT64, conn.data()), MIO_OK);

    const std::int64_t selected[] = {4};
    mio_refine_opts opts;
    mio_refine_opts_init(&opts);
    opts.cells = selected;
    opts.num_cells = 1;
    opts.record_hierarchy = 1;
    opts.record_levels = 1;
    mio_refine_result* r = mio_refine_ex(coarse, &opts);
    ASSERT_NE(r, nullptr);
    const mio_mesh* fine = mio_refine_result_mesh(r);

    std::int64_t fine_cells = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(fine, 0, &fine_cells, nullptr, nullptr), MIO_OK);

    std::int64_t num_groups_undone = -1, num_cells_removed = -1;
    mio_mesh* undone = mio_undo_green(coarse, fine, &num_groups_undone, &num_cells_removed);
    ASSERT_NE(undone, nullptr);
    EXPECT_GT(num_groups_undone, 0);
    EXPECT_GT(num_cells_removed, 0);

    std::int64_t undone_cells = 0, coarse_cells = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(undone, 0, &undone_cells, nullptr, nullptr), MIO_OK);
    ASSERT_EQ(mio_mesh_cell_block_info(coarse, 0, &coarse_cells, nullptr, nullptr), MIO_OK);
    EXPECT_LT(undone_cells, fine_cells);
    EXPECT_GT(undone_cells, coarse_cells);

    // The reserved refine:* arrays are dropped entirely.
    const void* data = nullptr;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[8] = {};
    EXPECT_NE(mio_mesh_get_cell_data(undone, "refine:cell_id", 0, &data, &dt, &ndim, shape),
              MIO_OK);
    EXPECT_NE(mio_mesh_get_point_data(undone, "refine:entity", &data, &dt, &ndim, shape), MIO_OK);

    // Nullable counters really are nullable.
    mio_mesh* undone2 = mio_undo_green(coarse, fine, nullptr, nullptr);
    ASSERT_NE(undone2, nullptr);
    mio_mesh_free(undone2);

    mio_mesh_free(undone);
    mio_refine_result_free(r);
    mio_mesh_free(coarse);
}

TEST(CApi, UndoGreenErrorsAreGuardedNotThrown) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3};
    mio_mesh* coarse = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(coarse, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(coarse, "quad", 1, 4, MIO_INT64, conn.data()), MIO_OK);

    // No hierarchy at all on "fine" (here just the coarse mesh itself).
    EXPECT_EQ(mio_undo_green(coarse, coarse, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    // NULL meshes never let an exception cross the ABI.
    EXPECT_EQ(mio_undo_green(nullptr, coarse, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_undo_green(coarse, nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    mio_mesh_free(coarse);
}

TEST(CApi, RefineChainsOnABorrowedMesh) {
    mio_mesh* m = build_tet_mesh();  // 2 tetrahedra

    mio_refine_result* one = mio_refine(m, 1, 0);
    ASSERT_NE(one, nullptr);
    const mio_mesh* once = mio_refine_result_mesh(one);
    EXPECT_EQ(block_type(once, 0), "tetra");

    // Feed the borrowed mesh straight back in: levels=1 twice == levels=2.
    mio_refine_result* two = mio_refine(once, 1, 0);
    ASSERT_NE(two, nullptr);
    mio_refine_result* direct = mio_refine(m, 2, 0);
    ASSERT_NE(direct, nullptr);
    EXPECT_EQ(mio_mesh_num_points(mio_refine_result_mesh(two)),
              mio_mesh_num_points(mio_refine_result_mesh(direct)));

    mio_refine_result_free(direct);
    mio_refine_result_free(two);
    mio_refine_result_free(one);
    mio_mesh_free(m);
}

TEST(CApi, PartitionPiecesAndMaps) {
    // A 2x2 quad grid split into 2 parts: exactly nparts pieces, blocks kept
    // 1:1, cell maps assign every input cell to exactly one piece.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1,
                                     0, 2, 1, 0, 0, 2, 0, 1, 2, 0, 2, 2, 0};
    const std::vector<std::int64_t> conn = {0, 1, 4, 3, 1, 2, 5, 4, 3, 4, 7, 6, 4, 5, 8, 7};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 9, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 4, 4, MIO_INT64, conn.data()), MIO_OK);

    mio_partition_result* r =
        mio_partition(m, 2, "sfc", 0.03, "eco", 0, /*record_ids=*/1, /*ghost_layers=*/0, "");
    ASSERT_NE(r, nullptr);
    ASSERT_EQ(mio_partition_result_num_pieces(r), 2);
    ASSERT_EQ(mio_partition_result_num_cell_maps(r), 1);

    std::int64_t total_out = 0;
    std::vector<int> owners(4, 0);
    for (std::int64_t i = 0; i < 2; ++i) {
        EXPECT_EQ(mio_partition_result_part_id(r, i), static_cast<int>(i));
        const mio_mesh* piece = mio_partition_result_mesh(r, i);
        ASSERT_NE(piece, nullptr);
        ASSERT_EQ(mio_mesh_num_cell_blocks(piece), 1);  // blocks kept 1:1
        std::int64_t num_cells = 0;
        ASSERT_EQ(mio_mesh_cell_block_info(piece, 0, &num_cells, nullptr, nullptr), MIO_OK);
        total_out += num_cells;

        const void* data = nullptr;
        mio_dtype dtype = MIO_FLOAT64;
        std::int64_t n = -1;
        ASSERT_EQ(mio_partition_result_point_map(r, i, &data, &dtype, &n), MIO_OK);
        EXPECT_EQ(dtype, MIO_INT64);
        EXPECT_EQ(n, 9);
        ASSERT_EQ(mio_partition_result_cell_map(r, i, 0, &data, &dtype, &n), MIO_OK);
        EXPECT_EQ(n, 4);
        for (std::int64_t c = 0; c < n; ++c)
            if (static_cast<const std::int64_t*>(data)[c] >= 0)
                ++owners[static_cast<std::size_t>(c)];
    }
    EXPECT_EQ(total_out, 4);
    for (int o : owners)
        EXPECT_EQ(o, 1);  // partition of unity

    mio_mesh* owned = mio_partition_result_take_mesh(r, 0);
    ASSERT_NE(owned, nullptr);
    mio_mesh_free(owned);
    mio_partition_result_free(r);
    mio_mesh_free(m);
}

TEST(CApi, PartitionLabelsFillsTheCallerBuffer) {
    mio_mesh* m = build_tet_mesh();  // 2 tetrahedra
    std::vector<std::int64_t> labels(2, -1);
    ASSERT_EQ(mio_partition_labels(m, 2, "sfc", 0.03, "eco", 0, "", labels.data(),
                                   static_cast<std::int64_t>(labels.size())),
              MIO_OK);
    for (std::int64_t v : labels) {
        EXPECT_GE(v, 0);
        EXPECT_LT(v, 2);
    }
    // A wrong buffer size is rejected with a named message, not written past.
    EXPECT_EQ(mio_partition_labels(m, 2, "sfc", 0.03, "eco", 0, "", labels.data(), 5),
              MIO_ERR_INVALID_ARG);
    EXPECT_NE(std::string(mio_last_error()).find("labels_size"), std::string::npos);
    mio_mesh_free(m);
}

TEST(CApi, PartitionGhostLayers) {
    mio_mesh* m = build_tet_mesh();
    ASSERT_NE(m, nullptr);
    mio_partition_result* plain = mio_partition(m, 2, "sfc", 0.03, "eco", 0, 0, 0, "");
    mio_partition_result* halo = mio_partition(m, 2, "sfc", 0.03, "eco", 0, 0, 1, "");
    ASSERT_NE(plain, nullptr) << mio_last_error();
    ASSERT_NE(halo, nullptr) << mio_last_error();
    ASSERT_EQ(mio_partition_result_num_pieces(halo), mio_partition_result_num_pieces(plain));

    // With only two tetras sharing a face, one ghost layer makes each piece the
    // whole mesh -- and the tag must come along to say which cell is owned.
    const mio_mesh* piece = mio_partition_result_mesh(halo, 0);
    ASSERT_NE(piece, nullptr);
    EXPECT_GE(mio_mesh_num_cell_blocks(piece), 1);
    const void* data = nullptr;
    mio_dtype dt{};
    int32_t ndim = 0;
    int64_t shape[8]{};
    ASSERT_EQ(mio_mesh_get_cell_data(piece, "partition:ghost", 0, &data, &dt, &ndim, shape), MIO_OK)
        << mio_last_error();
    EXPECT_EQ(dt, MIO_INT64);

    mio_partition_result_free(plain);
    mio_partition_result_free(halo);
    mio_mesh_free(m);
}

TEST(CApi, PartitionErrorsAreGuardedNotThrown) {
    mio_mesh* m = build_tet_mesh();
    // Bad method name, bad nparts, negative ghost layers: NULL + last_error,
    // no throw. (A POSITIVE ghost_layers is a supported request since v9.0.0
    // and is exercised in PartitionGhostLayers below.)
    EXPECT_EQ(mio_partition(m, 2, "metis", 0.03, "eco", 0, 0, 0, ""), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_partition(m, 0, "sfc", 0.03, "eco", 0, 0, 0, ""), nullptr);
    EXPECT_EQ(mio_partition(m, 2, "sfc", 0.03, "eco", 0, 0, /*ghost_layers=*/-1, ""), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("ghost_layers"), std::string::npos);
    EXPECT_EQ(mio_partition(nullptr, 2, "sfc", 0.03, "eco", 0, 0, 0, ""), nullptr);
    EXPECT_EQ(mio_partition_result_num_pieces(nullptr), -1);
    EXPECT_EQ(mio_partition_result_part_id(nullptr, 0), -1);
    EXPECT_EQ(mio_partition_result_mesh(nullptr, 0), nullptr);
    EXPECT_EQ(mio_partition_result_point_map(nullptr, 0, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_partition_result_cell_map(nullptr, 0, 0, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_partition_result_free(nullptr);  // NULL-safe
    EXPECT_EQ(mio_partition_labels(nullptr, 2, "sfc", 0.03, "eco", 0, "", nullptr, 0),
              MIO_ERR_INVALID_ARG);
    mio_mesh_free(m);
}

TEST(CApi, RefineErrorsAreGuardedNotThrown) {
    // A pyramid has no same-type subdivision: the C++ exception must be turned
    // into NULL + last_error, never allowed to escape the ABI.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "pyramid", 1, 5, MIO_INT64, conn.data()), MIO_OK);

    EXPECT_EQ(mio_refine(m, 1, 0), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_refine(nullptr, 1, 0), nullptr);
    EXPECT_EQ(mio_refine_result_mesh(nullptr), nullptr);
    EXPECT_EQ(mio_refine_result_num_cell_maps(nullptr), -1);
    EXPECT_EQ(mio_refine_result_point_map(nullptr, nullptr, nullptr, nullptr), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_refine_result_cell_map(nullptr, 0, nullptr, nullptr, nullptr),
              MIO_ERR_INVALID_ARG);
    mio_refine_result_free(nullptr);  // NULL-safe
    mio_mesh_free(m);
}

TEST(CApi, IsosurfaceContoursAScalarPointField) {
    // The unit cube with f = x on its corners: the 0.5 level set is the unit
    // square at x = 0.5, on which f reads back as exactly the isovalue.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    const std::vector<double> fx = {0, 1, 1, 0, 0, 1, 1, 0};
    const std::int64_t shape[1] = {8};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "fx", MIO_FLOAT64, 1, shape, fx.data()), MIO_OK);

    const double isovalues[2] = {0.25, 0.75};
    mio_mesh* out = mio_isosurface(m, "fx", isovalues, 1, /*component=*/-1,
                                   /*record_parent_ids=*/1);
    ASSERT_NE(out, nullptr);
    ASSERT_GT(mio_mesh_num_points(out), 0);
    ASSERT_GT(mio_mesh_num_cell_blocks(out), 0);

    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int32_t ndim = 0;
    std::int64_t got_shape[MIO_MAX_NDIM] = {0};
    ASSERT_EQ(mio_mesh_get_point_data(out, "fx", &data, &dt, &ndim, got_shape), MIO_OK);
    ASSERT_EQ(dt, MIO_FLOAT64);
    const double* v = static_cast<const double*>(data);
    for (std::int64_t i = 0; i < got_shape[0]; ++i)
        EXPECT_EQ(v[i], 0.25) << "the contoured field must be exact";

    // Every contour cell is tagged with its value and its ordinal.
    EXPECT_GT(mio_mesh_cell_data_num_blocks(out, "iso:value"), 0);
    EXPECT_GT(mio_mesh_cell_data_num_blocks(out, "iso:index"), 0);
    EXPECT_GT(mio_mesh_cell_data_num_blocks(out, "iso:parent_cell"), 0);
    mio_mesh_free(out);

    // Two isovalues land in one mesh.
    mio_mesh* both = mio_isosurface(m, "fx", isovalues, 2, -1, 0);
    ASSERT_NE(both, nullptr);
    ASSERT_EQ(mio_mesh_get_cell_data(both, "iso:index", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_EQ(dt, MIO_INT64);
    mio_mesh_free(both);

    mio_mesh_free(m);
}

TEST(CApi, GradientCarriesTheOperatorsAndCounters) {
    // A frustum, not a cube: on a cube every face is a parallelogram whose
    // corner average IS its area centroid, so the exactness assertion below
    // would pass even with a broken quadrature.
    const std::vector<double> pts = {0,   0,   0, 2,   0,   0, 2,   2,   0, 0,   2,   0,
                                     0.5, 0.5, 1, 1.5, 0.5, 1, 1.5, 1.5, 1, 0.5, 1.5, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> f(8), u(24);
    for (std::size_t i = 0; i < 8; ++i) {
        const double x = pts[i * 3 + 0], y = pts[i * 3 + 1], z = pts[i * 3 + 2];
        f[i] = 3.0 * x - 2.0 * y + 5.0 * z + 7.0;
        u[i * 3 + 0] = 7.0 * z;
        u[i * 3 + 1] = 11.0 * x;
        u[i * 3 + 2] = 13.0 * y;
    }
    const std::int64_t sshape[1] = {8};
    const std::int64_t vshape[2] = {8, 3};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "f", MIO_FLOAT64, 1, sshape, f.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "u", MIO_FLOAT64, 2, vshape, u.data()), MIO_OK);

    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int32_t ndim = 0;
    std::int64_t got_shape[MIO_MAX_NDIM] = {0};

    std::int64_t skipped = -1, fallback = -1;
    mio_mesh* g =
        mio_gradient(m, "f", nullptr, nullptr, nullptr, nullptr, -1, 0, &skipped, &fallback);
    ASSERT_NE(g, nullptr) << mio_last_error();
    EXPECT_EQ(skipped, 0);
    EXPECT_EQ(fallback, 0);
    ASSERT_EQ(mio_mesh_get_cell_data(g, "f:gradient", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(got_shape[1], 3);
    const double* v = static_cast<const double*>(data);
    EXPECT_NEAR(v[0], 3.0, 1e-12);
    EXPECT_NEAR(v[1], -2.0, 1e-12);
    EXPECT_NEAR(v[2], 5.0, 1e-12);
    mio_mesh_free(g);

    // curl of (7z, 11x, 13y) is (13, 7, 11): three distinct nonzero components,
    // so any index permutation or sign flip fails.
    mio_mesh* c = mio_gradient(m, "u", "curl", "green-gauss", "cell", "w", -1, 0, nullptr, nullptr);
    ASSERT_NE(c, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_get_cell_data(c, "w", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    v = static_cast<const double*>(data);
    EXPECT_NEAR(v[0], 13.0, 1e-12);
    EXPECT_NEAR(v[1], 7.0, 1e-12);
    EXPECT_NEAR(v[2], 11.0, 1e-12);
    mio_mesh_free(c);

    // Least squares on a lone cell has no neighbourhood at all, so the fallback
    // counter must fire -- asserting it stays 0 on a nice mesh would be inert.
    skipped = fallback = -1;
    mio_mesh* l = mio_gradient(m, "f", "gradient", "least-squares", "cell", nullptr, -1, 0,
                               &skipped, &fallback);
    ASSERT_NE(l, nullptr) << mio_last_error();
    EXPECT_EQ(fallback, 1);
    EXPECT_EQ(skipped, 0);
    mio_mesh_free(l);

    // Point location moves the array to point_data and drops the intermediate.
    mio_mesh* p =
        mio_gradient(m, "f", "gradient", "green-gauss", "point", nullptr, -1, 0, nullptr, nullptr);
    ASSERT_NE(p, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_get_point_data(p, "f:gradient", &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_LE(mio_mesh_cell_data_num_blocks(p, "f:gradient"), 0)
        << "the intermediate cell array must be dropped";
    mio_mesh_free(p);

    mio_mesh_free(m);
}

TEST(CApi, GradientErrorsAreGuardedNotThrown) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3};
    const std::vector<double> h = {0, 1, 1, 0};
    const std::int64_t shape[1] = {4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 1, 4, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "h", MIO_FLOAT64, 1, shape, h.data()), MIO_OK);

    // Unknown array, unknown operator/method, a scalar divergence, an
    // out-of-range component and NULL arguments: NULL + last_error, never an
    // exception across the ABI.
    EXPECT_EQ(mio_gradient(m, "nope", nullptr, nullptr, nullptr, nullptr, -1, 0, nullptr, nullptr),
              nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_gradient(m, "h", "laplacian", nullptr, nullptr, nullptr, -1, 0, nullptr, nullptr),
              nullptr);
    EXPECT_EQ(mio_gradient(m, "h", nullptr, "magic", nullptr, nullptr, -1, 0, nullptr, nullptr),
              nullptr);
    EXPECT_EQ(
        mio_gradient(m, "h", "divergence", nullptr, nullptr, nullptr, -1, 0, nullptr, nullptr),
        nullptr)
        << "a scalar has no divergence";
    EXPECT_EQ(mio_gradient(m, "h", nullptr, nullptr, nullptr, nullptr, 7, 0, nullptr, nullptr),
              nullptr);
    EXPECT_EQ(
        mio_gradient(nullptr, "h", nullptr, nullptr, nullptr, nullptr, -1, 0, nullptr, nullptr),
        nullptr);
    EXPECT_EQ(mio_gradient(m, nullptr, nullptr, nullptr, nullptr, nullptr, -1, 0, nullptr, nullptr),
              nullptr);
    mio_mesh_free(m);
}

TEST(CApi, HessianCarriesCountersAndShape) {
    // A single hexahedron; a linear field has an exactly zero Hessian
    // regardless of mesh shape (see hessian.hpp), so this needs no
    // structured-grid fixture the way the C++ gtest suite's own
    // interior-exactness test does.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> f(8);
    for (std::size_t i = 0; i < 8; ++i) {
        const double x = pts[i * 3 + 0], y = pts[i * 3 + 1], z = pts[i * 3 + 2];
        f[i] = 3.0 * x - 2.0 * y + 5.0 * z + 7.0;
    }
    const std::int64_t sshape[1] = {8};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "f", MIO_FLOAT64, 1, sshape, f.data()), MIO_OK);

    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int32_t ndim = 0;
    std::int64_t got_shape[MIO_MAX_NDIM] = {0};

    std::int64_t skipped = -1, fallback = -1;
    mio_mesh* h = mio_hessian(m, "f", nullptr, nullptr, nullptr, 0, &skipped, &fallback);
    ASSERT_NE(h, nullptr) << mio_last_error();
    EXPECT_EQ(skipped, 0);
    EXPECT_EQ(fallback, 0);
    ASSERT_EQ(mio_mesh_get_cell_data(h, "f:hessian", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(got_shape[1], 9);
    const double* v = static_cast<const double*>(data);
    for (int i = 0; i < 9; ++i)
        EXPECT_NEAR(v[i], 0.0, 1e-9) << "entry " << i;
    mio_mesh_free(h);

    // Point location moves the array to point_data and drops the intermediate.
    mio_mesh* p = mio_hessian(m, "f", "green-gauss", "point", nullptr, 0, nullptr, nullptr);
    ASSERT_NE(p, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_get_point_data(p, "f:hessian", &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_LE(mio_mesh_cell_data_num_blocks(p, "f:hessian"), 0)
        << "the intermediate cell array must be dropped";
    mio_mesh_free(p);

    // Unknown array, a vector (multi-component) input, and NULL arguments:
    // NULL + last_error, never an exception across the ABI.
    EXPECT_EQ(mio_hessian(m, "nope", nullptr, nullptr, nullptr, 0, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    const std::vector<double> u(24, 0.0);
    const std::int64_t vshape[2] = {8, 3};
    ASSERT_EQ(mio_mesh_add_point_data(m, "u", MIO_FLOAT64, 2, vshape, u.data()), MIO_OK);
    EXPECT_EQ(mio_hessian(m, "u", nullptr, nullptr, nullptr, 0, nullptr, nullptr), nullptr)
        << "hessian is scalar-only";
    EXPECT_EQ(mio_hessian(nullptr, "f", nullptr, nullptr, nullptr, 0, nullptr, nullptr), nullptr);
    EXPECT_EQ(mio_hessian(m, nullptr, nullptr, nullptr, nullptr, 0, nullptr, nullptr), nullptr);

    mio_mesh_free(m);
}

TEST(CApi, EstimateError) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> f(8);
    for (std::size_t i = 0; i < 8; ++i) {
        const double x = pts[i * 3 + 0], y = pts[i * 3 + 1], z = pts[i * 3 + 2];
        f[i] = 3.0 * x - 2.0 * y + 5.0 * z + 7.0;
    }
    const std::int64_t sshape[1] = {8};

    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "f", MIO_FLOAT64, 1, sshape, f.data()), MIO_OK);

    double global_error = -1.0;
    std::int64_t skipped = -1, marked = -1;
    mio_mesh* out = mio_estimate_error(m, "f", nullptr, nullptr, 0.0, nullptr, nullptr, 0,
                                       &global_error, &skipped, &marked);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(skipped, 0);
    EXPECT_EQ(marked, 0);
    EXPECT_NEAR(global_error, 0.0, 1e-9);

    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int32_t ndim = 0;
    std::int64_t got_shape[MIO_MAX_NDIM] = {0};
    ASSERT_EQ(mio_mesh_get_cell_data(out, "error:zz", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    EXPECT_LE(mio_mesh_cell_data_num_blocks(out, "error:marked"), 0)
        << "no marking requested, so no marked array";
    mio_mesh_free(out);

    // Custom names + marking.
    double ge2 = -1.0;
    std::int64_t sk2 = -1, mk2 = -1;
    mio_mesh* out2 =
        mio_estimate_error(m, "f", "zz", "absolute", 1e-6, "ind", "flag", 1, &ge2, &sk2, &mk2);
    ASSERT_NE(out2, nullptr) << mio_last_error();
    EXPECT_EQ(mk2, 0);
    ASSERT_EQ(mio_mesh_get_cell_data(out2, "ind", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    ASSERT_EQ(mio_mesh_get_cell_data(out2, "flag", 0, &data, &dt, &ndim, got_shape), MIO_OK);
    EXPECT_EQ(dt, MIO_INT64);
    mio_mesh_free(out2);

    mio_mesh_free(m);
}

TEST(CApi, CheckConformity) {
    // A coarse quad beside two fine ones: the fine pair's shared node sits on
    // the coarse quad's edge, a T-junction, and the edges overlap partly.
    const std::vector<double> pts = {0, 0, 0, 1, 0,   0, 1, 1, 0, 0, 1,   0,
                                     2, 0, 0, 2, 0.5, 0, 2, 1, 0, 1, 0.5, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 1, 4, 5, 7, 7, 5, 6, 2};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 3, 4, MIO_INT64, conn.data()), MIO_OK);

    mio_conformity_opts o;
    mio_conformity_opts_init(&o);
    EXPECT_EQ(o.tolerance, 1e-9);
    EXPECT_EQ(o.attach_flags, 1);
    mio_conformity_report report{};
    report.reserved[0] = 7;
    mio_mesh* out = mio_check_conformity(m, &o, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.dimension, 2);
    EXPECT_EQ(report.num_t_junctions, 1);
    EXPECT_EQ(report.num_hanging_nodes, 0);
    EXPECT_EQ(report.num_partial_overlaps, 2);
    EXPECT_EQ(report.num_unwelded_facets, 0);
    EXPECT_EQ(report.conforming, 0);
    EXPECT_EQ(report.interface_mode, MIO_CONFORMITY_PROXIMITY);
    EXPECT_EQ(report.num_warnings, 0);
    EXPECT_EQ(report.reserved[0], 0);
    const void* data = nullptr;
    mio_dtype dt = MIO_FLOAT64;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {0};
    ASSERT_EQ(mio_mesh_get_point_data(out, "conformity:hanging", &data, &dt, &ndim, shape), MIO_OK);
    ASSERT_EQ(dt, MIO_INT64);
    const std::int64_t* flags = static_cast<const std::int64_t*>(data);
    for (std::int64_t p = 0; p < 8; ++p)
        EXPECT_EQ(flags[p], p == 7 ? 2 : 0) << p;
    mio_mesh_free(out);

    // NULL options are the defaults; a negative tolerance is refused.
    out = mio_check_conformity(m, nullptr, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    mio_mesh_free(out);
    o.tolerance = -1.0;
    EXPECT_EQ(mio_check_conformity(m, &o, nullptr), nullptr);
    EXPECT_EQ(mio_check_conformity_warnings(nullptr, 0), 0);
    mio_mesh_free(m);
}

TEST(CApi, ErrorSize) {
    // Two unit right triangles of the unit square (diameter sqrt 2), indicator
    // 1 and 4: with eta = 1 and N = 2 they ask 1 and a quarter (Zienkiewicz
    // and Zhu 1987); a point takes the least of its cells'.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 0, 2, 3};
    const std::vector<double> eta = {1.0, 4.0};
    const std::int64_t eshape[1] = {2};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 2, 3, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_append_cell_data(m, "error:zz", MIO_FLOAT64, 1, eshape, eta.data()), MIO_OK);

    mio_error_size_opts o;
    mio_error_size_opts_init(&o);
    EXPECT_EQ(o.error_array, nullptr);
    EXPECT_EQ(o.output, nullptr);
    EXPECT_EQ(o.target_error, 0.0);
    EXPECT_EQ(o.order, 1.0);
    EXPECT_EQ(o.overwrite, 0);
    o.target_error = 1.0;
    mio_error_size_report report{};
    report.reserved[0] = 7;
    mio_mesh* out = mio_error_size(m, &o, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_cells, 2);
    EXPECT_EQ(report.num_skipped, 0);
    EXPECT_EQ(report.num_clamped, 0);
    EXPECT_NEAR(report.min_size, 0.25, 1e-14);
    EXPECT_NEAR(report.max_size, 1.0, 1e-14);
    EXPECT_EQ(report.reserved[0], 0);
    const void* data = nullptr;
    mio_dtype dt = MIO_INT64;
    std::int32_t ndim = 0;
    std::int64_t shape[MIO_MAX_NDIM] = {0};
    ASSERT_EQ(mio_mesh_get_point_data(out, "size", &data, &dt, &ndim, shape), MIO_OK);
    ASSERT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(shape[0], 4);
    const double* s = static_cast<const double*>(data);
    EXPECT_NEAR(s[0], 0.25, 1e-14);
    EXPECT_NEAR(s[1], 1.0, 1e-14);
    EXPECT_NEAR(s[3], 0.25, 1e-14);

    // Named arrays, bounds (counted), and a NULL report.
    o.output = "h";
    o.min_size = 0.5;
    mio_mesh* bounded = mio_error_size(out, &o, &report);
    ASSERT_NE(bounded, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_clamped, 1);
    EXPECT_NEAR(report.min_size, 0.5, 1e-14);
    ASSERT_EQ(mio_mesh_get_point_data(bounded, "h", &data, &dt, &ndim, shape), MIO_OK);
    mio_mesh_free(bounded);
    mio_mesh* quiet = mio_error_size(m, &o, nullptr);
    ASSERT_NE(quiet, nullptr) << mio_last_error();
    mio_mesh_free(quiet);

    // Refusals come back as NULL with the reason, never thrown.
    o.output = nullptr;
    EXPECT_EQ(mio_error_size(out, &o, nullptr), nullptr) << "size exists, no overwrite";
    EXPECT_NE(std::string(mio_last_error()).find("overwrite"), std::string::npos);
    o.overwrite = 1;
    mio_mesh* again = mio_error_size(out, &o, nullptr);
    ASSERT_NE(again, nullptr) << mio_last_error();
    mio_mesh_free(again);
    o.overwrite = 0;
    o.target_error = 0.0;
    EXPECT_EQ(mio_error_size(m, &o, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("target_error"), std::string::npos);
    o.target_error = 1.0;
    o.min_size = 2.0;
    o.max_size = 1.0;
    EXPECT_EQ(mio_error_size(m, &o, nullptr), nullptr);
    o.min_size = o.max_size = 0.0;
    o.error_array = "nope";
    EXPECT_EQ(mio_error_size(m, &o, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("nope"), std::string::npos);
    EXPECT_EQ(mio_error_size(m, nullptr, nullptr), nullptr);
    EXPECT_EQ(mio_error_size(nullptr, &o, nullptr), nullptr);
    mio_error_size_opts_init(nullptr);

    mio_mesh_free(out);
    mio_mesh_free(m);
}

TEST(CApi, EstimateErrorErrorsAreGuardedNotThrown) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<double> f(8, 1.0);
    const std::int64_t sshape[1] = {8};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 1, 8, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "f", MIO_FLOAT64, 1, sshape, f.data()), MIO_OK);

    EXPECT_EQ(mio_estimate_error(m, "nope", nullptr, nullptr, 0.0, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_estimate_error(m, "f", "kelly", nullptr, 0.0, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr);
    EXPECT_EQ(mio_estimate_error(m, "f", nullptr, "median", 0.0, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr);
    EXPECT_EQ(mio_estimate_error(m, "f", nullptr, "fraction", 1.5, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr)
        << "out-of-range marking_value";
    EXPECT_EQ(mio_estimate_error(nullptr, "f", nullptr, nullptr, 0.0, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr);
    EXPECT_EQ(mio_estimate_error(m, nullptr, nullptr, nullptr, 0.0, nullptr, nullptr, 0, nullptr,
                                 nullptr, nullptr),
              nullptr);
    mio_mesh_free(m);
}

namespace {
// A regular octahedron -- the smallest closed 2-manifold triangle mesh with
// no coplanar faces, small enough to hand-write for a C API smoke test.
mio_mesh* c_api_octahedron() {
    const std::vector<double> pts = {1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
    const std::vector<std::int64_t> conn = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
                                            2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
    mio_mesh* m = mio_mesh_create();
    mio_mesh_set_points(m, MIO_FLOAT64, 6, 3, pts.data());
    mio_mesh_add_cell_block(m, "triangle", 8, 3, MIO_INT64, conn.data());
    return m;
}
}  // namespace

TEST(CApi, RemeshProducesTheRequestedNumberOfClusters) {
    mio_mesh* m = c_api_octahedron();

    std::int64_t num_clusters = -1, num_iterations = -1, num_isolated = -1, num_nonmanifold = -1;
    int subdivide_applied = -1;
    mio_mesh* out =
        mio_remesh(m, 30, -1, 10.0, 4, 100, 10, "isotropic", 0.0, 1, &num_clusters, &num_iterations,
                   &subdivide_applied, &num_isolated, &num_nonmanifold);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(num_clusters, 30);
    EXPECT_GT(subdivide_applied, 0) << "6 vertices cannot support 30 clusters without subdividing";
    EXPECT_EQ(mio_mesh_num_points(out), 30);
    EXPECT_GE(num_nonmanifold, 0);
    mio_mesh_free(out);

    // The quadric ("feature-preserving") metric, gradation and boundary
    // preservation are all accepted too.
    mio_mesh* out2 = mio_remesh(m, 30, -1, 10.0, 4, 100, 10, "quadric", 1.5, 0, nullptr, nullptr,
                                nullptr, nullptr, nullptr);
    ASSERT_NE(out2, nullptr) << mio_last_error();
    mio_mesh_free(out2);

    mio_mesh_free(m);
}

TEST(CApi, RemeshErrorsAreGuardedNotThrown) {
    mio_mesh* m = c_api_octahedron();

    EXPECT_EQ(mio_remesh(m, 3, -1, 10.0, 4, 100, 10, nullptr, 0.0, 1, nullptr, nullptr, nullptr,
                         nullptr, nullptr),
              nullptr)
        << "too few clusters";
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_remesh(m, 10, -1, 10.0, 4, 100, 10, "bogus", 0.0, 1, nullptr, nullptr, nullptr,
                         nullptr, nullptr),
              nullptr)
        << "unknown metric";
    EXPECT_EQ(mio_remesh(nullptr, 10, -1, 10.0, 4, 100, 10, nullptr, 0.0, 1, nullptr, nullptr,
                         nullptr, nullptr, nullptr),
              nullptr);

    mio_mesh_free(m);
}

TEST(CApi, RemeshExAcceptsAnisotropicAndReportsCounters) {
    mio_mesh* m = c_api_octahedron();

    mio_remesh_opts opts;
    mio_remesh_opts_init(&opts);
    opts.num_clusters = 30;
    opts.metric = "anisotropic";
    opts.max_anisotropy = 3.0;
    mio_remesh_report report;
    mio_mesh* out = mio_remesh_ex(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_GT(report.num_clusters, 0);
    EXPECT_EQ(mio_mesh_num_points(out), report.num_clusters);
    EXPECT_GE(report.num_isolated_clusters, 0);
    EXPECT_GE(report.num_non_manifold_vertices, 0);
    mio_mesh_free(out);

    // NULL report is fine.
    mio_mesh* out2 = mio_remesh_ex(m, &opts, nullptr);
    ASSERT_NE(out2, nullptr) << mio_last_error();
    mio_mesh_free(out2);

    // NULL opts falls back to mio_remesh_opts_init's own defaults, whose
    // num_clusters == 0 fails validation exactly like the flat mio_remesh's
    // own "too few clusters" case above.
    EXPECT_EQ(mio_remesh_ex(m, nullptr, nullptr), nullptr);

    // max_anisotropy away from the default under a non-anisotropic metric
    // is guarded, not silently ignored, exactly as in the C++ core.
    mio_remesh_opts bad;
    mio_remesh_opts_init(&bad);
    bad.num_clusters = 30;
    bad.metric = "isotropic";
    bad.max_anisotropy = 2.0;
    EXPECT_EQ(mio_remesh_ex(m, &bad, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    mio_mesh_free(m);
}

TEST(CApi, RemeshVolumeOptsInitDefaultsAndRemeshVolume) {
    mio_mesh* octa = c_api_octahedron();

    mio_remesh_volume_opts opts;
    mio_remesh_volume_opts_init(&opts);
    EXPECT_EQ(opts.cell_size, 0.0);
    EXPECT_DOUBLE_EQ(opts.padding_relative, 0.1);
    EXPECT_EQ(opts.max_cells, 20000000);
    EXPECT_EQ(opts.max_tets, 20000000);
    EXPECT_DOUBLE_EQ(opts.warp_fraction, 0.35);
    // mio_sdf_opts_init's own defaults are honoured for the embedded distance
    // options -- the mio_compute_sdf_opts nesting precedent.
    EXPECT_EQ(opts.distance.sign, MIO_SDF_PSEUDONORMAL);
    opts.cell_size = 0.4;
    opts.distance.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    mio_remesh_volume_report report;
    mio_mesh* out = mio_remesh_volume_ex(octa, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_GT(report.num_tets, 0);
    EXPECT_EQ(mio_mesh_num_cell_blocks(out), 1);
    EXPECT_EQ(block_type(out, 0), "tetra");
    EXPECT_GE(report.num_vertices_warped, 0);
    EXPECT_GE(report.num_tets_rejected, 0);
    EXPECT_GE(report.num_non_manifold_edges, 0);
    mio_mesh_free(out);

    // NULL report is fine.
    mio_mesh* out2 = mio_remesh_volume_ex(octa, &opts, nullptr);
    ASSERT_NE(out2, nullptr) << mio_last_error();
    mio_mesh_free(out2);

    // NULL opts is an error, mio_voxelize's own rule: a resolution or cell
    // size must be given, and there is no sensible all-defaults fallback.
    EXPECT_EQ(mio_remesh_volume_ex(octa, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    mio_mesh_free(octa);
}

TEST(CApi, RemeshVolumeErrorsAreGuardedNotThrown) {
    mio_mesh* octa = c_api_octahedron();
    mio_remesh_volume_opts opts;
    mio_remesh_volume_opts_init(&opts);
    opts.cell_size = 0.4;
    opts.warp_fraction = -0.1;  // negative warp fraction

    EXPECT_EQ(mio_remesh_volume_ex(octa, &opts, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");

    mio_mesh_free(octa);
}

TEST(CApi, IsosurfaceErrorsAreGuardedNotThrown) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::vector<std::int64_t> conn = {0, 1, 2, 3};
    const std::vector<double> h = {0, 1, 1, 0};
    const std::int64_t shape[1] = {4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 1, 4, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "h", MIO_FLOAT64, 1, shape, h.data()), MIO_OK);

    const double iso[1] = {0.5};
    // Unknown array, no isovalues, and NULL arguments: NULL + last_error, never
    // an exception across the ABI.
    EXPECT_EQ(mio_isosurface(m, "nope", iso, 1, -1, 0), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_isosurface(m, "h", iso, 0, -1, 0), nullptr);
    EXPECT_EQ(mio_isosurface(m, "h", iso, 1, 7, 0), nullptr);
    EXPECT_EQ(mio_isosurface(nullptr, "h", iso, 1, -1, 0), nullptr);
    EXPECT_EQ(mio_isosurface(m, nullptr, iso, 1, -1, 0), nullptr);
    EXPECT_EQ(mio_isosurface(m, "h", nullptr, 1, -1, 0), nullptr);
    mio_mesh_free(m);
}

// --- named regions (doc/regions.md) ------------------------------------------
// Regions are the first thing on this ABI that carries named *groups* of
// entities; before meshio++ 8.1 they never left the Python layer.

TEST(CApi, RegionsRoundTripThroughTheOpaqueHandle) {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 4, 0, 2, 3, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, conn.data()), MIO_OK);

    // Entries are given unsorted and with a duplicate: canonicalization is part
    // of the contract, so what comes back is sorted and de-duplicated.
    const std::int64_t nodes[4] = {3, 0, 3, 1};
    ASSERT_EQ(mio_mesh_add_region(m, "fixed", MIO_REGION_POINT, -1, -1, nodes, 4), MIO_OK);
    const std::int64_t cells[2] = {1, 0};
    ASSERT_EQ(mio_mesh_add_region(m, "solid", MIO_REGION_CELL, 3, 42, cells, 2), MIO_OK);
    const std::int64_t sides[4] = {1, 3, 0, 1};  // (cell, facet) pairs
    ASSERT_EQ(mio_mesh_add_region(m, "wall", MIO_REGION_SIDE, 2, -1, sides, 4), MIO_OK);

    mio_regions* r = mio_regions_create(m);
    ASSERT_NE(r, nullptr);
    ASSERT_EQ(mio_regions_count(r), 3);

    // Order is (kind, name, dim, tag): Point, then Cell, then Side.
    char buf[64];
    ASSERT_GT(mio_regions_name(r, 0, buf, sizeof(buf)), 0);
    EXPECT_EQ(std::string(buf), "fixed");
    ASSERT_GT(mio_regions_name(r, 2, buf, sizeof(buf)), 0);
    EXPECT_EQ(std::string(buf), "wall");

    mio_region_info info{};
    ASSERT_EQ(mio_regions_info(r, 0, &info), MIO_OK);
    EXPECT_EQ(info.kind, MIO_REGION_POINT);
    EXPECT_EQ(info.stride, 1);
    EXPECT_EQ(info.num_entries, 3);  // the duplicate 3 was dropped

    ASSERT_EQ(mio_regions_info(r, 1, &info), MIO_OK);
    EXPECT_EQ(info.kind, MIO_REGION_CELL);
    EXPECT_EQ(info.dim, 3);
    EXPECT_EQ(info.tag, 42);

    ASSERT_EQ(mio_regions_info(r, 2, &info), MIO_OK);
    EXPECT_EQ(info.kind, MIO_REGION_SIDE);
    EXPECT_EQ(info.stride, 2);
    EXPECT_EQ(info.num_entries, 2);

    std::int64_t count = 0;
    const std::int64_t* e = mio_regions_entries(r, 0, &count);
    ASSERT_NE(e, nullptr);
    ASSERT_EQ(count, 3);
    EXPECT_EQ(e[0], 0);
    EXPECT_EQ(e[1], 1);
    EXPECT_EQ(e[2], 3);

    e = mio_regions_entries(r, 2, &count);
    ASSERT_NE(e, nullptr);
    ASSERT_EQ(count, 4);  // 2 pairs, lexicographically sorted
    EXPECT_EQ(e[0], 0);
    EXPECT_EQ(e[1], 1);
    EXPECT_EQ(e[2], 1);
    EXPECT_EQ(e[3], 3);

    mio_regions_free(r);
    mio_mesh_free(m);
}

TEST(CApi, RegionErrorsAreGuardedNotThrown) {
    mio_mesh* m = mio_mesh_create();
    const std::int64_t one[1] = {0};

    // No exception may cross the ABI: every one of these is a status/NULL.
    EXPECT_EQ(mio_mesh_add_region(nullptr, "x", MIO_REGION_POINT, -1, -1, one, 1),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_add_region(m, nullptr, MIO_REGION_POINT, -1, -1, one, 1),
              MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_mesh_add_region(m, "x", MIO_REGION_POINT, -1, -1, nullptr, 1),
              MIO_ERR_INVALID_ARG);
    // A side region needs (cell, facet) pairs.
    EXPECT_EQ(mio_mesh_add_region(m, "x", MIO_REGION_SIDE, -1, -1, one, 1), MIO_ERR_INVALID_ARG);

    EXPECT_EQ(mio_regions_create(nullptr), nullptr);
    mio_regions* r = mio_regions_create(m);
    ASSERT_NE(r, nullptr);
    EXPECT_EQ(mio_regions_count(r), 0);
    mio_region_info info{};
    EXPECT_EQ(mio_regions_info(r, 0, &info), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_regions_entries(r, 5, nullptr), nullptr);
    EXPECT_EQ(mio_regions_name(r, 5, nullptr, 0), -1);
    mio_regions_free(r);
    mio_regions_free(nullptr);  // NULL is safe
    mio_mesh_free(m);
}

TEST(CApi, RegionsSurviveAnOperation) {
    // crop keeps both tetra, so the cell region survives whole; the point
    // region is remapped through the pruned point numbering.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 1};
    const std::vector<std::int64_t> conn = {0, 1, 2, 4, 0, 2, 3, 4};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, conn.data()), MIO_OK);
    const std::int64_t cells[1] = {0};
    ASSERT_EQ(mio_mesh_add_region(m, "solid", MIO_REGION_CELL, 3, 7, cells, 1), MIO_OK);

    const double lo[3] = {-1, -1, -1};
    const double hi[3] = {2, 2, 2};
    mio_mesh* out = mio_crop_bbox(m, lo, hi, 0, 0);
    ASSERT_NE(out, nullptr);

    mio_regions* r = mio_regions_create(out);
    ASSERT_NE(r, nullptr);
    ASSERT_EQ(mio_regions_count(r), 1);
    mio_region_info info{};
    ASSERT_EQ(mio_regions_info(r, 0, &info), MIO_OK);
    EXPECT_EQ(info.num_entries, 1);
    EXPECT_EQ(info.tag, 7);  // the format-native id rides along
    mio_regions_free(r);

    mio_mesh_free(out);
    mio_mesh_free(m);
}

// ---- transient XDMF ------------------------------------------------------
//
// The one writer the flat API exposes as a handle rather than a (path, mesh)
// call. Written here, then read back through the C reader one step at a time.

TEST(CApi, XdmfTimeSeries) {
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    const double pts[15] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 0.5};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, pts), MIO_OK) << mio_last_error();
    const int64_t conn[8] = {0, 1, 2, 4, 0, 2, 3, 4};
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 2, 4, MIO_INT64, conn), MIO_OK)
        << mio_last_error();

    const std::string path = mt::temp_path(".xdmf");

    const int64_t anchors[2] = {0, 4};
    const int64_t wall[2] = {0, 1};
    ASSERT_EQ(mio_mesh_add_region(m, "anchors", MIO_REGION_POINT, 0, 7, anchors, 2), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(m, "wall", MIO_REGION_SIDE, 2, 9, wall, 2), MIO_OK);

    mio_xdmf_series* s = mio_xdmf_series_create(path.c_str(), "XML", -1);
    ASSERT_NE(s, nullptr) << mio_last_error();
    ASSERT_EQ(mio_xdmf_series_write_points_cells(s, m), MIO_OK) << mio_last_error();
    const int64_t shape[1] = {5};
    for (int k = 0; k < 3; ++k) {
        const double t[5] = {100.0 * k, 100.0 * k + 1, 100.0 * k + 2, 100.0 * k + 3, 100.0 * k + 4};
        ASSERT_EQ(mio_mesh_add_point_data(m, "T", MIO_FLOAT64, 1, shape, t), MIO_OK)
            << mio_last_error();
        ASSERT_EQ(mio_xdmf_series_write_data(s, 0.25 * k, m), MIO_OK) << mio_last_error();
    }
    EXPECT_EQ(mio_xdmf_series_num_steps(s), 3);
    ASSERT_EQ(mio_xdmf_series_finalize(s), MIO_OK) << mio_last_error();
    mio_xdmf_series_free(s);

    // Every step's time value is reachable without decoding a payload.
    mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), nullptr);
    ASSERT_NE(meta, nullptr) << mio_last_error();
    ASSERT_EQ(mio_read_metadata_num_time_values(meta), 3);
    double times[3] = {0, 0, 0};
    ASSERT_EQ(mio_read_metadata_time_values(meta, times, 3), 3);
    for (int k = 0; k < 3; ++k)
        EXPECT_DOUBLE_EQ(times[k], 0.25 * k);
    mio_read_metadata_free(meta);

    // ... and each step reads back with its own values.
    for (int k = 0; k < 3; ++k) {
        mio_read_opts opts;
        mio_read_opts_init(&opts);
        opts.time_step = k;
        mio_mesh* out = mio_read_ex(path.c_str(), nullptr, &opts);
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_num_points(out), 5);
        mio_regions* regions = mio_regions_create(out);
        ASSERT_NE(regions, nullptr);
        ASSERT_EQ(mio_regions_count(regions), 2);
        mio_region_info region_info{};
        ASSERT_EQ(mio_regions_info(regions, 1, &region_info), MIO_OK);
        EXPECT_EQ(region_info.kind, MIO_REGION_SIDE);
        EXPECT_EQ(region_info.tag, 9);
        int64_t entry_count = 0;
        const auto* entries = mio_regions_entries(regions, 1, &entry_count);
        ASSERT_NE(entries, nullptr);
        EXPECT_EQ(entry_count, 2);
        EXPECT_EQ(entries[0], 0);
        EXPECT_EQ(entries[1], 1);
        mio_regions_free(regions);
        const void* data = nullptr;
        mio_dtype dt = MIO_FLOAT64;
        int32_t ndim = 0;
        int64_t got_shape[MIO_MAX_NDIM] = {0};
        ASSERT_EQ(mio_mesh_get_point_data(out, "T", &data, &dt, &ndim, got_shape), MIO_OK)
            << mio_last_error();
        ASSERT_EQ(dt, MIO_FLOAT64);
        ASSERT_EQ(got_shape[0], 5);
        const double* values = static_cast<const double*>(data);
        for (int i = 0; i < 5; ++i)
            EXPECT_DOUBLE_EQ(values[i], 100.0 * k + i) << "step " << k << " entry " << i;
        mio_mesh_free(out);
    }

    mio_mesh_free(m);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(CApi, XdmfTimeSeriesUnknownFormat) {
    const std::string path = mt::temp_path(".xdmf");
    EXPECT_EQ(mio_xdmf_series_create(path.c_str(), "Zarr", -1), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("Zarr"), std::string::npos);
}

// ---------------------------------------------------------------------------
// Settings pipeline (v9.11.0): JSON text only across the ABI.
// ---------------------------------------------------------------------------

#ifdef MESHIOPLUSPLUS_HAS_JSON

TEST(CApi, PipelineRunsASettingsFile) {
    ASSERT_EQ(mio_pipeline_has_json(), 1);
    const std::string in_path = mt::temp_path("_capi_pipe.vtk");
    const std::string out_path = mt::temp_path("_capi_pipe_out.vtk");
    const std::string settings_path = mt::temp_path("_capi_settings.json");
    {
        mio_mesh* m = build_tet_mesh();
        ASSERT_EQ(mio_write(in_path.c_str(), m, "vtk"), MIO_OK) << mio_last_error();
        mio_mesh_free(m);
    }
    {
        std::ofstream settings(settings_path);
        settings << R"({"Input": {"Path": ")" << in_path << R"("},
                       "Operations": [{"Op": "Quality"}],
                       "Output": {"Path": ")"
                 << out_path << R"("}})";
    }
    ASSERT_EQ(mio_pipeline_run_file(settings_path.c_str()), MIO_OK) << mio_last_error();
    mio_mesh* out = mio_read(out_path.c_str(), "vtk");
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_GE(mio_mesh_cell_data_num_blocks(out, "quality:scaled_jacobian"), 1);
    mio_mesh_free(out);
    std::error_code ec;
    std::filesystem::remove(in_path, ec);
    std::filesystem::remove(out_path, ec);
    std::filesystem::remove(settings_path, ec);
}

TEST(CApi, PipelineJsonSchemaErrorsNameTheOffender) {
    EXPECT_NE(mio_pipeline_run_json(R"({"Input": {"Path": "a"}, "Output": {"Path": "b"},
                                        "Operations": [{"Op": "Nope"}]})"),
              MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("Nope"), std::string::npos);
    EXPECT_NE(mio_pipeline_run_file(nullptr), MIO_OK);
    EXPECT_NE(mio_pipeline_run_json(nullptr), MIO_OK);
}

TEST(CApi, SequencePipelineRunsASettingsFile) {
    // The sequence document shares the pipeline's parser and guard.
    const std::string dir = mt::temp_path("_capi_seqdoc");
    std::filesystem::create_directories(dir);
    for (int i = 0; i < 3; ++i) {
        mio_mesh* m = build_tet_mesh();
        ASSERT_EQ(mio_write((dir + "/in_" + std::to_string(i) + ".vtu").c_str(), m, "vtu"), MIO_OK);
        mio_mesh_free(m);
    }
    const std::string settings = dir + "/settings.json";
    {
        std::ofstream f(settings);
        f << R"({"Version": 1, "Input": {"Pattern": ")" << dir
          << R"(/in_*.vtu"}, "Operations": [{"Op": "Quality"}], "Output": {"Path": ")" << dir
          << R"(/out_{step}.vtu"}})";
    }
    ASSERT_EQ(mio_sequence_pipeline_run_file(settings.c_str()), MIO_OK) << mio_last_error();
    EXPECT_TRUE(std::filesystem::exists(dir + "/out_0002.vtu"));
    EXPECT_NE(mio_sequence_pipeline_run_file(nullptr), MIO_OK);
    EXPECT_NE(mio_sequence_pipeline_run_json(nullptr), MIO_OK);
    std::filesystem::remove_all(dir);
}

#else  // !MESHIOPLUSPLUS_HAS_JSON

TEST(CApi, PipelineCompiledOutFailsNamingTheFlag) {
    // Never a link error, never a silent no-op: the entry points exist and
    // the error names the CMake option (fires on the WITH_JSON=OFF CI leg).
    EXPECT_EQ(mio_pipeline_has_json(), 0);
    EXPECT_NE(mio_pipeline_run_json("{}"), MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_JSON"), std::string::npos);
}

#endif  // MESHIOPLUSPLUS_HAS_JSON

// ---------------------------------------------------------------------------
// Sequences. The handle is a PLAN, not a cache: mio_sequence_read hands back an
// OWNED mesh precisely so nothing accumulates, which is the C ABI's expression
// of the streaming guarantee.
// ---------------------------------------------------------------------------

namespace {

/// Write `count` single-step .vtu files named `<dir>/in_<i>.vtu`.
std::string capi_seq_dir(int count) {
    const std::string dir = mt::temp_path("_capi_seq");
    std::filesystem::create_directories(dir);
    for (int i = 0; i < count; ++i) {
        mio_mesh* m = build_tet_mesh();
        mio_write((dir + "/in_" + std::to_string(i) + ".vtu").c_str(), m, "vtu");
        mio_mesh_free(m);
    }
    return dir;
}

}  // namespace

TEST(CApi, SequenceOpensAPatternInNaturalOrder) {
    const std::string dir = capi_seq_dir(12);
    mio_sequence* seq = mio_sequence_open((dir + "/in_*.vtu").c_str());
    ASSERT_NE(seq, nullptr) << mio_last_error();
    ASSERT_EQ(mio_sequence_count(seq), 12);

    // in_9 must precede in_10 -- the whole point of the ordering rule.
    char buf[512];
    ASSERT_GT(mio_sequence_path(seq, 9, buf, sizeof(buf)), 0);
    EXPECT_NE(std::string(buf).find("in_9.vtu"), std::string::npos);
    ASSERT_GT(mio_sequence_path(seq, 10, buf, sizeof(buf)), 0);
    EXPECT_NE(std::string(buf).find("in_10.vtu"), std::string::npos);

    double t = -1.0;
    ASSERT_EQ(mio_sequence_time(seq, 3, &t), MIO_OK);
    EXPECT_DOUBLE_EQ(t, 3.0);
    EXPECT_EQ(mio_sequence_time_source(seq, 3), 2);  // from the filename
    EXPECT_EQ(mio_sequence_step(seq, 3), 0);

    // The string-getter protocol: a short buffer reports the needed length.
    char tiny[2];
    EXPECT_GT(mio_sequence_path(seq, 0, tiny, sizeof(tiny)), 2);

    mio_sequence_free(seq);
    std::filesystem::remove_all(dir);
}

TEST(CApi, SequenceReadHandsBackAnOwnedMesh) {
    const std::string dir = capi_seq_dir(3);
    mio_sequence* seq = mio_sequence_open((dir + "/in_*.vtu").c_str());
    ASSERT_NE(seq, nullptr);
    for (int64_t i = 0; i < 3; ++i) {
        mio_mesh* m = mio_sequence_read(seq, i);
        ASSERT_NE(m, nullptr) << mio_last_error();
        EXPECT_EQ(mio_mesh_num_points(m), 5);
        // Owned: freeing it here must be correct, and the sequence must stay
        // usable afterwards (it holds no reference to what it produced).
        mio_mesh_free(m);
    }
    EXPECT_EQ(mio_sequence_count(seq), 3);
    mio_sequence_free(seq);
    std::filesystem::remove_all(dir);
}

TEST(CApi, SequenceExplicitListAndOptions) {
    const std::string dir = capi_seq_dir(3);
    const std::string a = dir + "/in_0.vtu";
    const std::string b = dir + "/in_1.vtu";
    const char* paths[2] = {b.c_str(), a.c_str()};  // a caller's stated order

    mio_sequence_opts opts;
    mio_sequence_opts_init(&opts);
    const double times[2] = {10.0, 20.0};
    opts.times = times;
    opts.num_times = 2;

    mio_sequence* seq = mio_sequence_open_list(paths, 2, &opts);
    ASSERT_NE(seq, nullptr) << mio_last_error();
    char buf[512];
    ASSERT_GT(mio_sequence_path(seq, 0, buf, sizeof(buf)), 0);
    EXPECT_NE(std::string(buf).find("in_1.vtu"), std::string::npos)
        << "an explicit list is a stated order and must not be re-sorted";
    double t = 0.0;
    ASSERT_EQ(mio_sequence_time(seq, 1, &t), MIO_OK);
    EXPECT_DOUBLE_EQ(t, 20.0);
    EXPECT_EQ(mio_sequence_time_source(seq, 1), 0);  // explicit
    mio_sequence_free(seq);

    // sort=1 reorders it.
    opts.sort = 1;
    seq = mio_sequence_open_list(paths, 2, &opts);
    ASSERT_NE(seq, nullptr);
    ASSERT_GT(mio_sequence_path(seq, 0, buf, sizeof(buf)), 0);
    EXPECT_NE(std::string(buf).find("in_0.vtu"), std::string::npos);
    mio_sequence_free(seq);
    std::filesystem::remove_all(dir);
}

TEST(CApi, SequenceFanInAndFanOut) {
    const std::string dir = capi_seq_dir(4);
    mio_sequence* seq = mio_sequence_open((dir + "/in_*.vtu").c_str());
    ASSERT_NE(seq, nullptr);

    const std::string series = dir + "/series.xdmf";
    ASSERT_EQ(mio_sequence_to_timeseries(seq, series.c_str(), nullptr), MIO_OK) << mio_last_error();
    mio_sequence_free(seq);
    ASSERT_TRUE(std::filesystem::exists(series));

    const std::string pattern = dir + "/back_{step}.vtu";
    ASSERT_EQ(mio_timeseries_to_sequence(series.c_str(), nullptr, pattern.c_str(), nullptr), MIO_OK)
        << mio_last_error();
    EXPECT_TRUE(std::filesystem::exists(dir + "/back_0003.vtu"));
    std::filesystem::remove_all(dir);
}

TEST(CApi, SequenceToTimeseriesExSelectsEncodingAndRejectsWhatItCannotHonour) {
    // The transient writer bypasses mio_write_ex's registry path entirely, so
    // it needs its own opt-in encoding selection -- and the same "reject an
    // option you cannot honour" rule (this is what a build without HDF5, e.g.
    // the Julia/R notebook environments, needs to select the XML data format).
    const std::string dir = capi_seq_dir(2);
    mio_sequence* seq = mio_sequence_open((dir + "/in_*.vtu").c_str());
    ASSERT_NE(seq, nullptr);

    mio_write_opts ascii;
    mio_write_opts_init(&ascii);
    ascii.encoding = MIO_ENCODING_ASCII;
    const std::string xml = dir + "/xml_series.xdmf";
    EXPECT_EQ(mio_sequence_to_timeseries_ex(seq, xml.c_str(), nullptr, &ascii), MIO_OK)
        << mio_last_error();
    EXPECT_TRUE(std::filesystem::exists(xml));
    // "XML" means the light data is inline -- no sibling .h5 companion.
    EXPECT_FALSE(std::filesystem::exists(dir + "/xml_series.h5"));

    mio_write_opts bad;
    mio_write_opts_init(&bad);
    bad.codec = MIO_CODEC_ZLIB;
    EXPECT_NE(mio_sequence_to_timeseries_ex(seq, (dir + "/s2.xdmf").c_str(), nullptr, &bad),
              MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("Codec"), std::string::npos);

    mio_sequence_free(seq);
    std::filesystem::remove_all(dir);
}

TEST(CApi, SequenceErrorsAreReportedNotCrashed) {
    EXPECT_EQ(mio_sequence_open(nullptr), nullptr);
    EXPECT_EQ(mio_sequence_count(nullptr), -1);
    EXPECT_EQ(mio_sequence_step(nullptr, 0), -1);
    EXPECT_EQ(mio_sequence_time_source(nullptr, 0), -1);
    EXPECT_EQ(mio_sequence_read(nullptr, 0), nullptr);
    EXPECT_NE(mio_sequence_time(nullptr, 0, nullptr), MIO_OK);
    EXPECT_EQ(mio_sequence_open_list(nullptr, 0, nullptr), nullptr);
    // A pattern matching nothing is an error, never an empty sequence.
    EXPECT_EQ(mio_sequence_open("/nonexistent-dir-for-meshio/out_*.vtu"), nullptr);
    mio_sequence_free(nullptr);  // must be a no-op

    const std::string dir = capi_seq_dir(2);
    mio_sequence* seq = mio_sequence_open((dir + "/in_*.vtu").c_str());
    ASSERT_NE(seq, nullptr);
    EXPECT_EQ(mio_sequence_read(seq, 99), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("out of range"), std::string::npos);
    // A target that cannot hold a series fails by name.
    EXPECT_NE(mio_sequence_to_timeseries(seq, (dir + "/no.vtu").c_str(), nullptr), MIO_OK);
    EXPECT_NE(std::string(mio_last_error()).find("{step}"), std::string::npos);
    mio_sequence_free(seq);
    std::filesystem::remove_all(dir);
}

// --- regular grids and signed distance --------------------------------------

namespace {

/// The unit cube [0,1]^3 as a closed, outward-wound triangle surface.
mio_mesh* capi_cube_surface() {
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> tris = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
                                            1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7};
    mio_mesh* m = mio_mesh_create();
    EXPECT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_add_cell_block(m, "triangle", 12, 3, MIO_INT64, tris.data()), MIO_OK);
    return m;
}

}  // namespace

TEST(CApi, GridBuildsALattice) {
    const std::int64_t dims[3] = {2, 3, 4};
    const double origin[3] = {-1.0, 0.0, 0.5};
    const double spacing[3] = {0.5, 2.0, 0.25};
    mio_mesh* g = mio_grid(dims, origin, spacing, 20000000);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(mio_mesh_num_points(g), 3 * 4 * 5);
    ASSERT_EQ(mio_mesh_num_cell_blocks(g), 1);
    EXPECT_EQ(block_type(g, 0), "hexahedron");
    std::int64_t num_cells = 0, npc = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(g, 0, &num_cells, &npc, nullptr), MIO_OK);
    EXPECT_EQ(num_cells, 24);
    EXPECT_EQ(npc, 8);
    mio_mesh_free(g);

    // NULL dims and an exceeded budget both fail rather than crash.
    EXPECT_EQ(mio_grid(nullptr, origin, spacing, 0), nullptr);
    const std::int64_t big[3] = {100, 100, 100};
    EXPECT_EQ(mio_grid(big, nullptr, nullptr, 1000), nullptr);
}

TEST(CApi, GridDefaultsOriginAndSpacingWhenNull) {
    const std::int64_t dims[3] = {1, 1, 1};
    mio_mesh* g = mio_grid(dims, nullptr, nullptr, 0);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(mio_mesh_num_points(g), 8);
    mio_mesh_free(g);
}

TEST(CApi, GenOptsInitDefaultsAndLayout) {
    mio_gen_opts opts;
    std::memset(&opts, 0xAB, sizeof opts);
    mio_gen_opts_init(&opts);
    EXPECT_EQ(sizeof(mio_gen_opts), 160u);
    EXPECT_EQ(opts.counts[0], 8);
    EXPECT_EQ(opts.size[2], 1.0);
    EXPECT_EQ(opts.radius, 1.0);
    EXPECT_EQ(opts.order, 1);
    EXPECT_EQ(opts.max_cells, 20000000);
    EXPECT_EQ(opts.pad, 0);
    EXPECT_EQ(opts.radius2, 0.5);
    EXPECT_EQ(opts.angle, 6.283185307179586);
    EXPECT_EQ(opts.num_twists, 1);
    EXPECT_EQ(opts.complete, 0);
    EXPECT_EQ(opts.ratio, 1.0);
    EXPECT_EQ(opts.seed, 0u);
    EXPECT_EQ(offsetof(mio_gen_opts, num_twists), 128u);
    EXPECT_EQ(offsetof(mio_gen_opts, seed), 144u);
    for (int64_t r : opts.reserved)
        EXPECT_EQ(r, 0);
    mio_gen_opts_init(nullptr);  // a NULL options pointer is a no-op, not a crash
}

namespace {

// Points, cells and the first block's type name of a generated mesh (freed).
struct CapiGen {
    int64_t mPoints = -1;
    int64_t mCells = 0;
    std::string mType;
};

CapiGen capi_gen(mio_mesh* pMesh) {
    CapiGen out;
    if (!pMesh)
        return out;
    out.mPoints = mio_mesh_num_points(pMesh);
    for (int64_t b = 0; b < mio_mesh_num_cell_blocks(pMesh); ++b) {
        mio_cell_block_info info{};
        EXPECT_EQ(mio_mesh_cell_block_info_ex(pMesh, b, &info), MIO_OK);
        out.mCells += info.num_cells;
    }
    char buf[64] = {0};
    if (mio_mesh_num_cell_blocks(pMesh) > 0)
        mio_mesh_cell_block_type(pMesh, 0, buf, sizeof buf);
    out.mType = buf;
    mio_mesh_free(pMesh);
    return out;
}

}  // namespace

TEST(CApi, SecondTierGeneratorsReadTheirFields) {
    mio_gen_opts o;
    mio_gen_opts_init(&o);
    o.counts[0] = 4;
    o.size[0] = 0.0;
    o.size[1] = 0.0;
    o.size[2] = 2.0;  // the end point
    CapiGen g = capi_gen(mio_line(&o));
    EXPECT_EQ(g.mPoints, 5);
    EXPECT_EQ(g.mType, "line");

    mio_gen_opts_init(&o);
    o.counts[0] = 6;
    o.order = 2;
    g = capi_gen(mio_circle(&o));
    EXPECT_EQ(g.mPoints, 12);
    EXPECT_EQ(g.mType, "line3");

    mio_gen_opts_init(&o);
    o.counts[0] = 10;
    o.angle = 0.0;  // two turns
    g = capi_gen(mio_helix(&o));
    EXPECT_EQ(g.mCells, 10);

    mio_gen_opts_init(&o);
    o.counts[0] = 7;
    o.seed = 9;
    EXPECT_EQ(capi_gen(mio_sphere_points(&o)).mType, "vertex");
    EXPECT_EQ(capi_gen(mio_box_points(&o)).mPoints, 7);
    o.order = 2;  // a point has no quadratic form
    EXPECT_EQ(mio_box_points(&o), nullptr);

    mio_gen_opts_init(&o);
    o.counts[0] = 1;
    EXPECT_EQ(capi_gen(mio_tetra_sphere(&o)).mCells, 16);
    EXPECT_EQ(capi_gen(mio_octa_sphere(&o)).mCells, 32);
    o.counts[0] = 2;
    EXPECT_EQ(capi_gen(mio_geo_sphere(&o)).mCells, 80);
    o.order = 2;
    o.complete = 1;
    g = capi_gen(mio_cube_sphere(&o));
    EXPECT_EQ(g.mCells, 24);
    EXPECT_EQ(g.mType, "quad9");

    mio_gen_opts_init(&o);
    o.counts[0] = 6;
    o.counts[1] = 4;
    o.radius2 = 0.0;  // the default minor radius, 0.25
    o.cell_type = MIO_GEN_TRIANGLE;
    EXPECT_EQ(capi_gen(mio_torus(&o)).mCells, 48);
    o.cell_type = MIO_GEN_DEFAULT;
    o.counts[1] = 2;
    o.num_twists = 3;
    EXPECT_EQ(capi_gen(mio_moebius(&o)).mCells, 12);
    o.counts[2] = 2;
    EXPECT_EQ(capi_gen(mio_tube(&o)).mCells, 12);
    EXPECT_EQ(capi_gen(mio_cone(&o)).mCells, 6 * 4);
    o.radius = 0.5;
    EXPECT_EQ(capi_gen(mio_capsule(&o)).mCells, 6 * (2 * 2 + 2 - 1 + 1));

    mio_gen_opts_init(&o);
    o.counts[0] = 1;
    o.counts[1] = 2;
    o.counts[2] = 3;
    EXPECT_EQ(capi_gen(mio_box_surface(&o)).mCells, 22);

    mio_gen_opts_init(&o);
    for (const auto& [faces, corners] :
         std::vector<std::pair<int64_t, int64_t>>{{4, 4}, {6, 8}, {8, 6}, {12, 20}, {20, 12}}) {
        o.counts[0] = faces;
        g = capi_gen(mio_platonic(&o));
        EXPECT_EQ(g.mCells, faces);
        EXPECT_EQ(g.mPoints, corners);
    }
    o.counts[0] = 7;  // not a platonic face count
    EXPECT_EQ(mio_platonic(&o), nullptr);

    mio_gen_opts_init(&o);
    o.counts[0] = 1;
    o.cell_type = MIO_GEN_TETRA;
    EXPECT_EQ(capi_gen(mio_ball(&o)).mCells, 48);
    EXPECT_EQ(capi_gen(mio_fichera(&o)).mCells, 7 * 6);
    o.cell_type = MIO_GEN_DEFAULT;
    o.counts[0] = 2;
    o.counts[1] = 1;
    o.radius2 = 0.0;  // the default inner radius, 0.5
    o.ratio = 0.0;    // a uniform wall
    EXPECT_EQ(capi_gen(mio_spherical_shell(&o)).mCells, 24);
    o.counts[0] = 6;
    o.counts[1] = 1;
    o.counts[2] = 2;
    EXPECT_EQ(capi_gen(mio_pipe(&o)).mCells, 12);
    o.counts[0] = 2;
    o.counts[1] = 5;
    o.order = 3;
    g = capi_gen(mio_solid_torus(&o));
    EXPECT_EQ(g.mCells, 20);
    EXPECT_EQ(g.mType, "VTK_LAGRANGE_HEXAHEDRON");
    o.cell_type = MIO_GEN_WEDGE;  // not built by this shape
    EXPECT_EQ(mio_solid_torus(&o), nullptr);
}

TEST(CApi, GeneratorsBuildMeshesAndRefuseBadRequests) {
    // NULL options mean mio_gen_opts_init's defaults: 8 x 8 x 8 hexahedra.
    mio_mesh* b = mio_box(nullptr);
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(mio_mesh_num_points(b), 9 * 9 * 9);
    EXPECT_EQ(block_type(b, 0), "hexahedron");
    mio_mesh_free(b);

    mio_gen_opts o;
    mio_gen_opts_init(&o);
    o.counts[0] = o.counts[1] = o.counts[2] = 2;
    o.cell_type = MIO_GEN_TETRA;
    o.order = 2;
    b = mio_box(&o);
    ASSERT_NE(b, nullptr);
    EXPECT_EQ(block_type(b, 0), "tetra10");
    mio_mesh_free(b);

    mio_gen_opts_init(&o);
    o.counts[0] = 3;
    o.counts[1] = 2;
    o.cell_type = MIO_GEN_TRIANGLE;
    o.split = MIO_GEN_SPLIT_CENTER;
    mio_mesh* r = mio_rectangle(&o);
    ASSERT_NE(r, nullptr);
    EXPECT_EQ(mio_mesh_num_points(r), 12 + 6);
    mio_mesh_free(r);

    mio_gen_opts_init(&o);
    o.counts[0] = 8;
    o.counts[1] = 4;
    o.cell_type = MIO_GEN_QUAD;
    mio_mesh* s = mio_uv_sphere(&o);
    ASSERT_NE(s, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(s), 2);
    mio_mesh_free(s);

    mio_gen_opts_init(&o);
    o.counts[0] = 1;
    mio_mesh* ico = mio_icosa_sphere(&o);
    ASSERT_NE(ico, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ico), 42);
    mio_mesh_free(ico);

    for (auto fn : {mio_disk, mio_cylinder}) {
        mio_gen_opts_init(&o);
        o.counts[0] = 2;
        o.counts[1] = 2;
        mio_mesh* m = fn(&o);
        ASSERT_NE(m, nullptr);
        mio_mesh_free(m);
    }

    // The planar domains. A struct filled by an older caller (radius2 and angle
    // both zero) reads as the defaults.
    mio_gen_opts_init(&o);
    o.counts[0] = 3;
    mio_mesh* tri = mio_reference_triangle(&o);
    ASSERT_NE(tri, nullptr);
    EXPECT_EQ(mio_mesh_num_points(tri), 10);
    mio_mesh_free(tri);

    mio_gen_opts_init(&o);
    o.counts[0] = 5;
    o.counts[1] = 2;
    mio_mesh* ng = mio_ngon(&o);
    ASSERT_NE(ng, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ng), 16);
    mio_mesh_free(ng);

    mio_gen_opts_init(&o);
    o.counts[0] = 12;
    o.counts[1] = 3;
    mio_mesh* ring = mio_annulus(&o);
    ASSERT_NE(ring, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ring), 48);  // a full ring welds its seam
    mio_mesh_free(ring);
    o.radius2 = 0.0;  // zeros from an older caller: the defaults
    o.angle = 0.0;
    ring = mio_annulus(&o);
    ASSERT_NE(ring, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ring), 48);
    mio_mesh_free(ring);
    mio_gen_opts_init(&o);
    o.counts[0] = 6;
    o.counts[1] = 2;
    o.angle = 1.5707963267948966;
    ring = mio_annulus(&o);
    ASSERT_NE(ring, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ring), 21);  // a sector keeps both end rays
    mio_mesh_free(ring);

    mio_gen_opts_init(&o);
    o.counts[0] = 3;
    o.size[0] = 2.0;
    mio_mesh* ls = mio_lshape(&o);
    ASSERT_NE(ls, nullptr);
    EXPECT_EQ(mio_mesh_num_points(ls), 40);
    mio_mesh_free(ls);

    mio_gen_opts_init(&o);
    o.counts[0] = 3;
    o.counts[1] = 2;
    o.radius = 0.25;
    mio_mesh* plate = mio_plate_with_hole(&o);
    ASSERT_NE(plate, nullptr);
    EXPECT_EQ(mio_mesh_num_points(plate), 36);
    mio_mesh_free(plate);
    mio_gen_opts_init(&o);  // radius 1 is not below half of the unit size
    EXPECT_EQ(mio_plate_with_hole(&o), nullptr);

    mio_gen_opts_init(&o);
    o.counts[0] = 3;
    o.counts[1] = 2;
    mio_mesh* hc = mio_honeycomb(&o);
    ASSERT_NE(hc, nullptr);
    EXPECT_EQ(mio_mesh_num_points(hc), 22);
    mio_mesh_free(hc);
    o.order = 2;  // a hexagon has no quadratic form
    EXPECT_EQ(mio_honeycomb(&o), nullptr);

    mio_mesh* ref = mio_reference_cells("tetra,wedge", 0);
    ASSERT_NE(ref, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(ref), 2);
    mio_mesh_free(ref);
    ref = mio_reference_cells(nullptr, 1);
    ASSERT_NE(ref, nullptr);
    EXPECT_EQ(mio_mesh_num_cell_blocks(ref), 7);
    mio_mesh_free(ref);

    // A cell type the shape does not build, a nonzero reserved tail, a bad
    // count, an exceeded budget and an unknown type all fail rather than crash.
    mio_gen_opts_init(&o);
    o.cell_type = MIO_GEN_WEDGE;
    EXPECT_EQ(mio_disk(&o), nullptr);
    mio_gen_opts_init(&o);
    o.reserved[0] = 1;
    EXPECT_EQ(mio_box(&o), nullptr);
    mio_gen_opts_init(&o);
    o.counts[0] = -1;
    EXPECT_EQ(mio_box(&o), nullptr);
    mio_gen_opts_init(&o);
    o.counts[0] = o.counts[1] = o.counts[2] = 100;
    o.max_cells = 1000;
    EXPECT_EQ(mio_box(&o), nullptr);
    EXPECT_EQ(mio_reference_cells("blob", 0), nullptr);
}

TEST(CApi, VoxelOptsInitDefaultsAndVoxelize) {
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);
    EXPECT_EQ(opts.fill, MIO_VOXEL_ALL);
    EXPECT_EQ(opts.sign, MIO_SDF_PSEUDONORMAL);
    EXPECT_EQ(opts.max_cells, 20000000);
    EXPECT_EQ(opts.resolution, nullptr);

    mio_mesh* cube = capi_cube_surface();
    const std::int64_t res[3] = {4, 4, 4};
    opts.resolution = res;

    std::int64_t dims[3] = {0, 0, 0};
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    std::int64_t occupied = -1;
    mio_mesh* g = mio_voxelize(cube, &opts, dims, origin, spacing, &occupied);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(occupied, 64);
    EXPECT_EQ(dims[0], 4);
    EXPECT_DOUBLE_EQ(spacing[0], 0.25);
    EXPECT_DOUBLE_EQ(origin[0], 0.0);
    mio_mesh_free(g);

    // Every out-param is nullable.
    mio_mesh* g2 = mio_voxelize(cube, &opts, nullptr, nullptr, nullptr, nullptr);
    ASSERT_NE(g2, nullptr);
    mio_mesh_free(g2);

    // NULL options is an error, not a default: a resolution must be chosen.
    EXPECT_EQ(mio_voxelize(cube, nullptr, nullptr, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, ComputeSdfOptsInitDefaultsAndVoxelStructure) {
    mio_compute_sdf_opts opts;
    mio_compute_sdf_opts_init(&opts);
    EXPECT_EQ(opts.structure, MIO_SDF_VOXEL);
    EXPECT_EQ(opts.root_resolution, 8);
    EXPECT_EQ(opts.max_depth, 4);
    EXPECT_DOUBLE_EQ(opts.band_cells, 1.0);
    EXPECT_DOUBLE_EQ(opts.padding_relative, 0.1);
    EXPECT_EQ(opts.max_cells, 20000000);
    EXPECT_EQ(opts.resolution, nullptr);
    // The embedded distance options are initialized too, or a caller who only
    // touched the outer struct would get sign = Unsigned by accident.
    EXPECT_EQ(opts.distance.sign, MIO_SDF_PSEUDONORMAL);
    EXPECT_DOUBLE_EQ(opts.distance.max_winding_work, 2.0e9);

    mio_mesh* cube = capi_cube_surface();
    const std::int64_t res[3] = {4, 4, 4};
    opts.resolution = res;
    opts.distance.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    std::int64_t dims[3] = {0, 0, 0}, depth = -1, banded = -1;
    double origin[3] = {0, 0, 0}, spacing[3] = {0, 0, 0};
    mio_surface_quality q{};
    mio_mesh* g = mio_compute_sdf(cube, &opts, dims, origin, spacing, &depth, &banded, &q);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(dims[0], 4);
    EXPECT_EQ(depth, 0);
    EXPECT_EQ(banded, 0);
    EXPECT_NE(q.watertight, 0);
    EXPECT_EQ(mio_mesh_num_point_data(g), 1);
    mio_mesh_free(g);

    // Every out-param is nullable.
    mio_mesh* g2 =
        mio_compute_sdf(cube, &opts, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
    ASSERT_NE(g2, nullptr);
    mio_mesh_free(g2);

    // NULL options is an error, not a default: a sizing must be chosen.
    EXPECT_EQ(mio_compute_sdf(cube, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr),
              nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, ComputeSdfOctreeRefinesAndRejectsAVoxelSizing) {
    mio_mesh* cube = capi_cube_surface();
    mio_compute_sdf_opts opts;
    mio_compute_sdf_opts_init(&opts);
    opts.structure = MIO_SDF_OCTREE;
    opts.root_resolution = 4;
    opts.max_depth = 2;
    opts.distance.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    std::int64_t depth = -1;
    mio_mesh* g = mio_compute_sdf(cube, &opts, nullptr, nullptr, nullptr, &depth, nullptr, nullptr);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(depth, 2);
    std::int64_t nc = 0, npc = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(g, 0, &nc, &npc, nullptr), MIO_OK);
    EXPECT_GT(nc, 64);
    EXPECT_LT(nc, 16 * 16 * 16);
    mio_mesh_free(g);

    // resolution/cell_size size a voxel grid; an octree's finest cell is already
    // determined, so accepting either would silently ignore one of the two.
    const std::int64_t res[3] = {8, 8, 8};
    opts.resolution = res;
    EXPECT_EQ(mio_compute_sdf(cube, &opts, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr),
              nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, CropPredicate) {
    mio_mesh* m = mio_mesh_create();
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0,
                                     0, 1, 0, 1, 1, 0, 2, 1, 0, 3, 1, 0};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    const std::vector<std::int64_t> conn = {0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6};
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 3, 4, MIO_INT64, conn.data()), MIO_OK);
    const std::vector<double> t = {0.0, 1.0, 2.0};
    const std::int64_t shape[1] = {3};
    ASSERT_EQ(mio_mesh_append_cell_data(m, "t", MIO_FLOAT64, 1, shape, t.data()), MIO_OK);

    mio_mesh* kept = mio_crop_predicate(m, "t", MIO_REFINE_LT, 1.5, 0);
    ASSERT_NE(kept, nullptr);
    std::int64_t nc = 0, npc = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(kept, 0, &nc, &npc, nullptr), MIO_OK);
    EXPECT_EQ(nc, 2);
    EXPECT_EQ(mio_mesh_num_points(kept), 6);
    mio_mesh_free(kept);

    // An unknown array and an out-of-range comparison both fail rather than
    // defaulting.
    EXPECT_EQ(mio_crop_predicate(m, "nope", MIO_REFINE_LT, 0.0, 0), nullptr);
    EXPECT_EQ(mio_crop_predicate(m, "t", 99, 0.0, 0), nullptr);
    EXPECT_EQ(mio_crop_predicate(nullptr, "t", MIO_REFINE_LT, 0.0, 0), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, VoxelizeInsideKeepsTheInterior) {
    mio_mesh* cube = capi_cube_surface();
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);
    const std::int64_t res[3] = {5, 5, 5};
    const double bounds[6] = {-0.5, -0.5, -0.5, 1.5, 1.5, 1.5};
    opts.resolution = res;
    opts.bounds = bounds;
    opts.fill = MIO_VOXEL_INSIDE;
    opts.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    std::int64_t occupied = -1;
    mio_mesh* g = mio_voxelize(cube, &opts, nullptr, nullptr, nullptr, &occupied);
    ASSERT_NE(g, nullptr);
    EXPECT_EQ(occupied, 27);
    mio_mesh_free(g);
    mio_mesh_free(cube);
}

TEST(CApi, OctreeAroundASurfaceRefinesOnlyWhereItCuts) {
    mio_mesh* cube = capi_cube_surface();
    mio_octree_opts opts;
    mio_octree_opts_init(&opts);
    const double bounds[6] = {-0.5, -0.5, -0.5, 1.5, 1.5, 1.5};
    opts.bounds = bounds;
    opts.max_level = 3;
    opts.fill = MIO_VOXEL_INSIDE;
    opts.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    mio_octree_report report;
    mio_mesh* tree = mio_octree(cube, &opts, &report);
    ASSERT_NE(tree, nullptr);
    EXPECT_EQ(report.dimension, 3);
    EXPECT_GT(report.num_leaves, 0);
    EXPECT_EQ(report.finest_level, 3);
    mio_mesh_free(tree);
    mio_mesh_free(cube);
}

TEST(CApi, OctreeConformingTransitionsHaveNoHangingNodes) {
    mio_mesh* cube = capi_cube_surface();
    for (int transition : {MIO_OCTREE_POLYHEDRA, MIO_OCTREE_TEMPLATES}) {
        mio_octree_opts opts;
        mio_octree_opts_init(&opts);
        const double bounds[6] = {-0.5, -0.5, -0.5, 1.5, 1.5, 1.5};
        opts.bounds = bounds;
        opts.max_level = 3;
        opts.transition = transition;
        opts.watertight_check = MIO_SDF_WATERTIGHT_OFF;
        mio_octree_report report;
        mio_mesh* tree = mio_octree(cube, &opts, &report);
        ASSERT_NE(tree, nullptr);
        EXPECT_EQ(report.num_hanging, 0);
        mio_mesh_free(tree);
    }
    mio_octree_opts bad;
    mio_octree_opts_init(&bad);
    bad.transition = 99;
    EXPECT_EQ(mio_octree(cube, &bad, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, OctreeMeshTetrahedraAndCastellatedCells) {
    mio_mesh* cube = capi_cube_surface();
    mio_octree_mesh_opts opts;
    mio_octree_mesh_opts_init(&opts);
    EXPECT_EQ(opts.output, MIO_OCTREE_MESH_TET);
    EXPECT_EQ(opts.preserve_features, 1);
    EXPECT_DOUBLE_EQ(opts.alpha_long, 0.24799);
    opts.max_size = 0.5;
    opts.min_size = 0.125;
    mio_octree_mesh_report report;
    mio_mesh* tets = mio_octree_mesh(cube, &opts, &report);
    ASSERT_NE(tets, nullptr) << mio_last_error();
    EXPECT_EQ(report.output, MIO_OCTREE_MESH_TET);
    EXPECT_GT(report.num_tetrahedra, 0);
    EXPECT_EQ(report.num_inverted, 0);
    EXPECT_EQ(report.num_corners, 8);
    EXPECT_EQ(report.num_corners_snapped, 8);
    EXPECT_EQ(report.num_boundary_edges_open, 0);
    std::int64_t nc = 0, npc = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(tets, 0, &nc, &npc, nullptr), MIO_OK);
    EXPECT_EQ(nc, report.num_tetrahedra);
    EXPECT_EQ(npc, 4);
    // The warnings accessor follows string rule 5.
    const std::int64_t need = mio_octree_mesh_warnings(nullptr, 0);
    EXPECT_GE(need, 0);
    EXPECT_EQ(need == 0, report.num_warnings == 0);
    mio_mesh_free(tets);

    opts.output = MIO_OCTREE_MESH_CASTELLATED;
    mio_mesh* cells = mio_octree_mesh(cube, &opts, &report);
    ASSERT_NE(cells, nullptr) << mio_last_error();
    EXPECT_GT(report.num_hexahedra, 0);
    EXPECT_GT(report.num_surrogate_faces, 0);
    mio_mesh_free(cells);

    opts.output = MIO_OCTREE_MESH_HEX;
    opts.surface_size = 0.25;
    mio_mesh* hexes = mio_octree_mesh(cube, &opts, &report);
    ASSERT_NE(hexes, nullptr) << mio_last_error();
    EXPECT_EQ(report.output, MIO_OCTREE_MESH_HEX);
    EXPECT_GT(report.num_hexahedra, 0);
    EXPECT_EQ(report.num_tetrahedra, 0);
    EXPECT_EQ(report.num_inverted, 0);
    EXPECT_EQ(report.num_boundary_edges_open, 0);
    mio_mesh_free(hexes);

    EXPECT_EQ(opts.hex_grading, MIO_OCTREE_MESH_HEX_UNIFORM);
    opts.hex_grading = MIO_OCTREE_MESH_HEX_GRADED;
    mio_mesh* graded = mio_octree_mesh(cube, &opts, &report);
    ASSERT_NE(graded, nullptr) << mio_last_error();
    EXPECT_GT(report.num_hexahedra, 0);
    EXPECT_EQ(report.num_inverted, 0);
    EXPECT_EQ(report.num_boundary_edges_open, 0);
    mio_mesh_free(graded);
    opts.hex_grading = 99;
    EXPECT_EQ(mio_octree_mesh(cube, &opts, nullptr), nullptr);
    opts.hex_grading = MIO_OCTREE_MESH_HEX_UNIFORM;

    // Out-of-range enums fail.
    opts.output = 99;
    EXPECT_EQ(mio_octree_mesh(cube, &opts, nullptr), nullptr);
    EXPECT_EQ(mio_octree_mesh(nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, OctreeRejectsAnUnknownBalance) {
    mio_mesh* cube = capi_cube_surface();
    mio_octree_opts opts;
    mio_octree_opts_init(&opts);
    opts.balance = 99;  // drift guard: an out-of-range enum fails, never reinterprets
    EXPECT_EQ(mio_octree(cube, &opts, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, VoxelizeRejectsAnUnknownFill) {
    mio_mesh* cube = capi_cube_surface();
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);
    const std::int64_t res[3] = {2, 2, 2};
    opts.resolution = res;
    opts.fill = 99;  // drift guard: an out-of-range enum fails, never reinterprets
    EXPECT_EQ(mio_voxelize(cube, &opts, nullptr, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, VoxelizeTakesAPlane) {
    // The cube surface's boundary-free skin in 2-D: `fill` all needs only its
    // box, and `plane` chooses which two axes it spans.
    mio_mesh* cube = capi_cube_surface();
    mio_voxel_opts opts;
    mio_voxel_opts_init(&opts);
    EXPECT_EQ(opts.plane, MIO_PLANE_AUTO);
    const std::int64_t res[3] = {4, 2, 0};
    opts.resolution = res;
    opts.plane = MIO_PLANE_YZ;
    std::int64_t dims[3] = {0, 0, 0};
    mio_mesh* out = mio_voxelize(cube, &opts, dims, nullptr, nullptr, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(dims[0], 4);
    EXPECT_EQ(dims[1], 2);
    EXPECT_EQ(dims[2], 1);
    mio_mesh_free(out);
    opts.plane = 7;  // drift guard, as for fill
    EXPECT_EQ(mio_voxelize(cube, &opts, nullptr, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, SurfaceWatertightCheckCountsDefects) {
    mio_mesh* cube = capi_cube_surface();
    mio_surface_quality q;
    ASSERT_EQ(mio_surface_watertight_check(cube, &q), MIO_OK);
    EXPECT_NE(q.watertight, 0);
    EXPECT_EQ(q.boundary_edges, 0);
    EXPECT_EQ(q.degenerate_triangles, 0);
    mio_mesh_free(cube);

    // A lone triangle is a sheet: three boundary edges, reported as a number.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    const std::vector<std::int64_t> tri = {0, 1, 2};
    mio_mesh* sheet = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(sheet, MIO_FLOAT64, 3, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(sheet, "triangle", 1, 3, MIO_INT64, tri.data()), MIO_OK);
    ASSERT_EQ(mio_surface_watertight_check(sheet, &q), MIO_OK);
    EXPECT_EQ(q.watertight, 0);
    EXPECT_EQ(q.boundary_edges, 3);
    mio_mesh_free(sheet);

    EXPECT_EQ(mio_surface_watertight_check(nullptr, &q), MIO_ERR_INVALID_ARG);
}

TEST(CApi, SampleDistanceWritesIntoTheCallersBuffer) {
    mio_mesh* cube = capi_cube_surface();
    mio_sdf_opts opts;
    mio_sdf_opts_init(&opts);
    opts.watertight_check = MIO_SDF_WATERTIGHT_OFF;

    const double points[9] = {0.5, 0.5, 0.5, 2.0, 0.5, 0.5, -1.0, 0.5, 0.5};
    double out[3] = {0, 0, 0};
    ASSERT_EQ(mio_sample_distance(cube, points, 3, &opts, out), MIO_OK);
    EXPECT_NEAR(out[0], -0.5, 1e-12);  // the centre, inside
    EXPECT_NEAR(out[1], 1.0, 1e-12);
    EXPECT_NEAR(out[2], 1.0, 1e-12);

    // NULL options means the defaults, not a failure.
    ASSERT_EQ(mio_sample_distance(cube, points, 1, nullptr, out), MIO_OK);
    EXPECT_NEAR(out[0], -0.5, 1e-12);

    // Zero queries is a no-op, not an error; NULL buffers are refused.
    EXPECT_EQ(mio_sample_distance(cube, points, 0, &opts, out), MIO_OK);
    EXPECT_EQ(mio_sample_distance(cube, nullptr, 1, &opts, out), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_sample_distance(cube, points, 1, &opts, nullptr), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_sample_distance(nullptr, points, 1, &opts, out), MIO_ERR_INVALID_ARG);
    mio_mesh_free(cube);
}

TEST(CApi, DistanceToSurfaceAttachesPointData) {
    mio_mesh* cube = capi_cube_surface();
    const std::int64_t dims[3] = {2, 2, 2};
    const double origin[3] = {-0.5, -0.5, -0.5};
    const double spacing[3] = {1.0, 1.0, 1.0};
    mio_mesh* q = mio_grid(dims, origin, spacing, 0);
    ASSERT_NE(q, nullptr);

    mio_sdf_opts opts;
    mio_sdf_opts_init(&opts);
    opts.watertight_check = MIO_SDF_WATERTIGHT_OFF;
    opts.record_inside = 1;

    std::int64_t banded = -1;
    mio_surface_quality quality;
    mio_mesh* out = mio_distance_to_surface(q, cube, &opts, &banded, &quality);
    ASSERT_NE(out, nullptr);
    EXPECT_EQ(banded, 0);  // no band in force
    EXPECT_NE(quality.watertight, 0);
    EXPECT_EQ(mio_mesh_num_point_data(out), 2);  // sdf:distance and sdf:inside
    mio_mesh_free(out);

    EXPECT_EQ(mio_distance_to_surface(nullptr, cube, &opts, nullptr, nullptr), nullptr);
    mio_mesh_free(q);
    mio_mesh_free(cube);
}

// ---------------------------------------------------------------------------
// Provenance (v10.16.0): the C ABI bridge over detail::ProvenanceScope.
// ---------------------------------------------------------------------------

TEST(CApi, ProvenanceOffByDefaultIsByteIdenticalToNoScope) {
    mio_mesh* m = capi_cube_surface();
    const std::string a = mt::temp_path("_capi_prov_a.obj");
    const std::string b = mt::temp_path("_capi_prov_b.obj");
    ASSERT_EQ(mio_write(a.c_str(), m, "obj"), MIO_OK);

    ASSERT_EQ(mio_provenance_scope_begin(MIO_PROVENANCE_OFF), MIO_OK);
    ASSERT_EQ(mio_write(b.c_str(), m, "obj"), MIO_OK);
    mio_provenance_scope_end();

    std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
    std::string sa((std::istreambuf_iterator<char>(fa)), std::istreambuf_iterator<char>());
    std::string sb((std::istreambuf_iterator<char>(fb)), std::istreambuf_iterator<char>());
    EXPECT_EQ(sa, sb);

    mio_mesh_free(m);
    std::error_code ec;
    std::filesystem::remove(a, ec);
    std::filesystem::remove(b, ec);
}

TEST(CApi, ProvenanceBestEffortRecordsSourceAndNotes) {
    mio_mesh* m = capi_cube_surface();
    const std::string path = mt::temp_path("_capi_prov.obj");

    ASSERT_EQ(mio_provenance_scope_begin(MIO_PROVENANCE_BEST_EFFORT), MIO_OK);
    mio_provenance_set_source("in.vtu", "vtu");
    mio_provenance_note("regions-dropped", "1 named region(s) dropped");
    ASSERT_EQ(mio_write(path.c_str(), m, "obj"), MIO_OK);
    mio_provenance_scope_end();

    std::ifstream f(path);
    std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("in.vtu"), std::string::npos);
    EXPECT_NE(text.find("regions-dropped"), std::string::npos);

    // Scope closed: the next write must not see the stale note.
    const std::string path2 = mt::temp_path("_capi_prov2.obj");
    ASSERT_EQ(mio_write(path2.c_str(), m, "obj"), MIO_OK);
    std::ifstream f2(path2);
    std::string text2((std::istreambuf_iterator<char>(f2)), std::istreambuf_iterator<char>());
    EXPECT_EQ(text2.find("in.vtu"), std::string::npos);

    mio_mesh_free(m);
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(path2, ec);
}

TEST(CApi, ProvenanceScopeBeginRejectsAnUnknownMode) {
    EXPECT_EQ(mio_provenance_scope_begin(99), MIO_ERR_INVALID_ARG);
    EXPECT_EQ(mio_provenance_scope_begin(-1), MIO_ERR_INVALID_ARG);
}

TEST(CApi, ProvenanceScopeEndWithNoneOpenIsANoOp) {
    mio_provenance_scope_end();  // must not crash
    SUCCEED();
}

TEST(CApi, ProvenanceReadBackThroughMetadata) {
    mio_mesh* m = capi_cube_surface();
    const std::string path = mt::temp_path("_capi_rb.obj");

    ASSERT_EQ(mio_provenance_scope_begin(MIO_PROVENANCE_BEST_EFFORT), MIO_OK);
    mio_provenance_set_source("in.msh", "gmsh");
    mio_provenance_note("regions-dropped", "Side regions dropped");
    ASSERT_EQ(mio_write(path.c_str(), m, "obj"), MIO_OK);
    mio_provenance_scope_end();

    mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), "obj");
    ASSERT_NE(meta, nullptr);
    EXPECT_EQ(mio_read_metadata_provenance_recognised(meta), 1);
    const int64_t n = mio_read_metadata_num_provenance_lines(meta);
    ASSERT_GE(n, 3);

    bool saw_source = false, saw_note = false;
    char buf[512];
    for (int64_t i = 0; i < n; ++i) {
        ASSERT_GT(mio_read_metadata_provenance_line(meta, i, buf, sizeof(buf)), 0);
        std::string line(buf);
        if (line.find("in.msh") != std::string::npos)
            saw_source = true;
        if (line.find("regions-dropped") != std::string::npos)
            saw_note = true;
    }
    EXPECT_TRUE(saw_source);
    EXPECT_TRUE(saw_note);

    // Out-of-range and NULL are errors, not silent zeros.
    EXPECT_EQ(mio_read_metadata_provenance_line(meta, n, buf, sizeof(buf)), -1);
    EXPECT_EQ(mio_read_metadata_num_provenance_lines(nullptr), -1);
    EXPECT_EQ(mio_read_metadata_provenance_recognised(nullptr), -1);

    mio_read_metadata_free(meta);
    mio_mesh_free(m);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(CApi, ProvenanceReadBackIsHonestAboutAForeignFile) {
    const std::string path = mt::temp_path("_capi_foreign.obj");
    {
        std::ofstream f(path);
        f << "# Created by SomeOtherTool 3.2\nv 0 0 0\nv 1 0 0\nv 1 1 0\nf 1 2 3\n";
    }
    mio_read_metadata* meta = mio_read_metadata_create(path.c_str(), "obj");
    ASSERT_NE(meta, nullptr);
    EXPECT_EQ(mio_read_metadata_provenance_recognised(meta), 0);
    EXPECT_EQ(mio_read_metadata_num_provenance_lines(meta), 0);
    mio_read_metadata_free(meta);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

namespace {

/// A regular octahedron: the smallest closed, consistently wound triangle
/// surface, so Gauss-Bonnet applies and every edge is used exactly twice.
mio_mesh* capi_curv_octahedron() {
    const std::vector<double> pts = {1.0, 0.0,  0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                     0.0, -1.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0, -1.0};
    const std::vector<std::int64_t> conn = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
                                            2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
    mio_mesh* m = mio_mesh_create();
    EXPECT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 6, 3, pts.data()), MIO_OK);
    EXPECT_EQ(mio_mesh_add_cell_block(m, "triangle", 8, 3, MIO_INT64, conn.data()), MIO_OK);
    return m;
}

}  // namespace

TEST(CApi, RepairRewindsAndClosesThroughTheFlatAbi) {
    // The octahedron with two facets flipped: repair must rewind exactly
    // those two, report a watertight output, and leave the point count alone.
    // Then the same octahedron minus one facet: the hole is filled by three
    // triangles in a trailing block.
    const std::vector<double> pts = {1.0, 0.0,  0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                     0.0, -1.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0, -1.0};
    std::vector<std::int64_t> conn = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
                                      2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
    std::swap(conn[1], conn[2]);    // flip facet 0
    std::swap(conn[13], conn[14]);  // flip facet 4
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 6, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 8, 3, MIO_INT64, conn.data()), MIO_OK);

    mio_repair_opts opts;
    mio_repair_opts_init(&opts);
    EXPECT_EQ(opts.fix_orientation, 1);  // ON by default -- all-zero is NOT the default
    EXPECT_EQ(opts.fill_holes, 1);
    EXPECT_EQ(opts.max_hole_edges, 10);
    mio_repair_report report;
    mio_mesh* out = mio_repair(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.quality_before.inconsistent_pairs, 4);  // 0 and 4 share edge (0,2)
    EXPECT_EQ(report.quality_after.inconsistent_pairs, 0);
    EXPECT_NE(report.quality_after.watertight, 0);
    EXPECT_EQ(report.num_flipped, 2);
    EXPECT_EQ(report.num_components, 1);
    EXPECT_EQ(report.num_holes_detected, 0);
    EXPECT_EQ(mio_mesh_num_points(out), 6);
    EXPECT_EQ(mio_mesh_num_cell_blocks(out), 1);
    mio_mesh_free(out);
    mio_mesh_free(m);

    mio_mesh* holed = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(holed, MIO_FLOAT64, 6, 3, pts.data()), MIO_OK);
    const std::vector<std::int64_t> seven = {2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0,
                                             5, 1, 2, 5, 3, 1, 5, 0, 3, 5};
    ASSERT_EQ(mio_mesh_add_cell_block(holed, "triangle", 7, 3, MIO_INT64, seven.data()), MIO_OK);
    mio_mesh* filled = mio_repair(holed, nullptr, &report);
    ASSERT_NE(filled, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_holes_filled, 1);
    EXPECT_EQ(report.num_faces_added, 3);
    EXPECT_EQ(report.num_points_added, 1);
    EXPECT_NE(report.quality_after.watertight, 0);
    EXPECT_EQ(mio_mesh_num_points(filled), 7);
    EXPECT_EQ(mio_mesh_num_cell_blocks(filled), 2);
    mio_mesh_free(filled);
    mio_mesh_free(holed);
}

TEST(CApi, ShrinkwrapProjectsOntoAPlaneWithAnOffset) {
    // A z = 0 square as the target; three points above and below it land at
    // exactly z = offset. Weights by array name select the second one out.
    const std::vector<double> tpts = {0, 0, 0, 4, 0, 0, 4, 4, 0, 0, 4, 0};
    const std::vector<std::int64_t> tconn = {0, 1, 2, 0, 2, 3};
    mio_mesh* target = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(target, MIO_FLOAT64, 4, 3, tpts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(target, "triangle", 2, 3, MIO_INT64, tconn.data()), MIO_OK);

    const std::vector<double> pts = {1, 1, 0.5, 2, 2, -0.7, 3, 1, 1.5};
    const std::vector<std::int64_t> vconn = {0, 1, 2};
    const std::vector<std::int32_t> pick = {1, 0, 1};
    const int64_t shape[1] = {3};
    mio_mesh* src = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(src, MIO_FLOAT64, 3, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(src, "vertex", 3, 1, MIO_INT64, vconn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(src, "pick", MIO_INT32, 1, shape, pick.data()), MIO_OK);

    mio_shrinkwrap_opts opts;
    mio_shrinkwrap_opts_init(&opts);
    EXPECT_EQ(opts.normal_weight, MIO_SDF_WEIGHT_ANGLE);
    opts.offset = 0.25;
    opts.weights = "pick";
    opts.record_distance = 1;
    mio_shrinkwrap_report report;
    mio_mesh* out = mio_shrinkwrap(src, target, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_projected, 2);
    EXPECT_EQ(report.num_skipped, 1);
    EXPECT_EQ(report.num_missed, 0);
    EXPECT_NEAR(report.max_displacement, 1.25, 1e-15);
    const void* data = nullptr;
    mio_dtype dt;
    ASSERT_EQ(mio_mesh_get_points(out, &data, &dt), MIO_OK);
    const double* p = static_cast<const double*>(data);
    EXPECT_NEAR(p[2], 0.25, 1e-15);
    EXPECT_DOUBLE_EQ(p[5], -0.7);  // unselected: untouched
    EXPECT_NEAR(p[8], 0.25, 1e-15);
    mio_mesh_free(out);

    // A volume target is refused by name.
    opts.weights = nullptr;
    const std::vector<double> vpts = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    const std::vector<std::int64_t> tet = {0, 1, 2, 3};
    mio_mesh* vol = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(vol, MIO_FLOAT64, 4, 3, vpts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(vol, "tetra", 1, 4, MIO_INT64, tet.data()), MIO_OK);
    EXPECT_EQ(mio_shrinkwrap(src, vol, &opts, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("extract_surface"), std::string::npos);
    mio_mesh_free(vol);
    mio_mesh_free(src);
    mio_mesh_free(target);
}

TEST(CApi, SobolevDeformFiltersADisplacementAndPinsByArray) {
    // A 3x3 triangle grid with a checkerboard z-displacement: the filter
    // damps it, a pinned corner does not move, and a constant field with
    // nothing pinned is preserved exactly with zero iterations.
    std::vector<double> pts;
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i) {
            pts.push_back(i);
            pts.push_back(j);
            pts.push_back(0.0);
        }
    std::vector<std::int64_t> conn;
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i) {
            const std::int64_t a = j * 4 + i, b = a + 1, c = a + 4, d = c + 1;
            conn.insert(conn.end(), {a, b, d, a, d, c});
        }
    std::vector<double> disp(16 * 3, 0.0);
    for (int k = 0; k < 16; ++k)
        disp[k * 3 + 2] = ((k % 4 + k / 4) % 2) ? 0.1 : -0.1;
    std::vector<std::int32_t> pin(16, 0);
    pin[0] = 1;
    const int64_t dshape[2] = {16, 3};
    const int64_t pshape[1] = {16};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 16, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "triangle", 18, 3, MIO_INT64, conn.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "d", MIO_FLOAT64, 2, dshape, disp.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(m, "pin", MIO_INT32, 1, pshape, pin.data()), MIO_OK);

    mio_sobolev_opts opts;
    mio_sobolev_opts_init(&opts);
    EXPECT_EQ(opts.max_iterations, 128);
    EXPECT_DOUBLE_EQ(opts.tolerance, 1e-10);
    // The array name is required.
    EXPECT_EQ(mio_sobolev_deform(m, &opts, nullptr), nullptr);
    opts.array = "d";
    opts.length_scale = 2.0;
    opts.fixed_points_array = "pin";
    opts.record_filtered = 1;
    mio_sobolev_report report;
    mio_mesh* out = mio_sobolev_deform(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_NE(report.converged, 0);
    EXPECT_GT(report.num_iterations, 0);
    EXPECT_EQ(report.num_fixed, 1);
    EXPECT_EQ(report.num_isolated, 0);
    EXPECT_LT(report.residual, 1e-10);
    EXPECT_GT(report.max_displacement, 0.0);
    EXPECT_LT(report.max_displacement, 0.1);  // damped
    const void* data = nullptr;
    mio_dtype dt;
    ASSERT_EQ(mio_mesh_get_points(out, &data, &dt), MIO_OK);
    const double* p = static_cast<const double*>(data);
    EXPECT_DOUBLE_EQ(p[2], 0.0);                 // pinned
    EXPECT_EQ(mio_mesh_num_point_data(out), 3);  // d, pin, sobolev:displacement
    mio_mesh_free(out);

    // A constant field with nothing pinned: exact, zero iterations.
    std::vector<double> cst(16 * 3, 0.0);
    for (int k = 0; k < 16; ++k)
        cst[k * 3 + 1] = 0.4;
    ASSERT_EQ(mio_mesh_add_point_data(m, "c", MIO_FLOAT64, 2, dshape, cst.data()), MIO_OK);
    opts.array = "c";
    opts.fixed_points_array = nullptr;
    mio_mesh* moved = mio_sobolev_deform(m, &opts, &report);
    ASSERT_NE(moved, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_iterations, 0);
    EXPECT_EQ(report.num_fixed, 0);
    ASSERT_EQ(mio_mesh_get_points(moved, &data, &dt), MIO_OK);
    p = static_cast<const double*>(data);
    for (int k = 0; k < 16; ++k)
        EXPECT_DOUBLE_EQ(p[k * 3 + 1], pts[k * 3 + 1] + 0.4);
    mio_mesh_free(moved);
    mio_mesh_free(m);
}

TEST(CApi, ComputeNormalsSplitsACubeAtItsCreases) {
    // A unit cube of six outward-wound quads: 8 points smooth, 24 once split at
    // a 30 degree crease, with every corner of the original on three faces.
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                                     0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    const std::vector<std::int64_t> quads = {0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4,
                                             3, 7, 6, 2, 0, 4, 7, 3, 1, 2, 6, 5};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 8, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 6, 4, MIO_INT64, quads.data()), MIO_OK);

    mio_normals_opts opts;
    mio_normals_opts_init(&opts);
    EXPECT_EQ(opts.point_normals, 1);  // ON by default -- an all-zero struct is NOT the default
    EXPECT_EQ(opts.split, 0);
    EXPECT_DOUBLE_EQ(opts.split_angle, 30.0);
    EXPECT_EQ(opts.weight, MIO_SDF_WEIGHT_ANGLE);

    mio_normals_report report;
    mio_mesh* smooth = mio_compute_normals(m, &opts, &report);
    ASSERT_NE(smooth, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(smooth), 8);
    EXPECT_EQ(report.num_added_points, 0);
    EXPECT_NE(report.quality.watertight, 0);
    mio_mesh_free(smooth);

    opts.split = 1;
    opts.cell_normals = 1;
    opts.record_parent_ids = 1;
    mio_mesh* split = mio_compute_normals(m, &opts, &report);
    ASSERT_NE(split, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(split), 24);
    EXPECT_EQ(report.num_added_points, 16);
    EXPECT_EQ(report.num_split_points, 8);
    EXPECT_EQ(report.num_isolated, 0);
    mio_mesh_free(split);

    // NULL options means every default, and a NULL report is not written.
    mio_mesh* plain = mio_compute_normals(m, nullptr, nullptr);
    ASSERT_NE(plain, nullptr) << mio_last_error();
    mio_mesh_free(plain);

    // An out-of-range weight or split angle is refused, never clamped, and no
    // exception crosses the ABI.
    opts.weight = 9;
    EXPECT_EQ(mio_compute_normals(m, &opts, nullptr), nullptr);
    opts.weight = MIO_SDF_WEIGHT_AREA;
    opts.split_angle = 270.0;
    EXPECT_EQ(mio_compute_normals(m, &opts, nullptr), nullptr);
    EXPECT_EQ(mio_compute_normals(nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, ComputeCurvatureSatisfiesGaussBonnet) {
    // The tessellation-independent oracle: on a closed surface the angle
    // defects sum to 2*pi*chi, which is 4*pi for anything sphere-like. It is
    // reported through the C ABI precisely so a flat-binding consumer can check
    // a result without reimplementing the estimator.
    mio_mesh* m = capi_curv_octahedron();

    mio_curvature_opts opts;
    mio_curvature_opts_init(&opts);
    EXPECT_EQ(opts.mean, 1);  // ON by default -- an all-zero struct is NOT the default
    EXPECT_EQ(opts.gaussian, 1);
    EXPECT_EQ(opts.dual_area, MIO_CURVATURE_MIXED_VORONOI);
    opts.record_area = 1;
    opts.record_principal = 1;

    mio_curvature_report report;
    mio_mesh* out = mio_compute_curvature(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_NEAR(report.total_angle_defect, 4.0 * 3.14159265358979323846, 1e-12);
    EXPECT_EQ(report.num_boundary, 0);
    EXPECT_EQ(report.num_isolated, 0);
    EXPECT_EQ(report.num_degenerate, 0);
    EXPECT_NE(report.quality.watertight, 0);
    EXPECT_EQ(report.quality.inconsistent_pairs, 0);
    EXPECT_EQ(mio_mesh_num_points(out), 6);  // a pure data step
    mio_mesh_free(out);

    // NULL options means every mio_curvature_opts_init default, and a NULL
    // report is simply not written.
    mio_mesh* plain = mio_compute_curvature(m, nullptr, nullptr);
    ASSERT_NE(plain, nullptr) << mio_last_error();
    mio_mesh_free(plain);

    // Barycentric gives the same defect: it changes the dual area, and the
    // defect never touches it.
    opts.dual_area = MIO_CURVATURE_BARYCENTRIC;
    mio_curvature_report bary;
    mio_mesh* out2 = mio_compute_curvature(m, &opts, &bary);
    ASSERT_NE(out2, nullptr) << mio_last_error();
    EXPECT_EQ(bary.total_angle_defect, report.total_angle_defect);
    mio_mesh_free(out2);

    // An out-of-range dual area is refused rather than silently clamped.
    opts.dual_area = 7;
    EXPECT_EQ(mio_compute_curvature(m, &opts, nullptr), nullptr);
    mio_mesh_free(m);

    // A volume mesh is refused by name; no exception crosses the ABI.
    const std::vector<double> tp = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    const std::vector<std::int64_t> tc = {0, 1, 2, 3};
    mio_mesh* vol = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(vol, MIO_FLOAT64, 4, 3, tp.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(vol, "tetra", 1, 4, MIO_INT64, tc.data()), MIO_OK);
    EXPECT_EQ(mio_compute_curvature(vol, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("extract_surface"), std::string::npos);
    mio_mesh_free(vol);
}

TEST(CApi, WriteToADatNeverLetsTheOldFileChoose) {
    // A Marc deck at the path reads as Marc, but mio_write resolves the write
    // from the extension alone: Tecplot's (until v16.17.0: "no writer for
    // format 'marc'").
    const std::string dat = std::string(::testing::TempDir()) + "capi_over_marc.dat";
    {
        std::ofstream out(dat);
        out << "title x\nend\nconnectivity\n    1    0    1\n    1  134    1    2    3    4\n"
               "coordinates\n    3    4    0    1\n1, 0., 0., 0.\n2, 1., 0., 0.\n3, 0., 1., 0.\n"
               "4, 0., 0., 1.\nend option\n";
    }
    mio_mesh* deck = mio_read(dat.c_str(), nullptr);  // read as Marc
    ASSERT_NE(deck, nullptr) << mio_last_error();
    mio_mesh_free(deck);
    mio_mesh* m = build_tet_mesh();
    ASSERT_EQ(mio_write(dat.c_str(), m, nullptr), MIO_OK) << mio_last_error();
    mio_mesh* back = mio_read(dat.c_str(), "tecplot");
    ASSERT_NE(back, nullptr) << mio_last_error();
    mio_mesh_free(back);
    mio_mesh_free(m);
    std::remove(dat.c_str());
}

// --- analysis and editing (v16.23.0) ----------------------------------------

TEST(CApi, FeatureEdgesOfACube) {
    mio_mesh* m = capi_cube_surface();
    mio_feature_edges_opts opts;
    mio_feature_edges_opts_init(&opts);
    EXPECT_DOUBLE_EQ(opts.feature_angle, 30.0);
    EXPECT_EQ(opts.boundary, 1);  // ON by default -- an all-zero struct is NOT the default
    mio_feature_edges_report report;
    mio_mesh* edges = mio_feature_edges(m, &opts, &report);
    ASSERT_NE(edges, nullptr) << mio_last_error();
    // The fan diagonals of each face are flat; the 12 cube edges are sharp.
    EXPECT_EQ(report.num_feature, 12);
    EXPECT_EQ(report.num_boundary, 0);
    ASSERT_EQ(mio_mesh_num_cell_blocks(edges), 1);
    EXPECT_EQ(block_type(edges, 0), "line");
    mio_mesh_free(edges);
    opts.feature_angle = 200.0;
    EXPECT_EQ(mio_feature_edges(m, &opts, nullptr), nullptr);
    EXPECT_EQ(mio_feature_edges(nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, TessellateATetra10) {
    // One straight-sided tetra10 carrying T = x: the default two levels give
    // the 35-point, 64-tetra lattice, and T interpolates exactly.
    // clang-format off
    const std::vector<double> pts = {0,   0,   0,   1,   0,   0,   0, 1,   0,   0, 0,   1,
                                     0.5, 0,   0,   0.5, 0.5, 0,   0, 0.5, 0,   0, 0,   0.5,
                                     0.5, 0,   0.5, 0,   0.5, 0.5};
    // clang-format on
    const std::vector<int64_t> conn = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 10, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra10", 1, 10, MIO_INT64, conn.data()), MIO_OK);
    std::vector<double> temperature(10);
    for (std::size_t i = 0; i < 10; ++i)
        temperature[i] = pts[3 * i];
    const int64_t tshape[1] = {10};
    ASSERT_EQ(mio_mesh_add_point_data(m, "T", MIO_FLOAT64, 1, tshape, temperature.data()), MIO_OK);

    mio_tessellate_opts opts;
    mio_tessellate_opts_init(&opts);
    EXPECT_EQ(opts.levels, 2);  // an all-zero struct is NOT the default
    EXPECT_EQ(opts.curved, 1);
    EXPECT_EQ(opts.fields, 1);
    EXPECT_EQ(opts.record_stencil, 0);
    mio_tessellate_report report;
    mio_mesh* out = mio_tessellate(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_source_points, 10);
    EXPECT_EQ(report.num_source_cells, 1);
    EXPECT_EQ(report.num_curved_source_cells, 1);
    EXPECT_EQ(report.num_pass_through_source_cells, 0);
    EXPECT_EQ(report.num_points, 35);
    EXPECT_EQ(report.num_cells, 64);
    EXPECT_EQ(report.num_boundary_facets, 64);
    EXPECT_EQ(report.num_interior_facets, 96);
    EXPECT_EQ(report.num_facets_with_bad_count, 0);
    EXPECT_EQ(mio_mesh_num_points(out), 35);

    const void* data = nullptr;
    mio_dtype dtype = MIO_FLOAT64;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mesh_get_point_data(out, "tessellate:source_point", &data, &dtype, &ndim, shape),
              MIO_OK)
        << mio_last_error();
    EXPECT_EQ(dtype, MIO_INT64);
    ASSERT_EQ(ndim, 1);
    ASSERT_EQ(shape[0], 35);
    const auto* source = static_cast<const int64_t*>(data);
    for (int64_t i = 0; i < 4; ++i)
        EXPECT_EQ(source[i], i);  // the kept corners first, ascending
    EXPECT_EQ(source[4], -1);
    EXPECT_EQ(mio_mesh_cell_data_num_blocks(out, "tessellate:source_cell"), 1);
    EXPECT_EQ(mio_mesh_cell_data_num_blocks(out, "tessellate:sub_index"), 1);
    // T = x interpolates exactly onto the lattice.
    ASSERT_EQ(mio_mesh_get_point_data(out, "T", &data, &dtype, &ndim, shape), MIO_OK);
    ASSERT_EQ(dtype, MIO_FLOAT64);
    const void* points = nullptr;
    mio_dtype pdtype = MIO_FLOAT64;
    const int64_t np = mio_mesh_num_points(out);
    const int64_t dim = mio_mesh_point_dim(out);
    ASSERT_EQ(mio_mesh_get_points(out, &points, &pdtype), MIO_OK);
    ASSERT_EQ(pdtype, MIO_FLOAT64);
    for (int64_t i = 0; i < np; ++i)
        EXPECT_NEAR(static_cast<const double*>(data)[i],
                    static_cast<const double*>(points)[i * dim], 1e-12);
    mio_mesh_free(out);

    // NULL opts means the defaults; the report is optional.
    out = mio_tessellate(m, nullptr, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(mio_mesh_num_points(out), 35);
    EXPECT_EQ(
        mio_mesh_get_point_data(out, "tessellate:stencil", nullptr, nullptr, nullptr, nullptr),
        MIO_ERR_NOT_FOUND);
    mio_mesh_free(out);

    // One level with the stencil recorded: (P, 27) stencil and weights.
    opts.levels = 1;
    opts.record_stencil = 1;
    out = mio_tessellate(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_cells, 8);
    ASSERT_EQ(mio_mesh_get_point_data(out, "tessellate:weights", &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(dtype, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 10);
    EXPECT_EQ(shape[1], 27);
    ASSERT_EQ(mio_mesh_get_point_data(out, "tessellate:stencil", &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    mio_mesh_free(out);

    // curved = 0 passes the cell through.
    opts.curved = 0;
    out = mio_tessellate(m, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_pass_through_source_cells, 1);
    EXPECT_EQ(report.num_cells, 1);
    mio_mesh_free(out);

    // Negative levels and a NULL mesh fail with an error.
    opts.curved = 1;
    opts.levels = -1;
    EXPECT_EQ(mio_tessellate(m, &opts, &report), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("levels"), std::string::npos) << mio_last_error();
    EXPECT_EQ(mio_tessellate(nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, SegmentAndFitACube) {
    mio_mesh* m = capi_cube_surface();
    mio_segment_opts sopts;
    mio_segment_opts_init(&sopts);
    EXPECT_DOUBLE_EQ(sopts.feature_angle, 30.0);
    EXPECT_EQ(sopts.min_patch_cells, 1);  // an all-zero struct is NOT the default
    EXPECT_EQ(sopts.method, MIO_SEGMENT_ANGLE);
    EXPECT_DOUBLE_EQ(sopts.tolerance, 1e-3);
    {
        // The primitive method (v17.9.0) finds the cube's six planes too.
        mio_segment_opts popts = sopts;
        popts.method = MIO_SEGMENT_PRIMITIVE;
        mio_segment_report preport;
        mio_mesh* pseg = mio_segment(m, &popts, &preport);
        ASSERT_NE(pseg, nullptr) << mio_last_error();
        EXPECT_EQ(preport.num_patches, 6);
        mio_mesh_free(pseg);
        popts.method = 4;
        EXPECT_EQ(mio_segment(m, &popts, nullptr), nullptr);
        // A struct zero-filled before the fields existed is the angle method.
        popts.method = MIO_SEGMENT_ANGLE;
        popts.tolerance = 0.0;
        pseg = mio_segment(m, &popts, nullptr);
        ASSERT_NE(pseg, nullptr) << mio_last_error();
        mio_mesh_free(pseg);
    }
    mio_segment_report sreport;
    mio_mesh* seg = mio_segment(m, &sopts, &sreport);
    ASSERT_NE(seg, nullptr) << mio_last_error();
    EXPECT_EQ(sreport.num_patches, 6);  // the fan diagonals are flat
    EXPECT_EQ(sreport.num_feature_edges, 12);
    mio_regions* regions = mio_regions_create(seg);
    EXPECT_EQ(mio_regions_count(regions), 6);
    mio_regions_free(regions);

    mio_fit_primitives_opts fopts;
    mio_fit_primitives_opts_init(&fopts);
    EXPECT_EQ(fopts.types, MIO_FIT_ALL_TYPES);
    mio_fit_primitives_report freport;
    mio_mesh* fit = mio_fit_primitives(seg, &fopts, &freport);
    ASSERT_NE(fit, nullptr) << mio_last_error();
    EXPECT_EQ(freport.num_patches, 6);
    EXPECT_EQ(freport.num_plane, 6);
    mio_mesh_free(fit);

    // No label array on the raw cube; a bad tolerance; NULL meshes.
    EXPECT_EQ(mio_fit_primitives(m, &fopts, nullptr), nullptr);
    fopts.patch_data = "";
    fopts.tolerance = 0.0;
    EXPECT_EQ(mio_fit_primitives(m, &fopts, nullptr), nullptr);
    sopts.feature_angle = 200.0;
    EXPECT_EQ(mio_segment(m, &sopts, nullptr), nullptr);
    EXPECT_EQ(mio_segment(nullptr, nullptr, nullptr), nullptr);
    EXPECT_EQ(mio_fit_primitives(nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(seg);
    mio_mesh_free(m);
}

TEST(CApi, FitCadBuildsOrNamesTheBackend) {
    mio_mesh* m = capi_cube_surface();
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_init(&opts);
    EXPECT_EQ(opts.max_demotion_rounds, 8);
    EXPECT_EQ(opts.exact_edges, 1);
    EXPECT_EQ(opts.types, MIO_FIT_ALL_TYPES);
    EXPECT_EQ(opts.freeform, 1);
    EXPECT_EQ(opts.freeform_split, 1);
    EXPECT_EQ(opts.fillets, 1);
    EXPECT_EQ(opts.threads, 1);
    EXPECT_EQ(opts.merge_faces, 1);
    EXPECT_EQ(opts.backend, MIO_FIT_CAD_AUTO);
    EXPECT_EQ(opts.reserved32, 0);
    mio_fit_cad_report report;
    mio_mesh* solid = mio_fit_cad(m, nullptr, &opts, &report);
    {
        // The native backend runs in every build (v17.9.0).
        mio_fit_cad_opts native = opts;
        native.backend = MIO_FIT_CAD_NATIVE;
        mio_fit_cad_report nr;
        mio_fit_cad_report_ex nex;
        mio_mesh* skin = mio_fit_cad_ex(m, nullptr, &native, &nr, &nex);
        ASSERT_NE(skin, nullptr) << mio_last_error();
        EXPECT_EQ(nex.backend, MIO_FIT_CAD_NATIVE);
        EXPECT_EQ(nr.valid, 1);
        EXPECT_EQ(nr.num_plane, 6);
        EXPECT_EQ(nr.num_faceted_faces, 0);
        mio_mesh_free(skin);
        native.backend = 9;
        EXPECT_EQ(mio_fit_cad(m, nullptr, &native, nullptr), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("backend"), std::string::npos);
    }
    if (mio_has_occt()) {
        ASSERT_NE(solid, nullptr) << mio_last_error();
        EXPECT_EQ(report.valid, 1);
        EXPECT_EQ(report.num_plane, 6);
        EXPECT_EQ(report.num_faces, 6);
        EXPECT_EQ(report.num_exact_edges, 12);  // each a single mesh edge
        EXPECT_NEAR(report.volume_ratio, 1.0, 1e-9);
        EXPECT_EQ(report.num_bspline, 0);
        EXPECT_EQ(report.num_split_patches, 0);
        EXPECT_EQ(report.num_fillet, 0);
        EXPECT_EQ(report.num_helical, 0);
        EXPECT_EQ(report.num_variable_fillet, 0);
        EXPECT_EQ(report.num_blend, 0);
        EXPECT_EQ(report.num_rounded_thread, 0);
        EXPECT_EQ(mio_mesh_num_cell_blocks(solid), 1);
        mio_mesh_free(solid);
        // The extended report (v17.2.0): a box has none of its features.
        mio_fit_cad_report_ex ex;
        std::memset(&ex, 0xff, sizeof ex);
        solid = mio_fit_cad_ex(m, nullptr, &opts, nullptr, &ex);
        ASSERT_NE(solid, nullptr) << mio_last_error();
        EXPECT_EQ(ex.num_vertex_blend, 0);
        EXPECT_EQ(ex.num_chamfer, 0);
        EXPECT_EQ(ex.num_conic_blend, 0);
        EXPECT_EQ(ex.num_tapered_thread, 0);
        EXPECT_EQ(ex.num_variable_pitch_thread, 0);
        EXPECT_EQ(ex.backend, MIO_FIT_CAD_OCCT);
        for (int64_t r : ex.reserved)
            EXPECT_EQ(r, 0);
        mio_mesh_free(solid);
        // Fillets, threads and splitting off still build the six planes.
        opts.freeform_split = 0;
        opts.fillets = 0;
        opts.threads = 0;
        solid = mio_fit_cad(m, nullptr, &opts, &report);
        ASSERT_NE(solid, nullptr) << mio_last_error();
        EXPECT_EQ(report.num_plane, 6);
        mio_mesh_free(solid);
        opts.max_demotion_rounds = -1;
        EXPECT_EQ(mio_fit_cad(m, nullptr, &opts, nullptr), nullptr);
    } else {
        // Without OpenCASCADE the automatic backend is the native one, and
        // asking for OpenCASCADE names the flag.
        ASSERT_NE(solid, nullptr) << mio_last_error();
        EXPECT_EQ(report.num_plane, 6);
        mio_mesh_free(solid);
        opts.backend = MIO_FIT_CAD_OCCT;
        EXPECT_EQ(mio_fit_cad(m, nullptr, &opts, nullptr), nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_OCCT"),
                  std::string::npos);
    }
    EXPECT_EQ(mio_fit_cad(nullptr, nullptr, nullptr, nullptr), nullptr);
    EXPECT_EQ(mio_fit_cad_ex(nullptr, nullptr, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
}

namespace {

/// The unit cube's skin as n x n quads per side, outward-wound (the core
/// tests' analytic::box_skin, through the C surface).
mio_mesh* capi_box_skin(int n) {
    const int m = n + 1;
    std::vector<std::int64_t> ids(static_cast<std::size_t>(m) * m * m, -1);
    std::vector<double> pts;
    const auto pid = [&](const int x[3]) {
        std::int64_t& id = ids[(static_cast<std::size_t>(x[0]) * m + x[1]) * m + x[2]];
        if (id < 0) {
            id = static_cast<std::int64_t>(pts.size() / 3);
            for (int d = 0; d < 3; ++d)
                pts.push_back(double(x[d]) / n);
        }
        return id;
    };
    std::vector<std::int64_t> quads;
    for (int axis = 0; axis < 3; ++axis)
        for (int side = 0; side < 2; ++side) {
            const int u = (axis + 1) % 3, v = (axis + 2) % 3;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) {
                    const int corners[4][2] = {{a, b}, {a + 1, b}, {a + 1, b + 1}, {a, b + 1}};
                    std::int64_t ring[4];
                    for (int c = 0; c < 4; ++c) {
                        int x[3];
                        x[axis] = side * n;
                        x[u] = corners[c][0];
                        x[v] = corners[c][1];
                        ring[c] = pid(x);
                    }
                    // (u, v, axis) is right-handed, so this ring faces +axis.
                    if (side == 0)
                        std::swap(ring[1], ring[3]);
                    quads.insert(quads.end(), ring, ring + 4);
                }
        }
    mio_mesh* mesh = mio_mesh_create();
    EXPECT_EQ(mio_mesh_set_points(mesh, MIO_FLOAT64, static_cast<std::int64_t>(pts.size() / 3), 3,
                                  pts.data()),
              MIO_OK);
    EXPECT_EQ(mio_mesh_add_cell_block(mesh, "quad", static_cast<std::int64_t>(quads.size() / 4), 4,
                                      MIO_INT64, quads.data()),
              MIO_OK);
    return mesh;
}

}  // namespace

TEST(CApi, RefinePatches) {
    mio_refine_patches_opts opts;
    std::memset(&opts, 0xff, sizeof opts);
    mio_refine_patches_opts_init(&opts);
    EXPECT_DOUBLE_EQ(opts.tolerance, 1e-3);
    EXPECT_DOUBLE_EQ(opts.normal_weight, 1.0);
    EXPECT_DOUBLE_EQ(opts.epsilon, 1e-3);
    EXPECT_DOUBLE_EQ(opts.feature_angle, 30.0);
    EXPECT_EQ(opts.max_iterations, 10);
    EXPECT_EQ(opts.min_patch_cells, 3);
    EXPECT_EQ(opts.ransac_seed, MIO_REFINE_DEFAULT_SEED);
    EXPECT_EQ(opts.ransac_seed, 0x6d657368696f2b2bull);
    EXPECT_EQ(opts.types, MIO_FIT_ALL_TYPES);
    EXPECT_EQ(opts.relative, 1);
    EXPECT_EQ(opts.allow_split, 1);
    EXPECT_EQ(opts.allow_merge, 1);
    EXPECT_EQ(opts.prefix, nullptr);
    for (int64_t r : opts.reserved)
        EXPECT_EQ(r, 0);

    // With a feature angle above the box's 90 degrees the whole box is one
    // patch; refinement splits the six planes out of it.
    mio_mesh* box = capi_box_skin(3);
    mio_segment_opts sopts;
    mio_segment_opts_init(&sopts);
    sopts.feature_angle = 100.0;
    mio_segment_report sreport;
    mio_mesh* seg = mio_segment(box, &sopts, &sreport);
    ASSERT_NE(seg, nullptr) << mio_last_error();
    EXPECT_EQ(sreport.num_patches, 1);

    opts.feature_angle = 100.0;
    opts.max_iterations = 20;
    mio_refine_patches_report report;
    std::memset(&report, 0xff, sizeof report);
    mio_mesh* out = mio_refine_patches(seg, &opts, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_patches, 6);
    EXPECT_EQ(report.num_primitive, 6);
    EXPECT_GT(report.num_splits, 0);
    EXPECT_NEAR(report.explained_area_fraction, 1.0, 1e-12);
    EXPECT_EQ(report.initial_explained_area_fraction, 0.0);
    EXPECT_GT(report.tolerance, 0.0);
    EXPECT_LE(report.hausdorff_max, report.tolerance * (1.0 + 1e-12));
    EXPECT_GE(report.stop, MIO_REFINE_UNCHANGED);
    EXPECT_LE(report.stop, MIO_REFINE_ENERGY_ROSE);
    EXPECT_EQ(report.reserved32, 0);
    for (int64_t r : report.reserved)
        EXPECT_EQ(r, 0);
    EXPECT_EQ(mio_mesh_get_cell_data(out, "segment:patch", 0, nullptr, nullptr, nullptr, nullptr),
              MIO_OK);
    mio_regions* regions = mio_regions_create(out);
    EXPECT_EQ(mio_regions_count(regions), 6);
    mio_regions_free(regions);
    mio_mesh_free(out);

    // A partition segmented at the creases is already six planes; NULL
    // options are the defaults.
    mio_mesh_free(seg);
    mio_segment_opts_init(&sopts);
    seg = mio_segment(box, &sopts, nullptr);
    ASSERT_NE(seg, nullptr) << mio_last_error();
    out = mio_refine_patches(seg, nullptr, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.stop, MIO_REFINE_UNCHANGED);
    EXPECT_EQ(report.num_patches, 6);
    EXPECT_EQ(report.num_splits + report.num_merges + report.num_absorbed, 0);
    mio_mesh_free(out);

    // A NULL mesh, a mesh without "segment:patch" and an option out of range
    // fail with an error.
    EXPECT_EQ(mio_refine_patches(nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()), "");
    EXPECT_EQ(mio_refine_patches(box, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("segment"), std::string::npos) << mio_last_error();
    mio_refine_patches_opts_init(&opts);
    opts.tolerance = 0.0;
    EXPECT_EQ(mio_refine_patches(seg, &opts, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("tolerance"), std::string::npos)
        << mio_last_error();
    mio_mesh_free(seg);
    mio_mesh_free(box);
}

TEST(CApi, FitCad2Refines) {
    mio_fit_cad_opts_ex ex;
    std::memset(&ex, 0xff, sizeof ex);
    mio_fit_cad_opts_ex_init(&ex);
    EXPECT_EQ(ex.refine_segmentation, 0);
    EXPECT_EQ(ex.freeform_general, 0);
    EXPECT_EQ(ex.freeform_tolerance, 0.0);
    EXPECT_DOUBLE_EQ(ex.freeform_smoothing, 1e-12);
    EXPECT_EQ(ex.max_freeform_poles, 40);
    EXPECT_EQ(ex.max_freeform_iterations, 30);
    EXPECT_EQ(ex.max_freeform_splits, 2);
    EXPECT_EQ(ex.sew_tolerance, 0.0);
    EXPECT_EQ(ex.extend_surfaces, 0);
    EXPECT_EQ(ex.snap_junctions, 0);
    EXPECT_EQ(ex.face_deviation, 0);
    EXPECT_EQ(ex.reserved32b, 0);
    EXPECT_DOUBLE_EQ(ex.refine.tolerance, 1e-3);
    EXPECT_EQ(ex.refine.max_iterations, 10);
    EXPECT_EQ(ex.refine.ransac_seed, MIO_REFINE_DEFAULT_SEED);
    EXPECT_EQ(ex.refine.prefix, nullptr);
    for (int64_t r : ex.reserved)
        EXPECT_EQ(r, 0);

    // The native backend runs in every build.
    mio_mesh* m = capi_cube_surface();
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_init(&opts);
    opts.backend = MIO_FIT_CAD_NATIVE;
    ex.refine_segmentation = 1;
    mio_fit_cad_report report;
    mio_fit_cad_report_ex report_ex;
    mio_fit_cad_report_ex2 report_ex2;
    std::memset(&report_ex2, 0xff, sizeof report_ex2);
    mio_mesh* skin = mio_fit_cad2(m, nullptr, &opts, &ex, &report, &report_ex, &report_ex2);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report_ex2.refined, 1);
    EXPECT_TRUE(report_ex2.refinement_kept == 0 || report_ex2.refinement_kept == 1);
    EXPECT_GE(report_ex2.refine_iterations, 0);
    EXPECT_NEAR(report_ex2.explained_area_fraction, 1.0, 1e-12);
    EXPECT_EQ(report_ex.backend, MIO_FIT_CAD_NATIVE);
    EXPECT_EQ(report.num_plane, 6);
    for (int64_t r : report_ex2.reserved)
        EXPECT_EQ(r, 0);
    mio_mesh_free(skin);

    // NULL opts_ex leaves refinement off, and every report is optional.
    std::memset(&report_ex2, 0xff, sizeof report_ex2);
    skin = mio_fit_cad2(m, nullptr, &opts, nullptr, nullptr, nullptr, &report_ex2);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report_ex2.refined, 0);
    EXPECT_EQ(report_ex2.refine_iterations, 0);
    EXPECT_EQ(report_ex2.num_refine_splits, 0);
    EXPECT_EQ(report_ex2.num_refine_merges, 0);
    EXPECT_EQ(report_ex2.explained_area_fraction, 0.0);
    // The Phase 3 fields are filled on every run; nothing measured is NaN.
    EXPECT_EQ(report_ex2.closed, 1);
    EXPECT_EQ(report_ex2.refinement_kept, 0);
    EXPECT_EQ(report_ex2.num_freeform_general, 0);
    EXPECT_EQ(report_ex2.num_snapped_junctions, 0);
    EXPECT_TRUE(std::isnan(report_ex2.hausdorff_max));
    EXPECT_TRUE(std::isnan(report_ex2.hausdorff_rms));
    for (int64_t r : report_ex2.reserved)
        EXPECT_EQ(r, 0);
    mio_mesh_free(skin);

    // mio_fit_cad_ex forwards to mio_fit_cad2 unchanged.
    skin = mio_fit_cad_ex(m, nullptr, &opts, &report, &report_ex);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_plane, 6);
    EXPECT_EQ(report_ex.backend, MIO_FIT_CAD_NATIVE);
    mio_mesh_free(skin);

    // A refinement option out of range fails only when refinement is on.
    ex.refine.min_patch_cells = 0;
    EXPECT_EQ(mio_fit_cad2(m, nullptr, &opts, &ex, nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("min_patch_cells"), std::string::npos)
        << mio_last_error();
    ex.refine_segmentation = 0;
    skin = mio_fit_cad2(m, nullptr, &opts, &ex, nullptr, nullptr, nullptr);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    mio_mesh_free(skin);
    EXPECT_EQ(mio_fit_cad2(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, FitCad2PhaseThreeOptions) {
    mio_mesh* m = capi_cube_surface();
    mio_fit_cad_opts_ex ex;
    mio_fit_cad_opts_ex_init(&ex);
    // A struct zero-filled before the fields existed: zero poles and
    // iterations read as their defaults, so the call runs.
    ex.max_freeform_poles = 0;
    ex.max_freeform_iterations = 0;
    ex.freeform_smoothing = 0.0;
    ex.max_freeform_splits = 0;
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_init(&opts);
    opts.backend = MIO_FIT_CAD_NATIVE;
    mio_mesh* skin = mio_fit_cad2(m, nullptr, &opts, &ex, nullptr, nullptr, nullptr);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    mio_mesh_free(skin);
    // Any other value out of range is refused, as given.
    ex.max_freeform_poles = 3;
    EXPECT_EQ(mio_fit_cad2(m, nullptr, &opts, &ex, nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("max_freeform_poles"), std::string::npos)
        << mio_last_error();
    ex.max_freeform_poles = 40;
    ex.sew_tolerance = -1.0;
    EXPECT_EQ(mio_fit_cad2(m, nullptr, &opts, &ex, nullptr, nullptr, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("sew_tolerance"), std::string::npos)
        << mio_last_error();
    ex.sew_tolerance = 0.0;

    if (!mio_has_occt()) {
        mio_mesh_free(m);
        GTEST_SKIP() << "built without OpenCASCADE";
    }
    // Snapped corners and the measured deviation, on the box's six planes.
    opts.backend = MIO_FIT_CAD_OCCT;
    ex.snap_junctions = 1;
    ex.face_deviation = 1;
    mio_fit_cad_report report;
    mio_fit_cad_report_ex2 report_ex2;
    std::memset(&report_ex2, 0xff, sizeof report_ex2);
    skin = mio_fit_cad2(m, nullptr, &opts, &ex, &report, nullptr, &report_ex2);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report.valid, 1);
    EXPECT_EQ(report_ex2.closed, 1);
    EXPECT_GT(report_ex2.num_snapped_junctions, 0);
    EXPECT_TRUE(std::isfinite(report_ex2.hausdorff_max));
    EXPECT_TRUE(std::isfinite(report_ex2.hausdorff_rms));
    EXPECT_LE(report_ex2.hausdorff_rms, report_ex2.hausdorff_max);
    EXPECT_EQ(report_ex2.num_extended_edges, 0);
    EXPECT_EQ(report_ex2.num_freeform_general, 0);
    for (int64_t r : report_ex2.reserved)
        EXPECT_EQ(r, 0);
    // The per-face deviation travels as field data, (faces, 2) and (faces).
    const void* data = nullptr;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[8] = {};
    ASSERT_EQ(mio_mesh_get_field_data(skin, "fit_cad:face_hausdorff", &data, &dt, &ndim, shape),
              MIO_OK)
        << mio_last_error();
    EXPECT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], report.num_faces);
    EXPECT_EQ(shape[1], 2);
    ASSERT_EQ(mio_mesh_get_field_data(skin, "fit_cad:face_rms", &data, &dt, &ndim, shape), MIO_OK)
        << mio_last_error();
    EXPECT_EQ(shape[0], report.num_faces);
    mio_mesh_free(skin);
    mio_mesh_free(m);
}

TEST(CApi, BrepSideChannel) {
    // A degree-(1, 1) face with one outer loop of one trim, and an edge on it.
    const auto path = mt::temp_path(".cad.json");
    std::ofstream(path, std::ios::binary)
        << R"({"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "surface": {
           "is_trimmed": true, "degrees": [1, 1], "knot_vectors": [[0, 1], [0, 1]],
           "control_points": [[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]},
           "boundary_loops": [{"loop_type": "outer", "trimming_curves": [
             {"trim_index": 0, "curve_direction": false, "parameter_curve": {"degree": 1,
              "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
              "control_points": [[[0, 0, 0, 1]], [[1, 0, 0, 1]]]}}]}]}],
           "edges": [{"brep_id": 3, "topology": [{"brep_id": 2, "trim_index": 0}]}]}]})";
    mio_format_info* info = nullptr;
    mio_mesh* mesh = mio_read_with_info(path.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    ASSERT_NE(info, nullptr);
    char fmt[16];
    ASSERT_EQ(mio_format_info_format(info, fmt, sizeof(fmt)), 8);
    EXPECT_STREQ(fmt, "cad_json");
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_FILE), 1);
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_FACES), 1);
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_LOOPS), 1);
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_TRIMS), 1);
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_EDGES), 1);
    EXPECT_EQ(mio_brep_info_count(info, MIO_BREP_REFS), 1);
    EXPECT_EQ(mio_brep_info_count(info, 99), -1);
    EXPECT_EQ(mio_nurbs_info_count(info, MIO_NURBS_PATCHES), -1);  // not an mfem handle
    int64_t v = -7;
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_FACES, 0, 0, &v), MIO_OK);
    EXPECT_EQ(v, 2);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_FACES, 0, 5, &v), MIO_OK);  // reduced knots
    EXPECT_EQ(v, 1);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_FACES, 0, 9, &v), MIO_OK);  // trims
    EXPECT_EQ(v, 1);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_TRIMS, 0, 3, &v), MIO_OK);  // direction
    EXPECT_EQ(v, 0);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_TRIMS, 0, 10, &v), MIO_OK);  // curve degree
    EXPECT_EQ(v, 1);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_EDGES, 0, 1, &v), MIO_OK);  // no 3-D curve
    EXPECT_EQ(v, 0);
    EXPECT_EQ(mio_brep_info_int(info, MIO_BREP_EDGES, 0, 10, &v), MIO_ERR_INVALID_ARG);
    ASSERT_EQ(mio_brep_info_int(info, MIO_BREP_REFS, 0, 2, &v), MIO_OK);  // direction absent
    EXPECT_EQ(v, -1);
    EXPECT_EQ(mio_brep_info_int(info, MIO_BREP_FACES, 1, 0, &v), MIO_ERR_INVALID_ARG);
    const void* data = nullptr;
    mio_dtype dtype;
    int32_t ndim;
    int64_t shape[MIO_MAX_NDIM]{};
    ASSERT_EQ(mio_brep_info_array(info, MIO_BREP_PATCHES, 0, 0, &data, &dtype, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(dtype, MIO_FLOAT64);
    EXPECT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 4);
    EXPECT_EQ(shape[1], 3);
    EXPECT_EQ(static_cast<const double*>(data)[3], 2.0);
    ASSERT_EQ(mio_brep_info_array(info, MIO_BREP_PATCHES, 0, 2, &data, &dtype, &ndim, shape),
              MIO_OK);  // the full knot vector in u
    EXPECT_EQ(shape[0], 4);
    ASSERT_EQ(mio_brep_info_array(info, MIO_BREP_TRIMS, 0, 14, &data, &dtype, &ndim, shape),
              MIO_OK);  // the active range
    EXPECT_EQ(shape[0], 2);
    ASSERT_EQ(mio_brep_info_array(info, MIO_BREP_BREPS, 0, 1, &data, &dtype, &ndim, shape), MIO_OK);
    EXPECT_EQ(dtype, MIO_INT64);
    EXPECT_EQ(static_cast<const int64_t*>(data)[0], 0);
    char buf[16];
    EXPECT_EQ(mio_brep_info_string(info, MIO_BREP_LOOPS, 0, 0, buf, sizeof(buf)), 5);
    EXPECT_STREQ(buf, "outer");
    EXPECT_EQ(mio_brep_info_string(info, MIO_BREP_FACES, 0, 7, buf, sizeof(buf)), -1);

    const auto output = mt::temp_path(".cad.json");
    ASSERT_EQ(mio_write_with_info(output.c_str(), mesh, nullptr, info), MIO_OK) << mio_last_error();
    mio_format_info* back = nullptr;
    mio_mesh* again = mio_read_with_info(output.c_str(), nullptr, nullptr, &back);
    ASSERT_NE(again, nullptr) << mio_last_error();
    EXPECT_EQ(mio_brep_info_count(back, MIO_BREP_TRIMS), 1);
    EXPECT_EQ(mio_write_with_info(mt::temp_path(".mesh").c_str(), mesh, "mfem", info),
              MIO_ERR_INVALID_ARG);  // a cad_json side channel is not an mfem one
    mio_format_info_free(back);
    mio_mesh_free(again);
    mio_format_info_free(info);
    mio_mesh_free(mesh);
}

TEST(CApi, StepAndIgesReadOrNameTheBackend) {
    // A STEP written by the faceted writer reads back through OpenCASCADE.
    mio_mesh* m = capi_cube_surface();
    const auto step = mt::temp_path(".step");
    ASSERT_EQ(mio_write(step.c_str(), m, nullptr), MIO_OK) << mio_last_error();
    mio_format_info* info = nullptr;
    mio_mesh* back = mio_read_with_info(step.c_str(), nullptr, nullptr, &info);
    if (mio_has_occt()) {
        ASSERT_NE(back, nullptr) << mio_last_error();
        EXPECT_GE(mio_brep_info_count(info, MIO_BREP_FACES), 6);
        EXPECT_GT(mio_brep_info_count(info, MIO_BREP_EDGES), 0);
        mio_format_info_free(info);
        mio_mesh_free(back);
    } else {
        EXPECT_EQ(back, nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_OCCT"),
                  std::string::npos);
    }
    mio_mesh_free(m);
}

TEST(CApi, WritesAStepSolid) {
    mio_mesh* m = capi_cube_surface();
    const std::string path =
        (std::filesystem::temp_directory_path() / "meshioplusplus_capi_cube.step").string();
    ASSERT_EQ(mio_write(path.c_str(), m, nullptr), MIO_OK) << mio_last_error();
    std::ifstream in(path);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_EQ(text.rfind("ISO-10303-21;", 0), 0u);
    EXPECT_NE(text.find("MANIFOLD_SOLID_BREP"), std::string::npos);
    std::filesystem::remove(path);
    mio_mesh_free(m);
}

TEST(CApi, RegionAdjacencyBetweenCellBlocks) {
    mio_mesh* m = mio_mesh_create();
    const double points[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1};
    const std::int64_t upper[] = {0, 1, 2, 3};
    const std::int64_t lower[] = {0, 2, 1, 4};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 5, 3, points), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 1, 4, MIO_INT64, upper), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "tetra", 1, 4, MIO_INT64, lower), MIO_OK);
    mio_mesh* adjacency = mio_region_adjacency(m, nullptr, 0);
    ASSERT_NE(adjacency, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_num_cell_blocks(adjacency), 1);
    EXPECT_EQ(block_type(adjacency, 0), "triangle");
    int64_t num_cells = 0, nodes_per_cell = 0;
    int32_t is_ragged = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(adjacency, 0, &num_cells, &nodes_per_cell, &is_ragged),
              MIO_OK);
    EXPECT_EQ(num_cells, 1);
    mio_mesh_free(adjacency);
    EXPECT_EQ(mio_region_adjacency(m, nullptr, 1), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, FindContactAndSplitInterface) {
    mio_mesh* mesh = mio_mesh_create();
    const double points[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1};
    const std::int64_t cells[] = {0, 1, 2, 3, 0, 2, 1, 4};
    const std::int64_t cell_a[] = {0};
    const std::int64_t cell_b[] = {1};
    const std::int64_t slave_node[] = {0};
    ASSERT_EQ(mio_mesh_set_points(mesh, MIO_FLOAT64, 5, 3, points), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(mesh, "tetra", 2, 4, MIO_INT64, cells), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(mesh, "upper", MIO_REGION_CELL, 3, -1, cell_a, 1), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(mesh, "lower", MIO_REGION_CELL, 3, -1, cell_b, 1), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(mesh, "slave", MIO_REGION_POINT, 0, -1, slave_node, 1), MIO_OK);

    mio_region_selector upper, lower, slave;
    mio_region_selector_init(&upper, "upper");
    mio_region_selector_init(&lower, "lower");
    mio_region_selector_init(&slave, "slave");
    upper.kind = lower.kind = MIO_REGION_CELL;
    slave.kind = MIO_REGION_POINT;

    mio_find_interface_opts find_opts;
    mio_find_interface_opts_init(&find_opts);
    EXPECT_DOUBLE_EQ(find_opts.angle_tolerance, 30.0);
    mio_find_interface_result* found =
        mio_find_interface(mesh, &upper, nullptr, &lower, &find_opts);
    ASSERT_NE(found, nullptr) << mio_last_error();
    std::int64_t pairs = 0, side_count = 0;
    double area = 0.0;
    ASSERT_EQ(mio_find_interface_result_report(found, &pairs, &area, nullptr, nullptr, nullptr),
              MIO_OK);
    EXPECT_EQ(pairs, 1);
    EXPECT_GT(area, 0.0);
    const std::int64_t* side_a = mio_find_interface_result_side_a(found, &side_count);
    ASSERT_NE(side_a, nullptr);
    ASSERT_EQ(side_count, 1);
    const std::int64_t side_entry[] = {side_a[0], side_a[1]};
    ASSERT_EQ(mio_mesh_add_region(mesh, "cut", MIO_REGION_SIDE, 2, -1, side_entry, 2), MIO_OK);

    mio_contact_pairs_opts contact_opts;
    mio_contact_pairs_opts_init(&contact_opts);
    contact_opts.tolerance = 10.0;
    mio_contact_pairs_result* projected =
        mio_contact_pairs(mesh, &slave, nullptr, &lower, &contact_opts);
    ASSERT_NE(projected, nullptr) << mio_last_error();
    std::int64_t count = 0, unmatched = 0;
    ASSERT_EQ(mio_contact_pairs_result_info(projected, &count, &unmatched), MIO_OK);
    EXPECT_EQ(count, 1);
    EXPECT_EQ(unmatched, 0);
    ASSERT_NE(mio_contact_pairs_master_cell(projected), nullptr);
    EXPECT_EQ(mio_contact_pairs_slave_point(projected)[0], 0);
    EXPECT_NE(mio_contact_pairs_closest_point(projected), nullptr);
    mio_contact_pairs_result_free(projected);

    mio_split_interface_opts split_opts;
    mio_split_interface_opts_init(&split_opts);
    split_opts.add_cohesive = 1;
    std::int64_t duplicated = 0, cohesive = 0;
    mio_mesh* split = mio_split_interface(mesh, &upper, &split_opts, &duplicated, &cohesive);
    EXPECT_EQ(split, nullptr);  // a Cell selector is not a Side selector
    mio_region_selector cut;
    mio_region_selector_init(&cut, "cut");
    cut.kind = MIO_REGION_SIDE;
    split = mio_split_interface(mesh, &cut, &split_opts, &duplicated, &cohesive);
    ASSERT_NE(split, nullptr) << mio_last_error();
    EXPECT_EQ(duplicated, 3);
    EXPECT_EQ(cohesive, 1);
    mio_mesh_free(split);
    mio_find_interface_result_free(found);
    mio_mesh_free(mesh);
}

TEST(CApi, HausdorffOfAMeshWithItselfIsZero) {
    mio_mesh* a = capi_cube_surface();
    mio_mesh* b = capi_cube_surface();
    mio_hausdorff_opts opts;
    mio_hausdorff_opts_init(&opts);
    opts.face_samples = 2;
    mio_hausdorff_report report;
    ASSERT_EQ(mio_hausdorff_distance(a, b, &opts, &report), MIO_OK) << mio_last_error();
    EXPECT_NEAR(report.distance, 0.0, 1e-12);  // face samples: rounding, not exactly 0
    EXPECT_EQ(report.num_samples_a, 8 + 12 * 4);
    EXPECT_NE(mio_hausdorff_distance(a, b, &opts, nullptr), MIO_OK);
    mio_mesh_free(a);
    mio_mesh_free(b);
}

TEST(CApi, EditRegionsAndPeriodicPairs) {
    mio_mesh* m = capi_cube_surface();
    const std::int64_t bottom[] = {0, 1, 2, 3};
    const std::int64_t top[] = {4, 5, 6, 7};
    ASSERT_EQ(mio_mesh_add_region(m, "bottom", MIO_REGION_POINT, 0, -1, bottom, 4), MIO_OK);
    ASSERT_EQ(mio_mesh_add_region(m, "top", MIO_REGION_POINT, 0, -1, top, 4), MIO_OK);

    mio_region_selector in[2];
    mio_region_selector_init(&in[0], "bottom");
    mio_region_selector_init(&in[1], "top");
    EXPECT_EQ(in[0].kind, -1);
    EXPECT_EQ(in[0].dim, MIO_REGION_ANY);
    mio_mesh* u =
        mio_edit_regions(m, MIO_REGION_UNION, in, 2, "all", MIO_REGION_ANY, MIO_REGION_ANY, 1);
    ASSERT_NE(u, nullptr) << mio_last_error();
    mio_regions* regs = mio_regions_create(u);
    EXPECT_EQ(mio_regions_count(regs), 3);
    mio_regions_free(regs);
    ASSERT_EQ(mio_mesh_remove_region(u, 0), MIO_OK);
    regs = mio_regions_create(u);
    EXPECT_EQ(mio_regions_count(regs), 2);
    mio_regions_free(regs);
    EXPECT_NE(mio_mesh_remove_region(u, 9), MIO_OK);
    mio_mesh_free(u);
    EXPECT_EQ(mio_edit_regions(m, 42, in, 2, "x", MIO_REGION_ANY, MIO_REGION_ANY, 1), nullptr);

    mio_periodic_opts opts;
    mio_periodic_opts_init(&opts);
    EXPECT_EQ(opts.matrix[0], 1.0);
    EXPECT_EQ(opts.matrix[15], 1.0);
    opts.matrix[11] = 1.0;  // translate +1 in z
    mio_periodic_pairs* pairs = mio_match_periodic_nodes(m, &in[0], &in[1], &opts);
    ASSERT_NE(pairs, nullptr) << mio_last_error();
    std::int64_t n = 0, unmatched = -1, fixed = -1;
    double residual = -1.0;
    ASSERT_EQ(mio_periodic_pairs_info(pairs, &n, &unmatched, &fixed, &residual), MIO_OK);
    EXPECT_EQ(n, 4);
    EXPECT_EQ(unmatched, 0);
    EXPECT_EQ(fixed, 0);
    EXPECT_EQ(residual, 0.0);
    std::int64_t count = 0;
    const std::int64_t* s = mio_periodic_pairs_slave(pairs, &count);
    const std::int64_t* t = mio_periodic_pairs_master(pairs, nullptr);
    ASSERT_EQ(count, 4);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(s[i], i);
        EXPECT_EQ(t[i], i + 4);
    }
    mio_periodic_pairs_free(pairs);
    opts.matrix[11] = 0.5;  // nothing lands on a master
    EXPECT_EQ(mio_match_periodic_nodes(m, &in[0], &in[1], &opts), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, CheckQualityGates) {
    mio_mesh* m = capi_cube_surface();
    mio_quality_gate_report report;
    char summary[512];
    ASSERT_EQ(mio_check_quality(m, "min_angle >= 30", 0, 0, &report, summary, sizeof(summary)),
              MIO_OK)
        << mio_last_error();
    EXPECT_EQ(report.passed, 1);
    EXPECT_EQ(report.num_checks, 3);
    EXPECT_EQ(report.num_failed, 0);
    EXPECT_EQ(report.num_cells, 12);
    EXPECT_GT(report.summary_length, 0);
    EXPECT_NE(std::string(summary).find("PASS"), std::string::npos);
    ASSERT_EQ(mio_check_quality(m, "min_angle >= 50", -1, -1, &report, nullptr, 0), MIO_OK);
    EXPECT_EQ(report.passed, 0);
    EXPECT_EQ(report.num_checks, 1);
    EXPECT_NE(mio_check_quality(m, "bogus >= 1", 0, 0, &report, nullptr, 0), MIO_OK);
    EXPECT_NE(mio_check_quality(m, nullptr, 0, 0, nullptr, nullptr, 0), MIO_OK);
    mio_mesh_free(m);
}

TEST(CApi, AgglomerateExAndBlendAndResample) {
    // Two unit hexes in a row: one group, and with coplanar merging the four
    // side faces pairs fuse (4 faces removed).
    const std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1,
                                     1, 1, 1, 0, 1, 1, 2, 0, 0, 2, 1, 0, 2, 0, 1, 2, 1, 1};
    const std::vector<std::int64_t> hexes = {0, 1, 2, 3, 4, 5, 6, 7, 1, 8, 9, 2, 5, 10, 11, 6};
    mio_mesh* m = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 12, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "hexahedron", 2, 8, MIO_INT64, hexes.data()), MIO_OK);
    mio_agglomerate_opts opts;
    mio_agglomerate_opts_init(&opts);
    EXPECT_EQ(opts.target_group_size, 8);
    EXPECT_EQ(opts.coplanar_angle, 1.0);
    opts.merge_coplanar_faces = 1;
    std::int64_t merged = -1, rejected = -1;
    mio_agglomerate_result* r = mio_agglomerate_ex(m, &opts, &merged, &rejected);
    ASSERT_NE(r, nullptr) << mio_last_error();
    EXPECT_EQ(merged, 4);
    EXPECT_EQ(rejected, 0);
    mio_agglomerate_result_free(r);
    opts.min_sphericity = 2.0;
    EXPECT_EQ(mio_agglomerate_ex(m, &opts, nullptr, nullptr), nullptr);

    // blend_steps: a point field halfway between two steps.
    const double u0[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const double u1[] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
    mio_mesh* b = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(b, MIO_FLOAT64, 12, 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(b, "hexahedron", 2, 8, MIO_INT64, hexes.data()), MIO_OK);
    const std::int64_t shape[] = {12};
    ASSERT_EQ(mio_mesh_add_point_data(m, "u", MIO_FLOAT64, 1, shape, u0), MIO_OK);
    ASSERT_EQ(mio_mesh_add_point_data(b, "u", MIO_FLOAT64, 1, shape, u1), MIO_OK);
    mio_mesh* half = mio_blend_steps(m, b, 0.25, 0);
    ASSERT_NE(half, nullptr) << mio_last_error();
    mio_mesh_free(half);
    mio_resample_opts ro;
    mio_resample_opts_init(&ro);
    EXPECT_NE(mio_sequence_resample(nullptr, "x_{index}.vtu", nullptr, &ro), MIO_OK);
    mio_mesh_free(b);
    mio_mesh_free(m);
}

// ---------------------------------------------------------------------------
// mio_format_info: the MDPA side channel on the flat ABI (v16.27.0)
// ---------------------------------------------------------------------------

namespace {

const char* const kCapiMdpaDeck = R"(Begin ModelPartData
    DOMAIN_SIZE 3
    SOLVER_TYPE static
End ModelPartData
Begin Table 1 TIME VALUE
    0.0 1.0
    1.0 2.0
End Table
Begin Properties 1
    DENSITY 7850
    CONSTITUTIVE_LAW LinearElastic3DLaw
End Properties
Begin Nodes
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
End Nodes
Begin Elements Element3D4N
1 1 1 2 3 4
End Elements
Begin Geometries Triangle3D3
7 1 2 3
End Geometries
Begin SubModelPart Inlet
    Begin SubModelPartData
        VELOCITY_X 2.5
    End SubModelPartData
    Begin SubModelPartNodes
        1
    End SubModelPartNodes
    Begin SubModelPartGeometries
        7
    End SubModelPartGeometries
    Begin SubModelPartConstraints
        1
    End SubModelPartConstraints
End SubModelPart
Begin Mesh 3
    Begin MeshNodes
        2
        4
    End MeshNodes
End Mesh
Begin Constraints LinearMasterSlaveConstraint
    1 1 DISPLACEMENT_X 2 DISPLACEMENT_X 1.0 0.0
End Constraints
)";

std::string capi_info_string(const mio_format_info* pInfo, int32_t Section, int64_t Index,
                             int32_t Field) {
    const int64_t n = mio_mdpa_info_string(pInfo, Section, Index, Field, nullptr, 0);
    if (n < 0)
        return "<error>";
    std::string out(static_cast<std::size_t>(n) + 1, '\0');
    mio_mdpa_info_string(pInfo, Section, Index, Field, out.data(), n + 1);
    out.resize(static_cast<std::size_t>(n));
    return out;
}

}  // namespace

TEST(CApi, ReadWithInfoKeepsTheMdpaSideChannel) {
    const std::string in = mt::temp_path("_capi_info.mdpa");
    {
        std::ofstream f(in);
        f << kCapiMdpaDeck;
    }
    // mio_read refuses what only the side channel can hold, by name.
    EXPECT_EQ(mio_read(in.c_str(), nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("mio_read_with_info"), std::string::npos)
        << mio_last_error();

    mio_format_info* info = nullptr;
    mio_mesh* mesh = mio_read_with_info(in.c_str(), nullptr, nullptr, &info);
    ASSERT_NE(mesh, nullptr) << mio_last_error();
    ASSERT_NE(info, nullptr);
    char fmt[16];
    EXPECT_EQ(mio_format_info_format(info, fmt, sizeof(fmt)), 4);
    EXPECT_STREQ(fmt, "mdpa");

    EXPECT_EQ(mio_mdpa_info_count(info, MIO_MDPA_MODEL_PART_DATA), 1);
    ASSERT_EQ(mio_mdpa_info_data_count(info, MIO_MDPA_MODEL_PART_DATA, 0), 1);
    EXPECT_EQ(mio_mdpa_info_data_kind(info, MIO_MDPA_MODEL_PART_DATA, 0, 0), MIO_MDPA_VALUE_TEXT);

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_TABLES), 1);
    EXPECT_EQ(capi_info_string(info, MIO_MDPA_TABLES, 0, 0), "1 TIME VALUE");
    const void* data = nullptr;
    mio_dtype dt = MIO_INT8;
    int32_t ndim = 0;
    int64_t shape[MIO_MAX_NDIM] = {};
    ASSERT_EQ(mio_mdpa_info_array(info, MIO_MDPA_TABLES, 0, 0, &data, &dt, &ndim, shape), MIO_OK);
    EXPECT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 2);
    EXPECT_DOUBLE_EQ(static_cast<const double*>(data)[3], 2.0);

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_GEOMETRIES), 1);
    EXPECT_EQ(capi_info_string(info, MIO_MDPA_GEOMETRIES, 0, 0), "Triangle3D3");
    EXPECT_EQ(capi_info_string(info, MIO_MDPA_GEOMETRIES, 0, 1), "triangle");
    ASSERT_EQ(mio_mdpa_info_array(info, MIO_MDPA_GEOMETRIES, 0, 1, &data, &dt, &ndim, shape),
              MIO_OK);
    EXPECT_EQ(dt, MIO_INT64);
    EXPECT_EQ(static_cast<const int64_t*>(data)[0], 7);

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_MESH_BLOCKS), 1);
    int64_t id = 0;
    ASSERT_EQ(mio_mdpa_info_int(info, MIO_MDPA_MESH_BLOCKS, 0, 0, &id), MIO_OK);
    EXPECT_EQ(id, 3);
    ASSERT_EQ(mio_mdpa_info_array(info, MIO_MDPA_MESH_BLOCKS, 0, 0, &data, &dt, &ndim, shape),
              MIO_OK);
    ASSERT_EQ(shape[0], 2);
    EXPECT_EQ(static_cast<const int64_t*>(data)[1], 3);  // node 4 -> row 3

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_SUBMODELPARTS), 1);
    EXPECT_EQ(capi_info_string(info, MIO_MDPA_SUBMODELPARTS, 0, 0), "Inlet");
    EXPECT_EQ(mio_mdpa_info_data_kind(info, MIO_MDPA_SUBMODELPARTS, 0, 0), MIO_MDPA_VALUE_NUMBER);
    for (int32_t field : {1, 2}) {
        ASSERT_EQ(
            mio_mdpa_info_array(info, MIO_MDPA_SUBMODELPARTS, 0, field, &data, &dt, &ndim, shape),
            MIO_OK);
        EXPECT_EQ(dt, MIO_INT64);
        EXPECT_EQ(ndim, 1);
        ASSERT_EQ(shape[0], 1);
        EXPECT_EQ(static_cast<const int64_t*>(data)[0], field == 1 ? 7 : 1);
    }

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_PROPERTIES), 1);
    ASSERT_EQ(mio_mdpa_info_data_count(info, MIO_MDPA_PROPERTIES, 0), 2);
    char key[32];
    mio_mdpa_info_data_string(info, MIO_MDPA_PROPERTIES, 0, 1, 0, key, sizeof(key));
    EXPECT_STREQ(key, "CONSTITUTIVE_LAW");

    ASSERT_EQ(mio_mdpa_info_count(info, MIO_MDPA_RAW_BLOCKS), 1);
    EXPECT_EQ(capi_info_string(info, MIO_MDPA_RAW_BLOCKS, 0, 2), "End Constraints");

    // Bad section/field/index fail without crossing an exception.
    EXPECT_EQ(mio_mdpa_info_count(info, 99), -1);
    EXPECT_EQ(mio_mdpa_info_string(info, MIO_MDPA_TABLES, 5, 0, nullptr, 0), -1);
    EXPECT_NE(mio_mdpa_info_int(info, MIO_MDPA_TABLES, 0, 0, &id), MIO_OK);

    // Round trip: the side channel survives mio_write_with_info.
    const std::string out = mt::temp_path("_capi_info_out.mdpa");
    ASSERT_EQ(mio_write_with_info(out.c_str(), mesh, nullptr, info), MIO_OK) << mio_last_error();
    // ...but not as another format.
    const std::string vtu = mt::temp_path("_capi_info_out.vtu");
    EXPECT_EQ(mio_write_with_info(vtu.c_str(), mesh, nullptr, info), MIO_ERR_INVALID_ARG);
    mio_format_info* info2 = nullptr;
    mio_mesh* mesh2 = mio_read_with_info(out.c_str(), nullptr, nullptr, &info2);
    ASSERT_NE(mesh2, nullptr) << mio_last_error();
    for (int32_t field : {1, 2}) {
        ASSERT_EQ(
            mio_mdpa_info_array(info2, MIO_MDPA_SUBMODELPARTS, 0, field, &data, &dt, &ndim, shape),
            MIO_OK);
        ASSERT_EQ(shape[0], 1);
        EXPECT_EQ(static_cast<const int64_t*>(data)[0], field == 1 ? 7 : 1);
    }
    for (int32_t section : {MIO_MDPA_TABLES, MIO_MDPA_GEOMETRIES, MIO_MDPA_MESH_BLOCKS,
                            MIO_MDPA_SUBMODELPARTS, MIO_MDPA_RAW_BLOCKS, MIO_MDPA_PROPERTIES})
        EXPECT_EQ(mio_mdpa_info_count(info2, section), mio_mdpa_info_count(info, section))
            << section;
    EXPECT_EQ(capi_info_string(info2, MIO_MDPA_RAW_BLOCKS, 0, 1),
              capi_info_string(info, MIO_MDPA_RAW_BLOCKS, 0, 1));

    // A format with no side channel reads normally and hands back no handle.
    mio_format_info* none = reinterpret_cast<mio_format_info*>(0x1);
    ASSERT_EQ(mio_write(vtu.c_str(), mesh, nullptr), MIO_OK);
    mio_mesh* mesh3 = mio_read_with_info(vtu.c_str(), nullptr, nullptr, &none);
    ASSERT_NE(mesh3, nullptr) << mio_last_error();
    EXPECT_EQ(none, nullptr);

    mio_mesh_free(mesh3);
    mio_mesh_free(mesh2);
    mio_mesh_free(mesh);
    mio_format_info_free(info2);
    mio_format_info_free(info);
    mio_format_info_free(nullptr);
    std::error_code ec;
    for (const std::string& p : {in, out, vtu})
        std::filesystem::remove(p, ec);
}

// --- software rendering (v16.33.0, field rendering v16.34.0) ------------------

TEST(CApi, RenderDefaultsAreInitializedNotZero) {
    mio_render_opts o;
    mio_render_opts_init(&o);
    EXPECT_EQ(o.width, 320);
    EXPECT_EQ(o.height, 240);
    EXPECT_EQ(o.shading, MIO_SHADING_FLAT);
    EXPECT_EQ(o.two_sided, 1);
    EXPECT_EQ(o.component, -1);
    EXPECT_EQ(o.fill_color, 0xC8C5BDFFu);
    EXPECT_DOUBLE_EQ(o.elevation, 35.264389682754654);
    mio_text_opts t;
    mio_text_opts_init(&t);
    EXPECT_EQ(t.cols, 80);
    EXPECT_EQ(t.notes, 1);
    // An all-zero struct is not a default: its zero width is refused by name.
    mio_render_opts zero{};
    mio_mesh* cube = capi_cube_surface();
    EXPECT_EQ(mio_render(cube, &zero), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("width and height"), std::string::npos);
    mio_mesh_free(cube);
}

TEST(CApi, RenderAFrameAndEncodeIt) {
    mio_mesh* cube = capi_cube_surface();
    mio_render_opts o;
    mio_render_opts_init(&o);
    o.width = 48;
    o.height = 32;
    o.background = 0xFFFFFFFFu;
    mio_frame* f = mio_render(cube, &o);
    ASSERT_NE(f, nullptr) << mio_last_error();
    EXPECT_EQ(mio_frame_width(f), 48);
    EXPECT_EQ(mio_frame_height(f), 32);
    ASSERT_NE(mio_frame_rgba(f), nullptr);
    ASSERT_NE(mio_frame_cell_ids(f), nullptr);
    EXPECT_EQ(mio_frame_range(f, nullptr, nullptr), 0);
    std::size_t covered = 0;
    for (int i = 0; i < 48 * 32; ++i)
        covered += mio_frame_cell_ids(f)[i] >= 0 ? 1 : 0;
    EXPECT_GT(covered, 100u);

    // PNG: length query, then the bytes.
    const std::int64_t need = mio_frame_png(f, 0, nullptr, 0);
    ASSERT_GT(need, 8);
    std::vector<std::uint8_t> png(static_cast<std::size_t>(need));
    EXPECT_EQ(mio_frame_png(f, 0, png.data(), need), need);
    EXPECT_EQ(std::string(png.begin(), png.begin() + 8), std::string("\x89PNG\r\n\x1a\n", 8));
    EXPECT_EQ(mio_frame_png(f, 10, nullptr, 0), -1);

    // Text: a cell encoding needs a frame sized to its grid (48 x 32 suits half blocks).
    mio_text_opts t;
    mio_text_opts_init(&t);
    t.color_depth = MIO_COLOR_MONO;
    const std::int64_t tn = mio_frame_text(f, &t, nullptr, 0);
    ASSERT_GT(tn, 0);
    std::string text(static_cast<std::size_t>(tn) + 1, '\0');
    EXPECT_EQ(mio_frame_text(f, &t, text.data(), tn + 1), tn);
    text.resize(static_cast<std::size_t>(tn));
    EXPECT_EQ(std::count(text.begin(), text.end(), '\n'), 16);  // 32 rows of pixels = 16 cells
    mio_frame_free(f);
    mio_frame_free(nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, RenderTextIsSizedToTheCells) {
    mio_mesh* cube = capi_cube_surface();
    mio_text_opts t;
    mio_text_opts_init(&t);
    t.cols = 20;
    t.rows = 6;
    t.color_depth = MIO_COLOR_MONO;
    const std::int64_t n = mio_render_text(cube, nullptr, &t, nullptr, 0);
    ASSERT_GT(n, 0) << mio_last_error();
    std::string out(static_cast<std::size_t>(n) + 1, '\0');
    EXPECT_EQ(mio_render_text(cube, nullptr, &t, out.data(), n + 1), n);
    out.resize(static_cast<std::size_t>(n));
    EXPECT_EQ(std::count(out.begin(), out.end(), '\n'), 6);
    EXPECT_EQ(mio_render_text(nullptr, nullptr, &t, nullptr, 0), -1);
    mio_mesh_free(cube);
}

TEST(CApi, RenderFieldOptionsAndNotes) {
    mio_mesh* m = mio_mesh_create();
    const double points[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0};
    const std::int64_t quad[] = {0, 1, 2, 3};
    const double u[] = {0, 1, 2, 3};
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, 4, 3, points), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", 1, 4, MIO_INT64, quad), MIO_OK);
    const std::int64_t u_shape[] = {4};
    ASSERT_EQ(mio_mesh_add_point_data(m, "u", MIO_FLOAT64, 1, u_shape, u), MIO_OK);
    mio_render_opts o;
    mio_render_opts_init(&o);
    o.width = 40;
    o.height = 40;
    o.view = "+z";
    o.color_by = "u";
    o.colorbar = 1;
    o.symmetric = 1;
    const double levels[] = {1.0, 2.0};
    o.iso_levels = levels;
    o.num_iso_levels = 2;
    mio_frame* f = mio_render(m, &o);
    ASSERT_NE(f, nullptr) << mio_last_error();
    double lo = 0, hi = 0;
    EXPECT_EQ(mio_frame_range(f, &lo, &hi), 1);
    EXPECT_DOUBLE_EQ(lo, -hi);
    const std::int64_t notes = mio_frame_num_notes(f);
    ASSERT_GE(notes, 2);
    char buf[256];
    EXPECT_GT(mio_frame_note(f, 0, buf, sizeof buf), 0);
    EXPECT_NE(std::string(buf).find("u: "), std::string::npos);
    EXPECT_EQ(mio_frame_note(f, notes, buf, sizeof buf), -1);
    mio_frame_free(f);

    // A bad enum, a reserved word set and a misnamed array are refused by name.
    o.num_iso_levels = 0;
    o.iso_levels = nullptr;
    o.scale = 9;
    EXPECT_EQ(mio_render(m, &o), nullptr);
    o.scale = MIO_SCALE_LINEAR;
    o.reserved[0] = 1;
    EXPECT_EQ(mio_render(m, &o), nullptr);
    o.reserved[0] = 0;
    o.color_by = "nope";
    EXPECT_EQ(mio_render(m, &o), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("nope"), std::string::npos);
    mio_mesh_free(m);
}

TEST(CApi, RenderCutawaysClipAndAreValidated) {
    mio_render_opts o;
    mio_render_opts_init(&o);
    EXPECT_EQ(o.num_cutaways, 0);
    EXPECT_EQ(o.cutaways, nullptr);
    EXPECT_EQ(o.cutaway_tint, 0xE88034FFu);
    mio_mesh* cube = capi_cube_surface();
    o.width = 32;
    o.height = 32;
    o.view = "-x";  // the camera sits on the -x side
    o.shading = MIO_SHADING_NONE;
    o.background = 0x000000FFu;
    mio_frame* whole = mio_render(cube, &o);
    ASSERT_NE(whole, nullptr) << mio_last_error();
    const std::int64_t whole_centre = mio_frame_cell_ids(whole)[16 * 32 + 16];
    const std::vector<std::uint8_t> whole_rgba(mio_frame_rgba(whole),
                                               mio_frame_rgba(whole) + 32 * 32 * 4);
    mio_frame_free(whole);

    // One plane keeping x >= 0.5: the near face goes, the far one shows from inside.
    const double plane[6] = {0.5, 0, 0, 1, 0, 0};
    o.cutaways = plane;
    o.num_cutaways = 1;
    o.cutaway_tint = 0x0AC81EFFu;
    mio_frame* cut = mio_render(cube, &o);
    ASSERT_NE(cut, nullptr) << mio_last_error();
    EXPECT_NE(mio_frame_cell_ids(cut)[16 * 32 + 16], whole_centre);
    const std::uint8_t* centre = mio_frame_rgba(cut) + (16 * 32 + 16) * 4;
    EXPECT_EQ(centre[0], 0x0A);
    EXPECT_EQ(centre[1], 0xC8);
    EXPECT_EQ(centre[2], 0x1E);
    EXPECT_NE(std::vector<std::uint8_t>(mio_frame_rgba(cut), mio_frame_rgba(cut) + 32 * 32 * 4),
              whole_rgba);
    mio_frame_free(cut);

    // A pointer without a count, a count over two, a negative one, a zero normal,
    // and a non-zero reserved word are refused by name.
    o.cutaways = nullptr;
    EXPECT_EQ(mio_render(cube, &o), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("cutaways"), std::string::npos);
    const double three[18] = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
    o.cutaways = three;
    o.num_cutaways = 3;
    EXPECT_EQ(mio_render(cube, &o), nullptr);
    o.num_cutaways = -1;
    EXPECT_EQ(mio_render(cube, &o), nullptr);
    const double flat[6] = {0, 0, 0, 0, 0, 0};
    o.cutaways = flat;
    o.num_cutaways = 1;
    EXPECT_EQ(mio_render(cube, &o), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("normal"), std::string::npos);
    o.cutaways = plane;
    o.reserved[0] = 1;
    EXPECT_EQ(mio_render(cube, &o), nullptr);
    mio_mesh_free(cube);
}

TEST(CApi, RenderStreamlinesDrawAndAreValidated) {
    mio_render_opts o;
    mio_render_opts_init(&o);
    EXPECT_EQ(o.streamlines, nullptr);
    EXPECT_EQ(o.stream_seeds, 40);
    EXPECT_EQ(o.stream_length, 0.5);
    EXPECT_EQ(o.stream_color, 0xF050A0FFu);
    EXPECT_EQ(sizeof(mio_render_opts), 432u);
    // A 6 x 6 grid of quads with a swirling vector point array.
    const int n = 6;
    std::vector<double> pts;
    std::vector<double> swirl;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) {
            const double x = double(i) / n;
            const double y = double(j) / n;
            pts.insert(pts.end(), {x, y, 0.0});
            swirl.insert(swirl.end(), {-(y - 0.5), x - 0.5, 0.0});
        }
    std::vector<std::int64_t> quads;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const std::int64_t p = j * (n + 1) + i;
            quads.insert(quads.end(), {p, p + 1, p + n + 2, p + n + 1});
        }
    mio_mesh* m = mio_mesh_create();
    ASSERT_NE(m, nullptr);
    ASSERT_EQ(mio_mesh_set_points(m, MIO_FLOAT64, (n + 1) * (n + 1), 3, pts.data()), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(m, "quad", n * n, 4, MIO_INT64, quads.data()), MIO_OK);
    const std::int64_t shape[2] = {(n + 1) * (n + 1), 3};
    ASSERT_EQ(mio_mesh_add_point_data(m, "swirl", MIO_FLOAT64, 2, shape, swirl.data()), MIO_OK);
    o.width = 96;
    o.height = 96;
    o.view = "+z";
    o.shading = MIO_SHADING_NONE;
    mio_frame* plain = mio_render(m, &o);
    ASSERT_NE(plain, nullptr) << mio_last_error();
    auto count = [](mio_frame* f, std::uint32_t rgba) {
        std::size_t c = 0;
        const std::uint8_t* px = mio_frame_rgba(f);
        for (std::size_t i = 0; i < 96u * 96u; ++i)
            c += px[4 * i] == (rgba >> 24) && px[4 * i + 1] == ((rgba >> 16) & 255) &&
                 px[4 * i + 2] == ((rgba >> 8) & 255);
        return c;
    };
    EXPECT_EQ(count(plain, 0xFF0000FFu), 0u);
    mio_frame_free(plain);
    o.streamlines = "swirl";
    o.stream_color = 0xFF0000FFu;
    mio_frame* lines = mio_render(m, &o);
    ASSERT_NE(lines, nullptr) << mio_last_error();
    EXPECT_GT(count(lines, 0xFF0000FFu), 100u);
    mio_frame_free(lines);
    o.streamlines = "nope";
    EXPECT_EQ(mio_render(m, &o), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("nope"), std::string::npos);
    o.streamlines = "swirl";
    o.stream_seeds = 0;
    EXPECT_EQ(mio_render(m, &o), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("seeds"), std::string::npos);
    mio_mesh_free(m);
}

TEST(CApi, WriteSnapshotChoosesTheFormByExtension) {
    mio_mesh* cube = capi_cube_surface();
    const std::string png = mt::temp_path(".png");
    const std::string txt = mt::temp_path(".txt");
    mio_render_opts o;
    mio_render_opts_init(&o);
    o.width = 32;
    o.height = 24;
    EXPECT_EQ(mio_write_snapshot(png.c_str(), cube, &o, nullptr, nullptr), MIO_OK)
        << mio_last_error();
    mio_text_opts t;
    mio_text_opts_init(&t);
    t.cols = 12;
    t.rows = 4;
    EXPECT_EQ(mio_write_snapshot(txt.c_str(), cube, nullptr, &t, nullptr), MIO_OK);
    std::ifstream in(png, std::ios::binary);
    char head[8];
    in.read(head, 8);
    EXPECT_EQ(std::string(head, 8), std::string("\x89PNG\r\n\x1a\n", 8));
    EXPECT_NE(mio_write_snapshot("frame.vtu", cube, nullptr, nullptr, nullptr), MIO_OK);
    mio_snapshot_opts s;
    mio_snapshot_opts_init(&s);
    EXPECT_EQ(s.cast_frames, 36);
    EXPECT_NE(mio_write_snapshot(nullptr, cube, nullptr, nullptr, nullptr), MIO_OK);
    mio_mesh_free(cube);
}

TEST(CApi, DetectColorDepth) {
    EXPECT_EQ(mio_detect_color_depth("1", "truecolor", "xterm-256color"), MIO_COLOR_MONO);
    EXPECT_EQ(mio_detect_color_depth(nullptr, "24bit", nullptr), MIO_COLOR_TRUECOLOR);
    EXPECT_EQ(mio_detect_color_depth(nullptr, nullptr, "xterm-256color"), MIO_COLOR_256);
    EXPECT_EQ(mio_detect_color_depth(nullptr, nullptr, "xterm"), MIO_COLOR_16);
}

// mio_fit_cad3 (v17.12.0): the freeform-mode options and report. The legacy
// default reports a facet area of 0 on the cube's six planes and no table;
// another mode adds the "fit_cad:face_table" field data; an option out of
// range is refused, naming it. mio_fit_cad2 forwards with the defaults.
TEST(CApi, FitCad3FreeformModes) {
    mio_mesh* m = capi_cube_surface();
    mio_fit_cad_opts opts;
    mio_fit_cad_opts_init(&opts);
    opts.backend = MIO_FIT_CAD_NATIVE;
    mio_fit_cad_opts_ex2 ex2;
    mio_fit_cad_opts_ex2_init(&ex2);
    EXPECT_EQ(ex2.freeform_mode, MIO_FREEFORM_LEGACY);
    EXPECT_EQ(ex2.extraordinary_vertex, MIO_EV_G1X);
    EXPECT_EQ(ex2.freeform_max_control_points, 400);
    EXPECT_EQ(ex2.layout_snap_to_features, 1);
    EXPECT_EQ(ex2.constrain_primitives, 0);
    mio_fit_cad_report_ex3 report_ex3;
    std::memset(&report_ex3, 0xff, sizeof report_ex3);
    mio_mesh* skin =
        mio_fit_cad3(m, nullptr, &opts, nullptr, &ex2, nullptr, nullptr, nullptr, &report_ex3);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report_ex3.brep_checked, 0);
    EXPECT_EQ(report_ex3.freeform_mode, MIO_FREEFORM_LEGACY);
    EXPECT_EQ(report_ex3.facet_area_fraction, 0.0);
    EXPECT_EQ(report_ex3.num_face_records, 6);
    EXPECT_EQ(report_ex3.num_constrained_primitives, 0);
    EXPECT_EQ(report_ex3.num_constraint_groups, 0);
    for (int64_t r : report_ex3.reserved)
        EXPECT_EQ(r, 0);
    const void* data = nullptr;
    mio_dtype dt;
    int32_t ndim = 0;
    int64_t shape[8] = {};
    EXPECT_NE(mio_mesh_get_field_data(skin, "fit_cad:face_table", &data, &dt, &ndim, shape),
              MIO_OK);
    mio_mesh_free(skin);

    ex2.freeform_mode = MIO_FREEFORM_ADAPTIVE;
    skin = mio_fit_cad3(m, nullptr, &opts, nullptr, &ex2, nullptr, nullptr, nullptr, &report_ex3);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report_ex3.freeform_mode, MIO_FREEFORM_ADAPTIVE);
    ASSERT_EQ(mio_mesh_get_field_data(skin, "fit_cad:face_table", &data, &dt, &ndim, shape), MIO_OK)
        << mio_last_error();
    EXPECT_EQ(dt, MIO_FLOAT64);
    ASSERT_EQ(ndim, 2);
    EXPECT_EQ(shape[0], 6);
    EXPECT_EQ(shape[1], 13);
    mio_mesh_free(skin);

    ex2.freeform_mode = 9;
    EXPECT_EQ(mio_fit_cad3(m, nullptr, &opts, nullptr, &ex2, nullptr, nullptr, nullptr, nullptr),
              nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("freeform_mode"), std::string::npos)
        << mio_last_error();
    // Zero reads as the default for the tolerances and budgets.
    mio_fit_cad_opts_ex2_init(&ex2);
    ex2.freeform_g1_tolerance_deg = 0.0;
    ex2.freeform_max_control_points = 0;
    ex2.layout_target_patch_size = 0.0;
    ex2.layout_max_valence_deviation = 0;
    skin = mio_fit_cad3(m, nullptr, &opts, nullptr, &ex2, nullptr, nullptr, nullptr, nullptr);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    mio_mesh_free(skin);
    // The constrained primitives: the cube's opposite planes share a normal
    // and its neighbours are made exactly perpendicular, every face kept.
    mio_fit_cad_opts_ex2_init(&ex2);
    ex2.constrain_primitives = 1;
    mio_fit_cad_report report;
    skin = mio_fit_cad3(m, nullptr, &opts, nullptr, &ex2, &report, nullptr, nullptr, &report_ex3);
    ASSERT_NE(skin, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_demoted, 0);
    EXPECT_GE(report_ex3.num_constrained_primitives, 0);
    EXPECT_LE(report_ex3.num_constrained_primitives, 6);
    for (int64_t r : report_ex3.reserved)
        EXPECT_EQ(r, 0);
    mio_mesh_free(skin);
    mio_mesh_free(m);
}

TEST(CApi, ExtrudeRevolveAndMirror) {
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = 3;
    g.counts[1] = 2;
    g.cell_type = MIO_GEN_QUAD;
    mio_mesh* rect = mio_rectangle(&g);
    ASSERT_NE(rect, nullptr) << mio_last_error();

    mio_extrude_opts e;
    mio_extrude_opts_init(&e);
    EXPECT_EQ(e.direction[2], 1.0);
    EXPECT_EQ(e.num_layers, 1);
    EXPECT_EQ(e.ratio, 1.0);
    e.num_layers = 4;
    mio_sweep_report report{};
    report.reserved[0] = 7;
    mio_mesh* out = mio_extrude(rect, &e, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_hexahedra, 24);
    EXPECT_EQ(report.num_cells, 24);
    EXPECT_EQ(report.full, 0);
    EXPECT_EQ(report.reserved[0], 0);
    EXPECT_EQ(mio_mesh_num_points(out), 4 * 3 * 5);
    mio_mesh_free(out);
    // Explicit offsets; NULL options are the defaults.
    const double offsets[2] = {0.25, 1.0};
    mio_extrude_opts_init(&e);
    e.offsets = offsets;
    e.num_offsets = 2;
    out = mio_extrude(rect, &e, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.num_cells, 12);
    mio_mesh_free(out);
    out = mio_extrude(rect, nullptr, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    mio_mesh_free(out);
    mio_extrude_opts_init(&e);
    e.direction[1] = 1.0;
    e.direction[2] = 0.0;
    EXPECT_EQ(mio_extrude(rect, &e, nullptr), nullptr);
    EXPECT_NE(std::string(mio_last_error()).find("parallel"), std::string::npos);

    mio_revolve_opts r;
    mio_revolve_opts_init(&r);
    EXPECT_EQ(r.angle, 360.0);
    EXPECT_EQ(r.num_layers, 16);
    EXPECT_EQ(r.axis[1], 1.0);
    r.num_layers = 8;
    out = mio_revolve(rect, &r, &report);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(report.full, 1);
    EXPECT_EQ(report.num_on_axis, 3);
    EXPECT_EQ(report.num_wedges, 2 * 8);
    EXPECT_EQ(report.num_warnings, 0);
    EXPECT_EQ(mio_revolve_warnings(nullptr, 0), 0);
    mio_mesh_free(out);
    r.num_layers = 2;
    EXPECT_EQ(mio_revolve(rect, &r, nullptr), nullptr);
    mio_mesh_free(rect);

    mio_gen_opts_init(&g);
    g.counts[0] = g.counts[1] = g.counts[2] = 2;
    g.cell_type = MIO_GEN_TETRA;
    mio_mesh* box = mio_box(&g);
    ASSERT_NE(box, nullptr) << mio_last_error();
    mio_mirror_opts mo;
    mio_mirror_opts_init(&mo);
    EXPECT_EQ(mo.normal[0], 1.0);
    EXPECT_EQ(mo.weld, 1);
    EXPECT_EQ(mo.keep_original, 0);
    mio_mirror_report mr{};
    out = mio_mirror(box, &mo, &mr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(mr.num_flipped_cells, 48);
    std::int64_t nc = 0, ninv = 0, ndeg = 0;
    ASSERT_EQ(mio_quality_counts(out, &nc, &ninv, &ndeg), MIO_OK);
    EXPECT_EQ(ninv, 0);
    mio_mesh_free(out);
    mo.keep_original = 1;
    out = mio_mirror(box, &mo, &mr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    EXPECT_EQ(mr.num_points_on_plane, 9);
    EXPECT_EQ(mio_mesh_num_points(out), 2 * 27 - 9);
    EXPECT_EQ(mr.num_warnings, 1);
    const std::int64_t len = mio_mirror_warnings(nullptr, 0);
    ASSERT_GT(len, 0);
    std::string text(static_cast<std::size_t>(len) + 1, '\0');
    EXPECT_EQ(mio_mirror_warnings(text.data(), len + 1), len);
    EXPECT_NE(text.find("xmin"), std::string::npos);
    mio_mesh_free(out);
    mo.normal[0] = 0.0;
    EXPECT_EQ(mio_mirror(box, &mo, nullptr), nullptr);
    mio_mesh_free(box);
}

// The MMG backend (roadmap 3.2): layouts, the capability queries, mio_adapt
// with the backend or its refusal by name without it, and the MMG-free
// mio_metric_from_hessian.
TEST(CApi, MmgOptionLayoutsAndDefaults) {
    EXPECT_EQ(sizeof(mio_mmg_local), 32u);
    EXPECT_EQ(sizeof(mio_mmg_params), 168u);
    EXPECT_EQ(sizeof(mio_mmg_report), 192u);
    EXPECT_EQ(sizeof(mio_adapt_opts), 224u);
    EXPECT_EQ(offsetof(mio_adapt_opts, hsiz), 168u);
    mio_adapt_opts o;
    mio_adapt_opts_init(&o);
    EXPECT_EQ(o.params.hmin, -1.0);
    EXPECT_EQ(o.params.angle_detection, 45.0);
    EXPECT_EQ(o.params.verbose, -1);
    EXPECT_EQ(o.hsiz, -1.0);
    EXPECT_EQ(mio_mmg_index_bits() != 0, mio_has_mmg() != 0);
}

TEST(CApi, AdaptWithOrWithoutMmg) {
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = 4;
    g.counts[1] = 4;
    g.cell_type = MIO_GEN_TRIANGLE;
    mio_mesh* square = mio_rectangle(&g);
    ASSERT_NE(square, nullptr) << mio_last_error();
    mio_adapt_opts o;
    mio_adapt_opts_init(&o);
    o.hsiz = 0.1;
    const int64_t frozen[] = {12};
    o.params.frozen = frozen;
    o.params.num_frozen = 1;
    mio_mmg_report report;
    mio_mesh* out = mio_adapt(square, &o, &report);
    if (!mio_has_mmg()) {
        EXPECT_EQ(out, nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_MMG"), std::string::npos);
    } else {
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_EQ(report.status, 0);
        EXPECT_STREQ(report.library, "mmg2d");
        EXPECT_EQ(report.num_points_in, 25);
        EXPECT_EQ(report.num_points_out, mio_mesh_num_points(out));
        EXPECT_GT(report.num_elements_out, report.num_elements_in);
        EXPECT_GE(report.quality_min, 1.0);
        mio_mesh_free(out);
        // NULL options and report are fine (MMG sizes from the mesh).
        mio_mesh* plain = mio_adapt(square, nullptr, nullptr);
        ASSERT_NE(plain, nullptr) << mio_last_error();
        mio_mesh_free(plain);
        // A frozen id out of range is refused.
        const int64_t bad[] = {999};
        o.params.frozen = bad;
        EXPECT_EQ(mio_adapt(square, &o, nullptr), nullptr);
    }
    mio_mesh_free(square);
}

TEST(CApi, MetricFromHessian) {
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = 2;
    g.counts[1] = 2;
    g.cell_type = MIO_GEN_TRIANGLE;
    mio_mesh* m = mio_rectangle(&g);
    ASSERT_NE(m, nullptr) << mio_last_error();
    const int64_t n = mio_mesh_num_points(m);
    std::vector<double> h(static_cast<std::size_t>(n) * 4, 0.0);
    for (int64_t i = 0; i < n; ++i)
        h[static_cast<std::size_t>(i) * 4] = 1e6;
    const int64_t shape[] = {n, 4};
    ASSERT_EQ(mio_mesh_add_point_data(m, "h", MIO_FLOAT64, 2, shape, h.data()), MIO_OK);
    mio_mesh* out = mio_metric_from_hessian(m, "h", 1.0, 0.1, 1.0, nullptr);
    ASSERT_NE(out, nullptr) << mio_last_error();
    const void* data = nullptr;
    int64_t got[MIO_MAX_NDIM] = {0};
    int32_t ndim = 0;
    ASSERT_EQ(mio_mesh_get_point_data(out, "metric", &data, nullptr, &ndim, got), MIO_OK);
    EXPECT_EQ(ndim, 2);
    EXPECT_EQ(got[1], 3);
    const double* metric = static_cast<const double*>(data);
    EXPECT_DOUBLE_EQ(metric[0], 100.0);  // clamped to 1/hmin^2
    EXPECT_DOUBLE_EQ(metric[2], 1.0);    // clamped to 1/hmax^2
    mio_mesh_free(out);
    EXPECT_EQ(mio_metric_from_hessian(m, "nope", 1.0, 0.1, 1.0, nullptr), nullptr);
    mio_mesh_free(m);
}

TEST(CApi, DiscretizeLevelsetAndMeshDomain) {
    EXPECT_EQ(sizeof(mio_mmg_material), 32u);
    EXPECT_EQ(sizeof(mio_discretize_levelset_opts), 240u);
    EXPECT_EQ(sizeof(mio_mesh_domain_opts), 232u);
    // A disc's level set in a square.
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = 8;
    g.counts[1] = 8;
    g.size[0] = g.size[1] = 2.0;
    g.origin[0] = g.origin[1] = -1.0;
    g.cell_type = MIO_GEN_TRIANGLE;
    mio_mesh* square = mio_rectangle(&g);
    ASSERT_NE(square, nullptr) << mio_last_error();
    const int64_t n = mio_mesh_num_points(square);
    const void* pts = nullptr;
    ASSERT_EQ(mio_mesh_get_points(square, &pts, nullptr), MIO_OK);
    std::vector<double> ls(static_cast<std::size_t>(n));
    for (int64_t i = 0; i < n; ++i) {
        const double* p = static_cast<const double*>(pts) + i * 3;
        ls[static_cast<std::size_t>(i)] = std::sqrt(p[0] * p[0] + p[1] * p[1]) - 0.5;
    }
    const int64_t shape[] = {n};
    ASSERT_EQ(mio_mesh_add_point_data(square, "ls", MIO_FLOAT64, 1, shape, ls.data()), MIO_OK);
    mio_discretize_levelset_opts o;
    mio_discretize_levelset_opts_init(&o);
    o.level_set = "ls";
    mio_mmg_report report;
    mio_mesh* cut = mio_discretize_levelset(square, &o, &report);
    if (!mio_has_mmg()) {
        EXPECT_EQ(cut, nullptr);
        EXPECT_NE(std::string(mio_last_error()).find("MESHIOPLUSPLUS_WITH_MMG"), std::string::npos);
    } else {
        ASSERT_NE(cut, nullptr) << mio_last_error();
        EXPECT_STREQ(report.library, "mmg2d");
        EXPECT_GT(report.num_elements_out, report.num_elements_in);
        mio_mesh_free(cut);
    }
    EXPECT_EQ(mio_discretize_levelset(square, nullptr, nullptr), nullptr);
    mio_mesh_free(square);

    // A square's boundary, meshed.
    const double xy[] = {0, 0, 1, 0, 1, 1, 0, 1};
    const int64_t lines[] = {0, 1, 1, 2, 2, 3, 3, 0};
    mio_mesh* boundary = mio_mesh_create();
    ASSERT_EQ(mio_mesh_set_points(boundary, MIO_FLOAT64, 4, 2, xy), MIO_OK);
    ASSERT_EQ(mio_mesh_add_cell_block(boundary, "line", 4, 2, MIO_INT64, lines), MIO_OK);
    mio_mesh_domain_opts d;
    mio_mesh_domain_opts_init(&d);
    d.hsiz = 0.1;
    int64_t subdomains = -1, holes = -1;
    mio_mesh* domain = mio_mesh_domain(boundary, &d, &report, &subdomains, &holes);
    if (!mio_has_mmg()) {
        EXPECT_EQ(domain, nullptr);
    } else {
        ASSERT_NE(domain, nullptr) << mio_last_error();
        EXPECT_EQ(subdomains, 1);
        EXPECT_EQ(holes, 0);
        EXPECT_GT(report.num_elements_out, 100);
        mio_mesh_free(domain);
    }
    mio_mesh_free(boundary);
}

TEST(CApi, MoveTranslatesAFaceOrNamesWhatIsMissing) {
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = g.counts[1] = g.counts[2] = 4;
    g.cell_type = MIO_GEN_TETRA;
    mio_mesh* box = mio_box(&g);
    ASSERT_NE(box, nullptr) << mio_last_error();
    const int64_t n = mio_mesh_num_points(box);
    const void* raw = nullptr;
    mio_dtype dt;
    ASSERT_EQ(mio_mesh_get_points(box, &raw, &dt), MIO_OK);
    ASSERT_EQ(dt, MIO_FLOAT64);
    const double* p = static_cast<const double*>(raw);
    std::vector<double> d(static_cast<std::size_t>(n) * 3, 0.0);
    for (int64_t i = 0; i < n; ++i)
        if (std::abs(p[i * 3] - 1.0) < 1e-12)
            d[static_cast<std::size_t>(i) * 3] = 0.1;
    const int64_t shape[2] = {n, 3};
    ASSERT_EQ(mio_mesh_add_point_data(box, "d", MIO_FLOAT64, 2, shape, d.data()), MIO_OK);
    mio_move_opts o;
    mio_move_opts_init(&o);
    EXPECT_EQ(o.lagrangian, 1);
    const char* moving[] = {"xmax"};
    o.displacement = "d";
    o.moving = moving;
    o.num_moving = 1;
    mio_mmg_report report;
    mio_mesh* out = mio_move(box, &o, &report);
    if (!mio_has_mmg() || !mio_has_mmg_elas()) {
        EXPECT_EQ(out, nullptr);
        EXPECT_NE(
            std::string(mio_last_error()).find(mio_has_mmg() ? "ELAS" : "MESHIOPLUSPLUS_WITH_MMG"),
            std::string::npos);
    } else {
        ASSERT_NE(out, nullptr) << mio_last_error();
        EXPECT_STREQ(report.library, "mmg3d");
        const void* q = nullptr;
        ASSERT_EQ(mio_mesh_get_points(out, &q, &dt), MIO_OK);
        double xmax = 0.0;
        for (int64_t i = 0; i < mio_mesh_num_points(out); ++i)
            xmax = std::max(xmax, static_cast<const double*>(q)[i * 3]);
        EXPECT_NEAR(xmax, 1.1, 1e-9);
        mio_mesh_free(out);
    }
    EXPECT_EQ(mio_move(box, nullptr, nullptr), nullptr);
    mio_mesh_free(box);
}
