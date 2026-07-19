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
 * @brief Tests for the C API (bindings_c/). Compiled into the gtest suite
 *        only when MESHIOPLUSPLUS_BUILD_C_API=ON; written purely against the
 *        public C surface (plus mt:: fixtures for reference meshes), so it
 *        runs identically under every mesh backend.
 */

// System includes
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/meshioplusplus.h"

#include "mesh_fixtures.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/med.hpp"
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

TEST(CApi, FormatAvailability) {
    EXPECT_EQ(mio_format_readable("vtu"), 1);
    EXPECT_EQ(mio_format_writable("vtu"), 1);
    EXPECT_EQ(mio_format_readable("openfoam"), 1);
    EXPECT_EQ(mio_format_writable("openfoam"), 0);  // read-only format
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
    // openfoam is read-only: resolvable format, no writer.
    EXPECT_EQ(mio_write("mesh.foam", m, "openfoam"), MIO_ERR_NOT_FOUND);

#ifndef MESHIOPLUSPLUS_HAS_HDF5
    // Compiled-out formats name the missing dependency.
    EXPECT_EQ(mio_write("mesh.med", m, nullptr), MIO_ERR_NOT_FOUND);
    EXPECT_NE(std::strstr(mio_last_error(), "HDF5"), nullptr) << mio_last_error();
#endif

    mio_mesh_free(m);
    mio_mesh_free(nullptr);  // NULL-safe
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

#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(CApi, RaggedBlocksAreReportedButNotAccessible) {
    // Ragged blocks cannot be constructed through the C API, and MED is the
    // one C++ writer that serializes them (POG polygons) -- build the mesh
    // through the C++ API, round-trip through .med, inspect via C.
    meshioplusplus::Mesh cpp_mesh;
    cpp_mesh.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0.5, 0}}));
    cpp_mesh.AddPolygonBlock("polygon", {{0, 1, 2, 3}, {1, 4, 2}});
    const std::string med = mt::temp_path("_capi_ragged.med");
    meshioplusplus::write_med(med, cpp_mesh, meshioplusplus::MedInfo{});

    mio_mesh* m = mio_read(med.c_str(), nullptr);
    ASSERT_NE(m, nullptr) << mio_last_error();
    ASSERT_EQ(mio_mesh_num_cell_blocks(m), 1);
    std::int64_t num_cells = 0, npc = -1;
    std::int32_t ragged = 0;
    ASSERT_EQ(mio_mesh_cell_block_info(m, 0, &num_cells, &npc, &ragged), MIO_OK);
    EXPECT_EQ(num_cells, 2);
    EXPECT_EQ(ragged, 1);
    const void* conn = nullptr;
    mio_dtype dt;
    EXPECT_EQ(mio_mesh_cell_block_conn(m, 0, &conn, &dt), MIO_ERR_UNSUPPORTED);
    EXPECT_STRNE(mio_last_error(), "");
    mio_mesh_free(m);
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

}  // namespace
