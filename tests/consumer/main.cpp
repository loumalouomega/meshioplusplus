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
 * @file tests/consumer/main.cpp
 * @brief Smoke test for the INSTALLED meshio++, built as a separate CMake
 * project against an install prefix with the source tree unreachable.
 *
 * This is the regression test for packaging, not for the library's behaviour --
 * the gtest suite covers behaviour. What it pins is that a consumer really can:
 *
 *  1. write and read a file through the registry (`registry_read` / the writer
 *     tables), i.e. the format layer is reachable and linked;
 *  2. call the operations layer (`clean`);
 *  3. include `kratos_bridge.hpp` and instantiate `to_model_part` /
 *     `from_model_part` against a Kratos-shaped model part -- the header exists
 *     precisely to serve an out-of-tree Kratos application, and it is
 *     unreachable through the C ABI, so nothing else would catch it regressing;
 *  4. link the C API into the SAME binary, proving the two components coexist
 *     under one find_package() without symbol clashes.
 *
 * The "Kratos-shaped model part" is meshio++'s own `meshioplusplus::ModelPart`
 * (backends/model_part.hpp), exactly as tests/cpp/test_kratos_bridge.cpp uses it:
 * the bridge is templated on a consumer-supplied type, so no real Kratos (or
 * CoSimIO) is needed to exercise it.
 */

// System includes
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// Project includes -- the full C++ API, from the install prefix
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/kratos_bridge.hpp"
#include "meshioplusplus/backends/model_part.hpp"
#include "meshioplusplus/operations/clean.hpp"

// ...and the C API, in the same translation unit
#include "meshioplusplus/meshioplusplus.h"

namespace {

int g_failures = 0;

void check(bool cond, const char* what) {
    std::printf("%s  %s\n", cond ? "  ok  " : "  FAIL", what);
    if (!cond) {
        ++g_failures;
    }
}

/// Two triangles sharing an edge, plus one deliberately duplicated cell for clean() to drop.
meshioplusplus::Mesh make_mesh() {
    using meshioplusplus::DType;
    using meshioplusplus::NDArray;

    const std::vector<double> pts{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0};
    NDArray points(DType::Float64, {4, 3});
    for (std::size_t i = 0; i < pts.size(); ++i) {
        points.As<double>()[i] = pts[i];
    }

    // Third row repeats the first: clean()'s duplicate-cell pass must remove it.
    const std::vector<std::int64_t> rows{0, 1, 2, 0, 2, 3, 0, 1, 2};
    NDArray conn(DType::Int64, {3, 3});
    for (std::size_t i = 0; i < rows.size(); ++i) {
        conn.As<std::int64_t>()[i] = rows[i];
    }

    meshioplusplus::Mesh m;
    m.AssignPoints(std::move(points));
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

}  // namespace

int main() {
    std::printf("meshio++ installed-consumer smoke test\n");
    std::printf("  mesh backend in use: %s\n",
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
                "NATIVE"
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
                "KRATOS"
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO)
                "MESHIO"
#else
                "<none defined -- the backend macro did NOT propagate>"
#endif
    );

    // The macro must arrive from the imported target's INTERFACE_COMPILE_DEFINITIONS,
    // never from this file. Without it mesh.hpp silently assumes MESHIO and the
    // link-time sentinel in detail/mesh_backend_check.hpp is what saves us.
#if !defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO) && \
    !defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE) && !defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
#error "no MESHIOPLUSPLUS_MESH_BACKEND_* macro: the exported target failed to propagate it"
#endif

    // ---------------------------------------------------------------- 1. registry
    const std::string path = "consumer_smoketest.vtu";
    {
        const meshioplusplus::Mesh m = make_mesh();
        const auto& writers = meshioplusplus::registry_writers();
        const auto it = writers.find("vtu");
        check(it != writers.end(), "registry_writers() knows 'vtu'");
        if (it != writers.end()) {
            it->second(path, m);
        }
    }
    {
        const meshioplusplus::Mesh m =
            meshioplusplus::registry_read(path, "vtu", meshioplusplus::ReadOptions{});
        check(m.NumPoints() == 4, "registry_read() round-trips 4 points");
        check(m.NumCellBlocks() == 1, "registry_read() round-trips 1 cell block");
    }
    check(!meshioplusplus::resolve_format(path, "").empty(),
          "resolve_format() infers a format from the extension");
    std::remove(path.c_str());

    // -------------------------------------------------------------- 2. operations
    {
        const meshioplusplus::Mesh m = make_mesh();
        meshioplusplus::CleanOptions opts;  // defaults drop degenerate + duplicate cells
        const meshioplusplus::CleanResult r = meshioplusplus::clean(m, opts);
        check(r.mCellsDroppedDuplicate == 1, "clean() drops the duplicated triangle");
        check(r.mMesh.NumPoints() == 4, "clean() keeps all four points");
    }

    // ----------------------------------------------------------- 3. kratos_bridge
    {
        meshioplusplus::ModelPart source("Main");
        source.CreateNewNode(1, 0.0, 0.0, 0.0);
        source.CreateNewNode(2, 1.0, 0.0, 0.0);
        source.CreateNewNode(3, 0.0, 1.0, 0.0);
        source.CreateNewElement(meshioplusplus::CellType::Triangle, 1, {1, 2, 3}, 0);

        // to_model_part() into a second Kratos-shaped part (the template
        // parameter is where a real Kratos::ModelPart would go).
        meshioplusplus::ModelPart dest("Dest");
        meshioplusplus::to_model_part(source, dest);
        check(dest.NumberOfNodes() == 3, "to_model_part() copies 3 nodes");
        check(dest.NumberOfElements() == 1, "to_model_part() copies 1 element");

        // ...and back out through the traits-based reader.
        const meshioplusplus::ModelPart round = meshioplusplus::from_model_part(dest, "Round");
        check(round.NumberOfNodes() == 3, "from_model_part() reads 3 nodes back");
        check(round.NumberOfElements() == 1, "from_model_part() reads 1 element back");
    }

    // ---------------------------------------------------------------- 4. C API
    {
        const char* v = mio_version();
        check(v != nullptr && v[0] != '\0', "mio_version() links from the same binary");

        mio_mesh* cm = mio_mesh_create();
        check(cm != nullptr, "mio_mesh_create() works alongside the C++ API");
        if (cm != nullptr) {
            const double coords[9] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
            check(mio_mesh_set_points(cm, MIO_FLOAT64, 3, 3, coords) == MIO_OK,
                  "mio_mesh_set_points() succeeds");
            check(mio_mesh_num_points(cm) == 3, "mio_mesh_num_points() reports 3");
            mio_mesh_free(cm);
        }
    }

    std::printf(g_failures == 0 ? "\nALL CHECKS PASSED\n" : "\n%d CHECK(S) FAILED\n", g_failures);
    return g_failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
