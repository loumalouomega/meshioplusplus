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
 * @file bench_backends.cpp
 * @brief Pure-C++ micro-benchmark comparing the mesh backends
 * (MESHIO / NATIVE / KRATOS) on the same synthetic workload.
 *
 * Because the mesh backend is an exclusive compile-time choice, one binary
 * measures one backend; `benchmark/bench_backends.sh` builds and runs all
 * three and collates the CSV. Method mirrors `benchmark/bench.py`: warmup
 * run + median of N timed runs (`std::chrono::steady_clock`), no external
 * benchmark framework.
 *
 * Workload: a structured tetrahedral cube (the C++ analogue of
 * `benchmark/inputs.py`'s `synthetic_tet_grid`; 6 tets per hex cell, shared
 * vertices). Timed operations, one CSV row each
 * (`backend,op,format,cells,median_s,runs` to stdout):
 *
 *  - `ingest`      — building the mesh through the uniform ingestion API
 *                    (points + connectivity + one point/cell data array),
 *                    i.e. what every format reader pays.
 *  - `traverse`    — a full writer-side accessor sweep (points, per-cell
 *                    connectivity, data arrays; checksummed so it can't be
 *                    optimized away), i.e. what every format writer pays.
 *  - `to_modelpart`— KRATOS only: `GetModelPart()` materialization on a
 *                    freshly ingested mesh (Nodes/Elements/Conditions +
 *                    variables + tag SubModelParts).
 *  - `write`/`read` per format — full file round-trips for gmsh (4.1
 *                    binary), vtu (binary+zlib when available), vtk
 *                    (binary), medit (ASCII) and su2 (ASCII).
 *
 * Usage: `meshioplusplus_bench [n]` — n is the grid size per edge
 * (default 35 -> 6*35^3 = 257k tets, 46k points).
 */

// System includes
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

constexpr int kRuns = 5;

struct RawGrid {
    NDArray mPoints;     // Float64 (npts, 3)
    NDArray mConn;       // Int64 (ntets, 4)
    NDArray mPointData;  // Float64 (npts,)
    NDArray mCellTags;   // Int64 (ntets,)
    std::size_t mNumCells = 0;
};

/** @brief Structured tet cube: (n+1)^3 shared vertices, 6 tets per hex cell. */
RawGrid make_tet_grid(std::size_t n) {
    const std::size_t np = n + 1;
    const std::size_t npts = np * np * np;
    RawGrid g;
    g.mPoints = NDArray::Uninit(DType::Float64, {npts, 3});
    double* p = g.mPoints.As<double>();
    for (std::size_t k = 0; k < np; ++k)
        for (std::size_t j = 0; j < np; ++j)
            for (std::size_t i = 0; i < np; ++i) {
                const std::size_t idx = (k * np + j) * np + i;
                p[idx * 3 + 0] = static_cast<double>(i) / static_cast<double>(n);
                p[idx * 3 + 1] = static_cast<double>(j) / static_cast<double>(n);
                p[idx * 3 + 2] = static_cast<double>(k) / static_cast<double>(n);
            }
    const std::size_t ntets = 6 * n * n * n;
    g.mNumCells = ntets;
    g.mConn = NDArray::Uninit(DType::Int64, {ntets, 4});
    std::int64_t* c = g.mConn.As<std::int64_t>();
    auto vid = [np](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * np + j) * np + i);
    };
    // Kuhn 6-tet decomposition of each cube.
    static const int tets[6][4][3] = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 1, 1}},
    };
    std::size_t t = 0;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                for (int s = 0; s < 6; ++s, ++t)
                    for (int v = 0; v < 4; ++v)
                        c[t * 4 + v] = vid(i + tets[s][v][0], j + tets[s][v][1], k + tets[s][v][2]);
    g.mPointData = NDArray::Uninit(DType::Float64, {npts});
    double* pd = g.mPointData.As<double>();
    for (std::size_t i = 0; i < npts; ++i)
        pd[i] = p[i * 3] + 2.0 * p[i * 3 + 1];
    g.mCellTags = NDArray::Uninit(DType::Int64, {ntets});
    std::int64_t* tag = g.mCellTags.As<std::int64_t>();
    for (std::size_t i = 0; i < ntets; ++i)
        tag[i] = static_cast<std::int64_t>(i % 4 + 1);
    return g;
}

NDArray copy_array(const NDArray& rA) {
    NDArray out = NDArray::Uninit(rA.Dtype(), rA.Shape());
    std::copy(rA.Data(), rA.Data() + rA.Nbytes(), out.Data());
    return out;
}

/** @brief Ingest the raw buffers through the uniform API (arrays copied first
 *  so every run pays the same allocation, then moved in like a reader). */
Mesh ingest(const RawGrid& rGrid) {
    Mesh m;
    m.AssignPoints(copy_array(rGrid.mPoints));
    m.AddCellBlock("tetra", copy_array(rGrid.mConn));
    m.AddPointData("field", copy_array(rGrid.mPointData));
    m.AddCellData("gmsh:physical", {copy_array(rGrid.mCellTags)});
    return m;
}

/** @brief Writer-side sweep: touch all points, connectivity and data. */
double traverse(const Mesh& rMesh) {
    double sum = 0.0;
    const NDArray& points = rMesh.Points();
    const double* p = points.As<double>();
    const std::size_t np = rMesh.NumPoints() * rMesh.PointDim();
    for (std::size_t i = 0; i < np; ++i)
        sum += p[i];
    std::int64_t csum = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& conn = cb.Conn();
        const std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < conn.Size(); ++i)
            csum += c[i];
    }
    for (const auto& r_name : rMesh.PointDataNames())
        sum += rMesh.PointData(r_name).As<double>()[0];
    return sum + static_cast<double>(csum % 1000);
}

double median_seconds(const std::function<void()>& rF) {
    rF();  // warmup
    std::vector<double> times;
    times.reserve(kRuns);
    for (int r = 0; r < kRuns; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        rF();
        const auto t1 = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration<double>(t1 - t0).count());
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
}

void row(const char* pOp, const char* pFormat, std::size_t cells, double median) {
    std::printf("%s,%s,%s,%zu,%.6f,%d\n", meshioplusplus::mesh_backend_name(), pOp, pFormat, cells,
                median, kRuns);
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? static_cast<std::size_t>(std::stoul(argv[1])) : 35;
    const RawGrid grid = make_tet_grid(n);
    std::fprintf(stderr, "backend=%s n=%zu cells=%zu points=%zu\n",
                 meshioplusplus::mesh_backend_name(), n, grid.mNumCells, grid.mPoints.Shape()[0]);
    std::printf("backend,op,format,cells,median_s,runs\n");

    row("ingest", "-", grid.mNumCells, median_seconds([&] {
            Mesh m = ingest(grid);
            (void)m;
        }));

    const Mesh mesh = ingest(grid);
    volatile double sink = 0.0;
    row("traverse", "-", grid.mNumCells, median_seconds([&] { sink = traverse(mesh); }));
    (void)sink;

#ifdef MESHIOPLUSPLUS_MESH_BACKEND_KRATOS
    row("to_modelpart", "-", grid.mNumCells, median_seconds([&] {
            Mesh m = ingest(grid);
            (void)m.GetModelPart();  // Nodes/Elements/variables/SubModelParts
        }));
#endif

    struct Fmt {
        const char* mName;
        const char* mExt;
        std::function<void(const std::string&, const Mesh&)> mWrite;
        std::function<Mesh(const std::string&)> mRead;
    };
    const std::vector<Fmt> formats = {
        {"gmsh41-bin", ".msh",
         [](const std::string& rP, const Mesh& rM) { meshioplusplus::write_gmsh41(rP, rM, true); },
         [](const std::string& rP) { return meshioplusplus::read_gmsh(rP); }},
        {"vtu-bin", ".vtu",
         [](const std::string& rP, const Mesh& rM) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
             meshioplusplus::write_vtu(rP, rM, true, true);
#else
             meshioplusplus::write_vtu(rP, rM, true, false);
#endif
         },
         [](const std::string& rP) { return meshioplusplus::read_vtu(rP); }},
        {"vtk-bin", ".vtk",
         [](const std::string& rP, const Mesh& rM) {
             meshioplusplus::write_vtk(rP, rM, true, false);
         },
         [](const std::string& rP) { return meshioplusplus::read_vtk(rP); }},
        {"medit-ascii", ".mesh",
         [](const std::string& rP, const Mesh& rM) { meshioplusplus::write_medit_ascii(rP, rM); },
         [](const std::string& rP) { return meshioplusplus::read_medit_ascii(rP); }},
        {"su2", ".su2",
         [](const std::string& rP, const Mesh& rM) { meshioplusplus::write_su2(rP, rM); },
         [](const std::string& rP) { return meshioplusplus::read_su2(rP); }},
    };

    const auto tmp = std::filesystem::temp_directory_path();
    for (const auto& r_fmt : formats) {
        const std::string path = (tmp / (std::string("meshio_bench") + r_fmt.mExt)).string();
        row("write", r_fmt.mName, grid.mNumCells,
            median_seconds([&] { r_fmt.mWrite(path, mesh); }));
        row("read", r_fmt.mName, grid.mNumCells, median_seconds([&] {
                Mesh m = r_fmt.mRead(path);
                (void)m;
            }));
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
    return 0;
}
