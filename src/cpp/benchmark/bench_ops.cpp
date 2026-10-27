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
 * @file bench_ops.cpp
 * @brief Timings of the mesh operations over a size sweep (roadmap §2/§3).
 *
 * The companion of `bench_backends.cpp` for operations rather than I/O: the
 * operations below on the same structured tetrahedral cube (or its surface),
 * at several sizes. The parallel backend is a compile-time choice, so a binary
 * measures one backend; `tools/bench_ops.sh` builds SEQ, OpenMP and TBB trees
 * and sweeps `OMP_NUM_THREADS`. Warmup plus the median of N runs,
 * `std::chrono::steady_clock`, no framework.
 *
 * Output, one CSV row per (tier, op) on stdout:
 * `backend,threads,op,cells,median_s,runs,digest`.
 *
 * `digest` is empty unless `--hash` is given; then it is a 64-bit FNV-1a
 * digest of the warmup run's result -- points, every block's connectivity,
 * point/cell/field data and regions, in block and sorted-name order, plus the
 * result's index maps and counters -- so `tools/bench_ops.sh --check` can
 * prove a row's output is byte-identical across parallel backends and thread
 * counts (roadmap §3, "The determinism check the section assumes"). The warmup
 * run is the one hashed, so hashing never enters a timing.
 *
 * Usage: `meshioplusplus_bench_ops [--tier S|M|L|XL]... [--ops a,b,...] [--runs N] [--hash]`.
 * The default tiers are S, M and L (about 10k, 160k and 750k tetrahedra);
 * XL (about 10M) is opt-in, since it needs several GB.
 */

// System includes
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/project_to_cad.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/sobolev_deform.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/render.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/operations/voxelize.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/skin.hpp"
#include "mesh_digest.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::bench::MeshDigest;

namespace {

/** @brief Structured tet cube: (n+1)^3 shared vertices, 6 tets per hex (Kuhn). */
Mesh bench_ops_tet_cube(std::size_t n) {
    const std::size_t np = n + 1;
    NDArray pts = NDArray::Uninit(DType::Float64, {np * np * np, 3});
    double* p = pts.As<double>();
    for (std::size_t k = 0; k < np; ++k)
        for (std::size_t j = 0; j < np; ++j)
            for (std::size_t i = 0; i < np; ++i) {
                const std::size_t idx = (k * np + j) * np + i;
                p[idx * 3 + 0] = static_cast<double>(i) / static_cast<double>(n);
                p[idx * 3 + 1] = static_cast<double>(j) / static_cast<double>(n);
                p[idx * 3 + 2] = static_cast<double>(k) / static_cast<double>(n);
            }
    static const int tets[6][4][3] = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 1, 1}},
    };
    // Half the Kuhn tets are negatively oriented as listed; swapping their
    // last two corners makes every tet positive, so the extracted surface is
    // consistently wound (decimate and compute_sdf need that to do real work).
    int order[6][4];
    for (int s = 0; s < 6; ++s) {
        const int* a = tets[s][0];
        int e[3][3];
        for (int v = 1; v < 4; ++v)
            for (int d = 0; d < 3; ++d)
                e[v - 1][d] = tets[s][v][d] - a[d];
        const int det = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                        e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0]) +
                        e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
        const bool flip = det < 0;
        order[s][0] = 0;
        order[s][1] = 1;
        order[s][2] = flip ? 3 : 2;
        order[s][3] = flip ? 2 : 3;
    }
    NDArray conn = NDArray::Uninit(DType::Int64, {6 * n * n * n, 4});
    std::int64_t* c = conn.As<std::int64_t>();
    std::size_t t = 0;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                for (int s = 0; s < 6; ++s, ++t)
                    for (int v = 0; v < 4; ++v)
                        c[t * 4 + v] = static_cast<std::int64_t>(
                            ((k + tets[s][order[s][v]][2]) * np + (j + tets[s][order[s][v]][1])) *
                                np +
                            (i + tets[s][order[s][v]][0]));
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

/**
 * @brief A planar triangle lattice of `m x m` squares (two triangles each) at
 * unit size, with two-component points: the 2-D input the Triangle reader
 * needs, which the 3-D cube cannot supply.
 */
Mesh bench_ops_tri_plate(std::size_t m) {
    const std::size_t np = m + 1;
    NDArray pts = NDArray::Uninit(DType::Float64, {np * np, 2});
    double* p = pts.As<double>();
    for (std::size_t j = 0; j < np; ++j)
        for (std::size_t i = 0; i < np; ++i) {
            p[(j * np + i) * 2 + 0] = static_cast<double>(i) / static_cast<double>(m);
            p[(j * np + i) * 2 + 1] = static_cast<double>(j) / static_cast<double>(m);
        }
    NDArray conn = NDArray::Uninit(DType::Int64, {2 * m * m, 3});
    std::int64_t* c = conn.As<std::int64_t>();
    std::size_t t = 0;
    for (std::size_t j = 0; j < m; ++j)
        for (std::size_t i = 0; i < m; ++i) {
            const std::int64_t a = static_cast<std::int64_t>(j * np + i);
            const std::int64_t b = a + 1;
            const std::int64_t d = a + static_cast<std::int64_t>(np);
            const std::int64_t e = d + 1;
            c[t * 3 + 0] = a;
            c[t * 3 + 1] = b;
            c[t * 3 + 2] = e;
            ++t;
            c[t * 3 + 0] = a;
            c[t * 3 + 1] = e;
            c[t * 3 + 2] = d;
            ++t;
        }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    mesh.AddCellBlock("triangle", std::move(conn));
    return mesh;
}

/**
 * @brief A copy of `rMesh` whose points are moved by `f(point index, xyz)`.
 * Used for the jittered cube (`optimize_volume` has nothing to flip on the
 * regular one) and the inflated surface `shrinkwrap` projects back.
 */
Mesh bench_ops_moved(const Mesh& rMesh, const std::function<void(std::size_t, double*)>& rF) {
    NDArray pts = NDArray::Uninit(DType::Float64, {rMesh.NumPoints(), 3});
    std::memcpy(pts.Data(), rMesh.Points().Data(), pts.Nbytes());
    double* p = pts.As<double>();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        rF(i, p + 3 * i);
    Mesh m;
    m.AssignPoints(std::move(pts));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        NDArray conn = NDArray::Uninit(cb.Conn().Dtype(), cb.Conn().Shape());
        std::memcpy(conn.Data(), cb.Conn().Data(), conn.Nbytes());
        m.AddCellBlock(cb.Type(), std::move(conn));
    }
    return m;
}

/**
 * @brief `rVolume`'s points with its tetrahedra recast as a polyhedron block
 * (four triangular faces each) and one triangle per tetrahedron as a polygon
 * block: the ragged input the `reorder` CSR rows run on.
 */
Mesh bench_ops_ragged(const Mesh& rVolume) {
    const auto cb = rVolume.Cells(0);
    const std::size_t nc = cb.NumCells();
    static const int faces[4][3] = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
    std::vector<std::int64_t> flat, rows{0}, cells{0}, ring, poly_rows{0};
    flat.reserve(nc * 12);
    rows.reserve(nc * 4 + 1);
    cells.reserve(nc + 1);
    ring.reserve(nc * 3);
    poly_rows.reserve(nc + 1);
    for (std::size_t c = 0; c < nc; ++c) {
        const std::int64_t* v = cb.Conn().As<std::int64_t>() + 4 * c;
        for (const auto& f : faces) {
            for (int k = 0; k < 3; ++k)
                flat.push_back(v[f[k]]);
            rows.push_back(static_cast<std::int64_t>(flat.size()));
        }
        cells.push_back(static_cast<std::int64_t>(rows.size() - 1));
        for (int k = 0; k < 3; ++k)
            ring.push_back(v[faces[0][k]]);
        poly_rows.push_back(static_cast<std::int64_t>(ring.size()));
    }
    NDArray pts = NDArray::Uninit(DType::Float64, {rVolume.NumPoints(), 3});
    std::memcpy(pts.Data(), rVolume.Points().Data(), pts.Nbytes());
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddPolyhedronBlock("polyhedron4", std::move(flat), std::move(rows), std::move(cells));
    m.AddPolygonBlock("polygon", std::move(ring), std::move(poly_rows));
    return m;
}

/**
 * @brief `rMesh` with non-canonical storage: float32 points and point data,
 * int32 connectivity. The `_narrow` rows run on it, because the dtype-hoisted
 * loops are zero-copy on float64/int64 and only a narrower input exercises the
 * one converted copy (roadmap §3.3.2).
 */
Mesh bench_ops_narrowed(const Mesh& rMesh) {
    const auto narrow = [](const NDArray& rA, DType To) {
        NDArray out = NDArray::Uninit(To, rA.Shape());
        const std::size_t n = rA.Size();
        if (To == DType::Float32) {
            float* d = out.As<float>();
            for (std::size_t i = 0; i < n; ++i)
                d[i] = static_cast<float>(meshioplusplus::detail::read_double(rA, i));
        } else {
            std::int32_t* d = out.As<std::int32_t>();
            for (std::size_t i = 0; i < n; ++i)
                d[i] = static_cast<std::int32_t>(meshioplusplus::detail::read_int(rA, i));
        }
        return out;
    };
    Mesh m;
    m.AssignPoints(narrow(rMesh.Points(), DType::Float32));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        m.AddCellBlock(cb.Type(), narrow(cb.Conn(), DType::Int32));
    }
    for (const std::string& name : rMesh.PointDataNames())
        m.AddPointData(name, narrow(rMesh.PointData(name), DType::Float32));
    return m;
}

double bench_ops_median(const std::function<void()>& rFn, int Runs) {
    std::vector<double> ts;
    for (int r = 0; r < Runs; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        rFn();
        ts.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    std::sort(ts.begin(), ts.end());
    return ts[ts.size() / 2];
}

// A derived input built on first use. At tier XL (10M tetrahedra) every
// derived mesh held at once is more than a 15 GiB machine has, so a run that
// selects rows with `--ops` builds only the inputs those rows read. What a row
// times and hashes is unchanged: the mesh is built before its warmup run.
class BenchLazyMesh {
public:
    explicit BenchLazyMesh(std::function<Mesh()> Make) : mMake(std::move(Make)) {}
    const Mesh& Get() const {
        if (!mMesh)
            mMesh.emplace(mMake());
        return *mMesh;
    }
    operator const Mesh&() const { return Get(); }
    const Mesh* operator->() const { return &Get(); }

private:
    std::function<Mesh()> mMake;
    mutable std::optional<Mesh> mMesh;
};

}  // namespace

int main(int argc, char** argv) {
    namespace mio = meshioplusplus;
    const std::map<std::string, std::size_t> tiers = {{"S", 12}, {"M", 30}, {"L", 50}, {"XL", 120}};
    std::vector<std::string> chosen;
    std::vector<std::string> ops;
    int runs = 3;
    bool hash = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--tier" && i + 1 < argc) {
            chosen.push_back(argv[++i]);
        } else if (a == "--ops" && i + 1 < argc) {
            const std::string list = argv[++i];
            for (std::size_t at = 0; at <= list.size();) {
                const std::size_t comma = std::min(list.find(',', at), list.size());
                if (comma > at)
                    ops.push_back(list.substr(at, comma - at));
                at = comma + 1;
            }
        } else if (a == "--runs" && i + 1 < argc) {
            runs = std::max(1, std::atoi(argv[++i]));
        } else if (a == "--hash") {
            hash = true;
        } else {
            std::fprintf(stderr,
                         "usage: %s [--tier S|M|L|XL]... [--ops a,b,...] [--runs N] [--hash]\n",
                         argv[0]);
            return 2;
        }
    }
    if (chosen.empty())
        chosen = {"S", "M", "L"};
    // OMP_NUM_THREADS is honoured by OpenMP itself; TBB has no such variable,
    // so the same one caps it here, keeping tools/bench_ops.sh's sweep uniform.
    const char* env_threads = std::getenv("OMP_NUM_THREADS");
    const std::string threads =
        env_threads ? env_threads : std::to_string(std::thread::hardware_concurrency());
#if defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    std::optional<tbb::global_control> cap;
    if (env_threads && std::atoi(env_threads) > 0)
        cap.emplace(tbb::global_control::max_allowed_parallelism,
                    static_cast<std::size_t>(std::atoi(env_threads)));
#endif
    const auto wanted = [&](const char* pOp) {
        return ops.empty() || std::find(ops.begin(), ops.end(), pOp) != ops.end();
    };

    std::printf("backend,threads,op,cells,median_s,runs,digest\n");
    for (const std::string& tier : chosen) {
        auto it = tiers.find(tier);
        if (it == tiers.end()) {
            std::fprintf(stderr, "unknown tier '%s'\n", tier.c_str());
            return 2;
        }
        const std::size_t n = it->second;
        const Mesh volume = bench_ops_tet_cube(n);
        const Mesh surface = mio::extract_surface(volume);
        const std::size_t ncells = volume.Cells(0).NumCells();
        // A deterministic jitter of up to 0.3 of a lattice step on every
        // point (boundary points slide in their face plane only, which keeps
        // the cube's boundary intact for `optimize_volume`).
        const double h = 1.0 / static_cast<double>(n);
        const BenchLazyMesh jittered([&] {
            return bench_ops_moved(volume, [&](std::size_t i, double* pP) {
                for (int d = 0; d < 3; ++d) {
                    if (pP[d] <= 0.0 || pP[d] >= 1.0)
                        continue;
                    pP[d] += 0.3 * h * std::sin(1.7 * static_cast<double>(3 * i + d) + 0.3);
                }
            });
        });
        const BenchLazyMesh inflated([&] {
            return bench_ops_moved(surface, [](std::size_t, double* pP) {
                for (int d = 0; d < 3; ++d)
                    pP[d] = 0.5 + 1.25 * (pP[d] - 0.5);
            });
        });
        // Two copies of the cube, unwelded: the welding rows' input.
        const BenchLazyMesh doubled([&] {
            mio::MergeOptions o;
            o.source_tag = false;
            return mio::merge({&volume, &volume}, o).mMesh;
        });
        // The volume with a smooth point field, for `hessian`.
        const BenchLazyMesh with_field([&] {
            Mesh m = bench_ops_moved(volume, [](std::size_t, double*) {});
            NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
            const double* p = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i)
                u.As<double>()[i] = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 2];
            m.AddPointData("u", std::move(u));
            return m;
        });
        // The cube elevated to quadratic tetrahedra, for `linearize`.
        const BenchLazyMesh quadratic([&] {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Elevate;
            return mio::convert_cells(volume, o).mMesh;
        });
        // A green-closed refinement of every 7th cell, for `undo_green`.
        const BenchLazyMesh green([&] {
            mio::RefineOptions o;
            for (std::size_t c = 0; c < ncells; c += 7)
                o.mCells.push_back(static_cast<std::int64_t>(c));
            o.mRecordHierarchy = true;
            o.mRecordLevels = true;
            return mio::refine(volume, o).mMesh;
        });
        // The same meshes with float32 points/data and int32 connectivity.
        const BenchLazyMesh volume_narrow([&] { return bench_ops_narrowed(volume); });
        const BenchLazyMesh with_field_narrow([&] { return bench_ops_narrowed(with_field); });
        const BenchLazyMesh jittered_narrow([&] { return bench_ops_narrowed(jittered); });
        const BenchLazyMesh quadratic_narrow([&] { return bench_ops_narrowed(quadratic); });
        using Hashed = std::function<void(MeshDigest*)>;
        const auto row = [&](const char* pOp, const Hashed& rFn) {
            if (!wanted(pOp))
                return;
            MeshDigest digest;
            rFn(hash ? &digest : nullptr);  // warmup, hashed with --hash
            const double s = bench_ops_median([&] { rFn(nullptr); }, runs);
            char dig[24] = "";
            if (hash)
                std::snprintf(dig, sizeof dig, "%016llx",
                              static_cast<unsigned long long>(digest.Value()));
            std::printf("%s,%s,%s,%zu,%.6f,%d,%s\n", mio::parallel_backend_name(), threads.c_str(),
                        pOp, ncells, s, runs, dig);
            std::fflush(stdout);
        };
        // Wraps a call returning a mesh (or a result carrying one) so the
        // digest covers what the row produces.
        const auto of = [](MeshDigest* pD, const Mesh& rM) {
            if (pD)
                pD->Of(rM);
        };
        row("extract_surface", [&](MeshDigest* pD) { of(pD, mio::extract_surface(volume)); });
        row("extract_skin", [&](MeshDigest* pD) { of(pD, mio::extract_skin(volume)); });
        row("smooth", [&](MeshDigest* pD) {
            mio::SmoothOptions o;
            o.mIterations = 10;
            auto r = mio::smooth(surface, o);
            of(pD, r.mMesh);
            if (pD)
                pD->U64(static_cast<std::uint64_t>(r.mNumNodesMoved));
        });
        row("smooth_volume", [&](MeshDigest* pD) {
            mio::SmoothOptions o;
            o.mIterations = 10;
            auto r = mio::smooth(jittered, o);
            of(pD, r.mMesh);
            if (pD)
                pD->U64(static_cast<std::uint64_t>(r.mNumNodesMoved));
        });
        // A topology-preserving step feeding a facet reader: both build their own
        // facet table over the same cells. The row times the chain, so the cost of
        // building that table twice is visible next to the `smooth_volume` and
        // `extract_surface` rows it is made of (doc/benchmarks.md, "Pipeline row:
        // what a shared facet table could save").
        row("pipeline_smooth_surface", [&](MeshDigest* pD) {
            mio::PipelineStep smooth_step;
            smooth_step.mOp = "Smooth";
            smooth_step.mParams["Iterations"] = std::int64_t{10};
            mio::PipelineStep surface_step;
            surface_step.mOp = "ExtractSurface";
            mio::PipelineReport report;
            of(pD, mio::run_pipeline_steps(jittered, {smooth_step, surface_step}, report));
        });
        row("refine", [&](MeshDigest* pD) {
            auto r = mio::refine(volume);
            of(pD, r.mMesh);
            if (pD)
                pD->Array(r.mPointMap);
        });
        row("elevate", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Elevate;
            auto r = mio::convert_cells(volume, o);
            of(pD, r.mMesh);
            if (pD) {
                pD->Array(r.mPointMap);
                pD->Arrays(r.mCellMaps);
            }
        });
        row("merge", [&](MeshDigest* pD) {
            mio::MergeOptions o;
            o.weld = true;
            auto r = mio::merge({&volume, &volume}, o);
            of(pD, r.mMesh);
            if (pD) {
                pD->Arrays(r.mPointMaps);
                pD->Arrays(r.mCellMaps);
            }
        });
        const auto clean_row = [&](const char* pOp, const Mesh& rIn, bool weld) {
            row(pOp, [&](MeshDigest* pD) {
                mio::CleanOptions o;
                o.weld = weld;
                auto r = mio::clean(rIn, o);
                of(pD, r.mMesh);
                if (pD) {
                    pD->Array(r.mPointMap);
                    pD->Arrays(r.mCellMaps);
                    pD->U64(static_cast<std::uint64_t>(r.mPointsWelded));
                    pD->U64(static_cast<std::uint64_t>(r.mCellsDroppedDuplicate));
                }
            });
        };
        clean_row("clean", volume, false);
        clean_row("clean_weld", doubled, true);
        row("compute_sdf", [&](MeshDigest* pD) {
            mio::SdfOptions o;
            o.mResolution = std::array<std::int64_t, 3>{48, 48, 48};
            of(pD, mio::compute_sdf(surface, o).mMesh);
        });
        row("decimate", [&](MeshDigest* pD) {
            mio::DecimateOptions o;
            o.mTargetRatio = 0.5;
            auto r = mio::decimate(surface, o);
            of(pD, r.mMesh);
            if (pD)
                pD->Array(r.mPointMap);
        });
        row("partition", [&](MeshDigest* pD) {
            mio::PartitionOptions o;
            o.mNParts = 8;
            auto r = mio::partition(volume, o);
            if (pD)
                for (const auto& piece : r.mPieces) {
                    pD->U64(static_cast<std::uint64_t>(piece.mPartId));
                    pD->Of(piece.mMesh);
                    pD->Array(piece.mPointMap);
                    pD->Arrays(piece.mCellMaps);
                }
        });
        const auto reorder_row = [&](const char* pOp, mio::ReorderMethod method) {
            row(pOp, [&, method](MeshDigest* pD) {
                auto r = mio::reorder(volume, method);
                of(pD, r.mMesh);
                if (pD) {
                    pD->Array(r.mNodePermutation);
                    pD->Arrays(r.mCellPermutations);
                }
            });
        };
        reorder_row("reorder", mio::ReorderMethod::RCM);
        reorder_row("reorder_hilbert", mio::ReorderMethod::Hilbert);
        const BenchLazyMesh ragged([&] { return bench_ops_ragged(volume); });
        row("reorder_ragged", [&](MeshDigest* pD) {
            auto r = mio::reorder(ragged, mio::ReorderMethod::RCM);
            of(pD, r.mMesh);
            if (pD) {
                pD->Array(r.mNodePermutation);
                pD->Arrays(r.mCellPermutations);
            }
        });
        // Simplexify of the same ragged mesh: every polyhedron fans into
        // tetrahedra and every polygon into triangles, which is where the
        // fan paths' output vectors grow (roadmap §3.2.2).
        row("simplexify_ragged", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Simplexify;
            auto r = mio::convert_cells(ragged, o);
            of(pD, r.mMesh);
            if (pD) {
                pD->Array(r.mPointMap);
                pD->Arrays(r.mCellMaps);
            }
        });
        // Linearize of the same ragged mesh changes nothing in either block, so
        // both pass through convert_cells' staging: the copy of every polyhedron
        // and polygon into the block that is emitted.
        row("linearize_ragged", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Linearize;
            auto r = mio::convert_cells(ragged, o);
            of(pD, r.mMesh);
            if (pD) {
                pD->Array(r.mPointMap);
                pD->Arrays(r.mCellMaps);
            }
        });
        // Subdivide of the same ragged mesh: one apex point per polyhedron,
        // whose source-node list is the allocation of roadmap §3.2.2.
        row("subdivide_ragged", [&](MeshDigest* pD) {
            auto r = mio::subdivide(ragged);
            of(pD, r.mMesh);
            if (pD)
                pD->Arrays(r.mCellMaps);
        });
        // The ragged mesh with a smooth point field and one cell field per
        // block: the per-cell and per-polyhedron scratch rows of roadmap
        // §3.2.3 (quality, data_average, interpolate's cell_data and the SFC
        // partition's centroids) run on it.
        const BenchLazyMesh ragged_field([&] {
            Mesh m = bench_ops_ragged(volume);
            NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
            const double* p = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i)
                u.As<double>()[i] = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 2];
            m.AddPointData("u", std::move(u));
            std::vector<NDArray> blocks;
            for (const auto cb : m.CellRange()) {
                NDArray c = NDArray::Uninit(DType::Float64, {cb.NumCells()});
                for (std::size_t i = 0; i < cb.NumCells(); ++i)
                    c.As<double>()[i] = std::sin(0.01 * static_cast<double>(i));
                blocks.push_back(std::move(c));
            }
            m.AddCellData("c", std::move(blocks));
            return m;
        });
        row("quality", [&](MeshDigest* pD) { of(pD, mio::attach_quality(volume)); });
        row("quality_ragged", [&](MeshDigest* pD) { of(pD, mio::attach_quality(ragged)); });
        row("data_average",
            [&](MeshDigest* pD) { of(pD, mio::point_data_to_cell_data(with_field)); });
        row("data_average_ragged",
            [&](MeshDigest* pD) { of(pD, mio::point_data_to_cell_data(ragged_field)); });
        row("interpolate_cells", [&](MeshDigest* pD) {
            mio::InterpolateOptions o;
            o.mArrays = {"c"};
            o.mOnConflict = mio::InterpolateConflict::Overwrite;
            of(pD, mio::interpolate(ragged_field, ragged_field, o));
        });
        row("partition_ragged", [&](MeshDigest* pD) {
            mio::PartitionOptions o;
            o.mNParts = 8;
            o.mMethod = mio::PartitionMethod::SFC;
            auto r = mio::partition(ragged, o);
            if (pD)
                for (const auto& piece : r.mPieces) {
                    pD->U64(static_cast<std::uint64_t>(piece.mPartId));
                    pD->Of(piece.mMesh);
                    pD->Array(piece.mPointMap);
                    pD->Arrays(piece.mCellMaps);
                }
        });
        row("feature_edges_ragged",
            [&](MeshDigest* pD) { of(pD, mio::feature_edges(ragged).mMesh); });
        // Reads through the registry: the file is written once, outside the
        // timed region. A format whose writer declines the mesh is skipped, so
        // the rows never fail the sweep.
        const auto read_row = [&](const char* pOp, const char* pFormat, const char* pExt,
                                  const Mesh& rSource) {
            if (!wanted(pOp))
                return;
            const std::filesystem::path dir = std::filesystem::temp_directory_path() /
                                              (std::string("meshioplusplus_bench_read_") + pFormat);
            std::filesystem::create_directories(dir);
            const std::string path = (dir / (std::string("mesh") + pExt)).string();
            try {
                mio::registry_writers().at(pFormat)(path, rSource);
                row(pOp, [&](MeshDigest* pD) {
                    of(pD, mio::registry_read(path, pFormat, mio::ReadOptions{}));
                });
            } catch (const std::exception& rExc) {
                std::fprintf(stderr, "skip %s: %s\n", pOp, rExc.what());
            }
            std::filesystem::remove_all(dir);
        };
        read_row("read_ragged_ensight", "ensight", ".case", ragged);
        read_row("read_ragged_tecplot", "tecplot", ".dat", ragged);
        // The registry's STL writer is ASCII, so this times the ASCII reader
        // over a triangle surface (roadmap §3.1.1.1).
        read_row("read_stl", "stl", ".stl", surface);
        // OBJ over the cube's triangle surface with a normal and a texture
        // coordinate per point, so the `vn` and `vt` lines are parsed too; the
        // write row times the same mesh (roadmap §3.1.1.1 keeps an OBJ writer
        // observation open, and the two rows isolate the read from it).
        const Mesh obj_surface = [&] {
            Mesh m = bench_ops_moved(surface, [](std::size_t, double*) {});
            const std::size_t np = m.NumPoints();
            const double* p = m.Points().As<double>();
            NDArray vn = NDArray::Uninit(DType::Float64, {np, std::size_t(3)});
            NDArray vt = NDArray::Uninit(DType::Float64, {np, std::size_t(2)});
            for (std::size_t i = 0; i < np; ++i) {
                for (std::size_t d = 0; d < 3; ++d)
                    vn.As<double>()[3 * i + d] = 0.5 - p[3 * i + d];
                vt.As<double>()[2 * i] = p[3 * i];
                vt.As<double>()[2 * i + 1] = p[3 * i + 1];
            }
            m.AddPointData("obj:vn", std::move(vn));
            m.AddPointData("obj:vt", std::move(vt));
            return m;
        }();
        read_row("read_obj", "obj", ".obj", obj_surface);
        if (wanted("write_obj")) {
            const std::filesystem::path dir =
                std::filesystem::temp_directory_path() / "meshioplusplus_bench_write_obj";
            std::filesystem::create_directories(dir);
            const std::string path = (dir / "mesh.obj").string();
            row("write_obj", [&](MeshDigest* pD) {
                mio::registry_writers().at("obj")(path, obj_surface);
                if (pD)
                    of(pD, mio::registry_read(path, "obj", mio::ReadOptions{}));
            });
            std::filesystem::remove_all(dir);
        }
        // The text-token readers of roadmap §3.1.1.1: TetGen's `.node`/`.ele`
        // pair over the cube with its point field (so an attribute column is
        // parsed too), and Triangle's over a planar lattice of about as many
        // triangles as the cube has tetrahedra.
        read_row("read_tetgen", "tetgen", ".node", with_field);
        const Mesh tri_plate = bench_ops_tri_plate(
            static_cast<std::size_t>(std::sqrt(3.0 * static_cast<double>(n * n * n))));
        read_row("read_triangle", "triangle", ".node", tri_plate);
        // FreeFEM's `.msh` (a name only: the registry is given the format) and
        // UGRID's ASCII flavour (`.ugrid` carries no binary key), both over the
        // tetrahedral cube with its point field.
        read_row("read_freefem", "freefem", ".msh", with_field);
        read_row("read_ugrid", "ugrid", ".ugrid", with_field);
        // OpenFOAM's ASCII polyMesh (a `.foam` marker; the case directory is
        // written beside it) over the tetrahedral cube with its point field.
        read_row("read_openfoam", "openfoam", ".foam", with_field);
        // Elmer's text mesh (a directory of `mesh.nodes`, `mesh.elements` and
        // `mesh.boundary`) over the cube's tetrahedra plus the triangles of its
        // surface, so a boundary file is parsed too.
        const Mesh elmer_mesh = [&] {
            Mesh m = bench_ops_moved(volume, [](std::size_t, double*) {});
            m.AddCellBlock("triangle", surface.Cells(0).Conn());
            return m;
        }();
        read_row("read_elmer", "elmer", ".elmer", elmer_mesh);
        row("optimize_volume", [&](MeshDigest* pD) {
            auto r = mio::optimize_volume(jittered);
            of(pD, r.mMesh);
            if (pD)
                pD->U64(static_cast<std::uint64_t>(r.mNumFlips));
        });
        row("agglomerate", [&](MeshDigest* pD) {
            auto r = mio::agglomerate(volume);
            of(pD, r.mMesh);
            if (pD)
                pD->Array(r.mCellMap);
        });
        // The tetrahedral cube cut into two Cell regions (the lower and upper
        // halves of its cell list, a planar interface), for the interface rows.
        const Mesh with_halves = [&] {
            Mesh m = bench_ops_moved(volume, [](std::size_t, double*) {});
            const std::size_t half = ncells / 2;
            for (int side = 0; side < 2; ++side) {
                mio::Region region;
                region.mName = side == 0 ? "lower" : "upper";
                region.mKind = mio::RegionKind::Cell;
                region.mDim = 3;
                region.mTag = -1;
                const std::size_t first = side == 0 ? 0 : half;
                const std::size_t count = side == 0 ? half : ncells - half;
                region.mEntries = NDArray::Uninit(DType::Int64, {count});
                std::int64_t* e = region.mEntries.As<std::int64_t>();
                for (std::size_t i = 0; i < count; ++i)
                    e[i] = static_cast<std::int64_t>(first + i);
                m.AddRegion(std::move(region));
            }
            return m;
        }();
        row("region_adjacency",
            [&](MeshDigest* pD) { of(pD, mio::region_adjacency(with_halves)); });
        row("find_interface", [&](MeshDigest* pD) {
            auto r = mio::find_interface(with_halves, mio::RegionSelector{"lower"},
                                         mio::RegionSelector{"upper"});
            of(pD, r.mMesh);
            if (pD) {
                pD->U64(static_cast<std::uint64_t>(r.mReport.mNumPairs));
                pD->Bytes(&r.mReport.mArea, sizeof r.mReport.mArea);
            }
        });
        row("split", [&](MeshDigest* pD) {
            auto r = mio::split(doubled, mio::SplitBy::Component);
            if (pD)
                for (const auto& piece : r.mPieces) {
                    pD->Str(piece.mKey);
                    pD->Of(piece.mMesh);
                    pD->Array(piece.mPointMap);
                    pD->Arrays(piece.mCellMaps);
                }
        });
        row("undo_green", [&](MeshDigest* pD) {
            auto r = mio::undo_green(volume, green);
            of(pD, r.mMesh);
            if (pD)
                pD->Arrays(r.mCellMaps);
        });
        row("hessian", [&](MeshDigest* pD) {
            mio::HessianOptions o;
            o.mArrayName = "u";
            o.mLocation = mio::DataLocation::Point;
            of(pD, mio::hessian(with_field, o).mMesh);
        });
        row("compute_normals",
            [&](MeshDigest* pD) { of(pD, mio::compute_normals(surface).mMesh); });
        row("compute_curvature",
            [&](MeshDigest* pD) { of(pD, mio::compute_curvature(surface).mMesh); });
        // Rows for roadmap §5's analysis operations (v16.23.0).
        row("feature_edges", [&](MeshDigest* pD) { of(pD, mio::feature_edges(volume).mMesh); });
        // Rows for roadmap §6.3's skin to CAD: the cube's skin segments into its
        // six faces (the growth is the serial phase), and each face fits a plane.
        const BenchLazyMesh segmented([&] { return mio::segment(surface).mMesh; });
        row("segment", [&](MeshDigest* pD) { of(pD, mio::segment(surface).mMesh); });
        row("fit_primitives",
            [&](MeshDigest* pD) { of(pD, mio::fit_primitives(segmented).mMesh); });
        // Rows for fit_cad's freeform modes (v17.12.0): the same skin with a
        // monkey saddle lifted onto its top, fitted kernel-free, the saddle as
        // one adaptive B-spline face and as a patch network.
        const BenchLazyMesh saddle_fitted([&] {
            const double* q = surface.Points().As<double>();
            double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
            for (std::size_t i = 0; i < surface.NumPoints(); ++i)
                for (int j = 0; j < 3; ++j) {
                    lo[j] = std::min(lo[j], q[3 * i + j]);
                    hi[j] = std::max(hi[j], q[3 * i + j]);
                }
            const Mesh lifted = bench_ops_moved(surface, [&](std::size_t, double* pP) {
                if (pP[2] < hi[2])
                    return;
                const double u = 2.0 * (pP[0] - lo[0]) / (hi[0] - lo[0]) - 1.0;
                const double v = 2.0 * (pP[1] - lo[1]) / (hi[1] - lo[1]) - 1.0;
                pP[2] += 0.15 * (hi[2] - lo[2]) * (u * u * u - 3.0 * u * v * v);
            });
            return mio::fit_primitives(mio::segment(lifted).mMesh).mMesh;
        });
        for (const auto mode :
             {mio::FitCadFreeformMode::Adaptive, mio::FitCadFreeformMode::Network}) {
            row(mode == mio::FitCadFreeformMode::Adaptive ? "fit_cad_adaptive" : "fit_cad_network",
                [&, mode](MeshDigest* pD) {
                    mio::FitCadOptions o;
                    o.mBackend = mio::FitCadBackend::Native;
                    o.mFreeformMode = mode;
                    of(pD, mio::fit_cad(saddle_fitted, "", o).mMesh);
                });
        }
        // Rows for roadmap §7's software rasterizer (v16.33.0): a terminal-sized
        // frame of the volume (its skin extraction is part of the cost), and a
        // 4x-supersampled image of the field mesh with smooth shading, a mapped
        // field and every edge. The digest covers the pixels and the id buffer,
        // so the determinism check holds the raster to its contract.
        const auto frame_of = [](MeshDigest* pD, const mio::Frame& rF) {
            if (!pD)
                return;
            pD->Bytes(rF.mRgba.data(), rF.mRgba.size());
            pD->Bytes(rF.mCellIds.data(), rF.mCellIds.size() * sizeof(std::int64_t));
        };
        row("render_frame_160x96", [&](MeshDigest* pD) {
            mio::RenderOptions o;
            o.mWidth = 160;
            o.mHeight = 96;
            frame_of(pD, mio::render(volume, o));
        });
        row("render_ss4_800x480", [&](MeshDigest* pD) {
            mio::RenderOptions o;
            o.mWidth = 800;
            o.mHeight = 480;
            o.mSupersample = 4;
            o.mShading = mio::RenderShading::Smooth;
            o.mEdges = mio::RenderEdges::All;
            o.mColorBy = "u";
            frame_of(pD, mio::render(with_field, o));
        });
        // v16.37.0: streamlines through a cut-away volume. The digest covers the
        // pixels and the id buffer, so it holds the tracer (its point locator, its
        // seeds and its parallel per-seed buffers) to the same contract.
        const Mesh with_flow = [&] {
            Mesh m = bench_ops_moved(volume, [](std::size_t, double*) {});
            NDArray v = NDArray::Uninit(DType::Float64, {m.NumPoints(), std::size_t(3)});
            const double* p = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i) {
                v.As<double>()[3 * i] = -(p[3 * i + 1] - 0.5);
                v.As<double>()[3 * i + 1] = p[3 * i] - 0.5;
                v.As<double>()[3 * i + 2] = 0.1;
            }
            m.AddPointData("flow", std::move(v));
            return m;
        }();
        row("render_streamlines_320x192", [&](MeshDigest* pD) {
            mio::RenderOptions o;
            o.mWidth = 320;
            o.mHeight = 192;
            o.mSupersample = 2;
            o.mStreamlines = "flow";
            o.mStreamSeeds = 200;
            mio::RenderCutaway cut;
            cut.mPoint = {0.0, 0.0, 0.5};
            cut.mNormal = {0.0, 0.0, -1.0};
            o.mCutaways = {cut};
            frame_of(pD, mio::render(with_flow, o));
        });
        // v16.38.0: the synthwave theme (a sunset, the grid floor and the three
        // post-processes) on the cut-away volume.
        row("render_synthwave_320x192", [&](MeshDigest* pD) {
            mio::RenderOptions o;
            o.mWidth = 320;
            o.mHeight = 192;
            o.mSupersample = 2;
            o.mTheme = mio::RenderTheme::Synthwave;
            o.mGridFloor = true;
            o.mBloom = true;
            o.mFringe = true;
            o.mScanlines = true;
            o.mThemePhase = 3;
            frame_of(pD, mio::render(with_flow, o));
        });
        row("hausdorff", [&](MeshDigest* pD) {
            mio::HausdorffOptions o;
            o.mFaceSamples = 2;
            const mio::HausdorffResult r = mio::hausdorff_distance(inflated, surface, o);
            if (!pD)
                return;
            pD->Bytes(&r.mDistance, sizeof r.mDistance);
            pD->Bytes(&r.mMeanAtoB, sizeof r.mMeanAtoB);
            pD->Bytes(&r.mRmsBtoA, sizeof r.mRmsBtoA);
            pD->Bytes(r.mWorstPointA.data(), sizeof r.mWorstPointA);
        });
        row("shrinkwrap",
            [&](MeshDigest* pD) { of(pD, mio::shrinkwrap(inflated, surface).mMesh); });
        row("voxelize", [&](MeshDigest* pD) {
            mio::VoxelOptions o;
            o.mResolution = std::array<std::int64_t, 3>{48, 48, 48};
            of(pD, mio::voxelize(volume, o).mMesh);
        });
        // The parametric generators read no input mesh, so each row measures
        // construction alone, sized by the tier's lattice resolution. The budget
        // is lifted: the XL tier is far above the default limit by design.
        row("box", [&](MeshDigest* pD) {
            const auto m = static_cast<std::int64_t>(n);
            of(pD, mio::box({{m, m, m}}, {{1.0, 1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "hexahedron", 1, 0));
        });
        row("rectangle", [&](MeshDigest* pD) {
            const auto m = static_cast<std::int64_t>(4 * n);
            of(pD, mio::rectangle({{m, m}}, {{1.0, 1.0}}, {{0.0, 0.0, 0.0}}, "triangle", "zigzag",
                                  1, 0));
        });
        row("uv_sphere", [&](MeshDigest* pD) {
            const auto m = static_cast<std::int64_t>(4 * n);
            of(pD, mio::uv_sphere(1.0, 2 * m, m, {{0.0, 0.0, 0.0}}, "triangle", 1, 0));
        });
        row("icosa_sphere",
            [&](MeshDigest* pD) { of(pD, mio::icosa_sphere(1.0, 6, {{0.0, 0.0, 0.0}}, 1, 0)); });
        row("cylinder", [&](MeshDigest* pD) {
            const auto m = static_cast<std::int64_t>(n);
            of(pD, mio::cylinder(1.0, 1.0, m, m, {{0.0, 0.0, 0.0}}, "hexahedron", 1, 0));
        });
        // Rows for roadmap §3's remaining operation items (v16.19.0).
        row("remesh", [&](MeshDigest* pD) {
            mio::RemeshOptions o;
            o.mNumClusters = static_cast<std::int64_t>(surface.NumPoints() / 4);
            of(pD, mio::remesh(surface, o).mMesh);
        });
        row("remesh_volume", [&](MeshDigest* pD) {
            mio::RemeshVolumeOptions o;
            o.mCellSize = 1.5 / static_cast<double>(n);
            of(pD, mio::remesh_volume(surface, o).mMesh);
        });
        // The octree around the cube's surface (the baseline of the next two
        // rows), and `octree_mesh`'s castellated and tetrahedral outputs on a
        // size field from one to eight lattice steps.
        const double octree_step = 1.0 / static_cast<double>(n);
        row("octree", [&](MeshDigest* pD) {
            mio::OctreeOptions o;
            o.mMaxLevel =
                static_cast<std::int64_t>(std::ceil(std::log2(static_cast<double>(n)))) + 1;
            o.mFill = mio::VoxelFill::Inside;
            of(pD, mio::octree(surface, o).mMesh);
        });
        row("octree_mesh_castellated", [&](MeshDigest* pD) {
            mio::OctreeMeshOptions o;
            o.mOutput = mio::OctreeMeshOutput::Castellated;
            o.mMinSize = octree_step;
            o.mMaxSize = 8.0 * octree_step;
            of(pD, mio::octree_mesh(surface, o).mMesh);
        });
        row("octree_mesh_tet", [&](MeshDigest* pD) {
            mio::OctreeMeshOptions o;
            o.mMinSize = octree_step;
            o.mMaxSize = 8.0 * octree_step;
            of(pD, mio::octree_mesh(surface, o).mMesh);
        });
        row("octree_mesh_hex", [&](MeshDigest* pD) {
            mio::OctreeMeshOptions o;
            o.mOutput = mio::OctreeMeshOutput::Hexa;
            o.mSurfaceSize = 2.0 * octree_step;
            o.mMaxSize = 8.0 * octree_step;
            of(pD, mio::octree_mesh(surface, o).mMesh);
        });
        row("octree_mesh_hex_graded", [&](MeshDigest* pD) {
            mio::OctreeMeshOptions o;
            o.mOutput = mio::OctreeMeshOutput::Hexa;
            o.mHexGrading = mio::OctreeMeshHexGrading::Graded;
            o.mMinSize = 2.0 * octree_step;
            o.mMaxSize = 8.0 * octree_step;
            of(pD, mio::octree_mesh(surface, o).mMesh);
        });
        // The CAD models of the mesh_cad_surface rows, built inside each row's
        // timed body as before.
        const auto cad_box = [&]() {
            mio::BrepInfo box;
            const double sides[6][3][3] = {
                {{0, 0, 0}, {0, 1, 0}, {1, 0, 0}}, {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
                {{0, 0, 0}, {1, 0, 0}, {0, 0, 1}}, {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
                {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
            for (int f = 0; f < 6; ++f) {
                mio::NurbsPatch p;
                p.mDegrees = {1, 1};
                p.mKnots = {{0, 0, 1, 1}, {0, 0, 1, 1}};
                p.mControlPoints = NDArray(mio::DType::Float64, {4, 3});
                p.mWeights = NDArray(mio::DType::Float64, {4});
                for (int j = 0; j < 2; ++j)
                    for (int i = 0; i < 2; ++i)
                        for (int k = 0; k < 3; ++k)
                            p.mControlPoints.As<double>()[3 * (2 * j + i) + k] =
                                sides[f][0][k] + i * sides[f][1][k] + j * sides[f][2][k];
                for (int i = 0; i < 4; ++i)
                    p.mWeights.As<double>()[i] = 1.0;
                box.mPatches.push_back(p);
                mio::BrepFace face;
                face.mId = f + 1;
                face.mPatch = f;
                box.mFaces.push_back(face);
            }
            return box;
        };
        const auto cad_plate = [&]() {
            constexpr double kTwoPi = 6.283185307179586;
            // A full circle in the plane's parameters: four rational quarters
            // from angle A0, run counter-clockwise or not.
            const auto circle = [&](double Cx, double Cy, double R, double A0, bool Ccw) {
                mio::BrepTrim t;
                t.mDirection = Ccw;
                mio::BrepCurve& c = t.mCurve;
                c.mDegree = 2;
                c.mRational = true;
                const double q = 0.25 * kTwoPi;
                c.mKnots = {0, 0, 0, q, q, 2 * q, 2 * q, 3 * q, 3 * q, 4 * q, 4 * q, 4 * q};
                c.mControlPoints = NDArray(mio::DType::Float64, {9, 3});
                c.mWeights = NDArray(mio::DType::Float64, {9});
                const double w = std::sqrt(0.5);
                for (int k = 0; k < 9; ++k) {
                    const double a = A0 + 0.5 * q * k;
                    const double r = k % 2 ? R / w : R;
                    double* p = c.mControlPoints.As<double>() + 3 * k;
                    p[0] = Cx + r * std::cos(a);
                    p[1] = Cy + r * std::sin(a);
                    p[2] = 0.0;
                    c.mWeights.As<double>()[k] = k % 2 ? w : 1.0;
                }
                return t;
            };
            mio::BrepInfo plate;
            mio::NurbsPatch p;
            p.mDegrees = {1, 1};
            p.mKnots = {{-2.5, -2.5, 2.5, 2.5}, {-2.5, -2.5, 2.5, 2.5}};
            p.mControlPoints = NDArray(mio::DType::Float64, {4, 3});
            p.mWeights = NDArray(mio::DType::Float64, {4});
            for (int j = 0; j < 2; ++j)
                for (int i = 0; i < 2; ++i) {
                    double* x = p.mControlPoints.As<double>() + 3 * (2 * j + i);
                    x[0] = i ? 2.5 : -2.5;
                    x[1] = j ? 2.5 : -2.5;
                    x[2] = 0.0;
                    p.mWeights.As<double>()[2 * j + i] = 1.0;
                }
            plate.mPatches.push_back(p);
            mio::BrepFace face;
            face.mId = 1;
            face.mPatch = 0;
            face.mTrimmed = true;
            face.mRational = true;
            face.mLoops.resize(2);
            face.mLoops[0].mTrims.push_back(circle(0.0, 0.0, 2.0, 0.625 * kTwoPi, true));
            face.mLoops[1].mType = "inner";
            face.mLoops[1].mTrims.push_back(circle(0.98, 0.0, 1.0, 0.0, false));
            mio::BrepTrim line;
            line.mCurve.mDegree = 1;
            line.mCurve.mKnots = {0, 0, 1, 1};
            line.mCurve.mControlPoints = NDArray(mio::DType::Float64, {2, 3});
            line.mCurve.mWeights = NDArray(mio::DType::Float64, {2});
            const double ends[6] = {-1.6, -0.6, 0.0, -0.4, 1.2, 0.0};
            std::copy(ends, ends + 6, line.mCurve.mControlPoints.As<double>());
            line.mCurve.mWeights.As<double>()[0] = 1.0;
            line.mCurve.mWeights.As<double>()[1] = 1.0;
            face.mEmbeddedEdges.push_back(line);
            plate.mFaces.push_back(face);
            return plate;
        };
        // mesh_cad_surface (v17.13.0): the unit box as six untrimmed bilinear
        // patches with no edge topology (the trims are paired by their
        // images), meshed at one lattice step.
        row("mesh_cad_surface", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            of(pD, mio::mesh_cad_surface(box, o).mMesh);
        });
        // mesh_cad_surface on trimmed geometry (v17.14.0): a disc plate of
        // radius 2 with a hole of radius 1 whose edge comes within 0.02 of the
        // rim, and an embedded line across the plate. At h = 1.1 the coarse
        // loops cross and are split until they do not; at the tier's step
        // the plate is filled around the embedded edge.
        row("mesh_cad_surface_trimmed", [&](MeshDigest* pD) {
            const mio::BrepInfo plate = cad_plate();
            for (const double h : {1.1, 4.0 * octree_step}) {
                mio::MeshCadSurfaceOptions o;
                o.mTargetSize = h;
                of(pD, mio::mesh_cad_surface(plate, o).mMesh);
            }
        });
        // Delaunay refinement (v17.15.0): the same two models refined to 25
        // degrees on the surface (the plate also smoothed twice), against
        // the two rows above at the same sizes for refinement's cost.
        row("mesh_cad_surface_refined", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            o.mMinAngleDeg = 25.0;
            of(pD, mio::mesh_cad_surface(box, o).mMesh);
            for (const double h : {1.1, 4.0 * octree_step}) {
                o.mTargetSize = h;
                o.mSmoothIterations = 2;
                of(pD, mio::mesh_cad_surface(plate, o).mMesh);
            }
        });
        // Curved elements (v17.16.0): the same two models at orders 2 and 3,
        // their nodes placed on the geometry and every element's Jacobian
        // bounded, against the unrefined rows above at the same sizes.
        row("mesh_cad_surface_curved", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            const mio::BrepInfo plate = cad_plate();
            for (const std::int32_t order : {2, 3}) {
                mio::MeshCadSurfaceOptions o;
                o.mOrder = order;
                o.mTargetSize = octree_step;
                of(pD, mio::mesh_cad_surface(box, o).mMesh);
                for (const double h : {1.1, 4.0 * octree_step}) {
                    o.mTargetSize = h;
                    of(pD, mio::mesh_cad_surface(plate, o).mMesh);
                }
            }
        });
        // A size field (v17.18.0): the unit box squashed to a slab two
        // lattice steps thick, with two cells across its walls, so the first
        // uniform pass, the gap rays, the graded field and the sized edges
        // and faces all run, against mesh_cad_surface at the same step.
        row("mesh_cad_surface_sized", [&](MeshDigest* pD) {
            mio::BrepInfo slab = cad_box();
            for (mio::NurbsPatch& r_p : slab.mPatches)
                for (std::size_t i = 0; i < 4; ++i)
                    r_p.mControlPoints.As<double>()[3 * i + 2] *= 2.0 * octree_step;
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            o.mCellsAcrossGaps = 2.0;
            o.mGrading = 1.3;
            of(pD, mio::mesh_cad_surface(slab, o).mMesh);
        });
        // Given sizes (v17.19.0): the unit box at four lattice steps with one
        // face given one step, graded into the rest.
        row("mesh_cad_surface_given_sizes", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = 4.0 * octree_step;
            mio::MeshCadSizeSources src;
            src.mFaceSizes = {{1, octree_step}};
            of(pD, mio::mesh_cad_surface(box, o, src).mMesh);
        });
        // A metric source (roadmap 6.10.2): the unit box at four lattice
        // steps, a constant Hessian asking for one step along x alone, so the
        // triangles stretch along y and z.
        row("mesh_cad_surface_metric", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = 4.0 * octree_step;
            o.mEdges = false;
            mio::Mesh bg = mio::mesh_cad_surface(box, o).mMesh;
            NDArray hess(mio::DType::Float64, {bg.NumPoints(), 9});
            for (std::size_t i = 0; i < bg.NumPoints(); ++i)
                for (std::size_t k = 0; k < 9; ++k)
                    hess.As<double>()[9 * i + k] = k == 0 ? 8.0 / (octree_step * octree_step) : 0.0;
            bg.AddPointData("hessian", std::move(hess));
            mio::MeshCadSizeSources src;
            src.pMetricBackground = &bg;
            src.mMetricArray = "hessian";
            src.mMetricError = 1.0;
            of(pD, mio::mesh_cad_surface(box, o, src).mMesh);
        });
        // Bilinear patches over z = 0, each from its four corners (u fastest).
        const auto bilinear = [&](const std::vector<std::array<double, 12>>& rCorners) {
            mio::BrepInfo b;
            for (std::size_t f = 0; f < rCorners.size(); ++f) {
                mio::NurbsPatch p;
                p.mDegrees = {1, 1};
                p.mKnots = {{0, 0, 1, 1}, {0, 0, 1, 1}};
                p.mControlPoints = NDArray(mio::DType::Float64, {4, 3});
                p.mWeights = NDArray(mio::DType::Float64, {4});
                for (int i = 0; i < 12; ++i)
                    p.mControlPoints.As<double>()[i] = rCorners[f][static_cast<std::size_t>(i)];
                for (int i = 0; i < 4; ++i)
                    p.mWeights.As<double>()[i] = 1.0;
                b.mPatches.push_back(p);
                mio::BrepFace face;
                face.mId = static_cast<std::int64_t>(f) + 1;
                face.mPatch = static_cast<std::int64_t>(f);
                b.mFaces.push_back(face);
            }
            return b;
        };
        // The discrete reparametrisation (v17.20.0): a patch whose two corner
        // control points coincide, flattened by its conformal map, refined to
        // 25 degrees.
        row("mesh_cad_surface_reparam", [&](MeshDigest* pD) {
            const mio::BrepInfo patch = bilinear({{0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 2, 0}});
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = 2.0 * octree_step;
            o.mMinAngleDeg = 25.0;
            o.mReparamConditionNumber = 4.0;
            of(pD, mio::mesh_cad_surface(patch, o).mMesh);
        });
        // Virtual faces (v17.20.0): three plates in a row, the middle one a
        // twentieth of the step wide, merged into a neighbour.
        row("mesh_cad_surface_virtual_faces", [&](MeshDigest* pD) {
            const double w = 0.05 * octree_step;
            const mio::BrepInfo strip =
                bilinear({{0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0},
                          {0, 1, 0, 1, 1, 0, 0, 1 + w, 0, 1, 1 + w, 0},
                          {0, 1 + w, 0, 1, 1 + w, 0, 0, 2 + w, 0, 1, 2 + w, 0}});
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            o.mMinAngleDeg = 25.0;
            mio::MeshCadSurfaceExtras extras;
            extras.mVirtualFaces = true;
            of(pD, mio::mesh_cad_surface(strip, o, mio::MeshCadSizeSources{}, extras).mMesh);
        });
        // Corner lopping (v17.20.0): a parallelogram with two 30-degree
        // corners, lopped and refined to 20 degrees.
        row("mesh_cad_surface_lop_corners", [&](MeshDigest* pD) {
            const double c = std::sqrt(3.0), sh = 1.0;
            const mio::BrepInfo para = bilinear({{0, 0, 0, 2, 0, 0, c, sh, 0, 2 + c, sh, 0}});
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            o.mMinAngleDeg = 20.0;
            mio::MeshCadSurfaceExtras extras;
            extras.mLopCorners = true;
            of(pD, mio::mesh_cad_surface(para, o, mio::MeshCadSizeSources{}, extras).mMesh);
        });
        // Untangling (roadmap 6.6.3): the plate's curved rows again, its
        // folded elements' own nodes moved on the surface before any is
        // straightened.
        row("mesh_cad_surface_untangle", [&](MeshDigest* pD) {
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceExtras extras;
            extras.mUntangle = true;
            for (const std::int32_t order : {2, 3})
                for (const double h : {1.1, 4.0 * octree_step}) {
                    mio::MeshCadSurfaceOptions o;
                    o.mOrder = order;
                    o.mTargetSize = h;
                    of(pD,
                       mio::mesh_cad_surface(plate, o, mio::MeshCadSizeSources{}, extras).mMesh);
                }
        });
        // Quadrilaterals (roadmap 6.8): the box and the plate refined to 25
        // degrees, recombined quad-dominant, and the box all-quad.
        row("mesh_cad_surface_quads", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceOptions o;
            o.mMinAngleDeg = 25.0;
            mio::MeshCadSurfaceExtras extras;
            extras.mQuads = mio::MeshCadQuads::Dominant;
            o.mTargetSize = octree_step;
            of(pD, mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, extras).mMesh);
            o.mTargetSize = 4.0 * octree_step;
            of(pD, mio::mesh_cad_surface(plate, o, mio::MeshCadSizeSources{}, extras).mMesh);
            extras.mQuads = mio::MeshCadQuads::All;
            o.mTargetSize = octree_step;
            of(pD, mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, extras).mMesh);
        });
        // The same quad-dominant meshes paired greedily (roadmap 6.8.1: the
        // matching is the default).
        row("mesh_cad_surface_quads_greedy", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceOptions o;
            o.mMinAngleDeg = 25.0;
            mio::MeshCadSurfaceExtras extras;
            extras.mQuads = mio::MeshCadQuads::Dominant;
            extras.mQuadPairs = mio::MeshCadQuadPairs::Greedy;
            o.mTargetSize = octree_step;
            of(pD, mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, extras).mMesh);
            o.mTargetSize = 4.0 * octree_step;
            of(pD, mio::mesh_cad_surface(plate, o, mio::MeshCadSizeSources{}, extras).mMesh);
        });
        // Curved all-quad meshes (roadmap 6.8.3): the box and the plate as
        // quad9, and the plate sized by its gaps at twice the size.
        row("mesh_cad_surface_quad9", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceOptions o;
            o.mOrder = 2;
            mio::MeshCadSurfaceExtras extras;
            extras.mQuads = mio::MeshCadQuads::All;
            o.mTargetSize = octree_step;
            of(pD, mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, extras).mMesh);
            o.mTargetSize = 4.0 * octree_step;
            of(pD, mio::mesh_cad_surface(plate, o, mio::MeshCadSizeSources{}, extras).mMesh);
            o.mOrder = 1;
            o.mCellsAcrossGaps = 2.0;
            of(pD, mio::mesh_cad_surface(plate, o, mio::MeshCadSizeSources{}, extras).mMesh);
        });
        // Back onto the CAD (roadmap 6.10.1): the plate's mesh lifted off its
        // plane and moved in it, projected back onto its face and edges.
        row("project_to_cad", [&](MeshDigest* pD) {
            const mio::BrepInfo plate = cad_plate();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            const mio::Mesh m = mio::mesh_cad_surface(plate, o).mMesh;
            mio::Mesh off = mio::detail::clone_mesh(m);
            mio::NDArray p =
                mio::NDArray::Uninit(mio::DType::Float64, {m.NumPoints(), std::size_t{3}});
            const double* src = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i) {
                p.As<double>()[3 * i] = src[3 * i] + 0.01 * std::sin(double(i));
                p.As<double>()[3 * i + 1] = src[3 * i + 1] + 0.01 * std::cos(double(i));
                p.As<double>()[3 * i + 2] = src[3 * i + 2] + 0.05 * std::sin(3.0 * double(i));
            }
            off.AssignPoints(std::move(p));
            of(pD, mio::project_to_cad(off, plate).mMesh);
        });
        // Boundary layers (roadmap 6.9): four layers grown out of the box's
        // surface and into its all-quad one.
        row("boundary_layer", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            mio::BoundaryLayerOptions b;
            b.mFirstHeight = 0.05 * octree_step;
            b.mNumLayers = 4;
            of(pD, mio::boundary_layer(mio::mesh_cad_surface(box, o).mMesh, b).mMesh);
            mio::MeshCadSurfaceExtras ex;
            ex.mQuads = mio::MeshCadQuads::All;
            b.mReverse = true;
            of(pD, mio::boundary_layer(
                       mio::mesh_cad_surface(box, o, mio::MeshCadSizeSources{}, ex).mMesh, b)
                       .mMesh);
        });
        // Tetrahedra filling the box's skin (roadmap 6.7), only where the
        // Netgen backend is built in; its calls are serialised, so the digest
        // is the same at every thread count.
        if (mio::tetrahedralize_has_netgen())
            row("tetrahedralize", [&](MeshDigest* pD) {
                const mio::BrepInfo box = cad_box();
                mio::MeshCadSurfaceOptions o;
                o.mTargetSize = octree_step;
                mio::TetrahedralizeOptions t;
                t.mBackend = mio::TetrahedralizeBackend::Netgen;
                t.mMaxSize = octree_step;
                of(pD, mio::tetrahedralize(mio::mesh_cad_surface(box, o).mMesh, t).mMesh);
            });
        // The same skin through the native octree backend, always built.
        row("tetrahedralize_octree", [&](MeshDigest* pD) {
            const mio::BrepInfo box = cad_box();
            mio::MeshCadSurfaceOptions o;
            o.mTargetSize = octree_step;
            mio::TetrahedralizeOptions t;
            t.mBackend = mio::TetrahedralizeBackend::Octree;
            of(pD, mio::tetrahedralize(mio::mesh_cad_surface(box, o).mMesh, t).mMesh);
        });
        row("sample_distance", [&](MeshDigest* pD) {
            const NDArray d = mio::sample_distance(surface, jittered->Points());
            if (pD)
                pD->Array(d);
        });
        row("distance_to_surface",
            [&](MeshDigest* pD) { of(pD, mio::distance_to_surface(jittered, surface).mMesh); });
        // The accelerator alone (roadmap §3.3.1.2): the bucket grid and the
        // normal tables of a surface's triangle soup, per call. The soup is
        // made once, by the hashed warm-up call, so it is outside the timed
        // runs. In hash mode the digest also walks the grid -- every bucket of
        // the occupied key box, in the traversal order queries use -- so a
        // change of its layout cannot hide behind the tables.
        std::optional<mio::detail::TriangleSoup> distance_soup;
        row("distance_build", [&](MeshDigest* pD) {
            if (!distance_soup)
                distance_soup = mio::detail::build_triangle_soup(surface, "");
            const mio::detail::DistanceQuery q =
                mio::detail::build_distance_query(*distance_soup, mio::SurfaceDistanceOptions{});
            if (!pD)
                return;
            pD->Bytes(&q.mCellSize, sizeof q.mCellSize);
            pD->Bytes(q.mFaceNormal.data(), q.mFaceNormal.size() * sizeof(mio::detail::Vec3));
            pD->Bytes(q.mVertexNormal.data(), q.mVertexNormal.size() * sizeof(mio::detail::Vec3));
            pD->Bytes(q.mEdgeOfCorner.data(), q.mEdgeOfCorner.size() * sizeof(std::int64_t));
            pD->Bytes(q.mEdgeNormals.data(), q.mEdgeNormals.size() * sizeof(mio::detail::Vec3));
            const mio::detail::GridKey lo = q.mGrid.OccupiedLo();
            const mio::detail::GridKey hi = q.mGrid.OccupiedHi();
            pD->Bytes(&lo, sizeof lo);
            pD->Bytes(&hi, sizeof hi);
            q.mGrid.ForEachInBox(lo, hi, [&](const auto& rIds) {
                pD->U64(rIds.size());
                for (const std::int64_t t : rIds)
                    pD->U64(static_cast<std::uint64_t>(t));
            });
        });
        row("isosurface", [&](MeshDigest* pD) {
            mio::IsosurfaceOptions o;
            o.mArrayName = "u";
            o.mIsovalues = {0.3, 0.6};
            of(pD, mio::isosurface(with_field, o));
        });
        row("slice", [&](MeshDigest* pD) {
            mio::SliceOptions o;
            o.mOrigin = {0.5, 0.45, 0.47};
            o.mNormal = {0.3, 0.2, 1.0};
            of(pD, mio::slice(with_field, o));
        });
        row("linearize", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Linearize;
            of(pD, mio::convert_cells(quadratic, o).mMesh);
        });
        row("interpolate", [&](MeshDigest* pD) {
            mio::InterpolateOptions o;
            o.mMethod = mio::InterpolateMethod::Barycentric;
            o.mExtrapolate = true;
            of(pD, mio::interpolate(with_field, jittered, o));
        });
        row("repair", [&](MeshDigest* pD) { of(pD, mio::repair(surface).mMesh); });
        row("sobolev_deform", [&](MeshDigest* pD) {
            Mesh m = bench_ops_moved(surface, [](std::size_t, double*) {});
            NDArray disp = NDArray::Uninit(DType::Float64, {m.NumPoints(), 3});
            for (std::size_t i = 0; i < m.NumPoints() * 3; ++i)
                disp.As<double>()[i] = 0.01 * static_cast<double>((i * 7) % 13) - 0.06;
            m.AddPointData("d", std::move(disp));
            mio::SobolevOptions o;
            o.mArrayName = "d";
            o.mLengthScale = 2.0 / static_cast<double>(n);
            of(pD, mio::sobolev_deform(m, o).mMesh);
        });
        row("gradient", [&](MeshDigest* pD) {
            mio::GradientOptions o;
            o.mArrayName = "u";
            o.mMethod = mio::GradientMethod::LeastSquares;
            o.mLocation = mio::DataLocation::Point;
            of(pD, mio::gradient(with_field, o).mMesh);
        });
        // Green-Gauss reads each cell's corner connectivity and field values
        // through grad_load_cell, which the least-squares row above does not.
        row("gradient_gg", [&](MeshDigest* pD) {
            mio::GradientOptions o;
            o.mArrayName = "u";
            o.mMethod = mio::GradientMethod::GreenGauss;
            o.mLocation = mio::DataLocation::Point;
            of(pD, mio::gradient(with_field, o).mMesh);
        });
        row("gradient_narrow", [&](MeshDigest* pD) {
            mio::GradientOptions o;
            o.mArrayName = "u";
            o.mMethod = mio::GradientMethod::GreenGauss;
            o.mLocation = mio::DataLocation::Point;
            of(pD, mio::gradient(with_field_narrow, o).mMesh);
        });
        // Green-Gauss on the triangulated surface: the 2-D ring path.
        row("gradient_surface", [&](MeshDigest* pD) {
            Mesh m = bench_ops_moved(surface, [](std::size_t, double*) {});
            NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
            const double* p = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i)
                u.As<double>()[i] = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 2];
            m.AddPointData("u", std::move(u));
            mio::GradientOptions o;
            o.mArrayName = "u";
            o.mMethod = mio::GradientMethod::GreenGauss;
            o.mLocation = mio::DataLocation::Point;
            of(pD, mio::gradient(m, o).mMesh);
        });
        // The ordered diff of the cube against its jittered twin: every cell
        // block is compared through diff_cell_nodes.
        const auto diff_digest = [](MeshDigest* pD, const mio::DiffReport& rR) {
            if (!pD)
                return;
            const auto f64 = [&](double v) {
                std::uint64_t bits;
                std::memcpy(&bits, &v, sizeof bits);
                pD->U64(bits);
            };
            pD->U64(static_cast<std::uint64_t>(rR.mVerdict));
            f64(rR.mPoints.mMaxAbsError);
            pD->U64(static_cast<std::uint64_t>(rR.mPoints.mNumExceeding));
            for (const mio::BlockDiff& b : rR.mBlocks)
                pD->U64(static_cast<std::uint64_t>(b.mConnMismatchCount));
        };
        row("diff", [&](MeshDigest* pD) { diff_digest(pD, mio::diff(volume, jittered)); });
        row("diff_narrow",
            [&](MeshDigest* pD) { diff_digest(pD, mio::diff(volume_narrow, jittered_narrow)); });
        row("refine_narrow", [&](MeshDigest* pD) {
            auto r = mio::refine(volume_narrow);
            of(pD, r.mMesh);
            if (pD)
                pD->Array(r.mPointMap);
        });
        row("linearize_narrow", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Linearize;
            of(pD, mio::convert_cells(quadratic_narrow, o).mMesh);
        });
        row("elevate_narrow", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Elevate;
            auto r = mio::convert_cells(volume_narrow, o);
            of(pD, r.mMesh);
            if (pD) {
                pD->Array(r.mPointMap);
                pD->Arrays(r.mCellMaps);
            }
        });
    }
    return 0;
}
