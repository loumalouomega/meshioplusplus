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
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/normals.hpp"
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
    ring.reserve(nc * 3);
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
        const Mesh jittered = bench_ops_moved(volume, [&](std::size_t i, double* pP) {
            for (int d = 0; d < 3; ++d) {
                if (pP[d] <= 0.0 || pP[d] >= 1.0)
                    continue;
                pP[d] += 0.3 * h * std::sin(1.7 * static_cast<double>(3 * i + d) + 0.3);
            }
        });
        const Mesh inflated = bench_ops_moved(surface, [](std::size_t, double* pP) {
            for (int d = 0; d < 3; ++d)
                pP[d] = 0.5 + 1.25 * (pP[d] - 0.5);
        });
        // Two copies of the cube, unwelded: the welding rows' input.
        const Mesh doubled = [&] {
            mio::MergeOptions o;
            o.source_tag = false;
            return mio::merge({&volume, &volume}, o).mMesh;
        }();
        // The volume with a smooth point field, for `hessian`.
        const Mesh with_field = [&] {
            Mesh m = bench_ops_moved(volume, [](std::size_t, double*) {});
            NDArray u = NDArray::Uninit(DType::Float64, {m.NumPoints()});
            const double* p = m.Points().As<double>();
            for (std::size_t i = 0; i < m.NumPoints(); ++i)
                u.As<double>()[i] = p[3 * i] * p[3 * i] + p[3 * i + 1] * p[3 * i + 2];
            m.AddPointData("u", std::move(u));
            return m;
        }();
        // The cube elevated to quadratic tetrahedra, for `linearize`.
        const Mesh quadratic = [&] {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Elevate;
            return mio::convert_cells(volume, o).mMesh;
        }();
        // A green-closed refinement of every 7th cell, for `undo_green`.
        const Mesh green = [&] {
            mio::RefineOptions o;
            for (std::size_t c = 0; c < ncells; c += 7)
                o.mCells.push_back(static_cast<std::int64_t>(c));
            o.mRecordHierarchy = true;
            o.mRecordLevels = true;
            return mio::refine(volume, o).mMesh;
        }();
        // The same meshes with float32 points/data and int32 connectivity.
        const Mesh volume_narrow = bench_ops_narrowed(volume);
        const Mesh with_field_narrow = bench_ops_narrowed(with_field);
        const Mesh jittered_narrow = bench_ops_narrowed(jittered);
        const Mesh quadratic_narrow = bench_ops_narrowed(quadratic);
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
        const Mesh ragged = bench_ops_ragged(volume);
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
        const Mesh ragged_field = [&] {
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
        }();
        row("quality", [&](MeshDigest* pD) { of(pD, mio::attach_quality(volume)); });
        row("quality_ragged", [&](MeshDigest* pD) { of(pD, mio::attach_quality(ragged)); });
        row("data_average", [&](MeshDigest* pD) { of(pD, mio::point_data_to_cell_data(with_field)); });
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
            const std::filesystem::path dir =
                std::filesystem::temp_directory_path() /
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
        row("feature_edges",
            [&](MeshDigest* pD) { of(pD, mio::feature_edges(volume).mMesh); });
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
        row("sample_distance", [&](MeshDigest* pD) {
            const NDArray d = mio::sample_distance(surface, jittered.Points());
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
