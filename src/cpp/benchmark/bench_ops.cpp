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
 * @brief Timings of the mesh operations over a size sweep (roadmap §3/§4).
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
 * counts (roadmap §4, "The determinism check the section assumes"). The warmup
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
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/partition.hpp"
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
#include "meshioplusplus/operations/surface.hpp"
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
        // A green-closed refinement of every 7th cell, for `undo_green`.
        const Mesh green = [&] {
            mio::RefineOptions o;
            for (std::size_t c = 0; c < ncells; c += 7)
                o.mCells.push_back(static_cast<std::int64_t>(c));
            o.mRecordHierarchy = true;
            o.mRecordLevels = true;
            return mio::refine(volume, o).mMesh;
        }();
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
        row("shrinkwrap",
            [&](MeshDigest* pD) { of(pD, mio::shrinkwrap(inflated, surface).mMesh); });
        row("voxelize", [&](MeshDigest* pD) {
            mio::VoxelOptions o;
            o.mResolution = std::array<std::int64_t, 3>{48, 48, 48};
            of(pD, mio::voxelize(volume, o).mMesh);
        });
        // Rows for roadmap §4's remaining operation items (v16.18.0).
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
        row("convert_cells", [&](MeshDigest* pD) {
            mio::ConvertCellsOptions o;
            o.mMode = mio::ConvertCellsMode::Elevate;
            of(pD, mio::convert_cells(volume, o).mMesh);
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
    }
    return 0;
}
